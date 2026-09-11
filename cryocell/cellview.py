"""
Interactive cell view.

Zoom and pan, hover any compartment for its live state and its Human Protein
Atlas annotation, toggle compartments on and off, and drag the cell or its
organelles to perturb the mechanics directly. Ice closes in from every
direction, not from top and bottom.

Level of detail is tied to zoom: past ~3.5x the plasma membrane resolves into
two leaflets with CPA partitioned into the headgroup region, pores opening
where the MD regime says they should, and the actin cortex underneath.
"""
from __future__ import annotations
import math
import numpy as np
from PyQt6.QtCore import Qt, QPointF, QRectF, pyqtSignal, QTimer
from PyQt6.QtGui import (QPainter, QPen, QBrush, QColor, QPainterPath,
                         QRadialGradient, QFont, QPolygonF, QPixmap)
from PyQt6.QtWidgets import QWidget

from .softbody import SoftBody
from .icefield import IceField
from .hpa import HPA, GEOMETRY
from .model import dT_from_osm, osm_from_dT   # non-ideal liquidus for the phase-diagram inset

# isotonic cell radius (um) of the calibration cell (Viso = 1800 um^3); the view
# draws at a fixed px/um referenced to this, so different cell sizes are visible.
R_ISO_REF = (3 * 1800.0 / (4 * math.pi)) ** (1 / 3)

# Characteristic organelle dimensions in MICRONS (mid-range of published
# mammalian values), so every organelle is drawn at its true physical size
# relative to the whole cell via the fixed px/um scale — a mitochondrion is
# ~1.5 um long in a ~15 um cell, not an eyeballed pixel blob. Mitochondrion is
# (half-length, half-width). In cultured mesenchymal / most somatic cells the
# network is largely punctate-to-short-tubular, ~1-1.5 um long x ~0.4 um wide,
# smaller than the classic 2 um textbook figure. Refs: Cell Biology by the
# Numbers (Milo & Phillips); typical mito diameter 0.25-0.5 um radius.
MITO_UM = (0.72, 0.21)          # ~1.4 x 0.4 um punctate/short-tubular mitochondrion
ORG_UM = {"lysosomes": 0.25,    # ~0.5 um diameter
          "peroxisomes": 0.19,  # ~0.4 um
          "endosomes": 0.22,    # ~0.45 um
          "lipid_drop": 0.32,   # ~0.65 um
          "vesicles": 0.11}     # ~0.2 um transport/secretory vesicle

# ---- palette (matches the HTML build) --------------------------------------
C = dict(
    surface0="#f4f4f1", surface1="#fcfcfb", surface2="#eceae4",
    border="#d9d6cd", borderStrong="#b8b4a8",
    text="#0b0b0b", text2="#52514e", muted="#7a7873",
    nucleus="#2a78d6", mito="#eb6834", er="#1baf7a", golgi="#c06fb8",
    lyso="#e34948", perox="#eda100", endo="#8a6bd1", lipid="#d9a441",
    vesicle="#e87ba4", centro="#eda100", mtub="#5b6470", actin="#a8985f",
    ifil="#9aa3ad", cpa="#4a3aa7", ice="#9fc4e8",
    good="#008300", warn="#eda100", serious="#eb6834", crit="#e34948",
)
# active render mode: "illustrative" | "fluor" (confocal fluorescence) |
# "phase" (label-free phase-contrast greyscale). Set by CellView before paint.
_RENDER = {"mode": "illustrative"}
_SURF = ("surface0", "surface1", "surface2", "surface3")
_LINEK = ("text", "text2", "muted", "faint", "border", "borderStrong")

def col(k, a=255):
    c = QColor(C.get(k, "#808080"))
    m = _RENDER["mode"]
    if m == "fluor":
        if k in _SURF:                         # dark field
            c = QColor(7, 8, 12)
        elif k == "text":
            c = QColor(232, 238, 250)
        elif k in _LINEK:
            c = QColor(150, 160, 178)
        else:                                  # channel colours glow brighter
            c = c.lighter(140)
    elif m == "dark":                          # AIDO-style: dark ground, green cytoplasm
        if k == "surface0":
            c = QColor(9, 13, 11)              # near-black background
        elif k == "surface1":
            c = QColor(46, 120, 78)           # cytoplasm (bright green)
        elif k in ("surface2", "surface3"):
            c = QColor(22, 70, 46)            # cytoplasm shadow
        elif k == "text":
            c = QColor(226, 236, 228)
        elif k in _LINEK:
            c = QColor(150, 178, 160)
        elif k == "er":
            c = QColor(120, 232, 182)         # mint tubules, stand out from cytoplasm
        elif k == "actin":
            c = QColor(168, 210, 150)         # free-ribosome speckle on green
        elif k == "nucleus":
            c = QColor(96, 156, 232)
        elif k == "ice":
            c = QColor(150, 190, 225)
        else:
            c = c.lighter(118)                # organelles pop a little on dark
    elif m == "phase":
        lum = int(0.30 * c.red() + 0.59 * c.green() + 0.11 * c.blue())
        if k in _SURF:                         # light-grey ground / cell body
            c = QColor(198, 199, 196) if k == "surface0" else QColor(206, 207, 205)
        elif k == "text":
            c = QColor(34, 34, 36)
        elif k in _LINEK:
            c = QColor(96, 98, 100)
        else:                                  # structures as mid-grey
            g = int(70 + lum * 0.5)
            c = QColor(g, g, g)
    c.setAlpha(a); return c

LAYERS = [
    ("membrane",     "Plasma membrane",       "plasma_mem"),
    ("cortex",       "Actin cortex",          "actin"),
    ("microtubules", "Microtubules",          "microtubules"),
    ("interm_fil",   "Intermediate filaments","interm_fil"),
    ("centrosome",   "Centrosome",            "centrosome"),
    ("nucleus",      "Nucleus",               "nucleoplasm"),
    ("nucleoli",     "Nucleoli",              "nucleoli"),
    ("nuclear_mem",  "Nuclear envelope",      "nuclear_mem"),
    ("mitochondria", "Mitochondria",          "mitochondria"),
    ("er",           "Endoplasmic reticulum", "er"),
    ("golgi",        "Golgi apparatus",       "golgi"),
    ("lysosomes",    "Lysosomes",             "lysosomes"),
    ("peroxisomes",  "Peroxisomes",           "peroxisomes"),
    ("endosomes",    "Endosomes",             "endosomes"),
    ("lipid_drop",   "Lipid droplets",        "lipid_drop"),
    ("vesicles",     "Vesicles",              "vesicles"),
    ("cpa",          "CPA molecules",         None),
    ("ice",          "Extracellular ice",     None),
    ("iif",          "Intracellular ice",     None),
    ("waterflux",    "Water flux",            None),
    ("stress",       "Stress pathways",       None),
]

# ---- cryo cell-stress pathways (BioGPU-style pathway readout) ---------------
# Each pathway's ACTIVITY comes from a state variable the engine actually
# computes; HPA supplies only the compartment it acts in. Pathways the model does
# not simulate are flagged (simulated=False) and never given a faked activity.
#   key, label, HPA compartment key, glow target, severity colour, simulated
STRESS_PATHWAYS = [
    ("oxid", "Oxidative stress (ROS)",              "mitochondria", "mito",     "crit", True),
    ("apop", "Intrinsic apoptosis (MPT→caspase-3)", "mitochondria", "mito",     "crit", True),
    ("mech", "Mechanotransduction / anoikis",       "focal_adh",    "membrane", "warn", True),
    ("prot", "Unfolded-protein / proteostasis",     "er",           "er",       "warn", True),
    ("memb", "Membrane integrity stress",           "plasma_mem",   "membrane", "warn", True),
    ("ca",   "Ca²⁺ / ionic stress",                 "cytosol",      "cytosol",  "warn", True),
    ("cold", "Cold-shock RNA (CIRBP/RBM3)",         "nucleoplasm",  None,       "muted", False),
    ("dna",  "DNA-damage response",                 "nucleoplasm",  None,       "muted", False),
]

# which compartments a cell type actually has (None = all, for nucleated cells).
# A mature RBC is anucleate and organelle-free (cytosol + membrane only); a
# platelet is anucleate but keeps mitochondria, granules, actin and a marginal
# microtubule band.
CELL_COMPARTMENTS = {
    "rbc": {"cytosol", "plasma_mem"},
    "platelet": {"cytosol", "plasma_mem", "mitochondria", "vesicles", "actin", "microtubules"},
    # Avian RBC: nucleated, keeps mitochondria and a cytoskeleton and does limited
    # protein synthesis, but lacks the full secretory apparatus of a somatic cell.
    "avian_rbc": {"cytosol", "plasma_mem", "nucleoplasm", "nuclear_mem", "nucleoli",
                  "mitochondria", "actin", "microtubules"},
}
def cell_has_compartment(cell_type, comp):
    s = CELL_COMPARTMENTS.get(cell_type)
    return True if s is None else comp in s

def pathway_applicable(key, cell_type):
    """A stress pathway only applies if the cell has the machinery for it."""
    if cell_type == "rbc":        # no nucleus, mitochondria or ER
        return key in ("oxid", "memb", "ca")
    if cell_type == "platelet":   # anucleate: no nuclear DNA/RNA responses
        return key not in ("cold", "dna")
    if cell_type == "avian_rbc":  # nucleated + mitochondria: no ER-stress machinery
        return key != "prot"
    return True

def stress_activity(f):
    """Per-frame activity (0–1) of each stress pathway, from real model state."""
    c = lambda v: float(np.clip(v, 0, 1))
    return {
        "oxid": c(f.get("ros", 0.0)),
        "apop": c(0.5 * f.get("mpt", 0.0) + 0.5 * f.get("casp3", 0.0)),
        "mech": c(max(f.get("rock", 0.0), f.get("bleb", 0.0),
                      f.get("pMLC", 0.0) * (0.3 + 0.7 * f.get("actin", 0.0)))),
        "prot": c(1.0 - f.get("prot", 1.0)),
        "memb": c(f.get("pore", 0.0) + 0.3 * f.get("bleb", 0.0)),
        "ca":   c(f.get("caCyt", 0.0) / 0.4),
        "cold": 0.0, "dna": 0.0,
    }


class Organelle:
    __slots__ = ("kind", "a", "r", "rot", "scale", "x", "y", "px", "py", "drag")
    def __init__(self, kind, a, r, rot, scale):
        self.kind, self.a, self.r, self.rot, self.scale = kind, a, r, rot, scale
        self.x = self.y = self.px = self.py = 0.0
        self.drag = np.zeros(2)


class CellView(QWidget):
    hovered  = pyqtSignal(str)        # compartment key under the cursor ("" if none)
    selected = pyqtSignal(str)
    poked    = pyqtSignal(float)      # magnitude of a user-applied deformation

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMouseTracking(True)
        self.setMinimumSize(560, 460)
        self.setAutoFillBackground(True)

        self.sb = SoftBody(n=160, radius=110.0)
        self.nuc = SoftBody(n=64, radius=44.0, seed=11)
        self.ice = IceField(n_grains=24, seed=7)
        self.R0 = 110.0

        self.zoom, self.pan = 1.5, QPointF(0, 0)
        self.visible = {k for k, _, _ in LAYERS}
        self.frame = None            # dict of scalars for the current instant
        self.sci = True              # scientific overlay: scale bar, labels, cryo-stage ice
        self.render_mode = "dark"           # dark (default) | illustrative | fluor | phase
        self.sel = ""
        self.hover = ""
        self._drag_mode = None
        self._drag_last = None
        self._drag_org = None
        self._ext = None
        self._settle = 200

        rs = np.random.RandomState(23)
        # denser organelle population for a realistic, crowded cytoplasm
        DENS = {"mitochondria": 1.9, "vesicles": 1.3, "lysosomes": 1.4,
                "peroxisomes": 1.3, "endosomes": 1.3, "lipid_drop": 1.2}
        self._ribo_pix = {}          # cached cytoplasm-texture pixmaps, keyed by mode
        self.org = []
        for kind, (vf, count, _src) in GEOMETRY.items():
            if kind in ("nucleus", "er", "golgi", "centrosome"):
                continue
            for i in range(int(count * DENS.get(kind, 1.0))):
                a = rs.uniform(0, 2 * np.pi)
                r = math.sqrt(rs.uniform(0.14, 0.90))
                self.org.append(Organelle(kind, a, r, rs.uniform(0, np.pi),
                                          rs.uniform(0.70, 1.35)))
        self.golgi_a = 2.35
        self.centro_a = 2.05
        self.cpa_seed = np.stack([np.arange(220) * 2.39996,
                                  np.sqrt((np.arange(220) + 0.5) / 220)], axis=1)
        # extracellular solute field (CPA + salt) placed in the unfrozen channel,
        # so freeze-concentration during cooling and washout during dilution show
        self.ecpa_seed = np.stack([rs.uniform(0, 2 * np.pi, 700),
                                   rs.uniform(0, 1, 700),
                                   rs.uniform(0, 1, 700)], axis=1)
        self.mt_seed = rs.uniform(-0.25, 0.25, 64)
        self.if_seed = rs.uniform(0, 2 * np.pi, 40)
        # cytoplasmic crowding / ribosome texture: seeded polar field baked once
        # into a cached pixmap (so density is nearly free per frame). ~1e7
        # ribosomes fill a real cell; this reads as the dense granular cytoplasm
        # of the reference rather than an empty interior. (angle, radius, size, shade)
        nrib = 2400
        self.ribo = np.stack([rs.uniform(0, 2 * np.pi, nrib),
                              np.sqrt(rs.uniform(0.02, 0.99, nrib)),
                              rs.uniform(0.5, 1.4, nrib),
                              rs.uniform(0.35, 0.85, nrib)], axis=1)
        # chromatin granularity inside the nucleus
        nchr = 130
        self.chrom_seed = np.stack([rs.uniform(0, 2 * np.pi, nchr),
                                    np.sqrt(rs.uniform(0.0, 0.94, nchr)),
                                    rs.uniform(0.5, 1.6, nchr)], axis=1)

        self.timer = QTimer(self); self.timer.timeout.connect(self._tick); self.timer.start(33)

    # ------------------------------------------------------------ state feed
    def set_frame(self, f: dict, jumped: bool = False):
        self.frame = f
        if jumped: self._settle = 160
        self.update()

    def set_layer(self, key, on):
        (self.visible.add(key) if on else self.visible.discard(key)); self.update()

    def set_sci(self, on):
        self.sci = bool(on); self.update()

    def set_render_mode(self, mode):
        self.render_mode = mode if mode in ("illustrative", "fluor", "phase", "dark") else "illustrative"
        self.update()

    def _ribo_pixmap(self, mode):
        """Cytoplasm-crowding texture, rendered once per mode into a pixmap."""
        pix = self._ribo_pix.get(mode)
        if pix is not None:
            return pix
        R0 = self.R0; D = int(2 * R0) + 4
        pix = QPixmap(D, D); pix.fill(Qt.GlobalColor.transparent)
        p = QPainter(pix); p.setRenderHint(QPainter.RenderHint.Antialiasing)
        base = col("er") if mode == "fluor" else col("actin")
        p.setPen(Qt.PenStyle.NoPen); c0 = D / 2
        for a, rf, sz, sh in self.ribo:
            c = QColor(base); c.setAlpha(int((58 if mode == "fluor" else 55) * sh))
            p.setBrush(c)
            p.drawEllipse(QPointF(c0 + math.cos(a) * R0 * rf, c0 + math.sin(a) * R0 * rf),
                          sz * 0.85, sz * 0.85)
        p.end()
        self._ribo_pix[mode] = pix
        return pix

    # ------------------------------------------------------------- mechanics
    def _tick(self):
        f = self.frame
        if not f: return
        # FIXED pixels-per-micron (referenced to the calibration cell, r_iso 7.55
        # um at Viso 1800), so absolute cell size is visible: a smaller cell line
        # (HeLa, A549) draws smaller than a large mesenchymal cell (hMSC), instead
        # of every cell being renormalised to fill the view.
        px_um = self.R0 / R_ISO_REF
        sizef = f["r_iso_um"] / R_ISO_REF
        R = self.R0 * sizef * max(f["Vn"], 0.05) ** (1 / 3)

        self.ice.update(f["fIce"] if f["frozen"] else 0.0, f["grain"],
                        px_um, f["chanW"], R)

        mem_per = 2 * math.pi * self.R0 * sizef * f["APL"] * (f["Vn"] ** (1/3) if f["Vn"] > 1 else 1)
        L0 = mem_per / self.sb.n
        slack = float(np.clip(mem_per / (2 * math.pi * max(R, 1)) - 1, 0, 1.5))
        kB = float(np.clip(0.008 + 0.014 * f["gel"] + 0.00025 * f["sterol"]
                           - 0.004 * f["pore"], 0.004, 0.030))
        kT = 0.075 * (1 + 0.6 * f["tension"])
        kC = float(np.clip(0.085 * (0.5 + 0.9 * f["pMLC"]) * (1 - 0.55 * f["pore"])
                           * (0.35 + 0.65 * f["actin"]), 0.012, 0.16))

        crystals = []
        if f["Piif"] > 0.02:
            k = np.arange(max(1, int(f["Piif"] * 11)))
            aa = k * 2.39996; rr = np.sqrt((k + 0.5) / 12) * R * 0.7
            crystals = list(zip(np.cos(aa) * rr, np.sin(aa) * rr,
                                np.full(len(k), 3 + f["Piif"] * 9)))

        # Gel phase = 2D solid: bending rigidity jumps ~10x below Tm and thermal
        # fluctuations are quenched (amplitude ~ sqrt(kT/kappa); Dimova 2014;
        # Steltenkamp 2006). Freeze the contour by damping out velocity and
        # killing the jiggle, so a cold membrane visibly holds its (crenated)
        # shape rather than wobbling. Existing folds are preserved because the
        # freeze acts on motion, not on the bending set-point.
        gel = f["gel"]
        # Below Tm the bilayer is in the ordered gel (Lβ') phase — a 2D solid:
        # undulations are quenched (amplitude ~ sqrt(kT/kappa)), lateral lipid
        # remodelling and cortex flow cease, and the contour locks into the
        # (crenated) shape it froze in. Ramp the freeze with the gel fraction, so
        # cooling through Tm visibly stiffens the membrane and a stored cell holds
        # a rigid shape instead of wobbling like a fluid vesicle.
        damp_g  = 0.88 * (1 - 0.95 * gel)              # velocity killed in gel
        noise_g = (0.055 * slack + 0.005) * (1 - 0.98 * gel)   # thermal jiggle quenched
        kT_g = kT * (1 - 0.70 * gel)                   # active tension remodelling stops
        kC_g = kC * (1 - 0.85 * gel)                   # actin cortex frozen out
        self.sb.update_blebs(f["bleb"])
        pocket = self.ice.pocket_radius if (f["frozen"] and self.ice.active) else None
        if self._settle > 0:
            steps = 22; self._settle -= 22
        elif gel > 0.80:
            steps = 0                                  # solid gel: hold shape, no dynamics
        else:
            steps = 7
        for _ in range(steps):
            self.sb.step(L0, kT_g, kB, 0.60, math.pi * R * R, damp_g, slack,
                         k_cortex=kC_g, pocket_r=pocket, crystals=crystals,
                         noise=noise_g, external=self._ext)
            _ct = f.get("cell_type")                                       # N:C ratio by cell type
            nuc_ratio = 0.62 if _ct == "tcell" else 0.55 if _ct == "avian_rbc" else 0.40
            Rn = self.R0 * nuc_ratio * max(f["Vnuc"], 0.1) ** (1 / 3)
            self.nuc.step(2 * math.pi * Rn / self.nuc.n, 0.11, 0.30, 0.42,
                          math.pi * Rn * Rn, 0.88, 0.02, k_cortex=0.09)
            # nucleus stays inside the plasma membrane
            ang = np.arctan2(self.nuc.p[:, 1], self.nuc.p[:, 0])
            lim = self.sb.radius_at(ang) * 0.80
            d = np.hypot(self.nuc.p[:, 0], self.nuc.p[:, 1]); d[d < 1e-6] = 1e-6
            over = d > lim
            if over.any():
                self.nuc.p[over] *= (lim[over] / d[over])[:, None]
                self.nuc.v[over] *= 0.4
        self._ext = None
        self._place_organelles()
        self.update()

    def _place_organelles(self):
        f = self.frame
        # cytoplasm vitrification: a glass has no long-range diffusion, so as the
        # cytosol vitrifies (WLF mobility -> 0 near Tg') the organelles lock in
        # place and stop tracking the membrane. track = 1 fluid .. 0 glass.
        glass = float(np.clip(f.get("glass", 0.0), 0, 1))
        track = 1.0 - min(glass / 0.85, 1.0)
        angs = np.array([o.a for o in self.org])
        rads = self.sb.radius_at(angs)
        nucR = np.sqrt(self.nuc.area() / np.pi)
        for o, ra in zip(self.org, rads):
            rr = ra * o.r * 0.92
            # push organelles out of the nucleus
            if rr < nucR * 1.12: rr = nucR * 1.12 + (ra - nucR * 1.12) * 0.25
            o.px, o.py = o.x, o.y
            tx = math.cos(o.a) * rr + o.drag[0]
            ty = math.sin(o.a) * rr + o.drag[1]
            o.x += track * (tx - o.x)
            o.y += track * (ty - o.y)
            o.drag *= 0.90

    # ------------------------------------------------------------- transform
    def _to_screen(self, x, y):
        return QPointF(self.width() / 2 + self.pan.x() + x * self.zoom,
                       self.height() / 2 + self.pan.y() + y * self.zoom)
    def _to_world(self, pt):
        return ((pt.x() - self.width() / 2 - self.pan.x()) / self.zoom,
                (pt.y() - self.height() / 2 - self.pan.y()) / self.zoom)

    # ----------------------------------------------------------------- input
    def wheelEvent(self, e):
        wx, wy = self._to_world(e.position())
        self.zoom = float(np.clip(self.zoom * (1.0015 ** e.angleDelta().y()), 0.35, 14.0))
        p = self._to_screen(wx, wy)
        self.pan += e.position() - p
        self.update()

    def mousePressEvent(self, e):
        self._drag_last = e.position()
        if e.button() == Qt.MouseButton.MiddleButton or \
           (e.button() == Qt.MouseButton.LeftButton and e.modifiers() & Qt.KeyboardModifier.ShiftModifier):
            self._drag_mode = "pan"
        elif e.button() == Qt.MouseButton.LeftButton:
            wx, wy = self._to_world(e.position())
            hit, org = self._hit(wx, wy)
            self.sel = hit; self.selected.emit(hit)
            self._drag_mode = "org" if org is not None else "poke"
            self._drag_org = org
        self.update()

    def mouseMoveEvent(self, e):
        if self._drag_mode == "pan":
            self.pan += e.position() - self._drag_last
            self._drag_last = e.position(); self.update(); return
        if self._drag_mode in ("poke", "org"):
            d = e.position() - self._drag_last
            self._drag_last = e.position()
            wdx, wdy = d.x() / self.zoom, d.y() / self.zoom
            if self._drag_mode == "org" and self._drag_org is not None:
                self._drag_org.drag += np.array([wdx, wdy]) * 0.9
                self.poked.emit(float(np.hypot(wdx, wdy)))
            else:
                wx, wy = self._to_world(e.position())
                d2 = self.sb.p - np.array([wx, wy])
                dist = np.hypot(d2[:, 0], d2[:, 1])
                w = np.exp(-(dist / 42.0) ** 2)
                ext = np.zeros_like(self.sb.p)
                ext[:, 0] = w * wdx * 0.28; ext[:, 1] = w * wdy * 0.28
                self._ext = ext
                self.poked.emit(float(np.hypot(wdx, wdy)))
            return
        wx, wy = self._to_world(e.position())
        h, _ = self._hit(wx, wy)
        if h != self.hover:
            self.hover = h; self.hovered.emit(h); self.update()

    def mouseReleaseEvent(self, e):
        self._drag_mode = None; self._drag_org = None; self._ext = None

    def _hit(self, wx, wy):
        d = math.hypot(wx, wy)
        best, org, bd = "", None, 1e9
        for o in self.org:
            if o.kind not in self.visible: continue
            dd = math.hypot(wx - o.x, wy - o.y)
            if dd < 13 * o.scale and dd < bd:
                bd, best, org = dd, o.kind, o
        if best: return best, org
        nucR = np.sqrt(self.nuc.area() / np.pi)
        if "nucleus" in self.visible and d < nucR: return "nucleus", None
        r_here = float(self.sb.radius_at(np.array([math.atan2(wy, wx)]))[0])
        if d > r_here * 0.90 and d < r_here * 1.18 and "membrane" in self.visible:
            return "membrane", None
        if d < r_here: return "cytosol", None
        if self.ice.active and "ice" in self.visible: return "ice", None
        return "", None

    # --------------------------------------------------------------- drawing
    def paintEvent(self, _):
        _RENDER["mode"] = self.render_mode          # make col() mode-aware for this paint
        q = QPainter(self)
        q.setRenderHint(QPainter.RenderHint.Antialiasing)
        q.fillRect(self.rect(), col("surface0"))
        f = self.frame
        if not f:
            q.setPen(col("muted")); q.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter,
                                               "Run a protocol")
            return
        Z, S = self.zoom, self._to_screen
        px_um = self.R0 / R_ISO_REF

        # ---- extracellular ice, closing in from every direction
        if self.ice.active and "ice" in self.visible:
            q.save()
            q.setPen(Qt.PenStyle.NoPen)
            ice_base = col("ice")
            shade0 = QColor(ice_base); shade0.setAlpha(int(55 + 120 * f["fIce"]))
            q.fillRect(self.rect(), shade0)
            extent = max(self.width(), self.height()) / max(Z, 1e-3)
            cells = self.ice.grain_cells(extent)
            for poly, sh in cells:
                cc = QColor(ice_base)
                cc.setAlpha(int(50 + 90 * f["fIce"] * (0.45 + 0.55 * sh)))
                q.setBrush(cc); q.setPen(Qt.PenStyle.NoPen)
                q.drawPolygon(QPolygonF([S(x, y) for x, y in poly]))
            q.setPen(QPen(col("surface2", 230), max(1.0, 1.3 * min(Z, 2.0))))
            for v1, v2 in self.ice.grain_boundaries(extent):
                q.drawLine(S(v1[0], v1[1]), S(v2[0], v2[1]))
            # punch out the unfrozen pocket
            pts, _th, _r = self.ice.polygon(240)
            poly = QPolygonF([S(x, y) for x, y in pts])
            q.setBrush(col("surface1")); q.setPen(QPen(col("ice"), 2.0)); q.drawPolygon(poly)
            # cryo-stage look: dendritic ice fingers growing inward from the
            # freezing front into the unfrozen channel (directional solidification)
            if self.sci:
                q.setPen(QPen(col("ice", 210), max(0.8, 1.0 * min(Z, 2.0))))
                npts = len(pts)
                for i in range(0, npts, 5):
                    x0, y0 = pts[i]
                    rr = math.hypot(x0, y0) or 1e-6
                    ix, iy = -x0 / rr, -y0 / rr                    # inward normal
                    tx, ty = -iy, ix                               # tangent
                    spike = (6 + 10 * f["fIce"]) * (0.6 + 0.4 * math.sin(i * 1.7))
                    tipx, tipy = x0 + ix * spike, y0 + iy * spike
                    q.drawLine(S(x0, y0), S(tipx, tipy))           # primary dendrite
                    for s in (-1, 1):                              # side branches
                        bx = x0 + ix * spike * 0.5 + tx * s * spike * 0.28
                        by = y0 + iy * spike * 0.5 + ty * s * spike * 0.28
                        q.drawLine(S(x0 + ix * spike * 0.5, y0 + iy * spike * 0.5), S(bx, by))
            q.restore()
        else:
            # subtle cell shadow — kept dark in AIDO mode so it blends into the ground
            q.save(); q.setPen(Qt.PenStyle.NoPen)
            q.setBrush(col("surface0" if self.render_mode == "dark" else "surface1"))
            r = self.sb.radii().max() * 1.5
            q.drawEllipse(S(0, 0), r * Z, r * Z); q.restore()

        # extracellular CPA / solute in the unfrozen channel (freeze-concentration)
        self._draw_extracellular_cpa(q, f)

        # ---- cell body
        # Avian erythrocytes are ELLIPSOIDAL (an oval cell with an oval nucleus),
        # unlike the round mammalian RBC. Draw the whole cell body through an
        # anisotropic scale about the cell centre (rendering only ~2:1 long axis;
        # area preserved) so the membrane, nucleus and organelles all read oval.
        avian = f.get("cell_type") == "avian_rbc"
        q.save()
        if avian:
            _c0 = S(0.0, 0.0)
            q.translate(_c0.x(), _c0.y()); q.scale(1.42, 0.70); q.translate(-_c0.x(), -_c0.y())
        path = self._path(self.sb.p)
        rmean = math.sqrt(self.sb.area() / math.pi)
        fluor = self.render_mode == "fluor"
        phase = self.render_mode == "phase"
        if not fluor:                                  # dark-field: no cytoplasm fill
            g = QRadialGradient(S(-rmean * .3, -rmean * .35), rmean * 1.25 * Z)
            g.setColorAt(0, col("surface1")); g.setColorAt(1, col("surface2"))
            q.setBrush(QBrush(g)); q.setPen(Qt.PenStyle.NoPen); q.drawPath(path)
        if phase:                                      # phase-contrast bright halo
            q.setPen(QPen(QColor(250, 250, 250, 200), max(2.0, 3.0 * min(Z, 2.2))))
            q.setBrush(Qt.BrushStyle.NoBrush); q.drawPath(self._path(self.sb.p * 1.02))
            q.setPen(QPen(QColor(60, 60, 62, 160), max(1.0, 1.4 * Z)))
            q.drawPath(self._path(self.sb.p * 0.98))

        q.save(); q.setClipPath(path)
        if fluor:                                      # structures glow additively
            q.setCompositionMode(QPainter.CompositionMode.CompositionMode_Plus)
        self._draw_interior(q, f, rmean, px_um)
        # cytoplasm vitrification: a cool glassy sheen over the interior once the
        # cytosol passes into the glass (WLF mobility -> 0). Signals that the
        # whole cell, not just the membrane, is now a rigid solid.
        glass = f.get("glass", 0.0)
        if glass > 0.25:
            gg = QRadialGradient(S(-rmean * .3, -rmean * .35), rmean * 1.3 * Z)
            gg.setColorAt(0, col("ice", int(28 * glass)))
            gg.setColorAt(1, col("nucleus", int(46 * glass)))
            q.setBrush(QBrush(gg)); q.setPen(Qt.PenStyle.NoPen); q.drawPath(path)
        q.restore()

        self._draw_membrane(q, f, path, rmean)
        self._draw_stress_glow(q, f, rmean)
        self._draw_ice_penetration(q, f, rmean)
        self._draw_water_flux(q, f, rmean)
        q.restore()                                    # end avian elliptical transform
        self._draw_overlay(q, f, px_um, rmean)
        q.end()

    def _path(self, pts):
        p = QPainterPath()
        S = self._to_screen; n = len(pts)
        mid = lambda i, j: S((pts[i][0] + pts[j][0]) / 2, (pts[i][1] + pts[j][1]) / 2)
        p.moveTo(mid(n - 1, 0))
        for i in range(n):
            p.quadTo(S(pts[i][0], pts[i][1]), mid(i, (i + 1) % n))
        p.closeSubpath(); return p

    def _draw_extracellular_cpa(self, q, f):
        """Extracellular CPA + salt in the space around the cell. Density tracks
        the extracellular osmolality (Osme): sparse in isotonic medium, crowding
        as the unfrozen channel freeze-concentrates during cooling, and thinning
        again on warming/dilution — so the transmembrane CPA gradient is visible
        against the intracellular CPA drawn inside the cell."""
        if "cpa" not in self.visible:
            return
        S, Z = self._to_screen, self.zoom
        osme = float(f.get("Osme", 1.0))
        dens = float(np.clip((osme - 0.8) / 9.0, 0.03, 1.0))    # 0 isotonic → 1 freeze-conc.
        n = int(dens * len(self.ecpa_seed))
        if n <= 0:
            return
        thetas = self.ecpa_seed[:n, 0]
        memb = self.sb.radius_at(thetas)
        active = self.ice.active and f.get("frozen", False)
        if active:                                              # confined to the unfrozen channel
            outer = np.minimum(self.ice.pocket_radius(thetas), memb * 3.0)
        else:                                                   # surrounding medium
            outer = memb * 1.6
        cpa_c = col("cpa", 165); salt_c = col("ice", 150)
        q.setPen(Qt.PenStyle.NoPen)
        for k in range(n):
            m = float(memb[k]); o = float(outer[k])
            if o <= m + 1.0:
                continue
            a = float(thetas[k])
            rr = m + (o - m) * (0.06 + 0.9 * float(self.ecpa_seed[k, 1]))
            grow = 0.85 + 0.5 * dens                           # heavier dots when concentrated
            if float(self.ecpa_seed[k, 2]) < 0.68:             # CPA molecule
                q.setBrush(cpa_c); rad = 1.7 * Z * grow
            else:                                              # salt ion (rest of the solution)
                q.setBrush(salt_c); rad = 1.05 * Z * grow
            q.drawEllipse(S(math.cos(a) * rr, math.sin(a) * rr), rad, rad)

    def _draw_ice_penetration(self, q, f, rmean):
        """Without CPA, the extracellular ice front breaches the plasma membrane
        and seeds intracellular ice (surface-catalysed nucleation through membrane
        defects; Mazur 1965, Toner 1990). Jagged ice fingers cross the membrane
        where P_iif is high. A cryoprotectant keeps the membrane a barrier and
        depresses nucleation, so the fingers vanish (P_iif → 0)."""
        if not ("iif" in self.visible and self.ice.active and f.get("frozen", False)):
            return
        pf = float(np.clip(f.get("Piif", 0.0), 0, 1))
        if pf < 0.25:                                   # protected: no penetration
            return
        S, Z = self._to_screen, self.zoom
        ic = col("ice", int(200 + 55 * pf))
        n = int(3 + pf * 7)
        for j in range(n):
            a = j / n * 2 * math.pi + 0.4
            memb = float(self.sb.radius_at(np.array([a]))[0])
            front = float(self.ice.pocket_radius(np.array([a]))[0])
            r_out = max(front, memb * 1.05)             # start at the ice front, outside
            r_in = memb * (0.45 - 0.20 * pf)            # penetrate into the cytoplasm
            ca, sa = math.cos(a), math.sin(a); ta, tb = -sa, ca
            q.setPen(QPen(ic, max(1.2, 2.0 * Z * pf))); q.setBrush(Qt.BrushStyle.NoBrush)
            p = QPainterPath(S(ca * r_out, sa * r_out)); steps = 7
            for s in range(1, steps + 1):
                t = s / steps; rr = r_out + (r_in - r_out) * t
                off = math.sin(t * math.pi * 3.0) * memb * 0.05 * (1 - t)
                p.lineTo(S(ca * rr + ta * off, sa * rr + tb * off))
            q.drawPath(p)
            tipr = r_in + (r_out - r_in) * 0.22         # side dendrites near the tip
            for sgn in (-1, 1):
                q.drawLine(S(ca * tipr, sa * tipr),
                           S(ca * tipr + ta * sgn * memb * 0.08, sa * tipr + tb * sgn * memb * 0.08))
            q.setBrush(col("crit", 170)); q.setPen(Qt.PenStyle.NoPen)  # membrane breach
            q.drawEllipse(S(ca * memb, sa * memb), 2.6 * Z, 2.6 * Z)

    def _draw_stress_glow(self, q, f, rmean):
        """BioGPU-style: cryo cell-stress pathways light up on the compartments
        they act in, intensity from the model's live state (see stress_activity).
        Only simulated pathways glow; cold-shock/DNA-damage are never faked."""
        if "stress" not in self.visible:
            return
        S, Z = self._to_screen, self.zoom
        act = stress_activity(f)
        def glow(pt, r, key, a):
            if a < 0.12:
                return
            g = QRadialGradient(pt, r)
            c0 = QColor(col(key)); c0.setAlpha(int(165 * min(a, 1)))
            c1 = QColor(col(key)); c1.setAlpha(0)
            g.setColorAt(0, c0); g.setColorAt(1, c1)
            q.setBrush(QBrush(g)); q.setPen(Qt.PenStyle.NoPen)
            q.drawEllipse(pt, r, r)
        q.save()
        # mitochondria: oxidative stress + intrinsic apoptosis (red)
        mv = max(act["oxid"], act["apop"])
        if mv > 0.12 and "mitochondria" in self.visible:
            for o in self.org:
                if o.kind == "mitochondria":
                    glow(S(o.x, o.y), 13 * Z, "crit", mv)
        # proteostasis: unfolded-protein stress across the ER/cytosol (amber ring)
        if act["prot"] > 0.12:
            g = QRadialGradient(S(0, 0), rmean * 0.9 * Z)
            c0 = QColor(col("warn")); c0.setAlpha(int(70 * act["prot"]))
            c1 = QColor(col("warn")); c1.setAlpha(0)
            g.setColorAt(0.35, c1); g.setColorAt(0.7, c0); g.setColorAt(1, c1)
            q.setBrush(QBrush(g)); q.setPen(Qt.PenStyle.NoPen)
            q.drawEllipse(S(0, 0), rmean * 0.9 * Z, rmean * 0.9 * Z)
        # membrane: mechanotransduction + integrity stress (amber ring on the cortex)
        mm = max(act["mech"], act["memb"])
        if mm > 0.12:
            q.setBrush(Qt.BrushStyle.NoBrush)
            pen = QPen(QColor(col("warn"))); pen.setWidthF(max(3.0, 8.0 * Z * mm))
            cc = QColor(col("warn")); cc.setAlpha(int(150 * mm)); pen.setColor(cc)
            q.setPen(pen); q.drawPath(self._path(self.sb.p * 0.99))
        q.restore()

    def _draw_water_flux(self, q, f, rmean):
        """Osmotic water flux across the membrane — the engine that drives volume
        change. Efflux (arrows out) as the cell dehydrates during CPA loading and
        freeze-concentration; influx (arrows in) as it rehydrates on dilution.
        Magnitude from the model's volume rate of change (dVw)."""
        if "waterflux" not in self.visible:
            return
        dvw = float(f.get("dVw", 0.0))
        if abs(dvw) < 0.004:
            return
        S, Z = self._to_screen, self.zoom
        mag = float(np.clip(abs(dvw) / 0.06, 0.25, 1.0))
        out = dvw < 0                                  # water leaving the cell
        wcol = QColor(64, 150, 226)                    # water blue
        wcol.setAlpha(int(150 + 90 * mag))
        q.setPen(QPen(wcol, max(1.4, 2.2 * Z * mag)))
        L = rmean * (0.16 + 0.20 * mag)                # arrow length
        for k in range(10):
            a = k / 10 * 2 * math.pi + 0.31
            rm = float(self.sb.radius_at(np.array([a]))[0])
            ca, sa = math.cos(a), math.sin(a)
            if out:                                    # tail just inside → tip outside
                x0, y0 = rm * 0.86, rm * 0.86
                x1, y1 = (rm + L), (rm + L)
                tail = S(ca * x0, sa * y0); tip = S(ca * x1, sa * y1)
            else:                                      # tail outside → tip just inside
                x0, y0 = (rm + L), (rm + L)
                x1, y1 = rm * 0.86, rm * 0.86
                tail = S(ca * x0, sa * y0); tip = S(ca * x1, sa * y1)
            q.drawLine(tail, tip)
            ang = math.atan2(tip.y() - tail.y(), tip.x() - tail.x()); hs = 4 + 4 * mag
            q.setBrush(wcol); q.setPen(Qt.PenStyle.NoPen)
            q.drawPolygon(QPolygonF([tip,
                QPointF(tip.x() - hs * math.cos(ang - 0.5), tip.y() - hs * math.sin(ang - 0.5)),
                QPointF(tip.x() - hs * math.cos(ang + 0.5), tip.y() - hs * math.sin(ang + 0.5))]))
            q.setPen(QPen(wcol, max(1.4, 2.2 * Z * mag)))

    def _draw_phase_diagram(self, q, f):
        """Inset phase diagram: the extracellular solution rides the liquidus as it
        freeze-concentrates. Plots the non-ideal freezing-point-depression curve
        and the current (osmolality, temperature) operating point."""
        w, h = self.width(), self.height()
        bw, bh = 148, 104
        x0, y0 = 12, h - bh - 34            # bottom-left, above the scale bar
        q.setBrush(col("surface1", 235)); q.setPen(QPen(col("border"), 1))
        q.drawRoundedRect(QRectF(x0, y0, bw, bh), 6, 6)
        pl, pt, pr, pb = x0 + 30, y0 + 16, x0 + bw - 8, y0 + bh - 20
        OSM_MAX, T_MIN = 12.0, -22.0                    # axes ranges (osmol/L, °C)
        def PX(osm): return pl + (pr - pl) * min(osm, OSM_MAX) / OSM_MAX
        def PY(T):   return pt + (pb - pt) * min(max(-T, 0), -T_MIN) / (-T_MIN)
        q.setPen(QPen(col("muted"), 1))
        q.drawLine(QPointF(pl, pt), QPointF(pl, pb)); q.drawLine(QPointF(pl, pb), QPointF(pr, pb))
        # liquidus curve T = -dT_from_osm(osm): below it, ice + freeze-concentrated brine
        q.setPen(QPen(col("ice"), 1.8)); path = None
        pts = []
        for j in range(41):
            osm = OSM_MAX * j / 40; T = -dT_from_osm(osm)
            if T < T_MIN: break
            pts.append(QPointF(PX(osm), PY(T)))
        for j in range(1, len(pts)):
            q.drawLine(pts[j - 1], pts[j])
        # extracellular point (on the liquidus — in equilibrium with external ice)
        osm, T = float(f.get("Osme", 1.0)), float(f.get("T", 22.0))
        q.setBrush(col("ice")); q.setPen(QPen(col("ice").darker(140), 1))
        q.drawEllipse(QPointF(PX(osm), PY(T)), 3.5, 3.5)
        # intracellular supercooling — the SCN driver (Toner 1990; Li et al. 2020).
        # The cell lags below its own liquidus by dTsc; when that gap is large the
        # surface-catalysed nucleation rate J ~ exp(-K/(T^3 dTsc^2)) explodes → IIF.
        dTsc = float(f.get("dTsc", 0.0))
        if dTsc > 0.2 and T < 0:
            osm_i = osm_from_dT(max(0.0, -T - dTsc))    # intracellular osmolality
            Tm_i = -dT_from_osm(osm_i)                   # its equilibrium liquidus point
            xi = PX(osm_i)
            hot = dTsc > 4.0                              # SCN firing regime
            gc = col("crit") if hot else col("warn")
            q.setPen(QPen(gc, 1.6))                       # supercooling gap up to the liquidus
            q.drawLine(QPointF(xi, PY(T)), QPointF(xi, PY(Tm_i)))
            q.setBrush(col("cpa")); q.setPen(QPen(col("cpa").darker(140), 1))
            q.drawEllipse(QPointF(xi, PY(T)), 3.5, 3.5)   # supercooled intracellular point
            q.setPen(gc); q.setFont(QFont("", 6, QFont.Weight.Bold))
            q.drawText(QRectF(x0, pt - 2, bw - 6, 10), Qt.AlignmentFlag.AlignRight,
                       f"ΔTsc {dTsc:.0f}°" + (" · SCN→IIF" if hot else ""))
        q.setPen(col("text2")); q.setFont(QFont("", 7))
        q.drawText(QRectF(x0, y0 + 1, bw, 12), Qt.AlignmentFlag.AlignHCenter, "Liquidus · supercooling (SCN)")
        q.setFont(QFont("", 6)); q.setPen(col("muted"))
        q.drawText(QPointF(x0 + 2, pb + 9), "osmolality →")
        q.save(); q.translate(x0 + 9, pb); q.rotate(-90)
        q.drawText(QRectF(0, 0, pb - pt, 10), Qt.AlignmentFlag.AlignLeft, "T (°C)"); q.restore()

    def _draw_interior(self, q, f, rmean, px_um):
        S, Z = self._to_screen, self.zoom
        vis = self.visible
        # resting cytoskeletal density of this cell line: mesenchymal cells
        # (hMSC) carry a denser network than epithelial lines (HeLa, A549).
        cyto = float(np.clip(f.get("cyto", 1.0), 0.3, 2.0))
        # cell-type morphology (rendering only): a mature red blood cell is
        # anucleate and organelle-free, its interior packed with haemoglobin; a
        # T lymphocyte is small with a high nucleus-to-cytoplasm ratio.
        ctype = f.get("cell_type", "msc")
        # A mature MAMMALIAN red cell is anucleate and organelle-free. An AVIAN
        # red cell (Bissoyi et al. 2025, ACS Polym Au 6:366) keeps its nucleus and
        # mitochondria and performs limited protein synthesis, so it is drawn
        # haemoglobin-red BUT with a nucleus and organelles layered on top.
        is_rbc_like = ctype in ("rbc", "avian_rbc")
        has_org = ctype != "rbc"                     # mammalian RBC: none; platelet/avian: yes
        has_nucleus = ctype not in ("rbc", "platelet")  # mammalian RBC + platelet anucleate

        # ---- cytoplasmic interior ----
        if self.render_mode != "phase":
            St = self._to_screen(0, 0); rpx = rmean * Z
            if is_rbc_like:
                # interior packed with haemoglobin (red radial fill); an avian
                # RBC keeps its nucleus + organelles, drawn over this fill below.
                hb = QColor(196, 60, 52)
                g = QRadialGradient(St, rpx)
                g.setColorAt(0.0, QColor(232, 150, 140)); g.setColorAt(0.42, hb)
                g.setColorAt(1.0, hb.darker(118))
                q.setBrush(QBrush(g)); q.setPen(Qt.PenStyle.NoPen)
                q.drawEllipse(St, rpx, rpx)
            else:
                gl = f.get("glass", 0.0)
                pix = self._ribo_pixmap(self.render_mode)
                q.save(); q.setOpacity(max(0.0, 1.0 - 0.5 * gl))
                q.drawPixmap(QRectF(St.x() - rpx, St.y() - rpx, 2 * rpx, 2 * rpx),
                             pix, QRectF(0, 0, pix.width(), pix.height()))
                q.restore()

        # intermediate filaments — a loose basket under the cortex
        if has_org and "interm_fil" in vis:
            q.setPen(QPen(col("ifil", int(90 * min(cyto, 1.4))), max(0.6, 0.9 * Z * cyto)))
            n_if = max(6, int(len(self.if_seed) * cyto))
            for a in self.if_seed[:n_if]:
                r1 = float(self.sb.radius_at(np.array([a]))[0]) * 0.92
                r2 = float(self.sb.radius_at(np.array([a + 1.9]))[0]) * 0.55
                p = QPainterPath(S(math.cos(a) * r1, math.sin(a) * r1))
                p.quadTo(S(0, 0), S(math.cos(a + 1.9) * r2, math.sin(a + 1.9) * r2))
                q.drawPath(p)

        # microtubules radiating from the centrosome (cold-labile); count and
        # weight scale with the line's cytoskeletal density
        if has_org and "microtubules" in vis and f["mt"] > 0.03:
            cr = float(self.sb.radius_at(np.array([self.centro_a]))[0]) * 0.30
            cxp, cyp = math.cos(self.centro_a) * cr, math.sin(self.centro_a) * cr
            q.setPen(QPen(col("mtub", int((40 + 150 * f["mt"]) * min(cyto, 1.3))),
                          max(0.7, 1.1 * Z * min(cyto, 1.4))))
            n_mt = max(8, int(len(self.mt_seed) * cyto))
            for j, jit in enumerate(self.mt_seed[:n_mt]):
                a = j / n_mt * 2 * math.pi + jit
                rr = float(self.sb.radius_at(np.array([a]))[0]) * (0.30 + 0.66 * f["mt"])
                p = QPainterPath(S(cxp, cyp))
                p.quadTo(S(math.cos(a) * rr * 0.5 + jit * 26, math.sin(a) * rr * 0.5),
                         S(math.cos(a) * rr, math.sin(a) * rr))
                q.drawPath(p)

        # ER — perinuclear rough-ER cisternae (membrane sheets studded with
        # ribosomes) continuous with the nuclear envelope, plus peripheral smooth
        # tubules reaching toward the cortex. Ca-store depletion fades it.
        if has_org and "er" in vis and self.render_mode != "phase":
            alpha = int(np.clip(120 + 110 * f["caER"], 70, 240))
            ercol = col("er", alpha); ribo = col("er", min(255, alpha + 25))
            # nested perinuclear sheets (partial arcs, so it reads as stacks not rings)
            for k in range(5):
                a0 = k * 1.15; span = 4.4 + 0.3 * k; base_r = 0.44 + k * 0.043
                p = QPainterPath(); studs = []
                N = 60
                for s in range(N + 1):
                    th = a0 + span * s / N
                    env = float(self.sb.radius_at(np.array([th]))[0])
                    rr = env * (base_r + 0.017 * math.sin(s * 0.7 + k * 1.3))
                    pt = S(math.cos(th) * rr, math.sin(th) * rr)
                    p.moveTo(pt) if s == 0 else p.lineTo(pt)
                    if s % 3 == 0: studs.append(pt)
                q.setBrush(Qt.BrushStyle.NoBrush)
                q.setPen(QPen(ercol, max(0.9, 1.5 * Z))); q.drawPath(p)
                if Z > 0.6:                              # ribosomes studding the rough ER
                    q.setPen(Qt.PenStyle.NoPen); q.setBrush(ribo)
                    rdot = max(0.8, 1.05 * Z)
                    for pt in studs: q.drawEllipse(pt, rdot, rdot)
            # peripheral smooth-ER tubules
            q.setBrush(Qt.BrushStyle.NoBrush)
            q.setPen(QPen(col("er", int(alpha * 0.65)), max(0.7, 1.0 * Z)))
            for k in range(8):
                a0 = k * 0.79
                p = QPainterPath()
                for s in range(18):
                    th = a0 + s * 0.045
                    env = float(self.sb.radius_at(np.array([th]))[0])
                    rr = env * (0.60 + s * 0.019 + 0.03 * math.sin(s * 1.1 + k))
                    pt = S(math.cos(th) * rr, math.sin(th) * rr)
                    p.moveTo(pt) if s == 0 else p.lineTo(pt)
                q.drawPath(p)

        # Golgi — stacked cisternae beside the nucleus
        if has_org and "golgi" in vis:
            gr = float(self.sb.radius_at(np.array([self.golgi_a]))[0]) * 0.52
            gx, gy = math.cos(self.golgi_a) * gr, math.sin(self.golgi_a) * gr
            q.setPen(QPen(col("golgi", 210), max(1.2, 2.4 * Z)))
            for k in range(5):
                p = QPainterPath()
                for s in range(11):
                    t = (s / 10 - 0.5) * 46
                    off = k * 7 - 14
                    p_ = S(gx + t * 0.95 - off * 0.25, gy + off + 0.010 * t * t)
                    p.moveTo(p_) if s == 0 else p.lineTo(p_)
                q.drawPath(p)

        # centrosome
        if has_org and "centrosome" in vis:
            cr = float(self.sb.radius_at(np.array([self.centro_a]))[0]) * 0.30
            cxp, cyp = math.cos(self.centro_a) * cr, math.sin(self.centro_a) * cr
            q.setBrush(col("centro", 200)); q.setPen(QPen(col("centro"), 1.4))
            for dd in (-4, 4):
                q.drawEllipse(S(cxp + dd, cyp + dd * 0.4), 4.5 * Z, 2.2 * Z)

        # small organelles — sized in true microns via px_um, so each keeps its
        # real-world aspect ratio to the whole cell. Organelles lose less water
        # than the cell, so they shrink only mildly under dehydration (Vn^(1/6)).
        vm, psi = f["Vmito"], f["dPsi"]
        dehyd = max(f.get("Vn", 1.0), 0.05) ** (1 / 6)
        for o in self.org:
            if not has_org or o.kind not in vis: continue
            pt = S(o.x, o.y)
            if o.kind == "mitochondria":
                L = MITO_UM[0] * px_um * vm ** (1 / 3) * o.scale * dehyd
                W = MITO_UM[1] * px_um * vm ** (1 / 3) * o.scale * dehyd * (0.7 + 0.5 * min(vm, 2.4))
                cc = col("mito", int(90 + 140 * psi)) if psi > .5 else col("muted", 150)
                q.save(); q.translate(pt); q.rotate(math.degrees(o.rot))
                q.setBrush(cc); q.setPen(QPen(cc.darker(130), max(0.8, 1.3 * Z)))
                q.drawEllipse(QPointF(0, 0), L * Z, W * Z)
                if vm < 1.6 and psi > 0.4 and Z > 0.7:
                    # cristae drawn as inner zig-zag folds (more realistic than bars)
                    q.setPen(QPen(cc.darker(150), max(0.5, 0.8 * Z)))
                    q.setBrush(Qt.BrushStyle.NoBrush)
                    ncr = 5; span = L * 1.5 * Z; amp = W * 0.62 * Z
                    path = QPainterPath(QPointF(-span / 2, 0))
                    for c2 in range(ncr + 1):
                        xx = -span / 2 + span * c2 / ncr
                        path.lineTo(QPointF(xx, (amp if c2 % 2 else -amp)))
                    q.drawPath(path)
                elif psi <= 0.4 and Z > 0.7:
                    # depolarised / swollen: rounder, cristae unfolded (blank)
                    pass
                q.restore()
            else:
                key = dict(lysosomes="lyso", peroxisomes="perox", endosomes="endo",
                           lipid_drop="lipid", vesicles="vesicle").get(o.kind, "muted")
                rr = ORG_UM[o.kind] * px_um * o.scale * dehyd
                q.setBrush(col(key, 190)); q.setPen(QPen(col(key), max(0.7, 1.1 * Z)))
                q.drawEllipse(pt, rr * Z, rr * Z)

        # CPA molecules. Density is normalised against the loaded target
        # concentration, not the freeze-concentrated peak. The old /14 was
        # tuned for the ~15 M freeze-concentration limit, so at a normal
        # loading of ~1-2 M only a dozen faint dots drew and the cytosol read
        # as empty even though a permeant CPA has fully equilibrated inside.
        # sqrt keeps both regimes on screen: a fully loaded cell sits near
        # half density, freeze-concentration saturates, and the dilution phase
        # (CPA trapped inside while the wash strips it outside) still shows.
        if "cpa" in vis and f["Cin"] > 0.02:
            ref = max(f.get("cpaLoad", 0.0), 0.5)
            nd = int(np.clip(math.sqrt(f["Cin"] / ref) / 2.0, 0, 1) * 220)
            q.setPen(Qt.PenStyle.NoPen); q.setBrush(col("cpa", 130))
            for k in range(nd):
                a, rr0 = self.cpa_seed[k]
                rr = float(self.sb.radius_at(np.array([a]))[0]) * rr0 * 0.94
                q.drawEllipse(S(math.cos(a) * rr, math.sin(a) * rr), 1.7 * Z, 1.7 * Z)

        # nucleus (absent in a mature red blood cell)
        if has_nucleus and "nucleus" in vis:
            np_ = self._path(self.nuc.p)
            rn = math.sqrt(self.nuc.area() / math.pi)
            g = QRadialGradient(self._to_screen(-rn * .3, -rn * .3), rn * 1.3 * Z)
            g.setColorAt(0, col("nucleus", 66)); g.setColorAt(1, col("nucleus", 26))
            q.setBrush(QBrush(g))
            if "nuclear_mem" in vis:
                q.setPen(QPen(col("nucleus"), max(1.2, 2.2 * Z)))
            else:
                q.setPen(Qt.PenStyle.NoPen)
            q.drawPath(np_)
            # nuclear pores
            if "nuclear_mem" in vis and Z > 1.3:
                q.setPen(QPen(col("nucleus", 220), max(1.0, 2.6 * Z)))
                for k in range(26):
                    a = k / 26 * 2 * math.pi
                    r2 = float(self.nuc.radius_at(np.array([a]))[0])
                    q.drawPoint(self._to_screen(math.cos(a) * r2, math.sin(a) * r2))
            cond = float(np.clip((f["Osme"] - 1) / 10, 0, 1))
            q.save(); q.setClipPath(np_)
            q.setPen(Qt.PenStyle.NoPen)
            # granular chromatin — denser and speckled, condensing with osmolality
            for a, rf, sz in self.chrom_seed:
                r2 = rf * rn * (0.94 - 0.30 * cond)
                q.setBrush(col("nucleus", int(45 + 120 * cond)))
                rr = (0.9 + 1.9 * cond) * sz * Z
                q.drawEllipse(self._to_screen(math.cos(a) * r2, math.sin(a) * r2), rr, rr)
            if "nucleoli" in vis:
                q.setBrush(col("nucleus", 120)); q.setPen(QPen(col("nucleus", 190), 1.2))
                for k, (aa, rr) in enumerate([(0.8, .38), (2.6, .30), (4.7, .34)]):
                    q.drawEllipse(self._to_screen(math.cos(aa) * rn * rr, math.sin(aa) * rn * rr),
                                  rn * 0.20 * Z, rn * 0.18 * Z)

            # ---- intracellular ice in the nucleus (freeze-substitution look) ----
            # Yu, Marquez-Curtis & Elliott 2026 (npj Imaging): visible IIF displaces
            # nuclear material into dark angular edges around light "hole"-like voids.
            #   small ice = many voids sectioned by dark edges (crystals >~1.1 um)
            #   big ice   = one large void with a single dark rim, material at border
            iif = float(np.clip(f.get("Piif", 0.0), 0, 1))
            iceph = f.get("frozen", False) or f.get("phase", "") in ("store", "warm", "melt")
            if iif > 0.12 and iceph and "iif" in vis:
                void = QColor(250, 246, 248)               # unstained ice-crystal void
                edge = col("nucleus").darker(160)          # condensed material at borders
                rsq = np.random.RandomState(97)
                if iif > 0.62:                             # single large crystal
                    rv = rn * (0.30 + 0.36 * iif)
                    poly = QPolygonF([self._to_screen(
                        math.cos(t) * rv * (0.82 + 0.34 * rsq.rand()),
                        math.sin(t) * rv * (0.82 + 0.34 * rsq.rand()))
                        for t in np.linspace(0, 2 * math.pi, 11)])
                    q.setBrush(void); q.setPen(QPen(edge, max(1.4, 2.3 * Z)))
                    q.drawPolygon(poly)
                else:                                      # many small crystals
                    q.setPen(QPen(edge, max(0.9, 1.4 * Z)))
                    for _ in range(int(4 + 16 * iif)):
                        a = rsq.uniform(0, 2 * math.pi); r0 = math.sqrt(rsq.rand()) * rn * 0.80
                        cx0, cy0 = math.cos(a) * r0, math.sin(a) * r0
                        sz = rn * (0.06 + 0.11 * rsq.rand()) * (0.8 + 0.6 * iif)
                        m = int(rsq.randint(4, 7))         # angular ice facets
                        poly = QPolygonF([self._to_screen(
                            cx0 + math.cos(t) * sz * (0.7 + 0.5 * rsq.rand()),
                            cy0 + math.sin(t) * sz * (0.7 + 0.5 * rsq.rand()))
                            for t in np.linspace(0, 2 * math.pi, m + 1)])
                        q.setBrush(void); q.drawPolygon(poly)
            q.restore()

            # ---- LINC complex & nuclear mechanotransduction ----
            # Actomyosin force (p-MLC), transmitted through the LINC complex
            # (nesprin-SUN) wherever the actin cortex is intact, loads the nuclear
            # lamina and drives mechanosensitive signalling. Kirby & Lammerding
            # 2018; Swift 2013 (lamin A/C as a mechanical set-point).
            force = f.get("pMLC", 0.0) * (0.30 + 0.70 * f.get("actin", 0.0))
            linc = f.get("actin", 0.0)
            # lamin A/C ring just inside the envelope; tints amber as force loads it
            lamc = col("path") if force < 0.4 else col("warn")
            q.setPen(QPen(lamc, max(1.0, 1.5 * Z))); q.setBrush(Qt.BrushStyle.NoBrush)
            q.drawPath(self._path(self.nuc.p * 0.88))
            # LINC tethers (nesprin-SUN) spanning the envelope to the cytoskeleton
            if linc > 0.15:
                q.setPen(QPen(col("actin", int(110 + 130 * linc)), max(0.7, 1.0 * Z)))
                for k in range(16):
                    a = k / 16 * 2 * math.pi
                    r0 = float(self.nuc.radius_at(np.array([a]))[0])
                    q.drawLine(self._to_screen(math.cos(a) * r0 * 0.98, math.sin(a) * r0 * 0.98),
                               self._to_screen(math.cos(a) * r0 * (1.18 + 0.10 * linc),
                                               math.sin(a) * r0 * (1.18 + 0.10 * linc)))
            # actomyosin force arrows pulling the nucleus inward via the LINC
            if force > 0.25:
                q.setPen(QPen(col("mito", int(130 + 120 * min(force, 1))), max(1.3, 1.9 * Z)))
                for k in range(6):
                    a = k / 6 * 2 * math.pi + 0.35
                    r0 = float(self.nuc.radius_at(np.array([a]))[0])
                    tip = self._to_screen(math.cos(a) * r0 * 1.04, math.sin(a) * r0 * 1.04)
                    tail = self._to_screen(math.cos(a) * r0 * (1.32 + 0.45 * force),
                                           math.sin(a) * r0 * (1.32 + 0.45 * force))
                    q.drawLine(tail, tip)
                    ox, oy, px, py = math.cos(a), math.sin(a), -math.sin(a), math.cos(a)
                    for s in (1, -1):
                        q.drawLine(tip, QPointF(tip.x() + ox * 7 + px * s * 4,
                                                tip.y() + oy * 7 + py * s * 4))
                # mechanosensitive-transcription glow when force reaches the nucleus
                q.setPen(Qt.PenStyle.NoPen); q.setBrush(col("path", int(55 * min(force, 1))))
                q.drawEllipse(self._to_screen(0, 0), rn * 0.42 * Z, rn * 0.42 * Z)

        # intracellular ice — angular crystals nucleated through the cytoplasm when
        # the cell is UNPROTECTED (high P_iif). CPA drives P_iif → 0, so this whole
        # field is the visual signature of the "no protection" freezing regime and
        # vanishes when a cryoprotectant is present. While frozen (cool→store→warm)
        # they read as ice (blue); after melt they persist as the pale unstained
        # "void / freeze-damage" left behind in fixed, thawed cells (the cytoplasmic
        # counterpart of the nuclear freeze-substitution look above, and of the
        # interior perforation seen in thawed spheroids, Gao/Bissoyi 2024).
        if f["Piif"] > 0.05 and "iif" in vis:
            pf = float(np.clip(f["Piif"], 0, 1))
            nX = int(6 + pf * 46)
            rsq = np.random.RandomState(23)
            post_thaw = f.get("phase", "") in ("melt", "dilute", "recover", "end")
            if post_thaw:
                void_fill = QColor(250, 246, 248); void_fill.setAlpha(int(70 + 110 * pf))
                void_edge = QColor(150, 152, 160, int(150 + 90 * pf))
                q.setPen(QPen(void_edge, max(1.0, 1.5 * Z)))
            else:
                q.setPen(QPen(col("ice", int(170 + 80 * pf)), max(1.0, 1.5 * Z)))
            for k in range(nX):
                a = k * 2.39996
                # surface-catalysed nucleation (Li et al. 2020; Toner 1990): ice
                # nucleates at the inner membrane surface and grows inward, so it
                # is densest at the periphery and thins toward the centre.
                r2 = rmean * (0.30 + 0.66 * rsq.rand() ** 0.45)
                cx, cy = math.cos(a) * r2, math.sin(a) * r2
                cr = (2.2 + pf * 6.5) * (0.6 + 0.8 * rsq.rand())
                m = int(rsq.randint(5, 7))                 # angular ice facets
                q.setBrush(void_fill if post_thaw else col("ice", int(55 + 95 * pf)))
                q.drawPolygon(QPolygonF([self._to_screen(
                    cx + math.cos(t) * cr * (0.7 + 0.5 * rsq.rand()),
                    cy + math.sin(t) * cr * (0.7 + 0.5 * rsq.rand()))
                    for t in np.arange(m) * 2 * math.pi / m + 0.3]))

    def _draw_membrane(self, q, f, path, rmean):
        Z = self.zoom
        pore, gel, tens = f["pore"], f["gel"], f["tension"]
        lod = Z > 3.5
        base_w = (3.2 - 1.3 * pore) * (1 + .45 * gel)
        pen = QPen(col("crit") if pore > .55 else col("serious") if pore > .2
                   else col("nucleus") if gel > .5 else col("text"),
                   max(0.8, base_w * (0.55 if lod else min(Z, 2.2))))
        if pore > .12:
            pen.setStyle(Qt.PenStyle.DashLine)
        elif gel > .5:
            pen.setStyle(Qt.PenStyle.DotLine)
        q.setBrush(Qt.BrushStyle.NoBrush); q.setPen(pen); q.drawPath(path)

        # ---- gel-phase rigidity: a cold membrane reads as a solid crystalline
        #      shell, not a wobbly line. Below Tm the bilayer is a 2D solid
        #      (Dimova 2014). Overdraw a thick, cool-toned stroke and short
        #      radial "frost" ticks that make the rigidity legible.
        if gel > 0.35:
            q.setPen(QPen(col("nucleus", int(110 + 130 * gel)),
                          max(1.6, (2.4 + 2.2 * gel) * min(Z, 2.0))))
            q.setBrush(Qt.BrushStyle.NoBrush); q.drawPath(path)
            if gel > 0.55 and Z > 0.6:
                n = self.sb.n; stp = max(1, n // 44)
                q.setPen(QPen(col("nucleus", int(80 + 110 * gel)), max(0.8, 1.0 * Z)))
                o = 3.2
                for i in range(0, n, stp):
                    p0 = self.sb.p[i]; ang = math.atan2(p0[1], p0[0])
                    cx, cy = math.cos(ang) * o, math.sin(ang) * o
                    q.drawLine(self._to_screen(p0[0] - cx, p0[1] - cy),
                               self._to_screen(p0[0] + cx, p0[1] + cy))

        # ---- membrane damage: visible at ANY zoom, scaled by the instantaneous
        #      defect index and the accumulated membrane damage Dmem. Once the
        #      bilayer enters the MD pore regime / phase-transition leak
        #      (Leontiadou 2004; the model's bilayer_regime), the contour breaks
        #      into discrete lesions that leak intracellular solute outward.
        dmg = float(np.clip(pore + 0.7 * f.get("Dmem", 0.0), 0, 1))
        if dmg > 0.05:
            n = self.sb.n
            nles = int(dmg * 24)
            for k in range(nles):
                p0 = self.sb.p[(k * 23) % n]
                rr = (1.8 + 4.5 * dmg) * max(Z, 0.6)
                q.setPen(QPen(col("crit", int(120 + 130 * dmg)), max(1.0, 1.4 * min(Z, 2.2))))
                q.setBrush(col("surface1", 200))
                q.drawEllipse(self._to_screen(p0[0], p0[1]), rr, rr)
                if dmg > 0.28:                       # frank leak of contents
                    ang = math.atan2(p0[1], p0[0])
                    q.setPen(Qt.PenStyle.NoPen)
                    for j in (1, 2, 3):
                        lp = (p0[0] + math.cos(ang) * j * 4.5, p0[1] + math.sin(ang) * j * 4.5)
                        q.setBrush(col("cpa", max(0, int(95 - j * 26))))
                        q.drawEllipse(self._to_screen(lp[0], lp[1]), 1.4 * max(Z, 0.6), 1.4 * max(Z, 0.6))
            q.setBrush(Qt.BrushStyle.NoBrush)

        if f["pMLC"] > .15:
            q.setPen(QPen(col("mito", int(np.clip((f["pMLC"] - .15) * 240, 0, 210))),
                          max(1.0, (1.4 + 2.2 * min(f["pMLC"], 1.2))
                              * (0.5 if lod else min(Z, 2.0)))))
            q.drawPath(path)

        # high zoom: resolve the bilayer into leaflets, with CPA in the headgroups
        if lod:
            n = self.sb.n
            step = max(1, n // 96)
            # bilayer thickness is ~5 nm, invisible at cell scale, so it is drawn
            # exaggerated at a fixed world size and the two leaflets separated
            th = 3.0 * f["thick"]
            for i in range(0, n, step):
                p0 = self.sb.p[i]; p1 = self.sb.p[(i + step) % n]
                d = p1 - p0; L = math.hypot(*d) or 1e-6
                nx, ny = d[1] / L, -d[0] / L
                if p0[0] * nx + p0[1] * ny < 0: nx, ny = -nx, -ny
                for sgn, alpha in ((1, 215), (-1, 165)):
                    q.setPen(QPen(col("text2", alpha), max(0.9, 0.55 * Z)))
                    q.drawLine(self._to_screen(p0[0] + nx * sgn * th, p0[1] + ny * sgn * th),
                               self._to_screen(p1[0] + nx * sgn * th, p1[1] + ny * sgn * th))
                if f["Cin"] > 0.05 and (i // step) % 2 == 0:
                    q.setPen(Qt.PenStyle.NoPen); q.setBrush(col("cpa", 170))
                    o = th * 0.55 * (1 if (i // step) % 4 else -1)
                    q.drawEllipse(self._to_screen(p0[0] + nx * o, p0[1] + ny * o),
                                  0.9 * Z, 0.9 * Z)
            if pore > 0.05:
                q.setPen(QPen(col("crit", 210), max(1.0, 0.7 * Z)))
                q.setBrush(Qt.BrushStyle.NoBrush)
                for k in range(int(pore * 14)):
                    i = (k * 17) % n
                    p0 = self.sb.p[i]
                    q.drawEllipse(self._to_screen(*p0), th * 0.9 * Z, th * 0.9 * Z)

        if "cortex" in self.visible and f["actin"] > 0.05:
            cyto = float(np.clip(f.get("cyto", 1.0), 0.3, 2.0))
            q.setPen(QPen(col("actin", int((60 + 120 * f["actin"]) * min(cyto, 1.4))),
                          max(0.8, 1.6 * Z * min(cyto, 1.4))))
            for k in range(max(1, int(3 * cyto))):
                sc = 0.955 - k * 0.022
                pts = self.sb.p * sc
                q.drawPath(self._path(pts))

        if tens > .5:
            q.setPen(QPen(col("crit", int(np.clip((tens - .5) * 300, 0, 180))), max(1.0, 1.2 * Z)))
            q.setBrush(Qt.BrushStyle.NoBrush)
            q.drawEllipse(self._to_screen(0, 0), (rmean + 6) * Z, (rmean + 6) * Z)

        ap, ne = f["apop"], f["necr"]
        if ap + ne > .02:
            q.setPen(Qt.PenStyle.NoPen)
            q.setBrush(col("crit" if ne > ap else "warn", int(np.clip((ap + ne) * 78, 0, 78))))
            q.drawPath(path)

    def _draw_overlay(self, q, f, px_um, rmean):
        h, w = self.height(), self.width()
        ppum = px_um * self.zoom                       # screen pixels per micron

        # ---- calibrated scale bar: choose a round micron length ~90 px long
        target = 90.0 / max(ppum, 1e-6)
        nice = min([1, 2, 5, 10, 20, 25, 50, 100, 200], key=lambda v: abs(v - target))
        barpx = nice * ppum
        q.setPen(QPen(col("text"), 2)); q.setFont(QFont("", 9))
        y0 = h - 24
        q.drawLine(QPointF(18, y0), QPointF(18 + barpx, y0))
        q.drawLine(QPointF(18, y0 - 4), QPointF(18, y0 + 4))
        q.drawLine(QPointF(18 + barpx, y0 - 4), QPointF(18 + barpx, y0 + 4))
        q.setPen(col("text2"))
        q.drawText(QRectF(18, y0 + 3, barpx, 16), Qt.AlignmentFlag.AlignHCenter, f"{nice} µm")

        if not self.sci:
            q.setPen(col("muted")); q.setFont(QFont("", 9))
            mem_per = 2 * math.pi * self.R0 * f["APL"]
            slack = max(0.0, mem_per / (2 * math.pi * max(rmean, 1)) - 1) * 100
            bits = [f"zoom {self.zoom:.1f}×", f"membrane slack {slack:.0f}%"]
            if self.sb.blebs: bits.append(f"{len(self.sb.blebs)} blebs")
            if self.ice.active:
                bits.append(f"pocket {self.ice.pocket_mean / max(px_um,1e-6):.1f} µm")
            q.drawText(QPointF(w - 260, h - 12), "  ·  ".join(bits))
            if self.hover:
                q.setPen(col("text2")); q.setFont(QFont("", 10, QFont.Weight.Bold))
                q.drawText(QPointF(18, 22), self.hover.replace("_", " "))
            return

        # ============================ SCIENTIFIC OVERLAY ============================
        S = self._to_screen
        cx, cy = S(0, 0).x(), S(0, 0).y()
        Rs = rmean * self.zoom
        dia = 2 * rmean / max(px_um, 1e-6)              # cell diameter, microns

        # title strip: phase, temperature, time
        names = dict(dry1="Primary drying", dry2="Secondary drying", drystore="Dry storage",
                     rehydrate="Rehydration",
                     load="CPA loading", cool="Cooling", seed="Seeding", store="Storage",
                     warm="Warming", melt="Melting", dilute="Dilution", recover="Recovery", end="End")
        q.setPen(col("text")); q.setFont(QFont("", 12, QFont.Weight.Bold))
        q.drawText(QPointF(18, 26), names.get(f.get("phase",""), f.get("phase","")))
        q.setPen(col("text2")); q.setFont(QFont("", 11))
        q.drawText(QPointF(18, 44), f"T = {f['T']:.1f} °C")

        # leader-lined labels: (anchor world/screen point, corner text)
        def label(ax, ay, tx, ty, lines, ckey=None):
            # leader line from structure anchor to the label block
            q.setPen(QPen(col(ckey, 230) if ckey else col("borderStrong", 200), 1.3))
            q.drawLine(QPointF(ax, ay), QPointF(tx, ty))
            q.setBrush(col("surface1", 235)); q.setPen(QPen(col("border"), 1))
            bw = 4 + max(len(s) for s in lines) * 6.4
            bh = 4 + len(lines) * 13
            ox = tx if tx < w / 2 else tx - bw
            q.drawRoundedRect(QRectF(ox, ty - 10, bw, bh), 4, 4)
            q.setFont(QFont("", 8, QFont.Weight.Bold)); q.setPen(col("text"))
            q.drawText(QRectF(ox + 4, ty - 9, bw, 12), Qt.AlignmentFlag.AlignLeft, lines[0])
            q.setFont(QFont("", 8)); q.setPen(col("text2"))
            for i, s in enumerate(lines[1:], 1):
                q.drawText(QRectF(ox + 4, ty - 9 + i * 13, bw, 12), Qt.AlignmentFlag.AlignLeft, s)

        d = 0.7071
        # cell-type morphology gates: only label structures the cell actually has
        # (mammalian RBC is anucleate + organelle-free; platelet is anucleate).
        lc = f.get("cell_type", "msc")
        lbl_nucleus = lc not in ("rbc", "platelet")
        lbl_org = lc != "rbc"
        # cell / membrane  (top)
        label(cx, cy - Rs, w * 0.62, 74,
              ["Whole cell", f"{f['Vn']*100:.0f}% Viso · Ø {dia:.1f} µm",
               f"membrane: gel {f['gel']*100:.0f}%"])
        # nucleus (centre) — only for nucleated cells
        if lbl_nucleus:
            label(cx, cy, 24, h * 0.40,
                  ["Nucleus", f"{f['Vnuc']*100:.0f}% resting vol"], "nucleus")
        # mitochondria — pick a drawn one; report its true physical size
        mpt = next((o for o in self.org if o.kind == "mitochondria"), None)
        if lbl_org and mpt is not None:
            mp = S(mpt.x, mpt.y); ms = f["Vmito"] ** (1 / 3)
            label(mp.x(), mp.y(), w - 4, h * 0.30,
                  ["Mitochondria", f"{2*MITO_UM[0]*ms:.1f} × {2*MITO_UM[1]*ms:.1f} µm",
                   f"ΔΨm {f['dPsi']*100:.0f}% · MPT {f['mpt']*100:.0f}%"], "mito")
        # cytoskeleton — ties the cell view to the FA·LINC mechanics (cold-labile
        # actin/MT depolymerise; vimentin persists). Skip for the mammalian RBC,
        # whose skeleton is a spectrin membrane mesh, not actin/MT/vimentin.
        if lbl_org:
            label(cx - Rs * 0.45, cy - Rs * 0.55, 24, h * 0.22,
                  ["Cytoskeleton", f"actin {f.get('actin',0)*100:.0f}% · MT {f.get('mt',0)*100:.0f}%",
                   f"vimentin {f.get('intf',0)*100:.0f}% (cold-stable)"], "ifil")
        # cytosol / CPA (interior) — with the transmembrane CPA gradient (in vs out)
        label(cx - Rs * 0.35, cy + Rs * 0.35, 24, h * 0.62,
              ["Cytosol", f"CPA in {f['Cin']:.1f} M · out {f.get('Cout',0):.1f} M",
               f"glass {f.get('glass',0)*100:.0f}% · protein {f.get('prot',1)*100:.0f}%"], "cpa")
        # extracellular ice (outside, when frozen)
        if self.ice.active:
            ax, ay = cx + Rs * d + 30, cy - Rs * d - 20
            label(ax, ay, w - 4, 74,
                  ["Extracellular ice", f"{f['fIce']*100:.0f}% frozen · {f['Osme']:.1f} osmol/L",
                   f"grain {f['grain']:.0f} µm"], "ice")

        # ---- compact legend (bottom-right)
        leg = [("Nucleus", "nucleus"), ("Mitochondria", "mito"), ("ER", "er"),
               ("CPA", "cpa"), ("Ice / freeze-damage voids", "ice")]
        q.setFont(QFont("", 8)); lx, ly = w - 186, h - 92
        q.setBrush(col("surface1", 220)); q.setPen(QPen(col("border"), 1))
        q.drawRoundedRect(QRectF(lx - 8, ly - 12, 186, len(leg) * 15 + 10), 5, 5)
        for i, (nm, ck) in enumerate(leg):
            yy = ly + i * 15
            q.setBrush(col(ck)); q.setPen(Qt.PenStyle.NoPen)
            q.drawEllipse(QPointF(lx, yy), 4, 4)
            q.setPen(col("text2")); q.drawText(QPointF(lx + 10, yy + 4), nm)

        # phase-diagram inset last, so it sits on top of any leader labels
        if self.ice.active or f.get("T", 22) < 5:
            self._draw_phase_diagram(q, f)
