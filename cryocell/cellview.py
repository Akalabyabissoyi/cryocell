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

# isotonic cell radius (um) of the calibration cell (Viso = 1800 um^3); the view
# draws at a fixed px/um referenced to this, so different cell sizes are visible.
R_ISO_REF = (3 * 1800.0 / (4 * math.pi)) ** (1 / 3)

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
]


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
        self.sci = False             # scientific overlay: scale bar, labels, cryo-stage ice
        self.render_mode = "illustrative"   # illustrative | fluor | phase
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
        self.mt_seed = rs.uniform(-0.25, 0.25, 64)
        self.if_seed = rs.uniform(0, 2 * np.pi, 40)
        # cytoplasmic crowding / ribosome texture: seeded polar field mapped into
        # the deforming cell each frame (angle, radius, size, shade)
        nrib = 620
        self.ribo = np.stack([rs.uniform(0, 2 * np.pi, nrib),
                              np.sqrt(rs.uniform(0.02, 0.95, nrib)),
                              rs.uniform(0.6, 1.9, nrib),
                              rs.uniform(0.25, 0.75, nrib)], axis=1)
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
        self.render_mode = mode if mode in ("illustrative", "fluor", "phase") else "illustrative"
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
            c = QColor(base); c.setAlpha(int((46 if mode == "fluor" else 34) * sh))
            p.setBrush(c)
            p.drawEllipse(QPointF(c0 + math.cos(a) * R0 * rf, c0 + math.sin(a) * R0 * rf),
                          sz * 0.9, sz * 0.9)
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
        damp_g  = 0.88 * (1 - 0.80 * gel)
        noise_g = (0.055 * slack + 0.005) * (1 - 0.92 * gel)
        self.sb.update_blebs(f["bleb"])
        pocket = self.ice.pocket_radius if (f["frozen"] and self.ice.active) else None
        steps = 22 if self._settle > 0 else 7
        if self._settle > 0: self._settle -= steps
        for _ in range(steps):
            self.sb.step(L0, kT, kB, 0.60, math.pi * R * R, damp_g, slack,
                         k_cortex=kC, pocket_r=pocket, crystals=crystals,
                         noise=noise_g, external=self._ext)
            Rn = self.R0 * 0.40 * max(f["Vnuc"], 0.1) ** (1 / 3)
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
            q.save(); q.setPen(Qt.PenStyle.NoPen); q.setBrush(col("surface1"))
            r = self.sb.radii().max() * 1.5
            q.drawEllipse(S(0, 0), r * Z, r * Z); q.restore()

        # ---- cell body
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

    def _draw_interior(self, q, f, rmean, px_um):
        S, Z = self._to_screen, self.zoom
        vis = self.visible
        # resting cytoskeletal density of this cell line: mesenchymal cells
        # (hMSC) carry a denser network than epithelial lines (HeLa, A549).
        cyto = float(np.clip(f.get("cyto", 1.0), 0.3, 2.0))

        # ---- cytoplasmic crowding / ribosome texture (behind organelles) ----
        # rendered once into a cached pixmap and blitted, scaled to the current
        # cell radius and clipped to the cell path — O(1) per frame. Quenched as
        # the cytosol vitrifies (glass); skipped in clean phase-contrast.
        if self.render_mode != "phase":
            gl = f.get("glass", 0.0)
            pix = self._ribo_pixmap(self.render_mode)
            St = self._to_screen(0, 0); rpx = rmean * Z
            q.save(); q.setOpacity(max(0.0, 1.0 - 0.5 * gl))
            q.drawPixmap(QRectF(St.x() - rpx, St.y() - rpx, 2 * rpx, 2 * rpx),
                         pix, QRectF(0, 0, pix.width(), pix.height()))
            q.restore()

        # intermediate filaments — a loose basket under the cortex
        if "interm_fil" in vis:
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
        if "microtubules" in vis and f["mt"] > 0.03:
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

        # ER — a tubular network around the nucleus
        if "er" in vis:
            alpha = int(np.clip(70 + 150 * f["caER"], 40, 235))
            q.setPen(QPen(col("er", alpha), max(1.0, 1.9 * Z)))
            for k in range(9):
                a0 = k * 0.70
                p = QPainterPath()
                for s in range(24):
                    th = a0 + s * 0.085
                    rr = float(self.sb.radius_at(np.array([th]))[0]) * (0.40 + 0.16 * math.sin(s * .8 + k))
                    pt = S(math.cos(th) * rr, math.sin(th) * rr)
                    p.moveTo(pt) if s == 0 else p.lineTo(pt)
                q.drawPath(p)

        # Golgi — stacked cisternae beside the nucleus
        if "golgi" in vis:
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
        if "centrosome" in vis:
            cr = float(self.sb.radius_at(np.array([self.centro_a]))[0]) * 0.30
            cxp, cyp = math.cos(self.centro_a) * cr, math.sin(self.centro_a) * cr
            q.setBrush(col("centro", 200)); q.setPen(QPen(col("centro"), 1.4))
            for dd in (-4, 4):
                q.drawEllipse(S(cxp + dd, cyp + dd * 0.4), 4.5 * Z, 2.2 * Z)

        # small organelles
        vm, psi = f["Vmito"], f["dPsi"]
        for o in self.org:
            if o.kind not in vis: continue
            pt = S(o.x, o.y)
            if o.kind == "mitochondria":
                L = 11.5 * vm ** (1 / 3) * o.scale
                W = 5.0 * vm ** (1 / 3) * o.scale * (0.6 + 0.55 * min(vm, 2.4))
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
                rr = dict(lysosomes=6.0, peroxisomes=4.2, endosomes=5.0,
                          lipid_drop=5.4, vesicles=3.0)[o.kind] * o.scale
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

        # nucleus
        if "nucleus" in vis:
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

        # intracellular ice
        if f["Piif"] > 0.02:
            q.setPen(QPen(col("ice", int(255 * min(f["Piif"], 1))), max(1.2, 1.8 * Z)))
            q.setBrush(col("ice", int(70 * f["Piif"])))
            nX = int(f["Piif"] * 11)
            for k in range(nX):
                a = k * 2.39996; r2 = math.sqrt((k + .5) / 12) * rmean * .7
                cx, cy = math.cos(a) * r2, math.sin(a) * r2
                cr = 3 + f["Piif"] * 9
                q.drawPolygon(QPolygonF([self._to_screen(cx + math.cos(t) * cr,
                                                         cy + math.sin(t) * cr)
                                         for t in np.arange(6) * math.pi / 3 + .3]))

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
        names = dict(load="CPA loading", cool="Cooling", seed="Seeding", store="Storage",
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
        # cell / membrane  (top)
        label(cx, cy - Rs, w * 0.62, 74,
              ["Whole cell", f"{f['Vn']*100:.0f}% Viso · Ø {dia:.1f} µm",
               f"membrane: gel {f['gel']*100:.0f}%"])
        # nucleus (centre)
        label(cx, cy, 24, h * 0.40,
              ["Nucleus", f"{f['Vnuc']*100:.0f}% resting vol"], "nucleus")
        # mitochondria — pick a drawn one
        mpt = next((o for o in self.org if o.kind == "mitochondria"), None)
        if mpt is not None:
            mp = S(mpt.x, mpt.y)
            label(mp.x(), mp.y(), w - 4, h * 0.34,
                  ["Mitochondria", f"ΔΨm {f['dPsi']*100:.0f}% · MPT {f['mpt']*100:.0f}%"], "mito")
        # cytosol / CPA (interior)
        label(cx - Rs * 0.35, cy + Rs * 0.35, 24, h * 0.62,
              ["Cytosol", f"CPA {f['Cin']:.1f} M · glass {f.get('glass',0)*100:.0f}%",
               f"protein {f.get('prot',1)*100:.0f}% native"], "cpa")
        # extracellular ice (outside, when frozen)
        if self.ice.active:
            ax, ay = cx + Rs * d + 30, cy - Rs * d - 20
            label(ax, ay, w - 4, 74,
                  ["Extracellular ice", f"{f['fIce']*100:.0f}% frozen · {f['Osme']:.1f} osmol/L",
                   f"grain {f['grain']:.0f} µm"], "ice")

        # ---- compact legend (bottom-right)
        leg = [("Nucleus", "nucleus"), ("Mitochondria", "mito"), ("ER", "er"),
               ("CPA", "cpa"), ("Ice", "ice")]
        q.setFont(QFont("", 8)); lx, ly = w - 96, h - 92
        q.setBrush(col("surface1", 220)); q.setPen(QPen(col("border"), 1))
        q.drawRoundedRect(QRectF(lx - 8, ly - 12, 96, len(leg) * 15 + 10), 5, 5)
        for i, (nm, ck) in enumerate(leg):
            yy = ly + i * 15
            q.setBrush(col(ck)); q.setPen(Qt.PenStyle.NoPen)
            q.drawEllipse(QPointF(lx, yy), 4, 4)
            q.setPen(col("text2")); q.drawText(QPointF(lx + 10, yy + 4), nm)
