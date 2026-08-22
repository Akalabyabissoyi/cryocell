"""CryoCell desktop application."""
from __future__ import annotations
import sys, math, traceback, os, shutil, subprocess
import numpy as np
from PyQt6.QtCore import Qt, QThread, pyqtSignal, QTimer, QPointF, QRectF
from PyQt6.QtGui import QFont, QColor, QAction, QImage, QPainter, QPen, QPainterPath, QPolygonF
from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout,
    QHBoxLayout, QGridLayout, QLabel, QSlider, QComboBox, QPushButton, QGroupBox,
    QScrollArea, QCheckBox, QTabWidget, QTextEdit, QSplitter, QSizePolicy,
    QTableWidget, QTableWidgetItem, QHeaderView, QProgressBar, QFileDialog)
import pyqtgraph as pg

from .model import Params, simulate, CPAS, ADHESION_STATES, ADDITIVES, osm_from_dT
from .analysis import (simulate_population, knockout_screen, stability_check,
                       next_experiment)
from .hpa import HPA, HPA_TOTALS, GEOMETRY, summary_line
from .cellview import CellView, LAYERS, C

pg.setConfigOptions(antialias=True, background=C["surface1"], foreground=C["text2"])
SER = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#4a3aa7"]
def clamp(v, a, b): return max(a, min(b, v))

PRESETS = {
    "MSC bank standard":            dict(cpa_key="dmso", conc_pct=10, T_add=22, hold_min=10, add_steps=3, T_seed=-6, CR=1, T_store=-196, WR=400, dilution="step", adhesion="suspension", rock_drug="none"),
    "MSC + ROCK inhibitor":         dict(cpa_key="dmso", conc_pct=10, T_add=22, hold_min=10, add_steps=3, T_seed=-6, CR=1, T_store=-196, WR=400, dilution="step", adhesion="suspension", rock_drug="y27632", rock_conc=10, rock_when="both"),
    "Frozen attached, sheared":     dict(cpa_key="dmso", conc_pct=10, T_add=22, hold_min=10, add_steps=3, T_seed=-6, CR=1, T_store=-196, WR=400, dilution="step", adhesion="sheared"),
    "Clinical HSC, infused":        dict(cpa_key="dmso", conc_pct=10, T_add=4,  hold_min=5,  add_steps=1, T_seed=-4, CR=1, T_store=-196, WR=400, dilution="none", adhesion="suspension"),
    "Rapid cool / rapid warm":      dict(cpa_key="dmso", conc_pct=10, T_add=22, hold_min=10, add_steps=1, T_seed=-4, CR=100, T_store=-196, WR=3000, dilution="step"),
    "Vitrification-like":           dict(cpa_key="dmso", conc_pct=25, T_add=4,  hold_min=8,  add_steps=5, T_seed=-20, CR=500, T_store=-196, WR=10000, dilution="step", iri=60),
    "DMSO-free: trehalose + IRI":   dict(cpa_key="tre",  conc_pct=10, T_add=22, hold_min=15, add_steps=3, T_seed=-6, CR=1, T_store=-196, WR=400, dilution="none", iri=85),
    "Combo: reduced DMSO + trehalose": dict(cpa_key="dmso", conc_pct=7, additive="tre", add_conc=5, T_add=22, hold_min=10, add_steps=3, T_seed=-6, CR=1, T_store=-196, WR=400, dilution="step", adhesion="suspension"),
    "Combo: DMSO + polyampholyte":  dict(cpa_key="dmso", conc_pct=7, additive="pll", add_conc=6, T_add=22, hold_min=10, add_steps=3, T_seed=-6, CR=1, T_store=-196, WR=400, dilution="step", adhesion="suspension"),
    "Combo: DMSO + HES (clinical)": dict(cpa_key="dmso", conc_pct=5, additive="hes", add_conc=6, T_add=4, hold_min=5, add_steps=1, T_seed=-4, CR=1, T_store=-196, WR=400, dilution="none", adhesion="suspension"),
    "Self-deactivating CPA kit":    dict(cpa_key="sd",   conc_pct=12, T_add=4,  hold_min=8,  add_steps=3, T_seed=-6, CR=1, T_store=-196, WR=400, dilution="none", thalf=8),
    "Everything wrong":             dict(cpa_key="dmso", conc_pct=10, T_add=37, hold_min=45, add_steps=1, T_seed=-1, CR=0.2, T_store=-80, WR=3, dilution="direct"),
    # cell-line profiles — now with each line's characteristic SIZE (Viso) and
    # cytoskeletal density (cyto) as well as its death phenotype, so switching
    # line visibly changes the cell, not just the survival number. hMSC is the
    # large, cytoskeleton-rich calibrated primary baseline (Viso 1800, all axes
    # 0); HeLa and A549 are smaller epithelial cancer lines that resist death.
    "Cell line — hMSC (primary)":   dict(Viso=1800, cyto=1.35, sterol=25, adhesion="suspension",
                                         apop_resist=0.0,  anoikis_resist=0.0,  glycolytic=0.0,  antioxidant=0.0),
    "Cell line — HeLa":             dict(Viso=1400, cyto=1.00, sterol=30, adhesion="suspension",
                                         apop_resist=0.75, anoikis_resist=0.80, glycolytic=0.60, antioxidant=0.35),
    "Cell line — A549":             dict(Viso=1100, cyto=0.85, sterol=28, adhesion="suspension",
                                         apop_resist=0.40, anoikis_resist=0.55, glycolytic=0.55, antioxidant=0.85),
}

# (attr, label, min, max, step, decimals, log)
SLIDERS = [
    ("conc_pct",  "CPA concentration",         0,   30,  0.5, 1, False, "%"),
    ("add_conc",  "Additive concentration",    0,   20,  0.5, 1, False, "% w/v"),
    ("sterol",    "Membrane sterol",           0,   45,  1,   0, False, " mol%"),
    ("iri",       "Added IRI activity",        0,   95,  1,   0, False, "%"),
    ("T_add",     "CPA addition temperature",  0,   37,  1,   0, False, " °C"),
    ("hold_min",  "Equilibration hold",        0,   60,  1,   0, False, " min"),
    ("T_seed",    "Ice seeding temperature",  -20,  -1,  0.5, 1, False, " °C"),
    ("CR",        "Cooling rate",             -1,    3,  0.02,2, True,  " °C/min"),
    ("T_store",   "Storage temperature",     -196, -20,  1,   0, False, " °C"),
    ("days",      "Storage time",             -1,    3,  0.05,2, True,  " d"),
    ("WR",        "Warming rate",              0,    4,  0.02,1, True,  " °C/min"),
    ("rock_conc", "ROCK inhibitor conc.",      0,  100,  1,   0, False, " µM"),
    ("zvad",      "z-VAD-fmk",                 0,   95,  5,   0, False, "% inhib"),
    ("calpain_i", "Calpain inhibitor",         0,   95,  5,   0, False, "% inhib"),
    ("gsmtx",     "GsMTx4 (Piezo1)",           0,   95,  5,   0, False, "% block"),
    ("T_recover", "Post-thaw hold temperature",4,   37,  1,   0, False, " °C"),
    ("recover_h", "Recovery window",           1,   48,  1,   0, False, " h"),
    ("lp",        "Lp at 25 °C",               0.02, 0.8,0.01,2, False, " µm/min/atm"),
    ("ps",        "Ps at 25 °C",               0.005,0.4,0.005,3,False, " µm/s"),
    ("vb",        "Inactive volume Vb",        0.10, 0.45,0.01,2,False, ""),
    ("Viso",      "Isotonic cell volume",      800, 4000, 50, 0, False, " µm³"),
    ("thalf",     "SD linker t½ at 37 °C",     0.5, 180, 0.5, 1, False, " min"),
    ("frag_tox",  "Fragment toxicity vs parent",0,   0.5, 0.01,2, False, ""),
    ("frag_n",    "Fragments per molecule",     1,   3,   1,   0, False, "×"),
    # cell-line robustness phenotype (0 = primary/normal; higher = transformed)
    ("apop_resist",   "Apoptosis resistance",    0, 1, 0.05, 2, False, ""),
    ("anoikis_resist","Anoikis resistance",      0, 1, 0.05, 2, False, ""),
    ("antioxidant",   "Antioxidant capacity",    0, 1, 0.05, 2, False, ""),
    ("glycolytic",    "Glycolytic (Warburg)",    0, 1, 0.05, 2, False, ""),
    ("cyto",          "Cytoskeletal density",    0.3, 2.0, 0.05, 2, False, "×"),
]

COMBOS = [
    ("cpa_key",  "Cryoprotectant (penetrating)", [(k, v["name"]) for k, v in CPAS.items()]),
    ("additive", "Extracellular additive", [(k, v["name"]) for k, v in ADDITIVES.items()]),
    ("adhesion", "Adhesion state", [(k, v[0]) for k, v in ADHESION_STATES.items()]),
    ("add_steps","CPA addition",   [(1, "Single step (bolus)"), (3, "3-step ramp"), (5, "5-step ramp")]),
    ("dilution", "Post-thaw dilution", [("direct", "Direct 1:10 into isotonic"),
                                        ("step", "Stepwise with 0.25 M sucrose"),
                                        ("none", "None — infuse in situ")]),
    ("rock_drug","ROCK inhibitor", [("none", "None"), ("y27632", "Y-27632"),
                                    ("fasudil", "Fasudil"), ("y39983", "Y-39983")]),
    ("rock_when","Inhibitor present in", [("thaw", "Thaw / recovery medium"),
                                          ("freeze", "Freezing medium"), ("both", "Both")]),
]


class Worker(QThread):
    done = pyqtSignal(object)
    fail = pyqtSignal(str)
    def __init__(self, fn, *a, **k):
        super().__init__(); self.fn, self.a, self.k = fn, a, k
    def run(self):
        try: self.done.emit(self.fn(*self.a, **self.k))
        except Exception:
            self.fail.emit(traceback.format_exc())


UC_RATES = [0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 30.0, 100.0]

def _ucurve_job(pdict):
    """Re-run the model across cooling rates for the two-factor survival curve."""
    surv = []
    for cr in UC_RATES:
        d = dict(pdict); d["CR"] = cr; d["ko"] = dict(d.get("ko") or {})
        surv.append(simulate(Params(**d))[1]["S_24"] * 100)
    return UC_RATES, surv


class MechanoView(QWidget):
    """Live mechanotransduction / cell-death signalling network, drawn from the
    model state — membrane tension -> Piezo1 -> Ca; microtubule loss/buckling ->
    RhoA -> ROCK -> p-MLC -> blebbing; focal-adhesion loss -> anoikis; the
    caspase cascade -> apoptosis; with the protective YAP / Akt / IF / dPsi arms."""
    # id, label, x, y (0..1), kind: death|protect|input, state key
    NODES = [
        ("ecm",    "ECM · stiffness",   0.50, 0.06, "input",   None),
        ("tension","Membrane tension",  0.50, 0.16, "death",   "tension"),
        ("piezo",  "Piezo1",            0.17, 0.27, "death",   "piezo"),
        ("ca",     "Ca²⁺ cytosol",      0.17, 0.42, "death",   "caCyt"),
        ("mt",     "Microtubules",      0.50, 0.29, "protect", "mt"),
        ("rock",   "RhoA → ROCK",       0.50, 0.42, "death",   "rock"),
        ("pmlc",   "phospho-MLC",       0.50, 0.55, "death",   "pMLC"),
        ("bleb",   "Blebbing",          0.50, 0.68, "death",   "bleb"),
        ("fa",     "Focal adhesions",   0.84, 0.29, "protect", "FA"),
        ("yap",    "YAP nuclear",       0.84, 0.44, "protect", "yapN"),
        ("akt",    "Akt survival",      0.84, 0.59, "protect", "akt"),
        ("intf",   "Interm. filaments", 0.15, 0.61, "protect", "intf"),
        ("dpsi",   "ΔΨm",               0.33, 0.60, "protect", "dPsi"),
        ("mpt",    "Mito MPT",          0.33, 0.73, "death",   "mpt"),
        ("casp",   "Caspase-3",         0.50, 0.81, "death",   "casp3"),
        ("apop",   "Apoptosis",         0.50, 0.93, "death",   "apop"),
    ]
    EDGES = [
        ("ecm","tension","mech"), ("tension","piezo","act"), ("piezo","ca","act"),
        ("ca","mpt","act"), ("tension","mt","mech"), ("mt","rock","act"),
        ("rock","pmlc","act"), ("pmlc","bleb","act"), ("fa","yap","act"),
        ("yap","akt","act"), ("akt","casp","inh"), ("fa","casp","inh"),
        ("rock","casp","act"), ("ca","casp","act"), ("mpt","casp","act"),
        ("bleb","apop","act"), ("casp","apop","act"), ("intf","bleb","inh"),
        ("mpt","dpsi","inh"),
    ]
    NW, NH = 118, 38

    def __init__(self, parent=None):
        super().__init__(parent)
        self.frame = None
        self.setMinimumSize(420, 460)
        self.setAutoFillBackground(True)

    def set_frame(self, f):
        self.frame = f; self.update()

    def _risk(self, kind, v):
        if kind == "protect":
            return float(max(0.0, min(1.0, 1.0 - v)))
        return float(max(0.0, min(1.0, v / 1.05)))          # death / input

    def _riskcol(self, r, a=255):
        # green (0) -> amber (0.5) -> red (1)
        if r < 0.5:
            t = r / 0.5; R, G, B = 47 + t * 177, 192 - t * 26, 136 - t * 78
        else:
            t = (r - 0.5) / 0.5; R, G, B = 224, 166 - t * 79, 58 - t * 21
        return QColor(int(R), int(G), int(B), a)

    def paintEvent(self, _):
        q = QPainter(self); q.setRenderHint(QPainter.RenderHint.Antialiasing)
        q.fillRect(self.rect(), QColor(C["surface0"]))
        w, hh = self.width(), self.height()
        pad = 14
        def cxy(nx, ny):
            return (pad + nx * (w - 2 * pad), pad + ny * (hh - 2 * pad))
        pos = {nid: cxy(nx, ny) for nid, _l, nx, ny, _k, _key in self.NODES}
        f = self.frame

        if f is None:
            q.setPen(QColor(C["muted"])); q.drawText(self.rect(),
                Qt.AlignmentFlag.AlignCenter, "Run a protocol"); return

        val = {nid: (float(f.get(key, 0.0)) if key else 0.0)
               for nid, _l, _x, _y, _k, key in self.NODES}
        kind = {nid: k for nid, _l, _x, _y, k, _key in self.NODES}

        # ---- edges ----
        for a, b, kind_e in self.EDGES:
            (x0, y0), (x1, y1) = pos[a], pos[b]
            # exit/enter node boxes vertically-ish
            y0 += self.NH / 2 if y1 > y0 else -self.NH / 2
            y1 += -self.NH / 2 if y1 > y0 else self.NH / 2
            active = self._risk(kind[a], val[a])
            base = 60 + int(150 * (active if kind_e != "inh" else 1 - active))
            col = QColor(C["borderStrong"]); col.setAlpha(base)
            pen = QPen(col, 1.4 + 1.6 * active)
            if kind_e == "mech":
                pen.setStyle(Qt.PenStyle.DashLine)
            q.setPen(pen); q.setBrush(Qt.BrushStyle.NoBrush)
            mx, my = (x0 + x1) / 2, (y0 + y1) / 2
            path = QPainterPath(QPointF(x0, y0))
            path.quadTo(QPointF(mx, my), QPointF(x1, y1))
            q.drawPath(path)
            # arrow / bar head
            import math
            ang = math.atan2(y1 - my, x1 - mx)
            if kind_e == "inh":
                bx, by = math.cos(ang + math.pi / 2), math.sin(ang + math.pi / 2)
                q.drawLine(QPointF(x1 - bx * 6, y1 - by * 6), QPointF(x1 + bx * 6, y1 + by * 6))
            else:
                for s in (2.6, -2.6):
                    hx = x1 - math.cos(ang) * 9 + math.cos(ang + math.pi / 2) * s * 2.2
                    hy = y1 - math.sin(ang) * 9 + math.sin(ang + math.pi / 2) * s * 2.2
                    q.drawLine(QPointF(x1, y1), QPointF(hx, hy))

        # ---- nodes ----
        q.setFont(QFont("", 8, QFont.Weight.Bold))
        for nid, label, _x, _y, k, key in self.NODES:
            x, y = pos[nid]
            r = self._risk(k, val[nid])
            rc = self._riskcol(r)
            box = QRectF(x - self.NW / 2, y - self.NH / 2, self.NW, self.NH)
            q.setBrush(QColor(C["surface1"])); q.setPen(QPen(rc, 2.0))
            q.drawRoundedRect(box, 8, 8)
            # activation bar along the bottom edge
            q.setBrush(self._riskcol(r, 210)); q.setPen(Qt.PenStyle.NoPen)
            q.drawRoundedRect(QRectF(x - self.NW / 2 + 3, y + self.NH / 2 - 6,
                                     (self.NW - 6) * (val[nid] if k != "protect" else val[nid]),
                                     3), 2, 2)
            q.setPen(QColor(C["text"])); q.setFont(QFont("", 8, QFont.Weight.Bold))
            q.drawText(QRectF(box.x(), box.y() + 3, self.NW, 14),
                       Qt.AlignmentFlag.AlignHCenter, label)
            q.setPen(rc); q.setFont(QFont("", 9))
            vtxt = "—" if key is None else f"{val[nid]*100:.0f}%"
            q.drawText(QRectF(box.x(), box.y() + 17, self.NW, 14),
                       Qt.AlignmentFlag.AlignHCenter, vtxt)

        # legend
        q.setFont(QFont("", 8)); q.setPen(QColor(C["muted"]))
        q.drawText(QPointF(pad, hh - 6),
                   "green = quiescent/protected · red = active/lost   ·   → promotes  ⊣ inhibits  ⇢ mechanical")


_PHCOL = {"load":"#22d3ee","cool":"#3b82f6","seed":"#8b7ff0","store":"#6b7a99",
          "warm":"#eb6834","melt":"#eb6834","dilute":"#22d3ee","recover":"#1baf7a","end":"#1baf7a"}

class TimelineBar(QWidget):
    """Full-width freeze-thaw timeline: phase bands + temperature curve + scrubber."""
    seek = pyqtSignal(int)
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setFixedHeight(74); self.setAutoFillBackground(True); self.setMouseTracking(True)
        self.T = []; self.phase = []; self.pos = 0; self._drag = False
    def set_data(self, T, phase):
        self.T = list(T); self.phase = list(phase); self.update()
    def set_pos(self, i):
        self.pos = i; self.update()
    def _xf(self, i, w, pad):
        n = max(1, len(self.T) - 1); return pad + (i / n) * (w - 2 * pad)
    def paintEvent(self, _):
        q = QPainter(self); q.setRenderHint(QPainter.RenderHint.Antialiasing)
        q.fillRect(self.rect(), QColor(C["surface2"]))
        if not self.T: return
        w, h, pad, n = self.width(), self.height(), 12, len(self.T)
        for i in range(n - 1):
            x0, x1 = self._xf(i, w, pad), self._xf(i + 1, w, pad)
            c = QColor(_PHCOL.get(self.phase[i], "#888")); c.setAlpha(46)
            q.fillRect(QRectF(x0, 6, x1 - x0, h - 24), c)
        seen = set(); q.setFont(QFont("", 8, QFont.Weight.Bold))
        for i, ph in enumerate(self.phase):
            if ph not in seen and ph != "end":
                seen.add(ph); q.setPen(QColor(_PHCOL.get(ph, "#888")))
                q.drawText(QPointF(self._xf(i, w, pad) + 3, h - 6), ph)
        Tmin, Tmax = -200, 40
        yT = lambda v: 6 + (1 - (v - Tmin) / (Tmax - Tmin)) * (h - 24)
        q.setPen(QPen(QColor("#e0a63a"), 2)); q.setBrush(Qt.BrushStyle.NoBrush)
        path = QPainterPath()
        for i, t in enumerate(self.T):
            x, y = self._xf(i, w, pad), yT(t)
            path.moveTo(x, y) if i == 0 else path.lineTo(x, y)
        q.drawPath(path)
        xs = self._xf(min(self.pos, n - 1), w, pad)
        q.setPen(QPen(QColor(C["nucleus"]), 2)); q.drawLine(QPointF(xs, 3), QPointF(xs, h - 16))
        q.setBrush(QColor(C["nucleus"])); q.setPen(Qt.PenStyle.NoPen)
        q.drawEllipse(QPointF(xs, yT(self.T[min(self.pos, n - 1)])), 4, 4)
    def _seekx(self, x):
        w, pad, n = self.width(), 12, max(1, len(self.T) - 1)
        self.seek.emit(max(0, min(n, int(round((x - pad) / (w - 2 * pad) * n)))))
    def mousePressEvent(self, e): self._drag = True; self._seekx(e.position().x())
    def mouseMoveEvent(self, e):
        if self._drag: self._seekx(e.position().x())
    def mouseReleaseEvent(self, e): self._drag = False


class MolecularView(QWidget):
    """Molecular-detail focal-adhesion & LINC schematic — ECM → integrins →
    talin/vinculin/FAK plaque → actin stress fibre → RhoA-ROCK-LIMK-cofilin →
    nesprin-SUN-lamin A/C → nucleus / mechanosensitive genes. A schematic driven
    by the model's aggregate focal-adhesion, tension and cytoskeletal state
    (per-molecule occupancy is a documented visual proxy, not a modelled species)."""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.frame = None; self.setMinimumSize(360, 520); self.setAutoFillBackground(True)
    def set_frame(self, f): self.frame = f; self.update()

    def paintEvent(self, _):
        q = QPainter(self); q.setRenderHint(QPainter.RenderHint.Antialiasing)
        q.fillRect(self.rect(), QColor(C["surface0"]))
        f = self.frame
        if f is None:
            q.setPen(QColor(C["muted"])); q.drawText(self.rect(),
                Qt.AlignmentFlag.AlignCenter, "Run a protocol"); return
        w, h = self.width(), self.height()
        X = lambda fx: fx * w; Y = lambda fy: fy * h
        FA = f.get("FA", 0); ten = f.get("tension", 0); act = f.get("actin", 0)
        rock = f.get("rock", 0); yap = f.get("yapN", 0); force = f.get("pMLC", 0) * (0.3 + 0.7 * act)
        engaged = clamp(FA / 0.5, 0, 1)              # adhesion engagement (suspension FA0 is low)
        talin_stretch = clamp(ten / 0.5, 0, 1) * (0.3 + 0.7 * engaged)
        good = QColor(47, 192, 136); bad = QColor(224, 87, 79); amb = QColor(224, 166, 58)
        def mix(a, b, t): return QColor(int(a.red()+(b.red()-a.red())*t),
            int(a.green()+(b.green()-a.green())*t), int(a.blue()+(b.blue()-a.blue())*t))
        def blob(cx, cy, r, col, label, sub=None, pP=False):
            q.setBrush(col); q.setPen(QPen(col.darker(140), 1.5))
            q.drawEllipse(QPointF(cx, cy), r, r)
            q.setPen(QColor(C["text"])); q.setFont(QFont("", 8, QFont.Weight.Bold))
            q.drawText(QRectF(cx-40, cy-7, 80, 14), Qt.AlignmentFlag.AlignHCenter, label)
            if sub: q.setFont(QFont("", 7)); q.setPen(QColor(C["muted"]))
            if sub: q.drawText(QRectF(cx-45, cy+r+1, 90, 12), Qt.AlignmentFlag.AlignHCenter, sub)
            if pP:
                q.setBrush(QColor(224,166,58)); q.setPen(Qt.PenStyle.NoPen)
                q.drawEllipse(QPointF(cx+r*0.7, cy-r*0.7), 5, 5)
                q.setPen(QColor("#3a2a00")); q.setFont(QFont("",7,QFont.Weight.Bold))
                q.drawText(QRectF(cx+r*0.7-5, cy-r*0.7-6, 10, 12), Qt.AlignmentFlag.AlignCenter, "P")

        # region labels (right margin) — matching the reference figure
        q.setPen(QColor(C["muted"])); q.setFont(QFont("", 8, QFont.Weight.Bold))
        for fy, lab in [(0.80,"Focal-adhesion\nmechanosensing"),(0.50,"Cytoskeletal\nmechanotransduction"),
                        (0.18,"Nuclear\nmechanotransduction")]:
            for j,ln in enumerate(lab.split("\n")):
                q.drawText(QRectF(w-118, Y(fy)+j*11, 112, 12), Qt.AlignmentFlag.AlignRight, ln)

        # ---- ECM substrate (bottom) ----
        q.setPen(QPen(QColor(C["muted"]), 1))
        for k in range(0, int(w*0.7), 9):
            q.drawLine(QPointF(10+k, h-8), QPointF(2+k, h-2))
        q.setPen(QColor(C["muted"])); q.setFont(QFont("",8))
        q.drawText(QPointF(10, h-14), f"ECM · substrate stiffness")

        # ---- membrane + integrins ----
        my = Y(0.84); q.setPen(QPen(QColor(C["borderStrong"]),1))
        q.setBrush(QColor(C["surface2"])); q.drawRoundedRect(QRectF(X(0.10), my-7, X(0.66), 14), 6, 6)
        for ix in (0.30, 0.52):
            blob(X(ix), my, 8, mix(bad,good,engaged), "ITα", None);
            blob(X(ix)+16, my, 8, mix(bad,good,engaged), "ITβ", None)
            q.setPen(QPen(mix(bad,good,engaged),3)); q.drawLine(QPointF(X(ix)+8, my+8), QPointF(X(ix)+8, h-8))
        # ---- talin (stretched vs folded) ----
        tx = X(0.41); ty0 = my-8; ty1 = Y(0.66)
        q.setPen(QPen(mix(QColor(C["muted"]), good, talin_stretch), 2.4))
        coils = 3 if talin_stretch>0.5 else 6
        path = QPainterPath(QPointF(tx, ty0))
        for c in range(1, coils+1):
            path.lineTo(QPointF(tx + (10 if c%2 else -10), ty0 + (ty1-ty0)*c/coils))
        q.drawPath(path)
        q.setPen(QColor(C["text"])); q.setFont(QFont("",7,QFont.Weight.Bold))
        q.drawText(QRectF(tx-52, (ty0+ty1)/2-6, 40, 12), Qt.AlignmentFlag.AlignRight,
                   "talin" + (" ⟺" if talin_stretch>0.5 else ""))
        q.setFont(QFont("",7)); q.setPen(QColor(C["muted"]))
        q.drawText(QRectF(tx-72, (ty0+ty1)/2+5, 60, 12), Qt.AlignmentFlag.AlignRight,
                   "stretched" if talin_stretch>0.5 else "folded")
        # ---- plaque proteins ----
        py = Y(0.70)
        blob(X(0.55), py, 13, mix(QColor(C["muted"]), good, talin_stretch), "VCL", "vinculin")
        blob(X(0.30), py, 13, mix(bad,good,engaged), "FAK", None, pP=engaged>0.4)
        blob(X(0.30), Y(0.78), 11, mix(QColor(C["muted"]),good,engaged), "PAX", "paxillin")
        blob(X(0.55), Y(0.78), 10, mix(QColor(C["muted"]),good,engaged*0.8), "p130Cas", None)
        if engaged < 0.25:
            q.setPen(bad); q.setFont(QFont("",8,QFont.Weight.Bold))
            q.drawText(QRectF(X(0.20), Y(0.63), X(0.5), 14), Qt.AlignmentFlag.AlignHCenter,
                       "adhesion lost → anoikis")

        # ---- actin stress fibre ----
        ax = X(0.41)
        q.setPen(QPen(QColor(168,152,96, int(120+130*act)), 2+5*act))
        q.drawLine(QPointF(ax, ty1), QPointF(ax, Y(0.34)))
        for m in range(4):  # myosin
            yy = Y(0.66) - (Y(0.66)-Y(0.34))*(m+.5)/4
            q.setBrush(mix(good,bad,clamp(f.get("pMLC",0),0,1))); q.setPen(Qt.PenStyle.NoPen)
            q.drawEllipse(QPointF(ax, yy), 3.5, 2)
        q.setPen(QColor(C["muted"])); q.setFont(QFont("",7))
        q.drawText(QPointF(ax+8, Y(0.50)), "actin + myosin II")

        # ---- RhoA→ROCK→LIMK→cofilin (left cascade) ----
        casc=[("RhoA",0.44),("ROCK",0.40),("LIMK",0.44),("CFL-P",0.40)]
        for i,(nm,_)in enumerate([("RhoA",0),("ROCK",0),("LIMK",0),("cofilin",0)]):
            cyy=Y(0.62-i*0.075); c=mix(good,bad,clamp(rock,0,1)) if nm!="cofilin" else mix(bad,good,clamp(rock,0,1))
            blob(X(0.13), cyy, 12, c, nm, "P" if nm=="cofilin" else None)
            if i>0: q.setPen(QPen(QColor(C["borderStrong"]),1.4));q.drawLine(QPointF(X(0.13),Y(0.62-(i-1)*0.075)+12),QPointF(X(0.13),cyy-12))
        q.setPen(QPen(QColor(C["borderStrong"]),1.4,Qt.PenStyle.DashLine))
        q.drawLine(QPointF(X(0.13),Y(0.335)+2),QPointF(ax-6,Y(0.36)))  # ROCK→actin stabilisation

        # ---- nuclear envelope: nesprin-SUN-lamin ----
        ny=Y(0.30)
        q.setBrush(QColor(C["surface2"]));q.setPen(QPen(QColor(C["borderStrong"]),1))
        q.drawRoundedRect(QRectF(X(0.10),ny-6,X(0.66),12),5,5)  # ONM
        q.drawRoundedRect(QRectF(X(0.10),ny+8,X(0.66),12),5,5)  # INM
        # nesprin (outer) + SUN (inner) linking actin to lamina
        q.setPen(QPen(mix(QColor(C["muted"]),good,act),2.2))
        q.drawLine(QPointF(ax,Y(0.34)),QPointF(ax,ny-6))       # actin→nesprin
        blob(ax, ny, 8, mix(QColor(C["muted"]),good,act),"nesprin",None)
        blob(ax, ny+14, 8, mix(QColor(C["muted"]),good,act),"SUN",None)
        q.setPen(QColor(C["muted"]));q.setFont(QFont("",7))
        q.drawText(QPointF(X(0.60),ny-8),"outer NM");q.drawText(QPointF(X(0.60),ny+26),"inner NM")
        # lamin A/C mesh (loads amber under force)
        lam=mix(QColor(139,127,240),amb,clamp(force,0,1))
        q.setPen(QPen(lam,2));q.drawLine(QPointF(X(0.12),ny+22),QPointF(X(0.74),ny+22))
        for k in range(6):q.drawLine(QPointF(X(0.14+k*0.11),ny+20),QPointF(X(0.18+k*0.11),ny+26))
        q.setPen(lam);q.setFont(QFont("",7,QFont.Weight.Bold));q.drawText(QPointF(X(0.12),ny+38),"lamin A/C + emerin")

        # ---- nucleus / mechanosensitive genes ----
        q.setBrush(QColor(42,120,214,40));q.setPen(QPen(QColor(C["nucleus"]),1.5))
        q.drawEllipse(QPointF(X(0.42),Y(0.10)),X(0.30),Y(0.09))
        mech=mix(bad,good,clamp(yap,0,1))
        q.setBrush(mech);q.setPen(Qt.PenStyle.NoPen)
        q.drawEllipse(QPointF(X(0.42),Y(0.10)),6,6)
        q.setPen(QColor(C["text"]));q.setFont(QFont("",8,QFont.Weight.Bold))
        q.drawText(QRectF(X(0.20),Y(0.10)+8,X(0.44),14),Qt.AlignmentFlag.AlignHCenter,
                   "mechanosensitive genes")
        q.setFont(QFont("",7));q.setPen(QColor(C["muted"]))
        q.drawText(QRectF(X(0.20),Y(0.10)+21,X(0.44),12),Qt.AlignmentFlag.AlignHCenter,
                   f"YAP/TAZ nuclear {yap*100:.0f}%")


class Main(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("CryoCell — a virtual cell for cryoprotectant testing")
        self.resize(1680, 1000)
        self.P = Params()
        self.S = self.R = None
        self.idx = 0
        self.playing = False
        self._workers = []

        split = QSplitter(Qt.Orientation.Horizontal)
        split.addWidget(self._controls())
        split.addWidget(self._centre())
        split.addWidget(self._right())
        split.setSizes([330, 800, 430])
        # AIDO-style shell: three panels above a full-width freeze-thaw timeline
        self.timelinebar = TimelineBar()
        self.timelinebar.seek.connect(self._seek_frame)
        central = QWidget(); cvl = QVBoxLayout(central)
        cvl.setContentsMargins(0, 0, 0, 0); cvl.setSpacing(0)
        cvl.addWidget(split, 1); cvl.addWidget(self.timelinebar)
        self.setCentralWidget(central)

        self.play_timer = QTimer(self); self.play_timer.timeout.connect(self._advance)
        self.run_timer = QTimer(self); self.run_timer.setSingleShot(True)
        self.run_timer.timeout.connect(self.run)
        self.statusBar().showMessage("Ready")
        self.run()

    # ------------------------------------------------------------- controls
    def _controls(self):
        w = QWidget(); v = QVBoxLayout(w); v.setContentsMargins(8, 8, 8, 8); v.setSpacing(6)
        top = QHBoxLayout()
        self.preset = QComboBox(); self.preset.addItem("Load a preset…")
        for k in PRESETS: self.preset.addItem(k)
        self.preset.currentTextChanged.connect(self._preset)
        top.addWidget(self.preset)
        b = QPushButton("Run"); b.clicked.connect(self.run); top.addWidget(b)
        v.addLayout(top)

        self.widgets = {}
        box = QWidget(); g = QVBoxLayout(box); g.setSpacing(4)
        for attr, lab, options in COMBOS:
            g.addWidget(QLabel(lab))
            cb = QComboBox()
            for val, txt in options: cb.addItem(txt, val)
            cb.setCurrentIndex(max(0, [o[0] for o in options].index(getattr(self.P, attr))
                                   if getattr(self.P, attr) in [o[0] for o in options] else 0))
            cb.currentIndexChanged.connect(lambda _i, a=attr, c=cb: self._set(a, c.currentData()))
            self.widgets[attr] = cb; g.addWidget(cb)
        for attr, lab, lo, hi, stp, dec, log, unit in SLIDERS:
            row = QHBoxLayout()
            l1 = QLabel(lab); l1.setStyleSheet("color:#52514e")
            l2 = QLabel(); l2.setAlignment(Qt.AlignmentFlag.AlignRight)
            l2.setStyleSheet("font-family:monospace;color:#0b0b0b")
            row.addWidget(l1); row.addStretch(); row.addWidget(l2)
            g.addLayout(row)
            s = QSlider(Qt.Orientation.Horizontal)
            n = int(round((hi - lo) / stp))
            s.setRange(0, n)
            cur = getattr(self.P, attr)
            raw = math.log10(cur) if log else cur
            s.setValue(int(round((float(np.clip(raw, lo, hi)) - lo) / stp)))
            s.valueChanged.connect(lambda val, a=attr, L=lo, S=stp, lg=log, d=dec,
                                          u=unit, out=l2: self._slide(a, L + val * S, lg, d, u, out))
            self.widgets[attr] = s
            self._slide(attr, lo + s.value() * stp, log, dec, unit, l2, apply=False)
            g.addWidget(s)
        g.addStretch()
        sc = QScrollArea(); sc.setWidget(box); sc.setWidgetResizable(True)
        sc.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        v.addWidget(sc, 1)
        return w

    def _slide(self, attr, raw, log, dec, unit, out, apply=True):
        val = 10 ** raw if log else raw
        if attr in ("add_steps",): val = int(val)
        out.setText(f"{val:.{dec}f}{unit}")
        if apply: self._set(attr, val)

    def _set(self, attr, val):
        setattr(self.P, attr, val)
        self.run_timer.start(220)

    def _preset(self, name):
        p = PRESETS.get(name)
        if not p: return
        self.P = Params(**{**{k: v for k, v in vars(self.P).items() if k in Params.__annotations__}, **p})
        for attr, w in self.widgets.items():
            val = getattr(self.P, attr, None)
            if val is None: continue
            if isinstance(w, QComboBox):
                i = w.findData(val)
                if i >= 0: w.blockSignals(True); w.setCurrentIndex(i); w.blockSignals(False)
            else:
                spec = next(s for s in SLIDERS if s[0] == attr)
                _a, _l, lo, hi, stp, dec, log, unit = spec
                raw = math.log10(max(val, 1e-6)) if log else val
                w.blockSignals(True)
                w.setValue(int(round((float(np.clip(raw, lo, hi)) - lo) / stp)))
                w.blockSignals(False)
        self.run()

    # --------------------------------------------------------------- centre
    def _centre(self):
        w = QWidget(); v = QVBoxLayout(w); v.setContentsMargins(6, 6, 6, 6); v.setSpacing(6)
        # experiment thread — the protocol as an ordered, clickable phase sequence
        thr = QHBoxLayout(); thr.setSpacing(3)
        thr.addWidget(QLabel("Experiment:"))
        self.threadbtns = {}
        for ph, lab in [("load", "Load"), ("cool", "Cool"), ("seed", "Seed"),
                        ("store", "Store"), ("warm", "Warm"), ("dilute", "Dilute"),
                        ("recover", "Recover")]:
            b = QPushButton(lab); b.setCheckable(True)
            b.setStyleSheet("QPushButton{font-size:11px;padding:3px 9px;border:1px solid #cfd8e6;"
                            "border-radius:12px;background:palette(base)}"
                            "QPushButton:checked{background:#2a78d6;color:white;border-color:#2a78d6}")
            b.clicked.connect(lambda _c, p=ph: self._jump_phase(p))
            self.threadbtns[ph] = b; thr.addWidget(b)
        thr.addStretch(); v.addLayout(thr)
        self.view = CellView()
        self.view.hovered.connect(self._hover)
        self.view.selected.connect(self._select)
        v.addWidget(self.view, 1)

        row = QHBoxLayout()
        self.playbtn = QPushButton("▶"); self.playbtn.setFixedWidth(38)
        self.playbtn.clicked.connect(self._play); row.addWidget(self.playbtn)
        self.scrub = QSlider(Qt.Orientation.Horizontal)
        self.scrub.valueChanged.connect(self._scrub); row.addWidget(self.scrub, 1)
        self.phase = QLabel(); self.phase.setMinimumWidth(210)
        self.phase.setStyleSheet("font-family:monospace;color:#52514e")
        self.phase.setAlignment(Qt.AlignmentFlag.AlignRight); row.addWidget(self.phase)
        self.rendcombo = QComboBox()
        for key, lbl in [("illustrative", "Illustrative"), ("fluor", "Fluorescence"),
                         ("phase", "Phase-contrast")]:
            self.rendcombo.addItem(lbl, key)
        self.rendcombo.setToolTip("Rendering style: illustrative colour, confocal-fluorescence "
                                  "(dark field, glowing channels), or label-free phase-contrast greyscale")
        self.rendcombo.currentIndexChanged.connect(
            lambda _i: self.view.set_render_mode(self.rendcombo.currentData()))
        row.addWidget(self.rendcombo)
        self.scibox = QCheckBox("Scientific labels")
        self.scibox.setToolTip("Overlay scale bar, leader-line labels with live quantitative "
                               "state, legend, and a cryo-stage dendritic ice front (any rendering)")
        self.scibox.toggled.connect(self.view.set_sci); row.addWidget(self.scibox)
        self.savebtn = QPushButton("Save video…")
        self.savebtn.setToolTip("Render the cell-view animation across the whole protocol "
                                "to an MP4 or animated GIF")
        self.savebtn.clicked.connect(self._export_animation); row.addWidget(self.savebtn)
        v.addLayout(row)

        lay = QGroupBox("Compartments"); lg = QGridLayout(lay); lg.setSpacing(2)
        self.layer_boxes = {}
        for i, (key, label, hpa_key) in enumerate(LAYERS):
            cb = QCheckBox(label); cb.setChecked(True)
            cb.setStyleSheet("font-size:11px")
            if hpa_key and HPA.get(hpa_key) and HPA[hpa_key][1]:
                cb.setToolTip(summary_line(hpa_key))
            cb.toggled.connect(lambda on, k=key: self.view.set_layer(k, on))
            self.layer_boxes[key] = cb
            lg.addWidget(cb, i // 4, i % 4)
        v.addWidget(lay)
        hint = QLabel("scroll to zoom · shift-drag to pan · drag the cell or an organelle to deform it · "
                      "hover for live state and Human Protein Atlas data")
        hint.setStyleSheet("color:#7a7873;font-size:11px"); hint.setWordWrap(True)
        v.addWidget(hint)
        return w

    # ---------------------------------------------------------------- right
    def _right(self):
        tabs = QTabWidget()
        # --- outcome + inspector
        w = QWidget(); v = QVBoxLayout(w)
        self.hero = QLabel("—"); self.hero.setFont(QFont("", 34, QFont.Weight.Bold))
        v.addWidget(self.hero)
        self.hero2 = QLabel(); self.hero2.setStyleSheet("color:#52514e"); self.hero2.setWordWrap(True)
        v.addWidget(self.hero2)
        self.quad = QLabel(); self.quad.setStyleSheet("font-family:monospace"); v.addWidget(self.quad)
        self.dmg = QTableWidget(0, 2); self.dmg.horizontalHeader().setVisible(False)
        self.dmg.verticalHeader().setVisible(False)
        self.dmg.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        self.dmg.setMaximumHeight(210); v.addWidget(self.dmg)
        insp_l = QLabel("Inspector"); insp_l.setStyleSheet("font-weight:600;color:#52514e")
        v.addWidget(insp_l)
        self.insp = QTextEdit(); self.insp.setReadOnly(True); v.addWidget(self.insp, 1)
        tabs.addTab(w, "Outcome")

        # --- plots
        pw = QWidget(); pv = QVBoxLayout(pw)
        self.plots, self.curves = [], []
        specs = [("Compartment volumes", [("Whole cell","Vn"),("Mito matrix","Vmito"),("Nucleus","Vnuc")]),
                 ("Solute load", [("CPA in (M)","Cin"),("Osm out","Osme"),("mol% at bilayer","x")]),
                 ("Ice & mitochondria", [("P(IIF)","Piif"),("ΔΨm","dPsi"),("Ca²⁺ cyt","caCyt")]),
                 ("Membrane", [("Gel","gel"),("Defects","pore"),("Tension","tension")]),
                 ("Mechanotransduction", [("ROCK","rock"),("Blebbing","bleb"),("Caspase-3","casp3")]),
                 ("Death", [("Focal adhesions","FA"),("Nuclear YAP","yapN"),
                            ("Apoptotic","apop"),("Necrotic","necr")])]
        for title, series in specs:
            p = pg.PlotWidget(title=title); p.showGrid(x=True, y=True, alpha=0.15)
            p.addLegend(offset=(-8, 8), labelTextSize="7pt")
            cs = [p.plot([], [], pen=pg.mkPen(SER[i % len(SER)], width=2), name=nm)
                  for i, (nm, _k) in enumerate(series)]
            p.setMinimumHeight(150)
            self.vline = None
            self.plots.append((p, series)); self.curves.append(cs)
            pv.addWidget(p)
        sc = QScrollArea(); sc.setWidget(pw); sc.setWidgetResizable(True)
        tabs.addTab(sc, "Plots")

        # --- cryo-signatures: the three classic quantitative cryobiology plots
        gw = QWidget(); gv = QVBoxLayout(gw)
        self.sig_uc = pg.PlotWidget(title="Two-factor survival vs cooling rate (Mazur inverted-U)")
        self.sig_uc.showGrid(x=True, y=True, alpha=0.15); self.sig_uc.setLogMode(x=True, y=False)
        self.sig_uc.setLabel("bottom", "cooling rate", "°C/min"); self.sig_uc.setLabel("left", "24 h survival", "%")
        self.sig_uc.setMinimumHeight(190); gv.addWidget(self.sig_uc)
        self.sig_uc_curve = self.sig_uc.plot([], [], pen=pg.mkPen(SER[0], width=2.5))
        self.sig_uc_pt = self.sig_uc.plot([], [], pen=None, symbol='o', symbolBrush=SER[1], symbolSize=12)

        self.sig_pd = pg.PlotWidget(title="Phase diagram — extracellular solution on the liquidus")
        self.sig_pd.showGrid(x=True, y=True, alpha=0.15)
        self.sig_pd.setLabel("bottom", "temperature", "°C"); self.sig_pd.setLabel("left", "extracellular osmolality", "osmol/L")
        self.sig_pd.setMinimumHeight(190); gv.addWidget(self.sig_pd)
        self.sig_pd_liq = self.sig_pd.plot([], [], pen=pg.mkPen(SER[3], width=2, style=Qt.PenStyle.DashLine), name="liquidus")
        self.sig_pd_traj = self.sig_pd.plot([], [], pen=pg.mkPen(SER[0], width=2))
        self.sig_pd_pt = self.sig_pd.plot([], [], pen=None, symbol='o', symbolBrush=SER[1], symbolSize=11)

        self.sig_vt = pg.PlotWidget(title="Cell volume vs temperature (osmotic dehydration response)")
        self.sig_vt.showGrid(x=True, y=True, alpha=0.15)
        self.sig_vt.setLabel("bottom", "temperature", "°C"); self.sig_vt.setLabel("left", "V / Viso", "")
        self.sig_vt.setMinimumHeight(190); gv.addWidget(self.sig_vt)
        self.sig_vt_curve = self.sig_vt.plot([], [], pen=pg.mkPen(SER[0], width=2))
        self.sig_vt_pt = self.sig_vt.plot([], [], pen=None, symbol='o', symbolBrush=SER[1], symbolSize=11)
        note = QLabel("The U-curve is computed by re-running the model across cooling rates; "
                      "the dot marks your protocol. The phase diagram and V–T trace the current run.")
        note.setStyleSheet("color:#7a7873;font-size:11px"); note.setWordWrap(True); gv.addWidget(note)
        gsc = QScrollArea(); gsc.setWidget(gw); gsc.setWidgetResizable(True)
        tabs.addTab(gsc, "Cryo-signatures")

        # --- mechanotransduction / cell-death signalling network (live)
        self.mechano = MechanoView()
        msc = QScrollArea(); msc.setWidget(self.mechano); msc.setWidgetResizable(True)
        tabs.addTab(msc, "Signalling")

        # --- molecular-detail focal-adhesion & LINC inset (live)
        self.molec = MolecularView()
        molsc = QScrollArea(); molsc.setWidget(self.molec); molsc.setWidgetResizable(True)
        tabs.addTab(molsc, "FA · LINC")

        # --- analysis
        aw = QWidget(); av = QVBoxLayout(aw)
        row = QHBoxLayout()
        for lab, fn in [("Population", self._pop), ("Perturbation screen", self._screen),
                        ("Measure next?", self._next), ("Stability", self._stab)]:
            b = QPushButton(lab); b.clicked.connect(fn); row.addWidget(b)
        av.addLayout(row)
        self.prog = QProgressBar(); self.prog.setRange(0, 0); self.prog.hide(); av.addWidget(self.prog)
        self.pop_plot = pg.PlotWidget(title="Per-cell 24 h survival")
        self.pop_plot.setMinimumHeight(190); av.addWidget(self.pop_plot)
        self.anal = QTextEdit(); self.anal.setReadOnly(True); av.addWidget(self.anal, 1)
        tabs.addTab(aw, "Analysis")

        # --- HPA
        hw = QWidget(); hv = QVBoxLayout(hw)
        t = QTextEdit(); t.setReadOnly(True)
        rows = ["<b>Human Protein Atlas — Cell Atlas</b>",
                f"{HPA_TOTALS['genes_experimental']:,} genes ({HPA_TOTALS['pct_experimental']}% of "
                f"protein-coding) localised experimentally; {HPA_TOTALS['genes_with_prediction']:,} "
                f"({HPA_TOTALS['pct_with_prediction']}%) including predictions, across "
                f"{HPA_TOTALS['locations_annotated']} annotated subcellular locations.<br>",
                "<table cellpadding=4><tr><th align=left>Compartment</th><th align=right>Genes</th>"
                "<th align=right>% of proteome</th><th align=right>also elsewhere</th></tr>"]
        for k, (nm, genes, pct, multi, note) in HPA.items():
            if genes is None: continue
            rows.append(f"<tr><td>{nm}</td><td align=right>{genes:,}</td>"
                        f"<td align=right>{pct}%</td><td align=right>{multi*100:.0f}%</td></tr>")
        rows.append("</table><br><i>Volume fractions used for the drawing are conventional "
                    "mammalian values, not HPA data — the HPA counts proteins, not volume.</i>")
        t.setHtml("\n".join(rows)); hv.addWidget(t)
        tabs.addTab(hw, "HPA reference")
        return tabs

    # ------------------------------------------------------------------ run
    def run(self):
        try:
            self.S, self.R, self.E = simulate(self.P)
        except Exception:
            self.statusBar().showMessage("simulation failed"); traceback.print_exc(); return
        n = len(self.S)
        self.scrub.blockSignals(True); self.scrub.setRange(0, n - 1)
        self.idx = min(self.idx, n - 1); self.scrub.setValue(self.idx)
        self.scrub.blockSignals(False)
        R = self.R
        self.hero.setText(f"{R['S_24']*100:.0f}%")
        colr = C["good"] if R['S_24'] > .7 else C["warn"] if R['S_24'] > .4 else C["crit"]
        self.hero.setStyleSheet(f"color:{colr}")
        self.hero2.setText(f"viable at {self.P.recover_h:.0f} h post-thaw · membrane-integrity "
                           f"{R['S_imm']*100:.0f}% · settles to {R['S_72']*100:.0f}% by 72 h "
                           f"(delayed-onset death) · functional recovery {R['F_rec']*100:.0f}%")
        self.quad.setText(f"live {R['liveFrac']*100:.0f}%   apoptotic {R['apopFrac']*100:.0f}%   "
                          f"necrotic {R['necrFrac']*100:.0f}%")
        rows = [("Intracellular ice (+ propagation)", R['D_iif']), ("Osmotic / solution effects", R['D_osm']),
                ("Membrane (MD pore + phase)", R['D_mem']), ("CPA toxicity", R['D_tox']),
                ("Mechanical: ice squeeze", R['D_mech']), ("Regulated death (0-24 h)", R['D_apop']),
                ("Delayed-onset death (ROS, 1-3 d)", R['D_delayed'])]
        self.dmg.setRowCount(len(rows))
        for i, (lab, val) in enumerate(rows):
            self.dmg.setItem(i, 0, QTableWidgetItem(lab))
            it = QTableWidgetItem(f"{val*100:.0f}%")
            it.setTextAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
            it.setForeground(QColor(C["crit"] if val > .3 else C["serious"] if val > .1 else C["muted"]))
            self.dmg.setItem(i, 1, it)
        self.dmg.setColumnWidth(1, 60)
        x = np.arange(n)
        for (p, series), cs in zip(self.plots, self.curves):
            for c, (_nm, key) in zip(cs, series):
                y = np.array(getattr(self.S, key), dtype=float)
                if key == "x": y = y / 30.0
                if key == "Osme": y = y / 25.0
                if key == "Cin": y = y / 14.0
                c.setData(x, y)
        self._update_signatures()
        self._compute_ucurve()
        self.timelinebar.set_data(self.S.T, self.S.phase)
        self._show(self.idx, jumped=True)
        self.statusBar().showMessage(
            f"{n} frames · peak ROCK {R['rockPeak']:.2f} · grain {R['grainMax']:.0f} µm · "
            f"min channel/cell {R['squeezeMin']:.2f} · ATP min {R['atpMin']:.2f}")

    # ------------------------------------------------------- cryo-signatures
    def _update_signatures(self):
        if not self.S: return
        S = self.S
        T = np.array(S.T, dtype=float); Osme = np.array(S.Osme, dtype=float)
        Vn = np.array(S.Vn, dtype=float)
        tt = np.linspace(-50, -0.5, 120)
        self.sig_pd_liq.setData(tt, [osm_from_dT(-x, 0.033) for x in tt])
        m = T < 0.5
        self.sig_pd_traj.setData(T[m], Osme[m])
        self.sig_vt_curve.setData(T, Vn)
        self._update_sig_point()

    def _update_sig_point(self):
        if not self.S: return
        i = min(self.idx, len(self.S) - 1)
        self.sig_pd_pt.setData([self.S.T[i]], [self.S.Osme[i]])
        self.sig_vt_pt.setData([self.S.T[i]], [self.S.Vn[i]])

    def _compute_ucurve(self):
        pdict = {k: v for k, v in vars(self.P).items()}
        pdict["ko"] = dict(pdict.get("ko") or {})
        w = Worker(_ucurve_job, pdict); self._workers.append(w)
        w.done.connect(self._ucurve_done)
        w.fail.connect(lambda e: None)
        w.start()

    def _ucurve_done(self, res):
        rates, surv = res
        self.sig_uc_curve.setData(list(rates), list(surv))
        if self.R:
            self.sig_uc_pt.setData([self.P.CR], [self.R["S_24"] * 100])

    # ---------------------------------------------------------------- frames
    def _frame_dict(self, i):
        S, R = self.S, self.R
        ph = S.phase[i]
        frozen = S.fIce[i] > 0.02 and S.T[i] < 0 and ph not in ("dilute", "recover", "end")
        return dict(i=i, T=S.T[i], Vn=S.Vn[i], Vmito=S.Vmito[i], Vnuc=S.Vnuc[i],
                    Cin=S.Cin[i], Osme=S.Osme[i], x=S.x[i], pore=S.pore[i],
                    Dmem=S.Dmem[i], fluid=S.fluid[i], glass=S.glass[i],
                    mcpa=S.mcpa[i], Pmito=S.Pmito[i], prot=S.prot[i],
                    gel=S.gel[i], APL=S.APL[i], thick=S.thick[i], tension=S.tension[i],
                    mt=S.mt[i], actin=S.actin[i], pMLC=S.pMLC[i], bleb=S.bleb[i],
                    intf=S.intf[i], sigMT=S.sigMT[i], sigIF=S.sigIF[i],
                    Piif=S.Piif[i], dPsi=S.dPsi[i], caCyt=S.caCyt[i], caER=S.caER[i],
                    FA=S.FA[i], yapN=S.yapN[i], apop=S.apop[i], necr=S.necr[i],
                    piezo=S.piezo[i], akt=S.akt[i],
                    atp=S.atp[i], rock=S.rock[i], casp3=S.casp3[i], mpt=S.mpt[i],
                    fIce=S.fIce[i], grain=S.grain[i], chanW=S.chanW[i],
                    squeeze=S.squeeze[i], frozen=frozen, phase=ph,
                    sterol=self.P.sterol, cpaLoad=self.P.molar(), cyto=self.P.cyto,
                    r_iso_um=(3 * self.P.Viso / (4 * math.pi)) ** (1 / 3))

    def _show(self, i, jumped=False):
        f = self._frame_dict(i)
        self.view.set_frame(f, jumped)
        self.mechano.set_frame(f)
        self.molec.set_frame(f)
        names = dict(load="CPA loading", cool="Cooling", seed="Seeding", store="Storage",
                     warm="Warming", melt="Melting", dilute="Dilution", recover="Recovery", end="End")
        tt = self.S.t[i]
        tstr = f"{tt/60:.1f} min" if tt < 7200 else f"{tt/3600:.1f} h"
        self.phase.setText(f"{names.get(f['phase'], f['phase'])} · {f['T']:.1f} °C · {tstr}")
        self._update_sig_point()
        self.timelinebar.set_pos(i)
        cur = f['phase']
        for ph, b in self.threadbtns.items():
            b.blockSignals(True); b.setChecked(ph == cur); b.blockSignals(False)
        if self.view.sel or self.view.hover:
            self._hover(self.view.hover or self.view.sel)

    def _scrub(self, v):
        jumped = abs(v - self.idx) > 2
        self.idx = v; self._show(v, jumped)

    def _seek_frame(self, i):
        i = max(0, min(i, len(self.S) - 1)) if self.S else 0
        self.scrub.blockSignals(True); self.scrub.setValue(i); self.scrub.blockSignals(False)
        self.idx = i; self._show(i, jumped=True)

    def _jump_phase(self, ph):
        if not self.S: return
        for i, p in enumerate(self.S.phase):
            if p == ph:
                self._seek_frame(i); return

    def _play(self):
        self.playing = not self.playing
        self.playbtn.setText("❚❚" if self.playing else "▶")
        self.play_timer.start(40) if self.playing else self.play_timer.stop()

    def _advance(self):
        if not self.S: return
        self.idx = (self.idx + 2) % len(self.S)
        self.scrub.blockSignals(True); self.scrub.setValue(self.idx); self.scrub.blockSignals(False)
        self._show(self.idx)

    # ------------------------------------------------------- animation export
    def _export_animation(self, _checked=False, path=None, W=900, H=760, fps=25):
        """Render the cell view across the whole trajectory to an MP4 or GIF.

        Uses ffmpeg (piped raw frames) when it is on PATH; otherwise falls back
        to writing a numbered PNG sequence next to the chosen file. Rendering is
        done on a private off-screen CellView so the live view is undisturbed.
        The `path` argument lets this be driven headlessly in tests.
        """
        if not self.S:
            return
        if path is None:
            default = os.path.join(os.path.expanduser("~"), "cryocell.mp4")
            path, _ = QFileDialog.getSaveFileName(
                self, "Save cell-view animation", default,
                "MP4 video (*.mp4);;Animated GIF (*.gif)")
            if not path:
                return
        ext = os.path.splitext(path)[1].lower()
        if ext not in (".mp4", ".gif"):
            path += ".mp4"; ext = ".mp4"
        ffmpeg = shutil.which("ffmpeg")

        n = len(self.S)
        step = max(1, n // 600)                 # cap the number of rendered frames
        idxs = list(range(0, n, step))

        # private off-screen renderer, seeded to the current view settings
        ev = CellView(); ev.timer.stop(); ev.resize(W, H)
        ev.visible = set(self.view.visible)
        ev.zoom = self.view.zoom; ev.pan = QPointF(0, 0)
        ev.set_frame(self._frame_dict(idxs[0]), jumped=True)
        for _ in range(60):                     # warm up the soft body from rest
            ev._tick()

        self.savebtn.setEnabled(False)
        self._busy(True, "Rendering animation…")
        self.prog.show(); self.prog.setRange(0, len(idxs)); self.prog.setValue(0)

        proc = None; framedir = None; msg = ""
        try:
            if ffmpeg:
                base = [ffmpeg, "-y", "-loglevel", "error", "-f", "rawvideo",
                        "-pix_fmt", "rgba", "-s", f"{W}x{H}", "-r", str(fps), "-i", "-"]
                if ext == ".mp4":
                    cmd = base + ["-c:v", "libx264", "-pix_fmt", "yuv420p",
                                  "-crf", "18", "-movflags", "+faststart", path]
                else:
                    cmd = base + ["-vf", "split[a][b];[a]palettegen[p];[b][p]paletteuse", path]
                proc = subprocess.Popen(cmd, stdin=subprocess.PIPE,
                                        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            else:
                framedir = os.path.splitext(path)[0] + "_frames"
                os.makedirs(framedir, exist_ok=True)

            for k, i in enumerate(idxs):
                ev.set_frame(self._frame_dict(i))
                for _ in range(3):
                    ev._tick()
                img = QImage(W, H, QImage.Format.Format_RGBA8888)
                p = QPainter(img); ev.render(p); p.end()
                if proc:
                    buf = img.constBits(); buf.setsize(img.sizeInBytes())
                    proc.stdin.write(bytes(buf))
                else:
                    img.save(os.path.join(framedir, f"frame_{k:04d}.png"))
                self.prog.setValue(k + 1)
                if k % 4 == 0:
                    QApplication.processEvents()

            if proc:
                proc.stdin.close(); proc.wait()
                secs = len(idxs) / fps
                msg = f"Saved {os.path.basename(path)} — {len(idxs)} frames, {secs:.0f}s"
            else:
                msg = (f"ffmpeg not found — wrote {len(idxs)} PNG frames to "
                       f"{os.path.basename(framedir)}/ (assemble with any encoder)")
        except Exception:
            traceback.print_exc(); msg = "Animation export failed — see console"
            if proc and proc.stdin:
                try: proc.stdin.close(); proc.wait()
                except Exception: pass
        finally:
            self.prog.hide(); self._busy(False); self.savebtn.setEnabled(True)
        self.statusBar().showMessage(msg, 9000)
        return path

    # ------------------------------------------------------------- inspector
    def _select(self, key): self._hover(key)

    def _hover(self, key):
        if not self.S: return
        f = self._frame_dict(self.idx); R = self.R
        HK = {"nucleus": "nucleoplasm", "nucleoli": "nucleoli", "nuclear_mem": "nuclear_mem",
              "mitochondria": "mitochondria", "er": "er", "golgi": "golgi",
              "lysosomes": "lysosomes", "peroxisomes": "peroxisomes",
              "endosomes": "endosomes", "lipid_drop": "lipid_drop",
              "vesicles": "vesicles", "membrane": "plasma_mem", "cytosol": "cytosol",
              "centrosome": "centrosome", "microtubules": "microtubules",
              "cortex": "actin", "interm_fil": "interm_fil"}
        h = []
        title = key.replace("_", " ").title() if key else "Extracellular ice field"
        h.append(f"<b style='font-size:13px'>{title}</b>")
        hk = HK.get(key)
        if hk and HPA.get(hk):
            nm, genes, pct, multi, note = HPA[hk]
            if genes:
                h.append(f"<span style='color:#52514e'>Human Protein Atlas — <b>{genes:,}</b> genes "
                         f"({pct}% of protein-coding), <b>{multi*100:.0f}%</b> also localise "
                         f"elsewhere.<br>{note}</span>")
            else:
                h.append(f"<span style='color:#52514e'>{note}</span>")
        h.append("<table cellpadding=3>")
        def r(k, v): h.append(f"<tr><td style='color:#52514e'>{k}</td>"
                              f"<td align=right style='font-family:monospace'>{v}</td></tr>")
        if key == "mitochondria":
            r("Matrix volume", f"{f['Vmito']*100:.0f}% of resting")
            r("Membrane potential ΔΨm", f"{f['dPsi']*100:.0f}%")
            r("Permeability transition", f"{f['mpt']*100:.1f}% open")
            r("Cristae", "unfolded" if f['Vmito'] > 1.6 else "intact")
            r("Osmometer", "aquaporin-8 fast water exchange; matrix ~ cytosol")
            r("Matrix CPA", f"{f.get('mcpa',0):.2f} M (loads across inner membrane)")
            r("Intramitochondrial ice", f"{f.get('Pmito',0)*100:.0f}% (cardiolipin nucleator)")
            if f.get('Pmito', 0) > 0.02:
                r("→ seeds cytosolic ice + MPT", f"{R['D_mitoice']*100:.0f}% lethal")
        elif key in ("nucleus", "nucleoli", "nuclear_mem"):
            r("Nuclear volume", f"{f['Vnuc']*100:.0f}% of resting")
            r("Chromatin condensation", f"{np.clip((f['Osme']-1)/10,0,1)*100:.0f}%")
            r("Nuclear YAP", f"{f['yapN']*100:.0f}%")
            r("Envelope", "pores make it leaky to small solutes")
        elif key == "er":
            r("Lumenal Ca²⁺ remaining", f"{f['caER']*100:.0f}%")
            r("Cytosolic Ca²⁺", f"{f['caCyt']:.2f} (1 = overload)")
            r("SERCA", "active" if f['T'] > 10 else "cold-arrested")
        elif key == "membrane" or key == "cortex":
            r("CPA at the bilayer", f"{f['x']:.1f} mol%")
            r("Reflection coeff σ", f"{CPAS[self.P.cpa_key]['sigma']:.2f} (KK 3-parameter)")
            r("Area per lipid", f"{f['APL']*100:.1f}% of neat")
            r("Bilayer thickness", f"{f['thick']*100:.1f}% of neat")
            r("Gel-phase fraction", f"{f['gel']*100:.0f}%")
            rig = ("rigid — gel solid (~10× stiffer)" if f['gel'] > 0.6
                   else "stiffening" if f['gel'] > 0.25 else "fluid")
            r("Bending rigidity", rig)
            r("Fluidity", f"{f.get('fluid', 0)*100:.0f}% of neat")
            r("Defect / pore index", f"{f['pore']:.3f}")
            r("Accumulated damage", f"{f.get('Dmem', 0)*100:.0f}%")
            r("Membrane tension", f"{f['tension']:.2f} (1 = lytic)")
            r("Actomyosin pMLC", f"{f['pMLC']:.2f}")
            r("Focal adhesions", f"{f['FA']*100:.0f}%")
        elif key == "microtubules":
            r("Microtubules (compression struts)", f"{f['mt']*100:.0f}% · load σ {f.get('sigMT',0):.2f}")
            r("F-actin (tension cortex)", f"{f['actin']*100:.0f}% (cold + shear labile)")
            r("Intermediate filaments (cables)", f"{f.get('intf',1)*100:.0f}% · load σ {f.get('sigIF',0):.2f}")
            r("Tensegrity", "MT buckling under compression → GEF-H1 → RhoA → ROCK")
            r("Consequence", "mechanical stress also depresses ΔΨm (see mitochondria)")
        elif key == "cytosol":
            r("Cell volume", f"{f['Vn']*100:.0f}% of isotonic")
            r("CPA inside", f"{f['Cin']:.2f} M")
            gl = f.get('glass', 0)
            state = ("vitrified glass (arrested)" if gl > 0.85 else
                     "vitrifying" if gl > 0.4 else "viscous" if gl > 0.1 else "fluid")
            r("Cytoplasm state", f"{state}")
            r("Molecular mobility", f"{(1-gl)*100:.0f}% (WLF; 0 below Tg′)")
            pr = f.get('prot', 1.0)
            pstate = ("denatured" if pr < 0.6 else "stressed" if pr < 0.85 else "native")
            r("Protein integrity", f"{pr*100:.0f}% ({pstate})")
            r("CPA water-replacement", f"σ-excl / H-bond: {CPAS[self.P.cpa_key]['wRepl']:.2f}")
            r("ATP pool", f"{f['atp']*100:.0f}% of resting")
            r("Caspase-3", f"{f['casp3']:.2f}")
        elif key == "ice" or not key:
            r("Ice fraction", f"{f['fIce']*100:.1f}% of volume")
            r("Mean grain size", f"{f['grain']:.1f} µm")
            r("Unfrozen pocket", "—" if f['chanW'] > 500 else f"{f['chanW']:.2f} µm")
            r("Pocket / cell diameter", f"{f['squeeze']:.2f}")
            r("Extracellular osmolality", f"{f['Osme']:.2f} osmol/L")
            if self.P.additive != "none" and self.P.add_conc > 0:
                r("Extracellular additive", f"{ADDITIVES[self.P.additive]['name']}, {self.P.add_conc:.1f}% w/v")
            if R.get('junction', 0) > 0:
                r("Junction coupling", f"{R['junction']:.2f} (intercellular propagation)")
            r("P(intracellular ice)", f"{R['P_iif']*100:.0f}%")
        else:
            gk = GEOMETRY.get(key)
            if gk: r("Volume fraction (conventional)", f"{gk[0]*100:.1f}%  ({gk[1]} objects)")
            r("Cell volume", f"{f['Vn']*100:.0f}% of isotonic")
        h.append("</table>")
        self.insp.setHtml("\n".join(h))

    # ------------------------------------------------------------- analysis
    def _busy(self, on, msg=""):
        self.prog.setVisible(on)
        self.statusBar().showMessage(msg if on else "Ready")

    def _spawn(self, fn, cb, msg):
        self._busy(True, msg); self.anal.setHtml(f"<i>{msg}</i>")
        w = Worker(fn); self._workers.append(w)
        w.done.connect(lambda res: (self._busy(False), cb(res)))
        w.fail.connect(lambda e: (self._busy(False), self.anal.setPlainText(e)))
        w.start()

    def _pop(self):
        P = self.P
        self._spawn(lambda: simulate_population(P, n=120), self._pop_done,
                    "Simulating 120 individual cells across cores…")

    def _pop_done(self, pop):
        self.pop_plot.clear()
        b = np.array(pop["bins"], dtype=float)
        xs = (np.arange(20) + 0.5) * 5
        bg = pg.BarGraphItem(x=xs, height=b, width=4.2, brush=C["nucleus"])
        self.pop_plot.addItem(bg)
        self.pop_plot.setLabel("bottom", "24 h survival (%)")
        self.pop_plot.setLabel("left", "cells")
        rows = "".join(f"<tr><td>{lab}</td><td align=right style='font-family:monospace'>"
                       f"{rho:+.2f}</td></tr>" for _k, lab, rho in pop["attribution"])
        tail = sum(1 for c in pop["cells"] if c["s24"] < pop["median"] * 0.6) / pop["n"]
        note = ("The distribution is bimodal because nucleation is binary: a cell either forms "
                "intracellular ice or it does not, and the graded parameters barely move the "
                "survivors. That is what a real post-thaw viability histogram looks like."
                if tail > 0.02 else
                "The distribution is tight — under this protocol, cell-to-cell parameter variation "
                "does not change the outcome much. Push the protocol towards its failure boundary "
                "to see the population separate.")
        self.anal.setHtml(
            f"<b>Population of {pop['n']} cells</b><br>"
            f"mean <b>{pop['mean']*100:.0f}%</b> · median {pop['median']*100:.0f}% · "
            f"IQR {pop['p25']*100:.0f}–{pop['p75']*100:.0f}% · "
            f"10–90th {pop['p10']*100:.0f}–{pop['p90']*100:.0f}%<br>"
            f"{pop['nucleated']*100:.0f}% of cells nucleated intracellular ice · "
            f"mean functional recovery {pop['frec_mean']*100:.0f}%<br><br>"
            f"<b>What drives the spread</b> (Spearman ρ with survival)<table cellpadding=3>{rows}</table>"
            f"<br><i>{note}</i>")

    def _screen(self):
        P = self.P
        self._spawn(lambda: knockout_screen(P), self._screen_done,
                    "Disabling each model component in turn…")

    def _screen_done(self, S):
        rows = "".join(
            f"<tr><td>{r['name']}</td>"
            f"<td align=right style='font-family:monospace;color:"
            f"{C['crit'] if r['delta']<-0.01 else C['good'] if r['delta']>0.01 else C['muted']}'>"
            f"{r['delta']*100:+.0f} pts</td>"
            f"<td align=right style='font-family:monospace'>{r['ko']*100:.0f}%</td>"
            f"<td>{r['cls']}</td><td style='color:#7a7873'>{r['bench']}</td></tr>"
            for r in S["rows"])
        warn = ("" if S["headroom"] else
                f"<br><b style='color:{C['serious']}'>No headroom.</b> Baseline is "
                f"{S['wt']*100:.0f}%, so most perturbations cannot move it.")
        self.anal.setHtml(
            f"<b>In silico perturbation screen</b> — baseline {S['wt']*100:.0f}% at 24 h.{warn}"
            "<table cellpadding=3><tr><th align=left>Component disabled</th><th>Δ 24 h</th>"
            "<th>Survival</th><th align=left>Class</th><th align=left>Bench analogue</th></tr>"
            f"{rows}</table>")

    def _next(self):
        P = self.P
        self._spawn(lambda: next_experiment(P), self._next_done,
                    "Probing the model's own sensitivity…")

    def _next_done(self, N):
        rows = "".join(
            f"<tr><td align=right style='color:#7a7873'>{i+1}</td><td>{r['label']}</td>"
            f"<td>{r['tier']}</td>"
            f"<td align=right style='font-family:monospace'>{r['influence']*100:.0f} pts</td>"
            f"<td style='color:#7a7873'>{r['how']}</td></tr>"
            for i, r in enumerate(N["rows"]))
        colr = C["good"] if N["confident"] else C["warn"] if N["band"] < .2 else C["crit"]
        self.anal.setHtml(
            f"<b>Lab in the loop</b><br><span style='font-size:20px'>{N['base']*100:.0f}% "
            f"± {N['band']*100:.0f} points</span><br>"
            f"<b style='color:{colr}'>{N['verdict']}</b><br><br>"
            "<table cellpadding=3><tr><th></th><th align=left>Measure this</th><th>Currently</th>"
            f"<th>Swing</th><th align=left>How</th></tr>{rows}</table>")

    def _stab(self):
        P = self.P
        self._spawn(lambda: stability_check(P), self._stab_done,
                    "Re-integrating at half and double step size…")

    def _stab_done(self, S):
        colr = C["good"] if S["verdict"] == "stable" else C["warn"] if S["verdict"] == "marginal" else C["crit"]
        self.anal.setHtml(
            "<b>Numerical stability</b><br>Cui et al. name mathematical instability as the core "
            "weakness of whole-cell models built this way. That is testable, so it gets tested."
            "<table cellpadding=4>"
            f"<tr><td>Half step</td><td align=right style='font-family:monospace'>{S['half']*100:.2f}%</td></tr>"
            f"<tr><td>Standard step</td><td align=right style='font-family:monospace'><b>{S['base']*100:.2f}%</b></td></tr>"
            f"<tr><td>Double step</td><td align=right style='font-family:monospace'>{S['dbl']*100:.2f}%</td></tr>"
            f"<tr><td>Maximum drift</td><td align=right style='font-family:monospace;color:{colr}'>"
            f"<b>{S['drift']*100:.2f} points</b></td></tr>"
            f"<tr><td>Verdict</td><td align=right style='color:{colr}'><b>{S['verdict'].upper()}</b></td></tr>"
            "</table>")


def main():
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    w = Main(); w.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
