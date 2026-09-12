"""CryoCell desktop application."""
from __future__ import annotations
import sys, math, traceback, os, shutil, subprocess
import numpy as np
from PyQt6.QtCore import Qt, QThread, pyqtSignal, QTimer, QPointF, QRectF
from PyQt6.QtGui import (QFont, QColor, QAction, QImage, QPainter, QPen, QPainterPath,
    QPolygonF, QRadialGradient, QLinearGradient, QStandardItemModel, QStandardItem)
from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout,
    QHBoxLayout, QGridLayout, QLabel, QSlider, QComboBox, QPushButton, QGroupBox,
    QScrollArea, QCheckBox, QTabWidget, QTextEdit, QSplitter, QSizePolicy,
    QTableWidget, QTableWidgetItem, QHeaderView, QProgressBar, QFileDialog,
    QStackedWidget)
import pyqtgraph as pg

from .model import Params, simulate, CPAS, ADHESION_STATES, ADDITIVES, osm_from_dT
from .analysis import (simulate_population, knockout_screen, stability_check,
                       next_experiment)
from .hpa import HPA, HPA_TOTALS, GEOMETRY, summary_line
from .cellview import (CellView, LAYERS, C, STRESS_PATHWAYS, stress_activity,
                       pathway_applicable, cell_has_compartment)
from .pathways import REACTOME, REACTOME_RETRIEVED
try:
    from .pathways import HPA_COMPARTMENTS, HPA_RETRIEVED
except Exception:
    HPA_COMPARTMENTS, HPA_RETRIEVED = {}, ""
_COMP_SHORT = {"cytosol":"cytosol","mitochondria":"mito","nucleoplasm":"nucleus","er":"ER",
               "golgi":"Golgi","plasma_mem":"membrane","vesicles":"vesicles","nucleoli":"nucleoli",
               "actin":"actin","focal_adh":"FA","interm_fil":"IF","microtubules":"MT",
               "peroxisomes":"perox","lysosomes":"lyso","endosomes":"endo","centrosome":"centro",
               "nuclear_mem":"nuc.env","lipid_drop":"lipid"}

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
    # Adipose-derived MSC (hADSC), the cell type characterised by Li et al. 2020
    # (J Therm Biol 93:102689, doi:10.1016/j.jtherbio.2020.102689) for water
    # transport + surface-catalysed nucleation. Encoded as a PRIMARY-MSC
    # LITERATURE PRIOR (size / permeability / primary-cell phenotype), NOT the
    # paper's exact fitted Lp, Ea, Omega, kappa — those numeric values were not
    # available to encode. Provenance flag: literature_prior, not cell-type-validated.
    "Cell line — hADSC (Li 2020, prior)": dict(Viso=1700, cyto=1.30, sterol=22, lp=0.25,
                                         adhesion="suspension",
                                         apop_resist=0.0, anoikis_resist=0.0,
                                         glycolytic=0.0, antioxidant=0.0),
    # Red blood cell — the extreme case: tiny, ANUCLEATE, no organelles, a very
    # high water permeability (AQP1) and a cholesterol-rich membrane. Its huge
    # surface-to-volume and high Lp let it dehydrate almost instantly, so it
    # rarely forms intracellular ice even at fast cooling — but it is osmotically
    # fragile. No caspase apoptosis (no nucleus/mito), so apop_resist = 1.
    # Values are LITERATURE PRIORS (not calibrated). Lp ~1.6-3 um/min/atm (AQP1).
    # RBC uses GLYCEROL, not DMSO — the standard clinical cryoprotectant for red
    # cells (low-glycerol ~20% here; the high-glycerol method uses ~40%). Step
    # dilution mimics deglycerolisation. The RBC solute permeability is HIGH
    # (ps = 0.25), because human red cells carry the AQP3 aquaglyceroporin and are
    # unusually glycerol-permeable — the reason glycerol works for them. With a
    # low ps the cell cannot shed glycerol fast enough on dilution and lyses
    # (osmotic swelling past the lytic limit); ps = 0.25 gives a realistic ~64%
    # recovery for the 20% low-glycerol method. Literature-prior parameters.
    "Cell — Red blood cell (prior)":     dict(cell_type="rbc", Viso=90, lp=1.8, ps=0.25,
                                         sterol=40, cyto=0.55, nuc_scale=0.4, adhesion="suspension",
                                         cpa_key="glycerol", conc_pct=20, T_add=22, dilution="step",
                                         apop_resist=1.0, anoikis_resist=1.0, glycolytic=1.0, antioxidant=0.4),
    # Avian red blood cell — UNLIKE the mammalian RBC it is NUCLEATED, retains
    # mitochondria and does limited protein synthesis, and is larger. Bissoyi,
    # Oral, Choudhary & Gibson 2025 (ACS Polym Au 6:366, doi:10.1021/acspolymersau.5c00151):
    # 10% DMSO + 60 mg/mL polyampholyte (~6% w/v, COOH-PLL here) gives ~90%
    # recovery, outperforming glycerol, with reduced apoptosis and recovered actin
    # polymerisation post-thaw. Nucleated -> apoptosis IS possible (apop_resist 0).
    # Literature-prior parameters, not calibrated.
    "Cell — Avian red blood cell (Bissoyi 2025)": dict(cell_type="avian_rbc", Viso=140, lp=1.2,
                                         ps=0.05, sterol=35, cyto=0.60, nuc_scale=0.55,
                                         adhesion="suspension", cpa_key="dmso", conc_pct=10,
                                         additive="pll", add_conc=6, T_add=22, dilution="step",
                                         apop_resist=0.0, anoikis_resist=1.0, glycolytic=0.5, antioxidant=0.4),
    # T-lymphocyte — small round cell with a high nucleus-to-cytoplasm ratio and
    # only a thin cytoplasmic rim. Moderate Lp; a common DMSO cryopreservation
    # target. Literature-prior parameters, not calibrated.
    "Cell — T lymphocyte (prior)":       dict(cell_type="tcell", Viso=180, lp=0.30, sterol=25,
                                         cyto=0.90, adhesion="suspension",
                                         apop_resist=0.0, anoikis_resist=0.0, glycolytic=0.0, antioxidant=0.0),
    # Platelet — a tiny ANUCLEATE cell fragment (~8 fL) full of granules. Poorly
    # served by conventional freezing; trehalose-loaded freeze-drying is the
    # active research route (see the freeze-dry option). Literature-prior params.
    "Cell — Platelet (prior)":           dict(cell_type="platelet", Viso=8, lp=0.40, sterol=25,
                                         cyto=0.70, adhesion="suspension",
                                         apop_resist=1.0, anoikis_resist=1.0, glycolytic=0.5, antioxidant=0.3),
    # Freeze-drying (lyophilisation) — freeze, sublimate the water away, store dry
    # at room temperature, rehydrate before use. Trehalose + a permeant water-
    # replacement CPA (glycerol) is the best the model can do; full cytoplasmic
    # protection needs INTRACELLULAR trehalose (an unsolved delivery problem), so
    # recovery is honestly low. Active research for platelets and RBC.
    "Freeze-dry — Platelet (trehalose)": dict(cell_type="platelet", Viso=8, lp=0.40, sterol=25, cyto=0.70,
                                         adhesion="suspension", apop_resist=1.0, anoikis_resist=1.0,
                                         glycolytic=0.5, antioxidant=0.3,
                                         freeze_dry=True, cpa_key="glycerol", conc_pct=10,
                                         additive="tre", add_conc=6, dry_residual=0.05),
    "Freeze-dry — RBC (trehalose)":      dict(cell_type="rbc", Viso=90, lp=1.8, ps=0.05, sterol=40,
                                         cyto=0.55, nuc_scale=0.4, adhesion="suspension",
                                         apop_resist=1.0, anoikis_resist=1.0, glycolytic=1.0, antioxidant=0.4,
                                         freeze_dry=True, cpa_key="glycerol", conc_pct=15,
                                         additive="tre", add_conc=6, dry_residual=0.05),
    # 3D spheroid / organoid construct (Gao, Bissoyi, Guo & Gibson 2024, ACS
    # Biomater Sci Eng 11:208, doi:10.1021/acsbiomaterials.4c00958). 10% DMSO
    # alone supercools (~-16 C), which sheds surface cells and perforates the
    # interior; open the "3D spheroid" tab, tick the extracellular ice nucleator
    # (IN+) and change the diameter to see the protection and the size effect.
    # Cells are junction-coupled (adhesion = spheroid): Irimia & Karlsson 2002.
    "Construct — 3D spheroid (Bissoyi 2024)": dict(Viso=1800, cyto=1.35, sterol=25,
                                         adhesion="spheroid", cpa_key="dmso", conc_pct=10,
                                         apop_resist=0.0, anoikis_resist=0.0,
                                         glycolytic=0.0, antioxidant=0.0),
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

# One-line explanation for each control (what it does + an honesty caveat where
# the value is a teaching/literature prior rather than a measured constant).
# Shown as muted microcopy under the control, in the spirit of an honest,
# teaching-first scientific simulator.
HELP = {
    # combos
    "cpa_key":  "The permeating cryoprotectant: sets toxicity, glass-forming ability and permeation speed.",
    "additive": "A non-permeating co-solute that stays outside the cell (ice-recrystallisation inhibition, membrane stabilisation).",
    "adhesion": "Cell–cell / matrix state. Junction-coupled states let intracellular ice propagate between cells.",
    "add_steps":"How the CPA is introduced; more steps reduce osmotic shock during loading.",
    "dilution": "How the CPA is washed out after thaw. Stepwise sucrose limits osmotic swelling; direct risks lysis.",
    "rock_drug":"ROCK-pathway inhibitor that blunts cold-induced blebbing and anoikis.",
    "rock_when":"When the ROCK inhibitor is present in the medium.",
    # sliders
    "conc_pct": "Concentration of the penetrating CPA. More protects against ice but raises toxicity and osmotic stress.",
    "add_conc": "Concentration of the extracellular additive (full membrane coverage near ~6% w/v).",
    "sterol":   "Membrane cholesterol. Stiffens the bilayer and shifts the gel-phase transition.",
    "iri":      "Ice-recrystallisation inhibition added on top of the CPA (slows grain growth during warming).",
    "T_add":    "Temperature at which the CPA is added; warmer speeds permeation but raises toxicity exposure.",
    "hold_min": "Equilibration time before cooling; longer lets the CPA reach the cell interior.",
    "T_seed":   "Temperature at which extracellular ice is induced. Warmer seeding = less supercooling = gentler freezing.",
    "CR":       "Cooling rate. The Mazur two-factor optimum: too fast → intracellular ice; too slow → solution injury.",
    "T_store":  "Long-term storage temperature (−196 °C = liquid nitrogen).",
    "days":     "Time held in storage before warming.",
    "WR":       "Warming rate. Fast warming limits ice recrystallisation on the way out.",
    "rock_conc":"Dose of the selected ROCK inhibitor.",
    "zvad":     "Pan-caspase inhibitor (z-VAD-fmk): blocks apoptotic execution. % is an efficacy prior.",
    "calpain_i":"Calpain protease inhibitor: blunts calcium-driven necrosis. % is an efficacy prior.",
    "gsmtx":    "GsMTx4 blocks the Piezo1 mechanosensitive channel. % is an efficacy prior.",
    "T_recover":"Post-thaw culture temperature during the recovery window.",
    "recover_h":"How long recovery is scored after thaw (delayed-onset death still accrues).",
    "lp":       "Water permeability of the membrane. Literature-typical prior; higher = faster dehydration.",
    "ps":       "Solute (CPA) permeability. Literature-typical prior; higher = faster CPA in/out.",
    "vb":       "Osmotically inactive volume fraction (Boyle–van 't Hoff intercept).",
    "Viso":     "Isotonic cell volume; sets the surface-to-volume ratio that governs dehydration.",
    "thalf":    "Half-life of the self-deactivating CPA linker at 37 °C (only for the SD-CPA).",
    "frag_tox": "Toxicity of SD-CPA breakdown fragments relative to the parent molecule.",
    "frag_n":   "Number of fragments each SD-CPA molecule breaks into.",
    "apop_resist":   "Cell-line phenotype: resistance to apoptosis (0 = primary/normal, 1 = transformed).",
    "anoikis_resist":"Resistance to detachment-induced death (0 = anchorage-dependent).",
    "antioxidant":   "Capacity to buffer reactive oxygen species during recovery.",
    "glycolytic":    "Warburg metabolism: how much ATP is made without mitochondria.",
    "cyto":          "Cytoskeletal density; scales mechanical stiffness and cold-labile depolymerisation.",
}


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
          "warm":"#eb6834","melt":"#eb6834","dilute":"#22d3ee","recover":"#1baf7a","end":"#1baf7a",
          "dry1":"#c99a4a","dry2":"#d9a441","drystore":"#a0895c","rehydrate":"#22d3ee"}

class TimelineBar(QWidget):
    """Full-width freeze-thaw timeline: phase bands + temperature curve + scrubber."""
    seek = pyqtSignal(int)
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setFixedHeight(74); self.setAutoFillBackground(True); self.setMouseTracking(True)
        self.T = []; self.phase = []; self.Vn = []; self.pos = 0; self._drag = False
    def set_data(self, T, phase, Vn=None):
        self.T = list(T); self.phase = list(phase)
        self.Vn = list(Vn) if Vn is not None else []; self.update()
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
        # cell-volume trace (shrink-then-swell on loading; deep dehydration on freezing)
        if self.Vn:
            yV = lambda v: 6 + (1 - min(max(v, 0), 1.5) / 1.5) * (h - 24)
            q.setPen(QPen(QColor(64, 150, 226), 1.6)); vpath = QPainterPath()
            for i, v in enumerate(self.Vn):
                x, y = self._xf(i, w, pad), yV(v)
                vpath.moveTo(x, y) if i == 0 else vpath.lineTo(x, y)
            q.drawPath(vpath)
            q.setFont(QFont("", 7, QFont.Weight.Bold))
            q.setPen(QColor("#e0a63a")); q.drawText(QPointF(w - 74, 14), "T")
            q.setPen(QColor(64, 150, 226)); q.drawText(QPointF(w - 60, 14), "volume")
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

        # The FA-LINC axis is focal-adhesion-to-nucleus mechanotransduction; it does
        # not exist in an anucleate cell (RBC, platelet), nor in a non-adherent
        # nucleated cell that has no focal adhesions (avian RBC).
        _ct = f.get("cell_type")
        if _ct in ("rbc", "platelet", "avian_rbc"):
            reason = {
                "rbc": "A red blood cell is anucleate — no nucleus, nuclear lamina or "
                       "LINC complex, so there is no focal-adhesion-to-nucleus "
                       "mechanotransduction axis.",
                "platelet": "A platelet is anucleate — no nucleus, nuclear lamina or "
                            "LINC complex, so there is no focal-adhesion-to-nucleus "
                            "mechanotransduction axis.",
                "avian_rbc": "An avian red blood cell is nucleated but non-adherent — a "
                             "spectrin membrane skeleton, no focal adhesions — so there "
                             "is no focal-adhesion-to-nucleus mechanotransduction axis.",
            }[_ct]
            q.setPen(QColor(C["text2"])); q.setFont(QFont("", 12, QFont.Weight.Bold))
            q.drawText(QRectF(20, h/2 - 30, w - 40, 24), Qt.AlignmentFlag.AlignHCenter,
                       "Not applicable")
            q.setPen(QColor(C["muted"])); q.setFont(QFont("", 9))
            q.drawText(QRectF(20, h/2 - 2, w - 40, 60), Qt.AlignmentFlag.AlignHCenter | Qt.TextFlag.TextWordWrap,
                       reason)
            return

        # ---- model state → mechanical drivers along the FA–LINC axis ----
        FA = f.get("FA", 0); ten = f.get("tension", 0); act = f.get("actin", 0)
        mt = f.get("mt", 0); intf = f.get("intf", 0); rock = f.get("rock", 0)
        yap = f.get("yapN", 0); pMLC = f.get("pMLC", 0)
        force = clamp(pMLC * (0.3 + 0.7 * act), 0, 1)
        engaged = clamp(FA / 0.5, 0, 1)                 # adhesion engagement
        Vn = f.get("Vn", 1.0)
        shrink = clamp((1.0 - Vn) / 0.55, 0, 1)         # osmotic dehydration
        iif = clamp(f.get("Piif", 0), 0, 1)             # intracellular-ice probability
        frozen = f.get("frozen", False)
        # Mazur (1963, 1972) two-factor hypothesis, expressed mechanically on this axis:
        #   slow cooling  → extracellular ice → dehydration/"solution effects" → shrink
        #   fast cooling  → intracellular ice formation (IIF) → filament/envelope rupture
        sol_arm = shrink if iif < 0.6 else shrink * 0.4
        iif_arm = iif if frozen else iif * 0.6
        shatter = clamp(iif_arm, 0, 1)                  # network fragmentation
        buckle = clamp(sol_arm, 0, 1)                   # slack / buckled filaments

        good = QColor(47, 192, 136); bad = QColor(224, 87, 79); amb = QColor(224, 166, 58)
        gold = QColor(198, 156, 74); violet = QColor(139, 127, 240)
        ifcol = QColor(58, 58, 66); mtcol = QColor(86, 138, 190)
        nuc_fill = QColor(247, 232, 150)               # euchromatin ground (as in the reference)
        def mix(a, b, t):
            t = clamp(t, 0, 1)
            return QColor(int(a.red()+(b.red()-a.red())*t), int(a.green()+(b.green()-a.green())*t),
                          int(a.blue()+(b.blue()-a.blue())*t))
        def label(s, cx, cy, wd=90, col=None, sz=7, bold=True):
            q.setPen(QColor(col or C["text"])); q.setFont(QFont("", sz, QFont.Weight.Bold if bold else QFont.Weight.Normal))
            q.drawText(QRectF(cx-wd/2, cy-7, wd, 14), Qt.AlignmentFlag.AlignHCenter, s)
        def blob(cx, cy, r, col, lab=None, pP=False):
            q.setBrush(col); q.setPen(QPen(col.darker(150), 1.4))
            q.drawEllipse(QPointF(cx, cy), r, r)
            if lab: label(lab, cx, cy, 70, C["text"], 7)
            if pP:
                q.setBrush(amb); q.setPen(Qt.PenStyle.NoPen)
                q.drawEllipse(QPointF(cx+r*0.7, cy-r*0.7), 4.5, 4.5)
        def rng(x0, y0, x1, y1):  # length of a segment
            return math.hypot(x1-x0, y1-y0)
        def perp(x0, y0, x1, y1):
            dx, dy = -(y1-y0), (x1-x0); L = math.hypot(dx, dy) or 1; return dx/L, dy/L

        # ---- actin filament: two beaded strands. `poly` is the polymerised
        # fraction (S.actin) — as it falls (cold depolymerisation) the cable
        # dissolves into scattered free G-actin monomers; buckles when the cell
        # dehydrates (slack); gaps open where intracellular ice shatters it.
        def actin(x0, y0, x1, y1, base, poly, weight=1.0):
            poly = clamp(poly, 0, 1)
            n = max(6, int(rng(x0, y0, x1, y1) / 6)); dx, dy = perp(x0, y0, x1, y1)
            amp = 1.6 + 4.5 * buckle * weight
            col = QColor(base); col.setAlpha(int(90 + 150 * poly))
            mono = QColor(base); mono.setAlpha(75)
            rs = np.random.RandomState(int(abs(x0 + y1 * 3)) % 9999)
            q.setPen(Qt.PenStyle.NoPen)
            for k in range(n):
                t = k/(n-1)
                if shatter*weight > 0.25 and (k % max(2, int(2+(1-shatter)*7)) == 0): continue
                bx = x0+(x1-x0)*t + dx*math.sin(t*math.pi*3)*amp
                by = y0+(y1-y0)*t + dy*math.sin(t*math.pi*3)*amp
                if ((k*0.61803) % 1.0) <= poly:              # segment still polymerised
                    q.setBrush(col)
                    for s in (-1.6, 1.6): q.drawEllipse(QPointF(bx+dx*s, by+dy*s), 2.2, 2.2)
                else:                                        # depolymerised → free monomer
                    q.setBrush(mono)
                    q.drawEllipse(QPointF(bx+rs.uniform(-10, 10), by+rs.uniform(-6, 6)), 1.8, 1.8)

        # ---- intermediate filament: thin dark wavy line
        def if_fil(x0, y0, x1, y1, col):
            dx, dy = perp(x0, y0, x1, y1); p = QPainterPath(QPointF(x0, y0)); n = 14
            for k in range(1, n+1):
                t = k/n; a = math.sin(t*math.pi*4)*(3.0+2*buckle)
                p.lineTo(QPointF(x0+(x1-x0)*t+dx*a, y0+(y1-y0)*t+dy*a))
            q.setPen(QPen(col, 1.6)); q.setBrush(Qt.BrushStyle.NoBrush); q.drawPath(p)

        # ---- microtubule: hollow dashed tube. `poly` (S.mt) is the polymerised
        # fraction — the tube catastrophes from the plus-end (membrane side) inward
        # as it falls, leaving scattered free tubulin dimers (cold-labile, like actin).
        def mtub(x0, y0, x1, y1, base, poly):
            poly = clamp(poly, 0, 1); dx, dy = perp(x0, y0, x1, y1)
            L = rng(x0, y0, x1, y1); n = max(4, int(L/10))
            col = QColor(base); col.setAlpha(int(70 + 160 * poly))
            if poly > 0.12:                                  # intact segment from nucleus/MTOC end
                xe, ye = x0+(x1-x0)*poly, y0+(y1-y0)*poly
                for s in (-2.2, 2.2):
                    q.setPen(QPen(col, 1.3)); q.drawLine(QPointF(x0+dx*s, y0+dy*s), QPointF(xe+dx*s, ye+dy*s))
                q.setPen(QPen(col, 0.8))
                for k in range(int(n*poly)):
                    t = k/max(1, n-1); xx = x0+(x1-x0)*t; yy = y0+(y1-y0)*t
                    q.drawLine(QPointF(xx+dx*2.2, yy+dy*2.2), QPointF(xx-dx*2.2, yy-dy*2.2))
            if poly < 0.85:                                  # free tubulin where depolymerised
                mc = QColor(base); mc.setAlpha(70); q.setBrush(mc); q.setPen(Qt.PenStyle.NoPen)
                rs = np.random.RandomState(int(abs(x1 + y0 * 3)) % 9999)
                for _ in range(int(9*(1-poly))):
                    t = rs.uniform(min(poly, 0.95), 1.0); xx = x0+(x1-x0)*t; yy = y0+(y1-y0)*t
                    q.drawEllipse(QPointF(xx+rs.uniform(-8, 8), yy+rs.uniform(-6, 6)), 1.8, 1.8)

        def arrow(x0, y0, x1, y1, wgt, col):  # thick force arrow
            q.setPen(QPen(col, wgt)); q.drawLine(QPointF(x0, y0), QPointF(x1, y1))
            ang = math.atan2(y1-y0, x1-x0); s = 7+wgt
            poly = QPolygonF([QPointF(x1, y1),
                QPointF(x1-s*math.cos(ang-0.5), y1-s*math.sin(ang-0.5)),
                QPointF(x1-s*math.cos(ang+0.5), y1-s*math.sin(ang+0.5))])
            q.setBrush(col); q.setPen(Qt.PenStyle.NoPen); q.drawPolygon(poly)

        # ================= NUCLEUS INTERIOR (top) =================
        q.setPen(QColor(C["muted"])); q.setFont(QFont("", 8, QFont.Weight.Bold))
        q.drawText(QRectF(0, 2, w, 14), Qt.AlignmentFlag.AlignHCenter, "Nucleus interior")
        nx, nyc, nrx, nry = X(0.50), Y(0.155), X(0.44)*(0.85+0.15*(1-shrink)), Y(0.125)*(0.8+0.2*(1-shrink))
        q.setBrush(nuc_fill); q.setPen(QPen(violet.darker(120), 1.4))
        q.drawEllipse(QPointF(nx, nyc), nrx, nry)
        # euchromatin (open pink loops) + G-actin/F-actin inside
        rs = np.random.RandomState(7)
        q.setPen(QPen(QColor(214, 120, 170, 150), 1.3)); q.setBrush(Qt.BrushStyle.NoBrush)
        for _ in range(9):
            a = rs.uniform(0, 2*math.pi); rr = rs.uniform(0.15, 0.6)
            cxp, cyp = nx+math.cos(a)*nrx*rr, nyc+math.sin(a)*nry*rr
            q.drawArc(QRectF(cxp-7, cyp-4, 14, 8), 0, 300*16)
        q.setBrush(QColor(150, 90, 200)); q.setPen(Qt.PenStyle.NoPen)   # G-actin monomers
        for _ in range(7):
            a = rs.uniform(0, 2*math.pi); rr = rs.uniform(0.1, 0.5)
            q.drawEllipse(QPointF(nx+math.cos(a)*nrx*rr, nyc+math.sin(a)*nry*rr), 2.3, 2.3)
        label("euchromatin", nx, nyc-nry*0.5, 90, QColor(150,70,120), 7)
        label("G-actin", X(0.66), Y(0.075), 60, QColor(120,70,170), 7)
        # heterochromatin clumps tethered to the lamina (LAD) around the rim
        q.setBrush(QColor(120, 96, 40)); q.setPen(Qt.PenStyle.NoPen)
        for k in range(12):
            a = math.pi + k/11*math.pi   # lower rim, facing envelope
            cxp, cyp = nx+math.cos(a)*nrx*0.86, nyc+math.sin(a)*nry*0.86
            q.drawEllipse(QPointF(cxp, cyp), 4.5, 3.2)
        label("heterochromatin", nx, nyc+nry*0.55, 110, QColor(90,70,30), 7)

        # ================= NUCLEAR ENVELOPE + LINC =================
        onm = Y(0.315); inm = Y(0.345)                 # outer / inner nuclear membrane
        for my, nm in ((onm, "ONM"), (inm, "INM")):
            q.setBrush(QColor(236, 206, 120)); q.setPen(QPen(QColor(180, 150, 70), 1))
            q.drawRoundedRect(QRectF(X(0.10), my-4, X(0.80), 8), 4, 4)
            label(nm, X(0.955)*w/w*0.0+X(0.955), my, 40, QColor(150,120,60), 7)
        label("Perinuclear space", X(0.86), (onm+inm)/2+16, 120, C["muted"], 7)
        # nuclear lamina (mesh lining the INM), loads amber → red under mechanical force
        lam = mix(violet, bad, force)
        q.setPen(QPen(lam, 2)); q.drawLine(QPointF(X(0.12), inm+9), QPointF(X(0.88), inm+9))
        for k in range(9): q.drawLine(QPointF(X(0.13+k*0.083), inm+6), QPointF(X(0.17+k*0.083), inm+12))
        label("nuclear lamina (lamin A/C)", nx, inm+20, 200, lam, 7)
        # LAD tethers (heterochromatin → lamina)
        q.setPen(QPen(QColor(120, 96, 40, 170), 1));
        for k in range(5):
            xx = X(0.30+k*0.10); q.drawLine(QPointF(xx, nyc+nry*0.8), QPointF(xx, inm+6))
        label("LAD", X(0.20), inm+2, 40, QColor(90,70,30), 7)
        # NPC — nuclear pore spanning both membranes
        npx = X(0.20)
        q.setBrush(QColor(150, 110, 180)); q.setPen(QPen(QColor(110,80,140),1))
        q.drawRoundedRect(QRectF(npx-9, onm-6, 18, (inm-onm)+12), 4, 4)
        q.setBrush(QColor(C["surface0"])); q.setPen(Qt.PenStyle.NoPen)
        q.drawRect(QRectF(npx-3, onm-4, 6, (inm-onm)+8))
        label("NPC", npx, onm-12, 50, QColor(110,80,140), 7)
        # LINC: nesprin (ONM, binds cytoskeleton) + SUN (INM, binds lamina)
        for lx in (0.42, 0.62):
            cx = X(lx); lc = mix(QColor(C["muted"]), good, act)
            blob(cx, onm, 7, lc); blob(cx, inm, 7, mix(QColor(C["muted"]), violet, 0.6))
            q.setPen(QPen(lc, 2)); q.drawLine(QPointF(cx, onm+5), QPointF(cx, inm-5))  # SUN-nesprin bridge
        label("LINC (nesprin–SUN)", X(0.52), onm-12, 150, good.darker(110), 7)
        # formin at the perinuclear area, nucleating an actin cable
        fmx = X(0.76)
        q.setBrush(QColor(120, 200, 150)); q.setPen(QPen(QColor(70,150,100),1.2))
        q.drawEllipse(QPointF(fmx, onm+2), 9, 6)
        label("formin", fmx, onm-11, 55, QColor(60,140,95), 7)

        # ================= CYTOSKELETON (cytoplasm) =================
        q.setPen(QColor(C["muted"])); q.setFont(QFont("", 8, QFont.Weight.Bold))
        q.drawText(QRectF(0, Y(0.50)-7, X(0.30), 14), Qt.AlignmentFlag.AlignLeft, "  Cytoplasm")
        memY = Y(0.86)
        # actin cables: LINC/formin → focal adhesions (the force-bearing lines).
        # Driven by S.actin, so they dissolve into free monomers in the cold.
        abase = mix(QColor(150,120,60), gold, 0.7)
        actin(X(0.42), onm+6, X(0.28), memY-6, abase, act, 1.0)   # to FA-left
        actin(X(0.62), onm+6, X(0.68), memY-6, abase, act, 1.0)   # to FA-right
        actin(fmx, onm+6, X(0.50), memY-6, abase, act, 0.8)       # formin cable → membrane
        # myosin II mini-filaments — only where an actin cable exists
        for m in range(4):
            t = (m+0.5)/4; xx = X(0.42)+(X(0.50)-X(0.42))*t; yy = (onm+6)+(memY-6-(onm+6))*t
            mc = mix(good, bad, force); mc.setAlpha(int(255*clamp(act, 0, 1)))
            q.setBrush(mc); q.setPen(Qt.PenStyle.NoPen); q.drawEllipse(QPointF(xx, yy), 3.2, 2)
        acol = mix(QColor(150,120,60), gold, act)
        label(f"Actin filament ({act*100:.0f}%)", X(0.30), Y(0.58), 130, acol.darker(130), 7)
        # intermediate filaments (vimentin) — cold-stable, persist when actin/MT are gone
        if_fil(X(0.30), inm+14, X(0.20), memY-6, mix(QColor(90,90,98), ifcol, intf))
        if_fil(X(0.70), inm+14, X(0.80), memY-6, mix(QColor(90,90,98), ifcol, intf))
        label("Intermediate filament (cold-stable)", X(0.18), Y(0.66), 190, ifcol, 7)
        # microtubules radiating from the perinuclear MTOC region (also cold-labile)
        mtub(X(0.52), inm+16, X(0.44), memY-6, mtcol, mt)
        mtub(X(0.56), inm+16, X(0.60), memY-6, mtcol, mt)
        label(f"Microtubule ({mt*100:.0f}%)", X(0.74), Y(0.62), 120, mtcol, 7)

        # ================= PLASMA MEMBRANE + FOCAL ADHESIONS =================
        # lipid bilayer (two rows of heads with tails)
        q.setPen(Qt.PenStyle.NoPen)
        for row, yy in ((0, memY-4), (1, memY+4)):
            for k in range(int(X(0.05)), int(X(0.95)), 9):
                q.setBrush(QColor(210, 205, 150)); q.drawEllipse(QPointF(k, yy), 3, 3)
                q.setPen(QPen(QColor(190, 185, 130), 1))
                q.drawLine(QPointF(k, yy+(4 if row else -4)), QPointF(k, memY))
                q.setPen(Qt.PenStyle.NoPen)
        label("Plasma membrane", X(0.50), memY-14, 140, C["muted"], 7)
        # focal adhesions: integrin clusters + plaque, engaged→green / lost→red
        ecol = mix(bad, good, engaged)
        for fx in (0.28, 0.68):
            cx = X(fx)
            for dxi in (-8, 0, 8):                      # integrin cluster
                blob(cx+dxi, memY, 5, ecol)
            q.setBrush(mix(QColor(C["muted"]), good, engaged)); q.setPen(Qt.PenStyle.NoPen)
            q.drawRoundedRect(QRectF(cx-16, memY+7, 32, 7), 3, 3)   # talin/vinculin/FAK plaque
            label("FA", cx, memY+22, 40, ecol.darker(120), 8)
        # ECM hatching + substrate force arrows (the big arrows in the reference)
        q.setPen(QPen(QColor(C["muted"]), 1))
        for k in range(int(X(0.05)), int(X(0.95)), 10):
            q.drawLine(QPointF(k+6, memY+30), QPointF(k, memY+36))
        label("ECM · substrate", X(0.50), memY+40, 140, C["muted"], 7)
        traction = 0.3 + 0.7*engaged*force
        arrow(X(0.30), memY+30, X(0.10), memY+30, 3+5*traction, mix(QColor(C["muted"]), bad, traction))
        arrow(X(0.66), memY+30, X(0.90), memY+30, 3+5*traction, mix(QColor(C["muted"]), bad, traction))
        # force transmitted up the axis to the nucleus (curved implied by two arrows)
        arrow(X(0.90), memY-10, X(0.90), Y(0.36), 2+4*force, mix(QColor(C["muted"]), amb, force))
        label("force → nucleus", X(0.90), Y(0.50), 90, mix(QColor(C["muted"]), amb, force), 6)

        # ================= MAZUR TWO-FACTOR OVERLAY =================
        # Only name a Mazur injury regime in a freezing context; CPA-loading
        # shrinkage is intentional equilibration, and suspension FA is low by design.
        phase = f.get("phase", "")
        freezing_ctx = frozen or phase in ("cool", "seed", "store", "hold", "warm")
        banner = None
        if freezing_ctx and shatter >= buckle and shatter > 0.15:
            banner = (f"IIF regime — intracellular ice fragments the network (P_iif {iif*100:.0f}%)", bad)
        elif freezing_ctx and buckle > 0.15:
            banner = (f"Solution-effects regime — dehydration collapses the cytoskeleton (V {Vn*100:.0f}%)", amb)
        elif phase == "load" and shrink > 0.25:
            banner = (f"CPA loading — osmotic shrinkage (V {Vn*100:.0f}%)", amb)
        elif phase == "recover" and engaged < 0.20:
            banner = ("Adhesion not re-established → anoikis risk", bad)
        if banner:
            msg, mc = banner
            q.setBrush(QColor(mc.red(), mc.green(), mc.blue(), 30)); q.setPen(QPen(mc, 1))
            q.drawRoundedRect(QRectF(X(0.04), h-20, X(0.92), 16), 4, 4)
            q.setPen(mc); q.setFont(QFont("", 7, QFont.Weight.Bold))
            q.drawText(QRectF(X(0.04), h-20, X(0.92), 16), Qt.AlignmentFlag.AlignCenter, msg)
        # YAP/TAZ readout (nuclear mechanotransduction outcome) — top-right, clear of chromatin
        label(f"YAP/TAZ nuclear {yap*100:.0f}%", X(0.80), Y(0.035), 150,
              mix(bad, good, yap).darker(120), 7)


class StressView(QWidget):
    """BioGPU-style cryo stress-pathway readout. Each pathway lights up by an
    activity the engine actually computes (see cellview.stress_activity), anchored
    to its Human Protein Atlas compartment. Pathways the model does not simulate
    are shown greyed and flagged, never given a faked value."""
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
        w = self.width()
        act = stress_activity(f)
        ctype = f.get("cell_type", "msc")
        cname = {"msc": "MSC", "rbc": "red blood cell", "tcell": "T cell",
                 "platelet": "platelet"}.get(ctype, ctype)
        q.setPen(QColor(C["text"])); q.setFont(QFont("", 12, QFont.Weight.Bold))
        q.drawText(QRectF(16, 10, w - 32, 22), Qt.AlignmentFlag.AlignLeft, "Cell-stress pathways")
        q.setPen(QColor(C["muted"])); q.setFont(QFont("", 8))
        q.drawText(QRectF(16, 32, w - 32, 16), Qt.AlignmentFlag.AlignLeft,
                   f"live model state · Reactome + HPA · cell type: {cname}")
        y = 60; rowh = 52; barx = 20; barw = w - 40
        for key, label, hpa_key, _tgt, sev, sim in STRESS_PATHWAYS:
            a = float(act.get(key, 0.0))
            hp = HPA.get(hpa_key); rc = REACTOME.get(key)
            applic = bool(sim) and pathway_applicable(key, ctype)
            if not applic and sim:              # cell type lacks the machinery
                comp = f"not applicable — {cname} lacks this compartment"
            elif rc:                            # membership + HPA distribution, cell-filtered
                dist = [(k, n) for k, n in HPA_COMPARTMENTS.get(key, [])
                        if cell_has_compartment(ctype, k)][:3]
                dtxt = " · ".join(f"{_COMP_SHORT.get(k, k)} {n}" for k, n in dist)
                comp = f"{rc['count']} genes → {dtxt}" if dtxt else f"{rc['count']} genes · {rc['id']}"
            elif sim:
                comp = f"{hp[0] if hp else ''} · derived from model state"
            else:
                comp = f"{hp[0]}" if hp else ""
            q.setPen(QColor(C["text"] if applic else C["muted"]))
            q.setFont(QFont("", 10, QFont.Weight.Bold))
            q.drawText(QRectF(barx, y, barw - 60, 16), Qt.AlignmentFlag.AlignLeft, label)
            q.setPen(QColor(C["muted"])); q.setFont(QFont("", 8))
            q.drawText(QRectF(barx, y + 17, barw - 60, 14), Qt.AlignmentFlag.AlignLeft, comp)
            by = y + 34
            q.setBrush(QColor(C["surface2"])); q.setPen(QColor(C["border"]))
            q.drawRoundedRect(QRectF(barx, by, barw, 8), 4, 4)
            if applic:
                col = QColor(C[sev]) if a > 0.45 else QColor(C["warn"]) if a > 0.2 else QColor(C["good"])
                q.setBrush(col); q.setPen(Qt.PenStyle.NoPen)
                q.drawRoundedRect(QRectF(barx, by, max(2.0, barw * a), 8), 4, 4)
                q.setPen(QColor(C["text2"])); q.setFont(QFont("", 9, QFont.Weight.Bold))
                q.drawText(QRectF(barx, y, barw, 16), Qt.AlignmentFlag.AlignRight, f"{a*100:.0f}%")
            else:
                q.setPen(QColor(C["muted"])); q.setFont(QFont("", 8, QFont.Weight.Bold))
                q.drawText(QRectF(barx, y, barw, 16), Qt.AlignmentFlag.AlignRight,
                           "not simulated" if not sim else "n/a")
            y += rowh
        q.setPen(QColor(C["muted"])); q.setFont(QFont("", 7))
        q.drawText(QRectF(16, y + 4, w - 32, 42), Qt.AlignmentFlag.AlignLeft | Qt.TextFlag.TextWordWrap,
                   "Activity = simulated state (ROS, ΔΨm/MPT, caspase, RhoA-ROCK, protein "
                   "denaturation, pore, Ca²⁺). Gene membership from Reactome (CC0); the arrow "
                   "shows how those genes distribute across compartments (Human Protein Atlas, "
                   f"CC BY-SA, per-gene main location, retrieved {HPA_RETRIEVED}).")


class CompartmentAtlas(QWidget):
    """HPA-style isolate view: each subcellular compartment shown on its own in a
    fluorescence-micrograph style (target in green, nucleus in blue as a DAPI
    reference), like the Human Protein Atlas Cell Atlas gallery. A reference view
    of the compartments the model resolves, with their HPA gene counts."""
    # (key, style) in a sensible reading order; matches the drawn compartments
    ATLAS = [
        ("plasma_mem", "ring"), ("nucleoplasm", "nucfill"), ("nucleoli", "nucspots"),
        ("nuclear_mem", "nucring"), ("mitochondria", "tubules"), ("er", "network"),
        ("golgi", "ribbon"), ("microtubules", "radial"), ("actin", "cortex"),
        ("interm_fil", "wavy"), ("centrosome", "dot"), ("focal_adh", "edge"),
        ("lysosomes", "puncta"), ("peroxisomes", "puncta_s"), ("endosomes", "puncta"),
        ("vesicles", "puncta_m"), ("lipid_drop", "droplets"), ("cytosol", "diffuse"),
    ]
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumWidth(360); self.setAutoFillBackground(True)

    def paintEvent(self, _):
        q = QPainter(self); q.setRenderHint(QPainter.RenderHint.Antialiasing)
        q.fillRect(self.rect(), QColor(C["surface0"]))
        w = self.width()
        q.setPen(QColor(C["text"])); q.setFont(QFont("", 12, QFont.Weight.Bold))
        q.drawText(QRectF(14, 8, w - 28, 22), Qt.AlignmentFlag.AlignLeft, "Compartment atlas")
        q.setPen(QColor(C["muted"])); q.setFont(QFont("", 8))
        q.drawText(QRectF(14, 30, w - 28, 14), Qt.AlignmentFlag.AlignLeft,
                   "each compartment isolated · HPA immunofluorescence style · gene counts from HPA")
        cols = max(2, int(w // 190)); tw = (w - 20) / cols; th = tw * 0.82 + 34
        top = 50
        for i, (key, style) in enumerate(self.ATLAS):
            cx = 10 + (i % cols) * tw; cy = top + (i // cols) * th
            self._tile(q, cx, cy, tw - 8, th - 10, key, style)
        rows = (len(self.ATLAS) + cols - 1) // cols
        self.setMinimumHeight(int(top + rows * th + 12))

    def _tile(self, q, x, y, w, h, key, style):
        img_h = h - 26
        BG = QColor(8, 12, 18); NUC = QColor(74, 120, 226); GRN = QColor(72, 224, 132)
        REF = QColor(150, 160, 172, 70)
        q.setBrush(BG); q.setPen(QPen(QColor(C["border"]), 1)); q.drawRoundedRect(QRectF(x, y, w, img_h), 6, 6)
        q.save(); q.setClipRect(QRectF(x, y, w, img_h))
        cx, cy = x + w / 2, y + img_h / 2; R = min(w, img_h) * 0.40
        rs = np.random.RandomState(hash(key) & 0xffff)
        nuc_target = style in ("nucfill", "nucspots", "nucring")
        # membrane reference
        q.setPen(QPen(REF, 1.4)); q.setBrush(Qt.BrushStyle.NoBrush)
        q.drawEllipse(QPointF(cx, cy), R, R)
        # nucleus (blue DAPI reference), unless the nucleus IS the target
        if not nuc_target:
            q.setBrush(QColor(NUC.red(), NUC.green(), NUC.blue(), 90)); q.setPen(Qt.PenStyle.NoPen)
            q.drawEllipse(QPointF(cx, cy), R * 0.42, R * 0.42)
        def pts(n, rmax, rmin=0.0):
            for _ in range(n):
                a = rs.uniform(0, 2 * math.pi); rr = R * (rmin + (rmax - rmin) * math.sqrt(rs.uniform(0, 1)))
                yield cx + math.cos(a) * rr, cy + math.sin(a) * rr
        q.setPen(Qt.PenStyle.NoPen); q.setBrush(GRN)
        if style == "ring":
            q.setBrush(Qt.BrushStyle.NoBrush); q.setPen(QPen(GRN, 3)); q.drawEllipse(QPointF(cx, cy), R, R)
        elif style == "nucfill":
            q.setBrush(GRN); q.drawEllipse(QPointF(cx, cy), R * 0.5, R * 0.5)
        elif style == "nucspots":
            q.setBrush(QColor(NUC.red(), NUC.green(), NUC.blue(), 80)); q.drawEllipse(QPointF(cx, cy), R * 0.5, R * 0.5)
            q.setBrush(GRN)
            for px, py in pts(3, 0.32): q.drawEllipse(QPointF(px, py), R * 0.10, R * 0.10)
        elif style == "nucring":
            q.setBrush(QColor(NUC.red(), NUC.green(), NUC.blue(), 70)); q.drawEllipse(QPointF(cx, cy), R * 0.5, R * 0.5)
            q.setBrush(Qt.BrushStyle.NoBrush); q.setPen(QPen(GRN, 2.6)); q.drawEllipse(QPointF(cx, cy), R * 0.5, R * 0.5)
        elif style == "tubules":
            for px, py in pts(9, 0.9, 0.5):
                q.save(); q.translate(px, py); q.rotate(rs.uniform(0, 180))
                q.setBrush(GRN); q.drawEllipse(QPointF(0, 0), R * 0.16, R * 0.05); q.restore()
        elif style == "network":
            q.setPen(QPen(GRN, 1.6)); q.setBrush(Qt.BrushStyle.NoBrush)
            for k in range(5):
                p = QPainterPath()
                for s in range(20):
                    th_ = k * 0.7 + s * 0.09; rr = R * (0.5 + 0.28 * math.sin(s * 0.8 + k))
                    pt = QPointF(cx + math.cos(th_) * rr, cy + math.sin(th_) * rr)
                    p.moveTo(pt) if s == 0 else p.lineTo(pt)
                q.drawPath(p)
        elif style == "ribbon":
            gx, gy = cx + R * 0.5, cy - R * 0.15
            q.setPen(QPen(GRN, 2)); q.setBrush(Qt.BrushStyle.NoBrush)
            for k in range(4):
                p = QPainterPath(QPointF(gx - R * 0.28, gy + k * 4 - 6))
                p.quadTo(QPointF(gx, gy + k * 4 - 14), QPointF(gx + R * 0.28, gy + k * 4 - 6)); q.drawPath(p)
        elif style == "radial":
            q.setPen(QPen(GRN, 1.5))
            for k in range(14):
                a = k / 14 * 2 * math.pi
                q.drawLine(QPointF(cx + math.cos(a) * R * 0.18, cy + math.sin(a) * R * 0.18),
                           QPointF(cx + math.cos(a) * R * 0.95, cy + math.sin(a) * R * 0.95))
        elif style == "cortex":
            q.setPen(QPen(GRN, 3)); q.setBrush(Qt.BrushStyle.NoBrush); q.drawEllipse(QPointF(cx, cy), R * 0.92, R * 0.92)
            q.setPen(QPen(GRN, 1.4))
            for k in range(5):
                a = rs.uniform(0, 2 * math.pi)
                q.drawLine(QPointF(cx + math.cos(a) * R * 0.3, cy + math.sin(a) * R * 0.3),
                           QPointF(cx + math.cos(a + 0.4) * R * 0.85, cy + math.sin(a + 0.4) * R * 0.85))
        elif style == "wavy":
            q.setPen(QPen(GRN, 1.5)); q.setBrush(Qt.BrushStyle.NoBrush)
            for k in range(6):
                a0 = k * 1.0; p = QPainterPath(QPointF(cx + math.cos(a0) * R * 0.9, cy + math.sin(a0) * R * 0.9))
                for s in range(1, 12):
                    t = s / 11; rr = R * (0.9 - 0.8 * t); aa = a0 + math.sin(t * 6) * 0.4
                    p.lineTo(QPointF(cx + math.cos(aa) * rr, cy + math.sin(aa) * rr))
                q.drawPath(p)
        elif style == "dot":
            q.setBrush(GRN); q.drawEllipse(QPointF(cx + R * 0.25, cy - R * 0.2), R * 0.12, R * 0.12)
        elif style == "edge":
            q.setBrush(GRN)
            for k in range(6):
                a = math.pi * (0.25 + k / 5 * 0.5)
                q.drawEllipse(QPointF(cx + math.cos(a) * R * 0.95, cy + math.sin(a) * R * 0.95), R * 0.08, R * 0.05)
        elif style in ("puncta", "puncta_s", "puncta_m", "droplets"):
            n = {"puncta": 12, "puncta_s": 16, "puncta_m": 22, "droplets": 7}[style]
            sz = {"puncta": 0.09, "puncta_s": 0.06, "puncta_m": 0.07, "droplets": 0.12}[style]
            for px, py in pts(n, 0.9, 0.45):
                if style == "droplets":
                    q.setBrush(QColor(230, 200, 90)); q.setPen(QPen(GRN, 1.4))
                else:
                    q.setBrush(GRN); q.setPen(Qt.PenStyle.NoPen)
                q.drawEllipse(QPointF(px, py), R * sz, R * sz)
        elif style == "diffuse":
            q.setBrush(QColor(GRN.red(), GRN.green(), GRN.blue(), 70)); q.setPen(Qt.PenStyle.NoPen)
            q.drawEllipse(QPointF(cx, cy), R * 0.95, R * 0.95)
            q.setBrush(BG); q.drawEllipse(QPointF(cx, cy), R * 0.42, R * 0.42)
        q.restore()
        # label
        hp = HPA.get(key)
        nm = hp[0] if hp else key
        genes = f"{hp[1]:,} genes · {hp[2]}%" if hp and hp[1] is not None else ""
        q.setPen(QColor(C["text"])); q.setFont(QFont("", 8, QFont.Weight.Bold))
        q.drawText(QRectF(x, y + img_h + 1, w, 12), Qt.AlignmentFlag.AlignHCenter, nm)
        q.setPen(QColor(C["muted"])); q.setFont(QFont("", 7))
        q.drawText(QRectF(x, y + img_h + 13, w, 11), Qt.AlignmentFlag.AlignHCenter, genes)


# ---------------------------------------------------------------------------
# 3D spheroid / multicellular-construct cryopreservation view.
#
# The core single-cell model is run once; this view overlays the physics that
# make a 3D construct behave differently from a cell in suspension, resolved
# ACROSS THE RADIUS of the spheroid:
#   1. CPA loading is diffusion-limited -> the core stays under-protected in a
#      fixed loading hold, and the bigger the spheroid the worse it gets
#      (penetration depth delta ~ sqrt(D_eff * t_hold)).
#   2. Intracellular ice, once nucleated, propagates cell-to-cell through gap
#      junctions -> a single event sweeps the coupled cluster (Irimia & Karlsson
#      2002; Acker & McGann 2000). The model's junction-coupled P_iif drives the
#      invasion front here.
#   3. Large spheroids carry a pre-existing hypoxic / necrotic core before any
#      freezing (a standard 3D-culture feature, not a cryo effect) -> flagged
#      separately so it is not confused with freezing damage.
#
# This is an illustrative radial overlay on the single-cell solve, not a full 3D
# reaction-diffusion solve. The diffusion coefficient is an order-of-magnitude
# tissue value and is labelled as such. Each drawn circle is a representative
# cell for its shell, not a literal count.
_D_EFF_CPA = 3.0e-11        # effective CPA diffusivity in packed tissue, m^2/s
                            # (tortuous; ~10x lower than free solution). Order of
                            # magnitude only -- Xu 2014; Devireddy tissue reviews.
# Measured spheroid nucleation temperatures, Gao, Bissoyi, Guo & Gibson 2024
# (ACS Biomater Sci Eng 11:208, doi:10.1021/acsbiomaterials.4c00958):
TN_INP = -9.25              # 10% DMSO + extracellular ice nucleator (IN+)
TN_DMSO = -15.77            # 10% DMSO alone (supercools before nucleating)

def _hex_axial(rings):
    """Axial (q, r) coords for a hex-packed disk of `rings` rings around 0,0."""
    cells = [(0, 0)]
    dirs = [(1, 0), (1, -1), (0, -1), (-1, 0), (-1, 1), (0, 1)]
    for rad in range(1, rings + 1):
        q, r = -rad, rad                    # start corner
        for i in range(6):
            for _ in range(rad):
                cells.append((q, r))
                q += dirs[i][0]; r += dirs[i][1]
    return cells

_HEX_NEI = [(1, 0), (1, -1), (0, -1), (-1, 0), (-1, 1), (0, 1)]


class _SpheroidCanvas(QWidget):
    """Custom-painted radial cross-section of a cryopreserved spheroid."""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.frame = None
        self.diam_um = 300.0
        self.in_plus = False                                     # extracellular ice nucleator?
        self.setMinimumSize(420, 460)
        self.setAutoFillBackground(True)
        self._layout_cache = {}

    def set_frame(self, f):
        self.frame = f; self.update()

    def set_diam(self, d):
        self.diam_um = float(d); self.update()

    def set_in_plus(self, on):
        self.in_plus = bool(on); self.update()

    # ---- physics helpers ---------------------------------------------------
    def _penetration_um(self, f):
        """CPA penetration depth from the surface, microns.

        Equilibration is set during the warm loading hold: diffusion is
        Arrhenius-slow once cooling starts and stops when the construct freezes,
        so whatever the core failed to take up during the hold stays missing
        into the freeze. During the hold itself the front is still advancing, so
        it grows with elapsed time up to the full hold; after that it is fixed.
        """
        t_hold = max(f.get("hold_min", 10.0), 0.1) * 60.0        # s
        if f.get("phase") == "load":
            t_eff = clamp(f.get("t", t_hold), 0.0, t_hold)       # still loading
        else:
            t_eff = t_hold                                       # frozen at hold value
        return math.sqrt(_D_EFF_CPA * t_eff) * 1e6              # m -> um

    def _cpa_frac(self, u, f):
        """Local CPA fraction (0..1) at fractional radius u (0 core .. 1 rim)."""
        loaded = clamp(f.get("cpaLoad", 0.0) / 0.6, 0.0, 1.0)    # any CPA on board?
        if loaded <= 0.01:
            return 0.0
        R = self.diam_um / 2.0
        depth = R * (1.0 - u)                                    # um from surface
        delta = self._penetration_um(f)
        if depth <= delta:
            local = 1.0
        else:
            local = math.exp(-(depth - delta) / max(0.6 * delta, 1e-6))
        return clamp(loaded * local, 0.0, 1.0)

    def _necrotic(self, u, f):
        """Pre-existing hypoxic/necrotic core fraction (0..1), size-dependent.
        A standard 3D-culture feature (before any freezing), so it is shown
        separately from cryo-damage. Onset in large spheroids only."""
        big = clamp((self.diam_um - 350.0) / 250.0, 0.0, 1.0)    # onset well above 200 um
        core = clamp((0.30 - u) / 0.30, 0.0, 1.0)
        return big * core

    def _nucleation_T(self):
        """Effective nucleation temperature, degC. Measured values from Gao,
        Bissoyi, Guo & Gibson 2024 (ACS Biomater Sci Eng 11:208): 10% DMSO alone
        supercools and nucleates at -15.77 degC; adding an extracellular ice
        nucleator (IN+) raises it to -9.25 degC."""
        return TN_INP if self.in_plus else TN_DMSO

    def _supercool_sev(self):
        """Supercooling severity 0..1. IN+ (-9.25 C) ~ 0; DMSO-only (-15.77 C) ~ 1.
        Deeper supercooling before nucleation is the dominant damage driver for
        spheroids (delayed nucleation -> longer CPA exposure, less dehydration,
        violent freezing) -- Gao/Bissoyi 2024; Gao/Bissoyi 2023 Chem Commun."""
        Tn = self._nucleation_T()
        return clamp((-Tn - (-TN_INP)) / ((-TN_DMSO) - (-TN_INP)), 0.0, 1.0)

    def _sizefac(self):
        """Size penalty. 200 um spheroids recover best; 400 um markedly worse
        (Gao/Bissoyi 2024, measured WST-1 recovery)."""
        return clamp((self.diam_um - 180.0) / 220.0, 0.0, 1.0)

    def _shed(self, u, f):
        """Surface-shedding damage at shell u. Without induced nucleation the
        outermost cells detach from the spheroid -- the paper's key finding, and
        the same detachment seen in supercooled 2D monolayers. Confined to the
        outer ~two shells; scales with supercooling and (weakly) size."""
        ss = self._supercool_sev()
        rim = clamp((u - 0.72) / 0.28, 0.0, 1.0)                 # outer shells only
        return clamp(ss * rim * (0.7 + 0.3 * self._sizefac()), 0.0, 1.0)

    def _perforate(self, u, f):
        """Interior perforation at shell u. Without induced nucleation the core
        shows voids/perforation (Gao/Bissoyi 2024) -- pronounced in large
        spheroids, minimal in small ones. Confined to the inner shells."""
        ss = self._supercool_sev()
        core = clamp((0.62 - u) / 0.62, 0.0, 1.0)               # inner shells
        return clamp(ss * core * (0.35 + 0.65 * self._sizefac()), 0.0, 1.0)

    def _dead(self, u, f):
        """Post-thaw lethality at shell u: shed rim OR perforated interior."""
        return clamp(max(self._shed(u, f), self._perforate(u, f)), 0.0, 1.0)

    # ---- painting ----------------------------------------------------------
    def paintEvent(self, _):
        q = QPainter(self); q.setRenderHint(QPainter.RenderHint.Antialiasing)
        q.fillRect(self.rect(), QColor(9, 12, 17))
        w, h = self.width(), self.height()
        f = self.frame
        if not f:
            q.setPen(QColor(150, 160, 172)); q.setFont(QFont("", 10))
            q.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, "run a protocol to view the spheroid")
            return

        rings = 4
        cells = self._layout_cache.get(rings)
        if cells is None:
            cells = _hex_axial(rings); self._layout_cache[rings] = cells
        cellset = set(cells)
        post_thaw = f.get("phase") in ("melt", "dilute", "recover", "end") and f.get("T", 20) > 0.5
        # pixel geometry — keep the cluster in the left ~58% so the readout column
        # on the right never overlaps the cells. A flat-top hex disk of `rings`
        # rings spans ~2*sqrt(3)*rings*size, so normalise `size` by that.
        cx, cy = w * 0.31, h * 0.50
        Rpix = min(w * 0.28, h * 0.42)
        size = Rpix / (math.sqrt(3) * rings + 0.8)              # hex spacing (fits Rpix)
        crad = size * 0.92
        def px(qc, rc):
            x = cx + size * 1.5 * qc
            y = cy + size * math.sqrt(3) * (rc + qc / 2.0)
            return x, y
        def udist(qc, rc):                                       # cube distance / rings
            xx, zz = qc, rc; yy = -xx - zz
            return (abs(xx) + abs(yy) + abs(zz)) / 2.0 / rings

        # protocol state: is ice present yet, and (if supercooling) has the medium
        # reached its nucleation temperature?
        T = f.get("T", 20.0); phase = f.get("phase", "")
        Tn = self._nucleation_T()
        cold = T < -0.5
        nucleated = cold and (T <= Tn + 0.5 or phase in ("seed", "store", "warm", "melt"))
        supercooled = cold and not nucleated                     # below 0, not yet frozen
        post_thaw = phase in ("melt", "dilute", "recover", "end") and T > 0.5
        damaged = nucleated or post_thaw

        # junction lines first (behind cells)
        q.setPen(QPen(QColor(120, 130, 145, 90), 1.2))
        for (qc, rc) in cells:
            x0, y0 = px(qc, rc)
            for dq, dr in _HEX_NEI:
                nb = (qc + dq, rc + dr)
                if nb in cellset and (nb > (qc, rc)):
                    x1, y1 = px(*nb); q.drawLine(QPointF(x0, y0), QPointF(x1, y1))

        R = self.diam_um / 2.0
        # per-cell draw
        alive_n = under_n = shed_n = perf_n = necr_n = 0
        for (qc, rc) in cells:
            x, y = px(qc, rc); u = udist(qc, rc)
            cpa = self._cpa_frac(u, f); nec = self._necrotic(u, f)
            shed = self._shed(u, f); perf = self._perforate(u, f)
            g = QRadialGradient(QPointF(x - crad * 0.3, y - crad * 0.3), crad * 1.4)
            has_nuc = True; dx = dy = 0.0; alpha = 255
            if nec > 0.5:                                        # pre-existing hypoxic core
                g.setColorAt(0, QColor(90, 74, 60)); g.setColorAt(1, QColor(54, 44, 36))
                ring_c = QColor(120, 96, 70); has_nuc = False; necr_n += 1
            elif damaged and shed > 0.5:                         # surface cell shed / detached
                g.setColorAt(0, QColor(180, 80, 80)); g.setColorAt(1, QColor(110, 40, 40))
                ring_c = QColor(200, 100, 100); has_nuc = False; shed_n += 1
                ang = math.atan2(y - cy, x - cx); off = crad * (0.25 + 0.45 * shed)
                dx, dy = math.cos(ang) * off, math.sin(ang) * off  # detaching outward
                alpha = 150
            elif damaged and perf > 0.5:                         # interior perforation / void
                g.setColorAt(0, QColor(70, 40, 44)); g.setColorAt(1, QColor(34, 20, 22))
                ring_c = QColor(120, 70, 74); has_nuc = False; perf_n += 1
            elif supercooled:                                    # below 0 C, not yet nucleated
                g.setColorAt(0, QColor(150, 190, 220)); g.setColorAt(1, QColor(80, 120, 160))
                ring_c = QColor(170, 205, 230)
            elif nucleated and self.in_plus:                     # gently frozen, protected
                g.setColorAt(0, QColor(120, 210, 210)); g.setColorAt(1, QColor(60, 150, 160))
                ring_c = QColor(150, 220, 220); alive_n += 1
            elif cpa >= 0.66:                                    # alive; hue by CPA protection
                g.setColorAt(0, QColor(90, 220, 150)); g.setColorAt(1, QColor(40, 150, 110))
                ring_c = QColor(120, 230, 170); alive_n += 1
            elif cpa >= 0.33:
                g.setColorAt(0, QColor(235, 200, 90)); g.setColorAt(1, QColor(180, 140, 50))
                ring_c = QColor(240, 210, 120); under_n += 1
            else:
                g.setColorAt(0, QColor(232, 130, 70)); g.setColorAt(1, QColor(170, 80, 45))
                ring_c = QColor(240, 150, 100); under_n += 1
            if alpha < 255:
                ring_c = QColor(ring_c.red(), ring_c.green(), ring_c.blue(), alpha)
            q.setOpacity(alpha / 255.0)
            q.setBrush(g); q.setPen(QPen(ring_c, 1.4))
            q.drawEllipse(QPointF(x + dx, y + dy), crad, crad)
            if has_nuc:                                          # nucleus dot for living cells
                q.setBrush(QColor(60, 90, 170, 150)); q.setPen(Qt.PenStyle.NoPen)
                q.drawEllipse(QPointF(x + dx, y + dy), crad * 0.34, crad * 0.34)
            q.setOpacity(1.0)

        # CPA penetration ring marker (how deep the front reached) — during loading
        delta = self._penetration_um(f)
        if f.get("cpaLoad", 0.0) > 0.02 and delta < R and not damaged and not supercooled:
            rr = Rpix * clamp(1.0 - delta / R, 0.0, 1.0)
            q.setBrush(Qt.BrushStyle.NoBrush)
            q.setPen(QPen(QColor(90, 220, 150, 150), 1.6, Qt.PenStyle.DashLine))
            q.drawEllipse(QPointF(cx, cy), rr, rr)
            q.setPen(QColor(120, 230, 170)); q.setFont(QFont("", 7))
            q.drawText(QRectF(cx - 60, cy - rr - 13, 120, 12),
                       Qt.AlignmentFlag.AlignHCenter, "CPA front")

        # supercooling warning banner (the invisible danger)
        if supercooled and not self.in_plus:
            q.setPen(QColor(150, 190, 220)); q.setFont(QFont("", 8, QFont.Weight.Bold))
            q.drawText(QRectF(cx - Rpix, cy - Rpix - 20, 2 * Rpix, 14),
                       Qt.AlignmentFlag.AlignHCenter,
                       f"supercooled to {T:.1f} °C — ice not yet nucleated")

        # scale bar (true selected diameter)
        bar_um = 50.0 if self.diam_um <= 250 else 100.0
        blen = bar_um * (Rpix / R)
        bx, by = cx - Rpix, cy + Rpix + 18
        q.setPen(QPen(QColor(210, 214, 220), 2)); q.drawLine(QPointF(bx, by), QPointF(bx + blen, by))
        q.setFont(QFont("", 8)); q.setPen(QColor(210, 214, 220))
        q.drawText(QRectF(bx, by + 3, blen + 60, 14), Qt.AlignmentFlag.AlignLeft, f"{bar_um:.0f} µm")

        self._draw_readout(q, w, h, f, Tn, alive_n, shed_n, perf_n, necr_n, len(cells))

    def _draw_readout(self, q, w, h, f, Tn, alive, shed, perf, necr, total):
        x0 = min(w * 0.60, w - 220.0)
        y = 14
        def line(txt, col=QColor(210, 214, 220), sz=8, bold=False, dy=13):
            nonlocal y
            q.setPen(col); q.setFont(QFont("", sz, QFont.Weight.Bold if bold else QFont.Weight.Normal))
            q.drawText(QRectF(x0, y, w - x0 - 6, 16), Qt.AlignmentFlag.AlignLeft, txt)
            y += dy
        line(f"Spheroid {self.diam_um:.0f} µm", QColor(240, 240, 240), 10, True, 17)
        line(f"{f.get('phase','')} · {f.get('T',0):.1f} °C", QColor(160, 170, 182), 8, False, 16)
        # nucleation regime
        if self.in_plus:
            line("IN+ ice nucleator", QColor(150, 220, 220), 9, True)
            line(f"  nucleates at {Tn:.1f} °C", QColor(150, 200, 210))
        else:
            line("DMSO only", QColor(235, 170, 120), 9, True)
            line(f"  supercools to {Tn:.1f} °C", QColor(235, 170, 120))
        y += 3
        line(f"Cryo-damage (of {total} cells)", QColor(240, 240, 240), 9, True)
        line(f"  surface shed: {shed}", QColor(210, 110, 110) if shed else QColor(150, 160, 172))
        line(f"  interior voids: {perf}", QColor(180, 110, 110) if perf else QColor(150, 160, 172))
        line(f"  intact: {alive}", QColor(120, 230, 170))
        if necr:
            line(f"  hypoxic core: {necr}", QColor(180, 150, 120))
        y += 3
        line("Radial CPA load", QColor(240, 240, 240), 9, True)
        cc, rc = self._cpa_frac(0.0, f), self._cpa_frac(1.0, f)
        line(f"  core: {cc*100:.0f}%   rim: {rc*100:.0f}%",
             QColor(120, 230, 170) if cc > 0.6 else QColor(240, 200, 110))
        y += 5
        line("Legend", QColor(240, 240, 240), 9, True)
        for c, t in [(QColor(60, 185, 130), "intact / protected"),
                     (QColor(235, 200, 90), "under-loaded core"),
                     (QColor(180, 80, 80), "shed surface cell"),
                     (QColor(70, 40, 44), "interior void"),
                     (QColor(150, 190, 220), "supercooled"),
                     (QColor(90, 74, 60), "hypoxic core")]:
            q.setBrush(c); q.setPen(Qt.PenStyle.NoPen)
            q.drawEllipse(QPointF(x0 + 6, y + 5), 5, 5)
            q.setPen(QColor(200, 205, 212)); q.setFont(QFont("", 8))
            q.drawText(QRectF(x0 + 16, y, w - x0 - 20, 14), Qt.AlignmentFlag.AlignLeft, t)
            y += 15
        y += 3
        q.setPen(QColor(130, 138, 150)); q.setFont(QFont("", 7))
        q.drawText(QRectF(x0, y, w - x0 - 4, 120), Qt.TextFlag.TextWordWrap,
                   "Damage pattern from Gao, Bissoyi, Guo & Gibson 2024 (ACS Biomater "
                   "Sci Eng 11:208): DMSO-only supercooling shears surface cells and "
                   "perforates the interior; an extracellular nucleator (IN+) raises the "
                   "nucleation temperature and protects both. Illustrative radial overlay "
                   "on the single-cell solve, not a 3D reaction-diffusion solve.")


class SpheroidView(QWidget):
    """3D spheroid cryopreservation view with a live diameter selector and an
    extracellular-ice-nucleator (IN+) toggle."""
    def __init__(self, parent=None):
        super().__init__(parent)
        v = QVBoxLayout(self); v.setContentsMargins(8, 8, 8, 8); v.setSpacing(6)
        head = QHBoxLayout()
        t = QLabel("3D spheroid cryopreservation"); t.setStyleSheet("font-weight:600;")
        head.addWidget(t); head.addStretch(1)
        head.addWidget(QLabel("diameter"))
        self.diam = QComboBox()
        for d in (120, 200, 300, 400, 500):
            self.diam.addItem(f"{d} µm", d)
        self.diam.setCurrentIndex(3)                            # 400 um (paper's large case)
        self.diam.currentIndexChanged.connect(self._diam)
        head.addWidget(self.diam)
        self.inp = QCheckBox("extracellular ice nucleator (IN+)")
        self.inp.toggled.connect(self._inp)
        head.addWidget(self.inp)
        v.addLayout(head)
        sub = QLabel("Follows Gao, Bissoyi, Guo & Gibson 2024 (ACS Biomater Sci Eng 11:208): "
                     "10% DMSO alone supercools to ~-16 °C, shedding surface cells and "
                     "perforating the interior; an extracellular ice nucleator (IN+) raises "
                     "nucleation to ~-9 °C and protects the whole spheroid. Smaller spheroids "
                     "recover better. Toggle IN+ and change the diameter to show it live.")
        sub.setWordWrap(True); sub.setStyleSheet("color:#8a8873; font-size:11px;")
        v.addWidget(sub)
        self.canvas = _SpheroidCanvas()
        self.canvas.set_diam(self.diam.currentData())
        v.addWidget(self.canvas, 1)

    def _diam(self, _):
        self.canvas.set_diam(self.diam.currentData())

    def _inp(self, on):
        self.canvas.set_in_plus(on)

    def set_frame(self, f):
        self.canvas.set_frame(f)


class _DropdownPanel(QWidget):
    """A tab-like container that picks pages from a dropdown instead of a tab
    strip — more ergonomic when there are many views (no overflow / scroll
    arrows). Views can be grouped under non-selectable section headers. Mimics
    the small QTabWidget API used here (addTab / setCurrentIndex / currentIndex);
    indices are stack-page indices, not combo rows."""
    def __init__(self, parent=None):
        super().__init__(parent)
        lay = QVBoxLayout(self); lay.setContentsMargins(0, 0, 0, 0); lay.setSpacing(6)
        row = QHBoxLayout()
        lbl = QLabel("View"); lbl.setStyleSheet("font-weight:600; color:#52514e;")
        self.sel = QComboBox(); self.sel.setMinimumWidth(200)
        self._model = QStandardItemModel(self.sel); self.sel.setModel(self._model)
        row.addWidget(lbl); row.addWidget(self.sel, 1)
        lay.addLayout(row)
        self.stack = QStackedWidget()
        lay.addWidget(self.stack, 1)
        self._last_group = None; self._first_row = None
        self.sel.currentIndexChanged.connect(self._on_sel)

    def addTab(self, w, name, group=None):
        page = self.stack.count(); self.stack.addWidget(w)
        if group and group != self._last_group:            # non-selectable header
            hdr = QStandardItem(group.upper())
            hdr.setFlags(Qt.ItemFlag.NoItemFlags)
            hf = QFont(); hf.setBold(True); hf.setPointSize(8)
            hdr.setData(hf, Qt.ItemDataRole.FontRole)
            hdr.setData(QColor("#8a8873"), Qt.ItemDataRole.ForegroundRole)
            self._model.appendRow(hdr); self._last_group = group
        it = QStandardItem(("   " + name) if group else name)
        it.setData(page, Qt.ItemDataRole.UserRole)
        self._model.appendRow(it)
        if self._first_row is None:
            self._first_row = it.row(); self.sel.setCurrentIndex(it.row())
        return page

    def _on_sel(self, row):
        it = self._model.item(row)
        if it is None: return
        pg = it.data(Qt.ItemDataRole.UserRole)
        if pg is not None: self.stack.setCurrentIndex(pg)

    def setCurrentIndex(self, page):                        # page = stack index
        for r in range(self._model.rowCount()):
            it = self._model.item(r)
            if it and it.data(Qt.ItemDataRole.UserRole) == page:
                self.sel.setCurrentIndex(r); return

    def currentIndex(self):
        it = self._model.item(self.sel.currentIndex())
        return it.data(Qt.ItemDataRole.UserRole) if it else 0


class Main(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("CryoCell — a virtual cell for cryoprotectant testing")
        self.resize(1680, 1000)
        self.P = Params()
        self.S = self.R = None
        self.ref = None                 # saved reference run for A/B comparison
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
        def caption(attr):
            txt = HELP.get(attr)
            if not txt: return
            c = QLabel(txt); c.setWordWrap(True)
            c.setStyleSheet("color:#8a8873; font-size:10px; margin-bottom:4px;")
            g.addWidget(c)
        for attr, lab, options in COMBOS:
            g.addWidget(QLabel(lab))
            cb = QComboBox()
            for val, txt in options: cb.addItem(txt, val)
            cb.setCurrentIndex(max(0, [o[0] for o in options].index(getattr(self.P, attr))
                                   if getattr(self.P, attr) in [o[0] for o in options] else 0))
            cb.currentIndexChanged.connect(lambda _i, a=attr, c=cb: self._set(a, c.currentData()))
            self.widgets[attr] = cb; g.addWidget(cb)
            caption(attr)
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
            caption(attr)
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

    # ------------------------------------------------- A/B reference comparison
    def _save_ref(self):
        if not self.R: return
        R = self.R
        lab = self.preset.currentText()
        if lab in ("", "Load a preset…"):
            lab = f"{CPAS[self.P.cpa_key]['name']} {self.P.conc_pct:.0f}%"
        self.ref = dict(label=lab, S_24=R['S_24'], F_rec=R['F_rec'],
                        P_iif=R['P_iif'], maxV=R['maxV'])
        self.refclr.setEnabled(True)
        self._update_compare()

    def _clear_ref(self):
        self.ref = None; self.refclr.setEnabled(False); self.cmp.setText("")

    def _update_compare(self):
        if not self.ref or not self.R:
            self.cmp.setText(""); return
        R, ref = self.R, self.ref
        def chip(cur, base, better_up, unit, mul=1.0, dec=0):
            cv, bv = cur * mul, base * mul; dv = cv - bv
            if abs(dv) < 10 ** (-dec) / 2:      col = C['muted']
            elif (dv > 0) == better_up:         col = C['good']
            else:                               col = C['crit']
            sign = '+' if dv >= 0 else '−'
            return (f"{cv:.{dec}f}{unit} <span style=\"color:{col}\">"
                    f"(Δ{sign}{abs(dv):.{dec}f}{unit})</span>")
        self.cmp.setText(
            f"<b>vs reference</b> · {ref['label']}<br>"
            f"24 h viable {chip(R['S_24'], ref['S_24'], True, '%', 100, 0)} &nbsp; "
            f"recovery {chip(R['F_rec'], ref['F_rec'], True, '%', 100, 0)}<br>"
            f"P(IIF) {chip(R['P_iif'], ref['P_iif'], False, '', 1, 2)} &nbsp; "
            f"max swell {chip(R['maxV'], ref['maxV'], False, '×', 1, 2)}")

    def _preset(self, name):
        p = PRESETS.get(name)
        if not p: return
        # Cell-identity and mode fields must NOT persist across presets, or a cell
        # picked after RBC/platelet would stay anucleate (or freeze-dried). Reset
        # them to defaults unless the chosen preset sets them explicitly; protocol
        # parameters (CPA, temperatures, rates) still carry over.
        d = Params()
        RESET = ("cell_type", "freeze_dry", "Viso", "lp", "ps", "sterol", "cyto", "nuc_scale",
                 "apop_resist", "anoikis_resist", "antioxidant", "glycolytic",
                 "dry_residual", "dry_hours", "drystore_days", "adhesion")
        cur = {k: v for k, v in vars(self.P).items() if k in Params.__annotations__}
        for k in RESET:
            cur[k] = getattr(d, k)
        self.P = Params(**{**cur, **p})
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
        for key, lbl in [("dark", "Dark"), ("illustrative", "Illustrative"),
                         ("fluor", "Fluorescence"), ("phase", "Phase-contrast")]:
            self.rendcombo.addItem(lbl, key)
        self.rendcombo.setToolTip("Rendering style: dark green-cytoplasm view (default), "
                                  "illustrative colour, confocal-fluorescence (glowing channels), "
                                  "or label-free phase-contrast greyscale")
        self.rendcombo.currentIndexChanged.connect(
            lambda _i: self.view.set_render_mode(self.rendcombo.currentData()))
        row.addWidget(self.rendcombo)
        self.scibox = QCheckBox("Scientific labels")
        self.scibox.setToolTip("Overlay scale bar, leader-line labels with live quantitative "
                               "state, legend, and a cryo-stage dendritic ice front (any rendering)")
        self.scibox.setChecked(True)  # scientific overlay on by default
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
        tabs = _DropdownPanel()
        # --- outcome + inspector (PINNED — always visible above the dropdown views)
        w = QWidget(); v = QVBoxLayout(w); v.setContentsMargins(0, 0, 0, 0)
        oc_l = QLabel("Outcome"); oc_l.setStyleSheet("font-weight:600; color:#52514e;")
        v.addWidget(oc_l)
        self.hero = QLabel("—"); self.hero.setFont(QFont("", 34, QFont.Weight.Bold))
        v.addWidget(self.hero)
        self.hero2 = QLabel(); self.hero2.setStyleSheet("color:#52514e"); self.hero2.setWordWrap(True)
        v.addWidget(self.hero2)
        self.quad = QLabel(); self.quad.setStyleSheet("font-family:monospace"); v.addWidget(self.quad)
        # A/B comparison: save the current run as a reference, then see the delta
        # of every later run against it — the fastest way to show cause and effect.
        cmp_row = QHBoxLayout()
        self.refbtn = QPushButton("Save as reference"); self.refbtn.clicked.connect(self._save_ref)
        self.refbtn.setToolTip("Freeze the current run as a baseline. Later runs show the change "
                               "in survival, ice risk and swelling against it.")
        self.refclr = QPushButton("Clear"); self.refclr.clicked.connect(self._clear_ref)
        self.refclr.setMaximumWidth(60); self.refclr.setEnabled(False)
        cmp_row.addWidget(self.refbtn); cmp_row.addWidget(self.refclr); cmp_row.addStretch()
        v.addLayout(cmp_row)
        self.cmp = QLabel(); self.cmp.setWordWrap(True)
        self.cmp.setStyleSheet("font-size:11px; color:#52514e;")
        self.cmp.setTextFormat(Qt.TextFormat.RichText)
        v.addWidget(self.cmp)
        self.dmg = QTableWidget(0, 2); self.dmg.horizontalHeader().setVisible(False)
        self.dmg.verticalHeader().setVisible(False)
        self.dmg.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        self.dmg.setMaximumHeight(210); v.addWidget(self.dmg)
        insp_l = QLabel("Inspector"); insp_l.setStyleSheet("font-weight:600;color:#52514e")
        v.addWidget(insp_l)
        self.insp = QTextEdit(); self.insp.setReadOnly(True); v.addWidget(self.insp, 1)
        outcome_w = w                                  # pinned, not added to the dropdown

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
        tabs.addTab(sc, "Plots", "Curves")

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
        tabs.addTab(gsc, "Cryo-signatures", "Curves")

        # --- mechanotransduction / cell-death signalling network (live)
        self.mechano = MechanoView()
        msc = QScrollArea(); msc.setWidget(self.mechano); msc.setWidgetResizable(True)
        tabs.addTab(msc, "Signalling", "Cell mechanisms")

        # --- molecular-detail focal-adhesion & LINC inset (live)
        self.molec = MolecularView()
        molsc = QScrollArea(); molsc.setWidget(self.molec); molsc.setWidgetResizable(True)
        tabs.addTab(molsc, "FA · LINC", "Cell mechanisms")

        # --- BioGPU-style cell-stress pathway readout (live)
        self.stress = StressView()
        strsc = QScrollArea(); strsc.setWidget(self.stress); strsc.setWidgetResizable(True)
        tabs.addTab(strsc, "Stress pathways", "Cell mechanisms")

        # --- HPA-style compartment atlas (each compartment isolated)
        atsc = QScrollArea(); atsc.setWidget(CompartmentAtlas()); atsc.setWidgetResizable(True)
        tabs.addTab(atsc, "Atlas", "3D & atlas")

        # --- 3D spheroid / multicellular construct cryopreservation (live)
        self.spheroid = SpheroidView()
        tabs.addTab(self.spheroid, "3D spheroid", "3D & atlas")

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
        tabs.addTab(aw, "Analysis", "Data & reference")

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
        tabs.addTab(hw, "HPA reference", "Data & reference")

        # Outcome stays pinned at the top; the grouped dropdown views sit below,
        # in a vertical splitter so the presenter can size them.
        right = QSplitter(Qt.Orientation.Vertical)
        right.addWidget(outcome_w); right.addWidget(tabs)
        right.setStretchFactor(0, 3); right.setStretchFactor(1, 4)
        right.setSizes([440, 560])
        return right

    # ------------------------------------------------------------------ run
    def run(self):
        try:
            self.S, self.R, self.E = simulate(self.P)
        except Exception:
            self.statusBar().showMessage("simulation failed"); traceback.print_exc(); return
        self.view.ice.reseed()          # fresh polycrystal each run (no two freezes alike)
        n = len(self.S)
        self.scrub.blockSignals(True); self.scrub.setRange(0, n - 1)
        self.idx = min(self.idx, n - 1); self.scrub.setValue(self.idx)
        self.scrub.blockSignals(False)
        R = self.R
        self.hero.setText(f"{R['S_24']*100:.0f}%")
        colr = C["good"] if R['S_24'] > .7 else C["warn"] if R['S_24'] > .4 else C["crit"]
        self.hero.setStyleSheet(f"color:{colr}")
        methb = R.get('metHb', 0.0)
        methb_txt = (f" · methaemoglobin {methb*100:.0f}% (O₂-carrying lost)"
                     if self.P.cell_type == "rbc" and methb > 0.01 else "")
        self.hero2.setText(f"viable at {self.P.recover_h:.0f} h post-thaw · membrane-integrity "
                           f"{R['S_imm']*100:.0f}% · settles to {R['S_72']*100:.0f}% by 72 h "
                           f"(delayed-onset death) · functional recovery {R['F_rec']*100:.0f}%{methb_txt}")
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
        self.timelinebar.set_data(self.S.T, self.S.phase, self.S.Vn)
        self._show(self.idx, jumped=True)
        self.statusBar().showMessage(
            f"{n} frames · peak ROCK {R['rockPeak']:.2f} · grain {R['grainMax']:.0f} µm · "
            f"min channel/cell {R['squeezeMin']:.2f} · ATP min {R['atpMin']:.2f}")
        self._update_compare()

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
                    Cout=S.Cout[i], dVw=(S.Vn[i] - S.Vn[i - 1]) if i > 0 else 0.0,
                    dTsc=S.dTsc[i], ros=S.ros[i],
                    sterol=self.P.sterol, cpaLoad=self.P.molar(), cyto=self.P.cyto,
                    cell_type=self.P.cell_type,
                    t=float(self.S.t[i]),
                    hold_min=float(getattr(self.P, "hold_min", 10.0)),
                    adhesion=getattr(self.P, "adhesion", "suspension"),
                    r_iso_um=(3 * self.P.Viso / (4 * math.pi)) ** (1 / 3))

    def _show(self, i, jumped=False):
        f = self._frame_dict(i)
        self.view.set_frame(f, jumped)
        self.mechano.set_frame(f)
        self.molec.set_frame(f)
        self.stress.set_frame(f)
        self.spheroid.set_frame(f)
        names = dict(dry1="Primary drying", dry2="Secondary drying", drystore="Dry storage",
                     rehydrate="Rehydration",
                     load="CPA loading", cool="Cooling", seed="Seeding", store="Storage",
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
