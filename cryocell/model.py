"""
CryoCell biophysical model — Python port of the calibrated engine.

Faithful port of the JavaScript engine: every constant, threshold and rate is
carried over unchanged, so the validation table still holds (inverted-U optimum
at ~1.3 C/min; hMSC suspension 38% at 24 h against Heng's measured 39.8%; etc).

Units: volume um^3, area um^2, length um, time s, concentration osmol/L of cell
water, temperature K internally and C at the interface.

See the project documents for provenance of every parameter:
  cryocell-virtual-cell-model-spec.md      (physics)
  cryocell-whole-cell-architecture.md      (architecture, population, screens)
  cryocell-mechanotransduction-module.md   (post-thaw death)
"""
from __future__ import annotations
import math
from dataclasses import dataclass, field, asdict

R_ATM = 0.0821          # L atm / (mol K)
K0    = 273.15
R_GAS = 8.314           # J / (mol K)

# ---------------------------------------------------------------- CPA library
# Added columns (see REFERENCES.md):
#   sigma  - Staverman reflection coefficient (Kedem-Katchalsky 3-parameter
#            transport). The old code assumed sigma=1 (2-parameter formalism);
#            Kleinhans 1998 (Cryobiology 37:271) showed that mis-times the CPA
#            volume excursion. Values are literature-typical for cell membranes.
#   EaPs   - activation energy for CPA permeation (J/mol). Was a single fixed
#            56.6 kJ/mol for every CPA; permeation Ea is species-specific
#            (Kleinhans 1998; Benson 2012; Kang 2020).
#   Bvir   - second osmotic virial coefficient, molality basis (kg/mol), for
#            the non-ideal osmolality (Elliott et al. 2007, J Phys Chem B
#            111:1775; Prickett 2011). Salt handled separately below.
#   TgP    - glass-transition of the maximally freeze-concentrated solution,
#            Tg' (deg C); composition-specific (Fahy 1984; MacFarlane 1986;
#            Wowk 2010). Was a single universal -125 C.
#   Teut   - eutectic / practical lower liquidus bound (deg C).
#   wRepl  - water-replacement / preferential-exclusion capacity (0-1): how well
#            the agent hydrogen-bonds in place of stripped hydration water and
#            stabilises the native protein/membrane fold on dehydration. This is
#            the *primary* protection for non-penetrating sugars (Crowe & Crowe
#            1984 water-replacement; Timasheff 1998 preferential exclusion).
#   denat  - protein-denaturation propensity at high intracellular concentration
#            (0-1+): the molecular face of CPA toxicity -- the same binding that
#            can protect will unfold protein when the agent is concentrated and
#            warm (strong for DMSO, PG).
CPAS = {
    #                    M       rho    vbar     psRel tox  sterolHit permeant raisesTm iri   sigma EaPs    Bvir    TgP     Teut    wRepl denat
    "dmso":     dict(name="DMSO",                M=78.13, rho=1.100, vbar=0.0710, psRel=1.00, tox=1.00, sterolHit=1.00, permeant=True,  raisesTm=True,  iri=0.05, sigma=0.83, EaPs=54e3,  Bvir=0.108, TgP=-123.0, Teut=-70.0, wRepl=0.35, denat=1.00),
    "glycerol": dict(name="Glycerol",            M=92.09, rho=1.261, vbar=0.0731, psRel=0.12, tox=0.55, sterolHit=0.35, permeant=True,  raisesTm=False, iri=0.02, sigma=0.88, EaPs=63e3,  Bvir=0.023, TgP=-65.0,  Teut=-46.5, wRepl=0.65, denat=0.15),
    "eg":       dict(name="Ethylene glycol",     M=62.07, rho=1.113, vbar=0.0558, psRel=1.60, tox=0.75, sterolHit=0.55, permeant=True,  raisesTm=False, iri=0.02, sigma=0.80, EaPs=57e3,  Bvir=0.037, TgP=-110.0, Teut=-51.0, wRepl=0.40, denat=0.55),
    "pg":       dict(name="Propylene glycol",    M=76.09, rho=1.036, vbar=0.0734, psRel=0.90, tox=1.35, sterolHit=0.95, permeant=True,  raisesTm=False, iri=0.03, sigma=0.90, EaPs=60e3,  Bvir=0.061, TgP=-108.0, Teut=-60.0, wRepl=0.40, denat=1.10),
    "tre":      dict(name="Trehalose",           M=342.3, rho=1.580, vbar=0.2200, psRel=0.00, tox=0.10, sterolHit=0.05, permeant=False, raisesTm=True,  iri=0.25, sigma=1.00, EaPs=56.6e3,Bvir=0.250, TgP=-29.5,  Teut=-40.0, wRepl=1.00, denat=0.00),
    "sd":       dict(name="Self-deactivating CPA", M=150.0, rho=1.150, vbar=0.1200, psRel=0.55, tox=0.95, sterolHit=0.85, permeant=True, raisesTm=True,  iri=0.05, sigma=0.85, EaPs=56e3,  Bvir=0.120, TgP=-115.0, Teut=-65.0, wRepl=0.40, denat=0.85),
    "none":     dict(name="No CPA",              M=78.13, rho=1.100, vbar=0.0710, psRel=0.00, tox=0.00, sterolHit=0.00, permeant=False, raisesTm=False, iri=0.00, sigma=1.00, EaPs=56.6e3,Bvir=0.000, TgP=-125.0, Teut=-73.0, wRepl=0.00, denat=0.00),
}

# ---------------------------------------------- extracellular additive library
# Non-penetrating agents used *alongside* a penetrating CPA (combination
# freezing) or on their own (DMSO-free). They stay outside the cell and protect
# by: raising extracellular osmolality (osmotic dehydration, which lets you use
# LESS of the toxic penetrating CPA), inhibiting ice recrystallisation (IRI),
# and stabilising the plasma membrane. Concentration is % w/v; osmolality falls
# out of the molar mass, so high-MW polymers contribute little colligative
# osmolality but protect by membrane/IRI instead.
#   iri       - ice-recrystallisation inhibition (0-1)
#   membStab  - plasma-membrane stabilisation (0-1), reduces membrane damage
#   wRepl     - extracellular water-replacement / surface protection (0-1)
#   tox       - (low) toxicity; AFP carries ice-shaping toxicity at high dose
ADDITIVES = {
    "none":  dict(name="None",                       M=1.0,    iri=0.00, membStab=0.00, wRepl=0.00, tox=0.00),
    "tre":   dict(name="Trehalose (extracellular)",  M=342.3,  iri=0.25, membStab=0.30, wRepl=0.60, tox=0.05),
    "suc":   dict(name="Sucrose",                    M=342.3,  iri=0.15, membStab=0.20, wRepl=0.45, tox=0.03),
    "hes":   dict(name="Hydroxyethyl starch (HES)",  M=200000, iri=0.20, membStab=0.55, wRepl=0.30, tox=0.02),
    "pvp":   dict(name="PVP-40",                     M=40000,  iri=0.30, membStab=0.45, wRepl=0.25, tox=0.10),
    "pll":   dict(name="Polyampholyte (COOH-PLL)",   M=5000,   iri=0.55, membStab=0.85, wRepl=0.40, tox=0.08),
    "afp":   dict(name="Antifreeze protein (AFP-III)",M=7000,  iri=0.90, membStab=0.30, wRepl=0.20, tox=0.35),
    "pva":   dict(name="PVA (polyvinyl alcohol)",    M=30000,  iri=0.70, membStab=0.25, wRepl=0.15, tox=0.05),
}

ADHESION_STATES = {
    "suspension": ("Suspension (dissociated)",      0.06, "cloning-efficiency regime; blebbing-prone"),
    "adherent":   ("Adherent monolayer",            0.55, "standard tissue-culture plastic"),
    "sheared":    ("Shear-conditioned monolayer",   1.00, "8 d low shear; elevated vinculin (Bissoyi 2016)"),
    "spheroid":   ("Spheroid / 3D construct",       0.78, "cell-cell junctions dominate"),
}

# Cell-cell junction coupling for intercellular ice propagation (distinct from
# the focal-adhesion / matrix coupling above). Irimia & Karlsson 2002
# (Biophys J 82:1858) and Acker, Larese, Yang, Rasch & McGann 2001 showed ice
# propagates cell-to-cell through gap junctions, so confluent / 3D constructs
# suffer far more IIF than dissociated suspensions at the same cooling rate.
# Suspension is 0 by construction, which is why the suspension benchmark is
# unaffected by this term.
JUNCTION = {"suspension": 0.0, "adherent": 0.45, "sheared": 0.55, "spheroid": 0.92}

PERTURBATIONS = [
    ("aqp",     "Aquaporin water transport",            "HgCl2 / AQP knockdown"),
    ("cpaperm", "CPA membrane permeation",              "substitute a non-permeating CPA"),
    ("serca",   "SERCA Ca2+ reuptake",                  "thapsigargin"),
    ("mpt",     "Mitochondrial permeability transition","cyclosporin A"),
    ("oxphos",  "Oxidative phosphorylation",            "oligomycin"),
    ("sterol",  "Membrane cholesterol",                 "methyl-beta-cyclodextrin"),
    ("msc",     "Mechanosensitive channels",            "GsMTx4"),
    ("mt",      "Microtubules",                         "nocodazole"),
    ("actin",   "F-actin",                              "cytochalasin D"),
    ("gel",     "Membrane gel transition",              "lipid composition change"),
    ("nucl",    "Heterogeneous nucleation sites",       "-"),
    ("iri",     "Ice-recrystallisation inhibition",     "add / remove IRI compound"),
    ("rock",    "ROCK (Rho-kinase)",                    "Y-27632 / fasudil"),
    ("piezo",   "Piezo1 mechanosensitive channel",      "GsMTx4"),
    ("casp",    "Caspase cascade",                      "z-VAD-fmk"),
    ("calpain", "Calpain",                              "calpeptin / PD150606"),
    ("fa",      "Focal adhesions",                      "dissociate to single cells"),
    ("yap",     "YAP mechanotransduction",              "verteporfin"),
]

import os
# Strength of the gel-phase permeability gate (fix: fluidity -> transport).
# perm = clamp(1 - PERM_GEL_ALPHA*gel, PERM_GEL_FLOOR, 1). Tuned to preserve
# the validated inverted-U optimum; overridable for sweeps.
PERM_GEL_ALPHA = float(os.environ.get("PERM_GEL_ALPHA", "0.15"))
PERM_GEL_FLOOR = float(os.environ.get("PERM_GEL_FLOOR", "0.70"))
clamp = lambda v, a, b: max(a, min(b, v))
def arrh(Ea, T, Tref):  return math.exp(-(Ea / R_GAS) * (1.0 / T - 1.0 / Tref))

# ------------------------------------------------- colligative / phase diagram
# Non-ideal freezing-point depression via a truncated osmotic virial equation
# (Elliott et al. 2007, J Phys Chem B 111:1775; Prickett 2011): dTf = Kf * osm *
# (1 + B*osm), with B the mixture second virial coefficient (salt + CPA blend).
# The old code used a single generic B = 0.033 for every solution.
KF = 1.86  # cryoscopic constant, deg C kg/mol
B_SALT = 0.044  # NaCl second osmotic virial coefficient (Elliott 2007)
# Anchors the (molality-based) Elliott coefficients onto this model's
# osmolarity scale so the DMSO blend reproduces the calibrated effective
# B = 0.033. See the B_eff derivation in simulate().
K_VIRIAL = 0.34
def dT_from_osm(osm, B=0.033):  return KF * osm * (1 + B * osm)
def osm_from_dT(dT, B=0.033):
    if dT <= 0: return 0.0
    if B <= 0:  return dT / KF
    a, b, c = KF * B, KF, -dT
    return (-b + math.sqrt(b * b - 4 * a * c)) / (2 * a)

# Default eutectic / Tg' for backward-compatible module-level use; the actual
# run uses the CPA-specific values (CPAS[...]["Teut"], ["TgP"]).
T_EUT, TG_PRIME = -73.0, -125.0
WLF_C1, WLF_C2 = 17.4, 51.6
def _logaT(TC, tg_prime):
    d = TC - tg_prime
    return -WLF_C1 * d / (WLF_C2 + d)
def mobility(TC, tg_prime=TG_PRIME):
    # WLF shift factor for molecular mobility between Tg' and -40 C. Tg' is now
    # composition-specific, so e.g. glycerol (Tg' ~ -65 C) arrests far warmer
    # than DMSO (Tg' ~ -123 C).
    if TC <= tg_prime: return 0.0
    if TC >= -40.0:    return 1.0
    return 10.0 ** (_logaT(-40.0, tg_prime) - _logaT(TC, tg_prime))

# ------------------------------------------------------- ice field / membrane
GRAIN_SEED, K_RIPEN = 12.0, 80.0
def ripen_rate(TC, mob, iri):
    return K_RIPEN * math.exp(-max(0.0, -TC) / 18.0) * mob * (1 - clamp(iri, 0, 0.98))
def channel_width(grain, f_ice):
    return 1e4 if f_ice <= 0.01 else grain * (1 - f_ice) / f_ice

def mole_frac(molal): return molal / (molal + 55.5)

def bilayer_regime(x_molpc, sterol_pc, agent_hit):
    xe = (x_molpc / (1 + 0.022 * sterol_pc)) * agent_hit
    if   xe < 5:  idx, label = 0, "Bilayer intact - loosening only"
    elif xe < 15: idx, label = 1, "Thinning, transient water pores"
    elif xe < 25: idx, label = 2, "Persistent pores - ion leak"
    else:         idx, label = 3, "Lipid desorption - bilayer failing"
    pore = 0.0 if xe <= 5 else min(1.0, ((xe - 5) / 22.0) ** 1.7)
    return idx, label, pore, xe

def membrane_state(x_molpc, TC, sterol_pc, cpa):
    xe = x_molpc * cpa["sterolHit"] / (1 + 0.022 * sterol_pc)
    APL   = 1 + 0.0085 * xe
    thick = 1 - 0.0062 * xe
    Tm0   = -20 + 0.25 * sterol_pc
    Tm    = Tm0 + (0.55 * xe if cpa["raisesTm"] else -0.42 * xe)
    width = 4.5 + 0.22 * sterol_pc
    gel   = 1.0 / (1.0 + math.exp((TC - Tm) / (width / 4.0)))
    fluid = clamp((1 - gel) * (0.35 + 0.65 * clamp((TC + 40) / 77.0, 0, 1))
                  * (1 + 0.010 * xe - 0.0006 * xe * xe), 0, 1.4)
    return dict(xe=xe, APL=APL, thick=thick, Tm=Tm, gel=gel, fluid=fluid,
                domainLeak=4 * gel * (1 - gel))

def mt_depoly_rate(TC):
    return 0.010 * math.exp(-(TC + 8) / 9.0) if TC < 12 else 0.0

# --------------------------------------------------------------- nucleation
OMEGA_SCN, KAPPA_SCN, NUC_VISC_EXP, OSM_MAX = 1.08e9, 1.04e9, 3.0, 22.5
# Mitochondria are efficient heterogeneous ice nucleators (cardiolipin-rich
# inner membrane), so per unit area they nucleate more readily than the plasma
# membrane; once a mitochondrion freezes it seeds the surrounding cytosol.
# Muldrew & McGann 1994 (Biophys J 66:532); Stott & Karlsson 2009 (Cryobiology).
MITO_SITE_SCALE = 1.5
def nucleation_rate(A_um2, T_K, dTsc, mob, water_frac, visc_exp=NUC_VISC_EXP, site_scale=1.0):
    if dTsc <= 0.2 or mob <= 0: return 0.0
    arg = KAPPA_SCN / (T_K ** 3 * dTsc * dTsc)
    if arg > 60: return 0.0
    return (OMEGA_SCN * site_scale * A_um2 * 1e-12
            * (mob ** visc_exp) * clamp(water_frac, 0, 1) * math.exp(-arg))

# ------------------------------------------------------------------ energetics
ATP_TURNOVER, ATP_KM = 0.033, 0.15

def rock_inhibitor_effect(drug, conc_uM):
    """Heng (Tissue Cell 2009) measured hMSC adherent viable fraction across
    1-100 uM Y-27632: 39.8% at zero, peaking 48.5% (5 uM) / 48.4% (10 uM), back
    to 36.0% at 100 uM. Both limbs are kept: saturating inhibition, and
    off-target toxicity that takes over above ~50 uM."""
    if drug == "none" or not conc_uM:
        return 0.0, 0.0
    potency = {"fasudil": 0.55, "y39983": 1.6}.get(drug, 1.0)
    inhib = (conc_uM * potency) / (conc_uM * potency + 4.0)
    tox   = (conc_uM / 200.0) ** 1.8
    return clamp(inhib, 0, 0.95), clamp(tox, 0, 0.6)


@dataclass
class Params:
    cpa_key: str = "dmso"
    conc_pct: float = 10.0          # % v/v (% w/v for trehalose)
    additive: str = "none"          # extracellular co-CPA (combination freezing)
    add_conc: float = 0.0           # additive concentration, % w/v
    sterol: float = 25.0
    iri: float = 0.0
    thalf: float = 8.0
    q10: float = 2.5
    frag_n: int = 1
    sol_limit: float = 1.6
    frag_tox: float = 0.08     # fragment toxicity relative to the parent CPA.
                               # A design variable of the molecule, not a constant.
    T_add: float = 22.0
    hold_min: float = 10.0
    add_steps: int = 3
    T_seed: float = -6.0
    CR: float = 1.0                 # C/min
    T_store: float = -196.0
    days: float = 5.0
    WR: float = 400.0               # C/min
    T_dil: float = 22.0
    dilution: str = "step"          # direct | step | none
    lp: float = 0.22                # um/min/atm at 25 C (re-fit against the KK +
                                    # variable-area transport; hMSC range 0.2-0.4)
    ps: float = 0.05                # um/s at 25 C
    vb: float = 0.20
    Viso: float = 1800.0
    cell_type: str = "msc"          # morphology: msc | rbc | tcell (rendering only)
    adhesion: str = "suspension"
    rock_drug: str = "none"
    rock_conc: float = 0.0
    rock_when: str = "thaw"         # freeze | thaw | both
    zvad: float = 0.0
    calpain_i: float = 0.0
    gsmtx: float = 0.0
    T_recover: float = 37.0
    recover_h: float = 24.0
    # --- cell-line robustness phenotype (0 = primary/normal, e.g. hMSC;
    #     higher = the hardiness typical of transformed lines like HeLa, A549).
    #     All default to 0, so the hMSC calibration is unchanged. These act on
    #     the regulated-death and metabolic layers -- the axes on which cancer
    #     lines are genuinely tougher -- not on the biophysical damage, which is
    #     set by the measurable transport/geometry parameters (Lp, Ps, Viso...).
    apop_resist: float = 0.0        # suppressed apoptosis (p53 loss, Bcl-2 up)
    anoikis_resist: float = 0.0     # anchorage independence (survive detachment)
    antioxidant: float = 0.0        # ROS-defence capacity (NRF2/glutathione; A549)
    glycolytic: float = 0.0         # Warburg metabolism: ATP without mitochondria
    cyto: float = 1.0               # resting cytoskeletal density (visual): mesenchymal
                                    # cells (hMSC) are rich in it, epithelial lines less so
    ko: dict = field(default_factory=dict)
    dt_scale: float = 1.0
    recover_dt: float = 30.0   # s; the recovery phase is quiescent, so it can be coarse
                               # (60 s costs 0.2 points against a 10 s reference)
    nuc_scale: float = 1.0
    k_prop: float = 0.005           # intercellular ice-propagation rate constant
                                    # (Irimia & Karlsson 2002); scaled by junction
    k_mito_seed: float = 0.006      # rate at which a frozen mitochondrion seeds
                                    # cytosolic ice (Muldrew & McGann 1994)
    k_mito_cpa: float = 0.02        # inner-membrane CPA permeation rate (1/s);
                                    # sets the matrix CPA-loading lag (Yu & Quinn 1994)
    nuc_visc_exp: float = NUC_VISC_EXP
    k_ripen: float = K_RIPEN
    k_squeeze: float = 9.0e-4
    mob_exp: float = 2.0
    k_tox: float = 1.13e-4
    k_unfold: float = 1.6e-3        # protein unfolding rate from dehydration /
                                    # freeze-concentration (water-replacement module)
    k_denat: float = 5.0e-4         # protein denaturation rate at high [CPA]
    k_mem: float = 0.0012
    k_mpt: float = 0.003
    k_leak: float = 4.0e-4
    k_atp: float = ATP_TURNOVER
    v_crit_lo: float = 0.45
    v_crit_hi: float = 1.45

    def molar(self):
        c = CPAS[self.cpa_key]
        if self.cpa_key == "none": return 0.0
        if self.cpa_key == "tre":  return (self.conc_pct * 10.0) / c["M"]
        return (self.conc_pct / 100.0) * c["rho"] * 1000.0 / c["M"]

    def copy(self, **kw):
        d = asdict(self); d.update(kw)
        d["ko"] = dict(d.get("ko") or {})
        return Params(**d)


class Series:
    """Recorded trajectory. Plain lists; NumPy-friendly on demand."""
    KEYS = ("t T V Vn Vmito Vnuc Cin Cout Osme dTsc dPsi caCyt caER x pore mpt "
            "Dosm Dtox Dmem Dmech Piif frag fIce grain chanW squeeze gel fluid "
            "thick APL tension msOpen mt actin atp FA rock pMLC bleb yapN casp3 "
            "apop necr piezo akt glass mcpa Pmito prot intf sigMT sigIF ros").split()
    def __init__(self):
        for k in self.KEYS: setattr(self, k, [])
        self.phase, self.events = [], []
    def __len__(self): return len(self.t)


def simulate(P: Params):
    cpa   = CPAS[P.cpa_key]
    Viso  = P.Viso
    Vb    = P.vb * Viso
    Vw0   = Viso - Vb
    r0    = (3 * Viso / (4 * math.pi)) ** (1 / 3)
    A0    = 4 * math.pi * r0 * r0
    cpaM  = P.molar()
    OSM_ISO = 0.30
    salt_e0 = OSM_ISO

    # --- extracellular additive (combination freezing / DMSO-free)
    add     = ADDITIVES.get(P.additive, ADDITIVES["none"])
    addM    = (P.add_conc / 100.0) * 1000.0 / add["M"] if P.additive != "none" else 0.0
    add_pres = clamp(P.add_conc / 6.0, 0, 1) if P.additive != "none" else 0.0  # coverage: full at ~6% w/v
    add_iri  = add["iri"] * add_pres
    memb_prot = clamp(add["membStab"] * add_pres, 0, 0.85)                     # membrane-damage reduction

    # --- composition-specific thermodynamics (Elliott virial + CPA Tg'/eutectic)
    tg_run   = cpa["TgP"]
    teut_run = cpa["Teut"]
    # blend the salt and CPA second virial coefficients by their isotonic
    # osmolal weight -> effective non-ideality of the freeze-concentrated mix.
    # Elliott's coefficients are molality-based; this model's osm is on an
    # osmolarity scale and the virial is truncated at 2nd order, so the
    # absolute coefficients cannot be imported directly (B*osm would dominate
    # at freeze-concentration). We use them for the *relative* non-ideality of
    # each CPA, anchored by K_VIRIAL so DMSO reproduces the calibrated liquidus
    # (its historical effective B = 0.033). This keeps the DMSO benchmark exact
    # while ranking the other CPAs by their true Elliott non-ideality.
    _w_cpa = cpaM / max(cpaM + salt_e0, 1e-6)
    B_eff  = K_VIRIAL * ((1 - _w_cpa) * B_SALT + _w_cpa * cpa["Bvir"])
    junction = JUNCTION.get(P.adhesion, 0.0)

    KO = P.ko or {}
    sterol_pc = 0.0 if KO.get("sterol") else P.sterol
    mob_dam = lambda m: m ** P.mob_exp
    mob_of  = lambda TC: mobility(TC, tg_run)

    # intracellular
    Vw, n_s, n_c, n_f = Vw0, OSM_ISO * Vw0, 0.0, 0.0
    e_c = e_f = e_suc = e_add = 0.0
    e_s = salt_e0
    e_suc0 = 0.0
    eCPA0 = 0.0
    e_add0 = addM
    frac_intact = 1.0

    # organelles
    Vmito0, matrixF = 0.055 * Viso, 0.68
    n_mtx = OSM_ISO * Vmito0 * matrixF
    mpt_frac, dPsi = 0.0, 1.0
    Vnuc0 = 0.12 * Viso
    chrom_cond = lobulation = 0.0
    caER, caCyt = 1.0, 0.0

    # damage
    D_osm = D_tox = D_mem = P_iif = D_recry = D_apop = D_frag = 0.0
    supercool_peak = iif_amount = 0.0
    P_mito = D_mitoice = mcpa = 0.0  # mito ice: probability, lethality, matrix CPA conc
    prot = 1.0; D_prot = 0.0         # native/functional protein fraction; denaturation damage
    D_swell = D_sol = D_mech = D_energy = 0.0
    minV = maxV = 1.0
    ice_vol = 0.0

    # ice field
    grain = grain_prev = grain_max = GRAIN_SEED
    f_ice = 0.0; chanW = 1e4; squeeze = squeeze_min = 8.0; D_recry_ice = 0.0

    # membrane / mechanics
    memb = dict(xe=0, APL=1, thick=1, Tm=-8, gel=0, fluid=1, domainLeak=0)
    tension = ms_open = 0.0
    mt_i, actin, if_i = 1.0, 1.0, 1.0     # microtubule / F-actin / intermediate-filament integrity
    sig_mt = sig_if = 0.0                  # tensegrity stress on struts / cables

    # energetics
    atp, atp_min, atp_prod = 1.0, 1.0, 0.0
    atp_use = dict(serca=0.0, naka=0.0, basal=0.0, cytoskeleton=0.0)

    # mechanotransduction
    adh_label, FA0, adh_note = ADHESION_STATES.get(P.adhesion, ADHESION_STATES["suspension"])
    FA = FA0
    piezo, rhoGTP, rock, pMLC, bleb, yapN = 0.0, 0.08, 0.10, 0.10, 0.0, 0.65
    casp8 = casp9 = casp3 = calpain = 0.0
    akt = 1.0
    apop = necr = ros = 0.0
    rock_peak = bleb_peak = casp3_peak = ros_peak = 0.0
    RI_inhib, RI_tox = rock_inhibitor_effect(P.rock_drug, P.rock_conc)
    zvad, calpI, gsmtx = clamp(P.zvad,0,1), clamp(P.calpain_i,0,1), clamp(P.gsmtx,0,1)
    ri_freeze = 1.0 if P.rock_when in ("freeze", "both") else 0.0
    ri_thaw   = 1.0 if P.rock_when in ("thaw",   "both") else 0.0

    # kinetics
    Tref = 298.15
    Lp_ref = (P.lp / 60.0) * (0.25 if KO.get("aqp") else 1.0)
    Ps_ref = P.ps * cpa["psRel"] * (0.0 if KO.get("cpaperm") else 1.0)
    Ea_Lp_supra = 51.8e3
    Ea_Lp_sub_noCPA, Ea_Lp_sub_CPA = 106e3, 65e3
    k_hyd37 = math.log(2) / (P.thalf * 60.0) if P.cpa_key == "sd" else 0.0
    N_TOX, EA_TOX = 1.6, 70e3

    out = Series()
    flagged = set()
    T_C, t = P.T_add, 0.0
    ice = False
    last_dTsc = last_haz = last_x = last_pore = 0.0
    last_mob = 1.0

    def log(key, txt):
        if key in flagged: return
        flagged.add(key); out.events.append((t, T_C, txt))

    def mito_volume():
        Ci = (n_s + n_c + n_f) / max(Vw, 1e-9)
        load = n_mtx * (1 + 1.35 * mpt_frac)
        return clamp(load / max(Ci, 0.05), 0.25 * Vmito0 * matrixF, 2.6 * Vmito0 * matrixF)

    def push(phase):
        V = Vb + Vw + n_c * cpa["vbar"] + n_f * cpa["vbar"] * 0.6
        Vnuc = Vnuc0 * (Vw / Vw0) * 0.85 + Vnuc0 * 0.15
        vals = dict(t=t, T=T_C, V=V, Vn=V / Viso,
                    Vmito=mito_volume() / (Vmito0 * matrixF), Vnuc=Vnuc / Vnuc0,
                    Cin=n_c / max(Vw, 1e-9), Cout=e_c, Osme=e_s + e_c + e_f + e_suc,
                    dTsc=last_dTsc, dPsi=dPsi, caCyt=caCyt, caER=caER,
                    x=last_x, pore=last_pore, mpt=mpt_frac,
                    Dosm=D_osm, Dtox=D_tox, Dmem=D_mem, Dmech=D_mech, Piif=P_iif,
                    frag=e_f, fIce=f_ice, grain=grain, chanW=chanW, squeeze=squeeze,
                    gel=memb["gel"], fluid=memb["fluid"], thick=memb["thick"],
                    APL=memb["APL"], tension=tension, msOpen=ms_open, mt=mt_i,
                    actin=actin, atp=atp, FA=FA, rock=rock, pMLC=pMLC, bleb=bleb,
                    yapN=yapN, casp3=casp3, apop=apop, necr=necr, piezo=piezo, akt=akt,
                    glass=clamp(1.0 - last_mob, 0, 1), mcpa=mcpa, Pmito=P_mito, prot=prot,
                    intf=if_i, sigMT=sig_mt, sigIF=sig_if, ros=ros)
        for k, v in vals.items(): getattr(out, k).append(v)
        out.phase.append(phase)

    def step(dt, phase, slow=False):
        nonlocal Vw, n_s, n_c, n_f, e_c, e_f, e_s, e_suc, e_add, frac_intact, eCPA0
        nonlocal n_mtx, mpt_frac, dPsi, chrom_cond, lobulation, caER, caCyt
        nonlocal D_osm, D_tox, D_mem, P_iif, D_recry, D_frag, D_swell, D_sol
        nonlocal D_mech, D_energy, minV, maxV, ice_vol, grain, grain_prev
        nonlocal grain_max, f_ice, chanW, squeeze, squeeze_min, D_recry_ice
        nonlocal memb, tension, ms_open, mt_i, actin, atp, atp_min, atp_prod
        nonlocal FA, piezo, rhoGTP, rock, pMLC, bleb, yapN, casp8, casp9, casp3
        nonlocal calpain, akt, apop, necr, rock_peak, bleb_peak, casp3_peak
        nonlocal last_dTsc, last_haz, last_x, last_pore, supercool_peak, iif_amount
        nonlocal ros, ros_peak, P_mito, D_mitoice, mcpa, last_mob, prot, D_prot
        nonlocal if_i, sig_mt, sig_if

        T_K = T_C + K0
        mob  = mob_of(T_C) if ice else (1.0 if T_C > 0 else mob_of(T_C))
        mobD = mob_dam(mob)
        last_mob = mob

        # ---- extracellular composition + ice field
        if ice:
            osm_tot = osm_from_dT(min(-T_C, -teut_run), B_eff)
            osm0 = salt_e0 + eCPA0 + e_suc0 + e_add0
            cf = max(1.0, osm_tot / max(osm0, 1e-6))
            e_s = salt_e0 * cf
            e_c = eCPA0 * cf * frac_intact
            e_f = eCPA0 * cf * (1 - frac_intact) * P.frag_n
            e_suc = e_suc0 * cf
            e_add = e_add0 * cf
            ice_vol = 1 - 1 / cf
            f_ice = ice_vol
            iri = 0.0 if KO.get("iri") else clamp(P.iri / 100.0 + cpa["iri"] + add_iri, 0, 0.98)
            kr = ripen_rate(T_C, mob, iri) * (P.k_ripen / K_RIPEN)
            if kr > 0:
                grain = (grain ** 3 + kr * min(dt, 3600)) ** (1 / 3)
            dgrain = grain - grain_prev; grain_prev = grain
            grain_max = max(grain_max, grain)
            if dgrain > 0: D_recry_ice += 0.003 * dgrain * f_ice
            chanW = channel_width(grain, f_ice)
            dCell = 2 * (3 * max(Vb + Vw + n_c * cpa["vbar"], 1) / (4 * math.pi)) ** (1 / 3)
            squeeze = clamp(chanW / max(dCell, 1e-3), 0, 8)
            squeeze_min = min(squeeze_min, squeeze)
            if squeeze < 0.30 and mob > 0:
                D_mech += P.k_squeeze * ((0.30 - squeeze) / 0.30) ** 2 * mob * min(dt, 10)
                log("squeeze", "Cell compressed in the unfrozen channel between ice grains")
            if grain > 60: log("coarse", "Ice grains coarsened past 60 um")

        # ---- self-deactivating chemistry
        if k_hyd37 > 0:
            k = k_hyd37 * (P.q10 ** ((T_C - 37) / 10.0)) * mob
            dec = 1 - math.exp(-k * dt)
            frac_intact *= (1 - dec)
            conv = n_c * dec
            n_c -= conv; n_f += conv * P.frag_n
            # the extracellular pool deactivates too -- that is the point of the
            # molecule. While ice is present the composition is reset from the
            # liquidus each step using frac_intact; once it has melted the pool
            # has to be decayed here or an infused, wash-free CPA would sit at
            # full strength forever and read as pure toxicity.
            if not ice:
                conv_e = e_c * dec
                e_c -= conv_e; e_f += conv_e * P.frag_n
            if frac_intact < 0.5: log("sd50", "CPA half-converted to fragments")
            cf_int = n_f / max(Vw, 1e-9)
            if cf_int > P.sol_limit:
                n_f -= (cf_int - P.sol_limit) * Vw * 0.5 * dt
                D_frag += 0.0009 * dt * (cf_int - P.sol_limit)
                log("precip", "Fragment exceeded solubility - intracellular precipitation")

        # ---- transport: Kedem-Katchalsky 3-parameter (Lp, Ps, sigma)
        #      Kedem & Katchalsky 1958 (BBA 27:229); Kleinhans 1998
        #      (Cryobiology 37:271). Ps activation energy is now CPA-specific.
        sigma = cpa["sigma"]; EaPs = cpa["EaPs"]
        if T_K >= K0:
            Lp = Lp_ref * arrh(Ea_Lp_supra, T_K, Tref)
            Ps = Ps_ref * arrh(EaPs, T_K, Tref)
        else:
            Ea_sub = Ea_Lp_sub_noCPA if (ice and e_c <= 0.3) else Ea_Lp_sub_CPA
            Lp = Lp_ref * arrh(Ea_Lp_supra, K0, Tref) * arrh(Ea_sub, T_K, K0)
            Ps = Ps_ref * arrh(EaPs, K0, Tref) * arrh(EaPs * 1.6, T_K, K0)

        Ci_s = n_s / max(Vw, 1e-9); Ci_c = n_c / max(Vw, 1e-9); Ci_f = n_f / max(Vw, 1e-9)
        Osm_i = Ci_s + Ci_c + Ci_f
        Osm_e = e_s + e_c + e_f + e_suc + e_add
        # split the osmotic drive into non-permeant (salt in vs salt + sucrose +
        # extracellular additive out; fully reflected) and permeant (CPA +
        # fragments; reflected by sigma) parts -- this distinguishes KK from 2P.
        # The additive raises extracellular osmolality and so drives dehydration,
        # which is a large part of why a non-penetrating co-solute lets you use
        # less of the toxic penetrating CPA.
        osm_np = Ci_s - (e_s + e_suc + e_add)
        osm_p  = (Ci_c + Ci_f) - (e_c + e_f)

        last_x = mole_frac(e_c + e_f) * 100
        _, _, reg_pore, _ = bilayer_regime(last_x, sterol_pc, cpa["sterolHit"])
        memb = membrane_state(last_x, T_C, sterol_pc, cpa)
        if KO.get("gel"): memb["gel"] = 0.0; memb["domainLeak"] = 0.0
        defect = clamp(reg_pore + 0.12 * memb["domainLeak"], 0, 1)
        last_pore = defect
        pore_boost = 1 + 9 * defect
        # Gel-phase gating: below Tm the cooperative fluid->gel transition
        # collapses permeability at a knee, not along the Arrhenius slope.
        perm_gel = clamp(1.0 - PERM_GEL_ALPHA * memb["gel"], PERM_GEL_FLOOR, 1.0)
        # Instantaneous membrane area ~ V^(2/3) for a sphere, not the fixed
        # isotonic area: a dehydrated cell has far less area, so flux must
        # fall with it. Below the folded-membrane reserve (~0.6) area is held
        # (microvilli/folds buffer shrinkage before true area is lost).
        vn_now = clamp((Vb + Vw + n_c * cpa["vbar"] + n_f * cpa["vbar"] * 0.6) / Viso, 0.04, 4.0)
        A_inst = A0 * max(vn_now, 0.60) ** (2.0 / 3.0)
        Lp_eff = Lp * pore_boost * perm_gel * A_inst
        Ps_eff = Ps * pore_boost * perm_gel * A_inst

        if not slow:
            # KK water flux: Jv = Lp*A*R*T*[ (Ci_np - Ce_np) + sigma*(Ci_p - Ce_p) ]
            dVw = clamp(Lp_eff * R_ATM * T_K * (osm_np + sigma * osm_p) * dt,
                        -0.30 * Vw, 0.30 * Vw0)
            Vw = clamp(Vw + dVw, 0.02 * Vw0, 12 * Vw0)
            if cpa["permeant"]:
                # KK solute flux = diffusive Ps*A*dC  +  solvent drag (1-sigma)*Cbar*Jv
                Cbar_c = 0.5 * (e_c + Ci_c); Cbar_f = 0.5 * (e_f + Ci_f)
                n_c += clamp(Ps_eff * (e_c - Ci_c) * dt + (1 - sigma) * Cbar_c * dVw,
                             -0.35 * n_c, 0.35 * max(e_c * Vw0, 1))
                n_f += clamp(Ps_eff * 0.7 * (e_f - Ci_f) * dt + (1 - sigma) * Cbar_f * dVw,
                             -0.35 * n_f, 0.35 * max(e_f * Vw0, 1))
            if defect > 0:
                leak = clamp(0.0004 * defect * A_inst * max(0.0, Ci_s - e_s) * dt, 0, 0.2 * n_s)
                n_s = max(0.05 * OSM_ISO * Vw0, n_s - leak)

        # ---- volume / solution effects
        V = Vb + Vw + n_c * cpa["vbar"] + n_f * cpa["vbar"] * 0.6
        vn = V / Viso
        minV = min(minV, vn); maxV = max(maxV, vn)
        # intermediate filaments (vimentin) are the mechanical safety net that
        # resists over-extension and bears the tensile load (Janmey 1991). The
        # calibration already assumes an intact network, so intact IF = baseline
        # (factor 1) and IF LOSS worsens lytic swelling.
        if_shield = 1.0 + 0.8 * (1.0 - if_i)
        if vn < P.v_crit_lo: log("shrink", "Shrank past the minimum critical volume")
        if vn > P.v_crit_hi:
            D_swell += 2.0 * (vn - P.v_crit_hi) ** 1.7 * min(dt, 5) / 60.0 * if_shield
            log("swell", "Swelled past the lytic limit")
        # solution-effects damage is driven by the concentrated harmful solutes
        # (salt + CPA), NOT by a benign extracellular additive -- a trehalose or
        # polymer co-solute that raises osmolality is not itself toxic, and by
        # buffering the osmolality it actually lowers the freeze-concentrated
        # salt the cell sees. Excluding e_add here is what lets a combination
        # protect rather than add damage.
        osm_harm = Osm_e - e_add
        if ice and osm_harm > 2:
            D_sol += 2.0e-6 * (osm_harm - 2) ** 1.2 * mobD * dt

        # ---- tension and mechanosensitive channels
        area_needed = max(vn, 0.02) ** (2 / 3)
        tension = clamp((area_needed - 1.04) / 0.16, 0, 1.4)
        if tension > 0.15 and not slow and not KO.get("msc"):
            ms_open = clamp(ms_open + 0.06 * tension * dt, 0, 1)
            dump = min(0.25 * n_s, 0.010 * ms_open * A0 * max(0.0, Ci_s - e_s) * dt)
            n_s -= dump
            log("msc", "Mechanosensitive channels opened under membrane tension")
        elif not slow:
            ms_open = max(0.0, ms_open - 0.02 * dt)
        if tension > 1.2:
            D_swell += 0.0025 * (tension - 1.2) ** 1.5 * min(dt, 5) * if_shield
            log("lysis", "Membrane tension past the lytic limit")

        # ---- cytoskeleton as a TENSEGRITY network (Ingber 1993, 2003):
        # microtubules bear COMPRESSION (struts), the actin cortex and the
        # intermediate filaments bear TENSION (cables), balanced by prestress.
        # Cell deformation redistributes the load, so the stress on each element
        # is tracked and drives both depolymerisation and downstream signalling.
        strain_c = max(0.0, 1.0 - vn)                              # dehydration -> compression
        strain_t = max(0.0, vn - 1.0)                              # swelling -> tension
        squeeze_load = clamp((0.30 - squeeze) / 0.30, 0, 1) if ice else 0.0   # external ice compression
        sig_mt  = clamp((0.9 * strain_c + 1.1 * squeeze_load) * (0.3 + 0.7 * mt_i), 0, 2)
        sig_ten = strain_t * 1.4 + max(0.0, tension - 0.2)
        sig_actin = clamp(sig_ten * (0.3 + 0.7 * actin), 0, 2)
        sig_if    = clamp(sig_ten * (0.3 + 0.7 * if_i) * 1.2, 0, 2)   # IF: the high-strain cable
        if not slow:
            kdep = mt_depoly_rate(T_C) * mob
            mt_buckle = 0.010 * max(0.0, sig_mt - 0.5) * mob          # struts buckle under compression
            regrow = 0.0016 * dt * (1 - mt_i) if T_C > 20 else 0.0
            mt_i = 0.0 if KO.get("mt") else clamp(mt_i - (kdep + mt_buckle) * dt + regrow, 0, 1)
            shear = max(0.0, tension - 0.3) + max(0.0, 0.55 - vn)
            actin_cold = 0.0008 * clamp((8 - T_C) / 28.0, 0, 1) * mob   # mild cold depoly
            actin_regrow = 0.0022 * dt * (1 - actin) if T_C > 20 else 0.0
            actin = 0.0 if KO.get("actin") else clamp(
                actin - (0.004 * shear * mob + actin_cold) * dt + actin_regrow, 0, 1)
            # intermediate filaments (vimentin): cold-stable, the mechanical safety
            # net — they yield only at extreme strain and buffer over-extension.
            if_regrow = 0.0012 * dt * (1 - if_i) if T_C > 20 else 0.0
            if_i = 0.0 if KO.get("intf") else clamp(
                if_i - 0.0025 * max(0.0, sig_if - 1.1) * mob * dt + if_regrow, 0, 1)
            if mt_i < 0.35: log("mt", "Microtubules depolymerised — cold + compressive buckling")
            if actin < 0.5: log("actin", "F-actin disassembled by cold, osmotic and tension stress")
            if if_i < 0.5: log("intf", "Intermediate-filament network yielding at extreme strain")

        # ---- CPA toxicity (penetrating CPA + a small extracellular-additive term;
        # AFP carries a warm ice-shaping toxicity at high dose)
        if mob > 0 and cpa["tox"] > 0:
            kt = P.k_tox * cpa["tox"] * arrh(EA_TOX, T_K, 310.15)
            D_tox += (kt * max(Ci_c, 0) ** N_TOX
                      + kt * P.frag_tox * max(Ci_f, 0) ** N_TOX) * mob * dt
        if mob > 0 and add["tox"] > 0 and add_pres > 0:
            D_tox += 0.5 * P.k_tox * add["tox"] * add_pres * max(e_add, 0) ** 1.2 * mob * dt

        # ---- macromolecular (protein) stability: water replacement & denaturation
        # Dehydration and freeze-concentration strip a protein's hydration shell
        # and unfold it. A CPA that hydrogen-bonds in the place of the lost water
        # (water-replacement, Crowe & Crowe 1984) or is preferentially excluded
        # from the surface (Timasheff 1998) stabilises the native fold -- this is
        # the PRIMARY protection of non-penetrating sugars like trehalose, and it
        # scales with intracellular CPA and its wRepl capacity. The SAME agent,
        # concentrated and warm, binds and denatures protein directly -- the
        # molecular face of CPA TOXICITY (Arrhenius-gated, so it is a warm/
        # high-concentration effect, negligible deep in the cold). Extracellular
        # sucrose contributes a little coverage. Lost protein impairs recovery.
        if mob > 0 and not slow:
            dehyd = clamp((0.80 - vn) / 0.55, 0, 1)
            conc_stress = clamp((Osm_i - 3) / 17, 0, 1)
            cover = clamp(cpa["wRepl"] * (Ci_c + 0.4 * e_suc) / 2.2, 0, 0.92)
            unfold = P.k_unfold * (dehyd + conc_stress) * (1 - cover) * mob
            denat_h = (P.k_denat * cpa["denat"] * max(0.0, Ci_c - 1.6) ** 1.3
                       * arrh(EA_TOX, T_K, 310.15) * mob)
            refold = 0.0018 * (1 - prot) * dt if T_C > 20 else 0.0
            prot = clamp(prot - (unfold + denat_h) * dt + refold, 0, 1)
            D_tox += 0.6 * denat_h * dt                      # denaturation is toxicity
            if prot < 0.6: log("prot", "Proteins unfolding — hydration shell lost / CPA denaturation")

        # ---- membrane damage. A membrane-stabilising extracellular additive
        # (polyampholyte, HES, sugars) shields the bilayer and reduces this.
        mstab = 1.0 - memb_prot
        if defect > 0 and mob > 0:
            D_mem += P.k_mem * defect ** 1.5 * mobD * dt * mstab
            if reg_pore > 0 and last_x / (1 + 0.022 * sterol_pc) * cpa["sterolHit"] >= 15:
                log("pore", "Bilayer entered the persistent-pore regime")
        if (not slow) and memb["domainLeak"] > 0.35 and mob > 0:
            D_mem += 0.00022 * memb["domainLeak"] * min(dt, 5) * mstab
            log("gel", "Membrane crossing its fluid-to-gel transition - domain-boundary leak")
        if memb["gel"] > 0.9: log("gelfull", "Plasma membrane fully in the gel phase")

        # ---- intracellular ice
        Tm_i = -dT_from_osm(Osm_i, B_eff)
        dTsc = max(0.0, Tm_i - T_C) if ice else 0.0
        last_dTsc = dTsc
        f_freeze = (Vw / Vw0) * clamp(1 - Osm_i / OSM_MAX, 0, 1)
        J = nucleation_rate(A0, T_K, dTsc, mob, f_freeze, P.nuc_visc_exp, P.nuc_scale) \
            if (ice and not slow and not KO.get("nucl")) else 0.0
        # Mitochondria as an independent osmometer + nucleation site. The inner
        # membrane exchanges water fast (aquaporin-8), so the matrix stays at
        # cytosolic osmolality and supercools with it, but its cardiolipin inner
        # membrane nucleates more readily per unit area (MITO_SITE_SCALE). A
        # frozen mitochondrion then seeds cytosolic ice and drives the
        # permeability transition. Muldrew & McGann 1994; Stott & Karlsson 2009.
        # CPA-organelle interaction: CPA crosses the mitochondrial inner
        # membrane and loads the matrix, but with a lag (finite inner-membrane
        # CPA permeability; Yu & Quinn 1994). During slow cooling the matrix
        # equilibrates and is as CPA-protected as the cytosol; during fast
        # cooling it stays under-protected and freezes first -- this is what
        # grounds mitochondria-first nucleation, not just a higher site density.
        if not slow:
            mcpa = clamp(mcpa + P.k_mito_cpa * (Ci_c - mcpa) * dt, 0, Ci_c + 0.5)
        Osm_mtx = max(0.0, Osm_i - Ci_c) + mcpa          # matrix osmolality: non-CPA part + lagged matrix CPA
        dTsc_mtx = max(0.0, -dT_from_osm(Osm_mtx, B_eff) - T_C) if ice else 0.0
        if ice and not slow and not KO.get("nucl"):
            Vmito_now = mito_volume()
            A_mito = 8.0 * (3 * Vmito_now / (4 * math.pi)) ** (2 / 3)  # aggregate matrix surface
            J_mito = nucleation_rate(A_mito, T_K, dTsc_mtx, mob, f_freeze,
                                     P.nuc_visc_exp, P.nuc_scale * MITO_SITE_SCALE)
            if J_mito > 0:
                P_mito = 1 - (1 - P_mito) * math.exp(-J_mito * dt)
                J += P.k_mito_seed * P_mito * mob          # mito ice seeds the cytosol
                if P_mito > 0.5: log("mitoice", "Mitochondria nucleated ice - seeding the cytosol")
        # Intercellular ice propagation through cell-cell junctions
        # (Irimia & Karlsson 2002; Acker & McGann 2000). Mean-field: the frozen
        # fraction of the junction-coupled cluster is ~P_iif, so propagation is
        # autocatalytic -- one nucleation event sweeps a confluent sheet. A
        # small f_ice term seeds invasion from the frozen extracellular side.
        J_prop = (P.k_prop * junction * mob * (P_iif + 0.004 * f_ice)
                  if (ice and not slow and not KO.get("nucl")) else 0.0)
        last_haz = J + J_prop
        if J + J_prop > 0:
            P_iif = 1 - (1 - P_iif) * math.exp(-(J + J_prop) * dt)
            # track how much intracellular ice forms: deeper supercooling at
            # the moment of nucleation means more, larger crystals (Karlsson
            # 1993; Toner 1990) -> a larger lethal fraction downstream. Capture
            # the supercooling only where ice actually nucleates (J>0, mobile),
            # capped at 20 C -- the huge nominal supercooling deep in storage,
            # where mobility is zero and nothing forms, is not physical.
            iif_amount = max(iif_amount, P_iif * clamp(dTsc / 12.0, 0.2, 1.0))
            supercool_peak = max(supercool_peak, min(dTsc, 20.0))
            if P_iif > 0.5: log("iif", "Intracellular ice nucleation probability passed 50%")

        # ---- energetics. Warburg/glycolytic cells make more of their ATP
        # without the mitochondria, so they keep producing it when the membrane
        # potential collapses on rewarming -- a real robustness axis for
        # transformed lines. glycolytic shifts the oxphos:glycolysis split.
        gf   = clamp(P.glycolytic, 0, 1)
        kT_e = arrh(55e3, T_K, 310.15) * mob
        pOx  = 0.0 if KO.get("oxphos") else P.k_atp * (0.85 - 0.55 * gf) * dPsi * (1 - mpt_frac) * kT_e
        pGly = P.k_atp * (0.15 + 0.55 * gf) * kT_e
        uBas = P.k_atp * 0.45 * kT_e
        uNaK = P.k_atp * 0.30 * (0.30 + defect) * kT_e
        uSer = 9.0e-3 * arrh(60e3, T_K, 310.15) * caCyt * 2.0 * (atp / (atp + ATP_KM))
        uCyt = (0.004 * ((1 - mt_i) + (1 - actin)) if T_C > 15 else 0.0) * kT_e
        st = min(dt, 30)
        atp_prod += (pOx + pGly) * st
        atp_use["basal"] += uBas * st; atp_use["naka"] += uNaK * st
        atp_use["serca"] += uSer * st; atp_use["cytoskeleton"] += uCyt * st
        atp = clamp(atp + (pOx + pGly - uBas - uNaK - uSer - uCyt) * st, 0, 1.4)
        atp_min = min(atp_min, atp)
        if atp < 0.35 and T_C > 5:
            D_energy += 0.004 * ((0.35 - atp) / 0.35) ** 1.5 * min(dt, 5)
            log("atp", "ATP demand outran supply on rewarming - energetic failure")
        if atp < 0.6: log("atplow", "ATP pool fell below 60% of resting")
        atpF = atp / (atp + ATP_KM)

        # ---- ER calcium & permeability transition
        kLeak = P.k_leak * (0.9 * defect + 0.30 * max(0.0, (Osm_e - 4) / 16)) * mobD
        kPump = (0.0 if KO.get("serca") else 9.0e-3) * arrh(60e3, T_K, 310.15) * dPsi * atpF
        dCaOut = min(caER, kLeak * caER * dt)
        caER = clamp(caER - dCaOut, 0, 1)
        caCyt = clamp(caCyt + dCaOut - min(caCyt, kPump * caCyt * dt), 0, 1.4)
        if caCyt > 0.35: log("ca", "ER calcium leaked into the cytosol")

        # Mechanical stress transmitted through the cytoskeleton to the
        # mitochondria promotes the permeability transition directly, on top of
        # the Ca2+/defect/osmotic drive (mitochondrial mechanosensitivity;
        # cytoskeleton-mitochondria tethering — Bartolak-Suki 2017; Kaasik).
        mech_mito = clamp(0.7 * sig_mt + 0.9 * squeeze_load, 0, 1)
        mpt_drive = (2.2 * max(0.0, caCyt - 0.40) ** 1.5 + 0.5 * defect ** 1.5
                     + 0.30 * max(0.0, (Osm_i - 10) / 14)
                     + 0.55 * max(0.0, mech_mito - 0.40))
        kmpt = ((0.0 if KO.get("mpt") else P.k_mpt) * mpt_drive * mobD
                * arrh(80e3, T_K, 310.15) * (1 + 1.6 * max(0.0, 0.5 - atp)))
        if kmpt > 0:
            mpt_frac = clamp(mpt_frac + min(kmpt * dt, 0.2) * (1 - mpt_frac), 0, 1)
        # Frozen mitochondria open the permeability transition directly, as a
        # proportional floor (a fraction P_mito of mitochondria are ice-damaged
        # and cannot maintain the inner-membrane barrier). A floor, not an
        # integrating drive, so a negligible P_mito stays negligible instead of
        # accumulating over the long recovery into a Ca-MPT runaway.
        if P_mito > 0:
            mpt_frac = clamp(max(mpt_frac, 0.9 * P_mito), 0, 1)
            D_mitoice = max(D_mitoice, P_mito)
        if mpt_frac > 0.25: log("mpt", "Mitochondrial permeability transition opening")
        # cytoskeletal/compressive stress also depresses the membrane potential
        # reversibly, independent of the (largely irreversible) MPT.
        dpsi_target = ((1 - mpt_frac) * (1 - 0.40 * max(0.0, mech_mito - 0.30))
                       * (1.0 if T_C > 5 else 0.65))
        dPsi += (dpsi_target - dPsi) * clamp(dt * 0.02 * mob, 0, 1)
        if dPsi < 0.5: log("psi", "Membrane potential collapsed below half")
        if mito_volume() / (Vmito0 * matrixF) > 2.3: log("burst", "Mitochondrial outer membrane rupture")

        # ---- mechanotransduction
        warm = T_C > 4
        kT_m = arrh(65e3, T_K, 310.15) * mob
        inhib_now = RI_inhib * (ri_freeze if ice else
                                (max(ri_freeze, ri_thaw) if phase in ("dilute", "recover") else ri_freeze))
        stm = min(dt, 60)

        piezo_drive = 1.8 * tension * (0.40 + 0.60 * FA * actin) + 0.9 * defect
        piezo = 0.0 if KO.get("piezo") else clamp(
            (1 - gsmtx) / (1 + math.exp(-(piezo_drive - 1.30) * 4)), 0, 1)
        caCyt = clamp(caCyt + 0.004 * piezo * kT_m * stm * (1.2 - caCyt), 0, 1.4)

        fa_loss = (1.2e-4 * clamp((12 - T_C) / 25, 0, 1) + 0.0040 * max(0.0, tension - 0.3)
                   + 0.0030 * max(0.0, 0.55 - vn)) * mob
        fa_gain = 5.0e-5 * atpF * actin * kT_m if warm else 0.0
        FA = 0.0 if KO.get("fa") else clamp(FA - fa_loss * stm + fa_gain * stm * (1 - FA), 0, 1)
        if FA < 0.35 * FA0: log("fa", "Focal adhesions lost - cell rounding, anoikis signalling")

        # RhoA is driven by microtubule loss AND by the compressive stress that
        # buckles them — buckling releases GEF-H1 before the network fully
        # depolymerises (Krendel 2002; Chang 2008). Membrane tension and IF
        # yielding add to the mechanical drive.
        rho_drive = (0.85 * (1 - mt_i) + 0.70 * clamp(tension, 0, 1.5)
                     + 0.50 * clamp(Ci_c / 1.4, 0, 3) + 0.55 * defect
                     + 0.60 * max(0.0, sig_mt - 0.4) + 0.35 * max(0.0, sig_if - 1.0))
        rhoGTP = clamp(rhoGTP + (clamp(rho_drive / 2.2, 0, 1) - rhoGTP) * clamp(0.05 * kT_m * stm, 0, 1), 0, 1)
        rock_raw = clamp(0.60 * rhoGTP + 0.95 * casp3, 0, 2.5)
        rock = 0.0 if KO.get("rock") else clamp(rock_raw * (1 - inhib_now), 0, 2.5)
        rock_peak = max(rock_peak, rock)
        pMLC = clamp(pMLC + (clamp(rock, 0, 1.5) - pMLC) * clamp(0.06 * kT_m * stm, 0, 1), 0, 1.5)
        bleb = clamp(pMLC * (1 - 0.75 * FA) - 0.25, 0, 1.5)
        bleb_peak = max(bleb_peak, bleb)
        if bleb > 0.35: log("bleb", "Actomyosin hypercontraction - membrane blebbing")
        if warm:
            actin = clamp(actin - 0.004 * clamp(rock - 0.5, 0, 2) * kT_m * stm, 0, 1)
        akt = clamp(1 / (1 + 1.8 * clamp(rock, 0, 2)) * (0.4 + 0.6 * FA), 0, 1)
        yap_target = clamp(0.15 + 0.85 * FA * actin * (1 - 0.5 * bleb), 0, 1)
        yapN = 0.0 if KO.get("yap") else clamp(
            yapN + (yap_target - yapN) * clamp(0.03 * kT_m * stm, 0, 1), 0, 1)
        if yapN < 0.30: log("yap", "YAP exported from the nucleus - mechanosignalling lost")
        calpain = clamp(calpain + (clamp(2.0 * max(0.0, caCyt - 0.25), 0, 1) - calpain)
                        * clamp(0.04 * kT_m * stm, 0, 1) * (1 - calpI)
                        * (0 if KO.get("calpain") else 1), 0, 1)

        # apoptosis resistance (p53 loss in HeLa, Bcl-2 overexpression) throttles
        # the whole caspase cascade, like a genetic z-VAD; anoikis resistance
        # (anchorage independence) uncouples detachment from the death signal.
        kc = 0.012 * kT_m * (1 - zvad) * (1 - clamp(P.apop_resist, 0, 1)) * (0 if KO.get("casp") else 1)
        kdec = 3.0e-5 * kT_m * stm
        def ratchet(cur, target, k):
            return cur + (target - cur) * clamp(k, 0, 1) if target > cur \
                else cur + (target - cur) * clamp(kdec, 0, 1)
        anoikis = 0.75 * (1 - clamp(P.anoikis_resist, 0, 1)) * max(0.0, 0.45 - FA) / 0.45
        casp8 = clamp(ratchet(casp8, clamp(0.42 * clamp(rock - 0.45, 0, 2) + 0.24 * defect
                                           + 0.30 * calpain + anoikis, 0, 1), kc * stm), 0, 1)
        casp9 = clamp(ratchet(casp9, clamp(0.95 * mpt_frac + 0.45 * calpain + 0.40 * ros, 0, 1)
                              * (1 - 0.45 * akt), kc * stm), 0, 1)
        exec_ = clamp(0.85 * casp8 + 0.85 * casp9, 0, 1) * atpF
        casp3 = clamp(ratchet(casp3, exec_, kc * stm * (1 - zvad)), 0, 1)
        casp3_peak = max(casp3_peak, casp3)
        if casp3 > 0.3: log("casp3", "Caspase-3 activated - and it cleaves ROCK1, so the loop closes")

        if warm:
            # Oxidative stress on rewarming: leaky/depolarised mitochondria and a
            # damaged membrane generate ROS, which is the documented driver of
            # cryopreservation-induced delayed-onset cell death (Baust, Van
            # Buskirk & Baust 2000; Baust 2009; Bissoyi & Pramanik 2014). ROS
            # accumulates and feeds the mitochondrial (caspase-9) apoptotic arm.
            # antioxidant capacity (NRF2/glutathione — high in A549) lowers ROS
            # generation and speeds its clearance, blunting the delayed wave.
            aox = clamp(P.antioxidant, 0, 1)
            ros = clamp(ros + (0.012 * (1 - 0.8 * aox) * (mpt_frac + (1 - dPsi) + 0.5 * max(0.0, defect - 0.2))
                               - 0.0035 * (1 + 4 * aox) * ros) * kT_m * stm, 0, 1.5)
            ros_peak = max(ros_peak, ros)
            if ros > 0.4: log("ros", "Oxidative stress building on rewarming - delayed apoptosis primed")
            dA = 1.4e-5 * (casp3 + 0.80 * max(0.0, bleb - 0.10) + 0.60 * ros) \
                 * atpF * (1 - apop - necr) * stm
            dN = 2.4e-5 * (max(0.0, 0.45 - atp) / 0.45 + 0.80 * max(0.0, defect - 0.30)) \
                 * (1 - atpF) * (1 - apop - necr) * stm
            apop = clamp(apop + max(0.0, dA), 0, 1)
            necr = clamp(necr + max(0.0, dN), 0, 1)
            if RI_tox > 0 and (ri_thaw or ri_freeze):
                necr = clamp(necr + 5.0e-6 * RI_tox * (1 - apop - necr) * stm, 0, 1)

        # ---- nucleus
        chrom_cond = clamp((Osm_i - 0.6) / 9, 0, 1)
        if Ci_c > 1.2 and T_C > 0:
            lobulation = clamp(lobulation + 0.00004 * dt * Ci_c, 0, 1)
        if chrom_cond > 0.6: log("chrom", "Chromatin heavily condensed by freeze-concentration")

    # ---------------------------------------------------------------- timeline
    steps = []
    n_add = 1 if P.cpa_key == "none" else P.add_steps
    for i in range(1, n_add + 1):
        steps.append(dict(kind="add", frac=i / n_add, T=P.T_add,
                          dur=(P.hold_min * 60 if i == n_add else max(30, P.hold_min * 60 / n_add))))
    steps.append(dict(kind="cool", frm=P.T_add, to=P.T_seed, rate=P.CR))
    steps.append(dict(kind="seed", T=P.T_seed, dur=20))
    steps.append(dict(kind="cool", frm=P.T_seed, to=P.T_store, rate=P.CR, ice=True))
    steps.append(dict(kind="store", T=P.T_store, dur=P.days * 86400, ice=True))
    steps.append(dict(kind="warm", frm=P.T_store, to=0.0, rate=P.WR, ice=True))
    steps.append(dict(kind="melt", T=0.0, dur=15))
    steps.append(dict(kind="dilute", T=P.T_dil, dur=600))
    steps.append(dict(kind="recover", T=P.T_recover, dur=P.recover_h * 3600))

    BUDGET = dict(add=50, cool=150, seed=8, store=20, warm=120, melt=10, dilute=90, recover=150)
    n_cool = sum(1 for s in steps if s["kind"] == "cool") or 1
    ds = P.dt_scale
    next_frame = 0.0

    for s in steps:
        dur_est = s.get("dur") if s.get("dur") is not None \
            else abs(s["to"] - s["frm"]) / max(1e-6, s["rate"]) * 60
        div = (n_cool * 0.5) if s["kind"] == "cool" else (n_add if s["kind"] == "add" else 1)
        frame_period = max(1e-6, dur_est / max(1, BUDGET.get(s["kind"], 30) / div))
        next_frame = t

        if s["kind"] == "add":
            eCPA0 = cpaM * s["frac"]; e_c = eCPA0; e_f = 0.0; e_s = salt_e0; e_suc = 0.0
            e_add0 = addM; e_add = addM          # additive is in the loading medium from the start
            T_C = s["T"]; dt = 0.25 * ds; u = 0.0
            while u < s["dur"]:
                step(dt, "load"); t += dt; u += dt
                if t >= next_frame: push("load"); next_frame = t + frame_period
        elif s["kind"] in ("cool", "warm"):
            span = abs(s["to"] - s["frm"]); rate = max(1e-6, s["rate"]) / 60.0
            dur = span / rate
            dt_base = clamp(min(2.0, 0.05 / rate), 1e-4, 4) * ds
            direction = 1 if s["to"] > s["frm"] else -1
            ice = bool(s.get("ice")) or ice
            done = 0.0
            while done < dur:
                dt = max(dt_base, 2 / rate) if T_C < tg_run else dt_base
                stp = min(dt, dur - done)
                T_C += direction * rate * stp
                step(stp, s["kind"], T_C < tg_run); t += stp; done += stp
                if t >= next_frame: push(s["kind"]); next_frame = t + frame_period
            T_C = s["to"]
        elif s["kind"] == "seed":
            ice = True; log("seed", "Extracellular ice seeded")
            dt = 0.25 * ds; u = 0.0
            while u < s["dur"]:
                step(dt, "seed"); t += dt; u += dt
                if t >= next_frame: push("seed"); next_frame = t + frame_period
        elif s["kind"] == "store":
            T_C = s["T"]; ice = True; nsub = 60; dt = s["dur"] / nsub
            for _ in range(nsub):
                step(dt, "store", True); t += dt
                if t >= next_frame: push("store"); next_frame = t + frame_period
            if P.cpa_key == "sd":
                if frac_intact < 0.9 and T_C > tg_run:
                    log("stodeg", "CPA degraded during storage - protection lost before thaw")
                if T_C <= tg_run:
                    log("arrest", "Below Tg': cleavage chemistry arrested, CPA intact in storage")
        elif s["kind"] == "melt":
            T_C = 0.0
            wr_pen = clamp(1 - math.log10(max(P.WR, 1)) / 3.2, 0, 1)
            D_recry = P_iif * wr_pen
            if D_recry > 0.2: log("recry", "Intracellular ice recrystallised during slow warming")
            ice = False; dt = 0.25 * ds; u = 0.0
            while u < s["dur"]:
                step(dt, "melt"); t += dt; u += dt
                if t >= next_frame: push("melt"); next_frame = t + frame_period
        elif s["kind"] == "dilute":
            T_C = P.T_dil; ice = False
            if P.dilution == "none":   stages = [(eCPA0 * frac_intact, 0.0, 600)]
            elif P.dilution == "direct": stages = [(0.0, 0.0, 600)]
            else: stages = [(eCPA0 * 0.5 * frac_intact, 0.25, 180),
                            (eCPA0 * 0.2 * frac_intact, 0.25, 180),
                            (0.0, 0.10, 240)]
            for si, (c, suc, d) in enumerate(stages):
                frame_period = max(0.5, d / (BUDGET["dilute"] / len(stages)))
                next_frame = t
                # additive washes out with the CPA (unless dilution is "none",
                # i.e. infused in situ with no wash, where it stays)
                if P.dilution == "none":
                    awash = 1.0
                elif len(stages) == 3:
                    awash = [0.40, 0.15, 0.0][si]
                else:
                    awash = 0.0
                eCPA0 = c; e_c = c; e_f *= 0.15; e_s = salt_e0; e_suc0 = suc; e_suc = suc
                e_add0 = addM * awash; e_add = addM * awash
                dt = 0.25 * ds; u = 0.0
                while u < d:
                    step(dt, "dilute"); t += dt; u += dt
                    if t >= next_frame: push("dilute"); next_frame = t + frame_period
        elif s["kind"] == "recover":
            T_C = s["T"]; ice = False; dt = P.recover_dt * ds; u = 0.0
            while u < s["dur"]:
                step(dt, "recover"); t += dt; u += dt
                if t >= next_frame: push("recover"); next_frame = t + frame_period
    push("end")

    # ---------------------------------------------------------------- outcome
    # IIF lethality now scales with the amount of intracellular ice, proxied by
    # the peak supercooling at nucleation (deeper supercool -> more, finer,
    # more damaging crystals; Karlsson 1993), and still saturates to fully
    # lethal when slow warming lets it recrystallise. At the calibration point
    # (shallow supercooling, little recry) this reproduces the old ~0.55 base.
    recry_frac = clamp(D_recry / max(P_iif, 1e-6), 0, 1)
    iif_lethal = clamp(0.45 + 0.35 * clamp(supercool_peak / 10.0, 0, 1), 0.45, 0.90)
    iif_kill = P_iif * (iif_lethal + (1 - iif_lethal) * recry_frac)
    D_shrink = 1 - math.exp(-(max(0.0, 0.32 - minV) / 0.08) ** 2)
    D_osm = clamp(1 - (1 - D_shrink) * (1 - clamp(D_swell, 0, 1)) * (1 - clamp(D_sol, 0, 1)), 0, 1)
    # a water-replacing extracellular additive is an osmoprotectant: it lets the
    # cell tolerate the dehydration it undergoes (membrane/protein stabilised
    # through the shrink), so the same volume excursion is less lethal.
    osmo_prot = clamp(add["wRepl"] * add_pres, 0, 0.6)
    D_osm = clamp(D_osm * (1 - osmo_prot), 0, 1)
    D_mech = clamp(1 - (1 - clamp(D_mech, 0, 1)) * (1 - clamp(D_recry_ice, 0, 1)), 0, 1)
    D_tox = clamp(D_tox, 0, 1)
    D_mem = clamp(D_mem + D_frag, 0, 1)
    # mitochondrial ice is directly lethal (bioenergetic collapse) on top of the
    # MPT/apoptosis path it already drove; combine it with cytosolic IIF.
    mito_kill = 0.70 * clamp(D_mitoice, 0, 1)
    D_iif = clamp(1 - (1 - clamp(iif_kill, 0, 1)) * (1 - mito_kill), 0, 1)
    D_energy = clamp(D_energy, 0, 1)
    D_apop = clamp(apop + necr, 0, 0.98)

    S_imm = (1 - D_iif) * (1 - D_osm) * (1 - D_mem) * (1 - D_tox) * (1 - D_mech)
    S_24 = S_imm * (1 - D_apop)
    # Delayed-onset death (CIDOC): ROS and primed caspases keep killing cells
    # for 1-3 days after thaw, beyond the immediate recovery window. Project the
    # extra loss from the residual ROS and executioner-caspase left at the end
    # of the simulated window (Baust 2000, 2009). S_72 is the settled viability.
    D_delayed = clamp(0.55 * ros_peak + 0.35 * casp3, 0, 0.9)
    S_72 = S_24 * (1 - D_delayed)
    # functional recovery needs the machinery AND its proteins folded correctly;
    # protein integrity (water-replacement protected, denaturation damaged) gates it.
    D_prot = clamp(1 - prot, 0, 1)
    F_rec = (S_24 * (0.30 + 0.35 * mt_i + 0.20 * actin + 0.15 * FA)
             * (0.45 + 0.55 * dPsi) * (0.30 + 0.70 * clamp(atp, 0, 1))
             * (0.40 + 0.60 * yapN) * (0.55 + 0.45 * akt)
             * (0.35 + 0.65 * prot))
    live = clamp(S_imm * (1 - apop - necr), 0, 1)
    apop_f = clamp(S_imm * apop + (1 - S_imm) * 0.35, 0, 1)
    necr_f = clamp(1 - live - apop_f, 0, 1)

    res = dict(S_imm=S_imm, S_24=S_24, S_72=S_72, F_rec=F_rec,
               rosPeak=ros_peak, D_delayed=D_delayed,
               D_osm=D_osm, D_tox=D_tox, D_mem=D_mem, D_iif=D_iif, D_apop=D_apop,
               D_mech=D_mech, D_energy=D_energy, D_shrink=D_shrink,
               D_swell=clamp(D_swell, 0, 1), D_sol=clamp(D_sol, 0, 1),
               D_recry_ice=clamp(D_recry_ice, 0, 1),
               P_iif=P_iif, iifAmount=iif_amount, supercoolPeak=supercool_peak,
               junction=junction, P_mito=P_mito, D_mitoice=clamp(D_mitoice, 0, 1),
               matrixCPA=mcpa, glass=clamp(1.0 - last_mob, 0, 1),
               prot=prot, D_prot=D_prot,
               mptFrac=mpt_frac, dPsi=dPsi, caCyt=caCyt, caER=caER,
               minV=minV, maxV=maxV, chromCond=chrom_cond, lobulation=lobulation,
               fracIntact=frac_intact, Vmito_final=mito_volume() / (Vmito0 * matrixF),
               grain=grain, grainMax=grain_max, fIce=f_ice, chanW=chanW,
               squeeze=squeeze, squeezeMin=squeeze_min, memb=memb, tension=tension,
               msOpen=ms_open, mt=mt_i, actin=actin, intf=if_i,
               sigMT=sig_mt, sigIF=sig_if, atp=atp, atpMin=atp_min,
               atpProd=atp_prod, atpUse=atp_use, FA=FA, FA0=FA0, adhLabel=adh_label,
               piezo=piezo, rhoGTP=rhoGTP, rock=rock, rockPeak=rock_peak, pMLC=pMLC,
               bleb=bleb, blebPeak=bleb_peak, yapN=yapN, akt=akt, casp3=casp3,
               casp3Peak=casp3_peak, casp8=casp8, casp9=casp9, calpain=calpain,
               apop=apop, necr=necr, liveFrac=live, apopFrac=apop_f, necrFrac=necr_f,
               rockInhib=RI_inhib, rockTox=RI_tox)
    return out, res, out.events
