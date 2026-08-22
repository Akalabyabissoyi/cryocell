"""
Population statistics, in silico perturbation screens, self-assessment.

Structure follows Karr et al. (Cell 2012) for the population runs and the
disruption screen, and Cui et al. (Nature 2025) for uncertainty quantification
and the lab-in-the-loop next-experiment recommender.

Everything here parallelises across cores, because a 150-cell population is
150 independent simulations and there is no reason to run them one at a time.
"""
from __future__ import annotations
import math, os
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
from .model import Params, simulate, PERTURBATIONS, clamp

# cell-to-cell heterogeneity: key, label, CV, distribution
HETEROGENEITY = [
    ("Viso",      "Cell volume",             0.22, "lognormal"),
    ("lp",        "Water permeability Lp",   0.25, "lognormal"),
    ("ps",        "CPA permeability Ps",     0.25, "lognormal"),
    ("vb",        "Inactive volume Vb",      0.12, "normal"),
    ("sterol",    "Membrane sterol",         0.18, "normal"),
    ("nuc_scale", "Nucleation site density", 0.90, "lognormal"),
]

CALIBRATION_TARGETS = [
    ("lp",           "Lp of your MSC line",            0.60, "estimated",
     "Coulter shrinkage in hypertonic non-permeant, 3 temperatures"),
    ("ps",           "Ps for the candidate CPA",       0.60, "estimated",
     "shrink-swell kinetics on CPA addition, 3 temperatures"),
    ("vb",           "Osmotically inactive volume Vb", 0.25, "measured",
     "Boyle-van't Hoff plot"),
    ("Viso",         "Isotonic cell volume",           0.20, "measured",
     "Coulter counter, isotonic"),
    ("sterol",       "Membrane sterol content",        0.40, "estimated",
     "lipidomics or filipin / Amplex cholesterol assay"),
    ("nuc_visc_exp", "Nucleation viscosity exponent",  0.35, "fitted",
     "cryomicroscopy: IIF fraction vs cooling rate"),
    ("nuc_scale",    "Nucleation site density",        1.50, "fitted",
     "cryomicroscopy: nucleation temperature distribution across cells"),
    ("k_ripen",      "Ice grain-growth constant",      0.80, "fitted",
     "splat / sucrose-sandwich assay, grain size vs annealing time"),
    ("k_squeeze",    "Mechanical squeeze coefficient", 1.00, "fitted",
     "directional-solidification stage: survival vs interface velocity"),
    ("mob_exp",      "Mobility exponent for damage",   0.50, "fitted",
     "storage-stability series at -80 C over months"),
    ("k_tox",        "CPA toxicity coefficient",       0.70, "fitted",
     "viability vs CPA exposure time at 2-3 temperatures, no freezing"),
    ("k_mpt",        "Permeability-transition rate",   1.00, "estimated",
     "post-thaw dPsi (TMRM) time course +/- cyclosporin A"),
    ("k_leak",       "ER calcium leak rate",           1.00, "estimated",
     "Fura-2 / GCaMP cytosolic Ca2+ through freeze-thaw"),
    ("k_atp",        "ATP turnover rate",              0.60, "estimated",
     "luciferase ATP assay through rewarming"),
]


def _run(kwargs):
    """Top-level so it can be pickled for the process pool."""
    return simulate(Params(**kwargs))[1]


def _pool_map(fn, items, workers=None):
    n = workers or max(1, (os.cpu_count() or 4) - 1)
    if n <= 1 or len(items) < 4:
        return [fn(i) for i in items]
    with ProcessPoolExecutor(max_workers=n) as ex:
        return list(ex.map(fn, items, chunksize=max(1, len(items) // (n * 3))))


# --------------------------------------------------------------- population
class _RNG:
    def __init__(self, seed): self.a = seed & 0xFFFFFFFF
    def next(self):
        self.a = (self.a + 0x6D2B79F5) & 0xFFFFFFFF
        t = (self.a ^ (self.a >> 15)) * (1 | self.a) & 0xFFFFFFFF
        t = (t + ((t ^ (t >> 7)) * (61 | t) & 0xFFFFFFFF)) ^ t & 0xFFFFFFFF
        return ((t ^ (t >> 14)) & 0xFFFFFFFF) / 4294967296.0
    def gauss(self):
        u = self.next() or 1e-9; v = self.next() or 1e-9
        return math.sqrt(-2 * math.log(u)) * math.cos(2 * math.pi * v)


def simulate_population(P: Params, n=150, seed=20260816, progress=None):
    r = _RNG(seed)
    base = asdict(P); base["ko"] = dict(P.ko or {})
    base["recover_dt"] = max(base.get("recover_dt", 30), 30)
    jobs, samples = [], {k: [] for k, _, _, _ in HETEROGENEITY}
    draws = []
    for _ in range(n):
        q = dict(base); q["ko"] = dict(base["ko"])
        d = {}
        for key, _lab, cv, kind in HETEROGENEITY:
            b = 1.0 if key == "nuc_scale" else base[key]
            if kind == "lognormal":
                v = b * math.exp(r.gauss() * math.sqrt(math.log(1 + cv * cv)))
            else:
                v = max(0.0, b * (1 + r.gauss() * cv))
            if key == "vb":     v = clamp(v, 0.05, 0.6)
            if key == "sterol": v = clamp(v, 0.0, 60.0)
            q[key] = v; samples[key].append(v); d[key] = v
        draws.append(r.next())
        jobs.append(q)

    results = _pool_map(_run, jobs)

    cells = []
    for R, u in zip(results, draws):
        lethality = R["D_iif"] / R["P_iif"] if R["P_iif"] > 1e-9 else 0.55
        b = (1 - R["D_osm"]) * (1 - R["D_mem"]) * (1 - R["D_tox"]) \
            * (1 - R["D_mech"]) * (1 - R["D_apop"])
        nucleated = u < R["P_iif"]
        cells.append(dict(s24=clamp(b * ((1 - lethality) if nucleated else 1.0), 0, 1),
                          frec=R["F_rec"], nucleated=nucleated,
                          minV=R["minV"], atpMin=R["atpMin"], dPsi=R["dPsi"]))

    s = sorted(c["s24"] for c in cells)
    q = lambda f: s[clamp(int(f * (len(s) - 1)), 0, len(s) - 1)]

    def rank(a):
        order = sorted(range(len(a)), key=lambda i: a[i])
        rk = [0] * len(a)
        for j, i in enumerate(order): rk[i] = j
        return rk
    ys = rank([c["s24"] for c in cells])
    attribution = []
    for key, lab, _cv, _k in HETEROGENEITY:
        xs = rank(samples[key])
        mx = (len(xs) - 1) / 2; my = (len(ys) - 1) / 2
        num = dx = dy = 0.0
        for i in range(len(xs)):
            a_ = xs[i] - mx; b_ = ys[i] - my
            num += a_ * b_; dx += a_ * a_; dy += b_ * b_
        attribution.append((key, lab, num / math.sqrt(dx * dy) if dx and dy else 0.0))
    attribution.sort(key=lambda z: -abs(z[2]))

    bins = [0] * 20
    for c in cells: bins[clamp(int(c["s24"] * 20), 0, 19)] += 1
    return dict(n=n, cells=cells, bins=bins, samples=samples,
                mean=sum(s) / len(s), median=q(0.5), p10=q(0.10), p25=q(0.25),
                p75=q(0.75), p90=q(0.90),
                nucleated=sum(1 for c in cells if c["nucleated"]) / n,
                frec_mean=sum(c["frec"] for c in cells) / n,
                attribution=attribution)


# ----------------------------------------------------------------- screens
def knockout_screen(P: Params):
    base = asdict(P); base["ko"] = dict(P.ko or {})
    jobs = [base]
    for pid, _n, _b in PERTURBATIONS:
        q = dict(base); q["ko"] = dict(base["ko"]); q["ko"][pid] = True
        jobs.append(q)
    res = _pool_map(_run, jobs)
    wt = res[0]
    rows = []
    for (pid, name, bench), R in zip(PERTURBATIONS, res[1:]):
        d = R["S_24"] - wt["S_24"]
        cls = ("essential" if d < -0.20 else "buffering" if d < -0.05
               else "protective target" if d > 0.05 else "neutral")
        rows.append(dict(id=pid, name=name, bench=bench, ko=R["S_24"], delta=d,
                         cls=cls, dfrec=R["F_rec"] - wt["F_rec"]))
    rows.sort(key=lambda r: r["delta"])
    return dict(wt=wt["S_24"], rows=rows,
                headroom=0.05 < wt["S_24"] < 0.95)


def stability_check(P: Params):
    base = asdict(P); base["ko"] = dict(P.ko or {})
    jobs = []
    for ds in (1.0, 0.5, 2.0):
        q = dict(base); q["ko"] = dict(base["ko"]); q["dt_scale"] = ds; jobs.append(q)
    b, h, d = [r["S_24"] for r in _pool_map(_run, jobs)]
    drift = max(abs(h - b), abs(d - b))
    return dict(base=b, half=h, dbl=d, drift=drift,
                verdict="stable" if drift < 0.02 else "marginal" if drift < 0.05 else "unstable")


_DEFAULTS = dict(nuc_scale=1.0, nuc_visc_exp=3.0, k_ripen=80.0, k_squeeze=9.0e-4,
                 mob_exp=2.0, k_tox=1.13e-4, k_mpt=0.003, k_leak=4.0e-4, k_atp=0.033)


def next_experiment(P: Params):
    base = asdict(P); base["ko"] = dict(P.ko or {})
    jobs = [base]
    for key, _lab, u, _t, _how in CALIBRATION_TARGETS:
        cur = _DEFAULTS.get(key, base.get(key))
        for f in (1 - min(u, 0.95), 1 + u):
            q = dict(base); q["ko"] = dict(base["ko"]); q[key] = cur * f; jobs.append(q)
    res = _pool_map(_run, jobs)
    base_s = res[0]["S_24"]
    rows = []
    for i, (key, lab, u, tier, how) in enumerate(CALIBRATION_TARGETS):
        lo = res[1 + 2 * i]["S_24"]; hi = res[2 + 2 * i]["S_24"]
        rows.append(dict(key=key, label=lab, tier=tier, how=how, u=u,
                         lo=lo, hi=hi, influence=abs(hi - lo)))
    rows.sort(key=lambda r: -r["influence"])
    tot = sum(r["influence"] for r in rows) or 1.0
    for r in rows: r["share"] = r["influence"] / tot
    band = math.sqrt(sum((r["influence"] / 2) ** 2 for r in rows))
    return dict(base=base_s, rows=rows, band=band, confident=band < 0.08,
                verdict=("The prediction is tight enough to act on." if band < 0.08 else
                         "Treat the prediction as a ranking, not a number." if band < 0.20 else
                         "Too uncertain to quote. Measure the top rows before believing this."))
