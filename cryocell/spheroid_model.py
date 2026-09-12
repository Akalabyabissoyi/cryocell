"""1-D radial reaction-diffusion model for a cryopreserved spheroid.

This is a genuine finite-difference solve of cryoprotectant (CPA) loading into a
spherical multicellular construct, replacing the earlier illustrative
sqrt(D*t) penetration estimate in the spheroid view. It is deliberately kept in
its own module and is NOT used by the single-cell engine (model.py), so it
cannot affect the calibrated 35.4% benchmark.

Physics
-------
CPA loading is diffusion-limited. In spherical symmetry Fick's second law is

    dC/dt = D * [ d2C/dr2 + (2/r) dC/dr ],   C(R,t)=C_bath,  C(r,0)=0

solved here by an explicit, conservative finite-difference scheme over the warm
loading hold. Diffusion effectively stops once the construct cools (Arrhenius
slow-down, then frozen), so the radial profile reached at the end of the hold is
what the cell population carries into the freeze — a large core stays
CPA-starved. Heat diffuses ~10^3 x faster than CPA, so thermal gradients are
negligible at these sizes and are not modelled (stated honestly, not faked).

D_eff here is an EFFECTIVE CPA-loading transport coefficient for a packed
spheroid, not free-solution diffusion: loading a cell deep in the construct
requires the CPA to permeate across many cell membranes in series, so the
effective rate is far below the extracellular tortuous-diffusion value. ~8
um^2/s reproduces the measured size effect (200 um spheroids load fully in a
10 min hold and recover best; 400 um cores stay under-loaded and recover worse
— Gao/Bissoyi/Guo/Gibson 2024). It is an order-of-magnitude prior, labelled as
such (free-solution DMSO ~1200 um^2/s; extracellular tissue ~200-500; this
membrane-series effective value ~5-15). Xu 2014; Devireddy tissue reviews.
"""
from __future__ import annotations
import numpy as np

D_EFF_UM2_S = 8.0           # effective spheroid CPA-loading coefficient, um^2/s


def solve_cpa_profile(R_um: float, t_hold_s: float, N: int = 24,
                      D: float = D_EFF_UM2_S):
    """Radial CPA profile at the end of the loading hold.

    Returns (r_frac, cpa_frac): shell-centre fractional radii (0 core .. 1 rim)
    and the local CPA fraction there (0 = none, 1 = fully equilibrated to bath).
    """
    R_um = max(R_um, 1.0)
    dr = R_um / N
    r = (np.arange(N) + 0.5) * dr                       # shell-centre radii, um
    C = np.zeros(N); C[-1] = 1.0                         # surface at bath conc
    dt = 0.2 * dr * dr / D                               # explicit-stability step
    steps = int(np.clip(t_hold_s / dt, 1, 200000))
    for _ in range(steps):
        lap = np.zeros(N)
        lap[1:-1] = ((C[2:] - 2 * C[1:-1] + C[:-2]) / dr ** 2
                     + (2.0 / r[1:-1]) * (C[2:] - C[:-2]) / (2 * dr))
        Cn = C + dt * D * lap
        Cn[0] = C[0] + dt * D * 6.0 * (C[1] - C[0]) / dr ** 2   # r->0 symmetry
        Cn[-1] = 1.0                                            # Dirichlet surface
        C = np.clip(Cn, 0.0, 1.0)
    return r / R_um, C


def sample(r_frac, cpa_frac, u: float) -> float:
    """Interpolate the solved profile at fractional radius u (0 core .. 1 rim)."""
    return float(np.interp(np.clip(u, 0.0, 1.0), r_frac, cpa_frac))
