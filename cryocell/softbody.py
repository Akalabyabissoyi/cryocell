"""
Deformable cell mechanics.

A closed contour of point masses. The physics that matters:

  * tension springs carry a REST LENGTH taken from the membrane's own contour
    length, which is conserved when the cell dehydrates while the enclosed area
    is not. The contour has nowhere to go but buckle -- crenation is a
    consequence of the osmotic module, not a drawing decision. DMSO makes it
    worse because area per lipid grows +0.85%/mol%.
  * bending stiffness sets the fold wavelength. It has a narrow useful window:
    swept numerically, above ~0.030 the membrane stops expressing slack as folds
    at all, so it is clamped into [0.004, 0.030].
  * an actin cortex with a DEAD-BAND. Pointwise restoring irons the folds flat;
    purely smoothed restoring lets spikes grow into a star. Folds are free up to
    an amplitude set by the slack, the cortex bites past it, and a weak pull on
    the smoothed profile keeps the cell globally round.
  * an area constraint driven by the simulated cell volume.
  * confinement by the unfrozen pocket, applied radially in every direction.
"""
from __future__ import annotations
import numpy as np


class SoftBody:
    def __init__(self, n=160, radius=100.0, seed=3):
        rs = np.random.RandomState(seed)
        a = np.linspace(0, 2 * np.pi, n, endpoint=False)
        # buckling is an instability: a perfect circle under symmetric forces
        # stays a perfect circle, so the radii carry a broadband seed
        j = radius * (1 + 0.012 * rs.uniform(-1, 1, n))
        self.n = n
        self.p = np.stack([np.cos(a) * j, np.sin(a) * j], axis=1)
        self.v = np.zeros((n, 2))
        self.blebs = []          # dicts: i, w, amp, amp_max, t, life
        self.rs = rs

    # ------------------------------------------------------------------ core
    def area(self):
        x, y = self.p[:, 0], self.p[:, 1]
        return 0.5 * abs(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))

    def perimeter(self):
        d = np.roll(self.p, -1, axis=0) - self.p
        return float(np.sum(np.hypot(d[:, 0], d[:, 1])))

    def radii(self):
        return np.hypot(self.p[:, 0], self.p[:, 1])

    def radius_at(self, ang):
        """Outline radius in given directions, by nearest-vertex angle."""
        ang = np.atleast_1d(ang)
        va = np.arctan2(self.p[:, 1], self.p[:, 0])
        d = np.abs((va[None, :] - ang[:, None] + 3 * np.pi) % (2 * np.pi) - np.pi)
        return self.radii()[np.argmin(d, axis=1)]

    def step(self, L0, kT, kB, kP, A_target, damp, slack,
             k_cortex=0.0, pocket_r=None, crystals=None, noise=0.0,
             external=None):
        n = self.n
        p = self.p
        F = np.zeros_like(p)
        A = self.area()
        press = kP * (A_target - A) / max(A_target, 1.0)

        prev = np.roll(p, 1, axis=0); nxt = np.roll(p, -1, axis=0)

        # cortex softening under blebs
        stiff = np.full(n, kT)
        for b in self.blebs:
            idx = (np.arange(n) - b["i"] + n) % n
            dd = np.minimum(idx, n - idx)
            stiff[dd < b["w"]] *= 0.18

        for nb in (prev, nxt):
            d = nb - p
            L = np.hypot(d[:, 0], d[:, 1]); L[L < 1e-6] = 1e-6
            f = (stiff * (L - L0) / L)[:, None]
            F += f * d

        # bending: discrete Laplacian
        F += kB * (prev + nxt - 2 * p)

        # outward normals from the local tangent
        t = nxt - prev
        tl = np.hypot(t[:, 0], t[:, 1]); tl[tl < 1e-6] = 1e-6
        nrm = np.stack([t[:, 1] / tl, -t[:, 0] / tl], axis=1)
        flip = (p[:, 0] * nrm[:, 0] + p[:, 1] * nrm[:, 1]) < 0
        nrm[flip] *= -1
        F += press * nrm

        # blebs push outward against the detached cortex
        for b in self.blebs:
            idx = (np.arange(n) - b["i"] + n) % n
            dd = np.minimum(idx, n - idx)
            m = dd < b["w"]
            g = b["amp"] * np.cos(dd[m] / b["w"] * np.pi / 2)
            F[m] += g[:, None] * nrm[m]

        # actin cortex: dead-band pointwise + weak pull on the smoothed profile
        if k_cortex:
            Rt = np.sqrt(A_target / np.pi)
            rad = self.radii(); rad[rad < 1e-6] = 1e-6
            w = max(3, n // 12)
            k = np.ones(2 * w + 1) / (2 * w + 1)
            sm = np.convolve(np.r_[rad[-w:], rad, rad[:w]], k, mode="valid")[:n]
            band = Rt * (0.012 + 0.60 * min(max(slack, 0.0), 0.6))
            dev = rad - Rt
            f = np.zeros(n)
            f[dev > band] = -k_cortex * (dev[dev > band] - band) / rad[dev > band]
            f[dev < -band] = -k_cortex * (dev[dev < -band] + band) / rad[dev < -band]
            f += 0.35 * k_cortex * (Rt - sm) / rad
            F += f[:, None] * p

        # intracellular ice pushes the membrane out from within
        if crystals is not None and len(crystals):
            for cx, cy, cr in crystals:
                d = p - np.array([cx, cy])
                dist = np.hypot(d[:, 0], d[:, 1]); dist[dist < 1e-6] = 1e-6
                m = dist < cr
                if m.any():
                    F[m] += ((cr - dist[m]) * 1.4 / dist[m])[:, None] * d[m]

        # confinement by the unfrozen pocket, from every direction
        if pocket_r is not None:
            ang = np.arctan2(p[:, 1], p[:, 0])
            lim = pocket_r(ang)
            rad = self.radii(); rad[rad < 1e-6] = 1e-6
            over = rad - lim
            m = over > 0
            if m.any():
                F[m] -= (0.55 * over[m] / rad[m])[:, None] * p[m]

        if external is not None:
            F += external

        if noise:
            F += (self.rs.random_sample((n, 2)) - 0.5) * noise

        self.v = (self.v + F) * damp
        sp = np.hypot(self.v[:, 0], self.v[:, 1])
        big = sp > 6.0
        if big.any(): self.v[big] *= (6.0 / sp[big])[:, None]
        self.p += self.v
        # keep registered on the origin so mapped organelles do not drift out
        self.p -= self.p.mean(axis=0)

    # ----------------------------------------------------------------- blebs
    def update_blebs(self, bleb_index, max_blebs=10):
        if bleb_index > 0.12 and self.rs.random_sample() < min(bleb_index * 0.22, 0.6) \
                and len(self.blebs) < max_blebs:
            self.blebs.append(dict(i=int(self.rs.randint(self.n)),
                                   w=int(5 + self.rs.randint(6)),
                                   amp=0.0, amp_max=0.30 + 0.9 * bleb_index,
                                   t=0, life=int(34 + self.rs.randint(46))))
        for b in self.blebs:
            b["t"] += 1
            b["amp"] = b["amp_max"] * np.sin(np.pi * min(b["t"] / b["life"], 1.0))
        self.blebs = [b for b in self.blebs if b["t"] < b["life"]]
