"""
Extracellular ice as a field of grains surrounding the cell on every side.

The cell does not sit in a slab between two ice sheets. It sits in an irregular
unfrozen pocket bounded by the growth fronts of many grains, and that pocket
closes in from all directions as the ice fraction rises. This module builds the
pocket geometrically from the grain fronts, so:

  * the boundary in any direction is set by whichever grain front is nearest,
  * the pocket is irregular because grains sit at different distances,
  * coarse ice gives a wide, angular pocket and fine ice a narrow, polygonal one,
  * the mean pocket radius tracks the unfrozen channel width the model computes.

The membrane soft body is then confined by pocket_radius(theta) in every
direction, which is what produces realistic all-round compression.
"""
from __future__ import annotations
import numpy as np


class IceField:
    def __init__(self, n_grains: int = 22, seed: int = 7):
        rs = np.random.RandomState(seed)
        self.n = n_grains
        # angles spread around the cell with jitter, so no direction is special
        base = np.linspace(0, 2 * np.pi, n_grains, endpoint=False)
        self.theta = base + rs.uniform(-0.5, 0.5, n_grains) * (2 * np.pi / n_grains)
        self.size_jit = rs.uniform(0.65, 1.45, n_grains)     # grain size variation
        self.dist_jit = rs.uniform(0.72, 1.34, n_grains)     # how close each front sits
        self.rot = rs.uniform(0, np.pi, n_grains)            # facet orientation
        self.shade = rs.uniform(0.0, 1.0, n_grains)
        self.centres = np.zeros((n_grains, 2))
        self.radii = np.zeros(n_grains)
        self.pocket_mean = 1e4
        self.long = 1e4
        self.aniso = 1.0
        self.active = False

    def update(self, f_ice: float, grain_um: float, px_per_um: float,
               chan_um: float, r_cell_px: float):
        """Place the grain fronts so the pocket has the model's channel width in
        its narrow direction, while keeping enough enclosed area for the cell.

        A cell is not compressible. When the intergranular channel narrows below
        the cell diameter the cell does not shrink to fit -- the pocket it
        occupies is elongated, and the cell deforms into it. So the pocket is
        made anisotropic, with the anisotropy growing as the channel narrows,
        and its area is held at or above the cell's own. A pocket that simply
        closed in isotropically would crush the cell, which is not what happens.
        """
        self.active = f_ice > 0.02
        if not self.active:
            self.pocket_mean = 1e4
            self.long = 1e4
            self.aniso = 1.0
            return
        grain_px = max(8.0, grain_um * px_per_um)
        half = chan_um * px_per_um * 0.5
        # required area to hold the cell with a little clearance
        need = np.pi * (r_cell_px ** 2) * 1.25
        narrow = float(np.clip(half, r_cell_px * 0.16, r_cell_px * 4.5))
        # elongate along one axis until the ellipse of half-width `narrow` has
        # at least the area the cell needs
        self.aniso = float(np.clip(need / (np.pi * narrow * narrow), 1.0, 9.0))
        self.pocket_mean = narrow
        self.long = narrow * self.aniso
        self.radii = 0.5 * grain_px * self.size_jit
        # grain fronts sit on that ellipse, jittered, so the pocket is irregular
        e = np.hypot(np.cos(self.theta) * self.long, np.sin(self.theta) * narrow)
        d = e * self.dist_jit + self.radii
        self.centres = np.stack([np.cos(self.theta) * d, np.sin(self.theta) * d], axis=1)

    def pocket_radius(self, theta):
        """Distance from the cell centre to the nearest grain front along each ray."""
        theta = np.atleast_1d(np.asarray(theta, dtype=float))
        if not self.active:
            return np.full(theta.shape, 1e4)
        dx = np.cos(theta)[:, None]; dy = np.sin(theta)[:, None]   # (T,1)
        cx = self.centres[None, :, 0]; cy = self.centres[None, :, 1]  # (1,G)
        r = self.radii[None, :]
        # ray from origin: |s*d - c|^2 = r^2  ->  s^2 - 2 s (d.c) + |c|^2 - r^2 = 0
        b = dx * cx + dy * cy
        c2 = cx * cx + cy * cy - r * r
        disc = b * b - c2
        hit = disc > 0
        s = np.where(hit, b - np.sqrt(np.maximum(disc, 0.0)), np.inf)
        s = np.where(s > 0, s, np.inf)
        out = np.min(s, axis=1)
        # cap by the pocket envelope itself, so directions with no grain in the
        # way do not spike out to infinity
        env = 1.0 / np.sqrt((np.cos(theta) / self.long) ** 2
                            + (np.sin(theta) / self.pocket_mean) ** 2)
        out = np.minimum(np.where(np.isfinite(out), out, 1e4), env * 1.18)
        return out

    def polygon(self, n: int = 220):
        th = np.linspace(0, 2 * np.pi, n, endpoint=False)
        r = np.minimum(self.pocket_radius(th), 4000.0)
        return np.stack([np.cos(th) * r, np.sin(th) * r], axis=1), th, r

    def grain_boundaries(self, extent: float):
        """Grain boundaries as the Voronoi tessellation of the grain centres.

        Drawing each grain as a disc gives overlapping blobs; a polycrystal is
        a tessellation, so the boundaries are where two fronts met. Ridges are
        clipped to the drawing extent.
        """
        try:
            from scipy.spatial import Voronoi
        except Exception:
            return []
        # ring of far points so the outer cells are bounded
        far = np.stack([np.cos(np.linspace(0, 2*np.pi, 24, endpoint=False)) * extent * 2.2,
                        np.sin(np.linspace(0, 2*np.pi, 24, endpoint=False)) * extent * 2.2], axis=1)
        pts = np.vstack([self.centres, far])
        try:
            vor = Voronoi(pts)
        except Exception:
            return []
        segs = []
        for (a, b), (p, q) in zip(vor.ridge_points, vor.ridge_vertices):
            if p < 0 or q < 0: continue
            v1, v2 = vor.vertices[p], vor.vertices[q]
            if max(abs(v1).max(), abs(v2).max()) > extent * 2.0: continue
            segs.append((v1, v2))
        return segs

    def grain_cells(self, extent: float):
        """Filled Voronoi cells, so the ice reads as a polycrystal."""
        try:
            from scipy.spatial import Voronoi
        except Exception:
            return []
        ang = np.linspace(0, 2*np.pi, 24, endpoint=False)
        far = np.stack([np.cos(ang) * extent * 2.2, np.sin(ang) * extent * 2.2], axis=1)
        pts = np.vstack([self.centres, far])
        try:
            vor = Voronoi(pts)
        except Exception:
            return []
        out = []
        for i in range(self.n):
            reg = vor.regions[vor.point_region[i]]
            if not reg or -1 in reg: continue
            poly = vor.vertices[reg]
            if np.abs(poly).max() > extent * 2.0: continue
            out.append((poly, self.shade[i]))
        return out
