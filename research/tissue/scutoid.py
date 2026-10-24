"""Scutoids in the outer cell layer (L1, epidermis) of Arabidopsis ovules (PlantSeg ground truth).

A scutoid (Gómez-Gálvez et al., Nat Commun 2018) is a cell of a curved epithelium whose neighbours on the apical
surface differ from its neighbours on the basal surface. Here the epithelium is the L1 layer: cells touching the
background (outside of the ovule). For every voxel of the layer the relative depth is t = d_out / (d_out + d_in),
d_out = distance to the background, d_in = distance to the inner cells (both in µm). A wall between two L1 cells
counts on the apical side if its area at t < BAND is at least A_MIN µm², on the basal side if its area at
t > 1 − BAND is at least A_MIN. A cell is a scutoid when its apical and basal sets of L1 neighbours differ.
Cells touching the volume border or ignore regions, or lacking an apical or basal surface, are not counted.

Also geometric null models (column_model): straight columns along the surface normals, or along the axis through
the real centres of each cell's outer and inner faces; Voronoi of the axis points at each depth.

Usage: python scutoid.py file.h5 [vz vy vx] → data/scutoid/<name>.json
"""
import json
import pathlib
import sys

import h5py
import numpy as np
from scipy import ndimage as ndi

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "data" / "scutoid"
BAND, A_MIN_LIST = 0.2, (0.5, 1.0, 2.0)
IGN = 65000


def wall_pairs(L, t, isL1, vox):
    """For every pair of touching L1 cells: wall area at t < BAND and at t > 1 − BAND (µm²)."""
    area = {0: vox[1] * vox[2], 1: vox[0] * vox[2], 2: vox[0] * vox[1]}
    acc = {}
    for ax in range(3):
        sl0 = [slice(None)] * 3
        sl1 = [slice(None)] * 3
        sl0[ax], sl1[ax] = slice(None, -1), slice(1, None)
        a, b = L[tuple(sl0)], L[tuple(sl1)]
        m = (a != b) & isL1[a] & isL1[b]
        aa, bb = a[m].astype(np.int64), b[m].astype(np.int64)
        tt = (t[tuple(sl0)][m] + t[tuple(sl1)][m]) / 2
        lo, hi = np.minimum(aa, bb), np.maximum(aa, bb)
        key = lo * 65536 + hi
        u, inv = np.unique(key, return_inverse=True)
        ap = np.bincount(inv, weights=(tt < BAND) * area[ax], minlength=len(u))
        ba = np.bincount(inv, weights=(tt > 1 - BAND) * area[ax], minlength=len(u))
        tot = np.bincount(inv, weights=np.full(len(tt), area[ax]), minlength=len(u))
        for k, x, y, z in zip(u.tolist(), ap, ba, tot):
            p = acc.setdefault(k, [0.0, 0.0, 0.0])
            p[0] += x; p[1] += y; p[2] += z
    return acc


def fields(path, vox):
    """Label volume (ignore regions = IGN), relative depth t in the L1 layer, L1 flags, distance to the outside."""
    h = h5py.File(path, "r")
    L = h["label"][:].astype(np.int32)
    if "label_with_ignore" in h:
        L[h["label_with_ignore"][:] == -1] = IGN
    M = IGN + 1
    bg = L == 0
    # L1: cells with a substantial contact with the background
    near_bg = ndi.binary_dilation(bg) & ~bg
    cnt_bg = np.bincount(L[near_bg].ravel(), minlength=M)
    vol = np.bincount(L.ravel(), minlength=M)
    isL1 = np.zeros(M, bool)
    isL1[(cnt_bg >= 20) & (np.arange(M) > 0) & (np.arange(M) != IGN)] = True
    # the epithelium proper: L1 cells that span the layer, touching the outside and also cells that never touch it
    # (in ovules several L1 cells can be stacked, and integuments are sheets with background on both sides)
    deep_cells = (L > 0) & (L != IGN) & ~isL1[L]
    near_deep = ndi.binary_dilation(deep_cells) & isL1[L]
    span = isL1 & (np.bincount(L[near_deep].ravel(), minlength=M) >= 20)
    inner = (L > 0) & (L != IGN) & ~span[L]  # everything below the spanning layer, stacked L1 cells included
    d_out = ndi.distance_transform_edt(~bg, sampling=vox).astype(np.float32)
    d_in = ndi.distance_transform_edt(~inner, sampling=vox).astype(np.float32)
    t = d_out / np.maximum(d_out + d_in, 1e-6)
    del d_in
    return L, t, span, inner, d_out, cnt_bg, vol


def detect(L, t, isL1, good, vox, extra):
    """Apico-basal neighbour changes among layer cells of label volume L (see the module docstring)."""
    acc = wall_pairs(L, t, isL1, vox)
    anyc = {}  # all layer neighbours of a cell (any depth)
    for k in acc:
        a, b = k // 65536, k % 65536
        anyc.setdefault(a, set()).add(b)
        anyc.setdefault(b, set()).add(a)
    out = {}
    for amin in A_MIN_LIST:
        # strict: a pair is apical-only if its apical wall area >= amin and it has no wall at all in the basal band
        # (and vice versa); loose: presence in each band judged by area >= amin alone
        ap_only, ba_only, loose = set(), set(), set()
        for k, (x, y, _) in acc.items():
            if x >= amin and y == 0:
                ap_only.add(k)
            if y >= amin and x == 0:
                ba_only.add(k)
            if (x >= amin) != (y >= amin):
                loose.add(k)

        def quartet(k, other):
            """k = (a, b) touching on one side only; is there (c, d) touching on the other side only, with c and d
            both neighbours of a and b? That is an apico-basal T1 exchange."""
            a, b = k // 65536, k % 65536
            common = anyc.get(a, set()) & anyc.get(b, set())
            return any((min(c, d) * 65536 + max(c, d)) in other for c in common for d in common if c < d)

        t1 = {k for k in ap_only if quartet(k, ba_only)} | {k for k in ba_only if quartet(k, ap_only)}
        quads = [[k // 65536, k % 65536, "apical" if k in ap_only else "basal"] for k in sorted(t1)]
        cells = []
        for i in np.flatnonzero(good):
            i = int(i)
            if i not in anyc:
                continue
            mine = lambda S: sum(1 for k in S if i in (k // 65536, k % 65536))  # noqa: E731
            cells.append({"id": i, "loose": mine(loose), "strict": mine(ap_only | ba_only), "t1": mine(t1), **extra.get(i, {})})
        out[str(amin)] = {"cells": cells, "n": len(cells), "t1_pairs": quads,
                          **{f"share_{m}": float(np.mean([c[m] > 0 for c in cells])) if cells else 0.0
                             for m in ("loose", "strict", "t1")}}
    return out


def column_model(L, t, isL1, d_out, vox, levels=20, real_axes=False):
    """Geometric null model: every layer cell is a straight column along the outward normal at the centre of its
    outer face, as thick as the layer there; at each depth a layer voxel goes to the nearest column axis
    (a Voronoi tessellation of the axis points at that depth, as in the Voronoi tubes of Gómez-Gálvez et al.)."""
    from scipy.spatial import cKDTree
    layer = isL1[L]
    Z = np.argwhere(layer)
    lab, tt = L[layer], t[layer]
    X = Z * np.array(vox)
    g = np.stack(np.gradient(d_out, *vox), axis=-1)[layer]
    ids, axes = [], []
    for i in np.unique(lab):
        m = lab == i
        ap, ba = m & (tt < BAND), m & (tt > 1 - BAND)
        if ap.sum() < 5 or ba.sum() < 5:
            continue
        n = g[ap].mean(0)
        n /= np.linalg.norm(n) + 1e-12
        a0, t0 = X[ap].mean(0), tt[ap].mean()
        h = (X[ba].mean(0) - a0) @ n / max(tt[ba].mean() - t0, 1e-3)  # depth per unit t along the normal
        if real_axes:  # axis through the real centres of the outer and inner faces (keeps the real tilt)
            d = X[ba].mean(0) - a0
            h = np.linalg.norm(d) / max(tt[ba].mean() - t0, 1e-3)
            n = d / (np.linalg.norm(d) + 1e-12)
        ids.append(i)
        axes.append((a0, t0, n, h))
    ids = np.array(ids)
    Lm = L.copy()
    newlab = np.zeros(len(lab), L.dtype)
    edges = np.linspace(0, 1, levels + 1)
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = (tt >= lo) & (tt <= hi)
        if not sel.any():
            continue
        tc = (lo + hi) / 2
        P = np.array([a0 + (tc - t0) * h * n for a0, t0, n, h in axes])
        _, j = cKDTree(P).query(X[sel])
        newlab[sel] = ids[j]
    Lm[tuple(Z.T)] = newlab
    return Lm


def analyse(path, vox):
    M = IGN + 1
    L, t, isL1, inner, d_out, cnt_bg, vol = fields(path, vox)
    # cells to judge: layer cells with a basal surface, not touching the border or ignore regions
    near_in = ndi.binary_dilation(inner) & isL1[L]
    has_basal = np.bincount(L[near_in].ravel(), minlength=M) >= 20
    bad = np.zeros(M, bool)
    for ax in range(3):
        for idx in (0, L.shape[ax] - 1):
            bad[np.unique(np.take(L, idx, axis=ax))] = True
    near_ign = ndi.binary_dilation(L == IGN)
    bad[np.unique(L[near_ign])] = True
    good = isL1 & has_basal & ~bad
    thick = ndi.maximum(d_out, L, index=np.arange(M))
    ap_area = cnt_bg * vox[1] * vox[2]
    # basal area: layer voxels of the cell next to the inner tissue
    near_in_count = np.bincount(L[near_in].ravel(), minlength=M)
    extra = {int(i): {"volume_um3": float(vol[i] * np.prod(vox)), "thickness_um": float(thick[i]),
                      "apical_area_um2": float(ap_area[i]), "basal_area_um2": float(near_in_count[i] * vox[1] * vox[2])}
             for i in np.flatnonzero(good)}
    res = {"file": pathlib.Path(path).name, "n_L1": int(isL1.sum()), "n_judged": int(good.sum()),
           "by_threshold": detect(L, t, isL1, good, vox, extra)}
    Lm = column_model(L, t, isL1, d_out, vox)
    res["model"] = detect(Lm, t, isL1, good, vox, {})
    Lm = column_model(L, t, isL1, d_out, vox, real_axes=True)
    res["model_real_axes"] = detect(Lm, t, isL1, good, vox, {})
    return res


if __name__ == "__main__":
    path = sys.argv[1]
    vox = tuple(float(v) for v in sys.argv[2:5]) if len(sys.argv) > 4 else (0.235, 0.150, 0.150)
    OUT.mkdir(parents=True, exist_ok=True)
    r = analyse(path, vox)
    (OUT / (pathlib.Path(path).stem + ".json")).write_text(json.dumps(r))
    print(r["file"], "L1", r["n_L1"], "judged", r["n_judged"],
          {k: (v["n"], round(v["share_loose"], 3), round(v["share_strict"], 3), round(v["share_t1"], 3)) for k, v in r["by_threshold"].items()},
          "model", {k: (v["n"], round(v["share_strict"], 3), round(v["share_t1"], 3)) for k, v in r["model"].items()},
          "real axes", {k: (round(v["share_strict"], 3), round(v["share_t1"], 3)) for k, v in r["model_real_axes"].items()}, flush=True)
