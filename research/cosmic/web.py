"""Космическая паутина: ячейки Вороного галактик SDSS против случайных точек.

Объёмно-ограниченная выборка главной выборки галактик SDSS (z < Z_MAX, абсолютная величина ярче предела,
при котором выборка полна до Z_MAX): одинаковая плотность по всему объёму. Точные 3D-ячейки Вороного
(триангуляция Делоне), считаются только ячейки глубоко внутри обзора. Для каждой: число граней F,
объём, полный комбинаторный тип (код Вайнберга, research/tissue/topology.py). То же для случайных точек
в том же объёме и того же числа. Картинка: тонкий слой с галактиками и 2D-ячейками, цвет — плотность.
Выход: data/web_summary.json, docs/img/cosmic-web.png
"""
import json
import pathlib
import sys
from collections import Counter

import numpy as np
from scipy.integrate import quad
from scipy.spatial import ConvexHull, Delaunay, Voronoi

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "tissue"))
from topology import classify, face_graphs_from_delaunay  # noqa: E402

H0, OM = 70.0, 0.3
C = 299792.458
Z_MIN, Z_MAX, M_LIM_APP = 0.02, 0.10, 17.77
RA, DEC = (130.0, 235.0), (-2.0, 60.0)
CORE_RA, CORE_DEC, CORE_Z = (145.0, 220.0), (8.0, 50.0), (0.035, 0.088)


def dc(z):
    return C / H0 * quad(lambda x: 1 / np.sqrt(OM * (1 + x) ** 3 + 1 - OM), 0, z)[0]


zg = np.linspace(0, 0.13, 2000)
dg = np.array([dc(z) for z in zg])
D = lambda z: np.interp(z, zg, dg)  # noqa: E731  Мпк
Zof = lambda d: np.interp(d, dg, zg)  # noqa: E731


def xyz(ra, dec, d):
    ra, dec = np.radians(ra), np.radians(dec)
    return np.c_[d * np.cos(dec) * np.cos(ra), d * np.cos(dec) * np.sin(ra), d * np.sin(dec)]


def in_core(ra, dec, z):
    return (ra > CORE_RA[0]) & (ra < CORE_RA[1]) & (dec > CORE_DEC[0]) & (dec < CORE_DEC[1]) & (z > CORE_Z[0]) & (z < CORE_Z[1])


def analyse(P, core_idx, tag):
    vor = Voronoi(P)
    vols = []
    for i in core_idx:
        reg = vor.regions[vor.point_region[i]]
        if -1 in reg or not reg:
            vols.append(np.nan)
            continue
        vols.append(ConvexHull(vor.vertices[reg]).volume)
    cells = classify(face_graphs_from_delaunay(P, set(core_idx.tolist())))
    ok = [c for c in cells if c]
    F = np.array([c["F"] for c in ok])
    lv = np.log(np.array(vols)[np.isfinite(vols)])
    byF = {}
    rng = np.random.default_rng(0)
    for f in range(10, 21):
        c = [x for x in ok if x["F"] == f]
        if len(c) >= 100:
            sub = [c[k] for k in rng.choice(len(c), 100, replace=False)]
            byF[f] = {"n": len(c), "sym": float(np.mean([x["sym"] > 1 for x in c])), "distinct_100": len({x["code"] for x in sub})}
    out = {"cells": len(ok), "mean_F": float(F.mean()), "sd_F": float(F.std()), "sym": float(np.mean([c["sym"] > 1 for c in ok])),
           "sd_log_volume": float(lv.std()), "F_hist": sorted(Counter(F.tolist()).items()), "byF": byF}
    print(f"{tag:8s} cells {len(ok):6d}  F {out['mean_F']:.2f} ± {out['sd_F']:.2f}  sym {out['sym']:.3f}  sd(log V) {out['sd_log_volume']:.2f}")
    for f, b in byF.items():
        print(f"   F={f:2d}: sym {b['sym']:.2f}  distinct/100 {b['distinct_100']:3d}  n {b['n']}")
    return out


def picture(P, ra, dec, z, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    sel = np.abs(dec - 30) < 2.0
    d = D(z[sel])
    a = np.radians(ra[sel] - 182.5)
    Q = np.c_[d * np.sin(a), d * np.cos(a)]
    vor = Voronoi(Q)
    ang = np.arctan2(Q[:, 0], Q[:, 1])
    amin, amax = np.percentile(ang, 0.5) + 0.03, np.percentile(ang, 99.5) - 0.03
    polys, vals = [], []
    for i, reg in enumerate(vor.point_region):
        R = vor.regions[reg]
        if -1 in R or not R:
            continue
        V = vor.vertices[R]
        r = np.linalg.norm(V, axis=1)
        av = np.arctan2(V[:, 0], V[:, 1])
        if r.max() > D(Z_MAX) * 0.99 or r.min() < D(Z_MIN) * 1.3 or av.min() < amin or av.max() > amax:
            continue
        polys.append(V)
        vals.append(np.log(ConvexHull(V).volume))
    fig, ax = plt.subplots(figsize=(10, 5.6), dpi=140)
    pc = PolyCollection(polys, array=-np.array(vals), cmap="magma", edgecolors=(0.6, 0.6, 0.7, 0.35), linewidths=0.3)
    ax.add_collection(pc)
    inside = (np.linalg.norm(Q, axis=1) < D(Z_MAX) * 0.99) & (ang > amin) & (ang < amax)
    ax.scatter(Q[inside, 0], Q[inside, 1], s=0.6, c="white", alpha=0.7, linewidths=0)
    allv = np.vstack(polys)
    ax.set_xlim(allv[:, 0].min(), allv[:, 0].max()); ax.set_ylim(allv[:, 1].min(), allv[:, 1].max())
    ax.set_aspect("equal"); ax.axis("off")
    fig.patch.set_alpha(0)
    fig.savefig(out, bbox_inches="tight", transparent=True)
    print("picture:", out, "galaxies in slice", int(sel.sum()))


if __name__ == "__main__":
    raw = np.genfromtxt(HERE / "data" / "sdss_mgs.csv", delimiter=",", names=True, skip_header=1)
    ra, dec, z = raw["ra"], raw["dec"], raw["z"]
    m = raw["petroMag_r"] - raw["extinction_r"]
    DL = D(z) * (1 + z)
    M = m - 5 * np.log10(DL * 1e5)  # D_L в Мпк → модуль расстояния
    M_lim = M_LIM_APP - 5 * np.log10(D(Z_MAX) * (1 + Z_MAX) * 1e5)
    keep = (z > Z_MIN) & (z < Z_MAX) & (M < M_lim)
    ra, dec, z = ra[keep], dec[keep], z[keep]
    P = xyz(ra, dec, D(z))
    core = np.where(in_core(ra, dec, z))[0]
    print(f"volume-limited sample: {len(P)} galaxies (M_r < {M_lim:.2f}), core {len(core)}")
    res = {"n_sample": int(len(P)), "M_lim": float(M_lim), "galaxies": analyse(P, core, "galaxies")}
    # случайные точки в той же области (равномерно по объёму), то же число
    rng = np.random.default_rng(1)
    n = len(P)
    u = rng.random(n)
    d_r = np.cbrt(D(Z_MIN) ** 3 + u * (D(Z_MAX) ** 3 - D(Z_MIN) ** 3))
    ra_r = rng.uniform(*RA, n)
    sd = np.sin(np.radians(DEC))
    dec_r = np.degrees(np.arcsin(rng.uniform(*sd, n)))
    Pr = xyz(ra_r, dec_r, d_r)
    core_r = np.where(in_core(ra_r, dec_r, Zof(d_r)))[0]
    res["random"] = analyse(Pr, core_r, "random")
    (HERE / "data" / "web_summary.json").write_text(json.dumps(res, indent=1))
    picture(P, ra, dec, z, HERE.parents[1] / "docs" / "img" / "cosmic-web.png")
