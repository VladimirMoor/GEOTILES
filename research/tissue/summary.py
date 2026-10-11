"""Сводка по тканям: распределение числа граней, сравнение с Вороным центров и с Вороным случайных точек.

Базовая линия — разбиение Вороного равномерно случайных точек (Пуассона) в кубе, внутренние ячейки,
с тем же относительным порогом площади грани, что и для клеток. Результат: data/summary.json.
"""
import json
import pathlib
from collections import Counter

import numpy as np
from scipy.spatial import Voronoi

HERE = pathlib.Path(__file__).resolve().parent
FRACS = ("0.0", "0.01", "0.02", "0.05")


def poly_area(P):
    c = P.mean(axis=0)
    return 0.5 * np.linalg.norm(sum(np.cross(P[i] - c, P[(i + 1) % len(P)] - c) for i in range(len(P))))


def poisson(n=6000, seed=1):
    rng = np.random.default_rng(seed)
    X = rng.random((n, 3))
    vor = Voronoi(X)
    area = {}
    for (i, j), rv in zip(vor.ridge_points, vor.ridge_vertices):
        if -1 not in rv:
            area[(i, j)] = poly_area(vor.vertices[rv])
    surf = np.zeros(n)
    for (i, j), a in area.items():
        surf[i] += a; surf[j] += a
    inner = np.where(np.all((X > 0.15) & (X < 0.85), axis=1))[0]
    out = {}
    for fr in FRACS:
        f = float(fr); nb = np.zeros(n, int)
        for (i, j), a in area.items():
            nb[i] += a >= f * surf[i]; nb[j] += a >= f * surf[j]
        out[fr] = nb[inner].tolist()
    return out


def kind(name):
    return "ovules" if name.startswith("N_") else "lrp"


res = {}
pois = poisson()
for k in ("ovules", "lrp"):
    files = [json.loads(p.read_text()) for p in sorted((HERE / "data").glob("*.json"))
             if p.name not in ("summary.json",) and kind(p.name) == k]
    if not files:
        continue
    d = {"files": len(files), "cells": sum(f["cells"] for f in files), "inner": sum(f["inner"] for f in files), "by_frac": {}}
    for fr in FRACS:
        F = sum((f["by_frac"][fr]["faces"] for f in files), [])
        VF = sum((f["by_frac"][fr]["vor_faces"] for f in files), [])
        vol = sum((f["by_frac"][fr]["volume"] for f in files), [])
        ok = [f for f in files if f["inner"] > 0]
        jac = np.average([f["by_frac"][fr]["jaccard_mean"] for f in ok], weights=[f["inner"] for f in ok])
        lewis = {}
        for n, v in zip(F, vol):
            lewis.setdefault(n, []).append(v)
        mv = np.mean(vol)
        d["by_frac"][fr] = {"mean_faces": float(np.mean(F)), "sd_faces": float(np.std(F)), "mean_vor_faces": float(np.mean(VF)),
                            "mean_poisson": float(np.mean(pois[fr])), "jaccard": float(jac),
                            "hist": sorted(Counter(F).items()), "hist_vor": sorted(Counter(VF).items()),
                            "hist_poisson": sorted(Counter(pois[fr]).items()),
                            "lewis": {int(n): float(np.mean(v) / mv) for n, v in sorted(lewis.items()) if len(v) >= 20}}
    res[k] = d
    print(k, "files", d["files"], "cells", d["cells"], "inner", d["inner"])
    for fr in FRACS:
        b = d["by_frac"][fr]
        print(f"  frac {fr:4s}: faces {b['mean_faces']:.2f}±{b['sd_faces']:.2f}  voronoi(centres) {b['mean_vor_faces']:.2f}  "
              f"poisson {b['mean_poisson']:.2f}  jaccard {b['jaccard']:.3f}")
    print("  Lewis (frac 0.01, V/<V> by faces):", {n: round(v, 2) for n, v in d["by_frac"]["0.01"]["lewis"].items()})
(HERE / "data" / "summary.json").write_text(json.dumps(res))
