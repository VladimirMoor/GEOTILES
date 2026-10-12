"""Независимая проверка найденных плиток (data/nondv/*.npy).

1. Покрытие: случайные точки в шаре радиуса 3·R(P) — каждая должна лежать ровно в одной копии γP
   (перебор всех γ = операция × сдвиг решётки, копии которых достают до шара).
2. Соседи по граням: для каждой грани — точки чуть снаружи (центр и около вершин); в какой копии они лежат.
   Если у грани больше одного соседа — разбиение не «грань к грани».
3. Дирихле: если разбиение грань к грани, ищем x, для которого каждая грань — серединный перпендикуляр
   к x и γ'x, где γ' — настоящий сосед по этой грани (наименьшие квадраты, линейно по x).
"""
import itertools
import json
import pathlib
import sys

import numpy as np
from scipy.spatial import ConvexHull

from nondv import Search


def copies_near(S, center, radius):
    B, Binv = S.B, np.linalg.inv(S.B)
    ops = [(B @ R @ Binv, B @ t) for R, t in zip(S.G.R, S.G.t)]
    out = []
    for (R, t), n in itertools.product(ops, itertools.product(range(-3, 4), repeat=3)):
        tt = t + B @ np.array(n, float)
        if np.linalg.norm(R @ center + tt - center) < radius:
            out.append((R, tt))
    return out


def inside(H, y, tol=1e-9):
    return np.all(H[:, :3] @ y + H[:, 3] <= tol)


def verify(S, p, n_samples=3000, seed=0):
    H = S.halfspaces(p)
    c, _ = S.cheb(H)
    _, V = S.vol(H, c)
    hull = ConvexHull(V)
    cen = V.mean(axis=0)
    Rp = np.max(np.linalg.norm(V - cen, axis=1))
    G = copies_near(S, cen, 8 * Rp)
    Hs = [S.transform(H, R, t) for R, t in G]
    # 1. покрытие
    rng = np.random.default_rng(seed)
    bad = 0
    for _ in range(n_samples):
        d = rng.normal(size=3); d /= np.linalg.norm(d)
        y = cen + d * 2.5 * Rp * rng.random() ** (1 / 3)
        k = sum(inside(h, y, 1e-12) for h in Hs)
        bad += int(k != 1)
    # 2. грани и их соседи
    eq = np.unique(np.round(hull.equations, 8), axis=0)
    facets = []
    for e in eq:
        on = np.abs(V @ e[:3] + e[3]) < 1e-7
        P = V[on]
        if len(P) < 3:
            continue
        fc = P.mean(axis=0)
        probes = [fc] + [0.6 * fc + 0.4 * q for q in P]
        nbr = set()
        for q in probes:
            y = q + 1e-6 * e[:3]
            for k, h in enumerate(Hs):
                if inside(h, y, 1e-12):
                    nbr.add(k)
        facets.append((e, P, nbr))
    f2f = all(len(nb) == 1 for _, _, nb in facets)
    # 3. Дирихле по настоящим соседям
    dv = None
    if f2f:
        rows, rhs = [], []
        for e, P, nb in facets:
            R, t = G[next(iter(nb))]
            for y in P:
                rows.append(2 * (R.T @ y - y) - 2 * R.T @ t)
                rhs.append(t @ t - 2 * y @ t)
        A, b = np.array(rows), np.array(rhs)
        x, *_ = np.linalg.lstsq(A, b, rcond=None)
        dv = float(np.sqrt(np.mean((A @ x - b) ** 2)))
    return {"facets": len(facets), "coverage_bad": int(bad), "samples": n_samples, "face_to_face": f2f,
            "dv_residual_true": dv, "neighbours_per_facet": sorted(len(nb) for _, _, nb in facets)}


if __name__ == "__main__":
    for f in sorted(pathlib.Path("data/nondv").glob(sys.argv[1] if len(sys.argv) > 1 else "*.npy")):
        hm_tag, sig, seed = f.stem.split("_")
        hm = {"P213": "P 21 3", "I213": "I 21 3", "P4132": "P 41 3 2", "P4332": "P 43 3 2", "I4132": "I 41 3 2"}.get(hm_tag, hm_tag)
        S = Search(hm, seed=int(seed))
        r = verify(S, np.load(f))
        print(f.name, json.dumps(r), flush=True)
