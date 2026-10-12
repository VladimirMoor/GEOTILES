"""Изолированные точки остатка для прямых равносторонних шестиугольных призм.

По неравенству 720 = Σγ ≤ Σ T_r (см. RESULTS) открытыми могут быть только точки, где ≥ 3 соотношения
покрывают все вершины. Поиск: сетка в пространстве формы (внешние углы e1, e2, e3), для каждого соотношения
(Σa ≤ K) — клетки, где оно меняет знак; в каждой клетке перебираем тройки «своих» соотношений с покрытием
всех вершин и решаем систему 3×3 (Левенберг–Марквардт). Найденные точки проверяем полным фильтром и
отбрасываем лежащие на «хороших» тройках (сумма 360°: такие шестиугольники замощают плоскость).
Выход: data/hexprisms_points.json
"""
import itertools
import json
import pathlib
import sys

import numpy as np
from scipy.optimize import least_squares

from hexprisms import ROOT, TARGETS, TRIPLES, hexagon, on_good_triple

K = int(sys.argv[1]) if len(sys.argv) > 1 else 5
N = int(sys.argv[2]) if len(sys.argv) > 2 else 60


def relations(K):
    out = []
    for total in range(1, K + 1):
        for combo in itertools.combinations_with_replacement(range(6), total):
            a = np.bincount(combo, minlength=6)
            for T in TARGETS:
                if T < 180 * a.sum() - 1e-9:
                    out.append((a, T))
    return out


def passes(g, K=8, tol=1e-6):
    cov = set()
    for total in range(1, K + 1):
        for combo in itertools.combinations_with_replacement(range(6), total):
            a = np.bincount(combo, minlength=6)
            if any(abs(a @ g - T) < tol for T in TARGETS):
                cov |= set(combo)
        if len(cov) == 6:
            return True
    return len(cov) == 6


def main():
    ax = np.linspace(0.5, 179.5, N + 1)
    G = np.full((N + 1,) * 3 + (6,), np.nan)
    for i, a in enumerate(ax):
        for j, b in enumerate(ax):
            for k, c in enumerate(ax):
                g = hexagon([a, b, c])
                if g is not None:
                    G[i, j, k] = g
    valid = ~np.isnan(G[..., 0])
    cellok = np.ones((N, N, N), bool)
    for d in itertools.product((0, 1), repeat=3):
        cellok &= valid[d[0]:N + d[0], d[1]:N + d[1], d[2]:N + d[2]]
    print("valid cells:", int(cellok.sum()), flush=True)
    rels = relations(K)
    A = np.array([a for a, _ in rels], float)
    Tv = np.array([T for _, T in rels])
    supp = [frozenset(np.nonzero(a)[0].tolist()) for a, _ in rels]
    cell_rels = {}
    Gf = np.nan_to_num(G, nan=0.0)
    for r in range(len(rels)):
        V = Gf @ A[r] - Tv[r]
        lo = np.full((N, N, N), np.inf); hi = np.full((N, N, N), -np.inf)
        for d in itertools.product((0, 1), repeat=3):
            s = V[d[0]:N + d[0], d[1]:N + d[1], d[2]:N + d[2]]
            lo = np.minimum(lo, s); hi = np.maximum(hi, s)
        m = (lo <= 0) & (hi >= 0) & cellok
        for c in map(tuple, np.argwhere(m)):
            cell_rels.setdefault(c, []).append(r)
    print("cells with ≥3 relations:", sum(len(v) >= 3 for v in cell_rels.values()), flush=True)
    found = []
    tried = set()
    for c, rs in cell_rels.items():
        if len(rs) < 3:
            continue
        for t in itertools.combinations(rs, 3):
            if len(supp[t[0]] | supp[t[1]] | supp[t[2]]) < 6:
                continue
            if t in tried:
                continue
            tried.add(t)
            M = np.vstack([np.ones(6)] + [A[r] for r in t])
            if np.linalg.matrix_rank(M) < 4:
                continue
            e0 = np.array([ax[c[0]], ax[c[1]], ax[c[2]]]) + (ax[1] - ax[0]) / 2

            def F(e, t=t):
                g = hexagon(e)
                if g is None:
                    return np.full(3, 1e3)
                return np.array([A[r] @ g - Tv[r] for r in t])
            res = least_squares(F, e0, method="lm", xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=200)
            if np.max(np.abs(F(res.x))) > 1e-9:
                continue
            g = hexagon(res.x)
            if g is None or np.any(g < 0.5) or np.any(g > 179.5) or on_good_triple(g):
                continue
            if not passes(g):
                continue
            # канонизация по D6
            keys = []
            for sh in range(6):
                for rev in (1, -1):
                    keys.append(tuple(np.round(np.roll(g[::rev], sh), 5)))
            key = min(keys)
            if any(np.allclose(key, f["key"], atol=1e-4) for f in found):
                continue
            found.append({"key": key, "angles": g.round(6).tolist(), "relations": [[A[r].astype(int).tolist(), Tv[r]] for r in t]})
    print(f"isolated open points (K={K}, grid {N}):", len(found))
    for f in found:
        print("  ", f["angles"], f["relations"])
    for f in found:
        f["key"] = list(f["key"])
    (ROOT / "data" / f"hexprisms_points_K{K}.json").write_text(json.dumps(found, indent=1))


if __name__ == "__main__":
    main()
