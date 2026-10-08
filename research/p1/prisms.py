"""Шаг 6. Прямые равносторонние призмы: какие основания проходят 3D-фильтр углов.

Двугранные углы прямой призмы над многоугольником P: углы γᵢ многоугольника (боковые рёбра) и 90°
(рёбра оснований). Фильтр: для каждого i есть целые aⱼ ≥ 0, aᵢ ≥ 1, m ≥ 0 с
    Σ aⱼ γⱼ + 90·m = 360°   (180° от грани соседа = 2·90°, отдельно не нужен),
то есть Σ aⱼ γⱼ ∈ {360, 270, 180, 90}.

Страты (для пятиугольников; Σγ = 540°):
- кривые, на которых ОДНО 0/1-соотношение вместе с дополнением покрывает все вершины:
  γᵢ + γⱼ = 180° (дополнение 360°) и γᵢ + γⱼ = 270° (дополнение 270°);
- точки — пересечения двух соотношений, покрывающих все вершины.
Кривая «180°» — ровно пятиугольники, замощающие плоскость (Hirschhorn–Hunt, кроме P7).
Всё остальное — «остаток»: прямые призмы, проходящие фильтр, но с основанием, не замощающим плоскость
(частные случаи вопроса Куперберга 1993).

Запуск: python prisms.py  → data/prisms_right.json
"""
import itertools
import json
import math
import pathlib

import numpy as np
from scipy.optimize import fsolve

ROOT = pathlib.Path(__file__).resolve().parent
TARGETS = (360.0, 270.0, 180.0, 90.0)


# ---------- равносторонний пятиугольник по углам при A и B ----------

def pentagon(alpha, beta):
    """A=(0,0), B=(1,0); углы при A и B в градусах. Возвращает вершины A..E или None (не выпуклый)."""
    a, b = math.radians(alpha), math.radians(beta)
    A, B = np.array([0.0, 0]), np.array([1.0, 0])
    E = np.array([math.cos(a), math.sin(a)])
    C = B + np.array([-math.cos(b), math.sin(b)])
    M = (C + E) / 2
    s = np.linalg.norm(C - E) / 2
    if s >= 1 or s < 1e-9:
        return None
    n = np.array([-(E - C)[1], (E - C)[0]]) / (2 * s)  # ⊥ CE
    if n @ (M - (A + B) / 2) < 0:
        n = -n
    D = M + math.sqrt(1 - s * s) * n
    P = np.array([A, B, C, D, E])
    return P


def angles(P):
    k = len(P)
    out = []
    for i in range(k):
        u, w = P[i - 1] - P[i], P[(i + 1) % k] - P[i]
        cross = (P[i] - P[i - 1])[0] * (P[(i + 1) % k] - P[i])[1] - (P[i] - P[i - 1])[1] * (P[(i + 1) % k] - P[i])[0]
        ang = math.degrees(math.acos(np.clip(u @ w / (np.linalg.norm(u) * np.linalg.norm(w)), -1, 1)))
        out.append(ang if cross > 0 else 360 - ang)
    return np.array(out)


def convex_ok(P, margin=0.5):
    g = angles(P)
    return bool(np.all(g < 180 - margin) and np.all(g > margin) and abs(g.sum() - 180 * (len(P) - 2)) < 1e-6)


def gamma(x):
    P = pentagon(*x)
    return None if P is None else angles(P)


# ---------- соотношения ----------

def relations(n, K):
    """Все векторы a ∈ Z≥0^n с 1 ≤ Σa ≤ K."""
    out = []
    for total in range(1, K + 1):
        for combo in itertools.combinations_with_replacement(range(n), total):
            out.append(np.bincount(combo, minlength=n))
    return out


def covered(g, tol=1e-6, K=12):
    """Множество вершин, входящих в какое-нибудь соотношение при углах g."""
    cov = set()
    for a in relations(len(g), K):
        s = a @ g
        if any(abs(s - T) < tol for T in TARGETS):
            cov |= set(np.nonzero(a)[0].tolist())
    return cov


def passes(g, K=12):
    return len(covered(g, K=K)) == len(g)


# ---------- страты ----------

def curve_strata():
    """0/1-соотношения γᵢ + γⱼ = T, T ∈ {180, 270}: покрывают все вершины вместе с дополнением."""
    out = []
    for i, j in itertools.combinations(range(5), 2):
        for T in (180.0, 270.0):
            out.append({"pair": [i, j], "T": T, "adjacent": (j - i) % 5 in (1, 4)})
    return out


def point_strata(K=7, grid=120):
    """Пересечения двух соотношений, после которых (с учётом 0/1-дополнений) покрыты все 5 вершин.
    Кандидатные пары: соотношения, меняющие знак в одной и той же клетке сетки (ищутся матричным
    произведением масок); затем точное решение fsolve и полная проверка фильтра."""
    al = np.linspace(2, 178, grid)
    be = np.linspace(2, 178, grid)
    Gr = np.full((grid, grid, 5), np.nan)
    for i, x in enumerate(al):
        for j, y in enumerate(be):
            P = pentagon(x, y)
            if P is not None and convex_ok(P, margin=0.0):
                Gr[i, j] = angles(P)
    rels = relations(5, K)
    labels, masks = [], []
    valid = ~np.isnan(Gr[..., 0])
    cellok = valid[:-1, :-1] & valid[1:, :-1] & valid[:-1, 1:] & valid[1:, 1:]
    for a in rels:
        V = np.nan_to_num(Gr @ a, nan=1e9)
        for T in TARGETS:
            W = V - T
            c = [W[:-1, :-1], W[1:, :-1], W[:-1, 1:], W[1:, 1:]]
            lo = np.minimum.reduce(c)
            hi = np.maximum.reduce(c)
            m = (lo <= 0) & (hi >= 0) & cellok
            if m.any():
                labels.append((a, T))
                masks.append(m.ravel())
    M = np.array(masks, dtype=np.float32)
    print(f"  {len(labels)} relation/target pairs change sign somewhere; intersecting masks...")
    C = M @ M.T
    supports = [set(np.nonzero(a)[0].tolist()) for a, _ in labels]

    def cover(k):
        a, T = labels[k]
        s = set(supports[k])
        if np.all(a <= 1) and (540 - T) in TARGETS:
            s |= set(range(5)) - s  # 0/1-дополнение — тоже соотношение
        return s

    covers = [cover(k) for k in range(len(labels))]
    cand = [(p, q) for p, q in zip(*np.nonzero(np.triu(C, 1) > 0)) if len(covers[p] | covers[q]) >= 4]
    print(f"  {len(cand)} candidate pairs")
    cells = np.argwhere(cellok)
    found = []
    pairs_done = 0
    for p, q in cand:
        (a, Ta), (b, Tb) = labels[p], labels[q]
        common = np.nonzero(M[p].astype(bool) & M[q].astype(bool))[0]
        pairs_done += 1
        for cidx in common[:2]:
            i, j = divmod(int(cidx), grid - 1)
            start = np.array([al[i] + (al[1] - al[0]) / 2, be[j] + (be[1] - be[0]) / 2])
            def F(x):
                g = gamma(x)
                if g is None:
                    return [1e3, 1e3]
                return [a @ g - Ta, b @ g - Tb]
            x, info, ier, _ = fsolve(F, start, full_output=True, xtol=1e-13)
            if ier != 1 or np.max(np.abs(F(x))) > 1e-9:
                continue
            P = pentagon(*x)
            if P is None or not convex_ok(P, margin=0.01):
                continue
            g = angles(P)
            if not passes(g):
                continue
            key = tuple(np.round(np.sort(g), 5))
            if any(np.allclose(key, f["key"], atol=1e-4) for f in found):
                continue
            on_curve = [(u, w, round(g[u] + g[w], 6)) for u, w in itertools.combinations(range(5), 2)
                        if abs(g[u] + g[w] - 180) < 1e-6 or abs(g[u] + g[w] - 270) < 1e-6]
            found.append({"key": key, "alpha_beta": x.tolist(), "angles": g.round(6).tolist(),
                          "onCurve": on_curve})
    return found, pairs_done


def main():
    # кривые: проверим, что они непусты среди выпуклых равносторонних пятиугольников
    curves = curve_strata()
    probe = []
    for al in np.linspace(5, 175, 400):
        for be in np.linspace(5, 175, 400):
            P = pentagon(al, be)
            if P is not None and convex_ok(P):
                probe.append(angles(P))
    probe = np.array(probe)
    for c in curves:
        i, j = c["pair"]
        v = probe[:, i] + probe[:, j] - c["T"]
        c["nonEmpty"] = bool(v.min() < 0 < v.max())
    print("Curve strata (each covers all vertices on its own):")
    for c in curves:
        print(f"  γ{'ABCDE'[c['pair'][0]]} + γ{'ABCDE'[c['pair'][1]]} = {c['T']:.0f}°  adjacent={c['adjacent']}  realizable={c['nonEmpty']}")
    pts, n = point_strata()
    off = [p for p in pts if not p["onCurve"]]
    print(f"\nPoint strata: {len(pts)} filter-passing pentagons from {n} relation pairs; "
          f"{len(off)} lie on neither the 180° nor the 270° curve:")
    for p in off:
        print("  angles", p["angles"])
    (ROOT / "data" / "prisms_right.json").write_text(json.dumps({"curves": curves, "points": pts}, default=float, indent=1))


if __name__ == "__main__":
    main()
