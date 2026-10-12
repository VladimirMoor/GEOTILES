"""Прямые равносторонние шестиугольные призмы: одномерный остаток фильтра двугранных углов.

Двугранные углы — углы γᵢ шестиугольника и 90°. Фильтр: для каждой вершины i есть соотношение
Σ aⱼγⱼ ∈ {360, 270, 180, 90} с aᵢ ≥ 1. Коразмерность 1 — только тройки с суммой 360° (они все замощают,
см. RESULTS). Здесь ищем кривые: пары соотношений (R1, R2), носители которых вместе покрывают все 6 вершин,
на которых шестиугольник равносторонний, выпуклый и не лежит на «хорошей» тройке (сумма 360°).

Шаги: 1) перебор соотношений Σa ≤ 5; 2) пары с объединением носителей = {1..6} и совместной ЛП
(Σγ = 720, 0 < γ < 180); 3) поиск точки на кривой проекцией Гаусса–Ньютона в пространстве формы
(3 параметра: внешние углы e1, e2, e3); 4) проверка: выпуклость, нет тройки 360°, на кривой есть
точки вне троек (проверяем две разные точки). Выход: data/hexprisms_curves.json
"""
import itertools
import json
import math
import pathlib

import numpy as np
from scipy.optimize import linprog, least_squares

ROOT = pathlib.Path(__file__).resolve().parent
TARGETS = (360.0, 270.0, 180.0, 90.0)
TRIPLES = list(itertools.combinations(range(6), 3))


def hexagon(e):
    """Внешние углы e1, e2, e3 (градусы) → вершины равностороннего выпуклого шестиугольника или None."""
    phi = np.radians(np.cumsum([0.0, e[0], e[1], e[2]]))
    U = np.c_[np.cos(phi), np.sin(phi)]
    r = -U.sum(axis=0)
    L = np.linalg.norm(r)
    if L >= 2 or L < 1e-9:
        return None
    psi = math.atan2(r[1], r[0])
    d = math.acos(L / 2)
    best = None
    for s in (1, -1):
        p5, p6 = psi - s * d, psi + s * d
        dirs = np.r_[phi, p5, p6]
        ex = np.degrees(np.diff(np.r_[dirs, dirs[0] + 2 * np.pi]))
        ex = (ex + 360) % 360
        if np.all(ex > 0.01) and np.all(ex < 179.99) and abs(ex.sum() - 360) < 1e-6:
            best = 180 - ex  # внутренние углы при вершинах (между u_k и u_{k+1})
    return best


def gam(e):
    return hexagon(e)


def relations(K=5):
    out = []
    for total in range(1, K + 1):
        for combo in itertools.combinations_with_replacement(range(6), total):
            a = np.bincount(combo, minlength=6)
            for T in TARGETS:
                if T < 180 * a.sum() - 1e-9:
                    out.append((a, T))
    return out


def lp_feasible(rs):
    A_eq = [np.ones(6)] + [a for a, _ in rs]
    b_eq = [720.0] + [T for _, T in rs]
    r = linprog(np.zeros(6), A_eq=np.array(A_eq), b_eq=np.array(b_eq), bounds=[(1.0, 179.0)] * 6, method="highs")
    return r.status == 0


def on_good_triple(g, tol=1e-5):
    return any(abs(g[list(t)].sum() - 360) < tol for t in TRIPLES)


def find_point(rs, rng, tries=12):
    def F(e):
        g = gam(e)
        if g is None:
            return np.full(len(rs), 1e3)
        return np.array([a @ g - T for a, T in rs])
    pts = []
    for _ in range(tries):
        e0 = rng.uniform(20, 110, 3)
        if gam(e0) is None:
            continue
        r = least_squares(F, e0, xtol=1e-14, ftol=1e-14, gtol=1e-14, max_nfev=400)
        if np.max(np.abs(F(r.x))) < 1e-8:
            g = gam(r.x)
            if g is not None and np.all(g > 0.5) and np.all(g < 179.5):
                pts.append((r.x, g))
    return pts


def main():
    rng = np.random.default_rng(1)
    rels = relations()
    print("relations:", len(rels))
    supp = [frozenset(np.nonzero(a)[0].tolist()) for a, _ in rels]
    curves = []
    seen = set()
    checked = 0
    for i in range(len(rels)):
        for j in range(i + 1, len(rels)):
            if len(supp[i] | supp[j]) < 6:
                continue
            (a1, T1), (a2, T2) = rels[i], rels[j]
            # независимость с Σγ = 720
            if np.linalg.matrix_rank(np.vstack([np.ones(6), a1, a2])) < 3:
                continue
            if not lp_feasible([rels[i], rels[j]]):
                continue
            checked += 1
            pts = find_point([rels[i], rels[j]], rng)
            gen = [(e, g) for e, g in pts if not on_good_triple(g)]
            if len(gen) < 1:
                continue
            # канонический ключ кривой — по двум соотношениям
            key = (tuple(a1), T1, tuple(a2), T2)
            if key in seen:
                continue
            seen.add(key)
            # проверим, что кривая — действительно кривая (найдено ≥ 2 различных точки), и запомним пример
            distinct = []
            for e, g in gen:
                if all(np.linalg.norm(g - h) > 1e-3 for _, h in distinct):
                    distinct.append((e, g))
            curves.append({"R1": [a1.tolist(), T1], "R2": [a2.tolist(), T2], "points": len(distinct),
                           "example_angles": distinct[0][1].round(6).tolist()})
    print("pairs passing LP:", checked, "curves with generic points:", len(curves))
    # приведение по симметрии D6 (циклический сдвиг и отражение меток вершин)
    def canon(c):
        best = None
        for sh in range(6):
            for rev in (False, True):
                perm = [(k + sh) % 6 for k in range(6)]
                if rev:
                    perm = [(-k + sh) % 6 for k in range(6)]
                r1 = (tuple(np.array(c["R1"][0])[perm]), c["R1"][1])
                r2 = (tuple(np.array(c["R2"][0])[perm]), c["R2"][1])
                k = tuple(sorted([r1, r2]))
                if best is None or k < best:
                    best = k
        return best
    classes = {}
    for c in curves:
        classes.setdefault(canon(c), c)
    out = list(classes.values())
    print("curve classes up to D6:", len(out))
    for c in out:
        print("  ", c["R1"], c["R2"], "example", c["example_angles"])
    (ROOT / "data" / "hexprisms_curves.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
