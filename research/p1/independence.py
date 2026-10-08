"""Лемма о независимости (точная проверка над ℚ).

Для суммы Минковского P = Δ₁ ⊕ … ⊕ Δ_k (отрезки и правильные треугольники) двугранный угол ребра из зоны
направления d — это π − ∠_d(a, b), где ∠_d(a, b) — ориентированный угол между проекциями направлений a, b
на плоскость ⊥ d. Поэтому все углы — разности «аргументов» A_d(x) = ∠_d(x₀, x).

Лемма: функции {A_d(x) : d — направление ребра, x ∈ L_d \\ {x₀}} линейно независимы по модулю констант
на пространстве поворотов слагаемых. Здесь L_d — прямые, видимые из d: рёбра остальных слагаемых плюс
след плоскости собственного треугольника.
Следствие: все тождества Σ kᵢθᵢ ≡ const между углами P — циркуляции в графах зон.

Проверка: правильный треугольник conv(e₁, e₂, e₃), отрезок (1,0,0); повороты — Кэли с рациональными
параметрами. Производная atan2(Y, X) равна (X dY − Y dX)/(X² + Y²); здесь Y = (d·(x₀×x))/|d|, где |d| —
константа (повороты сохраняют длину). Поэтому столбец градиента — рациональное число, умноженное на
постоянный множитель, что ранг не меняет. Строим строки (точка, параметр) в нескольких рациональных
точках и вычисляем точный ранг над ℚ. Если ранг равен числу функций, лемма доказана для этого типа.
"""
import itertools
import json
import pathlib
import random
import sys
from fractions import Fraction

import sympy as sp

ROOT = pathlib.Path(__file__).resolve().parent


def cayley_sym(p):
    a, b, c = p
    A = sp.Matrix([[0, -c, b], [c, 0, -a], [-b, a, 0]])
    I = sp.eye(3)
    return (I - A).inv() * (I + A)


def summand_dirs(name):
    """Направления рёбер (и для треугольника — какие пары лежат в одной плоскости)."""
    if name == "seg":
        return [sp.Matrix([1, 0, 0])]
    e1, e2, e3 = sp.Matrix([1, 0, 0]), sp.Matrix([0, 1, 0]), sp.Matrix([0, 0, 1])
    return [e2 - e1, e3 - e2, e1 - e3]


def build(kind):
    names = kind.split("+")
    params = []
    dirs = []  # (summand index, symbolic direction)
    for i, nm in enumerate(names):
        D = summand_dirs(nm)
        if i > 0:
            ps = sp.symbols(f"p{i}_0:3")
            params += list(ps)
            R = cayley_sym(ps)
            D = [R * d for d in D]
        dirs += [(i, d) for d in D]
    funcs = []  # (описание, выражение для градиента: список по параметрам — считаем позже)
    for di, (si, d) in enumerate(dirs):
        # видимые прямые: одно другое ребро своего треугольника (след плоскости) + все рёбра других слагаемых
        own = [x for (sj, x) in dirs if sj == si and x is not d]
        seen = ([own[0]] if own else []) + [x for (sj, x) in dirs if sj != si]
        x0 = seen[0]
        dd = (d.T * d)[0]
        for x in seen[1:]:
            X = (x0.T * x)[0] - (d.T * x0)[0] * (d.T * x)[0] / dd
            Yt = (d.T * x0.cross(x))[0]  # Y = Yt/|d|
            funcs.append((di, X, Yt, dd))
    return params, funcs, len(dirs)


def rank_fraction(rows):
    """Точный ранг матрицы из Fraction (метод Гаусса)."""
    M = [list(r) for r in rows]
    rank, ncol = 0, len(M[0]) if M else 0
    for col in range(ncol):
        piv = next((i for i in range(rank, len(M)) if M[i][col] != 0), None)
        if piv is None:
            continue
        M[rank], M[piv] = M[piv], M[rank]
        inv = 1 / M[rank][col]
        M[rank] = [v * inv for v in M[rank]]
        for i in range(len(M)):
            if i != rank and M[i][col] != 0:
                f = M[i][col]
                M[i] = [a - f * b for a, b in zip(M[i], M[rank])]
        rank += 1
    return rank


def certify(kind, npts=None, seed=1):
    params, funcs, ndirs = build(kind)
    n = len(funcs)
    npts = npts or (n // max(1, len(params)) + 3)
    # градиент atan2 ∝ (X dYt − Yt dX)/(X² + Yt²/dd)  (постоянный множитель 1/|d| отброшен)
    grads = []
    for (di, X, Yt, dd) in funcs:
        den = X ** 2 + Yt ** 2 / dd
        grads.append([sp.together((X * sp.diff(Yt, q) - Yt * sp.diff(X, q)) / den) for q in params])
    lam = [[sp.lambdify(params, g, modules="sympy") for g in col] for col in grads]
    rnd = random.Random(seed)
    rows = []
    got = 0
    while got < npts:
        pt = [sp.Rational(rnd.randint(-9, 9), rnd.randint(1, 7)) for _ in params]
        block = []
        try:
            for j in range(len(params)):
                row = []
                for col in lam:
                    v = sp.sympify(col[j](*pt))
                    if not v.is_finite:
                        raise ValueError
                    v = sp.Rational(v)
                    row.append(Fraction(int(v.p), int(v.q)))
                block.append(row)
        except (ValueError, TypeError, ZeroDivisionError):
            continue  # вырожденная точка (параллельные рёбра) — берём другую
        rows += block
        got += 1
    r = rank_fraction(rows)
    return {"functions": n, "params": len(params), "points": npts, "rank": r, "independent": r == n}


if __name__ == "__main__":
    kinds = sys.argv[1:] or ["tri+tri", "tri+seg+seg", "tri+tri+seg", "tri+seg+seg+seg", "seg+seg+seg+seg+seg",
                             "tri+tri+tri", "tri+tri+seg+seg", "tri+seg+seg+seg+seg", "seg+seg+seg+seg+seg+seg"]
    out = {}
    for k in kinds:
        out[k] = certify(k)
        print(k, out[k], flush=True)
    (ROOT / "data" / "independence.json").write_text(json.dumps(out, indent=1))
