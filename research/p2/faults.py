"""Плоскости разлома: плоские звёзды вершин.

Пусть разбиение не face-to-face, Π — плоскость разлома, S — множество точек Π, лежащих в гранях плиток
с обеих сторон. Грани в Π с «+»-стороны образуют плоское разбиение 𝒫₊ области S треугольниками,
конгруэнтными граням T. То же верно для «−»-стороны (𝒫₋), и 𝒫₊ ≠ 𝒫₋.

(E) Вдоль внутреннего отрезка σ разбиения 𝒫₊ (треугольники G, G' по разные стороны, их рёбра g, g' ⊃ σ)
    плитки над Π образуют полузвезду суммы π. Отсюда π − α_g − α_{g'} ∈ ℕ⟨α⟩, то есть это
    неотрицательная целая комбинация двугранных углов, возможно пустая.
(V) Во внутренней вершине q разбиения 𝒫₊ углы треугольников при q дают в сумме 2π (полная звезда).
    Если q лежит внутри ребра другого треугольника 𝒫₊, то сумма с одной стороны равна π (T-стык).

Здесь перечисляются все плоские звёзды, удовлетворяющие (E) и (V). Если их нет, у разбиений 𝒫± нет
внутренних вершин.
Углы граней иррациональны и сравниваются в 60-значной арифметике с порогом 1e-40. Настоящее
соотношение даёт расхождение около 1e-60, поэтому ни одно истинное соотношение не теряется, и
утверждения о несуществовании надёжны.
"""
import itertools
import sys
from fractions import Fraction as Fr
from functools import lru_cache

import mpmath as mp

from tetra import EDGES, PRINTED_A, EXCLUDED_BY_LP, Tet, angle_of
from obstruct import combos
from faceangles import face_angles

mp.mp.dps = 60
TOL = mp.mpf(10) ** -40
KNOWN = {1, 3}


def nat_sums(T, cap=Fr(1)):
    """Все значения Σ k_e α_e ≤ cap (в долях π), k ≥ 0, включая 0."""
    a = [angle_of(T.angles, *e) for e in EDGES]
    vals = {Fr(0)}
    frontier = [Fr(0)]
    while frontier:
        x = frontier.pop()
        for y in a:
            z = x + y
            if z <= cap and z not in vals:
                vals.add(z)
                frontier.append(z)
    return vals


class Planar:
    def __init__(self, T):
        self.T = T
        self.beta = face_angles(T)  # (k, v) -> угол
        self.sums = nat_sums(T)
        self.alpha = {e: angle_of(T.angles, *e) for e in EDGES}
        # уголки: (k, v, w1, w2): треугольник — грань k, вершина v в точке q, первый луч (v,w1), второй (v,w2)
        self.corners = []
        for k in (1, 2, 3, 4):
            F = [x for x in (1, 2, 3, 4) if x != k]
            for v in F:
                w1, w2 = [x for x in F if x != v]
                self.corners.append((k, v, w1, w2))
                self.corners.append((k, v, w2, w1))

    def compat(self, e1, e2):
        a1 = self.alpha[tuple(sorted(e1))]
        a2 = self.alpha[tuple(sorted(e2))]
        return (Fr(1) - a1 - a2) in self.sums

    def ang(self, c):
        return self.beta[(c[0], c[1])]


def value_classes(P):
    vals = []
    for c in P.corners:
        x = P.ang(c)
        for i, v in enumerate(vals):
            if abs(v - x) < TOL:
                break
        else:
            vals.append(x)
    cls = {c: next(i for i, v in enumerate(vals) if abs(v - P.ang(c)) < TOL) for c in P.corners}
    return vals, cls


def multisets(vals, target):
    out = []
    n = len(vals)

    def rec(i, acc, cur):
        if i == n:
            if abs(acc - target) < TOL:
                out.append(tuple(cur))
            return
        k = 0
        while acc + k * vals[i] <= target + TOL:
            rec(i + 1, acc + k * vals[i], cur + [k])
            k += 1
    rec(0, mp.mpf(0), [])
    return [m for m in out if any(m)]


def full_stars(P, limit=None):
    """Циклические последовательности уголков с суммой 2π и совместимыми лучами."""
    vals, cls = value_classes(P)
    sols = []
    for m in multisets(vals, mp.mpf(2)):
        # ищем цикл: фиксируем первый уголок (минимальный класс, чтобы сократить симметрию)
        first_cls = next(i for i, k in enumerate(m) if k)
        for c0 in [c for c in P.corners if cls[c] == first_cls]:
            rem = list(m)
            rem[first_cls] -= 1
            path = [c0]

            def dfs(rem, last):
                if not any(rem):
                    if P.compat((last[1], last[3]), (c0[1], c0[2])):
                        sols.append(tuple(path))
                        return True
                    return False
                found = False
                for c in P.corners:
                    i = cls[c]
                    if rem[i] and P.compat((last[1], last[3]), (c[1], c[2])):
                        rem[i] -= 1
                        path.append(c)
                        if dfs(rem, c):
                            found = True
                        path.pop()
                        rem[i] += 1
                        if found and limit:
                            return True
                return found
            dfs(rem, c0)
            if limit and len(sols) >= limit:
                return sols
    return sols


def t_stars(P, limit=None):
    """T-стык: q внутри ребра e0 = (a, b) треугольника; с другой стороны уголки суммы π."""
    vals, cls = value_classes(P)
    sols = []
    ms = multisets(vals, mp.mpf(1))
    for e0 in EDGES:
        for m in ms:
            def dfs(rem, last, path):
                if not any(rem):
                    if P.compat((last[1], last[3]), e0):
                        sols.append((e0, tuple(path)))
                        return True
                    return False
                found = False
                for c in P.corners:
                    i = cls[c]
                    prev = e0 if last is None else (last[1], last[3])
                    if rem[i] and P.compat(prev, (c[1], c[2])):
                        rem[i] -= 1
                        path.append(c)
                        if dfs(rem, c, path):
                            found = True
                        path.pop()
                        rem[i] += 1
                        if found and limit:
                            return True
                return found
            dfs(list(m), None, [])
            if limit and len(sols) >= limit:
                return sols
    return sols


def fmt_corner(c):
    return f"F{c[0]}@{c[1]}[{c[1]}{c[2]}|{c[1]}{c[3]}]"


if __name__ == "__main__":
    which = [int(x) for x in sys.argv[1:]] or None
    for idx, a in enumerate(PRINTED_A, 1):
        if a in EXCLUDED_BY_LP or (which and idx not in which):
            continue
        T = Tet(a)
        P = Planar(T)
        fs = full_stars(P, limit=1)
        ts = t_stars(P, limit=1)
        tag = " (known tiler)" if idx in KNOWN else ""
        print(f"#{idx:2d}{tag}: full vertex star {'YES' if fs else 'none'}; T-junction {'YES' if ts else 'none'}"
              + (f"  e.g. {[fmt_corner(c) for c in fs[0]]}" if fs else "")
              + (f"  T e.g. {ts[0][0]} {[fmt_corner(c) for c in ts[0][1]]}" if ts else ""), flush=True)
