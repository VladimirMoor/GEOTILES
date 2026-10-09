"""Сравнение 40 спорадических тетраэдров с известными замостителями.

Известные (по Goldberg 1974 и Bongiovanni–Diaz–Kakkar–Sothanaphan, arXiv:1709.04139, §3, §5):
  * Sommerville №1–4; №1 = T_{π/3}, №2 = симплекс Кокстера [4,3^{1,1}], №3 — из табл. 1 у BDKS;
  * первое семейство Goldberg = семейство Хилла 𝓕₁;
  * второе семейство Goldberg: углы 12:α, 13:π/2, 14:π/3, 23:π/2−α, 24:π/2+β, 34:π/2−β, sin β = cot α / 2;
  * третье семейство Goldberg: 12:α, 13:π/2, 14:π/6, 23:π−2α, 24:π/2−γ, 34:π/2+γ, sin γ = cos α / √3.
Тетраэдр однозначно (с точностью до подобия) задаётся упорядоченными двугранными углами, поэтому
сравниваем углы с учётом 24 перестановок вершин.
"""
import itertools
from fractions import Fraction as Fr

import mpmath as mp

from tetra import EDGES, PRINTED_A, EXCLUDED_BY_LP, SPORADIC_A, Tet

mp.mp.dps = 50


def relabel(angles, perm):
    """perm: старая вершина → новая. Возвращает углы в порядке EDGES для новой разметки."""
    d = {tuple(sorted((perm[i], perm[j]))): a for (i, j), a in zip(EDGES, angles)}
    return tuple(d[e] for e in EDGES)


def canon(angles):
    return min(relabel(angles, dict(zip((1, 2, 3, 4), p))) for p in itertools.permutations((1, 2, 3, 4)))


def from_dict(d):
    return tuple(Fr(d[e]) for e in EDGES)


def rational(x, qmax=720):
    """Рациональное приближение x (в долях π), если оно точное до 1e-40."""
    f = Fr(str(mp.nstr(x, 45))).limit_denominator(qmax)
    return f if abs(x - mp.mpf(f.numerator) / f.denominator) < mp.mpf(10) ** -40 else None


def goldberg2(qmax=240):
    out = []
    for q in range(2, qmax + 1):
        for p in range(1, q):
            a = Fr(p, q)
            if a.denominator != q or not (0 < a < Fr(1, 2)):
                continue
            s = mp.cot(mp.pi * a.numerator / a.denominator) / 2
            if not (0 < s < 1):
                continue
            b = rational(mp.asin(s) / mp.pi)
            if b is None:
                continue
            ang = {(1, 2): a, (1, 3): Fr(1, 2), (1, 4): Fr(1, 3), (2, 3): Fr(1, 2) - a,
                   (2, 4): Fr(1, 2) + b, (3, 4): Fr(1, 2) - b}
            out.append(("Goldberg II", a, b, from_dict(ang)))
    return out


def goldberg3(qmax=240):
    out = []
    for q in range(2, qmax + 1):
        for p in range(1, q):
            a = Fr(p, q)
            if a.denominator != q or not (0 < a < Fr(1, 2)):
                continue
            s = mp.cos(mp.pi * a.numerator / a.denominator) / mp.sqrt(3)
            g = rational(mp.asin(s) / mp.pi)
            if g is None:
                continue
            ang = {(1, 2): a, (1, 3): Fr(1, 2), (1, 4): Fr(1, 6), (2, 3): 1 - 2 * a,
                   (2, 4): Fr(1, 2) - g, (3, 4): Fr(1, 2) + g}
            out.append(("Goldberg III", a, g, from_dict(ang)))
    return out


SOMMERVILLE = {
    # из табл. 1 BDKS: n12 n13 n14 n23 n24 n34, угол 2π/n
    "Sommerville No.1": (4, 6, 6, 6, 6, 4),
    "Sommerville No.2": (4, 4, 4, 6, 6, 8),
    "Sommerville No.3": (3, 6, 6, 8, 8, 4),
}


def sommerville():
    out = {}
    order = [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
    for name, ns in SOMMERVILLE.items():
        d = {e: Fr(2, n) for e, n in zip(order, ns)}
        out[name] = from_dict(d)
    return out


if __name__ == "__main__":
    known = {}
    for name, a in sommerville().items():
        Tet(a)  # проверка, что это тетраэдр
        known[canon(a)] = name
    fams = goldberg2() + goldberg3()
    for name, x, y, a in fams:
        try:
            Tet(a)
        except AssertionError:
            continue
        known.setdefault(canon(a), f"{name} (α={x}π, second angle {y}π)")
    print("rational members of Goldberg II/III found:", [(n, str(x), str(y)) for n, x, y, _ in fams])
    for idx, a in enumerate(PRINTED_A, 1):
        c = canon(a)
        tag = " (excluded by LP)" if a in EXCLUDED_BY_LP else ""
        if c in known:
            print(f"#{idx}{tag}: {tuple(map(str, a))} = {known[c]}")
    # пары тетраэдров из 𝒜, конгруэнтных друг другу (дубликаты)
    seen = {}
    for idx, a in enumerate(PRINTED_A, 1):
        seen.setdefault(canon(a), []).append(idx)
    print("duplicates:", [v for v in seen.values() if len(v) > 1])
