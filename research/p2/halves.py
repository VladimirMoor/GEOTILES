"""Половинки симметричных тетраэдров.

Если у T есть зеркальная симметрия (перестановка i ↔ j, сохраняющая все длины), плоскость симметрии
проходит через k, l и середину m ребра ij и режет T на две конгруэнтные половинки T' = (i, k, l, m).
Любое разбиение копиями T даёт разбиение копиями T' (каждую плитку разрезаем). Поэтому:
    T' не замощает  ⟹  T не замощает.
Двугранные углы T' рациональны: на рёбрах в плоскости симметрии — половины углов T, на рёбрах,
перпендикулярных зеркалу, — те же, на новых рёбрах — π/2 (грани перпендикулярны зеркалу).
К T' применяем все критерии: рёберный LP, предложение 2.3 (не f2f) и звёзды рёбер (f2f).
"""
import itertools
from fractions import Fraction as Fr

import mpmath as mp

from tetra import EDGES, PRINTED_A, EXCLUDED_BY_LP, Tet, angle_of
from obstruct import combos, edge_lp, cone_contains_float
from stars import Lengths, edge_star_exists
from nonf2f import status as nonf2f_status

mp.mp.dps = 60


def mirror_symmetries(a):
    out = []
    for i, j in itertools.combinations((1, 2, 3, 4), 2):
        perm = {1: 1, 2: 2, 3: 3, 4: 4}
        perm[i], perm[j] = j, i
        if all(angle_of(a, *e) == angle_of(a, perm[e[0]], perm[e[1]]) for e in EDGES):
            out.append((i, j))
    return out


def dihedral(P, i, j):
    """Двугранный угол тетраэдра P (словарь 1..4 → точка) вдоль ребра ij, в долях π."""
    k, l = [x for x in (1, 2, 3, 4) if x not in (i, j)]
    a, b = P[i], P[j]
    d = (b - a) / mp.norm(b - a)

    def perp(p):
        w = p - a
        return w - d * (w.T * d)[0]
    u, v = perp(P[k]), perp(P[l])
    return mp.acos((u.T * v)[0] / (mp.norm(u) * mp.norm(v))) / mp.pi


def half(a, i, j):
    T = Tet(a)
    k, l = [x for x in (1, 2, 3, 4) if x not in (i, j)]
    m = (T.V[i] + T.V[j]) / 2
    P = {1: T.V[i], 2: T.V[k], 3: T.V[l], 4: m}
    ang = []
    for e in EDGES:
        x = dihedral(P, *e)
        f = Fr(str(mp.nstr(x, 40))).limit_denominator(1000)
        assert abs(x - mp.mpf(f.numerator) / f.denominator) < mp.mpf(10) ** -35, (x, f)
        ang.append(f)
    return tuple(ang)


def verdict(h):
    Th = Tet(h)
    cs, lens = edge_lp(Th)
    lp_ok, _ = cone_contains_float([list(c) for c in cs], lens)
    s23, _ = nonf2f_status(h)
    Lc = Lengths(Th)
    bad = [f"{e[0]}{e[1]}" for e in EDGES if not edge_star_exists(Th, e, Lc)[0]]
    no_f2f = bool(bad)
    no_nonf2f = s23.startswith("f2f only")
    excluded = (not lp_ok) or (no_f2f and no_nonf2f)
    return {"LP": "pass" if lp_ok else "FAIL", "prop2.3": s23, "edges_without_star": bad, "excluded": excluded}


if __name__ == "__main__":
    for idx, a in enumerate(PRINTED_A, 1):
        if a in EXCLUDED_BY_LP or idx in (1, 3):
            continue
        syms = mirror_symmetries(a)
        if not syms:
            continue
        for (i, j) in syms:
            h = half(a, i, j)
            v = verdict(h)
            print(f"#{idx} mirror ({i}{j}) half T' = {tuple(map(str, h))}: {v}", flush=True)
