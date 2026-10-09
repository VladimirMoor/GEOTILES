"""Звёзды рёбер для разбиений «лицом к лицу» (face-to-face, f2f).

В f2f-разбиении тетраэдрами вокруг каждого ребра [P, Q] плитки стоят по кругу. Соседние плитки делят
треугольную грань с основанием PQ, которая задаётся парой расстояний (|P·apex|, |Q·apex|). Плитка,
у которой ребро (a, b) лежит на PQ (a в P, b в Q), — это клин между гранями (a, b, c) и (a, b, d) с
углом α_ab. Значит, звезда ребра — замкнутый маршрут в графе «грань → грань» с суммой углов ровно 2π.
Зеркальные копии допустимы автоматически, потому что используются только метрические данные.

Если у какого-то ребра e тетраэдра T нет звезды, проходящей через саму плитку T (с ребром e на PQ),
то f2f-разбиения нет. Предложение 2.7 у Chentouf–Sun — частный случай этого критерия.

Равенство длин определяется в 60-значной арифметике (порог 1e-40). Длины — алгебраические числа,
и все совпадения, которые мы используем, можно проверить точно (см. exact_lengths).
"""
import itertools
from fractions import Fraction as Fr

import mpmath as mp

from tetra import EDGES, PRINTED_A, EXCLUDED_BY_LP, Tet, angle_of

TOL = mp.mpf(10) ** -40


class Lengths:
    """Классы равных длин рёбер: номер класса для каждого ребра."""

    def __init__(self, T):
        vals = []
        self.cls = {}
        for e in EDGES:
            L = T.len[e]
            for k, v in enumerate(vals):
                if abs(v - L) < TOL:
                    self.cls[e] = k
                    break
            else:
                self.cls[e] = len(vals)
                vals.append(L)
        self.vals = vals

    def __call__(self, i, j):
        return self.cls[tuple(sorted((i, j)))]


def wedges(T, Lc, base_cls):
    """Все клинья (переходы) для основания класса base_cls: (вход, выход, угол, описание)."""
    out = []
    for (a, b) in EDGES:
        if Lc(a, b) != base_cls:
            continue
        c, d = [x for x in (1, 2, 3, 4) if x not in (a, b)]
        ang = angle_of(T.angles, a, b)
        for (p, q) in ((a, b), (b, a)):  # p в точке P, q в точке Q
            for (s, t) in ((c, d), (d, c)):  # входим через грань (p,q,s), выходим через (p,q,t)
                fin = (Lc(p, s), Lc(q, s))
                fout = (Lc(p, t), Lc(q, t))
                out.append((fin, fout, ang, (p, q, s, t)))
    return out


def edge_star_exists(T, e, Lc=None, want_all=False):
    """Есть ли замкнутый маршрут суммы 2, содержащий клин самой T с ребром e."""
    Lc = Lc or Lengths(T)
    a, b = e
    W = wedges(T, Lc, Lc(a, b))
    starts = [w for w in W if {w[3][0], w[3][1]} == {a, b}]
    sols = []
    for w0 in starts:
        # ищем путь от выхода w0 к входу w0 с суммой 2 − α
        target_face = w0[0]
        need = Fr(2) - w0[2]
        # DP по (грань, накопленный угол)
        frontier = {(w0[1], Fr(0)): [w0]}
        seen = set(frontier)
        while frontier:
            new = {}
            for (face, acc), path in frontier.items():
                if acc == need and face == target_face:
                    sols.append(path)
                    if not want_all:
                        return True, sols
                    continue
                for w in W:
                    if w[0] == face and acc + w[2] <= need:
                        st = (w[1], acc + w[2])
                        if st not in seen:
                            seen.add(st)
                            new[st] = path + [w]
            frontier = new
    return bool(sols), sols


def report(T):
    Lc = Lengths(T)
    res = {}
    for e in EDGES:
        ok, sols = edge_star_exists(T, e, Lc)
        res[e] = (ok, sols[0] if sols else None)
    return Lc, res


if __name__ == "__main__":
    for idx, a in enumerate(PRINTED_A, 1):
        if a in EXCLUDED_BY_LP:
            continue
        T = Tet(a)
        Lc, res = report(T)
        bad = [f"{i}{j}" for (i, j), (ok, _) in res.items() if not ok]
        ncls = len(Lc.vals)
        verdict = "NO f2f tiling (edge " + ",".join(bad) + " has no star)" if bad else "edge stars exist for all edges"
        print(f"#{idx:2d} {tuple(map(str, a))} length classes={ncls}: {verdict}")
