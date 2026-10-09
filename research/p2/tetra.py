"""Рациональные тетраэдры: данные и геометрия.

Обозначения как у Kedlaya–Kolpakov–Poonen–Rubinstein и Chentouf–Sun (arXiv:2312.01654):
вершины 1..4, двугранные углы (α12, α34, α13, α24, α14, α23) в долях π, αij — угол вдоль ребра ij.

SPORADIC_A — множество 𝒜 из теоремы 1.2 у Chentouf–Sun: спорадические тетраэдры, для которых
не доказано, что они НЕ замощают пространство. В напечатанном списке 42 кортежа; два из них
(EXCLUDED_BY_LP) авторы в §3.1 сами исключают по LP-критерию, так что настоящее 𝒜 — это 40 остальных.
"""
from fractions import Fraction as Fr
import itertools

import mpmath as mp

mp.mp.dps = 60

EDGES = [(1, 2), (3, 4), (1, 3), (2, 4), (1, 4), (2, 3)]


def _t(*xs):
    return tuple(Fr(x) for x in xs)


PRINTED_A = [
    _t("1/4", "1/3", "1/3", "1/4", "1/2", "2/3"),
    _t("5/24", "3/8", "1/3", "1/4", "13/24", "5/8"),
    _t("1/4", "1/2", "1/2", "1/3", "1/3", "1/2"),
    _t("7/24", "11/24", "13/24", "7/24", "1/3", "1/2"),
    _t("1/5", "1/5", "1/3", "1/5", "2/3", "2/3"),
    _t("1/5", "1/5", "4/15", "4/15", "3/5", "11/15"),
    _t("1/5", "1/5", "1/3", "1/3", "3/5", "3/5"),
    _t("1/3", "1/3", "3/5", "1/3", "2/5", "2/5"),
    _t("4/15", "8/15", "1/3", "1/3", "7/15", "7/15"),
    _t("1/7", "3/7", "1/3", "1/3", "4/7", "4/7"),
    _t("2/7", "2/7", "1/3", "1/3", "3/7", "5/7"),
    _t("1/5", "2/5", "1/2", "1/3", "1/3", "2/3"),
    _t("2/15", "7/15", "1/2", "1/3", "2/5", "3/5"),
    _t("2/15", "7/15", "31/60", "19/60", "5/12", "7/12"),
    _t("1/5", "2/5", "7/12", "1/4", "5/12", "7/12"),
    _t("13/60", "23/60", "7/12", "1/4", "2/5", "3/5"),
    _t("1/5", "3/5", "1/3", "1/3", "1/2", "1/2"),
    _t("3/10", "7/10", "1/3", "1/3", "2/5", "2/5"),
    _t("1/5", "1/5", "2/5", "1/3", "1/2", "2/3"),
    _t("1/6", "7/30", "11/30", "11/30", "1/2", "2/3"),
    _t("1/5", "1/5", "9/20", "17/60", "11/20", "37/60"),
    _t("1/5", "1/3", "1/2", "1/3", "2/5", "3/5"),
    _t("1/6", "11/30", "1/2", "1/3", "13/30", "17/30"),
    _t("1/6", "11/30", "29/60", "7/20", "5/12", "7/12"),
    _t("1/5", "1/3", "1/3", "1/5", "1/2", "4/5"),
    _t("7/60", "5/12", "1/3", "1/5", "7/12", "43/60"),
    _t("1/5", "2/5", "2/5", "1/5", "1/2", "2/3"),
    _t("1/5", "2/5", "1/3", "1/3", "1/2", "3/5"),
    _t("7/30", "13/30", "3/10", "3/10", "1/2", "3/5"),
    _t("1/4", "7/20", "1/3", "1/3", "9/20", "13/20"),
    _t("1/5", "1/2", "3/5", "1/5", "1/3", "2/3"),
    _t("1/5", "1/2", "17/30", "7/30", "3/10", "7/10"),
    _t("3/20", "11/20", "17/30", "7/30", "7/20", "13/20"),
    _t("3/20", "11/20", "11/20", "1/4", "1/3", "2/3"),
    _t("11/60", "31/60", "11/20", "1/4", "3/10", "7/10"),
    _t("1/5", "1/2", "1/2", "1/3", "2/5", "1/2"),
    _t("1/5", "1/2", "7/15", "11/30", "11/30", "8/15"),
    _t("4/15", "13/30", "17/30", "4/15", "2/5", "1/2"),
    _t("4/15", "13/30", "3/5", "3/10", "11/30", "7/15"),
    _t("4/15", "17/30", "2/5", "3/10", "11/30", "8/15"),
    _t("3/10", "2/5", "3/5", "3/10", "1/3", "1/2"),
    _t("1/3", "2/5", "2/5", "1/3", "1/2", "2/5"),
]
EXCLUDED_BY_LP = [_t("3/20", "11/20", "11/20", "1/4", "1/3", "2/3"),
                  _t("11/60", "31/60", "11/20", "1/4", "3/10", "7/10")]
SPORADIC_A = [t for t in PRINTED_A if t not in EXCLUDED_BY_LP]


def angle_of(angles, i, j):
    return angles[EDGES.index(tuple(sorted((i, j))))]


class Tet:
    """Тетраэдр по двугранным углам (в долях π): нормали, вершины, длины рёбер."""

    def __init__(self, angles, name=""):
        self.angles = tuple(angles)
        self.name = name
        # грань k противоположна вершине k; ребро ij лежит на гранях {k, l} = дополнение
        G = mp.matrix(4, 4)
        for k in range(4):
            G[k, k] = 1
        for (i, j), a in zip(EDGES, angles):
            k, l = [x - 1 for x in (1, 2, 3, 4) if x not in (i, j)]
            G[k, l] = G[l, k] = -mp.cos(mp.pi * mp.mpf(a.numerator) / a.denominator)
        self.gram = G
        E, Q = mp.eigsy(G)
        self.gram_eigs = [E[i] for i in range(4)]
        # ранг 3, положительно полуопределённая: одно нулевое собственное значение
        order = sorted(range(4), key=lambda i: E[i])
        assert abs(E[order[0]]) < mp.mpf(10) ** -40, f"det≠0: {E[order[0]]}"
        assert E[order[1]] > 0
        N = []
        for k in range(4):
            N.append(mp.matrix([Q[k, i] * mp.sqrt(E[i]) for i in order[1:]]))
        self.normals = N  # внешние единичные нормали, n_k·n_l = −cos
        ker = [Q[k, order[0]] for k in range(4)]
        if ker[0] < 0:
            ker = [-x for x in ker]
        assert all(x > 0 for x in ker), "не тетраэдр (ядро не положительно)"
        V = {}
        for i in range(4):
            rows = [k for k in range(4) if k != i]
            A = mp.matrix([[N[k][c] for c in range(3)] for k in rows])
            V[i + 1] = mp.lu_solve(A, mp.matrix([1, 1, 1]))
        self.V = V
        vol = abs(mp.det(mp.matrix([[V[j][c] - V[1][c] for c in range(3)] for j in (2, 3, 4)]))) / 6
        s = vol ** (-mp.mpf(1) / 3)  # нормировка на единичный объём
        for i in V:
            V[i] = V[i] * s
        self.volume = mp.mpf(1)
        self.len = {e: mp.norm(V[e[0]] - V[e[1]]) for e in EDGES}

    def edge_lengths(self):
        return [self.len[e] for e in EDGES]

    def dihedral_check(self):
        """Пересчёт двугранных углов по вершинам (для контроля)."""
        out = []
        for (i, j) in EDGES:
            k, l = [x for x in (1, 2, 3, 4) if x not in (i, j)]
            a, b = self.V[i], self.V[j]
            d = (b - a) / mp.norm(b - a)

            def perp(p):
                w = p - a
                return w - d * (w.T * d)[0]
            u, v = perp(self.V[k]), perp(self.V[l])
            out.append(mp.acos((u.T * v)[0] / (mp.norm(u) * mp.norm(v))) / mp.pi)
        return out

    def solid_angles(self):
        """Телесный угол при вершине v: сумма трёх двугранных углов − π (в долях π)."""
        res = {}
        for v in (1, 2, 3, 4):
            res[v] = sum(angle_of(self.angles, v, w) for w in (1, 2, 3, 4) if w != v) - 1
        return res


if __name__ == "__main__":
    print(len(PRINTED_A), "printed;", len(SPORADIC_A), "in 𝒜")
    for idx, a in enumerate(PRINTED_A, 1):
        T = Tet(a)
        err = max(abs(x - mp.mpf(y.numerator) / y.denominator) for x, y in zip(T.dihedral_check(), a))
        L = [mp.nstr(x, 6) for x in T.edge_lengths()]
        tag = " (excluded by LP)" if a in EXCLUDED_BY_LP else ""
        print(f"#{idx:2d}{tag} err={mp.nstr(err, 3)} lengths={L} solid={T.solid_angles()}")
