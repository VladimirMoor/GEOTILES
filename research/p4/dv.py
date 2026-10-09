"""Ячейки Дирихле–Вороного (DV-стереоэдры) орбит пространственных групп.

Для группы Γ (операции из gemmi, кубическая решётка с a = 1) и точки x строим орбиту Γx в окрестности
x, полупространства биссекторов {z : (y − x)·z ≤ (|y|² − |x|²)/2} и их пересечение — ячейку x.
Число граней — число биссекторов, которые реально высекают двумерную грань.
"""
import itertools
from fractions import Fraction as Fr

import gemmi
import numpy as np
from scipy.spatial import HalfspaceIntersection, ConvexHull


class Group:
    def __init__(self, hm, cell=np.eye(3)):
        sg = gemmi.SpaceGroup(hm)
        self.number = sg.number
        self.ops = []
        for op in sg.operations():
            R = np.array(op.rot, dtype=float) / op.DEN
            t = np.array(op.tran, dtype=float) / op.DEN
            self.ops.append((R, t))
        self.cell = cell  # столбцы — базисные векторы решётки (декартовы)
        shifts = np.array(list(itertools.product(range(-2, 3), repeat=3)), dtype=float)
        self.shifts = shifts

    def orbit(self, x, radius):
        """Декартовы точки орбиты (кроме самой x) в шаре радиуса radius вокруг x (x — во фракционных)."""
        pts = []
        for R, t in self.ops:
            y = R @ x + t
            for n in self.shifts:
                z = y + n
                d = self.cell @ (z - x)
                if 1e-9 < np.linalg.norm(d) < radius:
                    pts.append(self.cell @ z)
        return np.array(pts)


def dv_cell(G, x, radius=1.0, tol=1e-9):
    """Возвращает (число граней, f-вектор, индексы соседей) для DV-ячейки точки x (фракционные координаты)."""
    X = G.cell @ x
    Y = G.orbit(x, radius)
    # уникальные соседи (на спецпозициях орбита может повторяться)
    Y = np.unique(np.round(Y, 12), axis=0)
    A = Y - X
    b = (np.sum(Y * Y, axis=1) - X @ X) / 2
    hs = np.hstack([A, -b[:, None]])  # A z − b ≤ 0
    H = HalfspaceIntersection(hs, X)
    V = H.intersections
    # для каждой плоскости — вершины на ней
    facets = []
    for i in range(len(A)):
        on = np.abs(V @ A[i] - b[i]) < tol * max(1.0, np.linalg.norm(A[i]))
        if on.sum() >= 3:
            P = V[on]
            if np.linalg.matrix_rank(P - P[0], tol=1e-7) >= 2:
                facets.append(i)
    # f-вектор по выпуклой оболочке вершин
    Vu = np.unique(np.round(V, 10), axis=0)
    hull = ConvexHull(Vu)
    nV = len(Vu)
    nF = len(facets)
    nE = nV + nF - 2
    return nF, (nV, nE, nF), Y[facets]


if __name__ == "__main__":
    G = Group("I 41 3 2")
    for name, x in [("Schmitt 38", (427 / 6984, 761 / 6984, 1421 / 6984)),
                    ("Schmitt 37", (443 / 6984, 259 / 2328, 1387 / 6984)),
                    ("Schmitt 36", (445 / 6984, 259 / 2328, 347 / 1746))]:
        nF, f, _ = dv_cell(G, np.array(x))
        print(name, "→ facets", nF, "f-vector", f)
