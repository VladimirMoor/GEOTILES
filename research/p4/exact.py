"""Точный сертификат нижней оценки числа граней DV-ячейки (рациональная арифметика).

Работаем во фракционных координатах с метрическим тензором M = BᵀB, где M рационален, если рационально c²
(для гексагональной решётки M = [[1, −1/2, 0], [−1/2, 1, 0], [0, 0, c²]]). Ячейка точки x:
    { z : 2 (y − x)ᵀ M z ≤ yᵀMy − xᵀMx  для всех y ∈ Γx \\ {x} }.
Сертификат «не меньше k граней»: k различных плоскостей, на каждой — три точно допустимые неколлинеарные
точки ячейки (вершины, найденные в плавающей точке и пересчитанные точно как пересечения трёх плоскостей).
"""
import itertools
from fractions import Fraction as Fr

import numpy as np

from dvfast import Group, Evaluator


def frac_vec(v):
    return [Fr(x) for x in v]


def metric(G, c2):
    if G.system == "cubic":
        return [[Fr(1), Fr(0), Fr(0)], [Fr(0), Fr(1), Fr(0)], [Fr(0), Fr(0), Fr(1)]]
    if G.system == "tetragonal":
        return [[Fr(1), Fr(0), Fr(0)], [Fr(0), Fr(1), Fr(0)], [Fr(0), Fr(0), c2]]
    return [[Fr(1), Fr(-1, 2), Fr(0)], [Fr(-1, 2), Fr(1), Fr(0)], [Fr(0), Fr(0), c2]]


def mv(M, v):
    return [sum(M[i][j] * v[j] for j in range(3)) for i in range(3)]


def dot(u, v):
    return sum(a * b for a, b in zip(u, v))


def solve3(A, b):
    """Точное решение 3×3 (Крамер)."""
    def det(m):
        return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
                + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))
    D = det(A)
    if D == 0:
        return None
    out = []
    for k in range(3):
        m = [row[:] for row in A]
        for i in range(3):
            m[i][k] = b[i]
        out.append(det(m) / D)
    return out


def certify(hm, x, c2=Fr(1), radius_factor=1.0):
    """x — тройка Fraction (фракционные координаты), c2 — (c/a)² как Fraction. Возвращает число граней с сертификатом."""
    G = Group(hm)
    c = float(c2) ** 0.5
    E = Evaluator(G, c)
    xf = np.array([float(t) for t in x])
    X, A, b, V = E.cell_of(xf)  # плавающая точка: кандидаты
    Binv = np.linalg.inv(E.B)
    M = metric(G, c2)
    # точные соседи: те же операции, точные сдвиги
    R = [[[Fr(int(round(v * 24)), 24) for v in row] for row in Rm] for Rm in G.R]
    T = [[Fr(int(round(v * 24)), 24) for v in tv] for tv in G.t]
    Yf = (Binv @ (A + X).T).T  # фракционные координаты соседей (float)
    planes = []
    for yf in Yf:
        # восстановим точную точку орбиты: найдём операцию и сдвиг
        best = None
        for Rm, tv in zip(R, T):
            base = [sum(Rm[i][j] * x[j] for j in range(3)) + tv[i] for i in range(3)]
            n = [round(yf[i] - float(base[i])) for i in range(3)]
            y = [base[i] + n[i] for i in range(3)]
            if max(abs(float(y[i]) - yf[i]) for i in range(3)) < 1e-9:
                best = y
                break
        assert best is not None
        y = best
        a = [2 * t for t in mv(M, [y[i] - x[i] for i in range(3)])]
        rhs = dot(y, mv(M, y)) - dot(x, mv(M, x))
        planes.append((a, rhs))
    # кандидатные вершины: тройки плоскостей, активных в плавающей точке
    Vf = (Binv @ V.T).T
    Af = np.array([[float(t) for t in a] for a, _ in planes])
    bf = np.array([float(r) for _, r in planes])
    exact_vertices = []
    for v in Vf:
        act = np.where(np.abs(Af @ v - bf) < 1e-7 * np.maximum(1, np.abs(bf)))[0]
        found = None
        for i, j, k in itertools.combinations(act, 3):
            sol = solve3([planes[i][0], planes[j][0], planes[k][0]], [planes[i][1], planes[j][1], planes[k][1]])
            if sol is not None:
                found = sol
                break
        if found is None:
            continue
        if all(dot(a, found) <= r for a, r in planes):  # точно допустима
            exact_vertices.append(found)
    # грани: плоскости, на которых ≥ 3 неколлинеарные точные вершины
    nf = 0
    for a, r in planes:
        on = [v for v in exact_vertices if dot(a, v) == r]
        ok = False
        for p, q, s in itertools.combinations(on, 3):
            u = [q[i] - p[i] for i in range(3)]
            w = [s[i] - p[i] for i in range(3)]
            cr = [u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0]]
            if any(cr):
                ok = True
                break
        nf += ok
    return nf, len(exact_vertices), len(planes)


if __name__ == "__main__":
    print("214 Schmitt 38:", certify("I 41 3 2", (Fr(427, 6984), Fr(761, 6984), Fr(1421, 6984))))
    print("98 Schmitt 35:", certify("I 41 2 2", (Fr(62, 125), Fr(41, 125), Fr(79, 1000)), Fr(727, 500) ** 2))
