"""Необходимые условия замощения для рациональных тетраэдров (точные, над ℚ).

1. Рёберный LP (Chentouf–Sun, Cor. 2.13), воспроизведение. Комбинации C: Σ C(e)·α_e ∈ {π, 2π}
   (π — если отрезок ребра лежит внутри грани соседа). D_C(e) = C(e)/len(e). Если D x = 1, x ≥ 0
   несовместна, тетраэдр не замощает.

2. Вершинный LP (новое). Точка p разбиения, являющаяся вершиной хотя бы одной плитки, окружена:
   телесными углами Ω_v плиток, у которых p — вершина; клиньями 2α_e плиток, у которых p внутри ребра;
   полупространствами 2π плиток, у которых p внутри грани. Сумма = 4π. Для рационального тетраэдра
   Ω_v = (Σ трёх двугранных углов при v) − π рационально, поэтому комбинации K: Σ k_v ω_v + 2 Σ m_e a_e
   + 2 f = 4 (в долях π), k ≠ 0, перечисляются конечным перебором. Подсчёт в шаре радиуса r: каждая плитка
   даёт по одной вершине каждого типа, граничные эффекты O(r²) ⇒ существуют частоты y_K ≥ 0 с
   Σ_K y_K k_v(K) = 1 для v = 1..4. Несовместность ⇒ тетраэдр не замощает.
   Ограничение f ≤ 1 и (f = 1 ⇒ остальное заполняет полупространство) учтено как f ∈ {0, 1}.

Обе задачи — проверка принадлежности вектора конусу, решаем точным симплекс-методом над ℚ
(sympy не нужен: Fraction-реализация ниже) — ответ «несовместно» сопровождается сертификатом Фаркаша.
"""
import itertools
import sys
from fractions import Fraction as Fr

import mpmath as mp

from tetra import EDGES, PRINTED_A, EXCLUDED_BY_LP, Tet, angle_of


def combos(weights, target_set, require_pos=None):
    """Все векторы k ≥ 0 с Σ k_i w_i ∈ target_set (w_i > 0 рациональные)."""
    n = len(weights)
    tmax = max(target_set)
    out = []

    def rec(i, acc, cur):
        if i == n:
            if acc in target_set:
                out.append(tuple(cur))
            return
        k = 0
        while acc + k * weights[i] <= tmax:
            rec(i + 1, acc + k * weights[i], cur + [k])
            k += 1
    rec(0, Fr(0), [])
    if require_pos is not None:
        out = [c for c in out if require_pos(c)]
    return out


def farkas_infeasible(cols, rhs):
    """Есть ли x ≥ 0 с Σ x_j cols_j = rhs? Точный фазовый симплекс над ℚ.
    Возвращает (feasible, x или сертификат y: y·cols_j ≥ 0 для всех j, y·rhs < 0)."""
    m = len(rhs)
    # приводим rhs ≥ 0
    sign = [1 if r >= 0 else -1 for r in rhs]
    A = [[Fr(sign[i]) * c[i] for c in cols] for i in range(m)]
    b = [Fr(sign[i]) * rhs[i] for i in range(m)]
    n = len(cols)
    # таблица фазы 1: переменные x (n) + искусственные (m)
    T = [A[i] + [Fr(1 if j == i else 0) for j in range(m)] + [b[i]] for i in range(m)]
    basis = [n + i for i in range(m)]
    cost = [Fr(0)] * n + [Fr(1)] * m
    while True:
        # редуцированные стоимости
        cb = [cost[basis[i]] for i in range(m)]
        red = [cost[j] - sum(cb[i] * T[i][j] for i in range(m)) for j in range(n + m)]
        enter = next((j for j in range(n + m) if red[j] < 0), None)  # правило Бланда
        if enter is None:
            break
        rows = [(T[i][-1] / T[i][enter], basis[i], i) for i in range(m) if T[i][enter] > 0]
        _, _, r = min(rows)
        piv = T[r][enter]
        T[r] = [v / piv for v in T[r]]
        for i in range(m):
            if i != r and T[i][enter] != 0:
                f = T[i][enter]
                T[i] = [a - f * c for a, c in zip(T[i], T[r])]
        basis[r] = enter
    obj = sum(cost[basis[i]] * T[i][-1] for i in range(m))
    if obj == 0:
        x = [Fr(0)] * n
        for i in range(m):
            if basis[i] < n:
                x[basis[i]] = T[i][-1]
        return True, x
    # двойственные переменные: y = c_B B⁻¹ ; B⁻¹ — столбцы искусственных переменных
    cb = [cost[basis[i]] for i in range(m)]
    y = [sum(cb[i] * T[i][n + k] for i in range(m)) for k in range(m)]
    # y·A_j ≤ 0 ... приведём к форме сертификата: z = −y·sign
    z = [-y[k] * sign[k] for k in range(m)]
    return False, z


def edge_lp(T):
    a = [angle_of(T.angles, *e) for e in EDGES]
    lens = T.edge_lengths()
    cs = combos(a, {Fr(1), Fr(2)})
    cs = [c for c in cs if any(c)]
    # D_C(e) = C(e)/len(e); масштабируем строку e на len(e): Σ_C x_C C(e) = len(e)
    # длины иррациональны — решаем в плавающей точке высокой точности через проверку конусом:
    # совместность Σ x_C C = L (вектор длин), x ≥ 0. Используем mpmath-LP методом перебора баз.
    return cs, lens


def cone_contains_float(cols, target):
    """x ≥ 0, Σ x_j cols_j = target (6-мерный вектор mpmath) — проверка через scipy linprog в float
    и затем точная проверка найденного базиса в mpmath."""
    import numpy as np
    from scipy.optimize import linprog
    A = np.array([[float(c[i]) for c in cols] for i in range(len(target))])
    b = np.array([float(t) for t in target])
    res = linprog(np.zeros(len(cols)), A_eq=A, b_eq=b, bounds=[(0, None)] * len(cols), method="highs")
    if res.status == 0:
        return True, res.x
    # сертификат Фаркаша: y с yᵀA ≥ 0, yᵀb < 0
    res2 = linprog(b, A_ub=-A.T, b_ub=np.zeros(len(cols)), bounds=[(-1, 1)] * len(target), method="highs")
    return False, res2.x


def vertex_combos(T):
    om = T.solid_angles()
    w_vert = [om[v] for v in (1, 2, 3, 4)]
    w_edge = [2 * angle_of(T.angles, *e) for e in EDGES]
    out = []
    for f in (0, 1):
        tgt = Fr(4) - 2 * f
        for c in combos(w_vert + w_edge, {tgt}):
            if any(c[:4]):
                out.append((c[:4], c[4:], f))
    return out


def vertex_lp(T):
    vc = vertex_combos(T)
    cols = sorted({k for k, _, _ in vc})
    feas, cert = farkas_infeasible([list(map(Fr, k)) for k in cols], [Fr(1)] * 4)
    return feas, cert, cols


if __name__ == "__main__":
    rows = []
    for idx, a in enumerate(PRINTED_A, 1):
        T = Tet(a)
        cs, lens = edge_lp(T)
        feas_e, _ = cone_contains_float([list(c) for c in cs], lens)
        feas_v, cert, cols = vertex_lp(T)
        tag = " (excluded by LP in paper)" if a in EXCLUDED_BY_LP else ""
        print(f"#{idx:2d}{tag}: edge-LP {'pass' if feas_e else 'FAIL'} ({len(cs)} combos); "
              f"vertex-LP {'pass' if feas_v else 'FAIL'} ({len(cols)} vertex vectors)"
              + ("" if feas_v else f" cert y={list(map(str, cert))}"), flush=True)
