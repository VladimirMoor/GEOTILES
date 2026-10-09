"""Пространство деформаций изоэдрального разбиения с фиксированной комбинаторикой (через вершины).

Неизвестные — координаты всех вершин плитки P. Условия разбиения при фиксированной комбинаторике:
  (L) спаривание граней: для каждого соседа γ элемент γ переводит вершины грани F_{γ⁻¹} в вершины F_γ
      (соответствие вершин берём из DV-ячейки) — линейные уравнения γ v_i − v_j = 0;
  (P) плоскостность каждой грани: (w₂−w₁)×(w₃−w₁)·(w_k−w₁) = 0, k ≥ 4.
Вблизи DV-ячейки любое решение (L)+(P) — снова плитка изоэдрального разбиения той же группы: копии
сходятся по граням, а суммы углов вокруг рёбер непрерывны и кратны 2π, значит, остаются 2π.
Размерность касательного пространства (ядро якобиана) минус 3 (сдвиг порождающей точки x) —
число независимых не-DV деформаций первого порядка.
"""
import sys

import numpy as np

from dvfast import Group, Evaluator


def dv_structure(E, x, tol=1e-9):
    """Вершины DV-ячейки и для каждой грани — (элемент γ = (Rc, tc) в декартовых, индексы вершин)."""
    X, A, b, V = E.cell_of(np.asarray(x, float))
    # уникальные вершины
    Vu = []
    for v in V:
        if not any(np.linalg.norm(v - u) < 1e-8 for u in Vu):
            Vu.append(v)
    Vu = np.array(Vu)
    B = E.B
    Binv = np.linalg.inv(B)
    faces = []
    for i in range(len(A)):
        on = np.where(np.abs(Vu @ A[i] - b[i]) / np.linalg.norm(A[i]) < tol)[0]
        if len(on) >= 3 and np.linalg.matrix_rank(Vu[on] - Vu[on][0], tol=1e-8) >= 2:
            Y = A[i] + X  # соседняя точка орбиты (декартовы)
            # найдём γ: γ X = Y
            yf = Binv @ Y
            xf = Binv @ X
            g = None
            for R, t in zip(E.G.R, E.G.t):
                s = yf - (R @ xf + t)
                if np.allclose(s, np.round(s), atol=1e-8):
                    g = (B @ R @ Binv, B @ (t + np.round(s)))
                    break
            faces.append((g, list(on), Y))
    return X, Vu, faces


def equations(Vu, faces):
    """Линейные уравнения (матрица Lmat: Lmat · vec(V) = 0) и список граней для плоскостности."""
    nV = len(Vu)
    rows = []
    nmatch = 0
    for (Rc, tc), on, Y in faces:
        # грань F_γ = on; грань F_{γ⁻¹}: та, чей сосед — γ⁻¹X. Соответствие: γ v (v ∈ F_{γ⁻¹}) ∈ F_γ.
        pass
    # проще: для каждой грани F_γ и каждой её вершины w: w' = γ⁻¹ w должна быть вершиной P
    for (Rc, tc), on, Y in faces:
        Rinv = Rc.T
        for j in on:
            w = Vu[j]
            wp = Rinv @ (w - tc)
            d = np.linalg.norm(Vu - wp, axis=1)
            i = int(np.argmin(d))
            if d[i] > 1e-7:
                raise RuntimeError("vertex correspondence not found")
            # уравнение: Rc v_i + tc − v_j = 0 → однородно по (V, 1): добавим столбец свободных членов
            for c in range(3):
                row = np.zeros(3 * nV + 1)
                row[3 * i: 3 * i + 3] += Rc[c]
                row[3 * j + c] -= 1
                row[-1] = tc[c]
                rows.append(row)
            nmatch += 1
    return np.array(rows), nmatch


def planarity_jacobian(Vu, faces):
    nV = len(Vu)
    rows = []
    for _, on, _ in faces:
        if len(on) <= 3:
            continue
        w1, w2, w3 = (Vu[on[0]], Vu[on[1]], Vu[on[2]])
        for k in on[3:]:
            wk = Vu[k]
            # f = ((w2−w1)×(w3−w1))·(wk−w1)
            a, b, c = w2 - w1, w3 - w1, wk - w1
            n = np.cross(a, b)
            row = np.zeros(3 * nV + 1)
            # производные
            d_wk = n
            d_w2 = np.cross(b, c)          # ∂/∂a of (a×b)·c = b×c
            d_w3 = np.cross(c, a)          # ∂/∂b = c×a
            d_w1 = -(d_wk + d_w2 + d_w3)
            for idx, dv in ((on[0], d_w1), (on[1], d_w2), (on[2], d_w3), (k, d_wk)):
                row[3 * idx: 3 * idx + 3] += dv
            rows.append(row)
    return np.array(rows)


def analyse(hm, x, c=1.0):
    G = Group(hm)
    E = Evaluator(G, c)
    X, Vu, faces = dv_structure(E, x)
    L, nmatch = equations(Vu, faces)
    Pj = planarity_jacobian(Vu, faces)
    J = np.vstack([L[:, :-1], Pj[:, :-1]]) if len(Pj) else L[:, :-1]
    # проверка: DV-вершины удовлетворяют линейным уравнениям
    resid = np.abs(L[:, :-1] @ Vu.ravel() + L[:, -1]).max()
    s = np.linalg.svd(J, compute_uv=False)
    rank = int(np.sum(s > 1e-8 * s[0]))
    null = J.shape[1] - rank
    print(f"{hm} c={c}: facets {len(faces)}, vertices {len(Vu)}, unknowns {J.shape[1]}, "
          f"linear eqs {len(L)} (residual {resid:.1e}), planarity eqs {len(Pj)}")
    print(f"  tangent dimension of the tiling family = {null}; minus 3 (moving x) → non-DV directions: {null - 3}")
    gap = s[rank - 1] / s[0] if rank > 0 else 0
    print(f"  smallest kept / next singular values: {s[rank-1]:.2e} / {s[rank]:.2e}" if rank < len(s) else "")
    return null


if __name__ == "__main__":
    analyse("I 41 3 2", (427 / 6984, 761 / 6984, 1421 / 6984))
    analyse("I 41 2 2", (62 / 125, 41 / 125, 79 / 1000), 727 / 500)
