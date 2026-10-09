"""Деформации стереоэдров: пространство изоэдральных разбиений вблизи DV-ячейки.

Плитка P с группой Γ (тривиальный стабилизатор) даёт разбиение ⇔ копии γP не перекрываются и vol P = V₀
(объём фундаментальной области). Для каждого соседа γ (пары {γ, γ⁻¹}) задаём разделяющую плоскость
π_γ = {n·z = d}, P ⊂ {n·z ≤ d}; тогда автоматически P ⊂ γ⁻¹{n·z ≥ d} = {(−Rᵀn)·z ≤ n·t − d}, и P, γP
разделены. Разбиения — это точки, где vol P достигает максимума V₀; DV-ячейка — одна из них
(n = γx − x, d = (|γx|² − |x|²)/2).

Здесь считаем гессиан vol по параметрам плоскостей (n, d) в DV-точке. Ядро гессиана — касательные
направления, вдоль которых плитка остаётся разбиением (до второго порядка). Тривиальные направления:
масштаб каждой пары (n, d) → λ(n, d) и сдвиг точки x (3 направления). Если ядро больше — рядом есть
не-DV стереоэдры.
"""
import itertools
import sys

import numpy as np
from scipy.spatial import HalfspaceIntersection, ConvexHull

from dvfast import Group, Evaluator


def dv_pairs(E, x):
    """Соседи-грани DV-ячейки: список пар (R, t, n, d) в декартовых координатах (кубические/тетр. решётки)."""
    G = E.G
    B = E.B
    X = B @ x
    nf, _, _ = E.facets(x)
    pairs = []
    seen = []
    for (R, t) in zip(G.R, G.t):
        for s in E.shifts:
            y = R @ x + t + s
            Y = B @ y
            if np.linalg.norm(Y - X) < 1e-9 or np.linalg.norm(Y - X) > E.radius:
                continue
            pairs.append((R, t + s, Y))
    # оставляем только тех, чьи биссекторы — грани
    _, A, b, V = E.cell_of(x)
    keep = []
    for (R, ts, Y) in pairs:
        n = Y - X
        d = (Y @ Y - X @ X) / 2
        on = np.abs(V @ n - d) / np.linalg.norm(n) < 1e-9
        if on.sum() >= 3 and np.linalg.matrix_rank(V[on] - V[on][0], tol=1e-8) >= 2:
            keep.append((R, ts, n, d))
    # группируем по парам {γ, γ⁻¹}: оставляем одного представителя
    reps = []
    for (R, ts, n, d) in keep:
        Rinv = R.T  # ортогональные (кубическая/тетрагональная решётка с B = diag)
        tinv = -Rinv @ ts
        dup = False
        for (R2, t2, _, _) in reps:
            if np.allclose(R2, Rinv) and np.allclose(t2, tinv):
                dup = True
                break
        if not dup:
            reps.append((R, ts, n, d))
    return reps, nf


def involution_frame(R, ts, B):
    """Для инволюции γ (γ² = id) — базис Q допустимых нормалей (⊥ неподвижному подпространству) и точка c."""
    Rc = B @ R @ np.linalg.inv(B)
    tc = B @ ts
    if not np.allclose(Rc @ (Rc @ np.zeros(3) + tc) + tc, 0, atol=1e-9) or not np.allclose(Rc @ Rc, np.eye(3)):
        return None
    c = np.linalg.lstsq(Rc - np.eye(3), -tc, rcond=None)[0]
    w, V = np.linalg.eig(Rc)
    U = np.real(V[:, np.abs(w - 1) < 1e-9])  # неподвижные направления
    if U.size == 0:
        Q = np.eye(3)
    else:
        Uo, _ = np.linalg.qr(U)
        P = np.eye(3) - Uo @ Uo.T
        Q = np.linalg.svd(P)[0][:, :3 - Uo.shape[1]]
    return Q, c


def make_param(reps, B):
    """Возвращает (p0, decode): p0 — плоский вектор параметров, decode(v) -> список (n, d)."""
    blocks = []
    p0 = []
    for (R, ts, n, d) in reps:
        fr = involution_frame(R, ts, B)
        if fr is None:
            blocks.append(("gen", len(p0), 4, None))
            p0 += list(n) + [d]
        else:
            Q, c = fr
            a = Q.T @ n
            blocks.append(("inv", len(p0), Q.shape[1], (Q, c)))
            p0 += list(a)
    def decode(v):
        out = []
        for kind, i, m, data in blocks:
            if kind == "gen":
                out.append((v[i:i + 3], v[i + 3]))
            else:
                Q, c = data
                n = Q @ v[i:i + m]
                out.append((n, n @ c))
        return out
    return np.array(p0, dtype=float), decode, blocks


def halfspaces(reps, params, B):
    """Полупространства A z ≤ b для параметров плоскостей. params: (k, 4) — (n, d) для каждого представителя."""
    rows = []
    for (R, ts, _, _), (n, d) in zip(reps, params):
        rows.append(np.append(n, -d))
        # γ в декартовых координатах: z ↦ Rc z + tc, где Rc = B R B⁻¹, tc = B ts
        Rc = B @ R @ np.linalg.inv(B)
        tc = B @ ts
        m = -Rc.T @ n
        rows.append(np.append(m, -(n @ tc - d)))
    return np.array(rows)


def volume(reps, params, B, X):
    hs = halfspaces(reps, params, B)
    try:
        H = HalfspaceIntersection(hs, X)
        return ConvexHull(H.intersections).volume
    except Exception:
        return np.nan


def analyse(hm, x, c=1.0, h=1e-4):
    G = Group(hm)
    E = Evaluator(G, c)
    x = np.array(x, dtype=float)
    reps, nf = dv_pairs(E, x)
    X = E.B @ x
    flat, decode, blocks = make_param(reps, E.B)
    k = len(reps)
    n_inv = sum(1 for b in blocks if b[0] == "inv")
    V0 = volume(reps, decode(flat), E.B, X)
    Vfd = abs(np.linalg.det(E.B)) / len(G.R)
    m = len(flat)

    def f(v):
        return volume(reps, decode(v), E.B, X)
    # гессиан конечными разностями
    Hs = np.zeros((m, m))
    f0 = f(flat)
    for i in range(m):
        for j in range(i, m):
            e_i = np.zeros(m); e_i[i] = h
            e_j = np.zeros(m); e_j[j] = h
            val = (f(flat + e_i + e_j) - f(flat + e_i - e_j) - f(flat - e_i + e_j) + f(flat - e_i - e_j)) / (4 * h * h)
            Hs[i, j] = Hs[j, i] = val
    g = np.array([(f(flat + np.eye(m)[i] * h) - f(flat - np.eye(m)[i] * h)) / (2 * h) for i in range(m)])
    w = np.linalg.eigvalsh(Hs)
    scale = np.max(np.abs(w))
    null = np.sum(np.abs(w) < 1e-6 * scale)
    print(f"{hm} x={np.round(x, 5)} c={c}: facets {nf}, pairs {k} (involutions {n_inv}), params {m}")
    print(f"  vol(P) = {V0:.10f}, fundamental domain = {Vfd:.10f}  (tiling: {abs(V0 - Vfd) < 1e-9})")
    print(f"  |grad| = {np.linalg.norm(g):.2e}  (≈0 at a maximum)")
    print(f"  Hessian eigenvalues: max {w.max():.2e}, min {w.min():.2e}; near-zero: {null}")
    print(f"  expected trivial null directions: {k} (scales) + 3 (moving x) = {k + 3}")
    print(f"  → extra tiling directions (non-DV deformations, 2nd order): {null - k - 3}")
    return w, null, k


if __name__ == "__main__":
    analyse("I 41 3 2", (427 / 6984, 761 / 6984, 1421 / 6984))
