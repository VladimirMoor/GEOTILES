"""Точная размерность пространства изоэдральных деформаций DV-ячейки (рациональная арифметика).

Всё во фракционных координатах: группа действует рационально (z ↦ R z + t), уравнения спаривания граней
и плоскостности аффинно-инвариантны, поэтому метрика нужна только для построения исходной DV-ячейки.

  * Точные вершины DV-ячейки — пересечения троек точных плоскостей (как в exact.py).
  * Уравнения: (L) γ⁻¹ w = v для каждой вершины w грани F_γ (v — вершина P); (P) плоскостность граней.
  * Ранг якобиана — по модулю двух больших простых (ранг над ℚ ≥ ранга mod p, равенство почти наверное;
    берём максимум по двум простым).
  * Направления «сдвиг x» (DV-семейство) — точная производная вершин по x; проверяем, что они лежат в ядре.
Ответ: dim ker J − rank(DV-направлений) = число не-DV направлений первого порядка.
"""
import itertools
import sys
from fractions import Fraction as Fr

import numpy as np

from dvfast import Group, Evaluator
from exact import metric, mv, dot, solve3

PRIMES = [2147483629, 2147483587]


def exact_cell(hm, x, c2):
    G = Group(hm)
    E = Evaluator(G, float(c2) ** 0.5)
    xf = np.array([float(t) for t in x])
    X, A, b, V = E.cell_of(xf)
    Binv = np.linalg.inv(E.B)
    M = metric(G, c2)
    R = [[[Fr(int(round(v * 24)), 24) for v in row] for row in Rm] for Rm in G.R]
    T = [[Fr(int(round(v * 24)), 24) for v in tv] for tv in G.t]
    Yf = (Binv @ (A + X).T).T
    planes = []  # (a, rhs, Rm, t(full, with shift), y)
    for yf in Yf:
        for Rm, tv in zip(R, T):
            base = [sum(Rm[i][j] * x[j] for j in range(3)) + tv[i] for i in range(3)]
            n = [round(yf[i] - float(base[i])) for i in range(3)]
            y = [base[i] + n[i] for i in range(3)]
            if max(abs(float(y[i]) - yf[i]) for i in range(3)) < 1e-9:
                tt = [tv[i] + n[i] for i in range(3)]
                a = [2 * s for s in mv(M, [y[i] - x[i] for i in range(3)])]
                rhs = dot(y, mv(M, y)) - dot(x, mv(M, x))
                planes.append((a, rhs, Rm, tt, y))
                break
    Vf = (Binv @ V.T).T
    Af = np.array([[float(t) for t in p[0]] for p in planes])
    bf = np.array([float(p[1]) for p in planes])
    verts = []
    vdef = []  # индексы трёх определяющих плоскостей
    for v in Vf:
        act = np.where(np.abs(Af @ v - bf) < 1e-7 * np.maximum(1, np.abs(bf)))[0]
        for i, j, k in itertools.combinations(act, 3):
            sol = solve3([planes[i][0], planes[j][0], planes[k][0]], [planes[i][1], planes[j][1], planes[k][1]])
            if sol is not None and all(dot(p[0], sol) <= p[1] for p in planes):
                if sol not in verts:
                    verts.append(sol)
                    vdef.append((i, j, k))
                break
    faces = []
    for pi, p in enumerate(planes):
        on = [vi for vi, v in enumerate(verts) if dot(p[0], v) == p[1]]
        if len(on) >= 3:
            faces.append((pi, on))
    return G, M, planes, verts, vdef, faces


def rank_mod(rows, ncols, p):
    A = np.array([[int((r.numerator % p) * pow(r.denominator % p, p - 2, p) % p) for r in row] for row in rows],
                 dtype=np.int64)
    A %= p
    rank = 0
    nr = A.shape[0]
    for c in range(ncols):
        piv = None
        for r in range(rank, nr):
            if A[r, c] != 0:
                piv = r
                break
        if piv is None:
            continue
        A[[rank, piv]] = A[[piv, rank]]
        inv = pow(int(A[rank, c]), p - 2, p)
        A[rank] = (A[rank] * inv) % p
        nz = np.nonzero(A[:, c])[0]
        for r in nz:
            if r != rank:
                A[r] = (A[r] - A[r, c] * A[rank]) % p
        rank += 1
    return rank


def analyse(hm, x, c2=Fr(1)):
    G, M, planes, verts, vdef, faces = exact_cell(hm, x, c2)
    nV = len(verts)
    ncols = 3 * nV
    rows = []
    # (L) спаривание: для грани F_γ (плоскость между x и y = γx) и её вершины w: γ⁻¹ w = R⁻¹(w − t) — вершина P
    idx = {tuple(v): i for i, v in enumerate(verts)}
    for pi, on in faces:
        _, _, Rm, tt, _ = planes[pi]
        Rinv = [[r for r in row] for row in np.array(Rm, dtype=object).T.tolist()]  # R ортогональна? для фракц. нет
        # общий обратный для рациональной матрицы
        Rf = np.array([[float(v) for v in row] for row in Rm])
        Rinv = [[Fr(round(v * 24), 24) for v in row] for row in np.linalg.inv(Rf)]
        for j in on:
            w = verts[j]
            wp = [sum(Rinv[a][b] * (w[b] - tt[b]) for b in range(3)) for a in range(3)]
            i = idx.get(tuple(wp))
            if i is None:
                raise RuntimeError("no vertex correspondence (degenerate cell?)")
            # уравнение R v_i + t − v_j = 0 (однородная часть по вариациям: R δv_i − δv_j = 0)
            for cc in range(3):
                row = [Fr(0)] * ncols
                for b in range(3):
                    row[3 * i + b] += Rm[cc][b]
                row[3 * j + cc] -= 1
                rows.append(row)
    nL = len(rows)
    # (P) плоскостность
    for pi, on in faces:
        if len(on) <= 3:
            continue
        w1, w2, w3 = verts[on[0]], verts[on[1]], verts[on[2]]
        for k in on[3:]:
            wk = verts[k]
            a = [w2[t] - w1[t] for t in range(3)]
            b = [w3[t] - w1[t] for t in range(3)]
            c = [wk[t] - w1[t] for t in range(3)]
            def cross(u, v):
                return [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]
            d_wk = cross(a, b); d_w2 = cross(b, c); d_w3 = cross(c, a)
            d_w1 = [-(d_wk[t] + d_w2[t] + d_w3[t]) for t in range(3)]
            row = [Fr(0)] * ncols
            for vi, dv in ((on[0], d_w1), (on[1], d_w2), (on[2], d_w3), (k, d_wk)):
                for t in range(3):
                    row[3 * vi + t] += dv[t]
            rows.append(row)
    rk = max(rank_mod(rows, ncols, p) for p in PRIMES)
    null = ncols - rk
    # DV-направления: производная вершин по x_k
    dirs = []
    for kx in range(3):
        e = [Fr(int(t == kx)) for t in range(3)]
        dv_all = []
        for v, (i, j, k) in zip(verts, vdef):
            Arows, rhs = [], []
            for pi in (i, j, k):
                a, r, Rm, tt, y = planes[pi]
                Re = [sum(Rm[s][t] * e[t] for t in range(3)) for s in range(3)]
                da = [2 * s for s in mv(M, [Re[t] - e[t] for t in range(3)])]
                dr = 2 * dot(y, mv(M, Re)) - 2 * dot(x, mv(M, e))
                Arows.append(a)
                rhs.append(dr - dot(da, v))
            dv_all += solve3(Arows, rhs)
        dirs.append(dv_all)
    # для некубических групп — ещё направление «изменение метрики» (c²): dM/dc² = e3 e3ᵀ
    if G.system != "cubic":
        D = [[Fr(0)] * 3 for _ in range(3)]; D[2][2] = Fr(1)
        dv_all = []
        for v, (i, j, k) in zip(verts, vdef):
            Arows, rhs = [], []
            for pi in (i, j, k):
                a, r, Rm, tt, y = planes[pi]
                da = [2 * s for s in mv(D, [y[t] - x[t] for t in range(3)])]
                dr = dot(y, mv(D, y)) - dot(x, mv(D, x))
                Arows.append(a); rhs.append(dr - dot(da, v))
            dv_all += solve3(Arows, rhs)
        dirs.append(dv_all)
    # проверка: DV-направления в ядре
    resid = max(abs(sum(rw[t] * d[t] for t in range(ncols))) for rw in rows for d in dirs)
    rk_dirs = max(rank_mod(dirs, ncols, p) for p in PRIMES)
    print(f"{hm} x={tuple(map(str, x))} c²={c2}: facets {len(faces)}, vertices {nV}, unknowns {ncols}")
    print(f"  equations: pairing {nL}, planarity {len(rows) - nL}; exact rank {rk}; kernel dim {null}")
    print(f"  DV directions (moving x{'' if G.system == 'cubic' else ', changing c/a'}): rank {rk_dirs}, max residual in J·d = {float(resid):.1e} (must be 0)")
    print(f"  → non-DV first-order deformation directions: {null - rk_dirs}")
    return null - rk_dirs


if __name__ == "__main__":
    analyse("I 41 3 2", (Fr(427, 6984), Fr(761, 6984), Fr(1421, 6984)))
    analyse("I 41 2 2", (Fr(62, 125), Fr(41, 125), Fr(79, 1000)), Fr(727, 500) ** 2)
