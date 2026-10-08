"""Обратное утверждение (общий случай): суммы отрезков и треугольников суммарной размерности ≥ 5
почти всегда не замощают ℝ³.

Для члена p семейства и ребра e: если ни одно соотношение Σ kᵢθᵢ = π или 2π с kₑ ≥ 1 не выполнено
тождественно на окрестности p в камере, то множество членов камеры, проходящих фильтр для e, —
счётное объединение собственных аналитических подмножеств, то есть множество меры ноль.
Тождество Σ kᵢθᵢ ≡ const ⇒ Σ kᵢ ∇θᵢ ≡ 0. Поэтому ребро e — «свидетель», если единичный вектор eₑ
лежит в пространстве строк матрицы градиентов (по точкам камеры): тогда любая зависимость
градиентов имеет kₑ = 0. Градиенты — конечные разности в 60-значной арифметике; точки — малые
сдвиги p внутри той же камеры (комбинаторика проверяется).

Запуск: python converse.py [--members 20]
"""
import argparse
import itertools
import json
import pathlib
import sys

import mpmath as mp
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from tilecheck import Tile  # noqa: E402

mp.mp.dps = 60


def rot(q):
    """Поворот из (ненормированного) кватерниона q (mpmath)."""
    n = mp.sqrt(sum(x * x for x in q))
    w, x, y, z = (t / n for t in q)
    return mp.matrix([[1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
                      [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
                      [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)]])


BASE = {
    "seg": [[0, 0, 0], [1, 0, 0]],
    "tri": [[0, 0, 0], [1, 0, 0], [mp.mpf(1) / 2, mp.sqrt(3) / 2, 0]],
}


def member(kind, params):
    """params — плоский список: по кватерниону на каждое слагаемое, кроме первого (оно фиксировано)."""
    parts = kind.split("+")
    out = []
    for i, nm in enumerate(parts):
        V = [mp.matrix(v) for v in BASE[nm]]
        if i > 0:
            R = rot(params[4 * (i - 1): 4 * i])
            V = [R * v for v in V]
        out.append(V)
    return out


def labels_and_faces(parts):
    """Комбинаторика по float-оболочке: вершины-метки и грани (метки в циклическом порядке)."""
    labs = list(itertools.product(*[range(len(P)) for P in parts]))
    pts = np.array([[float(sum(parts[k][lab[k]][j] for k in range(len(parts)))) for j in range(3)] for lab in labs])
    T = Tile(pts)
    vlab = [labs[int(np.argmin(np.linalg.norm(pts - v, axis=1)))] for v in T.vertices]
    faces = []
    for f in T.faces():
        Q = T.vertices[f]
        cen = Q.mean(axis=0)
        nrm = np.linalg.svd(Q - cen)[2][-1]
        u = (Q[0] - cen) / np.linalg.norm(Q[0] - cen)
        w = np.cross(nrm, u)
        order = np.argsort(np.arctan2((Q - cen) @ w, (Q - cen) @ u))
        faces.append(tuple(vlab[f[i]] for i in order))
    return tuple(sorted(vlab)), frozenset(frozenset(F) for F in faces), faces


def angles(kind, params, faces_ordered):
    parts = member(kind, params)
    pt = lambda lab: sum((parts[k][lab[k]] for k in range(len(parts))), mp.matrix([0, 0, 0]))
    allv = {lab for F in faces_ordered for lab in F}
    cen = sum((pt(l) for l in allv), mp.matrix([0, 0, 0])) / len(allv)
    normals = []
    for F in faces_ordered:
        p0, p1, p2 = pt(F[0]), pt(F[1]), pt(F[2])
        n = mp.matrix([(p1 - p0)[1] * (p2 - p0)[2] - (p1 - p0)[2] * (p2 - p0)[1],
                       (p1 - p0)[2] * (p2 - p0)[0] - (p1 - p0)[0] * (p2 - p0)[2],
                       (p1 - p0)[0] * (p2 - p0)[1] - (p1 - p0)[1] * (p2 - p0)[0]])
        if sum(n[i] * (cen - p0)[i] for i in range(3)) > 0:
            n = -n
        normals.append(n / mp.norm(n))
    edges = {}
    for fi, F in enumerate(faces_ordered):
        for i in range(len(F)):
            e = frozenset((F[i], F[(i + 1) % len(F)]))
            edges.setdefault(e, []).append(fi)
    keys = sorted(edges, key=lambda e: sorted(e))
    th = []
    for e in keys:
        f1, f2 = edges[e]
        cosv = sum(normals[f1][i] * normals[f2][i] for i in range(3))
        th.append(mp.pi - mp.acos(max(-1, min(1, cosv))))
    return keys, th


def analyse(kind, rng, npts=6, p0=None):
    nparams = 4 * (len(kind.split("+")) - 1)
    if p0 is None:
        p0 = [mp.mpf(float(x)) for x in rng.normal(size=nparams)]
    vk, fk, faces = labels_and_faces(member(kind, p0))
    keys, th0 = angles(kind, p0, faces)
    # классы рёбер с тождественно равными углами: сравним значения в нескольких точках
    pts = [p0] + [[x + mp.mpf(float(d)) for x, d in zip(p0, rng.normal(scale=1e-3, size=nparams))] for _ in range(npts)]
    vals, grads = [], []
    h = mp.mpf(10) ** -25
    for p in pts:
        vk2, fk2, faces2 = labels_and_faces(member(kind, p))
        if (vk2, fk2) != (vk, fk):
            return None  # вышли из камеры — пропускаем член
        _, th = angles(kind, p, faces)
        vals.append(th)
        G = []
        for j in range(nparams):
            q = list(p)
            q[j] += h
            _, thp = angles(kind, q, faces)
            G.append([(a - b) / h for a, b in zip(thp, th)])
        grads.append(G)
    E = len(keys)
    classes = []
    for e in range(E):
        for cl in classes:
            r = cl[0]
            if all(abs(v[e] - v[r]) < mp.mpf(10) ** -40 for v in vals):
                cl.append(e)
                break
        else:
            classes.append([e])
    C = len(classes)
    # матрица: строки — (точка, параметр), столбцы — классы
    M = mp.matrix(len(pts) * nparams, C)
    r = 0
    for G in grads:
        for j in range(nparams):
            for ci, cl in enumerate(classes):
                M[r, ci] = G[j][cl[0]]
            r += 1
    # тождества фильтра: k ≥ 0 целые, kₑ ≥ 1, Σ kᵢθᵢ ≡ π или 2π на камере.
    # Перебор по значениям в p0 (полон: Σk ≤ 2π/θ_min), затем проверка градиента (тождество ⇒ Σk∇θ = 0).
    th0 = [vals[0][cl[0]] for cl in classes]
    order = sorted(range(C), key=lambda i: -th0[i])
    A = [th0[i] for i in order]
    tol_v = mp.mpf(10) ** -35
    found = []

    def rec(i, rest, k):
        if len(found) > 20000:
            return
        for T in (mp.pi, 2 * mp.pi):
            pass
        if i == C:
            return
        # можно взять 0..m копий A[i]
        m = int(mp.floor((rest + tol_v) / A[i]))
        for t in range(m, -1, -1):
            k[i] = t
            r2 = rest - t * A[i]
            if t > 0:
                for T in (mp.pi, 2 * mp.pi):
                    used = 2 * mp.pi - r2  # Σ уже выбранных
                    if abs(used - T) < tol_v:
                        found.append((list(k[:i + 1]) + [0] * (C - i - 1), T))
            rec(i + 1, r2, k)
        k[i] = 0

    rec(0, 2 * mp.pi, [0] * C)
    covered = set()
    identities = 0
    for k, T in found:
        kk = [0] * C
        for pos, cnt in enumerate(k):
            kk[order[pos]] = cnt
        res = max(abs(sum(kk[ci] * M[r, ci] for ci in range(C))) for r in range(M.rows))
        if res < mp.mpf(10) ** -20:
            identities += 1
            covered |= {ci for ci in range(C) if kk[ci]}
    witnesses = [ci for ci in range(C) if ci not in covered]
    thmin = min(th0)
    return {"faces": len(fk), "vertices": len(vk), "edgeClasses": C, "paramDim": nparams,
            "valueMatches": len(found), "identities": identities, "witnesses": len(witnesses),
            "minAngleDeg": float(mp.degrees(thmin)), "combinatorics": hash((vk, fk))}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--members", type=int, default=10)
    ap.add_argument("--kinds", default="tri+tri+seg,tri+seg+seg+seg,seg+seg+seg+seg+seg,tri+tri,tri+seg+seg")
    args = ap.parse_args()
    rng = np.random.default_rng(5)
    out = {}
    for kind in args.kinds.split(","):
        rows = []
        while len(rows) < args.members:
            r = analyse(kind, rng)
            if r:
                rows.append(r)
                print(kind, r, flush=True)
        out[kind] = rows
    (ROOT / "data" / "converse.json").write_text(json.dumps(out, indent=1, default=str))


if __name__ == "__main__":
    main()
