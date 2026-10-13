"""Стержни: замощения тетраэдром бесконечной призмы копиями g^k T (g — изометрия, склеивающая T с gT по грани).

Если орбита {g^k T} (для несобственного g — орбита g² с мотивом из двух копий) заполняет бесконечную
выпуклую призму без щелей, то сечение призмы — выпуклый многоугольник C. Если C замощает плоскость
(любой треугольник и четырёхугольник), то стержни, уложенные по замощению плоскости и сдвинутые вдоль осей
как угодно, дают замощение пространства тетраэдром — как правило, не «лицом к лицу». Такие разбиения
не видны перебору конечных выпуклых кластеров (clusters.py): у куска стержня концы не плоские.

Критерий стержня: площадь проекции орбиты на плоскость ⊥ оси × шаг вдоль оси = (копий на шаг) × объём T.
Запуск: python rods.py  → data/rods.json
"""
import itertools
import json
import pathlib

import numpy as np
from scipy.spatial import ConvexHull

from clusters import tet_float, place_on_face, congruences
from tetra import PRINTED_A, EXCLUDED_BY_LP, SPORADIC_A

HERE = pathlib.Path(__file__).resolve().parent


def affine(src, dst):
    """Аффинное отображение, переводящее 4 точки src в dst: (R, t)."""
    A = np.c_[src, np.ones(4)]
    X = np.linalg.solve(A, dst)  # 4×3
    return X[:3].T, X[3]


def rod_test(V, R, t, K=12):
    vol = abs(np.linalg.det(np.array([V[1] - V[0], V[2] - V[0], V[3] - V[0]]))) / 6
    proper = np.linalg.det(R) > 0
    if proper:
        Rh, th, m = R, t, 1
    else:
        Rh, th, m = R @ R, R @ t + t, 2
    # ось винта: собственный вектор Rh с собственным числом 1
    w, U = np.linalg.eig(Rh)
    k = int(np.argmin(np.abs(w - 1)))
    d = np.real(U[:, k]); d /= np.linalg.norm(d)
    tau = abs(th @ d)
    if tau < 1e-7:
        return None
    # орбита
    pts = []
    Rk, tk = np.eye(3), np.zeros(3)
    Rinv, tinv = np.linalg.inv(R), -np.linalg.inv(R) @ t
    for direction in (1, -1):
        Rk, tk = np.eye(3), np.zeros(3)
        for _ in range(K):
            pts.append(V @ Rk.T + tk)
            if direction == 1:
                Rk, tk = R @ Rk, R @ tk + t
            else:
                Rk, tk = Rinv @ Rk, Rinv @ tk + tinv
    P = np.vstack(pts)
    # проекция на плоскость ⊥ d
    e1 = np.cross(d, [1, 0, 0]) if abs(d[0]) < 0.9 else np.cross(d, [0, 1, 0])
    e1 /= np.linalg.norm(e1); e2 = np.cross(d, e1)
    Q = np.c_[P @ e1, P @ e2]
    try:
        h = ConvexHull(Q)
    except Exception:
        return None
    ratio = h.volume * tau / (m * vol)
    return {"ratio": float(ratio), "tau": float(tau), "proper": bool(proper), "polygon": Q[h.vertices].tolist(),
            "n_sides": int(len(h.vertices))}


def simplify_polygon(poly, tol=1e-6):
    P = np.array(poly)
    keep = []
    n = len(P)
    for i in range(n):
        a, b, c = P[i - 1], P[i], P[(i + 1) % n]
        u, w = b - a, c - b
        if abs(u[0] * w[1] - u[1] * w[0]) > tol * max(np.linalg.norm(u), 1):
            keep.append(b)
    return np.array(keep)


def rods_for(a):
    T, V = tet_float(a)
    found = []
    for fa, fb in itertools.product(range(4), repeat=2):
        A = [v for v in range(4) if v != fa]  # грань fa (противоположна вершине fa)
        Bf = [v for v in range(4) if v != fb]
        for perm in congruences(V, Bf, V[A]):
            # g: вершины грани Bf → вершины грани A в порядке perm; противоположная вершина fb → по другую сторону от fa
            dst_face = V[A][list(perm)]
            new_apex = place_on_face(None, dst_face, V[fa], V[Bf], V[fb])
            src = np.vstack([V[Bf], V[fb]])
            dst = np.vstack([dst_face, new_apex])
            R, t = affine(src, dst)
            if np.abs(R @ R.T - np.eye(3)).max() > 1e-6:
                continue
            r = rod_test(V, R, t)
            if r and abs(r["ratio"] - 1) < 1e-6:
                C = simplify_polygon(r["polygon"])
                found.append({"faces": [fa + 1, fb + 1], "perm": list(perm), "proper": r["proper"], "tau": r["tau"],
                              "section_sides": int(len(C))})
    return found


def main():
    labels = {i + 1: a for i, a in enumerate(PRINTED_A)}
    out = {}
    for i, a in labels.items():
        tag = "excluded_by_LP" if a in EXCLUDED_BY_LP else "A"
        f = rods_for(a)
        out[i] = {"angles": [str(x) for x in a], "set": tag, "rods": f}
        sides = sorted({r["section_sides"] for r in f})
        print(f"#{i:2d} {tag:15s} rods: {len(f):3d}  section sides: {sides}", flush=True)
    (HERE / "data" / "rods.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()


# ---------- стержни с мотивом из нескольких копий ----------
from clusters import overlap  # noqa: E402


def glue_all(V0, tets):
    """Все способы приклеить новую копию (вершины в порядке V0) к свободной грани кластера, без перекрытий."""
    out = []
    for ti, T in enumerate(tets):
        for f in range(4):
            face = [v for v in range(4) if v != f]
            for g in range(4):  # какая грань новой копии ляжет на face
                src = [v for v in range(4) if v != g]
                for perm in congruences(V0, src, T[face]):
                    dst_face = T[face][list(perm)]
                    apex = place_on_face(None, dst_face, T[f], V0[src], V0[g])
                    new = np.empty((4, 3))
                    for s, d in zip(src, dst_face):
                        new[s] = d
                    new[g] = apex
                    if any(overlap(new, U) for U in tets):
                        continue
                    out.append(new)
    return out


def boundary_faces(tets):
    """Свободные грани кластера: (индекс копии, номер противоположной вершины)."""
    faces = []
    for ti, T in enumerate(tets):
        for f in range(4):
            pts = T[[v for v in range(4) if v != f]]
            shared = False
            for tj, U in enumerate(tets):
                if tj == ti:
                    continue
                for h in range(4):
                    q = U[[v for v in range(4) if v != h]]
                    if all(min(np.linalg.norm(p - x) for x in q) < 1e-7 for p in pts):
                        shared = True
            if not shared:
                faces.append((ti, f))
    return faces


def rods_motif(a, k_max=3):
    T0, V = tet_float(a)
    clusters = [[V.copy()]]
    level = [[V.copy()]]
    for k in range(2, k_max + 1):
        nxt = []
        for cl in level:
            for new in glue_all(V, cl):
                nxt.append(cl + [new])
        level = nxt
        clusters += nxt
    found = []
    for cl in clusters[1:]:
        m = len(cl)
        P = np.vstack(cl)
        bf = boundary_faces(cl)
        for (ta, fa), (tb, fb) in itertools.product(bf, repeat=2):
            if (ta, fa) == (tb, fb):
                continue
            A = cl[ta][[v for v in range(4) if v != fa]]
            Bv = cl[tb][[v for v in range(4) if v != fb]]
            for perm in itertools.permutations(range(3)):
                dst = A[list(perm)]
                if any(abs(np.linalg.norm(Bv[i] - Bv[j]) - np.linalg.norm(dst[i] - dst[j])) > 1e-7 for i, j in itertools.combinations(range(3), 2)):
                    continue
                apex = place_on_face(None, dst, cl[ta][fa], Bv, cl[tb][fb])
                R, t = affine(np.vstack([Bv, cl[tb][fb]]), np.vstack([dst, apex]))
                if np.abs(R @ R.T - np.eye(3)).max() > 1e-6:
                    continue
                img = [U @ R.T + t for U in cl]
                if any(overlap(x, y) for x in img for y in cl):
                    continue
                r = rod_test_cluster(cl, R, t)
                if r and abs(r - 1) < 1e-6:
                    found.append({"motif": m, "faces": [(int(ta), int(fa) + 1), (int(tb), int(fb) + 1)], "proper": bool(np.linalg.det(R) > 0)})
    return found


def rod_test_cluster(cl, R, t, K=10):
    m = len(cl)
    V = np.vstack(cl)
    proper = np.linalg.det(R) > 0
    Rh, th, mm = (R, t, m) if proper else (R @ R, R @ t + t, 2 * m)
    w, U = np.linalg.eig(Rh)
    k = int(np.argmin(np.abs(w - 1)))
    d = np.real(U[:, k]); d /= np.linalg.norm(d)
    tau = abs(th @ d)
    if tau < 1e-7:
        return None
    vol = abs(np.linalg.det(np.array([cl[0][1] - cl[0][0], cl[0][2] - cl[0][0], cl[0][3] - cl[0][0]]))) / 6
    pts = []
    Rinv, tinv = np.linalg.inv(R), -np.linalg.inv(R) @ t
    for direction in (1, -1):
        Rk, tk = np.eye(3), np.zeros(3)
        for _ in range(K):
            pts.append(V @ Rk.T + tk)
            Rk, tk = (R @ Rk, R @ tk + t) if direction == 1 else (Rinv @ Rk, Rinv @ tk + tinv)
    P = np.vstack(pts)
    e1 = np.cross(d, [1, 0, 0]) if abs(d[0]) < 0.9 else np.cross(d, [0, 1, 0])
    e1 /= np.linalg.norm(e1); e2 = np.cross(d, e1)
    try:
        h = ConvexHull(np.c_[P @ e1, P @ e2])
    except Exception:
        return None
    return h.volume * tau / (mm * vol)


def main_motif(k_max=3):
    out = {}
    for i, a in enumerate(PRINTED_A, start=1):
        f = rods_motif(a, k_max)
        out[i] = f
        print(f"#{i:2d} {'excl' if a in EXCLUDED_BY_LP else 'A   '} motif rods (k≤{k_max}): {len(f)}", flush=True)
    (HERE / "data" / f"rods_motif_k{k_max}.json").write_text(json.dumps(out, indent=1))
