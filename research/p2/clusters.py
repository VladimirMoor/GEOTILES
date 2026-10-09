"""Кластеры: выпуклые многогранники Q, склеенные face-to-face из k копий T (зеркальные разрешены).

Если Q замощает пространство, то T тоже замощает. Внутренние триангуляции соседних граней Q могут не
совпадать (например, грань-параллелограмм разрезана разными диагоналями), и именно так возникают
плоскости разлома. Для 38 открытых тетраэдров face-to-face разбиений нет, поэтому любое разбиение
должно устроено примерно так.

Шаг 1 (здесь): перечисляем кластеры до K копий с точностью до конгруэнтности и отбираем выпуклые
(объём выпуклой оболочки равен k·vol T).
Шаг 2 (poincare_q.py): ищем изоэдральные разбиения Q (теорема Пуанкаре со спариванием граней Q).
"""
import itertools
import sys

import numpy as np
from scipy.spatial import ConvexHull

from tetra import EDGES, PRINTED_A, EXCLUDED_BY_LP, Tet

EPS = 1e-7


def tet_float(a):
    T = Tet(a)
    V = np.array([[float(T.V[i][c]) for c in range(3)] for i in (1, 2, 3, 4)])
    return T, V


def face_ids(k):
    return [v for v in range(4) if v != k]


def place_on_face(P, Q3, other, src_face_pts, src_opp):
    """Изометрия, переводящая треугольник src_face_pts в Q3 (по порядку), а src_opp — в сторону,
    противоположную точке other. Возвращает 4 вершины новой копии в порядке исходных индексов."""
    a, b, c = src_face_pts
    A, B, C = Q3
    def frame(a, b, c):
        u = b - a; w = c - a; n = np.cross(u, w)
        return np.column_stack([u, w, n])
    M1 = frame(a, b, c)
    coef = np.linalg.solve(M1, src_opp - a)
    u2 = B - A; w2 = C - A; n2 = np.cross(u2, w2)
    cand = [A + coef[0] * u2 + coef[1] * w2 + s * coef[2] * n2 for s in (1, -1)]
    side_other = np.dot(other - A, n2)
    return [p for p in cand if np.dot(p - A, n2) * side_other < 0][0]


def overlap(T1, T2):
    """Пересекаются ли внутренности двух тетраэдров (SAT)."""
    axes = []
    for T in (T1, T2):
        for k in range(4):
            f = [T[i] for i in range(4) if i != k]
            axes.append(np.cross(f[1] - f[0], f[2] - f[0]))
    e1 = [T1[j] - T1[i] for i, j in itertools.combinations(range(4), 2)]
    e2 = [T2[j] - T2[i] for i, j in itertools.combinations(range(4), 2)]
    for u in e1:
        for v in e2:
            axes.append(np.cross(u, v))
    for ax in axes:
        n = np.linalg.norm(ax)
        if n < 1e-9:
            continue
        ax = ax / n
        p1 = T1 @ ax; p2 = T2 @ ax
        if p1.max() <= p2.min() + EPS or p2.max() <= p1.min() + EPS:
            return False
    return True


def lengths(V):
    return {(i, j): np.linalg.norm(V[i] - V[j]) for i in range(4) for j in range(4) if i != j}


def congruences(V, face_src, face_dst_pts):
    """Биекции вершин грани face_src (индексы в V) на точки face_dst_pts, сохраняющие длины."""
    out = []
    for perm in itertools.permutations(range(3)):
        ok = True
        for x, y in itertools.combinations(range(3), 2):
            d1 = np.linalg.norm(V[face_src[x]] - V[face_src[y]])
            d2 = np.linalg.norm(face_dst_pts[perm[x]] - face_dst_pts[perm[y]])
            if abs(d1 - d2) > 1e-7:
                ok = False
                break
        if ok:
            out.append(perm)
    return out


def canon_key(tets):
    pts = np.vstack(tets)
    uniq = []
    for p in pts:
        if not any(np.linalg.norm(p - q) < 1e-6 for q in uniq):
            uniq.append(p)
    D = sorted(round(float(np.linalg.norm(p - q)), 5) for p, q in itertools.combinations(uniq, 2))
    return (len(tets), len(uniq), tuple(D))


def enumerate_clusters(V, K):
    vol = abs(np.linalg.det(np.array([V[1] - V[0], V[2] - V[0], V[3] - V[0]]))) / 6
    start = [V.copy()]
    seen = {canon_key(start)}
    layer = [start]
    convex = []
    for k in range(2, K + 1):
        nxt = []
        for cl in layer:
            # свободные грани
            faces = []
            for ti, Tt in enumerate(cl):
                for kk in range(4):
                    fpts = [Tt[i] for i in range(4) if i != kk]
                    shared = False
                    for tj, Tu in enumerate(cl):
                        if tj == ti:
                            continue
                        if all(any(np.linalg.norm(p - q) < 1e-6 for q in Tu) for p in fpts):
                            shared = True
                            break
                    if not shared:
                        faces.append((fpts, Tt[kk]))
            for fpts, other in faces:
                for ks in range(4):
                    src = face_ids(ks)
                    for perm in congruences(V, src, fpts):
                        dst = [fpts[perm[x]] for x in range(3)]
                        newopp = place_on_face(None, dst, other, [V[s] for s in src], V[ks])
                        newT = np.zeros((4, 3))
                        for x, s in enumerate(src):
                            newT[s] = dst[x]
                        newT[ks] = newopp
                        if any(overlap(newT, Tu) for Tu in cl):
                            continue
                        ncl = cl + [newT]
                        key = canon_key(ncl)
                        if key in seen:
                            continue
                        seen.add(key)
                        nxt.append(ncl)
                        pts = np.vstack(ncl)
                        try:
                            hv = ConvexHull(pts).volume
                        except Exception:
                            continue
                        if abs(hv - k * vol) < 1e-6 * k * vol:
                            convex.append(ncl)
        layer = nxt
        print(f"   k={k}: {len(nxt)} clusters, convex so far {len(convex)}", flush=True)
    return convex


if __name__ == "__main__":
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    which = [int(x) for x in sys.argv[2:]] or None
    for idx, a in enumerate(PRINTED_A, 1):
        if a in EXCLUDED_BY_LP or (which and idx not in which):
            continue
        T, V = tet_float(a)
        print(f"#{idx}", flush=True)
        cv = enumerate_clusters(V, K)
        print(f"#{idx}: convex clusters up to {K}: {len(cv)}", flush=True)
