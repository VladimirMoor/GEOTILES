"""Невыпуклые кластеры как фундаментальные области (расширенный поиск конструкций).

Кластер Q — объединение k копий T, склеенных face-to-face. Его граница — треугольники копий, не общие
для двух копий. Грани Q — максимальные плоские связные области; рёбра Q — отрезки, по которым
встречаются две грани Q; двугранный угол может быть больше π.

Отсечение по лемме о кластерах: каждое плохое ребро каждой копии должно лежать на ребре Q (не внутри Q
и не внутри грани Q).
Проверка Пуанкаре: спаривание граней Q изометриями и замыкание циклов рёбер с суммой 2π. Кандидаты
затем проверяются численно: строим образы Q в шаре и считаем покрытие случайных точек.
"""
import itertools
import sys
from fractions import Fraction as Fr

import numpy as np

from tetra import EDGES, PRINTED_A, EXCLUDED_BY_LP, Tet
from stars import Lengths, edge_star_exists
from clusters import tet_float, face_ids, congruences, place_on_face, overlap, canon_key
from poincare_q import affine_from, is_isometry, partitions

TOL = 1e-6


def all_clusters(V, K):
    start = [V.copy()]
    seen = {canon_key(start)}
    layer = [start]
    out = [start]
    for k in range(2, K + 1):
        nxt = []
        for cl in layer:
            for ti, Tt in enumerate(cl):
                for kk in range(4):
                    fpts = [Tt[i] for i in range(4) if i != kk]
                    if any(tj != ti and all(any(np.linalg.norm(p - q) < 1e-6 for q in Tu) for p in fpts)
                           for tj, Tu in enumerate(cl)):
                        continue
                    for ks in range(4):
                        src = face_ids(ks)
                        for perm in congruences(V, src, fpts):
                            dst = [fpts[perm[x]] for x in range(3)]
                            newT = np.zeros((4, 3))
                            for x, s in enumerate(src):
                                newT[s] = dst[x]
                            newT[ks] = place_on_face(None, dst, Tt[kk], [V[s] for s in src], V[ks])
                            if any(overlap(newT, Tu) for Tu in cl):
                                continue
                            ncl = cl + [newT]
                            key = canon_key(ncl)
                            if key in seen:
                                continue
                            seen.add(key)
                            nxt.append(ncl)
        layer = nxt
        out += nxt
    return out


class Boundary:
    """Граничный комплекс кластера: грани-многоугольники, рёбра, двугранные углы."""

    def __init__(self, cl):
        pts = []

        def vid(p):
            for i, q in enumerate(pts):
                if np.linalg.norm(p - q) < 1e-6:
                    return i
            pts.append(p)
            return len(pts) - 1
        tets = [[vid(p) for p in Tt] for Tt in cl]
        self.P = np.array(pts)
        cnt = {}
        tri_owner = {}
        for ti, t in enumerate(tets):
            for k in range(4):
                f = frozenset(t[i] for i in range(4) if i != k)
                cnt[f] = cnt.get(f, 0) + 1
                tri_owner.setdefault(f, (ti, t[k]))
        tris = [f for f, c in cnt.items() if c == 1]
        self.ok = True
        cen = self.P.mean(axis=0)
        # внешние нормали треугольников
        info = {}
        for f in tris:
            a, b, c = sorted(f)
            n = np.cross(self.P[b] - self.P[a], self.P[c] - self.P[a])
            n /= np.linalg.norm(n)
            opp = tri_owner[f][1]
            if np.dot(self.P[opp] - self.P[a], n) > 0:
                n = -n
            info[f] = (n, float(np.dot(n, self.P[a])))
        # слияние копланарных соседних треугольников
        edge_tris = {}
        for f in tris:
            for e in itertools.combinations(sorted(f), 2):
                edge_tris.setdefault(frozenset(e), []).append(f)
        if any(len(v) != 2 for v in edge_tris.values()):
            self.ok = False  # не многообразие
            return
        parent = {f: f for f in tris}

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x
        for e, (f1, f2) in edge_tris.items():
            n1, d1 = info[f1]
            n2, d2 = info[f2]
            if np.allclose(n1, n2, atol=1e-7) and abs(d1 - d2) < 1e-7:
                parent[find(f1)] = find(f2)
        groups = {}
        for f in tris:
            groups.setdefault(find(f), []).append(f)
        self.faces, self.normals = [], []
        tri_face = {}
        for gi, (_, fs) in enumerate(groups.items()):
            for f in fs:
                tri_face[f] = gi
            self.normals.append(info[fs[0]][0])
        # граничные отрезки каждой грани
        seg_face = {}
        for e, (f1, f2) in edge_tris.items():
            g1, g2 = tri_face[f1], tri_face[f2]
            if g1 != g2:
                seg_face[e] = (g1, g2)
        for gi in range(len(groups)):
            segs = [e for e, gg in seg_face.items() if gi in gg]
            # обход цикла (ориентация против часовой стрелки снаружи)
            adj = {}
            for e in segs:
                a, b = tuple(e)
                adj.setdefault(a, []).append(b)
                adj.setdefault(b, []).append(a)
            if any(len(v) != 2 for v in adj.values()):
                self.ok = False  # грань с дыркой или касанием
                return
            start = min(adj)
            cyc = [start]
            prev, cur = None, start
            while True:
                nx = [x for x in adj[cur] if x != prev]
                nxt = nx[0] if prev is not None else adj[cur][0]
                if nxt == start:
                    break
                cyc.append(nxt)
                prev, cur = cur, nxt
            n = self.normals[gi]
            area = sum(np.dot(np.cross(self.P[cyc[i]], self.P[cyc[(i + 1) % len(cyc)]]), n) for i in range(len(cyc)))
            if area < 0:
                cyc.reverse()
            # убираем вершины, где обе стороны граничат с одной и той же гранью и точка коллинеарна
            def neighbor(a, b):
                g1, g2 = seg_face[frozenset((a, b))]
                return g2 if g1 == gi else g1
            changed = True
            while changed and len(cyc) > 3:
                changed = False
                for i in range(len(cyc)):
                    a, b, c = cyc[i - 1], cyc[i], cyc[(i + 1) % len(cyc)]
                    u, w = self.P[b] - self.P[a], self.P[c] - self.P[b]
                    if np.linalg.norm(np.cross(u, w)) < 1e-9 and np.dot(u, w) > 0 and neighbor(a, b) == neighbor(b, c):
                        # слить: запомнить отрезок (a, c)
                        g = neighbor(a, b)
                        seg_face[frozenset((a, c))] = (gi, g)
                        cyc.pop(i)
                        changed = True
                        break
            self.faces.append(cyc)
        # согласуем удаления: вершина удаляется, если она удалена в обеих гранях — иначе оставим
        self.edges = {}
        for fi, f in enumerate(self.faces):
            for a, b in zip(f, f[1:] + f[:1]):
                self.edges.setdefault(frozenset((a, b)), []).append(fi)
        if any(len(v) != 2 for v in self.edges.values()):
            self.ok = False
            return
        self.dihedral = {}
        for e, (f1, f2) in self.edges.items():
            a, b = tuple(e)
            n1, n2 = self.normals[f1], self.normals[f2]
            cosv = float(np.clip(-np.dot(n1, n2), -1, 1))
            ang = np.arccos(cosv)
            # выпуклое или вогнутое ребро: проверяем, смотрит ли грань f2 «внутрь» относительно f1
            mid = (self.P[a] + self.P[b]) / 2
            c2 = self.P[[v for v in self.faces[f2] if v not in (a, b)]].mean(axis=0)
            if np.dot(c2 - mid, n1) > 1e-9:
                ang = 2 * np.pi - ang
            self.dihedral[e] = Fr(float(ang / np.pi)).limit_denominator(720)


def bad_edges_ok(cl, V, bad):
    """Все плохие рёбра всех копий лежат на рёбрах Q (а не внутри граней или внутри Q)."""
    B = Boundary(cl)
    if not B.ok:
        return False, B
    qedges = [(B.P[list(e)[0]], B.P[list(e)[1]]) for e in B.edges]
    for Tt in cl:
        for (i, j) in bad:
            p, q = Tt[i - 1], Tt[j - 1]
            on = False
            for a, b in qedges:
                d = b - a
                L = np.linalg.norm(d)
                if np.linalg.norm(np.cross(p - a, d)) / L < 1e-6 and np.linalg.norm(np.cross(q - a, d)) / L < 1e-6:
                    t1, t2 = np.dot(p - a, d) / L ** 2, np.dot(q - a, d) / L ** 2
                    if -1e-6 <= min(t1, t2) and max(t1, t2) <= 1 + 1e-6:
                        on = True
                        break
            if not on:
                return False, B
    return True, B


def pairings(B, fG, fF):
    G, F = B.faces[fG], B.faces[fF]
    m = len(G)
    if len(F) != m:
        return []
    out = []
    nG, nF = B.normals[fG], B.normals[fF]
    for s in range(m):
        for d in (1, -1):
            vmap = {G[i]: F[(s + d * i) % m] for i in range(m)}
            # три неколлинеарные вершины
            idx = None
            for a, b, c in itertools.combinations(range(m), 3):
                if np.linalg.norm(np.cross(B.P[G[b]] - B.P[G[a]], B.P[G[c]] - B.P[G[a]])) > 1e-6:
                    idx = (a, b, c)
                    break
            a, b, c = idx
            src = [B.P[G[a]], B.P[G[b]], B.P[G[c]], B.P[G[a]] - nG]
            dst = [B.P[vmap[G[a]]], B.P[vmap[G[b]]], B.P[vmap[G[c]]], B.P[vmap[G[a]]] + nF]
            M = affine_from(src, dst)
            if not is_isometry(M):
                continue
            if all(np.allclose((M @ np.append(B.P[g], 1))[:3], B.P[vmap[g]], atol=1e-6) for g in G):
                if fG == fF and not np.allclose(M @ M, np.eye(4), atol=1e-6):
                    continue
                out.append((M, vmap))
    return out


def edge_check(B, gam):
    for e0, (fa, fb) in B.edges.items():
        g = np.eye(4)
        edge, leave = e0, fa
        total = Fr(0)
        while True:
            total += B.dihedral[edge]
            M, vmap, src = gam[leave]
            g = g @ M
            inv = {v: k for k, v in vmap.items()}
            if not all(x in inv for x in edge):
                return False
            edge2 = frozenset(inv[x] for x in edge)
            if edge2 not in B.edges:
                return False
            f1, f2 = B.edges[edge2]
            nxt = f2 if f1 == src else f1
            edge, leave = edge2, nxt
            if total > 2:
                return False
            if total == 2:
                if np.allclose(g, np.eye(4), atol=1e-6) and edge == e0 and leave == fa:
                    break
                return False
            if np.allclose(g, np.eye(4), atol=1e-6):
                return False
    return True


def search(B, cap=100000):
    nf = len(B.faces)
    cache = {}

    def opts(a, b):
        if (a, b) not in cache:
            cache[(a, b)] = pairings(B, a, b)
        return cache[(a, b)]
    count = 0
    for part in partitions(nf):
        choices = []
        for blk in part:
            o = ([(blk[0], blk[0], p) for p in opts(blk[0], blk[0])] if len(blk) == 1
                 else [(blk[0], blk[1], p) for p in opts(blk[1], blk[0])])
            if not o:
                break
            choices.append(o)
        else:
            for combo in itertools.product(*choices):
                count += 1
                if count > cap:
                    return None
                gam = {}
                for f, gsrc, (M, vmap) in combo:
                    gam[f] = (M, vmap, gsrc)
                    if gsrc != f:
                        gam[gsrc] = (np.linalg.inv(M), {v: k for k, v in vmap.items()}, f)
                if edge_check(B, gam):
                    return gam
    return False


def coverage(B, gam, cl, R=None, samples=400, rng=np.random.default_rng(0)):
    """Численная проверка: образы Q под группой, порождённой спариваниями, покрывают шар ровно один раз."""
    gens = [M for (M, _, _) in gam.values()]
    cen = B.P.mean(axis=0)
    rad = max(np.linalg.norm(p - cen) for p in B.P)
    R = R or 2.5 * rad
    elems = [np.eye(4)]
    seen = [tuple(np.round(np.eye(4), 5).ravel())]
    frontier = [np.eye(4)]
    while frontier:
        nf = []
        for g in frontier:
            for M in gens:
                h = g @ M
                c = (h @ np.append(cen, 1))[:3]
                if np.linalg.norm(c - cen) > R + 2 * rad:
                    continue
                key = tuple(np.round(h, 5).ravel())
                if key in seen:
                    continue
                seen.append(key)
                elems.append(h)
                nf.append(h)
        frontier = nf
        if len(elems) > 4000:
            break
    # точки: в шаре радиуса R/2 вокруг центра; проверка принадлежности через копии тетраэдров
    def inside_tet(x, Tt):
        A = np.column_stack([Tt[1] - Tt[0], Tt[2] - Tt[0], Tt[3] - Tt[0]])
        l = np.linalg.solve(A, x - Tt[0])
        return l.min() > 1e-9 and l.sum() < 1 - 1e-9
    bad = 0
    for _ in range(samples):
        x = cen + (rng.random(3) - 0.5) * R
        hits = 0
        for h in elems:
            for Tt in cl:
                Ti = np.array([(h @ np.append(p, 1))[:3] for p in Tt])
                if inside_tet(x, Ti):
                    hits += 1
        if hits != 1:
            bad += 1
    return bad, len(elems)


if __name__ == "__main__":
    K = int(sys.argv[1])
    which = [int(x) for x in sys.argv[2:]] or None
    for idx, a in enumerate(PRINTED_A, 1):
        if a in EXCLUDED_BY_LP or (which and idx not in which):
            continue
        T, V = tet_float(a)
        Lc = Lengths(T)
        bad = [e for e in EDGES if not edge_star_exists(T, e, Lc)[0]]
        cls = all_clusters(V, K)
        cand = tiled = 0
        hits = []
        for cl in cls:
            ok, B = bad_edges_ok(cl, V, bad)
            if not ok:
                continue
            cand += 1
            g = search(B)
            if g:
                bad_pts, n = coverage(B, g, cl)
                hits.append((len(cl), len(B.faces), bad_pts, n))
        print(f"#{idx}: clusters {len(cls)}, admissible (bad edges on Q-edges) {cand}, Poincaré hits {hits}", flush=True)
