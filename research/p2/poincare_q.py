"""Теорема Пуанкаре для выпуклого многогранника Q (кластера из копий T).

Грани Q — многоугольники (копланарные треугольники выпуклой оболочки сливаются). Спаривание граней:
G → F изометрией γ, при этом γQ лежит по другую сторону F, а γ_{F'} = γ_F⁻¹. Самоспаривание —
инволюция: отражение в плоскости грани, полуоборот вокруг оси симметрии грани или центральная
симметрия грани, скомпонованная с отражением. Условие на рёбрах: обход замыкается при сумме
двугранных углов ровно 2π с тождественным итоговым элементом.
Если условие выполнено, образы Q разбивают пространство, а значит, T замощает пространство.
"""
import itertools
import sys
from fractions import Fraction as Fr

import numpy as np
from scipy.spatial import ConvexHull

TOL = 1e-6


class Poly:
    def __init__(self, pts):
        H = ConvexHull(pts)
        self.P = pts[H.vertices]
        idx = {int(v): i for i, v in enumerate(H.vertices)}
        c = self.P.mean(axis=0)
        planes = []
        for eq in H.equations:
            n, d = eq[:3], eq[3]
            if not any(np.allclose(n, m, atol=1e-7) and abs(d - e) < 1e-7 for m, e in planes):
                planes.append((n, d))
        self.faces = []
        self.normals = []
        for n, d in planes:
            vs = [i for i in range(len(self.P)) if abs(self.P[i] @ n + d) < 1e-7]
            fc = self.P[vs].mean(axis=0)
            u = self.P[vs[0]] - fc; u /= np.linalg.norm(u)
            w = np.cross(n, u)
            vs.sort(key=lambda i: np.arctan2((self.P[i] - fc) @ w, (self.P[i] - fc) @ u))
            self.faces.append(vs)  # против часовой стрелки, если смотреть снаружи
            self.normals.append(n)
        self.edges = {}
        for fi, f in enumerate(self.faces):
            for a, b in zip(f, f[1:] + f[:1]):
                self.edges.setdefault(frozenset((a, b)), []).append(fi)
        assert all(len(v) == 2 for v in self.edges.values())
        self.dihedral = {}
        for e, (f1, f2) in self.edges.items():
            cosv = -self.normals[f1] @ self.normals[f2]
            self.dihedral[e] = Fr(float(np.arccos(np.clip(cosv, -1, 1)) / np.pi)).limit_denominator(720)


def affine_from(src, dst):
    A = np.vstack([np.array(src).T, np.ones(4)])
    B = np.vstack([np.array(dst).T, np.ones(4)])
    return B @ np.linalg.inv(A)


def is_isometry(M):
    L = M[:3, :3]
    return np.allclose(L.T @ L, np.eye(3), atol=1e-7)


def pairings(Q, fG, fF):
    """Изометрии γ: грань fG → грань fF, переводящие Q на другую сторону fF. Возвращает (M, vmap)."""
    G, F = Q.faces[fG], Q.faces[fF]
    m = len(G)
    if len(F) != m:
        return []
    out = []
    nG, nF = Q.normals[fG], Q.normals[fF]
    for s in range(m):
        for d in (1, -1):
            vmap = {G[i]: F[(s + d * i) % m] for i in range(m)}
            src = [Q.P[G[0]], Q.P[G[1]], Q.P[G[2]], Q.P[G[0]] - nG]
            dst = [Q.P[vmap[G[0]]], Q.P[vmap[G[1]]], Q.P[vmap[G[2]]], Q.P[vmap[G[0]]] + nF]
            M = affine_from(src, dst)
            if not is_isometry(M):
                continue
            if all(np.allclose((M @ np.append(Q.P[g], 1))[:3], Q.P[vmap[g]], atol=1e-6) for g in G):
                if fG == fF:
                    if not np.allclose(M @ M, np.eye(4), atol=1e-6):
                        continue
                out.append((M, vmap))
    return out


def partitions(n):
    def rec(rest):
        if not rest:
            yield []
            return
        a = rest[0]
        for tail in rec(rest[1:]):
            yield [(a,)] + tail
        for b in rest[1:]:
            r2 = [x for x in rest[1:] if x != b]
            for tail in rec(r2):
                yield [(a, b)] + tail
    yield from rec(list(range(n)))


def edge_check(Q, gam):
    """gam[f] = (M, vmap, g): изометрия переводит грань g → f. Проверка циклов рёбер."""
    for e0, (fa, fb) in Q.edges.items():
        g = np.eye(4)
        edge, leave = e0, fa
        total = Fr(0)
        while True:
            total += Q.dihedral[edge]
            M, vmap, src = gam[leave]
            g = g @ M
            inv = {v: k for k, v in vmap.items()}
            edge2 = frozenset(inv[x] for x in edge)
            f1, f2 = Q.edges[edge2]
            nxt = f2 if f1 == src else f1
            if f1 == f2:
                return False
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


def search(Q, max_schemes=200000):
    nf = len(Q.faces)
    cache = {}
    found = []
    count = 0

    def opts(a, b):
        if (a, b) not in cache:
            cache[(a, b)] = pairings(Q, a, b)
        return cache[(a, b)]

    for part in partitions(nf):
        choices = []
        ok = True
        for blk in part:
            if len(blk) == 1:
                o = [(blk[0], blk[0], p) for p in opts(blk[0], blk[0])]
            else:
                o = [(blk[0], blk[1], p) for p in opts(blk[1], blk[0])]
            if not o:
                ok = False
                break
            choices.append(o)
        if not ok:
            continue
        for combo in itertools.product(*choices):
            count += 1
            if count > max_schemes:
                return found, count
            gam = {}
            for f, gsrc, (M, vmap) in combo:
                gam[f] = (M, vmap, gsrc)
                if gsrc != f:
                    Minv = np.linalg.inv(M)
                    gam[gsrc] = (Minv, {v: k for k, v in vmap.items()}, f)
            if edge_check(Q, gam):
                found.append(combo)
                return found, count
    return found, count
