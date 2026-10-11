"""Полная топология клеток (по Лазару и др., PRL 2012): канонический код Вайнберга, порядок симметрии, p-вектор.

Клетка-многогранник задаётся своим «графом граней»: вершины — соседние клетки (грани), рёбра — пары граней,
сходящихся по общему ребру (тройной линии). Для простого многогранника это планарная триангуляция
(двойственный граф), и классы изоморфизма многогранников и их двойственных триангуляций совпадают.

Источники графа граней:
- воксельная сегментация: тройки меток в кубиках 2×2×2 (тройные линии); ребро (j, k) у клетки i есть,
  если тройка {i, j, k} встретилась хотя бы в MIN_TRIPLE кубиках;
- разбиение Вороного: треугольники триангуляции Делоне, содержащие точку i.
Граф принимается, только если это 3-связная планарная триангуляция (E = 3V − 6); иначе клетка «не разрешена».

Код Вайнберга считается на триангуляции; минимум по всем стартам (ориентированное ребро, ориентация),
число стартов с минимальным кодом = порядок группы автоморфизмов (вместе с отражениями).
"""
import json
import pathlib
import sys
from collections import Counter

import h5py
import networkx as nx
import numpy as np
from scipy.spatial import Delaunay

HERE = pathlib.Path(__file__).resolve().parent
MIN_TRIPLE = 2
QUADS = True
PRUNE = True


def triples_and_contacts(L, bg_ign):
    """Счётчики троек меток в кубиках 2×2×2 (ключ a·2³² + b·2¹⁶ + c, a < b < c)."""
    tri = Counter()
    nz = L.shape[0]
    for z0 in range(0, nz - 1, 12):
        z1 = min(z0 + 13, nz)
        B = L[z0:z1]
        V = np.stack([B[dz:B.shape[0] - 1 + dz, dy:B.shape[1] - 1 + dy, dx:B.shape[2] - 1 + dx].ravel()
                      for dz in (0, 1) for dy in (0, 1) for dx in (0, 1)], axis=1)
        V.sort(axis=1)
        d = V[:, 1:] != V[:, :-1]
        cnt = 1 + d.sum(axis=1)
        sel = cnt == 3
        if sel.any():
            S = V[sel]; D = d[sel]
            other = S[:, 1:][D].reshape(-1, 2)
            t = np.concatenate([S[:, :1], other], axis=1).astype(np.int64)
            keys = (t[:, 0] << 32) | (t[:, 1] << 16) | t[:, 2]
            u, c = np.unique(keys, return_counts=True)
            for k, n in zip(u.tolist(), c.tolist()):
                tri[k] += n
        if QUADS:
            sel = cnt == 4  # четверные точки (вершины): все 4 тройки, сходящиеся в вершине
            if sel.any():
                S = V[sel]; D = d[sel]
                q = np.concatenate([S[:, :1], S[:, 1:][D].reshape(-1, 3)], axis=1).astype(np.int64)
                for a, b, c in ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)):
                    keys = (q[:, a] << 32) | (q[:, b] << 16) | q[:, c]
                    u, cc = np.unique(keys, return_counts=True)
                    for k, n in zip(u.tolist(), cc.tolist()):
                        tri[k] += n
        del V
    return tri


def prune(g):
    """Узкие «щепки»: грань со степенью 2 стягиваем в ребро между её соседями, степень < 2 — удаляем."""
    g = g.copy()
    changed = True
    while changed and g.number_of_nodes() > 4:
        changed = False
        for v, d in list(g.degree()):
            if d < 3:
                nb = list(g.neighbors(v))
                g.remove_node(v)
                if d == 2:
                    g.add_edge(*nb)
                changed = True
                break
    return g


def face_graphs_from_voxels(L, bg_ign):
    tri = triples_and_contacts(L, bg_ign)
    G = {}
    border = set(bg_ign)
    for k, n in tri.items():
        a, b, c = k >> 32, (k >> 16) & 0xFFFF, k & 0xFFFF
        trip = (a, b, c)
        if any(x in bg_ign for x in trip):
            border.update(trip)
            continue
        if n < MIN_TRIPLE:
            continue
        for i, j, m in ((a, b, c), (b, a, c), (c, a, b)):
            G.setdefault(i, nx.Graph()).add_edge(j, m, w=n)
    # клетки на краю объёма
    for ax in range(3):
        for idx in (0, L.shape[ax] - 1):
            border.update(np.unique(np.take(L, idx, axis=ax)).tolist())
    return {i: (prune(g) if PRUNE else g) for i, g in G.items() if i not in border}


def face_graphs_from_delaunay(P, inner):
    tri = Delaunay(P)
    G = {i: nx.Graph() for i in inner}
    for s in tri.simplices:
        for a in range(4):
            i = s[a]
            if i in G:
                o = [s[b] for b in range(4) if b != a]
                G[i].add_edges_from(((o[0], o[1]), (o[0], o[2]), (o[1], o[2])))
    return G


def weinberg(g):
    """(канонический код, порядок симметрии, p-вектор) для 3-связной планарной триангуляции, иначе None."""
    n, m = g.number_of_nodes(), g.number_of_edges()
    if n < 4 or m != 3 * n - 6 or min(d for _, d in g.degree()) < 3:
        return None
    ok, emb = nx.check_planarity(g)
    if not ok or nx.node_connectivity(g) < 3:
        return None
    rot = {v: list(emb.neighbors_cw_order(v)) for v in g}
    pos = {v: {u: k for k, u in enumerate(rot[v])} for v in g}
    deg = dict(g.degree())
    dmin = min(deg.values())
    best, mult = None, 0
    for u in g:
        if deg[u] != dmin:
            continue
        for v in rot[u]:
            for sgn in (1, -1):
                lab, code, used = {u: 1}, [1], set()
                a, b = u, v
                while True:
                    used.add((a, b))
                    if b not in lab:
                        lab[b] = len(lab) + 1
                        code.append(lab[b])
                        nxt = rot[b][(pos[b][a] + sgn) % deg[b]]
                    else:
                        code.append(lab[b])
                        if (b, a) not in used:
                            nxt = a
                        else:
                            k = pos[b][a]
                            nxt = None
                            for s in range(1, deg[b] + 1):
                                c = rot[b][(k + sgn * s) % deg[b]]
                                if (b, c) not in used:
                                    nxt = c
                                    break
                            if nxt is None:
                                break
                    a, b = b, nxt
                    if best is not None and tuple(code) > best[:len(code)]:
                        break
                t = tuple(code)
                if best is None or t < best:
                    best, mult = t, 1
                elif t == best:
                    mult += 1
    pv = Counter(deg.values())
    pvec = tuple(pv.get(k, 0) for k in range(3, max(pv) + 1))
    return best, mult, pvec


def repair(g):
    """Лишние рёбра (почти вырожденные вершины, где сходятся 4 клетки): убираем самые слабые, пока E > 3V − 6."""
    g = g.copy()
    while g.number_of_edges() > 3 * g.number_of_nodes() - 6:
        e = min(g.edges(data="w"), key=lambda t: t[2] or 0)
        g.remove_edge(e[0], e[1])
    return prune(g)


def classify(graphs, fix=False):
    out = []
    for i, g in graphs.items():
        r = weinberg(g)
        if r is None and fix:
            r = weinberg(repair(g))
        out.append(None if r is None else {"F": g.number_of_nodes(), "code": ",".join(map(str, r[0])), "sym": r[1], "p": list(r[2])})
    return out


if __name__ == "__main__":
    path = sys.argv[1]
    vox = np.array([float(v) for v in sys.argv[2:5]])
    h = h5py.File(path, "r")
    L = h["label"][:]
    bg_ign = {0}
    if "label_with_ignore" in h:
        Iw = h["label_with_ignore"]
        for z0 in range(0, L.shape[0], 16):
            I = Iw[z0:z0 + 16]
            L[z0:z0 + 16][I == -1] = 65000
        bg_ign.add(65000)
    vg = face_graphs_from_voxels(L, bg_ign)
    cells = classify(vg)
    cells_fix = classify(vg, fix=True)
    nfaces = [g.number_of_nodes() for g in vg.values()]
    # Вороной центров тех же (внутренних) клеток
    M = int(L.max()) + 1
    cnt = np.zeros(M); s = np.zeros((M, 3))
    yy, xx = np.indices(L.shape[1:])
    for z in range(L.shape[0]):
        f = L[z].ravel()
        c = np.bincount(f, minlength=M)
        cnt += c; s[:, 0] += c * z
        s[:, 1] += np.bincount(f, weights=yy.ravel(), minlength=M)
        s[:, 2] += np.bincount(f, weights=xx.ravel(), minlength=M)
    del L
    labels = [l for l in range(1, M) if cnt[l] > 0 and l not in bg_ign]
    idx = {l: k for k, l in enumerate(labels)}
    P = np.array([s[l] / cnt[l] for l in labels]) * vox
    vor = classify(face_graphs_from_delaunay(P, {idx[l] for l in vg}))
    res = {"file": pathlib.Path(path).name, "cells": cells, "cells_repaired": cells_fix, "faces_all": nfaces, "voronoi": vor}
    out = HERE / "data" / "topo" / (pathlib.Path(path).stem + ".json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(res))
    ok = [c for c in cells if c]
    print(res["file"], "inner", len(cells), "resolved", len(ok), "repaired", sum(1 for c in cells_fix if c), "voronoi resolved", sum(1 for c in vor if c), "sym>1:", sum(1 for c in ok if c["sym"] > 1), "/", len(ok), " voronoi sym>1:", sum(1 for c in vor if c and c["sym"] > 1))
