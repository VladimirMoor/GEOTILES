"""Комбинаторное предсказание свидетелей через графы зон и сравнение с прямым поиском тождеств.

Зона направления d: грани P, содержащие рёбра ∥ d, идут цепочкой (для ребра треугольника — от его копии
«сверху» к копии «снизу», для отрезка — замкнутым поясом). Каждой грани соответствует «видимая прямая»:
x для параллелограмма d ⊕ x и след плоскости треугольника для его собственной грани. Ребро зоны — дуга
(прямая грани до → прямая грани после). Граф формул: вершины — прямые, дуги — различные пары.
По лемме о независимости тождества между углами — ровно циркуляции в этих графах. Тождество фильтра —
неотрицательная целая циркуляция со значением Σ(π − φ) ∈ {π, 2π}; она раскладывается на ориентированные
циклы, каждый со значением ≥ π. Значит, ребро покрыто ⇔ оно лежит на ориентированном цикле значения ≤ 2π.
"""
import itertools
import json
import pathlib
import sys

import mpmath as mp
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
import converse as cv  # noqa: E402


def zone_structure(kind, p):
    parts = cv.member(kind, p)
    vk, fk, faces = cv.labels_and_faces(parts)
    keys, th = cv.angles(kind, p, faces)
    names = kind.split("+")
    ang = {e: float(t) for e, t in zip(keys, th)}

    def varying(labels):
        return tuple(i for i in range(len(names)) if len({lab[i] for lab in labels}) > 1)

    def face_line(F, zone_summand, zone_pair):
        vs = varying(F)
        if vs == (zone_summand,):
            return ("own", zone_summand)
        j = next(i for i in vs if i != zone_summand)
        vals = sorted({lab[j] for lab in F})
        return ("edge", j, tuple(vals[:2]) if names[j] == "tri" else (0, 1))

    # направления зон: (слагаемое i, неупорядоченная пара вершин ребра)
    zones = {}
    for e in keys:
        a, b = sorted(e)
        i = varying([a, b])[0]
        pair = tuple(sorted((a[i], b[i])))
        zones.setdefault((i, pair), []).append(e)
    covered, total = set(), 0
    edge_face = {}
    for fi, F in enumerate(faces):
        for k in range(len(F)):
            edge_face.setdefault(frozenset((F[k], F[(k + 1) % len(F)])), []).append(fi)
    for (i, pair), edges in zones.items():
        # цепочка граней: соседство через рёбра зоны
        adj = {}
        for e in edges:
            f1, f2 = edge_face[e]
            adj.setdefault(f1, []).append((f2, e))
            adj.setdefault(f2, []).append((f1, e))
        # начало: грань-треугольник (конец цепочки) или любая грань для замкнутого пояса
        start = next((f for f in adj if len(adj[f]) == 1), next(iter(adj)))
        order, used, cur = [start], set(), start
        seq_edges = []
        while True:
            nxt = [(g, e) for g, e in adj[cur] if e not in used]
            if not nxt:
                break
            g, e = nxt[0]
            used.add(e)
            seq_edges.append(e)
            order.append(g)
            cur = g
            if g == start:
                break
        lines = [face_line(faces[f], i, pair) for f in order]
        # дуги формул: (прямая_до, прямая_после) → угол φ = π − θ
        arcs = {}
        for k, e in enumerate(seq_edges):
            arc = (lines[k], lines[k + 1])
            arcs.setdefault(arc, []).append(e)
        # ориентированные простые циклы в графе формул
        nodes = sorted({x for arc in arcs for x in arc}, key=str)
        out = {}
        for (u, v) in arcs:
            out.setdefault(u, []).append(v)
        phi = {arc: np.pi - ang[es[0]] for arc, es in arcs.items()}
        good_arcs = set()

        def dfs(start_node, node, path, visited):
            for w in out.get(node, []):
                arc = (node, w)
                if w == start_node:
                    cyc = path + [arc]
                    val = sum(np.pi - phi[a] for a in cyc)
                    if val <= 2 * np.pi + 1e-9:
                        good_arcs.update(cyc)
                elif w not in visited and len(path) < 8:
                    dfs(start_node, w, path + [arc], visited | {w})

        for s in nodes:
            dfs(s, s, [], {s})
        for arc, es in arcs.items():
            total += 1
            if arc in good_arcs:
                covered.add((i, pair, arc))
    return total, total - len(covered)


def compare(kind, members, seed=5):
    rng = np.random.default_rng(seed)
    rows = []
    while len(rows) < members:
        nparams = 4 * (len(kind.split("+")) - 1)
        p0 = [mp.mpf(float(x)) for x in rng.normal(size=nparams)]
        try:
            total, predicted = zone_structure(kind, p0)
        except Exception as ex:  # noqa: BLE001
            continue
        rows.append((total, predicted))
    return rows


if __name__ == "__main__":
    for kind in (sys.argv[1:] or ["tri+tri", "tri+seg+seg", "seg+seg+seg+seg", "tri+tri+seg", "tri+seg+seg+seg",
                                  "seg+seg+seg+seg+seg", "tri+tri+tri"]):
        rows = compare(kind, 15)
        print(kind, "arc classes / predicted witnesses:", rows[:8], flush=True)
