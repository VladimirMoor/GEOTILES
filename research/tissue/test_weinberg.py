"""Проверка кода Вайнберга на многогранниках с известным порядком симметрии (вместе с отражениями).
Графы — двойственные триангуляции (вершины = грани многогранника)."""
import itertools
import random

import networkx as nx

from topology import weinberg


def tetrakis():  # двойственная к усечённому октаэдру (ячейка Кельвина, 14 граней)
    g = nx.Graph()
    sq = [(i, s) for i in range(3) for s in (1, -1)]
    hx = list(itertools.product((1, -1), repeat=3))
    for (i, s) in sq:
        for h in hx:
            if h[i] == s:
                g.add_edge(("q", i, s), ("h", h))
    for a, b in itertools.combinations(hx, 2):
        if sum(x != y for x, y in zip(a, b)) == 1:
            g.add_edge(("h", a), ("h", b))
    return g


CASES = {
    "tetrahedron (24)": (nx.complete_graph(4), 24),
    "cube (48)": (nx.octahedral_graph(), 48),
    "dodecahedron (120)": (nx.icosahedral_graph(), 120),
    "triangular prism (12)": (nx.complete_graph(5).copy(), 12),
    "truncated octahedron (48)": (tetrakis(), 48),
}
g = CASES["triangular prism (12)"][0]
g.remove_edge(3, 4)  # треугольная бипирамида = двойственная к треугольной призме
ok = True
for name, (g, sym) in CASES.items():
    code, mult, p = weinberg(g)
    # переименование вершин не меняет код
    perm = list(g.nodes()); random.shuffle(perm)
    code2, _, _ = weinberg(nx.relabel_nodes(g, dict(zip(g.nodes(), perm))))
    good = mult == sym and code == code2
    ok &= good
    print(f"{name:28s} sym {mult:3d} p-vector {p}  {'OK' if good else 'FAIL'}")
# разные многогранники с одинаковым числом граней получают разные коды
a = weinberg(nx.octahedral_graph())[0]
pent = nx.Graph([(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)] + [(5, i) for i in range(5)] + [(6, i) for i in range(5)])  # пентагональная призма? нет: двойственная — пентагональная бипирамида
print("distinct codes for 6-face vs 7-face:", a != weinberg(pent)[0])
print("ALL OK" if ok else "SOMETHING FAILED")
