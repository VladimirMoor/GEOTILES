"""Суммы метрических цепочек через плохие рёбра (кандидаты на двугранные углы кластера < π)."""
from fractions import Fraction as Fr
from tetra import EDGES, PRINTED_A, EXCLUDED_BY_LP, Tet
from stars import Lengths, edge_star_exists
from halfstars import chain_angles

for idx, a in enumerate(PRINTED_A, 1):
    if a in EXCLUDED_BY_LP or idx in (1, 3):
        continue
    T = Tet(a); Lc = Lengths(T)
    parts = []
    for e in EDGES:
        if edge_star_exists(T, e, Lc)[0]:
            continue
        R = sorted(x for x in chain_angles(T, e, Lc) if x < 1)
        parts.append(f"{e[0]}{e[1]}:{{{', '.join(map(str, R))}}}")
    print(f"#{idx}: " + "  ".join(parts))
