"""Метрические звёзды частичного угла: для каждого ребра e — множество сумм углов θ ≤ 2π, при которых
существует цепочка клиньев (соседние плитки делят одинаковую грань), начинающаяся с самой T на ребре e.
2π с замыканием — полная звезда; π — полузвезда (обе крайние грани лежат в одной плоскости)."""
from fractions import Fraction as Fr
from tetra import EDGES, PRINTED_A, EXCLUDED_BY_LP, Tet
from stars import Lengths, wedges, edge_star_exists

KNOWN = {1, 3}


def chain_angles(T, e, Lc):
    a, b = e
    W = wedges(T, Lc, Lc(a, b))
    starts = [w for w in W if {w[3][0], w[3][1]} == {a, b}]
    reach = set()
    for w0 in starts:
        # цепочка может продолжаться в обе стороны; сумма углов всех клиньев
        # вперёд от выхода w0 и назад от входа w0
        fwd = {(w0[1], Fr(0))}
        stack = list(fwd)
        while stack:
            f, acc = stack.pop()
            for w in W:
                if w[0] == f and acc + w[2] + w0[2] <= 2:
                    st = (w[1], acc + w[2])
                    if st not in fwd:
                        fwd.add(st); stack.append(st)
        bwd = {(w0[0], Fr(0))}
        stack = list(bwd)
        while stack:
            f, acc = stack.pop()
            for w in W:
                if w[1] == f and acc + w[2] + w0[2] <= 2:
                    st = (w[0], acc + w[2])
                    if st not in bwd:
                        bwd.add(st); stack.append(st)
        for _, x in fwd:
            for _, y in bwd:
                if x + y + w0[2] <= 2:
                    reach.add(x + y + w0[2])
    return reach


if __name__ == "__main__":
    for idx, a in enumerate(PRINTED_A, 1):
        if a in EXCLUDED_BY_LP:
            continue
        T = Tet(a)
        Lc = Lengths(T)
        row = []
        for e in EDGES:
            full, _ = edge_star_exists(T, e, Lc)
            R = chain_angles(T, e, Lc)
            half = Fr(1) in R
            row.append(f"{e[0]}{e[1]}:" + ("2π" if full else ("π" if half else "—")))
        print(f"#{idx:2d}{' (known tiler)' if idx in KNOWN else ''}: " + "  ".join(row))
