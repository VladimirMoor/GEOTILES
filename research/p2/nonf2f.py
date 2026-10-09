"""Какие из 40 тетраэдров вообще допускают не «лицом к лицу» разбиения (Chentouf–Sun, Prop. 2.3).
π-комбинация: Σ k_e α_e = π, k ≥ 0. Не-f2f исключено, если π-комбинаций нет, или она одна и
затрагивает ровно пару противоположных рёбер."""
from fractions import Fraction as Fr
from tetra import EDGES, PRINTED_A, EXCLUDED_BY_LP, angle_of
from obstruct import combos

OPP = [{(1, 2), (3, 4)}, {(1, 3), (2, 4)}, {(1, 4), (2, 3)}]


def pi_combos(a):
    return [c for c in combos(list(a), {Fr(1)}) if any(c)]


def status(a):
    pcs = pi_combos(a)
    if not pcs:
        return "f2f only (no π-combination)", pcs
    if len(pcs) == 1:
        supp = {EDGES[i] for i, k in enumerate(pcs[0]) if k}
        if supp in OPP:
            return "f2f only (single opposite-pair π-combination)", pcs
    return "non-f2f possible", pcs


if __name__ == "__main__":
    for idx, a in enumerate(PRINTED_A, 1):
        if a in EXCLUDED_BY_LP:
            continue
        s, pcs = status(a)
        show = [" + ".join(f"{k}·α{EDGES[i][0]}{EDGES[i][1]}" for i, k in enumerate(c) if k) for c in pcs]
        print(f"#{idx:2d} {tuple(map(str, a))}: {s}; π-combos: {show}")
