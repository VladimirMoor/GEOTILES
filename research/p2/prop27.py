"""Покрытие предложением 2.7 Chentouf–Sun: f2f исключено, если T параллелограммный или все длины
различны, и какой-то двугранный угол не делит 2π."""
from fractions import Fraction as Fr
from tetra import EDGES, PRINTED_A, EXCLUDED_BY_LP, Tet
from stars import Lengths

KNOWN = {1: "Sommerville No.3", 3: "Sommerville No.2"}


def pattern(Lc):
    c = [Lc.cls[e] for e in EDGES]  # порядок 12,34,13,24,14,23 — пары противоположных рёбер подряд
    return c


def parallelogram_like(c):
    pairs = [(c[0], c[1]), (c[2], c[3]), (c[4], c[5])]
    eq = [p for p in pairs if p[0] == p[1]]
    return len(eq) == 2 and eq[0][0] != eq[1][0] and len(set(c)) == 4


if __name__ == "__main__":
    new = []
    for idx, a in enumerate(PRINTED_A, 1):
        if a in EXCLUDED_BY_LP or idx in KNOWN:
            continue
        T = Tet(a)
        c = pattern(Lengths(T))
        distinct = len(set(c)) == 6
        par = parallelogram_like(c)
        not_div = any((Fr(2) / x).denominator != 1 for x in a)
        covered = (distinct or par) and not_div
        if not covered:
            new.append(idx)
        print(f"#{idx:2d} lengths pattern {c} distinct={distinct} parallelogram-like={par} "
              f"-> Prop2.7 {'covers' if covered else 'does NOT cover'}")
    print("f2f exclusion new (not covered by Prop 2.7):", new, len(new))
