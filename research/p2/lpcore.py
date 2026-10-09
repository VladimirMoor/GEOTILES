"""Ядро рёберного LP: какие комбинации обязаны иметь положительный вес в любом допустимом решении."""
import sys
import numpy as np
from scipy.optimize import linprog
from fractions import Fraction as Fr
from tetra import EDGES, PRINTED_A, Tet, angle_of
from obstruct import combos

idx = int(sys.argv[1])
a = PRINTED_A[idx - 1]
T = Tet(a)
al = [angle_of(a, *e) for e in EDGES]
L = np.array([float(x) for x in T.edge_lengths()])
cs = [c for c in combos(al, {Fr(1), Fr(2)}) if any(c)]
tgt = [sum(k * x for k, x in zip(c, al)) for c in cs]
A = np.array([[c[i] for c in cs] for i in range(6)], dtype=float)
names = [f"{e[0]}{e[1]}" for e in EDGES]
print("angles", dict(zip(names, map(str, al))), "lengths", dict(zip(names, np.round(L, 4))))
def show(j):
    return ("π:" if tgt[j] == 1 else "2π:") + "+".join(f"{k}·{names[i]}" if k > 1 else names[i] for i, k in enumerate(cs[j]) if k)
must = []
for j in range(len(cs)):
    c = np.zeros(len(cs)); c[j] = 1
    r = linprog(c, A_eq=A, b_eq=L, bounds=[(0, None)] * len(cs), method="highs")
    if r.status == 0 and r.fun > 1e-9:
        must.append((show(j), round(r.fun, 5)))
r = linprog(np.zeros(len(cs)), A_eq=A, b_eq=L, bounds=[(0, None)] * len(cs), method="highs")
print("feasible:", r.status == 0)
print("combos forced in every solution (min weight):", must)
print("an example solution:", [(show(j), round(x, 4)) for j, x in enumerate(r.x) if x > 1e-9])
