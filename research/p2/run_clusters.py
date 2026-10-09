"""Перебор выпуклых кластеров и проверка Пуанкаре."""
import sys
import numpy as np
from tetra import PRINTED_A, EXCLUDED_BY_LP
from clusters import tet_float, enumerate_clusters
from poincare_q import Poly, search

K = int(sys.argv[1])
which = [int(x) for x in sys.argv[2:]] or None
for idx, a in enumerate(PRINTED_A, 1):
    if a in EXCLUDED_BY_LP or (which and idx not in which):
        continue
    T, V = tet_float(a)
    cv = enumerate_clusters(V, K)
    hits = []
    for cl in cv:
        Q = Poly(np.vstack(cl))
        f, n = search(Q)
        if f:
            hits.append((len(cl), len(Q.faces)))
    print(f"#{idx}: convex clusters {len(cv)}, tiling clusters {hits}", flush=True)
