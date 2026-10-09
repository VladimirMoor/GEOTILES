"""Граф разрезаний: какие тетраэдры списка складываются из копий других (или собственных) тетраэдров.

Если тетраэдр U склеен face-to-face из k копий T, то:
    U замощает ⟹ T замощает;   T не замощает ⟹ U не замощает.
Ищем кластеры из ≤ K копий T, выпуклая оболочка которых — тетраэдр (4 грани), и опознаём его
по упорядоченным двугранным углам с точностью до перестановки вершин (и подобия).
"""
import sys
from fractions import Fraction as Fr

import mpmath as mp
import numpy as np
from scipy.spatial import ConvexHull

from tetra import EDGES, PRINTED_A, EXCLUDED_BY_LP
from known import canon, sommerville
from clusters import tet_float, enumerate_clusters
from halves import dihedral

mp.mp.dps = 60


def tet_angles(pts):
    P = {i + 1: mp.matrix([mp.mpf(float(x)) for x in p]) for i, p in enumerate(pts)}
    out = []
    for e in EDGES:
        x = dihedral(P, *e)
        f = Fr(float(x)).limit_denominator(720)
        if abs(float(x) - float(f)) > 1e-9:
            return None
        out.append(f)
    return tuple(out)


def hill(x):
    return (Fr(1, 2), Fr(1, 2), 1 - 2 * x, Fr(1, 3), x, x)


if __name__ == "__main__":
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    names = {}
    for idx, a in enumerate(PRINTED_A, 1):
        if a not in EXCLUDED_BY_LP:
            names.setdefault(canon(a), f"#{idx}")
    for nm, a in sommerville().items():
        names.setdefault(canon(a), nm)
    for idx, a in enumerate(PRINTED_A, 1):
        if a in EXCLUDED_BY_LP:
            continue
        T, V = tet_float(a)
        found = {}
        for cl in enumerate_clusters(V, K):
            pts = np.vstack(cl)
            H = ConvexHull(pts)
            verts = pts[H.vertices]
            if len(verts) != 4:
                continue
            ang = tet_angles(verts)
            if ang is None:
                continue
            c = canon(ang)
            nm = names.get(c)
            if nm is None:
                # Hill family?
                for q in range(2, 61):
                    for p in range(1, q):
                        x = Fr(p, q)
                        if Fr(1, 6) < x < Fr(1, 2) and canon(hill(x)) == c:
                            nm = f"Hill(x={x})"
                nm = nm or "other " + str(tuple(map(str, ang)))
            found.setdefault(nm, len(cl))
        print(f"#{idx}: tetrahedra made of ≤{K} copies: {found}", flush=True)
