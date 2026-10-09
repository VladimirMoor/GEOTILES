import numpy as np
from tetra import PRINTED_A, EDGES, Tet
from stars import Lengths, edge_star_exists
from clusters import tet_float
from nonconvex import all_clusters, bad_edges_ok, search, coverage
a = PRINTED_A[41]; T, V = tet_float(a); Lc = Lengths(T)
bad = [e for e in EDGES if not edge_star_exists(T, e, Lc)[0]]
for cl in all_clusters(V, 5):
    ok, B = bad_edges_ok(cl, V, bad)
    if not ok: continue
    g = search(B)
    if not g: continue
    print("cluster size", len(cl), "faces", [len(f) for f in B.faces])
    print("dihedral angles of Q:", sorted(set(str(x) for x in B.dihedral.values())))
    # точная проверка: сколько раз покрыта точка внутри самого Q
    cen = np.vstack(cl).mean(axis=0)
    elems = [np.eye(4)] + [M for (M, _, _) in g.values()]
    def inside(x, Tt):
        A = np.column_stack([Tt[1]-Tt[0], Tt[2]-Tt[0], Tt[3]-Tt[0]]); l = np.linalg.solve(A, x-Tt[0])
        return l.min() > 1e-9 and l.sum() < 1-1e-9
    rng = np.random.default_rng(1)
    # точки внутри Q: центроиды тетраэдров кластера
    for Tt in cl:
        x = Tt.mean(axis=0)
        hits = [k for k, h in enumerate(elems) for U in cl if inside(x, np.array([(h @ np.append(p,1))[:3] for p in U]))]
        print("point in Q covered by images:", hits)
    print("generators: det of linear parts", [round(float(np.linalg.det(M[:3,:3])),3) for (M,_,_) in g.values()])
    print("coverage:", coverage(B, g, cl, samples=100))
    break
