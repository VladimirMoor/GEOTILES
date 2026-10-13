"""Слои: замощения тетраэдром плоского слоя, периодические по двум сдвигам.

Если кластер M из k копий (склеенных грань к грани) и решётка сдвигов L = ⟨t1, t2⟩ дают M + L, заполняющее
слой между двумя параллельными плоскостями, то слои, уложенные друг на друга с любым сдвигом, замощают
пространство (как правило, не «лицом к лицу»). Сдвиги-кандидаты: t, для которых свободная грань F2 + t
совпадает со свободной гранью F1 и внешние нормали противоположны (M и M + t склеены по грани).
Критерий слоя: k·vol(T) = толщина(M вдоль нормали t1×t2) · |t1×t2| и нет перекрытий M с соседними M + a t1 + b t2.
Слоистые группы с поворотами/скольжениями покрываются за счёт мотива: подрешётка сдвигов с мотивом ≤ K копий.
Запуск: python slabs.py K [номера тетраэдров]  → data/slabs_K.json
"""
import itertools
import json
import pathlib
import sys

import numpy as np

from clusters import canon_key, congruences, face_ids, overlap, place_on_face, tet_float
from tetra import EXCLUDED_BY_LP, PRINTED_A

HERE = pathlib.Path(__file__).resolve().parent


def free_faces(cl):
    out = []
    for ti, Tt in enumerate(cl):
        for kk in range(4):
            fpts = np.array([Tt[i] for i in range(4) if i != kk])
            shared = any(all(any(np.linalg.norm(p - q) < 1e-6 for q in Tu) for p in fpts)
                         for tj, Tu in enumerate(cl) if tj != ti)
            if not shared:
                n = np.cross(fpts[1] - fpts[0], fpts[2] - fpts[0])
                n /= np.linalg.norm(n)
                if n @ (fpts[0] - Tt[kk]) < 0:
                    n = -n
                out.append((fpts, Tt[kk], n))
    return out


def all_clusters(V, K):
    start = [V.copy()]
    seen = {canon_key(start)}
    level, every = [start], [start]
    for k in range(2, K + 1):
        nxt = []
        for cl in level:
            for fpts, other, _ in free_faces(cl):
                for ks in range(4):
                    src = face_ids(ks)
                    for perm in congruences(V, src, fpts):
                        dst = [fpts[perm[x]] for x in range(3)]
                        newopp = place_on_face(None, dst, other, [V[s] for s in src], V[ks])
                        newT = np.zeros((4, 3))
                        for x, s in enumerate(src):
                            newT[s] = dst[x]
                        newT[ks] = newopp
                        if any(overlap(newT, Tu) for Tu in cl):
                            continue
                        ncl = cl + [newT]
                        key = canon_key(ncl)
                        if key in seen:
                            continue
                        seen.add(key)
                        nxt.append(ncl)
        level = nxt
        every += nxt
    return every


def translations(cl):
    F = free_faces(cl)
    ts = []
    for (A, _, nA), (B, _, nB) in itertools.permutations(F, 2):
        if nA @ nB > -1 + 1e-7:
            continue
        for perm in itertools.permutations(range(3)):
            t = A[0] - B[perm[0]]
            if all(np.linalg.norm(A[i] - B[perm[i]] - t) < 1e-6 for i in range(3)):
                if not any(np.linalg.norm(t - s) < 1e-6 for s in ts):
                    ts.append(t)
    return ts


def slab_test(cl, t1, t2, vol):
    c = np.cross(t1, t2)
    A = np.linalg.norm(c)
    if A < 1e-7:
        return False
    n = c / A
    P = np.vstack(cl)
    H = (P @ n).max() - (P @ n).min()
    if abs(len(cl) * vol - H * A) > 1e-6 * len(cl) * vol:
        return False
    for a, b in itertools.product(range(-2, 3), repeat=2):
        if (a, b) == (0, 0):
            continue
        s = a * t1 + b * t2
        if any(overlap(X + s, Y) for X in cl for Y in cl):
            return False
    return True


def slabs_for(V, K):
    vol = abs(np.linalg.det(np.array([V[1] - V[0], V[2] - V[0], V[3] - V[0]]))) / 6
    found = []
    for cl in all_clusters(V, K):
        ts = translations(cl)
        for t1, t2 in itertools.combinations(ts, 2):
            if slab_test(cl, t1, t2, vol):
                found.append({"k": len(cl), "t1": t1.tolist(), "t2": t2.tolist()})
                break
    return found


if __name__ == "__main__":
    K = int(sys.argv[1])
    which = [int(x) for x in sys.argv[2:]] or list(range(1, len(PRINTED_A) + 1))
    out = {}
    for i in which:
        a = PRINTED_A[i - 1]
        _, V = tet_float(a)
        f = slabs_for(V, K)
        out[i] = f
        print(f"#{i:2d} {'excl' if a in EXCLUDED_BY_LP else 'A   '} slabs (k≤{K}): {len(f)}"
              + (f"  smallest motif k={min(x['k'] for x in f)}" if f else ""), flush=True)
    (HERE / "data" / f"slabs_K{K}{'_' + '_'.join(map(str, which)) if len(which) < 10 else ''}.json").write_text(json.dumps(out, indent=1))
