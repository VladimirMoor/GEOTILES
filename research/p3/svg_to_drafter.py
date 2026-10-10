import re, sys, collections
import numpy as np
from scipy.spatial import cKDTree
s = open("basic6.svg").read()
polys = [(re.search(r'class="([^"]+)"', p).group(1),
          np.array([[float(v) for v in q.split(',')] for q in re.search(r'points="([^"]+)"', p).group(1).split()]))
         for p in re.findall(r'<polygon([^>]*)>', s)]
cells = [(int(v[0]), int(v[1]), np.array(list(map(float, v[2:]))).reshape(3, 2)) for v in (l.split() for l in open('/tmp/drafter_cells.txt'))]
C = np.array([t.mean(axis=0) for _, _, t in cells]); keys = [(x, y) for x, y, _ in cells]; tree = cKDTree(C)

def drafters_of(group):
    out = []
    for p in group:
        if len(p) == 3:
            out.append(p)
        else:
            c = p.mean(axis=0)
            for i in range(6):
                tri = [c, p[i], p[(i + 1) % 6]]; g = sum(tri) / 3
                for j in range(3):
                    v = tri[j]
                    for w in (tri[(j + 1) % 3], tri[(j + 2) % 3]):
                        out.append(np.array([g, v, (v + w) / 2]))
    return out

def align(dr):
    D = np.array([d.mean(axis=0) for d in dr]); Z = (D - D.mean(axis=0)) @ np.array([1, 1j])
    best = None
    for refl in (False, True):
        z0 = np.conj(Z) if refl else Z
        for r in range(12):
            for sc in np.linspace(0.2270, 0.2310, 17):
                w = np.exp(1j * r * np.pi / 6) * sc * z0; P = np.stack([w.real, w.imag], 1)
                Q = P - P.mean(axis=0)
                for _ in range(10):
                    d, idx = tree.query(Q); Q = Q + (C[idx] - Q).mean(axis=0)
                d, idx = tree.query(Q)
                if best is None or np.median(d) < best[1]:
                    best = (d.max(), np.median(d), len(set(idx)), idx.copy())
    return best

# копии плитки: группируем многоугольники одного класса в связные кластеры по 75
by_class = collections.defaultdict(list)
for c, p in polys:
    by_class[c].append(p)
done = 0
for c in sorted(by_class, key=lambda k: len(by_class[k]))[:30]:
    ps = by_class[c]
    if len(ps) % 75: continue
    # разбиение на копии: кластеризация центров (односвязная по расстоянию 40)
    cen = np.array([p.mean(axis=0) for p in ps])
    lab = -np.ones(len(ps), int); k = 0
    for i in range(len(ps)):
        if lab[i] >= 0: continue
        stack = [i]; lab[i] = k
        while stack:
            a = stack.pop()
            for b in np.where((np.linalg.norm(cen - cen[a], axis=1) < 40) & (lab < 0))[0]:
                lab[b] = k; stack.append(b)
        k += 1
    for t in range(k):
        grp = [ps[i] for i in range(len(ps)) if lab[i] == t]
        if len(grp) != 75: continue
        dr = drafters_of(grp)
        dmax, dmed, uniq, idx = align(dr)
        print(f"class {c} copy {t}: drafters {len(dr)}, max {dmax:.3f}, median {dmed:.3f}, unique {uniq}", flush=True)
        if dmax < 0.4 and uniq == len(dr):
            open("basic6_drafter.txt", "w").write("D? " + " ".join(f"{keys[i][0]} {keys[i][1]}" for i in idx) + "\n")
            print("written from", c); sys.exit(0)
        done += 1
        if done > 12: sys.exit(1)
