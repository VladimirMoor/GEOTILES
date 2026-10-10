"""SVG плитки на сетке драфтеров по координатам heesch-sat (клетки + внешний контур)."""
import sys
from collections import Counter
import numpy as np

cells = {}
for line in open("data/drafter_cells.txt"):
    v = line.split()
    cells[(int(v[0]), int(v[1]))] = np.array(list(map(float, v[2:]))).reshape(3, 2)
tile_file, out = sys.argv[1], sys.argv[2]
c = list(map(int, open(tile_file).read().split()[1:]))
tris = [cells[(c[i], c[i + 1])] for i in range(0, len(c), 2)]
P = np.vstack(tris)
lo, hi = P.min(axis=0), P.max(axis=0)
pad = 2.0
W, H = hi - lo + 2 * pad
key = lambda p: (round(p[0], 3), round(p[1], 3))
edges = Counter()
for t in tris:
    for i in range(3):
        a, b = key(t[i]), key(t[(i + 1) % 3])
        edges[tuple(sorted((a, b)))] += 1
def xy(p):
    return f"{p[0] - lo[0] + pad:.3f},{hi[1] - p[1] + pad:.3f}"
svg = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W:.2f} {H:.2f}" width="{W*6:.0f}" height="{H*6:.0f}">']
for t in tris:
    svg.append(f'<polygon points="{" ".join(xy(p) for p in t)}" fill="#4f8fd9" fill-opacity="0.55" stroke="#4f8fd9" stroke-opacity="0.35" stroke-width="0.08"/>')
for (a, b), n in edges.items():
    if n == 1:
        svg.append(f'<line x1="{a[0]-lo[0]+pad:.3f}" y1="{hi[1]-a[1]+pad:.3f}" x2="{b[0]-lo[0]+pad:.3f}" y2="{hi[1]-b[1]+pad:.3f}" stroke="#1b2a3a" stroke-width="0.35" stroke-linecap="round"/>')
svg.append("</svg>")
open(out, "w").write("\n".join(svg))
print("cells", len(tris), "boundary edges", sum(1 for n in edges.values() if n == 1))
