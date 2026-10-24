"""Picture of one apico-basal neighbour exchange: L1 voxels near a T1 quartet, in a local tangent frame, at the
apical band (t < 0.2) and at the basal band (t > 0.8), coloured by cell.
Usage: python scutoid_figure.py file.h5 a b out.png [c d]  (with c, d: label only the four cells of the exchange)"""
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from scutoid import fields  # noqa: E402

path, a, b, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
quad = [a, b] + [int(x) for x in sys.argv[5:7]]
vox = np.array((0.235, 0.150, 0.150))
L, t, isL1, inner, d_out, _, _ = fields(path, tuple(vox))
P = np.argwhere((L == a) | (L == b)) * vox
c0 = P.mean(0)
Z = np.argwhere(isL1[L]) * vox
near = np.linalg.norm(Z - c0, axis=1) < 9.0
Z = Z[near]
idx = np.round(Z / vox).astype(int)
lab = L[idx[:, 0], idx[:, 1], idx[:, 2]]
tt = t[idx[:, 0], idx[:, 1], idx[:, 2]]
g = np.array(np.gradient(d_out))[:, idx[:, 0], idx[:, 1], idx[:, 2]].T / vox
n = g.mean(0); n /= np.linalg.norm(n)
e1 = np.cross(n, [1, 0, 0]); e1 /= np.linalg.norm(e1); e2 = np.cross(n, e1)
X = np.c_[(Z - c0) @ e1, (Z - c0) @ e2]
ids = np.unique(lab)
rng = np.random.default_rng(3)
col = {i: rng.random(3) * 0.7 + 0.2 for i in ids}
col = {i: np.array([0.82, 0.82, 0.82]) if len(quad) == 4 else c for i, c in col.items()}
for i, c in zip(quad, ([0.85, 0.15, 0.15], [0.15, 0.35, 0.9], [0.95, 0.7, 0.1], [0.2, 0.65, 0.3])):
    col[i] = np.array(c)
fig, ax = plt.subplots(1, 2, figsize=(12, 6))
for k, (sel, name) in enumerate(((tt < 0.2, "outer (apical) surface"), (tt > 0.8, "inner (basal) surface"))):
    ax[k].scatter(X[sel, 0], X[sel, 1], c=[col[i] for i in lab[sel]], s=3, marker="s")
    for i in ids:
        m = sel & (lab == i)
        if m.sum() > 50 and (len(quad) < 4 or i in quad):
            ax[k].text(*X[m].mean(0), "ABCD"[quad.index(i)] if len(quad) == 4 else str(i), fontsize=14 if len(quad) == 4 else 8,
                       ha="center", va="center", weight="bold")
    ax[k].set_title(name)
    ax[k].axis("off")
    ax[k].set_aspect("equal")
plt.tight_layout()
plt.savefig(out, dpi=80)
