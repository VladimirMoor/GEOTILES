"""Figure for the site: a slice of the real foam with segmented bubbles, and the share of symmetric cells by
number of faces for real foam, Voronoi of the foam's bubble centres, coarsening foam (Potts), tissue, random cells.
Inputs: downloaded slices (foam_fetch.py), data/compare_systems.json. Output: docs/img/foam-real.png
"""
import json
import pathlib

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from skimage.segmentation import find_boundaries  # noqa: E402

from foam_seg import load, segment  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE.parents[1] / "docs" / "img" / "foam-real.png"

V, origin = load(260, 340, (700, 1300, 700, 1300), 1)
L, _ = segment(V, origin, 1, 940.0, 12, 1.0)
z = L.shape[0] // 2
raw = V[z]
lo, hi = np.percentile(raw, [1, 99.5])

C = json.load(open(HERE / "data" / "compare_systems.json"))
series = [("foam_real", "real foam (bubbles)", "#d62728", "-"), ("foam_real_PV", "Voronoi of bubble centres", "#d62728", "--"),
          ("potts_small", "coarsening foam (Potts)", "#ff7f0e", "-"), ("tissue", "plant tissue (ovules)", "#2ca02c", "-"),
          ("exact_PV", "random cells (Poisson–Voronoi)", "#7f7f7f", "-")]

fig, ax = plt.subplots(1, 2, figsize=(13, 6), gridspec_kw={"width_ratios": [1, 1.15]})
ax[0].imshow(raw, cmap="gray", vmin=lo, vmax=hi)
ax[0].imshow(np.ma.masked_where(~find_boundaries(L[z], mode="inner"), np.ones_like(raw)), cmap="autumn", alpha=0.7)
ax[0].set_title("Real foam, X-ray tomography (1.65 × 1.65 mm)")
ax[0].axis("off")
for key, label, col, ls in series:
    byF = C[key]["byF"]
    F = sorted(int(f) for f in byF)
    ax[1].plot(F, [byF[str(f)]["sym"] for f in F], ls, color=col, marker="o", label=label)
ax[1].set_xlabel("number of faces F")
ax[1].set_ylabel("share of cells with a symmetry")
ax[1].set_xlim(7.5, 20.5)
ax[1].set_ylim(0, 1)
ax[1].grid(alpha=0.3)
ax[1].legend(frameon=False)
plt.tight_layout()
plt.savefig(OUT, dpi=90)
print(OUT)
