"""Resolution check for the foam segmentation on a central block: full resolution vs 2x binning.

Segments z 100..228, y,x 600..1400 at bin 1 and bin 2, runs the topology pipeline on both, prints
face-count and symmetry statistics, and draws label boundaries over a slice (data/foam/test_overlay.png).
"""
import pathlib

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from skimage.segmentation import find_boundaries  # noqa: E402

from foam_seg import load, segment  # noqa: E402
from topology import classify, face_graphs_from_voxels  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "data" / "foam"
Z, CROP = (100, 228), (600, 1400, 600, 1400)

res = {}
fig, ax = plt.subplots(1, 2, figsize=(14, 7))
for k, b in enumerate((1, 2)):
    V, origin = load(*Z, CROP, b)
    L, info = segment(V, origin, b, 940.0, max(2, round(12 / b)), 1.0)
    cells = classify(face_graphs_from_voxels(L, {0}))
    ok = [c for c in cells if c]
    F = np.array([c["F"] for c in ok])
    res[b] = dict(info, inner=len(cells), resolved=len(ok), meanF=F.mean(), sdF=F.std(),
                  sym=np.mean([c["sym"] > 1 for c in ok]))
    print(b, res[b], flush=True)
    zs = L.shape[0] // 2
    raw = V[zs]
    ax[k].imshow(raw, cmap="gray")
    ax[k].imshow(np.ma.masked_where(~find_boundaries(L[zs]), np.ones_like(raw)), cmap="autumn", alpha=0.8)
    ax[k].set_title(f"bin {b}: {info['bubbles']} markers, F = {F.mean():.2f} ± {F.std():.2f}")
    np.savez_compressed(OUT / f"test_bin{b}.npz", L=L)
plt.tight_layout()
plt.savefig(OUT / "test_overlay.png", dpi=80)
