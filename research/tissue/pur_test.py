"""Choice of the marker depth h for the PUR foam: segment a central block for several h, report cell counts,
volume spread, topology statistics, and draw label boundaries over a slice (data/pur/test_overlay.png)."""
import pathlib
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from skimage.segmentation import find_boundaries  # noqa: E402

from pur_seg import load, segment  # noqa: E402
from topology import classify, face_graphs_from_voxels  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
Z, CROP = (150, 450), (175, 575, 175, 575)
HS = [float(x) for x in sys.argv[1:]] or [1.0, 2.0, 4.0]

S = load(*Z, CROP)
fig, ax = plt.subplots(1, len(HS), figsize=(6 * len(HS), 6))
for k, h in enumerate(HS):
    L, info = segment(S, (CROP[0], CROP[2]), 360.0, h)
    border = set()
    for a in range(3):
        for i in (0, L.shape[a] - 1):
            border |= set(np.unique(np.take(L, i, axis=a)).tolist())
    v = np.bincount(L.ravel())
    inner = [l for l in range(1, len(v)) if v[l] and l not in border]
    w = v[inner]
    ok = [c for c in classify(face_graphs_from_voxels(L, {0})) if c]
    F = np.array([c["F"] for c in ok])
    print(f"h={h}: {info}, inner {len(inner)}, vol cv {w.std() / w.mean():.2f}, D_eq {(6 * np.median(w) / np.pi) ** (1 / 3):.1f} vox, "
          f"resolved {len(ok)}, F {F.mean():.2f}±{F.std():.2f}, sym {np.mean([c['sym'] > 1 for c in ok]):.3f}", flush=True)
    z = L.shape[0] // 2
    ax[k].imshow(S[z], cmap="gray_r")
    ax[k].imshow(np.ma.masked_where(~find_boundaries(L[z]), np.ones(L[z].shape)), cmap="autumn", alpha=0.9)
    ax[k].set_title(f"h = {h}: F = {F.mean():.2f} ± {F.std():.2f}")
plt.tight_layout()
plt.savefig(HERE / "data" / "pur" / "test_overlay.png", dpi=70)
