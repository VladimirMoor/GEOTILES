"""Robustness of the iron result to reconstruction quality: on a central block, compare the symmetric share at
fixed F between grains with high and low mean LabDCT completeness, and between large and small grains.
If reconstruction noise were what destroys symmetry, reliable grains would be much more symmetric.
Usage: python fe_check.py"""
import h5py
import numpy as np

import topology
from topology import classify, face_graphs_from_voxels

f = h5py.File(topology.HERE / "data" / "fe" / "T0_AGG100_crop_gm3d.h5", "r")
Z = slice(1000, 1600)
L = f["LabDCT/Data/GrainId"][Z]
L[f["LabDCT/Data/Mask"][Z] == 0] = 0
Cm = f["LabDCT/Data/Completeness"][Z]
M = int(L.max()) + 1
cnt = np.bincount(L.ravel(), minlength=M)
comp = np.bincount(L.ravel(), weights=Cm.ravel(), minlength=M) / np.maximum(cnt, 1)
vg = face_graphs_from_voxels(L, {0})
ids = list(vg)
cl = classify(vg)
ok = [(i, c) for i, c in zip(ids, cl) if c]
cv = np.array([comp[i] for i, _ in ok])
vv = np.array([cnt[i] for i, _ in ok])
print(f"resolved {len(ok)}; mean completeness: median {np.median(cv):.3f}, quartiles {np.percentile(cv, [25, 75]).round(3)}")
for name, sel in (("completeness top half", cv >= np.median(cv)), ("completeness bottom half", cv < np.median(cv)),
                  ("volume top half", vv >= np.median(vv)), ("volume bottom half", vv < np.median(vv))):
    cs = [c for (_, c), s in zip(ok, sel) if s]
    row = {F: (round(float(np.mean([c["sym"] > 1 for c in cs if c["F"] == F])), 2), sum(c["F"] == F for c in cs)) for F in range(10, 16)}
    print(f"{name:26s} mean F {np.mean([c['F'] for c in cs]):.2f}  sym (share, n) by F {row}", flush=True)
