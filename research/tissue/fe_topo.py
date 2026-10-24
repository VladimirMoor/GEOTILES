"""Complete topology of real metal grains: pure iron (alpha-Fe), laboratory diffraction contrast tomography.

Data: Ånes et al., Zenodo 20323210 (CC BY 4.0), state T0 of a grain-growth series, GrainMapper3D HDF5.
GrainId is a space-filling label volume (3292 x 371 x 370 voxels, 3 µm), with 0 outside the sample mask, so the
tissue pipeline (topology.py) applies directly. A grain touching the mask or the ends of the volume is not inner.
Also: grain volumes and centroids, and the exact Voronoi tessellation of all grain centroids for grains at least
two median diameters (32 voxels) inside the sample.

Output: data/topo_fe/fe_T0.json (same keys as data/topo_foam/foam_step102.json). Usage: python fe_topo.py
"""
import json
import pathlib

import h5py
import numpy as np
from scipy import ndimage as ndi

from topology import classify, face_graphs_from_delaunay, face_graphs_from_voxels

HERE = pathlib.Path(__file__).resolve().parent
SRC = HERE / "data" / "fe" / "T0_AGG100_crop_gm3d.h5"
OUT = HERE / "data" / "topo_fe"
DEEP = 32

if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    f = h5py.File(SRC, "r")
    L = f["LabDCT/Data/GrainId"][:]
    L[f["LabDCT/Data/Mask"][:] == 0] = 0
    assert L.max() < 65000
    M = int(L.max()) + 1
    cnt = np.bincount(L.ravel(), minlength=M).astype(float)
    sums = np.zeros((M, 3))
    for z in range(L.shape[0]):  # centroids slice by slice to keep memory low
        s = L[z]
        sums[:, 0] += np.bincount(s.ravel(), minlength=M) * z
        yy, xx = np.indices(s.shape)
        sums[:, 1] += np.bincount(s.ravel(), weights=yy.ravel(), minlength=M)
        sums[:, 2] += np.bincount(s.ravel(), weights=xx.ravel(), minlength=M)
    with np.errstate(invalid="ignore", divide="ignore"):
        C = sums / cnt[:, None]
    vg = face_graphs_from_voxels(L, {0})
    print("inner grains", len(vg), flush=True)
    E = ndi.distance_transform_edt(L > 0).astype(np.float32)  # distance to the outside of the sample
    del L
    labels = np.flatnonzero(cnt[1:] > 0) + 1
    inner = sorted(vg)
    ci = np.round(C[inner]).astype(int)
    depth = E[ci[:, 0], ci[:, 1], ci[:, 2]]
    deep = (depth > DEEP).tolist()
    del E
    cells, fixed = classify({i: vg[i] for i in inner}), classify({i: vg[i] for i in inner}, fix=True)
    idx = {int(l): k for k, l in enumerate(labels)}
    P = C[labels]
    vor = classify(face_graphs_from_delaunay(P, {idx[i] for i, d in zip(inner, deep) if d}))
    vor_by = dict(zip(sorted(idx[i] for i, d in zip(inner, deep) if d), vor))
    res = {"file": SRC.name, "cells": cells, "cells_repaired": fixed, "faces_all": [vg[i].number_of_nodes() for i in inner],
           "voronoi": [vor_by.get(idx[i]) if d else None for i, d in zip(inner, deep)], "deep": deep,
           "volumes_inner": cnt[inner].tolist(), "centroids_inner": C[inner].tolist(), "depth_inner": depth.tolist(),
           "volumes_all": cnt[labels].tolist()}
    (OUT / "fe_T0.json").write_text(json.dumps(res))
    ok = [c for c in cells if c]
    print("inner", len(cells), "resolved", len(ok), "mean F", np.mean([c["F"] for c in ok]),
          "sym>1", np.mean([c["sym"] > 1 for c in ok]), "deep", sum(deep), flush=True)
