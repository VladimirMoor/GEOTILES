"""Complete topology of real foam bubbles (Zenodo 15564270, step 102) at full resolution, block by block.

The z-range of downloaded slices is cut into overlapping blocks (core + margin). Each block is segmented with
foam_seg.segment (bin 1, bubble diameter ≈ 42 voxels, close to the tissue resolution of ≈ 39), face graphs come
from topology.py, and each bubble is kept only from the block whose core contains its centroid; a bubble that
touches the block boundary or the sample mask is not "inner". Also: the exact Voronoi tessellation of all
bubble centroids (as for the tissue), to separate packing order from surface-tension order.

Output: data/topo_foam/foam_step102.json (same keys as data/topo/*.json). Usage: python foam_topo.py [z0 z1]
"""
import json
import pathlib
import sys
import time

import numpy as np

from foam_seg import N, load, segment
from topology import classify, face_graphs_from_delaunay, face_graphs_from_voxels

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "data" / "topo_foam"
CORE_Z, CORE_XY, MARGIN = 137, 627, 64
R_MASK, LO = 940.0, 68  # sample disk radius and its bounding box start (N/2 - R_MASK)


def blocks(z0, z1):
    zc = list(range(z0, z1, CORE_Z))
    xy = list(range(LO, N - LO, CORE_XY))
    for a in zc:
        for y in xy:
            for x in xy:
                core = (a, min(a + CORE_Z, z1), y, min(y + CORE_XY, N - LO), x, min(x + CORE_XY, N - LO))
                box = (max(core[0] - MARGIN, z0), min(core[1] + MARGIN, z1), max(core[2] - MARGIN, 0),
                       min(core[3] + MARGIN, N), max(core[4] - MARGIN, 0), min(core[5] + MARGIN, N))
                yield core, box


def centroids(L):
    M = int(L.max()) + 1
    cnt = np.bincount(L.ravel(), minlength=M).astype(float)
    zz = np.arange(L.shape[0])[:, None, None]
    yy = np.arange(L.shape[1])[None, :, None]
    xx = np.arange(L.shape[2])[None, None, :]
    c = np.stack([np.bincount(L.ravel(), weights=np.broadcast_to(g, L.shape).ravel(), minlength=M) for g in (zz, yy, xx)], 1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return c / cnt[:, None], cnt


if __name__ == "__main__":
    z0, z1 = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (100, 510)
    OUT.mkdir(parents=True, exist_ok=True)
    cells, fixed, faces, vols, P, inner_idx, infos = [], [], [], [], [], [], []
    for core, box in blocks(z0, z1):
        t = time.time()
        V, origin = load(box[0], box[1], (box[2], box[3], box[4], box[5]), 1)
        if V.size == 0:
            continue
        L, info = segment(V, origin, 1, R_MASK, 12, 1.0)
        del V
        C, cnt = centroids(L)
        G = C + np.array([box[0], box[2], box[4]])
        in_core = ((G[:, 0] >= core[0]) & (G[:, 0] < core[1]) & (G[:, 1] >= core[2]) & (G[:, 1] < core[3])
                   & (G[:, 2] >= core[4]) & (G[:, 2] < core[5]) & (cnt > 0))
        in_core[0] = False
        vg = face_graphs_from_voxels(L, {0})
        del L
        own = {i: g for i, g in vg.items() if in_core[i]}
        cl, cf = classify(own), classify(own, fix=True)
        base = len(P)
        lab = np.flatnonzero(in_core)
        pos = {int(l): base + k for k, l in enumerate(lab)}
        P += G[lab].tolist()
        vols += cnt[lab].tolist()
        for (i, g), a, b in zip(own.items(), cl, cf):
            cells.append(a); fixed.append(b); faces.append(g.number_of_nodes()); inner_idx.append(pos[i])
        info.update(core=core, box=box, inner=len(own), resolved=sum(1 for c in cl if c), sec=round(time.time() - t))
        infos.append(info)
        print(info, flush=True)
    P = np.array(P)
    vor = classify(face_graphs_from_delaunay(P, set(inner_idx)))
    vor_by = dict(zip(sorted(set(inner_idx)), vor))
    # Voronoi cells of centroids are reliable only away from the ends of the downloaded range and the sample edge
    Q = P[inner_idx]
    rad = np.hypot(Q[:, 1] - N / 2, Q[:, 2] - N / 2)
    deep = ((Q[:, 0] - z0 > 2 * 42) & (z1 - Q[:, 0] > 2 * 42) & (rad < R_MASK - 2 * 42)).tolist()
    res = {"file": f"tomo_rheoscopy_102 z{z0}-{z1}", "cells": cells, "cells_repaired": fixed, "faces_all": faces,
           "voronoi": [vor_by[i] if d else None for i, d in zip(inner_idx, deep)], "deep": deep,
           "volumes_inner": [vols[i] for i in inner_idx], "centroids_inner": Q.tolist(), "blocks": infos}
    (OUT / "foam_step102.json").write_text(json.dumps(res))
    ok = [c for c in cells if c]
    print("inner", len(cells), "resolved", len(ok), "mean F", np.mean([c["F"] for c in ok]),
          "sym>1", np.mean([c["sym"] > 1 for c in ok]), "voronoi resolved", sum(1 for c in vor if c), flush=True)
