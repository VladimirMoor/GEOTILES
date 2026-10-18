"""Complete topology of the cells of an open-cell polyurethane foam (DaRUS doi:10.18419/darus-3010, t0), block by block.

Same scheme as foam_topo.py: overlapping blocks, segmentation by pur_seg.segment (h = 1, chosen with pur_check.py),
face graphs from topology.py, each cell counted once from the block whose core holds its centroid, and the exact
Voronoi tessellation of all cell centroids for cells at least two diameters (110 voxels) from the data boundary.

Output: data/topo_pur/pur_t0.json (same keys as data/topo_foam/foam_step102.json). Usage: python pur_topo.py [h]
"""
import json
import pathlib
import sys
import time

import numpy as np

from foam_topo import centroids
from pur_seg import FILES, N, load, segment
from topology import classify, face_graphs_from_delaunay, face_graphs_from_voxels

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "data" / "topo_pur"
R_MASK, MARGIN, DEEP = 360.0, 80, 110
LO, HI = int(N / 2 - R_MASK), int(N / 2 + R_MASK) + 1


def blocks(nz):
    for a in range(0, nz, 200):
        for y in range(LO, HI, 250):
            for x in range(LO, HI, 250):
                core = (a, min(a + 200, nz), y, min(y + 250, HI), x, min(x + 250, HI))
                box = (max(core[0] - MARGIN, 0), min(core[1] + MARGIN, nz), max(core[2] - MARGIN, 0),
                       min(core[3] + MARGIN, N), max(core[4] - MARGIN, 0), min(core[5] + MARGIN, N))
                yield core, box


if __name__ == "__main__":
    h = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
    nz = len(FILES)
    OUT.mkdir(parents=True, exist_ok=True)
    cells, fixed, faces, vols, P, inner_idx, infos = [], [], [], [], [], [], []
    for core, box in blocks(nz):
        t = time.time()
        S = load(box[0], box[1], (box[2], box[3], box[4], box[5]))
        L, info = segment(S, (box[2], box[4]), R_MASK, h)
        del S
        if info["cells"] == 0:
            continue
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
    Q = P[inner_idx]
    rad = np.hypot(Q[:, 1] - (N - 1) / 2, Q[:, 2] - (N - 1) / 2)
    deep = ((Q[:, 0] > DEEP) & (nz - Q[:, 0] > DEEP) & (rad < R_MASK - DEEP)).tolist()
    res = {"file": f"darus-3010 t0_segmented, h={h}", "cells": cells, "cells_repaired": fixed, "faces_all": faces,
           "voronoi": [vor_by[i] if d else None for i, d in zip(inner_idx, deep)], "deep": deep,
           "volumes_inner": [vols[i] for i in inner_idx], "centroids_inner": Q.tolist(), "blocks": infos}
    (OUT / f"pur_t0_h{h:g}.json").write_text(json.dumps(res))
    ok = [c for c in cells if c]
    print("inner", len(cells), "resolved", len(ok), "mean F", np.mean([c["F"] for c in ok]),
          "sym>1", np.mean([c["sym"] > 1 for c in ok]), "deep", sum(deep), flush=True)
