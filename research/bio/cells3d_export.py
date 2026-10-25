"""3D scenes of protein-crystal Voronoi cells for the site (docs/data/cells/*.json, docs/js/cells3d.js).

- protein_chains.json: a P2₁2₁2₁ crystal with two chains in the asymmetric unit; the Voronoi cells of the eight
  chain centres of one unit cell, coloured by chain (A blue, B orange), shrunk a little so the cells separate.
- protein_parity.json: a single-chain P2₁ crystal whose molecule cell has 14 faces; faces are coloured in pairs:
  the neighbour reached by operation g and the one reached by g⁻¹ share a colour, which is why the count is even.
Entries are chosen from data/chains_*.json (cached downloads from pdb_chains.py).

Usage: python cells3d_export.py
"""
import json
import pathlib

import numpy as np
from scipy.spatial import ConvexHull, Voronoi

from pdb_chains import Group, analyse, ca_by_chain, cell_info, frac_matrix

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE.parents[1] / "docs" / "data" / "cells"
PAL = (np.array([[0.89, 0.35, 0.33], [0.31, 0.56, 0.85], [0.95, 0.70, 0.20], [0.31, 0.70, 0.54],
                 [0.55, 0.47, 0.85], [0.25, 0.70, 0.77], [0.85, 0.51, 0.69], [0.60, 0.67, 0.27]]) ** 2.2)


def setup(pid):
    info = cell_info(pid)
    cl = info["cell"]
    M = frac_matrix(cl["length_a"], cl["length_b"], cl["length_c"], cl["angle_alpha"], cl["angle_beta"], cl["angle_gamma"])
    G = Group(info["symmetry"]["space_group_name_H_M"])
    ch = {k: v for k, v in ca_by_chain(pid).items() if len(v) >= 20}
    return M, np.linalg.inv(M), G, [(M @ v.mean(0)) % 1.0 for v in ch.values()]


def images(points, G, reach=2):
    """All images: (frac position, chain index, operation index, lattice shift)."""
    out = []
    for k, x in enumerate(points):
        for o, (R, t) in enumerate(zip(G.R, G.t)):
            y = R @ x + t
            base = np.floor(y)
            for s in np.array([(i, j, l) for i in range(-reach, reach + 1) for j in range(-reach, reach + 1) for l in range(-reach, reach + 1)]):
                out.append((y - base + s, k, o, tuple((s - base).astype(int))))
    return out


def cell_faces(vor, i):
    """Facets of the Voronoi cell of point i: list of (neighbour index, vertex array)."""
    faces = []
    for (p, q), rv in zip(vor.ridge_points, vor.ridge_vertices):
        if i in (p, q) and -1 not in rv:
            faces.append((q if p == i else p, vor.vertices[rv]))
    return faces


def facet_mesh(faces, centre, shrink, colours):
    V, F, C = [], [], []
    for (nb, pts), col in zip(faces, colours):
        pts = centre + (pts - centre) * shrink
        c = pts.mean(0)
        nrm = np.cross(pts[1] - pts[0], pts[2] - pts[0])
        u = pts[0] - c
        w = np.cross(nrm, u)
        order = np.argsort(np.arctan2((pts - c) @ w, (pts - c) @ u))
        pts = pts[order]
        base = len(V)
        V += pts.tolist()
        C += [col] * len(pts)
        F += [[base, base + j, base + j + 1] for j in range(1, len(pts) - 1)]
    return np.array(V), np.array(F), np.array(C)


def pack(V, F, C, **kw):
    return {"v": np.round(V, 4).ravel().tolist(), "f": F.ravel().tolist(), "c": np.round(C, 3).ravel().tolist(), "flat": True, **kw}


def chains_scene(pid):
    M, O, G, pts = setup(pid)
    imgs = images(pts, G)
    P = np.array([y for y, *_ in imgs]) @ O.T
    vor = Voronoi(P)
    inside = [i for i, (y, *_) in enumerate(imgs) if (y >= 0).all() and (y < 1).all()]
    centre = np.array([0.5, 0.5, 0.5]) @ O.T
    meshes, allV = [], []
    for i in inside:
        k = imgs[i][1]
        col = PAL[1] if k == 0 else PAL[2]
        faces = cell_faces(vor, i)
        V, F, C = facet_mesh(faces, P[i], 0.9, [col] * len(faces))
        meshes.append((V, F, C))
        allV.append(V)
    allV = np.vstack(allV)
    r = np.linalg.norm(allV - centre, axis=1).max()
    out = [pack((V - centre) / r, F, C) for V, F, C in meshes]
    return out, len(inside)


def parity_scene(pid):
    M, O, G, pts = setup(pid)
    assert len(pts) == 1
    imgs = images(pts, G)
    P = np.array([y for y, *_ in imgs]) @ O.T
    vor = Voronoi(P)
    i0 = next(i for i, (y, k, o, s) in enumerate(imgs) if o == 0 and s == (0, 0, 0))
    faces = cell_faces(vor, i0)
    # neighbour j is the image of the molecule under the affine map A_j = (R_o, t_o + shift); its partner is A_j⁻¹
    def affine(j):
        _, _, o, s = imgs[j]
        return G.R[o], G.t[o] + np.array(s)
    keys = []
    for j, _ in faces:
        R, t = affine(j)
        keys.append((np.round(R, 6), np.round(t, 6)))
    colours, pair_of, used = [None] * len(faces), {}, 0
    for a in range(len(faces)):
        if colours[a] is not None:
            continue
        Ra, ta = keys[a]
        Ri = np.round(np.linalg.inv(Ra), 6)
        ti = np.round(-np.linalg.inv(Ra) @ ta, 6)
        b = next((b for b in range(len(faces)) if b != a and colours[b] is None
                  and np.allclose(keys[b][0], Ri) and np.allclose(keys[b][1], ti)), None)
        colours[a] = PAL[used % len(PAL)]
        if b is not None:
            colours[b] = colours[a]
        used += 1
    V, F, C = facet_mesh(faces, P[i0], 1.0, colours)
    c = V.mean(0)
    r = np.linalg.norm(V - c, axis=1).max()
    return [pack((V - c) / r, F, C)], len(faces), used


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    R2 = [r for r in json.load(open(HERE / "data" / "chains_P212121.json")) if r.get("n_chains") == 2 and r["chain_faces"] == [15, 15]]
    for r in R2:
        try:
            meshes, n = chains_scene(r["id"])
            break
        except Exception as e:  # noqa: BLE001
            print("skip", r["id"], e)
    (OUT / "protein_chains.json").write_text(json.dumps({"id": r["id"], "meshes": meshes}))
    print("protein_chains.json", r["id"], n, "cells")
    R1 = [r for r in json.load(open(HERE / "data" / "chains_P1211.json")) if r.get("n_chains") == 1 and r["asu_faces"] == 14]
    for r in R1:
        try:
            meshes, nf, pairs = parity_scene(r["id"])
            break
        except Exception as e:  # noqa: BLE001
            print("skip", r["id"], e)
    (OUT / "protein_parity.json").write_text(json.dumps({"id": r["id"], "meshes": meshes}))
    print("protein_parity.json", r["id"], nf, "faces in", pairs, "colour groups")
