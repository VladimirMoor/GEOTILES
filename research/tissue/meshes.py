"""3D surfaces of real cells for the site (docs/data/cells/*.json, shown by docs/js/cells3d.js).

- scutoid.json: the verified apico-basal exchange in ovule N_405_A (cells A = 33, B = 343, C = 31, D = 156),
  turned so that the outer surface of the layer is on top; cells are shrunk slightly so the walls show as gaps.
- gallery.json: one typical cell of each system: a liquid-foam bubble of Matzke's 1-10-2 type, a pentagonal-
  dodecahedron bubble, an open-cell polyurethane cell, an iron grain, an ovule tissue cell, and an exact
  Poisson–Voronoi cell, each with its number of faces and symmetry order.
Surfaces: marching cubes on a Gaussian-smoothed cell mask (voxel size respected), centred, scaled to radius 1.

Usage: python meshes.py [scutoid|gallery]   (needs the ovule file, the foam slices, the PUR stack, the iron map)
"""
import json
import pathlib
import sys
import urllib.request

import numpy as np
from scipy import ndimage as ndi
from scipy.spatial import ConvexHull, Voronoi
from skimage.measure import marching_cubes

from topology import classify, face_graphs_from_voxels

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE.parents[1] / "docs" / "data" / "cells"
OVULE = HERE / "data" / "raw" / "N_405_A_ds2x.h5"
OV_VOX = (0.235, 0.150, 0.150)


def surface(mask, vox, sigma=1.0, step=1):
    """Triangle mesh (vertices in µm, faces) of a binary mask."""
    m = np.pad(mask.astype(np.float32), 3)
    m = ndi.gaussian_filter(m, sigma)
    v, f, _, _ = marching_cubes(m, 0.5, spacing=vox, step_size=step)
    return v - 3 * np.array(vox), f


def pack(v, f, color, **kw):
    if "c" in kw:
        kw["c"] = np.round(kw["c"], 3).ravel().tolist()
    return {"v": np.round(v, 3).ravel().tolist(), "f": f.astype(int).ravel().tolist(), "color": color, **kw}


FACE_COLOURS = np.array([[0.89, 0.35, 0.33], [0.31, 0.56, 0.85], [0.95, 0.70, 0.20], [0.31, 0.70, 0.54],
                         [0.55, 0.47, 0.85], [0.25, 0.70, 0.77], [0.85, 0.51, 0.69], [0.60, 0.67, 0.27]]) ** 2.2  # sRGB → linear (three.js vertex colours)


def face_colours(Lbox, lab, v_vox, bg=(0,)):
    """Per-vertex colours: each face of the cell (its contact with one neighbour) gets one colour; neighbouring
    faces get different colours (greedy colouring of the face adjacency). v_vox: vertex positions in voxels of Lbox."""
    _, idx = ndi.distance_transform_edt(Lbox == lab, return_indices=True)
    near = Lbox[tuple(idx)]  # label of the nearest voxel outside the cell
    q = np.clip(np.round(v_vox).astype(int), 0, np.array(Lbox.shape) - 1)
    nb = near[q[:, 0], q[:, 1], q[:, 2]]
    # face adjacency: two neighbours whose patches meet at a mesh vertex neighbourhood
    labs = [x for x in np.unique(nb) if x not in bg and x != lab]
    adj = {x: set() for x in labs}
    pad = np.pad(near, 1, mode="edge")
    for d in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
        a = nb
        b = pad[q[:, 0] + 1 + d[0], q[:, 1] + 1 + d[1], q[:, 2] + 1 + d[2]]
        for x, y in set(zip(a.tolist(), b.tolist())):
            if x != y and x in adj and y in adj:
                adj[x].add(y); adj[y].add(x)
    col = {}
    for x in sorted(labs, key=lambda z: -len(adj[z])):
        used = {col[y] for y in adj[x] if y in col}
        col[x] = next(k for k in range(len(FACE_COLOURS)) if k not in used) if len(used) < len(FACE_COLOURS) else 0
    return np.array([FACE_COLOURS[col[x]] if x in col else [0.5, 0.5, 0.5] for x in nb.tolist()])


def normalise(parts, rot=None):
    """Centre all meshes at their common centroid, rotate, scale to radius 1."""
    c = np.vstack([v for v, _ in parts]).mean(0)
    out = []
    for v, f in parts:
        w = v - c
        if rot is not None:
            w = w @ rot.T
        out.append((w, f))
    r = max(np.linalg.norm(w, axis=1).max() for w, _ in out)
    return [(w / r, f) for w, f in out]


def rotation_to_z(n):
    """Rotation matrix taking unit vector n to +z."""
    n = n / np.linalg.norm(n)
    z = np.array([0.0, 0.0, 1.0])
    v, c = np.cross(n, z), n @ z
    if np.linalg.norm(v) < 1e-9:
        return np.eye(3) if c > 0 else np.diag([1, -1, -1.0])
    K = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + K + K @ K * (1 / (1 + c))


def ensure_ovule():
    if not OVULE.exists():
        OVULE.parent.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve("https://osf.io/download/kjhcw/", OVULE)


def to_xyz(v):
    """Volume axes (z, y, x) → (x, y, z) for display."""
    return v[:, ::-1].copy()


def scutoid_scene():
    from scutoid import fields
    ensure_ovule()
    L, t, span, inner, d_out, _, _ = fields(OVULE, OV_VOX)
    cells = {"A": (33, "#d62728"), "B": (343, "#2b59d9"), "C": (31, "#f0b419"), "D": (156, "#2ca25f")}
    sel = np.isin(L, [c for c, _ in cells.values()])
    zz, yy, xx = np.nonzero(sel)
    lo = np.maximum([zz.min() - 4, yy.min() - 4, xx.min() - 4], 0)
    hi = np.array([zz.max() + 5, yy.max() + 5, xx.max() + 5])
    box = tuple(slice(a, b) for a, b in zip(lo, hi))
    Lb = L[box]
    # outward normal = −grad(d_out), averaged over the four cells
    g = np.array(np.gradient(d_out[box], *OV_VOX))
    m = np.isin(Lb, [c for c, _ in cells.values()])
    n = -np.array([gi[m].mean() for gi in g])
    parts, meta = [], []
    for name, (lab, col) in cells.items():
        v, f = surface(Lb == lab, OV_VOX, sigma=1.2, step=2)
        v = to_xyz(v)
        cen = v.mean(0)
        v = cen + (v - cen) * 0.94  # shrink so that the walls between cells show as gaps
        parts.append((v, f))
        meta.append((name, col))
    rot = rotation_to_z(n[::-1])
    parts = normalise(parts, rot)
    meshes = []
    for (v, f), (nm, col) in zip(parts, meta):
        top = v[v[:, 2] > np.percentile(v[:, 2], 80)].mean(0) if nm in "AB" else v[v[:, 2] < np.percentile(v[:, 2], 20)].mean(0)
        meshes.append(pack(v, f, col, label=nm, labelPos=np.round(top + [0, 0, 0.12 if nm in "AB" else -0.12], 3).tolist()))
    allv = np.vstack([v for v, _ in parts])
    (x0, y0), (x1, y1) = allv[:, :2].min(0) - 0.15, allv[:, :2].max(0) + 0.15
    for z, txt in ((allv[:, 2].max() + 0.05, {"en": "outside", "ru": "снаружи", "es": "exterior"}),
                   (allv[:, 2].min() - 0.05, {"en": "inside", "ru": "внутри", "es": "interior"})):
        sq = np.array([[x0, y0, z], [x1, y0, z], [x1, y1, z], [x0, y1, z]])
        meshes.append(pack(sq, np.array([[0, 1, 2], [0, 2, 3]]), "#9fb3c8", opacity=0.18, label=txt,
                           labelPos=[round(float(x1), 3), round(float(y1), 3), round(float(z), 3)]))
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "scutoid.json").write_text(json.dumps({"meshes": meshes}))
    print("scutoid.json", [len(m["f"]) // 3 for m in meshes], "triangles")


def pick(L, bg, want, vox=(1, 1, 1)):
    """Label of a resolved inner cell matching want(cell), the one nearest the block centre."""
    vg = face_graphs_from_voxels(L, bg)
    ids = list(vg)
    cl = classify(vg)
    cen = np.array(L.shape) / 2
    best = None
    for i, c in zip(ids, cl):
        if c and want(c):
            p = np.argwhere(L == i).mean(0)
            d = np.linalg.norm((p - cen) * vox)
            if best is None or d < best[0]:
                best = (d, i, c)
    return best[1], best[2]


def cell_mesh(L, lab, vox, sigma=1.0, step=1, bg=(0,)):
    """Mesh of one cell, normalised, with per-vertex face colours."""
    zz, yy, xx = np.nonzero(L == lab)
    box = (slice(max(zz.min() - 3, 0), zz.max() + 4), slice(max(yy.min() - 3, 0), yy.max() + 4), slice(max(xx.min() - 3, 0), xx.max() + 4))
    Lb = L[box]
    v, f = surface(Lb == lab, vox, sigma, step)
    c = face_colours(Lb, lab, v / np.array(vox), bg)
    (w, f), = normalise([(to_xyz(v), f)])
    return w, f, c


def note(c):
    return {"en": f"{c['F']} faces · symmetry {c['sym']}", "ru": f"{c['F']} граней · симметрия {c['sym']}",
            "es": f"{c['F']} caras · simetría {c['sym']}"}


def gallery():
    items = []
    # liquid foam: Matzke's 1-10-2 bubble and a pentagonal dodecahedron
    from foam_seg import load as fload, segment as fseg
    V, origin = fload(280, 380, (800, 1200, 800, 1200), 1)
    Lf, _ = fseg(V, origin, 1, 940.0, 12, 1.0)
    for want, nm in ((lambda c: c["p"][:4] == [0, 1, 10, 2] and sum(c["p"][4:]) == 0,
                      {"en": "Soap-foam bubble, Matzke's 1-10-2", "ru": "Пузырь пены, тип Мацке 1-10-2", "es": "Burbuja de espuma, tipo 1-10-2 de Matzke"}),
                     (lambda c: c["p"][:3] == [0, 0, 12] and sum(c["p"][3:]) == 0,
                      {"en": "Soap-foam bubble, pentagonal dodecahedron", "ru": "Пузырь пены, пятиугольный додекаэдр", "es": "Burbuja de espuma, dodecaedro pentagonal"})):
        lab, c = pick(Lf, {0}, want)
        v, f, col = cell_mesh(Lf, lab, (1, 1, 1), sigma=1.5, step=2)
        items.append(pack(v, f, "#d62728", name=nm, note=note(c), c=col))
    # open-cell polyurethane foam
    from pur_seg import load as pload, segment as pseg
    S = pload(220, 380, (230, 520, 230, 520))
    Lp, _ = pseg(S, (230, 230), 360.0, 1.0)
    lab, c = pick(Lp, {0}, lambda c: c["F"] == 14)
    v, f, col = cell_mesh(Lp, lab, (1, 1, 1), sigma=1.5, step=2)
    items.append(pack(v, f, "#8d78d8", c=col, name={"en": "Polyurethane-foam cell", "ru": "Ячейка полиуретановой пены", "es": "Celda de espuma de poliuretano"}, note=note(c)))
    # iron grain
    import h5py
    h = h5py.File(HERE / "data" / "fe" / "T0_AGG100_crop_gm3d.h5", "r")
    Li = h["LabDCT/Data/GrainId"][1500:1620, 110:260, 110:260].astype(np.int32)
    lab, c = pick(Li, {0}, lambda c: c["F"] == 14)
    v, f, col = cell_mesh(Li, lab, (1, 1, 1), sigma=1.0, step=1)
    items.append(pack(v, f, "#1f77b4", c=col, name={"en": "Iron grain", "ru": "Зерно железа", "es": "Grano de hierro"}, note=note(c)))
    # ovule tissue cell
    ensure_ovule()
    hh = h5py.File(OVULE, "r")
    Lt = hh["label"][:].astype(np.int32)
    Lt[hh["label_with_ignore"][:] == -1] = 65000
    lab, c = pick(Lt, {0, 65000}, lambda c: c["F"] == 14, OV_VOX)
    v, f, col = cell_mesh(Lt, lab, OV_VOX, sigma=1.5, step=2, bg=(0, 65000))
    items.append(pack(v, f, "#2ca25f", c=col, name={"en": "Plant cell (ovule)", "ru": "Клетка растения (семязачаток)", "es": "Célula vegetal (óvulo)"}, note=note(c)))
    # exact Poisson–Voronoi cell with 15 faces
    rng = np.random.default_rng(4)
    P = rng.random((400, 3))
    vor = Voronoi(P)
    for i in np.argsort(np.linalg.norm(P - 0.5, axis=1)):
        reg = vor.regions[vor.point_region[i]]
        if -1 in reg or not reg:
            continue
        hull = ConvexHull(vor.vertices[reg])
        nf = sum(1 for (a, b) in vor.ridge_points if i in (a, b))
        if nf == 15:
            break
    pts = (hull.points - hull.points.mean(0)) / np.linalg.norm(hull.points - hull.points.mean(0), axis=1).max()
    eq = np.round(hull.equations, 6)
    facet = {tuple(e): k for k, e in enumerate(np.unique(eq, axis=0))}
    fid = np.array([facet[tuple(e)] for e in eq])
    fv = {k: set(hull.simplices[fid == k].ravel()) for k in set(fid)}
    colour = {}
    for k in sorted(fv):
        used = {colour[j] for j in colour if len(fv[k] & fv[j]) >= 2}
        colour[k] = next(c for c in range(len(FACE_COLOURS)) if c not in used)
    v = pts[hull.simplices].reshape(-1, 3)
    f = np.arange(len(v)).reshape(-1, 3)
    col = np.repeat(FACE_COLOURS[[colour[k] for k in fid]], 3, axis=0)
    items.append(pack(v, f, "#7f7f7f", flat=True, c=col,
                      name={"en": "Random Voronoi cell", "ru": "Случайная ячейка Вороного", "es": "Celda de Voronói aleatoria"},
                      note={"en": "15 faces · no symmetry", "ru": "15 граней · без симметрии", "es": "15 caras · sin simetría"}))
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "gallery.json").write_text(json.dumps({"meshes": items}))
    print("gallery.json", [(m["name"]["en"], len(m["f"]) // 3) for m in items])


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    if what in ("scutoid", "all"):
        scutoid_scene()
    if what in ("gallery", "all"):
        gallery()
