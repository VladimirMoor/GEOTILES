"""Dirichlet–Voronoi cells of protein chains in PDB crystals, for any space group and any cell.

pdb_dv.py takes one point per entry (the centroid of the asymmetric unit) and needs a lattice with one shape
parameter. Here every protein chain (label_asym_id with at least MIN_CA Cα atoms) gets its own centroid; all
centroids are expanded by the space-group operators and by lattice translations −2..2, and the Voronoi
tessellation of these points in Cartesian coordinates (real a, b, c, α, β, γ) gives the cell of every chain.
Faces = neighbouring chains. "asu" mode: one point per entry, as in pdb_dv.py, for comparison.

Usage: python pdb_chains.py "P 21 21 21" [sample size] [seed]  → data/chains_<group>.json
"""
import json
import pathlib
import random
import sys
import urllib.request
from collections import Counter

import numpy as np
from scipy.spatial import Voronoi

from pdb_dv import UA, cell_info, fetch, frac_matrix

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "p4"))
from dvfast import Group  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
MIN_CA = 20


def all_entries(sg):
    ids, start = [], 0
    while True:
        q = {"query": {"type": "group", "logical_operator": "and", "nodes": [
            {"type": "terminal", "service": "text", "parameters": {"attribute": "symmetry.space_group_name_H_M", "operator": "exact_match", "value": sg}},
            {"type": "terminal", "service": "text", "parameters": {"attribute": "exptl.method", "operator": "exact_match", "value": "X-RAY DIFFRACTION"}}]},
            "request_options": {"paginate": {"start": start, "rows": 10000}}, "return_type": "entry"}
        req = urllib.request.Request("https://search.rcsb.org/rcsbsearch/v2/query", data=json.dumps(q).encode(),
                                     headers={**UA, "Content-Type": "application/json"})
        d = json.loads(urllib.request.urlopen(req, timeout=120).read())
        page = [r["identifier"] for r in d.get("result_set", [])]
        ids += page
        start += len(page)
        if not page or start >= d.get("total_count", 0):
            return ids


def ca_by_chain(pid):
    url = f"https://models.rcsb.org/v1/{pid}/atoms?label_atom_id=CA&model_nums=1&encoding=cif"
    txt = fetch(url, f"{pid}_ca.cif")
    chains, cols, in_loop = {}, [], False
    for line in txt.splitlines():
        if line.startswith("_atom_site."):
            cols.append(line.split(".")[1].strip()); in_loop = True; continue
        if in_loop and line and not line.startswith(("_", "#", "loop_")):
            parts = line.split()
            if len(parts) == len(cols):
                d = dict(zip(cols, parts))
                chains.setdefault(d["label_asym_id"], []).append([float(d["Cartn_x"]), float(d["Cartn_y"]), float(d["Cartn_z"])])
        elif in_loop and line.startswith("#"):
            in_loop = False
    return {k: np.array(v) for k, v in chains.items()}


def voronoi_faces(points_frac, G, O, reach=2):
    """Faces of the Voronoi cells of the given fractional points under the group G, cell matrix O (frac→cart)."""
    imgs, owner = [], []
    for k, x in enumerate(points_frac):
        for R, t in zip(G.R, G.t):
            imgs.append((R @ x + t) % 1.0)
            owner.append(k)
    imgs = np.array(imgs)
    # merge images that coincide (points on special positions)
    keep, seen = [], []
    for i, p in enumerate(imgs):
        if not any(np.allclose((p - q + 0.5) % 1.0 - 0.5, 0, atol=1e-6) for q in seen):
            seen.append(p); keep.append(i)
    imgs, owner = imgs[keep], np.array(owner)[keep]
    shifts = np.array([(i, j, k) for i in range(-reach, reach + 1) for j in range(-reach, reach + 1) for k in range(-reach, reach + 1)])
    P = (imgs[None, :, :] + shifts[:, None, :]).reshape(-1, 3) @ O.T
    centre = np.flatnonzero((shifts == 0).all(1))[0] * len(imgs)
    vor = Voronoi(P)
    faces = np.zeros(len(imgs), int)
    for (p, q) in vor.ridge_points:
        for a in (p, q):
            if centre <= a < centre + len(imgs):
                faces[a - centre] += 1
    # one value per original point (its first image)
    first = {o: i for i, o in reversed(list(enumerate(owner)))}
    return [int(faces[first[k]]) for k in range(len(points_frac))]


def analyse(pid):
    info = cell_info(pid)
    cl = info["cell"]
    sg = info["symmetry"]["space_group_name_H_M"]
    M = frac_matrix(cl["length_a"], cl["length_b"], cl["length_c"], cl["angle_alpha"], cl["angle_beta"], cl["angle_gamma"])
    O = np.linalg.inv(M)
    G = Group(sg)
    ch = {k: v for k, v in ca_by_chain(pid).items() if len(v) >= MIN_CA}
    if not ch:
        raise ValueError("no protein chain")
    cents = np.array([v.mean(0) for v in ch.values()])
    allca = np.vstack(list(ch.values()))
    f_chain = voronoi_faces([(M @ c) % 1.0 for c in cents], G, O)
    f_asu = voronoi_faces([(M @ allca.mean(0)) % 1.0], G, O)[0]
    return {"id": pid, "sg": sg, "n_chains": len(ch), "chain_faces": f_chain, "asu_faces": f_asu,
            "chain_len": [len(v) for v in ch.values()], "cell": [cl[k] for k in ("length_a", "length_b", "length_c", "angle_alpha", "angle_beta", "angle_gamma")],
            "resolution": (info.get("rcsb_entry_info") or {}).get("resolution_combined")}


if __name__ == "__main__":
    sg = sys.argv[1]
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 1500
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    ids = all_entries(sg)
    random.Random(seed).shuffle(ids)
    ids = ids[:n]
    print(f"{sg}: sample of {len(ids)}", flush=True)
    out = []
    for i, pid in enumerate(ids):
        try:
            out.append(analyse(pid))
        except Exception as e:  # noqa: BLE001
            out.append({"id": pid, "error": str(e)[:80]})
        if (i + 1) % 100 == 0:
            print(f"  {i + 1}/{len(ids)}", flush=True)
    tag = sg.replace(" ", "")
    (HERE / "data" / f"chains_{tag}.json").write_text(json.dumps(out))
    ok = [r for r in out if "chain_faces" in r]
    print("asu faces:", sorted(Counter(r["asu_faces"] for r in ok).items()))
    print("chain faces:", sorted(Counter(f for r in ok for f in r["chain_faces"]).items()))
    print("errors:", sum(1 for r in out if "error" in r))
