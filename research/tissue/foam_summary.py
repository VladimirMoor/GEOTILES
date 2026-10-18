"""Summary of real-foam bubble topology (data/topo_foam/foam_step102.json) for RESULTS.md.

Face-count distribution, Matzke types (numbers of 4-, 5-, 6-gons), most common Weinberg types and the share of
the ten most common, symmetry at fixed F; real bubbles vs the Voronoi cells of their centroids (same bubbles,
deep ones only); dependence on the distance from the rotation axis (shear rate grows linearly with radius).
Output: console and data/topo_foam/summary.json
"""
import json
import pathlib
import sys
from collections import Counter

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
# usage: python foam_summary.py [topology json] [axis centre, voxels] [sample radius, voxels] [voxel, µm]
SRC = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "data" / "topo_foam" / "foam_step102.json"
D = json.load(open(SRC))
CEN, R_MASK, VOX = (float(x) for x in sys.argv[2:5]) if len(sys.argv) > 4 else (1008.0, 940.0, 2.75)


def stats(cells):
    ok = [c for c in cells if c]
    if not ok:
        return {"n": 0}
    F = np.array([c["F"] for c in ok])
    codes = Counter(c["code"] for c in ok)
    abc = Counter("-".join(map(str, (c["p"] + [0] * 4)[1:4])) if not c["p"][0] and sum(c["p"][4:]) == 0 else "other" for c in ok)
    top = codes.most_common(10)
    sym_code = {c["code"]: c["sym"] for c in ok}
    p_code = {c["code"]: c["p"] for c in ok}
    return {"n": len(ok), "mean_F": float(F.mean()), "sd_F": float(F.std()),
            "F_dist": {int(k): round(v / len(ok), 4) for k, v in sorted(Counter(F.tolist()).items())},
            "sym": float(np.mean([c["sym"] > 1 for c in ok])),
            "sym_byF": {int(f): round(float(np.mean([c["sym"] > 1 for c in ok if c["F"] == f])), 3)
                        for f in sorted(set(F.tolist())) if (F == f).sum() >= 30},
            "types": len(codes), "top10_share": sum(v for _, v in top) / len(ok),
            "top10": [{"p": p_code[k], "sym": sym_code[k], "share": round(v / len(ok), 4)} for k, v in top],
            "matzke_abc": {k: round(v / len(ok), 4) for k, v in abc.most_common(12)},
            "mean_edges_per_face": float(np.mean([sum((i + 3) * x for i, x in enumerate(c["p"])) / c["F"] for c in ok]))}


real = stats(D["cells"])
deep = [c for c, d in zip(D["cells"], D["deep"]) if d]
pairs = [(c, v) for c, v, d in zip(D["cells"], D["voronoi"], D["deep"]) if d and c and v]
out = {"real_all": real, "real_deep": stats(deep), "voronoi_deep": stats([v for c, v, d in zip(D["cells"], D["voronoi"], D["deep"]) if d]),
       "same_type_as_voronoi": float(np.mean([c["code"] == v["code"] for c, v in pairs])) if pairs else None,
       "same_F_as_voronoi": float(np.mean([c["F"] == v["F"] for c, v in pairs])) if pairs else None}
Q = np.array(D["centroids_inner"])
rad = np.hypot(Q[:, 1] - CEN, Q[:, 2] - CEN) / R_MASK
out["by_radius"] = {}
for lo, hi in ((0, .33), (.33, .66), (.66, 1.0)):
    sel = [c for c, r in zip(D["cells"], rad) if lo <= r < hi]
    s = stats(sel)
    out["by_radius"][f"{lo}-{hi}"] = {k: s.get(k) for k in ("n", "mean_F", "sd_F", "sym", "top10_share")}
V = np.array(D["volumes_inner"])
out["volume_cv"] = float(V.std() / V.mean())
out["R_eq_um"] = float((3 * np.median(V) / 4 / np.pi) ** (1 / 3) * VOX)
F_all = np.array([c["F"] if c else np.nan for c in D["cells"]])
out["volume_by_F"] = {int(f): round(float(np.median(V[F_all == f]) / np.median(V)), 3) for f in range(11, 18) if (F_all == f).sum() >= 30}
(SRC.parent / (SRC.stem + "_summary.json")).write_text(json.dumps(out, indent=1))
for k in ("real_all", "real_deep", "voronoi_deep"):
    s = out[k]
    if not s["n"]:
        print(k, "empty")
        continue
    print(f"{k:13s} n={s['n']:6d} F={s['mean_F']:.2f}±{s['sd_F']:.2f} sym={s['sym']:.3f} types={s['types']} top10={s['top10_share']:.3f} "
          f"edges/face={s['mean_edges_per_face']:.3f}")
    print("   F:", s["F_dist"])
    print("   sym by F:", s["sym_byF"])
    print("   Matzke a-b-c:", s["matzke_abc"])
    print("   top types:", [(t["p"], t["sym"], t["share"]) for t in s["top10"][:6]])
print("same type as Voronoi of centroids:", out["same_type_as_voronoi"], " same F:", out["same_F_as_voronoi"])
print("by radius:", out["by_radius"])
print("volume cv", round(out["volume_cv"], 3), "R_eq um", round(out["R_eq_um"], 1), "volume/median by F:", out["volume_by_F"])
