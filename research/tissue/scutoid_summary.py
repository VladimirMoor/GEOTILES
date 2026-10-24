"""Summary of the ovule scutoid census (data/scutoid/*.json): pooled shares of cells in an apico-basal T1 exchange
(and of cells with any one-sided neighbour change), real layer vs the two column models, per wall-area threshold;
per-volume paired comparison; surface ratio (outer / inner area of the judged cells). Output: console and
data/scutoid/summary.json"""
import glob
import json
import pathlib

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
R = [json.load(open(f)) for f in sorted(glob.glob(str(HERE / "data" / "scutoid" / "N_*.json")))]
out = {"volumes": len(R), "cells": int(sum(r["by_threshold"]["1.0"]["n"] for r in R)), "pooled": {}, "paired": {}}
for thr in ("0.5", "1.0", "2.0"):
    row = {}
    for key in ("by_threshold", "model_real_axes", "model"):
        cells = [c for r in R for c in r[key][thr]["cells"]]
        row[key] = {m: float(np.mean([c[m] > 0 for c in cells])) for m in ("t1", "strict")}
        row[key]["n"] = len(cells)
    out["pooled"][thr] = row
    real = np.array([r["by_threshold"][thr]["share_t1"] for r in R])
    ax = np.array([r["model_real_axes"][thr]["share_t1"] for r in R])
    nm = np.array([r["model"][thr]["share_t1"] for r in R])
    out["paired"][thr] = {"real_below_real_axes": int((real < ax).sum()), "real_above_real_axes": int((real > ax).sum()),
                          "real_below_normal": int((real < nm).sum()), "real_above_normal": int((real > nm).sum())}
ratios = []
for r in R:
    c = r["by_threshold"]["1.0"]["cells"]
    a, b = sum(x["apical_area_um2"] for x in c), sum(x["basal_area_um2"] for x in c)
    if b > 0:
        ratios.append(a / b)
out["surface_ratio"] = {"median": float(np.median(ratios)), "q25": float(np.percentile(ratios, 25)), "q75": float(np.percentile(ratios, 75))}
out["thickness_um_median"] = float(np.median([x["thickness_um"] for r in R for x in r["by_threshold"]["1.0"]["cells"]]))
(HERE / "data" / "scutoid" / "summary.json").write_text(json.dumps(out, indent=1))
print(f"{out['volumes']} ovules, {out['cells']} judged epidermal cells; surface ratio outer/inner {out['surface_ratio']}; "
      f"median layer thickness {out['thickness_um_median']:.1f} µm")
for thr, row in out["pooled"].items():
    print(f"threshold {thr} µm²: " + "  ".join(f"{k}: T1 {v['t1']:.3f}, one-sided {v['strict']:.3f}" for k, v in row.items()),
          " paired:", out["paired"][thr])
