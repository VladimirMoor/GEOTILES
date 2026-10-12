"""Полная топология зёрен из снимков модели Поттса и из воксельного Вороного (тем же конвейером, что ткань).

Для каждого data/potts/*.npz: граф граней каждого внутреннего зерна (не касается граней куба) по тройкам
и четвёркам меток в кубиках 2×2×2, отсечение «щепок», код Вайнберга (topology.py). Выход: data/topo_potts/*.json
"""
import json
import pathlib
import sys

import numpy as np

from topology import classify, face_graphs_from_voxels

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "data" / "topo_potts"
OUT.mkdir(parents=True, exist_ok=True)

if __name__ == "__main__":
    for f in sorted((HERE / "data" / "potts").glob(sys.argv[1] if len(sys.argv) > 1 else "*.npz")):
        o = OUT / (f.stem + ".json")
        if o.exists():
            continue
        L = np.load(f)["L"]
        vg = face_graphs_from_voxels(L, {0})
        cells = classify(vg)
        fixed = classify(vg, fix=True)
        o.write_text(json.dumps({"file": f.name, "cells": cells, "cells_repaired": fixed,
                                 "faces_all": [g.number_of_nodes() for g in vg.values()]}))
        ok = [c for c in cells if c]
        print(f"{f.name}: inner {len(cells)}, resolved {len(ok)}, mean F {np.mean([c['F'] for c in ok]):.2f}, "
              f"sym>1 {np.mean([c['sym'] > 1 for c in ok]):.3f}", flush=True)
