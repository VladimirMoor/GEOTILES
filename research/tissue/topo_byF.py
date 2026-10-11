"""Сравнение при фиксированном числе граней F: доля симметричных и разнообразие типов (энтропия, число типов
в подвыборке одного размера) для ткани, Вороного тех же центров и случайных точек → data/topo_byF.json."""
import json
import pathlib
from collections import Counter

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
files = [json.loads(p.read_text()) for p in sorted((HERE / "data" / "topo").glob("*.json"))]
SETS = {"tissue": sum((f["cells"] for f in files), []), "tissue_repaired": sum((f["cells_repaired"] for f in files), []),
        "voronoi_centres": sum((f["voronoi"] for f in files), []),
        "poisson": json.loads((HERE / "data" / "topo_poisson.json").read_text())["cells"]}
rng = np.random.default_rng(5)
K = 150  # размер подвыборки при каждом F
out = {}
print(f"{'F':>3} | " + " | ".join(f"{k:>24s}" for k in SETS))
for F in range(8, 21):
    row = {}
    for k, cells in SETS.items():
        c = [x for x in cells if x and x["F"] == F]
        if len(c) < K:
            row[k] = None
            continue
        sub = [c[i] for i in rng.choice(len(c), K, replace=False)]
        cnt = Counter(x["code"] for x in c)
        p = np.array(list(cnt.values())) / len(c)
        row[k] = {"n": len(c), "sym": float(np.mean([x["sym"] > 1 for x in c])), "distinct_K": len({x["code"] for x in sub}),
                  "top_share": float(p.max())}
    out[F] = row
    print(f"{F:>3} | " + " | ".join(
        (f"n={r['n']:5d} sym={r['sym']:.2f} d={r['distinct_K']:3d}" if r else f"{'—':>24s}") for r in row.values()))
(HERE / "data" / "topo_byF.json").write_text(json.dumps(out))
