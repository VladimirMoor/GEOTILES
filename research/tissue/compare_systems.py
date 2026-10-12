"""Сравнение полной топологии клеток трёх систем одним конвейером и при одинаковой детализации:
ткань (семязачатки PlantSeg), рост зёрен (модель Поттса — модель огрубления пены и зёрен металла),
воксельный Вороной случайных точек; плюс точный Вороной случайных точек для справки.

Для каждой системы: доля восстановленных клеток, среднее и разброс числа граней F; при фиксированном F —
доля клеток с нетривиальной симметрией и число разных типов в случайной подвыборке K клеток.
Выход: data/compare_systems.json и таблица в консоли.
"""
import glob
import json
import pathlib
from collections import Counter

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
D = HERE / "data"
K = 100


def load(pattern, key="cells"):
    cells = []
    for f in sorted(glob.glob(str(D / pattern))):
        cells += json.load(open(f))[key]
    return cells


SETS = {
    "tissue": load("topo/N_*.json"),
    "potts_small": load("topo_potts/potts_N180_*_t700.json") + load("topo_potts/potts_N180_*_t1000.json") + load("topo_potts/potts_N180_*_t400.json"),
    "potts_big": load("topo_potts/potts_N300_*_t1400.json") + load("topo_potts/potts_N300_*_t1800.json") + load("topo_potts/potts_N300_*_t2200.json"),
    "voxel_PV": load("topo_potts/pv_N300_n900_*.json"),
    "exact_PV": json.load(open(D / "topo_poisson.json"))["cells"],
}
rng = np.random.default_rng(1)
res = {}
for name, cells in SETS.items():
    ok = [c for c in cells if c]
    if not ok:
        continue
    F = np.array([c["F"] for c in ok])
    r = {"n": len(cells), "resolved": len(ok), "mean_F": float(F.mean()), "sd_F": float(F.std()),
         "sym": float(np.mean([c["sym"] > 1 for c in ok])), "byF": {}}
    for f in range(8, 21):
        c = [x for x in ok if x["F"] == f]
        if len(c) < K:
            continue
        sub = [c[i] for i in rng.choice(len(c), K, replace=False)]
        cnt = Counter(x["code"] for x in c)
        r["byF"][f] = {"n": len(c), "sym": float(np.mean([x["sym"] > 1 for x in c])),
                       "distinct_K": len({x["code"] for x in sub}), "top_share": max(cnt.values()) / len(c)}
    res[name] = r
(D / "compare_systems.json").write_text(json.dumps(res, indent=1))
names = list(res)
print(f"{'':10s}" + "".join(f"{n:>22s}" for n in names))
print(f"{'cells':10s}" + "".join(f"{res[n]['resolved']:>14d}/{res[n]['n']:<7d}" for n in names))
print(f"{'mean F':10s}" + "".join(f"{res[n]['mean_F']:>16.2f} ±{res[n]['sd_F']:4.2f}" for n in names))
print(f"{'sym all':10s}" + "".join(f"{res[n]['sym']:>22.3f}" for n in names))
for f in range(8, 21):
    row = []
    for n in names:
        b = res[n]["byF"].get(f) if res[n]["byF"] else None
        b = res[n]["byF"].get(f, None)
        row.append(f"{b['sym']:>12.2f} d={b['distinct_K']:3d} n={b['n']:<5d}"[-22:] if b else f"{'—':>22s}")
    print(f"F={f:<8d}" + "".join(f"{x:>22s}" for x in row))
