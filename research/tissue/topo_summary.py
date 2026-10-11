"""Сводка по полной топологии клеток: ткань (строго и с починкой), Вороной тех же центров, Вороной случайных точек.

Метрики:
- доля разрешённых клеток и смещение по числу граней (разрешённые vs все);
- доля клеток с нетривиальной симметрией — всего и с перевзвешиванием на распределение граней ткани;
- разнообразие типов: число разных типов в случайной подвыборке фиксированного размера;
- самые частые типы (p-вектор, симметрия) и их частота в каждом наборе.
Результат: data/topo_summary.json.
"""
import json
import pathlib
from collections import Counter

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
files = [json.loads(p.read_text()) for p in sorted((HERE / "data" / "topo").glob("*.json"))]
pois = json.loads((HERE / "data" / "topo_poisson.json").read_text())["cells"]

SETS = {
    "tissue": sum((f["cells"] for f in files), []),
    "tissue_repaired": sum((f["cells_repaired"] for f in files), []),
    "voronoi_centres": sum((f["voronoi"] for f in files), []),
    "poisson": pois,
}
faces_all = sum((f["faces_all"] for f in files), [])
N_SUB = 5000
rng = np.random.default_rng(3)


def stats(cells, ref_F=None):
    ok = [c for c in cells if c]
    F = Counter(c["F"] for c in ok)
    sym_by_F = {f: np.mean([c["sym"] > 1 for c in ok if c["F"] == f]) for f in F}
    out = {"n": len(cells), "resolved": len(ok), "mean_F": float(np.mean([c["F"] for c in ok])),
           "sym_frac": float(np.mean([c["sym"] > 1 for c in ok])),
           "sym_orders": sorted(Counter(c["sym"] for c in ok).items())}
    if ref_F:  # перевзвешивание на распределение числа граней ref_F
        tot = sum(v for f, v in ref_F.items() if f in sym_by_F)
        out["sym_frac_reweighted"] = float(sum(sym_by_F[f] * v for f, v in ref_F.items() if f in sym_by_F) / tot)
    sub = [ok[i] for i in rng.choice(len(ok), size=min(N_SUB, len(ok)), replace=False)]
    out["distinct_in_sub"] = len({c["code"] for c in sub})
    out["sub_size"] = len(sub)
    types = Counter(c["code"] for c in ok)
    out["top"] = [{"code": k, "count": v, "freq": v / len(ok),
                   "p": next(c["p"] for c in ok if c["code"] == k), "sym": next(c["sym"] for c in ok if c["code"] == k),
                   "F": next(c["F"] for c in ok if c["code"] == k)} for k, v in types.most_common(12)]
    out["freq_of"] = types
    return out


tF = Counter(c["F"] for c in SETS["tissue"] if c)
res = {k: stats(v, tF) for k, v in SETS.items()}
# частоты самых частых типов ткани в других наборах
for t in res["tissue"]["top"]:
    t["freq_in"] = {k: res[k]["freq_of"].get(t["code"], 0) / res[k]["resolved"] for k in SETS}
for k in res:
    del res[k]["freq_of"]
res["bias"] = {"mean_F_all_inner": float(np.mean(faces_all)), "n_all_inner": len(faces_all)}
(HERE / "data" / "topo_summary.json").write_text(json.dumps(res))

print(f"files {len(files)}; inner cells {len(faces_all)}, mean faces (graph, all inner) {np.mean(faces_all):.2f}")
for k in SETS:
    r = res[k]
    print(f"{k:16s} resolved {r['resolved']:6d}/{r['n']:6d}  mean F {r['mean_F']:.2f}  sym>1 {r['sym_frac']:.3f}  "
          f"(reweighted to tissue F: {r['sym_frac_reweighted']:.3f})  distinct in {r['sub_size']}: {r['distinct_in_sub']}")
print("top tissue types (freq in tissue / repaired / voronoi centres / poisson):")
for t in res["tissue"]["top"]:
    fi = t["freq_in"]
    print(f"  F={t['F']:2d} sym={t['sym']:2d} p={t['p']}  {fi['tissue']:.4f} {fi['tissue_repaired']:.4f} {fi['voronoi_centres']:.4f} {fi['poisson']:.4f}")
