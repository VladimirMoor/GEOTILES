"""Сводка по COD: максимумы граней по группам, сравнение со случайными точками и точная проверка лидеров.

Для каждой группы:
- позиции атомов делятся на общие (орбита полной длины) и частные;
- базовая линия: случайные точки x при c/a, взятых у тех же записей COD, — доля точек с ≥ k гранями,
  отсюда ожидаемое число общих позиций с ≥ k гранями при том же их количестве;
- до 3 лучших общих позиций проверяются точно (research/p4/exact.py).
Результат: data/summary.json.
"""
import json
import pathlib
import sys
from collections import Counter
from fractions import Fraction as Fr

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "p4"))
from dvfast import Group, Evaluator  # noqa: E402
from exact import certify  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
GROUPS = [("I 41 3 2", 38), ("I 41 2 2", 35), ("P 61 2 2", 34), ("P 65 2 2", 34), ("I -4 2 d", 33),
          ("P 43 3 2", 29), ("P 41 3 2", 29), ("I 41 c d", 24), ("I 21 3", 24), ("P 21 3", 24),
          ("P 41 21 2", 29), ("P 43 21 2", 29), ("I 41", 26), ("P 31 2 1", 25)]
N_RANDOM = 3000


def is_general(G, x):
    imgs = (np.einsum("nij,j->ni", G.R, x) + G.t) % 1.0
    imgs = np.round(imgs, 4) % 1.0
    return len({tuple(v) for v in imgs}) == len(G.R)


out = {}
for hm, record in GROUPS:
    tag = hm.replace(" ", "")
    G = Group(hm)
    data = [r for r in json.loads((HERE / "data" / f"cod_{tag}.json").read_text()) if r.get("sites")]
    gen, spec = [], []
    for r in data:
        for s in r["sites"]:
            if not s["facets"]:
                continue
            (gen if is_general(G, np.array(s["x"])) else spec).append((s["facets"], r["id"], s, r["c_over_a"]))
    # базовая линия
    rng = np.random.default_rng(7)
    cs = [r["c_over_a"] for r in data]
    picks = Counter(round(cs[i], 5) for i in rng.integers(0, len(cs), N_RANDOM))
    rand = Counter()
    for c, k in picks.items():
        E = Evaluator(G, c)
        for _ in range(k):
            try:
                rand[E.facets(rng.random(3))[0]] += 1
            except Exception:
                pass
    nr = sum(rand.values())
    kmax = max(f for f, *_ in gen)
    p_ge = sum(v for f, v in rand.items() if f >= kmax) / nr
    # точная проверка лучших общих позиций
    top = sorted(gen, key=lambda t: -t[0])[:3]
    cert = []
    for f, cid, s, c in top:
        nf = certify(hm, tuple(Fr(str(v)) for v in s["x"]), Fr(str(c)) ** 2)[0]
        cert.append({"id": cid, "label": s["label"], "el": s["el"], "float": f, "exact": nf, "occ": s["occ"], "c_over_a": c})
    out[hm] = {"record": record, "entries": len(data), "general_sites": len(gen), "special_sites": len(spec),
               "max_general": kmax, "max_special": max((f for f, *_ in spec), default=None),
               "hist_general": sorted(Counter(f for f, *_ in gen).items()),
               "random_hist": sorted(rand.items()), "random_n": nr, "p_random_ge_max": p_ge,
               "observed_ge_max": sum(1 for f, *_ in gen if f >= kmax), "expected_ge_max": p_ge * len(gen),
               "mean_general": float(np.mean([f for f, *_ in gen])),
               "mean_random": sum(f * v for f, v in rand.items()) / nr, "certified_top": cert}
    o = out[hm]
    print(f"{hm:10s} rec {record:2d} | entries {o['entries']:4d} gen {o['general_sites']:5d} | max {kmax:2d} "
          f"(exact {[c['exact'] for c in cert]}) | obs≥max {o['observed_ge_max']:3d} exp {o['expected_ge_max']:6.2f} "
          f"| mean COD {o['mean_general']:.2f} random {o['mean_random']:.2f}", flush=True)
(HERE / "data" / "summary.json").write_text(json.dumps(out, indent=1))
