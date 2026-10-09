"""Сводная таблица по множеству 𝒜 (40 тетраэдров) → data/a_summary.json."""
import json, re, pathlib
from tetra import PRINTED_A, EXCLUDED_BY_LP, Tet, EDGES
from stars import Lengths, edge_star_exists
from prop27 import pattern, parallelogram_like
from fractions import Fraction as Fr

clus = {}
for line in open("data_clusters7.log"):
    m = re.match(r"#(\d+): convex clusters (\d+)", line)
    if m:
        clus[int(m.group(1))] = int(m.group(2))
rows = []
for idx, a in enumerate(PRINTED_A, 1):
    if a in EXCLUDED_BY_LP:
        continue
    T = Tet(a); Lc = Lengths(T)
    bad = [f"{i}{j}" for (i, j) in EDGES if not edge_star_exists(T, (i, j), Lc)[0]]
    c = pattern(Lc)
    covered = (len(set(c)) == 6 or parallelogram_like(c)) and any((Fr(2) / x).denominator != 1 for x in a)
    status = {1: "Sommerville No. 3 (tiles)", 3: "Sommerville No. 2 (tiles)"}.get(idx, "no face-to-face tiling")
    rows.append({"no": idx, "angles": [str(x) for x in a], "status": status, "edges_without_star": bad,
                 "f2f_exclusion_new": (not covered) and idx not in (1, 3),
                 "convex_clusters_le7": clus.get(idx)})
pathlib.Path("data").mkdir(exist_ok=True)
json.dump(rows, open("data/a_summary.json", "w"), indent=1)
print(len(rows), "rows;", sum(r["f2f_exclusion_new"] for r in rows), "new f2f exclusions")
