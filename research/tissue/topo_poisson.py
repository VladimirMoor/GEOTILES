"""Базовая линия для топологии: ячейки Вороного случайных точек (Пуассона) → data/topo_poisson.json."""
import json

import numpy as np

from topology import HERE, classify, face_graphs_from_delaunay

rng = np.random.default_rng(11)
cells = []
for rep in range(4):
    X = rng.random((12000, 3))
    inner = set(np.where(np.all((X > 0.12) & (X < 0.88), axis=1))[0].tolist())
    cells += classify(face_graphs_from_delaunay(X, inner))
(HERE / "data" / "topo_poisson.json").write_text(json.dumps({"cells": cells}))
ok = [c for c in cells if c]
print("poisson cells", len(cells), "resolved", len(ok), "sym>1", sum(c["sym"] > 1 for c in ok), "distinct", len({c["code"] for c in ok}))
