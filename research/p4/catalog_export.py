"""Экспорт рекордных DV-стереоэдров в каталог сайта (docs/js/stereo-data.js).

Для каждого: вершины ячейки (декартовы, a = 1), решётка (столбцы B) и операции группы в декартовых
координатах (R, t). Число граней сверяется с точным сертификатом.
"""
import json
import pathlib
from fractions import Fraction as Fr

import numpy as np

from dvfast import Group, Evaluator
from exact import certify

ITEMS = [
    ("engel38", "I 41 3 2", (Fr(427, 6984), Fr(761, 6984), Fr(1421, 6984)), Fr(1)),
    ("schmitt35", "I 41 2 2", (Fr(62, 125), Fr(41, 125), Fr(79, 1000)), Fr(727, 500)),
]
cert = json.loads(pathlib.Path("certificates.json").read_text())
for n, key in ((122, "ours33"), (88, "ours29")):
    r = cert[str(n)]
    ITEMS.append((key, r["group"], tuple(Fr(v) for v in r["x"]), Fr(r["c_over_a"])))
# ячейки белковых молекул из PDB (research/bio): точка — центр асимметричной единицы, округлённый как в данных
for key, tag, pid in (("pdb4ux6", "P6122", "4UX6"), ("pdb8s97", "I4122", "8S97"), ("pdb1jky", "I4132", "1JKY")):
    r = next(r for r in json.loads(pathlib.Path(f"../bio/data/dv_{tag}.json").read_text()) if r["id"] == pid)
    ITEMS.append((key, r["sg"], tuple(Fr(str(v)) for v in r["x"]), Fr(str(r["c_over_a"]))))

out = {}
for key, hm, x, c in ITEMS:
    nf_exact, _, _ = certify(hm, x, c * c)
    G = Group(hm)
    E = Evaluator(G, float(c))
    X, A, b, V = E.cell_of(np.array([float(v) for v in x]))
    Vu = []
    for v in V:
        if not any(np.linalg.norm(v - u) < 1e-9 for u in Vu):
            Vu.append(v)
    Vu = np.array(Vu)
    B = E.B
    Binv = np.linalg.inv(B)
    ops = []
    for R, t in zip(G.R, G.t):
        Rc = B @ R @ Binv
        tc = B @ t
        ops.append([np.round(Rc, 12).tolist(), np.round(tc, 12).tolist()])
    nf_float, _, _ = E.facets(np.array([float(v) for v in x]))
    # точная структура граней: для каждой грани — вершины в циклическом порядке
    C0 = Vu.mean(axis=0)
    faces, normals, offs = [], [], []
    nA = np.linalg.norm(A, axis=1)
    for i in range(len(A)):
        on = np.where(np.abs(Vu @ A[i] - b[i]) / nA[i] < 1e-10)[0]
        if len(on) < 3 or np.linalg.matrix_rank(Vu[on] - Vu[on][0], tol=1e-11) < 2:
            continue
        nrm = A[i] / nA[i]
        c = Vu[on].mean(axis=0)
        u = Vu[on[0]] - c; u /= np.linalg.norm(u)
        w = np.cross(nrm, u)
        order = on[np.argsort(np.arctan2((Vu[on] - c) @ w, (Vu[on] - c) @ u))]
        faces.append([int(k) for k in order])
        normals.append(np.round(nrm, 12).tolist())
        offs.append(float(nrm @ Vu[order[0]]))
    assert len(faces) == nf_exact, (len(faces), nf_exact)
    vol_cell = abs(np.linalg.det(B)) / len(G.R)
    out[key] = {"group": hm, "number": G.number, "facets": nf_exact, "facets_float": nf_float,
                "x": [str(v) for v in x], "c_over_a": str(c),
                "points": np.round(Vu, 12).tolist(), "faces": faces, "normals": normals, "offsets": offs,
                "lattice": np.round(B.T, 12).tolist(), "ops": ops, "volume": vol_cell}
    print(key, hm, "facets exact", nf_exact, "float", nf_float, "vertices", len(Vu), "ops", len(ops))

js = "// Сгенерировано research/p4/catalog_export.py: рекордные DV-стереоэдры (вершины, решётка, операции группы).\n"
js += "export const STEREO = " + json.dumps(out) + ";\n"
pathlib.Path("../../docs/js/stereo-data.js").write_text(js)
