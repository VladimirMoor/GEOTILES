"""Шаг 4. Гибкие семейства: что меняется, что постоянно.

Для каждого семейства по выборке членов находим:
- постоянные двугранные углы (одинаковые у всех членов);
- тождественные соотношения «α_e + α_e' = 180°» (пары рёбер, углы которых дополнительны у всех членов);
- грани: какие из них правильные треугольники, ромбы, пятиугольники и т. д.;
- классы параллельных рёбер (зоны): рёбра, параллельные у всех членов.
"""
import json
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from realize import edges_of, verify  # noqa: E402


def main(src="realizations.json"):
    types = {t["id"]: t for t in json.loads((ROOT / "data" / "types.json").read_text())}
    R = json.loads((ROOT / "data" / src).read_text())
    report = {}
    for tid, rec in R.items():
        fam = [s for s in rec["solutions"] if s["flex"] > 0]
        if not fam:
            continue
        t = types[tid]
        E = edges_of(t["faces"])
        D = np.array([s["dihedrals"] for s in fam])  # члены × рёбра
        # отбрасываем почти вырожденные члены (угол < 2° или > 178°)
        keep = (D.min(axis=1) > 2) & (D.max(axis=1) < 178)
        D = D[keep]
        Xs = [np.array(s["X"]) for s, k in zip(fam, keep) if k]
        if len(D) < 2:
            report[tid] = {"dims": sorted({s["flex"] for s in fam}), "members": int(keep.sum()), "degenerate": True,
                           "minAngle": float(np.min([min(s["dihedrals"]) for s in fam]))}
            print(f"{tid:8} only near-degenerate samples ({len(fam)}), min dihedral {report[tid]['minAngle']:.3f}°")
            continue
        const = [i for i in range(len(E)) if np.ptp(D[:, i]) < 1e-4]
        supp = [(i, j) for i in range(len(E)) for j in range(i + 1, len(E))
                if np.ptp(D[:, i] + D[:, j]) < 1e-4 and abs(D[0, i] + D[0, j] - 180) < 1e-4]
        # зоны: рёбра, параллельные у всех членов
        par = []
        for i in range(len(E)):
            for j in range(i + 1, len(E)):
                ok = True
                for X in Xs:
                    a = X[E[i][1]] - X[E[i][0]]
                    b = X[E[j][1]] - X[E[j][0]]
                    if np.linalg.norm(np.cross(a, b)) > 1e-6:
                        ok = False
                        break
                if ok:
                    par.append((i, j))
        # компоненты параллельности
        parent = list(range(len(E)))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x
        for i, j in par:
            parent[find(i)] = find(j)
        zones = {}
        for i in range(len(E)):
            zones.setdefault(find(i), []).append(i)
        zone_sizes = sorted((len(z) for z in zones.values()), reverse=True)
        face_kinds = []
        for f, ang in zip(t["faces"], fam[0]["faceAngles"]):
            k = len(f)
            if k == 3:
                face_kinds.append("T")
            elif k == 4:
                face_kinds.append("R")
            else:
                face_kinds.append(f"{k}")
        report[tid] = {
            "dims": sorted({s["flex"] for s in fam}),
            "members": int(keep.sum()),
            "faces": "".join(sorted(face_kinds)),
            "vertexDegrees": t["vertexDegrees"],
            "constAngles": sorted({round(float(D[0, i]), 4) for i in const}),
            "nConstEdges": len(const), "nEdges": len(E),
            "supplementaryPairs": len(supp),
            "zones": zone_sizes,
        }
    (ROOT / "data" / "families.json").write_text(json.dumps(report, indent=1))
    for tid, r in report.items():
        if r.get("degenerate"):
            continue
        print(f"{tid:8} dim {r['dims']} members {r['members']:4}  faces {r['faces']:12} deg {r['vertexDegrees']}  "
              f"const {r['nConstEdges']}/{r['nEdges']} {r['constAngles']}  suppPairs {r['supplementaryPairs']}  zones {r['zones']}")


if __name__ == "__main__":
    main(*sys.argv[1:])
