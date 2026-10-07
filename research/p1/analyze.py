"""Шаг 3. Сводка реализаций: жёсткие тела + гибкие семейства, фильтр двугранных углов.

Фильтр (необходимое условие монотайлинга): каждый двугранный угол α входит в сумму
Σ kᵢαᵢ = 360° или Σ kᵢαᵢ = 180° (неотрицательные целые kᵢ, α с kα ≥ 1).
"""
import json
import pathlib
from collections import Counter

ROOT = pathlib.Path(__file__).resolve().parent


def distinct(angles, tol=1e-5):
    out = []
    for a in sorted(angles):
        if not out or abs(a - out[-1]) > tol:
            out.append(a)
    return out


def combos(A, target, tol=1e-5, limit=2000):
    A = sorted(A, reverse=True)
    sols, k = [], [0] * len(A)

    def rec(i, rest):
        if len(sols) >= limit:
            return
        if abs(rest) < tol:
            sols.append(list(k))
            return
        if i >= len(A) or rest < -tol:
            return
        for m in range(int((rest + tol) // A[i]), -1, -1):
            k[i] = m
            rec(i + 1, rest - m * A[i])
        k[i] = 0

    rec(0, target)
    return A, [s for s in sols if any(s)]


def dihedral_filter(angles):
    A = distinct(angles)
    _, s360 = combos(A, 360)
    A_sorted, s180 = combos(A, 180)
    bad = []
    for i, a in enumerate(A_sorted):
        if not any(s[i] for s in s360) and not any(s[i] for s in s180):
            bad.append(a)
    return {"ok": not bad, "failing": sorted(bad), "angles": A}


def face_word(t):
    c = Counter(len(f) for f in t["faces"])
    return " + ".join(f"{c[k]}×{k}" for k in sorted(c, reverse=True))


def main(src="realizations.json"):
    types = {t["id"]: t for t in json.loads((ROOT / "data" / "types.json").read_text())}
    R = json.loads((ROOT / "data" / src).read_text())
    rows = []
    for tid, rec in R.items():
        sols = rec["solutions"]
        if not sols:
            continue
        t = types[tid]
        flex = Counter(s["flex"] for s in sols)
        rigid = [s for s in sols if s["flex"] == 0]
        fam = [s for s in sols if s["flex"] > 0]
        entry = {"id": tid, "faces": face_word(t), "nV": t["nV"], "nE": t["nE"], "rigid": [], "families": {}}
        for s in rigid:
            f = dihedral_filter(s["dihedrals"])
            entry["rigid"].append({"hits": s["hits"], "angles": [round(a, 4) for a in f["angles"]], "filter": f["ok"], "failing": [round(a, 4) for a in f["failing"]], "X": s["X"]})
        if fam:
            # для семейств поэлементный фильтр не показателен (углы меняются непрерывно) — разбираем отдельно
            dims = Counter(s["flex"] for s in fam)
            entry["families"] = {"dims": dict(dims), "samples": len(fam), "example": fam[0]["X"],
                                 "faceAngleRange": [round(min(min(a) for a in s["faceAngles"]), 2) for s in fam[:1]]}
        entry["stats"] = rec["stats"]
        rows.append(entry)
    (ROOT / "data" / "summary.json").write_text(json.dumps(rows, indent=1))
    print(f"{len(rows)} of {len(R)} types have equilateral convex realizations\n")
    print(f"{'type':9} {'faces':18} {'V':>3} {'E':>3}  realizations")
    for e in rows:
        parts = []
        for r in e["rigid"]:
            mark = "PASS" if r["filter"] else "fail"
            parts.append(f"rigid[{r['hits']} hits] {mark} {r['angles']}")
        if e["families"]:
            f = e["families"]
            parts.append(f"family dim {f['dims']} ({f['samples']} samples)")
        print(f"{e['id']:9} {e['faces']:18} {e['nV']:>3} {e['nE']:>3}  " + " | ".join(parts))


if __name__ == "__main__":
    import sys
    main(*sys.argv[1:])
