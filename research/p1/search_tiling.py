"""Поиск разбиения вида «решётка» или «решётка + центральная симметрия» для конкретного тела.

Кандидаты решёточных векторов: разности вершин. Центр симметрии x ↦ w − x: w = X_i + X_j.
Отсев: |det L| = V (только сдвиги) или 2V (с симметрией); затем Монте-Карло проверка покрытия.

Запуск: python search_tiling.py F7-029 F8-249 F8-255 [--members 5]
"""
import argparse
import itertools
import json
import pathlib

import numpy as np

from tilecheck import Tile, coverage, crystal

ROOT = pathlib.Path(__file__).resolve().parent


def lattice_candidates(X):
    diffs = []
    for i, j in itertools.permutations(range(len(X)), 2):
        d = X[j] - X[i]
        if np.linalg.norm(d) > 1e-9 and not any(np.allclose(d, e, atol=1e-9) or np.allclose(d, -e, atol=1e-9) for e in diffs):
            diffs.append(d)
    return np.array(diffs)


def quick_check(tile, maps, samples=300, seed=0):
    try:
        c = coverage(tile, maps, samples=samples, R=2.6 * tile.radius + 1.5, seed=seed)
    except ValueError:
        return False
    return c["exactlyOne"] == c["samples"]


def search(X, max_found=3):
    X = np.asarray(X) - np.mean(X, axis=0)
    tile = Tile(X)
    V = tile.volume
    D = lattice_candidates(X)
    found = []
    I = np.eye(3)
    for (i, j, k) in itertools.combinations(range(len(D)), 3):
        det = abs(np.linalg.det(D[[i, j, k]]))
        for mode, need in (("translations", V), ("point-reflection", 2 * V)):
            if abs(det - need) > 1e-7 * need:
                continue
            L = D[[i, j, k]]
            if mode == "translations":
                options = [[(I, np.zeros(3))]]
            else:
                ws = {tuple(np.round(X[a] + X[b], 9)) for a in range(len(X)) for b in range(a, len(X))}
                options = [[(I, np.zeros(3)), (-I, np.array(w))] for w in ws]
            for motifs in options:
                maps = crystal(L, motifs, 1, R=2.6 * tile.radius + 1.5)
                if quick_check(tile, maps):
                    full = coverage(tile, crystal(L, motifs, 1, R=3 * tile.radius + 2), samples=4000, R=3 * tile.radius + 2)
                    if full["exactlyOne"] == full["samples"]:
                        found.append({"mode": mode, "lattice": L.tolist(),
                                      "w": motifs[1][1].tolist() if len(motifs) > 1 else None, "coverage": full})
                        if len(found) >= max_found:
                            return tile, found
    return tile, found


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ids", nargs="+")
    ap.add_argument("--members", type=int, default=4)
    args = ap.parse_args()
    R = json.loads((ROOT / "data" / "realizations.json").read_text())
    results = {}
    for tid in args.ids:
        fam = [s for s in R[tid]["solutions"] if min(s["dihedrals"]) > 8 and max(s["dihedrals"]) < 172]
        rng = np.random.default_rng(1)
        picks = rng.choice(len(fam), size=min(args.members, len(fam)), replace=False)
        results[tid] = []
        for p in picks:
            tile, found = search(fam[p]["X"], max_found=1)
            ok = bool(found)
            print(f"{tid} member {p}: V={tile.volume:.5f}  tiling found: {ok}" + (f"  ({found[0]['mode']}, coverage {found[0]['coverage']['exactlyOne']}/{found[0]['coverage']['samples']})" if ok else ""))
            results[tid].append({"member": int(p), "X": fam[p]["X"], "found": found})
    (ROOT / "data" / "tilings_found.json").write_text(json.dumps(results, default=float))


if __name__ == "__main__":
    main()
