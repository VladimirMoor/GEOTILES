"""Тождественные соотношения углов в семействе: Σ kᵢθᵢ = 360° или 180°, верные у всех членов."""
import itertools
import json
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent


def identities(D, K=6):
    nE = D.shape[1]
    groups = []
    for e in range(nE):
        for g in groups:
            if np.max(np.abs(D[:, g[0]] - D[:, e])) < 1e-6:
                g.append(e)
                break
        else:
            groups.append([e])
    F = np.array([D[:, g[0]] for g in groups]).T
    m = F.shape[1]
    covered, rels = set(), []
    for total in range(1, K + 1):
        for combo in itertools.combinations_with_replacement(range(m), total):
            k = np.bincount(combo, minlength=m)
            v = F @ k
            for tg in (360, 180):
                if np.max(np.abs(v - tg)) < 1e-6:
                    rels.append((k.tolist(), tg))
                    covered |= set(np.nonzero(k)[0].tolist())
    return groups, F, covered, rels


if __name__ == "__main__":
    R = json.loads((ROOT / "data" / "realizations.json").read_text())
    for tid in sys.argv[1:]:
        fam = [s for s in R[tid]["solutions"] if s["flex"] > 0 and min(s["dihedrals"]) > 5 and max(s["dihedrals"]) < 175]
        D = np.array([s["dihedrals"] for s in fam])
        groups, F, cov, rels = identities(D)
        free = [i for i in range(F.shape[1]) if i not in cov]
        print(f"{tid}: {len(fam)} members, {F.shape[1]} angle functions, edges per function {[len(g) for g in groups]}")
        print(f"   identities: {rels[:12]}")
        print(f"   functions NOT covered by any identity: {free}  (edge counts {[len(groups[i]) for i in free]})")
