"""Явные рациональные сертификаты новых нижних оценок: точка x, c/a и точное число граней.

Для каждой группы перебираем лучшие точки из логов прогонов (batch, batch_low, batch_high) и их
рациональные приближения/малые сдвиги; оставляем точку с наибольшим точно сертифицированным числом
граней. Результат — data/certificates.json и таблица в stdout.
"""
import json
import pathlib
import re
from fractions import Fraction as Fr

import gemmi
import numpy as np

from exact import certify

GROUPS = [88, 122, 110, 167, 86, 201, 228, 138, 184, 192]
LOGS = ["data/batch", "data/batch_low", "data/batch_high"]


def candidates(n):
    out = []
    for d in LOGS:
        f = pathlib.Path(d) / f"g{n:03d}.log"
        if not f.exists():
            continue
        for m in re.finditer(r"→ (\d+) \(gap [^)]*\) at x=\[([^\]]*)\](?:, c/a=([\d.]+))?", f.read_text()):
            out.append((int(m.group(1)), [float(v) for v in m.group(2).split()], float(m.group(3)) if m.group(3) else 1.0))
    out.sort(key=lambda r: -r[0])
    return out[:6]


def main():
    rng = np.random.default_rng(1)
    res = {}
    for n in GROUPS:
        hm = gemmi.find_spacegroup_by_number(n).hm
        best = None
        for _, x, c in candidates(n):
            for t in range(10):
                eps = 0 if t == 0 else 10.0 ** rng.uniform(-6, -4)
                xr = tuple(Fr(v + (rng.normal() * eps if t else 0)).limit_denominator(10 ** 6) for v in x)
                cr = Fr(c).limit_denominator(10 ** 4)
                try:
                    nf, nv, _ = certify(hm, xr, cr * cr)
                except Exception:
                    continue
                if best is None or nf > best[0]:
                    best = (nf, xr, cr, nv)
        nf, xr, cr, nv = best
        res[n] = {"group": hm, "facets": nf, "x": [str(v) for v in xr], "c_over_a": str(cr)}
        print(f"IT({n:3d}) {hm:12s} certified facets {nf:2d}  x = ({', '.join(map(str, xr))}),  c/a = {cr}", flush=True)
    pathlib.Path("data/certificates.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
