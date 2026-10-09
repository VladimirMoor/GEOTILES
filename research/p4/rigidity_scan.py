"""Скан жёсткости: у скольких DV-ячеек есть не-DV деформации первого порядка (точно)."""
import re, sys, io, contextlib
from fractions import Fraction as Fr
import numpy as np, gemmi
from deform_exact import analyse

rng = np.random.default_rng(7)
pts = []
for n in [214, 98, 122, 88, 92, 96, 80, 178, 152, 155, 212, 199, 198]:
    hm = gemmi.find_spacegroup_by_number(n).hm
    log = open(f"data/batch/g{n:03d}.log").read()
    best = sorted(((int(m.group(1)), m.group(2), m.group(3)) for m in re.finditer(r"→ (\d+) \(gap [^)]*\) at x=\[([^\]]*)\](?:, c/a=([\d.]+))?", log)), reverse=True)[:2]
    for nf, xs, cs in best:
        pts.append((hm, [float(v) for v in xs.split()], float(cs) if cs else 1.0))
    for _ in range(2):
        pts.append((hm, list(rng.random(3)), 1.0 if n >= 195 else float(rng.uniform(0.6, 2.5))))
tot = 0
for hm, x, c in pts:
    xr = tuple(Fr(v).limit_denominator(5000) for v in x); cr = Fr(c).limit_denominator(1000)
    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf):
            extra = analyse(hm, xr, cr * cr)
        line = [l for l in buf.getvalue().splitlines() if "facets" in l][0]
        nf = int(re.search(r"facets (\d+)", line).group(1))
        print(f"{hm:12s} facets {nf:2d}  non-DV directions: {extra}", flush=True)
        tot += extra > 0
    except Exception as e:
        print(f"{hm:12s} skipped ({type(e).__name__}: {str(e)[:50]})", flush=True)
print("cells with non-DV deformations:", tot)
