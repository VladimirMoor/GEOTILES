"""Точная перепроверка результатов пакетного прогона и сравнение с максимумами Шмитта (2016).

Для каждой группы берём лучшие локальные точки из лога (координаты округлены до 1e-6), строим рациональные
точки: саму округлённую и несколько случайных рациональных сдвигов 1e-6…1e-4, — и сертифицируем число граней
точно (exact.certify). Итог по группе — максимум сертифицированного числа граней.
"""
import pathlib
import re
import sys
from fractions import Fraction as Fr

import gemmi
import numpy as np

from exact import certify

LOG = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].isdigit() else "data/batch")
SCHMITT = pathlib.Path("lit/schmitt.txt")


def schmitt_max():
    best = {}
    cur = []
    for line in SCHMITT.read_text().split("\n"):
        if "Space group type" in line:
            cur = [int(n) for n in re.findall(r"IT\((\d+)\)", line)]
            continue
        for f in re.findall(r"\((\d+), (\d+), (\d+)\)", line):
            for n in cur:
                best[n] = max(best.get(n, 0), int(f[2]))
    return best


def points(text, k=3):
    rows = []
    for m in re.finditer(r"→ (\d+) \(gap [^)]*\) at x=\[([^\]]*)\](?:, c/a=([\d.]+))?", text):
        x = [float(v) for v in m.group(2).split()]
        c = float(m.group(3)) if m.group(3) else 1.0
        rows.append((int(m.group(1)), x, c))
    rows.sort(key=lambda r: -r[0])
    return rows[:k]


def certify_point(hm, x, c, rng, tries=6):
    best = 0
    for t in range(tries):
        eps = 0 if t == 0 else 10.0 ** rng.uniform(-6, -4)
        xr = tuple(Fr(v + (rng.normal() * eps if t else 0)).limit_denominator(10 ** 7) for v in x)
        cr = Fr(c).limit_denominator(10 ** 7)
        try:
            nf, _, _ = certify(hm, xr, cr * cr)
            best = max(best, nf)
        except Exception:
            pass
    return best


if __name__ == "__main__":
    S = schmitt_max()
    rng = np.random.default_rng(0)
    only = [int(a) for a in sys.argv[1:] if a.isdigit()]
    for f in sorted(LOG.glob("g*.log")):
        n = int(f.stem[1:])
        if only and n not in only:
            continue
        text = f.read_text()
        hm = gemmi.find_spacegroup_by_number(n).hm
        pts = points(text)
        if not pts:
            continue
        cert = max(certify_point(hm, x, c, rng) for _, x, c in pts)
        flt = pts[0][0]
        s = S.get(n)
        flag = ""
        if s is not None and cert > s:
            flag = "  ← above Schmitt"
        if cert >= 39:
            flag += "  ★★★ ≥ 39"
        print(f"IT({n:3d}) {hm:14s} float {flt:3d}  certified {cert:3d}  Schmitt {s}{flag}", flush=True)
