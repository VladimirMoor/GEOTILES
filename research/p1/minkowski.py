"""Суммы Минковского единичных многоугольников и отрезков: равносторонность, грани, замощение.

В общем положении каждое ребро суммы — сдвиг ребра одного из слагаемых, поэтому сумма единичных
многоугольников и отрезков равностороння. Инвариант Дена суммы равен нулю (рёбра одного направления
дают в сумме внешних углов π), так что фильтр Дена пройден автоматически.

Для случайных членов каждого семейства: строим сумму, проверяем рёбра и число граней, ищем разбиение
вида «решётка» или «решётка + центральная симметрия» (search_tiling.search) и проверяем покрытие.

Запуск: python minkowski.py [--members 8] [--kinds tri+tri,rhomb+tri,...]
"""
import argparse
import itertools
import json
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from tilecheck import Tile, coverage, crystal  # noqa: E402
from search_tiling import search  # noqa: E402


def rand_rot(rng):
    q = rng.normal(size=4)
    q /= np.linalg.norm(q)
    a, b, c, d = q
    return np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                     [2*(b*c+a*d), a*a-b*b+c*c-d*d, 2*(c*d-a*b)],
                     [2*(b*d-a*c), 2*(c*d+a*b), a*a-b*b-c*c+d*d]])


def regular(n):
    """Правильный n-угольник со стороной 1 в плоскости xy."""
    R = 1 / (2 * np.sin(np.pi / n))
    return np.array([[R * np.cos(2 * np.pi * k / n), R * np.sin(2 * np.pi * k / n), 0] for k in range(n)])


def polygon(kind, rng):
    if kind == "seg":
        return np.array([[0, 0, 0], [1, 0, 0.0]]) @ rand_rot(rng).T
    if kind == "tri":
        P = regular(3)
    elif kind == "rhomb":
        t = rng.uniform(np.radians(35), np.radians(145))
        P = np.array([[0, 0, 0], [1, 0, 0], [1 + np.cos(t), np.sin(t), 0], [np.cos(t), np.sin(t), 0]])
    elif kind == "hex":  # центрально-симметричный равносторонний шестиугольник
        a, b = rng.uniform(np.radians(30), np.radians(80)), rng.uniform(np.radians(90), np.radians(150))
        g = [np.array([1, 0, 0]), np.array([np.cos(a), np.sin(a), 0]), np.array([np.cos(b), np.sin(b), 0])]
        P = [np.zeros(3)]
        for v in g + [-x for x in g]:
            P.append(P[-1] + v)
        P = np.array(P[:-1])
    elif kind == "pent":  # равносторонний пятиугольник (случайный, выпуклый)
        from prisms import pentagon, convex_ok
        while True:
            P2 = pentagon(rng.uniform(60, 170), rng.uniform(60, 170))
            if P2 is not None and convex_ok(P2, margin=5):
                break
        P = np.c_[P2, np.zeros(5)]
    else:
        raise ValueError(kind)
    return (P - P.mean(axis=0)) @ rand_rot(rng).T


def msum(parts):
    pts = parts[0]
    for P in parts[1:]:
        pts = np.array([p + q for p in pts for q in P])
    return pts


def analyse(kind, members, rng):
    rows = []
    for m in range(members):
        parts = [polygon(k, rng) for k in kind.split("+")]
        pts = msum(parts)
        tile = Tile(pts)
        lens = tile.edge_lengths()
        nf = len(tile.planes)
        tile2, found = search(tile.vertices, max_found=1)
        rows.append({"faces": nf, "vertices": len(tile.vertices), "equilateral": bool(np.ptp(lens) < 1e-7),
                     "tiles": bool(found), "mode": found[0]["mode"] if found else None,
                     "coverage": found[0]["coverage"]["exactlyOne"] if found else None})
        print(f"  {kind} #{m}: F={nf} V={len(tile.vertices)} equilateral={rows[-1]['equilateral']} "
              f"tiling={'yes (' + found[0]['mode'] + ')' if found else 'not found'}", flush=True)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--members", type=int, default=8)
    ap.add_argument("--kinds", default="seg+tri,tri+tri,rhomb+tri,seg+seg+tri,pent+seg,tri+hex,tri+tri+seg,pent+tri")
    ap.add_argument("--seed", type=int, default=11)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    out = {}
    for kind in args.kinds.split(","):
        out[kind] = analyse(kind, args.members, rng)
    path = ROOT / "data" / "minkowski.json"
    old = json.loads(path.read_text()) if path.exists() else {}
    old.update(out)
    path.write_text(json.dumps(old, indent=1))


if __name__ == "__main__":
    main()
