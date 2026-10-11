"""Иллюстрация для страницы «Реальные кристаллы»: двумерный аналог стереоэдров.

Группа p3 (повороты на 120° и сдвиги шестиугольной решётки). Точка общего положения, её копии под группой
(маленькие «молекулы»-стрелки, повёрнутые вместе с точкой) и их ячейки Вороного: все конгруэнтны и
замощают плоскость. Одна ячейка выделена, её соседи по рёбрам подсвечены.
Выход: docs/img/orbit-voronoi.svg
"""
import pathlib

import numpy as np
from scipy.spatial import Voronoi

OUT = pathlib.Path(__file__).resolve().parents[2] / "docs" / "img" / "orbit-voronoi.svg"
W, H, S = 760, 380, 95.0  # размер картинки и масштаб (пикселей на период решётки)
a1, a2 = np.array([1.0, 0.0]), np.array([0.5, np.sqrt(3) / 2])
x0 = np.array([0.31, 0.13])
rot = lambda k: np.array([[np.cos(2 * np.pi * k / 3), -np.sin(2 * np.pi * k / 3)],  # noqa: E731
                          [np.sin(2 * np.pi * k / 3), np.cos(2 * np.pi * k / 3)]])

pts, angs = [], []
for i in range(-9, 10):
    for j in range(-9, 10):
        for k in range(3):
            pts.append(rot(k) @ x0 + i * a1 + j * a2)
            angs.append(120 * k)
pts, angs = np.array(pts), np.array(angs)
vor = Voronoi(pts)
cx, cy = W / 2, H / 2
P = lambda p: (cx + S * p[0], cy - S * p[1])  # noqa: E731

# выделенная ячейка — точка, ближайшая к центру
c0 = int(np.argmin(np.linalg.norm(pts, axis=1)))
nbrs = {int(b) for a, b in vor.ridge_points if a == c0} | {int(a) for a, b in vor.ridge_points if b == c0}

svg = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" font-family="Manrope, sans-serif">',
       f'<defs><clipPath id="c"><rect x="0" y="0" width="{W}" height="{H}" rx="14"/></clipPath>'
       '<path id="m" d="M -11 -6.5 L 11 0 L -11 6.5 L -6 0 Z"/></defs><g clip-path="url(#c)">']
for idx, reg in enumerate(vor.point_region):
    R = vor.regions[reg]
    if -1 in R or not R:
        continue
    poly = [P(vor.vertices[v]) for v in R]
    if all(x < -60 or x > W + 60 or y < -60 or y > H + 60 for x, y in poly):
        continue
    if idx == c0:
        fill, op = "#2f6fd0", 0.55
    elif idx in nbrs:
        fill, op = "#e3a33b", 0.42
    else:
        fill, op = "#8a909b", 0.06
    d = " ".join(f"{x:.1f},{y:.1f}" for x, y in poly)
    svg.append(f'<polygon points="{d}" fill="{fill}" fill-opacity="{op}" stroke="#8a909b" stroke-width="1.3"/>')
for p, ang in zip(pts, angs):
    x, y = P(p)
    if -20 < x < W + 20 and -20 < y < H + 20:
        svg.append(f'<use href="#m" transform="translate({x:.1f},{y:.1f}) rotate({-(ang + 20)})" fill="#9aa3b0" stroke="#5a6472" stroke-width="0.8"/>')
svg.append("</g></svg>")
OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text("\n".join(svg))
print(OUT, "neighbours of highlighted cell:", len(nbrs))
