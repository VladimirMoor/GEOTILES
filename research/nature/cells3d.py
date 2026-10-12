"""Замкнутые ячейки для каталога: пчелиная и ячейка Тота как параллелоэдры.

Если продолжить соты вверх и вниз слоями, сдвинутыми каждый раз на −d, все ячейки становятся параллельными
копиями одной: ячейка = призма над шестиугольником между поверхностью f (дно) и f(p + d) + H (крышка).
Пчелиная получается удлинённым ромбододекаэдром (12 граней), ячейка Тота — комбинаторно усечённым
октаэдром (14 граней). Выход: docs/js/nature-data.js.
"""
import json
import pathlib

import numpy as np
from scipy.spatial import ConvexHull

from honeycomb import A1, A2, HEX, S3, solve

H = 2.6  # длина ячейки в единицах стороны шестиугольника (для модели)


def height_fn(faces):
    """f(p) по кускам дна (плоскости), с учётом периодичности решётки."""
    planes = []
    for P, h in faces:
        P = np.array(P); h = np.array(h)
        A = np.c_[P, np.ones(len(P))]
        coef, *_ = np.linalg.lstsq(A, h, rcond=None)
        planes.append((P, coef))

    def inside(P, q):
        n = len(P)
        s = [(P[(i + 1) % n][0] - P[i][0]) * (q[1] - P[i][1]) - (P[(i + 1) % n][1] - P[i][1]) * (q[0] - P[i][0]) for i in range(n)]
        return all(x >= -1e-7 for x in s) or all(x <= 1e-7 for x in s)

    def f(q):
        for i in range(-3, 4):
            for j in range(-3, 4):
                r = q - i * A1 - j * A2
                for P, c in planes:
                    if inside(P, r):
                        return c[0] * r[0] + c[1] * r[1] + c[2]
        raise ValueError(q)
    return f


def cell(d, faces):
    f = height_fn(faces)
    pts = []
    # дно: вершины кусков (и вершины шестиугольника)
    for P, h in faces:
        for p, z in zip(P, h):
            pts.append([p[0], p[1], z])
    # крышка: g(p) = f(p + d) + H; её изломы — в точках, где p + d — вершина кусков дна
    cand = [np.array(p) for P, _ in faces for p in P] + list(HEX)
    for P, _ in faces:
        for p in P:
            for i in range(-2, 3):
                for j in range(-2, 3):
                    q = np.array(p) - d + i * A1 + j * A2
                    if np.all(np.abs(q) < 1.01) and max(np.linalg.norm(q - v) for v in HEX) >= 0 and \
                            all((HEX[(k + 1) % 6][0] - HEX[k][0]) * (q[1] - HEX[k][1]) - (HEX[(k + 1) % 6][1] - HEX[k][1]) * (q[0] - HEX[k][0]) >= -1e-7 for k in range(6)):
                        cand.append(q)
    for q in cand:
        pts.append([q[0], q[1], f(q + d) + H])
    for v in HEX:  # вертикальные рёбра призмы
        pts.append([v[0], v[1], f(v)])
    P = np.unique(np.round(np.array(pts), 9), axis=0)
    hull = ConvexHull(P)
    # грани: группируем треугольники оболочки по плоскостям
    eq = np.round(hull.equations, 7)
    nf = len({tuple(e) for e in eq})
    V = P[np.unique(hull.simplices)]
    V = V - V.mean(axis=0)
    lattice = [[A1[0], A1[1], 0.0], [A2[0], A2[1], 0.0], [-d[0], -d[1], H]]
    return V, nf, lattice


if __name__ == "__main__":
    out = {}
    for key, d in (("bee", np.array([0.0, 1.0])), ("toth", A1 / 2)):
        F, faces = solve(d, key)
        V, nf, lat = cell(d, faces)
        vol_check = ConvexHull(V).volume
        print(key, "faces", nf, "vertices", len(V), f"volume {vol_check:.6f} vs hex·H {1.5 * S3 * H:.6f}")
        out[key] = {"points": np.round(V, 9).tolist(), "lattice": lat, "faces": nf, "F": F}
    js = "// Сгенерировано research/nature/cells3d.py: пчелиная ячейка и ячейка Тота (параллелоэдры).\n"
    js += "export const NATURE = " + json.dumps(out) + ";\n"
    (pathlib.Path(__file__).resolve().parents[2] / "docs" / "js" / "nature-data.js").write_text(js)
