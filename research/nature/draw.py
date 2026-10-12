"""Иллюстрации для страницы «Замощения в природе».

1. docs/img/comb-defect.svg — соты с дислокацией: треугольная решётка центров ячеек, смещённая полем краевой
   дислокации (Вольтерра), и её ячейки Вороного. Пятиугольник и семиугольник в ядре дислокации выделены —
   такие пары пчёлы и осы ставят на стыке ячеек разного размера.
2. docs/img/capsids.svg — икосаэдрические замощения сферы по Каспару–Клугу для T = 1, 3, 4, 7, 13:
   12 пятиугольников и 10(T − 1) шестиугольников (сферический Вороной центров белков).
"""
import itertools
import pathlib

import numpy as np
from scipy.spatial import SphericalVoronoi, Voronoi

IMG = pathlib.Path(__file__).resolve().parents[2] / "docs" / "img"
COL = {5: ("#d8605d", 0.8), 6: ("#e3a33b", 0.28), 7: ("#4f8fd9", 0.8)}


def comb_defect():
    W, H, S = 760, 330, 40.0
    pts = []
    for j in range(-14, 15):
        for i in range(-22, 23):
            pts.append([i + 0.5 * (j % 2), j * np.sqrt(3) / 2])
    P = np.array(pts)
    b, nu = 1.0, 0.3  # вектор Бюргерса вдоль x
    x, y = P[:, 0] + 0.31, P[:, 1] + 0.27  # ядро дислокации между узлами
    r2 = x * x + y * y
    ux = b / (2 * np.pi) * (np.arctan2(y, x) + x * y / (2 * (1 - nu) * r2))
    uy = -b / (2 * np.pi) * ((1 - 2 * nu) / (4 * (1 - nu)) * np.log(r2) + (x * x - y * y) / (4 * (1 - nu) * r2))
    Q = P + np.c_[ux, uy]
    vor = Voronoi(Q)
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}"><defs><clipPath id="c"><rect width="{W}" height="{H}" rx="14"/></clipPath></defs><g clip-path="url(#c)">']
    count = {}
    for k, reg in enumerate(vor.point_region):
        R = vor.regions[reg]
        if -1 in R or not R:
            continue
        V = vor.vertices[R]
        if np.abs(V[:, 0]).max() * S > W / 2 + 40 or np.abs(V[:, 1]).max() * S > H / 2 + 40:
            continue
        n = len(R)
        if np.abs(Q[k]).max() < 6:
            count[n] = count.get(n, 0) + 1
        col, op = COL.get(n, ("#e3a33b", 0.7))
        d = " ".join(f"{W / 2 + S * v[0]:.1f},{H / 2 - S * v[1]:.1f}" for v in V)
        out.append(f'<polygon points="{d}" fill="{col}" fill-opacity="{op}" stroke="#8a909b" stroke-width="1.4"/>')
    out.append("</g></svg>")
    (IMG / "comb-defect.svg").write_text("\n".join(out))
    print("comb-defect.svg: polygon counts near the core", dict(sorted(count.items())))


def icosahedron():
    phi = (1 + 5 ** 0.5) / 2
    V = np.array([(0, s1, s2 * phi) for s1 in (-1, 1) for s2 in (-1, 1)] + [(s1, s2 * phi, 0) for s1 in (-1, 1) for s2 in (-1, 1)]
                 + [(s2 * phi, 0, s1) for s1 in (-1, 1) for s2 in (-1, 1)], float)
    V /= np.linalg.norm(V, axis=1)[:, None]
    F = [f for f in itertools.combinations(range(12), 3)
         if all(abs(np.dot(V[a], V[b]) - 1 / 5 ** 0.5) < 1e-6 for a, b in itertools.combinations(f, 2))]
    return V, F


def capsid_points(h, k):
    """Центры капсомеров для T = h² + hk + k²: узлы треугольной решётки в каждой грани икосаэдра."""
    V, F = icosahedron()
    e1, e2 = np.array([1.0, 0.0]), np.array([0.5, np.sqrt(3) / 2])
    P0 = h * e1 + k * e2
    rot = np.array([[0.5, -np.sqrt(3) / 2], [np.sqrt(3) / 2, 0.5]])
    Q0 = rot @ P0
    M = np.c_[P0, Q0]
    pts = []
    n = h + k + 1
    for a in range(-n, 2 * n):
        for c in range(-n, 2 * n):
            p = a * e1 + c * e2
            bc = np.linalg.solve(M, p)  # p = s·P0 + t·Q0
            s, t = bc
            if s >= -1e-9 and t >= -1e-9 and s + t <= 1 + 1e-9:
                for f in F:
                    A, B, C = V[list(f)]
                    # ориентация грани наружу
                    if np.dot(np.cross(B - A, C - A), A) < 0:
                        B, C = C, B
                    q = (1 - s - t) * A + s * B + t * C
                    pts.append(q / np.linalg.norm(q))
    P = np.unique(np.round(np.array(pts), 7), axis=0)
    return P


def capsids():
    cases = [(1, 0), (1, 1), (2, 0), (2, 1), (3, 1)]
    R, gap = 70, 152
    W, H = gap * len(cases), 200
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" font-family="JetBrains Mono, monospace">']
    view = np.array([0.35, -0.25, 1.0]); view /= np.linalg.norm(view)
    ex = np.cross([0, 1, 0], view); ex /= np.linalg.norm(ex); ey = np.cross(view, ex)
    for ci, (h, k) in enumerate(cases):
        T = h * h + h * k + k * k
        P = capsid_points(h, k)
        assert len(P) == 10 * T + 2, (h, k, len(P))
        sv = SphericalVoronoi(P)
        sv.sort_vertices_of_regions()
        cx, cy = gap * ci + gap / 2, 90
        polys = []
        for i, reg in enumerate(sv.regions):
            V = sv.vertices[reg]
            c = V.mean(axis=0)
            depth = np.dot(c, view)
            if depth <= 0.05:
                continue
            d = " ".join(f"{cx + R * np.dot(v, ex):.1f},{cy - R * np.dot(v, ey):.1f}" for v in V)
            col, op = COL[len(reg)]
            shade = 0.55 + 0.45 * depth
            polys.append((depth, f'<polygon points="{d}" fill="{col}" fill-opacity="{op * shade:.2f}" stroke="#8a909b" stroke-width="1"/>'))
        out += [p for _, p in sorted(polys)]
        out.append(f'<text x="{cx}" y="{H - 16}" text-anchor="middle" font-size="13" fill="#8a909b">T = {T}</text>')
        print(f"T={T} ({h},{k}): {len(P)} capsomers, pentagons {sum(len(r) == 5 for r in sv.regions)}, hexagons {sum(len(r) == 6 for r in sv.regions)}")
    out.append("</svg>")
    (IMG / "capsids.svg").write_text("\n".join(out))


if __name__ == "__main__":
    comb_defect()
    capsids()
