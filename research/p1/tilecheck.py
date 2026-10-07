"""Численная проверка разбиений: плитка = выпуклая оболочка, разбиение = набор изометрий.

Случайные точки шара, который заведомо покрыт построенным фрагментом, должны лежать строго внутри ровно
одной плитки. Плюс проверка объёма: V(плитки) · (число плиток на ячейку) = |det решётки|.
"""
import itertools
import json
import pathlib

import numpy as np
from scipy.spatial import ConvexHull

ROOT = pathlib.Path(__file__).resolve().parent


class Tile:
    def __init__(self, pts):
        self.pts = np.asarray(pts, float)
        h = ConvexHull(self.pts)
        self.volume = h.volume
        # плоскости гиперграней (qhull: n·x + c ≤ 0 внутри), склеиваем копланарные
        eqs = []
        for e in h.equations:
            if not any(np.allclose(e, q, atol=1e-9) for q in eqs):
                eqs.append(e)
        self.planes = np.array(eqs)
        self.centroid = self.pts[h.vertices].mean(axis=0)
        self.radius = np.max(np.linalg.norm(self.pts - self.centroid, axis=1))
        self.vertices = self.pts[h.vertices]

    def faces(self):
        out = []
        for e in self.planes:
            on = [i for i, p in enumerate(self.vertices) if abs(e[:3] @ p + e[3]) < 1e-9]
            out.append(on)
        return out

    def edge_lengths(self):
        """Рёбра = пары вершин, лежащие вместе ровно в двух гранях и образующие ребро обеих."""
        faces = self.faces()
        lens = []
        for f in faces:
            P = self.vertices[f]
            n = None
            # упорядочим вершины грани по углу
            c = P.mean(axis=0)
            nrm = np.linalg.svd(P - c)[2][-1]
            u = (P[0] - c) / np.linalg.norm(P[0] - c)
            w = np.cross(nrm, u)
            order = np.argsort(np.arctan2((P - c) @ w, (P - c) @ u))
            Q = P[order]
            lens.extend(np.linalg.norm(Q - np.roll(Q, -1, axis=0), axis=1))
        return np.array(lens)


def coverage(tile, maps, samples=3000, R=None, seed=0):
    """maps: список (M, t) — x ↦ M x + t. Точки в шаре радиуса R − 2·radius вокруг центра первой плитки."""
    rng = np.random.default_rng(seed)
    c0 = tile.centroid
    if R is None:
        R = max(np.linalg.norm(M @ c0 + t - c0) for M, t in maps)
    r = R - 2 * tile.radius
    if r <= 0:
        raise ValueError("фрагмент слишком мал")
    # обратные отображения: x ∈ M T + t  ⇔  Mᵀ(x − t) ∈ T (M ортогональна)
    inv = [(M.T, t) for M, t in maps]
    n, P = tile.planes[:, :3], tile.planes[:, 3]
    counts = np.zeros(samples, int)
    pts = rng.normal(size=(samples, 3))
    pts *= (r * rng.random(samples) ** (1 / 3) / np.linalg.norm(pts, axis=1))[:, None]
    pts += c0
    for Mi, t in inv:
        q = (pts - t) @ Mi.T
        inside = np.all(q @ n.T + P < -1e-10, axis=1)
        counts += inside
    return {"samples": samples, "exactlyOne": int(np.sum(counts == 1)), "gaps": int(np.sum(counts == 0)),
            "overlaps": int(np.sum(counts > 1)), "radius": r}


def lll(B, delta=0.75):
    """LLL-редукция базиса решётки (строки B)."""
    B = np.array(B, float)
    n = len(B)

    def gs(B):
        Q = np.zeros_like(B)
        mu = np.zeros((n, n))
        for i in range(n):
            Q[i] = B[i]
            for j in range(i):
                mu[i, j] = B[i] @ Q[j] / (Q[j] @ Q[j])
                Q[i] = Q[i] - mu[i, j] * Q[j]
        return Q, mu

    k = 1
    while k < n:
        for j in range(k - 1, -1, -1):
            Q, mu = gs(B)
            q = round(mu[k, j])
            if q:
                B[k] = B[k] - q * B[j]
        Q, mu = gs(B)
        if Q[k] @ Q[k] >= (delta - mu[k, k - 1] ** 2) * (Q[k - 1] @ Q[k - 1]):
            k += 1
        else:
            B[[k, k - 1]] = B[[k - 1, k]]
            k = max(k - 1, 1)
    return B


def crystal(lattice, motifs, N, R=None, center=None):
    """x ↦ M x + t + i a + j b + k c. Если задан R — диапазоны i, j, k подбираются так, чтобы
    сдвиги покрывали шар радиуса R (с запасом N шагов); иначе |i|,|j|,|k| ≤ N."""
    L = lll(lattice)
    if R is None:
        ranges = [range(-N, N + 1)] * 3
    else:
        Linv = np.linalg.inv(L)  # координаты: p = c @ L  ⇒  c = p @ Linv
        ext = [int(np.ceil(R * np.linalg.norm(Linv[:, i]))) + N for i in range(3)]
        ranges = [range(-e, e + 1) for e in ext]
    out = []
    for i, j, k in itertools.product(*ranges):
        s = i * L[0] + j * L[1] + k * L[2]
        if R is not None and np.linalg.norm(s) > R + 3 * max(np.linalg.norm(L, axis=1)):
            continue
        for M, t in motifs:
            out.append((np.asarray(M, float), np.asarray(t, float) + s))
    return out


# ---------- семейство F8-203: обобщённый гиробифастигиум ----------

def gyro_member(phi, psi1, psi2, side1=1, side2=1):
    """Ромб (a, b) с углом phi; призма 1 — стержни вдоль a над ромбом, треугольник на ребре b под углом psi1;
    призма 2 — стержни вдоль b под ромбом, треугольник на ребре a под углом psi2."""
    a = np.array([1.0, 0, 0])
    b = np.array([np.cos(phi), np.sin(phi), 0])
    z = np.array([0, 0, 1.0])
    mb = side1 * np.cross(z, b)          # в плоскости ромба, ⊥ b
    ma = side2 * np.cross(z, a)          # в плоскости ромба, ⊥ a
    u1 = np.cos(psi1) * mb + np.sin(psi1) * z
    u2 = np.cos(psi2) * ma - np.sin(psi2) * z
    apex1 = b / 2 + np.sqrt(3) / 2 * u1
    apex2 = a / 2 + np.sqrt(3) / 2 * u2
    T1 = [np.zeros(3), b, apex1]
    T2 = [np.zeros(3), a, apex2]
    pts = T1 + [p + a for p in T1] + T2 + [p + b for p in T2]
    return np.array(pts), a, b, apex1, apex2


def gyro_tiling(a, b, apex1, apex2, N=3):
    I, minusI = np.eye(3), -np.eye(3)
    motifs = [(I, np.zeros(3)), (minusI, a + b + apex1)]
    return crystal([a, b, apex1 - apex2], motifs, N)


if __name__ == "__main__":
    rng = np.random.default_rng(7)
    report = []
    for trial in range(12):
        phi = rng.uniform(40, 140) * np.pi / 180
        psi1 = rng.uniform(40, 140) * np.pi / 180
        psi2 = rng.uniform(40, 140) * np.pi / 180
        pts, a, b, ap1, ap2 = gyro_member(phi, psi1, psi2)
        try:
            T = Tile(pts)
        except Exception as e:
            print("hull failed", e)
            continue
        nf = len(T.planes)
        lens = T.edge_lengths()
        maps = gyro_tiling(a, b, ap1, ap2, N=4)
        cov = coverage(T, maps, samples=4000, R=3.2)
        cell = abs(np.linalg.det(np.array([a, b, ap1 - ap2]))) / 2
        print(f"phi={np.degrees(phi):6.1f} psi1={np.degrees(psi1):6.1f} psi2={np.degrees(psi2):6.1f}  faces={nf} "
              f"edges∈[{lens.min():.6f},{lens.max():.6f}]  V/cell={T.volume / cell:.6f}  coverage={cov}")


# ---------- семейство F8-255: обобщённый удлинённый гиробифастигиум ----------

def unit(v):
    return v / np.linalg.norm(v)


def elongated_member(a, b, c):
    """Параллелепипед (a, b, c); сверху призма — стержни вдоль a, треугольник на ребре b в плоскости (b, c);
    снизу призма — стержни вдоль b, треугольник на ребре a в плоскости (a, c)."""
    a, b, c = map(np.asarray, (a, b, c))
    u1 = unit(c - (c @ b) * b)
    u2 = unit(c - (c @ a) * a)
    apex1 = c + b / 2 + np.sqrt(3) / 2 * u1
    apex2 = a / 2 - np.sqrt(3) / 2 * u2
    O = np.zeros(3)
    box = [O, a, b, a + b, c, a + c, b + c, a + b + c]
    pts = box + [apex1, apex1 + a, apex2, apex2 + b]
    lattice = [a, b, c + apex1 - apex2]
    motifs = [(np.eye(3), O), (-np.eye(3), a + b + c + apex1)]
    return np.array(pts), lattice, motifs


def random_unit(rng):
    v = rng.normal(size=3)
    return v / np.linalg.norm(v)

