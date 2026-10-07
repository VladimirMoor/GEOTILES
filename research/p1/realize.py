"""Шаг 2. Равносторонние выпуклые реализации каждого комбинаторного типа.

Неизвестные: координаты вершин x_v, единичные нормали граней n_f и сдвиги d_f.
Уравнения: |x_u − x_v|² = 1 на рёбрах, n_f·x_v = d_f для вершин грани, |n_f|² = 1,
плюс 6 калибровочных (первая грань: x_a = 0, x_b на оси x, x_c в плоскости xy). Система квадратная:
3E + F + 6 = 3V + 4F по формуле Эйлера, поэтому в общем случае решений конечное число,
а у «гибких» типов есть непрерывные семейства. Ищем многократным запуском Левенберга–Марквардта,
каждое решение строго проверяем на выпуклость и совпадение комбинаторики.

Запуск: python realize.py [--starts N] [--only F7-012]  →  data/realizations.json
"""
import argparse
import json
import math
import pathlib
import time
from multiprocessing import Pool

import numpy as np
from scipy.optimize import least_squares

ROOT = pathlib.Path(__file__).resolve().parent
TOL_EQ = 1e-9      # точность уравнений
TOL_STRICT = 1e-6  # строгость выпуклости


def edges_of(faces):
    es = set()
    for f in faces:
        for i, a in enumerate(f):
            b = f[(i + 1) % len(f)]
            es.add((min(a, b), max(a, b)))
    return sorted(es)


class System:
    def __init__(self, t):
        self.faces = t["faces"]
        self.V, self.F = t["nV"], t["nF"]
        self.E = edges_of(self.faces)
        a, b, c = self.faces[0][:3]
        self.gauge = (a, b, c)
        self.n = 3 * self.V + 4 * self.F
        self.eu = np.array([u for u, _ in self.E])
        self.ev = np.array([v for _, v in self.E])
        self.inc_f = np.array([f for f, face in enumerate(self.faces) for _ in face])
        self.inc_v = np.array([v for face in self.faces for v in face])
        self.gauge_idx = np.array([3 * a, 3 * a + 1, 3 * a + 2, 3 * b + 1, 3 * b + 2, 3 * c + 2])

    def split(self, z):
        X = z[: 3 * self.V].reshape(self.V, 3)
        N = z[3 * self.V: 3 * self.V + 3 * self.F].reshape(self.F, 3)
        d = z[3 * self.V + 3 * self.F:]
        return X, N, d

    def residuals(self, z, gauge=True):
        X, N, d = self.split(z)
        parts = [
            np.sum((X[self.eu] - X[self.ev]) ** 2, axis=1) - 1,
            np.sum(X[self.inc_v] * N[self.inc_f], axis=1) - d[self.inc_f],
            np.sum(N ** 2, axis=1) - 1,
        ]
        if gauge:
            parts.append(z[self.gauge_idx])
        return np.concatenate(parts)

    def jacobian(self, z, gauge=True):
        X, N, d = self.split(z)
        V3, F3 = 3 * self.V, 3 * self.F
        nE, nI, nF = len(self.E), len(self.inc_v), self.F
        J = np.zeros((nE + nI + nF + (6 if gauge else 0), self.n))
        r = np.arange(nE)
        g = 2 * (X[self.eu] - X[self.ev])
        for k in range(3):
            J[r, 3 * self.eu + k] = g[:, k]
            J[r, 3 * self.ev + k] = -g[:, k]
        r = nE + np.arange(nI)
        for k in range(3):
            J[r, 3 * self.inc_v + k] = N[self.inc_f, k]
            J[r, V3 + 3 * self.inc_f + k] = X[self.inc_v, k]
        J[r, V3 + F3 + self.inc_f] = -1
        r = nE + nI + np.arange(nF)
        for k in range(3):
            J[r, V3 + 3 * np.arange(nF) + k] = 2 * N[:, k]
        if gauge:
            J[nE + nI + nF + np.arange(6), self.gauge_idx] = 1
        return J

    def pack(self, X):
        """Нормали и сдвиги по положениям вершин (плоскость МНК, нормаль наружу)."""
        c = X.mean(axis=0)
        N, d = [], []
        for face in self.faces:
            P = X[face]
            m = P.mean(axis=0)
            nrm = np.linalg.svd(P - m)[2][-1]
            if nrm @ (m - c) < 0:
                nrm = -nrm
            N.append(nrm)
            d.append(nrm @ m)
        return np.concatenate([X.ravel(), np.ravel(N), d])


def gauge_fix(X, sys):
    a, b, c = sys.gauge
    X = X - X[a]
    e1 = X[b] / np.linalg.norm(X[b])
    w = X[c] - (X[c] @ e1) * e1
    e2 = w / np.linalg.norm(w)
    e3 = np.cross(e1, e2)
    return X @ np.array([e1, e2, e3]).T


def spectral(t):
    A = np.zeros((t["nV"], t["nV"]))
    for u, nb in enumerate(t["adj"]):
        for v in nb:
            A[u, v] = 1
    L = np.diag(A.sum(1)) - A
    w, U = np.linalg.eigh(L)
    return U[:, 1:4]


def verify(X, faces, E):
    """None, если не выпуклая реализация именно этого типа; иначе словарь с геометрией."""
    lens = np.array([np.linalg.norm(X[u] - X[v]) for u, v in E])
    if np.max(np.abs(lens - 1)) > 1e-7:
        return None
    c = X.mean(axis=0)
    normals = []
    for face in faces:
        P = X[face]
        m = P.mean(axis=0)
        s, vt = np.linalg.svd(P - m)[1:]
        nrm = vt[-1]
        if s[-1] > 1e-7 if len(face) > 3 else False:
            return None  # грань не плоская
        if nrm @ (m - c) < 0:
            nrm = -nrm
        others = [v for v in range(len(X)) if v not in face]
        if np.max((X[others] - m) @ nrm) > -TOL_STRICT:
            return None  # плоскость грани не опорная или задевает другие вершины
        k = len(face)
        for i in range(k):
            p0, p1, p2 = X[face[i - 1]], X[face[i]], X[face[(i + 1) % k]]
            if np.cross(p1 - p0, p2 - p1) @ nrm < TOL_STRICT:
                return None  # угол грани ≥ 180° или обход не совпадает
        normals.append(nrm)
    # двугранные углы по рёбрам
    edge_faces = {}
    for f, face in enumerate(faces):
        for i, a in enumerate(face):
            b = face[(i + 1) % len(face)]
            edge_faces.setdefault((min(a, b), max(a, b)), []).append(f)
    dihedrals = []
    for e in E:
        f1, f2 = edge_faces[e]
        cosv = np.clip(normals[f1] @ normals[f2], -1, 1)
        dihedrals.append(180 - math.degrees(math.acos(cosv)))
    face_angles = []
    for face in faces:
        k = len(face)
        angs = []
        for i in range(k):
            p0, p1, p2 = X[face[i - 1]], X[face[i]], X[face[(i + 1) % k]]
            a, b = p0 - p1, p2 - p1
            angs.append(math.degrees(math.acos(np.clip(a @ b / (np.linalg.norm(a) * np.linalg.norm(b)), -1, 1))))
        face_angles.append(angs)
    dist = sorted(np.round([np.linalg.norm(X[i] - X[j]) for i in range(len(X)) for j in range(i + 1, len(X))], 6))
    return {"dihedrals": dihedrals, "faceAngles": face_angles, "signature": tuple(dist)}


def flex_dim(sys, z):
    J = sys.jacobian(z, gauge=False)
    s = np.linalg.svd(J, compute_uv=False)
    nullity = sys.n - int(np.sum(s > 1e-7 * s[0]))
    return nullity - 6


def solve_type(args):
    t, starts, seed = args
    rng = np.random.default_rng(seed)
    sys = System(t)
    base = spectral(t)
    found = {}
    stats = {"converged": 0, "convex": 0}
    for k in range(starts):
        mode = k % 4
        if mode == 0:
            X0 = base + rng.normal(scale=0.05, size=base.shape)
        elif mode == 1:
            X0 = base + rng.normal(scale=0.25, size=base.shape)
        elif mode == 2:
            X0 = base @ rng.normal(size=(3, 3)) + rng.normal(scale=0.1, size=base.shape)
        else:
            X0 = rng.normal(size=base.shape)
        edge_len = np.mean([np.linalg.norm(X0[u] - X0[v]) for u, v in sys.E])
        X0 = gauge_fix(X0 / edge_len, sys)
        z0 = sys.pack(X0)
        try:
            res = least_squares(sys.residuals, z0, jac=sys.jacobian, method="lm", xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=400)
        except Exception:
            continue
        if np.max(np.abs(res.fun)) > TOL_EQ:
            continue
        stats["converged"] += 1
        X = sys.split(res.x)[0]
        g = verify(X, t["faces"], sys.E)
        if g is None:
            continue
        stats["convex"] += 1
        key = g["signature"]
        if key in found:
            found[key]["hits"] += 1
            continue
        found[key] = {"X": X.tolist(), "dihedrals": g["dihedrals"], "faceAngles": g["faceAngles"],
                      "flex": flex_dim(sys, res.x), "hits": 1, "firstStart": k}
    return t["id"], list(found.values()), stats


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--starts", type=int, default=400)
    ap.add_argument("--only", default=None)
    ap.add_argument("--out", default="realizations.json")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    types = json.loads((ROOT / "data" / "types.json").read_text())
    if args.only:
        types = [t for t in types if t["id"] in args.only.split(",")]
    t0 = time.time()
    with Pool() as pool:
        results = pool.map(solve_type, [(t, args.starts, args.seed * 100003 + i) for i, t in enumerate(types)], chunksize=1)
    out = {}
    for tid, sols, stats in results:
        out[tid] = {"solutions": sols, "stats": stats}
    (ROOT / "data" / args.out).write_text(json.dumps(out))
    rigid = sum(1 for v in out.values() for s in v["solutions"] if s["flex"] == 0)
    flexible_types = sorted(k for k, v in out.items() if any(s["flex"] > 0 for s in v["solutions"]))
    with_sol = sum(1 for v in out.values() if v["solutions"])
    print(f"{len(types)} types, {with_sol} with equilateral convex realizations; rigid solutions: {rigid}; "
          f"flexible types: {len(flexible_types)}; {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
