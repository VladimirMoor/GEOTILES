"""Быстрые DV-ячейки с непрерывным метрическим параметром (c/a) и адаптивный поиск максимума граней.

Группа задаётся символом Германа–Могена (gemmi). Решётка: кубическая (a = 1), тетрагональная (1, 1, c)
или гексагональная/тригональная (a = 1, γ = 120°, c). Параметры поиска: точка x (фракционные координаты)
и c/a для некубических групп.

Проверка корректности: все вершины ячейки должны лежать ближе radius/2 к x. Тогда любой сосед, задающий
грань, находится ближе radius, и граней «от недостающих соседей» быть не может.
"""
import itertools
import sys
import time

import gemmi
import numpy as np
from scipy.spatial import HalfspaceIntersection


class Group:
    def __init__(self, hm):
        sg = gemmi.SpaceGroup(hm)
        self.number = sg.number
        self.hm = hm
        self.system = sg.crystal_system_str()
        ops = list(sg.operations())
        self.R = np.array([np.array(op.rot) / op.DEN for op in ops])  # (n,3,3)
        self.t = np.array([np.array(op.tran) / op.DEN for op in ops])  # (n,3)

    def cell(self, c):
        if self.system == "cubic":
            return np.eye(3)
        if self.system == "tetragonal":
            return np.diag([1.0, 1.0, c])
        if self.system in ("hexagonal", "trigonal"):
            return np.array([[1.0, -0.5, 0.0], [0.0, np.sqrt(3) / 2, 0.0], [0.0, 0.0, c]])
        raise ValueError(self.system)


class Evaluator:
    def __init__(self, G, c=1.0):
        self.G = G
        self.set_c(c)

    MAX_POINTS = 400_000  # жёсткий лимит точек орбиты на одно вычисление (≈ 10 МБ)

    def set_c(self, c):
        self.c = c
        self.B = self.G.cell(c)
        # Ячейка Дирихле точки x лежит внутри ячейки Вороного решётки (орбита содержит x + L),
        # поэтому её вершины не дальше радиуса покрытия решётки R_L ≤ половины длинной диагонали
        # элементарной ячейки. Соседи, задающие грани, лежат ближе 2·R_L. Радиус фиксирован и не растёт.
        corners = np.array(list(itertools.product((0, 1), repeat=3)), dtype=float) @ self.B.T
        diag = max(np.linalg.norm(u - v) for u in corners for v in corners)
        self.rcov = diag / 2
        self.radius = 2 * self.rcov * 1.02
        Binv = np.linalg.inv(self.B)
        rng = [int(np.ceil(self.radius * np.linalg.norm(Binv[i]))) + 1 for i in range(3)]
        n_shifts = np.prod([2 * r + 1 for r in rng])
        if n_shifts * len(self.G.R) > self.MAX_POINTS:
            raise MemoryError(f"too many orbit points: {n_shifts * len(self.G.R)}")
        grids = np.indices([2 * r + 1 for r in rng]).reshape(3, -1).T
        self.shifts = (grids - np.array(rng)).astype(float)

    def cell_of(self, x):
        X = self.B @ x
        base = np.einsum("nij,j->ni", self.G.R, x) + self.G.t  # (n,3)
        pts = (base[:, None, :] + self.shifts[None, :, :]).reshape(-1, 3) @ self.B.T
        d = np.linalg.norm(pts - X, axis=1)
        Y = pts[(d > 1e-9) & (d < self.radius)]
        Y = np.unique(np.round(Y, 11), axis=0)
        A = Y - X
        b = (np.sum(Y * Y, axis=1) - X @ X) / 2
        H = HalfspaceIntersection(np.hstack([A, -b[:, None]]), X)
        V = H.intersections
        if not np.all(np.isfinite(V)) or np.max(np.linalg.norm(V - X, axis=1)) > self.rcov * 1.0001:
            raise ArithmeticError("degenerate cell (numerical)")  # точку просто отбрасываем
        return X, A, b, V

    def facets(self, x, tol=1e-9):
        """(число граней, минимальный положительный зазор не-граней) — второе для навигации к стенам."""
        X, A, b, V = self.cell_of(x)
        nA = np.linalg.norm(A, axis=1)
        S = (V @ A.T - b) / nA  # знаковые расстояния вершин до плоскостей (≤ 0 внутри)
        nf = 0
        gaps = []
        for i in range(len(A)):
            on = np.abs(S[:, i]) < tol
            if on.sum() >= 3 and np.linalg.matrix_rank(V[on] - V[on][0], tol=1e-8) >= 2:
                nf += 1
            else:
                gaps.append(-S[:, i].max())  # насколько плоскость отстоит от ячейки
        g = min([v for v in gaps if v > tol], default=1.0)
        return nf, g, len(V)


def search(hm, n_random=20000, n_local=40, steps=600, c_range=(0.3, 4.0), seed=0, log=print):
    rng = np.random.default_rng(seed)
    G = Group(hm)
    cubic = G.system == "cubic"
    E = Evaluator(G, 1.0)

    def evaluate(p):
        x, c = p[:3] % 1.0, (1.0 if cubic else p[3])
        if not cubic:
            if not (c_range[0] <= c <= c_range[1]):
                return (-1, 1.0)
            E.set_c(c)
        try:
            nf, g, _ = E.facets(x)
        except Exception:
            return (-1, 1.0)
        return (nf, g)

    t0 = time.time()
    pts = []
    for _ in range(n_random):
        p = np.concatenate([rng.random(3), [np.exp(rng.uniform(np.log(c_range[0]), np.log(c_range[1])))]])
        pts.append((evaluate(p), p))
    pts.sort(key=lambda r: (-r[0][0], r[0][1]))
    hist = {}
    for (nf, _), _ in pts:
        hist[nf] = hist.get(nf, 0) + 1
    log(f"{hm}: random {n_random} in {time.time() - t0:.0f}s, top facet counts {sorted(hist.items())[-6:]}")
    best_overall = (0, None)
    for k in range(min(n_local, len(pts))):
        (nf, g), p = pts[k]
        cur, cp = (nf, g), p.copy()
        sigma = 1e-2
        for it in range(steps):
            q = cp + rng.normal(scale=sigma, size=4) * (np.array([1, 1, 1, 0]) if cubic else np.array([1, 1, 1, cp[3]]))
            val = evaluate(q)
            if val[0] > cur[0] or (val[0] == cur[0] and val[1] < cur[1]):
                cur, cp = val, q
            else:
                sigma *= 0.97 if sigma > 1e-10 else 1.0
                if sigma < 1e-10:
                    sigma = 1e-3
        if cur[0] > best_overall[0]:
            best_overall = (cur[0], cp.copy())
        log(f"  local {k}: start {nf} → {cur[0]} (gap {cur[1]:.2e}) at x={np.round(cp[:3] % 1, 6)}"
            + ("" if cubic else f", c/a={cp[3]:.6f}"))
    return best_overall


if __name__ == "__main__":
    hm = sys.argv[1]
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 20000
    nl = int(sys.argv[3]) if len(sys.argv) > 3 else 40
    st = int(sys.argv[4]) if len(sys.argv) > 4 else 600
    cr = tuple(float(v) for v in sys.argv[5].split(":")) if len(sys.argv) > 5 else (0.3, 4.0)
    best = search(hm, n_random=n, n_local=nl, steps=st, c_range=cr, log=lambda m: print(m, flush=True))
    print("BEST", hm, best[0], best[1])
