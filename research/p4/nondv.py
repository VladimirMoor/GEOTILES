"""Поиск стереоэдров вне класса Дирихле–Вороного: тайл как многогранник со свободными плоскостями граней.

Группа Γ (dvfast.Group, декартовы операции y → R y + t). Выбираем набор S элементов γ (кандидаты в соседи) и
для каждой пары {γ, γ⁻¹} одну плоскость H_γ = {n·y = d}; тогда плоскость пары γ⁻¹ — образ γ⁻¹(H_γ), а тайл
    P = ∩_γ {n_γ·y ≤ d_γ} ∩ {(−Rᵀ n_γ)·y ≤ −(d_γ − n_γ·t_γ)}.
Ячейка Дирихле точки x — частный случай (H_γ — серединный перпендикуляр к x, γx).

P — фундаментальная область Γ (значит, стереоэдр с тривиальным стабилизатором) тогда и только тогда, когда
vol(P) = V₀ = covol/|Γ:L| и копии γP не перекрываются. Минимизируем
    f = |vol P − V₀|/V₀ + Σ_{γ∈S} vol(P ∩ γP)/V₀
по параметрам плоскостей из случайных стартов. Для найденного P проверяем, ячейка ли это Дирихле:
ищем x, для которого все плоскости граней — серединные перпендикуляры (линейные наименьшие квадраты).
"""
import itertools
import sys

import numpy as np
from scipy.optimize import linprog, minimize
from scipy.spatial import ConvexHull, HalfspaceIntersection

from dvfast import Group


class Search:
    def __init__(self, hm, c=1.0, r_factor=1.35, seed=0):
        self.G = Group(hm)
        B = self.G.cell(c)
        self.B = B
        Binv = np.linalg.inv(B)
        self.V0 = abs(np.linalg.det(B)) / len(self.G.R)
        rng = np.random.default_rng(seed)
        self.x0 = B @ rng.random(3)
        # декартовы операции
        ops = [(B @ R @ Binv, B @ t) for R, t in zip(self.G.R, self.G.t)]
        # кандидаты: γ = (R, t + Bn), ближайшие к x0 образы
        cand = []
        for (R, t), n in itertools.product(ops, itertools.product(range(-2, 3), repeat=3)):
            tt = t + B @ np.array(n, float)
            if np.allclose(R, np.eye(3)) and np.allclose(tt, 0):
                continue
            cand.append((np.linalg.norm(R @ self.x0 + tt - self.x0), R, tt))
        cand.sort(key=lambda z: z[0])
        # радиус: несколько ближайших «оболочек» соседей
        rmin = cand[0][0]
        self.S = [(R, t) for dist, R, t in cand if dist < r_factor * 2.2 * rmin ** 0 * max(rmin, 1e-9) * 1.0 + 0 * dist and dist < 2.4 * rmin] \
            if False else [(R, t) for dist, R, t in cand if dist < r_factor * self._shell(cand)]
        # пары {γ, γ⁻¹}
        self.pairs = []
        used = set()
        for i, (R, t) in enumerate(self.S):
            if i in used:
                continue
            Ri, ti = R.T, -R.T @ t  # γ⁻¹ (R ортогональна)
            j = next((k for k, (R2, t2) in enumerate(self.S) if k not in used and np.allclose(R2, Ri, atol=1e-9) and np.allclose(t2, ti, atol=1e-9)), None)
            used.add(i)
            if j is not None:
                used.add(j)
            self.pairs.append((R, t))
        self.box = 4 * np.linalg.norm(B, axis=0).max()
        self.restrict()

    def restrict(self, k_pairs=1.4, k_over=2.6):
        p = self.dv_params(self.x0)
        H = self.halfspaces(p)
        c, _ = self.cheb(H)
        _, V = self.vol(H, c)
        dist = lambda R, t: np.linalg.norm(R @ self.x0 + t - self.x0)  # noqa: E731
        act = [k for k in range(len(self.pairs)) if any((np.abs(V @ H[2 * k + j][:3] + H[2 * k + j][3]) < 1e-7).sum() >= 3 for j in (0, 1))]
        rmax = max(dist(*self.pairs[k]) for k in act)
        rc = np.max(np.linalg.norm(V - self.x0, axis=1))
        self.pairs = [(R, t) for R, t in self.pairs if dist(R, t) <= k_pairs * rmax]
        self.S = [(R, t) for R, t in self.S if dist(R, t) <= k_over * rc]

    def _shell(self, cand):
        d = np.array([z[0] for z in cand])
        # радиус, внутри которого лежат соседи ячейки Дирихле x0 (с запасом) — грубая оценка: 2.2 × медиана 14 ближайших
        return 2.0 * np.median(d[:14])

    def dv_params(self, x):
        """Параметры плоскостей ячейки Дирихле точки x."""
        p = []
        for R, t in self.pairs:
            gx = R @ x + t
            n = gx - x
            p += list(n) + [(gx @ gx - x @ x) / 2]
        return np.array(p)

    def halfspaces(self, p):
        H = []
        for k, (R, t) in enumerate(self.pairs):
            n, d = p[4 * k:4 * k + 3], p[4 * k + 3]
            s = np.linalg.norm(n) + 1e-12
            H.append(np.append(n / s, -d / s))
            m = -R.T @ n
            H.append(np.append(m / s, (d - n @ t) / s))
        b = self.box
        for i in range(3):  # большая коробка, чтобы P был ограничен
            e = np.zeros(3); e[i] = 1
            H.append(np.append(e, -b - self.x0[i])); H.append(np.append(-e, -b + self.x0[i]))
        return np.array(H)

    @staticmethod
    def cheb(H):
        A, b = H[:, :3], -H[:, 3]
        nrm = np.linalg.norm(A, axis=1)
        r = linprog(np.r_[0, 0, 0, -1], A_ub=np.c_[A, nrm], b_ub=b, bounds=[(None, None)] * 3 + [(0, None)], method="highs")
        if r.status != 0:
            return None, 0.0
        return r.x[:3], r.x[3]

    @staticmethod
    def vol(H, c):
        try:
            hs = HalfspaceIntersection(H, c)
            return ConvexHull(hs.intersections).volume, hs.intersections
        except Exception:
            return 0.0, None

    def transform(self, H, R, t):
        """Полупространства для γP, если H — для P: γP = {y: γ⁻¹y ∈ P}."""
        A, b = H[:, :3], H[:, 3]
        A2 = A @ R.T
        return np.c_[A2, b - A2 @ t]

    def pool(self):
        """Все копии γ = операция × сдвиг решётки в пределах ±3 ячеек (для отбора перекрытий по текущему положению P)."""
        if not hasattr(self, "_pool"):
            B, Binv = self.B, np.linalg.inv(self.B)
            Rs, ts = [], []
            for R, t in zip(self.G.R, self.G.t):
                Rc, tc = B @ R @ Binv, B @ t
                for n in itertools.product(range(-3, 4), repeat=3):
                    tt = tc + B @ np.array(n, float)
                    if np.allclose(Rc, np.eye(3)) and np.allclose(tt, 0):
                        continue
                    Rs.append(Rc); ts.append(tt)
            self._pool = (np.array(Rs), np.array(ts))
        return self._pool

    def near_copies(self, c, radius):
        Rs, ts = self.pool()
        d = np.linalg.norm(np.einsum("nij,j->ni", Rs, c) + ts - c, axis=1)
        idx = np.where(d < radius)[0]
        return [(Rs[i], ts[i]) for i in idx]

    def f(self, p):
        H = self.halfspaces(p)
        c, r = self.cheb(H)
        if c is None or r < 1e-6:
            return 10.0
        v, Vt = self.vol(H, c)
        if Vt is None:
            return 10.0
        tot = abs(v - self.V0) / self.V0
        rad = 2.05 * np.max(np.linalg.norm(Vt - c, axis=1))
        for R, t in self.near_copies(c, rad):
            H2 = np.vstack([H, self.transform(H, R, t)])
            c2, r2 = self.cheb(H2)
            if c2 is not None and r2 > 1e-7:
                v2, _ = self.vol(H2, c2)
                tot += v2 / self.V0
        return tot

    def facets(self, p, tol=1e-7):
        H = self.halfspaces(p)
        c, _ = self.cheb(H)
        v, V = self.vol(H, c)
        n = 0
        planes = []
        for h in H[:-6]:
            on = np.abs(V @ h[:3] + h[3]) < tol
            if on.sum() >= 3 and np.linalg.matrix_rank(V[on] - V[on][0], tol=1e-6) == 2:
                n += 1
                planes.append(h)
        return n, v, np.array(planes), V

    def dv_residual(self, p):
        """Насколько P далёк от ячейки Дирихле: min_x Σ (плоскость грани − серединный перпендикуляр к x, γx)²."""
        H = self.halfspaces(p)
        c, _ = self.cheb(H)
        _, V = self.vol(H, c)
        rows, rhs = [], []
        for k, (R, t) in enumerate(self.pairs):
            for (Rg, tg), h in (((R, t), H[2 * k]), ((R.T, -R.T @ t), H[2 * k + 1])):
                on = np.abs(V @ h[:3] + h[3]) < 1e-7
                if on.sum() < 3 or np.linalg.matrix_rank(V[on] - V[on][0], tol=1e-6) < 2:
                    continue
                # бисектор x и γx: (γx − x)·y = (|γx|² − |x|²)/2 для y на грани; линейно по x не сразу —
                # используем: точки грани y равноудалены от x и γx ⇔ |y − x|² = |y − Rg x − tg|²
                # ⇔ 2 y·(Rg x + tg − x) = |Rg x + tg|² − |x|² = 2 x·Rgᵀtg + |tg|²  (Rg ортогональна)
                for y in V[on]:
                    rows.append(2 * (Rg.T @ y - y) - 2 * Rg.T @ tg)
                    rhs.append(tg @ tg - 2 * y @ tg)
        A, b = np.array(rows), np.array(rhs)
        x, *_ = np.linalg.lstsq(A, b, rcond=None)
        return float(np.sqrt(np.mean((A @ x - b) ** 2))), x


if __name__ == "__main__":
    hm = sys.argv[1] if len(sys.argv) > 1 else "P 21 3"
    S = Search(hm, seed=1)
    print(hm, "|Γ:L| =", len(S.G.R), "candidates", len(S.S), "pairs", len(S.pairs), "V0", S.V0)
    p = S.dv_params(S.x0)
    print("DV start: f =", S.f(p), "facets", S.facets(p)[0], "dv residual", S.dv_residual(p)[0])
