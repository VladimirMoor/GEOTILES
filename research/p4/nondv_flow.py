"""«Поток покрытия»: градиентный спуск к фундаментальной области по плоскостям граней.

J = Σ_{γ≠e} vol(P ∩ γP) + (V₀ − vol P) — это дисперсия числа покрывающих копий на период (с точностью до |Γ:L|),
J ≥ 0 и J = 0 ⇔ P замощает пространство копиями γP.

Сдвиг грани F наружу на s(y) меняет J на ∫_F s(y) (2 N(y) − 1) dA, где N(y) — число *других* копий,
покрывающих точку сразу за гранью. Плоскость пары k: H_a = {n·y = d} (грань к соседу γ) и H_b = γ⁻¹(H_a)
(грань к соседу γ⁻¹). Для y на H_a сдвиг s = δd − δn·y; для y на H_b сдвиг s = −δd + δn·(R y + t).
Отсюда градиент по (n, d); n держим единичным (градиент по n проецируем на касательную).
"""
import json
import sys
import time

import numpy as np

from nondv import Search
from nondv_verify import verify


def facet_samples(V, h, m=36):
    """Точки на грани с плоскостью h (равномерно по площади) и их веса."""
    on = np.abs(V @ h[:3] + h[3]) < 1e-8
    P = V[on]
    if len(P) < 3 or np.linalg.matrix_rank(P - P[0], tol=1e-9) < 2:
        return None, None
    c = P.mean(axis=0)
    n = h[:3]
    u = P[0] - c
    if np.linalg.norm(u) < 1e-12:
        u = P[1] - c
    u /= np.linalg.norm(u)
    w = np.cross(n, u)
    P = P[np.argsort(np.arctan2((P - c) @ w, (P - c) @ u))]
    pts, wts = [], []
    tris = [(c, P[i], P[(i + 1) % len(P)]) for i in range(len(P))]
    areas = np.array([0.5 * np.linalg.norm(np.cross(b - a, cc - a)) for a, b, cc in tris])
    A = areas.sum()
    if A < 1e-14:
        return None, None
    k = max(3, int(np.ceil(m / len(tris))))
    # барицентрическая сетка в каждом треугольнике
    bary = [(i, j) for i in range(k) for j in range(k - i)]
    for (a, b, cc), ar in zip(tris, areas):
        for i, j in bary:
            l1, l2 = (i + 1 / 3) / k, (j + 1 / 3) / k
            pts.append(a + l1 * (b - a) + l2 * (cc - a))
            wts.append(ar / len(bary))
    return np.array(pts), np.array(wts)


def coverage(S, H, c, rad, Y):
    """Число копий γP (γ ≠ e), содержащих каждую точку Y."""
    cnt = np.zeros(len(Y), int)
    for R, t in S.near_copies(c, rad):
        Hg = S.transform(H, R, t)
        cnt += np.all(Y @ Hg[:, :3].T + Hg[:, 3] <= 1e-12, axis=1)
    return cnt


def grad(S, p, eps=1e-7):
    H = S.halfspaces(p)
    c, r = S.cheb(H)
    if c is None or r < 1e-8:
        return None, None
    _, V = S.vol(H, c)
    rad = 2.05 * np.max(np.linalg.norm(V - c, axis=1))
    g = np.zeros_like(p)
    total_bad = 0.0
    for k, (R, t) in enumerate(S.pairs):
        for side in (0, 1):
            h = H[2 * k + side]
            Y, W = facet_samples(V, h)
            if Y is None:
                continue
            N = coverage(S, H, c, rad, Y + eps * h[:3])
            G = 2 * N - 1  # dJ на единицу сдвига наружу
            total_bad += np.sum(W * (N != 1))
            # параметры пары: n (единичный, хранится ненормированным), d; сдвиг наружу s
            nn = p[4 * k:4 * k + 3]
            s_ = np.linalg.norm(nn)
            if side == 0:
                g[4 * k + 3] += np.sum(W * G) / s_
                g[4 * k:4 * k + 3] += -(W * G) @ Y / s_
            else:
                g[4 * k + 3] += -np.sum(W * G) / s_
                g[4 * k:4 * k + 3] += (W * G) @ (Y @ R.T + t) / s_
    return g, total_bad


def J(S, p):
    return S.f_signed(p)


def flow(S, p, steps=2500, step0=0.02, decay=0.998, log=False):
    L = S.V0 ** (1 / 3)
    best = (np.inf, p)
    for it in range(steps):
        g, bad = grad(S, p)
        if g is None:
            break
        # проекция градиента по n на касательную к сфере |n| = const
        for k in range(len(S.pairs)):
            n = p[4 * k:4 * k + 3]
            u = n / np.linalg.norm(n)
            g[4 * k:4 * k + 3] -= (g[4 * k:4 * k + 3] @ u) * u
        gn = np.linalg.norm(g)
        if gn < 1e-14:
            break
        step = step0 * L * (decay ** it)
        p = p - step * g / gn * np.sqrt(len(S.pairs))
        for k in range(len(S.pairs)):  # нормировка: |n| = 1, d масштабируется вместе с n
            s_ = np.linalg.norm(p[4 * k:4 * k + 3])
            p[4 * k:4 * k + 4] /= s_
        if it % 100 == 0 or it == steps - 1:
            j = S.f(p)
            if j < best[0]:
                best = (j, p.copy())
            if log:
                print(f"  it {it:4d}  f {j:.3e}  bad-area {bad:.3e}  step {step:.2e}", flush=True)
            if j < 1e-9:
                break
    return best


# ---------- точная доводка: совпадение парных граней по вершинам ----------
from scipy.optimize import least_squares  # noqa: E402


def facet_vertices(S, p):
    H = S.halfspaces(p)
    c, r = S.cheb(H)
    if c is None or r < 1e-9:
        return None, None
    _, V = S.vol(H, c)
    F = {}
    for i, h in enumerate(H[:-6]):
        on = np.abs(V @ h[:3] + h[3]) < 1e-7
        if on.sum() >= 3 and np.linalg.matrix_rank(V[on] - V[on][0], tol=1e-7) == 2:
            F[i] = V[on]
    return F, V


def polish(S, p, iters=8):
    """Least squares: для пары k вершины γ_k(F_b) совпадают с вершинами F_a; непарные грани стягиваются в 0."""
    for _ in range(iters):
        F, _ = facet_vertices(S, p)
        if F is None:
            return p
        # сопоставление вершин фиксируем по текущему положению
        plan = []
        areas = {i: ConvexHull_area(P) for i, P in F.items()}
        thr = 0.03 * np.median(list(areas.values()))
        for k, (R, t) in enumerate(S.pairs):
            a, b = 2 * k, 2 * k + 1
            small = (a in F and areas[a] < thr) or (b in F and areas[b] < thr)
            if small:
                for i in (a, b):
                    if i in F:
                        plan.append(("kill", i))
                continue
            if a in F and b in F:
                img = F[b] @ R.T + t
                match = [int(np.argmin(np.linalg.norm(F[a] - q, axis=1))) for q in img]
                plan.append(("pair", k, len(F[b]), match))
            elif a in F or b in F:
                plan.append(("kill", a if a in F else b))

        MAXV = 14

        def res(q):
            F2, _ = facet_vertices(S, q)
            n_out = sum(2 * MAXV if x[0] == "pair" else 1 for x in plan)
            if F2 is None:
                return np.full(n_out, 1.0)
            out = []
            for item in plan:
                if item[0] == "pair":
                    k = item[1]
                    R, t = S.pairs[k]
                    a, b = 2 * k, 2 * k + 1
                    if a not in F2 or b not in F2:
                        out += [1.0] * (2 * MAXV)
                        continue
                    A_, B_ = F2[a], F2[b] @ R.T + t
                    d1 = [np.min(np.linalg.norm(A_ - q_, axis=1)) for q_ in B_][:MAXV]
                    d2 = [np.min(np.linalg.norm(B_ - q_, axis=1)) for q_ in A_][:MAXV]
                    out += d1 + [0.0] * (MAXV - len(d1)) + d2 + [0.0] * (MAXV - len(d2))
                else:
                    i_ = item[1]
                    out.append(np.sqrt(ConvexHull_area(F2[i_])) if i_ in F2 else 0.0)
            return np.array(out)
        r = least_squares(res, p, method="trf", max_nfev=600, xtol=1e-15, ftol=1e-15, gtol=1e-15)
        p = r.x
        if S.f(p) < 1e-11:
            break
    return p


def ConvexHull_area(P):
    c = P.mean(axis=0)
    u, s, vt = np.linalg.svd(P - c)
    Q = (P - c) @ vt[:2].T
    if len(Q) < 3:
        return 0.0
    from scipy.spatial import ConvexHull as CH
    try:
        return CH(Q).volume
    except Exception:
        return 0.0


def random_start(S, rng, spread=1.0):
    """Случайные плоскости, при которых P содержит точку x: для пары (γ) плоскость лежит между x и γx,
    d = n·x + u·n·(γx − x), u ∈ (0.15, 0.85); нормаль — направление на γx с поворотом порядка spread."""
    x = S.B @ rng.random(3)
    p = []
    for R, t in S.pairs:
        v = R @ x + t - x
        while True:
            n = v / np.linalg.norm(v) + spread * 0.5 * rng.normal(size=3)
            n /= np.linalg.norm(n)
            if n @ v > 0.15 * np.linalg.norm(v):
                break
        u = rng.uniform(0.15, 0.85)
        p += list(n) + [n @ x + u * (n @ v)]
    return np.array(p)


def dv_distance(S, p):
    """Расстояние от P до ближайшей ячейки Дирихле: x подбирается по настоящим соседям крупных граней,
    затем vol(P Δ DV(x)) / V₀."""
    from nondv_verify import copies_near, inside
    H = S.halfspaces(p)
    c, _ = S.cheb(H)
    _, V = S.vol(H, c)
    F, _ = facet_vertices(S, p)
    areas = {i: ConvexHull_area(P) for i, P in F.items()}
    thr = 0.03 * np.median(list(areas.values()))
    Rp = np.max(np.linalg.norm(V - c, axis=1))
    G = copies_near(S, c, 3 * Rp)
    Hs = [S.transform(H, R, t) for R, t in G]
    rows, rhs = [], []
    for i, P in F.items():
        if areas[i] < thr:
            continue
        y = P.mean(axis=0) + 1e-6 * H[i][:3]
        ks = [k for k, h in enumerate(Hs) if inside(h, y, 1e-12)]
        if len(ks) != 1:
            continue
        R, t = G[ks[0]]
        for y in P:
            rows.append(2 * (R.T @ y - y) - 2 * R.T @ t)
            rhs.append(t @ t - 2 * y @ t)
    x, *_ = np.linalg.lstsq(np.array(rows), np.array(rhs), rcond=None)
    q = S.dv_params(x)
    Hq = S.halfspaces(q)
    both = np.vstack([H, Hq])
    cb, rb = S.cheb(both)
    vi = S.vol(both, cb)[0] if cb is not None and rb > 1e-9 else 0.0
    vp = S.vol(H, c)[0]
    cq, _ = S.cheb(Hq)
    vq = S.vol(Hq, cq)[0]
    return float((vp + vq - 2 * vi) / S.V0), x


if __name__ == "__main__":
    hm, mode, sig, runs, seed0 = sys.argv[1], sys.argv[2], float(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
    for r in range(runs):
        seed = seed0 + 1000 * r
        S = Search(hm, seed=seed)
        rng = np.random.default_rng(seed + 7)
        if mode == "dv":
            x = S.B @ rng.random(3)
            p0 = S.dv_params(x)
            for k in range(len(S.pairs)):
                p0[4 * k:4 * k + 4] /= np.linalg.norm(p0[4 * k:4 * k + 3])
            p0 = p0 + rng.normal(0, sig * S.V0 ** (1 / 3), p0.size)
        else:
            p0 = random_start(S, rng, sig)
        t0 = time.time()
        fbest, pbest = flow(S, p0, log=("-v" in sys.argv))
        out = {"group": hm, "mode": mode, "sigma": sig, "seed": seed, "f_flow": float(fbest), "pairs": len(S.pairs)}
        if fbest < 1e-3:
            out["dv_distance"] = dv_distance(S, pbest)[0]
        if fbest < 1e-3 and out["dv_distance"] > 0.10:
            pp = polish(S, pbest)
            fp = S.f(pp)
            out["f_polished"] = float(fp)
            np.save(f"data/nondv/cand_{hm.replace(' ', '')}_{mode}{sig}_{seed}.npy", pp)
            if fp < 1e-9:
                out["verify"] = verify(S, pp, 1500)
        out["sec"] = round(time.time() - t0, 1)
        print(json.dumps(out), flush=True)
