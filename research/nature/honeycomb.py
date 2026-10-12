"""Задача о дне пчелиной ячейки: двухслойные соты из шестиугольных призм, слои разделены общей кусочно-плоской
поверхностью f. Сравниваем дно пчёл (3 ромба) и дно Фейеша Тота (2 шестиугольника + 2 ромба).

Постановка (сторона шестиугольника = 1). Слой A — призмы над f, слой B — под f, сдвинут в плоскости на d.
При фиксированном объёме ячейки количество воска на пару ячеек (A + B) равно
    W = F(f) + 6 V / A_hex,   F(f) = S(f) + ½ (∮_{∂H_B} f ds − ∮_{∂H_A} f ds),
где S — площадь поверхности f на одну ячейку. Второе слагаемое от формы дна не зависит, поэтому сравниваем F.
Поверхность f плоская на каждом куске H_A ∩ H_B (это и есть грани дна); высоты в вершинах — неизвестные,
связанные условиями периодичности и плоскостности. Минимизируем F по ним.
"""
import itertools
import json
import pathlib

import numpy as np
from scipy.optimize import minimize

S3 = np.sqrt(3)
A1, A2 = np.array([S3, 0.0]), np.array([S3 / 2, 1.5])  # решётка центров (сторона 1, «острые» вершины по y)
HEX = np.array([[np.cos(np.pi / 6 + k * np.pi / 3), np.sin(np.pi / 6 + k * np.pi / 3)] for k in range(6)])
A_HEX = 1.5 * S3


def clip(poly, a, b):
    """Отсечение многоугольника полуплоскостью слева от направленного отрезка a→b (Сазерленд–Ходжман)."""
    out = []
    n = len(poly)
    side = lambda p: (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0])  # noqa: E731
    for i in range(n):
        p, q = poly[i], poly[(i + 1) % n]
        sp, sq = side(p), side(q)
        if sp >= -1e-12:
            out.append(p)
        if (sp >= -1e-12) != (sq >= -1e-12):
            t = sp / (sp - sq)
            out.append(p + t * (q - p))
    return out


def regions(d):
    """Куски H_A0 ∩ (H_B + сдвиги); возвращает список (многоугольник, индекс ячейки B)."""
    res = []
    for i, j in itertools.product(range(-2, 3), repeat=2):
        cB = d + i * A1 + j * A2
        poly = [p.copy() for p in HEX]
        hb = HEX + cB
        for k in range(6):
            poly = clip(poly, hb[k], hb[(k + 1) % 6])
            if not poly:
                break
        # убираем совпадающие соседние точки
        clean = []
        for q in poly:
            if not clean or np.linalg.norm(q - clean[-1]) > 1e-9:
                clean.append(q)
        if len(clean) > 1 and np.linalg.norm(clean[0] - clean[-1]) < 1e-9:
            clean.pop()
        poly = clean
        if len(poly) >= 3:
            P = np.array(poly)
            area = 0.5 * abs(np.dot(P[:, 0], np.roll(P[:, 1], -1)) - np.dot(P[:, 1], np.roll(P[:, 0], -1)))
            if area > 1e-9:
                res.append((P, (i, j)))
    return res


def solve(d, label):
    R = regions(d)
    # вершины по классам периодичности
    keys, pts = {}, []

    def key(p):
        # приводим точку по модулю решётки к фундаментальной области (координаты в базисе A1, A2)
        M = np.array([A1, A2]).T
        c = np.linalg.solve(M, p)
        c = np.round(c - np.floor(c + 1e-9), 6) % 1.0
        return (round(c[0], 5) % 1.0, round(c[1], 5) % 1.0)
    for P, _ in R:
        for p in P:
            k = key(p)
            if k not in keys:
                keys[k] = len(keys)
    nv = len(keys)
    # неизвестные: высоты в классах вершин; плоскостность каждого куска — линейные условия
    rows = []
    for P, _ in R:
        ids = [keys[key(p)] for p in P]
        # высота плоская: h(p) = α·p + β; условия для вершин сверх трёх
        # три неколлинеарные опорные точки
        bi = [0, 1]
        for m in range(2, len(P)):
            u, w = P[1] - P[0], P[m] - P[0]
            if abs(u[0] * w[1] - u[1] * w[0]) > 1e-6:
                bi.append(m)
                break
        base = P[bi]
        for m in range(len(P)):
            if m in bi:
                continue
            # p_m = Σ λ_k base_k с Σλ = 1 → h_m − Σ λ_k h_k = 0
            Mx = np.vstack([base.T, np.ones(3)])
            lam = np.linalg.solve(Mx, np.append(P[m], 1.0))
            r = np.zeros(nv)
            r[ids[m]] += 1
            for kk in range(3):
                r[ids[bi[kk]]] -= lam[kk]
            rows.append(r)
    C = np.array(rows) if rows else np.zeros((0, nv))
    # параметризация решений: h = N z
    if len(C):
        u, s, vt = np.linalg.svd(C)
        rank = int((s > 1e-9).sum())
        N = vt[rank:].T
    else:
        N = np.eye(nv)

    def F(z):
        h = N @ z
        S = 0.0
        for P, _ in R:
            ids = [keys[key(p)] for p in P]
            H = h[ids]
            # нормаль плоскости через первые 3 точки
            n = np.zeros(3)
            for m in range(1, len(P) - 1):  # нормаль по сумме треугольников веера (плоскость одна и та же)
                n += np.cross(np.append(P[m] - P[0], H[m] - H[0]), np.append(P[m + 1] - P[0], H[m + 1] - H[0]))
            area2d = 0.5 * abs(np.dot(P[:, 0], np.roll(P[:, 1], -1)) - np.dot(P[:, 1], np.roll(P[:, 0], -1)))
            S += area2d * np.linalg.norm(n) / abs(n[2])
        # интегралы по границам: f линейна на каждом отрезке между вершинами кусков
        def edge_int(c0):
            tot = 0.0
            hexc = HEX + c0
            for k in range(6):
                a, b = hexc[k], hexc[(k + 1) % 6]
                # точки разбиения: вершины кусков, лежащие на отрезке (по модулю решётки — через h в классах)
                ts = [0.0, 1.0]
                for P, _ in R:
                    for p in P:
                        for sh in itertools.product(range(-2, 3), repeat=2):
                            q = p + sh[0] * A1 + sh[1] * A2
                            ab = b - a
                            t = np.dot(q - a, ab) / np.dot(ab, ab)
                            if -1e-9 < t < 1 + 1e-9 and np.linalg.norm(a + t * ab - q) < 1e-7:
                                ts.append(min(max(t, 0), 1))
                ts = sorted(set(round(t, 9) for t in ts))
                for t0, t1 in zip(ts[:-1], ts[1:]):
                    p0, p1 = a + t0 * (b - a), a + t1 * (b - a)
                    h0, h1 = h[keys[key(p0)]] if key(p0) in keys else None, h[keys[key(p1)]] if key(p1) in keys else None
                    if h0 is None or h1 is None:
                        raise RuntimeError("edge point is not a region vertex")
                    tot += 0.5 * (h0 + h1) * np.linalg.norm(p1 - p0)
            return tot
        return S + 0.5 * (edge_int(d) - edge_int(np.zeros(2)))

    z0 = np.random.default_rng(0).normal(0, 0.1, N.shape[1])
    best = None
    for trial in range(6):
        r = minimize(F, z0 if trial == 0 else np.random.default_rng(trial).normal(0, 0.3, N.shape[1]), method="BFGS",
                     options={"gtol": 1e-11, "maxiter": 5000})
        if best is None or r.fun < best.fun:
            best = r
    h = N @ best.x
    faces = [(P.tolist(), [float(h[keys[key(p)]]) for p in P]) for P, _ in R]
    print(f"{label}: bottom faces {len(R)}, F = {best.fun:.9f}")
    return best.fun, faces


if __name__ == "__main__":
    flat = 1.5 * S3  # плоское дно: S = площадь шестиугольника, граничные слагаемые нулевые
    bee_F, bee = solve(np.array([0.0, 1.0]), "bee (B centres at alternate vertices)")  # центр B — в вершине A
    toth_F, toth = solve((A1) / 2, "Toth (B centres at edge midpoints)")
    exact_bee = 3 / np.sqrt(2)
    print(f"flat bottom F = {flat:.9f}; bee exact 3/√2 = {exact_bee:.9f}")
    # доля экономии для пропорций рабочей ячейки: ширина 5.2 мм между гранями, глубина ~ 11 мм
    s_mm = 5.2 / S3
    L = 11.0 / s_mm  # глубина в единицах стороны
    V = A_HEX * L
    W_bee = bee_F + 6 * V / A_HEX
    print(f"Toth saves {bee_F - toth_F:.6f} (units s²) = {(bee_F - toth_F) / bee_F:.3%} of the bottom term, "
          f"{(bee_F - toth_F) / W_bee:.4%} of all wax per cell pair at worker-cell proportions")
    out = {"flat": flat, "bee": bee_F, "toth": toth_F, "bee_exact": exact_bee, "bee_faces": bee, "toth_faces": toth,
           "saving_bottom": (bee_F - toth_F) / bee_F, "saving_total_worker": (bee_F - toth_F) / W_bee}
    pathlib.Path(__file__).with_name("honeycomb.json").write_text(json.dumps(out))
