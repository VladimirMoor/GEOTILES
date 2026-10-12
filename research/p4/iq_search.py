"""Самый «круглый» стереоэдр: минимум A/V^(2/3) по ячейкам Дирихле–Вороного во всех группах 75–230.

Для параллелоэдров (замощения сдвигами) минимум даёт правильный усечённый октаэдр, 3(1+2√3)/4^(2/3) ≈ 5.3147
(Bezdek; доказательства Hales–Song и Cesaroni–Novaga, сентябрь 2026). Здесь разрешены повороты: ячейки
Дирихле точки x под пространственной группой Γ (с любым стабилизатором), метрика — c/a для некубических групп.

Для каждой группы: случайные (x, c/a), затем локальная оптимизация (Нелдер–Мид) лучших. Пишет JSON-строки.
Запуск: nice -n 15 python iq_search.py IT_lo IT_hi n_random n_local > data/iq/....log
"""
import json
import sys
import time

import gemmi
import numpy as np
from scipy.optimize import minimize
from scipy.spatial import ConvexHull

from dvfast import Group, Evaluator

TO = 3 * (1 + 2 * np.sqrt(3)) / 4 ** (2 / 3)


def score(G, x, c, cache):
    key = round(float(c), 6)
    if key not in cache:
        cache.clear()
        cache[key] = Evaluator(G, c)
    E = cache[key]
    try:
        _, _, _, V = E.cell_of(np.asarray(x) % 1.0)
        h = ConvexHull(V)
        return h.area / h.volume ** (2 / 3)
    except Exception:
        return 99.0


def special_subspaces(G):
    """Аффинные подпространства, неподвижные относительно элементов группы (и их пересечения):
    (x0, базис направлений) во фракционных координатах. Решаем (R − I) x = −t + n по модулю решётки."""
    subs = {}
    ops = [(R, t) for R, t in zip(G.R, G.t) if not (np.allclose(R, np.eye(3)) and np.allclose(t % 1.0, 0))]
    shifts = [np.array(n, float) for n in np.ndindex(3, 3, 3)]
    shifts = [n - 1 for n in shifts]

    def add(A, b):
        # решения A x = b: x0 + span(N)
        x0, *_ = np.linalg.lstsq(A, b, rcond=None)
        if np.linalg.norm(A @ x0 - b) > 1e-9:
            return
        u, sv, vt = np.linalg.svd(A)
        rank = int((sv > 1e-9).sum())
        N = vt[rank:]
        if N.shape[0] == 3:
            return
        # канонизация: x0 по модулю 1 с проекцией на ортогональное дополнение
        key = (N.shape[0], tuple(np.round(np.abs(N).sum(axis=0), 4)), tuple(np.round((x0 - N.T @ (N @ x0)) % 1.0, 4)))
        subs.setdefault(key, (x0 % 1.0, N))
    for R, t in ops:
        for n in shifts:
            add(R - np.eye(3), -t + n)
    # пересечения пар (точки и прямые высокой симметрии)
    base = list(subs.values())[:60]
    for (R1, t1), (R2, t2) in [(a, b) for a in ops[:24] for b in ops[:24]]:
        for n1 in shifts[::3]:
            for n2 in shifts[::3]:
                add(np.vstack([R1 - np.eye(3), R2 - np.eye(3)]), np.r_[-t1 + n1, -t2 + n2])
        if len(subs) > 400:
            break
    return list(subs.values())


def run(it, n_random, n_local, rng):
    hm = gemmi.find_spacegroup_by_number(it).hm
    G = Group(hm)
    cubic = G.system == "cubic"
    cache = {}
    cands = []
    for _ in range(n_random):
        x = rng.random(3)
        c = 1.0 if cubic else float(np.exp(rng.uniform(np.log(0.3), np.log(3.0))))
        cands.append((score(G, x, c, cache), x, c))
    # особые позиции: точки на подпространствах, неподвижных относительно элементов группы
    special = []
    for x0, N in special_subspaces(G):
        for _ in range(4 if N.shape[0] else 1):
            x = x0 + (rng.uniform(-0.5, 0.5, N.shape[0]) @ N if N.shape[0] else 0)
            c = 1.0 if cubic else float(np.exp(rng.uniform(np.log(0.3), np.log(3.0))))
            sc = score(G, x, c, cache)
            special.append((sc, x, c, x0, N))
    special.sort(key=lambda z: z[0])
    cands.sort(key=lambda z: z[0])
    best = min(cands[0], special[0][:3], key=lambda z: z[0]) if special else cands[0]
    # локальная оптимизация вдоль лучших особых подпространств
    for sc, x, c, x0, N in special[:n_local]:
        k = N.shape[0]
        if k == 0 and cubic:
            continue

        def fs(z, x0=x0, N=N, k=k):
            xx = x0 + (z[:k] @ N if k else 0)
            cc = 1.0 if cubic else float(np.exp(z[k]))
            return score(G, xx, cc, cache)
        z0 = np.r_[np.linalg.lstsq(N.T, x - x0, rcond=None)[0] if k else np.zeros(0), [] if cubic else [np.log(c)]]
        if z0.size == 0:
            continue
        r = minimize(fs, z0, method="Nelder-Mead", options={"xatol": 1e-8, "fatol": 1e-10, "maxfev": 1500})
        if r.fun < best[0]:
            xx = x0 + (r.x[:k] @ N if k else 0)
            best = (r.fun, xx % 1.0, 1.0 if cubic else float(np.exp(r.x[k])))
    for s0, x0, c0 in cands[:n_local]:
        if cubic:
            f = lambda z: score(G, z, 1.0, cache)  # noqa: E731
            r = minimize(f, x0, method="Nelder-Mead", options={"xatol": 1e-7, "fatol": 1e-9, "maxfev": 1500})
            cand = (r.fun, r.x % 1.0, 1.0)
        else:
            f = lambda z: score(G, z[:3], float(np.exp(z[3])), cache)  # noqa: E731
            r = minimize(f, np.r_[x0, np.log(c0)], method="Nelder-Mead", options={"xatol": 1e-7, "fatol": 1e-9, "maxfev": 2500})
            cand = (r.fun, r.x[:3] % 1.0, float(np.exp(r.x[3])))
        if cand[0] < best[0]:
            best = cand
    return {"IT": it, "group": hm, "best": round(best[0], 6), "ratio_to_TO": round(best[0] / TO, 6),
            "x": [round(float(v), 6) for v in best[1]], "c_over_a": round(best[2], 6)}


if __name__ == "__main__":
    lo, hi, n_random, n_local = (int(v) for v in sys.argv[1:5])
    rng = np.random.default_rng(lo)
    for it in range(lo, hi + 1):
        t0 = time.time()
        try:
            r = run(it, n_random, n_local, rng)
        except Exception as e:
            r = {"IT": it, "error": str(e)[:100]}
        r["sec"] = round(time.time() - t0, 1)
        print(json.dumps(r), flush=True)
