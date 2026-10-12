"""Рост зёрен (модель Поттса) в периодическом кубе: модель огрубления пены и роста зёрен металла.

Старт: разбиение Вороного случайных точек (каждое зерно — своя метка), затем Монте-Карло: воксель пробует
перекраситься в метку случайного из 26 соседей; энергия — число соседей с другой меткой; принимаем по
Метрополису с температурой kT. Обновления по 27 подрешёткам (шаг 3 по каждой оси) — соседи в одной подрешётке
не взаимодействуют. Снимки сохраняются в data/potts/ (не в git).
Запуск: python potts.py N n_seeds kT steps snapshot_steps(через запятую) seed
"""
import pathlib
import sys
import time

import numpy as np
from scipy.spatial import cKDTree

HERE = pathlib.Path(__file__).resolve().parent
OFFS = [(dx, dy, dz) for dx in (-1, 0, 1) for dy in (-1, 0, 1) for dz in (-1, 0, 1) if (dx, dy, dz) != (0, 0, 0)]


def voronoi_labels(N, n, rng):
    P = rng.random((n, 3)) * N
    tree = cKDTree(P, boxsize=N)
    g = (np.indices((N, N, N)).reshape(3, -1).T + 0.5) % N
    _, lab = tree.query(g, workers=1)
    return lab.reshape(N, N, N).astype(np.int32) + 1


def sweep(L, kT, rng):
    N = L.shape[0]
    for a in range(3):
        for b in range(3):
            for c in range(3):
                Lp = np.pad(L, 1, mode="wrap")
                C = L[a::3, b::3, c::3]
                sh = C.shape
                nb = np.empty((26,) + sh, np.int32)
                for k, (dx, dy, dz) in enumerate(OFFS):
                    nb[k] = Lp[1 + a + dx:1 + a + dx + 3 * sh[0]:3, 1 + b + dy:1 + b + dy + 3 * sh[1]:3,
                               1 + c + dz:1 + c + dz + 3 * sh[2]:3]
                pick = rng.integers(0, 26, sh)
                new = np.take_along_axis(nb, pick[None], axis=0)[0]
                dE = (nb == C[None]).sum(axis=0) - (nb == new[None]).sum(axis=0)
                acc = (dE <= 0) if kT == 0 else (rng.random(sh) < np.exp(-np.maximum(dE, 0) / kT))
                acc &= new != C
                C[acc] = new[acc]
                L[a::3, b::3, c::3] = C
    return L


if __name__ == "__main__":
    N, n_seeds, kT, steps = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
    snaps = {int(s) for s in sys.argv[5].split(",")}
    seed = int(sys.argv[6]) if len(sys.argv) > 6 else 0
    rng = np.random.default_rng(seed)
    out = HERE / "data" / "potts"
    out.mkdir(parents=True, exist_ok=True)
    L = voronoi_labels(N, n_seeds, rng)
    np.savez_compressed(out / f"pv_N{N}_n{n_seeds}_s{seed}.npz", L=L)  # контроль: Вороной случайных точек в вокселях
    t0 = time.time()
    for s in range(1, steps + 1):
        L = sweep(L, kT, rng)
        if s in snaps:
            ng = len(np.unique(L))
            np.savez_compressed(out / f"potts_N{N}_n{n_seeds}_kT{kT}_s{seed}_t{s}.npz", L=L)
            print(f"step {s}: grains {ng}, {time.time() - t0:.0f} s", flush=True)
