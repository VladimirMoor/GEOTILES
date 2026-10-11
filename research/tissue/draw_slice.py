"""Иллюстрация: срез настоящей ткани и срез разбиения Вороного центров тех же клеток.

Берём объём меток PlantSeg (семязачаток), считаем 3D-центры клеток (в мкм, с учётом анизотропии вокселя),
и в одной плоскости z красим каждый пиксель: слева — по настоящей клетке, справа — по ближайшему центру
(это сечение 3D-разбиения Вороного). Цвета у одной и той же клетки одинаковые, так что видно, где форма
совпадает, а где нет. Показываются только пиксели внутри ткани.
Запуск: python draw_slice.py путь.h5 z → docs/img/tissue-voronoi.png
"""
import pathlib
import sys

import h5py
import matplotlib
import numpy as np
from scipy.spatial import cKDTree

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

OUT = pathlib.Path(__file__).resolve().parents[2] / "docs" / "img" / "tissue-voronoi.png"
VOX = np.array([0.235, 0.150, 0.150])

path, zc = sys.argv[1], int(sys.argv[2])
h = h5py.File(path, "r")
L = h["label"][:]
if "label_with_ignore" in h:
    L[h["label_with_ignore"][:] == -1] = 0
labels = np.unique(L)
labels = labels[labels > 0]
M = int(L.max()) + 1
cnt = np.zeros(M); s = np.zeros((M, 3))
yy, xx = np.indices(L.shape[1:])
for z in range(L.shape[0]):
    f = L[z].ravel()
    cnt += np.bincount(f, minlength=M)
    s[:, 0] += np.bincount(f, minlength=M) * z
    s[:, 1] += np.bincount(f, weights=yy.ravel(), minlength=M)
    s[:, 2] += np.bincount(f, weights=xx.ravel(), minlength=M)
cen = s[labels] / cnt[labels][:, None] * VOX
S = L[zc]
mask = S > 0
pts = np.stack([np.full(S.size, zc), yy.ravel(), xx.ravel()], axis=1) * VOX
_, nn = cKDTree(cen).query(pts[mask.ravel()])
V = np.zeros_like(S)
V[mask] = labels[nn]

rng = np.random.default_rng(5)
pal = rng.uniform(0.35, 0.9, size=(M, 3))
pal[0] = (1, 1, 1)


def img(A):
    rgba = np.ones(A.shape + (4,))
    rgba[..., :3] = pal[A]
    rgba[A == 0, 3] = 0  # фон прозрачный
    # границы клеток тёмными линиями
    edge = np.zeros(A.shape, bool)
    edge[:-1] |= A[:-1] != A[1:]
    edge[:, :-1] |= A[:, :-1] != A[:, 1:]
    edge &= A > 0
    rgba[edge, :3] = 0.15
    rgba[edge, 3] = 1
    return rgba


ys, xs = np.where(mask)
y0, y1, x0, x1 = ys.min(), ys.max() + 1, xs.min(), xs.max() + 1
agree = (V[mask] == S[mask]).mean()
fig, ax = plt.subplots(1, 2, figsize=(10, 10 * (y1 - y0) / (2 * (x1 - x0)) + 0.2), dpi=130)
for a, A in zip(ax, (S, V)):
    a.imshow(img(A[y0:y1, x0:x1]), interpolation="nearest", aspect=VOX[1] / VOX[2])
    a.axis("off")
fig.subplots_adjust(0, 0, 1, 1, wspace=0.03)
fig.savefig(OUT, transparent=True)
print(OUT, "slice", zc, "cells in slice", len(np.unique(S[mask])), f"pixels in the same cell: {agree:.1%}")
