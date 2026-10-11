"""Топология клеток в 3D-тканях (сегментации PlantSeg) и сравнение с разбиением Вороного центров клеток.

Для объёма меток (HDF5, датасет 'label'):
1. Контакты: пары соседних вокселей с разными метками по трём осям → площадь контакта для каждой пары клеток
   (в мкм², с учётом размера вокселя). Считается слоями по z, чтобы не держать лишних копий массива.
2. Граничные клетки: касаются фона, края объёма или неразмеченной области — их статистику не берём
   (как соседи они учитываются).
3. Грани клетки: соседи с площадью контакта ≥ frac · (площадь поверхности клетки). Порог убирает шум сегментации;
   результат приводится для нескольких порогов.
4. Вороной: 3D-разбиение Вороного центров масс всех клеток; грани ячейки — соседи с площадью общей грани
   ≥ того же порога. Сравнение множеств соседей для внутренних клеток (совпадение, индекс Жаккара).
Запуск: python cells.py путь.h5 [фон] [vx vy vz]  →  data/<имя>.json
"""
import json
import pathlib
import sys

import h5py
import numpy as np
from scipy.spatial import Voronoi

HERE = pathlib.Path(__file__).resolve().parent
FRACS = (0.0, 0.01, 0.02, 0.05)


def contacts(L, vox, bg, ign):
    """Словарь (a, b) → площадь контакта, a < b; плюс множество граничных меток."""
    vz, vy, vx = vox
    area = {0: vy * vx, 1: vz * vx, 2: vz * vy}  # площадь грани вокселя, перпендикулярной оси
    acc = {}
    border = set()
    for ax in range(3):
        n = L.shape[ax]
        step = 32
        for s in range(0, n - 1, step):
            e = min(s + step + 1, n)
            sl = [slice(None)] * 3
            sl[ax] = slice(s, e)
            B = L[tuple(sl)]
            a = np.moveaxis(B, ax, 0)
            p, q = a[:-1].ravel(), a[1:].ravel()
            m = p != q
            p, q = p[m].astype(np.int64), q[m].astype(np.int64)
            lo, hi = np.minimum(p, q), np.maximum(p, q)
            key = lo * 70000 + hi
            u, c = np.unique(key, return_counts=True)
            for k, cnt in zip(u.tolist(), c.tolist()):
                acc[k] = acc.get(k, 0.0) + cnt * area[ax]
    pairs = {}
    for k, ar in acc.items():
        a, b = divmod(k, 70000)
        if a in bg or b in bg or a in ign or b in ign:
            border.update(x for x in (a, b) if x not in bg and x not in ign)
            continue
        pairs[(a, b)] = ar
    # клетки на краю объёма
    for ax in range(3):
        for idx in (0, L.shape[ax] - 1):
            border.update(np.unique(np.take(L, idx, axis=ax)).tolist())
    return pairs, border - bg - ign


def poly_area(P):
    c = P.mean(axis=0)
    return 0.5 * np.linalg.norm(sum(np.cross(P[i] - c, P[(i + 1) % len(P)] - c) for i in range(len(P))))


def analyse(path, bg_label, vox):
    h = h5py.File(path, "r")
    L = h["label"][:]
    ign = set()
    if "label_with_ignore" in h:  # неразмеченные клетки (-1) превращаем в отдельную метку; послойно, без копий
        Iw = h["label_with_ignore"]
        for z0 in range(0, L.shape[0], 16):
            I = Iw[z0:z0 + 16]
            L[z0:z0 + 16][I == -1] = 65000
        del I
        ign = {65000}
    bg = {0, bg_label}
    pairs, border = contacts(L, vox, bg, ign)
    labels = np.array([l for l in np.unique(L).tolist() if l not in bg and l not in ign])
    idx = {l: i for i, l in enumerate(labels)}
    # объёмы и центры масс — послойно по z, без больших временных массивов
    M = int(L.max()) + 1
    cnt = np.zeros(M); sz = np.zeros(M); sy = np.zeros(M); sx = np.zeros(M)
    yy, xx = np.indices(L.shape[1:])
    for z in range(L.shape[0]):
        f = L[z].ravel()
        c = np.bincount(f, minlength=M)
        cnt += c; sz += c * z
        sy += np.bincount(f, weights=yy.ravel(), minlength=M)
        sx += np.bincount(f, weights=xx.ravel(), minlength=M)
    vol = cnt[labels] * np.prod(vox)
    cen = np.stack([sz[labels], sy[labels], sx[labels]], axis=1) / cnt[labels][:, None] * np.array(vox)
    del L
    # полная поверхность клетки (вместе с контактом с фоном) — для относительного порога берём только межклеточную
    surf = np.zeros(len(labels))
    for (a, b), ar in pairs.items():
        surf[idx[a]] += ar
        surf[idx[b]] += ar
    # Вороной центров
    vor = Voronoi(cen)
    vpairs = {}
    for (i, j), rv in zip(vor.ridge_points, vor.ridge_vertices):
        if -1 in rv:
            continue
        vpairs[(min(i, j), max(i, j))] = poly_area(vor.vertices[rv])
    vsurf = np.zeros(len(labels))
    for (i, j), ar in vpairs.items():
        vsurf[i] += ar
        vsurf[j] += ar
    # «Внутренние» клетки: не на границе ткани, и все их соседи по Вороному тоже имеют конечные ячейки
    inner = [i for i, l in enumerate(labels) if l not in border]
    res = {"file": pathlib.Path(path).name, "cells": len(labels), "inner": len(inner), "voxel": vox, "by_frac": {}}
    for fr in FRACS:
        nb = [set() for _ in labels]
        for (a, b), ar in pairs.items():
            i, j = idx[a], idx[b]
            if ar >= fr * surf[i]:
                nb[i].add(j)
            if ar >= fr * surf[j]:
                nb[j].add(i)
        vb = [set() for _ in labels]
        for (i, j), ar in vpairs.items():
            if ar >= fr * vsurf[i]:
                vb[i].add(j)
            if ar >= fr * vsurf[j]:
                vb[j].add(i)
        F = [len(nb[i]) for i in inner]
        VF = [len(vb[i]) for i in inner]
        jac = [len(nb[i] & vb[i]) / max(1, len(nb[i] | vb[i])) for i in inner]
        res["by_frac"][str(fr)] = {
            "faces": F, "vor_faces": VF, "mean_faces": float(np.mean(F)), "mean_vor_faces": float(np.mean(VF)),
            "jaccard_mean": float(np.mean(jac)), "same_set": float(np.mean([nb[i] == vb[i] for i in inner])),
            "volume": [float(vol[i]) for i in inner]}
    return res


if __name__ == "__main__":
    path = sys.argv[1]
    bg_label = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    vox = tuple(float(v) for v in sys.argv[3:6]) if len(sys.argv) > 5 else (1.0, 1.0, 1.0)
    r = analyse(path, bg_label, vox)
    out = HERE / "data" / (pathlib.Path(path).stem + ".json")
    out.write_text(json.dumps(r))
    print(r["file"], "cells", r["cells"], "inner", r["inner"])
    for fr, d in r["by_frac"].items():
        print(f"  frac {fr:5s}: faces {d['mean_faces']:.2f}  voronoi {d['mean_vor_faces']:.2f}  "
              f"jaccard {d['jaccard_mean']:.3f}  identical sets {d['same_set']:.2%}")
