"""Cells of an open-cell polyurethane foam (DaRUS doi:10.18419/darus-3010, t0 = undeformed, segmented stack).

The data are binary (strut skeleton = 255, void = 0), 593 slices of 751 x 751, voxel ≈ 0.07 mm, cylinder 52.4 mm.
The cells have no membranes, only struts along their edges, so cells are separated where the void narrows at the
windows: markers at maxima of the distance to the struts (h-maxima, to merge plateaus in elongated cells),
watershed on −distance inside the cylinder. Struts are absorbed into the cells, so struts become triple lines,
as in the tissue and foam tessellations, and topology.py applies unchanged. Label 0 = outside (ignored).

Usage: python pur_seg.py OUT.npz [--z 0 593] [--crop y0 y1 x0 x1] [--h 2.0] [--r 360]
"""
import argparse
import pathlib

import numpy as np
from PIL import Image
from scipy import ndimage as ndi
from skimage.morphology import h_maxima
from skimage.segmentation import watershed

HERE = pathlib.Path(__file__).resolve().parent
SRC = HERE / "data" / "pur" / "raw" / "t0_segmented"
FILES = sorted(SRC.glob("slice_xy_*.tif"))
N = 751


def load(z0, z1, crop):
    y0, y1, x0, x1 = crop
    return np.stack([np.array(Image.open(f))[y0:y1, x0:x1] > 0 for f in FILES[z0:z1]])


def segment(S, origin, r_mask, h, sigma=1.0):
    yy, xx = np.ogrid[:S.shape[1], :S.shape[2]]
    cy, cx = (N - 1) / 2 - origin[0], (N - 1) / 2 - origin[1]
    disk = (yy - cy) ** 2 + (xx - cx) ** 2 < r_mask ** 2
    mask = np.broadcast_to(disk, S.shape)
    D = ndi.distance_transform_edt(~S).astype(np.float32)
    Ds = ndi.gaussian_filter(D, sigma)
    peak = h_maxima(Ds, h).astype(bool) & mask
    markers, n = ndi.label(peak, structure=np.ones((3, 3, 3)))
    L = watershed(-Ds, markers, mask=mask)
    return L, {"cells": int(n), "solid_fraction": float(S[mask].mean()), "mean_peak_D": float(Ds[peak].mean())}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--z", type=int, nargs=2, default=[0, len(FILES)])
    ap.add_argument("--crop", type=int, nargs=4, default=[0, N, 0, N])
    ap.add_argument("--h", type=float, default=2.0)
    ap.add_argument("--r", type=float, default=360.0)
    a = ap.parse_args()
    S = load(*a.z, a.crop)
    L, info = segment(S, (a.crop[0], a.crop[2]), a.r, a.h)
    info.update(z=a.z, crop=a.crop, h=a.h, r_mask=a.r, shape=list(L.shape))
    np.savez_compressed(a.out, L=L.astype(np.uint16), info=str(info))
    print(info, flush=True)
