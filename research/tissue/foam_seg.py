"""Bubble segmentation of the real-foam tomogram (slices from foam_fetch.py) into a space-filling label volume.

Same idea as FoamQuant (Schott et al. 2025): threshold the liquid phase, place one marker per bubble at a
maximum of the distance to the liquid, and grow the markers by watershed on -distance. The growth also fills
the liquid, so Plateau borders become triple lines between bubble labels, as in a cell tissue, and the
tissue topology pipeline (topology.py) applies unchanged. Label 0 = outside the sample (ignored).

Usage: python foam_seg.py OUT.npz [--bin 2] [--z 100 510] [--crop y0 y1 x0 x1] [--r 940]
"""
import argparse
import pathlib

import numpy as np
from scipy import ndimage as ndi
from skimage.filters import threshold_otsu
from skimage.segmentation import watershed

HERE = pathlib.Path(__file__).resolve().parent
SL = HERE / "data" / "foam" / "slices"
N = 2016


def load(z0, z1, crop, b):
    y0, y1, x0, x1 = crop
    y1 -= (y1 - y0) % b
    x1 -= (x1 - x0) % b
    z1 -= (z1 - z0) % b
    out = np.zeros(((z1 - z0) // b, (y1 - y0) // b, (x1 - x0) // b), np.float32)
    for k, z in enumerate(range(z0, z1)):
        s = np.load(SL / f"z{z:03d}.npy")[y0:y1, x0:x1].astype(np.float32)
        out[k // b] += s.reshape((y1 - y0) // b, b, (x1 - x0) // b, b).mean(axis=(1, 3))
    return out / b, (y0, x0)


def segment(V, origin, b, r_mask, min_dist, sigma):
    V = ndi.gaussian_filter(V, sigma)
    yy, xx = np.ogrid[:V.shape[1], :V.shape[2]]
    cy, cx = (N / 2 - origin[0]) / b, (N / 2 - origin[1]) / b
    disk = (yy - cy) ** 2 + (xx - cx) ** 2 < (r_mask / b) ** 2
    mask = np.broadcast_to(disk, V.shape)
    t = threshold_otsu(V[mask][::7])
    gas = (V < t) & mask
    D = ndi.distance_transform_edt(gas).astype(np.float32)
    # one marker per bubble: plateau maxima of D above a floor, separated by at least min_dist
    peak = (D == ndi.maximum_filter(D, size=2 * min_dist + 1)) & (D > min_dist / 2)
    markers, n = ndi.label(peak, structure=np.ones((3, 3, 3)))
    L = watershed(-D, markers, mask=mask)
    return L, {"threshold": float(t), "liquid_fraction": float(1 - gas.sum() / mask.sum()), "bubbles": int(n),
               "mean_bubble_radius_vox": float(np.mean(D[peak]))}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--bin", type=int, default=2)
    ap.add_argument("--z", type=int, nargs=2, default=[100, 510])
    ap.add_argument("--crop", type=int, nargs=4, default=[0, N, 0, N])
    ap.add_argument("--r", type=float, default=940.0, help="sample mask radius, full-resolution voxels")
    ap.add_argument("--min-dist", type=int, default=None, help="marker separation in binned voxels")
    ap.add_argument("--sigma", type=float, default=1.0)
    a = ap.parse_args()
    md = a.min_dist or max(2, round(12 / a.bin))
    V, origin = load(a.z[0], a.z[1], a.crop, a.bin)
    L, info = segment(V, origin, a.bin, a.r, md, a.sigma)
    assert L.max() < 65000, "labels must fit the 16-bit keys of topology.py"
    info.update(bin=a.bin, z=a.z, crop=a.crop, r_mask=a.r, min_dist=md, sigma=a.sigma, shape=list(L.shape))
    np.savez_compressed(a.out, L=L.astype(np.uint16), info=str(info))
    print(info, flush=True)
