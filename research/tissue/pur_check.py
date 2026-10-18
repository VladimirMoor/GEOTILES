"""Check of the PUR cell segmentation against the struts: in a real open-cell foam every cell edge is a strut, so
triple lines of the label volume (2x2x2 cubes with three or more labels) should lie on the struts. Reports, for
several marker depths h, the share of triple-line voxels within d voxels of the solid. Usage: python pur_check.py h..."""
import sys

import numpy as np
from scipy import ndimage as ndi

from pur_seg import load, segment

Z, CROP = (150, 450), (175, 575, 175, 575)
S = load(*Z, CROP)
Dsol = ndi.distance_transform_edt(~S)
for h in [float(x) for x in sys.argv[1:]]:
    L, info = segment(S, (CROP[0], CROP[2]), 360.0, h)
    V = np.stack([L[dz:L.shape[0] - 1 + dz, dy:L.shape[1] - 1 + dy, dx:L.shape[2] - 1 + dx]
                  for dz in (0, 1) for dy in (0, 1) for dx in (0, 1)], axis=-1)
    V.sort(axis=-1)
    nl = 1 + (V[..., 1:] != V[..., :-1]).sum(-1)
    tri = (nl >= 3) & (V[..., 0] > 0)
    d = Dsol[:-1, :-1, :-1][tri]
    # same test for faces (two labels): faces cross open windows, so most face voxels are away from struts
    fac = (nl == 2) & (V[..., 0] > 0)
    # converse: every strut should be an edge of the tessellation (merged cells leave struts inside a cell)
    T = np.zeros(S.shape, bool)
    T[:-1, :-1, :-1] = tri
    Dtri = ndi.distance_transform_edt(~T)
    inner = S & (Dsol == 0)
    cov = np.mean(Dtri[inner] <= 3)
    print(f"  struts within 3 vox of a triple line: {cov:.3f}")
    print(f"h={h}: cells {info['cells']}, triple-line voxels {tri.sum()}, within 2 vox of strut {np.mean(d <= 2):.3f}, "
          f"within 4 vox {np.mean(d <= 4):.3f}; face voxels within 2 vox {np.mean(Dsol[:-1, :-1, :-1][fac] <= 2):.3f}", flush=True)
