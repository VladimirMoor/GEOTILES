"""Пилот поиска не-вороновских стереоэдров: старт от ячейки Дирихле со случайным шумом в плоскостях,
минимизация f (Пауэлл). Пишет строки JSON: группа, шум, f, число граней, отклонение от Дирихле."""
import json
import sys
import time

import numpy as np
from scipy.optimize import minimize

from nondv import Search
from nondv_verify import verify

hm, sigmas, runs, seed0 = sys.argv[1], [float(v) for v in sys.argv[2].split(",")], int(sys.argv[3]), int(sys.argv[4])
for sig in sigmas:
    for r in range(runs):
        seed = seed0 + 1000 * r
        S = Search(hm, seed=seed)
        rng = np.random.default_rng(seed + 7)
        x = S.B @ rng.random(3)
        p0 = S.dv_params(x)
        scale = np.abs(p0).mean()
        p0 = p0 + rng.normal(0, sig * scale, p0.size)
        t0 = time.time()
        res = minimize(S.f, p0, method="Powell", options={"maxfev": 6000, "xtol": 1e-9, "ftol": 1e-12})
        out = {"group": hm, "sigma": sig, "seed": seed, "f": float(res.fun), "nfev": int(res.nfev), "sec": round(time.time() - t0, 1),
               "pairs": len(S.pairs)}
        if res.fun < 1e-5:
            n, v, _, _ = S.facets(res.x)
            out.update({"facets": n, "vol_ratio": v / S.V0, "verify": verify(S, res.x, 1500)})
            np.save(f"data/nondv/{hm.replace(' ', '')}_{sig}_{seed}.npy", res.x)
        print(json.dumps(out), flush=True)
