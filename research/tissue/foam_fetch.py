"""Download a z-range of the real-foam tomogram (Zenodo 15564270, Schott et al. 2025, series 2, step 102).

The record holds one tar with two BigTIFFs (steps 102 and 103), each 610 x 2016 x 2016 float64,
uncompressed, slices stored contiguously. Values are 0..255, so every slice is stored here as uint8.
HTTP range requests fetch only the slices we need. Usage: python foam_fetch.py Z0 Z1 [workers]
"""
import concurrent.futures as cf
import pathlib
import sys
import time
import urllib.request

import numpy as np

URL = "https://zenodo.org/api/records/15564270/files/tomo_rheoscopy_102_103.tar/content"
N = 2016
SLICE = N * N * 8
BASE = 512 + 448  # tar header + TIFF data offset of slice 0 in Tomo_rheoscopy_102.tif
OUT = pathlib.Path(__file__).resolve().parent / "data" / "foam" / "slices"


def fetch(z):
    out = OUT / f"z{z:03d}.npy"
    if out.exists():
        return z, 0.0
    a = BASE + z * SLICE
    buf = bytearray()  # a broken transfer resumes from the bytes already received
    t = time.time()
    for attempt in range(12):
        try:
            req = urllib.request.Request(URL, headers={"Range": f"bytes={a + len(buf)}-{a + SLICE - 1}"})
            with urllib.request.urlopen(req, timeout=600) as r:
                while chunk := r.read(1 << 20):
                    buf += chunk
            if len(buf) != SLICE:
                raise IOError(f"short read {len(buf)}")
            s = np.frombuffer(buf, dtype="<f8").reshape(N, N)
            np.save(out.with_suffix(".tmp.npy"), np.clip(np.rint(s), 0, 255).astype(np.uint8))
            out.with_suffix(".tmp.npy").rename(out)
            return z, time.time() - t
        except Exception as e:  # noqa: BLE001
            print(f"z={z} attempt {attempt}: {e} (have {len(buf)} bytes)", flush=True)
            time.sleep(5 * (attempt + 1))
    raise RuntimeError(f"z={z} failed")


if __name__ == "__main__":
    z0, z1 = int(sys.argv[1]), int(sys.argv[2])
    workers = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    OUT.mkdir(parents=True, exist_ok=True)
    done = 0
    with cf.ThreadPoolExecutor(workers) as ex:
        for z, dt in ex.map(fetch, range(z0, z1)):
            done += 1
            if done % 10 == 0:
                print(f"{done}/{z1 - z0} slices (last z={z}, {dt:.0f}s)", flush=True)
    print("DONE", flush=True)
