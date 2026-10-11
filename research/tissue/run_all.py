"""Прогон по всем объёмам PlantSeg: скачать один файл → посчитать (cells.py, с лимитом памяти) → удалить.

На диске одновременно один файл (≤ 250 МБ). Процесс расчёта убивается, если его память превышает LIMIT_MB.
Готовые результаты (data/<имя>.json) пропускаются, так что прогон можно перезапускать.
"""
import json
import pathlib
import subprocess
import sys
import time
import urllib.request

HERE = pathlib.Path(__file__).resolve().parent
RAW = HERE / "data" / "raw"
PY = str(HERE.parents[1] / ".venv" / "bin" / "python")
LIMIT_MB = 3000
# узел OSF → (набор, метка фона, размер вокселя z y x в мкм)
# семязачатки ds2x: 0.235 × 0.150 × 0.150 мкм (PlantSeg models_zoo.yaml, confocal_3D_unet_ovules_ds2x)
OV = (0.235, 0.150, 0.150)
NODES = {"x9yns": ("ovules", 0, OV), "xp5uf": ("ovules", 0, OV), "8jz7e": ("ovules", 0, OV),
         "9x3g2": ("lrp", 1, (0.25, 0.1625, 0.1625)), "vs6gb": ("lrp", 1, (0.25, 0.1625, 0.1625)),
         "tn4xj": ("lrp", 1, (0.25, 0.1625, 0.1625))}


def listing(node):
    u, out = f"https://api.osf.io/v2/nodes/{node}/files/osfstorage/", []
    while u:
        d = json.load(urllib.request.urlopen(u, timeout=60))
        out += [(x["attributes"]["name"], x["links"]["download"]) for x in d["data"] if x["attributes"]["kind"] == "file"]
        u = d["links"].get("next")
    return out


def download(url, dst):
    dst.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=600) as r, open(dst, "wb") as w:
                while chunk := r.read(1 << 20):
                    w.write(chunk)
            return
        except Exception:
            time.sleep(10 * (attempt + 1))
    raise RuntimeError(f"download failed: {dst.name}")


def run_guarded(args):
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    peak = 0
    while p.poll() is None:
        try:
            rss = int(subprocess.run(["ps", "-o", "rss=", "-p", str(p.pid)], capture_output=True, text=True).stdout.strip() or 0) // 1024
        except ValueError:
            rss = 0
        peak = max(peak, rss)
        if rss > LIMIT_MB:
            p.kill()
            p.wait()
            return f"KILLED: memory {rss} MB > {LIMIT_MB} MB", peak
        time.sleep(0.5)
    return p.stdout.read(), peak


if __name__ == "__main__":
    for node, (kind, bg, vox) in NODES.items():
        for name, url in listing(node):
            if (HERE / "data" / (pathlib.Path(name).stem + ".json")).exists():
                continue
            dst = RAW / kind / name
            t0 = time.time()
            download(url, dst)
            args = ["nice", "-n", "10", PY, str(HERE / "cells.py"), str(dst), str(bg)]
            if vox:
                args += [str(v) for v in vox]
            out, peak = run_guarded(args)
            dst.unlink()
            print(f"[{kind}] {name}  ({time.time() - t0:.0f} s, peak {peak} MB)\n{out.strip()}", flush=True)
    print("ALL DONE", flush=True)
