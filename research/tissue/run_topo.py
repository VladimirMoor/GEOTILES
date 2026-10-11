"""Прогон topology.py по всем семязачаткам PlantSeg: скачать → посчитать (с лимитом памяти) → удалить."""
import pathlib
import time

from run_all import HERE, RAW, PY, NODES, listing, download, run_guarded

if __name__ == "__main__":
  for attempt in range(4):  # проходы: неудачные загрузки пропускаем и повторяем в следующем проходе
    failed = 0
    for node, (kind, bg, vox) in NODES.items():
        if kind != "ovules":
            continue
        for name, url in listing(node):
            if (HERE / "data" / "topo" / (pathlib.Path(name).stem + ".json")).exists():
                continue
            dst = RAW / kind / name
            t0 = time.time()
            try:
                download(url, dst)
            except Exception as e:
                failed += 1
                dst.unlink(missing_ok=True)
                print(f"[skip] {name}: {e}", flush=True)
                continue
            out, peak = run_guarded(["nice", "-n", "10", PY, str(HERE / "topology.py"), str(dst)] + [str(v) for v in vox])
            dst.unlink()
            print(f"[{kind}] {name}  ({time.time() - t0:.0f} s, peak {peak} MB)\n{out.strip()}", flush=True)
    print(f"pass {attempt + 1}: failed {failed}", flush=True)
    if not failed:
        break
    time.sleep(120)
  print("ALL DONE", flush=True)
