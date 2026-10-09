"""Пакетный поиск по группам 75–230. Безопасный режим: 6 процессов, nice 15, сторож памяти (1 ГБ на процесс)."""
import os, subprocess, sys, pathlib, time, gemmi
from concurrent.futures import ThreadPoolExecutor

PY = "/Users/vladimirmuravev/Projects2026/GEOTILES/.venv/bin/python"
WORKERS = 4
MAX_RSS_KB = 1_000_000
# аргументы: n_random n_local steps [c_lo:c_hi] [outdir] [nomirror]
out = pathlib.Path(sys.argv[5] if len(sys.argv) > 5 else "data/batch"); out.mkdir(parents=True, exist_ok=True)
CR = sys.argv[4] if len(sys.argv) > 4 else "0.3:4"
groups = [(n, gemmi.find_spacegroup_by_number(n).hm) for n in range(75, 231)]
if len(sys.argv) > 6 and sys.argv[6] == "nomirror":
    groups = [(n, hm) for n, hm in groups if n < 195 and "m" not in hm]


def rss_kb(pid):
    r = subprocess.run(["ps", "-o", "rss=", "-p", str(pid)], capture_output=True, text=True)
    return int(r.stdout.strip() or 0)


def run(g):
    n, hm = g
    log = out / f"g{n:03d}.log"
    if log.exists() and "BEST" in log.read_text():
        return
    with open(log, "w") as f:
        p = subprocess.Popen(["nice", "-n", "15", PY, "dvfast.py", hm, *sys.argv[1:4], CR], stdout=f, stderr=subprocess.STDOUT)
        t0 = time.time()
        while p.poll() is None:
            time.sleep(2)
            if rss_kb(p.pid) > MAX_RSS_KB or time.time() - t0 > 3600:
                p.kill()
                f.write(f"KILLED (rss or time limit)\n")
                break


with ThreadPoolExecutor(WORKERS) as ex:
    list(ex.map(run, groups))
