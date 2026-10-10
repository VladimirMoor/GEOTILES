"""Запуск heesch-sat со сторожем памяти и времени. Использование: guarded_sat.py in.txt out.txt maxlevel max_gb max_hours [extra flags]"""
import subprocess, sys, time
SAT = "/Users/vladimirmuravev/Projects2026/GEOTILES/tools/heesch-sat/src/sat"
inp, out, lvl, gb, hours = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]), float(sys.argv[5])
extra = sys.argv[6:]
t0 = time.time(); peak = 0
with open(out, "w") as fo:
    p = subprocess.Popen(["nice", "-n", "10", SAT, *extra, "-maxlevel", lvl, inp], stdout=fo, stderr=subprocess.STDOUT)
    while p.poll() is None:
        time.sleep(0.5)
        r = subprocess.run(["ps", "-o", "rss=", "-p", str(p.pid)], capture_output=True, text=True)
        rss = int(r.stdout.strip() or 0); peak = max(peak, rss)
        if rss > gb * 1e6:
            p.kill(); print(f"KILLED: memory {rss/1e6:.1f} GB"); break
        if time.time() - t0 > hours * 3600:
            p.kill(); print("KILLED: time"); break
print(f"done in {time.time()-t0:.0f}s, peak {peak/1e6:.2f} GB, exit {p.returncode}")
print(open(out).read().splitlines()[1] if len(open(out).read().splitlines()) > 1 else open(out).read()[:300])
