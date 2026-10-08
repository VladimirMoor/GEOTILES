"""Строгая проверка в Singular (детерминированные базисы Грёбнера над ℚ).

Для типа строится идеал I (как в exact.py) и насыщается по многочленам невырожденности, которые не
обращаются в ноль ни на одной строго выпуклой реализации:
  |x_u − x_v|²  для несмежных вершин u, v (у выпуклого тела вершины различны и вещественны,
               поэтому |x_u − x_v|² > 0).
Насыщение I : D^∞ выкидывает все компоненты, лежащие в {D = 0}; вещественные выпуклые решения
остаются. Затем считается размерность насыщенного идеала:
  -1 (идеал = ⟨1⟩)  ⇒ строгое доказательство: выпуклых реализаций нет;
   0                ⇒ конечное число решений, вещественные изолирует msolve;
  >0                ⇒ семейство (или невыпуклые компоненты) — разбирается отдельно.

Запуск: python singular.py F6-001 F6-003 [--timeout 300]
"""
import argparse
import itertools
import json
import pathlib
import re
import subprocess
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from exact import build  # noqa: E402
from realize import edges_of  # noqa: E402


def singular_script(t, variables, eqs, chunk=6, basis_path="basis.txt"):
    E = set(edges_of(t["faces"]))
    a, b, c = t["faces"][0][:3]

    def X(v, k):
        if v == a:
            return "0"
        if v == b:
            return "1" if k == 0 else "0"
        if v == c and k == 2:
            return "0"
        return f"x{v}_{k}"

    nonadj = [(u, v) for u, v in itertools.combinations(range(t["nV"]), 2) if (u, v) not in E]
    dist = ["+".join(f"({X(u, k)}-{X(v, k)})^2" for k in range(3)) for u, v in nonadj]
    vs = ",".join(v.replace("_", "") for v in variables)
    fix = lambda s: re.sub(r"([xn]\d+)_(\d)", r"\1\2", s)
    lines = [
        'LIB "elim.lib";',
        f"ring r = 0, ({vs}), dp;",
        "option(redSB);",
        f"ideal I = {', '.join(fix(e) for e in eqs)};",
        "I = std(I);",
        'print("dim_before " + string(dim(I)));',
    ]
    # насыщаем порциями: произведение нескольких расстояний за раз
    for i in range(0, len(dist), chunk):
        prod = "*".join(f"({fix(d)})" for d in dist[i:i + chunk])
        lines.append(f"def S = sat(I, {prod}); if (typeof(S) == \"list\") {{ I = S[1]; }} else {{ I = S; }} kill S;")
        lines.append("if (dim(std(I)) == -1) { print(\"dim_after -1\"); quit; }")
    lines += ["I = std(I);", 'print("dim_after " + string(dim(I)));',
              'if (dim(I) == 0) { print("vdim " + string(vdim(I))); }',
              f'write(":w {basis_path}", I);', "quit;"]
    return "\n".join(lines)


def run(t, timeout):
    variables, eqs = build(t)
    (ROOT / "data" / "sing").mkdir(exist_ok=True)
    script = singular_script(t, variables, eqs, basis_path=str(ROOT / "data" / "sing" / f"{t['id']}.basis"))
    path = ROOT / "data" / "sing" / f"{t['id']}.sing"
    path.write_text(script)
    t0 = time.time()
    try:
        p = subprocess.run(["Singular", "-q", str(path)], capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return {"status": "timeout", "time": timeout}
    out = p.stdout
    res = {"time": round(time.time() - t0, 2)}
    for key in ("dim_before", "dim_after", "vdim"):
        m = re.search(key + r" (-?\d+)", out)
        if m:
            res[key] = int(m.group(1))
    if "dim_after" not in res:
        res["status"] = "error"
        res["stderr"] = (p.stderr or out)[-300:]
    return res


def _job(args):
    t, timeout = args
    return t["id"], run(t, timeout)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ids", nargs="*")
    ap.add_argument("--faces", default=None)
    ap.add_argument("--timeout", type=int, default=300)
    ap.add_argument("--out", default=None)
    ap.add_argument("--jobs", type=int, default=10)
    args = ap.parse_args()
    types = json.loads((ROOT / "data" / "types.json").read_text())
    if args.faces:
        lo, hi = map(int, args.faces.split("-"))
        types = [t for t in types if lo <= t["nF"] <= hi]
    if args.ids:
        types = [t for t in types if t["id"] in args.ids]
    from multiprocessing import Pool
    results = {}
    with Pool(args.jobs) as pool:
        for tid, r in pool.imap_unordered(_job, [(t, args.timeout) for t in types]):
            results[tid] = r
            print(tid, r, flush=True)
    if args.out:
        (ROOT / "data" / args.out).write_text(json.dumps(results, indent=1))


if __name__ == "__main__":
    main()
