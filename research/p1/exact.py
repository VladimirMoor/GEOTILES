"""Строгая часть: точная система уравнений над ℚ для каждого комбинаторного типа и решение msolve.

Переменные: координаты вершин x_v (после калибровки: x_a = 0, x_b = (1,0,0), z_c = 0 для трёх
последовательных вершин первой грани) и нормали n_f граней с ≥ 5 вершинами.
Уравнения:
  |x_u − x_v|² − 1 = 0                  для каждого ребра;
  x_a − x_b + x_c − x_d = 0             для каждой четырёхугольной грани (плоский равносторонний
                                        четырёхугольник — ромб, т. е. параллелограмм);
  n_f · (x_v − x_w) = 0, r · n_f = 1     для граней с ≥ 5 вершинами (r — фиксированный рациональный вектор).
msolve вычисляет базис Грёбнера над ℚ и изолирует все вещественные решения рациональными интервалами;
выпуклость затем проверяется на этих интервалах.

Запуск: python exact.py F6-005 [...]   или   python exact.py --faces 4-6
"""
import argparse
import json
import pathlib
import subprocess
import sys
import time
from fractions import Fraction

ROOT = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from realize import edges_of  # noqa: E402

R_VEC = (3, 5, 7)  # нормировка нормалей: 3 n1 + 5 n2 + 7 n3 = 1 (вторая серия — другой вектор для контроля)


def build(t, r=R_VEC, rabinowitsch=False):
    faces = t["faces"]
    V = t["nV"]
    a, b, c = faces[0][:3]
    E = edges_of(faces)

    def X(v, k):
        if v == a:
            return "0"
        if v == b:
            return "1" if k == 0 else "0"
        if v == c and k == 2:
            return "0"
        return f"x{v}_{k}"

    variables = []
    for v in range(V):
        for k in range(3):
            s = X(v, k)
            if s not in ("0", "1"):
                variables.append(s)
    eqs = []
    for u, w in E:
        terms = [f"({X(u, k)}-{X(w, k)})^2" for k in range(3)]
        eqs.append("+".join(terms) + "-1")
    for f, face in enumerate(faces):
        if len(face) == 4:
            p, q, s, u = face
            for k in range(3):
                eqs.append(f"{X(p, k)}-{X(q, k)}+{X(s, k)}-{X(u, k)}")
        elif len(face) >= 5:
            n = [f"n{f}_{k}" for k in range(3)]
            variables += n
            v0 = face[0]
            for w in face[1:]:
                eqs.append("+".join(f"{n[k]}*({X(w, k)}-{X(v0, k)})" for k in range(3)))
            if rabinowitsch:
                eqs.append(f"{n[0]}^2+{n[1]}^2+{n[2]}^2-1")  # без выделенного направления
            else:
                eqs.append(f"{r[0]}*{n[0]}+{r[1]}*{n[1]}+{r[2]}*{n[2]}-1")
        if rabinowitsch and len(face) >= 4:
            # s_f · |(x1 − x0) × (x2 − x0)|² = 1: для вещественных решений ровно «x0, x1, x2 не на одной прямой»
            p0, p1, p2 = face[0], face[1], face[2]
            u = [f"({X(p1, k)}-{X(p0, k)})" for k in range(3)]
            w = [f"({X(p2, k)}-{X(p0, k)})" for k in range(3)]
            cr = [f"({u[1]}*{w[2]}-{u[2]}*{w[1]})", f"({u[2]}*{w[0]}-{u[0]}*{w[2]})", f"({u[0]}*{w[1]}-{u[1]}*{w[0]})"]
            sv = f"s{f}"
            variables.append(sv)
            eqs.append(f"{sv}*({cr[0]}^2+{cr[1]}^2+{cr[2]}^2)-1")
    if rabinowitsch:
        # для каждого ребра (u, w): грани по обе стороны не копланарны, т. е.
        # det(x_w − x_u, x_p − x_u, x_q − x_u) ≠ 0, где p, q — третьи вершины этих граней
        ef = {}
        for f, face in enumerate(faces):
            for i, u0 in enumerate(face):
                w0 = face[(i + 1) % len(face)]
                ef.setdefault((min(u0, w0), max(u0, w0)), []).append(face)
        for ei, (u0, w0) in enumerate(E):
            f1, f2 = ef[(u0, w0)]
            p = next(v for v in f1 if v not in (u0, w0))
            q = next(v for v in f2 if v not in (u0, w0))
            A = [f"({X(w0, k)}-{X(u0, k)})" for k in range(3)]
            B = [f"({X(p, k)}-{X(u0, k)})" for k in range(3)]
            C = [f"({X(q, k)}-{X(u0, k)})" for k in range(3)]
            det = (f"({A[0]}*({B[1]}*{C[2]}-{B[2]}*{C[1]})-{A[1]}*({B[0]}*{C[2]}-{B[2]}*{C[0]})"
                   f"+{A[2]}*({B[0]}*{C[1]}-{B[1]}*{C[0]}))")
            variables.append(f"e{ei}")
            eqs.append(f"e{ei}*{det}-1")
    if rabinowitsch == "full":
        # любые три последовательные вершины грани не коллинеарны (угол грани ≠ 0°, 180°)
        for f, face in enumerate(faces):
            if len(face) < 4:
                continue
            k = len(face)
            for i in range(1, k):
                p0, p1, p2 = face[i], face[(i + 1) % k], face[(i + 2) % k]
                u = [f"({X(p1, j)}-{X(p0, j)})" for j in range(3)]
                w = [f"({X(p2, j)}-{X(p0, j)})" for j in range(3)]
                cr = [f"({u[1]}*{w[2]}-{u[2]}*{w[1]})", f"({u[2]}*{w[0]}-{u[0]}*{w[2]})", f"({u[0]}*{w[1]}-{u[1]}*{w[0]})"]
                variables.append(f"c{f}_{i}")
                eqs.append(f"c{f}_{i}*({cr[0]}^2+{cr[1]}^2+{cr[2]}^2)-1")
        # строгая выпуклость: никакая вершина не лежит в плоскости чужой грани
        for f, face in enumerate(faces):
            p0, p1, p2 = face[0], face[1], face[2]
            A = [f"({X(p1, k)}-{X(p0, k)})" for k in range(3)]
            B = [f"({X(p2, k)}-{X(p0, k)})" for k in range(3)]
            for v in range(V):
                if v in face:
                    continue
                C = [f"({X(v, k)}-{X(p0, k)})" for k in range(3)]
                det = (f"({A[0]}*({B[1]}*{C[2]}-{B[2]}*{C[1]})-{A[1]}*({B[0]}*{C[2]}-{B[2]}*{C[0]})"
                       f"+{A[2]}*({B[0]}*{C[1]}-{B[1]}*{C[0]}))")
                variables.append(f"h{f}_{v}")
                eqs.append(f"h{f}_{v}*{det}-1")
    # раскрываем скобки и выбрасываем тождественные нули (ребро между закреплёнными вершинами)
    import sympy as sp
    syms = sp.symbols(variables)
    loc = dict(zip(variables, syms))
    expanded = []
    for e in eqs:
        poly = sp.expand(sp.sympify(e.replace("^", "**"), locals=loc))
        if poly != 0:
            expanded.append(str(poly).replace("**", "^"))
    return variables, expanded


def run_msolve(variables, eqs, threads=8, timeout=3600, tag="msolve"):
    (ROOT / "data" / "ms").mkdir(exist_ok=True)
    inp = ROOT / "data" / "ms" / f"{tag}.in"
    out = ROOT / "data" / "ms" / f"{tag}.out"
    inp.write_text(",".join(variables) + "\n0\n" + ",\n".join(eqs) + "\n")
    t0 = time.time()
    proc = subprocess.run(["msolve", "-f", str(inp), "-o", str(out), "-t", str(threads)],
                          capture_output=True, text=True, timeout=timeout)
    return proc, out.read_text() if out.exists() else "", time.time() - t0


def _classify(job):
    t, rab, timeout, threads = job
    variables, eqs = build(t, rabinowitsch=rab)
    try:
        proc, out, dt = run_msolve(variables, eqs, timeout=timeout, tag=t["id"], threads=threads)
    except subprocess.TimeoutExpired:
        return f"{t['id']}: {len(variables)} vars, {len(eqs)} eqs, TIMEOUT>{timeout}s"
    s = out.strip()
    kind = "empty" if s.startswith("[-1]") else "zero-dim" if s.startswith("[0") else "positive-dim" if s.startswith("[1") else "?"
    return f"{t['id']}: {len(variables)} vars, {len(eqs)} eqs, {dt:.2f}s  {kind}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ids", nargs="*")
    ap.add_argument("--faces", default=None)
    ap.add_argument("--timeout", type=int, default=120)
    ap.add_argument("--rab", action="store_true", help="Rabinowitsch-переменные невырожденности граней")
    ap.add_argument("--full", action="store_true", help="плюс: ни одна вершина не в плоскости чужой грани")
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()
    types = json.loads((ROOT / "data" / "types.json").read_text())
    if args.faces:
        lo, hi = map(int, args.faces.split("-"))
        types = [t for t in types if lo <= t["nF"] <= hi]
    if args.ids:
        types = [t for t in types if t["id"] in args.ids]
    from multiprocessing import Pool
    jobs = [(t, "full" if args.full else args.rab, args.timeout, args.threads) for t in types]
    with Pool(args.jobs) as pool:
        for line in pool.imap_unordered(_classify, jobs):
            print(line, flush=True)


if __name__ == "__main__":
    main()
