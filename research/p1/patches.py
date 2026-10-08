"""Строгое доказательство нереализуемости через «лоскуты».

Лоскут — набор граней типа (звезда вершины, грань с соседями, две соседние звезды…). Для лоскута
строится та же система, что в exact.build(rabinowitsch=True), но только на его вершинах:
длины рёбер лоскута, ромбы как параллелограммы, плоскостность граней с ≥ 5 вершинами, невырожденность
граней лоскута и рёбер, обе грани которых лежат в лоскуте. Калибровка — по первой грани лоскута.
Каждое уравнение подсистемы входит в полную систему (после переобозначения калибровки), поэтому если
подсистема не имеет решений над ℂ (msolve: [-1]), то и полная система их не имеет.

Запуск: python patches.py F8-005 F8-006 ...   или   python patches.py --from-log data/exact_rab_F8.log
"""
import argparse
import itertools
import json
import pathlib
import re
import subprocess
import sys

import sympy as sp

ROOT = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from realize import edges_of  # noqa: E402


def patch_system(t, face_ids):
    faces = [t["faces"][i] for i in face_ids]
    verts = sorted({v for f in faces for v in f})
    a, b, c = faces[0][:3]

    def X(v, k):
        if v == a:
            return "0"
        if v == b:
            return "1" if k == 0 else "0"
        if v == c and k == 2:
            return "0"
        return f"x{v}_{k}"

    variables = [X(v, k) for v in verts for k in range(3) if X(v, k) not in ("0", "1")]
    eqs = []
    E = edges_of(faces)
    for u, w in E:
        eqs.append("+".join(f"({X(u, k)}-{X(w, k)})^2" for k in range(3)) + "-1")

    def cross_norm(p0, p1, p2):
        u = [f"({X(p1, j)}-{X(p0, j)})" for j in range(3)]
        w = [f"({X(p2, j)}-{X(p0, j)})" for j in range(3)]
        cr = [f"({u[1]}*{w[2]}-{u[2]}*{w[1]})", f"({u[2]}*{w[0]}-{u[0]}*{w[2]})", f"({u[0]}*{w[1]}-{u[1]}*{w[0]})"]
        return f"({cr[0]}^2+{cr[1]}^2+{cr[2]}^2)"

    for fi, face in zip(face_ids, faces):
        if len(face) == 4:
            p, q, s, u = face
            for k in range(3):
                eqs.append(f"{X(p, k)}-{X(q, k)}+{X(s, k)}-{X(u, k)}")
        elif len(face) >= 5:
            n = [f"n{fi}_{k}" for k in range(3)]
            variables += n
            for w in face[1:]:
                eqs.append("+".join(f"{n[k]}*({X(w, k)}-{X(face[0], k)})" for k in range(3)))
            eqs.append(f"{n[0]}^2+{n[1]}^2+{n[2]}^2-1")
        if len(face) >= 4:
            variables.append(f"s{fi}")
            eqs.append(f"s{fi}*{cross_norm(*face[:3])}-1")
    # невырожденность внутренних рёбер лоскута (обе грани в лоскуте)
    ef = {}
    for face in faces:
        for i, u0 in enumerate(face):
            w0 = face[(i + 1) % len(face)]
            ef.setdefault((min(u0, w0), max(u0, w0)), []).append(face)
    for ei, ((u0, w0), fs) in enumerate(sorted(ef.items())):
        if len(fs) != 2:
            continue
        p = next(v for v in fs[0] if v not in (u0, w0))
        q = next(v for v in fs[1] if v not in (u0, w0))
        A = [f"({X(w0, k)}-{X(u0, k)})" for k in range(3)]
        B = [f"({X(p, k)}-{X(u0, k)})" for k in range(3)]
        C = [f"({X(q, k)}-{X(u0, k)})" for k in range(3)]
        det = (f"({A[0]}*({B[1]}*{C[2]}-{B[2]}*{C[1]})-{A[1]}*({B[0]}*{C[2]}-{B[2]}*{C[0]})"
               f"+{A[2]}*({B[0]}*{C[1]}-{B[1]}*{C[0]}))")
        variables.append(f"e{ei}")
        eqs.append(f"e{ei}*{det}-1")
    syms = sp.symbols(variables)
    loc = dict(zip(variables, syms))
    out = []
    for e in eqs:
        p = sp.expand(sp.sympify(e.replace("^", "**"), locals=loc))
        if p != 0:
            out.append(str(p).replace("**", "^"))
    return variables, out


def patches(t):
    """Звёзды вершин, звёзды рёбер (объединение звёзд концов), грань + соседи."""
    F = t["faces"]
    star = {v: [i for i, f in enumerate(F) if v in f] for v in range(t["nV"])}
    out = []
    for v in range(t["nV"]):
        out.append(sorted(star[v]))
    for u, w in edges_of(F):
        out.append(sorted(set(star[u]) | set(star[w])))
    for i, f in enumerate(F):
        nb = {j for j, g in enumerate(F) if j != i and len(set(f) & set(g)) >= 2}
        out.append(sorted({i} | nb))
    # без повторов, от маленьких к большим; первая грань лоскута — с наибольшим числом вершин
    uniq = sorted({tuple(p) for p in out}, key=len)
    return [sorted(p, key=lambda i: -len(F[i])) for p in uniq]


def run(t, timeout=60):
    MS = ROOT / "data" / "ms"
    MS.mkdir(exist_ok=True)
    for k, p in enumerate(patches(t)):
        variables, eqs = patch_system(t, p)
        inp, out = MS / f"{t['id']}.patch{k}.in", MS / f"{t['id']}.patch{k}.out"
        inp.write_text(",".join(variables) + "\n0\n" + ",\n".join(eqs) + "\n")
        try:
            subprocess.run(["msolve", "-f", str(inp), "-o", str(out), "-t", "2"], check=True,
                           capture_output=True, timeout=timeout)
        except subprocess.TimeoutExpired:
            continue
        if out.read_text().strip().startswith("[-1]"):
            return {"status": "empty", "patch": p, "faces": [t["faces"][i] for i in p]}
    return {"status": "no empty patch"}


def _job(args):
    t, timeout = args
    return t["id"], run(t, timeout)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ids", nargs="*")
    ap.add_argument("--from-log", nargs="*", default=None)
    ap.add_argument("--timeout", type=int, default=60)
    ap.add_argument("--jobs", type=int, default=6)
    args = ap.parse_args()
    T = {t["id"]: t for t in json.loads((ROOT / "data" / "types.json").read_text())}
    ids = list(args.ids)
    if args.from_log:
        for log in args.from_log:
            for line in open(log):
                m = re.match(r"(F\d-\d+):.*TIMEOUT", line)
                if m:
                    ids.append(m.group(1))
    from multiprocessing import Pool
    res = {}
    with Pool(args.jobs) as pool:
        for tid, r in pool.imap_unordered(_job, [(T[i], args.timeout) for i in ids]):
            res[tid] = r
            print(tid, r["status"], r.get("faces", ""), flush=True)
    path = ROOT / "data" / "patches.json"
    old = json.loads(path.read_text()) if path.exists() else {}
    old.update(res)
    path.write_text(json.dumps(old, indent=1))


if __name__ == "__main__":
    main()
