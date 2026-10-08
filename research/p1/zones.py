"""Формулировка через векторы классов рёбер (зоны).

Ромб — параллелограмм, поэтому его противоположные рёбра — один и тот же вектор. Склеиваем рёбра
в классы (с учётом ориентации): вектор ребра (u→v, u < v) равен σ·w_z, где w_z — единичный вектор класса z.
Тогда положения вершин линейны по w (сумма вдоль остовного дерева), ромбы выполнены автоматически, и
неизвестных становится 3·(число классов) вместо 3·(число вершин).

Уравнения:
  |w_z|² = 1;
  замыкание каждой неромбической грани: Σ σ_e w_{z(e)} = 0 по её рёбрам в порядке обхода
  (границы граней порождают пространство циклов, поэтому этого достаточно для согласованности);
  плоскостность граней с ≥ 5 вершинами: n_f · w_{z(e)} = 0 для её рёбер, |n_f|² = 1;
  калибровка: два первых ребра первой грани — w = (1,0,0) и w = (p,q,0);
  Рабинович: |w_a × w_b|² ≠ 0 для соседних рёбер каждой грани с ≥ 4 вершинами;
             det(рёбра у вершины по разные стороны ребра) ≠ 0 — соседние грани не копланарны.
Если два ребра одного треугольника попали в один класс, система сразу противоречива (это правильно:
у треугольника нет параллельных сторон).

Запуск: python zones.py F7-002 F7-027 ... [--timeout 600 --jobs 3]  (те же исходы, что у exact.py)
"""
import argparse
import json
import pathlib
import subprocess
import sys
import time

import sympy as sp

ROOT = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from realize import edges_of  # noqa: E402


def zone_classes(t):
    """Классы рёбер с ориентацией: edge(u<v) → (класс, знак)."""
    E = edges_of(t["faces"])
    idx = {e: i for i, e in enumerate(E)}
    parent = list(range(len(E)))
    sign = [1] * len(E)  # знак относительно родителя

    def find(x):
        if parent[x] == x:
            return x, 1
        r, s = find(parent[x])
        parent[x], sign[x] = r, sign[x] * s
        return r, sign[x]

    def union(e1, e2, rel):  # vec(e1) = rel · vec(e2)
        (r1, s1), (r2, s2) = find(e1), find(e2)
        if r1 == r2:
            return s1 == rel * s2  # согласованность
        parent[r1] = r2
        sign[r1] = rel * s2 * s1
        return True

    def vec_sign(u, v):  # u→v относительно канонического (min,max)
        return 1 if u < v else -1

    consistent = True
    for face in t["faces"]:
        if len(face) == 4:
            a, b, c, d = face
            # x_b − x_a = x_c − x_d,  x_c − x_b = x_d − x_a
            for (p, q), (r, s) in (((a, b), (d, c)), ((b, c), (a, d))):
                e1, e2 = idx[(min(p, q), max(p, q))], idx[(min(r, s), max(r, s))]
                rel = vec_sign(p, q) * vec_sign(r, s)
                consistent &= union(e1, e2, rel)
    cls = {}
    roots = {}
    for e, i in idx.items():
        r, s = find(i)
        roots.setdefault(r, len(roots))
        cls[e] = (roots[r], s)
    return cls, len(roots), consistent


def build_zones(t):
    cls, Z, consistent = zone_classes(t)
    faces = t["faces"]
    f0 = faces[0]
    e0 = (min(f0[0], f0[1]), max(f0[0], f0[1]))
    e1 = (min(f0[1], f0[2]), max(f0[1], f0[2]))
    z0, z1 = cls[e0][0], cls[e1][0]

    def W(z, k):
        if z == z0:
            return "1" if k == 0 else "0"
        if z == z1 and k == 2:
            return "0"
        return f"w{z}_{k}"

    variables = [W(z, k) for z in range(Z) for k in range(3) if W(z, k) not in ("0", "1")]
    eqs = []
    for z in range(Z):
        if z != z0:
            eqs.append("+".join(f"({W(z, k)})^2" for k in range(3)) + "-1")

    def edge_vec(u, v):
        """Вектор u→v как список строк-выражений."""
        z, s = cls[(min(u, v), max(u, v))]
        s *= 1 if u < v else -1
        return [f"({s})*({W(z, k)})" for k in range(3)]

    for f, face in enumerate(faces):
        k = len(face)
        if k == 4:
            continue
        vecs = [edge_vec(face[i], face[(i + 1) % k]) for i in range(k)]
        for c in range(3):
            eqs.append("+".join(v[c] for v in vecs))  # замыкание
        if k >= 5:
            n = [f"n{f}_{c}" for c in range(3)]
            variables += n
            for v in vecs:
                eqs.append("+".join(f"{n[c]}*{v[c]}" for c in range(3)))
            eqs.append(f"{n[0]}^2+{n[1]}^2+{n[2]}^2-1")

    def cross_sq(a, b):
        cr = [f"({a[1]}*{b[2]}-{a[2]}*{b[1]})", f"({a[2]}*{b[0]}-{a[0]}*{b[2]})", f"({a[0]}*{b[1]}-{a[1]}*{b[0]})"]
        return f"({cr[0]}^2+{cr[1]}^2+{cr[2]}^2)"

    def det(a, b, c):
        return (f"(({a[0]})*(({b[1]})*({c[2]})-({b[2]})*({c[1]}))-({a[1]})*(({b[0]})*({c[2]})-({b[2]})*({c[0]}))"
                f"+({a[2]})*(({b[0]})*({c[1]})-({b[1]})*({c[0]})))")

    for f, face in enumerate(faces):
        k = len(face)
        if k >= 4:
            for i in range(1):  # минимально: одна пара соседних рёбер на грань
                a = edge_vec(face[i], face[(i + 1) % k])
                b = edge_vec(face[(i + 1) % k], face[(i + 2) % k])
                variables.append(f"s{f}_{i}")
                eqs.append(f"s{f}_{i}*{cross_sq(a, b)}-1")
    # соседние грани не копланарны: у ребра (u,w) берём следующее ребро каждой грани из вершины u
    ef = {}
    for face in faces:
        k = len(face)
        for i, u in enumerate(face):
            w = face[(i + 1) % k]
            ef.setdefault((min(u, w), max(u, w)), []).append((face, i))
    for ei, ((u, w), fl) in enumerate(sorted(ef.items())):
        (fa, ia), (fb, ib) = fl
        ka, kb = len(fa), len(fb)
        # в грани fa ребро идёт fa[ia]→fa[ia+1]; третья точка — fa[ia+2]; аналогично fb
        p = edge_vec(fa[(ia + 1) % ka], fa[(ia + 2) % ka])
        q = edge_vec(fb[(ib + 1) % kb], fb[(ib + 2) % kb])
        e = edge_vec(fa[ia], fa[(ia + 1) % ka])
        variables.append(f"r{ei}")
        eqs.append(f"r{ei}*{det(e, p, q)}-1")
    syms = sp.symbols(variables)
    loc = dict(zip(variables, syms))
    out = []
    for e in eqs:
        poly = sp.expand(sp.sympify(e.replace("^", "**"), locals=loc))
        if poly != 0:
            out.append(str(poly).replace("**", "^"))
    if not consistent:
        out.append("1")  # противоречивые ориентации ромбов: решений нет
    return variables, out, Z


def classify(job):
    t, timeout, threads = job
    variables, eqs, Z = build_zones(t)
    MS = ROOT / "data" / "ms"
    MS.mkdir(exist_ok=True)
    inp, out = MS / f"{t['id']}.zones.in", MS / f"{t['id']}.zones.out"
    inp.write_text(",".join(variables) + "\n0\n" + ",\n".join(eqs) + "\n")
    t0 = time.time()
    try:
        subprocess.run(["msolve", "-f", str(inp), "-o", str(out), "-t", str(threads)], check=True,
                       capture_output=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return f"{t['id']}: zones={Z} {len(variables)} vars, {len(eqs)} eqs, TIMEOUT>{timeout}s"
    s = out.read_text().strip()
    kind = "empty" if s.startswith("[-1]") else "zero-dim" if s.startswith("[0") else "positive-dim" if s.startswith("[1") else "?"
    return f"{t['id']}: zones={Z} {len(variables)} vars, {len(eqs)} eqs, {time.time() - t0:.2f}s  {kind}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ids", nargs="*")
    ap.add_argument("--timeout", type=int, default=600)
    ap.add_argument("--jobs", type=int, default=3)
    ap.add_argument("--threads", type=int, default=3)
    args = ap.parse_args()
    T = {t["id"]: t for t in json.loads((ROOT / "data" / "types.json").read_text())}
    from multiprocessing import Pool
    with Pool(args.jobs) as pool:
        for line in pool.imap_unordered(classify, [(T[i], args.timeout, args.threads) for i in args.ids]):
            print(line, flush=True)


if __name__ == "__main__":
    main()
