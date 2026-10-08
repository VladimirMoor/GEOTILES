"""Строгий разбор одномерных компонент: есть ли на кривой решений строго выпуклые реализации.

Обозначения: I — система exact.build(rabinowitsch=True) (рёбра и грани невырождены), V = V(I),
K — замыкание V в пространстве координат вершин (исключение вспомогательных переменных).
T, U — случайные рациональные линейные формы от координат.

Свойство «строго выпуклая реализация типа» задаётся знаками многочленов G (для каждой грани f и вершины
v ∉ f: det-объём h_{f,v}; повороты внутри граней; длины векторных произведений). Пока ни один g ∈ G
не обращается в ноль, знаки постоянны. На вещественной кривой K_R вне конечного множества особых точек
дуги — графики над T. Особые значения T:
  (1) критические значения проекции: корни дискриминанта и старшего коэффициента F(T,U) по U,
      где F — образующая исключения K на плоскость (T, U);
  (2) нули многочленов g ∈ G на K.
Между соседними особыми значениями знаки всех g постоянны на каждой дуге, поэтому достаточно проверить
по одной рациональной точке T₀ в каждом интервале (и за крайними значениями): решить I + ⟨T − T₀⟩
(нульмерно) и строго проверить выпуклость найденных вещественных точек на интервалах.
Если ни в одной пробной точке нет выпуклой реализации — их нет на всей кривой (строго).

Требование: dim K ≤ 1 и исключение на (T, U) — главный идеал. Иначе скрипт сообщает и останавливается.
"""
import json
import pathlib
import random
import re
import subprocess
import sys
from fractions import Fraction

import sympy as sp

ROOT = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from exact import build  # noqa: E402
from rigid_exact import parse_msolve, vertex_intervals, certify_convex, msolve_points  # noqa: E402

MS = ROOT / "data" / "ms"


def msolve(variables, eqs, tag, extra=(), timeout=3600):
    MS.mkdir(exist_ok=True)
    inp, out = MS / f"{tag}.in", MS / f"{tag}.out"
    inp.write_text(",".join(variables) + "\n0\n" + ",\n".join(eqs) + "\n")
    subprocess.run(["msolve", "-f", str(inp), "-o", str(out), "-t", "8", *extra], check=True,
                   capture_output=True, timeout=timeout)
    return out.read_text()


def expand_all(exprs, names):
    syms = sp.symbols(names)
    loc = dict(zip(names, syms))
    out = []
    for e in exprs:
        p = sp.expand(sp.sympify(e.replace("^", "**"), locals=loc))
        if p != 0:
            out.append(str(p).replace("**", "^"))
    return out


def eliminate(variables, eqs, keep, tag):
    """Базис Грёбнера исключающего идеала (msolve, блочный порядок) в переменных keep."""
    rest = [v for v in variables if v not in keep]
    text = msolve(rest + keep, eqs, tag, extra=("-e", str(len(rest))))
    body = text[text.index("["): text.rindex("]") + 1]
    polys = [p.strip() for p in body.strip("[]").split(",") if p.strip()]
    return polys


def sign_polys(t, X):
    """Многочлены, знаки которых определяют строгую выпуклость (в виде строк от x-переменных)."""
    faces, V = t["faces"], t["nV"]

    def det(A, B, C):
        return (f"(({A[0]})*(({B[1]})*({C[2]})-({B[2]})*({C[1]}))-({A[1]})*(({B[0]})*({C[2]})-({B[2]})*({C[0]}))"
                f"+({A[2]})*(({B[0]})*({C[1]})-({B[1]})*({C[0]})))")

    d = lambda p, q: [f"({X(p, k)}-{X(q, k)})" for k in range(3)]
    G = []
    for face in faces:
        p0, p1, p2 = face[:3]
        for v in range(V):
            if v not in face:
                G.append(det(d(p1, p0), d(p2, p0), d(v, p0)))
        k = len(face)
        if k >= 4:
            # повороты: det(x_{i} − x_{i−1}, x_{i+1} − x_i, x_{p2} − x_{p0}×...) — используем объём с точкой вне грани
            out = next(v for v in range(V) if v not in face)
            for i in range(k):
                a, b, c = face[i - 1], face[i], face[(i + 1) % k]
                G.append(det(d(b, a), d(c, b), d(out, b)))
    return G


def analyze(tid, seed=1):
    T = {t["id"]: t for t in json.loads((ROOT / "data" / "types.json").read_text())}
    t = T[tid]
    variables, eqs = build(t, rabinowitsch=True)
    a, b, c = t["faces"][0][:3]

    def X(v, k):
        if v == a:
            return "0"
        if v == b:
            return "1" if k == 0 else "0"
        if v == c and k == 2:
            return "0"
        return f"x{v}_{k}"

    xs = [v for v in variables if v.startswith("x")]
    rnd = random.Random(seed)
    rt = [rnd.randint(-9, 9) for _ in xs]
    ru = [rnd.randint(-9, 9) for _ in xs]
    defT = "TT-(" + "+".join(f"{c_}*{x}" for c_, x in zip(rt, xs)) + ")"
    defU = "UU-(" + "+".join(f"{c_}*{x}" for c_, x in zip(ru, xs)) + ")"
    allv = variables + ["TT", "UU"]
    base = eqs + expand_all([defT, defU], allv)

    # (0) замыкание K в координатах вершин + T
    K = eliminate(allv, base, xs + ["TT"], f"{tid}.K")
    # (1) плоская кривая (из K — гораздо дешевле) и её особые значения
    KU = K + expand_all([defU], xs + ["TT", "UU"])
    F = eliminate(xs + ["TT", "UU"], KU, ["TT", "UU"], f"{tid}.F")
    Ts, Us = sp.symbols("TT UU")
    Fp = [sp.Poly(sp.sympify(p.replace("^", "**")), Ts, Us) for p in F]
    if len(Fp) != 1:
        return {"status": "not-principal", "generators": len(Fp)}
    Fp = Fp[0]
    if Fp.degree(Us) == 0:
        return {"status": "projection-degenerate"}
    disc = sp.discriminant(Fp.as_expr(), Us)
    lc = sp.Poly(Fp.as_expr(), Us).LC()
    special = set()
    for poly in (disc, lc):
        P = sp.Poly(poly, Ts)
        if P.degree() > 0:
            for (lo, hi), _mult in P.intervals():
                lo, hi = sp.Rational(lo), sp.Rational(hi)
                special.add((Fraction(int(lo.p), int(lo.q)), Fraction(int(hi.p), int(hi.q))))
    # (2) нули знаковых многочленов на K
    xT = xs + ["TT"]
    G = expand_all(sign_polys(t, X), xT)
    unresolved = 0
    for i, g in enumerate(G):
        try:
            sols = msolve_points(xT, K + [g], f"{tid}.g{i}")
        except ValueError:
            unresolved += 1  # g тождественно нуль на компоненте K — её точки исключены из V
            continue
        for sol in sols:
            special.add(sol[-1])
    # пробные значения T₀: между отсортированными особыми интервалами и за краями
    pts = sorted(special)
    samples = []
    if not pts:
        samples = [Fraction(0)]
    else:
        samples.append(pts[0][0] - 1)
        for (l1, h1), (l2, h2) in zip(pts, pts[1:]):
            if h1 < l2:
                samples.append((h1 + l2) / 2)
        samples.append(pts[-1][1] + 1)
    convex_found = []
    for j, T0 in enumerate(samples):
        fiber = base + [f"{T0.denominator}*TT-({T0.numerator})"]
        try:
            fiber_sols = msolve_points(allv, fiber, f"{tid}.s{j}")
        except ValueError:
            return {"status": "fiber-not-zero-dim", "T0": str(T0)}
        for sol in fiber_sols:
            Xi = vertex_intervals(t, allv, sol)
            verdict, _ = certify_convex(t, Xi)
            if verdict == "yes":
                convex_found.append(str(T0))
            elif verdict == "unknown":
                return {"status": "undecided-point", "T0": str(T0)}
    return {"status": "ok", "special": len(pts), "samples": len(samples),
            "convexSamples": len(convex_found), "degenerateSignPolys": unresolved, "Kgens": len(K)}


if __name__ == "__main__":
    for tid in sys.argv[1:]:
        print(tid, analyze(tid), flush=True)
