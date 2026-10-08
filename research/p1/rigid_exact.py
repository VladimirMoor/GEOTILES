"""Строгая обработка нульмерных типов: все вещественные решения (msolve, рациональные интервалы),
интервальная проверка строгой выпуклости и интервальный фильтр двугранных углов.

Каждое утверждение здесь — доказательство (с точностью до корректности msolve/Singular):
  * «выпуклая реализация» — все проверки выпуклости выполнены строго на интервалах;
  * «не выпуклая» — найдено нарушение, строго вне интервала погрешности;
  * «фильтр не пройден» — для некоторого класса рёбер ни одна комбинация углов (при ограничении
    Σk ≤ ⌊360°/θ_min⌋, которое верно для любого замощения) не содержит 360° или 180° в своём интервале.

Запуск: python rigid_exact.py F4-001 F5-002 ...   (нужны data/sing/<id>.basis от singular.py)
"""
import itertools
import json
import pathlib
import re
import subprocess
import sys
from fractions import Fraction

import mpmath as mp
import sympy as sp

sys.set_int_max_str_digits(0)  # в выводе msolve встречаются числа на десятки тысяч цифр

ROOT = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from exact import build  # noqa: E402
from realize import edges_of  # noqa: E402

iv = mp.iv
iv.prec = 200
mp.mp.prec = 260  # глобально; переводы чисел ниже всё равно идут через интервалы с внешним округлением


def rat_iv(lo, hi):
    """Строгий интервал, содержащий [lo, hi] (Fraction): деление в интервальной арифметике."""
    a = iv.mpf(lo.numerator) / iv.mpf(lo.denominator)
    b = iv.mpf(hi.numerator) / iv.mpf(hi.denominator)
    return iv.mpf([a.a, b.b])


def sing_name(v):
    return v.replace("_", "")


def load_msolve(text):
    """Вывод msolve как структура Python: дроби → Fraction."""
    t = text.strip().rstrip(":")
    t = re.sub(r"(-?\d+)\s*/\s*(\d+)\^(\d+)", r"F(\1,\2**\3)", t)
    t = re.sub(r"(-?\d+)\s*/\s*(\d+)", r"F(\1,\2)", t)
    return eval(t, {"F": Fraction})


def msolve_points(variables, eqs, tag, precision=256, timeout=3600, seed=7):
    """Все вещественные решения нульмерной системы (msolve) в порядке variables.

    Уравнения раскрываются sympy (парсер msolve ошибается на выражениях вида «-(-35)»).
    Запуск с -P 1: вывод [0, [char, n, deg, [имена], ...], [1, [решения]]], решения идут в порядке
    перечисленных имён (msolve может переставлять переменные). После разбора — проверка невязок."""
    import sympy as sp
    syms = sp.symbols(variables)
    loc = dict(zip(variables, syms))
    clean = []
    for e in eqs:
        p = sp.expand(sp.sympify(e.replace("^", "**"), locals=loc))
        if p != 0:
            clean.append(str(p).replace("**", "^"))
    (ROOT / "data" / "ms").mkdir(exist_ok=True)
    inp = ROOT / "data" / "ms" / f"{tag}.pts.in"
    out = ROOT / "data" / "ms" / f"{tag}.pts.out"
    inp.write_text(",".join(variables) + "\n0\n" + ",\n".join(clean) + "\n")
    subprocess.run(["msolve", "-f", str(inp), "-o", str(out), "-p", str(precision), "-t", "8", "-P", "1"],
                   check=True, capture_output=True, timeout=timeout)
    text = out.read_text()
    if text.strip().startswith("[-1]"):
        return []
    d = load_msolve(text)
    if d[0] != 0:
        raise ValueError("not zero-dimensional")
    listed = [v for v in d[1][3] if v in variables]
    if sorted(listed) != sorted(variables):
        raise RuntimeError("msolve variable list does not match input")
    raw = d[2][1] if len(d) > 2 else []
    pos = {v: i for i, v in enumerate(listed)}
    sols = [[(Fraction(sol[pos[v]][0]), Fraction(sol[pos[v]][1])) for v in variables] for sol in raw]
    check_residuals(variables, clean, sols)
    return sols


def check_residuals(variables, eqs, sols, tol=mp.mpf(10) ** -25):
    """Невязки в 260-битной арифметике; каждый многочлен нормирован на максимальный |коэффициент|."""
    import sympy as sp
    syms = [sp.Symbol(n) for n in variables]
    loc = dict(zip(variables, syms))
    fs = []
    for e in eqs:
        P = sp.Poly(sp.sympify(e.replace("^", "**"), locals=loc), *syms)
        c = max(abs(x) for x in P.coeffs())
        fs.append(sp.lambdify(syms, P.as_expr() / c, modules="mpmath"))
    for sol in sols:
        mid = [mp.mpf(lo.numerator) / lo.denominator / 2 + mp.mpf(hi.numerator) / hi.denominator / 2 for lo, hi in sol]
        scale = max(1, max(abs(m) for m in mid))
        res = max(abs(f(*mid)) for f in fs)
        if res > tol * scale ** 6:
            raise RuntimeError(f"msolve output misaligned or inaccurate: residual {mp.nstr(res, 5)}")


def msolve_real(tid, variables, eqs, precision=256, timeout=3600):
    """Все вещественные решения системы с переменными Рабиновича (exact.build(rabinowitsch=True))."""
    return msolve_points(variables, eqs, f"{tid}.rigid", precision, timeout)


def frac(tok):
    tok = tok.strip()
    m = re.fullmatch(r"(-?\d+)\s*/\s*2\^(\d+)", tok)
    if m:
        return Fraction(int(m.group(1)), 2 ** int(m.group(2)))
    m = re.fullmatch(r"(-?\d+)\s*/\s*(\d+)", tok)
    if m:
        return Fraction(int(m.group(1)), int(m.group(2)))
    return Fraction(tok)


def parse_msolve(text, nvars):
    text = text.strip()
    if text.startswith("[-1]"):
        return []
    if not text.startswith("[0"):
        raise ValueError("not zero-dimensional: " + text[:40])
    nums = re.findall(r"-?\d+\s*/\s*2\^\d+|-?\d+\s*/\s*\d+|-?\d+", text[text.index("[[["):] if "[[[" in text else "")
    vals = [frac(n) for n in nums]
    per = 2 * nvars
    sols = []
    for i in range(0, len(vals) - per + 1, per):
        chunk = vals[i:i + per]
        sols.append([(chunk[2 * k], chunk[2 * k + 1]) for k in range(nvars)])
    return sols


def vertex_intervals(t, variables, sol):
    a, b, c = t["faces"][0][:3]
    vals = {v: rat_iv(lo, hi) for v, (lo, hi) in zip(variables, sol)}
    X = []
    for v in range(t["nV"]):
        row = []
        for k in range(3):
            if v == a:
                row.append(iv.mpf(0))
            elif v == b:
                row.append(iv.mpf(1 if k == 0 else 0))
            elif v == c and k == 2:
                row.append(iv.mpf(0))
            else:
                row.append(vals[f"x{v}_{k}"])
        X.append(row)
    return X


def sub(p, q):
    return [p[i] - q[i] for i in range(3)]


def dot(p, q):
    return p[0] * q[0] + p[1] * q[1] + p[2] * q[2]


def cross(p, q):
    return [p[1] * q[2] - p[2] * q[1], p[2] * q[0] - p[0] * q[2], p[0] * q[1] - p[1] * q[0]]


def certify_convex(t, X):
    """'yes' — строго выпуклая реализация этого типа; 'no' — строго не выпуклая; 'unknown'."""
    faces = t["faces"]
    centroid = [sum(X[v][k] for v in range(len(X))) / len(X) for k in range(3)]
    normals = []
    unknown = False
    for face in faces:
        p0, p1, p2 = X[face[0]], X[face[1]], X[face[2]]
        n = cross(sub(p1, p0), sub(p2, p1))
        s = dot(n, sub(centroid, p0))
        if s.b < 0:
            pass  # нормаль наружу
        elif s.a > 0:
            n = [-x for x in n]
        else:
            unknown = True
        normals.append(n)
        for v in range(len(X)):
            if v in face:
                continue
            h = dot(n, sub(X[v], p0))
            if h.a > 0:
                return "no", normals
            if not h.b < 0:
                unknown = True
        k = len(face)
        for i in range(k):
            q0, q1, q2 = X[face[i - 1]], X[face[i]], X[face[(i + 1) % k]]
            turn = dot(cross(sub(q1, q0), sub(q2, q1)), n)
            if turn.b < 0:
                return "no", normals
            if not turn.a > 0:
                unknown = True
    return ("unknown" if unknown else "yes"), normals


def dihedral_intervals(t, X, normals):
    E = edges_of(t["faces"])
    ef = {}
    for f, face in enumerate(t["faces"]):
        for i, u in enumerate(face):
            w = face[(i + 1) % len(face)]
            ef.setdefault((min(u, w), max(u, w)), []).append(f)
    out = []
    for e in E:
        f1, f2 = ef[e]
        n1, n2 = normals[f1], normals[f2]
        cosv = dot(n1, n2) / iv.sqrt(dot(n1, n1) * dot(n2, n2))
        out.append(acos_deg(cosv, complement=True))
    return out


def acos_deg(c, complement=False):
    """Строгий интервал для arccos (в градусах): arccos убывает, считаем на концах с запасом на округление."""
    mp.mp.prec = 260
    lo, hi = max(mp.mpf(c.a), -1), min(mp.mpf(c.b), 1)
    eps = mp.mpf(2) ** -200
    a = mp.degrees(mp.acos(hi)) - eps
    b = mp.degrees(mp.acos(lo)) + eps
    if complement:
        a, b = 180 - b, 180 - a
    return iv.mpf([a, b])


def certify_filter(angles):
    """Строго: True — фильтр не пройден (тело не замощает); False — пройти может."""
    groups = []
    for a in angles:
        for idx, g in enumerate(groups):
            if not (mp.mpf(a.b) < mp.mpf(g.a) or mp.mpf(a.a) > mp.mpf(g.b)):
                groups[idx] = iv.mpf([min(mp.mpf(g.a), mp.mpf(a.a)), max(mp.mpf(g.b), mp.mpf(a.b))])
                break
        else:
            groups.append(a)
    lo_min = min(mp.mpf(g.a) for g in groups)
    K = int(mp.floor(360 / lo_min))
    los = [mp.mpf(g.a) for g in groups]
    his = [mp.mpf(g.b) for g in groups]
    n = len(groups)

    def reachable(i):
        """Есть ли мультимножество с ≥1 копией groups[i], интервал суммы которого содержит 180 или 360.
        Поиск в глубину по индексам j ≥ 0 с отсечением по нижней границе (> 360 — дальше бессмысленно)."""
        def dfs(j, lo, hi):
            if lo > 360:
                return False
            if (lo <= 180 <= hi) or (lo <= 360 <= hi):
                return True
            if j == n:
                return False
            # берём ещё одну копию j (оставаясь на j) или переходим к j+1
            return dfs(j, lo + los[j], hi + his[j]) or dfs(j + 1, lo, hi)
        return dfs(0, los[i], his[i])

    import sys as _s
    _s.setrecursionlimit(100000)
    failing = [groups[i] for i in range(n) if not reachable(i)]
    return failing, groups, K


def main(ids):
    types = {t["id"]: t for t in json.loads((ROOT / "data" / "types.json").read_text())}
    report = {}
    for tid in ids:
        t = types[tid]
        variables, eqs = build(t, rabinowitsch=True)
        sols = msolve_real(tid, variables, eqs)
        convex = []
        for sol in sols:
            X = vertex_intervals(t, variables, sol)
            verdict, normals = certify_convex(t, X)
            if verdict == "yes":
                angs = dihedral_intervals(t, X, normals)
                convex.append(angs)
            elif verdict == "unknown":
                print(f"  {tid}: a real solution could not be decided at this precision")
        # конгруэнтные копии (зеркала калибровки) дают одинаковые наборы углов
        distinct = []
        for angs in convex:
            key = sorted(float(mp.mpf(a.mid)) for a in angs)
            if not any(max(abs(x - y) for x, y in zip(key, d)) < 1e-9 for d in distinct):
                distinct.append(key)
        verdicts = []
        for angs in convex:
            key = sorted(float(mp.mpf(a.mid)) for a in angs)
            if any(v[0] == key for v in verdicts):
                continue
            failing, groups, K = certify_filter(angs)
            verdicts.append((key, [str(mp.nstr(mp.mpf(g.mid), 8)) for g in failing], K))
        report[tid] = {"realSolutions": len(sols), "convexRealizations": len(distinct),
                       "filter": [{"angles": sorted({round(x, 6) for x in k}), "failingAngles": f, "K": K} for k, f, K in verdicts]}
        print(tid, json.dumps(report[tid]))
    path = ROOT / "data" / "rigid_exact.json"
    old = json.loads(path.read_text()) if path.exists() else {}
    old.update(report)
    path.write_text(json.dumps(old, indent=1))


if __name__ == "__main__":
    main(sys.argv[1:])
