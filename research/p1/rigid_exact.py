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

ROOT = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from exact import build  # noqa: E402
from realize import edges_of  # noqa: E402

iv = mp.iv
iv.prec = 200


def sing_name(v):
    return v.replace("_", "")


def msolve_real(tid, variables, eqs, precision=256, timeout=3600):
    """Все вещественные решения системы с переменными Рабиновича (exact.build(rabinowitsch=True))."""
    (ROOT / "data" / "ms").mkdir(exist_ok=True)
    inp = ROOT / "data" / "ms" / f"{tid}.rigid.in"
    out = ROOT / "data" / "ms" / f"{tid}.rigid.out"
    inp.write_text(",".join(variables) + "\n0\n" + ",\n".join(eqs) + "\n")
    subprocess.run(["msolve", "-f", str(inp), "-o", str(out), "-p", str(precision), "-t", "8"], check=True,
                   capture_output=True, timeout=timeout)
    return parse_msolve(out.read_text(), len(variables))


def frac(tok):
    tok = tok.strip()
    m = re.fullmatch(r"(-?\d+)\s*/\s*2\^(\d+)", tok)
    if m:
        return Fraction(int(m.group(1)), 2 ** int(m.group(2)))
    return Fraction(tok)


def parse_msolve(text, nvars):
    text = text.strip()
    if text.startswith("[-1]"):
        return []
    if not text.startswith("[0"):
        raise ValueError("not zero-dimensional: " + text[:40])
    nums = re.findall(r"-?\d+\s*/\s*2\^\d+|-?\d+", text[text.index("[[["):] if "[[[" in text else "")
    vals = [frac(n) for n in nums]
    per = 2 * nvars
    sols = []
    for i in range(0, len(vals) - per + 1, per):
        chunk = vals[i:i + per]
        sols.append([(chunk[2 * k], chunk[2 * k + 1]) for k in range(nvars)])
    return sols


def vertex_intervals(t, variables, sol):
    a, b, c = t["faces"][0][:3]
    vals = {v: iv.mpf([float(lo) if False else mp.mpf(lo.numerator) / lo.denominator, mp.mpf(hi.numerator) / hi.denominator])
            for v, (lo, hi) in zip(variables, sol)}
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
    failing = []
    for i, gi in enumerate(groups):
        ok = False
        for total in range(1, K + 1):
            for combo in itertools.combinations_with_replacement(range(len(groups)), total):
                if i not in combo:
                    continue
                s = sum((groups[j] for j in combo), iv.mpf(0))
                if (s.a <= 360 <= s.b) or (s.a <= 180 <= s.b):
                    ok = True
                    break
            if ok:
                break
        if not ok:
            failing.append(gi)
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
