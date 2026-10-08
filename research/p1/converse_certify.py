"""Строгая (интервальная) проверка: конкретные суммы отрезков и треугольников размерности ≥ 5
не проходят фильтр двугранных углов, а значит, не замощают ℝ³.

Повороты — рациональные (преобразование Кэли от целочисленной кососимметрической матрицы), треугольник —
правильный единичный (√3 в интервальной арифметике). Комбинаторика оболочки берётся из float-вычисления,
а все углы и проверка фильтра — в интервалах (rigid_exact.certify_filter, acos_deg).
Отказ фильтра на ребре строго доказывает, что тело не замощает пространство.
"""
import itertools
import json
import pathlib
import sys
from fractions import Fraction

import mpmath as mp
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from rigid_exact import iv, acos_deg, certify_filter  # noqa: E402
from converse import labels_and_faces  # noqa: E402


def cayley(x, y, z):
    """Рациональная матрица поворота (I − A)⁻¹(I + A), A — кососимметрическая (x, y, z)."""
    A = np.array([[0, -z, y], [z, 0, -x], [-y, x, 0]], dtype=object)
    I = np.eye(3, dtype=object)
    import sympy as sp
    M = sp.Matrix((I - A).tolist()).inv() * sp.Matrix((I + A).tolist())
    return [[Fraction(int(sp.Rational(M[i, j]).p), int(sp.Rational(M[i, j]).q)) for j in range(3)] for i in range(3)]


def summand(name):
    s3 = iv.sqrt(3)
    if name == "seg":
        return [[iv.mpf(0)] * 3, [iv.mpf(1), iv.mpf(0), iv.mpf(0)]]
    return [[iv.mpf(0)] * 3, [iv.mpf(1), iv.mpf(0), iv.mpf(0)], [iv.mpf(1) / 2, s3 / 2, iv.mpf(0)]]


def apply(R, v):
    return [sum((iv.mpf(R[i][j].numerator) / R[i][j].denominator) * v[j] for j in range(3)) for i in range(3)]


def certify(kind, rng):
    names = kind.split("+")
    parts = []
    for i, nm in enumerate(names):
        V = summand(nm)
        if i > 0:
            R = cayley(*[int(x) for x in rng.integers(-4, 5, size=3)])
            V = [apply(R, v) for v in V]
        parts.append(V)
    # комбинаторика по серединам интервалов
    mid = [[mp.matrix([float(mp.mpf(x.mid)) for x in v]) for v in P] for P in parts]
    _, _, faces = labels_and_faces(mid)
    pt = lambda lab: [sum((parts[k][lab[k]][j] for k in range(len(parts))), iv.mpf(0)) for j in range(3)]
    allv = {lab for F in faces for lab in F}
    cen = [sum((pt(l)[j] for l in allv), iv.mpf(0)) / len(allv) for j in range(3)]
    normals = []
    for F in faces:
        p0, p1, p2 = pt(F[0]), pt(F[1]), pt(F[2])
        u = [p1[j] - p0[j] for j in range(3)]
        w = [p2[j] - p0[j] for j in range(3)]
        n = [u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0]]
        sgn = sum(n[j] * (cen[j] - p0[j]) for j in range(3))
        if sgn.a > 0:
            n = [-x for x in n]
        elif not sgn.b < 0:
            return {"status": "orientation undecided"}
        normals.append(n)
    edges = {}
    for fi, F in enumerate(faces):
        for i in range(len(F)):
            edges.setdefault(frozenset((F[i], F[(i + 1) % len(F)])), []).append(fi)
    angs = []
    for e, (f1, f2) in edges.items():
        n1, n2 = normals[f1], normals[f2]
        cosv = sum(n1[j] * n2[j] for j in range(3)) / iv.sqrt(sum(x * x for x in n1) * sum(x * x for x in n2))
        angs.append(acos_deg(cosv, complement=True))
    failing, groups, K = certify_filter(angs)
    return {"status": "does not tile (certified)" if failing else "filter not excluded",
            "faces": len(faces), "failingAngleClasses": len(failing), "K": K}


if __name__ == "__main__":
    rng = np.random.default_rng(31)
    out = {}
    for kind in ["tri+tri+seg", "tri+seg+seg+seg", "seg+seg+seg+seg+seg", "tri+tri+tri", "tri+tri+seg+seg"]:
        res = [certify(kind, rng) for _ in range(5)]
        out[kind] = res
        print(kind, [r["status"] + f" ({r.get('failingAngleClasses')} classes, F={r.get('faces')})" for r in res], flush=True)
    (ROOT / "data" / "converse_certified.json").write_text(json.dumps(out, indent=1))
