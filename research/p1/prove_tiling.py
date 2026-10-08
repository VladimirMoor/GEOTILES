"""Доказательство (компьютерное, точное): суммы Минковского T₁ ⊕ T₂ (два треугольника) и Par ⊕ T
(параллелограмм и треугольник) общего положения в ℝ³ замощают пространство сдвигами и центральной
симметрией.

1. Аффинная нормализация. Разбиение сдвигами и центральными симметриями переходит в разбиение того же
   вида при любом аффинном преобразовании. Поэтому считаем
      T₁ ⊕ T₂:  T₁ = {0, e₁, e₂},  T₂ = {0, e₃, v},   v = (a, b, c);
      Par ⊕ T:  Par = {0, e₁, e₁+e₂, e₂},  T = {0, e₃, v}.
   Это все пары общего положения с точностью до аффинной эквивалентности.
2. Камеры. Общее положение — все определители троек векторов рёбер ненулевые. Это линейные формы
   a, b, a+b (только для T₁ ⊕ T₂), c, c−1. Камера — сектор плоскости (a, b) × интервал по c;
   её замыкание задаётся явными вершинами и лучами.
3. Внутри камеры:
   (i)  грани P (наборы меток вершин) — опорные: n_F·(x − x_F) = 0 тождественно на F и < 0 вне F;
   (ii) спаривания граней: для каждой грани F есть γ (сдвиг на вектор решётки, возможно с центральной
        симметрией) с γ(G) = F для грани G (тождества) и γP по другую сторону F (неравенства);
   (iii) циклы рёбер: обход вокруг каждого ребра по спариваниям возвращается к P с тождественным
        произведением; сумма двугранных углов равна 2π в пробной точке, а так как она непрерывна и
        кратна 2π, то и во всей камере.
   Все формы в (i)–(ii) линейны по (a, b, c) (вектор v входит в определитель не более одного раза),
   поэтому их знак на камере проверяется точно — в вершинах и на лучах замыкания.
4. По теореме Пуанкаре о многограннике (евклидов случай, компактная выпуклая область, спаривания
   граней изометриями, условие циклов рёбер) образы P под группой спариваний замощают ℝ³.

Запуск: python prove_tiling.py [tri+tri | par+tri]
"""
import itertools
import math
import sys
from fractions import Fraction

import numpy as np
import sympy as sp

sys.path.insert(0, __import__("pathlib").Path(__file__).resolve().parent.as_posix())
from tilecheck import Tile, coverage, crystal  # noqa: E402

a, b, c = sp.symbols("a b c")
E1, E2, E3 = sp.Matrix([1, 0, 0]), sp.Matrix([0, 1, 0]), sp.Matrix([0, 0, 1])
VV = sp.Matrix([a, b, c])
Z = sp.zeros(3, 1)


def summands(kind):
    if kind == "tri+tri":
        return [[Z, E1, E2], [Z, E3, VV]]
    return [[Z, E1, E1 + E2, E2], [Z, E3, VV]]


def chambers(kind):
    """Камеры: (лучи сектора в (a,b), интервал c, пробная точка)."""
    if kind == "tri+tri":
        # прямые a=0, b=0, a+b=0 — 6 секторов; лучи по часовой: (1,0),(1,-1),(0,-1),(-1,0),(-1,1),(0,1)
        rays = [(1, 0), (1, -1), (0, -1), (-1, 0), (-1, 1), (0, 1)]
    else:
        rays = [(1, 0), (0, -1), (-1, 0), (0, 1)]
    sectors = [(rays[i], rays[(i + 1) % len(rays)]) for i in range(len(rays))]
    intervals = [(None, 0), (0, 1), (1, None)]
    out = []
    for (r1, r2), (lo, hi) in itertools.product(sectors, intervals):
        mid_ab = (Fraction(r1[0] + r2[0]) * Fraction(7, 10) + Fraction(r1[0]) * Fraction(1, 10),
                  Fraction(r1[1] + r2[1]) * Fraction(7, 10) + Fraction(r1[1]) * Fraction(1, 10))
        cc = Fraction(-1, 2) if lo is None else Fraction(3, 2) if hi is None else Fraction(2, 5)
        out.append({"rays": (r1, r2), "c": (lo, hi), "sample": (mid_ab[0], mid_ab[1], cc)})
    return out


def sign_on_chamber(f, ch):
    """f — аффинная форма от a,b,c. Возвращает '+' (≥0 на замыкании, >0 в пробной точке), '-', '0' или '?'."""
    f = sp.expand(f)
    if f == 0:
        return "0"
    poly = sp.Poly(f, a, b, c)
    if poly.total_degree() > 1:
        return "?nonlinear"
    f0 = poly.coeff_monomial(1)
    fa, fb, fc = poly.coeff_monomial(a), poly.coeff_monomial(b), poly.coeff_monomial(c)
    lo, hi = ch["c"]
    verts = [v for v in (lo, hi) if v is not None]
    vals = [f0 + fc * v for v in verts]                         # вершины (0, 0, c)
    rays = [fa * r[0] + fb * r[1] for r in ch["rays"]]          # лучи сектора (по c = 0)
    if lo is None:
        rays.append(-fc)
    if hi is None:
        rays.append(fc)
    s = ch["sample"]
    at = f.subs({a: s[0], b: s[1], c: s[2]})
    if all(x >= 0 for x in vals) and all(r >= 0 for r in rays) and at > 0:
        return "+"
    if all(x <= 0 for x in vals) and all(r <= 0 for r in rays) and at < 0:
        return "-"
    return "?"


def num(vec, s):
    return np.array([float(x.subs({a: s[0], b: s[1], c: s[2]})) for x in vec], float)


def build_P(kind, s):
    parts = summands(kind)
    labels = list(itertools.product(range(len(parts[0])), range(len(parts[1]))))
    sym = {lab: parts[0][lab[0]] + parts[1][lab[1]] for lab in labels}
    pts = np.array([num(sym[lab], s) for lab in labels])
    T = Tile(pts)
    vlab = []
    for v in T.vertices:
        k = int(np.argmin(np.linalg.norm(pts - v, axis=1)))
        vlab.append(labels[k])
    faces = []
    for f in T.faces():
        Q = T.vertices[f]
        cen = Q.mean(axis=0)
        nrm = np.linalg.svd(Q - cen)[2][-1]
        u = (Q[0] - cen) / np.linalg.norm(Q[0] - cen)
        wv = np.cross(nrm, u)
        order = np.argsort(np.arctan2((Q - cen) @ wv, (Q - cen) @ u))  # циклический обход грани
        faces.append([vlab[f[i]] for i in order])
    return sym, vlab, faces, T


def tri_labels(V):
    for i in range(3):
        rest = [j for j in range(3) if j != i]
        for j1, j2 in (rest, rest[::-1]):
            yield V[i], V[j1] - V[i], V[j2] - V[i]


def candidates(kind):
    parts = summands(kind)
    out = []
    if kind == "tri+tri":
        for (o1, u1, u2), (o2, v1, v2) in itertools.product(tri_labels(parts[0]), tri_labels(parts[1])):
            cc = 2 * (o1 + o2)
            out.append(([v1, u1, u2 - v2], cc + v1 + u1 + u2))
            out.append(([u1 + v1, u1 + v2, u2 - u1], cc))
    else:
        R = parts[0]
        for i in range(4):
            o = R[i]
            nb = [R[(i + 1) % 4] - o, R[(i + 3) % 4] - o]
            for r1, r2 in (nb, nb[::-1]):
                for (ot, t1, t2) in tri_labels(parts[1]):
                    out.append(([r1 + r2 + t1, r1 + r2 + t2, r1 - r2 + t1 - t2], 2 * (o + ot) + r1 + t1))
    return out


def works_numeric(T, L, w):
    I = np.eye(3)
    if abs(abs(np.linalg.det(np.array(L))) - 2 * T.volume) > 1e-9:
        return False
    R = 3 * T.radius + 2
    cov = coverage(T, crystal(L, [(I, np.zeros(3)), (-I, w)], 1, R=R), samples=400, R=R)
    return cov["exactlyOne"] == cov["samples"]


class G:
    """Элемент группы: x ↦ eps·x + t (t — символьный вектор)."""
    def __init__(self, eps, t):
        self.eps, self.t = eps, sp.Matrix(t)

    def __call__(self, x):
        return self.eps * x + self.t

    def compose(self, other):  # self ∘ other
        return G(self.eps * other.eps, self.eps * other.t + self.t)

    def is_identity(self):
        return self.eps == 1 and all(sp.expand(x) == 0 for x in self.t)


def prove_chamber(kind, ch):
    s = ch["sample"]
    sym, vlab, faces, T = build_P(kind, s)
    # разбиение в пробной точке
    choice = None
    for L, w in candidates(kind):
        Ln = [num(x, s) for x in L]
        if works_numeric(T, Ln, num(w, s)):
            choice = (L, w)
            break
    if choice is None:
        return {"ok": False, "why": "no tiling formula at sample"}
    L, w = choice
    # (i) опорность граней
    normals = {}
    centroid = np.mean([num(sym[x], s) for x in vlab], axis=0)
    for F in faces:
        p0, p1, p2 = (sym[x] for x in F[:3])
        n = (p1 - p0).cross(p2 - p0)
        if num(n, s) @ (centroid - num(p0, s)) > 0:  # нормаль наружу
            n = -n
        normals[tuple(F)] = (n, p0)
        for lab in vlab:
            h = (n.T * (sym[lab] - p0))[0]
            sg = sign_on_chamber(h, ch)
            if lab in F and sg != "0":
                return {"ok": False, "why": f"face {F}: vertex {lab} not on plane identically"}
            if lab not in F and sg != "-":
                return {"ok": False, "why": f"face {F}: vertex {lab} side sign {sg}"}
    # (ii) спаривания граней: ищем γ среди сдвигов на n·L (|n|≤2) с/без центральной симметрии
    group = []
    for eps in (1, -1):
        for nvec in itertools.product(range(-2, 3), repeat=3):
            t = sum((k * v for k, v in zip(nvec, L)), Z)
            if eps == -1:
                t = t + w
            group.append(G(eps, t))
    pts_num = {lab: num(sym[lab], s) for lab in vlab}
    pairing = {}
    for F in faces:
        Fset = {tuple(np.round(pts_num[x], 9)) for x in F}
        found = None
        for g in group:
            if g.is_identity():
                continue
            gn = {lab: num(g(sym[lab]), s) for lab in vlab}
            for Gf in faces:
                if len(Gf) != len(F):
                    continue
                if {tuple(np.round(gn[x], 9)) for x in Gf} == Fset:
                    found = (g, Gf)
                    break
            if found:
                break
        if not found:
            return {"ok": False, "why": f"no face pairing for {F}"}
        g, Gf = found
        mp = {}
        for q in Gf:
            p = next(x for x in F if np.allclose(num(g(sym[q]), s), pts_num[x], atol=1e-9))
            if any(sp.expand(e) != 0 for e in (g(sym[q]) - sym[p])):
                return {"ok": False, "why": f"pairing of {F} not an identity"}
            mp[q] = p
        n, p0 = normals[tuple(F)]
        for lab in vlab:
            if lab in mp:
                continue
            h = (n.T * (g(sym[lab]) - p0))[0]
            if sign_on_chamber(h, ch) != "+":
                return {"ok": False, "why": f"neighbor across {F} not strictly outside"}
        pairing[tuple(F)] = (g, tuple(Gf), mp)
    # (ii') согласованность спариваний: спаривание грани G — ровно обратный элемент
    for F, (g, Gf, mp) in pairing.items():
        g2, F2, _ = pairing[Gf]
        if F2 != F or not g2.compose(g).is_identity():
            return {"ok": False, "why": f"pairings of {F} and {Gf} are not mutually inverse"}
    # (iii) циклы рёбер
    edges = set()
    for F in faces:
        for i in range(len(F)):
            edges.add(frozenset((F[i], F[(i + 1) % len(F)])))
    cycles = []
    for e in edges:
        F0 = next(F for F in faces if e <= set(F))
        g_total = G(1, Z)
        F_cur, e_cur = tuple(F0), e
        angle = 0.0
        for step in range(40):
            # двугранный угол P при ребре e_cur (в пробной точке)
            fa = [F for F in faces if e_cur <= set(F)]
            n1 = num(normals[tuple(fa[0])][0], s)
            n2 = num(normals[tuple(fa[1])][0], s)
            angle += math.pi - math.acos(np.clip(n1 @ n2 / (np.linalg.norm(n1) * np.linalg.norm(n2)), -1, 1))
            g, Gf, mp = pairing[F_cur]
            inv = {p: q for q, p in mp.items()}
            e_new = frozenset(inv[x] for x in e_cur)
            g_total = g_total.compose(g)
            H = next(F for F in faces if e_new <= set(F) and tuple(F) != Gf)
            F_cur, e_cur = tuple(H), e_new
            if g_total.is_identity() and F_cur == tuple(F0) and e_cur == e:
                break
        else:
            return {"ok": False, "why": f"edge cycle {set(e)} did not close"}
        if abs(angle - 2 * math.pi) > 1e-9:
            return {"ok": False, "why": f"edge cycle angle {angle}"}
        cycles.append(step + 1)
    return {"ok": True, "faces": len(faces), "vertices": len(vlab), "cycleLengths": sorted(set(cycles)),
            "formula": "A" if kind == "tri+tri" and choice[1] != 2 * Z and False else "ok"}


def main(kind):
    results = []
    for ch in chambers(kind):
        r = prove_chamber(kind, ch)
        results.append(r)
        print(kind, ch["rays"], ch["c"], r, flush=True)
    ok = sum(r["ok"] for r in results)
    print(f"\n{kind}: {ok}/{len(results)} chambers proved")
    return ok == len(results)


if __name__ == "__main__":
    for k in (sys.argv[1:] or ["tri+tri", "par+tri"]):
        main(k)
