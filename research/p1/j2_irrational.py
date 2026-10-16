"""F8-178: доказательство, что сумма двугранных углов S = 3π + 5δ + 15γ не лежит в ℚπ (значит, инвариант Дена ≠ 0
и семейство не замощает пространство, по теореме Лагариаса–Мюэса / Дебруннера).

δ = arccos(1/√5), γ = двугранный угол J2 между треугольником и пятиугольником, cos γ = √((5+2√5)/15).
S ∈ ℚπ ⇔ w = e^{i(5δ+15γ)} — корень из единицы. Корень из единицы — целое алгебраическое число; минимальный
многочлен w над ℚ считаем точно (sympy). Если у него нецелые коэффициенты, w не корень из единицы.
"""
import sympy as sp

s5 = sp.sqrt(5)
cd, sd = 1 / s5, 2 / s5
cg = sp.sqrt((5 + 2 * s5) / 15)
sg = sp.sqrt((10 - 2 * s5) / 15)
assert sp.simplify(cg**2 + sg**2 - 1) == 0
ed = cd + sp.I * sd
eg = cg + sp.I * sg
x = sp.Symbol("x")
# поэтапно: минимальные многочлены e^{iδ}, e^{iγ}, затем w
for name, z in (("e^{i delta}", ed), ("e^{i gamma}", eg)):
    print(name, sp.factor(sp.minimal_polynomial(z, x)))
w = sp.expand(ed**5 * eg**15)
p = sp.Poly(sp.minimal_polynomial(w, x), x)
coeffs = p.all_coeffs()
lc = coeffs[0]
monic = [sp.Rational(c, 1) / lc for c in coeffs]
print("degree of minimal polynomial of w:", p.degree())
print("monic coefficients:", monic)
print("all integer:", all(c.q == 1 for c in monic))


# ---------- проверка формулы S = 3π + 5δ + 15γ на 50 знаках для двух членов семейства ----------
def refine_and_check():
    import json
    import numpy as np
    import mpmath as mp
    from scipy.spatial import ConvexHull
    mp.mp.dps = 60
    d = [x for x in json.load(open("data/summary.json")) if x["id"] == "F8-178"][0]
    V0 = np.array(d["families"]["example"])
    h = ConvexHull(V0)
    eq = []
    for e in h.equations:
        if not any(np.allclose(e, f, atol=1e-8) for f in eq):
            eq.append(e)
    faces = [sorted(np.where(np.abs(V0 @ e[:3] + e[3]) < 1e-8)[0].tolist()) for e in eq]
    edges = sorted({tuple(sorted((a, b))) for f in faces for a in f for b in f if a < b and any(
        set((a, b)) <= set(g) for g in faces) and sum(set((a, b)) <= set(g) for g in faces) == 2})
    n = len(V0)
    delta, gamma = mp.acos(1 / mp.sqrt(5)), mp.acos(mp.sqrt((5 + 2 * mp.sqrt(5)) / 15))
    target = 3 * mp.pi + 5 * delta + 15 * gamma

    def S_of(X):
        P = [mp.matrix(X[3 * i:3 * i + 3]) for i in range(n)]
        normals = []
        for f in faces:
            a, b, c = P[f[0]], P[f[1]], P[f[2]]
            u, w = b - a, c - a
            nv = mp.matrix([u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0]])
            nv = nv / mp.norm(nv)
            cen = sum((P[k] for k in f), mp.matrix([0, 0, 0])) / len(f)
            inner = sum((P[k] for k in range(n)), mp.matrix([0, 0, 0])) / n
            if (nv.T * (cen - inner))[0] < 0:
                nv = -nv
            normals.append(nv)
        tot = mp.mpf(0)
        for i in range(len(faces)):
            for j in range(i + 1, len(faces)):
                if len(set(faces[i]) & set(faces[j])) == 2:
                    tot += mp.pi - mp.acos((normals[i].T * normals[j])[0])
        return tot

    def system(X, fix_val):
        P = [X[3 * i:3 * i + 3] for i in range(n)]
        res = []
        for a, b in edges:
            res.append(sum((P[a][k] - P[b][k]) ** 2 for k in range(3)) - 1)
        for f in faces:
            if len(f) > 3:
                a, b, c = f[:3]
                for k in f[3:]:
                    M = mp.matrix([[P[b][t] - P[a][t] for t in range(3)], [P[c][t] - P[a][t] for t in range(3)],
                                   [P[k][t] - P[a][t] for t in range(3)]])
                    res.append(mp.det(M))
        # калибровка: v0 = 0, v1 на оси x, v2 в плоскости xy; один параметр семейства фиксируем
        res += [P[0][0], P[0][1], P[0][2], P[1][1], P[1][2], P[2][2], P[3][0] - fix_val]
        return res

    def square(X, fix_val):
        r = system(X, fix_val)
        return r[:len(edges) + 2 + 3] + r[len(edges) + 2 + 4:]  # без одного (зависимого) условия плоскости семиугольника

    for shift in (0.0, 0.01):
        X0 = [mp.mpf(float(v)) for v in V0.ravel()]
        fix_val = X0[9] + shift
        f = lambda *X: square(list(X), fix_val)  # noqa: E731
        X = mp.findroot(f, X0, tol=mp.mpf(10) ** -55, maxsteps=50, method="mdnewton")
        full = system(list(X), fix_val)
        print("max residual of the full system (incl. dropped equation):", mp.nstr(max(abs(v) for v in full), 5))
        S = S_of(list(X))
        print(f"member (shift {shift}): |S - (3pi + 5delta + 15gamma)| = {mp.nstr(abs(S - target), 5)}")


if __name__ == "__main__":
    refine_and_check()
