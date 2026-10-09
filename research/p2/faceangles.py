"""Плоские углы граней и их целочисленные соотношения Σ k·β ∈ {π, 2π}."""
import itertools
from fractions import Fraction as Fr
import mpmath as mp
from tetra import PRINTED_A, EXCLUDED_BY_LP, Tet

mp.mp.dps = 60


def face_angles(T):
    """{(грань k, вершина v): угол при v в грани, противоположной k} в долях π."""
    out = {}
    for k in (1, 2, 3, 4):
        F = [v for v in (1, 2, 3, 4) if v != k]
        for v in F:
            a, b = [w for w in F if w != v]
            u, w = T.V[a] - T.V[v], T.V[b] - T.V[v]
            out[(k, v)] = mp.acos((u.T * w)[0] / (mp.norm(u) * mp.norm(w))) / mp.pi
    return out


def distinct(vals, tol=mp.mpf(10) ** -40):
    cls = []
    for key, x in vals.items():
        for c in cls:
            if abs(c[0] - x) < tol:
                c[1].append(key)
                break
        else:
            cls.append([x, [key]])
    return cls


def relations(xs, targets=(1, 2), kmax=None):
    """Все k ≥ 0 с Σ k_i x_i ∈ targets (численно, точность 1e-40)."""
    out = []
    n = len(xs)

    def rec(i, acc, cur):
        if i == n:
            for t in targets:
                if abs(acc - t) < mp.mpf(10) ** -40:
                    out.append((t, tuple(cur)))
            return
        k = 0
        while acc + k * xs[i] <= max(targets) + mp.mpf(10) ** -30:
            rec(i + 1, acc + k * xs[i], cur + [k])
            k += 1
    rec(0, mp.mpf(0), [])
    return out


if __name__ == "__main__":
    for idx, a in enumerate(PRINTED_A, 1):
        if a in EXCLUDED_BY_LP:
            continue
        T = Tet(a)
        fa = face_angles(T)
        cls = distinct(fa)
        xs = [c[0] for c in cls]
        rel = relations(xs)
        # «тривиальные» — кратные полным наборам углов одной грани; покажем все, но пометим рациональные углы
        rat = [mp.nstr(x, 8) for x in xs]
        print(f"#{idx:2d}: {len(xs)} distinct face angles {rat}; relations to π/2π: {len(rel)}")
        for t, k in rel[:12]:
            print("     ", t, "π =", " + ".join(f"{m}·β{j}" for j, m in enumerate(k) if m))
