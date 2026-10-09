"""Поиск изоэдральных замощений: тетраэдр T — фундаментальная область дискретной группы изометрий.

Теорема Пуанкаре о многограннике (евклидов компактный случай): пусть гранями T заданы спаривания
F ↦ F' изометриями γ_F с γ_F(F') = F, γ_F T по другую сторону F, γ_{F'} = γ_F⁻¹ (для F' = F —
инволюция: отражение в плоскости F или полуоборот вокруг оси симметрии треугольника F). Если для каждого
ребра обход «через грани» замыкается — сумма двугранных углов вдоль цикла ровно 2π и итоговый элемент
тождественный, — то образы T образуют разбиение пространства (лицом к лицу).

Углы рациональны, поэтому условие «ровно 2π» проверяется точно (дроби). Совпадение длин рёбер и
тождественность итогового элемента проверяются в 60-значной арифметике (порог 1e-35).
"""
import itertools
import sys
from fractions import Fraction as Fr

import mpmath as mp

from tetra import Tet, EDGES, angle_of, PRINTED_A, SPORADIC_A, EXCLUDED_BY_LP

TOL = mp.mpf(10) ** -35


def faces_of(k):
    """Грань, противоположная вершине k (как кортеж вершин)."""
    return tuple(v for v in (1, 2, 3, 4) if v != k)


def affine_from(P, Q):
    """Аффинное отображение, переводящее точки P[0..3] в Q[0..3] (4×4 матрица)."""
    A = mp.matrix(4, 4)
    B = mp.matrix(4, 4)
    for c in range(4):
        for r in range(3):
            A[r, c] = P[c][r]
            B[r, c] = Q[c][r]
        A[3, c] = 1
        B[3, c] = 1
    return B * mp.inverse(A)


def apply(M, p):
    v = M * mp.matrix([p[0], p[1], p[2], 1])
    return mp.matrix([v[0], v[1], v[2]])


def is_identity(M):
    return all(abs(M[i, j] - (1 if i == j else 0)) < TOL for i in range(4) for j in range(4))


def reflect_point(p, a, b, c):
    n = mp.matrix(list(_cross(b - a, c - a)))
    n = n / mp.norm(n)
    d = ((p - a).T * n)[0]
    return p - 2 * d * n


def _cross(u, v):
    return [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]


class Pairing:
    """γ: грань G → грань F по биекции вершин sigma (словарь вершина G → вершина F)."""

    def __init__(self, T, kF, kG, sigma):
        V = T.V
        F, G = faces_of(kF), faces_of(kG)
        src = [V[g] for g in G] + [V[kG]]
        # образ противоположной вершины G: по другую сторону F относительно T
        img3 = [V[sigma[g]] for g in G]
        # точка, симметричная V[kF] относительно плоскости F — это «другая сторона»
        far = reflect_point(V[kF], *img3)
        # вершина kG должна перейти в точку на расстоянии h_G от плоскости F с другой стороны.
        # Строим её: та изометрия, что переводит G→F по sigma, однозначна с точностью до отражения в F;
        # берём кандидата через решение по расстояниям и выбираем сторону.
        cand = _fourth_point(src, img3)
        sides = [((c - img3[0]).T * mp.matrix(_cross(img3[1] - img3[0], img3[2] - img3[0])))[0] for c in cand]
        ref = ((far - img3[0]).T * mp.matrix(_cross(img3[1] - img3[0], img3[2] - img3[0])))[0]
        pick = [c for c, s in zip(cand, sides) if s * ref > 0]
        assert len(pick) == 1
        self.M = affine_from(src, img3 + [pick[0]])
        self.kF, self.kG, self.sigma = kF, kG, sigma
        # контроль: это изометрия
        L = mp.matrix([[self.M[i, j] for j in range(3)] for i in range(3)])
        assert mp.norm(L.T * L - mp.eye(3)) < TOL


def _fourth_point(src, img3):
    """Две точки, у которых расстояния до img3 те же, что у src[3] до src[0..2]."""
    a, b, c = src[:3]
    p = src[3]
    # координаты p в базисе (b−a, c−a, n) и перенос в образ
    def frame(a, b, c):
        u = b - a
        w = c - a
        n = mp.matrix(_cross(u, w))
        return u, w, n
    u, w, n = frame(a, b, c)
    Mx = mp.matrix([[u[i], w[i], n[i]] for i in range(3)])
    coef = mp.lu_solve(Mx, p - a)
    A, B, C = img3
    u2, w2, n2 = frame(A, B, C)
    base = A + u2 * coef[0] + w2 * coef[1]
    return [base + n2 * coef[2], base - n2 * coef[2]]


def congruences(T, kF, kG, self_pair):
    """Биекции вершин G → F, сохраняющие длины рёбер."""
    F, G = faces_of(kF), faces_of(kG)
    out = []
    for perm in itertools.permutations(F):
        sigma = dict(zip(G, perm))
        if self_pair:
            # инволюция: σ² = id
            if any(sigma[sigma[g]] != g for g in G):
                continue
        ok = all(abs(T.len[tuple(sorted((g1, g2)))] - T.len[tuple(sorted((sigma[g1], sigma[g2])))]) < TOL
                 for g1, g2 in itertools.combinations(G, 2))
        if ok:
            out.append(sigma)
    return out


def face_partitions():
    """Разбиения {1,2,3,4} на пары и одиночки."""
    items = [1, 2, 3, 4]

    def rec(rest):
        if not rest:
            yield []
            return
        a = rest[0]
        for tail in rec(rest[1:]):
            yield [(a,)] + tail
        for b in rest[1:]:
            r2 = [x for x in rest[1:] if x != b]
            for tail in rec(r2):
                yield [(a, b)] + tail
    yield from rec(items)


def schemes(T):
    for part in face_partitions():
        choices = []
        for blk in part:
            if len(blk) == 1:
                k = blk[0]
                cs = congruences(T, k, k, True)
                choices.append([("self", k, s) for s in cs])
            else:
                kF, kG = blk
                cs = congruences(T, kF, kG, False)
                choices.append([("pair", kF, kG, s) for s in cs])
        for combo in itertools.product(*choices):
            pairs = {}
            for ch in combo:
                if ch[0] == "self":
                    _, k, s = ch
                    pairs[k] = Pairing(T, k, k, s)
                else:
                    _, kF, kG, s = ch
                    g = Pairing(T, kF, kG, s)
                    pairs[kF] = g
                    inv_sigma = {v: k for k, v in s.items()}
                    pairs[kG] = Pairing(T, kG, kF, inv_sigma)
            yield combo, pairs


def edge_cycles_ok(T, pairs):
    """Проверка условия Пуанкаре на всех рёбрах; возвращает (ok, причина)."""
    for (i, j) in EDGES:
        k1, k2 = [x for x in (1, 2, 3, 4) if x not in (i, j)]
        # состояние: текущее ребро (i, j) тайла g·T и грань, через которую уходим (противоположная k)
        g = mp.eye(4)
        edge = (i, j)
        face_k = k1  # уходим через грань, противоположную k1 (она содержит i, j)
        total = Fr(0)
        start = (edge, face_k)
        seen_id = 0
        while True:
            total += angle_of(T.angles, *edge)
            # переходим через грань face_k: сосед — g·γ·T, где γ = pairs[face_k] переводит F' в F
            gam = pairs[face_k]
            g = g * gam.M
            # в координатах соседа ребро — прообраз при σ
            inv = {v: k for k, v in gam.sigma.items()}
            edge2 = tuple(sorted((inv[edge[0]], inv[edge[1]])))
            # пришли через грань kG; уходим через другую грань соседа, содержащую edge2
            others = [x for x in (1, 2, 3, 4) if x not in edge2]
            nxt = [x for x in others if x != gam.kG]
            assert len(nxt) == 1
            edge, face_k = edge2, nxt[0]
            if total > 2:
                return False, f"edge {i}{j}: angle sum passes 2π"
            if total == 2:
                if is_identity(g) and (edge, face_k) == start:
                    break
                return False, f"edge {i}{j}: sum 2π but holonomy ≠ id"
            if is_identity(g):
                return False, f"edge {i}{j}: tile returns early (overlap)"
    return True, "ok"


def search(T):
    found = []
    n = 0
    for combo, pairs in schemes(T):
        n += 1
        ok, why = edge_cycles_ok(T, pairs)
        if ok:
            found.append(combo)
    return n, found


def describe(combo):
    parts = []
    for ch in combo:
        if ch[0] == "self":
            _, k, s = ch
            kind = "reflection" if all(s[g] == g for g in s) else "half-turn"
            parts.append(f"F{k}:{kind}")
        else:
            _, kF, kG, s = ch
            parts.append(f"F{kG}→F{kF} {s}")
    return "; ".join(parts)


if __name__ == "__main__":
    from fractions import Fraction
    tests = {
        "Sommerville No.1 (1/2,1/2,1/3,1/3,1/3,1/3)": tuple(map(Fraction, ["1/2", "1/2", "1/3", "1/3", "1/3", "1/3"])),
        "Hill x=1/4": tuple(map(Fraction, ["1/2", "1/2", "1/2", "1/3", "1/4", "1/4"])),
    }
    for name, a in tests.items():
        n, f = search(Tet(a))
        print(name, "schemes", n, "tilings", len(f), [describe(c) for c in f[:3]])
    for idx, a in enumerate(PRINTED_A, 1):
        n, f = search(Tet(a))
        tag = "  (excluded by LP)" if a in EXCLUDED_BY_LP else ""
        print(f"#{idx:2d}{tag} schemes={n} isohedral tilings={len(f)}", [describe(c) for c in f[:2]], flush=True)
