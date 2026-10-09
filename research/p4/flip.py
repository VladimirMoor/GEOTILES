"""Не-DV стереоэдры через флипы двойственного комплекса (эксперимент на разбиении Энгеля).

Изоэдральное face-to-face разбиение с группой Γ описывается Γ-инвариантным двойственным комплексом
на орбите Γx: вершины-тетраэдры τ = {x, a, b, c} (наборы точек орбиты), рёбра {x, y} ↔ грани клетки x.
Для DV-разбиения это триангуляция Делоне. Плитка P (клетка x): вершины — тетраэдры, содержащие x;
грань F_y — тетраэдры, содержащие x и y; спаривание: элемент g с g x = y переводит τ ∋ x, g⁻¹y… — в
терминах множеств: g⁻¹(τ) — снова тетраэдр, содержащий x, если τ ∋ y.

Флип 2→3 у x: тетраэдры τ₁ = {x, a, b, c} и τ₂ = {a, b, c, e} (соседние по треугольнику abc) заменяются
на {x, e, a, b}, {x, e, b, c}, {x, e, c, a} — вместе со всеми Γ-образами. Появляется ребро {x, e}, т.е.
новая треугольная грань P (и парная ей).

Реализация: неизвестные — координаты вершин P (по одной на тетраэдр, содержащий x), уравнения —
спаривание (через g⁻¹) и плоскостность граней. Ищем решение рядом с вырожденной точкой (новые вершины
в старых), где новая грань имеет площадь > 0, и проверяем выпуклость и объём (= объёму фундаментальной
области ⇒ разбиение).
"""
import itertools
import sys
from fractions import Fraction as Fr

import numpy as np
from scipy.optimize import least_squares
from scipy.spatial import ConvexHull, HalfspaceIntersection

from dvfast import Group, Evaluator


class Orbit:
    """Точки орбиты Γx (фракционные, float) и элементы g с g x = y."""

    def __init__(self, G, x):
        self.G = G
        self.x = np.asarray(x, float)
        self.cache = {}

    def key(self, y):
        return tuple(int(v) for v in np.round(np.asarray(y, float) * 1e6))

    def element(self, y):
        k = self.key(y)
        if k in self.cache:
            return self.cache[k]
        for R, t in zip(self.G.R, self.G.t):
            s = y - (R @ self.x + t)
            if np.allclose(s, np.round(s), atol=1e-7):
                g = (R, t + np.round(s))
                self.cache[k] = g
                return g
        raise KeyError("not in orbit")


def apply(g, p):
    R, t = g
    return R @ p + t


def inverse(g):
    R, t = g
    Ri = np.linalg.inv(R)
    return (Ri, -Ri @ t)


def tet_key(pts, orb):
    return frozenset(orb.key(p) for p in pts)


class Complex:
    """Γ-инвариантный двойственный комплекс, хранится как множество тетраэдров (каждый — frozenset
    ключей точек орбиты), содержащих x, и словарь ключ → координаты точки."""

    def __init__(self, orb):
        self.orb = orb
        self.pts = {}
        self.tets = set()

    def add_point(self, p):
        k = self.orb.key(p)
        self.pts[k] = np.asarray(p, float)
        return k

    def images_containing_x(self, tet):
        """Все Γ-образы тетраэдра tet (набор координат), содержащие x: g⁻¹(tet) для каждого члена y = g x."""
        out = []
        for y in tet:
            g = self.orb.element(y)
            gi = inverse(g)
            img = [apply(gi, p) for p in tet]
            out.append(img)
        return out

    def add_orbit(self, tet):
        for img in self.images_containing_x(tet):
            for p in img:
                self.add_point(p)
            self.tets.add(tet_key(img, self.orb))

    def remove_orbit(self, tet):
        for img in self.images_containing_x(tet):
            self.tets.discard(tet_key(img, self.orb))


def from_dv_exact(hm, x, c2):
    """Двойственный комплекс из точной DV-ячейки (рациональные вершины и инцидентности)."""
    from deform_exact import exact_cell
    from fractions import Fraction as Fr
    xr = tuple(Fr(v) for v in x) if not isinstance(x[0], Fr) else x
    G, M, planes, verts, vdef, faces = exact_cell(hm, xr, c2)
    xf = np.array([float(v) for v in xr])
    orb = Orbit(G, xf)
    C = Complex(orb)
    C.add_point(xf)
    facet_sites = {pi: np.array([float(v) for v in planes[pi][4]]) for pi, on in faces}
    inc = {vi: [] for vi in range(len(verts))}
    for pi, on in faces:
        for vi in on:
            inc[vi].append(pi)
    from exact import dot
    positions = {}
    for vi, v in enumerate(verts):
        # все точки орбиты, равноудалённые от вершины (а не только соседи по граням)
        eq = [pi for pi, pl in enumerate(planes) if dot(pl[0], v) == pl[1]]
        tet = [xf] + [np.array([float(t) for t in planes[pi][4]]) for pi in eq]
        for p in tet:
            C.add_point(p)
        k = tet_key(tet, orb)
        C.tets.add(k)
        positions[k] = np.array([float(t) for t in v])
    return C, positions


def from_dv(E, x):
    """Двойственный комплекс DV-разбиения: для каждой вершины клетки x — её тетраэдр {x и три соседа}."""
    G = E.G
    orb = Orbit(G, x)
    X, A, b, V = E.cell_of(np.asarray(x, float))
    Binv = np.linalg.inv(E.B)
    Vu = []
    for v in V:
        if not any(np.linalg.norm(v - u) < 1e-8 for u in Vu):
            Vu.append(v)
    Vu = np.array(Vu)
    nA = np.linalg.norm(A, axis=1)
    facet = np.zeros(len(A), bool)
    for i in range(len(A)):
        on = np.abs(Vu @ A[i] - b[i]) / nA[i] < 1e-9
        facet[i] = on.sum() >= 3 and np.linalg.matrix_rank(Vu[on] - Vu[on][0], tol=1e-8) >= 2
    C = Complex(orb)
    xk = C.add_point(np.asarray(x, float))
    positions = {}
    for v in Vu:
        act = np.where((np.abs(A @ v - b) / nA < 1e-9) & facet)[0]
        if len(act) < 3:
            raise RuntimeError("bad vertex")
        ys = [Binv @ (A[i] + X) for i in act]
        tet = [np.asarray(x, float)] + ys
        for p in tet:
            C.add_point(p)
        k = tet_key(tet, orb)
        C.tets.add(k)
        positions[k] = Binv @ v  # фракционные координаты вершины
    return C, positions


def structure(C):
    """Грани P: для каждого соседа y — вершины P (двойственные ячейки), содержащие y (порядок не нужен:
    плоскостность и спаривание от порядка не зависят)."""
    xk = C.orb.key(C.orb.x)
    faces = {}
    for t in C.tets:
        for y in t:
            if y != xk:
                faces.setdefault(y, []).append(t)
    # соседи, касающиеся P только вершиной или ребром (< 3 вершин), гранями не являются
    return {y: v for y, v in faces.items() if len(v) >= 3}


def flips_at_x(C):
    """Кандидаты флипа 2→3, добавляющие ребро {x, e}: (τ₁, τ₂)."""
    xk = C.orb.key(C.orb.x)
    out = []
    for t1 in C.tets:
        if len(t1) != 4:
            continue
        abc = [k for k in t1 if k != xk]
        a = abc[0]
        # τ₂ = {a, b, c, e}: его образ g_a⁻¹ τ₂ содержит x и g_a⁻¹(b), g_a⁻¹(c)
        ga = C.orb.element(C.pts[a])
        gi = inverse(ga)
        img_abc = [tet_key([apply(gi, C.pts[k])], C.orb) for k in abc]
        img_abc = set().union(*img_abc)  # = {x, g⁻¹b, g⁻¹c}
        img_t1 = tet_key([apply(gi, C.pts[k]) for k in t1], C.orb)
        for t in C.tets:
            if len(t) == 4 and img_abc <= t and t != img_t1:
                # τ₂ = g_a(t)
                e_img = [k for k in t if k not in img_abc][0]
                e = apply(ga, C.pts[e_img])
                out.append((t1, abc, e))
                break
    return out


def do_flip(C, t1, abc, e):
    D = Complex(C.orb)
    D.pts = dict(C.pts)
    D.tets = set(C.tets)
    pa, pb, pc = (C.pts[k] for k in abc)
    px = C.orb.x
    D.remove_orbit([px, pa, pb, pc])
    D.remove_orbit([pa, pb, pc, e])
    for tri in ((pa, pb), (pb, pc), (pc, pa)):
        D.add_orbit([px, e, tri[0], tri[1]])
    return D


def compose(g, h):
    R1, t1 = g
    R2, t2 = h
    return (R1 @ R2, R1 @ t2 + t1)


def realize(C, faces, init, target_area=1e-3, tri_key=None):
    """Неизвестные — по одной вершине на класс (v_t = h_t(rep)); спаривание выполняется автоматически,
    решаем плоскостность граней (+ площадь новой грани). Несогласованные циклы — дополнительные невязки."""
    tets = sorted(C.tets, key=lambda t: sorted(t))
    idx = {t: i for i, t in enumerate(tets)}
    edges = {}
    for y, cyc in faces.items():
        g = C.orb.element(C.pts[y])
        gi = inverse(g)
        for t in cyc:
            img = tet_key([apply(gi, C.pts[k]) for k in t], C.orb)
            if img not in idx:
                return None
            edges.setdefault(idx[img], []).append((idx[t], g))          # v_t = g(v_img)
            edges.setdefault(idx[t], []).append((idx[img], gi))         # v_img = g⁻¹(v_t)
    n = len(tets)
    h = [None] * n
    rep_of = [None] * n
    reps = []
    loops = []  # (i, map_i, j, map_j): h_i(rep) должно совпасть с map(h_j...) — проверяем как невязку
    I = (np.eye(3), np.zeros(3))
    for s0 in range(n):
        if h[s0] is not None:
            continue
        r = len(reps)
        reps.append(s0)
        h[s0] = I
        rep_of[s0] = r
        stack = [s0]
        while stack:
            a = stack.pop()
            for b, g in edges.get(a, []):
                hb = compose(g, h[a])
                if h[b] is None:
                    h[b] = hb
                    rep_of[b] = r
                    stack.append(b)
                elif not (np.allclose(h[b][0], hb[0]) and np.allclose(h[b][1], hb[1])):
                    loops.append((b, hb))
    v0 = np.array([init[tets[i]] for i in reps]).ravel()

    def positions(v):
        R = v.reshape(-1, 3)
        return np.array([apply(h[i], R[rep_of[i]]) for i in range(n)])

    def resid(v):
        V = positions(v)
        R = v.reshape(-1, 3)
        r = []
        for b, hb in loops:
            r.extend(V[b] - apply(hb, R[rep_of[b]]))
        for y, cyc in faces.items():
            P = V[[idx[t] for t in cyc]]
            if len(P) > 3:
                c = P.mean(axis=0)
                r.append(np.linalg.svd(P - c, compute_uv=False)[-1] * 10)
        if tri_key is not None and tri_key in faces:
            P = V[[idx[t] for t in faces[tri_key]]]
            area = 0.5 * np.linalg.norm(np.cross(P[1] - P[0], P[2] - P[0]))
            r.append((area - target_area) * 10)
        return np.array(r)

    sol = least_squares(resid, v0, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=3000)
    sol.x_full = positions(sol.x).ravel()
    return sol, tets, idx, len(reps), len(loops)


def check_tile(V, faces, idx, E, orb):
    """Выпуклость (все вершины по одну сторону каждой грани) и объём против объёма фундаментальной области."""
    Vc = (E.B @ V.T).T
    worst = 0.0
    for y, cyc in faces.items():
        P = Vc[[idx[t] for t in cyc]]
        c = P.mean(axis=0)
        n = np.linalg.svd(P - c)[2][-1]
        s = (Vc - c) @ n
        if abs(s.min()) > abs(s.max()):
            s = -s
        worst = max(worst, -s.min() if s.min() < 0 else 0)  # нарушение выпуклости
    vol = ConvexHull(Vc).volume
    Vfd = abs(np.linalg.det(E.B)) / len(E.G.R)
    return worst, vol, Vfd


def run(hm, x, c=1.0, areas=(1e-6, 1e-5, 1e-4, 1e-3)):
    G = Group(hm)
    E = Evaluator(G, c)
    from fractions import Fraction as Fr
    C, pos = from_dv_exact(hm, tuple(Fr(v).limit_denominator(10**6) for v in x), Fr(c).limit_denominator(10**6) ** 2)
    F0 = structure(C)
    print(f"{hm}: DV complex: {len(C.tets)} vertices of P, {len(F0)} facets", flush=True)
    cands = flips_at_x(C)
    # уникальные флипы по ключу нового соседа e
    seen = set()
    uniq = []
    for t1, abc, e in cands:
        k = C.orb.key(e)
        if k in seen:
            continue
        seen.add(k)
        uniq.append((t1, abc, e))
    print(f"  candidate flips adding a neighbour: {len(uniq)}", flush=True)
    results = []
    for n, (t1, abc, e) in enumerate(uniq):
        D = do_flip(C, t1, abc, e)
        F = structure(D)
        if F is None:
            print(f"  flip {n}: invalid complex")
            continue
        # начальные позиции: старые вершины на месте; новые — в старой вершине τ₁ (и образах)
        init = {}
        old_by_set = {t: pos[t] for t in pos}
        for t in D.tets:
            if t in old_by_set:
                init[t] = old_by_set[t]
            else:
                # ближайшая старая вершина с наибольшим пересечением множеств
                best = max(pos, key=lambda u: len(u & t))
                init[t] = pos[best]
        ek = C.orb.key(e)
        out = []
        for A_ in areas:
            r = realize(D, F, init, target_area=A_, tri_key=ek if ek in F else None)
            if r is None:
                out.append("no-pairing")
                break
            sol, tets, idx, nrep, nloop = r
            V = sol.x_full.reshape(-1, 3)
            worst, vol, Vfd = check_tile(V, F, idx, E, C.orb)
            out.append(f"[{nrep} classes, {nloop} loops] area {A_:.0e}: resid {np.abs(sol.fun).max():.1e}, convex-viol {worst:.1e}, vol/V0 {vol / Vfd:.6f}")
            init = {t: V[idx[t]] for t in tets}
        print(f"  flip {n} (new facets {len(F)}): " + " | ".join(out), flush=True)
        results.append((n, len(F), out))
    return results


if __name__ == "__main__":
    run("I 41 3 2", (427 / 6984, 761 / 6984, 1421 / 6984))
