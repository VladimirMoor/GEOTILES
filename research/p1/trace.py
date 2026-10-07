"""Шаг 5. Однопараметрические семейства: трассировка кривой и поиск членов, проходящих фильтр углов.

Кривая решений гауге-фиксированной системы (realize.System) прослеживается предиктором по ядру якобиана
и корректором Гаусса–Ньютона в обе стороны, пока тело остаётся строго выпуклым.

Фильтр для члена p: для каждого ребра e существует целый вектор k ≥ 0 с k_e ≥ 1 и
Σ kᵢ θᵢ(p) = 360° или 180°. Число слагаемых ограничено: Σ kᵢ ≤ ⌊360° / min θ(p)⌋, поэтому
перебор соотношений полон для каждого члена. Перебор идёт по различным функциям углов
(рёбра с тождественно равными углами склеены), корни ищутся по смене знака на сетке и уточняются.

Запуск: python trace.py F6-002 [F7-016 ...]   →  data/trace_<id>.json
"""
import itertools
import json
import pathlib
import sys

import numpy as np
from scipy.optimize import brentq

ROOT = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from realize import System, verify  # noqa: E402


def correct(sys_, z, tangent=None, z_pred=None):
    for _ in range(30):
        r = sys_.residuals(z)
        J = sys_.jacobian(z)
        if tangent is not None:
            r = np.append(r, tangent @ (z - z_pred))
            J = np.vstack([J, tangent])
        if np.max(np.abs(r)) < 1e-13:
            break
        dz = np.linalg.lstsq(J, -r, rcond=None)[0]
        z = z + dz
    return z, np.max(np.abs(sys_.residuals(z)))


def null_dir(sys_, z, prev=None):
    J = sys_.jacobian(z)
    _, s, vt = np.linalg.svd(J)
    v = vt[-1]
    if prev is not None and v @ prev < 0:
        v = -v
    return v / np.linalg.norm(v), s[-2] / s[0]


def trace(t, X0, h=0.01, max_steps=4000):
    sys_ = System(t)
    from realize import gauge_fix
    z0 = sys_.pack(gauge_fix(np.array(X0), sys_))
    z0, err = correct(sys_, z0)
    branch_pts = []
    for sign in (1, -1):
        z, v = z0.copy(), None
        v, gap = null_dir(sys_, z)
        v = sign * v
        pts = []
        for _ in range(max_steps):
            X = sys_.split(z)[0]
            g = verify(X, t["faces"], sys_.E)
            if g is None:
                break
            pts.append({"z": z.copy(), "dihedrals": np.array(g["dihedrals"]), "X": X.copy()})
            zp = z + h * v
            zn, err = correct(sys_, zp, v, zp)
            if err > 1e-9:
                break
            vn, gap = null_dir(sys_, zn, v)
            z, v = zn, vn
        branch_pts.append(pts)
    pts = branch_pts[1][::-1] + branch_pts[0][1:]
    return sys_, pts


def scan(D, tol_deg=1e-7):
    """D: точки × рёбра (градусы). Возвращает кандидатные индексы-интервалы и соотношения."""
    nP, nE = D.shape
    # группы тождественно равных функций углов
    groups = []
    for e in range(nE):
        for g in groups:
            if np.max(np.abs(D[:, g[0]] - D[:, e])) < 1e-6:
                g.append(e)
                break
        else:
            groups.append([e])
    F = np.array([D[:, g[0]] for g in groups]).T  # точки × функции
    m = F.shape[1]
    K = int(360 // max(F.min(), 1.0))
    # для каждой функции e: множество точек (интервалов сетки), где возможно соотношение с k_e ≥ 1
    cand = {}
    identically = set()
    for e in range(m):
        roots = []
        for total in range(1, K + 1):
            for combo in itertools.combinations_with_replacement(range(m), total):
                if e not in combo:
                    continue
                k = np.bincount(combo, minlength=m)
                val = F @ k
                for target in (360.0, 180.0):
                    g = val - target
                    if np.max(np.abs(g)) < 1e-6:
                        identically.add(e)
                        continue
                    sgn = np.sign(g)
                    for i in np.nonzero(sgn[:-1] * sgn[1:] <= 0)[0]:
                        roots.append((int(i), tuple(int(x) for x in k), target))
        cand[e] = roots
    return groups, F, K, cand, identically


def angles_of(sys_, t, z):
    g = verify(sys_.split(z)[0], t["faces"], sys_.E)
    return None if g is None else np.array(g["dihedrals"])


def refine(sys_, t, z, edges_k, target):
    """Точный член семейства: система + соотношение Σ k_e θ_e = target (по рёбрам-представителям)."""
    def rel(zz):
        th = angles_of(sys_, t, zz)
        return None if th is None else float(sum(k * th[e] for e, k in edges_k) - target)
    for _ in range(40):
        r = sys_.residuals(z)
        g = rel(z)
        if g is None:
            return None
        J = sys_.jacobian(z)
        h = 1e-7
        grad = np.array([(rel(z + h * ei) - g) / h if rel(z + h * ei) is not None else 0 for ei in np.eye(len(z))])
        rr = np.append(r, np.radians(g))
        JJ = np.vstack([J, np.radians(grad)])
        if np.max(np.abs(rr)) < 1e-13:
            break
        z = z + np.linalg.lstsq(JJ, -rr, rcond=None)[0]
    if np.max(np.abs(sys_.residuals(z))) > 1e-10 or abs(rel(z)) > 1e-8:
        return None
    return z


def full_filter(th, tol=1e-7):
    sys.path.insert(0, str(ROOT))
    from analyze import dihedral_filter
    return dihedral_filter(list(th))


def identify(total_deg):
    """PSLQ: Σθ как целая комбинация π и стандартных иррациональных углов."""
    import mpmath as mp
    mp.mp.dps = 15
    # arccos(1/3): тетраэдр, октаэдр, J1, J12; arccos(1/√5) = arctan 2: додекаэдр/J2-родственные углы;
    # arccos(√((5+2√5)/15)): угол пятиугольной пирамиды J2 при основании
    consts = {"pi": mp.pi, "acos(1/3)": mp.acos(mp.mpf(1) / 3), "acos(1/sqrt5)": mp.acos(1 / mp.sqrt(5)),
              "J2base": mp.acos(mp.sqrt((5 + 2 * mp.sqrt(5)) / 15))}
    x = mp.radians(total_deg)
    rel = mp.pslq([x] + list(consts.values()), maxcoeff=200, maxsteps=10 ** 5, tol=mp.mpf(10) ** -10)
    if rel is None:
        return None
    return {name: rel[i + 1] for i, name in enumerate(consts)} | {"sum": rel[0]}


def exact_candidates(t, sys_, pts, groups, F, cand, ident):
    need = [e for e in range(F.shape[1]) if e not in ident]
    found = []
    seen = []
    for e in need:
        for i, k, tg in cand[e]:
            edges_k = [(groups[j][0], kj) for j, kj in enumerate(k) if kj]
            for start in (i, i + 1):
                z = refine(sys_, t, pts[start]["z"].copy(), edges_k, tg)
                if z is None:
                    continue
                th = angles_of(sys_, t, z)
                key = tuple(np.round(np.sort(th), 6))
                if any(np.allclose(key, s, atol=1e-6) for s in seen):
                    continue
                seen.append(key)
                f = full_filter(th)
                if f["ok"]:
                    found.append({"z": z.tolist(), "X": sys_.split(z)[0].tolist(), "angles": sorted(np.round(th, 8).tolist()),
                                  "sumDeg": float(np.sum(th)), "dehn": identify(float(np.sum(th))), "pin": [edges_k, tg]})
    return found


def main(ids):
    types = {t["id"]: t for t in json.loads((ROOT / "data" / "types.json").read_text())}
    R = json.loads((ROOT / "data" / "realizations.json").read_text())
    for tid in ids:
        t = types[tid]
        fam = [s for s in R[tid]["solutions"] if s["flex"] == 1]
        fam.sort(key=lambda s: -min(s["dihedrals"]))
        sys_, pts = trace(t, fam[0]["X"])
        D = np.array([p["dihedrals"] for p in pts])
        print(f"{tid}: traced {len(pts)} convex members; angle ranges "
              f"min θ ∈ [{D.min(axis=1).min():.3f}, {D.min(axis=1).max():.3f}]°")
        groups, F, K, cand, ident = scan(D)
        print(f"  {F.shape[1]} distinct angle functions, max tiles around an edge K={K}; "
              f"identically satisfied functions: {sorted(ident)}")
        # точки сетки, где каждая функция имеет корень соотношения в том же интервале
        need = [e for e in range(F.shape[1]) if e not in ident]
        intervals = None
        for e in need:
            s = {i for i, _, _ in cand[e]}
            intervals = s if intervals is None else intervals & s
        intervals = sorted(intervals or [])
        print(f"  grid intervals where every non-trivial edge class has a relation: {len(intervals)}")
        out = []
        for i in intervals[:50]:
            rel = {e: [(k, tg) for j, k, tg in cand[e] if j == i][:3] for e in need}
            out.append({"interval": i, "angles": F[i].round(4).tolist(), "relations": {str(e): rel[e] for e in rel}})
            print(f"   interval {i}: angles {F[i].round(3).tolist()}")
            for e in need:
                print(f"     f{e}: {rel[e][:2]}")
        exact = exact_candidates(t, sys_, pts, groups, F, cand, ident)
        print(f"  exact members passing the full filter: {len(exact)}")
        for x in exact:
            print(f"   angles {np.round(x['angles'], 4).tolist()}\n   Σθ = {x['sumDeg']:.8f}°  PSLQ: {x['dehn']}")
        (ROOT / "data" / f"trace_{tid}.json").write_text(json.dumps({
            "members": len(pts), "groups": groups, "K": K, "identically": sorted(ident),
            "angles": F.round(6).tolist(), "candidates": out, "exact": exact}))


if __name__ == "__main__":
    main(sys.argv[1:])
