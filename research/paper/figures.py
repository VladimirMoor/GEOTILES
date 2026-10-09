"""Рисунки для препринта: тела-суммы Минковского и фрагмент разбиения."""
import itertools
import pathlib
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "p1"))
from tilecheck import Tile, crystal  # noqa: E402
from minkowski import polygon, msum  # noqa: E402

COLORS = {3: "#e3a33b", 4: "#4f8fd9", 5: "#4fb38a", 6: "#8d78d8"}


def ordered_faces(T):
    out = []
    for f in T.faces():
        Q = T.vertices[f]
        c = Q.mean(axis=0)
        n = np.linalg.svd(Q - c)[2][-1]
        u = (Q[0] - c) / np.linalg.norm(Q[0] - c)
        w = np.cross(n, u)
        order = np.argsort(np.arctan2((Q - c) @ w, (Q - c) @ u))
        out.append(Q[order])
    return out


def draw(ax, T, shift=np.zeros(3), M=np.eye(3), alpha=0.92, color=None):
    polys = [(F - T.centroid) @ M.T + shift for F in ordered_faces(T)]
    cols = [color or COLORS.get(len(F), "#999") for F in polys]
    pc = Poly3DCollection(polys, facecolors=cols, edgecolors="#1b1d22", linewidths=0.8, alpha=alpha)
    ax.add_collection3d(pc)


def setup(ax, r):
    ax.set_xlim(-r, r); ax.set_ylim(-r, r); ax.set_zlim(-r, r)
    ax.set_box_aspect((1, 1, 1))
    ax.set_axis_off()


def find_member(kind, faces, rng, tries=300):
    """Самый «круглый» член с нужным числом граней: максимум V / R³."""
    best, score = None, -1
    for _ in range(tries):
        P = msum([polygon(k, rng) for k in kind.split("+")])
        T = Tile(P)
        if len(T.planes) != faces:
            continue
        sc = T.volume / T.radius ** 3
        if sc > score:
            best, score = T, sc
    return best


def elongated():
    a = np.array([1.0, 0, 0]); b = np.array([np.cos(1.35), np.sin(1.35), 0])
    c = np.array([0.25, 0.2, 0.947]); c /= np.linalg.norm(c)
    u1 = c - (c @ b) * b; u1 /= np.linalg.norm(u1)
    u2 = c - (c @ a) * a; u2 /= np.linalg.norm(u2)
    h = np.sqrt(3) / 2
    ap1 = c + b / 2 + h * u1; ap2 = a / 2 - h * u2
    box = [np.zeros(3), a, b, a + b, c, a + c, b + c, a + b + c]
    return Tile(np.array(box + [ap1, ap1 + a, ap2, ap2 + b]))


def fig_bodies():
    rng = np.random.default_rng(4)
    bodies = [
        (find_member("tri+tri", 8, rng), "triangle ⊕ triangle\n8 faces"),
        (find_member("tri+tri", 9, rng), "triangle ⊕ triangle\n9 faces"),
        (find_member("rhomb+tri", 10, rng), "rhombus ⊕ triangle\n10 faces"),
        (elongated(), "tri ⊕ tri ⊕ seg (special)\n8 faces"),
    ]
    fig = plt.figure(figsize=(10, 3.0))
    for i, (T, title) in enumerate(bodies):
        ax = fig.add_subplot(1, 4, i + 1, projection="3d")
        draw(ax, T)
        setup(ax, T.radius * 0.72)
        ax.view_init(elev=22, azim=35 + 20 * i)
        ax.set_title(title, fontsize=9)
    fig.tight_layout()
    fig.savefig(HERE / "fig_bodies.pdf", bbox_inches="tight")
    fig.savefig(HERE / "fig_bodies.png", dpi=160, bbox_inches="tight")


def fig_tiling():
    import itertools as it
    rng = np.random.default_rng(9)
    from minkowski import polygon as pg

    def labelings(T):
        for i in range(3):
            o = T[i]; rest = [T[j] for j in range(3) if j != i]
            for a1, a2 in (rest, rest[::-1]):
                yield o, a1 - o, a2 - o

    from tilecheck import coverage
    I = np.eye(3)
    while True:
        T1, T2 = pg("tri", rng), pg("tri", rng)
        P = msum([T1, T2]); T = Tile(P)
        if len(T.planes) != 9 or T.volume / T.radius ** 3 < 0.55:
            continue
        found = None
        for (o1, u1, u2), (o2, v1, v2) in it.product(labelings(T1), labelings(T2)):
            L = [u1 + v1, u1 + v2, u2 - u1]; w = 2 * (o1 + o2)
            if abs(abs(np.linalg.det(np.array(L))) - 2 * T.volume) > 1e-9:
                continue
            R = 3 * T.radius + 2
            cov = coverage(T, crystal(L, [(I, np.zeros(3)), (-I, w)], 1, R=R), samples=300, R=R)
            if cov["exactlyOne"] == cov["samples"]:
                found = (L, w); break
        if found:
            break
    L, w = found
    maps = [(M, t) for M, t in crystal(L, [(I, np.zeros(3)), (-I, w)], 1, R=3 * T.radius)
            if np.linalg.norm(M @ T.centroid + t - T.centroid) < 1.25 * T.radius]
    fig = plt.figure(figsize=(5.2, 5.2))
    ax = fig.add_subplot(1, 1, 1, projection="3d")
    c0 = T.centroid
    for M, t in maps:
        cen = M @ c0 + t
        off = (cen - c0) * 0.55
        col = "#4f8fd9" if np.allclose(M, I) else "#e3a33b"
        draw(ax, T, shift=cen + off - c0 * 0, M=M, alpha=0.9, color=col)
    setup(ax, 1.9 * T.radius)
    ax.view_init(elev=20, azim=40)
    fig.savefig(HERE / "fig_tiling.pdf", bbox_inches="tight")
    fig.savefig(HERE / "fig_tiling.png", dpi=160, bbox_inches="tight")


if __name__ == "__main__":
    fig_bodies()
    fig_tiling()
    print("ok")
