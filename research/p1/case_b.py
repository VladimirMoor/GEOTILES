"""Случай (B) комбинаторной леммы: один треугольник и S ∈ {3, 4} отрезков.

Зона отрезка s₁ — окружность C = s₁^⊥. На ней лежат: по паре противоположных точек от каждого из
остальных S−1 отрезков (одна прямая, два вхождения) и три точки от треугольника (три разные прямые, по
одному вхождению). Перебираем ВСЕ абстрактные циклические расположения — это надмножество реализуемых:
пары точек отрезков стоят в позициях k·π/(S−1) и k·π/(S−1)+π (порядок прямых без ограничения общности),
а три точки треугольника раскладываются по 2(S−1) промежуткам.

Для каждого расположения строим граф формул (дуги — различные пары соседних прямых) и находим дуги,
лежащие на ориентированном цикле со значением Σ(π − φ) ≤ 2π. Значение цикла определяется комбинаторикой:
Σφ кратна π и вычисляется по конкретным углам точек.
Лемма (B): в каждом расположении есть дуга, не лежащая ни на одном таком цикле (свидетель).
"""
import itertools
import math


def arrangements(S):
    G = 2 * (S - 1)
    seg = [(k * math.pi / (S - 1), ("s", k)) for k in range(S - 1)] + \
          [(k * math.pi / (S - 1) + math.pi, ("s", k)) for k in range(S - 1)]
    gap = 2 * math.pi / G
    for slots in itertools.combinations_with_replacement(range(G), 3):
        pts = list(seg)
        count = {}
        for t, g in enumerate(slots):
            count[g] = count.get(g, 0) + 1
            pos = g * gap + gap * count[g] / 4.0  # внутри промежутка, без совпадений
            pts.append((pos, ("t", t)))
        pts.sort()
        yield slots, pts


def witnesses(pts):
    n = len(pts)
    arcs = {}
    for i in range(n):
        a, b = pts[i], pts[(i + 1) % n]
        phi = (b[0] - a[0]) % (2 * math.pi)
        arcs.setdefault((a[1], b[1]), phi)  # одинаковые пары прямых — одна формула (с тем же φ mod π)
    out = {}
    for (u, v) in arcs:
        out.setdefault(u, []).append(v)
    good = set()

    def dfs(start, node, path, visited):
        for w in out.get(node, []):
            arc = (node, w)
            if w == start:
                cyc = path + [arc]
                val = sum(math.pi - arcs[x] for x in cyc)
                if val <= 2 * math.pi + 1e-9:
                    good.update(cyc)
            elif w not in visited:
                dfs(start, w, path + [arc], visited | {w})

    for s in out:
        dfs(s, s, [], {s})
    return [a for a in arcs if a not in good], len(arcs)


if __name__ == "__main__":
    for S in (3, 4, 5):
        total, ok, worst = 0, 0, None
        for slots, pts in arrangements(S):
            total += 1
            w, narcs = witnesses(pts)
            if w:
                ok += 1
            if worst is None or len(w) < worst[0]:
                worst = (len(w), narcs, slots)
        print(f"S={S} (Σdim={2 + S}): {ok}/{total} arrangements have a witness arc; fewest witnesses {worst}")
