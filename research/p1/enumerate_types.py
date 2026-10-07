"""Шаг 1. Все комбинаторные типы выпуклых многогранников с 4..8 гранями.

plantri -pc3 n перечисляет 3-связные планарные графы с n вершинами (графы 3-многогранников, зеркальные
образы отождествлены); ключ -d выдаёт двойственные, то есть многогранники с n гранями.
Из системы вращений восстанавливаем грани. Результат: data/types.json.
"""
import json
import pathlib
import subprocess

ROOT = pathlib.Path(__file__).resolve().parent
PLANTRI = ROOT.parents[1] / "tools" / "plantri55" / "plantri"
EXPECTED = {4: 1, 5: 2, 6: 7, 7: 34, 8: 257}  # OEIS A000944


def parse(line):
    n, rest = line.split()
    adj = [[ord(c) - ord("a") for c in word] for word in rest.split(",")]
    assert len(adj) == int(n)
    return adj


def faces_of(adj):
    """Грани обходом направленных рёбер: после (u→v) идём в соседа v, следующего за u в вращении v."""
    seen, faces = set(), []
    for u in range(len(adj)):
        for v in adj[u]:
            if (u, v) in seen:
                continue
            face, a, b = [], u, v
            while (a, b) not in seen:
                seen.add((a, b))
                face.append(a)
                rot = adj[b]
                a, b = b, rot[(rot.index(a) + 1) % len(rot)]
            faces.append(face)
    return faces


def main():
    types = []
    for nf in range(4, 9):
        out = subprocess.run([str(PLANTRI), "-pc3", "-d", "-a", str(nf)], capture_output=True, text=True, check=True).stdout
        lines = [l for l in out.splitlines() if l.strip()]
        assert len(lines) == EXPECTED[nf], (nf, len(lines))
        for k, line in enumerate(lines):
            adj = parse(line)
            faces = faces_of(adj)
            nv, ne = len(adj), sum(map(len, adj)) // 2
            assert len(faces) == nf and nv - ne + nf == 2
            types.append({
                "id": f"F{nf}-{k + 1:03d}",
                "faces": faces,
                "adj": adj,
                "nV": nv, "nE": ne, "nF": nf,
                "faceSizes": sorted((len(f) for f in faces), reverse=True),
                "vertexDegrees": sorted((len(a) for a in adj), reverse=True),
            })
    (ROOT / "data").mkdir(exist_ok=True)
    (ROOT / "data" / "types.json").write_text(json.dumps(types))
    by = {}
    for t in types:
        by[t["nF"]] = by.get(t["nF"], 0) + 1
    print("types by face count:", by, "total", len(types))


if __name__ == "__main__":
    main()
