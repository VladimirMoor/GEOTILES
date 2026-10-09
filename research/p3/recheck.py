"""Перепроверка «неопределённых» и высоких фигур из быстрого прогона.

Выход sat: пары строк «D <координаты>» и строка класса: «~ Hc Hh» (не замощает), «I» (изоэдрально),
«!» (не определено на заданном уровне). Каждую «!»-фигуру и каждую с Hc ≥ MINH прогоняем отдельно с
большим -maxlevel и лимитом времени. Фигура с H ≥ 7 после такой проверки получит «~ 7 …».
"""
import pathlib
import subprocess
import sys
import tempfile

SAT = "/Users/vladimirmuravev/Projects2026/GEOTILES/tools/heesch-sat/src/sat"


def shapes(path):
    lines = pathlib.Path(path).read_text().splitlines()
    out = []
    for i in range(0, len(lines) - 1):
        if lines[i].startswith("D ") and not lines[i + 1].startswith("D"):
            out.append((lines[i][2:], lines[i + 1].split()))
    return out


def run_one(coords, maxlevel, timeout):
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
        f.write("D? " + coords + "\n")
        name = f.name
    try:
        r = subprocess.run([SAT, "-isohedral", "-periodic", "-maxlevel", str(maxlevel), name], capture_output=True, text=True,
                           timeout=timeout)
        cls = [l for l in r.stdout.splitlines() if not l.startswith("D")]
        return cls[0] if cls else "?"
    except subprocess.TimeoutExpired:
        return f"timeout>{timeout}s"


if __name__ == "__main__":
    path, maxlevel, timeout, minh = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    for coords, cls in shapes(path):
        hc = int(cls[1]) if cls[0] == "~" else None
        if cls[0] == "!" or (hc is not None and hc >= minh):
            res = run_one(coords, maxlevel, timeout)
            print(f"{path.split('/')[-1]}: {cls} -> {res} | {coords}", flush=True)
