"""Разложение идеала типа (с переменными Рабиновича для граней и рёбер) на минимальные простые компоненты
в Singular; для каждой компоненты — размерность и элиминированные соотношения на координаты вершин."""
import json
import pathlib
import re
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from exact import build  # noqa: E402

fix = lambda s: re.sub(r"([xnseh]\d+)_(\d+)", r"\1_\2", s)


def decompose(tid, timeout=1800):
    T = {t["id"]: t for t in json.loads((ROOT / "data" / "types.json").read_text())}
    t = T[tid]
    v, e = build(t, rabinowitsch=True)
    xs = [x for x in v if x.startswith("x")]
    others = [x for x in v if not x.startswith("x")]
    order = others + xs  # исключаемые переменные первыми
    script = f'''LIB "primdec.lib"; LIB "elim.lib";
ring r=0,({",".join(order)}),dp;
ideal I={",".join(e)};
list L=minAssGTZ(I);
int i;
for(i=1;i<=size(L);i++){{
  ideal J=std(L[i]);
  ideal K=eliminate(J, {"*".join(others) if others else "1"});
  print("COMPONENT "+string(i)+" dim "+string(dim(J)));
  print(K);
}}
quit;'''
    path = ROOT / "data" / f"decomp_{tid}.sing"
    path.write_text(script)
    out = subprocess.run(["Singular", "-q", str(path)], capture_output=True, text=True, timeout=timeout)
    return out.stdout + out.stderr[-500:]


if __name__ == "__main__":
    print(decompose(sys.argv[1]))
