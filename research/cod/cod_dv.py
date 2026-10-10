"""Ячейки Дирихле–Вороного орбит атомов в кристаллах Crystallography Open Database (COD).

Для каждой записи COD в данной группе: ячейка, операции симметрии и атомы из CIF (gemmi). Для каждой независимой
позиции атома (кроме H/D) строим DV-ячейку её орбиты под группой (dvfast из research/p4): это стереоэдр той же
природы, что рекорды задачи 4 (38 у Энгеля), только точка выбрана природой.

Запись пропускается, если операции из CIF не совпадают со стандартной установкой группы (другая установка или
начало координат) или метрика не подходит (a ≠ b для тетрагональной и т. п.). CIF кэшируются в data/cache/.
Запуск: nice -n 10 python cod_dv.py "I 41 3 2"
"""
import json
import pathlib
import sys
import time
import urllib.parse
import urllib.request

import gemmi
import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "p4"))
from dvfast import Group, Evaluator  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
CACHE = HERE / "data" / "cache"
CACHE.mkdir(parents=True, exist_ok=True)
UA = {"User-Agent": "GEOTILES-research/1.0 (academic; github.com/VladimirMoor/GEOTILES)"}


def fetch(url, name=None, tries=6):
    f = CACHE / name if name else None
    if f and f.exists():
        return f.read_bytes()
    for attempt in range(tries):
        try:
            data = urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=90).read()
            if f:
                f.write_bytes(data)
            time.sleep(0.3)
            return data
        except Exception:
            time.sleep(3 + 5 * attempt)
    raise RuntimeError(f"fetch failed: {url}")


def entries_in_group(sg):
    url = "https://www.crystallography.net/cod/result?spacegroup=" + urllib.parse.quote(sg) + "&format=lst"
    for _ in range(5):
        ids = [s.strip() for s in fetch(url).decode().split() if s.strip().isdigit()]
        if ids:
            return ids
        time.sleep(10)
    return []


def same_ops(small, G):
    """Совпадают ли операции из CIF с операциями стандартной установки (по модулю решётки)."""
    ops = [gemmi.Op(s) for s in small.symops] if small.symops else None
    if not ops:
        return small.spacegroup_hm.replace(" ", "") == G.hm.replace(" ", "")
    std = {op.wrap().triplet() for op in gemmi.SpaceGroup(G.hm).operations()}
    got = {op.wrap().triplet() for op in ops}
    return std == got


def metric_c(cell, system):
    a, b, c, al, be, ga = cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma
    ok = lambda u, v, tol=2e-3: abs(u - v) <= tol * max(abs(u), abs(v))  # noqa: E731
    if system == "cubic" and ok(a, b) and ok(a, c) and ok(al, 90) and ok(be, 90) and ok(ga, 90):
        return 1.0
    if system == "tetragonal" and ok(a, b) and ok(al, 90) and ok(be, 90) and ok(ga, 90):
        return c / a
    if system in ("hexagonal", "trigonal") and ok(a, b) and ok(al, 90) and ok(be, 90) and ok(ga, 120):
        return c / a
    return None


def analyse(cid, G, ev_cache):
    small = gemmi.read_small_structure(str(CACHE / f"{cid}.cif")) if (CACHE / f"{cid}.cif").exists() else None
    if small is None:
        fetch(f"https://www.crystallography.net/cod/{cid}.cif", f"{cid}.cif")
        small = gemmi.read_small_structure(str(CACHE / f"{cid}.cif"))
    if not same_ops(small, G):
        return {"id": cid, "skip": "setting"}
    c = metric_c(small.cell, G.system)
    if c is None:
        return {"id": cid, "skip": "metric"}
    key = round(c, 6)
    if key not in ev_cache:
        ev_cache.clear()
        ev_cache[key] = Evaluator(G, c)
    E = ev_cache[key]
    sites = []
    for s in small.sites:
        if s.element.name in ("H", "D"):
            continue
        x = np.array([s.fract.x, s.fract.y, s.fract.z]) % 1.0
        try:
            nf = E.facets(x)[0]
        except (ArithmeticError, MemoryError, ValueError) as e:
            nf = None
        sites.append({"label": s.label, "el": s.element.name, "x": [round(float(v), 5) for v in x], "facets": nf,
                      "occ": round(float(s.occ), 3)})
    good = [t["facets"] for t in sites if t["facets"]]
    try:
        formula = gemmi.cif.read(str(CACHE / f"{cid}.cif")).sole_block().find_value("_chemical_formula_sum")
        formula = gemmi.cif.as_string(formula) if formula else None
    except Exception:
        formula = None
    return {"id": cid, "formula": formula,
            "c_over_a": round(c, 5), "max": max(good) if good else None, "sites": sites}


if __name__ == "__main__":
    sg = sys.argv[1]
    G = Group(sg)
    ids = entries_in_group(sg)
    print(f"{sg}: {len(ids)} entries", flush=True)
    out, ev_cache = [], {}
    for i, cid in enumerate(ids):
        try:
            out.append(analyse(cid, G, ev_cache))
        except Exception as e:
            out.append({"id": cid, "error": str(e)[:100]})
        if (i + 1) % 50 == 0:
            print(f"  {i + 1}/{len(ids)}", flush=True)
    tag = sg.replace(" ", "")
    (HERE / "data" / f"cod_{tag}.json").write_text(json.dumps(out))
    from collections import Counter
    ok = [r for r in out if r.get("max")]
    print("entries with cells:", len(ok), " skipped:", Counter(r.get("skip") for r in out if "skip" in r),
          " errors:", sum("error" in r for r in out))
    print("max facets per entry:", sorted(Counter(r["max"] for r in ok).items()))
    print("top:", [(r["id"], r["max"]) for r in sorted(ok, key=lambda r: -r["max"])[:10]], flush=True)
