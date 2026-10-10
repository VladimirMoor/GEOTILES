"""Ячейки Дирихле–Вороного центров белковых молекул в кристаллах PDB.

Для записи PDB берём ячейку (a, b, c, α, β, γ), пространственную группу и Cα-атомы модели 1 (RCSB ModelServer),
считаем центр масс асимметричной единицы, переводим во фракционные координаты и строим DV-ячейку орбиты этой
точки под пространственной группой (dvfast из research/p4). Число граней — число соседних молекул по граням.
Данные кэшируются в data/cache/.
"""
import json
import math
import pathlib
import sys
import time
import urllib.request
import urllib.parse

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "p4"))
from dvfast import Group, Evaluator  # noqa: E402

CACHE = pathlib.Path(__file__).resolve().parent / "data" / "cache"
CACHE.mkdir(parents=True, exist_ok=True)
UA = {"User-Agent": "GEOTILES-research/1.0 (academic; github.com/VladimirMoor/GEOTILES)"}


def fetch(url, name, binary=False):
    f = CACHE / name
    if f.exists():
        return f.read_bytes() if binary else f.read_text()
    for attempt in range(4):
        try:
            req = urllib.request.Request(url, headers=UA)
            data = urllib.request.urlopen(req, timeout=60).read()
            f.write_bytes(data)
            time.sleep(0.15)
            return data if binary else data.decode()
        except Exception:
            time.sleep(2 + 3 * attempt)
    raise RuntimeError(f"fetch failed: {url}")


def entries_in_group(sg):
    q = {"query": {"type": "group", "logical_operator": "and", "nodes": [
        {"type": "terminal", "service": "text", "parameters": {"attribute": "symmetry.space_group_name_H_M", "operator": "exact_match", "value": sg}},
        {"type": "terminal", "service": "text", "parameters": {"attribute": "exptl.method", "operator": "exact_match", "value": "X-RAY DIFFRACTION"}}]},
        "request_options": {"paginate": {"start": 0, "rows": 10000}}, "return_type": "entry"}
    req = urllib.request.Request("https://search.rcsb.org/rcsbsearch/v2/query", data=json.dumps(q).encode(),
                                 headers={**UA, "Content-Type": "application/json"})
    d = json.loads(urllib.request.urlopen(req, timeout=120).read())
    return [r["identifier"] for r in d.get("result_set", [])]


def cell_info(pid):
    q = '{entry(entry_id:"%s"){cell{length_a length_b length_c angle_alpha angle_beta angle_gamma} symmetry{space_group_name_H_M} rcsb_entry_info{resolution_combined deposited_polymer_monomer_count polymer_entity_count_protein}}}' % pid
    url = "https://data.rcsb.org/graphql?query=" + urllib.parse.quote(q)
    return json.loads(fetch(url, f"{pid}_cell.json"))["data"]["entry"]


def ca_centroid(pid):
    url = f"https://models.rcsb.org/v1/{pid}/atoms?label_atom_id=CA&model_nums=1&encoding=cif"
    txt = fetch(url, f"{pid}_ca.cif")
    xyz = []
    cols, in_loop = [], False
    for line in txt.splitlines():
        if line.startswith("_atom_site."):
            cols.append(line.split(".")[1].strip()); in_loop = True; continue
        if in_loop and line and not line.startswith(("_", "#", "loop_")):
            parts = line.split()
            if len(parts) == len(cols):
                d = dict(zip(cols, parts))
                xyz.append([float(d["Cartn_x"]), float(d["Cartn_y"]), float(d["Cartn_z"])])
        elif in_loop and line.startswith("#"):
            in_loop = False
    xyz = np.array(xyz)
    return xyz.mean(axis=0), len(xyz)


def frac_matrix(a, b, c, al, be, ga):
    """PDB-ортогонализация: a ∥ x, b в плоскости xy. Возвращает M: frac = M · cart."""
    al, be, ga = (math.radians(v) for v in (al, be, ga))
    v = math.sqrt(1 - math.cos(al) ** 2 - math.cos(be) ** 2 - math.cos(ga) ** 2 + 2 * math.cos(al) * math.cos(be) * math.cos(ga))
    O = np.array([[a, b * math.cos(ga), c * math.cos(be)],
                  [0, b * math.sin(ga), c * (math.cos(al) - math.cos(be) * math.cos(ga)) / math.sin(ga)],
                  [0, 0, c * v / math.sin(ga)]])
    return np.linalg.inv(O)


def analyse(pid):
    info = cell_info(pid)
    cl = info["cell"]
    sg = info["symmetry"]["space_group_name_H_M"]
    cen, n = ca_centroid(pid)
    x = frac_matrix(cl["length_a"], cl["length_b"], cl["length_c"], cl["angle_alpha"], cl["angle_beta"], cl["angle_gamma"]) @ cen
    G = Group(sg)
    c_over_a = cl["length_c"] / cl["length_a"]
    E = Evaluator(G, 1.0 if G.system == "cubic" else c_over_a)
    nf, gap, nv = E.facets(x % 1.0)
    return {"id": pid, "sg": sg, "facets": nf, "c_over_a": round(c_over_a, 4), "x": [round(float(v), 5) for v in x % 1.0],
            "n_ca": n, "resolution": (info.get("rcsb_entry_info") or {}).get("resolution_combined")}


if __name__ == "__main__":
    sg = sys.argv[1]
    ids = entries_in_group(sg)
    print(f"{sg}: {len(ids)} entries", flush=True)
    out = []
    for i, pid in enumerate(ids):
        try:
            r = analyse(pid)
            out.append(r)
        except Exception as e:
            out.append({"id": pid, "error": str(e)[:80]})
        if (i + 1) % 25 == 0:
            print(f"  {i + 1}/{len(ids)}", flush=True)
    tag = sg.replace(" ", "")
    pathlib.Path(f"data/dv_{tag}.json").write_text(json.dumps(out, indent=0))
    ok = [r for r in out if "facets" in r]
    from collections import Counter
    print("facet histogram:", sorted(Counter(r["facets"] for r in ok).items()))
    print("top:", [(r["id"], r["facets"]) for r in sorted(ok, key=lambda r: -r["facets"])[:8]])
    print("errors:", sum(1 for r in out if "error" in r))
