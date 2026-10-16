"""Выгрузка точных систем (exact.build, с переменными Рабиновича) для непроверенных типов в JSON для Julia."""
import json
import pathlib

from exact import build

ROOT = pathlib.Path(__file__).resolve().parent
status = json.loads((ROOT / "data" / "rigor_status.json").read_text())
types = {t["id"]: t for t in json.loads((ROOT / "data" / "types.json").read_text())}
out = []
for group in ("uncertified_timeout", "uncertified_positive_dim"):
    for tid in status[group]:
        variables, eqs = build(types[tid], rabinowitsch=True)
        v0, e0 = build(types[tid], rabinowitsch=False)
        out.append({"id": tid, "group": group, "variables": variables, "eqs": eqs, "variables_plain": v0, "eqs_plain": e0,
                    "nF": types[tid]["nF"], "nV": types[tid]["nV"], "faces": types[tid]["faces"]})
(ROOT / "data" / "hc" / "systems.json").write_text(json.dumps(out))
print(len(out), "systems; rabinowitsch: max vars", max(len(o["variables"]) for o in out), "| plain: max vars",
      max(len(o["variables_plain"]) for o in out), "max eqs", max(len(o["eqs_plain"]) for o in out))
