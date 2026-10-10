"""Независимая проверка сертификатов из data/certificates.json (точная рациональная арифметика)."""
import json
from fractions import Fraction as Fr
from exact import certify

for n, r in sorted(json.load(open("data/certificates.json")).items(), key=lambda kv: int(kv[0])):
    x = tuple(Fr(v) for v in r["x"]); c = Fr(r["c_over_a"])
    nf, nv, npl = certify(r["group"], x, c * c)
    status = "OK" if nf >= r["facets"] else "FAIL"
    print(f"IT({int(n):3d}) {r['group']:12s} claimed {r['facets']:2d}, certified {nf:2d}  [{status}]")
