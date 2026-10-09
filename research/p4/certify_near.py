"""Точная проверка найденных точек: рациональные приближения разной точности."""
import sys
from fractions import Fraction as Fr
from exact import certify

def run(hm, x, c):
    best = None
    for den in (10**4, 10**5, 10**6, 10**7, 10**8):
        xr = tuple(Fr(v).limit_denominator(den) for v in x)
        c2 = Fr(c).limit_denominator(den) ** 2
        try:
            nf, nv, npl = certify(hm, xr, c2)
        except Exception as e:
            nf = None
        print(f"  den {den}: certified facets {nf}  x={tuple(map(str, xr))}, c/a={Fr(c).limit_denominator(den)}", flush=True)
        if nf is not None and (best is None or nf > best[0]):
            best = (nf, xr, Fr(c).limit_denominator(den))
    return best

if __name__ == "__main__":
    hm = sys.argv[1]; vals = [float(v) for v in sys.argv[2:]]
    print(run(hm, vals[:3], vals[3] if len(vals) > 3 else 1.0))
