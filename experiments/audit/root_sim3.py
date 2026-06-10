from fractions import Fraction
I128 = 2**127

def canon(raw):
    d = {}
    for e, c in raw:
        d[e] = d.get(e, Fraction(0)) + c
    return sorted(((e, c) for e, c in d.items() if c != 0), key=lambda p: -p[0])
def mul(a, b): return canon([(ea+eb, ca*cb) for ea, ca in a for eb, cb in b])
def add(a, b): return canon(a + b)
def neg(a): return [(e, -c) for e, c in a]
def sub(a, b): return add(a, neg(b))
def trunc(a, n): return a[:n] if len(a) > n else a
ONE = [(0, Fraction(1))]

def big(x):
    return any(abs(c.numerator) >= I128 or c.denominator >= I128 for _, c in x)

def nth_root(x, k, n):
    e0, c0 = x[0]
    root_m = [(Fraction(e0, k), Fraction(1))]
    m_inv = [(-e0, 1/c0)]
    r = sub(mul(m_inv, x), ONE)
    if not r: return root_m, None
    alpha = Fraction(1, k)
    w = 2*n + 8
    series = ONE[:]; power = ONE[:]; coeff = Fraction(1)
    ovf = None
    for j in range(1, 4*w + 16 + 1):
        coeff = coeff * (alpha - (j-1)) / j
        power = trunc(mul(power, r), w)
        if ovf is None and (big(power) or abs(coeff.numerator) >= I128 or coeff.denominator >= I128):
            ovf = j
        if not power: break
        if coeff == 0: continue
        contrib = trunc(mul([(0, coeff)], power), w)
        if ovf is None and big(contrib): ovf = j
        if len(series) >= w and contrib[0][0] < series[w-1][0]: break
        series = trunc(add(series, contrib), w)
    return trunc(mul(root_m, series), n), ovf

for gap in [22, 25, 28, 30, 35, 40, 50]:
    y = canon([(0, Fraction(1)), (-1, Fraction(-1)), (-gap, Fraction(1))])
    z = mul(y, y)
    got, ovf = nth_root(z, 2, 3)
    want = trunc(y, 3)
    print(f"sqrt gap={gap}: match={got==want} ovf_j={ovf} got={[(int(e),str(c)) for e,c in got]}")
