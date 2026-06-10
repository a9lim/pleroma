from fractions import Fraction
I128 = 2**127

def canon(raw):
    d = {}
    for e, c in raw:
        d[e] = d.get(e, Fraction(0)) + c
    return sorted(((e, c) for e, c in d.items() if c != 0), key=lambda p: -p[0])
def add(a, b): return canon(a + b)
def neg(a): return [(e, -c) for e, c in a]
def sub(a, b): return add(a, neg(b))
def mul(a, b): return canon([(ea+eb, ca*cb) for ea, ca in a for eb, cb in b])
def trunc(a, n): return a[:n] if len(a) > n else a
ONE = [(0, Fraction(1))]

def nth_root(x, k, n, log=False):
    e0, c0 = x[0]
    root_m = [(Fraction(e0, k), Fraction(1))]
    m_inv = [(-e0, 1/c0)]
    r = sub(mul(m_inv, x), ONE)
    if not r: return root_m, None
    alpha = Fraction(1, k)
    w = 2*n + 8
    series = ONE[:]; power = ONE[:]; coeff = Fraction(1)
    overflow_at = None
    for j in range(1, 4*w + 16 + 1):
        coeff = coeff * (alpha - (j-1)) / j
        if overflow_at is None and (abs(coeff.numerator) >= I128 or coeff.denominator >= I128):
            overflow_at = j
        power = trunc(mul(power, r), w)
        if not power: break
        if coeff == 0: continue
        contrib = trunc(mul([(0, coeff)], power), w)
        if len(series) >= w and contrib[0][0] < series[w-1][0]:
            if log: print(f"   break at j={j}")
            break
        series = trunc(add(series, contrib), w)
    else:
        if log: print(f"   ran all {4*w+16} iters")
    return trunc(mul(root_m, series), n), overflow_at

for gap in [14, 15, 16, 18, 21]:
    y = canon([(0, Fraction(1)), (-1, Fraction(-1)), (-gap, Fraction(1))])
    z = mul(mul(y, y), y)
    got, ovf = nth_root(z, 3, 3, log=False)
    want = trunc(y, 3)
    print(f"cbrt gap={gap}: match={got==want} overflow_j={ovf} got={[(int(e),str(c)) for e,c in got]}")

# sqrt variants: y = 1 - t + c*t^gap with coeff that breaks symmetry
for gap in [14, 16, 18, 21]:
    y = canon([(0, Fraction(1)), (-1, Fraction(-1)), (-gap, Fraction(1))])
    z = mul(y, y)
    got, ovf = nth_root(z, 2, 3)
    want = trunc(y, 3)
    print(f"sqrt gap={gap}: match={got==want} overflow_j={ovf} got={[(int(e),str(c)) for e,c in got]}")
