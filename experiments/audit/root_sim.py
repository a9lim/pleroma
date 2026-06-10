from fractions import Fraction

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

def binomial_series(r, alpha, n):
    w = 2*n + 8
    series = ONE[:]
    power = ONE[:]
    coeff = Fraction(1)
    for j in range(1, 4*w + 16 + 1):
        coeff = coeff * (alpha - (j-1)) / j
        power = trunc(mul(power, r), w)
        if not power: break
        if coeff == 0: continue
        contrib = trunc(mul([(0, coeff)], power), w)
        if len(series) >= w and contrib[0][0] < series[w-1][0]: break
        series = trunc(add(series, contrib), w)
    return series

def nth_root_to_terms(x, k, n):
    if not x: return []
    e0, c0 = x[0]
    # leading coeff must be perfect k-th power; here it's 1
    root_m = [(Fraction(e0, k), Fraction(1))]
    m_inv = [(-e0, 1/c0)]
    r = sub(mul(m_inv, x), ONE)
    if not r: return root_m
    series = binomial_series(r, Fraction(1, k), n)
    return trunc(mul(root_m, series), n)

# y = 1 - t + t^21  (t = omega^{-1});  z = y^2 finite support; sqrt(z) = y
y = [(0, Fraction(1)), (-1, Fraction(-1)), (-21, Fraction(1))]
z = mul(y, y)
print("z =", z)
got = nth_root_to_terms(z, 2, 3)
print("sqrt_to_terms(z,3) got :", got)
print("                  want :", trunc(y,3))
print("match:", got == trunc(y,3))

# also cube: y3 = 1 - t + t^21, z3 = y^3, cube root
z3 = mul(mul(y, y), y)
got3 = nth_root_to_terms(z3, 3, 3)
print("cbrt(y^3, 3) got:", got3, "match:", got3 == trunc(y,3))
