from fractions import Fraction
import random

# Surreal with integer exponents (powers of t = omega^-1): list[(exp, coeff)] desc by exp, coeff != 0
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

def inv_to_terms(x, n):
    if not x: return None
    if n == 0: return []
    e0, c0 = x[0]
    m_inv = [(-e0, 1/c0)]
    r = sub(mul(m_inv, x), ONE)
    if not r: return m_inv
    neg_r = neg(r)
    w = 2*n + 8
    series = ONE[:]
    power = ONE[:]
    for _ in range(4*w + 16):
        power = trunc(mul(power, neg_r), w)
        if not power: break
        if len(series) >= w and power[0][0] < series[w-1][0]: break
        series = trunc(add(series, power), w)
    return trunc(mul(m_inv, series), n)

# exact inverse first n nonzero terms, via long division on power series in t (exp = -t-degree)
def exact_inv_terms(x, n, depth=400):
    # x as poly in t: coeff of t^j is coeff at exp -j ; require x[0] exp == 0 (normalize first)
    e0, c0 = x[0]
    # shift so leading exp is 0: y = x * omega^{-e0}; 1/x = (1/y) * omega^{-e0}
    y = {-(e - e0): c for e, c in x}  # t-degree -> coeff, degrees >= 0
    inv = {}
    # series inversion: inv[0] = 1/y[0]; inv[j] = -(1/y0) * sum_{i=1..j} y[i]*inv[j-i]
    y0 = y[0]
    inv[0] = 1/y0
    for j in range(1, depth):
        s = sum(y.get(i, Fraction(0)) * inv[j-i] for i in range(1, j+1))
        inv[j] = -s / y0
    terms = [(-j - e0, c) for j, c in sorted(inv.items()) if c != 0]
    terms.sort(key=lambda p: -p[0])
    return terms[:n]

# Test 1: geometric x = 1 + t + ... + t^20 ; truth: 1/x = 1 - t + t^21 - t^22 + ...
x = [(-j, Fraction(1)) for j in range(21)]
got = inv_to_terms(x, 3)
want = exact_inv_terms(x, 3)
print("geometric m=20, n=3:")
print("  got :", got)
print("  want:", want)
print("  match:", got == want)

# Test 2: longer geometric
x = [(-j, Fraction(1)) for j in range(41)]
got = inv_to_terms(x, 3); want = exact_inv_terms(x, 3)
print("geometric m=40, n=3:", "match:" , got == want, "got", got, "want", want)

# Test 3: random dense polys, compare
random.seed(0)
bad = 0
for trial in range(300):
    deg = random.randint(1, 25)
    n = random.randint(1, 6)
    terms = [(0, Fraction(random.choice([1,2,3,-1,-2])))]
    for j in range(1, deg+1):
        c = random.randint(-3, 3)
        if c: terms.append((-j, Fraction(c)))
    x = canon(terms)
    got = inv_to_terms(x, n)
    want = exact_inv_terms(x, n)
    if got != want:
        bad += 1
        if bad <= 5:
            print(f"MISMATCH trial={trial} n={n} x={x}")
            print("   got :", got)
            print("   want:", want)
print("random trials mismatches:", bad, "/300")
