import random, itertools
from math import gcd
from fractions import Fraction

def div_euclid(a, b):
    q, r = divmod(a, b)
    return q + 1 if r < 0 else q

# ---- faithful port of smith_normal_form (integer.rs) ----
def ext_gcd(a, b):
    r0, r1 = a, b
    s0, s1 = 1, 0
    t0, t1 = 0, 1
    while r1 != 0:
        q = div_euclid(r0, r1)
        r0, r1 = r1, r0 - q*r1
        s0, s1 = s1, s0 - q*s1
        t0, t1 = t1, t0 - q*t1
    if r0 < 0:
        return -r0, -s0, -t0
    return r0, s0, t0

def trunc_div(a, b):  # Rust `/` is truncated
    q = abs(a)//abs(b)
    return q if (a<0)==(b<0) else -q

def smith_normal_form(m):
    m = [row[:] for row in m]
    rows = len(m)
    if rows == 0: return []
    cols = len(m[0])
    k = min(rows, cols)
    for t in range(k):
        while True:
            if m[t][t] == 0:
                piv = None
                for i in range(t, rows):
                    for j in range(t, cols):
                        if m[i][j] != 0:
                            piv = (i, j); break
                    if piv: break
                if piv is None:
                    break
                i, j = piv
                m[t], m[i] = m[i], m[t]
                for row in m:
                    row[t], row[j] = row[j], row[t]
            changed = False
            for i in range(t+1, rows):
                if m[i][t] == 0: continue
                if m[i][t] % m[t][t] == 0:
                    q = trunc_div(m[i][t], m[t][t])
                    for c in range(cols):
                        m[i][c] -= q*m[t][c]
                else:
                    g, x, y = ext_gcd(m[t][t], m[i][t])
                    u = trunc_div(-m[i][t], g); v = trunc_div(m[t][t], g)
                    for c in range(cols):
                        a0, b0 = m[t][c], m[i][c]
                        m[t][c] = x*a0 + y*b0
                        m[i][c] = u*a0 + v*b0
                    changed = True
            if changed: continue
            for j in range(t+1, cols):
                if m[t][j] == 0: continue
                if m[t][j] % m[t][t] == 0:
                    q = trunc_div(m[t][j], m[t][t])
                    for r in range(rows):
                        m[r][j] -= q*m[r][t]
                else:
                    g, x, y = ext_gcd(m[t][t], m[t][j])
                    u = trunc_div(-m[t][j], g); v = trunc_div(m[t][t], g)
                    for row in m:
                        a0, b0 = row[t], row[j]
                        row[t] = x*a0 + y*b0
                        row[j] = u*a0 + v*b0
                    changed = True
            if changed: continue
            p = m[t][t]
            violated = None
            for i in range(t+1, rows):
                for j in range(t+1, cols):
                    if m[i][j] % p != 0:
                        violated = i; break
                if violated is not None: break
            if violated is not None:
                for c in range(cols):
                    m[t][c] += m[violated][c]
            else:
                break
    return [abs(m[i][i]) for i in range(k)]

# ---- independent characterization: determinantal divisors ----
def minor_det(m, rsel, csel):
    sub = [[Fraction(m[r][c]) for c in csel] for r in rsel]
    n = len(sub)
    det = Fraction(1)
    for col in range(n):
        piv = next((r for r in range(col, n) if sub[r][col] != 0), None)
        if piv is None: return 0
        if piv != col:
            sub[col], sub[piv] = sub[piv], sub[col]
            det = -det
        det *= sub[col][col]
        inv = 1/sub[col][col]
        for r in range(col+1, n):
            f = sub[r][col]*inv
            if f:
                for c in range(col, n):
                    sub[r][c] -= f*sub[col][c]
    assert det.denominator == 1
    return int(det)

def determinantal_divisors_check(m, d):
    rows, cols = len(m), len(m[0])
    k = min(rows, cols)
    prod = 1
    for j in range(1, k+1):
        g = 0
        for rsel in itertools.combinations(range(rows), j):
            for csel in itertools.combinations(range(cols), j):
                g = gcd(g, abs(minor_det(m, rsel, csel)))
        # gcd of all j x j minors must equal d1*...*dj
        expect = prod * d[j-1]
        assert g == expect, (m, d, j, g, expect)
        prod = expect
        if prod == 0:
            # all larger minors must be 0 too; chain of zeros
            for jj in range(j, k):
                assert d[jj] == 0, (m, d)
            break

random.seed(20260609)
fails = 0
for trial in range(4000):
    rows = random.randint(1, 4)
    cols = random.randint(1, 4)
    m = [[random.randint(-9, 9) for _ in range(cols)] for _ in range(rows)]
    d = smith_normal_form(m)
    # chain
    for a, b in zip(d, d[1:]):
        if a == 0:
            assert b == 0, (m, d)
        else:
            assert b % a == 0, (m, d)
    determinantal_divisors_check(m, d)
print("SNF: 4000 random matrices OK (chain + determinantal divisors)")
