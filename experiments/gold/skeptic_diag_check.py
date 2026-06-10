#!/usr/bin/env python3
"""Adversarial re-verification of the gold-diagonal attack, independent code.

A. cross-check nim-mul against the mex/Turning-Corners DEFINITION on F_16
B. Fermat sesquimultiple facts
C. closed form lambda_1^(m) = XOR_{t=1..k-1} 2^(2^t-1), m up to 64 (PREDICTIVE:
   probe stopped at 32), verified as the unique trace dual directly
D. C(M) identity Tr_M(2^(M-1)(e_j^3 + e_j)) = 0, M up to 32
E. Fermat-coin witness q_{m/2} = a mod 2, m up to 64
F. COUNTEREXAMPLE CANDIDATE: lambda_1^(m) == w*w (+) w, w = XOR of coins at
   power-of-two indices >= 2 (m-independent index predicate + disjunctive sum
   + ONE game-squaring) -- a commutative, uniform, game-built source the
   attack's "commutative-local" class excludes only via the fixed-lambda clause
G. drift asserted for printed even-a lambdas (probe printed, never asserted)
H. cocycle properties of P(u,v) = Tr(u v^(2^a)) on F_256
I. e_0 = 1 in radical of Gold B
"""
import random
from functools import lru_cache

# my own nim-mul: Conway recursion, written independently of the probe
@lru_cache(maxsize=None)
def nmul(x, y):
    if x < 2 or y < 2:
        return x * y
    lev = 1
    while (max(x, y) >> (1 << lev)) > 0:
        lev += 1
    h = 1 << (lev - 1)      # half-exponent
    H = 1 << h              # Fermat point 2^(2^(lev-1))
    a, b = x >> h, x & (H - 1)
    c, d = y >> h, y & (H - 1)
    # (aH+b)(cH+d) = ac H^2 + (ad+bc) H + bd,  H^2 = H + H/2
    ac = nmul(a, c)
    ad_bc = nmul(a, d) ^ nmul(b, c)
    bd = nmul(b, d)
    return ((ac ^ ad_bc) << h) ^ bd ^ nmul(ac, H >> 1)

# A. mex-definition cross-check on F_16 (Turning Corners recurrence, from scratch)
tbl = [[0] * 16 for _ in range(16)]
for x in range(16):
    for y in range(16):
        seen = {tbl[i][y] ^ tbl[x][j] ^ tbl[i][j] for i in range(x) for j in range(y)}
        v = 0
        while v in seen:
            v += 1
        tbl[x][y] = v
for x in range(16):
    for y in range(16):
        assert tbl[x][y] == nmul(x, y), (x, y)
print("A. nim-mul == mex/Turning-Corners definition on all of F_16: OK")

# B. sesquimultiples
for t in range(1, 6):
    F = 1 << (1 << t)
    assert nmul(F, F) == F ^ (F >> 1)
print("B. Fermat sesquimultiple F^2 = F (+) F/2 for t=1..5: OK")

def sq(x):
    return nmul(x, x)

def frob(x, a):
    for _ in range(a):
        x = sq(x)
    return x

def tr(x, m):
    acc, t = x, x
    for _ in range(m - 1):
        t = sq(t)
        acc ^= t
    assert acc < 2
    return acc

def gold_q(i, a, m):
    e = 1 << i
    return tr(nmul(e, frob(e, a)), m)

# C. closed form as unique trace dual, m = 4..64 (m=64 is OUT-OF-SAMPLE for the probe)
print("C. closed form lambda_1^(m) = XOR_{t=1..k-1} 2^(2^t - 1):")
for m in (4, 8, 16, 32, 64):
    k = m.bit_length() - 1
    lam = 0
    for t in range(1, k):
        lam ^= 1 << ((1 << t) - 1)
    ok = all(tr(nmul(lam, 1 << i), m) == gold_q(i, 1, m) for i in range(m))
    print(f"   m={m:>2}: lambda={lam}  is-the-trace-dual: {ok}")
    assert ok
print("   (trace pairing nondegenerate => this IS the unique lambda_1^(m))")

# D. C(M): Tr_M(2^(M-1) (e_j^3 + e_j)) = 0 for all j < M
for M in (4, 8, 16, 32):
    u = 1 << (M - 1)
    ok = all(tr(nmul(u, nmul(1 << j, sq(1 << j)) ^ (1 << j)), M) == 0 for j in range(M))
    print(f"D. C({M}): {ok}")
    assert ok

# E. Fermat-coin witness
for m in (4, 8, 16, 32, 64):
    k = m.bit_length() - 1
    for a in range(1, k + 1):
        assert gold_q(m // 2, a, m) == a % 2, (m, a)
print("E. q_(m/2)^(m,a) = a mod 2 for all m<=64, 1<=a<=log2(m): OK")

# F. the counterexample candidate: lambda_1 = w^2 (+) w, w = XOR of Fermat coins
print("F. commutative game-built source for the a=1 diagonal:")
for m in (4, 8, 16, 32, 64):
    w = 0
    i = 2
    while i < m:                 # m-independent index predicate: i a 2-power >= 2
        w ^= 1 << i
        i <<= 1
    lam = sq(w) ^ w              # one game-squaring + disjunctive sum
    ok = all(tr(nmul(lam, 1 << i_), m) == gold_q(i_, 1, m) for i_ in range(m))
    print(f"   m={m:>2}: w={w}  wp(w)=w^2(+)w={lam}  equals lambda_1^(m): {ok}")
    assert ok

# G. assert even-a drift for the probe's printed duals (verify each IS the dual)
printed = {(4, 2): 0, (8, 2): 6, (16, 2): 102, (32, 2): 24582,
           (16, 4): 31, (32, 4): 8030}
for (m, a), lam in printed.items():
    assert all(tr(nmul(lam, 1 << i), m) == gold_q(i, a, m) for i in range(m)), (m, a)
assert len({printed[(m, 2)] for m in (4, 8, 16, 32)}) == 4
assert printed[(16, 4)] != printed[(32, 4)]
print("G. even-a duals verified, drift (all distinct) ASSERTED: OK")

# H. cocycle properties of P(u,v) = Tr(u v^(2^a)) on F_256, a=1
m, a = 8, 1
rng = random.Random(9)
def P(u, v):
    return tr(nmul(u, frob(v, a)), m)
for _ in range(300):
    u, v, w2 = (rng.randrange(256) for _ in range(3))
    assert P(u ^ w2, v) == P(u, v) ^ P(w2, v)
    assert P(u, v ^ w2) == P(u, v) ^ P(u, w2)
    assert P(u, u) == tr(nmul(u, frob(u, a)), m)
    Bv = P(u, v) ^ P(v, u)
    assert Bv == (P(u ^ v, u ^ v) ^ P(u, u) ^ P(v, v))
print("H. P bilinear, P(u,u)=Q, P+P^T=B (polarization) on F_256: OK")

# I. e_0 = 1 lies in radical of every Gold polar form (tested m<=32)
for m2 in (4, 8, 16, 32):
    k = m2.bit_length() - 1
    for a2 in range(1, k + 1):
        for _ in range(40):
            v = rng.randrange(1 << m2)
            Bv = tr(frob(v, a2), m2) ^ tr(v, m2)   # B(1,v) = Tr(v^{2^a}) + Tr(v)
            assert Bv == 0
print("I. 1 in R(B) for all Gold polar forms tested: OK")
print("\nALL SKEPTIC CHECKS PASSED")
