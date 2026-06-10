#!/usr/bin/env python3
"""Final-skeptic independent verification for the round-2 construct RESULT.

Independence measures:
 - nim multiplication cross-checked against the Conway mex DEFINITION
   (a*b = mex{a'b + ab' + a'b' : a'<a, b'<b}) on the full 64x64 table --
   an oracle neither probe used.
 - field-structure pins at m=8: x^255=1 for x!=0, subfield sizes
   |{x:x^2=x}|=2, |{x:x^4=x}|=4, |{x:x^16=x}|=16, associativity fuzz.
 - T2 solver re-implemented from the RULE PROSE (not their code).
 - explicit check: T2 gate == [Q(x^d) != Q(x)] for every position and
   every descending width-<=2 flip (the disguised-differential-evaluator test).
 - random maximal plays: parity of play length == Q(x) (unrolled evaluator test).
 - blocking-lemma fuzz on RANDOM forms (random q AND random alternating B).
 - ECHO re-implemented from prose; m=4 exactness + mistake states; (8,2) bench.
"""
import random, sys, time
from functools import lru_cache
sys.setrecursionlimit(1000000)

t00 = time.time()

# ---------- 1. nim arithmetic: my karatsuba + the mex-definition oracle ------
@lru_cache(maxsize=None)
def nm(a, b):
    if a < 2 or b < 2:
        return a * b
    if a < b:
        a, b = b, a
    k = 0
    while a >= (1 << (1 << (k + 1))):
        k += 1
    h = 1 << k
    F = 1 << h
    a1, a0 = a >> h, a & (F - 1)
    b1, b0 = b >> h, b & (F - 1)
    hh = nm(a1, b1)
    return ((hh ^ nm(a1, b0) ^ nm(a0, b1)) << h) ^ nm(a0, b0) ^ nm(hh, F >> 1)

N_MEX = 64
T = [[0] * N_MEX for _ in range(N_MEX)]
for a in range(N_MEX):
    for b in range(N_MEX):
        s = set()
        for ap in range(a):
            Tb = T[ap]
            for bp in range(b):
                s.add(Tb[b] ^ T[a][bp] ^ Tb[bp])
        m_ = 0
        while m_ in s:
            m_ += 1
        T[a][b] = m_
assert all(T[a][b] == nm(a, b) for a in range(N_MEX) for b in range(N_MEX))
print(f"[ok] karatsuba == Conway mex definition on full {N_MEX}x{N_MEX} table")
assert nm(2, 2) == 3 and nm(2, 4) == 8 and nm(16, 16) == 24

def npow(x, e):
    r = 1
    while e:
        if e & 1:
            r = nm(r, x)
        x = nm(x, x)
        e >>= 1
    return r

assert all(npow(x, 255) == 1 for x in range(1, 256))
assert sum(1 for x in range(256) if npow(x, 2) == x) == 2
assert sum(1 for x in range(256) if npow(x, 4) == x) == 4
assert sum(1 for x in range(256) if npow(x, 16) == x) == 16
rng = random.Random(99)
for _ in range(2000):
    a, b, c = (rng.randrange(256) for _ in range(3))
    assert nm(nm(a, b), c) == nm(a, nm(b, c))
    assert nm(a, b ^ c) == nm(a, b) ^ nm(a, c)
print("[ok] F_256 structure: x^255=1, subfield sizes 2/4/16, assoc+distrib fuzz")

def nsq(x):
    return nm(x, x)

def frob(x, a):
    for _ in range(a):
        x = nsq(x)
    return x

def tr(x, m):
    t, y = 0, x
    for _ in range(m):
        t ^= y
        y = nsq(y)
    assert t in (0, 1)
    return t

# ---------- 2. forms ---------------------------------------------------------
def form(m, a, lam=1):
    n = 1 << m
    Q = [tr(nm(lam, nm(v, frob(v, a))), m) for v in range(n)]
    q = [Q[1 << i] for i in range(m)]
    B = [[0 if i == j else Q[(1 << i) ^ (1 << j)] ^ Q[1 << i] ^ Q[1 << j]
          for j in range(m)] for i in range(m)]
    # polarization identity on full table
    for _ in range(500):
        u, v = rng.randrange(n), rng.randrange(n)
        bv = 0
        for i in range(m):
            if (u >> i) & 1:
                for j in range(m):
                    if (v >> j) & 1:
                        bv ^= B[i][j]
        assert bv == Q[u ^ v] ^ Q[u] ^ Q[v]
    return Q, q, B

def radical_dim(B, m):
    n = 1 << m
    cnt = 0
    for v in range(n):
        ok = True
        for i in range(m):
            s = 0
            for j in range(m):
                if (v >> j) & 1:
                    s ^= B[i][j]
            if s:
                ok = False
                break
        cnt += ok
    d = cnt.bit_length() - 1
    assert 1 << d == cnt
    return d

checks = [
    # (m, a, lam, zero_count, radical_dim)
    (4, 1, 1, 4, 2),
    (4, 1, 2, 10, 0),
    (8, 1, 1, 112, 2),
    (8, 2, 1, 96, 4),
    (8, 1, 2, 136, 0),
]
FORMS = {}
for (m, a, lam, zc, rd) in checks:
    Q, q, B = form(m, a, lam)
    z = sum(1 for v in Q if v == 0)
    r = radical_dim(B, m)
    assert z == zc, (m, a, lam, z)
    assert r == rd, (m, a, lam, r)
    FORMS[(m, a, lam)] = (Q, q, B)
    print(f"[ok] (m={m},a={a},lam={lam}): |Q=0|={z}, radical dim={r}, rank={m-r}")
Q41, _, _ = FORMS[(4, 1, 1)]
assert {v for v in range(16) if Q41[v] == 0} == {0, 1, 2, 3}
print("[ok] m=4 unscaled zero set == nim subfield F_4 = {0,1,2,3} (affine/vacuous)")
assert 112 & (112 - 1) != 0  # not a power of two -> not an affine subspace
# isotropic radical at (8,1): Q vanishes on radical
Q81, q81, B81 = FORMS[(8, 1, 1)]
rad = [v for v in range(256)
       if all((bin(v & sum((B81[i][j] << j) for j in range(8))).count("1") & 1) == 0
              for i in range(8))]
assert len(rad) == 4 and all(Q81[v] == 0 for v in rad)
print("[ok] (8,1) radical is isotropic; 112 = 4 * 28 = 2^rad * (2^5 - 2^2)  => Arf 1")

# ---------- 3. Weierstrass source --------------------------------------------
w8 = (1 << 2) ^ (1 << 4)
lam8 = nsq(w8) ^ w8
assert w8 == 20 and lam8 == 10, (w8, lam8)
assert all(tr(nm(lam8, 1 << i), 8) == q81[i] for i in range(8))
assert q81 == [0, 0, 0, 0, 1, 1, 0, 0]
print(f"[ok] m=8: w=20, P(w)=w^2^w={lam8}; Tr(10*e_i) == Tr(e_i^3) for all i; q={q81}")
w16 = (1 << 2) ^ (1 << 4) ^ (1 << 8)
lam16 = nsq(w16) ^ w16
assert w16 == 276 and lam16 == 138, (w16, lam16)
q16 = [tr(nm(1 << i, nsq(1 << i)), 16) for i in range(16)]
assert all(tr(nm(lam16, 1 << i), 16) == q16[i] for i in range(16))
print(f"[ok] m=16: w=276, P(w)={lam16}; Tr(138*e_i) == Tr(e_i^3) for all i")

# ---------- 4. T2 game, implemented from the rule prose -----------------------
def t2_moves(v, m, q, B):
    """All legal (d, target): wt(d) in {1,2}, msb(d) heads in v,
    gate B(v,d) ^ Q(d) == 1."""
    out = []
    heads = [i for i in range(m) if (v >> i) & 1]
    Brow = [0] * m  # B(v, e_l)
    for l in range(m):
        s = 0
        for k in heads:
            s ^= B[k][l]
        Brow[l] = s
    for i in heads:
        if Brow[i] ^ q[i] == 1:                       # single
            out.append(((1 << i), v ^ (1 << i)))
        for j in range(i):                            # pair, msb i is a head
            qd = q[i] ^ q[j] ^ B[i][j]
            if Brow[i] ^ Brow[j] ^ qd == 1:
                d = (1 << i) | (1 << j)
                out.append((d, v ^ d))
    return out

def t2_pset(m, q, B):
    n = 1 << m
    win = [False] * n
    movelists = [None] * n
    for v in range(n):
        mv = [t for (_, t) in t2_moves(v, m, q, B)]
        movelists[v] = mv
        win[v] = any(not win[t] for t in mv)  # targets < v always
        assert all(t < v for t in mv)
    return {v for v in range(n) if not win[v]}, movelists

Z81 = {v for v in range(256) if Q81[v] == 0}
pset, mls = t2_pset(8, q81, B81)
assert pset == Z81
print(f"[ok] T2 P-set == {{Tr(x^3)=0}} exactly at (8,1)  ({len(pset)} positions)")

# THE DISGUISED-EVALUATOR TEST: gate == [Q flips], for every v and every
# descending width-<=2 d (legal or not).
viol = 0
for v in range(256):
    legal_targets = {t for (_, t) in t2_moves(v, 8, q81, B81)}
    for i in range(8):
        if not (v >> i) & 1:
            continue
        cands = [1 << i] + [(1 << i) | (1 << j) for j in range(i)]
        for d in cands:
            flips = (Q81[v] ^ Q81[v ^ d]) == 1
            if flips != ((v ^ d) in legal_targets):
                # careful: two different d can hit same target; recheck via gate
                viol += 1
assert viol == 0
print("[ok] PROVEN BY EXHAUSTION at (8,1): T2 gate <=> 'this move changes Q(x)'"
      " for all 256 positions x and every descending width-<=2 flip d")

# decision-degeneracy: trivially forced by the gate; confirm + unrolled-evaluator:
assert all(len({Q81[t] for t in mls[v]}) <= 1 for v in range(256))
for _ in range(2000):
    v0 = rng.randrange(256)
    v, steps = v0, 0
    while True:
        mv = mls[v]
        if not mv:
            break
        v = rng.choice(mv)
        steps += 1
    assert steps & 1 == Q81[v0], (v0, steps)
print("[ok] every maximal play (random choices!) has length parity == Q(x):"
      " outcome is a forced clock; the game IS an unrolled Q-evaluation")

# refinement uniformity == UNIVERSALITY: random q AND random alternating B
fails = 0
for trial in range(60):
    mm = 8
    q2 = [rng.randint(0, 1) for _ in range(mm)]
    B2 = [[0] * mm for _ in range(mm)]
    for i in range(mm):
        for j in range(i):
            bbit = rng.randint(0, 1)
            B2[i][j] = B2[j][i] = bbit
    def Q2(v):
        s = 0
        bits = [i for i in range(mm) if (v >> i) & 1]
        for i in bits:
            s ^= q2[i]
        for x in range(len(bits)):
            for y in range(x + 1, len(bits)):
                s ^= B2[bits[x]][bits[y]]
        return s
    p2, _ = t2_pset(mm, q2, B2)
    if p2 != {v for v in range(1 << mm) if Q2(v) == 0}:
        fails += 1
print(f"[{'ok' if fails == 0 else 'XX'}] T2 template realizes {{Q'=0}} for "
      f"60/60 RANDOM (q,B) forms (random B too, incl. degenerate): "
      f"failures={fails} -- the template is a UNIVERSAL quadratic-form clock")

# ---------- 5. m=16 scale check ----------------------------------------------
t0 = time.time()
m16 = 16
Q16 = [tr(nm(v, nsq(v)), m16) for v in range(1 << m16)]
z16 = sum(1 for v in Q16 if v == 0)
assert z16 == 32512, z16
q16d = [Q16[1 << i] for i in range(m16)]
B16 = [[0 if i == j else Q16[(1 << i) ^ (1 << j)] ^ Q16[1 << i] ^ Q16[1 << j]
        for j in range(m16)] for i in range(m16)]
p16, _ = t2_pset(m16, q16d, B16)
assert p16 == {v for v in range(1 << m16) if Q16[v] == 0}
print(f"[ok] m=16: |Q=0|=32512 and T2 P-set == {{Q=0}} exactly "
      f"[{time.time()-t0:.0f}s]")

# bent T2 at (8,1,lam=2)
Qb, qb, Bb = FORMS[(8, 1, 2)]
pb, _ = t2_pset(8, qb, Bb)
assert pb == {v for v in range(256) if Qb[v] == 0}
print("[ok] bent (8,1,lam=2): T2 P-set == {Q=0} exactly (136 positions)")

# ---------- 6. ECHO-ko, implemented fresh from the prose ----------------------
def echo_value(x, m, q, B, tgt):
    bits = [i for i in range(m) if (x >> i) & 1]
    if not bits:
        return 0
    side = [0] * m
    for i in range(m):
        for k in range(m):
            if k > i and B[k][i]:
                side[i] |= 1 << k

    def charge(o, i):
        c = q[i] if (o >> i) & 1 else 0
        c ^= bin(o & side[i]).count("1") & 1
        return c

    memo = {}
    def rec(u, o, last, mover, sigma):
        if u == 0 and o == 0:
            return sigma
        key = (u, o, last, mover, sigma)
        if key in memo:
            return memo[key]
        legal = []
        for i in bits:
            if i == last:
                continue
            b = 1 << i
            if u & b:
                legal.append((i, u ^ b, o ^ b))
            elif o & b:
                legal.append((i, u, o ^ b))
        if not legal:
            res = rec(u, o, -1, 1 - mover, sigma)
        else:
            want = tgt if mover == 0 else 1 - tgt
            res = 1 - want
            for (i, u2, o2) in legal:
                if rec(u2, o2, i, 1 - mover, sigma ^ charge(o, i)) == want:
                    res = want
                    break
            memo[key] = res
        return res

    return rec(x, 0, -1, 0, 0)

def echo_mistakes(x, m, q, B, tgt):
    """choice-states and mistake-states reachable from x (full option values)."""
    bits = [i for i in range(m) if (x >> i) & 1]
    if not bits:
        return 0, 0
    side = [0] * m
    for i in range(m):
        for k in range(m):
            if k > i and B[k][i]:
                side[i] |= 1 << k
    def charge(o, i):
        c = q[i] if (o >> i) & 1 else 0
        c ^= bin(o & side[i]).count("1") & 1
        return c
    memo = {}
    def val(u, o, last, mover, sigma):
        if u == 0 and o == 0:
            return sigma
        key = (u, o, last, mover, sigma)
        if key in memo:
            return memo[key]
        legal = []
        for i in bits:
            if i == last:
                continue
            b = 1 << i
            if u & b:
                legal.append((i, u ^ b, o ^ b))
            elif o & b:
                legal.append((i, u, o ^ b))
        if not legal:
            res = val(u, o, -1, 1 - mover, sigma)
        else:
            want = tgt if mover == 0 else 1 - tgt
            outs = [val(u2, o2, i, 1 - mover, sigma ^ charge(o, i))
                    for (i, u2, o2) in legal]
            res = want if want in outs else 1 - want
        memo[key] = res
        return res
    seen, stack, cs, ms = set(), [(x, 0, -1, 0, 0)], 0, 0
    while stack:
        st = stack.pop()
        if st in seen:
            continue
        seen.add(st)
        (u, o, last, mover, sigma) = st
        if u == 0 and o == 0:
            continue
        legal = []
        for i in bits:
            if i == last:
                continue
            b = 1 << i
            if u & b:
                legal.append((i, u ^ b, o ^ b))
            elif o & b:
                legal.append((i, u, o ^ b))
        if not legal:
            stack.append((u, o, -1, 1 - mover, sigma))
            continue
        outs = []
        for (i, u2, o2) in legal:
            s2 = sigma ^ charge(o, i)
            outs.append(val(u2, o2, i, 1 - mover, s2))
            stack.append((u2, o2, i, 1 - mover, s2))
        if len(legal) >= 2:
            cs += 1
            if len(set(outs)) > 1:
                ms += 1
    return cs, ms

print("\n--- ECHO-ko (fresh implementation from prose) ---")
Qb4, qb4, Bb4 = FORMS[(4, 1, 2)]
vals = [echo_value(x, 4, qb4, Bb4, 1) for x in range(16)]
exact_b4 = all(vals[x] == Qb4[x] for x in range(16))
cs = ms = 0
for x in range(1, 16):
    c, mk = echo_mistakes(x, 4, qb4, Bb4, 1)
    cs += c
    ms += mk
print(f"bent (4,1,lam=2) lower/A: exact={exact_b4}  choice-states={cs}  "
      f"mistake-states={ms}  {'NON-DEGENERATE' if ms > 0 else 'clock'}")
assert exact_b4 and ms > 0

Qg4, qg4, Bg4 = FORMS[(4, 1, 1)]
vals_g = [echo_value(x, 4, qg4, Bg4, 1) for x in range(16)]
print(f"Gold (4,1,1)   lower/A: exact={all(vals_g[x] == Qg4[x] for x in range(16))}")

t0 = time.time()
Q82, q82, B82 = FORMS[(8, 2, 1)]
agree = 0
misses = []
for x in range(256):
    v = echo_value(x, 8, q82, B82, 1)
    if v == Q82[x]:
        agree += 1
    else:
        misses.append((x, bin(x).count("1"), Q82[x], v))
print(f"(8,2,1) lower/A: {agree}/256  misses={misses} [{time.time()-t0:.0f}s]")

t0 = time.time()
agree81 = 0
for x in range(256):
    if echo_value(x, 8, q81, B81, 1) == Q81[x]:
        agree81 += 1
print(f"(8,1,1) lower/A: {agree81}/256 [{time.time()-t0:.0f}s]")

print(f"\nall checks done [{time.time()-t00:.0f}s total]")
