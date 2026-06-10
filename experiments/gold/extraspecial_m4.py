"""Validate solver against an independent no-memo tree solver, then sweep m=4."""
import itertools, random
from extraspecial_core import *

validate()

# ---------- independent brute-force solver (lists, no memo, no bitmasks) ----------
def brute_value(x, Q, Brows, m, ko='self', maxfirst=True):
    S = [i for i in range(m) if (x >> i) & 1]
    k = len(S)
    if k == 0: return 0
    B = [[(Brows[S[i]] >> S[j]) & 1 for j in range(k)] for i in range(k)]
    q = [Q[1 << c] for c in S]
    def charge(t, i):
        # triangular cocycle on global coin order (S sorted ascending)
        c = q[i] if t[i] == 1 else 0
        for j in range(k):
            if t[j] == 1 and S[j] > S[i] and B[j][i]:
                c ^= 1
        return c
    def rec(t, last, mover, sigma):
        if all(v == 2 for v in t): return sigma
        legal = [i for i in range(k) if t[i] < 2 and i != last]
        if not legal:
            return rec(t, -1, 1 - mover, sigma)
        wantmax = (mover == 0) == maxfirst
        out = []
        for i in legal:
            t2 = list(t); t2[i] += 1
            out.append(rec(tuple(t2), i, 1 - mover, sigma ^ charge(t, i)))
        return max(out) if wantmax else min(out)
    return rec(tuple([0] * k), -1, 0, 0)

# cross-validate on all popcount<=3 positions of (8,1) and 30 random k=4 positions
Q81 = gold_q(8, 1); B81 = polar(Q81, 8)
mismatch = 0
tested = 0
for x in range(256):
    if bin(x).count('1') <= 3:
        for mf in (True, False):
            v1, _ = echo_value(x, Q81, B81, 8, maxfirst=mf)
            v2 = brute_value(x, Q81, B81, 8, maxfirst=mf)
            tested += 1
            if v1 != v2: mismatch += 1; print("MISMATCH", x, mf, v1, v2)
random.seed(1)
k4 = [x for x in range(256) if bin(x).count('1') == 4]
for x in random.sample(k4, 20):
    for mf in (True, False):
        v1, _ = echo_value(x, Q81, B81, 8, maxfirst=mf)
        v2 = brute_value(x, Q81, B81, 8, maxfirst=mf)
        tested += 1
        if v1 != v2: mismatch += 1; print("MISMATCH", x, mf, v1, v2)
print(f"solver cross-validation: {tested} checks, {mismatch} mismatches")

# section-independence sanity: reversed triangular convention must give same value
def echo_value_revtri(x, Q, Brows, m, maxfirst=True):
    S = [i for i in range(m) if (x >> i) & 1]
    k = len(S)
    if k == 0: return 0
    qb = [Q[1 << c] for c in S]
    Bh = []
    for li, ci in enumerate(S):
        mask = 0
        for lj, cj in enumerate(S):
            if cj < ci and ((Brows[ci] >> cj) & 1):   # reversed: k<j half
                mask |= 1 << lj
        Bh.append(mask)
    return Echo(qb, Bh).solve(maxfirst=maxfirst)[0]

ok = all(echo_value(x, Q81, B81, 8)[0] == echo_value_revtri(x, Q81, B81, 8)
         for x in range(256) if bin(x).count('1') <= 4)
print("cocycle-convention independence (k<=4, (8,1)):", "OK" if ok else "FAIL")

# ---------- m=4 sweep: unscaled Gold + all 15 lambda components ----------
print("\n=== m=4, a=1, ko=self ===")
for lam in range(1, 16):
    Q = gold_q(4, 1, lam)
    B = polar(Q, 4)
    rank_info = sum(1 for v in Q if v == 0)
    # genuinely quadratic? (skip purely affine-linear: B identically 0)
    Bzero = all(b == 0 for b in B)
    res = {}
    for mf in (True, False):
        agree = sum(1 for x in range(16)
                    if echo_value(x, Q, B, 4, maxfirst=mf)[0] == Q[x])
        res[mf] = agree
    tag = "LINEAR(B=0)" if Bzero else ""
    exact = [('P1max' if mf else 'P1min') for mf in (True, False) if res[mf] == 16]
    print(f"lam={lam:2d} |Z|={rank_info:2d} agree(P1max)={res[True]:2d}/16 "
          f"agree(P1min)={res[False]:2d}/16 exact={exact} {tag}")
