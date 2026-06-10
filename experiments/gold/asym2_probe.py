"""Round-2 [asymmetry] probe: ECHO-ko on the extraspecial group E, CORRECT solver.

Model: V = F_2^m, Q = scaled Gold form Tr(lam x^(1+2^a)), B polar.
Triangular cocycle c(u,v) = sum_i q_i u_i v_i XOR sum_{k>j} B_kj u_k v_j
(lower side; upper side flips the strict inequality). E = V x F_2 with
(s,u)(t,v) = (s+t+c(u,v), u+v): squaring map Q, commutator B.

Game ECHO(x): coins = bits(x), each touched exactly twice (open, then close).
Players alternate. Touching coin i with open set o charges sigma ^= c(o, e_i)
= right-multiplication of the running word by the lift (0, e_i). Ko: may not
touch the coin touched in the immediately preceding touch; stuck => pass
(clears ko). Complete word lies over 0 in V, equals 1 or z in E; readout =
central character. Orientation A: P1 wants sigma=1. Orientation B: P1 wants 0.

CORRECT solver: memo key includes accumulated sigma (the round-1 probe's bug
was omitting it: XOR payoff => odd prefix flips downstream objective).
Validated below against explicit no-memo tree enumeration.
"""
import sys
from functools import lru_cache

sys.setrecursionlimit(100000)


# ---------------------------------------------------------------- nim arithmetic
@lru_cache(maxsize=None)
def nim_mul(a, b):
    if a < b:
        a, b = b, a
    if b == 0:
        return 0
    if b == 1:
        return a
    k = 0
    while a >= (1 << (1 << (k + 1))):
        k += 1
    sh = 1 << k
    F = 1 << sh
    ah, al = a >> sh, a & (F - 1)
    bh, bl = b >> sh, b & (F - 1)
    t1 = nim_mul(ah, bh)
    t2 = nim_mul(ah, bl) ^ nim_mul(al, bh)
    t3 = nim_mul(al, bl)
    return ((t1 ^ t2) << sh) ^ nim_mul(t1, F >> 1) ^ t3


assert nim_mul(2, 2) == 3 and nim_mul(2, 4) == 8 and nim_mul(16, 16) == 24


def frob(x, a):
    for _ in range(a):
        x = nim_mul(x, x)
    return x


def trace(x, m):
    acc, t = x, x
    for _ in range(m - 1):
        t = nim_mul(t, t)
        acc ^= t
    assert acc in (0, 1)
    return acc


def make_form(m, a, lam):
    n = 1 << m
    Q = [trace(nim_mul(lam, nim_mul(v, frob(v, a))), m) for v in range(n)]
    qd = [Q[1 << i] for i in range(m)]
    B = [[Q[(1 << i) ^ (1 << j)] ^ Q[1 << i] ^ Q[1 << j] if i != j else 0
          for j in range(m)] for i in range(m)]
    return Q, qd, B


# pinned zero counts (goldarf.tex tables / round-1 cross-checks)
for (m, a, lam, zc) in [(4, 1, 1, 4), (8, 1, 1, 112), (8, 2, 1, 96), (8, 1, 2, 136)]:
    Q, _, _ = make_form(m, a, lam)
    assert sum(1 for v in Q if v == 0) == zc, (m, a, lam, sum(1 for v in Q if v == 0))
print("nim arithmetic + Gold zero counts: OK")


# ---------------------------------------------------------------- cocycle/charge
def cocycle(u, v, qd, B, m, lower=True):
    acc = 0
    for i in range(m):
        if (u >> i) & 1 and (v >> i) & 1:
            acc ^= qd[i]
    for k in range(m):
        if not (u >> k) & 1:
            continue
        for j in range(m):
            if (v >> j) & 1 and ((k > j) if lower else (k < j)):
                acc ^= B[k][j]
    return acc


def make_charge(qd, B, m, lower=True):
    # charge(o, i) = c(o, e_i): q_i if i open, plus B[k][i] over open k on the side
    side_mask = [0] * m
    for i in range(m):
        msk = 0
        for k in range(m):
            if ((k > i) if lower else (k < i)) and B[k][i]:
                msk |= 1 << k
        side_mask[i] = msk

    def charge(o, i):
        acc = qd[i] if (o >> i) & 1 else 0
        acc ^= bin(o & side_mask[i]).count("1") & 1
        return acc

    return charge


# identities
import random
rng = random.Random(2026)
for (m, a, lam) in [(4, 1, 1), (4, 1, 2), (8, 1, 1)]:
    Q, qd, B = make_form(m, a, lam)
    n = 1 << m
    pairs = [(u, v) for u in range(n) for v in range(n)] if m <= 4 else \
        [(rng.randrange(n), rng.randrange(n)) for _ in range(500)]
    for side in (True, False):
        for (u, v) in pairs:
            assert cocycle(v, v, qd, B, m, side) == Q[v]
            assert cocycle(u, v, qd, B, m, side) ^ cocycle(v, u, qd, B, m, side) == \
                (Q[u ^ v] ^ Q[u] ^ Q[v])
        ch = make_charge(qd, B, m, side)
        for _ in range(300):
            o, i = rng.randrange(n), rng.randrange(m)
            assert ch(o, i) == cocycle(o, 1 << i, qd, B, m, side)
print("cocycle identities (both sides): OK")

# chord-linking formula on random complete plays (no ko): sigma = l_diag + sum B_ij*linked
for (m, a, lam) in [(4, 1, 2), (8, 1, 1)]:
    Q, qd, B = make_form(m, a, lam)
    ch = make_charge(qd, B, m, True)
    for _ in range(300):
        x = rng.randrange(1, 1 << m)
        bits = [i for i in range(m) if (x >> i) & 1]
        seq = bits * 2
        rng.shuffle(seq)
        o, sig = 0, 0
        for i in seq:
            sig ^= ch(o, i)
            o ^= 1 << i
        pos = {}
        for t, i in enumerate(seq):
            pos.setdefault(i, []).append(t)
        pred = 0
        for i in bits:
            pred ^= qd[i]
        for ii in range(len(bits)):
            for jj in range(ii + 1, len(bits)):
                i, j = bits[ii], bits[jj]
                a1, a2 = pos[i]
                b1, b2 = pos[j]
                linked = (a1 < b1 < a2 < b2) or (b1 < a1 < b2 < a2)
                if linked:
                    pred ^= B[i][j]
        assert sig == pred, (x, seq)
print("chord-linking formula: OK")


# ---------------------------------------------------------------- CORRECT solver
def solve_position(x, m, charge, p1_target, prune=True):
    """Final sigma under optimal play; P1 (first mover) wants p1_target."""
    bits = [i for i in range(m) if (x >> i) & 1]
    if not bits:
        return 0
    memo = {}

    def rec(u, o, last, mover, sigma):
        if u == 0 and o == 0:
            return sigma
        key = (u, o, last, mover, sigma)
        r = memo.get(key)
        if r is not None:
            return r
        legal = []
        for i in bits:
            if i == last:
                continue
            bit = 1 << i
            if u & bit:
                legal.append((i, u ^ bit, o ^ bit))
            elif o & bit:
                legal.append((i, u, o ^ bit))
        if not legal:
            res = rec(u, o, -1, 1 - mover, sigma)
            memo[key] = res
            return res
        want = p1_target if mover == 0 else 1 - p1_target
        res = 1 - want
        for (i, u2, o2) in legal:
            r2 = rec(u2, o2, i, 1 - mover, sigma ^ charge(o, i))
            if r2 == want:
                res = want
                if prune:
                    break
        memo[key] = res
        return res

    return rec(x, 0, -1, 0, 0)


# explicit no-memo tree enumeration (validator)
def solve_explicit(x, m, charge, p1_target):
    bits = [i for i in range(m) if (x >> i) & 1]
    if not bits:
        return 0

    def rec(u, o, last, mover, sigma):
        if u == 0 and o == 0:
            return sigma
        legal = []
        for i in bits:
            if i == last:
                continue
            bit = 1 << i
            if u & bit:
                legal.append((i, u ^ bit, o ^ bit))
            elif o & bit:
                legal.append((i, u, o ^ bit))
        if not legal:
            return rec(u, o, -1, 1 - mover, sigma)
        want = p1_target if mover == 0 else 1 - p1_target
        outs = [rec(u2, o2, i, 1 - mover, sigma ^ charge(o, i)) for (i, u2, o2) in legal]
        return want if want in outs else 1 - want

    return rec(x, 0, -1, 0, 0)


# validate solver: all 16 positions of two m=4 forms, both orientations, both sides
for (m, a, lam) in [(4, 1, 1), (4, 1, 2), (4, 1, 7)]:
    Q, qd, B = make_form(m, a, lam)
    for side in (True, False):
        ch = make_charge(qd, B, m, side)
        for tgt in (1, 0):
            for x in range(1 << m):
                assert solve_position(x, m, ch, tgt) == solve_explicit(x, m, ch, tgt)
# random m=8 popcount<=4 spot checks
Q, qd, B = make_form(8, 2, 1)
ch = make_charge(qd, B, 8, True)
cnt = 0
for x in range(256):
    if bin(x).count("1") <= 3:
        assert solve_position(x, 8, ch, 1) == solve_explicit(x, 8, ch, 1)
        cnt += 1
print(f"solver == explicit tree enumeration: OK ({cnt} m=8 spot checks + full m=4)")
