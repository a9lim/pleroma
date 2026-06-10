"""Screen alternation-discipline variants on the abstract linking game.

Exactness of ANY echo-style rule on big fields REQUIRES: for every k-coin
B-pattern b (q=0), the forced value of sum b_ij*linked(i,j) equals the
all-linked parity sum b_ij, under ONE fixed orientation (since all patterns
occur as B-restrictions of Gold forms at m>=8).

k=3 has 8 patterns, k=4 has 64, k=5 has 1024.
"""
import sys
from itertools import product
sys.setrecursionlimit(100000)


def solve_banwin(k, B, w, p1_target, pass_parity=False):
    """ban-window w: may not touch any of the last w touched coins; pass clears.
    pass_parity: readout sigma ^= (number of passes mod 2)."""
    def charge(o, i):
        acc = 0
        for kk in range(k):
            if kk > i and (o >> kk) & 1 and B[kk][i]:
                acc ^= 1
        return acc
    memo = {}
    def rec(u, o, recent, mover, sigma, pp):
        if u == 0 and o == 0:
            return sigma ^ (pp if pass_parity else 0)
        key = (u, o, recent, mover, sigma, pp)
        r = memo.get(key)
        if r is not None:
            return r
        legal = []
        for i in range(k):
            if i in recent:
                continue
            bit = 1 << i
            if u & bit:
                legal.append((i, u ^ bit, o ^ bit))
            elif o & bit:
                legal.append((i, u, o ^ bit))
        if not legal:
            res = rec(u, o, (), 1 - mover, sigma, pp ^ 1)
        else:
            want = p1_target if mover == 0 else 1 - p1_target
            res = 1 - want
            for (i, u2, o2) in legal:
                r2 = rec(u2, o2, ((i,) + recent)[:w], 1 - mover,
                         sigma ^ charge(o, i), pp)
                if r2 == want:
                    res = want
                    break
        memo[key] = res
        return res
    full = (1 << k) - 1
    return rec(full, 0, (), 0, 0, 0)


def solve_closedelay(k, B, d, p1_target):
    """close-delay d: a coin may be closed only if >= d touches/passes since
    it was opened. Per-coin code: 0 untouched; 1+a open with age a (cap d); -1 done."""
    def charge(codes, i):
        acc = 0
        for kk in range(k):
            if kk > i and codes[kk] >= 1 and B[kk][i]:
                acc ^= 1
        return acc
    memo = {}
    def age(codes):
        return tuple(min(c + 1, 1 + d) if c >= 1 else c for c in codes)
    def rec(codes, mover, sigma):
        if all(c == -1 for c in codes):
            return sigma
        key = (codes, mover, sigma)
        r = memo.get(key)
        if r is not None:
            return r
        legal = []
        for i in range(k):
            c = codes[i]
            if c == 0:
                nc = list(codes)
                nc[i] = 1  # open, age 0
                legal.append((i, tuple(nc)))
            elif c >= 1 and (c - 1) >= d:
                nc = list(codes)
                nc[i] = -1
                legal.append((i, tuple(nc)))
        if not legal:
            res = rec(age(codes), 1 - mover, sigma)
        else:
            want = p1_target if mover == 0 else 1 - p1_target
            res = 1 - want
            for (i, nc) in legal:
                r2 = rec(age(nc), 1 - mover, sigma ^ charge(codes, i))
                # note: aging applies after each touch too (uniform tempo)
                if r2 == want:
                    res = want
                    break
        memo[key] = res
        return res
    return rec(tuple(0 for _ in range(k)), 0, 0)


def solve_fifo(k, B, discipline, ko1, p1_target):
    """discipline 'fifo': may close only the longest-open coin;
    'lifo': only the most recently opened. ko1: w=1 touch ban on top."""
    def charge(openseq, i):
        acc = 0
        for kk in openseq:
            if kk > i and B[kk][i]:
                acc ^= 1
        return acc
    memo = {}
    def rec(u, openseq, last, mover, sigma):
        if u == 0 and not openseq:
            return sigma
        key = (u, openseq, last, mover, sigma)
        r = memo.get(key)
        if r is not None:
            return r
        legal = []
        for i in range(k):
            if ko1 and i == last:
                continue
            if (u >> i) & 1:
                legal.append((i, u ^ (1 << i), openseq + (i,)))
        if openseq:
            c = openseq[0] if discipline == "fifo" else openseq[-1]
            if not (ko1 and c == last):
                rest = openseq[1:] if discipline == "fifo" else openseq[:-1]
                legal.append((c, u, rest))
        if not legal:
            res = rec(u, openseq, -1, 1 - mover, sigma)
        else:
            want = p1_target if mover == 0 else 1 - p1_target
            res = 1 - want
            for (i, u2, seq2) in legal:
                r2 = rec(u2, seq2, i, 1 - mover, sigma ^ charge(openseq, i))
                if r2 == want:
                    res = want
                    break
        memo[key] = res
        return res
    full = (1 << k) - 1
    return rec(full, (), -1, 0, 0)


def patterns(k):
    pairs = [(i, j) for i in range(k) for j in range(i + 1, k)]
    for bits in product((0, 1), repeat=len(pairs)):
        B = [[0] * k for _ in range(k)]
        for (b, (i, j)) in zip(bits, pairs):
            B[i][j] = B[j][i] = b
        yield bits, B, sum(bits) & 1


def screen(name, solver_for, ks=(3, 4)):
    for tgt in (1, 0):
        ok_all = True
        firstfail = None
        for k in ks:
            for bits, B, par in patterns(k):
                v = solver_for(k, B, tgt)
                if v != par:
                    ok_all = False
                    firstfail = (k, bits, par, v)
                    break
            if not ok_all:
                break
        status = "PASS" if ok_all else f"fail @ k={firstfail[0]} b={firstfail[1]} parity={firstfail[2]} forced={firstfail[3]}"
        print(f"{name:34s} tgt={'A' if tgt else 'B'}: {status}")


print("--- abstract pattern screen: forced(b) ?= parity(b), all patterns, k=3,4 ---")
screen("ban-window w=1 (current)", lambda k, B, t: solve_banwin(k, B, 1, t))
screen("ban-window w=2", lambda k, B, t: solve_banwin(k, B, 2, t))
screen("ban-window w=3", lambda k, B, t: solve_banwin(k, B, 3, t))
screen("ban-window w=1 + pass-parity", lambda k, B, t: solve_banwin(k, B, 1, t, True))
screen("ban-window w=2 + pass-parity", lambda k, B, t: solve_banwin(k, B, 2, t, True))
screen("close-delay d=1", lambda k, B, t: solve_closedelay(k, B, 1, t))
screen("close-delay d=2", lambda k, B, t: solve_closedelay(k, B, 2, t))
screen("close-delay d=3", lambda k, B, t: solve_closedelay(k, B, 3, t))
screen("FIFO-close", lambda k, B, t: solve_fifo(k, B, "fifo", False, t))
screen("FIFO-close + ko1", lambda k, B, t: solve_fifo(k, B, "fifo", True, t))
screen("LIFO-close (sanity: nested only)", lambda k, B, t: solve_fifo(k, B, "lifo", False, t))
screen("LIFO-close + ko1", lambda k, B, t: solve_fifo(k, B, "lifo", True, t))
