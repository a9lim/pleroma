"""ECHO-FIFO real Gold-form bench via decomposition, validated directly on m=4."""
import sys, time
from itertools import permutations
sys.path.insert(0, "/tmp")
from asym2_probe import make_form
from asym2_fifo import abstract_value, canon

sys.setrecursionlimit(100000)

_cache = {}


def support_value(S, qd, B, t):
    """Game value of position with support S under ECHO-FIFO, P1 wants sigma=t.
    value = l XOR V*(b_S, want = t XOR l)."""
    k = len(S)
    if k == 0:
        return 0
    l = 0
    for i in S:
        l ^= qd[i]
    edges = frozenset((a, b) for ai, a_ in enumerate(S) for bi, b_ in enumerate(S)
                      if ai < bi and B[a_][b_]
                      for (a, b) in [(ai, bi)])
    want = t ^ l
    if k <= 6:
        key = (k, canon(k, edges), want)
    else:
        key = (k, tuple(sorted(edges)), want)
    v = _cache.get(key)
    if v is None:
        v = abstract_value(k, edges, want)
        _cache[key] = v
    return l ^ v


# ------------------------------------------------ direct real-state solver (m<=4)
def direct_fifo_value(x, m, qd, B, t):
    bits = [i for i in range(m) if (x >> i) & 1]
    if not bits:
        return 0
    memo = {}

    def charge(omask, openset_has_i, i):
        acc = qd[i] if openset_has_i else 0
        for kk in bits:
            if kk > i and (omask >> kk) & 1 and B[kk][i]:
                acc ^= 1
        return acc

    def rec(u, openseq, last, mover, sigma):
        if u == 0 and not openseq:
            return sigma
        key = (u, openseq, last, mover, sigma)
        r = memo.get(key)
        if r is not None:
            return r
        omask = 0
        for c in openseq:
            omask |= 1 << c
        legal = []
        for i in bits:
            if i == last:
                continue
            if (u >> i) & 1:
                legal.append((i, charge(omask, False, i), u ^ (1 << i), openseq + (i,)))
        if openseq:
            c = openseq[0]
            if c != last:
                legal.append((c, charge(omask, True, c), u, openseq[1:]))
        if not legal:
            res = rec(u, openseq, -1, 1 - mover, sigma)
        else:
            want = t if mover == 0 else 1 - t
            res = 1 - want
            for (i, ch, u2, seq2) in legal:
                if rec(u2, seq2, i, 1 - mover, sigma ^ ch) == want:
                    res = want
                    break
        memo[key] = res
        return res

    return rec(x, 0, (), -1, 0) if False else rec(x, (), -1, 0, 0)


# fix call shape: rec(u, openseq, last, mover, sigma)
def direct_fifo_value2(x, m, qd, B, t):
    return direct_fifo_value(x, m, qd, B, t)


# validate decomposition against direct solver on all m=4 forms
print("validating decomposition vs direct real-state solver on m=4 ...")
for lam in range(1, 16):
    Q, qd, B = make_form(4, 1, lam)
    for t in (0, 1):
        for x in range(16):
            S = tuple(i for i in range(4) if (x >> i) & 1)
            assert support_value(S, qd, B, t) == direct_fifo_value(x, 4, qd, B, t), \
                (lam, t, x)
print("decomposition == direct solver (all 15 forms x 2 t x 16 x): OK")

print()
print("=" * 70)
print("ECHO-FIFO m=4, a=1: exactness sweep (t=0: P1 wants 0; t=1: P1 wants 1)")
print("=" * 70)
hits = {}
for lam in range(1, 16):
    Q, qd, B = make_form(4, 1, lam)
    row = []
    for t in (1, 0):
        agree = sum(1 for x in range(16)
                    if support_value(tuple(i for i in range(4) if (x >> i) & 1),
                                     qd, B, t) == Q[x])
        row.append(f"t={t}:{agree}{'*' if agree == 16 else ''}")
        if agree == 16:
            hits.setdefault(lam, []).append(t)
    print(f"lam={lam:2d} |Q=0|={sum(1 for v in Q if v==0):2d}  " + "  ".join(row))
print("EXACT hits m=4:", sorted(hits.items()))

print()
print("=" * 70)
print("ECHO-FIFO m=8 forms")
print("=" * 70)
for (m, a, lam) in [(8, 1, 1), (8, 2, 1), (8, 1, 2), (8, 1, 3)]:
    Q, qd, B = make_form(m, a, lam)
    for t in (1, 0):
        t0 = time.time()
        misses = []
        for x in range(256):
            S = tuple(i for i in range(m) if (x >> i) & 1)
            v = support_value(S, qd, B, t)
            if v != Q[x]:
                misses.append(x)
        agree = 256 - len(misses)
        mtxt = ""
        if 0 < len(misses) <= 10:
            mtxt = " misses=" + ",".join(
                f"{x}(pc{bin(x).count('1')})" for x in misses)
        print(f"(m={m},a={a},lam={lam}) t={t}: {agree}/256"
              f"{' EXACT' if agree == 256 else ''}{mtxt} [{time.time()-t0:.0f}s]")
