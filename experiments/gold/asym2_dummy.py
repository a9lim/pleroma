"""ECHO-FIFO + one dummy coin (q=0, B-isolated): the parity fix.

Even-size boards are perfect on both orientations (k=4, k=6 exhaustive up to
iso). Odd-k failures all occur on classes with no isolated vertex. Adjoining
one B-isolated, q=0 dummy coin to EVERY board (a uniform rule ingredient)
maps a k-support to a (k+1)-board containing an isolated vertex.
Test: do the boards-with-isolated-vertex tables stay perfect at 5,7,9?
Then: full m=4 and m=8 Gold benches.
"""
import sys, time
from itertools import permutations
sys.path.insert(0, "/tmp")
from asym2_probe import make_form
from asym2_fifo import abstract_value, canon
from asym2_fifo_bench import direct_fifo_value

sys.setrecursionlimit(1000000)

# 1) which k=5 classes fail want1? do they have isolated vertices?
print("k=5 want1 failing classes and isolated-vertex check:")
from itertools import product
pairs5 = [(i, j) for i in range(5) for j in range(i + 1, 5)]
seen = set()
for bits in product((0, 1), repeat=10):
    edges = frozenset(p for (b, p) in zip(bits, pairs5) if b)
    c = canon(5, edges)
    if c in seen:
        continue
    seen.add(c)
    par = len(edges) & 1
    for want, label in [(1, "want1"), (0, "want0")]:
        if abstract_value(5, edges, want) != par:
            deg = [0] * 5
            for (i, j) in edges:
                deg[i] += 1
                deg[j] += 1
            print(f"  {label} fail: {c}  degrees={deg}  isolated={0 in deg}")

# 2) boards-with-isolated-vertex screens at sizes 5 and 7:
#    size-5 boards = all k=4 classes + isolated vertex
#    size-7 boards = all k=6 classes + isolated vertex
def screen_iso_boards(kbase):
    pairs = [(i, j) for i in range(kbase) for j in range(i + 1, kbase)]
    classes = {}
    for bits in product((0, 1), repeat=len(pairs)):
        edges = frozenset(p for (b, p) in zip(bits, pairs) if b)
        c = canon(kbase, edges)
        if c not in classes:
            classes[c] = edges
    bad = []
    for c, edges in classes.items():
        par = len(edges) & 1
        for want in (0, 1):
            v = abstract_value(kbase + 1, edges, want)  # vertex kbase isolated
            if v != par:
                bad.append((c, want, v))
    return len(classes), bad

for kbase in (2, 3, 4, 5, 6):
    t0 = time.time()
    nc, bad = screen_iso_boards(kbase)
    print(f"boards size {kbase+1} = k={kbase} classes + dummy: {nc} classes, "
          f"failures={len(bad)} [{time.time()-t0:.0f}s]")
    for b in bad[:4]:
        print("   FAIL:", b)

# 3) real-form benches with the dummy
_cache = {}

def support_value_dummy(S, qd, B, t, m):
    k = len(S)
    l = 0
    for i in S:
        l ^= qd[i]
    edges = frozenset((ai, bi) for ai in range(k) for bi in range(ai + 1, k)
                      if B[S[ai]][S[bi]])
    want = t ^ l
    kk = k + 1  # + dummy as highest index (isolated)
    if kk <= 6:
        key = (kk, canon(kk, edges), want)
    else:
        key = (kk, tuple(sorted(edges)), want)
    v = _cache.get(key)
    if v is None:
        v = abstract_value(kk, edges, want)
        _cache[key] = v
    return l ^ v


# validate against direct solver with explicitly extended tables, all m=4 forms
print("\nvalidating dummy decomposition vs direct solver on m=4 ...")
for lam in range(1, 16):
    Q, qd, B = make_form(4, 1, lam)
    qd5 = qd + [0]
    B5 = [row + [0] for row in B] + [[0] * 5]
    for t in (0, 1):
        for x in range(16):
            S = tuple(i for i in range(4) if (x >> i) & 1)
            direct = direct_fifo_value(x | 16, 5, qd5, B5, t)
            assert support_value_dummy(S, qd, B, t, 4) == direct, (lam, t, x)
print("dummy decomposition == direct solver (15 forms x 2 t x 16 x): OK")

print()
print("ECHO-FIFO+dummy m=4, a=1 sweep:")
hits = {}
for lam in range(1, 16):
    Q, qd, B = make_form(4, 1, lam)
    for t in (1, 0):
        agree = sum(1 for x in range(16)
                    if support_value_dummy(tuple(i for i in range(4) if (x >> i) & 1),
                                           qd, B, t, 4) == Q[x])
        if agree == 16:
            hits.setdefault(lam, []).append(t)
        else:
            print(f"  lam={lam} t={t}: {agree}/16")
print("EXACT hits m=4:", sorted(hits.items()))

print()
print("ECHO-FIFO+dummy m=8 (popcount <= 7 first, then the full board):")
for (m, a, lam) in [(8, 1, 1), (8, 2, 1), (8, 1, 2), (8, 1, 3)]:
    Q, qd, B = make_form(m, a, lam)
    for t in (1, 0):
        t0 = time.time()
        misses = []
        for x in range(256):
            if bin(x).count("1") > 7:
                continue
            S = tuple(i for i in range(m) if (x >> i) & 1)
            if support_value_dummy(S, qd, B, t, m) != Q[x]:
                misses.append(x)
        print(f"(m={m},a={a},lam={lam}) t={t}: pc<=7: "
              f"{255-len(misses)}/255{' ALL' if not misses else ''} "
              f"misses={misses[:8]} [{time.time()-t0:.0f}s]")
