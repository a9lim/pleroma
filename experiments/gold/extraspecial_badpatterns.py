"""Inspect bad k=3 patterns for ko=self; test ko=self==ko=opp; try B-adaptive kos."""
from extraspecial_core import Echo
import itertools

def mk_Bh(k, Bdict):
    Bh = []
    for i in range(k):
        mask = 0
        for j in range(k):
            if j > i and Bdict.get((min(i, j), max(i, j)), 0):
                mask |= 1 << j
        Bh.append(mask)
    return Bh

def target(qbits, Bdict):
    t = 0
    for b in qbits: t ^= b
    for _, b in Bdict.items(): t ^= b
    return t

pairs3 = list(itertools.combinations(range(3), 2))
print("=== k=3 bad patterns for ko=self (P1max) ===")
bad3 = []
for qm in range(8):
    q = tuple((qm >> i) & 1 for i in range(3))
    for bm in range(8):
        B = {pairs3[p]: (bm >> p) & 1 for p in range(3)}
        v, ch = Echo(list(q), mk_Bh(3, B), ko='self').solve(True)
        t = target(q, B)
        if v != t:
            bad3.append((q, B))
            bpairs = tuple(p for p in pairs3 if B[p])
            print(f"q={q} B1={bpairs}  target={t} val={v}")
# classify: degree of B-graph, q-weight
print(f"\ntotal bad: {len(bad3)}")
from collections import Counter
cnt = Counter((sum(q), sum(B.values())) for q, B in bad3)
print("histogram (q-weight, B-edges):", dict(cnt))

# check ko=self vs ko=opp equality on all k=3,4 patterns (P1max)
same = True
for k in (3, 4):
    prs = list(itertools.combinations(range(k), 2))
    for qm in range(1 << k):
        q = tuple((qm >> i) & 1 for i in range(k))
        for bm in range(1 << len(prs)):
            B = {prs[p]: (bm >> p) & 1 for p in range(len(prs))}
            v1, _ = Echo(list(q), mk_Bh(k, B), ko='self').solve(True)
            v2, _ = Echo(list(q), mk_Bh(k, B), ko='opp').solve(True)
            if v1 != v2: same = False; print("DIFF", k, q, B, v1, v2); break
print("ko=self == ko=opp on all k=3,4 patterns:", same)
