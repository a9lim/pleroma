"""Pattern-level characterization: the ECHO game value depends only on
(q|S, B|S, ko, orientation).  Enumerate ALL k=3 and k=4 patterns and find
which ko variant (if any) solves every pattern, i.e. value == Q(pattern)
where Q = sum q_i + sum_{i<j} B_ij  (the all-linked target).
"""
from extraspecial_core import Echo
import itertools

def solve_pattern(k, qbits, Bdict, ko, maxfirst):
    # local coins 0..k-1 in increasing global order; Bhigh[i] = {j>i: B_ij=1}
    Bh = []
    for i in range(k):
        mask = 0
        for j in range(k):
            if j > i and Bdict.get((i, j), 0):
                mask |= 1 << j
        Bh.append(mask)
    v, ch = Echo(list(qbits), Bh, ko=ko).solve(maxfirst=maxfirst)
    return v, ch

def target(k, qbits, Bdict):
    t = 0
    for b in qbits: t ^= b
    for (i, j), b in Bdict.items(): t ^= b
    return t

def sweep(k, kos=('self', 'opp', 'w2', 'none')):
    pairs = list(itertools.combinations(range(k), 2))
    npairs = len(pairs)
    results = {}
    for ko in kos:
        for mf in (True, False):
            bad = []
            ndg = 0   # solved patterns that are decision-nondegenerate
            total = 0
            for qm in range(1 << k):
                qbits = tuple((qm >> i) & 1 for i in range(k))
                for bm in range(1 << npairs):
                    Bdict = {pairs[p]: (bm >> p) & 1 for p in range(npairs)}
                    total += 1
                    v, ch = solve_pattern(k, qbits, Bdict, ko, mf)
                    if v != target(k, qbits, Bdict):
                        bad.append((qbits, tuple(sorted((p, b) for p, b in Bdict.items() if b))))
                    elif ch > 0:
                        ndg += 1
            results[(ko, mf)] = (total - len(bad), total, bad, ndg)
            ori = 'P1max' if mf else 'P1min'
            print(f"k={k} ko={ko:4s} {ori}: solved {total-len(bad)}/{total} "
                  f"(nondeg among solved: {ndg})")
            if 0 < len(bad) <= 8:
                for b in bad: print(f"    BAD q={b[0]} B1pairs={b[1]}")
    return results

print("=== k=3 pattern sweep ===")
r3 = sweep(3)
print("\n=== k=4 pattern sweep ===")
r4 = sweep(4)
