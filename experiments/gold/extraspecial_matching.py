"""Matching-pattern exactness + the reduced (p,s) scheduling game.

If B|S is a partial matching (max degree 1), sigma = sum q_i + sum_matched linked.
Exactness for all matching patterns <=> the 'even-unlinked' player always forces
unlinked-parity 0 in the reduced game v(p, s).
"""
from extraspecial_core import Echo
import itertools

def mk_Bh(k, Bedges):
    Bh = []
    for i in range(k):
        mask = 0
        for (a, b) in Bedges:
            if a == i and b > i: mask |= 1 << b
            elif b == i and a > i: mask |= 1 << a
        Bh.append(mask)
    return Bh

def matchings(k):
    """all partial matchings on k vertices as edge lists"""
    verts = list(range(k))
    out = [[]]
    def rec(avail, cur):
        if len(avail) < 2: return
        a = avail[0]
        rest = avail[1:]
        # a unmatched
        rec(rest, cur)
        for i, b in enumerate(rest):
            e = cur + [(a, b)]
            out.append(e)
            rec(rest[:i] + rest[i+1:], e)
    rec(verts, [])
    # dedupe
    seen = set()
    uniq = []
    for e in out:
        key = tuple(sorted(e))
        if key not in seen:
            seen.add(key); uniq.append(e)
    return uniq

print("=== matching patterns, ko=self ===")
for k in range(3, 7):
    Ms = matchings(k)
    bad = 0; tot = 0; ndg = 0
    for edges in Ms:
        for qm in range(1 << k):
            q = [(qm >> i) & 1 for i in range(k)]
            t = (sum(q) + len(edges)) & 1
            for mf in (True, False):
                tot += 1
                v, ch = Echo(q, mk_Bh(k, edges), ko='self').solve(mf)
                if v != t: bad += 1
                elif ch: ndg += 1
    print(f"k={k}: {len(Ms)} matchings, {tot} (pattern,orient) cases, "
          f"bad={bad}, nondeg-solved={ndg}")

print("\n=== reduced game v(p, s): forced unlinked-parity "
      "(q=0, B=p disjoint edges, s isolated coins) ===")
print("rows p=1..4, cols s=0..4; entry (P1max-val, P1min-val); target = p&1")
for p in range(1, 5):
    row = []
    for s in range(0, 5):
        k = 2 * p + s
        if k > 9:
            row.append("  -  ")
            continue
        edges = [(2 * i, 2 * i + 1) for i in range(p)]
        q = [0] * k
        vmax, _ = Echo(q, mk_Bh(k, edges), ko='self').solve(True)
        vmin, _ = Echo(q, mk_Bh(k, edges), ko='self').solve(False)
        t = p & 1
        ok = "OK " if (vmax == t and vmin == t) else "BAD"
        row.append(f"{vmax}{vmin}{ok}")
    print(f"p={p}: " + " | ".join(str(c) for c in row))
