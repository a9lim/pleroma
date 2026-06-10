"""k=5: verify value depends only on (B-graph, target); census bad graphs.
Also extend reduced table v(p,s) to larger s."""
import itertools
from collections import defaultdict
from extraspecial_core import Echo

def mk_Bh(k, B):
    Bh = []
    for i in range(k):
        m = 0
        for j in range(k):
            if j > i and B.get((i, j), 0): m |= 1 << j
        Bh.append(m)
    return Bh

print("=== extended reduced table v(p,s) (q=0, P1max/P1min) ===")
for p in range(1, 4):
    row = []
    for s in range(0, 7):
        k = 2 * p + s
        if k > 10: row.append("-"); continue
        edges = [(2 * i, 2 * i + 1) for i in range(p)]
        Bh = mk_Bh(k, {e: 1 for e in edges})
        vmax, _ = Echo([0] * k, Bh, ko='self').solve(True)
        vmin, _ = Echo([0] * k, Bh, ko='self').solve(False)
        row.append(f"{vmax}{vmin}")
    print(f"p={p}: " + " ".join(row))

print("\n=== k=5: (graph,target)-dependence check ===")
k = 5
prs = list(itertools.combinations(range(k), 2))
fn = defaultdict(set)
badgraphs = defaultdict(set)   # bm -> set of failing targets
for qm in range(1 << k):
    q = [(qm >> i) & 1 for i in range(k)]
    sq = sum(q) & 1
    for bm in range(1 << len(prs)):
        B = {prs[x]: (bm >> x) & 1 for x in range(len(prs))}
        t = (sq + sum(B.values())) & 1
        v, _ = Echo(q, mk_Bh(k, B), ko='self').solve(True)
        fn[(bm, t)].add(v)
        if v != t: badgraphs[bm].add(t)
dep = all(len(s) == 1 for s in fn.values())
print("value determined by (B-graph, target):", dep)
nbad = len(badgraphs)
print(f"graphs with some failing target: {nbad}/{1 << len(prs)}")
# census by (edges, maxdeg, iso-vertices)
from collections import Counter
cen = Counter()
for bm, ts in badgraphs.items():
    edges = [prs[x] for x in range(len(prs)) if (bm >> x) & 1]
    deg = [0] * k
    for a, b in edges: deg[a] += 1; deg[b] += 1
    iso = sum(1 for d in deg if d == 0)
    cen[(len(edges), max(deg) if edges else 0, iso, tuple(sorted(ts)))] += 1
print("census (edges, maxdeg, iso, failing-targets): count")
for key in sorted(cen): print(f"  {key}: {cen[key]}")
