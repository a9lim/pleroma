"""Characterize the k=4 value function f(q, B|S) for ko=self.
Conjecture from k=3: failure iff target=0 and B-graph 'dense enough'.
Find the exact predicate."""
import itertools
from extraspecial_core import Echo
from collections import Counter

def mk_Bh(k, B):
    Bh = []
    for i in range(k):
        m = 0
        for j in range(k):
            if j > i and B.get((i, j), 0): m |= 1 << j
        Bh.append(m)
    return Bh

k = 4
prs = list(itertools.combinations(range(k), 2))
records = []
for qm in range(1 << k):
    q = [(qm >> i) & 1 for i in range(k)]
    for bm in range(1 << len(prs)):
        B = {prs[x]: (bm >> x) & 1 for x in range(len(prs))}
        t = (sum(q) + sum(B.values())) & 1
        v, ch = Echo(q, mk_Bh(k, B), ko='self').solve(True)
        records.append((tuple(q), bm, t, v))

bad = [(q, bm, t, v) for q, bm, t, v in records if t != v]
print(f"k=4 bad: {len(bad)}/1024")
# one-sidedness?
print("direction histogram (target, val):", Counter((t, v) for _, _, t, v in bad))

# graph-theoretic stats of bad B-graphs
def graph_stats(bm):
    edges = [prs[x] for x in range(len(prs)) if (bm >> x) & 1]
    deg = [0] * k
    for a, b in edges:
        deg[a] += 1; deg[b] += 1
    return len(edges), max(deg) if edges else 0, tuple(sorted(deg))

cnt_bad = Counter()
cnt_all = Counter()
for q, bm, t, v in records:
    e, mx, degs = graph_stats(bm)
    cnt_all[(e, mx, t)] += 1
for q, bm, t, v in bad:
    e, mx, degs = graph_stats(bm)
    cnt_bad[(e, mx, t)] += 1
print("\n(edges, maxdeg, target): bad/all")
for key in sorted(cnt_all):
    b = cnt_bad.get(key, 0)
    if b or key[0] >= 2:
        print(f"  {key}: {b}/{cnt_all[key]}")

# is value a function of (target, B-graph) only (q enters only via target)?
from collections import defaultdict
fn = defaultdict(set)
for q, bm, t, v in records:
    fn[(bm, t)].add(v)
print("\nvalue determined by (B-graph, target)?",
      all(len(s) == 1 for s in fn.values()))
# if yes, print the exceptional (B-graph -> which targets fail)
ex = {}
for (bm, t), s in fn.items():
    v = next(iter(s)) if len(s) == 1 else None
    if v is not None and v != t:
        ex.setdefault(bm, []).append(t)
if all(len(s) == 1 for s in fn.values()):
    print("B-graphs with failures (edge lists), failing targets:")
    for bm, ts in sorted(ex.items()):
        edges = [prs[x] for x in range(len(prs)) if (bm >> x) & 1]
        e, mx, degs = graph_stats(bm)
        print(f"  edges={edges} maxdeg={mx} fails targets {sorted(ts)}")
