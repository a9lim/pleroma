"""Verify the 'equivariance sandwich' at m=4 (dim-4 nondegenerate case).

Claims checked exhaustively over V = F_2^4, B = the standard symplectic form
(pairs (0,1),(2,3)), for BOTH Arf classes of quadratic refinement Q:

  (a) |Sp(B)| = 720, transitive on V\{0}  (the Tier-1 no-go input);
  (b) O(Q) = setwise stabilizer of {Q=0} in Sp(B), with |O^+|=72, |O^-|=120;
  (c) ORBITAL LEMMA: the O(Q)-orbits on ordered pairs V x V are EXACTLY the
      fibers of the invariant tuple (Q(v), Q(w), B(v,w), [v=0], [w=0], [v=w]);
      hence every O(Q)-invariant move relation is a Boolean function of
      Q-evaluations and one B-evaluation -- i.e. semantically a Q-evaluator;
  (d) MAXIMALITY: for every g in Sp(B) \ O(Q), <O(Q), g> = Sp(B);
      hence there is NO group strictly between O(Q) and Sp(B).
"""
import itertools
import random

M = 4
N = 1 << M
PAIRS = [(0, 1), (2, 3)]


def B(u, v):
    acc = 0
    for (i, j) in PAIRS:
        acc ^= ((u >> i) & (v >> j) & 1) ^ ((u >> j) & (v >> i) & 1)
    return acc


def make_Q(qd):
    def Q(v):
        acc = sum((qd[i] & (v >> i)) & 1 for i in range(M))
        for (i, j) in PAIRS:
            acc ^= (v >> i) & (v >> j) & 1
        return acc & 1
    return Q


def invertible(cols):
    rows = list(cols)
    r = 0
    for bit in range(M):
        piv = next((k for k in range(r, len(rows)) if (rows[k] >> bit) & 1), None)
        if piv is None:
            continue
        rows[r], rows[piv] = rows[piv], rows[r]
        for k in range(len(rows)):
            if k != r and (rows[k] >> bit) & 1:
                rows[k] ^= rows[r]
        r += 1
    return r == M


def to_perm(cols):
    out = []
    for v in range(N):
        x = 0
        for i in range(M):
            if (v >> i) & 1:
                x ^= cols[i]
        out.append(x)
    return tuple(out)


basis = [1 << i for i in range(M)]
sp = []
for cols in itertools.product(range(N), repeat=M):
    if not invertible(cols):
        continue
    if all(B(cols[i], cols[j]) == B(basis[i], basis[j])
           for i in range(M) for j in range(i + 1, M)):
        sp.append(to_perm(cols))
print(f"(a) |Sp(B)| = {len(sp)} (expect 720)")
orbit0 = set(g[1] for g in sp)
print(f"    orbit of e_0: {len(orbit0)} vectors (expect {N-1}: transitive on V\\0)")

sp_set = set(sp)


def closure_size(gens):
    seen = {tuple(range(N))}
    frontier = [tuple(range(N))]
    gens = [g for g in gens]
    while frontier:
        nxt = []
        for h in frontier:
            for g in gens:
                hg = tuple(h[g[v]] for v in range(N))
                if hg not in seen:
                    seen.add(hg)
                    nxt.append(hg)
        frontier = nxt
    return seen


def small_gens(group):
    rng = random.Random(0xA9)
    glist = list(group)
    for k in (2, 3, 4):
        for _ in range(60):
            gens = rng.sample(glist, k)
            if len(closure_size(gens)) == len(group):
                return gens
    raise RuntimeError("no small generating set found")


for name, qd, expect in [("Arf 0 (O+)", (0, 0, 0, 0), 72),
                         ("Arf 1 (O-)", (1, 1, 0, 0), 120)]:
    Q = make_Q(qd)
    OQ = [g for g in sp if all(Q(g[v]) == Q(v) for v in range(N))]
    # sanity: zero-set stabilizer equals function stabilizer over F_2
    zeros = frozenset(v for v in range(N) if Q(v) == 0)
    OQ_set = [g for g in sp if frozenset(g[v] for v in zeros) == zeros]
    print(f"\n(b) {name}: |O(Q)| = {len(OQ)} (expect {expect}); "
          f"zero-set stabilizer == form stabilizer: {OQ == OQ_set}")

    # (c) orbital lemma
    def invariant(v, w):
        return (Q(v), Q(w), B(v, w), v == 0, w == 0, v == w)

    pair_orbit = {}
    for v in range(N):
        for w in range(N):
            if (v, w) in pair_orbit:
                continue
            orb = set()
            stack = [(v, w)]
            while stack:
                p = stack.pop()
                if p in orb:
                    continue
                orb.add(p)
                for g in OQ:
                    stack.append((g[p[0]], g[p[1]]))
            for p in orb:
                pair_orbit[p] = (v, w)
    orbits = {}
    for p, rep in pair_orbit.items():
        orbits.setdefault(rep, set()).add(p)
    classes = {}
    for v in range(N):
        for w in range(N):
            classes.setdefault(invariant(v, w), set()).add((v, w))
    match = sorted(map(sorted, orbits.values())) == sorted(map(sorted, classes.values()))
    print(f"(c) orbital lemma: #orbits on VxV = {len(orbits)}, "
          f"#invariant classes = {len(classes)}, partitions equal: {match}")

    # (d) maximality
    gens = small_gens(OQ)
    OQ_frozen = set(OQ)
    bad = 0
    for g in sp:
        if g in OQ_frozen:
            continue
        if len(closure_size(gens + [g])) != len(sp):
            bad += 1
    print(f"(d) maximality: <O(Q),g> = Sp(B) for all {len(sp)-len(OQ)} "
          f"g outside O(Q): {bad == 0} (failures: {bad})")
