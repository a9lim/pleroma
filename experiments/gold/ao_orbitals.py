"""Obstruction arm, Lemma 1 + Lemma 2 verification at m=4.

Lemma 1 (AO(Q)-orbital coarsening): for nondegenerate Q on F_2^4 (both Arf
classes), the orbits of AO(Q) = Stab_AGL(Z) on V x V are exactly the fibers of
(Q(u), Q(w), [u==w]).  Hence an AO(Q)-invariant move relation's legality factors
through endpoint Q-values alone: a pure two-point evaluator.
Also: the canonical evaluator rule R_T (legal u->w iff (Q(u),Q(w)) in T) is
AO(Q)-invariant -- so ANY invariance requirement G <= AO(Q) admits the evaluator.

Lemma 2 (Inn(E)/conjugation vacuity): conjugation orbits on E x E leave the
full V x V data free; every pullback of every V-rule is conjugation-equivariant,
including the evaluator.  Predicted orbit count on E x E for m=4: 304.
"""
import itertools

# coordinate quadratic forms on F_2^4 (no nim arithmetic -- independent route)
def Q_arf0(x):  # x1x2 + x3x4
    b = [(x >> i) & 1 for i in range(4)]
    return (b[0] & b[1]) ^ (b[2] & b[3])
def Q_arf1(x):  # x1x2 + x1 + x2 + x3x4   (x^2 terms = diagonal)
    b = [(x >> i) & 1 for i in range(4)]
    return (b[0] & b[1]) ^ b[0] ^ b[1] ^ (b[2] & b[3])

def mat_apply(M, x):
    y = 0
    for i in range(4):
        if (x >> i) & 1:
            y ^= M[i]
    return y

# all invertible 4x4 matrices over F_2 (columns as masks)
def gl4():
    mats = []
    for c0 in range(1, 16):
        for c1 in range(1, 16):
            if c1 == c0: continue
            sp01 = {0, c0, c1, c0 ^ c1}
            for c2 in range(1, 16):
                if c2 in sp01: continue
                sp = {a ^ b for a in sp01 for b in (0, c2)}
                for c3 in range(1, 16):
                    if c3 in sp: continue
                    mats.append((c0, c1, c2, c3))
    return mats

GL = gl4()
print(f"|GL(4,2)| = {len(GL)} (expect 20160)")

for name, Qf in (("Arf0", Q_arf0), ("Arf1", Q_arf1)):
    Q = [Qf(x) for x in range(16)]
    Z = [x for x in range(16) if Q[x] == 0]
    AO = []
    for M in GL:
        # images under M for all x
        img = [mat_apply(M, x) for x in range(16)]
        for c in range(16):
            if all(Q[img[x] ^ c] == Q[x] for x in range(16)):
                AO.append((M, c))
    trans_parts = sorted({c for _, c in AO})
    print(f"{name}: |Z|={len(Z)} |AO(Q)|={len(AO)} translation parts == Z: "
          f"{trans_parts == sorted(Z)}")
    # orbit partition of V x V under AO(Q)
    # union-find
    parent = list(range(256))
    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]; a = parent[a]
        return a
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb: parent[ra] = rb
    for M, c in AO:
        img = [mat_apply(M, x) ^ c for x in range(16)]
        for u in range(16):
            for w in range(16):
                union(u * 16 + w, img[u] * 16 + img[w])
    orbits = {}
    for p in range(256):
        orbits.setdefault(find(p), []).append(p)
    # predicted fibers of (Q(u), Q(w), [u==w])
    fibers = {}
    for u in range(16):
        for w in range(16):
            fibers.setdefault((Q[u], Q[w], u == w), []).append(u * 16 + w)
    match = sorted(tuple(sorted(o)) for o in orbits.values()) == \
            sorted(tuple(sorted(f)) for f in fibers.values())
    print(f"   orbits on VxV: {len(orbits)}; fibers of (Q(u),Q(w),[u=w]): "
          f"{len(fibers)}; partitions match: {match}")
    # evaluator rule invariance: R = {(u,w): (Q[u],Q[w]) == (1,0)}
    R = {(u, w) for u in range(16) for w in range(16) if (Q[u], Q[w]) == (1, 0)}
    inv = all((img[u] ^ 0, ) for _ in ()) or True
    ok = True
    for M, c in AO:
        img = [mat_apply(M, x) ^ c for x in range(16)]
        if {(img[u], img[w]) for (u, w) in R} != R:
            ok = False; break
    print(f"   canonical evaluator rule AO(Q)-invariant: {ok}")

# ---------- Lemma 2: conjugation orbits on E x E ----------
# E via triangular cocycle for Q_arf1's data on F_2^4 (any nondeg form works)
Q = [Q_arf1(x) for x in range(16)]
q = [Q[1 << i] for i in range(4)]
Bbit = [[Q[(1 << i) ^ (1 << j)] ^ Q[1 << i] ^ Q[1 << j] if i != j else 0
         for j in range(4)] for i in range(4)]
def coc(u, v):  # c(u,v) = sum q_i u_i v_i + sum_{k>j} B_kj u_k v_j
    s = 0
    for i in range(4):
        if (u >> i) & 1 and (v >> i) & 1: s ^= q[i]
    for kk in range(4):
        for j in range(kk):
            if (u >> kk) & 1 and (v >> j) & 1: s ^= Bbit[kk][j]
    return s
def emul(g, h):  # g=(a,u), h=(b,v)
    return ((g[0] ^ h[0] ^ coc(g[1], h[1])), g[1] ^ h[1])
def einv(g):
    # g^{-1} = (a + c(u,u) ... ) solve: g*(b,u) = (0,0) -> b = a + c(u,u)? check
    a, u = g
    b = a ^ coc(u, u)
    assert emul(g, (b, u)) == (0, 0)
    return (b, u)
E = [(a, u) for a in range(2) for u in range(16)]
# squaring map check: g^2 = (Q(u), 0)
assert all(emul(g, g) == (Q[g[1]], 0) for g in E)
# conjugation orbits on E x E
idx = {g: i for i, g in enumerate(E)}
parent = list(range(32 * 32))
def find(a):
    while parent[a] != a:
        parent[a] = parent[parent[a]]; a = parent[a]
    return a
def union(a, b):
    ra, rb = find(a), find(b)
    if ra != rb: parent[ra] = rb
def conj(x, g):
    return emul(emul(x, g), einv(x))
for x in E:
    cmap = {g: conj(x, g) for g in E}
    for g in E:
        for h in E:
            union(idx[g] * 32 + idx[h], idx[cmap[g]] * 32 + idx[cmap[h]])
orbs = {}
for p in range(32 * 32):
    orbs.setdefault(find(p), 0)
    orbs[find(p)] += 1
print(f"\nInn(E) orbits on ExE: {len(orbs)} (predicted 304: "
      f"210 free fibers + 30+30 half-fibers + 30 diag-fiber pairs + 4 central)")
# verify: for independent (gbar,hbar) the full 4-element fiber is ONE orbit
free = 0
for gb in range(1, 16):
    for hb in range(1, 16):
        if hb == gb: continue
        reps = {find(idx[(a, gb)] * 32 + idx[(b, hb)]) for a in (0, 1) for b in (0, 1)}
        if len(reps) == 1: free += 1
print(f"independent-pair fibers that are single orbits: {free}/210 (predict 210)")
# pullback of evaluator rule is conjugation-equivariant
Rv = {(u, w) for u in range(16) for w in range(16) if (Q[u], Q[w]) == (1, 0)}
RE = {(g, h) for g in E for h in E if (g[1], h[1]) in Rv}
ok = all((conj(x, g), conj(x, h)) in RE for (g, h) in RE for x in E)
print(f"pullback of V-evaluator is Inn(E)-equivariant: {ok}")
