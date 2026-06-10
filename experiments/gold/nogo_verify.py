r"""Verification bench for the strengthened no-go (ogdoad OPEN.md problem 1).

Checks, on concrete game-built Gold instances:
  L1  affine stabilizer of {Q=0} is exactly AO(Q) = affine isometries of Q:
      translation parts = singular vectors compensated by Sp\O linear parts;
      pure translations and pure linear maps outside O(Q) are excluded;
      |AO(Q)| = |Z| * |O(Q)| = |Sp(B)|
  L2  transvection criterion: T_v in O(Q)  <=>  Q(v)=1
  L3  every 3-dim subspace contains a nonzero singular vector (Chevalley-Warning)
      and anisotropic planes exist (tightness of the t = 2r-2 escape hatch)
  L4  O(Q)-orbitals on V x V coincide with Gram classes (Q(u),Q(v),B(u,v)) + flags
  L5  refinement torsor: Arf(Q + B(c,.)) = Arf(Q) + Q(c); orbit sizes = Sp:O index
  L6  Theorem D class: every f(d,B(v,d))-gated flip rule is automatically
      undirected; loopy Loss-set is an affine subspace; |quadric| is never a
      power of 2 (r >= 2) so the class can never hit {Q=0}
  L7  commutative-monoid obstruction mechanism on the actual R8 quotient:
      squaring is an endomorphism (trivial polarization)
Cross-checks against repo-documented numbers: |Sp(4,2)|=720, Gold(4,1) zero
count 4 = R(B), bent counts on F_16/F_256, |O+(4,2)|=72 / |O-(4,2)|=120.
"""
from functools import lru_cache
from itertools import product, combinations
import random

# ----------------------------------------------------------------- nim arithmetic
@lru_cache(maxsize=None)
def nmul(a, b):
    if a < b:
        a, b = b, a
    if b == 0:
        return 0
    if b == 1:
        return a
    F = 2
    while F * F <= a:
        F = F * F
    a1, a0 = divmod(a, F)
    b1, b0 = divmod(b, F)
    a1b1 = nmul(a1, b1)
    cross = nmul(a1, b0) ^ nmul(a0, b1)
    return ((a1b1 ^ cross) * F) ^ nmul(a0, b0) ^ nmul(a1b1, F >> 1)

assert nmul(2, 2) == 3 and nmul(2, 3) == 1 and nmul(3, 3) == 2
assert nmul(4, 4) == 6 and nmul(2, 4) == 8

def trace(x, m):
    acc, t = x, x
    for _ in range(m - 1):
        t = nmul(t, t)
        acc ^= t
    assert acc in (0, 1)
    return acc

def frob(x, a):
    for _ in range(a):
        x = nmul(x, x)
    return x

# field sanity: F_16 = {0..15} closed, every nonzero invertible
F16 = list(range(16))
assert all(nmul(x, y) < 16 for x in F16 for y in F16)
assert all(any(nmul(x, y) == 1 for y in range(1, 16)) for x in range(1, 16))
print("[ok] nim arithmetic: F_16 is a field; small products match known table")

# ----------------------------------------------------------------- forms
def make_gold(lam, a, m):
    Q = [trace(nmul(lam, nmul(v, frob(v, a))), m) for v in range(1 << m)]
    def B(u, v):
        return Q[u ^ v] ^ Q[u] ^ Q[v]
    return Q, B

# repo cross-check: plain Gold (4,1): |{Q=0}|=4, rank B = 2, R(B) = {Q=0}
m = 4
Q41, B41 = make_gold(1, 1, 4)
Z41 = [v for v in range(16) if Q41[v] == 0]
RB41 = [v for v in range(16) if all(B41(v, u) == 0 for u in range(16))]
gram41 = [[B41(1 << i, 1 << j) for j in range(4)] for i in range(4)]
def f2rank(rows):
    rows = [int("".join(map(str, r)), 2) if isinstance(r, list) else r for r in rows]
    rk = 0
    for bit in range(16):
        piv = next((i for i, r in enumerate(rows) if (r >> bit) & 1), None)
        if piv is None:
            continue
        p = rows.pop(piv)
        rows = [r ^ p if (r >> bit) & 1 else r for r in rows]
        rk += 1
    return rk
assert len(Z41) == 4 and sorted(Z41) == sorted(RB41)
assert f2rank(gram41) == 2
print("[ok] repo cross-check Gold(4,1): |{Q=0}|=4 = R(B), rank B = 2")

# bent witness on F_16: Tr(lam * v^3), expect 2(2^4-1)/3 = 10 bent lambdas
bents = []
for lam in range(1, 16):
    Q, B = make_gold(lam, 1, 4)
    z = sum(1 for v in range(16) if Q[v] == 0)
    if z in (6, 10):
        bents.append((lam, z))
assert len(bents) == 10
lam, zc = bents[0]
Q, B = make_gold(lam, 1, 4)
ZSET = frozenset(v for v in range(16) if Q[v] == 0)
ARF = 0 if zc == 10 else 1   # 2^{2r-1} + (-1)^Arf 2^{r-1}, r=2: 10 / 6
print(f"[ok] bent witness lam={lam}: |Z|={zc}, Arf={ARF}; 10 bent lambdas (classical count)")

# ----------------------------------------------------------------- GL / Sp / O
def apply(cols, v):
    out = 0
    for i in range(4):
        if (v >> i) & 1:
            out ^= cols[i]
    return out

GL = []
for cols in product(range(16), repeat=4):
    if f2rank(list(cols)) == 4:
        GL.append(cols)
assert len(GL) == 20160
E = [1, 2, 4, 8]
SP = [g for g in GL
      if all(B(apply(g, E[i]), apply(g, E[j])) == B(E[i], E[j])
             for i in range(4) for j in range(i + 1, 4))]
assert len(SP) == 720
OQ = [g for g in GL if all(Q[apply(g, v)] == Q[v] for v in range(16))]
stab_set = [g for g in GL if frozenset(apply(g, v) for v in ZSET) == ZSET]
assert sorted(OQ) == sorted(stab_set), "setwise GL-stabilizer of {Q=0} != O(Q)"
assert all(g in SP for g in OQ), "O(Q) not inside Sp(B)?!"
expected_O = 72 if ARF == 0 else 120
assert len(OQ) == expected_O
print(f"[ok] L1a: Stab_GL({{Q=0}}) = O(Q), |O(Q)|={len(OQ)} (= {'O+' if ARF==0 else 'O-'}(4,2)), O(Q) <= Sp(B), |Sp|=720")

def transvect(v):
    return tuple(E[i] ^ (v if B(E[i], v) else 0) for i in range(4))

# L1b: affine stabilizer of {Q=0} is the affine isometry group AO(Q):
#   (g,c) stabilizes Z  <=>  Q(gx^c) == Q(x) for all x.
# Translation parts realize exactly Z (the singular vectors); pure translations
# (g=id, c!=0) never stabilize; pure linear stabilizers are exactly O(Q).
affine_stab = []
for g in GL:
    perm = [apply(g, v) for v in range(16)]
    for c in range(16):
        if frozenset(perm[v] ^ c for v in ZSET) == ZSET:
            affine_stab.append((g, c))
            # set-stabilizer <=> function isometry (F_2-valued, same zero set => equal)
            assert all(Q[perm[x] ^ c] == Q[x] for x in range(16))
            assert Q[c] == 0  # translation part is singular
            assert g in set(SP)  # linear part is symplectic
assert len(affine_stab) == len(ZSET) * len(OQ) == len(SP)  # |AO| = |Z||O| = |Sp|
assert all(c == 0 for (g, c) in affine_stab if g == (1, 2, 4, 8))  # no pure translation
assert sorted(g for (g, c) in affine_stab if c == 0) == sorted(OQ)  # pure linear = O(Q)
trans_parts = {c for (g, c) in affine_stab}
assert trans_parts == set(ZSET)  # translation parts = singular vectors exactly
# singular transvections enter AO only with their forced twist c = v, never linearly
for v in range(1, 16):
    if Q[v] == 0:
        Tv = transvect(v)
        assert (Tv, v) in set(affine_stab) and (Tv, 0) not in set(affine_stab)
print(f"[ok] L1b: Stab_AGL({{Q=0}}) = AO(Q), order {len(affine_stab)} = |Z|*|O(Q)| = |Sp(B)|;")
print("        pure translations excluded; pure linear stabilizers = O(Q) exactly;")
print("        singular transvections T_v enter only as the twisted x -> T_v x + v")

# ----------------------------------------------------------------- L2 transvections
for v in range(1, 16):
    Tv = transvect(v)
    assert Tv in set(SP), "transvection not symplectic"
    in_O = all(Q[apply(Tv, x)] == Q[x] for x in range(16))
    assert in_O == (Q[v] == 1), f"transvection criterion fails at v={v}"
print("[ok] L2: T_v in Sp(B) always; T_v in O(Q) <=> Q(v)=1  (all 15 v)")

# ----------------------------------------------------------------- L3 subspaces
def span(vs):
    s = {0}
    for v in vs:
        s |= {x ^ v for x in s}
    return frozenset(s)
dim3 = set()
for t in combinations(range(1, 16), 3):
    s = span(t)
    if len(s) == 8:
        dim3.add(s)
assert len(dim3) == 15  # number of 3-dim subspaces of F_2^4 = gaussian [4,3]_2 = 15
assert all(any(Q[v] == 0 for v in s if v) for s in dim3)
planes = set()
for t in combinations(range(1, 16), 2):
    s = span(t)
    if len(s) == 4:
        planes.add(s)
assert len(planes) == 35
aniso = [s for s in planes if all(Q[v] == 1 for v in s if v)]
print(f"[ok] L3: every 3-dim subspace has a nonzero singular vector; anisotropic planes exist ({len(aniso)} of 35)")

# ----------------------------------------------------------------- L4 pair orbitals
pairs = [(u, v) for u in range(16) for v in range(16)]
idx = {p: i for i, p in enumerate(pairs)}
parent = list(range(len(pairs)))
def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x
def union(x, y):
    x, y = find(x), find(y)
    if x != y:
        parent[x] = y
for g in OQ:
    perm = [apply(g, v) for v in range(16)]
    for (u, v) in pairs:
        union(idx[(u, v)], idx[(perm[u], perm[v])])
orbit_of = {p: find(idx[p]) for p in pairs}
def gram_class(u, v):
    return (Q[u], Q[v], B(u, v), u == 0, v == 0, u == v)
orbits = {}
for p in pairs:
    orbits.setdefault(orbit_of[p], set()).add(gram_class(*p))
classes = {}
for p in pairs:
    classes.setdefault(gram_class(*p), set()).add(orbit_of[p])
split = {c: o for c, o in classes.items() if len(o) > 1}
fused = {o: c for o, c in orbits.items() if len(c) > 1}
assert not fused, "one orbit with two Gram classes?!"
print(f"[ok] L4: O(Q)-orbitals on V x V == Gram classes (Q(u),Q(v),B(u,v))+flags: "
      f"{len(set(orbit_of.values()))} orbitals, {len(classes)} classes, splits={len(split)}")

# ----------------------------------------------------------------- L5 torsor
def arf_from_zero_count(z, n, r):
    if z == (1 << (n - 1)) + (1 << (r - 1)):
        return 0
    if z == (1 << (n - 1)) - (1 << (r - 1)):
        return 1
    return None
for c in range(16):
    Qc = [Q[x] ^ B(c, x) for x in range(16)]
    zc2 = sum(1 for x in range(16) if Qc[x] == 0)
    arf_c = arf_from_zero_count(zc2, 4, 2)
    assert arf_c == ARF ^ Q[c], f"torsor Arf shift fails at c={c}"
counts = {0: 0, 1: 0}
for c in range(16):
    counts[ARF ^ Q[c]] += 1
assert counts[0] in (10, 6) and counts[0] + counts[1] == 16
assert 720 // 72 == 10 and 720 // 120 == 6
print(f"[ok] L5: refinement torsor Arf(Q+B(c,.)) = Arf(Q)+Q(c); class sizes {counts} = Sp:O indices (10,6)")

# ----------------------------------------------------------------- L6 Theorem D class
# rule: flip direction d at v legal iff f(d, B(v,d)); behaviors per d:
#   0 never, 1 iff B=0, 2 iff B=1, 3 always
def loss_set_of(beh):
    def legal(v, d):
        b = beh[d]
        if b == 0:
            return False
        if b == 3:
            return True
        return (B(v, d) == 1) == (b == 2)
    # automatic symmetry check: B(v^d, d) == B(v, d)
    for v in range(16):
        for d in range(1, 16):
            assert B(v ^ d, d) == B(v, d)
            assert legal(v, d) == legal(v ^ d, d)
    # undirected loopy: Loss = isolated vertices, everything else Draw/Win;
    # standard fixpoint degenerates to: Loss = isolated (see loopy_quadric.rs)
    return frozenset(v for v in range(16)
                     if not any(legal(v, d) for d in range(1, 16)))

def is_affine(s):
    if not s:
        return True
    s0 = next(iter(s))
    sh = {x ^ s0 for x in s}
    return all((x ^ y) in sh for x in sh for y in sh)

rng = random.Random(0xA9)
tested = 0
for trial in range(400):
    beh = {d: rng.randrange(4) for d in range(1, 16)}
    L = loss_set_of(beh)
    assert is_affine(L), f"non-affine Loss set: {sorted(L)}"
    assert len(L) != zc or frozenset(L) != ZSET
    assert bin(len(L)).count("1") <= 1  # 0 or power of 2
    tested += 1
# uniform behaviors too (includes the repo's symmetric-B rule beh=2: Loss=R(B)={0})
for u in range(4):
    L = loss_set_of({d: u for d in range(1, 16)})
    assert is_affine(L)
assert loss_set_of({d: 2 for d in range(1, 16)}) == frozenset({0})  # bent => R(B)={0}
print(f"[ok] L6: {tested}+4 f(d,B(v,d))-gated rules: all auto-undirected, Loss-set always affine")
print(f"        (|{{Q=0}}| = {zc} = 2*odd, never a power of 2  =>  class can never hit the quadric)")

# ----------------------------------------------------------------- L6' on F_256
m8 = 8
lam8 = None
for cand in range(1, 256):
    Q8 = [trace(nmul(cand, nmul(v, frob(v, 1))), m8) for v in range(256)]
    z8 = sum(1 for v in range(256) if Q8[v] == 0)
    if z8 in (120, 136):
        lam8, Q8z, arf8 = cand, z8, (0 if z8 == 136 else 1)
        break
assert lam8 is not None
Q8 = [trace(nmul(lam8, nmul(v, frob(v, 1))), m8) for v in range(256)]
def B8(u, v):
    return Q8[u ^ v] ^ Q8[u] ^ Q8[v]
# transvection criterion at m=8 (no group enumeration needed)
for v in range(1, 256):
    def Tv(x):
        return x ^ (v if B8(x, v) else 0)
    in_O = all(Q8[Tv(x)] == Q8[x] for x in range(256))
    assert in_O == (Q8[v] == 1)
# every 3-dim subspace has a nonzero singular vector (random sample)
rng = random.Random(7)
for _ in range(500):
    vs = rng.sample(range(1, 256), 3)
    s = span(vs)
    if len(s) == 8:
        assert any(Q8[v] == 0 for v in s if v)
print(f"[ok] L2'/L3' on F_256 bent lam={lam8} (|Z|={Q8z}, Arf={arf8}): transvection criterion all 255 v; CW sample 500")

# ----------------------------------------------------------------- L7 R8 monoid
BASE = ["1", "b", "b2", "c"]
_MN = {
    ("1", "1"): (0, "1"), ("1", "b"): (0, "b"), ("1", "b2"): (0, "b2"), ("1", "c"): (0, "c"),
    ("b", "b"): (0, "b2"), ("b", "b2"): (0, "b"), ("b", "c"): (1, "b"),
    ("b2", "b2"): (0, "b2"), ("b2", "c"): (1, "b2"),
    ("c", "c"): (0, "b2"),
}
def _mn(x, y):
    return _MN.get((x, y)) or _MN[(y, x)]
ELEMS = [(i, mm) for i in (0, 1) for mm in BASE]
def mul(x, y):
    (i, mm), (j, nn) = x, y
    extra, base = _mn(mm, nn)
    return ((i + j + extra) % 2, base)
sq = {x: mul(x, x) for x in ELEMS}
assert all(sq[mul(x, y)] == mul(sq[x], sq[y]) for x in ELEMS for y in ELEMS)
print("[ok] L7: on R8 (smallest wild misere quotient) squaring IS an endomorphism")
print("        => its polarization B(x,y) := s(xy)s(x)^-1 s(y)^-1 (where defined) is trivial;")
print("        the commutative world can only carry the split refinement.")

print("\nALL CHECKS PASSED")
