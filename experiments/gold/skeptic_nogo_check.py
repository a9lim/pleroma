"""Independent adversarial re-verification of the no-go attack.

Different code path from /tmp/nogo_verify.py: coordinate quadratic forms (no nim
arithmetic), BOTH Arf classes exhaustively at dim 4, plus an end-to-end retrograde
P-set check of the Theorem C transvection step, and the Draw-set gap in Theorem D.
"""
from itertools import product, combinations
import random

N = 4
V = range(16)
bit = lambda x, i: (x >> i) & 1

def Q_arf0(x):  # x0x1 + x2x3, |Z| = 10
    return (bit(x,0)&bit(x,1)) ^ (bit(x,2)&bit(x,3))
def Q_arf1(x):  # x0^2 + x0x1 + x1^2 + x2x3, |Z| = 6
    return bit(x,0) ^ (bit(x,0)&bit(x,1)) ^ bit(x,1) ^ (bit(x,2)&bit(x,3))

def polar(Q):
    return lambda u, v: Q(u^v) ^ Q(u) ^ Q(v)

def apply(cols, v):
    out = 0
    for i in range(N):
        if (v >> i) & 1:
            out ^= cols[i]
    return out

def rank(rows):
    rows = list(rows); rk = 0
    for b in range(N):
        piv = next((i for i,r in enumerate(rows) if (r>>b)&1), None)
        if piv is None: continue
        p = rows.pop(piv); rows = [r^p if (r>>b)&1 else r for r in rows]; rk += 1
    return rk

GL = [c for c in product(V, repeat=N) if rank(c) == 4]
assert len(GL) == 20160

def spanset(vs):
    s={0}
    for v in vs: s|={x^v for x in s}
    return frozenset(s)

def loopy_solve(adj):
    out = {v: None for v in V}
    changed = True
    while changed:
        changed = False
        for v in V:
            if out[v] is not None: continue
            if all(out[w]=="W" for w in adj[v]):       # incl. terminal
                out[v]="L"; changed=True
            elif any(out[w]=="L" for w in adj[v]):
                out[v]="W"; changed=True
    return {v:(o if o else "D") for v,o in out.items()}

for tag, Q, zexp, oexp in (("Arf0", Q_arf0, 10, 72), ("Arf1", Q_arf1, 6, 120)):
    B = polar(Q)
    Z = frozenset(v for v in V if Q(v) == 0)
    assert len(Z) == zexp
    assert all(any(B(v,u) for u in V) for v in V if v), "polar form degenerate"
    E = [1,2,4,8]
    SP = [g for g in GL if all(B(apply(g,a),apply(g,b))==B(a,b) for a in E for b in E)]
    assert len(SP) == 720
    O = [g for g in GL if all(Q(apply(g,v))==Q(v) for v in V)]
    assert len(O) == oexp

    # --- Theorem A: Stab_AGL(Z) = AO(Q), exhaustive over all 322560 affine maps
    AO = []
    for g in GL:
        gi = [apply(g,v) for v in V]
        for c in V:
            if frozenset(gi[v]^c for v in Z) == Z:
                AO.append((g,c))
                assert all(Q(gi[x]^c)==Q(x) for x in V)   # set-stab => isometry
                assert Q(c)==0 and g in set(SP)
    assert len(AO) == len(Z)*len(O) == 720
    assert not any(c for g,c in AO if g == (1,2,4,8))            # no pure translation
    assert sorted(g for g,c in AO if c==0) == sorted(O)          # pure linear = O(Q)
    assert {c for g,c in AO} == set(Z)                           # translations = Z exactly

    # --- Theorem B core: O(Q)-orbitals on VxV = Gram classes + flags
    gram = lambda u,v: (Q(u),Q(v),B(u,v),u==0,v==0,u==v)
    orbit_reps = {}
    seen = set()
    for p in product(V,V):
        if p in seen: continue
        orb, stack = set(), [p]
        while stack:
            (u,v)=stack.pop()
            if (u,v) in orb: continue
            orb.add((u,v))
            for g in O:
                q2=(apply(g,u),apply(g,v))
                if q2 not in orb: stack.append(q2)
        seen |= orb
        gs = {gram(*x) for x in orb}
        assert len(gs)==1, "orbit spans two Gram classes"
        key = gs.pop()
        assert key not in orbit_reps, f"Gram class split into two orbits: {key}"
        orbit_reps[key]=len(orb)
    n_orb = len(orbit_reps)

    # --- Theorem C core: transvections + CW + exact escape boundary at t=2r-2
    T = lambda v: tuple(E[i]^(v if B(E[i],v) else 0) for i in range(N))
    for v in range(1,16):
        assert T(v) in set(SP)
        assert (T(v) in set(O)) == (Q(v)==1)
    d3 = {s for s in (spanset(t) for t in combinations(range(1,16),3)) if len(s)==8}
    assert len(d3)==15 and all(any(Q(v)==0 for v in s if v) for s in d3)
    planes = {s for s in (spanset(t) for t in combinations(range(1,16),2)) if len(s)==4}
    assert len(planes)==35
    escape, killed = 0, 0
    for W in planes:                       # W = span of the t=2r-2 constants
        Wp = frozenset(v for v in V if all(B(v,w)==0 for w in W))  # W-perp
        assert len(Wp)==4
        fixW = [g for g in SP if all(apply(g,w)==w for w in W)]
        if any(Q(v)==0 for v in Wp if v):
            v = next(v for v in Wp if v and Q(v)==0)
            assert T(v) in set(fixW) and T(v) not in set(O)        # bare transvection kills
            killed += 1
        else:
            assert all(g in set(O) for g in fixW), "escape claim fails"
            assert len(fixW)==6   # Sp(2,2)=GL(2,2)=S3
            escape += 1
    assert escape>0 and escape+killed==35

    # --- Theorem D incl. the Draw-set gap: undirected loopy outcomes
    rng = random.Random(1)
    for _ in range(200):
        beh = {d: rng.randrange(4) for d in range(1,16)}
        def legal(v,d):
            b=beh[d]
            return b==3 or (b in (1,2) and (B(v,d)==1)==(b==2))
        adj = {v:[v^d for d in range(1,16) if legal(v,d)] for v in V}
        assert all(v in adj[w] for v in V for w in adj[v])   # undirected
        out = loopy_solve(adj)
        L  = {v for v in V if out[v]=="L"}
        W_ = {v for v in V if out[v]=="W"}
        D  = {v for v in V if out[v]=="D"}
        assert L == {v for v in V if not adj[v]}    # Loss = isolated
        assert not W_                               # Win empty
        assert L != Z and D != Z and (L|D) != Z     # NO target set hits the quadric
    print(f"[ok] {tag}: |Z|={zexp} |O|={oexp} |AO|=720, translations=Z, "
          f"orbitals=Gram({n_orb}), transvection crit, CW, "
          f"escape-planes={escape}/35, loopy Loss=isolated & no target hits Z")

# --- end-to-end Theorem C spot check: t=1 B-oracle rules, retrograde-solved;
#     the singular transvection is a graph automorphism, outcome classes invariant
Q, B = Q_arf0, polar(Q_arf0)
Z = frozenset(v for v in V if Q(v)==0)
c1 = 3
rng = random.Random(9)
for trial in range(50):
    f = {bits: rng.randrange(2) for bits in product((0,1),repeat=3)}
    def legal(u,w):
        if u==w: return False
        return f[(B(u,w)&1, B(u,c1)&1, B(w,c1)&1)]==1
    adj = {v:[w for w in V if legal(v,w)] for v in V}
    v = next(v for v in V if v and Q(v)==0 and B(v,c1)==0)   # singular, perp to c1
    Tv = lambda x: x ^ (v if B(x,v) else 0)
    assert all(sorted(Tv(x) for x in adj[u]) == sorted(adj[Tv(u)]) for u in V)
    out = loopy_solve(adj)
    for lab in ("L","W","D"):
        S = {x for x in V if out[x]==lab}
        assert {Tv(x) for x in S} == S          # outcome classes T_v-invariant
        assert S != Z                            # and never the quadric
print("[ok] Theorem C end-to-end: 50 random t=1 B-oracle rules; singular transvection")
print("     is a graph automorphism, all outcome classes T_v-invariant, none equals Z")
print("\nINDEPENDENT CHECKS ALL PASS")
