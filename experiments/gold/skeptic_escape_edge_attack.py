"""SKEPTIC ATTACK on N3 (criterion angle, round 2 final-skeptic pass).

Construction: the pending-marker clock PLUS one q-blind 'escape' edge from
every clean position with x != 0 to a fixed always-N gadget (gamma -> t,
t terminal).  The escape move is strategically dead (it hands the opponent
an N-position), so the rule is morally still a clock.  Claim to verify:

  N1: Loss(iota(x)) == {Q_q(x)=0} for EVERY refinement q (P-set unchanged
      by the dead edge), every tested form.
  N2: (1,1)-local by construction (only the pending-completion edge queries
      the oracle, one weight-1 query; the escape edge is q-blind).
  N3: every loaded N-position is simultaneously outcome-critical (moves to
      both P pendings and the N gadget) and form-live (refinement flip).
  N3+: the critical-AND-live set is reached from a constant fraction of
      loaded positions (they ARE loaded positions).

If all pass, the sharp conjecture as stated is trivially TRUE, witnessed by
a clock wearing a dead edge -- N3 does not carry the content claimed.
"""
import random, sys
from functools import lru_cache
sys.setrecursionlimit(100000)

# nim arithmetic, identical to the artifact's (validated against repo pins)
@lru_cache(maxsize=None)
def nm(a, b):
    if a < b: a, b = b, a
    if b < 2: return a * b
    F = 2
    while F * F <= a: F *= F
    a1, a0 = divmod(a, F); b1, b0 = divmod(b, F)
    c  = nm(a1, b1)
    hi = c ^ nm(a1, b0) ^ nm(a0, b1)
    lo = nm(a0, b0) ^ nm(c, F // 2)
    return hi * F ^ lo
assert nm(2,2)==3 and nm(2,4)==8 and nm(16,16)==24

def frob(x,k,m):
    for _ in range(k): x = nm(x,x)
    return x
def tr(x,m):
    t,y = 0,x
    for _ in range(m): t ^= y; y = nm(y,y)
    assert t in (0,1); return t
def gold(x,a,m): return tr(nm(x, frob(x,a,m)), m)

def form_data(m,a):
    q = [gold(1<<i,a,m) for i in range(m)]
    B = [[0]*m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            if i != j: B[i][j] = gold((1<<i)^(1<<j),a,m) ^ q[i] ^ q[j]
    return q,B

def Qr(x,q,B,m):
    s=0; idx=[i for i in range(m) if (x>>i)&1]
    for i in idx: s ^= q[i]
    for ii in range(len(idx)):
        for jj in range(ii+1,len(idx)): s ^= B[idx[ii]][idx[jj]]
    return s

def Bv(v,d,B,m):
    s=0
    for i in range(m):
        if (v>>i)&1:
            for j in range(m):
                if (d>>j)&1: s ^= B[i][j]
    return s

GAMMA = ('gamma',); TEE = ('t',); BOT = ('bot',)

def attack_game(m, qvec, B):
    def succ(p):
        if p in (TEE, BOT): return []
        if p == GAMMA: return [TEE]                      # gamma is always N
        if p[0] == 'c':
            _, x, eps = p
            if x == 0:
                return [BOT] if eps == 1 else []
            ss = [('p', x, eps, i) for i in range(m) if (x>>i)&1]
            ss.append(GAMMA)                             # q-blind escape edge
            return ss
        _, x, eps, i = p
        # the ONLY oracle access: one weight-1 query q_i (plus public B)
        return [('c', x^(1<<i), eps ^ qvec[i] ^ Bv(x,1<<i,B,m))]
    out = {}
    def rec(p):
        if p in out: return out[p]
        out[p] = False
        out[p] = any(not rec(s) for s in succ(p))
        return out[p]
    return succ, rec

def check(m, a):
    qg, B = form_data(m, a)
    random.seed(1)
    refs = [qg] + [[random.randint(0,1) for _ in range(m)] for _ in range(5)]
    base_out = None
    live = set()
    for qv in refs:
        succ, rec = attack_game(m, qv, B)
        # N1: normal-play Loss slice == {Q_qv = 0}
        n1 = all((not rec(('c',x,0))) == (Qr(x,qv,B,m)==0) for x in range(1<<m))
        assert n1, ("N1 FAILS", m, a, qv)
        outs = {x: rec(('c',x,0)) for x in range(1<<m)}
        if base_out is None:
            base_out, base_succ, base_rec = outs, succ, rec
        else:
            live |= {x for x in range(1<<m) if outs[x] != base_out[x]}
    # N3 under the base (Gold) refinement: critical loaded positions
    critical = set()
    for x in range(1<<m):
        p = ('c', x, 0)
        cls = {base_rec(s) for s in base_succ(p)}
        if len(cls) == 2: critical.add(x)
    both = critical & live
    frac = len(both) / (1 << m)
    print(f"  ({m},{a}): N1 holds for Gold + 5 random refinements; "
          f"critical={len(critical)}, form-live={len(live)}, BOTH={len(both)} "
          f"({frac:.0%} of loaded positions)")
    assert both, "attack failed: conjunction empty"
    # sanity: every critical-and-live position is loaded (reachability trivial)
    return len(both)

print("escape-edge clock attack (clock + q-blind dead edge):")
for (m,a) in [(8,1),(8,2)]:
    check(m,a)
print("ATTACK SUCCEEDS: N1+N2(1,1)+N3+N3+ all pass; rule is morally a clock.")
print("=> the sharp conjecture as stated is trivially TRUE; N3 is gamed.")
