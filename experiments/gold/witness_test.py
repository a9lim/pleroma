"""The witness reduction test.

THEOREM (witness reduction, proved in the writeup): the width-k coin-turning
spin-flip rule T_k -- move v -> v^d for wt(d) <= k with msb(d) in supp(v),
legal iff the move flips Q (computable from q_i oracle bits + public B:
DeltaQ(v,d) = B(v,d) + Q(d), Q(d) = sum_{i in d} q_i + sum_{i<j in d} B_ij)
-- has P-set exactly {Q=0} IFF the instance has the k-LOCAL WITNESS property:

   every v with Q(v)=1 admits d, wt(d) <= k, msb(d) in supp(v), Q(v^d)=0.

(<=: ascending induction. Q(v)=0: every legal move flips Q, lands on Q=1 = Win
 by induction => v Loss. Q(v)=1: the witness move lands on Q=0 = Loss => Win.
 =>: a witness-less v with Q(v)=1 has NO legal move => terminal Loss in {Q=1}.)

So the Tier-2 conjecture REDUCES (sufficient direction) to: do bent Gold
components + all their refinements have k-local witnesses for constant k?

This script measures the minimal witness radius k:
  1. sanity: lambda=43, m=8 -- confirm zero diagonal framing (the S1 hits);
  2. m=8, a=1: all 170 bent lambdas, Gold framing: minimal k per instance;
  3. m=8: ALL 256 refinements for every bent lambda (full conjecture instance);
  4. m=16, a=1: first bent lambdas, Gold framing + random refinements.
"""
import random
import sys
import ogdoad as pl


def make_qtab(m, a, lam):
    N = 1 << m
    L = pl.Nimber(lam)
    tab = []
    for v in range(N):
        X = pl.Nimber(v)
        Y = X
        for _ in range(a):
            Y = Y * Y
        t = L * X * Y
        # trace
        acc, s = t, t
        for _ in range(m - 1):
            s = s * s
            acc = acc + s
        tab.append(acc.value)
    return tab


def refine(tab, ell, m):
    """Q' = Q + <ell, .> -- the refinement torsor action."""
    N = 1 << m
    return [tab[v] ^ (bin(v & ell).count("1") & 1) for v in range(N)]


def candidates(m, kmax):
    """cand[i][k] = flip sets d with msb(d)=i and wt(d)<=k."""
    cand = [[[] for _ in range(kmax + 1)] for _ in range(m)]
    for i in range(m):
        base = 1 << i
        cand[i][1] = [base]
        cand[i][2] = [base] + [base | (1 << j) for j in range(i)]
        if kmax >= 3:
            c3 = list(cand[i][2])
            for j in range(i):
                for l in range(j):
                    c3.append(base | (1 << j) | (1 << l))
            cand[i][3] = c3
    return cand


def min_witness_radius(tab, m, cand, kmax=3):
    """max over v with Q(v)=1 of the minimal k giving a witness; None if some
    v is blocked even at kmax."""
    N = 1 << m
    worst = 0
    blocked = []
    for v in range(N):
        if tab[v] == 0:
            continue
        found = None
        for k in range(1, kmax + 1):
            ok = False
            for i in range(m):
                if not (v >> i) & 1:
                    continue
                for d in cand[i][k]:
                    if tab[v ^ d] == 0:
                        ok = True
                        break
                if ok:
                    break
            if ok:
                found = k
                break
        if found is None:
            blocked.append(v)
        else:
            worst = max(worst, found)
    return worst, blocked


# ---- 1. lambda = 43 sanity -------------------------------------------------
M, A = 8, 1
tab43 = make_qtab(M, A, 43)
diag = [tab43[1 << i] for i in range(M)]
z = tab43.count(0)
print(f"lambda=43, m=8: diagonal framing q_i = {diag}, |{{Q=0}}| = {z} "
      f"(bent iff in {{120,136}})")

# ---- 2 & 3. m=8 full sweep ---------------------------------------------------
cand = candidates(M, 3)
N = 1 << M
bent = []
for lam in range(1, N):
    tab = make_qtab(M, A, lam)
    if tab.count(0) in (120, 136):
        bent.append((lam, tab))
print(f"\nm=8, a=1: {len(bent)} bent components")

radius_hist = {}
worst_overall = 0
any_blocked = []
for lam, tab in bent:
    w, blk = min_witness_radius(tab, M, cand)
    radius_hist[w] = radius_hist.get(w, 0) + 1
    worst_overall = max(worst_overall, w)
    if blk:
        any_blocked.append((lam, None, blk[:4]))
print(f"Gold framing only: minimal witness radius histogram {radius_hist}, "
      f"blocked instances: {len(any_blocked)}")

# all 256 refinements per lambda
worst_ref = 0
ref_hist = {}
blocked_refs = []
for lam, tab in bent:
    for ell in range(N):
        t2 = refine(tab, ell, M)
        w, blk = min_witness_radius(t2, M, cand)
        ref_hist[w] = ref_hist.get(w, 0) + 1
        worst_ref = max(worst_ref, w)
        if blk:
            blocked_refs.append((lam, ell, blk[:4]))
print(f"ALL {len(bent)}x256 refinements: radius histogram {ref_hist}, "
      f"blocked: {len(blocked_refs)}")
if blocked_refs[:5]:
    print("  blocked examples:", blocked_refs[:5])

# ---- 4. m=16 spot checks -----------------------------------------------------
M2 = 16
cand16 = candidates(M2, 3)
half, off = 1 << (M2 - 1), 1 << (M2 // 2 - 1)
rng = random.Random(0xA9)
checked = 0
print(f"\nm=16, a=1 spot checks (bent components, Gold framing + 4 random "
      f"refinements each):")
lam = 1
while checked < 3:
    tab = make_qtab(M2, A, lam)
    z = tab.count(0)
    if z in (half - off, half + off):
        w, blk = min_witness_radius(tab, M2, cand16)
        line = f"  lambda={lam}: bent, Gold framing radius={w}, blocked={len(blk)}"
        radii = []
        for _ in range(4):
            t2 = refine(tab, rng.randrange(1 << M2), M2)
            w2, blk2 = min_witness_radius(t2, M2, cand16)
            radii.append((w2, len(blk2)))
        print(line + f"; random refinements (radius, blocked): {radii}")
        checked += 1
    lam += 1
print("done")
