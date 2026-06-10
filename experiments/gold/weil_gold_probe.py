"""Verification probe: the Weil/discriminant-form route to the Gold-quadric game question.

Four claims, each falsifiable:

P1 (dictionary). The nonsingular core of a Gold form (rank 2r, Arf eps) is, as a
    finite quadratic module, the discriminant form of the even lattice
    U(2)^{r-1} (+) D4   (eps=1)   /   U(2)^r   (eps=0),
    and  Arf = milgram_signature_mod8 / 4 = -weil_s_prefactor/4 mod 2.
    Checked with the SHIPPED DiscriminantForm + arf_nimber (cross-pillar oracle).

P2 (torsor selection). For S = sigma * 2^{-r} ((-1)^{B(g,d)}) with sigma = (-1)^eps,
    among diagonal matrices T_f = diag((-1)^{f(g)}) built from quadratic refinements
    f of B, the metaplectic relation (S T)^3 = S^2 holds IFF Arf(f) = eps.
    Non-quadratic f fail. So the B-only Weil apparatus pins exactly the Arf CLASS.

P3 (difference-set no-go input). For the bent Gold component on F_2^8 (r=4):
    every nonzero v is a difference of two Q-zeros, with the predicted count
    N(v) = 2^{2r-2} + (-1)^Arf 2^{r-1}; hence I*I = E in the extraspecial group.

P4 (center inertness). The E-lift (standard cocycle) of a center-blind rule
    (bent_route Rule B) has outcomes o(a,v) = o(v): the extension is outcome-inert.
"""

import sys, itertools, cmath
sys.path.insert(0, "/Users/a9lim/Work/ogdoad/experiments")
import ogdoad as pl

# ---------------------------------------------------------------- nim/gold helpers
def frob(x, a):
    for _ in range(a):
        x = x * x
    return x

def nim_trace(x, m):
    acc = pl.Nimber(x); t = pl.Nimber(x)
    for _ in range(m - 1):
        t = t * t
        acc = acc + t
    assert acc.value in (0, 1)
    return acc.value

def bent_gold(v, lam, a, m):
    x = pl.Nimber(v)
    return nim_trace((pl.Nimber(lam) * x * frob(x, a)).value, m)

# ================================================================ P1: the dictionary
print("=" * 76)
print("P1 — Gold core = discriminant form of U(2)^{r-1} (+) D4 ; Arf = Milgram/4")
print("=" * 76)

U2 = pl.IntegralForm([[0, 2], [2, 0]])
D4 = pl.IntegralForm.d(4)

def dsum(a_gram, b_gram):
    n, k = len(a_gram), len(b_gram)
    g = [[0] * (n + k) for _ in range(n + k)]
    for i in range(n):
        for j in range(n):
            g[i][j] = a_gram[i][j]
    for i in range(k):
        for j in range(k):
            g[n + i][n + j] = b_gram[i][j]
    return pl.IntegralForm(g)

def gram_of(L):
    # py binding accessor
    for name in ("gram", "gram_matrix"):
        if hasattr(L, name):
            attr = getattr(L, name)
            return attr() if callable(attr) else attr
    raise AttributeError("no gram accessor")

def module_data(L):
    """Return (reps, Qvals in {0,1}, B matrix over F2) if the module is F2-valued."""
    disc = pl.DiscriminantForm.from_lattice(L)
    reps = disc.reps() if callable(disc.reps) else disc.reps
    Q = []
    for y in reps:
        qv = disc.quadratic_value_mod2(y)          # Rational in [0,2)
        num, den = qv.numer, qv.denom
        num = num() if callable(num) else num
        den = den() if callable(den) else den
        assert den == 1, f"q-value not integer: {num}/{den} — not F2-valued"
        Q.append(num % 2)
    n = len(reps)
    B = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            bv = disc.bilinear_value_mod1(reps[i], reps[j])
            num, den = bv.numer, bv.denom
            num = num() if callable(num) else num
            den = den() if callable(den) else den
            assert den in (1, 2)
            B[i][j] = (num * (2 // den)) % 2       # value in (1/2)Z/Z -> bit
    return disc, reps, Q, B

def arf_of_module(reps, Q, B):
    """Coordinatise A_L = (Z/2)^k by the basis of reps with a single 1-entry...
    safer: pick any F2-basis of reps under XOR-of-rep-vectors mod the lattice is
    awkward; instead use the rep index group structure implicitly by brute force:
    find basis among reps s.t. all reps are F2-combos (group is (Z/2)^k)."""
    # group is elementary abelian: addition of reps = vector add then reduce mod 2
    # (valid because HNF pivots here are all 2)
    n = len(reps)
    k = n.bit_length() - 1
    idx = {tuple(r): i for i, r in enumerate(reps)}
    def add(i, j):
        s = tuple((reps[i][t] + reps[j][t]) % 2 for t in range(len(reps[0])))
        return idx[s]
    # basis = reps that are "unit-like": greedily extend an independent set
    basis, span = [], {idx[tuple([0] * len(reps[0]))]: 0}
    for i in range(n):
        if i in span:
            continue
        new = {}
        for s, mask in span.items():
            new[add(s, i)] = mask | (1 << len(basis))
        basis.append(i); span.update(new)
        if len(basis) == k:
            break
    coord = [None] * n
    for s, mask in span.items():
        coord[s] = mask
    # build Metric over Nimber in this coordinatisation
    qvec = [pl.Nimber(0)] * k
    bmat = {}
    for t, bi in enumerate(basis):
        qvec[t] = pl.Nimber(Q[bi])
    for t in range(k):
        for u in range(t + 1, k):
            bval = B[basis[t]][basis[u]]
            if bval:
                bmat[(t, u)] = pl.Nimber(1)
    return pl.arf_nimber(pl.NimberAlgebra(q=qvec, b=bmat)), coord, span

for label, L, want_arf, r in [
    ("U(2)              ", U2, 0, 1),
    ("D4                 ", D4, 1, 1),
    ("U(2) (+) D4        ", dsum([[0, 2], [2, 0]], gram_of(D4)), 1, 2),
    ("U(2)^2 (+) D4      ", dsum([[0, 2], [2, 0]], dsum([[0, 2], [2, 0]], gram_of(D4)).gram()
                                if hasattr(dsum([[0,2],[2,0]], gram_of(D4)), 'gram') else None), 1, 3)
        if False else ("U(2) (+) U(2)      ", dsum([[0, 2], [2, 0]], [[0, 2], [2, 0]]), 0, 2),
]:
    disc, reps, Q, B = module_data(L)
    res, coord, span = arf_of_module(reps, Q, B)
    mil = disc.milgram_signature_mod8()
    pre = disc.weil_s_prefactor_phase_mod8()
    ok_weil = disc.verify_weil_relations()
    zeros = Q.count(0)
    print(f"  {label} |A|={len(reps):>3}  q-zeros={zeros:>2}  "
          f"Arf(shipped classifier)={res.arf} rank={res.rank}  "
          f"Milgram={mil}  S-prefactor={pre}  weil_ok={ok_weil}  "
          f"[expect Arf={want_arf}, Milgram={4*want_arf}, rank={2*r}]")
    assert res.arf == want_arf and res.rank == 2 * r and res.radical_dim == 0
    assert mil == 4 * want_arf and pre == (-4 * want_arf) % 8 and ok_weil

# Gold side: m=8, a=2 has rank 4 Arf 1 (goldarf.tex Table 1) -> same module as U(2)(+)D4.
g_res = None
qvec = [pl.Nimber(0)] * 8
from common import gold as gold_unscaled, polar as polar_unscaled
qv = [pl.Nimber(gold_unscaled(1 << i, 2, 8)) for i in range(8)]
bm = {}
for i in range(8):
    for j in range(i + 1, 8):
        if polar_unscaled(1 << i, 1 << j, 2, 8):
            bm[(i, j)] = pl.Nimber(1)
g_res = pl.arf_nimber(pl.NimberAlgebra(q=qv, b=bm))
print(f"  Gold m=8,a=2 shipped classifier: rank={g_res.rank}, Arf={g_res.arf} "
      f"-> core is THE rank-4 Arf-1 module above (Dickson: (rank,Arf) complete).")
assert (g_res.rank, g_res.arf) == (4, 1)

# ================================================================ P2: torsor selection
print()
print("=" * 76)
print("P2 — (ST)^3 = S^2 selects EXACTLY the Arf class of refinements of B")
print("=" * 76)

def mat_mul(A, Bm):
    n = len(A); m = len(Bm[0]); inner = len(Bm)
    out = [[0j] * m for _ in range(n)]
    for i in range(n):
        Ai = A[i]
        for kk in range(inner):
            a = Ai[kk]
            if a == 0:
                continue
            Bk = Bm[kk]
            oi = out[i]
            for j in range(m):
                oi[j] += a * Bk[j]
    return out

def mat_close(A, Bm, tol=1e-9):
    return all(abs(A[i][j] - Bm[i][j]) <= tol
               for i in range(len(A)) for j in range(len(A)))

def run_torsor(r, pairs):
    n2 = 2 * r
    N = 1 << n2
    def Bf(u, v):
        acc = 0
        for (i, j) in pairs:
            acc ^= ((u >> i) & (v >> j) & 1) ^ ((u >> j) & (v >> i) & 1)
        return acc
    F = [[((-1) ** Bf(g, d)) / (2 ** r) for d in range(N)] for g in range(N)]
    def Q0(v):  # split refinement
        acc = 0
        for (i, j) in pairs:
            acc ^= (v >> i) & (v >> j) & 1
        return acc
    def arf_lin(l):  # Arf of Q0 + <l,.> in the symplectic basis
        acc = 0
        for (i, j) in pairs:
            qi = Q0(1 << i) ^ ((l >> i) & 1)
            qj = Q0(1 << j) ^ ((l >> j) & 1)
            acc ^= qi & qj
        return acc
    results = {}
    for eps in (0, 1):
        sigma = (-1) ** eps
        S = [[sigma * F[i][j] for j in range(N)] for i in range(N)]
        S2 = mat_mul(S, S)
        sel = []
        for l in range(N):  # all 2^{2r} refinements: f = Q0 + <l,.>
            T = [(-1) ** (Q0(g) ^ bin(g & l).count("1") % 2) for g in range(N)]
            ST = [[S[i][j] * T[j] for j in range(N)] for i in range(N)]
            ST3 = mat_mul(mat_mul(ST, ST), ST)
            holds = mat_close(ST3, S2)
            assert holds == (arf_lin(l) == eps), (r, eps, l)
            sel.append(holds)
        results[eps] = sum(sel)
    # non-quadratic f sanity (only meaningful if non-refinement functions exist)
    import random
    rng = random.Random(7)
    bad_fail = 0; tried = 0
    for _ in range(8):
        f = [0] + [rng.randint(0, 1) for _ in range(N - 1)]
        # skip if f happens to be a refinement of B
        if all(f[u ^ v] == f[u] ^ f[v] ^ Bf(u, v) for u in range(N) for v in range(N)):
            continue
        tried += 1
        T = [(-1) ** f[g] for g in range(N)]
        for eps in (0, 1):
            S = [[((-1) ** eps) * F[i][j] for j in range(N)] for i in range(N)]
            ST = [[S[i][j] * T[j] for j in range(N)] for i in range(N)]
            ST3 = mat_mul(mat_mul(ST, ST), ST)
            if not mat_close(ST3, mat_mul(S, S)):
                bad_fail += 1
    return results, bad_fail, tried

for r, pairs in [(1, [(0, 1)]), (2, [(0, 1), (2, 3)])]:
    res, bad_fail, tried = run_torsor(r, pairs)
    tot = 1 << (2 * r)
    exp0 = 2 ** (2 * r - 1) + 2 ** (r - 1)
    exp1 = 2 ** (2 * r - 1) - 2 ** (r - 1)
    print(f"  r={r}: sigma=+1 selects {res[0]}/{tot} refinements (expect {exp0} = #Arf-0); "
          f"sigma=-1 selects {res[1]}/{tot} (expect {exp1} = #Arf-1); "
          f"non-quadratic f fail relation in {bad_fail}/{2*tried} trials")
    assert res[0] == exp0 and res[1] == exp1 and bad_fail == 2 * tried

# ================================================================ P3: I*I = E (bent, r=4)
print()
print("=" * 76)
print("P3 — bent Gold component on F_2^8: every nonzero v is a difference of zeros")
print("=" * 76)

m, a = 8, 1
half, off = 1 << (m - 1), 1 << (m // 2 - 1)
lam = None
for l in range(1, 1 << m):
    z = sum(1 for v in range(1 << m) if bent_gold(v, l, a, m) == 0)
    if z in (half + off, half - off):
        lam, zcount = l, z
        break
arf_bent = 0 if zcount == half + off else 1
Qb = [bent_gold(v, lam, a, m) for v in range(1 << m)]
zeros = [v for v in range(1 << m) if Qb[v] == 0]
zs = set(zeros)
pred = (1 << (2 * 4 - 2)) + ((-1) ** arf_bent) * (1 << 3)
counts = set()
ok = True
for v in range(1, 1 << m):
    Nv = sum(1 for x in zeros if (x ^ v) in zs)
    counts.add(Nv)
    if Nv == 0:
        ok = False
print(f"  lambda={lam}, Arf={arf_bent}, |zeros|={zcount}; "
      f"N(v) counts over nonzero v: {sorted(counts)} (predicted {pred}); all > 0: {ok}")
assert counts == {pred} and ok
print("  => with both central lifts (a,b free over each zero), I*I = E: the")
print("     translation-invariant (Cayley) kernel spec on E is impossible for r>=2.")

# ================================================================ P4: center inertness
print()
print("=" * 76)
print("P4 — E-lift of a center-blind rule is outcome-isomorphic to its projection")
print("=" * 76)

def Bb(u, v):
    return Qb[u ^ v] ^ Qb[u] ^ Qb[v]

qd = [Qb[1 << i] for i in range(m)]
def beta(x, y):  # standard cocycle: sum_i q_i x_i y_i + sum_{i>j} B_ij x_i y_j
    acc = 0
    for i in range(m):
        if not ((x >> i) & 1):
            continue
        if (y >> i) & 1:
            acc ^= qd[i]
        for j in range(i):
            if (y >> j) & 1:
                acc ^= Bb(1 << i, 1 << j)
    return acc

# downstairs: bent_route Rule B (turn off set bit i iff B(v, e_i) = 1)
n = 1 << m
succ_v = [[v ^ (1 << i) for i in range(m)
           if (v >> i) & 1 and Bb(v, 1 << i) == 1] for v in range(n)]
# upstairs: positions g = (cbit, v) = cbit * n + v ; move = left mult by (0, e_i)
succ_e = []
for g in range(2 * n):
    cb, v = divmod(g, n)
    opts = []
    for i in range(m):
        if (v >> i) & 1 and Bb(v, 1 << i) == 1:
            nc = cb ^ beta(1 << i, v)
            opts.append(nc * n + (v ^ (1 << i)))
    succ_e.append(opts)

def solve(succ):
    lab = [None] * len(succ)
    def go(u):
        if lab[u] is not None:
            return lab[u]
        lab[u] = "W"
        res = "L"
        for w in succ[u]:
            if go(w) == "L":
                res = "W"
                break
        lab[u] = res
        return res
    for u in range(len(succ)):
        go(u)
    return lab

lab_v = solve(succ_v)
lab_e = solve(succ_e)
inert = all(lab_e[cb * n + v] == lab_v[v] for cb in (0, 1) for v in range(n))
print(f"  o(a,v) == o(v) for all 512 lifted positions: {inert}")
assert inert
print("  => the central bit is outcome-inert unless the rule READS it; any")
print("     center-blind E-route reduces to the unsolved V-problem verbatim.")
print()
print("ALL CHECKS PASSED")
