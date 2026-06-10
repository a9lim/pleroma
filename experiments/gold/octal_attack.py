"""Attack probes for the octal/coin-turning angle on the Gold-quadric game rule.

Self-contained (no repo imports) so results are independent cross-checks:
  1. nim arithmetic validated against the repo's pinned small products and the
     arf_win_bias zero counts (F_16: 4, F_256 a=1: 112, a=2: 96).
  2. Theorem A check: normal-play coin-turning Grundy is XOR-additive in the
     configuration, hence P-sets are linear subspaces. Exhaustive over ALL
     2^15 coin-turning rules on 4 coins.
  3. Doubling lemma check: g(H)=2, options realize {0,1}, all proper followers
     g<=1  =>  H+H is a misere P-position.
  4. Misere coin-turning sweep: exhaustive over all 32768 rules on 4 coins:
     which misere P-sets are genuine quadrics? does any equal a (bent-)Gold
     zero/one set on F_16?
  5. Structured rules at m=8 (F_256) + Turning Corners misere on 3x3 and 4x4.
  6. The B-phase extension game (extraspecial-motivated candidate) on F_16.
"""
import sys
from functools import lru_cache

# ---------------------------------------------------------------- nim arithmetic
@lru_cache(maxsize=None)
def nim_mul(a, b):
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
    c2 = nim_mul(a1, b1)
    c1 = nim_mul(a1, b0) ^ nim_mul(a0, b1)
    c0 = nim_mul(a0, b0)
    return ((c1 ^ c2) * F) ^ c0 ^ nim_mul(c2, F >> 1)

def nim_sq(x):
    return nim_mul(x, x)

def frob(x, a):
    for _ in range(a):
        x = nim_sq(x)
    return x

def trace(x, m):
    acc, t = 0, x
    for _ in range(m):
        acc ^= t
        t = nim_sq(t)
    return acc

def gold(v, lam, a, m):
    return trace(nim_mul(lam, nim_mul(v, frob(v, a))), m)

def polar(u, v, lam, a, m):
    return gold(u ^ v, lam, a, m) ^ gold(u, lam, a, m) ^ gold(v, lam, a, m)

# validation
assert nim_mul(2, 2) == 3 and nim_mul(2, 3) == 1 and nim_mul(4, 4) == 6
assert nim_mul(2, 4) == 8 and nim_mul(16, 16) == 24
z16 = sum(1 for v in range(16) if gold(v, 1, 1, 4) == 0)
z256a1 = sum(1 for v in range(256) if gold(v, 1, 1, 8) == 0)
z256a2 = sum(1 for v in range(256) if gold(v, 1, 2, 8) == 0)
print(f"[nim ok] |Q=0| F16 a=1: {z16} (want 4); F256 a=1: {z256a1} (want 112); a=2: {z256a2} (want 96)")

# ------------------------------------------------------------- coin-turning solver
def solve_coin(n, comps, misere):
    """comps: list over r of list of companion masks (subsets of (1<<r)-1).
    Returns (pset, outcomes) with True = P-position."""
    N = 1 << n
    out = [False] * N
    for v in range(N):
        has_move = False
        to_p = False
        for r in range(n):
            if not (v >> r) & 1:
                continue
            for S in comps[r]:
                has_move = True
                if out[v ^ (1 << r) ^ S]:
                    to_p = True
                    break
            if to_p:
                break
        if not has_move:
            out[v] = not misere  # terminal: normal => P, misere => N
        else:
            out[v] = not to_p
    return [v for v in range(N) if out[v]], out

def grundy_coin(n, comps):
    N = 1 << n
    g = [0] * N
    for v in range(N):
        seen = set()
        for r in range(n):
            if not (v >> r) & 1:
                continue
            for S in comps[r]:
                seen.add(g[v ^ (1 << r) ^ S])
        x = 0
        while x in seen:
            x += 1
        g[v] = x
    return g

def is_linear_subspace(pts):
    s = set(pts)
    if 0 not in s:
        return False
    return all((x ^ y) in s for x in s for y in s)

def is_affine(pts):
    if not pts:
        return True
    s0 = pts[0]
    sh = {p ^ s0 for p in pts}
    return all((x ^ y) in sh for x in sh for y in sh)

# ------------------------------------------------------------------ quadric fit
def f2_rank(B, k):
    rows = [sum(B[i][j] << j for j in range(k)) for i in range(k)]
    rank = 0
    for col in range(k):
        piv = next((i for i in range(rank, k) if (rows[i] >> col) & 1), None)
        if piv is None:
            continue
        rows[rank], rows[piv] = rows[piv], rows[rank]
        for i in range(k):
            if i != rank and (rows[i] >> col) & 1:
                rows[i] ^= rows[rank]
        rank += 1
    return rank

def fit_quadratic(pset, k):
    N = 1 << k
    c = [1] * N
    for v in pset:
        c[v] = 0
    for i in range(k):
        bit = 1 << i
        for mask in range(N):
            if mask & bit:
                c[mask] ^= c[mask ^ bit]
    if any(c[m] and bin(m).count("1") > 2 for m in range(N)):
        return None
    B = [[0] * k for _ in range(k)]
    for i in range(k):
        for j in range(i + 1, k):
            if c[(1 << i) | (1 << j)]:
                B[i][j] = B[j][i] = 1
    return dict(const=c[0], qd=[c[1 << i] for i in range(k)], rank=f2_rank(B, k))

# =================================================================== experiment 2
print("\n== [2] Theorem A: normal-play coin-turning is XOR-additive (exhaustive m=4) ==")
def all_rules(n):
    """All companion-family choices per coin: families over subsets of lower coins."""
    import itertools
    per_coin = []
    for r in range(n):
        masks = list(range(1 << r))
        fams = []
        for bits in range(1 << len(masks)):
            fams.append([masks[i] for i in range(len(masks)) if (bits >> i) & 1])
        per_coin.append(fams)
    return per_coin

per_coin4 = all_rules(4)
total = bad_add = bad_lin = 0
for f0 in per_coin4[0]:
    for f1 in per_coin4[1]:
        for f2 in per_coin4[2]:
            for f3 in per_coin4[3]:
                comps = [f0, f1, f2, f3]
                g = grundy_coin(4, comps)
                # XOR-additivity of Grundy in the configuration
                ok = all(g[v] == (g[1] * 0 ^  # noqa
                          __import__("functools").reduce(lambda acc, r: acc ^ g[1 << r],
                                  [r for r in range(4) if (v >> r) & 1], 0))
                         for v in range(16))
                if not ok:
                    bad_add += 1
                pset = [v for v in range(16) if g[v] == 0]
                if not is_linear_subspace(pset):
                    bad_lin += 1
                total += 1
print(f"  rules checked: {total}; additivity failures: {bad_add}; nonlinear normal P-sets: {bad_lin}")

# =================================================================== experiment 3
print("\n== [3] Doubling lemma: g(H)=2, options {0,1}, followers g<=1 => H+H misere-P ==")
def mk(*opts):
    return frozenset(opts)
ZERO = mk()
STAR = mk(ZERO)
@lru_cache(maxsize=None)
def gsum(G, H):
    opts = {gsum(Gp, H) for Gp in G} | {gsum(G, Hp) for Hp in H}
    return frozenset(opts)
@lru_cache(maxsize=None)
def gr(G):
    seen = {gr(o) for o in G}
    x = 0
    while x in seen:
        x += 1
    return x
@lru_cache(maxsize=None)
def mis_n(G):
    if not G:
        return True  # terminal: mover wins (misere)
    return any(not mis_n(o) for o in G)

STAR2 = mk(ZERO, STAR)
B = mk(gsum(STAR, STAR))            # B = { *+* }, g(B) = 1, B != *
H1 = STAR2
H2 = mk(ZERO, STAR, B)              # custom H, options g in {0,1}
for name, H in [("*2", H1), ("{0,*,{*+*}}", H2)]:
    fol_ok = all(gr(o) <= 1 for o in H)  # options; deeper followers checked below
    def followers(G, acc):
        for o in G:
            if o not in acc:
                acc.add(o)
                followers(o, acc)
        return acc
    deep_ok = all(gr(x) <= 1 for x in followers(H, set()))
    hh_p = not mis_n(gsum(H, H))
    print(f"  H={name}: g(H)={gr(H)}, all proper followers g<=1: {deep_ok}, H+H misere-P: {hh_p}")

# =================================================================== experiment 4
print("\n== [4] Misere coin-turning, exhaustive over all 32768 rules on 4 coins ==")
gold_targets = {}
for lam in range(1, 16):
    zs = frozenset(v for v in range(16) if gold(v, lam, 1, 4) == 0)
    fitz = fit_quadratic(sorted(zs), 4)
    gold_targets[lam] = (zs, fitz["rank"] if fitz else None)
ranks = sorted({r for (_, r) in gold_targets.values()})
print(f"  Gold/bent-Gold targets on F16 (a=1): polar ranks present: {ranks}")
target_sets = set()
for lam, (zs, r) in gold_targets.items():
    if r and r >= 2:
        target_sets.add(zs)                       # {Q=0}
        target_sets.add(frozenset(range(16)) - zs)  # {Q=1}
print(f"  distinct genuine-quadric target sets (zero sets and one sets): {len(target_sets)}")

from collections import Counter
quad_psets = Counter()
quad_examples = {}
gold_hits = []
n_rules = n_quad = n_affine = 0
for f0 in per_coin4[0]:
    for f1 in per_coin4[1]:
        for f2 in per_coin4[2]:
            for f3 in per_coin4[3]:
                comps = [f0, f1, f2, f3]
                pset, _ = solve_coin(4, comps, misere=True)
                n_rules += 1
                fp = frozenset(pset)
                if is_affine(pset):
                    n_affine += 1
                fit = fit_quadratic(pset, 4)
                if fit and fit["rank"] >= 2:
                    n_quad += 1
                    quad_psets[fp] += 1
                    quad_examples.setdefault(fp, comps)
                    if fp in target_sets:
                        gold_hits.append((comps, sorted(fp)))
print(f"  rules: {n_rules}; affine misere P-sets: {n_affine}; genuine-quadric P-sets: {n_quad}")
print(f"  distinct genuine-quadric P-sets: {len(quad_psets)}")
for fp, cnt in quad_psets.most_common(10):
    fit = fit_quadratic(sorted(fp), 4)
    print(f"    P-set {sorted(fp)} (|{len(fp)}|, rank {fit['rank']}, const {fit['const']}): {cnt} rules"
          + ("  <-- GOLD TARGET MATCH" if fp in target_sets else ""))
if gold_hits:
    c, s = gold_hits[0]
    print(f"  *** GOLD HIT: companions={c} P-set={s} (total {len(gold_hits)} rules)")
else:
    print("  no misere coin-turning rule on 4 coins has a (bent-)Gold quadric P-set")

# =================================================================== experiment 5
print("\n== [5] Structured rules at m=8 (F_256) and Turning Corners misere ==")
def comps_ruler(n):      # turn r + any subset below: g(r)=2^r (field coordinates)
    return [list(range(1 << r)) for r in range(n)]
def comps_singleton(n):  # exactly one lower coin: g(r)=r
    return [[1 << i for i in range(r)] for r in range(n)]
def comps_turtles(n):    # at most one lower coin: g(r)=r+1
    return [[0] + [1 << i for i in range(r)] for r in range(n)]
def comps_mock(n):       # at most two lower coins (Mock Turtles): g odious
    out = []
    for r in range(n):
        fam = [0] + [1 << i for i in range(r)]
        fam += [(1 << i) | (1 << j) for i in range(r) for j in range(i + 1, r)]
        out.append(fam)
    return out
def comps_isotropic(n, lam, a):  # alphabet = turn sets with Q(chi_T)=0 (Q-referencing!)
    out = []
    for r in range(n):
        fam = [S for S in range(1 << r) if gold(S | (1 << r), lam, a, n) == 0]
        out.append(fam)
    return out

for name, mk_c in [("ruler", comps_ruler), ("singleton", comps_singleton),
                   ("turtles", comps_turtles), ("mock", comps_mock)]:
    for m in (4, 8):
        comps = mk_c(m)
        pset, _ = solve_coin(m, comps, misere=True)
        fit = fit_quadratic(pset, m)
        desc = ("affine" if is_affine(pset) else
                (f"quadric rank {fit['rank']}" if fit else "deg>2"))
        print(f"  misere {name:9s} m={m}: |P|={len(pset):4d}  {desc}")

# isotropic alphabet (Tier-3-flavored comparison)
for m, lam in [(4, 1), (8, 1)]:
    comps = comps_isotropic(m, lam, 1)
    pset, _ = solve_coin(m, comps, misere=True)
    fit = fit_quadratic(pset, m)
    target = frozenset(v for v in range(1 << m) if gold(v, lam, 1, m) == 0)
    desc = f"quadric rank {fit['rank']}" if fit else ("affine" if is_affine(pset) else "deg>2")
    print(f"  misere isotropic-alphabet m={m}: |P|={len(pset)} {desc}; equals "
          f"{{Q=0}}? {frozenset(pset)==target}; equals {{Q=1}}? {frozenset(pset)==frozenset(range(1<<m))-target}")
    # normal play too
    psetN, _ = solve_coin(m, comps, misere=False)
    print(f"  normal isotropic-alphabet m={m}: |P|={len(psetN)} linear? {is_linear_subspace(psetN)}")

# Turning Corners misere on k x k
def comps_corners(k):
    n = k * k
    idx = lambda x, y: x * k + y
    comps = [[] for _ in range(n)]
    for x in range(k):
        for y in range(k):
            fam = []
            for a_ in range(x):
                for b_ in range(y):
                    fam.append((1 << idx(a_, b_)) | (1 << idx(a_, y)) | (1 << idx(x, b_)))
            comps[idx(x, y)] = fam
    return comps
for k in (2, 3, 4):
    n = k * k
    pset, _ = solve_coin(n, comps_corners(k), misere=True)
    fit = fit_quadratic(pset, n)
    desc = ("affine" if is_affine(pset) else (f"quadric rank {fit['rank']}" if fit else "deg>2"))
    print(f"  misere TurningCorners {k}x{k}: |P|={len(pset)}  {desc}")

# random sparse rules at m=8
import random
random.seed(0)
n_q8 = 0
for trial in range(400):
    comps = []
    for r in range(8):
        kfam = random.randint(1, 4)
        comps.append([random.randrange(1 << r) for _ in range(kfam)])
    pset, _ = solve_coin(8, comps, misere=True)
    fit = fit_quadratic(pset, 8)
    if fit and fit["rank"] >= 2:
        n_q8 += 1
        target_hit = any(frozenset(pset) == frozenset(v for v in range(256) if gold(v, lam, 1, 8) == b)
                         for lam in range(1, 256) for b in (0, 1))
        print(f"  random m=8 rule #{trial}: GENUINE QUADRIC |P|={len(pset)} rank={fit['rank']} gold-match={target_hit}")
print(f"  random m=8 rules with genuine-quadric misere P-sets: {n_q8}/400")

# =================================================================== experiment 6
print("\n== [6] B-phase extension game (extraspecial candidate) on F_16, a=1 ==")
def solve_phase(m, lam, a, alphabet, win):
    """Positions (v, eps). Moves: r in supp(v), T = S|{r} in alphabet(r):
       (v,eps) -> (v^T, eps ^ B(T, v)).  win in {'phase','normal','misere'}.
       Returns P-set of (v,0) slice and full outcome map."""
    N = 1 << m
    out = {}
    for v in range(N):
        for eps in (0, 1):
            moves = []
            for r in range(m):
                if not (v >> r) & 1:
                    continue
                for S in alphabet[r]:
                    T = S | (1 << r)
                    moves.append((v ^ T, eps ^ polar(T, v, lam, a, m)))
            if not moves:
                if win == "phase":
                    out[(v, eps)] = (eps == 0)   # P iff phase 0 at the end
                elif win == "normal":
                    out[(v, eps)] = True
                else:
                    out[(v, eps)] = False
            else:
                out[(v, eps)] = not any(out[w] for w in moves)
    return [v for v in range(N) if out[(v, 0)]], out

m, lam, a = 4, 1, 1
targ0 = frozenset(v for v in range(16) if gold(v, lam, a, m) == 0)
alpha_full = [list(range(1 << r)) for r in range(m)]
alpha_iso = comps_isotropic(m, lam, a)
alpha_single = [[0] for _ in range(m)]
for an, alpha in [("full", alpha_full), ("isotropic", alpha_iso), ("singleton", alpha_single)]:
    for wn in ("phase", "normal", "misere"):
        pset, _ = solve_phase(m, lam, a, alpha, wn)
        fit = fit_quadratic(pset, m)
        desc = (f"quadric rank {fit['rank']}" if fit and fit["rank"] >= 2
                else ("affine" if is_affine(pset) else ("deg<=2 rank<2" if fit else "deg>2")))
        eq0 = frozenset(pset) == targ0
        eq1 = frozenset(pset) == frozenset(range(16)) - targ0
        print(f"  alphabet={an:9s} win={wn:6s}: |P(v,0)|={len(pset):2d} {desc}"
              f"{'  == {Q=0} !!' if eq0 else ''}{'  == {Q=1} !!' if eq1 else ''}")
print("done")
