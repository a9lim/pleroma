"""Independent skeptic verification of the T2 claim.

1. Blocking lemma + P-set(T2) = {Q=0} for RANDOM quadratic forms on F_2^10
   (random alternating B incl. degenerate, random diagonal q).
2. Same for Gold components Tr(lambda x^{1+2^a}) over the nimber field F_256,
   using my own nim multiplication (independent of the ogdoad crate and of the
   attacker's scripts). Includes unscaled degenerate (8,1),(8,2).
3. Spot-check the lambda=43 zero-diagonal claim.
"""
import random
from functools import lru_cache

# ------------- independent nim arithmetic (Conway recursion) -------------
@lru_cache(maxsize=None)
def nim_mul(a, b):
    if a < 2 or b < 2:
        return a * b
    F = 2
    while F * F <= max(a, b):
        F = F * F
    ah, al = divmod(a, F)
    bh, bl = divmod(b, F)
    hh = nim_mul(ah, bh)
    high = hh ^ nim_mul(ah, bl) ^ nim_mul(al, bh)
    low = nim_mul(al, bl) ^ nim_mul(hh, F >> 1)
    return high * F ^ low

assert nim_mul(2, 2) == 3 and nim_mul(2, 3) == 1 and nim_mul(4, 4) == 6
# field axiom spot-fuzz on F_16
rng = random.Random(7)
for _ in range(300):
    x, y, z = (rng.randrange(16) for _ in range(3))
    assert nim_mul(x, nim_mul(y, z)) == nim_mul(nim_mul(x, y), z)
    assert nim_mul(x, y ^ z) == nim_mul(x, y) ^ nim_mul(x, z)

def frob(x, a):
    for _ in range(a):
        x = nim_mul(x, x)
    return x

def nim_trace(x, m):
    acc, t = x, x
    for _ in range(m - 1):
        t = nim_mul(t, t)
        acc ^= t
    return acc

def gold_tab(lam, a, m):
    n = 1 << m
    return [nim_trace(nim_mul(lam, nim_mul(v, frob(v, a))), m) for v in range(n)]

# ------------- generic quadratic form machinery -------------
def q_from_tab(tab):
    """Return Q as the table; derive q_i and B_ij and check quadratic-ness."""
    return tab

def make_random_form(m, rng, rank_deficit=0):
    """Random alternating B (possibly degenerate) + random diagonal q; return Q table."""
    n = 1 << m
    B = [[0] * m for _ in range(m)]
    for i in range(m):
        for j in range(i + 1, m):
            B[i][j] = B[j][i] = rng.randrange(2)
    if rank_deficit:  # zero out some rows/cols to force degeneracy
        for k in range(rank_deficit):
            for j in range(m):
                B[k][j] = B[j][k] = 0
    q = [rng.randrange(2) for _ in range(m)]
    tab = []
    for v in range(n):
        s = sum(q[i] for i in range(m) if (v >> i) & 1) & 1
        for i in range(m):
            if not (v >> i) & 1:
                continue
            for j in range(i + 1, m):
                if (v >> j) & 1:
                    s ^= B[i][j]
        tab.append(s)
    return tab

def check_is_quadratic(tab, m):
    """Verify tab really is quadratic: third derivative vanishes (random checks)."""
    n = 1 << m
    rng = random.Random(1)
    for _ in range(200):
        v = rng.randrange(n); x = rng.randrange(n); y = rng.randrange(n); z = rng.randrange(n)
        d3 = 0
        for sx in (0, x):
            for sy in (0, y):
                for sz in (0, z):
                    d3 ^= tab[v ^ sx ^ sy ^ sz]
        if d3 != 0:
            return False
    return True

# ------------- T2 game: turn 1-2 coins, leading coin H->T, legal iff Q flips ----
def t2_pset(tab, m):
    n = 1 << m
    loss = [False] * n
    radius_fail = []
    for v in range(n):
        movs = []
        for j in range(m):
            if not (v >> j) & 1:
                continue
            d = 1 << j
            if tab[v ^ d] != tab[v]:
                movs.append(v ^ d)
            for i in range(j):           # i < j = msb(d); i may be outside supp(v)
                d2 = d | (1 << i)
                if tab[v ^ d2] != tab[v]:
                    movs.append(v ^ d2)
        loss[v] = not any(loss[w] for w in movs)   # all movs are < v
        if tab[v] == 1 and not movs:
            radius_fail.append(v)
    pset = frozenset(v for v in range(n) if loss[v])
    zset = frozenset(v for v in range(n) if tab[v] == 0)
    return pset == zset, radius_fail

# ------------- 1. random forms on F_2^10 -------------
m = 10
rng = random.Random(0xC1A0)
fails = 0
for trial in range(40):
    deficit = rng.choice([0, 0, 1, 3, m - 1, m])  # incl. fully degenerate B=0
    tab = make_random_form(m, rng, rank_deficit=deficit)
    ok, blocked = t2_pset(tab, m)
    if not ok or blocked:
        fails += 1
        print(f"  FAIL trial={trial} deficit={deficit} blocked={blocked[:5]}")
print(f"[1] random forms m={m}: 40 trials, failures = {fails}")

# ------------- 2. Gold components over F_256, my own nim arithmetic -------------
m, n = 8, 256
all_ok, n_bent, n_deg, blocked_any = True, 0, 0, 0
for a in (1, 2):
    for lam in range(1, n):
        tab = gold_tab(lam, a, m)
        z = tab.count(0)
        bent = z in (128 + 8, 128 - 8)
        n_bent += bent
        n_deg += (not bent)
        ok, blocked = t2_pset(tab, m)
        if not ok:
            all_ok = False
            print(f"  FAIL a={a} lam={lam}")
        if blocked:
            blocked_any += 1
assert check_is_quadratic(gold_tab(43, 1, 8), 8)
print(f"[2] Gold a=1,2 ALL lam in F_256: P-set == {{Q=0}} for all: {all_ok}; "
      f"bent={n_bent}, non-bent={n_deg}, blocked-Q=1-positions instances: {blocked_any}")

# unscaled degenerate Gold (8,1),(8,2),(4,1)
for (mm, aa) in ((8, 1), (8, 2), (4, 1)):
    tab = gold_tab(1, aa, mm)
    ok, blocked = t2_pset(tab, mm)
    print(f"    unscaled Gold ({mm},{aa}): zeros={tab.count(0)}/{1<<mm}, "
          f"P-set=={{Q=0}}: {ok}, blocked: {len(blocked)}")

# ------------- 3. lambda=43 diagonal claim -------------
tab43 = gold_tab(43, 1, 8)
diag = [tab43[1 << i] for i in range(8)]
z43 = tab43.count(0)
print(f"[3] lam=43 a=1 m=8: diagonal q_i = {diag} (claim: all zero), "
      f"|Q=0|={z43} (bent iff 120 or 136)")
# count bent lambdas and zero-diagonal bent lambdas (claims: 170 and 13)
bent_lams = []
zero_diag = []
for lam in range(1, 256):
    t = gold_tab(lam, 1, 8)
    if t.count(0) in (120, 136):
        bent_lams.append(lam)
        if all(t[1 << i] == 0 for i in range(8)):
            zero_diag.append(lam)
print(f"    bent lambdas (a=1): {len(bent_lams)} (claim 170); "
      f"zero-diagonal bent: {len(zero_diag)} (claim 13) -> {zero_diag}")
