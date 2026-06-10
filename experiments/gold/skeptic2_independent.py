"""Independent skeptic verification of ECHO-FIFO+dummy.

Everything from scratch: nim arithmetic, Gold forms, the game solver.
Checks:
 1. nim mul pinned products.
 2. m=4: outcome == Q for all 15 lambda x both orientations x all 16 x.
 3. Tautology tell: the same rule realizes ARBITRARY (q,B) - random
    non-Gold forms - i.e. the construction has no Gold-specific content;
    value(x) is identically the polynomial sum q_i + sum B_ij on support.
 4. Decision-nondegeneracy spot check (are there mistake states at all).
"""
import sys, random
sys.setrecursionlimit(1000000)

# ---------------- nim arithmetic (independent implementation) ----------------
_nm = {}
def nmul(a, b):
    if a < b:
        a, b = b, a
    if b == 0:
        return 0
    if b == 1:
        return a
    k = _nm.get((a, b))
    if k is not None:
        return k
    F = 2
    while F * F <= a:
        F *= F
    ah, al = a >> F.bit_length() - 1, a & (F - 1)
    bh, bl = b >> F.bit_length() - 1, b & (F - 1)
    t1 = nmul(ah, bh)
    t2 = nmul(ah, bl) ^ nmul(al, bh)
    t3 = nmul(al, bl)
    r = ((t1 ^ t2) * F) ^ t3 ^ nmul(t1, F >> 1)
    _nm[(a, b)] = r
    return r

assert nmul(2, 2) == 3 and nmul(2, 4) == 8 and nmul(16, 16) == 24
assert nmul(5, 9) == nmul(9, 5)

def npow2k(x, k):  # x^(2^k) by repeated nim squaring
    for _ in range(k):
        x = nmul(x, x)
    return x

def trace(x, m):
    t, y = 0, x
    for _ in range(m):
        t ^= y
        y = nmul(y, y)
    return t

def make_form(m, a, lam):
    n = 1 << m
    Q = [trace(nmul(lam, nmul(x, npow2k(x, a))), m) for x in range(n)]
    assert all(v in (0, 1) for v in Q), "trace not landing in F2"
    qd = [Q[1 << i] for i in range(m)]
    B = [[0] * m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            if i != j:
                B[i][j] = Q[(1 << i) ^ (1 << j)] ^ Q[1 << i] ^ Q[1 << j]
    # sanity: Q is quadratic-consistent on a few random triples
    rng = random.Random(1)
    for _ in range(50):
        u, v = rng.randrange(n), rng.randrange(n)
        bb = 0
        for i in range(m):
            for j in range(m):
                if j > i and (u >> i) & 1 and (v >> j) & 1:
                    bb ^= B[i][j]
                if j > i and (v >> i) & 1 and (u >> j) & 1:
                    bb ^= B[i][j]
        assert Q[u ^ v] == Q[u] ^ Q[v] ^ bb, "polar identity fails"
    return Q, qd, B

# ---------------- the game, written directly from the spec -------------------
def game_value(support, qd, B, t, dummy_idx):
    """ECHO-FIFO+dummy. Coins = support + dummy (q=0, B-isolated, highest idx).
    Open any untouched coin; close only openseq[0] (FIFO). Ko: may not touch
    the coin touched on the immediately preceding touch; stuck => pass (clears
    ko). Charge on touching c with open set o (lower-triangular cocycle):
    q_c if c in o, plus XOR_{k in o, k>c} B[k][c]. P1 wants final sigma == t."""
    coins = list(support) + [dummy_idx]
    q = {c: (0 if c == dummy_idx else qd[c]) for c in coins}
    def Bv(i, j):
        if i == dummy_idx or j == dummy_idx:
            return 0
        return B[i][j]
    memo = {}
    def val(u, openseq, last, mover, sigma):
        if not u and not openseq:
            return sigma
        key = (u, openseq, last, mover, sigma)
        r = memo.get(key)
        if r is not None:
            return r
        oset = set(openseq)
        legal = []
        for c in coins:
            if c == last:
                continue
            if c in u:  # open it
                ch = 0
                for k in oset:
                    if k > c:
                        ch ^= Bv(k, c)
                legal.append((c, ch, u - {c}, openseq + (c,)))
        if openseq:
            c = openseq[0]
            if c != last:
                ch = q[c]
                for k in oset:
                    if k > c:
                        ch ^= Bv(k, c)
                legal.append((c, ch, u, openseq[1:]))
        if not legal:
            res = val(u, openseq, -1, 1 - mover, sigma)
        else:
            want = t if mover == 0 else 1 - t
            res = 1 - want
            for (c, ch, u2, s2) in legal:
                if val(u2, s2, c, 1 - mover, sigma ^ ch) == want:
                    res = want
                    break
        memo[key] = res
        return res
    return val(frozenset(coins), (), -1, 0, 0)

# ---------------- check 2: m=4 Gold exactness --------------------------------
print("m=4 Gold exactness (independent implementation):")
allok = True
for lam in range(1, 16):
    Q, qd, B = make_form(4, 1, lam)
    for t in (0, 1):
        miss = [x for x in range(16)
                if game_value(tuple(i for i in range(4) if (x >> i) & 1),
                              qd, B, t, 4) != Q[x]]
        if miss:
            allok = False
            print(f"  lam={lam} t={t}: MISSES {miss}")
print("  ALL 15 lambda x 2 t x 16 x: EXACT" if allok else "  FAILURES FOUND")

# also m=4, a=... only a=1 valid for m=4 APN; try a=3 (gcd(3,4)=1)
ok3 = True
for lam in (1, 5, 9):
    Q, qd, B = make_form(4, 3, lam)
    for t in (0, 1):
        for x in range(16):
            S = tuple(i for i in range(4) if (x >> i) & 1)
            if game_value(S, qd, B, t, 4) != Q[x]:
                ok3 = False
print(f"  spot a=3 forms: {'EXACT' if ok3 else 'FAIL'}")

# ---------------- check 3: tautology tell - arbitrary (q,B) ------------------
print("\ntautology tell: random NON-Gold (q,B), m=5, all x, both t:")
rng = random.Random(99)
fails = expected_eq = 0
for trial in range(3):
    m = 5
    qd = [rng.randint(0, 1) for _ in range(m)]
    B = [[0] * m for _ in range(m)]
    for i in range(m):
        for j in range(i + 1, m):
            B[i][j] = B[j][i] = rng.randint(0, 1)
    for x in range(1 << m):
        S = tuple(i for i in range(m) if (x >> i) & 1)
        # the closed-form polynomial evaluation of Q(x):
        poly = 0
        for i in S:
            poly ^= qd[i]
        for ii in range(len(S)):
            for jj in range(ii + 1, len(S)):
                poly ^= B[S[ii]][S[jj]]
        for t in (0, 1):
            v = game_value(S, qd, B, t, m)
            if v == poly:
                expected_eq += 1
            else:
                fails += 1
print(f"  game value == polynomial Q(x) on support: {expected_eq} eq, {fails} neq")
print("  -> rule realizes EVERY (q,B); value is the closed-form polynomial"
      if fails == 0 else "  -> rule does NOT reduce to the polynomial")
