#!/usr/bin/env python3
"""Gold-diagonal source probe (ogdoad OPEN.md #1, [diagonal] angle).

Verifies, standalone (independent nim arithmetic, cross-checked against the
repo's published tables):

  L1  top-coin trace lemma: Tr_m(e_i) = [i == m-1]  on the bit frame
  L2  subfield vanishing:   q_i^{(m,a)} = 0 for i < m/2
  R   tower recursion:      q_{M+j}^{(2M,a)} = Tr_M((1 (+) u_a) e_j^{1+2^a}),
                            u_a = sum_{s<a} (2^{M-1})^{2^s}
  D   trace-dual lambda_a^{(m)} of the diagonal functional; lambda lives in the
      HALF FIELD F_{2^{m/2}}; drift lambda^{(2M)} != lambda^{(M)} (odd a forced)
  A   Arf classes of the canonical-source refinements
      {frame, frame+Tr, frame+ones, frame+ones+Tr} vs gold
"""
from functools import lru_cache

# ----------------------------------------------------------------- nim arithmetic
@lru_cache(maxsize=None)
def nm(x, y):
    if x > y:
        x, y = y, x
    if x < 2:
        return x * y
    t = 0
    while y >= (1 << (2 << t)):
        t += 1
    sh = 1 << t          # H = 2^(2^t)
    H = 1 << sh
    c, d = x >> sh, x & (H - 1)
    e, f = y >> sh, y & (H - 1)
    ce, cf, de, df = nm(c, e), nm(c, f), nm(d, e), nm(d, f)
    return ((ce ^ cf ^ de) << sh) ^ df ^ nm(ce, H >> 1)

assert nm(2, 2) == 3 and nm(2, 3) == 1 and nm(4, 4) == 6 and nm(2, 4) == 8

def frob(x, a):
    for _ in range(a):
        x = nm(x, x)
    return x

def tr(x, m):
    acc = t = x
    for _ in range(m - 1):
        t = nm(t, t)
        acc ^= t
    assert acc < 2, (x, m, acc)
    return acc

def gold(v, a, m):
    return tr(nm(v, frob(v, a)), m)

def pop(x):
    return bin(x).count("1") & 1

# --------------------------------------------------- L1: top-coin trace lemma
print("L1: Tr_m(e_i) on the bit frame (claim: indicator of the top coin)")
for m in (2, 4, 8, 16, 32):
    tvec = [tr(1 << i, m) for i in range(m)]
    ok = all(tvec[i] == (1 if i == m - 1 else 0) for i in range(m))
    print(f"  m={m:>2}: Tr(e_i) = {''.join(map(str,tvec))}  top-coin-indicator: {ok}")
    assert ok

# ------------------------------------------- L2 + diagonal tables + R recursion
print("\nL2/R: Gold diagonals q_i^(m,a) = Tr(e_i^(1+2^a)), low-half vanishing,")
print("      and the tower recursion q_{M+j}^(2M,a) = Tr_M((1+u_a) e_j^(1+2^a))")
Q = {}   # (m,a) -> q list
for m in (4, 8, 16, 32):
    k = m.bit_length() - 1
    for a in range(1, k + 1):
        q = [gold(1 << i, a, m) for i in range(m)]
        Q[(m, a)] = q
        lowzero = all(q[i] == 0 for i in range(m // 2))
        g = __import__("math").gcd(2 * a, m)
        print(f"  m={m:>2} a={a}: q={''.join(map(str,q))}  rank={m-g:>2}  "
              f"low-half-zero: {lowzero}")
        assert lowzero
print("\n  recursion check (each level 2M from level M):")
for M in (4, 8, 16):
    kM = M.bit_length() - 1
    for a in range(1, kM + 1):
        u = 1 << (M - 1)
        ua = 0
        x = u
        for s in range(a):
            ua ^= x
            x = nm(x, x)
        pred = [tr(nm(1 ^ ua, nm(1 << j, frob(1 << j, a))), M) for j in range(M)]
        got = Q[(2 * M, a)][M:]
        ok = pred == got
        print(f"  M={M:>2}->2M={2*M:>2} a={a}: u_a={ua:>6}  recursion holds: {ok}")
        assert ok

# -------------------------------------------------------- D: trace-dual lambda
def solve_f2(rows_in, rhs_in, n):
    rows, rhs = list(rows_in), list(rhs_in)
    lam, r, piv = 0, 0, []
    for col in range(n):
        p = next((kk for kk in range(r, n) if (rows[kk] >> col) & 1), None)
        if p is None:
            continue
        rows[r], rows[p] = rows[p], rows[r]
        rhs[r], rhs[p] = rhs[p], rhs[r]
        for kk in range(n):
            if kk != r and (rows[kk] >> col) & 1:
                rows[kk] ^= rows[r]
                rhs[kk] ^= rhs[r]
        piv.append(col)
        r += 1
    assert r == n
    for kk, col in enumerate(piv):
        if rhs[kk]:
            lam |= 1 << col
    return lam

print("\nD: trace-dual lambda_a^(m)  (unique nimber with Tr(lambda e_i) = q_i)")
LAM = {}
for m in (4, 8, 16, 32):
    k = m.bit_length() - 1
    T = [sum(tr(nm(1 << i, 1 << j), m) << j for j in range(m)) for i in range(m)]
    for a in range(1, k + 1):
        lam = solve_f2(T, Q[(m, a)], m)
        LAM[(m, a)] = lam
        # verify
        assert all(tr(nm(lam, 1 << i), m) == Q[(m, a)][i] for i in range(m))
        half = lam < (1 << (m // 2))
        # minimal subfield level
        d = 1
        while frob(lam, d) != lam:
            d *= 2
        print(f"  m={m:>2} a={a}: lambda = {lam:>10}  (binary {lam:b})  "
              f"in-half-field: {half}  min-subfield F_2^{d}")
        assert half
print("\n  tower drift (does one fixed nimber lambda work at every level?):")
for a in (1, 2, 3, 4):
    seq = [(m, LAM[(m, a)]) for m in (4, 8, 16, 32) if (m, a) in LAM]
    vals = [v for _, v in seq]
    coherent = len(set(vals)) == 1
    print(f"  a={a}: lambda^(m) = {seq}  coherent: {coherent}")

# ------------------------------------------------- A: Arf of canonical sources
print("\nA: Arf classes of canonical-source refinements (zero-count classifier)")

def classify(m, a, dvec):
    """Radical-adjusted class of Q_frame + sum d_i x_i (polar form = Gold B)."""
    q = Q[(m, a)]
    dmask = sum(b << i for i, b in enumerate(dvec))
    qmask = sum(b << i for i, b in enumerate(q))
    # gold values once, fast bilinear tables
    sq = [nm(1 << i, 1 << i) for i in range(m)]
    def frob_lin(v):
        out = 0
        for i in range(m):
            if (v >> i) & 1:
                out ^= sq[i]
        return out
    def frob_a(v):
        for _ in range(a):
            v = frob_lin(v)
        return v
    trmask = sum(tr(1 << i, m) << i for i in range(m))
    PROD = [[nm(1 << i, 1 << j) for j in range(m)] for i in range(m)]
    ROW = []
    for i in range(m):
        row = [0] * (1 << m)
        for w in range(1, 1 << m):
            lb = (w & -w).bit_length() - 1
            row[w] = row[w & (w - 1)] ^ PROD[i][lb]
        ROW.append(row)
    def mul_lin(v, w):
        out = 0
        for i in range(m):
            if (v >> i) & 1:
                out ^= ROW[i][w]
        return out
    N0 = 0
    goldset_check = 0
    for v in range(1 << m):
        gv = pop(mul_lin(v, frob_a(v)) & trmask)        # gold(v)
        fv = gv ^ pop(qmask & v)                        # frame(v)
        qv = fv ^ pop(dmask & v)                        # candidate
        N0 += (qv == 0)
        goldset_check += (qv == gv)
    exact = goldset_check == (1 << m)
    # radical of B and Q|radical
    Brow = []
    for i in range(m):
        bi = 0
        for j in range(m):
            if i != j:
                bij = gold((1 << i) ^ (1 << j), a, m) ^ q[i] ^ q[j]
                bi |= bij << j
        Brow.append(bi)
    # nullspace basis of Brow (vectors v with B(v, e_j)=0 for all j)
    # solve Brow^T v = 0 -> since B symmetric, Brow v = 0 rowwise
    basis = []
    rows = Brow[:]
    # gaussian elim to find nullspace of the m x m matrix over F2
    pivots = {}
    rr = []
    for row_i in range(m):
        row = rows[row_i]
        for c, prow in pivots.items():
            if (row >> c) & 1:
                row ^= prow
        if row:
            c = (row & -row).bit_length() - 1
            pivots[c] = row
            rr.append((c, row))
    free = [c for c in range(m) if c not in pivots]
    null = []
    for fcol in free:
        v = 1 << fcol
        # back-substitute
        changed = True
        while changed:
            changed = False
            for c, prow in pivots.items():
                if pop(prow & v) == 1 and not (v >> c) & 1:
                    v |= 1 << c
                    changed = True
                elif pop(prow & v) == 1 and (v >> c) & 1:
                    v &= ~(1 << c)
                    changed = True
        null.append(v)
    # verify nullspace + compute candidate Q on radical basis
    s = len(null)
    rad_aniso = False
    for v in null:
        assert all(pop(Brow[i] & v) == 0 for i in range(m))
        gv = gold(v, a, m)
        fv = gv ^ pop(qmask & v)
        if fv ^ pop(dmask & v):
            rad_aniso = True
    if rad_aniso:
        return ("balanced(aniso-radical)", N0, exact, s)
    twor = m - s
    bias = N0 - (1 << (m - 1))
    if bias == (1 << (s + twor // 2 - 1)):
        return ("Arf 0 (O+)", N0, exact, s)
    if bias == -(1 << (s + twor // 2 - 1)):
        return ("Arf 1 (O-)", N0, exact, s)
    return (f"?? N0={N0}", N0, exact, s)

published = {(4, 1): 4, (8, 1): 112, (8, 2): 96, (16, 1): 32512, (16, 4): 30720}
for m in (4, 8, 16):
    k = m.bit_length() - 1
    for a in range(1, k + 1):
        if m - __import__("math").gcd(2 * a, m) < 2:
            continue
        q = Q[(m, a)]
        trvec = [1 if i == m - 1 else 0 for i in range(m)]
        ones = [1] * m
        onestr = [1 ^ t for t in trvec]
        rows = []
        for label, d in [("gold (d=q)", q), ("frame (d=0)", [0] * m),
                         ("frame+Tr  (top coin)", trvec),
                         ("frame+ones (odious)", ones),
                         ("frame+ones+Tr", onestr)]:
            cls, N0, exact, s = classify(m, a, d)
            rows.append((label, cls, N0, exact))
        print(f"  m={m:>2} a={a} (rad dim {s}):")
        for label, cls, N0, exact in rows:
            star = "  <- = gold zero set" if exact else ""
            print(f"     {label:<22} {cls:<24} |Q=0|={N0:>6}{star}")
        if (m, a) in published:
            gN0 = next(N0 for label, cls, N0, e in rows if label.startswith("gold"))
            assert gN0 == published[(m, a)], (m, a, gN0)
            print(f"     [cross-check vs goldarf.tex Table 2: |Q=0|={gN0} OK]")
print("\nall assertions passed")
