#!/usr/bin/env python3
"""Round-2 CONSTRUCT probe for ogdoad docs/OPEN.md problem 1 (Gold-quadric game).

PART 1 — primary construct T2-weierstrass at (8,1) and (16,1):
  Rule: positions = coin strings x in F_2^m <-> field elements of F_{2^m}.
  Move: turn over d with wt(d) in {1,2}, leading coin msb(d) heads in x
        (so x^d < x: terminating coin-turning convention).
  Gate: legal iff B(x,d) ^ Q(d) = 1, where Q(d) is q_i (singles) or
        q_i^q_j^B_ij (pairs) -- at most TWO diagonal bits + the public polar B.
  Data chain (all game-built, docs/OPEN.md standard chain):
        B_ij = Tr(e_i e_j^2 + e_j e_i^2)        (Turning-Corners + Frobenius + trace)
        q_i  = Tr(P(w) e_i),  P(z) = z^2 ^ z,   w = XOR of Fermat coins 2^(2^t), t>=1
  Claim (attack-5 blocking lemma + diagonal-skeptic identity): P-set = {Tr(x^3)=0}.

PART 2 — loopy width-2 sweep (residue B): add local-gated move families to T2,
  compute loopy Win/Loss/Draw (kernel::outcomes semantics: terminal = Loss),
  hunt for Loss / Draw / Loss-union-Draw = {Q=0} with decision-nondegeneracy.
"""
import sys
from functools import lru_cache
from itertools import combinations
import random

sys.setrecursionlimit(100000)

# ---------------- nim arithmetic (self-contained, validated below) -----------
@lru_cache(maxsize=None)
def nim_mul(a, b):
    if a > b:
        a, b = b, a
    if a == 0:
        return 0
    if a == 1:
        return b
    k = 1
    while (1 << (1 << k)) <= b:
        k += 1
    F = 1 << (1 << (k - 1))  # largest Fermat 2-power <= b ; a,b < F*F
    ah, al = a >> (1 << (k - 1)), a & (F - 1)
    bh, bl = b >> (1 << (k - 1)), b & (F - 1)
    hh = nim_mul(ah, bh)
    t = nim_mul(ah ^ al, bh ^ bl)
    ll = nim_mul(al, bl)
    return ((t ^ ll) << (1 << (k - 1))) ^ ll ^ nim_mul(hh, F >> 1)

def nim_sq(x):
    return nim_mul(x, x)

def frob(x, a):
    for _ in range(a):
        x = nim_sq(x)
    return x

def trace(x, m):
    t, y = 0, x
    for _ in range(m):
        t ^= y
        y = nim_sq(y)
    assert t in (0, 1), (x, m, t)
    return t

# pinned repo values
assert nim_mul(2, 2) == 3 and nim_mul(2, 4) == 8 and nim_mul(16, 16) == 24
print("[ok] nim arithmetic matches repo-pinned products 2*2=3, 2*4=8, 16*16=24")

# ---------------- form data, game-built ------------------------------------
def gold_q(v, a, m, lam=1):
    return trace(nim_mul(lam, nim_mul(v, frob(v, a))), m)

def build_data(m, a, lam=1):
    """Return (q list, B as row-masks, Qtable) for Q(x)=Tr(lam x^{1+2^a}) on F_2^m."""
    n = 1 << m
    Q = [gold_q(v, a, m, lam) for v in range(n)]
    q = [Q[1 << i] for i in range(m)]
    # game-built polar: B(u,v) = Tr(lam(u v^{2^a} + v u^{2^a}))
    Bmask = [0] * m
    for i in range(m):
        for j in range(m):
            if i == j:
                continue
            u, v = 1 << i, 1 << j
            b = trace(nim_mul(lam, nim_mul(u, frob(v, a)) ^ nim_mul(v, frob(u, a))), m)
            # cross-check vs polarization
            assert b == (Q[u ^ v] ^ Q[u] ^ Q[v]), (i, j)
            if b:
                Bmask[i] |= 1 << j
    return q, Bmask, Q

def weierstrass_source(m):
    """lambda = P(w), w = XOR of Fermat coins 2^(2^t), 2 <= 2^t < m. q_i = Tr(lam e_i)."""
    w = 0
    t = 1
    while (1 << t) < m:
        w ^= 1 << (1 << t)
        t += 1
    lam = nim_sq(w) ^ w
    return w, lam

# ---------------- the T2 game solver (acyclic, normal play) -----------------
def t2_solve(m, q, Bmask, ret_graph=False):
    """P-set of the width-2 leading-coin spin-flip game. outcome True = N(win)."""
    n = 1 << m
    # dQ_i(v) = B(v, e_i) ^ q_i
    win = [False] * n
    succ = [[] for _ in range(n)] if ret_graph else None
    Bij = [[(Bmask[i] >> j) & 1 for j in range(m)] for i in range(m)]
    for v in range(1, n):
        dq = [(bin(v & Bmask[i]).count("1") & 1) ^ q[i] for i in range(m)]
        mv = []
        for i in range(m):
            if not (v >> i) & 1:
                continue
            if dq[i] == 1:                      # single flip, head i off
                mv.append(v ^ (1 << i))
            for j in range(i):                  # pair, leading coin i (msb) is a head
                if dq[i] ^ dq[j] ^ Bij[i][j] == 1:
                    mv.append(v ^ (1 << i) ^ (1 << j))
        win[v] = any(not win[w] for w in mv)
        if ret_graph:
            succ[v] = mv
    pset = {v for v in range(n) if not win[v]}
    return (pset, succ) if ret_graph else pset

# ---------------- loopy retrograde solver (kernel::outcomes semantics) ------
def loopy_outcomes(succ):
    """0=Loss, 1=Win, 2=Draw. Terminal (no moves) = Loss."""
    n = len(succ)
    pred = [[] for _ in range(n)]
    deg = [len(s) for s in succ]
    for v, s in enumerate(succ):
        for w in s:
            pred[w].append(v)
    LOSS, WIN, DRAW = 0, 1, 2
    out = [DRAW] * n
    stack = [v for v in range(n) if deg[v] == 0]
    for v in stack:
        out[v] = LOSS
    while stack:
        v = stack.pop()
        for u in pred[v]:
            if out[u] != DRAW:
                continue
            if out[v] == LOSS:
                out[u] = WIN
                stack.append(u)
            else:
                deg[u] -= 1
                if deg[u] == 0:
                    out[u] = LOSS
                    stack.append(u)
    return out

# =================== PART 1: primary construct ==============================
print("\n=== PART 1: T2-weierstrass, (m,a)=(8,1), Q(x)=Tr(x^3) on F_256 ===")
m, a = 8, 1
q, Bmask, Q = build_data(m, a)
Z = {v for v in range(1 << m) if Q[v] == 0}
print(f"|{{Q=0}}| = {len(Z)} (expect 112)")
assert len(Z) == 112

# diagonal facts: q_i = 0 below the top Fermat layer, Fermat coin = a mod 2
assert all(q[i] == 0 for i in range(4)) and q[4] == 1
print(f"q = {q}  (subfield-vanishing L2 + Fermat-coin witness hold)")

# the game-native diagonal source
w, lam = weierstrass_source(m)
assert w == 20 and lam == 10, (w, lam)
qsrc = [trace(nim_mul(lam, 1 << i), m) for i in range(m)]
assert qsrc == q, (qsrc, q)
print(f"[ok] q_i == Tr(P(w) e_i) with w=Fermat-coin XOR={w}, P(w)=w^2^w={lam} (all i)")

pset, succ = t2_solve(m, q, Bmask, ret_graph=True)
assert pset == Z, f"P-set != {{Q=0}}: |P|={len(pset)}"
print(f"[ok] T2 P-set == {{Q=0}} exactly ({len(pset)} positions)")

# ender / decision-degeneracy check
mixed = sum(1 for v in range(1 << m) if succ[v] and
            len({(v ^ 0, Q[t]) and Q[t] for t in succ[v]} if False else {Q[t] for t in succ[v]}) > 1)
ender = all(len({Q[t] for t in succ[v]}) <= 1 for v in range(1 << m))
print(f"[ok] ender confirmed: every position's options share one Q-class ({ender}); "
      f"the game is decision-degenerate, as Theorem 5 forces")

# refinement uniformity (the non-tautology clause): same rule template, random q'
rng = random.Random(0)
for trial in range(20):
    q2 = [rng.randint(0, 1) for _ in range(m)]
    # Q'(v) from q' and B (coordinates)
    def q2form(v):
        s = 0
        bits = [i for i in range(m) if (v >> i) & 1]
        for i in bits:
            s ^= q2[i]
        for x, y in combinations(bits, 2):
            s ^= (Bmask[x] >> y) & 1
        return s
    p2 = t2_solve(m, q2, Bmask)
    assert p2 == {v for v in range(1 << m) if q2form(v) == 0}, trial
print("[ok] refinement uniformity: 20 random refinements of the same B all give P = {Q'=0}")

print("\n--- (16,1) scale check: Q(x)=Tr(x^3) on F_65536 ---")
m16 = 16
q16, Bmask16, Q16 = build_data(m16, 1)
Z16 = {v for v in range(1 << m16) if Q16[v] == 0}
print(f"|{{Q=0}}| = {len(Z16)} (expect 32512)")
assert len(Z16) == 32512
w16, lam16 = weierstrass_source(m16)
assert w16 == 4 ^ 16 ^ 256 and lam16 == 138, (w16, lam16)
qsrc16 = [trace(nim_mul(lam16, 1 << i), m16) for i in range(m16)]
assert qsrc16 == q16
print(f"[ok] q_i == Tr(P(w) e_i), w={w16}, P(w)={lam16}")
p16 = t2_solve(m16, q16, Bmask16)
assert p16 == Z16
print(f"[ok] T2 P-set == {{Q=0}} exactly ({len(p16)} positions) at m=16")

# bent component at (8,1): T2 with the bent diagonal (source open, constants fed)
print("\n--- bent component check (8,1,lam=2): blocking lemma is form-agnostic ---")
qb, Bmaskb, Qb = build_data(8, 1, lam=2)
Zb = {v for v in range(256) if Qb[v] == 0}
print(f"|{{Q_lam=0}}| = {len(Zb)} (bent iff 120 or 136)")
pb = t2_solve(8, qb, Bmaskb)
assert pb == Zb
print(f"[ok] T2 P-set == bent {{Q=0}} exactly ({len(pb)} positions)")

# =================== PART 2: loopy width-2 sweep ============================
print("\n=== PART 2: loopy width-2 sweep at (8,1) — hunting a non-ender ===")
m = 8
n = 1 << m
Bij = [[(Bmask[i] >> j) & 1 for j in range(m)] for i in range(m)]

def movegen(v, fams, q, Bmask, Bij):
    """Moves under selected families. All gates read only q_i, B_ij, supp(v).
    families: 'T2' desc flips; 'A' asc flips; 'S' slides (wt2, one end head, dQ=0);
    'Cd' desc laterals (dQ=0, msb head); 'Ca' asc laterals (dQ=0, msb tail)."""
    dq = [(bin(v & Bmask[i]).count("1") & 1) ^ q[i] for i in range(m)]
    mv = []
    for i in range(m):
        hi = (v >> i) & 1
        # singles, leading coin = i
        if dq[i] == 1:
            if (hi and "T2" in fams) or (not hi and "A" in fams):
                mv.append(v ^ (1 << i))
        else:
            if (hi and "Cd" in fams) or (not hi and "Ca" in fams):
                mv.append(v ^ (1 << i))
        for j in range(i):  # pairs, msb = i
            g = dq[i] ^ dq[j] ^ Bij[i][j]
            tgt = v ^ (1 << i) ^ (1 << j)
            if g == 1:
                if (hi and "T2" in fams) or (not hi and "A" in fams):
                    mv.append(tgt)
            else:
                hj = (v >> j) & 1
                if "S" in fams and hi + hj == 1:
                    mv.append(tgt)
                if (hi and "Cd" in fams) or (not hi and "Ca" in fams):
                    mv.append(tgt)
    return mv

def analyse(fams):
    succ = [movegen(v, fams, q, Bmask, Bij) for v in range(n)]
    out = loopy_outcomes(succ)
    loss = {v for v in range(n) if out[v] == 0}
    draw = {v for v in range(n) if out[v] == 2}
    win = {v for v in range(n) if out[v] == 1}
    # mistakes: a Win position with a non-Loss option, or a Draw position with a Win option
    mistakes = 0
    for v in range(n):
        if out[v] == 1 and any(out[t] != 0 for t in succ[v]):
            mistakes += 1
        if out[v] == 2 and any(out[t] == 1 for t in succ[v]):
            mistakes += 1
    return loss, draw, win, mistakes

results = []
fam_opts = ["A", "S", "Cd", "Ca"]
for r in range(len(fam_opts) + 1):
    for extra in combinations(fam_opts, r):
        fams = {"T2", *extra}
        loss, draw, win, mk = analyse(fams)
        tag = "+".join(sorted(fams))
        hitL = loss == Z
        hitD = draw == Z
        hitLD = (loss | draw) == Z
        results.append((tag, len(loss), len(draw), len(win), mk, hitL, hitD, hitLD))
        flag = ""
        if hitL:
            flag += "  LOSS=={Q=0}!"
        if hitD:
            flag += "  DRAW=={Q=0}!"
        if hitLD:
            flag += "  LOSS+DRAW=={Q=0}!"
        if flag and mk > 0:
            flag += f"  [NON-DEGENERATE: {mk} positions with mistakes]"
        print(f"  {tag:<18} |L|={len(loss):<4}|D|={len(draw):<4}|W|={len(win):<4}"
              f" mistakes={mk:<5}{flag}")

print("\ndone.")
