"""Exhaustive sweep of the MINIMAL Tier-2 stratum over bent Gold components.

Position space V = F_2^8 = F_256 (nimber field), bent components
Q(v) = Tr(lambda * v^{1+2^a}), a = 1 (APN). The Tier-2 access model grants,
per candidate move, O(1) oracle bits: the per-coin framing q_i = Q(e_i) and
polar evaluations B(.,.). The minimal stratum is:

  S1 (acyclic single-bit): turn OFF a set bit i, legality = f(q_i, B(v,e_i))
     for a fixed gate f: F_2^2 -> F_2  (16 gates; includes bent_route's
     Rule A = XOR gate [spin-flip] and Rule B = projection-on-B gate).
  S2 (loopy single-bit): flip bit i either way, same gates; Loss- AND
     Draw-sets both tested (loopy_quadric's third-outcome escape).
  S3 (leading-coin descent): move to any w < v, legality = f(q_j, B(v, v^w))
     with j = the highest changed bit (the coin-turning 'leading turned coin'
     convention).  16 gates, sampled lambdas.

For every bent lambda and every gate: does the P-set (or Loss/Draw set)
equal {Q=0}?  Also: is it a quadric at all (ANF fit), and with which Arf?
"""
import sys
sys.path.insert(0, "/Users/a9lim/Work/ogdoad/experiments")
import ogdoad as pl

M, A = 8, 1
N = 1 << M


def nim_trace_val(x):
    acc = pl.Nimber(x)
    t = pl.Nimber(x)
    for _ in range(M - 1):
        t = t * t
        acc = acc + t
    return acc.value


def qtab_for(lam):
    L = pl.Nimber(lam)
    tab = []
    for v in range(N):
        X = pl.Nimber(v)
        tab.append(nim_trace_val((L * X * (X * X)).value))
    return tab


# ---------------------------------------------------------------- outcome solvers

def pset_acyclic(succ):
    """Loss-set of a DAG game where every move strictly decreases the position."""
    loss = [False] * N
    for v in range(N):
        loss[v] = not any(loss[w] for w in succ[v])
    return frozenset(v for v in range(N) if loss[v])


def outcomes_loopy(succ):
    """Retrograde Win/Loss/Draw on an arbitrary finite graph (kernel.rs port)."""
    from collections import deque
    pred = [[] for _ in range(N)]
    for u in range(N):
        for v in succ[u]:
            pred[v].append(u)
    remaining = [len(succ[v]) for v in range(N)]
    label = [None] * N
    q = deque()
    for v in range(N):
        if not succ[v]:
            label[v] = "L"
            q.append(v)
    while q:
        v = q.popleft()
        for u in pred[v]:
            if label[u] is not None:
                continue
            if label[v] == "L":
                label[u] = "W"
                q.append(u)
            else:
                remaining[u] -= 1
                if remaining[u] == 0:
                    label[u] = "L"
                    q.append(u)
    loss = frozenset(v for v in range(N) if label[v] == "L")
    draw = frozenset(v for v in range(N) if label[v] is None)
    return loss, draw


def fit_quadric(points):
    """ANF Mobius fit: is `points` the zero set of a quadratic form? -> (ok, deg>2?)
    Returns None if not quadratic, else (qd, pairs) of the fitted form."""
    coeffs = [1] * N
    for v in points:
        coeffs[v] = 0
    for i in range(M):
        bit = 1 << i
        for mask in range(N):
            if mask & bit:
                coeffs[mask] ^= coeffs[mask ^ bit]
    if any(c and bin(mask).count("1") > 2 for mask, c in enumerate(coeffs)):
        return None
    qd = [coeffs[1 << i] for i in range(M)]
    pairs = [(i, j) for i in range(M) for j in range(i + 1, M)
             if coeffs[(1 << i) | (1 << j)]]
    return qd, pairs


def arf_of(qd, pairs):
    q = [pl.Nimber(x) for x in qd]
    b = {(i, j): pl.Nimber(1) for (i, j) in pairs}
    return pl.arf_nimber(pl.NimberAlgebra(q=q, b=b))


# ---------------------------------------------------------------- the sweep

bent_lams = []
for lam in range(1, N):
    tab = qtab_for(lam)
    z = tab.count(0)
    if z in (N // 2 - (1 << (M // 2 - 1)), N // 2 + (1 << (M // 2 - 1))):
        bent_lams.append((lam, tab))
print(f"bent components of Tr(lambda x^3) over F_256: {len(bent_lams)} "
      f"(classical count 2(2^m-1)/3 = {2*(N-1)//3})")

GATES = list(range(16))           # f(q,b) = (g >> (2*q + b)) & 1


def gate(g, q, b):
    return (g >> ((q << 1) | b)) & 1


hits = {"S1": [], "S2L": [], "S2D": [], "S3": []}
s1_quadric_stats = {}
qblind_stats = {"quadric": 0, "right_arf": 0, "exact": 0, "tested": 0}

for lam, tab in bent_lams:
    zero = frozenset(v for v in range(N) if tab[v] == 0)
    qd = [tab[1 << i] for i in range(M)]
    target_fit = fit_quadric(zero)
    target_arf = arf_of(*target_fit).arf

    def Bve(v, i):
        return tab[v ^ (1 << i)] ^ tab[v] ^ qd[i]

    Btab = [[Bve(v, i) for i in range(M)] for v in range(N)]

    for g in GATES:
        # S1: acyclic single-bit turn-off
        succ = [[v ^ (1 << i) for i in range(M)
                 if (v >> i) & 1 and gate(g, qd[i], Btab[v][i])]
                for v in range(N)]
        P = pset_acyclic(succ)
        if P == zero:
            hits["S1"].append((lam, g))
        fit = fit_quadric(P)
        key = (g, "quadric" if fit and fit[1] else
               ("affine" if fit else "non-quadric"))
        s1_quadric_stats[key] = s1_quadric_stats.get(key, 0) + 1
        # q-blind gates: f(q,b) independent of q  <=>  gate rows equal
        if gate(g, 0, 0) == gate(g, 1, 0) and gate(g, 0, 1) == gate(g, 1, 1) \
           and g == 0b0100_0100 & 0xF or False:
            pass
        # S2: loopy single-bit (either direction)
        succ2 = [[v ^ (1 << i) for i in range(M)
                  if gate(g, qd[i], Btab[v][i])] for v in range(N)]
        loss, draw = outcomes_loopy(succ2)
        if loss == zero:
            hits["S2L"].append((lam, g))
        if draw == zero:
            hits["S2D"].append((lam, g))

    # q-blind baseline (bent_route Rule B): legality = B(v,e_i), turn-off only
    succB = [[v ^ (1 << i) for i in range(M)
              if (v >> i) & 1 and Btab[v][i]] for v in range(N)]
    PB = pset_acyclic(succB)
    qblind_stats["tested"] += 1
    fitB = fit_quadric(PB)
    if fitB and fitB[1]:
        qblind_stats["quadric"] += 1
        if arf_of(*fitB).arf == target_arf:
            qblind_stats["right_arf"] += 1
    if PB == zero:
        qblind_stats["exact"] += 1

print(f"\nS1 acyclic single-bit, all 16 gates x {len(bent_lams)} bent lambdas: "
      f"exact {{Q=0}} hits: {len(hits['S1'])}")
print(f"S2 loopy single-bit:  Loss-set hits: {len(hits['S2L'])}, "
      f"Draw-set hits: {len(hits['S2D'])}")

print("\nq-blind baseline (gate = B only, bent_route Rule B) across bent lambdas:")
print(f"  P-set genuinely quadratic: {qblind_stats['quadric']}/{qblind_stats['tested']}, "
      f"of those with the target's Arf: {qblind_stats['right_arf']}, "
      f"exact {{Q=0}}: {qblind_stats['exact']}")

# S3: leading-coin descent, sampled lambdas
sample = bent_lams[:24]
for lam, tab in sample:
    zero = frozenset(v for v in range(N) if tab[v] == 0)
    qd = [tab[1 << i] for i in range(M)]

    def Bvd(v, d):
        return tab[v ^ d] ^ tab[v] ^ tab[d]

    for g in GATES:
        succ = []
        for v in range(N):
            row = []
            for w in range(v):
                d = v ^ w
                j = d.bit_length() - 1
                if gate(g, qd[j], Bvd(v, d)):
                    row.append(w)
            succ.append(row)
        P = pset_acyclic(succ)
        if P == zero:
            hits["S3"].append((lam, g))
print(f"\nS3 leading-coin descent, 16 gates x {len(sample)} sampled bent lambdas: "
      f"exact hits: {len(hits['S3'])}")

print("\nhits detail:", {k: v for k, v in hits.items() if v})
