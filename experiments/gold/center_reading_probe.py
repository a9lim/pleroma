"""Exploratory sweep: center-READING translation rules on the extraspecial E.

The no-go theorems kill (i) pure-Cayley rules on E (I*I = E) and (ii) center-blind
rules (outcome-inert covering). The live regime is translation moves gated by
local center+frame data. This sweeps a principled small family on the bent Gold
component (m=8, F_2^8, the clean r=4 nonsingular case) and on the m=4 a=1 Gold
form (degenerate, rank 2), reporting whether any member's Loss-set is
I = pre-image of {Q=0}, or projects to {Q=0}.

Moves from g=(a,v): left-mult by (0,e_i) for set bits i (turn coin i off):
    (a, v) -> (a ^ beta(e_i, v), v ^ e_i),  beta(e_i,v) = q_i v_i + sum_{j<i} B_ij v_j
Gates (legality predicates, local in (a, coin i, B-row of i)):
    G0 always                      G1 a == 0
    G2 a == 1                      G3 a == q_i
    G4 a == B(v,e_i)               G5 B(v,e_i) == 1   (center-blind baseline)
    G6 a == q_i ^ B(v,e_i) (= local Ising Delta-Q gate vs center)
Optional central move Z: (a,v) -> (a^1, v) (multiply by z; introduces 2-cycles).
"""

import sys
sys.path.insert(0, "/Users/a9lim/Work/ogdoad/experiments")
import ogdoad as pl
from common import gold as gold_u, polar as polar_u

def frob(x, a):
    for _ in range(a):
        x = x * x
    return x

def nim_trace(x, m):
    acc = pl.Nimber(x); t = pl.Nimber(x)
    for _ in range(m - 1):
        t = t * t
        acc = acc + t
    return acc.value

def bent_gold(v, lam, a, m):
    x = pl.Nimber(v)
    return nim_trace((pl.Nimber(lam) * x * frob(x, a)).value, m)

def loopy_solve(succ):
    """Win/Loss/Draw retrograde solver (normal play): Loss = all opts Win;
    terminals Loss."""
    n = len(succ)
    pred = [[] for _ in range(n)]
    deg = [0] * n
    for u, opts in enumerate(succ):
        deg[u] = len(opts)
        for w in opts:
            pred[w].append(u)
    label = ["D"] * n
    from collections import deque
    queue = deque()
    remaining = deg[:]
    for u in range(n):
        if deg[u] == 0:
            label[u] = "L"
            queue.append(u)
    while queue:
        w = queue.popleft()
        for u in pred[w]:
            if label[u] != "D":
                continue
            if label[w] == "L":
                label[u] = "W"
                queue.append(u)
            else:
                remaining[u] -= 1
                if remaining[u] == 0:
                    label[u] = "L"
                    queue.append(u)
    return label

def run_world(name, m, Q):
    n = 1 << m
    Qv = [Q(v) for v in range(n)]
    zeros = set(v for v in range(n) if Qv[v] == 0)
    qd = [Qv[1 << i] for i in range(m)]
    Bbit = [[Qv[(1 << i) ^ (1 << j)] ^ qd[i] ^ qd[j] if i != j else 0
             for j in range(m)] for i in range(m)]
    def Brow(v, i):
        acc = 0
        for j in range(m):
            if (v >> j) & 1:
                acc ^= Bbit[i][j]
        return acc
    def beta_ei(i, v):  # beta(e_i, v), standard cocycle
        acc = qd[i] if (v >> i) & 1 else 0
        for j in range(i):
            if (v >> j) & 1:
                acc ^= Bbit[i][j]
        return acc

    I = set(v for v in zeros) # projection target; I itself = both lifts
    gates = {
        "G0 always":        lambda a, v, i: True,
        "G1 a==0":          lambda a, v, i: a == 0,
        "G2 a==1":          lambda a, v, i: a == 1,
        "G3 a==q_i":        lambda a, v, i: a == qd[i],
        "G4 a==B(v,e_i)":   lambda a, v, i: a == Brow(v, i),
        "G5 B(v,e_i)==1":   lambda a, v, i: Brow(v, i) == 1,
        "G6 a==q_i^B(v,e_i)": lambda a, v, i: a == (qd[i] ^ Brow(v, i)),
    }
    print(f"--- {name}: |V|={n}, |zeros|={len(zeros)} ---")
    print(f"{'gate':<22}{'Z':<3}{'|Loss|':>7} {'Loss==I':>8} {'pi(Loss)=={Q=0}':>16} "
          f"{'sec0 agree':>11} {'draws':>6}")
    for gname, gate in gates.items():
        for with_z in (False, True):
            succ = []
            for g in range(2 * n):
                a, v = divmod(g, n)
                opts = []
                for i in range(m):
                    if (v >> i) & 1 and gate(a, v, i):
                        opts.append((a ^ beta_ei(i, v)) * n + (v ^ (1 << i)))
                if with_z:
                    opts.append((a ^ 1) * n + v)
                succ.append(opts)
            lab = loopy_solve(succ)
            loss = [g for g in range(2 * n) if lab[g] == "L"]
            loss_set = set(loss)
            is_I = loss_set == set(av * n + v for av in (0, 1) for v in zeros)
            proj = set(g % n for g in loss)
            proj_ok = proj == zeros
            sec0 = set(g % n for g in loss if g < n)
            agree0 = sum(1 for v in range(n) if (v in sec0) == (v in zeros))
            draws = sum(1 for x in lab if x == "D")
            print(f"{gname:<22}{'+Z' if with_z else '  ':<3}{len(loss):>7} "
                  f"{str(is_I):>8} {str(proj_ok):>16} {agree0:>7}/{n:<3} {draws:>6}")
    print()

# bent Gold component, m=8 (lambda=2 found by the bent witness search; Arf 0, r=4)
run_world("bent Gold Tr(2*x^3), F_2^8", 8, lambda v: bent_gold(v, 2, 1, 8))
# the unscaled Gold m=4,a=1 (rank 2, radical dim 2) for contrast
run_world("Gold Q_1, F_2^4", 4, lambda v: gold_u(v, 1, 4))
