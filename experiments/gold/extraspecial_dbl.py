"""Double-touch variant: a move is EITHER a single touch (ko=self) OR a
double-touch of a fresh coin (t 0->2 in one move, charge q_i + B-parity twice
= q_i exactly).  E-natural: 'play a generator or its square'.
Tests: reduced (p,s) game, k=3/k=4 full pattern tables, m=4 sweep.
"""
import itertools
from extraspecial_core import nim_mul, frob, tr, gold_q, polar, validate

class EchoDbl:
    def __init__(self, qbits, Bhigh, allow_dbl=True):
        self.k = len(qbits); self.q = qbits; self.Bh = Bhigh
        self.dbl = allow_dbl
    def solve(self, maxfirst=True):
        k, q, Bh = self.k, self.q, self.Bh
        full = (1 << k) - 1
        memo = {}
        choice = [0]
        def val(open_m, done_m, last, mover, sigma):
            if done_m == full: return sigma
            key = (open_m, done_m, last, mover, sigma)
            if key in memo: return memo[key]
            wantmax = (mover == 0) == maxfirst
            outs = []
            for i in range(k):
                if (done_m >> i) & 1: continue
                opened = (open_m >> i) & 1
                if i != last:
                    # single touch
                    if opened:
                        ch = q[i] ^ (bin(open_m & Bh[i]).count('1') & 1)
                        outs.append(val(open_m & ~(1 << i), done_m | (1 << i),
                                        i, 1 - mover, sigma ^ ch))
                    else:
                        ch = bin(open_m & Bh[i]).count('1') & 1
                        outs.append(val(open_m | (1 << i), done_m,
                                        i, 1 - mover, sigma ^ ch))
                if self.dbl and not opened and i != last:
                    # double touch of fresh coin: charge q_i (B-parity cancels)
                    outs.append(val(open_m, done_m | (1 << i),
                                    i, 1 - mover, sigma ^ q[i]))
            if not outs:
                v = val(open_m, done_m, -1, 1 - mover, sigma)  # pass clears ko
                memo[key] = v
                return v
            v = max(outs) if wantmax else min(outs)
            if len(set(outs)) > 1: choice[0] += 1
            memo[key] = v
            return v
        return val(0, 0, -1, 0, 0), choice[0]

def mk_Bh_edges(k, edges):
    Bh = []
    for i in range(k):
        m = 0
        for (a, b) in edges:
            if a == i and b > i: m |= 1 << b
            elif b == i and a > i: m |= 1 << a
        Bh.append(m)
    return Bh

print("=== reduced game v(p,s), double-touch variant ===")
for p in range(1, 5):
    row = []
    for s in range(0, 4):
        k = 2 * p + s
        if k > 9: row.append(" - "); continue
        edges = [(2 * i, 2 * i + 1) for i in range(p)]
        vmax, _ = EchoDbl([0] * k, mk_Bh_edges(k, edges)).solve(True)
        vmin, _ = EchoDbl([0] * k, mk_Bh_edges(k, edges)).solve(False)
        t = p & 1
        row.append(f"{vmax}{vmin}{'OK ' if vmax == t == vmin else 'BAD'}")
    print(f"p={p}: " + " | ".join(row))

print("\n=== k=3, k=4 full pattern tables, double-touch ===")
for k in (3, 4):
    prs = list(itertools.combinations(range(k), 2))
    for mf in (True, False):
        bad = 0; tot = 0; ndg = 0
        for qm in range(1 << k):
            q = [(qm >> i) & 1 for i in range(k)]
            for bm in range(1 << len(prs)):
                B = {prs[x]: (bm >> x) & 1 for x in range(len(prs))}
                Bh = []
                for i in range(k):
                    m_ = 0
                    for j in range(k):
                        if j > i and B.get((i, j), 0): m_ |= 1 << j
                    Bh.append(m_)
                t = (sum(q) + sum(B.values())) & 1
                tot += 1
                v, ch = EchoDbl(q, Bh).solve(mf)
                if v != t: bad += 1
                elif ch: ndg += 1
        print(f"k={k} {'P1max' if mf else 'P1min'}: solved {tot-bad}/{tot} nondeg={ndg}")

print("\n=== m=4 sweep, double-touch, P1max & P1min ===")
validate()
for lam in range(1, 16):
    Q = gold_q(4, 1, lam); B = polar(Q, 4)
    res = {}
    for mf in (True, False):
        agree = 0
        for x in range(16):
            S = [i for i in range(4) if (x >> i) & 1]
            if not S:
                v = 0
            else:
                qb = [Q[1 << c] for c in S]
                Bh = []
                for li, ci in enumerate(S):
                    m_ = 0
                    for lj, cj in enumerate(S):
                        if cj > ci and ((B[ci] >> cj) & 1): m_ |= 1 << lj
                    Bh.append(m_)
                v, _ = EchoDbl(qb, Bh).solve(mf)
            agree += (v == Q[x])
        res[mf] = agree
    print(f"lam={lam:2d}: P1max {res[True]:2d}/16  P1min {res[False]:2d}/16")
