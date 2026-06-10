"""Core: standalone nim arithmetic + Gold forms + correct ECHO-ko solver.

Validated against ogdoad's pinned values:
  nim products 2*2=3, 2*4=8, 16*16=24
  Gold zero counts: (4,1):4  (8,1):112  (8,2):96  bent (8,1,lam=2):136
Solver: state includes accumulated sigma (the round-1 bug was omitting it).
"""
import sys
from functools import lru_cache

sys.setrecursionlimit(100000)

# ---------- nim arithmetic ----------
_nm = {}
def nim_mul(a, b):
    if a > b: a, b = b, a          # a <= b
    if a < 2: return a * b
    key = (a, b)
    v = _nm.get(key)
    if v is not None: return v
    # largest Fermat 2-power F = 2^(2^k) <= b
    k = 0
    while (1 << (2 << k)) <= b:    # 2^(2^(k+1)) <= b
        k += 1
    F = 1 << (1 << k)              # F <= b < F*F
    bh, bl = b >> (1 << k), b & (F - 1)
    if a < F:
        r = (nim_mul(a, bh) << (1 << k)) ^ nim_mul(a, bl)
    else:
        ah, al = a >> (1 << k), a & (F - 1)
        t1 = nim_mul(ah, bh)
        t2 = nim_mul(ah, bl) ^ nim_mul(al, bh)
        t3 = nim_mul(al, bl)
        r = ((t1 ^ t2) << (1 << k)) ^ t3 ^ nim_mul(t1, F >> 1)
    _nm[key] = r
    return r

def frob(x, a):
    for _ in range(a):
        x = nim_mul(x, x)
    return x

def tr(x, m):
    s, y = 0, x
    for _ in range(m):
        s ^= y
        y = nim_mul(y, y)
    assert s in (0, 1), (x, m, s)
    return s

def gold_q(m, a, lam=1):
    """Q(x) = Tr(lam * x^(1+2^a)) as a list over F_2^m."""
    return [tr(nim_mul(lam, nim_mul(x, frob(x, a))), m) for x in range(1 << m)]

def polar(Q, m):
    """B(u,v) = Q(u^v)+Q(u)+Q(v) as dict of row masks: Brow[i] = mask of j with B(e_i,e_j)=1."""
    rows = []
    for i in range(m):
        row = 0
        for j in range(m):
            if i == j: continue
            b = Q[(1 << i) ^ (1 << j)] ^ Q[1 << i] ^ Q[1 << j]
            row |= b << j
        rows.append(row)
    return rows

# ---------- validations ----------
def validate():
    assert nim_mul(2, 2) == 3 and nim_mul(2, 4) == 8 and nim_mul(16, 16) == 24
    assert sum(1 for v in gold_q(4, 1) if v == 0) == 4
    assert sum(1 for v in gold_q(8, 1) if v == 0) == 112
    assert sum(1 for v in gold_q(8, 2) if v == 0) == 96
    assert sum(1 for v in gold_q(8, 1, 2) if v == 0) == 136
    print("nim/form validations OK")

# ---------- ECHO-ko game, correct solver ----------
# Position x: support coins S (global indices). Each coin touched exactly twice.
# State: (open_mask, done_mask, last, mover, sigma) over local indices 0..k-1.
# Touch i: charge = q[i] if i open else 0, XOR parity(open & Bhigh[i])
#   (triangular cocycle c(u,v) = sum_i q_i u_i v_i + sum_{k>j} B_kj u_k v_j;
#    final value cocycle-choice independent: every coin touched twice.)
# ko variants: 'self' = may not touch the coin just touched (last);
#              'none' = no ko; 'opp' = may not touch the coin the OPPONENT
#              last touched (but may re-touch own); 'w2' = last two touches banned.
# stuck -> pass (clears ko memory); all done -> payoff sigma.
# orientation: maxfirst=True: P1 (mover 0) maximizes sigma.

class Echo:
    def __init__(self, qbits, Bhigh, ko='self'):
        self.k = len(qbits)
        self.q = qbits
        self.Bh = Bhigh   # Bhigh[i] = mask of j with B(i,j)=1 and S[j]>S[i] (local idx, S sorted)
        self.ko = ko

    def solve(self, maxfirst=True):
        k, q, Bh, ko = self.k, self.q, self.Bh, self.ko
        full = (1 << k) - 1
        memo = {}
        choice_states = [0]   # states with >=2 legal moves whose child values differ
        def banned(last, mover):
            # returns mask of banned coins given ko memory 'last' and current mover
            if ko == 'none': return 0
            if ko == 'self':
                return 0 if last[0] < 0 else (1 << last[0])
            if ko == 'opp':
                # last = (last_by_p0, last_by_p1); banned = opponent's last
                lb = last[1 - mover]
                return 0 if lb < 0 else (1 << lb)
            if ko == 'w2':
                m_ = 0
                for l in last:
                    if l >= 0: m_ |= 1 << l
                return m_
            raise ValueError(ko)
        def init_last():
            if ko == 'self': return (-1,)
            if ko == 'none': return ()
            if ko == 'opp': return (-1, -1)
            if ko == 'w2': return (-1, -1)
            raise ValueError(ko)
        def upd_last(last, i, mover):
            if ko == 'none': return ()
            if ko == 'self': return (i,)
            if ko == 'opp':
                l = list(last); l[mover] = i; return tuple(l)
            if ko == 'w2': return (last[1], i)
            raise ValueError(ko)
        def clear_last(last, mover):
            if ko == 'opp':
                # pass clears the ko against the passer's opponent? clear all.
                return (-1, -1)
            return init_last()
        def val(open_m, done_m, last, mover, sigma):
            if done_m == full:
                return sigma
            key = (open_m, done_m, last, mover, sigma)
            v = memo.get(key)
            if v is not None: return v
            avail = full & ~done_m & ~banned(last, mover)
            # legal: coin i with t_i < 2  <=> not done
            legal = [i for i in range(k) if (avail >> i) & 1]
            if not legal:
                # stuck: pass, clear ko
                v = val(open_m, done_m, clear_last(last, mover), 1 - mover, sigma)
                memo[key] = v
                return v
            wantmax = (mover == 0) == maxfirst
            best = None
            vals = set()
            for i in legal:
                if (open_m >> i) & 1:   # second touch
                    ch = q[i] ^ (bin(open_m & Bh[i]).count('1') & 1)
                    no, nd = open_m & ~(1 << i), done_m | (1 << i)
                else:                   # first touch
                    ch = bin(open_m & Bh[i]).count('1') & 1
                    no, nd = open_m | (1 << i), done_m
                cv = val(no, nd, upd_last(last, i, mover), 1 - mover, sigma ^ ch)
                vals.add(cv)
                if best is None: best = cv
                elif wantmax: best = max(best, cv)
                else: best = min(best, cv)
            if len(vals) > 1:
                choice_states[0] += 1
            memo[key] = best
            return best
        v0 = val(0, 0, init_last(), 0, 0)
        return v0, choice_states[0]

def echo_value(x, Q, Brows, m, ko='self', maxfirst=True, qover=None, Bover=None):
    """Game value for position x under form (Q,Brows) on the bit frame,
    or with overridden frame data (qover list, Bover matrix) for normal frames."""
    S = [i for i in range(m) if (x >> i) & 1]
    k = len(S)
    if k == 0: return 0, 0
    if qover is None:
        qb = [Q[1 << c] for c in S]
        Bh = []
        for li, ci in enumerate(S):
            mask = 0
            for lj, cj in enumerate(S):
                if cj > ci and ((Brows[ci] >> cj) & 1):
                    mask |= 1 << lj
            Bh.append(mask)
    else:
        qb = [qover[c] for c in S]
        Bh = []
        for li, ci in enumerate(S):
            mask = 0
            for lj, cj in enumerate(S):
                if cj > ci and ((Bover[ci] >> cj) & 1):
                    mask |= 1 << lj
            Bh.append(mask)
    return Echo(qb, Bh, ko=ko).solve(maxfirst=maxfirst)

if __name__ == '__main__':
    validate()
