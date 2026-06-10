"""Independent skeptic solver: ECHO-ko from the PROSE rules, written fresh.

Differences from extraspecial_core.Echo on purpose:
- charge q_i on the FIRST touch (not second); B-charge uses LOWER-triangular rows
  (B(i,j) with j < i) instead of upper. Same cocycle class; final sigma of a
  complete play must agree, so game values must agree.
- state = tuple of per-coin touch counts (not bitmasks), explicit last-touched.
- plain dict memo keyed on (counts, last, mover, sigma).
Compares against the original engine on the (8,2) adapted frame, all 256 supports,
and against the true Gold Q.
"""
import sys
sys.path.insert(0, '/tmp')
from extraspecial_core import gold_q, echo_value
from extraspecial_adapted import build_adapted_frame

def solve_indep(qb, Bfull, k):
    # Bfull[i] = set of j with B(i,j)=1 (local indices)
    from functools import lru_cache
    memo = {}
    def val(counts, last, mover, sigma):
        if all(c == 2 for c in counts):
            return sigma
        key = (counts, last, mover, sigma)
        if key in memo: return memo[key]
        legal = [i for i in range(k) if counts[i] < 2 and i != last]
        if not legal:
            v = val(counts, -1, 1 - mover, sigma)  # pass clears ko
            memo[key] = v
            return v
        # open set = coins touched exactly once
        openset = [i for i in range(k) if counts[i] == 1]
        best = None
        for i in legal:
            # lower-triangular cocycle: B-charge counts open j < i with B(i,j)=1
            ch = sum(1 for j in openset if j < i and j in Bfull[i]) & 1
            if counts[i] == 0:
                ch ^= qb[i]   # q on FIRST touch (other section)
            nc = list(counts); nc[i] += 1
            cv = val(tuple(nc), i, 1 - mover, sigma ^ ch)
            if best is None: best = cv
            elif mover == 0: best = max(best, cv)
            else: best = min(best, cv)
        memo[key] = best
        return best
    return val(tuple([0] * k), -1, 0, 0)

Q = gold_q(8, 2, 1)
m = 8
frame, pairs, radb, arf1 = build_adapted_frame(Q, m)
qover = [Q[v] for v in frame]
Bover = []
for i in range(m):
    row = 0
    for j in range(m):
        if i == j: continue
        row |= (Q[frame[i] ^ frame[j]] ^ Q[frame[i]] ^ Q[frame[j]]) << j
    Bover.append(row)

mismatch_engine = 0
mismatch_Q = 0
for cm in range(1 << m):
    S = [i for i in range(m) if (cm >> i) & 1]
    k = len(S)
    xf = 0
    for i in S: xf ^= frame[i]
    if k == 0:
        v_ind = 0
    else:
        qb = [qover[c] for c in S]
        Bfull = [set(lj for lj, cj in enumerate(S) if cj != ci and ((Bover[ci] >> cj) & 1))
                 for li, ci in enumerate(S)]
        v_ind = solve_indep(qb, Bfull, k)
    v_eng, _ = echo_value(cm, None, None, m, ko='self', maxfirst=True,
                          qover=qover, Bover=Bover)
    if v_ind != v_eng: mismatch_engine += 1
    if v_ind != Q[xf]: mismatch_Q += 1
print(f"(8,2) adapted frame: indep-vs-engine mismatches {mismatch_engine}/256, "
      f"indep-vs-Gold-Q mismatches {mismatch_Q}/256")
