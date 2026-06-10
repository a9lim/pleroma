"""Synthesis-round verification of the two load-bearing round-1 skeptic claims.

1. ECHO-ko with a PARITY-CORRECT solver (state includes accumulated charge):
   - validate solver vs no-memo explicit tree at m=4
   - (4,1,1) exact both orientations; bent (4,1,2) exact P1-max
   - F_16 solved set (skeptic claimed {1,2,12,13,14})
   - (8,1,1), (8,2,1), bent (8,1,*): agreement counts (skeptic: (8,2) 255/256)
   - decision-nondegeneracy: does some reachable state have moves to both values?
2. Diagonal game-native source (diag skeptic): q_i^(m,1) = Tr(P(w) * e_i) with
   w = XOR of Fermat coins 2^(2^t), t>=1, P(w) = w*w XOR w. Check m=8,16,32.

Independent nim arithmetic (Conway recursion); cross-checked against repo-pinned
values and goldarf.tex zero counts.
"""
import sys
from functools import lru_cache

sys.setrecursionlimit(100000)


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
    return (high * F) ^ low


assert nim_mul(2, 2) == 3 and nim_mul(2, 4) == 8 and nim_mul(16, 16) == 24


def frob(x, a):
    for _ in range(a):
        x = nim_mul(x, x)
    return x


def trace(x, m):
    acc, t = x, x
    for _ in range(m - 1):
        t = nim_mul(t, t)
        acc ^= t
    assert acc in (0, 1)
    return acc


def make_form(m, a, lam):
    qtab = [trace(nim_mul(lam, nim_mul(v, frob(v, a))), m) for v in range(1 << m)]
    qd = [qtab[1 << i] for i in range(m)]
    Bm = [[qtab[(1 << i) ^ (1 << j)] ^ qtab[1 << i] ^ qtab[1 << j] for j in range(m)]
          for i in range(m)]
    return qtab, qd, Bm


def charge_move(o, i, qd, Bm, lower=True):
    """c(o, e_i) with triangular cocycle (lower: B_kj for k>j contributes)."""
    acc = qd[i] if (o >> i) & 1 else 0
    for k in range(len(qd)):
        if k == i or not (o >> k) & 1:
            continue
        if (lower and k > i) or (not lower and k < i):
            acc ^= Bm[k][i]
    return acc


# ------------------------- parity-correct solver -------------------------

def solve(x, m, qd, Bm, p1_max, lower=True, memo_on=True, nondeg_out=None):
    """Final sigma under optimal play. P1 wants sigma=1 iff p1_max.
    State: (u untouched, o open, last touched coin (-1 after pass/start),
            mover_is_p1, accumulated sigma)."""
    bits = [i for i in range(m) if (x >> i) & 1]
    memo = {}

    def rec(u, o, last, p1, s):
        if u == 0 and o == 0:
            return s
        key = (u, o, last, p1, s)
        if memo_on and key in memo:
            return memo[key]
        legal = []
        for i in bits:
            if i == last:
                continue
            if (u >> i) & 1:
                legal.append((i, u ^ (1 << i), o ^ (1 << i)))
            elif (o >> i) & 1:
                legal.append((i, u, o ^ (1 << i)))
        if not legal:
            res = rec(u, o, -1, not p1, s)
        else:
            want = 1 if (p1 == p1_max) else 0
            vals = []
            for (i, u2, o2) in legal:
                ch = charge_move(o, i, qd, Bm, lower)
                vals.append(rec(u2, o2, i, not p1, s ^ ch))
                if vals[-1] == want:
                    break
            res = want if want in vals else 1 - want
            if nondeg_out is not None and len(set(vals)) > 1:
                nondeg_out[0] = True
        if memo_on:
            memo[key] = res
        return res

    return rec(x, 0, -1, True, 0)


def solve_all(m, qd, Bm, p1_max, lower=True):
    return [solve(x, m, qd, Bm, p1_max, lower) for x in range(1 << m)]


def nondeg_table(m, qd, Bm, p1_max, lower=True):
    """For each x: does ANY reachable state have legal moves to both game values?"""
    out = []
    for x in range(1 << m):
        flag = [False]
        solve(x, m, qd, Bm, p1_max, lower, nondeg_out=flag)
        out.append(flag[0])
    return out


print("== solver validation: memo vs explicit tree, m=4, all forms/orientations ==")
for lam in (1, 2, 3):
    qtab, qd, Bm = make_form(4, 1, lam)
    for p1m in (True, False):
        for x in range(16):
            v_memo = solve(x, 4, qd, Bm, p1m, memo_on=True)
            v_tree = solve(x, 4, qd, Bm, p1m, memo_on=False)
            assert v_memo == v_tree, (lam, p1m, x)
print("   memo == explicit tree on all 16 positions x 3 forms x 2 orientations: OK")

print("\n== goldarf.tex zero-count cross-checks ==")
for (m, a, lam, expect) in [(4, 1, 1, 4), (8, 1, 1, 112), (8, 2, 1, 96)]:
    qtab, _, _ = make_form(m, a, lam)
    z = sum(1 for t in qtab if t == 0)
    print(f"   (m={m},a={a},lam={lam}): |Q=0| = {z} (expect {expect})")
    assert z == expect

print("\n== F_16 sweep (4,1,lam), ko=1, lower cocycle, corrected solver ==")
solved = {}
for lam in range(1, 16):
    qtab, qd, Bm = make_form(4, 1, lam)
    res = {}
    for p1m in (True, False):
        tab = solve_all(4, qd, Bm, p1m)
        res["P1max" if p1m else "P1min"] = (tab == qtab,
                                            sum(t == q for t, q in zip(tab, qtab)))
    solved[lam] = res
    hit = [k for k, (ok, _) in res.items() if ok]
    z = sum(1 for t in qtab if t == 0)
    print(f"   lam={lam:2d} |Q=0|={z:2d}  exact:{hit if hit else '-'}  "
          f"agree: P1max {res['P1max'][1]}/16, P1min {res['P1min'][1]}/16")
exact_set = sorted(l for l, r in solved.items() if any(ok for ok, _ in r.values()))
print(f"   exact-solved lam set: {exact_set}  (skeptic claimed [1, 2, 12, 13, 14])")

print("\n== m=8 cases, corrected solver, ko=1 ==")
for (m, a, lam, label) in [(8, 1, 1, "Gold (8,1) rank6"), (8, 2, 1, "Gold (8,2) rank4")]:
    qtab, qd, Bm = make_form(m, a, lam)
    for lower in (True, False):
        for p1m in (True, False):
            tab = solve_all(m, qd, Bm, p1m, lower)
            ag = sum(t == q for t, q in zip(tab, qtab))
            miss = [x for x in range(256) if tab[x] != qtab[x]]
            pc = sorted(set(bin(x).count('1') for x in miss))
            print(f"   {label} cocycle={'lower' if lower else 'upper'} "
                  f"{'P1max' if p1m else 'P1min'}: agree {ag}/256"
                  + (f"  miss popcounts {pc} miss={miss[:6]}" if miss else "  == Q EXACTLY"))

# first bent component at (8,1)
for lam in range(1, 256):
    qtab, qd, Bm = make_form(8, 1, lam)
    z = sum(1 for t in qtab if t == 0)
    if z in (120, 136):
        for p1m in (True, False):
            tab = solve_all(8, qd, Bm, p1m)
            ag = sum(t == q for t, q in zip(tab, qtab))
            print(f"   bent (8,1,lam={lam}) |Q=0|={z} {'P1max' if p1m else 'P1min'}: "
                  f"agree {ag}/256" + ("  == Q EXACTLY" if ag == 256 else ""))
        break

print("\n== decision-nondegeneracy (mistakes exist?) on exact m=4 hits ==")
for lam in exact_set:
    qtab, qd, Bm = make_form(4, 1, lam)
    for p1m in (True, False):
        tab = solve_all(4, qd, Bm, p1m)
        if tab == qtab:
            nd = nondeg_table(4, qd, Bm, p1m)
            k = sum(nd)
            print(f"   lam={lam} {'P1max' if p1m else 'P1min'}: positions with a "
                  f"value-splitting reachable choice: {k}/16")

print("\n== diagonal game-native source: q_i^(m,1) == Tr(P(w) * e_i)? ==")
for m in (8, 16, 32):
    k = 1
    w = 0
    while (1 << (1 << k)) < (1 << m):
        w ^= 1 << (1 << k)  # Fermat coins 2^(2^t), t>=1
        k += 1
    Pw = nim_mul(w, w) ^ w
    lam_closed = 0
    t = 1
    while (1 << ((1 << t) - 1)) < (1 << m):
        lam_closed ^= 1 << ((1 << t) - 1)
        t += 1
    ok = all(trace(nim_mul(nim_mul(1 << i, frob(1 << i, 1)), 1), m) ==
             trace(nim_mul(Pw, 1 << i), m) for i in range(m))
    print(f"   m={m}: w={w}, P(w)={Pw}, closed-form lam={lam_closed}, "
          f"P(w)==lam: {Pw == lam_closed}, q_i match all i: {ok}")
print("\nall checks done")
