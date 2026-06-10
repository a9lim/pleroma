"""Final checks: k=9 full boards, board-8 class sample, m=16 spot checks,
decision-nondegeneracy + explicit position graph for the dummy game."""
import sys, time, random
from itertools import product
sys.path.insert(0, "/tmp")
from asym2_probe import make_form
from asym2_fifo import abstract_value, canon
from asym2_fifo_bench import direct_fifo_value
from asym2_dummy import support_value_dummy

sys.setrecursionlimit(1000000)

# ---- 1) full boards: popcount-8 supports -> board size 9, four forms, both t
print("full-board (pc=8 -> board 9) checks:")
for (m, a, lam) in [(8, 1, 1), (8, 2, 1), (8, 1, 2), (8, 1, 3)]:
    Q, qd, B = make_form(m, a, lam)
    S = tuple(range(8))
    for t in (1, 0):
        t0 = time.time()
        v = support_value_dummy(S, qd, B, t, m)
        ok = v == Q[255]
        print(f"  (m={m},a={a},lam={lam}) t={t}: value={v} Q={Q[255]} "
              f"{'OK' if ok else 'FAIL'} [{time.time()-t0:.0f}s]")

# ---- 2) board-8 sample: random k=7 patterns + dummy, both wants
print("\nboard-8 sample screen (random k=7 patterns + isolated dummy):")
rng = random.Random(2026)
pairs7 = [(i, j) for i in range(7) for j in range(i + 1, 7)]
fails = 0
t0 = time.time()
N = 40
for trial in range(N):
    bits = [rng.randint(0, 1) for _ in pairs7]
    edges = frozenset(p for (b, p) in zip(bits, pairs7) if b)
    par = len(edges) & 1
    for want in (0, 1):
        v = abstract_value(8, edges, want)  # vertex 7 = isolated dummy
        if v != par:
            fails += 1
            print(f"  FAIL: edges={sorted(edges)} want={want} v={v} par={par}")
print(f"  {N} random k=7 patterns x 2 wants: {fails} failures [{time.time()-t0:.0f}s]")

# ---- 3) m=16 spot checks (boards <= 7 are exhaustively-screened sizes)
print("\nm=16, a=1, lam=1 spot checks (random positions, popcount <= 6):")
Q16, qd16, B16 = make_form(16, 1, 1)
cnt = bad = 0
for _ in range(400):
    pc = rng.randint(1, 6)
    S = tuple(sorted(rng.sample(range(16), pc)))
    x = sum(1 << i for i in S)
    for t in (0, 1):
        v = support_value_dummy(S, qd16, B16, t, 16)
        cnt += 1
        if v != Q16[x]:
            bad += 1
            print(f"  MISS x={x} S={S} t={t}")
pc7 = 0
for _ in range(30):
    S = tuple(sorted(rng.sample(range(16), 7)))
    x = sum(1 << i for i in S)
    v = support_value_dummy(S, qd16, B16, 1, 16)
    cnt += 1
    pc7 += 1
    if v != Q16[x]:
        bad += 1
        print(f"  MISS x={x} S={S} (pc7)")
print(f"  {cnt} checks ({pc7} at pc=7): {bad} misses")

# ---- 4) decision-nondegeneracy of the dummy game (direct solver, m=4)
print("\ndecision-nondegeneracy, ECHO-FIFO+dummy, (4,1,7) t=0 (newly exact):")
Q, qd, B = make_form(4, 1, 7)
qd5 = qd + [0]
B5 = [row + [0] for row in B] + [[0] * 5]


def analyze_fifo(x5, m5, qd, B, t):
    bits = [i for i in range(m5) if (x5 >> i) & 1]
    memo = {}

    def moves_of(u, openseq, last):
        out = []
        omask = 0
        for c in openseq:
            omask |= 1 << c
        for i in bits:
            if i == last:
                continue
            if (u >> i) & 1:
                ch = 0
                for kk in bits:
                    if kk > i and (omask >> kk) & 1 and B[kk][i]:
                        ch ^= 1
                out.append((i, ch, u ^ (1 << i), openseq + (i,)))
        if openseq:
            c = openseq[0]
            if c != last:
                ch = qd[c]
                for kk in bits:
                    if kk > c and (omask >> kk) & 1 and B[kk][c]:
                        ch ^= 1
                out.append((c, ch, u, openseq[1:]))
        return out

    def val(u, openseq, last, mover, sigma):
        if u == 0 and not openseq:
            return sigma
        key = (u, openseq, last, mover, sigma)
        r = memo.get(key)
        if r is not None:
            return r
        legal = moves_of(u, openseq, last)
        if not legal:
            res = val(u, openseq, -1, 1 - mover, sigma)
        else:
            want = t if mover == 0 else 1 - t
            outs = [val(u2, s2, i, 1 - mover, sigma ^ ch) for (i, ch, u2, s2) in legal]
            res = want if want in outs else 1 - want
        memo[key] = res
        return res

    root = val(x5, (), -1, 0, 0)
    seen = {(x5, (), -1, 0, 0)}
    stack = [(x5, (), -1, 0, 0)]
    nonterm = dec = mist = 0
    wit = None
    while stack:
        st = stack.pop()
        (u, openseq, last, mover, sigma) = st
        if u == 0 and not openseq:
            continue
        nonterm += 1
        legal = moves_of(u, openseq, last)
        if not legal:
            nxt = (u, openseq, -1, 1 - mover, sigma)
            if nxt not in seen:
                seen.add(nxt)
                stack.append(nxt)
            continue
        kids = []
        for (i, ch, u2, s2) in legal:
            s = (u2, s2, i, 1 - mover, sigma ^ ch)
            kids.append((i, val(*s)))
            if s not in seen:
                seen.add(s)
                stack.append(s)
        if len(legal) >= 2:
            dec += 1
            if len({v for (_, v) in kids}) > 1:
                mist += 1
                if wit is None:
                    wit = (st, kids)
    return root, nonterm, dec, mist, wit


tot_dec = tot_mist = 0
for x in range(16):
    root, nt, dec, mist, wit = analyze_fifo(x | 16, 5, qd5, B5, 0)
    assert root == Q[x], (x, root, Q[x])
    tot_dec += dec
    tot_mist += mist
print(f"  all 16 outcomes == Q  |  decision states {tot_dec}, mistake states "
      f"{tot_mist} -> {'NON-DEGENERATE' if tot_mist else 'CLOCK'}")

# ---- 5) explicit position graph: x = e0+e1 (pc 2) + dummy, lam=7
print("\nEXPLICIT GRAPH: (4,1,7), x = 3 = e0+e1, board {0,1,d}, t=0:")
print(f"  q0,q1 = {qd[0]},{qd[1]}; B01 = {B[0][1]}; Q(3) = {Q[3]}")
root, nt, dec, mist, wit = analyze_fifo(3 | 16, 5, qd5, B5, 0)
print(f"  value={root} (=Q: {root==Q[3]}), {nt} nonterminal states, "
      f"{dec} decision states, {mist} mistake states")
if wit:
    (st, kids) = wit
    print(f"  sample mistake state: u={st[0]:05b} open={st[1]} ko={st[2]} "
          f"mover=P{st[3]+1} sigma={st[4]} -> options {kids}")
x7 = 7 | 16
root, nt, dec, mist, wit = analyze_fifo(x7, 5, qd5, B5, 0)
print(f"\n  x = 7 = e0+e1+e2 board {{0,1,2,d}}: value={root} Q={Q[7]}, "
      f"{nt} states, {dec} decision, {mist} mistake")
if wit:
    (st, kids) = wit
    print(f"  sample mistake state: u={st[0]:05b} open={st[1]} ko={st[2]} "
          f"mover=P{st[3]+1} sigma={st[4]} -> options {kids}")
