"""Decision-degeneracy analysis + miss dissection + m=8 lambda landscape."""
import sys, time
sys.path.insert(0, "/tmp")
from asym2_probe import make_form, make_charge, solve_position

def analyze(x, m, charge, p1_target):
    """Full (unpruned) evaluation of every reachable state.
    Returns root value, #reachable nonterminal states, #decision states
    (>=2 legal moves), #mistake states (children values differ),
    and one (state, moves->values) mistake witness."""
    bits = [i for i in range(m) if (x >> i) & 1]
    if not bits:
        return 0, 0, 0, 0, None
    memo = {}

    def legal_moves(u, o, last):
        out = []
        for i in bits:
            if i == last:
                continue
            bit = 1 << i
            if u & bit:
                out.append((i, u ^ bit, o ^ bit))
            elif o & bit:
                out.append((i, u, o ^ bit))
        return out

    def val(u, o, last, mover, sigma):
        if u == 0 and o == 0:
            return sigma
        key = (u, o, last, mover, sigma)
        r = memo.get(key)
        if r is not None:
            return r
        legal = legal_moves(u, o, last)
        if not legal:
            res = val(u, o, -1, 1 - mover, sigma)
        else:
            want = p1_target if mover == 0 else 1 - p1_target
            outs = [val(u2, o2, i, 1 - mover, sigma ^ charge(o, i))
                    for (i, u2, o2) in legal]
            res = want if want in outs else 1 - want
        memo[key] = res
        return res

    root = val(x, 0, -1, 0, 0)
    # BFS reachable, count decision/mistake states
    seen = {(x, 0, -1, 0, 0)}
    stack = [(x, 0, -1, 0, 0)]
    nonterm = dec = mist = 0
    witness = None
    while stack:
        (u, o, last, mover, sigma) = stack.pop()
        if u == 0 and o == 0:
            continue
        nonterm += 1
        legal = legal_moves(u, o, last)
        if not legal:
            nxt = (u, o, -1, 1 - mover, sigma)
            if nxt not in seen:
                seen.add(nxt)
                stack.append(nxt)
            continue
        kids = []
        for (i, u2, o2) in legal:
            s2 = (u2, o2, i, 1 - mover, sigma ^ charge(o, i))
            kids.append((i, val(*s2)))
            if s2 not in seen:
                seen.add(s2)
                stack.append(s2)
        if len(legal) >= 2:
            dec += 1
            vs = {v for (_, v) in kids}
            if len(vs) > 1:
                mist += 1
                if witness is None:
                    witness = ((u, o, last, mover, sigma), kids)
    return root, nonterm, dec, mist, witness


print("=" * 72)
print("DECISION-DEGENERACY at the m=4 EXACT hits (orientation A, lower side)")
print("round-1 rigidity question: is every exact realizer a clock?")
print("=" * 72)
for lam in [1, 2, 12, 13, 14]:
    Q, qd, B = make_form(4, 1, lam)
    ch = make_charge(qd, B, 4, True)
    tot_dec = tot_mist = 0
    rows = []
    for x in range(16):
        root, nt, dec, mist, wit = analyze(x, 4, ch, 1)
        assert root == Q[x], (lam, x)
        tot_dec += dec
        tot_mist += mist
        if mist and len(rows) < 3:
            rows.append((x, nt, dec, mist, wit))
    print(f"lam={lam:2d}: total decision states {tot_dec}, "
          f"MISTAKE states {tot_mist} -> {'NON-DEGENERATE' if tot_mist else 'CLOCK'}")
    for (x, nt, dec, mist, wit) in rows:
        (st, kids) = wit
        print(f"   x={x} (pc{bin(x).count('1')}): {nt} states, {dec} decision, "
              f"{mist} mistake; e.g. state(u={st[0]},o={st[1]},ko={st[2]},"
              f"mover=P{st[3]+1},sig={st[4]}) options {kids}")

print()
print("=" * 72)
print("MISS DISSECTION: (8,2,1) lower/A, x=224 = e5+e6+e7")
print("=" * 72)
Q, qd, B = make_form(8, 2, 1)
x = 224
S = [5, 6, 7]
l = qd[5] ^ qd[6] ^ qd[7]
print(f"q5,q6,q7 = {qd[5]},{qd[6]},{qd[7]}  l_diag = {l}")
print(f"B56,B57,B67 = {B[5][6]},{B[5][7]},{B[6][7]}")
print(f"Q(224) = l + B56+B57+B67 = {Q[224]}")
ch = make_charge(qd, B, 8, True)
root, nt, dec, mist, wit = analyze(224, 8, ch, 1)
print(f"value = {root} (P1 forces sigma=1 although Q=0): {nt} states, "
      f"{dec} decision, {mist} mistake states")

# abstract 3-coin linking game: which B-patterns are exact under w=1 ko?
print()
print("abstract popcount-3 classification: forced linked-sum vs all-linked sum")
print("(b56,b57,b67 as (b01,b02,b12) on 3 coins; q=0; both orientations)")
def linking_game_value(bpat, p1_target):
    # 3 coins {0,1,2}, B[i][j] = bpat[(i,j)], q=0
    Bm = [[0] * 3 for _ in range(3)]
    Bm[0][1] = Bm[1][0] = bpat[0]
    Bm[0][2] = Bm[2][0] = bpat[1]
    Bm[1][2] = Bm[2][1] = bpat[2]
    qd3 = [0, 0, 0]
    ch3 = make_charge(qd3, Bm, 3, True)
    return solve_position(7, 3, ch3, p1_target)

for bpat in [(a, b, c) for a in (0, 1) for b in (0, 1) for c in (0, 1)]:
    allk = bpat[0] ^ bpat[1] ^ bpat[2]
    vA = linking_game_value(bpat, 1)
    vB = linking_game_value(bpat, 0)
    flag = ""
    if vA != allk:
        flag += "  <-- A-orientation NOT all-linked"
    if vB != allk:
        flag += "  <-- B-orientation NOT all-linked"
    print(f"b={bpat} all-linked={allk} forcedA={vA} forcedB={vB}{flag}")

print()
print("=" * 72)
print("m=8,a=1 lambda landscape (lower/A): agreement with Q over all 255 lambda")
print("=" * 72)
t0 = time.time()
hist = {}
best = []
for lam in range(1, 256):
    Q, qd, B = make_form(8, 1, lam)
    ch = make_charge(qd, B, 8, True)
    agree = sum(1 for x in range(256)
                if solve_position(x, 8, ch, 1) == Q[x])
    hist[agree] = hist.get(agree, 0) + 1
    if agree >= 250:
        best.append((lam, agree))
print(f"[{time.time()-t0:.0f}s] agreement histogram: {dict(sorted(hist.items()))}")
print(f"lambda with agreement >= 250: {best}")
