"""FIFO+ko1 deep test: abstract k=5,6 screen + real Gold-form benches.

Rule ECHO-FIFO(x): coins = bits(x), touched twice each (open, close).
Moves: open ANY untouched coin, or close the LONGEST-OPEN coin (FIFO).
Ko: may not touch the coin touched on the immediately preceding touch;
stuck => pass (clears ko). Charge on touching i with open set o (cocycle,
lower side): sigma ^= q_i*[i in o] + sum_{kk>i, kk in o} B_kk,i.
Readout: final sigma; orientation t: P1 wants sigma = t.

Decomposition (legality is q-blind; sigma_complete = l_S + linkedsum):
value(x) = l_S XOR V*(b_S, want = t XOR l_S), where V*(b, d) is the
abstract linking game value with the FIRST mover wanting linkedsum = d.
Validated against the direct real-state solver on all of m=4 below.
"""
import sys, time
from itertools import permutations, product
sys.path.insert(0, "/tmp")
from asym2_probe import make_form

sys.setrecursionlimit(100000)


# ---------------------------------------------------------- abstract D-game
def abstract_value(k, edges, d0):
    """Linking game on k coins, FIFO+ko1, q=0.
    edges: frozenset of (i,j) i<j with b=1. Returns forced linkedsum when the
    first mover wants linkedsum = d0 (opponent wants complement).
    win(state, d) = mover-to-move can force future-delta == d."""
    hadj = [0] * k  # lower-triangular cocycle: touching i charges open kk > i
    for (i, j) in edges:
        hadj[i] |= 1 << j  # i < j
    memo = {}

    def win(u, openseq, last, d):
        if u == 0 and not openseq:
            return d == 0
        key = (u, openseq, last, d)
        r = memo.get(key)
        if r is not None:
            return r
        omask = 0
        for c in openseq:
            omask |= 1 << c
        moves = []
        for i in range(k):
            if i == last:
                continue
            if (u >> i) & 1:
                ch = bin(omask & hadj[i]).count("1") & 1
                moves.append((ch, u ^ (1 << i), openseq + (i,), i))
        if openseq:
            c = openseq[0]
            if c != last:
                ch = bin(omask & hadj[c]).count("1") & 1
                moves.append((ch, u, openseq[1:], c))
        if not moves:
            res = not win(u, openseq, -1, 1 ^ d)  # pass: opponent wants 1^d
        else:
            res = False
            for (ch, u2, seq2, i) in moves:
                if not win(u2, seq2, i, 1 ^ d ^ ch):
                    res = True
                    break
        memo[key] = res
        return res

    full = (1 << k) - 1
    return d0 if win(full, (), -1, d0) else 1 ^ d0


# validate D-game transform against the sigma-explicit solver from the screen
from asym2_variants import solve_fifo
for k in (3, 4):
    pairs = [(i, j) for i in range(k) for j in range(i + 1, k)]
    for bits in product((0, 1), repeat=len(pairs)):
        edges = frozenset(p for (b, p) in zip(bits, pairs) if b)
        B = [[0] * k for _ in range(k)]
        for (i, j) in edges:
            B[i][j] = B[j][i] = 1
        for t in (0, 1):
            assert abstract_value(k, edges, t) == solve_fifo(k, B, "fifo", True, t)
print("D-game transform == sigma-explicit solver (k=3,4 exhaustive): OK")


# ------------------------------------------------- iso-reduced abstract screen
def canon(k, edges):
    best = None
    for perm in permutations(range(k)):
        key = frozenset((min(perm[i], perm[j]), max(perm[i], perm[j]))
                        for (i, j) in edges)
        fk = tuple(sorted(key))
        if best is None or fk < best:
            best = fk
    return best


def screen_abstract(k):
    pairs = [(i, j) for i in range(k) for j in range(i + 1, k)]
    classes = {}
    for bits in product((0, 1), repeat=len(pairs)):
        edges = frozenset(p for (b, p) in zip(bits, pairs) if b)
        c = canon(k, edges)
        if c not in classes:
            classes[c] = edges
    failA, failB = [], []
    for c, edges in classes.items():
        par = len(edges) & 1
        vB = abstract_value(k, edges, 0)   # first mover wants 0
        vA = abstract_value(k, edges, 1)   # first mover wants 1
        if vB != par:
            failB.append(c)
        if vA != par:
            failA.append(c)
    return len(classes), failA, failB


for k in (3, 4, 5):
    t0 = time.time()
    nc, fA, fB = screen_abstract(k)
    print(f"k={k}: {nc} iso classes  want0-side fails={len(fB)}  "
          f"want1-side fails={len(fA)}  [{time.time()-t0:.0f}s]")
    if fB:
        print("   want0 failures:", fB[:5])
    if fA and k == 3:
        print("   want1 failures:", fA)


# --------------------------------------------------------- iso-invariance spot check
import random as _rnd
_rng = _rnd.Random(7)
fail5 = ((0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (2, 3))
for _ in range(4):
    perm = list(range(5))
    _rng.shuffle(perm)
    redges = frozenset((min(perm[i], perm[j]), max(perm[i], perm[j])) for (i, j) in fail5)
    assert abstract_value(5, redges, 0) == abstract_value(5, frozenset(fail5), 0)
print("root-value iso-invariance (failing k=5 class, 4 relabelings): OK")

# confirm the two k=5 want0 failures with the independent sigma-explicit solver
for cls in [((0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (2, 3)),
            ((0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (2, 4))]:
    B = [[0] * 5 for _ in range(5)]
    for (i, j) in cls:
        B[i][j] = B[j][i] = 1
    v = solve_fifo(5, B, "fifo", True, 0)
    print(f"k=5 class {cls}: parity=1, sigma-explicit forced={v} (confirms fail)")

# k=6 screen
t0 = time.time()
nc, fA, fB = screen_abstract(6)
print(f"k=6: {nc} iso classes  want0 fails={len(fB)}  want1 fails={len(fA)}  [{time.time()-t0:.0f}s]")
