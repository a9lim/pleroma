"""The abstract linking game: reductions, screens, and the verified strategy.

The general-m linking-theorem chase (2026-06-10) for the echo-fifo+dummy
realizer (writeups/goldarf.tex SS8.3): abstract the verified sigma-valued
FIFO+ko+pass+dummy rule away from Gold forms to arbitrary support graphs,
and reduce the m-uniform exactness claim to one combinatorial statement.

THE REDUCED GAME.  Board: a finite graph on "coins"; state (U, q, ko) with
U = untouched coins, q = FIFO queue of open coins, ko = last-touched coin.
Moves: OPEN any x in U (x != ko, push to back) or CLOSE the queue front f
(f != ko, pop); no legal move => forced pass (clears ko).  A close FLIPS a
bit iff deg_U(f) is odd at that moment.  One player wants total flips even,
the other odd.

Reduction lemmas (each a short whole-play identity, machine-validated here):
  R1. FIFO => coins close in opening order => no chord nesting; a graph
      edge is LINKED iff the two open-windows overlap.  sigma == overlap
      parity of the played interval graph restricted to E(G).
  R2. (overlap accounting) D := sigma ^ |undecided edges| is invariant
      under opens and passes and flips exactly on odd-deg_U(front) closes.
      So the sigma-game IS the odd-close parity game above, with
      sigma_forced = |E| ^ (forced flip parity).
  R3. Opens are never ko-blocked (ko is always a touched coin); the ko
      blocks a close only when the front was just opened onto an empty
      queue; forced passes occur only once U = 0, after which deg_U == 0
      and no flips are possible.  Passes are irrelevant to the flip fight.

THE LINKING THEOREM (target).  If the board contains an isolated coin (the
dummy), the flip count is forced even -- both seats, every graph.  Hence
sigma is forced = |E| mod 2, which on a Gold board is Q(x): m-uniform
exactness of the echo-fifo+dummy realizer.

STATUS (2026-06-10), machine-verified by this file:
  * Rigidity holds for ALL graph iso classes with k <= 7 real coins +
    dummy, both seats (k=7: 1044 classes) -- far beyond the Gold-arising
    boards of the m=8 sweep.
  * Without the dummy the failures ("Bad graphs") are exactly mover-
    controlled, census {3:1, 5:4, 7:34}; none contains an isolated vertex;
    33/34 at n=7 have a dominating vertex (one composite exception).
  * Core Lemma (the unique local obstruction; proof = 4-case check): with
    the queue empty, after opening v with R = U \\ {v}, the responder can
    re-even v before it becomes closable UNLESS R is a subset of N(v) with
    |R| even -- the "domination device" (ko-protected zugzwang, flip in 2
    plies).  An isolated coin in U defeats it at every root.  |R| odd
    explains the bonus even-n no-dummy rigidity.
  * A two-mode defender strategy (PREVENTION/DEBT menus, rule_R3/debt_D3
    below) beats an optimal unrestricted attacker on every class k <= 7,
    both seats, with NO fallback outside the menus.  NB: menu-EXISTENTIAL
    semantics -- the menus always contain a winning move; not every menu
    choice wins (Codex exhibited a losing poison choice on the star).
  * General-n proof: OPEN.  Architecture (after a Codex spar, thread
    019eb4ff-695b-7762-97e8-c0bea66c4e7e): segment the queue at firewall
    coins (deg_U == 0; the opened dummy is permanent, the untouched dummy
    virtual), mutual induction E (no debt) / O (one debt) per segment;
    certificates bounded by game-tree depth.  The hard obligation is the
    poison transition E -> O (recursive repair-potential), which is also
    exactly where parity-local invariants provably fail (the safe/unsafe
    label is NOT a function of 13 natural parity features; minimal
    distinguishing pairs differ in E(U) repair structure).

Stages: validate | screen [K] | strategy [K] | all   (default K = 5; the
k=7 screen ~45 s, k=7 strategy ~25 s).  Stdlib only, no venv needed.
Cross-validated against experiments/echo_solver.py (the adversarially
reviewed solver) through the SynthForm bridge in stage `validate`.
"""

import random
import sys
import time
from itertools import permutations

sys.setrecursionlimit(1000000)


# ---------------------------------------------------------------- solvers

def adj_of(n: int, edges) -> list:
    adj = [0] * n
    for (i, j) in edges:
        adj[i] |= 1 << j
        adj[j] |= 1 << i
    return adj


def legal_moves(n, adj, u, seq, last):
    """All legal moves as (kind, coin, flip, u2, seq2, touched)."""
    mv = []
    for i in range(n):
        if i != last and u >> i & 1:
            mv.append(("o", i, 0, u ^ (1 << i), seq + (i,), i))
    if seq and seq[0] != last:
        f = seq[0]
        fl = bin(adj[f] & u).count("1") & 1
        mv.append(("c", f, fl, u, seq[1:], f))
    return mv


def rigid_values(k: int, edges, dummy: bool) -> list:
    """[value(P1 wants 0), value(P1 wants 1)] of the flip-parity game,
    d-folded full-state solver (the same move semantics as
    echo_solver.fifo_value; charge convention differs only by bookkeeping
    -- totals agree, validated in stage `validate`)."""
    n = k + (1 if dummy else 0)
    adj = adj_of(n, edges)
    memo: dict = {}

    def win(u, seq, last, g):
        # mover can force future flip count == g (mod 2)
        if u == 0 and not seq:
            return g == 0
        key = (u, seq, last, g)
        r = memo.get(key)
        if r is not None:
            return r
        mv = legal_moves(n, adj, u, seq, last)
        if not mv:
            res = not win(u, seq, -1, 1 ^ g)  # forced pass clears ko
        else:
            res = any(not win(u2, s2, i, 1 ^ g ^ fl)
                      for (_t, _c, fl, u2, s2, i) in mv)
        memo[key] = res
        return res

    par = len(edges) & 1
    full = (1 << n) - 1
    out = []
    for t in (0, 1):
        g = t ^ par  # flips needed for sigma == t
        out.append(t if win(full, (), -1, g) else 1 ^ t)
    return out


def sigma_value(k: int, edges, dummy: bool, t: int) -> int:
    """Sigma-explicit oracle for the d-folding (lower-cocycle charges,
    byte-for-byte the echo_solver.fifo_value recursion shape with q = 0)."""
    n = k + (1 if dummy else 0)
    adj = adj_of(n, edges)
    hadj = [adj[i] & ~((2 << i) - 1) for i in range(n)]
    memo: dict = {}

    def rec(u, seq, last, mover, sigma):
        if u == 0 and not seq:
            return sigma
        key = (u, seq, last, mover, sigma)
        r = memo.get(key)
        if r is not None:
            return r
        omask = 0
        for c in seq:
            omask |= 1 << c
        legal = []
        for i in range(n):
            if i != last and u >> i & 1:
                ch = bin(omask & hadj[i]).count("1") & 1
                legal.append((i, ch, u ^ (1 << i), seq + (i,)))
        if seq and seq[0] != last:
            c = seq[0]
            ch = bin(omask & hadj[c]).count("1") & 1
            legal.append((c, ch, u, seq[1:]))
        if not legal:
            res = rec(u, seq, -1, 1 - mover, sigma)
        else:
            want = t if mover == 0 else 1 - t
            res = 1 - want
            for (i, ch, u2, s2) in legal:
                if rec(u2, s2, i, 1 - mover, sigma ^ ch) == want:
                    res = want
                    break
        memo[key] = res
        return res

    return rec((1 << n) - 1, (), -1, 0, 0)


# ---------------------------------------------------------------- iso classes

def iso_classes(k: int) -> list:
    """One representative edge-frozenset per isomorphism class on k labelled
    vertices (orbit marking; fine through k = 7)."""
    pairs = [(i, j) for i in range(k) for j in range(i + 1, k)]
    pidx = {p: ii for ii, p in enumerate(pairs)}
    perms = list(permutations(range(k)))
    seen = set()
    reps = []
    for gmask in range(1 << len(pairs)):
        if gmask in seen:
            continue
        edges = frozenset(p for ii, p in enumerate(pairs) if gmask >> ii & 1)
        reps.append(edges)
        for perm in perms:
            om = 0
            for (i, j) in edges:
                a, b = perm[i], perm[j]
                om |= 1 << pidx[(min(a, b), max(a, b))]
            seen.add(om)
    return reps


# ---------------------------------------------------------------- strategy

def rule_R3(n, adj, u, seq, last):
    """PREVENTION menu (debt 0).  P1 re-even / P2 safe opens + safe close /
    P3 poison-or-close trap branch / P4 endgame close."""
    front = seq[0] if seq else None
    allowed = set()
    if u == 0:
        if front is not None and front != last:
            allowed.add(("c", front))
        return allowed
    if front is not None and bin(adj[front] & u).count("1") & 1:
        for i in range(n):
            if i != last and (u >> i & 1) and (adj[front] >> i & 1):
                allowed.add(("o", i))
        return allowed
    nontog = {("o", i) for i in range(n)
              if i != last and (u >> i & 1)
              and (front is None or not adj[front] >> i & 1)}
    if nontog:
        allowed |= nontog
        if front is not None and front != last:
            nxt = seq[1] if len(seq) > 1 else None
            if nxt is None or bin(adj[nxt] & u).count("1") % 2 == 0:
                allowed.add(("c", front))
        return allowed
    for i in range(n):
        if i != last and (u >> i & 1):
            allowed.add(("o", i))
    if front is not None and front != last:
        allowed.add(("c", front))
    return allowed


def debt_D3(n, adj, u, seq, last):
    """DEBT menu (debt 1).  D1 counter-close / D2 ko stall / D3 toggle or
    advance / D4 bare opens."""
    front = seq[0] if seq else None
    allowed = set()
    if front is not None and bin(adj[front] & u).count("1") & 1:
        if front != last:
            return {("c", front)}
        for i in range(n):
            if i != last and u >> i & 1:
                allowed.add(("o", i))
        return allowed
    if front is not None:
        for i in range(n):
            if i != last and (u >> i & 1) and (adj[front] >> i & 1):
                allowed.add(("o", i))
        if front != last:
            allowed.add(("c", front))
        return allowed
    for i in range(n):
        if i != last and u >> i & 1:
            allowed.add(("o", i))
    return allowed


def strategy_holds(k: int, edges, seat: int) -> bool:
    """Defender (flips-even) restricted to the R3/D3 menus, attacker
    unrestricted optimal; STRICT (an empty/illegal menu = defender loss).
    Menu-existential: True means a winning move always exists IN the menu."""
    n = k + 1  # always with dummy
    adj = adj_of(n, edges)
    memo: dict = {}

    def W(u, seq, last, mover, g):
        if u == 0 and not seq:
            return g == 0
        key = (u, seq, last, mover, g)
        r = memo.get(key)
        if r is not None:
            return r
        lm = legal_moves(n, adj, u, seq, last)
        if mover == seat and lm:
            rule = rule_R3 if g == 0 else debt_D3
            allowed = rule(n, adj, u, seq, last)
            mv = [m for m in lm if (m[0], m[1]) in allowed]
            if not mv:
                memo[key] = False
                return False
        else:
            mv = lm
        if not mv:
            res = W(u, seq, -1, 1 - mover, g)
        elif mover == seat:
            res = any(W(u2, s2, i, 1 - mover, g ^ fl)
                      for (_t, _c, fl, u2, s2, i) in mv)
        else:
            res = all(W(u2, s2, i, 1 - mover, g ^ fl)
                      for (_t, _c, fl, u2, s2, i) in mv)
        memo[key] = res
        return res

    return W((1 << n) - 1, (), -1, 0, 0)


# ---------------------------------------------------------------- stages

def stage_validate() -> None:
    print("== d-folded vs sigma-explicit (k <= 4 exhaustive, dummy on/off) ==")
    cnt = 0
    for k in (2, 3, 4):
        pairs = [(i, j) for i in range(k) for j in range(i + 1, k)]
        for gmask in range(1 << len(pairs)):
            edges = frozenset(p for ii, p in enumerate(pairs)
                              if gmask >> ii & 1)
            for dummy in (True, False):
                vals = rigid_values(k, edges, dummy)
                for t in (0, 1):
                    assert vals[t] == sigma_value(k, edges, dummy, t), \
                        (k, edges, dummy, t)
                    cnt += 1
    print(f"   {cnt} agree")

    print("== sigma-explicit vs the verified echo_solver.fifo_value ==")
    import pathlib
    sys.path.insert(0, str(pathlib.Path(__file__).parent))
    from echo_solver import fifo_value, SynthForm
    rng = random.Random(2026)
    cnt = 0
    for _ in range(12):
        k = 5
        pairs = [(i, j) for i in range(k) for j in range(i + 1, k)]
        edges = frozenset(p for p in pairs if rng.random() < 0.5)
        B = [[0] * k for _ in range(k)]
        for (i, j) in edges:
            B[i][j] = B[j][i] = 1
        from typing import Any
        f: Any = SynthForm(k, [0] * k, B)  # duck-types Form for fifo_value
        for dummy in (True, False):
            for t in (0, 1):
                assert sigma_value(k, edges, dummy, t) == \
                    fifo_value(f, (1 << k) - 1, t, dummy=dummy), \
                    (edges, dummy, t)
                cnt += 1
    print(f"   {cnt} agree (SynthForm bridge, q = 0)")

    print("== reduction identities on random plays (R1/R2) ==")
    rng = random.Random(7)
    for _ in range(400):
        k = rng.randrange(2, 7)
        pairs = [(i, j) for i in range(k) for j in range(i + 1, k)]
        edges = frozenset(p for p in pairs if rng.random() < 0.5)
        adj = adj_of(k, edges)
        hadj = [adj[i] & ~((2 << i) - 1) for i in range(k)]
        u, sigma, tt = (1 << k) - 1, 0, 0
        seq: tuple = ()
        windows = {}
        flips = 0
        while u or seq:
            omask = 0
            for c in seq:
                omask |= 1 << c
            opens = [i for i in range(k) if u >> i & 1]
            if opens and (not seq or rng.random() < 0.6):
                i = rng.choice(opens)
                sigma ^= bin(omask & hadj[i]).count("1") & 1
                u ^= 1 << i
                seq = seq + (i,)
                windows[i] = [tt, None]
            else:
                c = seq[0]
                sigma ^= bin(omask & hadj[c]).count("1") & 1
                flips ^= bin(u & adj[c]).count("1") & 1
                seq = seq[1:]
                windows[c][1] = tt
            tt += 1
        overlap = 0
        for (i, j) in edges:
            (a1, b1), (a2, b2) = windows[i], windows[j]
            assert not (a1 < a2 and b2 < b1) and not (a2 < a1 and b1 < b2), \
                "FIFO nesting impossible (R1)"
            if a1 < a2 < b1 or a2 < a1 < b2:
                overlap ^= 1
        assert sigma == overlap, "sigma != overlap parity (R1)"
        assert flips == (len(edges) & 1) ^ sigma, "flips != |E| ^ sigma (R2)"
    print("   400 random plays: no nesting; sigma == overlap;"
          " odd-close flips == |E| ^ sigma")
    print("validate: PASS")


def stage_screen(kmax: int) -> None:
    for k in range(2, kmax + 1):
        t0 = time.time()
        reps = iso_classes(k)
        fails_d, bad = [], []
        for edges in reps:
            par = len(edges) & 1
            if rigid_values(k, edges, True) != [par, par]:
                fails_d.append(tuple(sorted(edges)))
            vn = rigid_values(k, edges, False)
            if vn != [par, par]:
                bad.append((tuple(sorted(edges)), vn))
        niso = sum(1 for (e, _v) in bad
                   if min(sum(1 for (a, b) in e if v in (a, b))
                          for v in range(k)) == 0) if bad else 0
        print(f"k={k}: {len(reps)} classes | WITH dummy fails: {len(fails_d)}"
              f" | no-dummy Bad: {len(bad)} (with isolated vertex: {niso},"
              f" all mover-controlled:"
              f" {all(v == [0, 1] for (_e, v) in bad)})"
              f"  [{time.time()-t0:.0f}s]", flush=True)
        for e in fails_d:
            print(f"   THEOREM COUNTEREXAMPLE {e}")


def stage_strategy(kmax: int) -> None:
    for k in range(2, kmax + 1):
        t0 = time.time()
        reps = iso_classes(k)
        fails = [(tuple(sorted(e)), seat)
                 for e in reps for seat in (0, 1)
                 if not strategy_holds(k, e, seat)]
        print(f"k={k}: {len(reps)} classes x 2 seats | R3/D3 strict fails:"
              f" {len(fails)}  [{time.time()-t0:.0f}s]", flush=True)
        for f in fails[:8]:
            print(f"   FAIL {f}")


def main() -> None:
    stage = sys.argv[1] if len(sys.argv) > 1 else "all"
    kmax = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    if stage in ("validate", "all"):
        stage_validate()
    if stage in ("screen", "all"):
        stage_screen(kmax)
    if stage in ("strategy", "all"):
        stage_strategy(kmax)


if __name__ == "__main__":
    main()
