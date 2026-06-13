#!/usr/bin/env python3
"""Fresh direct stateful solver for the echo rule family (adversarial review).

VERDICT (2026-06-10): CONFIRM.  The echo-fifo+dummy m=8 exactness claim is
re-derived in full -- 391,680/391,680 checks (765 scaled Gold forms x 256
positions x both stances), zero misses -- by the direct full-state solver
`fifo_value` below (stage `fifo2-all`), with no decomposition and no
isomorphism caching.  Record and corrected rule description / verification
record: `writeups/goldarf.tex` SS8.  NOTE: the EchoGame
class and the `pin-ko`/`fifo-m4`/`fifo-m8` stages implement the OLD SS8.3
prose readings, kept as the documented negative result (no reading of that
prose is m=8-exact -- the prose misdescribed the rule); the faithful
sigma-valued rule lives in `fifo_value`/`ko_value` and the `fifo2-*`/`ko2`
stages.

This is the decisive-experiment harness for the `echo-solver` task (successor
`echo-family-sweep` in `docs/COMPLETENESS.md`) and
`writeups/goldarf.tex` SS8-9, ranked move 1: an independent re-derivation of the
echo-fifo+dummy m=8 exactness claim (391,680 checks / 765 scaled Gold forms),
which was produced by a decomposition-plus-isomorphism-caching solver validated
only at m=4 and was recorded in the writeup as unverified.

Design constraints (from the pre-registered experiment):
  * direct stateful solver -- NO decomposition, NO isomorphism caching;
  * FULL state in the memoization key, including the accumulated charge sigma;
  * validated against explicit (unmemoized) tree enumeration before being
    trusted at m=8.

Clean-room provenance: the rules are implemented from the prose of
goldarf.tex SS8 alone; the original probes in experiments/gold/ were NOT read
before this harness produced numbers. Because the prose underdetermines a few
conventions, each is an explicit parameter here, and the conventions are pinned
by reproducing the *corrected* echo-ko results table of SS8.2 (16/16 at m=4 bent
lambda in {2,12}; 255/256 at (8,2,1) with the unique miss x=224; 228/256 at
(8,1,1); 212/256 at (8,1,2)).

Rule (echo-ko, goldarf SS8.2): positions are x in F_{2^m}; the coins of x (the
support bits) must each be touched twice (first touch opens, second closes);
players alternate single touches; each touch of e_i updates a charge
sigma ^= c(open-set, e_i), where c is the triangular cocycle of the
extraspecial extension (Lemma extdict(ii)):

    c(u, v) = sum_i u_i v_i q_i  +  sum_{i>j} u_i v_j B(e_i, e_j)

with q_i = Q(e_i) and B the polar form; a one-move ko forbids immediately
re-touching; the player completing the board wins iff the terminal charge is
sigma = 1. A player with no legal move loses (normal-play stall).

Rule (echo-fifo+dummy, goldarf SS8.3): the one-move ko is replaced by a FIFO
close-discipline (only the longest-open coin may be closed; opens are free) and
one neutral tempo coin (q = 0, B-row = 0) is adjoined.

Claim levels: the nim arithmetic and Gold forms are standard math, cross-checked
here against the independent mex/Turning-Corners recurrence; everything this
script prints about game outcomes is "implemented and tested" for THIS harness
and these conventions only.

Usage:
    python3 experiments/echo_solver.py selftest    # nim oracle + tree-vs-memo
    python3 experiments/echo_solver.py pin-ko      # convention sweep vs SS8.2 table
    python3 experiments/echo_solver.py fifo-m4     # full m=4 family, all conventions
    python3 experiments/echo_solver.py fifo-m8     # m=8 benchmarks (surviving convs)
    python3 experiments/echo_solver.py stratified  # >=20 stratified lambda at m=8
"""

from __future__ import annotations

import argparse
import sys
import time
from functools import lru_cache

sys.setrecursionlimit(10_000)

# ----------------------------------------------------------------------------
# Nim arithmetic (standard math; independent of src/scalar/)
# ----------------------------------------------------------------------------


@lru_cache(maxsize=None)
def nim_mul(a: int, b: int) -> int:
    """Nim multiplication on naturals (Conway), Fermat-power recursion."""
    if a > b:
        a, b = b, a
    if a == 0:
        return 0
    if a == 1:
        return b
    f = 2
    while f * f <= b:
        f *= f
    bh, bl = b // f, b % f
    if a < f:
        return nim_mul(a, bh) * f ^ nim_mul(a, bl)
    ah, al = a // f, a % f
    t = nim_mul(ah, bh)
    high = (t ^ nim_mul(ah, bl) ^ nim_mul(al, bh)) * f
    low = nim_mul(al, bl) ^ nim_mul(t, f // 2)
    return high ^ low


def nim_pow_2exp(x: int, a: int) -> int:
    """x^(2^a) by repeated nim squaring."""
    for _ in range(a):
        x = nim_mul(x, x)
    return x


def nim_trace(x: int, m: int) -> int:
    """Absolute trace F_{2^m} -> F_2 (XOR of Frobenius orbit)."""
    t, y = 0, x
    for _ in range(m):
        t ^= y
        y = nim_mul(y, y)
    assert t in (0, 1), f"trace escaped F_2: Tr({x}) = {t} (m={m})"
    return t


def mex_product_table(n: int) -> list[list[int]]:
    """Turning-Corners mex recurrence, the independent game-theoretic oracle."""
    tc = [[0] * n for _ in range(n)]
    for x in range(n):
        for y in range(n):
            seen = {tc[i][y] ^ tc[x][j] ^ tc[i][j] for i in range(x) for j in range(y)}
            v = 0
            while v in seen:
                v += 1
            tc[x][y] = v
    return tc


# ----------------------------------------------------------------------------
# Scaled Gold forms Q_{lambda,a}(x) = Tr(lambda * x^(1+2^a)) and their data
# ----------------------------------------------------------------------------


class Form:
    """Q, diagonal q_i, polar form B, and the triangular-cocycle masks."""

    def __init__(self, m: int, a: int, lam: int):
        self.m, self.a, self.lam = m, a, lam
        n = 1 << m
        self.Q = [nim_trace(nim_mul(lam, nim_mul(x, nim_pow_2exp(x, a))), m) for x in range(n)]
        self.q = [self.Q[1 << i] for i in range(m)]
        self.B = [[self.Q[(1 << i) ^ (1 << j)] ^ self.q[i] ^ self.q[j] if i != j else 0
                   for j in range(m)] for i in range(m)]
        # cocycle B-masks: orientation 'lower' uses c(v,e_i) += v_j B[j][i] for j>i,
        # 'upper' for j<i.  (q_i term: += v_i q_i in both.)
        self.mask_lower = [sum(1 << j for j in range(i + 1, m) if self.B[j][i]) for i in range(m)]
        self.mask_upper = [sum(1 << j for j in range(i) if self.B[j][i]) for i in range(m)]

    def rank(self) -> int:
        rows = [sum(self.B[i][j] << j for j in range(self.m)) for i in range(self.m)]
        r = 0
        for col in range(self.m):
            piv = next((k for k in range(r, self.m) if rows[k] >> col & 1), None)
            if piv is None:
                continue
            rows[r], rows[piv] = rows[piv], rows[r]
            for k in range(self.m):
                if k != r and rows[k] >> col & 1:
                    rows[k] ^= rows[r]
            r += 1
        return r

    def cocycle_inc(self, open_mask: int, i: int, orientation: str) -> int:
        """c(open_mask, e_i) -- the charge increment for touching coin i."""
        mask = self.mask_lower[i] if orientation == "lower" else self.mask_upper[i]
        inc = bin(open_mask & mask).count("1") & 1
        if open_mask >> i & 1:
            inc ^= self.q[i]
        return inc


# ----------------------------------------------------------------------------
# The games.  Coins are field-basis indices 0..m-1; the dummy coin is index m
# (charge-free: q = 0, B-row = 0, and it is excluded from real open-masks).
# ----------------------------------------------------------------------------

DUMMY = None  # set per-solve to the dummy index, or -1 when absent


class EchoGame:
    """Parameterized echo game on the support of x.

    discipline: 'ko'   -- close any open coin; touching the last-touched coin
                          is forbidden (one-move ko, w=1); stall loses.
                'fifo' -- only the longest-open coin may be closed; opens free.
    dummy:      'none' | 'required' (must be closed for completion)
                       | 'ignore'   (completion looks at real coins only)
                       | 'resolved' (real closed AND dummy not left open)
    timing:     'pre'  -- open-set BEFORE the touch feeds the cocycle
                'post' -- open-set AFTER the touch feeds the cocycle
    orientation:'lower' | 'upper' triangular cocycle.
    ko_mode:    'hard' -- the ko ban is absolute; a player whose only touchable
                          coin is ko-banned is stalled and loses (normal play).
                'soft' -- the ban yields when it is the ONLY touch available
                          (otherwise singletons can never complete).
    win_target: completer wins iff terminal sigma == win_target.  The SS8.2
                prose says sigma = 1; the corrected results table is only
                reproducible with win_target = 0 from a zero seed, i.e. the
                original run's charge is seeded at 1 (equivalent game).
    """

    def __init__(self, form: Form, x: int, discipline: str, dummy: str,
                 timing: str, orientation: str, ko_mode: str = "soft",
                 win_target: int = 0, dummy_touches: int = 2):
        self.form = form
        self.x = x
        self.discipline = discipline
        self.dummy = dummy
        self.timing = timing
        self.orientation = orientation
        self.ko_mode = ko_mode
        self.win_target = win_target
        self.dummy_touches = dummy_touches
        self.real = [i for i in range(form.m) if x >> i & 1]
        self.realmask = x
        self.dummy_idx = form.m if dummy != "none" else -1
        self.coins = self.real + ([self.dummy_idx] if dummy != "none" else [])

    # charge increment for touching coin i with open real-mask v (pre-touch)
    def _inc(self, v_pre: int, i: int, opening: bool) -> int:
        if i == self.dummy_idx:
            return 0
        if self.timing == "pre":
            v = v_pre
        else:  # post
            v = (v_pre | (1 << i)) if opening else (v_pre & ~(1 << i))
        return self.form.cocycle_inc(v, i, self.orientation)

    def _complete(self, untouched: frozenset, open_seq: tuple) -> bool:
        if self.dummy == "required" or self.dummy == "none":
            return not untouched and not open_seq
        real_open = any(c != self.dummy_idx for c in open_seq)
        real_untouched = any(c != self.dummy_idx for c in untouched)
        if self.dummy == "ignore":
            return not real_open and not real_untouched
        # 'resolved': real coins done and dummy not left hanging open
        return (not real_open and not real_untouched
                and self.dummy_idx not in open_seq)

    # --- move generation -----------------------------------------------------
    # A state is (untouched frozenset, open_seq tuple-in-open-order, last, sigma).
    # 'last' is the ko memory (last touched coin) for 'ko', and unused ('') for
    # 'fifo'.  sigma is the accumulated charge.  Side to move is implicit (no
    # passing in either rule), so it need not enter the key; everything else is
    # the FULL state, per the pre-registered spec.

    def moves(self, untouched: frozenset, open_seq: tuple, last):
        out = []
        if self.discipline == "ko":
            for c in untouched:
                if c != last:
                    out.append(("open", c))
            for c in open_seq:
                if c != last:
                    out.append(("close", c))
            if not out and self.ko_mode == "soft":
                # the ban yields when the ko-banned touch is the only one left
                if last in untouched:
                    out.append(("open", last))
                elif last in open_seq:
                    out.append(("close", last))
        else:  # fifo / lifo
            for c in untouched:
                if c == self.dummy_idx and self.dummy_touches == 1:
                    out.append(("token", c))
                else:
                    out.append(("open", c))
            if open_seq:
                out.append(("close", open_seq[0] if self.discipline == "fifo"
                            else open_seq[-1]))
        return out

    def apply(self, untouched, open_seq, sigma, move):
        kind, c = move
        if kind == "token":
            return untouched - {c}, open_seq, sigma
        v_pre = 0
        for o in open_seq:
            if o != self.dummy_idx:
                v_pre |= 1 << o
        if kind == "open":
            sigma ^= self._inc(v_pre, c, True)
            return untouched - {c}, open_seq + (c,), sigma
        sigma ^= self._inc(v_pre, c, False)
        return untouched, tuple(o for o in open_seq if o != c), sigma

    # --- solvers --------------------------------------------------------------

    def outcome(self) -> str:
        """'P' if the second player (previous player) wins, else 'N'."""
        untouched = frozenset(self.coins)
        if self._complete(untouched, ()):
            return "P"  # no game: mover never moves, second player wins
        memo: dict = {}

        def win(untouched, open_seq, last, sigma) -> bool:
            key = (untouched, open_seq, last, sigma)
            got = memo.get(key)
            if got is not None:
                return got
            result = False
            for mv in self.moves(untouched, open_seq, last):
                u2, o2, s2 = self.apply(untouched, open_seq, sigma, mv)
                if self._complete(u2, o2):
                    if s2 == self.win_target:  # completes the board on target: wins
                        result = True
                        break
                    continue  # completed off-target: mover loses this line
                if not win(u2, o2, mv[1] if self.discipline == "ko" else "", s2):
                    result = True
                    break
            memo[key] = result
            return result

        return "N" if win(untouched, (), "", 0) else "P"

    def outcome_tree(self) -> str:
        """Explicit tree enumeration -- NO memoization.  The validation oracle."""
        untouched = frozenset(self.coins)
        if self._complete(untouched, ()):
            return "P"

        def win(untouched, open_seq, last, sigma) -> bool:
            for mv in self.moves(untouched, open_seq, last):
                u2, o2, s2 = self.apply(untouched, open_seq, sigma, mv)
                if self._complete(u2, o2):
                    if s2 == self.win_target:
                        return True
                    continue
                if not win(u2, o2, mv[1] if self.discipline == "ko" else "", s2):
                    return True
            return False

        return "N" if win(untouched, (), "", 0) else "P"

    def liveness(self) -> tuple[int, int]:
        """(mistake_states, reachable_states): states whose mover has BOTH a
        winning and a losing option.  Full expansion -- no early exit."""
        untouched = frozenset(self.coins)
        if self._complete(untouched, ()):
            return 0, 0
        memo: dict = {}
        mistakes = 0

        def win(untouched, open_seq, last, sigma) -> bool:
            nonlocal mistakes
            key = (untouched, open_seq, last, sigma)
            got = memo.get(key)
            if got is not None:
                return got
            wins = loses = 0
            for mv in self.moves(untouched, open_seq, last):
                u2, o2, s2 = self.apply(untouched, open_seq, sigma, mv)
                if self._complete(u2, o2):
                    good = s2 == self.win_target
                else:
                    good = not win(u2, o2, mv[1] if self.discipline == "ko" else "", s2)
                if good:
                    wins += 1
                else:
                    loses += 1
            if wins and loses:
                mistakes += 1
            memo[key] = wins > 0
            return wins > 0

        win(untouched, (), "", 0)
        return mistakes, len(memo)


# ----------------------------------------------------------------------------
# Instances and conventions
# ----------------------------------------------------------------------------


def exactness(form: Form, conv: dict, positions=None, progress=False) -> tuple[int, list[int]]:
    """Agreement count of {outcome == P} vs {Q == 0} and the miss list."""
    n = 1 << form.m
    xs = range(n) if positions is None else positions
    misses = []
    for x in xs:
        g = EchoGame(form, x, **conv)
        if (g.outcome() == "P") != (form.Q[x] == 0):
            misses.append(x)
        if progress and x % 32 == 31:
            print(f"      ... {x + 1}/{n}, misses so far: {len(misses)}", flush=True)
    total = n if positions is None else len(list(xs))
    return total - len(misses), misses


KO_CONVS: list[dict] = [
    {"discipline": "ko", "dummy": "none", "timing": t, "orientation": o,
     "ko_mode": k, "win_target": w}
    for t in ("pre", "post") for o in ("lower", "upper")
    for k in ("soft", "hard") for w in (0, 1)
]

FIFO_CONVS: list[dict] = [
    {"discipline": disc, "dummy": d, "dummy_touches": k, "timing": t,
     "orientation": o, "win_target": w}
    for disc in ("fifo", "lifo")
    for (d, k) in (("required", 2), ("ignore", 2), ("resolved", 2),
                   ("required", 1), ("ignore", 1), ("none", 2))
    for t in ("pre", "post") for o in ("lower", "upper") for w in (0, 1)
]


def conv_name(conv: dict) -> str:
    bits = [conv["discipline"], conv["dummy"], conv["timing"], conv["orientation"]]
    if conv["discipline"] == "ko":
        bits.append(conv.get("ko_mode", "soft"))
    else:
        bits.append(f"x{conv.get('dummy_touches', 2)}")
    bits.append(f"tgt{conv.get('win_target', 0)}")
    return "/".join(map(str, bits))


# The corrected echo-ko results table (goldarf SS8.2), the convention-pinning
# fixture.  agreement is out of 2^m; the (8,2,1) row pins the miss set exactly.
KO_FIXTURES = [
    # (m, a, lam, expected_agreement, expected_miss_set_or_None)
    (4, 1, 2, 16, []),
    (4, 1, 12, 16, []),
    (8, 2, 1, 255, [224]),
    (8, 1, 1, 228, None),
    (8, 1, 2, 212, None),
]


# ----------------------------------------------------------------------------
# Stages
# ----------------------------------------------------------------------------


def stage_selftest() -> None:
    print("== nim arithmetic vs the Turning-Corners mex recurrence (0..31) ==")
    tc = mex_product_table(32)
    bad = [(x, y) for x in range(32) for y in range(32) if tc[x][y] != nim_mul(x, y)]
    assert not bad, f"nim_mul disagrees with mex recurrence at {bad[:5]}"
    print("   1024/1024 products agree (independent game-theoretic oracle)")

    print("== Gold rank spot checks (goldarf Table 1) ==")
    for m, a, lam, want in [(4, 1, 1, 2), (8, 1, 1, 6), (8, 2, 1, 4), (8, 1, 2, 8)]:
        f = Form(m, a, lam)
        got = f.rank()
        flag = "ok" if got == want else "MISMATCH"
        print(f"   ({m},{a},lam={lam}): rank {got} (expected {want}) {flag}")
        assert got == want

    print("== memoized solver vs explicit tree enumeration ==")
    # all m=4 positions, several lambdas, every convention shape
    checked = 0
    for lam in (1, 2, 12):
        f = Form(4, 1, lam)
        for conv in KO_CONVS + FIFO_CONVS:
            for x in range(16):
                g = EchoGame(f, x, **conv)
                a_, b_ = g.outcome(), g.outcome_tree()
                assert a_ == b_, f"memo/tree mismatch: m=4 lam={lam} x={x} {conv_name(conv)}: {a_} vs {b_}"
                checked += 1
    print(f"   m=4: {checked} (position, convention) solves agree with tree enumeration")
    # m=8 spot checks on small supports (tree enumeration is exponential)
    f8 = Form(8, 1, 1)
    small = [x for x in range(256) if bin(x).count("1") <= 3][:40]
    checked = 0
    for conv in (KO_CONVS[0], FIFO_CONVS[0], FIFO_CONVS[5], FIFO_CONVS[-1]):
        for x in small:
            g = EchoGame(f8, x, **conv)
            assert g.outcome() == g.outcome_tree(), f"m=8 mismatch x={x} {conv_name(conv)}"
            checked += 1
    print(f"   m=8 (|supp| <= 3): {checked} spot solves agree with tree enumeration")
    print("selftest: PASS")


def stage_pin_ko() -> None:
    print("== convention sweep vs the corrected echo-ko table (goldarf SS8.2) ==")
    for conv in KO_CONVS:
        print(f"-- {conv_name(conv)}")
        ok = True
        for m, a, lam, want, want_miss in KO_FIXTURES:
            t0 = time.time()
            f = Form(m, a, lam)
            agree, misses = exactness(f, conv)
            n = 1 << m
            status = "ok" if agree == want else "X"
            if want_miss is not None and sorted(misses) != sorted(want_miss):
                status = "X(miss-set)"
            ok &= status == "ok"
            print(f"   ({m},{a},lam={lam}): {agree}/{n}  expected {want}/{n}"
                  f"  misses={misses if len(misses) <= 6 else len(misses)}  [{status}]"
                  f"  ({time.time() - t0:.1f}s)")
        print(f"   => {'REPRODUCES the corrected table' if ok else 'does not match'}")


def stage_fifo_m4() -> list[dict]:
    print("== echo-fifo+dummy, full m=4 family (a=1, lambda in F_16^*), all conventions ==")
    forms = [Form(4, 1, lam) for lam in range(1, 16)]
    survivors = []
    for conv in FIFO_CONVS:
        total_miss = 0
        rows = []
        for f in forms:
            agree, misses = exactness(f, conv)
            total_miss += len(misses)
            rows.append((f.lam, agree, misses))
        tag = "EXACT on the full m=4 family" if total_miss == 0 else f"{total_miss} misses"
        print(f"-- {conv_name(conv)}: {tag}")
        if total_miss and total_miss <= 24:
            for lam, agree, misses in rows:
                if misses:
                    print(f"      lam={lam}: {agree}/16 misses={misses}")
        if total_miss == 0:
            survivors.append(conv)
    print(f"=> surviving conventions: {[conv_name(c) for c in survivors] or 'NONE'}")
    return survivors


M8_BENCHMARKS = [(8, 2, 1), (8, 1, 1), (8, 1, 2)]  # rank 4 (pin), rank 6, bent rank 8


def stage_fifo_m8(convs) -> list[dict]:
    """Benchmark sweep with rank-4 pruning: a convention that misses on the
    rank-4 (8,2,1) instance cannot satisfy the all-765-forms claim."""
    print("== echo-fifo+dummy m=8 benchmarks (the unverified claim) ==")
    survivors = []
    for conv in convs:
        line = [conv_name(conv)]
        alive = True
        for m, a, lam in M8_BENCHMARKS:
            f = Form(m, a, lam)
            agree, misses = exactness(f, conv)
            line.append(f"({m},{a},{lam})r{f.rank()}: {agree}/256"
                        + (f" miss{misses[:4]}{'+' if len(misses) > 4 else ''}"
                           if misses else ""))
            if agree < 256:
                alive = False
                break  # claim requires exactness on every form; this conv is dead
        print(("OK  " if alive else "X   ") + "  ".join(line), flush=True)
        if alive:
            survivors.append(conv)
    print(f"=> m=8-exact conventions on all three benchmarks: "
          f"{[conv_name(c) for c in survivors] or 'NONE'}")
    return survivors


def stage_stratified(convs, per_stratum: int = 7) -> None:
    print("== >=20 stratified lambda at m=8 (strata = polar rank of Q_{lambda,a}) ==")
    by_rank: dict[int, list] = {}
    for a in (1, 2, 3):
        for lam in range(1, 256):
            f = Form(8, a, lam)
            by_rank.setdefault(f.rank(), []).append((a, lam))
    print("   family census by rank:",
          {r: len(v) for r, v in sorted(by_rank.items())})
    picks = []
    for r in sorted(by_rank):
        stratum = by_rank[r]
        step = max(1, len(stratum) // per_stratum)
        picks += [(r, al) for al in stratum[::step][:per_stratum]]
    print(f"   testing {len(picks)} stratified (a, lambda) pairs")
    for conv in convs:
        print(f"-- {conv_name(conv)}")
        for r, (a, lam) in picks:
            f = Form(8, a, lam)
            t0 = time.time()
            agree, misses = exactness(f, conv)
            flag = "EXACT" if agree == 256 else f"misses={len(misses)}: {misses[:8]}"
            print(f"   (8,{a},lam={lam}) rank {r}: {agree}/256  {flag}"
                  f"  ({time.time() - t0:.0f}s)", flush=True)


# ----------------------------------------------------------------------------
# The ACTUAL implemented echo-fifo+dummy rule (reconciliation phase).
#
# Reading the original probes (rule definition only) shows the SS8.3 prose
# misdescribes the game in three load-bearing ways:
#   * it is a sigma-VALUED steering game: the value of a position is the
#     terminal charge under optimal play, P1 wanting sigma = t and P2 wanting
#     1 - t, with t swept over both values ("exact" = value(x) == Q(x) for all
#     x and both t; 765 forms x 256 x x 2 t = 391,680 checks);
#   * the FIFO close-discipline RETAINS the one-move ko (the queue-front close
#     is blocked right after the coin was opened);
#   * a stalled player PASSES and the pass clears the ko (no stall-loses).
# The dummy is an ordinary chargeless coin (q = 0, zero B row) adjoined to
# every board, fully FIFO-queued and ko-able.
# Charge conventions match the pinned echo-ko ones: pre-timing (open-set
# before the touch), lower-triangular cocycle, q_i collected at the close.
# ----------------------------------------------------------------------------


def fifo_value(form, x: int, t: int, dummy: bool = True,
               tree: bool = False) -> int:
    """Terminal charge of the echo-fifo(+dummy) game from x under optimal play.

    Direct stateful solver: full state (untouched, open-queue, ko, mover,
    sigma) in the memo key; no decomposition, no isomorphism caching.  With
    tree=True, plain enumeration with no memo at all (the validation oracle).
    """
    m = form.m
    bits = [i for i in range(m) if x >> i & 1] + ([m] if dummy else [])
    if not bits:
        return 0
    q = form.q + [0]
    B = [row + [0] for row in form.B] + [[0] * (m + 1)]
    memo: dict = {}

    def charge(omask: int, has_i: bool, i: int) -> int:
        acc = q[i] if has_i else 0
        for k in bits:
            if k > i and omask >> k & 1 and B[k][i]:
                acc ^= 1
        return acc

    def rec(u: int, openseq: tuple, last: int, mover: int, sigma: int) -> int:
        if u == 0 and not openseq:
            return sigma
        key = (u, openseq, last, mover, sigma)
        if not tree:
            r = memo.get(key)
            if r is not None:
                return r
        omask = 0
        for c in openseq:
            omask |= 1 << c
        legal = []
        for i in bits:
            if i != last and u >> i & 1:
                legal.append((i, charge(omask, False, i), u ^ (1 << i), openseq + (i,)))
        if openseq and openseq[0] != last:
            c = openseq[0]
            legal.append((c, charge(omask, True, c), u, openseq[1:]))
        if not legal:
            res = rec(u, openseq, -1, 1 - mover, sigma)  # forced pass clears ko
        else:
            want = t if mover == 0 else 1 - t
            res = 1 - want
            for (i, ch, u2, seq2) in legal:
                if rec(u2, seq2, i, 1 - mover, sigma ^ ch) == want:
                    res = want
                    break
        if not tree:
            memo[key] = res
        return res

    return rec(x | ((1 << m) if dummy else 0), (), -1, 0, 0)


def their_direct_solver():
    """Exec the ORIGINAL probe's direct m<=4 solver, verbatim, as a cross-oracle."""
    import pathlib
    src = pathlib.Path(__file__).parent / "gold" / "asym2_fifo_bench.py"
    text = src.read_text()
    start = text.index("def direct_fifo_value(")
    end = text.index("# fix call shape")
    ns: dict = {"sys": sys}
    exec(text[start:end], ns)  # noqa: S102 -- verbatim original, reconciliation oracle
    return ns["direct_fifo_value"]


def stage_fifo2_validate() -> None:
    print("== faithful rule: my direct solver vs tree enumeration (m=4, dummy on/off) ==")
    n = 0
    for lam in range(1, 16):
        f = Form(4, 1, lam)
        for t in (0, 1):
            for dummy in (True, False):
                for x in range(16):
                    assert fifo_value(f, x, t, dummy) == fifo_value(f, x, t, dummy, tree=True), \
                        (lam, t, dummy, x)
                    n += 1
    print(f"   {n} solves agree with explicit tree enumeration")
    f8 = Form(8, 1, 1)
    small = [x for x in range(256) if bin(x).count("1") <= 3][:40]
    n = 0
    for t in (0, 1):
        for x in small:
            assert fifo_value(f8, x, t) == fifo_value(f8, x, t, tree=True), (t, x)
            n += 1
    print(f"   m=8 (|supp| <= 3, dummy): {n} spot solves agree with tree enumeration")

    print("== my direct solver vs the ORIGINAL probe's direct solver (verbatim) ==")
    theirs = their_direct_solver()
    n = 0
    for lam in range(1, 16):
        f = Form(4, 1, lam)
        q5, B5 = f.q + [0], [row + [0] for row in f.B] + [[0] * 5]
        for t in (0, 1):
            for x in range(16):
                assert fifo_value(f, x, t, dummy=True) == theirs(x | 16, 5, q5, B5, t), \
                    ("dummy", lam, t, x)
                assert fifo_value(f, x, t, dummy=False) == theirs(x, 4, f.q, f.B, t), \
                    ("plain", lam, t, x)
                n += 2
    print(f"   {n} solves agree with the original direct_fifo_value")
    print("fifo2-validate: PASS")


def stage_fifo2_m4() -> None:
    print("== faithful echo-fifo+dummy, m=4 family exactness (value == Q, both t) ==")
    for dummy in (True, False):
        hits = []
        for lam in range(1, 16):
            f = Form(4, 1, lam)
            for t in (0, 1):
                miss = [x for x in range(16) if fifo_value(f, x, t, dummy) != f.Q[x]]
                if not miss:
                    hits.append((lam, t))
        print(f"   dummy={dummy}: EXACT (lam, t) pairs: {hits}"
              f"  [{len(hits)}/30]")


def stage_fifo2_m8(forms=None, dummy: bool = True) -> None:
    print(f"== faithful echo-fifo+dummy m=8 (direct, full-state memo; dummy={dummy}) ==")
    forms = forms or M8_BENCHMARKS + [(8, 1, 3)]
    for (m, a, lam) in forms:
        f = Form(m, a, lam)
        for t in (1, 0):
            t0 = time.time()
            misses = [x for x in range(256) if fifo_value(f, x, t, dummy) != f.Q[x]]
            agree = 256 - len(misses)
            mtxt = (" misses=" + ",".join(f"{x}(pc{bin(x).count('1')})" for x in misses[:8])
                    + ("..." if len(misses) > 8 else "")) if misses else ""
            print(f"   ({m},{a},lam={lam}) rank {f.rank()} t={t}: {agree}/256"
                  f"{' EXACT' if agree == 256 else ''}{mtxt} [{time.time()-t0:.0f}s]",
                  flush=True)


def stage_fifo2_stratified(per_stratum: int = 7) -> None:
    print("== faithful rule: >=20 stratified lambda at m=8 (strata = polar rank) ==")
    by_rank: dict[int, list] = {}
    for a in (1, 2, 3):
        for lam in range(1, 256):
            by_rank.setdefault(Form(8, a, lam).rank(), []).append((a, lam))
    print("   family census by rank:", {r: len(v) for r, v in sorted(by_rank.items())})
    picks = []
    for r in sorted(by_rank):
        stratum = by_rank[r]
        step = max(1, len(stratum) // per_stratum)
        picks += stratum[::step][:per_stratum]
    print(f"   testing {len(picks)} stratified (a, lambda) pairs, both t, dummy on")
    stage_fifo2_m8(forms=[(8, a, lam) for (a, lam) in picks])


class SynthForm:
    """A quadratic refinement (q', B) of a given polar form, for torsor sweeps.

    Q'(x) = sum_i q'_i x_i + sum_{i<j} B_ij x_i x_j over F_2 -- same B, any
    diagonal.  Duck-types the Form fields the solvers consume (m, Q, q, B).
    """

    def __init__(self, m: int, q: list[int], B: list[list[int]]):
        self.m, self.q, self.B = m, q, B
        self.Q = []
        for x in range(1 << m):
            v = 0
            for i in range(m):
                if x >> i & 1:
                    v ^= q[i]
                    for j in range(i + 1, m):
                        if x >> j & 1:
                            v ^= B[i][j]
            self.Q.append(v)


def stage_fifo2_torsor(per_form: int = 20) -> None:
    """Refinement uniformity: the rule must be exact for >=20 refinements of
    each benchmark polar form, not just the Gold member (single-refinement
    luck is the failure mode this rules out)."""
    import random
    rng = random.Random(0x06DD0AD)
    print(f"== torsor sweep: {per_form} refinements per benchmark B, both t, dummy on ==")
    for (m, a, lam) in M8_BENCHMARKS + [(8, 1, 3)]:
        base = Form(m, a, lam)
        diags = [base.q, [0] * m, [1] * m]
        seen = {tuple(d) for d in diags}
        while len(diags) < per_form:
            d = [rng.randint(0, 1) for _ in range(m)]
            if tuple(d) not in seen:
                seen.add(tuple(d))
                diags.append(d)
        worst = (256, None, None)
        all_exact = True
        t0 = time.time()
        for d in diags:
            f = SynthForm(m, d, base.B)
            for t in (0, 1):
                agree = sum(1 for x in range(1 << m) if fifo_value(f, x, t) == f.Q[x])
                if agree < 256:
                    all_exact = False
                if agree < worst[0]:
                    worst = (agree, d, t)
        tag = "ALL EXACT" if all_exact else f"worst {worst[0]}/256 at q'={worst[1]} t={worst[2]}"
        print(f"   ({m},{a},lam={lam}) rank {base.rank()}: {per_form} refinements x 2t: "
              f"{tag} [{time.time()-t0:.0f}s]", flush=True)


def fifo_liveness(form, x: int, t: int, dummy: bool = True) -> tuple[int, int]:
    """(decision_states, reachable_states): full-expansion count of states whose
    mover has both a want-achieving and a want-failing option."""
    m = form.m
    bits = [i for i in range(m) if x >> i & 1] + ([m] if dummy else [])
    if not bits:
        return 0, 0
    q = list(form.q) + [0]
    B = [list(row) + [0] for row in form.B] + [[0] * (m + 1)]
    memo: dict = {}
    decisions = 0

    def charge(omask, has_i, i):
        acc = q[i] if has_i else 0
        for k in bits:
            if k > i and omask >> k & 1 and B[k][i]:
                acc ^= 1
        return acc

    def rec(u, openseq, last, mover, sigma):
        nonlocal decisions
        if u == 0 and not openseq:
            return sigma
        key = (u, openseq, last, mover, sigma)
        r = memo.get(key)
        if r is not None:
            return r
        omask = 0
        for c in openseq:
            omask |= 1 << c
        legal = []
        for i in bits:
            if i != last and u >> i & 1:
                legal.append((i, charge(omask, False, i), u ^ (1 << i), openseq + (i,)))
        if openseq and openseq[0] != last:
            c = openseq[0]
            legal.append((c, charge(omask, True, c), u, openseq[1:]))
        if not legal:
            res = rec(u, openseq, -1, 1 - mover, sigma)
        else:
            want = t if mover == 0 else 1 - t
            good = bad = 0
            res = 1 - want
            for (i, ch, u2, seq2) in legal:
                sub = rec(u2, seq2, i, 1 - mover, sigma ^ ch)
                if sub == want:
                    good += 1
                    res = want
                else:
                    bad += 1
            if good and bad:
                decisions += 1
        memo[key] = res
        return res

    rec(x | ((1 << m) if dummy else 0), (), -1, 0, 0)
    return decisions, len(memo)


def stage_fifo2_liveness() -> None:
    print("== decision-liveness of the faithful rule on the m=8 benchmarks ==")
    for (m, a, lam) in M8_BENCHMARKS + [(8, 1, 3)]:
        f = Form(m, a, lam)
        t0 = time.time()
        for t in (0, 1):
            dec = reach = livepos = 0
            for x in range(1 << m):
                d, r = fifo_liveness(f, x, t)
                dec += d
                reach += r
                if d:
                    livepos += 1
            print(f"   ({m},{a},lam={lam}) t={t}: {dec} decision states / "
                  f"{reach} reachable, {livepos}/256 positions decision-live"
                  f" [{time.time()-t0:.0f}s]", flush=True)


def stage_fifo2_all() -> None:
    """The FULL claimed sweep, independently re-derived: all 765 scaled Gold
    forms (a in {1,2,3}, lam in F_256^*), all 256 positions, both wants --
    391,680 checks through the direct full-state solver."""
    checks = 0
    bad = []
    t_start = time.time()
    for a in (1, 2, 3):
        for lam in range(1, 256):
            f = Form(8, a, lam)
            for t in (0, 1):
                for x in range(256):
                    if fifo_value(f, x, t) != f.Q[x]:
                        bad.append((a, lam, t, x))
                    checks += 1
            if lam % 32 == 0:
                print(f"   a={a} lam={lam}: {checks} checks, {len(bad)} misses "
                      f"[{time.time()-t_start:.0f}s]", flush=True)
    print(f"TOTAL: {checks} checks, {len(bad)} misses [{time.time()-t_start:.0f}s]")
    if bad:
        print("   misses:", bad[:40])
    else:
        print("   FULL m=8 EXACTNESS REPRODUCED: 391,680/391,680")


def ko_value(form: Form, x: int, t: int, tree: bool = False) -> int:
    """Terminal charge of the echo-KO game (sigma-valued, faithful semantics:
    any open coin closable, one-move ko, forced pass clears ko, P1 wants t).
    Sigma-EXPLICIT full-state memo -- the independent check on the original
    probe's sigma-free fixed-stance recursion."""
    m = form.m
    bits = [i for i in range(m) if x >> i & 1]
    if not bits:
        return 0
    q, B = form.q, form.B
    memo: dict = {}

    def charge(omask: int, has_i: bool, i: int) -> int:
        acc = q[i] if has_i else 0
        for k in bits:
            if k > i and omask >> k & 1 and B[k][i]:
                acc ^= 1
        return acc

    def rec(u: int, o: int, last: int, mover: int, sigma: int) -> int:
        if u == 0 and o == 0:
            return sigma
        key = (u, o, last, mover, sigma)
        if not tree:
            r = memo.get(key)
            if r is not None:
                return r
        legal = []
        for i in bits:
            if i == last:
                continue
            if u >> i & 1:
                legal.append((i, charge(o, False, i), u ^ (1 << i), o | (1 << i)))
            elif o >> i & 1:
                legal.append((i, charge(o, True, i), u, o ^ (1 << i)))
        if not legal:
            res = rec(u, o, -1, 1 - mover, sigma)
        else:
            want = t if mover == 0 else 1 - t
            res = 1 - want
            for (i, ch, u2, o2) in legal:
                if rec(u2, o2, i, 1 - mover, sigma ^ ch) == want:
                    res = want
                    break
        if not tree:
            memo[key] = res
        return res

    return rec(x, 0, -1, 0, 0)


def stage_ko2() -> None:
    print("== sigma-explicit echo-ko vs tree enumeration (m=4) ==")
    n = 0
    for lam in (1, 2, 7, 12):
        f = Form(4, 1, lam)
        for t in (0, 1):
            for x in range(16):
                assert ko_value(f, x, t) == ko_value(f, x, t, tree=True), (lam, t, x)
                n += 1
    print(f"   {n} solves agree with tree enumeration")
    print("== sigma-explicit echo-ko vs the corrected SS8.2 table (value == Q) ==")
    for (m, a, lam, want_agree, want_miss) in KO_FIXTURES:
        f = Form(m, a, lam)
        for t in (1, 0):
            miss = [x for x in range(1 << m) if ko_value(f, x, t) != f.Q[x]]
            agree = (1 << m) - len(miss)
            note = f" misses={miss}" if 0 < len(miss) <= 8 else ""
            pin = ""
            if agree == want_agree and (want_miss is None or sorted(miss) == sorted(want_miss)):
                pin = "  [matches corrected table]"
            print(f"   ({m},{a},lam={lam}) t={t}: {agree}/{1 << m}{note}{pin}", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    ap.add_argument("stage", choices=["selftest", "pin-ko", "fifo-m4", "fifo-m8", "stratified",
                                      "fifo2-validate", "fifo2-m4", "fifo2-m8",
                                      "fifo2-stratified", "fifo2-torsor",
                                      "fifo2-liveness", "fifo2-all", "ko2"])
    ap.add_argument("--conv", action="append", default=None,
                    help="fifo convention as dummy/timing/orientation, e.g. ignore/pre/lower")
    args = ap.parse_args()

    if args.stage == "selftest":
        stage_selftest()
    elif args.stage == "pin-ko":
        stage_pin_ko()
    elif args.stage == "fifo-m4":
        stage_fifo_m4()
    elif args.stage == "fifo2-validate":
        stage_fifo2_validate()
    elif args.stage == "fifo2-m4":
        stage_fifo2_m4()
    elif args.stage == "fifo2-m8":
        stage_fifo2_m8()
    elif args.stage == "fifo2-stratified":
        stage_fifo2_stratified()
    elif args.stage == "fifo2-torsor":
        stage_fifo2_torsor()
    elif args.stage == "fifo2-liveness":
        stage_fifo2_liveness()
    elif args.stage == "fifo2-all":
        stage_fifo2_all()
    elif args.stage == "ko2":
        stage_ko2()
    else:
        if args.conv:
            convs = []
            for spec in args.conv:
                d, t, o = spec.split("/")[:3]
                convs.append({"discipline": "fifo", "dummy": d, "timing": t, "orientation": o})
        else:
            convs = FIFO_CONVS
        if args.stage == "fifo-m8":
            stage_fifo_m8(convs)
        else:
            stage_stratified(convs)


if __name__ == "__main__":
    main()
