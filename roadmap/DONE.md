# Cross-pillar bridges — DONE (the go-forward ledger)

The running ledger of cross-pillar work **completed from here on**.

The cross-pillar bridge-building era (bridges **A–O** plus **K** — lattice/Clifford/
Brauer–Wall, the char-2 Arf classifier, Frobenius outermorphisms, the transfinite
Clifford engine, theta/modular forms, Construction-A codes, the Weil representation, the
rational and full-`ℚ/ℤ` Brauer invariants, Newton polygons, the Brown invariant, the
unification pass, lexicodes) closed with every non-deferred bridge shipped. Its full
record — mathematical contracts, surfaces, oracles, boundaries, and the formal proofs —
lives in [`roadmap/CODA.md`](CODA.md); the terse working-notes summary is in the
`AGENTS.md` files (root + per-pillar).

What remains unbuilt is tracked in [`roadmap/TODO.md`](TODO.md) — the game-valued
ledger of buildable items plus the deferred stars `*1`/`*2`/`*4`; genuine open
problems stay in [`OPEN.md`](../OPEN.md), loopy-valued: `tis`/`tisn`, `on`/`off`,
`over`/`under` (the old numerals §1–§4 survive as aliases).

## How to use this ledger

When a new piece of cross-pillar work lands, add a short entry here:

```
## <name> — <one-line what-it-connects>   (<date>)
**Pillars:** … ↔ …    **Claim level:** standard math / implemented-and-tested / …
- surface: the functions/types that shipped
- oracles: the tests that pin it
- boundaries: the honest non-claims
```

Promote anything proof-heavy or too long for working notes into `roadmap/CODA.md`, and
fold the one-line structural fact into the relevant `AGENTS.md`.

## `brown-reduce` — Brown invariant by reduction, not enumeration (2026-06-11)
**Pillars:** forms    **Claim level:** standard math made computational + implemented and tested
- surface: `forms/char2/brown.rs::brown_f2` now splits off the polar radical and
  reduces the nonsingular `ℤ/4`-quadratic core into odd lines and even symplectic
  planes, adding their Brown phases in `ℤ/8`; `BROWN_MAX_ENUM_RANK` is gone, so the
  public route no longer panics above rank 26.
- oracles: the old exact Gauss-sum enumeration route is retained as a test-only
  checker; tests exhaust every four-dimensional input, compare the former rank-26
  budget edge against enumeration, preserve the Arf-doubling/additivity/radical
  checks, and verify a rank-40 form reduces past the old ceiling.
- boundaries: the `u128` bitmask surface still bounds represented dimensions to
  `n ≤ 128`; Brown's symmetric polar category remains distinct from the Clifford
  char-2 alternating polar, with `double_f2` still the explicit bridge.

## `echo-solver` — the echo-fifo+dummy adversarial review: CONFIRM (2026-06-10)
**Pillars:** games ↔ forms    **Claim level:** implemented and tested
- surface: `experiments/echo_solver.py` (stdlib-only, maintained) — direct full-state
  σ-explicit solvers for the echo family (`fifo_value`, `ko_value`), the parameterized
  prose-reading sweep that *failed* (all 80+ readings of the old §8.3 prose die at
  m=8 — the prose misdescribed the rule), nim/mex self-test, torsor and
  decision-liveness instrumentation, and the full-sweep stage.
- verdict: the formerly **unverified** echo-fifo+dummy `m = 8` exactness claim is
  **CONFIRMED** — 391,680/391,680 checks (765 scaled Gold forms × 256 positions ×
  both stances), zero misses, re-derived with no decomposition and no isomorphism
  caching. The m=4 family is 30/30 with the dummy and only 15/30 without it; the
  §8.2 echo-ko table was independently reproduced (σ-explicit, miss x=224 included)
  and shown to be the σ=1-stance face of a stance-asymmetric rule.
- oracles: explicit no-memo tree enumeration (m=4 exhaustive + m=8 small supports);
  the original `direct_fifo_value` executed verbatim (1,920 agreeing solves); the
  Turning-Corners mex recurrence for the nim product; a Codex cross-run of every
  stage including nim_mul vs the original probe on all 65,536 pairs.
- boundaries: the verified realizer is **σ-valued** (forced terminal charge, not a
  P-set); normal/misère/loopy recasting and the even-`a` diagonal lemma remain open
  (`OPEN.md` tis (§1)); the bounded-*window* blocker conjecture is untouched (FIFO
  memory is unbounded); goldarf §8.3's old prose rule description was corrected in
  the same pass.

## `linking-reduction` — the echo-fifo+dummy mechanism, reduced and screened (2026-06-10)
**Pillars:** games ↔ forms    **Claim level:** standard math (reductions) + implemented and tested (screens) + open (general n)
- surface: `experiments/linking_game.py` (stdlib-only, maintained) — the abstract
  odd-close parity game, validated against `echo_solver.fifo_value` through the
  `SynthForm` bridge; all-iso-classes rigidity/Bad screens; the strict menu
  verification of the two-mode defender strategy. Writeup: goldarf §8 "the linking
  reduction and the general-m theorem" (`\label{sec:linking}`).
- reductions (proven, machine-validated): FIFO ⇒ no nesting ⇒ linked = overlap;
  `D = σ ⊕ und` flips only on odd-`deg_U(front)` closes (the σ-game IS the
  odd-close parity game); ko localizes to fronts opened onto an empty queue and
  passes occur only after `U = ∅` (irrelevant to the flip fight).
- verified: the **linking theorem** (isolated coin ⇒ flips forced even ⇒
  `σ = |E| mod 2 = Q(x)` on Gold boards, i.e. exactness for all m) holds on every
  graph iso class `k ≤ 7` + dummy, both seats (1,044 classes at k=7) — beyond
  Gold-arising boards. No-dummy Bad census `{3:1, 5:4, 7:34}`: all
  mover-controlled, none with an isolated vertex, 33/34 dominated at n=7; all
  even-n boards rigid. The domination device (queue empty, `v` dominating an
  even nonempty remainder ⇒ forced flip in two plies) is the unique local
  obstruction; the dummy kills it at every root. Two-mode defender strategy
  (prevention/debt menus) strict-verified `k ≤ 7`, both seats —
  menu-existential semantics (a Codex review exhibited a losing in-menu choice).
- boundaries: the general-n proof is **open** — parity-local invariants provably
  insufficient (feature-mining inconsistency); the working architecture
  (Codex spar, thread 019eb4ff-695b-7762-97e8-c0bea66c4e7e) is firewall
  segmentation + mutual no-debt/one-debt induction + certificate-depth
  completeness, with the poison transition as the hard obligation.
