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

## `subfield-detect` — finite ordinal-nimber subfield detection (2026-06-11)
**Pillars:** scalar ↔ forms    **Claim level:** standard math made computational + implemented and tested
- surface: `Ordinal::finite_subfield_degree`, `scalar::ordinal_finite_subfield_degree`,
  and `scalar::ordinal_common_finite_subfield_degree` detect the minimal represented
  finite field `F_{2^m}` by generator support plus Frobenius minimization. The forms
  side exposes `forms::ordinal_metric_finite_subfield_degree` and routes
  `arf_ordinal_finite`, ordinal Witt classes, ordinal Brauer-Wall classes, and
  ordinal isometry through the detected/common degree.
- oracles: tests pin finite-nimber degrees, `ω` as degree 6, `ω^3` as degree 18,
  `ω^ω` as degree 20, `ω^(ω^2)` as degree 42, common-degree lcm behavior, inversion
  in a detected non-`F_64` field, ordinal Arf classification past the old `F_64`
  window, and rejection at `ω^(ω^ω)`.
- boundaries: detection is limited to the source-verified staged tower and the
  shipped Kummer excess table (`α_u` through `47`); genuinely transfinite
  ordinal-nimber metrics still return `None` for finite Arf/Witt/Brauer-Wall
  classification.

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

## `cyclic-trace-form` — cyclic algebra `Trd(z²)` trace form (2026-06-11)
**Pillars:** forms ↔ scalar    **Claim level:** standard math made computational + implemented and tested
- surface: `forms/trace_form.rs::cyclic_algebra_trace_form::<E>(&a)` builds the
  literal cyclic-algebra quadratic trace form `T_A(z)=Trd_A(z²)` for
  `A=(E/F,σ,a)=⊕ E·u^i`, ordered by the `E·u^i` lines. The `u^0` and, when present,
  `u^{n/2}` self-blocks reuse `assemble_twisted_form`; the `i`/`n-i` line pairs are
  pure polar blocks. The existing `MAX_BASIS_DIM=128` cap gives the boundary
  `[E:F]^2 ≤ 128`.
- oracles: over `Q(i)/Q`, the degree-2 trace form is pinned as
  `⟨2,-2,2a,2a⟩`, and a grid check verifies the honest reduced-norm relation
  `Trd(z²)=Trd(z)^2-2·Nrd(z)` against the shipped `⟨1,1,-a,-a⟩` norm form. Over
  `F_27/F_3`, the `u`/`u^2` cross block is checked to be Witt-hyperbolic.
- boundaries: this is not the reduced norm for general `n` (nor equal to it for
  quaternions); it is the quadratic trace companion named in CODA K §6(c).

## `milnor-d2` — dyadic cell of Milnor's exact sequence (2026-06-11)
**Pillars:** forms    **Claim level:** standard math made computational + implemented and tested
- surface: `forms/witt/milnor.rs::global_residues` now includes the `p=2` component
  of Milnor's residue map. The dyadic residue uses Milnor's hand convention, not the
  odd-prime Springer residue: a diagonal line contributes iff its `2`-adic valuation
  is odd, landing in the `W(F_2) ≅ Z/2` carrier
  `WittClassG::Char2 { field_degree: 1, arf }`.
- oracles: tests pin `⟨2⟩`, `⟨1,2⟩`, and `⟨-2⟩` as nonzero dyadic residue classes,
  verify `⟨2,-2⟩` cancels, check mixed support such as `⟨6⟩`, and cross-check
  reconstruction against `try_is_isotropic_q` for `⟨2⟩` vs `⟨8⟩` and `⟨2⟩` vs `⟨1⟩`.
- boundaries: this is only the dyadic `W(ℚ)` cell; the equal-characteristic
  `F_q(t)` twin shipped separately as `milnor-ff`.

## `milnor-ff` — split Milnor residues over `F_q(t)` (2026-06-11)
**Pillars:** forms    **Claim level:** standard math made computational + implemented and tested
- surface: `forms/witt/milnor.rs::global_residues_ff` implements the split
  equal-characteristic Milnor map for odd constant fields:
  `W(F_q(t)) ≅ W(F_q) ⊕ ⊕_π W(F_q[t]/π)`. The first component is the constant-field
  class selected by the even-valuation layer at the degree place `∞`;
  `FunctionFieldMilnorResidues<S>` records that class plus the finite vector of
  nonzero second residues at monic irreducible places, using the exact
  `local_global/function_field.rs` valuation, residue-unit, and residue-character
  helpers.
- oracles: tests pin constant forms, the `t` place, nonsquare constants, a degree-2
  irreducible place over `F_5`, square-multiple invariance, hyperbolic cancellation,
  and radical-entry rejection.
- boundaries: odd constant fields only (`FiniteOddField`); characteristic-2
  function fields stay on the separate Artin-Schreier/Aravire-Jacob layer, and this
  does not claim tame or wild norm-residue symbols beyond the second-residue Witt
  map.

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
