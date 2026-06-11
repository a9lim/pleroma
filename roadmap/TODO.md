# Cross-pillar work — TODO (the game-valued ledger)

Every bridge that was *explicitly on the build order* is done — the full record is in
[`roadmap/CODA.md`](CODA.md), and newly completed work goes in the
[`roadmap/DONE.md`](DONE.md) ledger. This file is the ledger of what is **buildable
but not built**: standard math made computational, verification harnesses, and elbow
grease. Nothing here is a genuine research question — those live in
[`OPEN.md`](../OPEN.md) (which carries the loopy-valued entries; open problems give
no termination guarantee).

Claim-level discipline (`AGENTS.md` → "Claim levels and non-claims") applies to every
item: each is **standard math** or **engineering** when built — not a new theorem.

## How items are valued

Natural numbers don't do roadmap items justice, so the ledger is a **game-valued
multivector**: each item is a term `g·e_B` — a game value `g` (its size and temper)
on a pillar blade `e_B` (which pillars it joins; the blade's grade is how
cross-cutting the item is). Blades: `e_s` scalar, `e_c` clifford, `e_f` forms,
`e_i` integral, `e_g` games, `e_y` py; pure-prose chores are scalar-grade (no blade).

| value | temper | meaning |
|---|---|---|
| `n` (numbers) | cold | buildable now; `n` ≈ focused days; `½` ≈ an afternoon |
| `±n` (switches) | hot | a real scope decision belongs to a9 first; size `n` either way |
| `↑` (ups) | infinitesimal | worth less than any number, still strictly positive |
| `*n` (stars) | confused with `0` | deferred not-yet-numbers: real, on-thesis, unscheduled |

Reference items by **slug**. The ledger's total value is the disjunctive sum; play it
in any order. (`echo-solver`, the formerly hottest cold item, was played 2026-06-10
with outcome **CONFIRM** — see `roadmap/DONE.md`; its successor move is the
σ-recasting target in `OPEN.md` tis (§1), which is loopy-valued, not a number.)

---

## numbers — forms & Witt (the classifier spine)

### `echo-family-sweep` — 1·(e_g∧e_f)
**The remaining pre-registered family axes** (`writeups/goldarf.tex` §§8–9, ranked
move 2), on the shipped harness `experiments/echo_solver.py`: ko-memory window
`w ∈ {1,2,3}`, pass semantics (clears-ko / forbidden / loses), single-coin plus pair
touches (the tartan-companion axis), and no-dummy controls — mapping which disciplines
besides fifo+dummy are exact. No longer decisive for existence (the fifo+dummy verdict
is in); it bounds the *mechanism* and finally puts the bounded-window blocker
conjecture on valid data. (Partially advanced by the 2026-06-10 `linking-reduction`
pass, `roadmap/DONE.md`: the no-dummy controls are fully mapped at the abstract-graph
level — the Bad census — and the fifo+dummy mechanism is identified
(`experiments/linking_game.py`, goldarf §8 `sec:linking`); the `w ≥ 2` ko-window and
pass/pair axes remain unswept, and the general-n linking *proof* is loopy-valued in
`OPEN.md` tis (§1), not a number here.)

### `bw-rational` — 2·e_f
**The graded rational Brauer–Wall class** — the lift Bridge F deliberately stopped
short of ("do not conflate `Brauer2Class` with the graded `BrauerWallClass` *until a
rational Brauer–Wall story is separately modeled*"). Model `BW(ℚ)` by Wall's exact
sequence `0 → Br(ℚ) → BW(ℚ) → Q(ℚ) → 0`, `Q(k) = ℤ/2 × k*/k*²` with the twisted
cocycle multiplication (Wall 1964; Lam GSM 67). Oracles: `BW(ℝ) ≅ ℤ/8` recovers the
Bott index `bw_class_real`; the graded class of `Cl(q)` over ℚ projects to Bridge F's
`c(q)` and to disc/dimension data.

### `tame-symbols` — 2·e_f
**Tamely ramified cyclic symbols for Bridge K.** K's local invariant is scoped
unramified-at-`v` (`inv = v(a)/n`). For tame ramification (`n | q−1`) the norm-residue
symbol is explicit (the tame symbol `(a,b)_v = (−1)^{v(a)v(b)} a^{v(b)}/b^{v(a)} mod 𝔪`
in `μ_n`); shipping it widens `BrauerClass` to ramified cyclic classes at tame places
on the `Qq`/`F_q(t)` legs. The **wild** symbol stays out — that is star `*4` below.

### `nikulin-existence` — 1·(e_i∧e_f)
**Nikulin's even-lattice existence theorem on top of `FqmWittClass`.** The
`fqm-witt` pass shipped the finite-quadratic-module normal form and native
constructor; it did **not** implement Nikulin 1.10.1 (which `(signature, FQM)` pairs
are realized by even lattices). Build the predicate with the rank/signature
congruence, local length/parity inequalities, and exact citations; it should decide
existence, not enumerate lattices.

## numbers — the integral wing

### `niemeier` — 3·e_i
**The Niemeier zoo and the non-degenerate Siegel–Weil.** Bridge E's Siegel–Weil check
is degenerate at `n = 16` (both classes share `θ = E₄²`). Build the 24 Niemeier
lattices (root systems + glue codes, Conway–Sloane Ch. 16/18 — curated tables per
`TABLES.md` discipline; Kneser 2-neighbors from Leech is the principled alternative)
with their `|Aut|` orders, then verify the genuine statement at weight 12:
`(Σ_L θ_L/|Aut L|) / mass(24) = E₁₂ = 1 + (65520/691)Σσ₁₁(m)qᵐ` — the Eisenstein
congruence prime **691** appearing in the codebase, and the first Siegel–Weil instance
where the classes genuinely differ. Free cross-checks: `Σ 1/|Aut| = mass_even_unimodular(24)`,
rootless-class uniqueness (Leech), and 24 new Nikulin/N.3 test points.

### `padic-symbols` — 3·e_i
**Conway–Sloane `p`-adic genus symbols** (Ch. 15, incl. the canonical 2-adic symbol
with trains/compartments/oddity fusion). Gives exact genus equality without the
budgeted `is_isomorphic` search — un-capping `ISO_GROUP_CAP` for genus decisions and
turning N.3 from "budgeted cross-check" into "exact symbol comparison". Oracle: full
agreement with the shipped `are_in_same_genus` + `DiscriminantForm::is_isomorphic`
route across the ADE zoo and the Milnor pair.

### `odd-lattices` — 2·e_i
**Type I (odd) lattices for the discriminant chain.** `DiscriminantForm::from_lattice`,
Milgram, Brown, and `theta_series` are all even-only. Ship the odd-lattice discriminant
form (`q_L` valued in `ℚ/ℤ` rather than `ℚ/2ℤ`), the oddity-corrected Milgram statement
(the `genus.rs` oddity is already the independent oracle), and odd Construction A
(`ℤⁿ` from the `[n,1]` repetition code's complement, `ℤ ⊕ E₈`, …). Document the theta
boundary honestly (odd `θ` lives at level 4).

### `constructions-bd` — 1·e_i
**Constructions B and D** (Conway–Sloane Ch. 5/8), extending Bridge H's code↔lattice
seam: B from doubly-even subcode data (oracle: `B(golay)` has the right det/min), D
from nested code towers. Keeps the same `Option`-on-non-integral-Gram boundary as
`construction_a`.

### `eichler` — ½·(e_i∧e_c)
**Eichler's theorem as a documented predicate** — the one cheap honest piece of star
`*1`: *indefinite, rank ≥ 3 ⇒ spinor genus = isometry class*, letting `Genus` upgrade
to a class statement in exactly that regime. No adelic machinery; just the predicate,
its citation (Eichler; Cassels), and tests on indefinite Grams. The full definite
computation stays `*1`.

## numbers — scalar worlds

### `hyperfield` — ½·e_s
**Viro's tropical hyperfield**, making Bridge J's lax tropicalization strict (Remark
J.2 names this exact repair): a small multivalued-addition type
(`x ⊞ y = {min}` off the vanishing locus, the interval/set on it) with the hyperfield
laws as tests and `tropicalize` factoring through it. A leaf, but it converts the one
"lax" asterisk in the J appendix into a theorem about a shipped type.

## numbers — games

### `lexicode-game` — 1·e_g
**The turning-game realization of lexicodes** — Bridge O cites the Conway–Sloane
game construction "for transcription in a formalization pass". Build the actual
turning-game move structure whose Grundy-0 positions are `L(n,d)`, so greedy = mex is
a `Game`-level theorem witnessed in code, not a comment. Subordinate to `OPEN.md` §1
(the solved degree-1 shadow), exactly as Bridge O says.

### `guy-smith` — 1·e_g
**Octal periodicity certificates.** Implement the Guy–Smith periodicity theorem (if
the Grundy sequence of an octal game repeats with period `p` over a window long enough
relative to the largest take, it is periodic forever — Winning Ways; Siegel CGT) as a
checked certificate, turning `octal_hunt`-style sweeps into proofs-of-periodicity
rather than bounded observations. The *conjecture* that every finite octal game is
ultimately periodic is famous, external, and not ours to claim — the checker is.

## numbers — engine & bindings

### `spinor-gauge` — 2·e_c
**Spinor reps and reversal through the antisymmetric gauge.** `spinor_rep` and
`reverse()` reject general-bilinear (`a ≠ 0`) metrics; in char ≠ 2 the general engine
is gauge-equivalent to the orthogonal one (the antisymmetric part is a "gauge", the
symmetric part fixes the iso class). First pin the gauge isomorphism against the
shipped `reduce_word` oracle on this engine's conventions, then transport the spinor
construction and the reversal anti-automorphism through it. Char 2 keeps its own
boundary.

---

## switches (a9's move first)

### `surreal-completion` — ±2·e_s
**The ω-place completion of No** — a capped Hahn-window backend (`PrecisionScalar`
discipline, finite window of CNF terms) that finally represents `1/(ω+1)`, `√2`-as-
series, and divisible-Γ Newton polygons, completing the (exact global, capped local)
pattern every other leg has. The decision: whether No gets an inexact leg at all —
Surreal is currently the *exact* char-0 home, and the precedent (`Rational` as an
engine-validation scalar) cuts both ways. Divisible-Γ polygons are the research-edged
corner (CODA J: "definable but not claimed or scheduled").

### `theta-level` — ±3·e_i
**Level-`N` theta identification** — `θ_L ∈ M_{n/2}(Γ₀(N), χ)` for non-unimodular
even lattices. The decision: how much modular-forms machinery this crate wants to own
(dimension formulas, level-`N` Eisenstein bases, Sturm bounds) versus keeping the
full-level `SL₂(ℤ)` story as the deliberate boundary tied to `level()`. Worth a
design conversation before any code.

### `mass-32` — ±1·e_i
**Mass past rank 24.** `mass_even_unimodular` caps at 24 because the `i128` rational
model overflows. Serre's "more than 80 million classes" at rank 32 is one
factored-rational representation away — but the repo's fixed-width-carrier policy is
deliberate. Decision: admit a factored/big-rational carrier for this one corner, or
keep the cap as the honest model boundary.

---

## ups (infinitesimal, strictly positive)

### `ps-regularity` — ↑
Verify the regularity hypothesis of Plambeck–Siegel Thm 6.4 against the published
JCTA 2008 paper — load-bearing for goldarf Theorem C, flagged there as the cheap gate
(ranked move 5a). Literature work, no code.

### `octal-hunt-reframe` — ↑
`examples/octal_hunt.rs` hunts `(ℤ/2)^k` misère quotients with `k ≥ 2` — a target
goldarf Theorem C proves **empty** (group misère quotients have order ≤ 2). Retarget
the probe at non-group monoids / kernels where the quadric framing can still apply,
and have `p_set_as_f2` check its labeling is a monoid homomorphism.

### `docs-experiments` — ↑
Root `AGENTS.md` and `README.md` don't mention the `experiments/{gold,excess,audit}`
subdirectories (the rescued 2026-06-10 research-run probes backing `goldarf.tex`,
`excess.tex`, and `AUDIT.md`) or their not-CI-tested status. One layout-table line
plus a sentence each.

---

## stars (deferred — the not-yet-numbers, confused with zero)

## `*1` — spinor genus (was Bridge G)

Refine `genus → spinor genus → isometry class` via the spinor norm (Eichler;
Cassels–Hall). `clifford/spinor_norm.rs` is the right primitive in spirit, but the full
bridge is **not buildable from the current surface**: `spinor_norm` computes one versor's
norm, whereas the spinor genus needs the local spinor-norm *images* `θ(O(L ⊗ ℤ_p))` at
every prime, adelic class-group bookkeeping, and the proper/improper class distinction.

The one cheap, honest piece is **Eichler's theorem** as a documented predicate —
*indefinite, rank ≥ 3* ⇒ spinor genus = isometry class — which would let `Genus` upgrade
to a class statement in exactly that regime (now filed as the buildable `eichler` above).
The full definite-lattice computation is the larger build; it sits adjacent to the
roadmap, not inside it.

## `*2` — the char-`p` Drinfeld/Carlitz mirror of the integral pillar (large)

The entire `integral/` wing — even-unimodular `ℤ`-lattices, `θ`-series,
`M_*(SL₂ℤ) = ℂ[E₄, E₆]`, Construction-A codes, Leech — is char 0. The project already
ships **exact** `F_q[t] ⊂ F_q(t)`, the char-`p` global field, whose arithmetic carries a
complete mirror:

- the **Carlitz module** `C_t(x) = t·x + x^q` is the char-`p` analogue of `exp` / the
  lattice exponential; the mirror of `E₄, E₆` are **Drinfeld modular forms** for
  `GL₂(F_q[t])`, with Goss `ζ`-values mirroring the Eisenstein constants;
- rank-`r` `F_q[t]`-lattices mirror even-unimodular `ℤ`-lattices and their reduction
  theory;
- **Goppa / algebraic-geometry codes** from function fields tie straight back into the
  existing `codes.rs` Construction-A machinery — the same code↔lattice seam in char `p`.

This is the `No ↔ On₂` / char-0 ↔ char-2 move applied to the richest pillar — the most
on-thesis possible "new structure." But it is a genuine new wing (Drinfeld modules, the
Carlitz exponential, rank-`r` reduction theory): weeks of specialized work, worth starting
only as a *second headline pillar* rather than a task. References: Goss, *Basic Structures
of Function Field Arithmetic*; Gekeler, Drinfeld modular forms; Goppa / AG codes.

## `*4` — the wild local symbol (full local class field theory)

Bridge K's invariant is unramified-only; `tame-symbols` (above) would add the tame
slice. The remainder — norm-residue symbols for **wildly ramified** cyclic extensions
(degree divisible by the residue characteristic: Lubin–Tate formal groups, or Dwork's
explicit formula; the dyadic Hilbert symbol's big siblings) — is a genuine wing of
machinery over the capped local models, and the precision-model honesty questions are
real (wild symbols read deep unit structure, not just `v(a)`). Deferred, not rejected.
Nimbered `*4` rather than `*3`, since `*3 = *1 + *2` is already spoken for as the sum
of the other two stars.
