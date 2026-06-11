# OPEN: Genuine Research Problems

This file is intentionally narrow. It lists directions from repo audits, roadmap
splits, and the draft notes that look like genuine new research rather than
implementation of known formulas, standard algorithms, or already-source-pinned
theory. Implemented mathematical facts and maintenance context live in
`README.md` and `AGENTS.md`; buildable work lives in `roadmap/TODO.md` (the
game-valued ledger — items there are referenced by slug from here).

Numbering: an open problem is a loopy game, played without a termination
guarantee, so every entry wears a value from the loopy-stopper lexicon — the
shipped catalogue (`games/loopy.rs`: `on`, `off`, `over`, `under`, `dud`) plus
Conway's swinging pair `tis`/`tisn`, whose arithmetic is itself still unbuilt
here (`roadmap/TODO.md` slug `loopy-partizan`); fittingly, the flagship problems
wear values the codebase cannot yet compute. The values come in dual pairs, and
so do the problems:

- **`tis`/`tisn`** (`{0|tisn}`/`{tis|0}` — "this is / this isn't") — the two
  game-native-quadratic-data questions: the outcome side (§1, where every round
  of constructions and no-gos swings the apparent answer) and the coefficient
  side (§2, where the obstructions lean *isn't*).
- **`on`/`off`** — the two transfinite-On₂ questions: the tower that climbs past
  every verified rung (§3), and the classifier that switches off beyond the
  finite windows (§4).
- **`over`/`under`** — the two mirror questions: the mod-8 spine above the Arf
  bit, and the MinPlus shadow beneath MaxPlus thermography.

The original numerals survive as aliases — the rest of the repo cites them.
`dud` stays unassigned: `dud + G = dud` for every `G`, and no problem has yet
earned absorbing the whole roadmap. May none ever.

## tis (§1) — natural Gold-quadric game rule

Find, or rule out under a precise naturality condition, a non-tautological game
rule whose P-positions are the zero set `{Q = 0}` of a game-built Gold quadratic
form.

The implemented bridge is already concrete. In a finite nimber field,

```text
x + y      = XOR = disjunctive sum of impartial game values
x * y      = nim product = Turning-Corners product value
x -> x^2   = Frobenius = diagonal product x*x
Tr(x)      = x + x^2 + ... + x^(2^(m-1))
Q_a(x)     = Tr(x * x^(2^a))
```

The Gold form `Q_a(x) = Tr(x^(1+2^a))` is therefore not just an abstract
characteristic-2 quadratic form; it is assembled from nim/game operations. The
Arf invariant then has the standard zero-count interpretation. For a nonsingular
quadratic form on `F_2^(2r)`,

```text
#{x : Q(x)=0} = 2^(2r-1) + (-1)^Arf * 2^(r-1).
```

For degenerate forms, the implementation uses the usual radical-adjusted count:
an anisotropic radical balances the values exactly, while an isotropic radical
scales the bias. So if a game had P-positions exactly `{x : Q(x)=0}`, Arf would
say which player wins from more starting positions and by what square-root-scale
margin. That interpretation is meaningful, but it is conditional; it does not
exhibit the game.

Why this is research:
- The repo already builds the Gold forms and tests several game routes. The
  missing datum is not code for `Q`; it is a play rule, or a definition of
  "natural" strong enough to make the question non-ad-hoc.
- Normal-play sums do not solve it. For impartial normal play the P-condition is
  `g_1 xor ... xor g_n = 0`, hence linear in Grundy coordinates, while
  characteristic-2 quadrics obey `Q(u+v) = Q(u) + Q(v) + B(u,v)`. The polar form is
  exactly the XOR-closure obstruction.
- Frame-blind rules are too symmetric, while rules that directly evaluate `Q`
  are too tautological. The open core is the middle: a fixed play rule that reads
  the bilinear/game structure as a quadratic outcome without being a disguised
  evaluator.

The lexicode shadow (standard math + interpretation; the solved linear case):

The degree-1 version of the question is classically solved, and it is rich.
Conway-Sloane lexicodes ("Lexicographic codes: error-correcting codes from game
theory", IEEE Trans. Inform. Theory 32 (1986) 337-348) are built by the greedy
lexicographic rule, which is the mex rule: the codewords are the Grundy-value-0
positions of an associated turning-game move structure, binary lexicodes are
linear *because of* Sprague-Grundy theory (XOR-closure is a game theorem, not a
coding theorem), and the length-24, d = 8 lexicode is the extended binary Golay
code. More generally, lexicodes over base `2^k` are closed under nim-addition
and are linear when the base is a Fermat power `2^(2^k)` — exactly the sizes at
which nim-multiplication makes the ordinals below the base a field. So natural,
fixed, non-tautological rules demonstrably realize rich *linear* codes as
P-sets; and the matching no-go (`writeups/goldarf.tex`, Theorem A:
every Winning Ways coin-turning P-set is the kernel of an `F_2`-linear map)
says linearity is also the ceiling for that architecture. Floor and ceiling
coincide at linear. Problem 1 is exactly whether the lexicode phenomenon admits
a quadratic refinement — a rule producing the XOR-closure failure that the
polar form `B` measures. Bridge O (built; see `roadmap/CODA.md`) makes the
lexicode chain executable (greedy = mex -> Golay -> Construction A -> theta);
that is context for this problem, not progress on it.

Current probe map:

- `forms::quadric_fit::fit_f2_quadratic` asks whether a subset of `F_2^k` is the
  zero set of a genuine quadratic polynomial rather than an affine set.
- `experiments/trace_form_arf.py` builds Gold forms and checks the Gold rank
  formula on the tested power-of-two fields.
- `experiments/gold_form_from_games.py` rebuilds the same form using literal
  Turning-Corners products on small fields.
- `experiments/tartan_bilinear.py` rebuilds the polar form from game products.
- `experiments/arf_win_bias.py` brute-forces value distributions and matches the
  Arf-predicted zero counts.
- `experiments/gold_family_survey.py` broadens from unscaled Gold forms to
  components `Tr(lambda*x^(1+2^a))`. Over `F_256`, for APN Gold exponents
  `gcd(a,m)=1`, 2/3 of nonzero `lambda` give bent components, reproducing the
  classical count. Bent forms are the cleanest target because `R(B) = {0}`.
- `experiments/framing_obstruction.py` shows that for tested Gold polar forms,
  the coordinate-frame quadratic refinement has Arf 0 and the diagonal term
  flips to the Gold form. The remaining problem is whether the diagonal framing
  `q_i = Q(e_i)` is itself game-natural.
- `experiments/misere_kernel.py` verifies the Plambeck-Siegel kernel obstruction
  concretely on `R8`: the kernel is `(Z/2)^2`, `P cap K = {0}` is linear, and the
  genuine misere P-element lies outside the group where a vector-space quadric
  framing applies.
- `examples/interactive_kernel.rs` confirms that arbitrary P-sets and direct
  `Q`-evaluators are easy, while the tested polar-form rules do not reproduce the
  Gold zero set.
- `examples/loopy_quadric.rs` adds Draw as a third route. The symmetric `B` rule
  has Loss-set equal to the radical `R(B)`, so it explains one small coincidence
  and then fails away from it.
- `examples/bent_route.rs` tests a bent Gold component. A `B` plus coordinate-frame
  rule reaches a bent quadric of the correct Arf class but not the specific Gold
  zero set; adding the naive per-coin Ising field leaves the quadric variety.

The program state (2026-06-10 — `writeups/goldarf.tex` §§5–9, backed by the
`experiments/gold/` probes):

- The naturality criterion asked for below now has a draft formalization — N1
  (decision-nondegeneracy), N2 (bounded framing access), N3 (strategic
  relevance / anti-clock). N3's exact formulation is itself an open definitional
  problem: the escape-edge construction passes N1–N3 while being morally a clock,
  and the natural repairs run into two-game criticality being unsatisfiable in
  two-class outcome semantics.
- A no-go ladder (Theorems B–H) kills Tier 1 outright and shows every known
  in-quarantine Tier-2 normal-play realizer is a clock. Five named escape hatches
  remain: loopy-Draw semantics, `t ≥ 2r−2` with anisotropic complement,
  Frobenius-aware access (where both the symmetry and oracle methods are provably
  silent), non-quarantined rules using the game-native `℘` diagonal source, and
  rank-1 / radical-anisotropic degenerate layers.
- The abelian obstruction conjectured here is now Lemma `abelian` in the draft:
  no commutative game monoid's intrinsic squaring realizes a nondegenerate polar
  form, so the quadratic datum must come from the move relation's directedness.
- The leading Tier-2 candidate was the `echo`-ko charge-counting family on the
  extraspecial cocycle, and its `echo`-`fifo`+dummy variant is now **verified**
  (2026-06-10, pre-registered adversarial review, `experiments/echo_solver.py`):
  full `m = 8` exactness across all 765 scaled Gold forms, both stances,
  391,680/391,680 checks re-derived by a fresh direct full-state solver — no
  decomposition, σ in the memo key, validated against tree enumeration and the
  original direct solver, with a second-model cross-run. Decision-live in bulk
  (1.5–4.4M decision states per benchmark instance), torsor-uniform across
  refinements of each `B`. Three honest boundaries: the realizer is
  **σ-valued** (it realizes `Q` as a forced terminal charge — the central
  character of the play word — not yet as a P-set in normal/misère/loopy
  semantics); the `echo`-ko table is stance-asymmetric (its exactness face is
  the σ=1 stance only, where `fifo`+dummy is exact at both); and the
  bounded-window blocker conjecture is untouched (the FIFO queue is unbounded
  memory). The recasting is now the load-bearing open step; the
  Plambeck–Siegel Thm 6.4 regularity gate is still slug `ps-regularity`.
- The mechanism behind the verified realizer is now reduced and largely
  explained (2026-06-10 second pass, goldarf §8 "linking reduction",
  `experiments/linking_game.py`): FIFO forces closes in opening order (no
  nesting, linked = overlap), the whole σ-game is equivalent to an
  **odd-close parity game** (only closing a queue front with an odd number
  of untouched neighbors flips the outcome bit), ko/passes localize away,
  and the **general-m linking theorem** — flips forced even on any board
  with an isolated coin, hence exactness for ALL m — is machine-verified
  for every graph isomorphism class through k = 7 (1,044 classes, both
  seats), far beyond Gold-arising boards. The dummy's role is identified
  (it defeats the unique local obstruction, the domination device, at
  every root — matching the no-dummy Bad-graph census 1/4/34 at n=3/5/7,
  all mover-controlled), and an explicit two-mode defender strategy
  (prevention/debt menus) is strictly verified through k = 7. What
  remains is the general-n induction (firewall segmentation
  architecture); parity-local invariants provably do not suffice.

The naturality dichotomy:

- **Tier 1: frame-blind, `G >= Sp(B)`: no.** If the move relation is invariant
  under the full symplectic group of the polar form, its P-set is a union of
  `Sp(B)`-orbits. In dimension at least 4, `Sp(B)` is transitive on `V \ {0}`, so
  invariant subsets are only `empty`, `{0}`, `V\{0}`, or `V`. These are not
  nondegenerate quadrics. Degenerate Gold forms require care because the no-go
  only constrains the nondegenerate core `V/R(B)`.
- **Tier 3: per-`x` evaluator circuit: yes, but tautological.** The circuit
  `Q_a(x) = Tr(x*x^(2^a))` is a fixed Galois-symmetric circuit of game operations,
  and Frobenius permutes its summands. Realized as a disjunctive sum of those
  subgames with inputs driven by `x`, its P-condition is exactly `{Q_a = 0}`.
  That is more structured than a lookup table, but the form is still fed in rather
  than produced by autonomous play.
- **Tier 2: fixed-rule middle: open.** Positions should be indexed by field
  elements, with one rule independent of the chosen `x`, and the single-position
  Grundy-zero / kernel / Loss / Draw set should be `{Q_a = 0}`. The rule may use
  the nim product, Frobenius, or coordinate-frame data if a naturality criterion
  justifies them, but it must not simply evaluate `Q_a(x)`.

The extraspecial-group reframing (interpretation; explains the misère obstruction):

A characteristic-2 quadratic form `Q` on `V = F_2^n` with polar form `B` is **the same
data as an extraspecial-type central extension**

```text
1 -> Z/2 -> E -> V -> 0,
```

whose commutator pairing is `B` and whose **squaring map** `x -> x^2` (landing in the
center `Z/2`) **is** `Q`, because `(xy)^2 = x^2 y^2 (-1)^{B(x,y)}` gives
`Q(x+y) = Q(x) + Q(y) + B(x,y)` for free. The Arf invariant is exactly what classifies
the two extraspecial 2-groups of order `2^{1+2n}` (the `D_8`-central-product "+" type
versus the `Q_8`-central-product "-" type). This is standard math — the Heisenberg /
Weil-representation picture, adjacent to the already-built Bridge I (`weil_s`/`weil_t`).

It bites on the misère probe. `experiments/misere_kernel.py` found that on `R8` the
kernel `K = (Z/2)^2` and `P cap K` is **linear** — the genuine misère P-element lies
outside the group where a vector-space quadric framing applies. The reframing **predicts
that obstruction**: a misère quotient is a *commutative* monoid, so its unit group is
abelian, hence its intrinsic commutator pairing is trivial, hence its squaring map can
realize only the **split** refinement (`B = 0`, `Q = 0` on that part). A *nondegenerate*
`B` — which a Gold form has on its nonsingular core — is the commutator pairing of a
**nonabelian** extraspecial group and therefore **cannot** arise from any abelian
structure's own multiplication. So the linear obstruction is forced, not unlucky, and the
quadratic datum `q_i = Q(e_i)` must enter from a genuinely **noncommutative** source —
which, in game terms, is the one structural noncommutativity normal/partizan play has and
the symmetric polar form `B` discards: the **first-/second-player asymmetry** (the
directedness of the move relation).

This yields a candidate **Tier-2 naturality criterion** strictly between the two solved
tiers: require the rule to realize the *extraspecial squaring map* of `B` — equivariant
under the extension `E`, **not** merely under `Sp(B)`. That sits properly between
frame-blind `Sp(B)` (Tier 1, the no-go) and direct `Q_a`-evaluation (Tier 3,
tautological), because `E` is a proper central extension of `V`: it carries the `q_i`
data structurally without being a `Q`-evaluator. Status: developed into the Tier-2
screen and no-go ladder of `writeups/goldarf.tex` §§5–6 (see the program-state block
above); it does not yet exhibit a game.

Concrete progress targets (aligned with the goldarf §9 ranked moves):
- ~~Adversarially verify or refute the `echo`-`fifo`+dummy `m = 8` exactness
  claim~~ — **done, CONFIRM** (2026-06-10; `experiments/echo_solver.py`, record
  in goldarf §8 and `roadmap/DONE.md`). The successor target: **recast the
  σ-valued charge readout into normal/misère/loopy outcome semantics**, or
  prove the recasting impossible — the step that converts the verified
  realizer into a Tier-2 witness in the original P-set sense. Alongside it:
  the family-boundary sweep (ko-window `w`, pass semantics, pair touches,
  no-dummy controls), which also puts the bounded-window blocker on valid data.
- Close the **general-n linking theorem** (the mechanism half, reduced
  2026-06-10): prove that the odd-close parity game on any graph with an
  isolated coin forces an even flip count from both seats. Verified for all
  classes k ≤ 7 with a strictly-verified two-mode strategy
  (`experiments/linking_game.py`); the open residue is the firewall-segmented
  no-debt/one-debt induction with certificate-depth completeness (goldarf
  §8). A proof upgrades the m∈{4,8} verification to exactness for all m.
- Repair or replace N3, the anti-clock axiom — the open definitional problem: the
  escape-edge construction passes N1–N3 while being morally a clock, and two-game
  criticality is unsatisfiable in two-class outcome semantics.
- Exhibit a fixed uniform rule satisfying N1, N2, and N3 simultaneously on a Gold
  quadric of core rank ≥ 6 — or close the remaining escape hatches (loopy-Draw,
  `t ≥ 2r−2` anisotropic, Frobenius-aware access, `℘`-sourced diagonals,
  rank-1/radical-anisotropic layers) with no-gos of their own.
- Enumerate the Frobenius-aware access window at `m = 4, 8` — the one hatch where
  both the symmetry-killing and oracle methods are provably silent.
- Decide whether the diagonal refinement `q_i = Q(e_i)` is game-native for all `a`:
  the `a = 1` case is answered affirmatively by the `℘`-construction
  (`Wp(w) = w·w + w`, verified at `m = 4..32`); the even-`a` analogue (the drifting
  dual `λ_a^{(m)}` tower) has no named preimage family beyond `m = 8`.
- Cheap gates: verify the Plambeck–Siegel Thm 6.4 regularity hypothesis (slug
  `ps-regularity`); enumerate conjugation-move rules on `E` (the left-translation
  kill of Theorem H does not apply to conjugation); exhaust the board-8 case of the
  `fifo` parity-pinning conjecture.

Relevant surfaces:
- `writeups/goldarf.tex`
- `experiments/open_question_probe.py`
- `experiments/framing_obstruction.py`
- `experiments/gold_family_survey.py`
- `experiments/misere_kernel.py`
- `examples/interactive_kernel.rs`
- `examples/loopy_quadric.rs`
- `examples/bent_route.rs`
- `src/forms/quadric_fit.rs`
- `src/games/kernel.rs`, `src/games/misere.rs`, `src/games/loopy.rs`

## tisn (§2) — quadratic deformation of the game exterior algebra

Decide whether the current `GameExterior` construction admits a genuinely
game-native quadratic deformation on torsion-carrying game subgroups, rather than
only the all-zero Grassmann metric.

What is implemented:
- `GameExterior` is deliberately the exterior algebra of the game group. It uses
  the `Z`-module structure of games under disjunctive sum and can include non-number
  games such as `*` and `up`.
- Relation propagation is quotient-aware. If the game group imposes a relation,
  the exterior ideal respects it; for example, torsion in grade 1 propagates to
  torsion constraints in higher grades.
- This does not pretend that arbitrary games form a scalar ring. The construction
  is an exterior algebra over an abelian group, not a Clifford algebra over games.

Why this is research:
- A Clifford deformation would require extra quadratic data compatible with the
  game-group relations. Over torsion-free integer coefficients, a relation such as
  `2* = 0` forces any bilinear pairing involving `*` to vanish, and also forces a
  `Z`-valued quadratic value on `*` to vanish.
- Supplying an arbitrary quotient-compatible bilinear/quadratic table is a bounded
  implementation exercise. The research question is whether there is a natural,
  non-tautological source of such data from game structure itself.
- Torsion and mixed torsion/free subgroups make this sharper than "add a metric":
  the coefficient target, polarization identity, and relation compatibility all
  matter.

Concrete progress targets:
- Formalize the algebraic object: a quadratic map on a game subgroup, its
  coefficient ring or module, its polar pairing, and the exact compatibility
  condition with integer game relations.
- Prove obstruction results for torsion generators and mixed torsion/free subgroups
  under `Z`-valued or torsion-free coefficient targets.
- Identify coefficient targets where torsion can support nonzero quadratic data,
  and decide whether those targets are game-native or merely chosen by hand.
- Exhibit a nonzero deformation on a restricted class of games, or prove that every
  natural relation-respecting deformation collapses to the Grassmann one.
- Separate any useful engineering artifact, such as a checked
  `GameClifford::with_quadratic_data`, from the stronger mathematical claim that
  the data is game-native.

Relevant surfaces:
- `src/games/game_exterior.rs`
- `src/games/AGENTS.md`
- `examples/tour.rs`
- `demo.py`

## on (§3) — ordinal nim multiplication beyond the verified excess table

Push transfinite nim multiplication beyond the source-verified Lenstra-DiMuro
excess table. Historically the first missing carry in this checkout was
`alpha_47`; a local fixed-base finite-field oracle now verifies that carry, but
the general closed-form problem remains open.

What is implemented:
- The algebraic closure of `F_2` is represented by ordinals `< omega^(omega^omega)`
  under nim-arithmetic.
- The prime-power generator tower is implemented in `src/scalar/big/ordinal/tower.rs`.
  Products are exact when every Kummer carry uses a verified excess `alpha_u` for
  an odd prime `u <= 47`: DiMuro Table 1 through `43`, plus the local
  `ordinal_excess_probe.py` verification for `47`.
- Stage 1 handles scalar excesses such as `alpha_3 = 2`, `alpha_5 = 4`, and
  `alpha_17 = 16`; Stage 2 handles nonscalar excesses such as `alpha_7 = omega+1`
  by branching the monomial and recursing to lower places.
- Rows through `43` are from DiMuro's source table; `47` is from the independent
  local fixed-base probe. Field-axiom sweeps test engine consistency, not the
  truth of the table values.

The verified rows currently used are:

| u | alpha_u | u | alpha_u | u | alpha_u |
|---|---|---|---|---|---|
| 3 | 2 | 13 | omega+4 | 29 | omega^(omega^2)+4 |
| 5 | 4 | 17 | 16 | 31 | omega^omega+1 |
| 7 | omega+1 | 19 | omega^3+4 | 37 | omega^3+4 |
| 11 | omega^omega+1 | 23 | omega^(omega^3)+1 | 41 | omega^omega+1 |
| | | 43 | omega^(omega^2)+1 | 47 | omega^(omega^7)+1 |

Current external state:
- The first OEIS unknown in the extended table is now `p = 719`, where
  `f(719) = 359` and `Q(359) = {359}`. The calculator notes the required finite
  exponent as `e_719 = 1258230380`, which is the practical wall for the direct
  Lenstra power test.
- A tempting pattern matches the checked OEIS/calculator records from this pass:
  `m_p = 0` when `Q(f(p))` is not a singleton odd prime-power; `m_p = 1` for a
  singleton odd `Q(f(p))`, except the observed `f(p) = 2*3^k` cases have
  `m_p = 4`. A local audit matched this rule against the 950 calculator records
  with known `Q`-sets, and against every OEIS-known row covered by those `Q`-sets.
  This is still only a candidate rule, not a theorem.
- The exact finite-field reformulation is sharper than root-search language. If
  `beta = kappa_{f(p)} + m` lies in the component field `F_{2^E}`, then `beta`
  has no `p`-th root exactly when `p` divides the multiplicative order of `beta`.
  Thus the excess is the least `m` such that
  `p | ord(kappa_{f(p)} + m)`.
- The local fixed-base probe uses that criterion to verify `m_47 = 1` from the
  lower verified rows. Since `f(47) = 23` and `Q(23) = {23}`, this gives the newly
  shipped carry `alpha_47 = omega^(omega^7)+1`.

Since the 2026-06 research pass (`writeups/excess.tex`, `experiments/excess/`,
`experiments/cyclotomic_3k_family.py`):

- The 3-power column is now structural: `C_k` — the exact formula
  `ord(kappa_{3^k} + 1) = 3^(k+1) * (2^(3^k) - 1)` with `gamma_k` primitive — is
  certified for `k <= 6` and consistent-but-uncertified for `k = 7, 8`, blocked
  only by the unfactored cofactors of `Phi_{3^7}(2)` and `Phi_{3^8}(2)` (FactorDB
  CF). Whether ECM/GNFS reaches those on a realistic budget is open.
- The `f(p) = 2*3^k` exception column is proved unconditionally to have
  `m_p >= 4` (including the new example `p = 87211`); whether `m_p = 4` *exactly*
  is open — the `kappa + 4` translate lives in the degree-`4*3^k` compositum, one
  level above the half-angle toolkit, and a failure of the splitting there is
  where any `m >= 5` counterexample would hide.
- Wieferich caveat: the order criterion `m_p = min m : p | ord(kappa_{f(p)} + m)`
  is valid only when `v_p(2^(f(p)) - 1) = 1`. The two known base-2 Wieferich
  primes `1093` and `3511` sit inside the extended range and need the full power
  criterion.
- Newly certified `m_r = 1` rows (`262657` at `f = 27`; `71119` and `97685839` at
  `f = 81`; representatives at `f = 243, 729, 2187, 6561`) keep the candidate
  `0/1/4` rule unbroken. Still no proof; boundedness outside the 3-power and
  `2*3^k` columns (the 11-chain, the 23/29/47 components) has no structural
  theory, and no `m_p >= 5` example is known.
- `p = 719` feasibility: the direct test needs ~3.5 million Frobenius steps in
  `F_{2^1258230380}`; tower-aware Frobenius arithmetic (De Feo–Randriam–Rousseau
  standard lattices) is the conjectured 10–100x lever — a cost model, not a
  theorem.

Why this is research:
- Rewriting the current table-driven code to compute the known shape
  `f(u)`, `Q(f(u))`, and the `chi`-sum, while hardcoding only the finite excess
  integer, is a useful implementation improvement but not new reach.
- Extending past the verified table is different. DiMuro's theorem proves that the
  excess has a formulaic transfinite shape plus a finite correction, but the finite
  correction has no closed form in the cited theorem.
- Weaker "closed forms" already fail: `Q(f(p))` alone does not determine the
  excess, since `Q = {9}` gives `m_19 = 4` but `m_73 = 1`; similarly
  `Q = {81}` gives `m_163 = 4` but `m_2593 = 1`, and `Q = {243}` gives
  `m_1459 = 4` but `m_487 = 1`.
- The candidate `0/1/4` rule above would imply a global bound `m_p <= 4`. Lenstra
  explicitly left absolute boundedness open after proving lower-bound rules such
  as singleton-odd `Q(f(p))` forcing positive excess and `f(p)=2*3^k` forcing
  excess at least `4`.
- The order formulation explains the first weak-formula failures without appealing
  to the production table. In the independent probe, `ord(kappa_9 + 1) =
  3^3*(2^9 - 1)`, so `73 | ord(kappa_9 + 1)` but `19` does not divide it; adding
  `4` changes the order and picks up `19`. This is why the same `Q = {9}` gives
  both `m_73 = 1` and `m_19 = 4`.
- Shipping new values would require an independent oracle, a root-search theorem,
  or a new algorithmic proof. Otherwise the project would be numerology with a
  pleasant API.

Concrete progress targets:
- Implement the principled same-coverage route: compute `f(u) = ord_u(2)`,
  compute `Q(f(u))`, construct the `chi`-sum, and hardcode only the finite excess
  integer. This should independently cross-check the published rows. (Filed as
  `roadmap/TODO.md` slug `ordinal-principled` — implementation, not research.)
- Decide whether to import more known OEIS/calculator values through `p <= 709` as
  cited data, or keep requiring a local finite-field oracle for each shipped row.
- Derive or certify finite excess terms beyond the published table.
- Prove or find a counterexample to the candidate `0/1/4` rule. The smallest
  pressure point is `p = 719`, where the rule predicts `m_719 = 1` but the direct
  calculator path is too large for ordinary local verification.
- Turn the order-divisibility criterion into an actual theorem about the prime
  divisors of `ord(kappa_q + m)`, especially for singleton odd `Q = {q}` and for
  the exceptional tower `q = 3^k`.
- Build a verified `u`-th-power/root-search oracle for the transfinite field.
- Prove enough about the search to avoid merely empirical extensions.
- Decide what evidence is acceptable for shipping `alpha_53` and beyond.

Relevant surfaces:
- `writeups/excess.tex`
- `experiments/ordinal_excess_probe.py`
- `src/scalar/big/ordinal/tower.rs`
- `src/scalar/big/ordinal/mod.rs`
- `src/scalar/AGENTS.md`
- `examples/tour.rs`

## off (§4) — transfinite Arf/Witt classification for ordinal-nimber coefficients

Decide what, if anything, should replace the finite-field Arf/Brauer-Wall bit for
`CliffordAlgebra<Ordinal>` metrics whose coefficients do not all lie in one finite
nim-subfield.

What is implementation, not research:
- `roadmap/CODA.md` Bridge D is the tractable engine bridge: make `Ordinal` usable as a
  checked Clifford coefficient domain on the source-verified tower, and test the
  Clifford relations for genuinely transfinite squares such as `omega`.
- If all metric entries lie in a common finite nim-subfield `F_{2^d} ⊂ On₂`,
  classification should route through the generic finite characteristic-2 Arf
  classifier from Bridge B after detecting that subfield.
- The finite-field answer is an `F₂` bit because the absolute trace
  `Tr_{F_{2^d}/F₂}` exists. That finite-subfield case should stay separated from
  the genuinely transfinite case.

Why this is research:
- For genuinely transfinite ordinal-nimber coefficients there is no finite degree,
  so the finite trace-to-`F₂` definition of the Arf bit does not apply as-is.
- General characteristic-2 quadratic form theory has invariants over the
  coefficient field, such as Artin-Schreier quotient data, but the repo's current
  finite-nimber facade is an `F₂`-valued Arf/BW classifier. Deciding the right
  computable invariant for the represented ordinal-nimber domain is not just
  genericizing `arf_nimber`.
- The implemented ordinal multiplication itself is partial outside the verified
  Kummer tower. Any classifier that needs Artin-Schreier solving, roots, or field
  closure must respect that same source-verified boundary.

Concrete progress targets:
- Define the classification domain exactly: common finite subfields, the
  source-verified transfinite tower, or the ideal full `On_2` nimber field.
- Implement and test common finite-subfield detection so Bridge D can honestly
  delegate those metrics to Bridge B. (Filed as `roadmap/TODO.md` slug
  `subfield-detect` — implementation, not research.)
- Decide whether genuinely transfinite metrics should expose no classifier, a
  coefficient-field Arf class, a direct-limit finite-subfield invariant, or some
  other replacement for the finite trace bit.
- If an Artin-Schreier quotient or root-search route is chosen, build a checked
  oracle and prove enough about its represented domain to avoid table-driven
  guesses.
- State separately whether a Brauer-Wall class exists on the same surface, and
  whether it agrees with any proposed Arf-like invariant.

Relevant surfaces:
- `roadmap/CODA.md` Bridge D
- `src/scalar/big/ordinal/`
- `src/forms/char2/`
- `src/forms/witt/brauer_wall.rs`
- `src/clifford/`

## over — the mod-8 spine in game semantics

Decide whether the Brown invariant — the char-2 cell of the mod-8 spine, shipped as
Bridge M — has a game-theoretic reading the way the Arf bit does, i.e. whether the
conditional win-bias interpretation of `tis` (§1) lifts from `ℤ/2` to `ℤ/8`.

What is implemented (Bridge M, `forms/char2/brown.rs`): a `ℤ/4`-valued quadratic
refinement `q : V -> Z/4` has Gauss sum

```text
Sum_{x in V} i^(q(x)) = 2^(n/2) * zeta_8^beta,
```

read off the integer value-census Gaussian integer `(n0 - n2) + i*(n1 - n3)`, where
`n_k = #{x : q(x) = k}`. Doubling a classical char-2 form gives `beta = 4*Arf` — the
shipped win-bias bit embeds as the 2-torsion `{0, 4}` of `ℤ/8`.

Why this is research:
- The Arf reading is a **two**-class census: P-positions versus N-positions, bias
  `2^(r-1)` with sign `(-1)^Arf`. The Brown phase is a **four**-class census with a
  complex bias. No shipped outcome semantics has four classes: normal play has two,
  loopy play three (W/L/D). The question is whether any natural four-way outcome
  partition — loopy outcomes crossed with a parity, normal/misère outcome pairs, a
  mod-4 scoring residue, or something not yet named — produces the `zeta_8` phase of
  a game-built `ℤ/4`-form as its census.
- Game-built doubled forms only ever reach `beta in {0, 4}`. A genuinely odd `beta`
  needs `b` symmetric-but-not-alternating with `b_ii = q_i mod 2` — diagonal data
  again, one level up: this is the diagonal-framing problem of `tis` (§1) with the
  diagonal now *forced* by `q mod 2` rather than vanishing. The two problems are
  entangled, not parallel.
- The extraspecial picture of `tis` (§1) lifts: `ℤ/4`-valued forms correspond to
  central extensions by `ℤ/4` (the Pauli/complex-extraspecial family) exactly as
  `F₂`-forms correspond to extensions by `ℤ/2`. If the abelian obstruction
  (Lemma `abelian`) survives the lift, the four-class census also cannot come from
  any commutative game structure's own multiplication — which would make the
  first-/second-player asymmetry carry *three* extra bits instead of one.

Conditional claim, same shape as `tis` (§1): if a game's positions admitted a natural
four-class outcome census matching `i^q` for a game-built `q`, then `beta` would be
the phase and magnitude of its outcome imbalance — `sign mod 8` as a win-bias octant.
That interpretation is meaningful but conditional; it does not exhibit the game.

Concrete progress targets:
- Census probe: tabulate `(n0, n1, n2, n3)` for `ℤ/4`-refinements of game-built
  polar forms (doubled Gold forms first) and check which Gaussian integers actually
  arise on the game-reachable slice.
- Decide whether any existing three-class route (loopy W/L/D, `examples/loopy_quadric.rs`)
  extends by one natural axis to a four-class census with nonvanishing phase.
- Formulate the `ℤ/4` analogue of the abelian obstruction and prove or refute it.
- Connect to the lattice side: on 2-elementary discriminant forms `beta ≡ sign mod 8`
  (shipped); a game realizing `beta` would be a game computing a lattice signature.

Relevant surfaces:
- `src/forms/char2/brown.rs`, `src/forms/integral/discriminant.rs` (Bridge M)
- `src/games/loopy.rs`, `src/games/misere.rs`
- `writeups/goldarf.tex` §5 (the extraspecial reframing this lifts)
- `tis` (§1) — the `ℤ/2` floor of this question

## under — thermography ↔ Newton polygons: one tropical object or two?

Decide whether the project's two tropical consumers — thermography (`MaxPlus`, the
games axis) and the valuation/Newton-polygon stack (`MinPlus`, the place axis,
Bridge J) — are connected by a substantive transport, or whether the mirror is
purely notational. Either answer is the contribution; today the duality is named
(`scalar/tropical.rs` enforces the two-type separation) but carries no theorem.

Why this is research:
- On the place axis, the valuation axiom `v(x+y) >= min(v(x), v(y))` makes Newton
  polygons additive under multiplication (Dumas), and passing to the graded ring
  `gr_v` (Lemma J.3) is what "freezes" leading terms. On the game axis, the
  candidate analogue fails in the most interesting way: **thermographs of
  disjunctive sums do not compose** — that failure is precisely why temperature
  theory needs sidling and why `t(G+H)` is only bounded by `max(t(G), t(H))`, not
  determined. The open question is to make the failure structural: exhibit the
  exact lax/hyperfield law that thermographs *do* satisfy under `+` (a Viro-style
  repair, as Remark J.2 does for the valuation's own laxness), or prove no such
  law with nontrivial content exists.
- Sharper sub-question: is cooling a residue map? Cooling by `t` and "freezing" to
  the mast value is formally a leading-term extraction; does
  `(mast value, temperature)` behave like `(ac(x), v(x))` — i.e. is there a graded
  object `gr_t(Games)` whose pieces are the frozen values, with a multiplicative
  (Norton/overheating?) structure making the analogy a homomorphism rather than a
  pun? Berlekamp's economist's dictionary is the informal version; the question is
  whether it survives being made exact.
- The sign mirror is suggestive but not content: `MinPlus ↔ MaxPlus` is a convention
  flip. Content would be a single statement instantiating to Theorem J.5 (slopes =
  root valuations) on one axis and to a thermographic fact (masts/temperatures of a
  one-parameter family) on the other.

Concrete progress targets:
- Formulate and test the lax law for `t(G+H)` as a hyperfield statement; locate
  exactly where sidling violates strictness (the game-side "vanishing locus").
- Build the one-object probe: a polynomial family of games (e.g. switches with
  parameterized stakes) whose thermograph IS a Newton polygon under an explicit
  change of axes; determine whether the dictionary extends beyond the family or is
  an artifact of one-parameter linearity.
- Decide the graded-ring question for cooling: does Norton multiplication /
  overheating give `gr_t` a product compatible with freezing, in any restricted
  class of games?
- If every transport trivializes, write the no-go: the precise sense in which
  temperature is not a valuation (which axiom fails, on which games, measured how).

Relevant surfaces:
- `src/scalar/tropical.rs`, `src/games/` thermography, `src/scalar/newton.rs`
- `roadmap/CODA.md` Bridge J (the formal appendix, esp. J.1–J.3, J.5–J.6)
- `examples/tropical.rs` (the shipped thermography = tropical identity)

## References For The Open Threads

- Conway, *On Numbers and Games*: surreal numbers and nimbers.
- Berlekamp-Conway-Guy, *Winning Ways*: coin-turning games, Turning-Corners/nim
  product theorem, and thermography.
- Siegel, *Combinatorial Game Theory*: temperature theory and thermography.
- Arf, *Untersuchungen uber quadratische Formen...*: quadratic forms in
  characteristic 2.
- Dickson, *Linear Groups*: binary quadratic forms and zero-count bias.
- Ovsienko, *Real Clifford algebras and quadratic forms over F_2*: useful
  char-0/char-2 analogy, not a blanket nim-field Clifford classification theorem.
- Lidl-Niederreiter, *Finite Fields*: finite-field trace/Frobenius background and
  Gold-rank checks.
- DiMuro, *On Onp*: source table and theorem for transfinite nim Kummer excesses.
- Brown, *Generalizations of the Kervaire invariant*, Ann. of Math. 95 (1972):
  `ℤ/4`-valued quadratic refinements and the `ℤ/8` invariant (for `over`).
- Wall, *Quadratic forms on finite groups*, Topology 2 (1963): the Witt group of
  finite quadratic forms (for `over`).
- Plambeck-Siegel, *Misere quotients for impartial games*, JCTA 115 (2008): the
  quotient/kernel theory behind the misère obstruction (for `tis`, §1).
- Berlekamp, *The economist's view of combinatorial games*, in Games of No Chance
  (1996): the informal cooling dictionary (for `under`).
- Maclagan-Sturmfels, *Introduction to Tropical Geometry*; Viro, *Hyperfields for
  tropical geometry I*: valuations as (lax) tropicalization and the strictness
  repair (for `under`).
- De Feo-Randriam-Rousseau, standard lattices of compatibly embedded finite
  fields: the conjectured tower-aware Frobenius lever (for `on`, §3).
