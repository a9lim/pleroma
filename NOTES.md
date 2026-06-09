# Draft notes: what is actually novel here

Status: working draft, not a theorem paper. These notes try to separate the
mathematical content from the repo-tour material and from the speculative
game-semantics question.

The short version:

Pleroma is a Clifford-algebra lab over the field-like cores of Conway game
worlds. The main thread is not "Clifford over all games"; that phrase
is false because games under disjunctive sum are an abelian group, not a scalar
ring. The new thread is:

1. In characteristic 2, the Clifford engine keeps the quadratic data `q` and
   the polar data `b` independent, which is necessary for faithful char-2
   Clifford theory.
2. The nimber backend supplies finite characteristic-2 fields that are also
   game-value fields: nim addition is XOR and nim multiplication is realized by
   coin-turning games.
3. Gold trace forms
   `Q_a(x) = Tr(x * x^(2^a)) = Tr(x^(1+2^a))`
   are genuine quadratic forms on those fields. Their polar forms and their
   values can be built from nim/game operations: nim product, Frobenius
   squaring, and XOR.
4. The Arf invariant of those forms has an exact zero-count interpretation. If
   some natural game had P-positions `{Q = 0}`, then Arf would be the sign of
   the second-player win bias.
5. The open question is therefore sharp and limited: does a natural game rule
   produce a Gold quadric as its P-set? Current probes say normal-play sums and
   frame-blind `B`-only rules do not.

Everything else in the repo is supporting infrastructure unless it directly
feeds that chain.

## Claim levels

Use these labels when rewriting or extending the paper.

- **Theorem / standard math**: external facts such as Sprague-Grundy, the
  Turning-Corners product theorem, Arf classification of nonsingular binary
  quadratic forms, the Gold rank formula, and the zero-count formula for
  quadratic forms over `F_2`.
- **Implemented and tested**: statements backed by `cargo test`, examples, or
  experiments in this checkout.
- **Interpretation**: the bridge from zero-count to "win bias" is conditional on
  having a game whose P-set is `{Q = 0}`.
- **Open**: existence of a natural game with that P-set. The repo provides test
  benches and obstructions, not a solution.

## What the project is not claiming

- Not a Clifford algebra over arbitrary partizan games. A Clifford algebra needs
  a commutative scalar ring. The game group as a whole is not a ring.
- Not a new classification theorem for all characteristic-2 Clifford algebras
  over arbitrary fields. The code computes Arf data for finite nimber subfields
  and keeps rank/radical data explicit.
- Not a solved game-semantics theorem. The Gold form is built from game
  operations, but it is not currently the Grundy value or P-set of a natural
  game.
- Not evidence that every `B + frame` quadratic form is split. The split-frame
  claim is only supported for the Gold polar forms tested here.
- Not an algebraically closed finite backend. `Nimber(u128)` is the finite
  field `F_{2^128}`. It contains the finite nimber subfields `F_{2^m}` for
  `m = 1, 2, 4, ..., 128`, not all finite fields of characteristic 2.

## The scalar landscape

Combinatorial games under disjunctive sum form an abelian group. The scalar
rings/fields used here are only the field-like cores and adjacent coefficient
systems:

| backend | role |
| --- | --- |
| `Nimber(u128)` | finite nim-field `F_{2^128}` with nim add/mul; main char-2 backend |
| `Surreal` | finite-support Hahn/CNF characteristic-0 backend; the real-closed form table is available only on represented exact square classes |
| `Surcomplex` | `Surreal[i]`; the algebraically-closed form table is available only on represented exact square classes |
| `Integer`, `Omnific` | coefficient rings for exterior/nilpotent structures |
| `Fp`, `Fpn`, `Zp`, `WittVec` | comparison scalar worlds for the characteristic trichotomy |
| `Qp`, `Qq`, `Laurent`, `Ramified`, `Gauss` | local-field-style backends/functors, mostly precision models; used for valuation and Springer/Hilbert-symbol experiments |
| `Adele`, `LocalQp` | a runtime-prime adelic precision model over `Q`; used by the local-global form layer |
| `Ordinal` in `big/ordinal/` | staged transfinite nimbers: nim-addition on represented CNF terms; nim-multiplication via the prime-power generator tower (Kummer carries through `α_u`, `u ≤ 43`) |

The writeup should focus on `Nimber`, and mention the others only as context.

## The char-2 Clifford point

In characteristic 2, the quadratic form is not determined by the polar form.
For basis vectors:

```text
e_i^2 = q_i
e_i e_j + e_j e_i = b_ij
```

The polar form is alternating, so `b_ii = 0`, but `q_i` can be nonzero. If the
engine collapses `q` and `b` into one symmetric bilinear form, the nimber
backend loses the actual char-2 content.

The repo's relevant implementation:

- `src/clifford/engine.rs`: product engine over a generic `Scalar`.
- `Metric { q, b, a }`: `q` is the quadratic diagonal, `b` is the
  anticommutator/polar form, and optional `a` is an asymmetric contraction used
  for general bilinear-form products.
- `src/forms/char2.rs`: Arf computation and quadratic-set fitting.

For nonsingular forms over `F_2`, Arf is the complete invariant. For finite
characteristic-2 fields, the Arf value is read in `F / wp(F)` and pushed to
`F_2` by the field trace. The implementation reports:

```text
ArfResult {
  arf,
  rank,
  radical_dim,
  radical_anisotropic,
  o_type,
}
```

That extra data matters. For degenerate forms, Arf of the nonsingular core alone
is not the whole form.

## The game-built Gold forms

The game bridge is specific and concrete.

Known game fact: the Grundy value of Conway's Turning-Corners product is the
nim-product of the factors. In code:

- `src/games/coin_turning.rs::nim_mul_mex` implements the mex recurrence.
- tests compare it against the algebraic nim product.
- the slow game recurrence is used only on small fields.

This makes the following operations game-realizable:

```text
x + y          = XOR = disjunctive sum of impartial game values
x * y          = nim product = Turning-Corners product value
x -> x^2       = Frobenius = diagonal product x*x
Tr(x)          = x + x^2 + ... + x^(2^(m-1))
Q_a(x)         = Tr(x * x^(2^a))
```

The experiments then check:

- `experiments/trace_form_arf.py`: builds `Q_a` over the bit basis of
  `F_{2^m}` and checks the Gold rank formula
  `rank = m - gcd(2a, m)` for `m = 2,4,8,16,32` (radical `= F_{2^{gcd(2a,m)}}`;
  the often-quoted `m - 2·gcd(a,m)` agrees only when `m/gcd(a,m)` is even — true
  for these power-of-two `m`, not in general).
- `experiments/gold_form_from_games.py`: rebuilds the form using literal
  Turning-Corners products on small fields.
- `experiments/tartan_bilinear.py`: rebuilds the polar form using game products.

The useful contribution is not that Gold forms are new. They are not. The repo
makes them into a concrete bridge object between:

- the nimber field of impartial game values,
- characteristic-2 quadratic form theory,
- Arf classification/counting,
- and candidate P-set game semantics.

## Broadening the form: a game-realizable quadratic trace family

The thread above fixes one form (Gold, coefficient 1) and hunts for a game. The
form side is much larger. A standard trace representation (e.g. Carlet;
arXiv:1305.3700) describes the quadratic part of Boolean functions on
`F_{2^m}` using trace monomials of the form

```text
Q_c(x) = Σ_{i=1}^{m/2-1} Tr_1^m(c_i · x^{1+2^i})   [ + a half-trace middle term ]
```

with `c_i ∈ F_{2^m}`, plus the usual affine terms and, in even degree, the
middle/half-trace term. Each displayed monomial
`c_i·x^{1+2^i} = c_i ⊗ x ⊗ x^{2^i}` is built from nim-products of `x` with its
`i`-fold Frobenius image. Thus a large trace-presented quadratic family is
available from coin-turning operations, not just the coefficient-1 Gold atom.
The current probe intentionally omits the middle term and affine bookkeeping;
that omission is harmless for the bent examples it tests, but it should not be
phrased as a complete implementation of every Boolean quadratic form.

Why this matters for the open question: `{Q=0}` is farthest from any XOR-subspace
exactly when `Q` is **bent** (nondegenerate polar form, rank `m`, trivial radical,
`m` even) — the maximal-nonlinearity case, hardest for a normal-play sum and the
cleanest Tier-2 target. **Implemented-and-tested** (`experiments/gold_family_survey.py`,
exhaustive over `F_256`):

- The **unscaled** Gold form `Tr(x^{1+2^a})` is not bent in the tested
  power-of-two fields — radical
  `F_{2^{gcd(2a,m)}}`, dim ≥ 1, rank `m − gcd(2a,m)`.
- But its **components** `Tr(λ·x^{1+2^a})` **are bent for 2/3 of `λ`** when
  `gcd(a,m)=1` (APN exponent) — exactly the classical count `2(2^m-1)/3` of bent
  components of a Gold power map, reproduced over `F_256` (170/255 for `a=1,3`).
  For `gcd(a,m)>1` (non-APN, e.g. `a=2` on `F_256`) the split differs (204 bent,
  51 of rank 4). **A single extra nim-multiplication — the coefficient `λ` —
  already unlocks nondegenerate game-realizable forms in these cases.** Random
  multi-term samples in the same script also produce bent examples.
- Bent witnesses validate the zero-count `#{Q=0} = 2^{m-1} + (−1)^Arf·2^{m/2-1}`
  exactly. **Observation (not yet a theorem):** all 170 bent components of
  `Tr(λ·x^{1+2})` over `F_256` carry **Arf 0** — single-component broadening
  reaches bent at only one win-bias sign; Arf-1 bent forms appear to need sums.

The route consequence sharpens the current probes. On a bent form `R(B) = {0}`, so
the symmetric-`B` loopy rule (Loss-set `= R(B)`, see `loopy_quadric.rs`) collapses
to `Loss = {0}` — the radical route is empty — and the frame-blind `Sp(B)` no-go
applies without a degenerate radical layer. Bent game-realizable forms are
therefore a cleaner Tier-2 test than the degenerate Gold examples: the `(m,a) =
(4,1)` radical coincidence from `loopy_quadric.rs` cannot recur for a bent form.

## Arf as a conditional win-bias

For a nonsingular quadratic form on `F_2^(2r)`:

```text
#{x : Q(x)=0} = 2^(2r-1) + (-1)^Arf * 2^(r-1)
```

For degenerate forms, the implementation uses the standard radical-adjusted
count: an anisotropic radical balances the values exactly; an isotropic radical
scales the bias.

`experiments/arf_win_bias.py` brute-forces the value distribution of the Gold
forms and matches the Arf-predicted zero counts.

Interpretation:

If a game had P-positions exactly `{x : Q(x)=0}`, then Arf would say which
player wins from more starting positions and by what square-root-scale margin.

This is meaningful, but conditional. It does not by itself exhibit such a game.

## Why normal play does not solve it

For a normal-play disjunctive sum of impartial games, the P-condition is:

```text
g_1 xor g_2 xor ... xor g_n = 0
```

So the P-set is linear in Grundy coordinates. A genuine quadratic zero set is not
linear.

For a quadratic form in characteristic 2:

```text
Q(u+v) = Q(u) + Q(v) + B(u,v)
```

If `Q(u)=Q(v)=0`, then:

```text
u+v in {Q=0}  iff  B(u,v)=0
```

So the polar form is exactly the obstruction to XOR-closure. The polar form is
already game-built from nim products, but the quadratic refinement is the missing
piece.

## The current probes

### Quadratic-set fitting

`fit_f2_quadratic` asks whether a subset `S <= F_2^k` is the zero set of some
quadratic polynomial

```text
c + sum q_i x_i + sum b_ij x_i x_j.
```

It returns whether the result has genuine quadratic content or is just affine
linear. This is the test bench for candidate games.

### Misere route

The misere quotient code is promising only because misere sums are not
XOR-linear. The current bounded results are negative:

- `star` gives `Z/2`, with a rank-0 linear P-set.
- small misere Nim quotients are not elementary 2-groups.
- the octal sweep over 292 codes, heap cutoff <= 4, found quotient orders
  `2, 6, 10, 12, 14`, no `(Z/2)^k` for `k >= 2`.

**Why — the kernel obstruction (a structural obstruction, not an empty search).** The
empirical negatives are explained by the Plambeck-Siegel structure theory
("Misere quotients for impartial games", JCTA 2008), and it blocks the current
misere-quotient route to *genuine* quadrics:

- Every finite misere quotient `Q` has a **kernel** `K` — the mutual-divisibility
  class of the product of all idempotents — which is the **maximal subgroup** of
  `Q`; the map `x -> zx` (`z` = kernel identity) surjects `Q ->> K`, and every
  homomorphism from `Q` to a group factors through it. So `K` is the canonical
  group shadow available inside this quotient formalism. (Tame
  `T_n = K_n ∪ {1,a}`, `K_n ~= (Z/2)^n`, `|T_n| = 2^n+2`
  — the genuine `(Z/2)^n` is the *kernel*, never the whole quotient, which is why
  the example's "is the quotient `(Z/2)^k`" test never fires for `k >= 2`.)
- **Theorem 6.4:** `z*Phi(G) = z*Phi(H) <=> G,H have the same *normal-play* Grundy
  value`. So `K` is (isomorphic to) the **normal-play nim-value group** `(Z/2)^k`
  under XOR, and `P ∩ K` is the normal-play `{Grundy = 0}` set — **XOR-linear**.

A genuine quadric is a nonlinear zero set on a vector space. In the standard
misere-quotient structure, the canonical vector-space candidate is `K`, and
there the P-structure is the linear normal-play one (Thm 6.4). The genuinely
misere non-linear behavior appears off the kernel, among non-group "fickle
units", where this vector-space/quadric framing no longer applies. This is a
real obstruction to the quotient-kernel route, not merely an empty search.
`experiments/misere_kernel.py` verifies the reading concretely on `R8`, the
smallest wild quotient: kernel `(Z/2)^2`, `P ∩ K = {0}` (linear), and the lone
genuine P-element outside `K`. Caveat: Thm 6.4 has a regularity hypothesis on
the closed game set; the script checks a regular finite quotient, not every
conceivable misere construction.

### Interactive route

`examples/interactive_kernel.rs` confirms three useful facts:

- Any subset can be made the P-set of some ad hoc acyclic game, so existence is
  trivial.
- A rule that directly references `Q` can reproduce `{Q=0}`, but that is
  tautological.
- Rules using the polar form `B` in the tested ways do not reproduce the Gold
  zero set.

### Loopy route

`examples/loopy_quadric.rs` (on `games/loopy.rs`, which lifts `kernel::outcomes`'
Win/Loss/**Draw** retrograde analysis to a first-class Loss-set *and* Draw-set)
adds the third non-normal-play escape beside interactive and misère. A cyclic move
graph has a Draw outcome — a position from which neither player forces a win — and
the Draw-set is not bound by the XOR-linearity that blocks normal-play sums, so it
is a new place to look for `{Q = 0}`.

What the symmetric B-coupling rule (move `v → v ⊕ d` whenever `B(v,d) = 1`) actually
produces is instructive: since `B` is symmetric the move graph is *undirected*, so
the only Losses are isolated vertices, and `v` is isolated exactly when `B(v,·) ≡ 0`
— i.e. `v ∈ R(B)`. So **Loss-set = R(B)** (the radical) regardless of `Q`. At
`(m,a) = (4,1)` this coincidentally equals `{Q=0}` (both 4 points), which breaks at
`m = 8` (`|R(B)| = 4` vs `|{Q=0}| = 112`). And `R(B)` is precisely the degenerate
part on which the Tier-1 `Sp(B)` no-go is *silent* (the no-go constrains the
nondegenerate core `V/R(B)`). So the loopy B-only rule reproduces the obstruction
from a new angle rather than escaping it. The instrument — a cyclic rule's Draw-set
fed through `fit_f2_quadratic` — is what's new; a genuine Tier-2 witness must hit
`{Q=0}` where it is not the radical.

### Bent route

`examples/bent_route.rs` runs the route probes on a **bent** game-realizable form
— a bent Gold component `Tr(λ x^{1+2^a})` (bent for 2/3 of `λ`, see "Broadening the
form"). Bent is the cleanest Tier-2 test: `R(B) = {0}`, so the symmetric-B loopy
Loss-set collapses to `{0}` (radical route dead, no `(m,a)=(4,1)`-style coincidence
possible) and the `Sp(B)` no-go applies in full. Reading the form as an **Ising
energy** `Q(v) = Σ_{i<j} B_ij v_i v_j + Σ_i q_i v_i` (couplings `B` + per-coin field
`q_i = Q(e_i)`, both game-realizable), two results stand out (`m=8`, `λ=2`, bent
Arf 0):

- **B + frame reaches the right quadric *class* in this probe.** A single-bit rule gated
  by `B` alone in the bit frame (no diagonal, no `Q`) produces a genuine **bent
  quadric of the correct Arf** — but a *different* member of the isometry class
  (agreement with `{Q=0}` exactly at chance, `128/256`). So the residual gap to the
  *specific* Gold `{Q=0}` is alignment within the `O(Q)`-orbit — i.e. the diagonal
  framing, sharpened to a nondegenerate form with no radical to muddy it.
- **The naive Ising completion fails.** Adding the per-coin field `q_i` as a local
  spin-flip gate (`ΔQ_i(v) = q_i ⊕ B(v,e_i)`) does **not** align `B`'s quadric to
  `{Q=0}`; it leaves the quadric variety entirely (P-set not a quadric). So the
  diagonal framing must enter some way *other than* a per-coin spin-flip gate — a
  concrete negative that narrows the search.

Net: on this clean bent case, a `B` plus coordinate-frame rule can land in the
right quadric class, but not on the specific Gold zero set. Aligning to that
specific quadric (the diagonal framing's naturality) remains open, and the
tested local-field/spin-flip assembly is a concrete negative example.

### Frame-blind no-go

For dimension at least 4, if a finite game on `V = F_2^(2r)` has a move relation
invariant under the full symplectic group `Sp(B)`, then the P-set is a union of
`Sp(B)`-orbits. Since `Sp(B)` is transitive on `V \ {0}`, the invariant subsets
are only:

```text
empty, {0}, V\{0}, V.
```

Those are not nondegenerate quadrics in dimension >= 4. This is a real
obstruction to frame-blind `B`-only rules.

Caveats:

- Dimension 2 has an exception: the anisotropic quadric has zero set `{0}`.
- The no-go does not explain coordinate/frame-dependent negative probes. Those
  have already broken `Sp(B)` symmetry.

## The frame and diagonal story

Once a coordinate frame is allowed, the quadratic refinement

```text
Q_frame(v) = sum_{i<j} B(e_i,e_j) v_i v_j
```

has polar form `B`. The experiment
`experiments/framing_obstruction.py` shows that, for the Gold polar forms tested,
the Gold form decomposes as:

```text
Q_gold = Q_frame + ell_diag
ell_diag(v) = sum_i Q_gold(e_i) v_i.
```

For the tested genuinely quadratic Gold forms up to `F_{2^16}`,
`Q_frame` has Arf 0 and the diagonal term flips to the Gold Arf 1 form.

This is a good way to state the remaining problem:

> Is the diagonal framing `q_i = Q_gold(e_i)` itself game-natural?

That is sharper than asking vaguely for a "quadratic refinement from games".

Caveat: do not state that every `B + frame` construction is split. Random
nondegenerate alternating forms can give `Q_frame` of Arf 1. The supported claim
is about the Gold polar forms tested here.

## The naturality dichotomy: the open problem depends on a definition

The probes (`open_question_probe.py`, `interactive_kernel.rs`,
`framing_obstruction.py`) keep hitting the same wall from different angles, and
together they put pressure on the word "natural": frame-blind rules are too
symmetrical, while constructions that directly evaluate `Q` are tautological.
What separates a legitimate game rule from an evaluator is now part of the
problem, not background prose.

Organize candidate games by the symmetry group of their **encoding** — the map
`x ↦ (initial configuration)` — as a subgroup `G ≤ GL(V)` under which the move
relation is equivariant. Three tiers:

**Tier 1 — frame-blind, `G ⊇ Sp(B)` — provably NO.** *(Theorem.)* If the move
relation is invariant under the full symplectic group of the polar form, its P-set
is a union of `Sp(B)`-orbits; `Sp(B)` is transitive on `V∖{0}`, so the only
invariant sets are `∅, {0}, V∖{0}, V` — no nondegenerate quadric in dim ≥ 4 (the
"Frame-blind no-go" above). Caveat (degeneracy): the Gold `B` has radical
`R(B) = F_{2^{gcd(2a,m)}}`, so this literally constrains only the nondegenerate
core `V/R(B)`; on `R(B)` the form is the linear `ℓ_diag` and the no-go is silent.

**Tier 3 — `x`-evaluator circuit, `⟨Frobenius⟩ ⊆ O(Q)` — YES, but tautological.**
*(Implemented-and-tested.)* Choosing the refinement `Q` drops the admissible
symmetry to `O(Q) ⊊ Sp(B)`, which is *not* transitive (it preserves `Q`), and
`{Q=0}` is an `O(Q)`-orbit union — so the Tier-1 engine is gone. Concretely
`Q_a(x) = Tr(x ⊗ x^{2^a}) = ⊕_i (x^{2^i} ⊗ x^{2^{i+a}})` is a fixed **circuit of
game operations** on `x` (`m−1` Frobenius squarings, `m` Turning-Corners products,
an XOR fold — `gold_form_from_games.py::gold_literal`, verified over `F_4, F_{16}`
against the algebraic product). Realized as the disjunctive sum of those `m`
Turning-Corners subgames with inputs driven by `x`, its P-condition is exactly
`{Q_a=0}`. The circuit is **Frobenius-symmetric** — the `m` summands form one
Galois orbit, so `x ↦ x²` merely permutes them — hence the encoding is
`⟨Frobenius⟩`-equivariant, and `⟨Frobenius⟩ ⊆ O(Q_a)` because Frobenius is itself
an `F_2`-linear isometry of `Q_a` (`Q_a(x²)=Q_a(x)`). So this is a *Galois-natural*
evaluator, not an arbitrary lookup table. What keeps it tautological: the inputs
are *driven by `x`* — the form's structure is fed in, not produced by autonomous
play (the same gap `open_question_probe.py` flags: "a rule that directly references
`Q` is tautological").

**Tier 2 — the genuine open core.** *(Open.)* Between the two: a *single
fixed-rule* game, positions indexed by field elements, whose **single-position**
Grundy-zero (or interactive-kernel) set is `{Q_a=0}` with **no per-`x`
scaffolding**. `open_question_probe.py` localizes the one missing ingredient — the
linear part is Grundy/XOR (game-realizable), the XOR-closure obstruction is exactly
`B` (game-realizable via coin-turning products), and what remains is a *play rule*
that reads the bilinear coupling `B` out as the quadratic outcome `Q`, necessarily
interactive or misère (normal-play sums give XOR-linear subspace P-sets).

So the honest open problem is a **dichotomy with a definitional gap**: frame-blind
rules provably cannot; per-`x` Galois-natural evaluators can; the question is
whether anything *in the gap* — a fixed-rule game more constrained than an
evaluator but not frame-blind — realizes the quadric. Resolving it requires first
*defining* the gap (an encoding-complexity / equivariance condition admitting the
Frobenius-diagonal symmetry but forbidding per-`x` lookup), at which point the
question becomes sharper than the current taste-level word "natural".

## What should be in the writeup

The draft paper should stay narrow:

1. Explain why arbitrary games are not a scalar ring.
2. Explain why char-2 Clifford needs independent `q` and `b`.
3. Build the Gold forms from nim/game operations.
4. Validate ranks and zero-counts.
5. State the conditional win-bias interpretation.
6. State the open problem as the naturality dichotomy: frame-blind (`Sp(B)`)
   no-go, Galois-natural `x`-evaluator yes, the fixed-rule middle open.

Do not make the paper a catalogue of every module. The odd-characteristic,
p-adic, Witt, Brauer-Wall, CGA, spinor, Hopf, and transfinite-ordinal modules are
useful infrastructure, but they belong in a separate implementation appendix or
README section unless they are needed for the Arf/game thread. The same goes for
the implementation round-out layer added later — `Qp` (p-adic field) + p-adic
Springer, surreal lazy inversion / real roots / Gonshor transfinite birthdays,
the `Fpn` Galois toolkit, field invariants (level/u/Pythagoras), Hermitian forms
over surcomplex, the Cayley bivector↔rotor transform + general multivector
inverse, and atomic weight (`atomic_weight.rs`, finishing thermography) — all
appendix material, none of it changes the Arf/game claims. The same goes for the
symmetry round-out that squared the "any number" table: the `Laurent<S,K>`
transcendental functor (`F_q((t))`, the equal-char local field) and the `Qq`
unramified field `Frac(W_N(F_q))`, the `HasFractionField`/`HasRingOfIntegers`
trait pair making the (field, ring-of-integers) pairing structural, the third
Springer sibling over `F_q((t))` (`springer_laurent.rs`), the divided-power
algebra `Γ` (`divided_power.rs`, the char-faithful symmetric mirror of the
exterior Hopf algebra), and the transfinite surreal↔game round trip via the sign
expansion (`number_game.rs` / `from_transfinite_sign_expansion`) — appendix
material too; none of it touches the Arf/game thread. The same is true of the
new adelic/global layer (`scalar/global`, `forms/local_global/adelic.rs`): useful context for
local-global form experiments, but not evidence for the game-semantics claim. Its
equal-characteristic mirror — the exact global function field `F_q(t)`
(`scalar/global/function_field.rs` = `RationalFunction` over the shared `Poly`
ring `F_q[t]`) and its local–global form layer (`forms/local_global/function_field.rs`: tame
Hilbert symbol, reciprocity `∏_v (a,b)_v = +1`, Hasse–Minkowski, cross-checked
against `springer_laurent`) — is the same kind of "any number" table round-out
(it completes local↔global into a 2-row char-0/char-`p` table and fills the
function-field row's `F_q[t]⊂F_q(t)` pairing), and likewise touches no Arf/game claim.
The analytic-layer unification (`scalar/analytic.rs`) is the same kind of
round-out: the scattered root-taking methods are now two traits — `ExactRoots`
(exact `sqrt`/`is_square` across ℚ, the nim/finite fields, the p-adics, `Surreal`,
`Laurent`) and `SeriesRoots` (the lazy `*_to_terms` surreal primitives) — with the
algebraic-closure surcomplex `√(a+bi)` falling out as a blanket
`Surcomplex<R: ExactRoots+Ordered>` impl (it used to be a private helper inside
`forms/char0`). It makes surcomplex root-taking first-class and the classifier call
a trait, but changes no Arf/game claim.
The latest symmetry round-out is five more trait/mirror closures of the same kind,
none touching the Arf/game thread: (1) `ResidueField` (`scalar/residue.rs`) now carries
not only `k = 𝒪/𝔪` and the angular component but also the multiplicative
Teichmuller section `τ : k → 𝒪`, so the local-field package `(K,𝒪,𝔪,k,Γ,ϖ)` is typed
all the way through; this is what lets the two discrete Springer siblings collapse
into one generic `springer_decompose_local` (`forms/springer_local.rs`) and gives the
unramified `Q_q` decomposition (residue `F_q`) the same generic path as Laurent. The
surreal Springer stays separate exactly because its value group is divisible.
(2) `NimberGame` (`games/nimber_game.rs`) is the char-2 mirror of `NumberGame` —
transfinite Nim heaps `⋆α` carried by their `Ordinal` Grundy value — closing the
`No ↔ On₂` symmetry at the *games* layer (it had only lived at the scalar layer via
the shared CNF core). (3) `FieldExtension` (`scalar/extension.rs`) gives one
relative-trace/norm interface across the algebraic-closure functor (`Surcomplex`,
deg 2), the finite-field tower (`Fpn/Fp`, delegating to the existing `FiniteField`
machinery), the unramified local extension (`Q_q/Q_p`, via the Witt Frobenius), and
the concrete nim-field (`Nimber/F_2`). (4) `CyclicGaloisExtension` now includes the
unramified `Q_q/Q_p` leg too, using a Teichmuller-lifted residue basis and the
Witt-Frobenius, so `forms::trace_form` sees `Surcomplex`, `Fpn`, `Qq`, and `Nimber`
through one `(basis,σ)` interface. (5) the functor corners carry their
valuation/residue signatures explicitly: `Ramified<S,E>` is `Valued` with the new
uniformizer `π` and residue field `k`, while `Gauss<S>` has residue field `k(tbar)`;
`ExactScalar`/`ExactFieldScalar`/`PrecisionScalar` separately mark exact represented
rings/fields versus capped-relative models. `Ramified` and `Gauss` remain excluded
from finite `FieldExtension`/`CyclicGaloisExtension` honestly (non-Galois/infinite
degree), the same boundary `analytic` draws. All appendix material; none of it
changes the Arf/game claims.

## Thermography as tropical arithmetic

A latent structure made explicit, **not a new theorem**. The thermograph
recursion (`games/thermography.rs`) already computes a tropical (min-plus /
max-plus) semiring without naming it: the option folds `combine(_, take_max)` are
the tropical `⊕` (pointwise `max`/`min`), and cooling shifts values by `±t` —
tropical `⊗` (tropical multiplication is ordinary `+`). The two scaffold walls
genuinely live in **dual** semirings: the left wall is a `(max, +)` ⊕-fold over the
Left options' right walls, the right wall a `(min, +)` ⊕-fold over the Right
options' left walls. That duality is why `scalar/tropical.rs` makes
`Tropical<MaxPlus>` and `Tropical<MinPlus>` *distinct types* (sharing one impl body,
the `Surcomplex<S>`/`Laurent<S,K>` move) — the type system then forbids mixing the
two walls.

- **Standard math**: the correspondence (thermography is tropical) — Berlekamp–
  Conway–Guy *Winning Ways*, Conway *ONAG*, Siegel *Combinatorial Game Theory*.
- **Implemented and tested**: `Semiring` + `Tropical<C>` (a sibling trait, *not* a
  `Scalar` — an idempotent `⊕` has no additive inverse, the same reason the game
  group can't carry Clifford), fuzzed in both conventions
  (`tests/tropical_axioms.rs`); and `games::thermograph_via_tropical`, a parallel
  recursion that routes the folds through the named `Pl::oplus_max`/`oplus_min` and
  is pinned **equal** to the golden `thermograph` (`cargo run --example tropical`).
  Naming, machine-checked — no claim beyond it.

## Ordinal nim-multiplication above ω^ω: the Lenstra–DiMuro excess elements

Claim level: **standard math, primary-source-verified** (J. DiMuro, *On Onp*,
arXiv:1108.0962, Thm 3.1 + Table 1; the On₂ table first appeared in H. W. Lenstra,
*Nim multiplication*, 1978). Verified *before* building, per the discipline that a
research-grade landmark must come from the source, not the implementation under test.

The algebraic closure of `F₂` is the ordinals `< ω^(ω^ω)` under nim-arithmetic. The
field generators are `χ_r` for prime powers `r`; for the `(k+1)`-th prime `u`,
`χ_u = ω^(ω^(k-1))` (so `χ_3 = ω`, `χ_5 = ω^ω`, `χ_7 = ω^(ω²)`, `χ_11 = ω^(ω³)`, …),
and `sup_n χ_{u^n}` is the next prime's generator. Two relation families:

- **odd primes `u` (Kummer):** `(χ_u)^u = α_u` and `(χ_{u^{n+1}})^u = χ_{u^n}`, where
  `α_u` is the *excess* — the smallest ordinal `< χ_u` with no `u`-th root below `χ_u`
  (DiMuro Thm 3.1.4). NOT an Artin–Schreier condition (that was a mischaracterization
  in an earlier draft of `ordinal/mod.rs`, now fixed).
- **the prime 2 (Fermat tower, Artin–Schreier):** `χ_2` is a root of `x²+x+1` (=
  nimber `2`), the `x^p−x−1` case `p = u` (DiMuro Thm 3.1.7 / Cor 3.11). This is the
  `F_{2^{2^k}}` tower the finite `Nimber` backend already implements.

Lenstra's effective rule for the excess: `f(u) = ord₂(u)` (multiplicative order of 2
mod `u` = degree of a primitive `u`-th root of unity), and `α_u = χ_{f(u)} + m` for a
finite excess `m`, with `χ_h = Σ_{r ∈ Q(h)} χ_r` over a unique prime-power set `Q(h)`
(Cor 4.4). Worked subtlety: `χ_{f(11)} = χ_10 = χ_5 = ω^ω` because `d(χ_5) = 20`
(`χ_5` is a 5th root of `4 ∈ F_16`, so degree `4·5`) and `2 | 20` ⇒ `Q(10) = {5}`.

**Verified On₂ table** (`α_u`, `u ≤ 43`; `[·]` denotes ordinary/Cantor exponentiation):

| u | α_u | u | α_u | u | α_u |
|---|---|---|---|---|---|
| 3 | 2 | 13 | ω+4 | 29 | ω^(ω²)+4 |
| 5 | 4 | 17 | 16 | 31 | ω^ω+1 |
| 7 | ω+1 | 19 | ω³+4 | 37 | ω³+4 |
| 11 | ω^ω+1 | 23 | ω^(ω³)+1 | 41 | ω^ω+1 |
| | | | | 43 | ω^(ω²)+1 |

`α_5 = 4` — the only new excess the quintic stage `< ω^(ω²)` needs — was *also*
verified independently from finite-field theory: the 5th powers of `F_16*` are
`F_4* = {1,2,3}` (`2 = 3^⊗5`, `3 = 2^⊗5`), so `4` is the smallest non-5th-power, and
it stays a non-5th-power throughout the 2,3-primary tower.

Implementation status: the general prime-power tower (`scalar/big/ordinal/tower.rs`)
is shipped. A monomial `ω^E` is keyed by `place m ↦ base-p(m) digit vector` of its
exponent `E`, and `ω^{E1} ⊗ ω^{E2}` adds the digit vectors per `(m,k)` and reduces
with the per-prime carries (`χ_{u^{k+1}}^u = χ_{u^k}`, `χ_u^u = α_u`).

- **Stage 1** — *scalar* `α_u` (`α_3=2`, `α_5=4`, `α_17=16`): the level-0 carry stays
  in the coefficient, so the product is one monomial. Closes every ordinal `< ω^(ω²)`
  (primes 3,5 only) plus all higher products that avoid a non-scalar Kummer carry.
- **Stage 2** (shipped) — *non-scalar* `α_u` (`α_7=ω+1`, `α_11=ω^ω+1`, …): the level-0
  carry **branches** the monomial into the sum `α_u`, which is nim-multiplied back in.
  The recursion **descends by place** — every `α_{p(m)}` is built from generators at
  places `< m` (verified from Table 1), bottoming out at `α_3=2` in the finite field,
  so it terminates in depth ≤ the largest place index. Carries the source-verified `α_u`
  for primes `u ≤ 43`.

Verified by a prime-3/prime-5 field-axiom sweep (Stage 1), a prime-7 field-axiom sweep
through the `α_7` branching (Stage 2), and source-derived landmarks (`(ω^ω)^⊗5 = 4`;
`(ω^(ω²))^⊗7 = α_7 = ω+1`; `(ω^(ω²))^⊗9 = ω^(ω²·2+1) + ω^(ω²·2)`, cross-checked two
ways so it also pins associativity through the carry). The α *values* come from the
source, never the engine (a wrong α still yields *an* associative ring, so the sweeps
test engine consistency, not the α — that division is deliberate).

**Operational boundary** (honest): a product is exact iff every Kummer carry it triggers
is at a prime `≤ 43`; a carry needing `α_47` or beyond returns `None`, as does anything
`≥ ω^(ω^ω)`. The `u ≤ 43` cap is a **math wall, not an engineering one** — see below.

### Why the table is hardcoded, and the route past it (future work)

DiMuro Thm 4.6 splits every excess as `α_u = (Σ_{q ∈ Q(f(u))} χ_q) + m`, finite `m`:

- the **shape** is formulaic — `f(u) = ord_u(2)` (verified against all 13 table rows),
  `Q(f(u))` = the prime-power factorization minimized under field containment (keep the
  odd parts; keep the 2-part `2^a` only if no `χ_{r^j}` already supplies degree `2^a` —
  the `d(χ_5)=20` absorption), and each `χ_q` is a known ordinal;
- the **finite excess `m` has no closed form**. Thm 4.6 proves only that it is finite and
  exists; DiMuro computes the actual values numerically (Mathematica). The `19 → m=4`
  outlier is why no clean rule fits (an earlier reverse-engineered rule matched values but
  did not self-cohere).

So two unbuilt routes, both worth exploring, neither attempted here:

- **(b) compute `f(u)` + `Q(f(u))` + the `χ`-sum by formula, hardcode only the small
  excess ints.** Same `u ≤ 43` coverage, but more principled — the magic ordinals become
  *computed* `χ_q`, and the `f`/`Q` computation becomes an independent cross-check against
  DiMuro's own columns. Cost: a recursive `d(χ_{r^j})` degree for `Q`'s 2-absorption.
- **(c) full automation** — also derive `m` by searching for the smallest finite
  correction with no `u`-th root below `χ_u`. Unbounded in `u`, but requires u-th-power
  testing in the transfinite field (research-grade) and is **unverifiable past `u = 43`**
  (no source table beyond it) — so it would ship numbers the "verify, don't claim"
  discipline can't sign off without an independent oracle.

Left at the hardcoded table (route "a") for now: minimal, fully source-pinned, explicit
boundary. (b)/(c) are the way to push the wall — (b) for self-verification, (c) for reach.

## Characteristic-2 local–global: the Artin–Schreier symbol over F_{2^m}(t)

Claim level: **standard math, source-verified** (the symbol + reciprocity, the local
Aravire–Jacob Witt decomposition + rank-by-rank local isotropy, **and** the global
Hasse–Minkowski isotropy over `F_{2^m}(t)` itself), **implemented-and-tested**
(`forms/local_global/function_field_char2.rs`, `forms/char2/field.rs`, `forms/springer_char2.rs`).
The local↔global table now commutes into char 2 end-to-end; it touches no Arf/game
claim. Appendix material.

The odd-`q` local–global layer (`forms/local_global/function_field.rs`) needs odd residue
characteristic: its tame Hilbert symbol uses the multiplicative square class
`u^{(|κ|−1)/2}`, undefined when `|κ| = 2^k`. The char-2 replacement is genuinely a
*different formula*, not the same one at `p = 2`:

- **The symbol is the asymmetric Artin–Schreier symbol** `[a, b)`: `a ∈ F_q(t)`
  additive (mod `℘(K)`, `℘(x)=x²+x`), `b ∈ F_q(t)*` multiplicative. The local
  invariant is the **Schmid formula** `s_v(a,b) = Tr_{κ/F₂}(Res_v(a·dlog b)) ∈ F₂`
  (Serre, *Local Fields* XIV §5; Gille–Szamuely §9.2). `[a,b)` splits at `v` iff
  `s_v = 0`; its reduced-norm form is the 2-Pfister `[1,a] ⊥ b·[1,a]`, isotropic iff
  the algebra splits. Reciprocity `∑_v s_v = 0` is the residue theorem on `P¹` — the
  additive analogue of `∏_v (a,b)_v = +1` — so ramification is even.
- **The square-class capability is additive too.** `FiniteChar2Field` mirrors
  `FiniteOddField` with `artin_schreier_class = Tr_{F_q/F₂}` where the odd one carries
  `is_square_value`: in a finite char-2 field *every* element is a square (Frobenius
  is onto), so the working local datum is the additive `F/℘(F) ≅ F₂`, read by the same
  trace the Arf reduction already uses.
- **The residue engine** is the Hensel-parametrization formula
  `Res_P(g dt) = [u^{m-1}]( B(T(u))·P'(T(u))⁻¹ )` with `T(u)` the local series
  `P(T) = u` over `κ = F_q[t]/(P)`; the place `∞` via `u = 1/t`, `dt = u⁻²du` (the
  char-2 sign vanishes). Verified against source residue oracles (`t/P²→1`, `1/P³→0`,
  `t/P³→1` at `P = t²+t+1` over F₄) **and** the reciprocity sweep (a wrong residue
  desyncs `∑_v s_v`).

**The wild-term finding (it saved a wrong build).** The natural next step — a char-2
Springer decomposition `W_q(F_{2^m}((π))) = W_q(k) ⊕ W_q(k)`, the odd-char story
verbatim — is **false**. In char 2 there is a third, *wild* summand: Aravire–Jacob
(*Quadratic forms over rational function fields in characteristic 2*, Thm 1.3) give
`φ = φ₀ ⊥ ψ ⊥ ⟨π⟩φ₁` with `ψ ∈ R_π = {[1,r] : r ∈ π⁻¹·k[π⁻²]}` (binary forms in
*negative odd* powers of `π`). Witness: `[1, π⁻¹]` has zero `W_q(k)` and zero
`π·W_q(k)` parts yet is anisotropic — `℘(x)` of a pole-order-`n` element has *even*
pole order `2n`, never `π⁻¹`'s odd order 1. The same odd/even obstruction surfaces in
the residue engine: Hermite reduction lowers a pole order `j → j−1` only when `j` is
even, so odd-order poles are irreducible — the differential-level shadow of `R_π`.

**The decomposition, delivered (`forms/springer_char2.rs`).** The three-layer local
decomposition `φ ↦ (φ₀, ψ, φ₁)` at a place of `F_{2^m}(t)` is now implemented: split
each block coefficient by Laurent-exponent parity (`K = K² ⊕ πK²`), apply the AJ
relation `[a,b] ≅ [1, a_ev·b] ⊥ ⟨π⟩[1, a_odd·b]`, and push each `[1,c]` to its
Artin–Schreier normal form (clear even/positive poles via `c_{n/2} += √c_n`, keep the
`κ`-constant Arf bit and the odd negative poles `R_π`). The rank-by-rank **local**
isotropy (`[a,b]` iso ⟺ `ab ∈ ℘(K_v)`; ranks 3/4 via the Part-A symbol `s_v`; `u(K_v)
= 4` ⇒ rank ≥ 5 isotropic) is pinned to ten source-derived oracles (independent of the
engine, cross-checked via Codex) spanning the `φ₀`/`ψ`/`φ₁` coordinates and every
rank-by-rank branch — including the genuinely anisotropic 4-dim class realising `u = 4`
(`[1,1] ⊥ [π, π⁻¹]`). The naive two-layer `W = W_q(k)²` version was correctly avoided.
Mixed one-class radical tails use the odd-dimensional Clifford invariant:
for `⟨c⟩ ⊥ ⊥[a_i,b_i]`, the local invariant is the Brauer sum
`Σ [a_i b_i, c/a_i)`; over `F_{2^m}((π))` its split/non-split value distinguishes
anisotropic dimension `1` from `3`, while the two-class quasilinear tail remains the
universal anisotropic plane.

**Global isotropy, delivered (`is_isotropic_global_char2`).** Hasse–Minkowski over
`F_q(t)` itself, on three source-pinned ingredients past the per-place symbol. (1) A
global `℘`-membership test: `f ∈ ℘(F_q(t))` ⟺ `f ∈ ℘(K_v)` at every place (the
local–global principle for `℘` holds over the rational function field — the map is
into `∏_v K_v/℘`, a *product*), and the only obstructions live at the poles of `f`
and at `∞` (which also reads the leftover constant's `Tr_{F_q/F₂}`), so a finite
sweep decides it. This settles **rank 2** (`[a,b]` iso ⟺ `ab ∈ ℘`) — *not* reducible
to a finite bad-place scan, since the constant-trace obstruction (`[1,1]/F₂(t)`,
anisotropic) sits at infinitely many odd-degree places. (2) The **totally-singular**
part is elementary because `[F_q(t):F_q(t)²] = 2`: `≥ 3` entries are `K²`-dependent
(isotropic), a binary `⟨c₁,c₂⟩` is iso ⟺ `c₁c₂ ∈ K²`, and an anisotropic binary
quasilinear part is *universal* (isotropises any form with a nonzero block) — so there
is **no** local–global failure to punt on here (the failure that does exist for
quasilinear forms is over more general base fields, not rational `F_q(t)`). (3) For
**rank 3/4 non-degenerate**, a finite bad-place Hasse–Minkowski sweep: at a good place
all coefficients are units, the rank-`≥3` reduction is a `>2`-variable form over the
finite residue field, isotropic by Chevalley–Warning and Hensel-liftable, so only the
finite bad set (places dividing a coefficient, plus `∞`) needs checking. `u(F_q(t)) =
4` (`F_q(t)` is a `C₂` field, Tsen–Lang) caps it — every `rank ≥ 5` form is isotropic.
Ten worked global oracles over `F₂(t)`/`F₄(t)` pin it (the `℘`-obstruction at both a
constant-trace and a pole case, the rank-3/4 branches, the `u=4` anisotropic 4-dim
class `[1,1]⊥[t,1/t]`, and the quasilinear cases). One Codex oracle was caught wrong in
review (`[1,1]⊥[t,t]` is **isotropic** — vector `(1,0,1,1)` — not the norm form of
`[1,t)`, which is `[1,1] ⊥ t·[1,1]`); the corrected witness is `[1,1]⊥[t,1/t]`.

References: Serre, *Local Fields* XIV; Gille–Szamuely, *Central Simple Algebras and
Galois Cohomology* §9; Elman–Karpenko–Merkurjev §§7, 14 (regular/non-degenerate and
totally-singular classification); Aravire–Jacob (the char-2 function-field Witt
theory); Csahók–Kutas–Montessinos–Zábrádi (arXiv 2203.04068, the explicit `℘`
minimization + Hasse–Minkowski over `F(t)`). Independent theory pass cross-checked via
Codex (including the oracle-7 correction above).

## Integral lattices and the genus (the arithmetic view of a form)

Claim level: **standard math** (classical lattice theory — Conway–Sloane *SPLAG*),
**implemented-and-tested** (`forms/integral/lattice.rs`; `linalg/integer.rs` SNF/ext-gcd).
Appendix material — it touches no Arf/game claim. This is the **adelic layer's
intended payload**: the local–global machinery (`forms/local_global/padic.rs`,
`forms/local_global/adelic.rs`) classifies forms over a field; the integral
lattice is the object whose *genus* is
"local equivalence at every place," so the genus computation (M3) is where those
p-adic primitives finally land on a ℤ-form rather than a square class.

Staged M1→M4, each shippable; M1 (delivered) is the geometry of one lattice:

- An `IntegralForm` is a symmetric integer Gram matrix `G`. **Exact** invariants:
  determinant (fraction-free Bareiss), even/unimodular, positive-definiteness
  (Sylvester), the discriminant group `L#/L ≅ ⨁ ℤ/dᵢ` (Smith normal form), and the
  **level** `N` = smallest positive integer with `N·G⁻¹` even-integral (the level of
  the theta modular form; `E₈→1`, `A_2→3`, `A_1=⟨2⟩→4`, and `⟨1⟩→2` because `ℤ` is
  odd — *not* 1, a standard subtlety).
- **Positive-definite geometry.** Short vectors by Fincke–Pohst (float LDLᵀ bounds
  the search box; the exact integer norm filters every leaf, so float error can
  neither add nor drop a vector), giving minimum / minimal vectors / kissing number,
  and `|Aut(L)|` by backtracking over the images of a basis (each basis vector maps
  to a lattice vector of the same norm with the same pairwise inner products; every
  complete assignment is an automorphism, so the count is exact). The automorphism
  search is exponential and **bounded by an explicit node budget** — past it the
  answer is an honest `None`, never a silently truncated count (`E₈`'s Weyl group of
  order ~7·10⁸ and the Leech lattice are deliberately out of brute-force reach).
- Oracles pin it to the classical numbers: `A_2`/`A_3`/`D_4`/`E_8` kissing numbers
  6/12/24/**240**, |Aut| 12/48/**1152**, det 3/4/4/**1**, and `Aut(ℤⁿ)=2ⁿ·n!`.

M2 (delivered) is the `A_n`/`D_n`/`E_{6,7,8}` catalogue (`forms/integral/root_lattices.rs`), with
the `E₈` ↔ mod-8 Bott / `brauer_wall` BW(ℝ)=ℤ/8 resonance: `E₈` is the unique rank-8
even unimodular lattice, so the periodicity that runs the real Clifford table and the
densest-packing lattice are the *same* mod-8 fact, read once on each pillar — a
resonance worth naming, not a theorem.

M3 (delivered, `forms/integral/genus.rs`) is the genus: signature + determinant + the per-prime
Conway–Sloane symbol, computed by an exact-rational **p-adic Jordan decomposition**
(odd `p` diagonalizes; `p = 2` peels type-I lines and type-II even planes). This is
the first time the adelic square-class machinery (`forms/local_global/padic.rs`) acts on a *lattice*
rather than a field square class. Claim level: **standard math** (Conway–Sloane
*SPLAG* Ch. 15), **implemented-and-tested**. `are_in_same_genus` is exact for odd
`p`, and at `p = 2` carries the full Conway–Sloane quintuple data `(scale, dim,
det mod 8, type, oddity)` through the Allcock-corrected fine-symbol reduction:
determinant residues are normalised, compartment oddities are fused, and
sign-walking along trains adds `4` to the affected compartment oddity. The Sage
`canonical_2_adic_reduction` examples and randomised `Uᵀ G U` isometry checks pin
the implementation, alongside the `ℤ⁸` (`1₀^{+8}`, type I) vs `E_8` (`1_{II}^{+8}`,
type II) distinction.

M4 (delivered, `forms/integral/mass_formula.rs`) is the Minkowski–Siegel mass and the Leech lattice.
The mass of the even-unimodular genus has the Bernoulli closed form `mass(n) =
|B_{n/2}|/n · ∏_{j=1}^{n/2-1} |B_{2j}|/(4j)`. At `n = 8` it is `1/696729600 =
1/|W(E_8)|` — and since `E_8` is alone in its genus, the formula **recovers** the
automorphism order `|Aut(E_8)|` that the brute-force counter deliberately refuses to
compute (a satisfying cross-check: the analytic mass hands back the group order the
combinatorial search can't reach). `n = 16` and `n = 24` (the Niemeier mass) match the
published values; the i128 rational model reaches exactly to rank 24, the honest
ceiling. The **Leech lattice** `Λ₂₄` is built from the binary Golay code `[24,12,8]`
(spanning set of `√8·Λ₂₄ ⊂ ℤ²⁴`, HNF basis, `Gram = B·Bᵀ/8`) and **validated, not
trusted**: a rank-24 even unimodular lattice with no roots is Leech by Niemeier's
classification, so the test checks `det = 1`, even, and `short_vectors(2)` empty.
`Aut(Λ₂₄) = Co₀` (order `8315553613086720000`) is recorded as a constant — far past
brute-force reach. Claim level: **standard math** (Conway–Sloane *SPLAG* Ch. 4/16,
Wilson Thm 1.4), **implemented-and-tested**; mass values and the Golay generator
cross-checked with an independent Codex pass.

The **Monster** stays a thematic line only: `Λ₂₄ → Co₀ → Co₁` is a quadratic-form
computation, but the Monster is monstrous moonshine (vertex operator algebras, the
Griess algebra, the `j`-function) — a research field, not a lattice invariant.

This closes Arc 4 and the four-arc expansion: tropical thermography, transfinite
ordinal nim-multiplication, the characteristic-2 local–global table, and now the
integral-lattice/genus/Leech chain — each an independent, shipped, tested layer.

## Useful commands

```sh
cargo test
.venv/bin/python experiments/trace_form_arf.py
.venv/bin/python experiments/gold_form_from_games.py
.venv/bin/python experiments/tartan_bilinear.py
.venv/bin/python experiments/arf_win_bias.py
.venv/bin/python experiments/open_question_probe.py
.venv/bin/python experiments/framing_obstruction.py
.venv/bin/python experiments/gold_family_survey.py
.venv/bin/python experiments/misere_kernel.py
cargo run --example misere_quotient
cargo run --example interactive_kernel
cargo run --example loopy_quadric
cargo run --example bent_route
cargo run --release --example octal_hunt
cargo run --example tropical
```

Do not treat the numbers in old notes as current. Re-run the relevant commands
when using a table or claim in the paper.

## References to keep close

- Conway, *On Numbers and Games*: surreal numbers, nimbers.
- Berlekamp-Conway-Guy, *Winning Ways*: coin-turning games and the
  Turning-Corners/nim-product theorem; thermography.
- Siegel, *Combinatorial Game Theory*: temperature theory / thermography (the
  tropical-semiring reading of the scaffold walls).
- Arf, *Untersuchungen uber quadratische Formen...*: quadratic forms in
  characteristic 2.
- Dickson, *Linear Groups*: binary quadratic forms and zero-count bias.
- Ovsienko, *Real Clifford algebras and quadratic forms over F_2*: useful
  analogy and classification bridge, but do not overstate it as a general
  nim-field Clifford classification theorem.
- Lidl-Niederreiter, *Finite Fields*: finite-field trace/Frobenius background
  and Gold-rank checks.
