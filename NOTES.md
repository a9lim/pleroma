# Draft notes: what is actually novel here

Status: working draft, not a theorem paper. These notes try to separate the
mathematical content from the repo-tour material and from the speculative
game-semantics question.

The short version:

Pleroma is a Clifford-algebra lab over the field-like cores of Conway game
worlds. The genuinely new thread is not "Clifford over all games"; that phrase
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
| `Qp` | the p-adic *field* `Q_p` (capped-relative precision model; Zp's field of fractions) — the empty cell in the "any number" table; feeds the p-adic Springer leg |
| `Ordinal` in `big/ordinal/` | staged transfinite nimbers below `omega^3`, including `F_4(omega) ~= F_64` |

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

The precise novelty is not that Gold forms are new. They are not. The novelty is
that this repo makes the Gold forms into a concrete bridge object between:

- the nimber field of impartial game values,
- characteristic-2 quadratic form theory,
- Arf classification/counting,
- and candidate P-set game semantics.

## Broadening the form: the game-realizable quadratic family

The thread above fixes one form (Gold, coefficient 1) and hunts for a game. The
form side is much larger, and all of it is game-realizable. **Theorem (trace
representation, e.g. Carlet; arXiv:1305.3700):** every quadratic Boolean function
on `F_{2^m}` is

```text
Q_c(x) = Σ_{i=1}^{m/2-1} Tr_1^m(c_i · x^{1+2^i})   [ + a half-trace middle term ]
```

with `c_i ∈ F_{2^m}`. Each monomial `c_i·x^{1+2^i} = c_i ⊗ x ⊗ x^{2^i}` is
nim-products of `x` with its `i`-fold Frobenius image, so **the whole space of
quadratic forms is built from coin-turning operations**, not just the
coefficient-1 Gold atom. Gold `Tr(x^{1+2^a})` is the single-monomial,
coefficient-1 case.

Why this matters for the open question: `{Q=0}` is farthest from any XOR-subspace
exactly when `Q` is **bent** (nondegenerate polar form, rank `m`, trivial radical,
`m` even) — the maximal-nonlinearity case, hardest for a normal-play sum and the
cleanest Tier-2 target. **Implemented-and-tested** (`experiments/gold_family_survey.py`,
exhaustive over `F_256`):

- The **unscaled** Gold form `Tr(x^{1+2^a})` is **never bent** — radical
  `F_{2^{gcd(2a,m)}}`, dim ≥ 1, rank `m − gcd(2a,m)`.
- But its **components** `Tr(λ·x^{1+2^a})` **are bent for 2/3 of `λ`** when
  `gcd(a,m)=1` (APN exponent) — exactly the classical count `2(2^m-1)/3` of bent
  components of a Gold power map, reproduced over `F_256` (170/255 for `a=1,3`).
  For `gcd(a,m)>1` (non-APN, e.g. `a=2` on `F_256`) the split differs (204 bent,
  51 of rank 4). **A single extra nim-multiplication — the coefficient `λ` —
  already unlocks nondegenerate game-realizable forms.** Multi-term sums reach
  bent generically too.
- Bent witnesses validate the zero-count `#{Q=0} = 2^{m-1} + (−1)^Arf·2^{m/2-1}`
  exactly. **Observation (not yet a theorem):** all 170 bent components of
  `Tr(λ·x^{1+2})` over `F_256` carry **Arf 0** — single-component broadening
  reaches bent at only one win-bias sign; Arf-1 bent forms appear to need sums.

The route consequence sharpens the open question. On a bent form `R(B) = {0}`, so
the symmetric-`B` loopy rule (Loss-set `= R(B)`, see `loopy_quadric.rs`) collapses
to `Loss = {0}` — the radical route is empty — and the frame-blind `Sp(B)` no-go
applies in full, with **no degenerate part for it to be silent on**. So bent
game-realizable forms are the *purest* Tier-2 statement: the `(m,a)=(4,1)` radical
coincidence that masqueraded as a hit in `loopy_quadric.rs` provably cannot recur.
The next route probe should feed a bent `{Q=0}` into the interactive/misère
instruments, where realizing it is hardest.

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

**Why — the kernel no-go (a structural obstruction, not an empty search).** The
empirical negatives are explained by the Plambeck-Siegel structure theory
("Misere quotients for impartial games", JCTA 2008), and it closes the misere
route for *genuine* quadrics:

- Every finite misere quotient `Q` has a **kernel** `K` — the mutual-divisibility
  class of the product of all idempotents — which is the **maximal subgroup** of
  `Q`; the map `x -> zx` (`z` = kernel identity) surjects `Q ->> K`, and every
  homomorphism `Q ->` group factors through it. So `K` is the **only** `F_2`-vector-
  space part of `Q`. (Tame `T_n = K_n u {1,a}`, `K_n ~= (Z/2)^n`, `|T_n| = 2^n+2`
  — the genuine `(Z/2)^n` is the *kernel*, never the whole quotient, which is why
  the example's "is the quotient `(Z/2)^k`" test never fires for `k >= 2`.)
- **Theorem 6.4:** `z*Phi(G) = z*Phi(H) <=> G,H have the same *normal-play* Grundy
  value`. So `K` is (isomorphic to) the **normal-play nim-value group** `(Z/2)^k`
  under XOR, and `P n K` is the normal-play `{Grundy = 0}` set — **XOR-linear**.

A genuine quadric is a *nonlinear* zero set *on a vector space*. The only vector-
space part of a misere quotient is `K`, and there the P-structure is the linear
normal-play one (Thm 6.4). The genuine misere non-linearity lives **off** the
kernel, among the non-group "fickle units", where there is no ambient vector space
and "quadric" is undefined. **So the misere quotient places its non-linearity
exactly where the quadratic-form framing cannot reach it** — the misere analog of
the frame-blind `Sp(B)` no-go. Verified concretely on `R8`, the smallest wild
quotient (`experiments/misere_kernel.py`): kernel `(Z/2)^2`, `P n K = {0}` (linear),
the lone genuine P-element a fickle unit outside `K`. (Caveat: Thm 6.4 carries a
regularity hypothesis on the closed game set, satisfied by the regular finite
quotients arising in practice.)

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

- **B + frame already reaches the right quadric *class*.** A single-bit rule gated
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

Net: on the cleanest case, "B + frame" is provably enough for the right quadric
*class*; aligning to the specific Gold quadric (the diagonal framing's naturality)
is the open core, and the obvious local-field assembly is now ruled out.

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

## The naturality dichotomy: the open problem is a definition problem

The probes (`open_question_probe.py`, `interactive_kernel.rs`,
`framing_obstruction.py`) keep hitting the same wall from different angles, and
together they show the wall is **definitional, not mathematical**: every concrete
obstruction to `{Q_a=0}` being a P-set dissolves, and every concrete construction
that reaches it is a tautological evaluator. What separates "tautological" from
"natural" is the only thing doing real work. So the open problem is, precisely, to
*define* that line.

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
question becomes decidable rather than a matter of taste.

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
the "completeness" round-out layer added later — `Qp` (p-adic field) + p-adic
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
material too; none of it touches the Arf/game thread.

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
```

Current checked state during the rewrite pass:

- `cargo test`: 208 passed.
- the listed experiments/examples matched the tables and claims used in the
  rewritten draft.

## References to keep close

- Conway, *On Numbers and Games*: surreal numbers, nimbers.
- Berlekamp-Conway-Guy, *Winning Ways*: coin-turning games and the
  Turning-Corners/nim-product theorem.
- Arf, *Untersuchungen uber quadratische Formen...*: quadratic forms in
  characteristic 2.
- Dickson, *Linear Groups*: binary quadratic forms and zero-count bias.
- Ovsienko, *Real Clifford algebras and quadratic forms over F_2*: useful
  analogy and classification bridge, but do not overstate it as a general
  nim-field Clifford classification theorem.
- Lidl-Niederreiter, *Finite Fields*: finite-field trace/Frobenius background
  and Gold-rank checks.
