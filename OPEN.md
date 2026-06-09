# TODO: Genuine Research Problems

This file is intentionally narrow. It lists directions from repo audits and the
Gold/Arf draft that look like genuine new research rather than implementation of
known formulas, standard algorithms, or already-source-pinned theory. Implemented
mathematical facts and maintenance context live in `README.md` and `AGENTS.md`.

## Natural Gold-quadric game rule

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

Concrete progress targets:
- Formalize a naturality criterion: equivariance, locality, encoding complexity,
  basis/framing access, or a combination of these.
- Prove no-go theorems for larger classes than the current frame-blind `Sp(B)`
  obstruction, especially for polar-form-only and low-complexity frame-dependent
  rules.
- Exhibit a fixed uniform rule, more constrained than an arbitrary lookup game,
  whose P-set, Loss-set, Draw-set, or canonical kernel set is a Gold quadric.
- Explain whether the diagonal refinement `q_i = Q(e_i)` has a game-native source,
  or prove that every acceptable source collapses to a split/incorrect refinement.

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

## Quadratic deformation of the game exterior algebra

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

## Ordinal nim multiplication beyond the verified excess table

Push transfinite nim multiplication beyond the source-verified Lenstra-DiMuro
excess table, especially at the first unverified Kummer carry such as `alpha_47`.

What is implemented:
- The algebraic closure of `F_2` is represented by ordinals `< omega^(omega^omega)`
  under nim-arithmetic.
- The prime-power generator tower is implemented in `src/scalar/big/ordinal/tower.rs`.
  Products are exact when every Kummer carry uses a source-verified excess
  `alpha_u` for an odd prime `u <= 43`.
- Stage 1 handles scalar excesses such as `alpha_3 = 2`, `alpha_5 = 4`, and
  `alpha_17 = 16`; Stage 2 handles nonscalar excesses such as `alpha_7 = omega+1`
  by branching the monomial and recursing to lower places.
- The table values are from DiMuro's source table, not inferred from the engine.
  Field-axiom sweeps test engine consistency, not the truth of the table values.

The source-verified rows currently used are:

| u | alpha_u | u | alpha_u | u | alpha_u |
|---|---|---|---|---|---|
| 3 | 2 | 13 | omega+4 | 29 | omega^(omega^2)+4 |
| 5 | 4 | 17 | 16 | 31 | omega^omega+1 |
| 7 | omega+1 | 19 | omega^3+4 | 37 | omega^3+4 |
| 11 | omega^omega+1 | 23 | omega^(omega^3)+1 | 41 | omega^omega+1 |
| | | | | 43 | omega^(omega^2)+1 |

Why this is research:
- Rewriting the current table-driven code to compute the known shape
  `f(u)`, `Q(f(u))`, and the `chi`-sum, while hardcoding only the finite excess
  integer, is a useful implementation improvement but not new reach.
- Extending past the verified table is different. DiMuro's theorem proves that the
  excess has a formulaic transfinite shape plus a finite correction, but the finite
  correction has no closed form in the cited theorem.
- Shipping new values would require an independent oracle, a root-search theorem,
  or a new algorithmic proof. Otherwise the project would be numerology with a
  pleasant API.

Concrete progress targets:
- Implement the principled same-coverage route: compute `f(u) = ord_u(2)`,
  compute `Q(f(u))`, construct the `chi`-sum, and hardcode only the finite excess
  integer. This should independently cross-check the published rows.
- Derive or certify finite excess terms beyond the published table.
- Build a verified `u`-th-power/root-search oracle for the transfinite field.
- Prove enough about the search to avoid merely empirical extensions.
- Decide what evidence is acceptable for shipping `alpha_47` and beyond.

Relevant surfaces:
- `src/scalar/big/ordinal/tower.rs`
- `src/scalar/big/ordinal/mod.rs`
- `src/scalar/AGENTS.md`
- `examples/tour.rs`

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
