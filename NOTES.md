# Notes: Clifford algebras over games, and the Arf thread

Why this project has the shape it does, and where the genuinely open question is.

## The setup

Conway's combinatorial games, under disjunctive sum, form a partially ordered
abelian group — but **not a ring**: the game product is only a congruence on the
*numbers* (the surreals). A Clifford algebra needs a commutative scalar ring, so
"Clifford over games" forces you onto the field-like cores of game-world:

| core | structure | Clifford flavour |
|---|---|---|
| surreals `No` | real-closed field, char 0 | the ℝ-Clifford classification (8-fold periodicity), with infinite/infinitesimal metric entries |
| surcomplex `No[i]` | algebraically closed, char 0 | ℂ-Clifford, 2-fold periodicity |
| nimbers `On₂` | algebraically closed, char **2** | the genuinely different one |

Each is a backend here. The char-2 nimber case is the only one where Clifford
gets a new flavour rather than reproducing ℝ/ℂ with exotic scalars.

## The char-2 / Arf story (solid)

In characteristic 2 the quadratic form `Q` is **not** determined by a symmetric
bilinear form, so the engine carries the squares `q[i]=Q(eᵢ)` and the
alternating polar form `b[(i,j)]={eᵢ,eⱼ}` as independent data.

The classifying invariant of a nonsingular quadratic form over F₂ is the **Arf
invariant** (Arf 1941): for a symplectic basis `{aₖ,bₖ}`,
`Arf(Q) = Σₖ Q(aₖ)Q(bₖ) ∈ F₂`. Two forms are equivalent iff their Arf
invariants agree.

The payoff: Bertram, Kervaire et al. / the survey *"Real Clifford algebras and
quadratic forms over F₂: two old problems become one"* (arXiv:1601.07664) prove
that the classification of (real) Clifford algebras **is** the classification of
F₂ quadratic forms, with the Arf invariant complete. So computing the Arf
invariant of a nim-Clifford form (see `src/arf.rs`, `pl.arf_invariant`) is not a
toy — it returns the isomorphism class of the char-2 Clifford algebra.

`A ⊕ A ≅ H ⊕ H` (two anisotropic planes ≅ two hyperbolic planes) is the
additivity of Arf made executable, and the tool confirms it.

The classifier works over any nim-subfield, not just F₂: a form with entries in
F_{2^{2^k}} is symplectically reduced over that field (pairs normalised with the
`nim_inv` from the versor layer), and the Arf sum is pushed down to F₂ by the
field trace `Tr_{F/F₂}(x) = x + x² + … + x^{2^{m-1}}`, realising the canonical
`k/℘(k) ≅ F₂`. Over F₄, e.g., `q=[*2,*2]` is anisotropic (O⁻) while `q=[*2,*3]`
is hyperbolic (O⁺).

## The games bridge (solid as far as it goes)

Games connect to the nim *field*, concretely:

> **Product Theorem for coin-turning games** (Berlekamp–Conway–Guy, *Winning
> Ways* vol. 3): the Grundy value of a product of coin-turning games (e.g.
> Turning Corners) is the **nim-product** of the factors' Grundy values.

So nim-addition (XOR) and nim-multiplication are literally the arithmetic of a
real class of impartial games. That is the bridge from games to `On₂`, the
scalar field underneath the nimber backend.

## The open question (where it stops)

The bridge above delivers **linear** structure (Grundy values are nim-sums of
single-coin values) and **bilinear** structure (coin-turning products are
nim-products). A quadratic form — the thing carrying an Arf invariant — is, in
char 2, *strictly more* than a bilinear form, and nothing in the standard theory
hands you one from a game:

- the nim-square map `x ↦ x⊗x` is the Frobenius, which is F₂-**linear**, so the
  "diagonal" of the coin-turning bilinear form carries no quadratic content;
- Welter's game has deep structure, but it lands in the representation theory of
  symmetric groups (Sato's conjecture, Irie 2018), not in Arf invariants.

So the genuinely open question this project points at:

> **Is there a natural quadratic refinement of the nim-bilinear form arising
> from a combinatorial game, whose Arf invariant is itself a game invariant?**

This is the char-2 game-theoretic analogue of a *quadratic refinement of a
symplectic form* — the same structure that produces the Arf–Kervaire invariant
from framings in topology. I have not found it in the literature, and I won't
claim it exists.

> **This framing is refined by the three sections below.** The probes show the
> answer is more than "open": the Arf-bearing forms *are* built from game
> operations, and the Arf invariant *is* a win-bias in the counting sense. What
> remains open narrows to one thing — a *natural game* whose P-positions are a
> form's zero set. Read on.

## Empirical probe: quadratic forms intrinsic to the nim-field

`experiments/trace_form_arf.py` runs the first probe — entirely on top of the
shipped library. The natural quadratic forms on a char-2 field are
`Q_a(x) = Tr(x·x^{2^a}) = Tr(x^{1+2^a})` (the Gold functions; `g = Frobenius^a`
is additive, so `Q_a` is genuinely quadratic). We build each over the bit-basis
of `F_{2^m}` and read off its Arf invariant.

Findings:

1. **Validation against known mathematics.** The classifier's polar-form rank
   reproduces the Gold-function rank formula `rank = m − 2·gcd(a,m)` exactly in
   all 15 cases tested (m up to 32) — independent confirmation that the nim
   arithmetic, trace, symplectic reduction, and Arf computation are all correct
   on nontrivial input, not just toy forms.
2. **The substrate carries real quadratic structure.** These forms are
   nondegenerate-of-positive-rank with nonzero Arf — the nim-field has genuine
   Arf-bearing structure beyond the linear (Grundy) and bilinear (coin-turning)
   game operations. Every positive-rank case in the family came out type O⁻.
3. **The bridge is still indirect.** This structure comes from the field's
   Frobenius/trace, not from a game's *move* structure. Connecting a specific
   game to a specific quadratic form remains the open step; the instrument to
   test candidates now exists and is validated.

## Chasing the bridge: the Gold form is game-built

`src/games.rs` implements nim-multiplication a second way — Conway's
Turning-Corners excludant recurrence,
`x ⊗ y = mex{(i⊗y) ⊕ (x⊗j) ⊕ (i⊗j) : i<x, j<y}` — which is the *game*
definition of the product. It agrees with the algebraic Fermat-power `nim_mul`
on every pair tested (x,y < 48). So "coin-turning = nim-multiplication" is now
realised and cross-validated in code, not just cited.

That promotes three operations to *game-realizable*:

- `⊗` nim-product = Turning-Corners Grundy value;
- `□` Frobenius `v ↦ v²` = `v ⊗ v`, the diagonal of Turning Corners;
- `⊕` XOR = disjunctive sum of single-coin positions;

and the trace `Tr(x) = x ⊕ x² ⊕ … ⊕ x^{2^{m-1}}` is iterated `□` and `⊕`. The
Gold form `Q_a(v) = Tr(v ⊗ v^{2^a})` is therefore a **composite of game
operations** on a position's nimber value (and under the 1-D game with
`g(n)=2ⁿ`, a position's value *is* a nimber). `experiments/gold_form_from_games.py`
rebuilds `Q_a` from literal Turning-Corners products and checks it equals the
algebraic form — so the Arf-bearing quadratic form really is *made of games*.

## So where does it actually stop now

The bridge closes at the level of **construction**: the Arf-bearing form is
game-built, and its Arf invariant is computed and validated. What is *not*
established is the **play-semantics**:

> `Q_a` is a derived quadratic quantity, not the Grundy value of a single
> position. Is its Arf invariant the answer to a *game* question — e.g. the
> outcome (or some misère/loopy invariant) of a game canonically attached to
> `Q_a`?

That is the sharp residue: the form is made of games; whether its Arf invariant
*means* something about play is open. The instrument to test a candidate game is
here and validated end to end.

### The Arf invariant *is* a win-bias (counting sense)

Dickson (1901): the Arf invariant of a quadratic form over F₂ is the value the
form takes most often. Quantitatively, for a nonsingular form on F₂^{2m'},
`#{v : Q(v)=0} = 2^{2m'-1} + (-1)^Arf · 2^{m'-1}`. So if a game had P-positions
(second-player wins) exactly `{v : Q(v)=0}`, the Arf invariant would be the
**sign of the win-bias**: which player wins from more starting positions, the
margin a fixed Gauss-sum `2^{m-rank/2-1}` (a square-root-scale fluctuation
around 50/50). `experiments/arf_win_bias.py` brute-forces the value
distribution of the game-built Gold forms and confirms the zero-count matches
the Arf-predicted bias exactly in every case.

So the Arf invariant *is* a win-bias — in the counting sense. The one missing
piece is a **natural game** whose P-positions are `{Q=0}`. Normal-play
disjunctive sums can't do it (their outcomes are XOR-linear, never quadratic),
so a candidate must be *interactive* (e.g. a coin-turning / lights-out style
game coupled through the polar form) or *misère* (where sums are genuinely
non-linear). Constructing or ruling out such a game is the open problem; the
win-bias check above is the target any candidate must hit.

### Sharpening the obstruction: it is exactly the polar form

`experiments/open_question_probe.py` pins down *why* normal play fails and what a
candidate must supply. The P-positions of a disjunctive sum of impartial games
are `{XOR of Grundy = 0}` — a **subspace**. The Gold zero set is a quadric, and in
char 2 the deviation is measured term-for-term by the polar form:

> `Q(u ⊕ v) = Q(u) ⊕ Q(v) ⊕ B(u,v)`, so for `u,v ∈ {Q=0}`:
> `u ⊕ v ∈ {Q=0} ⟺ B(u,v) = 0`.

The probe confirms this exactly: for the nondegenerate Gold forms (e.g. over F₂⁸)
`{Q=0}` is **not** a subspace, and its failure to be XOR-closed is governed
*precisely* by `B`. (Over F₂⁴ the low-rank/degenerate members collapse to a
subspace — there normal play is *not* excluded; the obstruction only bites once
the form is genuinely quadratic.)

This decomposes the problem into three layers, two already game-realizable:

- the **linear** part is Grundy/XOR — Sprague–Grundy;
- the **obstruction** to XOR-closure is exactly `B`, the coin-turning / nim-product
  bilinear form — the Product Theorem / Tartan games;
- the **only** genuinely missing ingredient is a *play rule* that couples
  positions through `B` and reads out the quadratic `Q` (not the bilinear `B`).

So the open question is now sharp and constructive: build (or rule out) an
interactive/misère game whose between-component coupling is the polar form `B`
and whose outcome is `Q`. The two game-realizable layers and the win-bias target
are all in place; what remains is the quadratic *play* rule. Concrete next steps:
misère coin-turning quotients (genuinely non-linear sums) and Tartan-product
couplings, both buildable on the shipped nim-product and the Arf/win-bias
instrument.

### Both next steps, built

*Tartan side — the bilinear layer is game-built.* `games.rs` now carries general
1-D coin-turning games (companion-set encoding, `grundy_1d`) and the 2-D Tartan
product (`tartan_grundy`), with the **Tartan/Product theorem verified**
(`tartan_grundy = nim-product of the component Grundy values`), recovering Turning
Corners as the tartan square of the game with `g(n)=n`. `experiments/tartan_bilinear.py`
then shows the Gold form's polar form `B(e_i,e_j) = Tr(e_i ⊗ e_j^{2^a} ⊕ …)` is
reproduced *exactly* by Turning-Corners products and the trace — so the obstruction
`B` identified above is, concretely, a composite of coin-turning games. Two of the
three layers (linear Grundy, bilinear `B`) are now realized in code from actual
games.

*Misère side — the non-linearity bar is cleared.* `misere.rs` is a memoised
misère-outcome evaluator for any finite impartial game, with misère Nim checked
against Bouton's theorem. The point it nails down: the misère P-set is **not**
`{⊕ = 0}` — `[1]` has nim-sum 1 yet is a P-position, `[1,1]` has nim-sum 0 yet is
an N-position — so it is neither a subspace nor a coset, and the outcome is not an
XOR-linear function of the position. That is exactly the property normal-play sums
*lack* and a quadratic `{Q=0}` P-set *requires*. So misère clears the bar that
ruled out normal play.

*Where it still stops.* Both prerequisites now hold in code — the coupling `B` is
game-built, and misère supplies genuine non-linearity — but neither yet exhibits a
game whose P-set is an actual Gold *quadric* `{Q=0}`. The misère Nim P-set is
non-linear but is not (a priori) a Gold quadric, and the Tartan layer realizes `B`
without yet a play rule that reads out `Q` rather than `B`. The remaining gap is
unchanged in kind — a *quadratic play rule* — but the surrounding scaffolding for
testing candidates (Tartan couplings + misère outcomes + the Arf/win-bias check)
is now all in place.

### Two probes, and a shared test bench

The instrument both probes feed into is `arf::fit_f2_quadratic`: given a set
`S ⊆ F₂^k`, it solves over F₂ for a quadratic form with `{Q=0}=S` (constant +
linear + quadratic terms by Gaussian elimination), returning `None` if `S` is no
quadric, and otherwise the form's Arf — and crucially whether it is *genuinely*
quadratic (nonzero polar-rank) or a mere affine flat. (Sanity: it finds exactly
`2^7` of the `2^8` subsets of F₂³ are quadrics.) So any candidate game's P-set can
be classified: not-a-quadric / linear / genuine quadric-with-Arf.

*Misère route (`misere.rs`, `examples/misere_quotient.rs`).* A bounded
indistinguishability-quotient computer (Plambeck–Siegel) over an abstract
impartial game; verified to give `⋆ ↦ ℤ/2`. Applied to small games it finds:
misère Nim heaps {1,2} has the order-6 quotient (matching the literature), {1,2,3}
likewise small — but these are *not* elementary-abelian 2-groups, so they do not
coordinatise as `F₂^k` and the quadric question doesn't even apply; and the one
that does, `⋆ ↦ ℤ/2`, has a rank-0 (linear) P-set. So no genuine quadric P-set
turns up among the tame games — consistent with tame ≈ linear. A *wild* quotient
of shape `(ℤ/2)^k`, `k ≥ 2`, with Arf-rank ≥ 2 is what would be needed; the
instrument is ready to test one.

*Interactive route (`kernel.rs`, `examples/interactive_kernel.rs`).* A
retrograde Win/Loss/Draw solver for any finite game graph (the P-positions are the
Loss positions). Two findings. (i) *Existence is trivial*: a hand-built acyclic
graph has P-set exactly `{Q=0}` (send every non-zero-of-`Q` to a fixed loss in the
set) — so the open question is never about existence, only about a *natural* rule.
(ii) Searching uniform downward (terminating) rules on `F₂^m`: the rule "move iff
you flip `Q`" reproduces `{Q=0}` exactly — but tautologically, since it references
`Q` itself in the move legality. The rules coupled only through `B` (the
game-legitimate, coin-turning ingredient) do *not*: B-coupled descent gives an
affine subspace, and a single-bit B-gated turn gives a *different* quadric (wrong
Arf). So the open problem reaches its sharpest form yet:

> a game whose moves are built from the combinatorial data (`B` / coin-turning)
> **alone — not from `Q` itself** — whose kernel is the Gold quadric `{Q=0}`.

Referencing `Q` is cheating; referencing only `B` does not (yet) integrate up to
`Q`. The kernel solver + `fit_f2_quadratic` are the bench any candidate runs on.

## The char-0 companion: a matrix-algebra classifier (`classify.rs`)

The Arf invariant returns the isomorphism class of a *char-2* Clifford algebra.
Until now char 0 had the engine but no classifier — an asymmetry. `classify.rs`
closes it: `Cl(p,q)` over a real-closed field follows the 8-fold Bott table
indexed by `s = (q − p) mod 8`, and over an algebraically closed field the 2-fold
table. Because the surreals are real-closed, this *is* the genuine ℝ-Clifford
classification — and the signature is read off the **signs** of the surreal
squares, which may be infinite (ω) or infinitesimal (ε); only the sign matters,
since a real-closed field has square roots of positives (`√ω = ω^{1/2}`).

Cross-checks worth keeping: it reproduces `Cl(1,3) ≅ M₂(ℍ)` but `Cl(3,1) ≅ M₄(ℝ)`
(the two spacetime conventions are genuinely different algebras), `Cl(4,1) ≅
M₄(ℂ)` (conformal GA), and — tying it to the even subalgebra — `Cl(3,0)⁰ ≅
Cl(0,2) ≅ ℍ`. So both characteristics now carry a real classifier: **Arf for
char 2, the matrix-algebra name for char 0.**

## Artin–Schreier ↔ Arf: one trace, two roles

The trace `Tr_{F_{2^m}/F₂}(x) = Σ x^{2^i}` that pushes the Arf invariant down to
F₂ (the canonical `k/℘(k) ≅ F₂`) is the *same* trace that obstructs the
Artin–Schreier equation `y² + y = c`: it is solvable iff `Tr(c) = 0`. So the two
halves of this repo — the Arf classifier and the field arithmetic of `On₂` — run
on one object. `nimber.rs` now exposes it directly: `nim_sqrt` (the inverse
Frobenius, `x^{2^{63}}` in F_{2^64}, always defined in char 2), `nim_trace`, and
`nim_solve_artin_schreier` (an exact F₂ linear solve, solvable exactly on the
trace-zero hyperplane — half the field).

This sharpens "the Arf invariant is a win-bias". Dickson's zero-count
`#{Q=0} = 2^{2m−1} + (−1)^Arf·2^{m−1}` is, term by term, counting how often the
form's value is Artin–Schreier-solvable; the win-bias sign *is* the trace
obstruction aggregated over the form. The field-level operation behind the bias
is now in the library, not just implied.

## Dickson: classifying `O(Q)`, not the form (`arf.rs`)

In char 2 the determinant of any `g ∈ O(Q)` is forced to 1, so it cannot tell a
rotation from a reflection. The **Dickson invariant** `D(g) = rank(g − I) mod 2`
is the replacement, with `SO(Q) = ker D`: a single reflection has `D = 1`, a
product of k reflections `D = k mod 2`. `dickson_matrix` computes it over any
nim-field; `dickson_of_versor` reads it off a Clifford versor as its grade
parity. This is the companion to Arf on the *other* side of the same geometry:
**Arf classifies the form, Dickson classifies the form's isometries.**

## The Witt group makes additivity a group law (`witt.rs`)

`A ⊕ A ≅ H ⊕ H` was checked pointwise via the Arf invariant. The Witt group
`W_q(F)` of nonsingular quadratic forms mod hyperbolics is the home of that fact:
over a finite nim-field it is `≅ ℤ/2`, classified completely by Arf, with `⊥` as
the group operation and the hyperbolic plane as identity. `WittClass` makes the
additivity a one-liner: `w(A) + w(A) = 0` *is* `A ⊕ A ≅ H ⊕ H`.

## General bilinear form: deforming the product (`clifford.rs`)

The engine now computes `Cl(V, B)` for an *arbitrary* (not necessarily symmetric)
bilinear form `B`, via the Chevalley product `e_i e_j = e_i∧e_j + B(e_i,e_j)` in
the wedge basis. `B` is stored factored as `(q, b, a)`: diagonal `q`, symmetric
polar `b` (the anticommutator), and the new strictly-upper / in-order contraction
`a`. With `a` empty this is the ordinary Clifford algebra (and the general
Chevalley product is cross-validated, blade for blade, against the original
swap-reduction now kept as a `#[cfg(test)]` oracle).

Honest scope: the antisymmetric part of `B` is a *gauge* — `Cl(V, B)` is
isomorphic as an algebra to `Cl(V, sym B)`, so `a` does not create new algebras;
it deforms the *product* and the identification between the geometric and
exterior structures (the quantum-Clifford / normal-ordering setting of
Fauser–Oziewicz, interpolating toward the Weyl side). It is the right amount of
generality, faithfully implemented, not a claim of new isomorphism classes.

## The exterior algebra of the game group (`partizan.rs`)

A Clifford algebra needs a commutative scalar *ring*, which is exactly why this
project only reaches the three field-like cores. An **exterior algebra** needs
only a commutative ring of *coefficients* (ℤ) and a *module* of generators — and
the group of partizan games under disjunctive sum is a ℤ-module. So `Λ(game
group)` is well defined on **all** of game-world, the one Clifford-adjacent
structure that does not require the (nonexistent) game product.

`partizan.rs` ships a small short-game engine (sum, negation, the recursive
order, birthday, the number test) and the bridge `Λ¹ → (game group)`,
`e_i ↦ g_i`, built on the shipped Grassmann engine over the new ℤ scalar. The
point is the generators may be **non-numbers** (`⋆`, `↑`) — precisely where the
Conway product, and hence the entire Clifford story, is undefined — yet the wedge
structure (antisymmetry, grading) is perfectly well defined on them. The
2-torsion of `⋆` even surfaces as a relation: `value(2·e_⋆) = ⋆ + ⋆ = 0`. This is
the structural answer to "what lives on the whole game group, not just its
numbers."

## References

- C. Arf, *Untersuchungen über quadratische Formen in Körpern der
  Charakteristik 2* (1941).
- Berlekamp, Conway, Guy, *Winning Ways for Your Mathematical Plays*, vol. 3
  (coin-turning games; Turning Turtles / Turning Corners; the Product Theorem).
- Conway, *On Numbers and Games* (the surreal/nimber fields; `On₂`).
- "Real Clifford algebras and quadratic forms over F₂: two old problems become
  one", arXiv:1601.07664.
- Y. Irie, *p-Saturations of Welter's Game and the Irreducible Representations
  of Symmetric Groups* (2018).
- P. Lounesto, *Clifford Algebras and Spinors* (2nd ed.), Table 16.4 — the
  `Cl(p,q)` classification by `(q−p) mod 8` used in `classify.rs`.
- C. Chevalley, *The Algebraic Theory of Spinors* (1954); B. Fauser & Z.
  Oziewicz, *Clifford Hopf gebra for two-dimensional space* / "Clifford algebra
  of an arbitrary bilinear form" — associativity of the deformed product.
- L. E. Dickson, *Linear Groups* (1901) — the Dickson invariant; the value a
  binary quadratic form takes most often (the Arf win-bias).
- E. Artin, *Geometric Algebra* (1957) — `SO(Q) = ker D` in characteristic 2.
