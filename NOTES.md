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
claim it exists. What this repo provides is the tooling to explore it
computationally: build a candidate form over the nimber backend and read off its
Arf invariant and orthogonal type.

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
