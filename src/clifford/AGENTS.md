# AGENTS.md — `src/clifford/`

The PILLAR holding the multivector engine and the geometric-algebra layer on top.
Everything here is generic over `S: Scalar` — the same code runs over nimbers
(char 2), surreals (char 0), surcomplex, integers, etc. Python wrappers (the
`backend!` macro) live in `src/py/engine.rs`.

`mod.rs` is a thin hub re-exporting the engine + versor + the structured-algebra
modules flat, so public paths stay shallow (`clifford::Metric`,
`clifford::coproduct`, … — note `sandwich` is an inherent `CliffordAlgebra` method,
called as `alg.sandwich(…)`, not a free path).

Fixed-width mathematical payloads here are `u128`/`i128`: blade masks,
divided-power exponents, spinor/Dickson parities, and Frobenius subfield data.
`usize` is for dimensions, basis indices, and matrix indexing.

## The engine (`engine.rs` + `engine/`)

`engine.rs` is a thin hub (+ product/regression tests). The associative-algebra
core is split by concept under `engine/`:

- **`basis.rs`** — `bits` / `grade` / `MAX_BASIS_DIM` / `wedge_sign`.
- **`metric.rs`** — `Metric {q, b, a}`, constructors, `direct_sum`, `q_val`/
  `has_upper`, `map` (coefficient base-change `Metric<S>→Metric<T>`, e.g. lifting an
  `F_2` trace form into `Metric<Nimber>` for the Arf classifier, or consuming
  `IntegralForm::clifford_metric*` from the integral-lattice bridge). **The metric
  carries `q` and `b` independently — see the rules.**
- **`product.rs`** — `geom_product_blades` (the general-bilinear Chevalley product)
  plus the `cfg(test)` `reduce_word` oracle it is cross-validated against.
- **`algebra.rs`** — `CliffordAlgebra<S>`: blade arithmetic, grade projection,
  wedge/reverse/graded_tensor/embeddings.
- **`multivector.rs`** — `Multivector<S>`: term store, zero/display helpers.
- **`inverse.rs`** — GENERAL `multivector_inverse` via the shared `linalg::field`
  solver (used when `1+B` is not a versor, e.g. in the Cayley transform).
- **`terms.rs`** — local term-map scale/merge helpers.

## The GA layer

- **`versor.rs`** — the layer on top of the associative core: `versor_inverse`,
  `sandwich`, `twisted_sandwich` (Pin action), `reflect`, left/right_contract,
  dual/undual, grade_involution, norm2, even_part / even_subalgebra. Plus the
  product/involution suite: `clifford_conjugate`, `scalar_product ⟨ab⟩₀`,
  commutator/anticommutator (½-free, char-faithful), the regressive meet `a∨b`. Plus
  the CAYLEY transform `cayley`/`cayley_inverse = (1−B)(1+B)⁻¹`: the exact RATIONAL
  bivector↔rotor map (Lie algebra ↔ Spin group, no cos/sin, char≠2).
- **`blade.rs`** — blade analysis: `blade_subspace {x : x∧A=0}`, `is_blade`,
  `factor_blade` (decompose a blade into grade-1 vectors). Char-faithful.
- **`outermorphism.rs`** — lift a grade-1 `LinearMap<S>` to all grades
  (`f(a∧b)=f(a)∧f(b)`); determinant as the pseudoscalar action `f(I)=det·I`; compose,
  `inverse_outermorphism`. Plus the char poly via exterior powers
  (`exterior_power_trace`, `trace`, `char_poly`). Char-faithful (the char-2
  determinant/permanent too).
- **`frobenius.rs`** — the scalar-Galois ↔ Clifford bridge: turns
  `CyclicGaloisExtension` bases/generators into `LinearMap<E::Base>` values via
  `galois_linear_map` / `frobenius_linear_map`, plus
  `nimber_subfield_frobenius_linear_map` for small represented nimber subfields. Its
  tests pin the outermorphism spectrum (`char_poly`, determinant, exterior traces)
  against Frobenius.
- **`hopf.rs`** — the exterior Hopf algebra: unshuffle coproduct (sign read off
  wedge), counit, antipode = grade involution. Hopf axioms tested over Rational AND
  Nimber.
- **`divided_power.rs`** — the CHAR-FAITHFUL symmetric mirror of `hopf.rs`: the
  divided power algebra Γ(V) (dual of Sym), with a BINOMIAL product and
  DECONCATENATION coproduct. Binomials reduce mod char: `(γᵢ⁽¹⁾)²=2γᵢ⁽²⁾=0` in char 2
  while `γᵢ⁽²⁾≠0` — the honest Γ≠Sym (mirror of exterior `eᵢ²=0`). Standalone (own
  monomials, not the blade engine); Python exposes it via the
  `<World>DividedPowerAlgebra` / `<World>DpVector` backend family.
- **`cga.rs`** — conformal (Cl(n+1,1) null basis: up/down/inner/sphere/plane/meet)
  + projective GA (`pga(n)` = `Cl(n,0,1)`, with the terminating `exp_nilpotent` motor
  exp). Char-0 (needs ½); surreal ∞/ε radii are exact.
- **`spinor.rs`** — concrete left-ideal spinor matrices. Char 0 uses the `∏½(1+w)`
  idempotent search and matches the real-table classifier when it reaches a minimal
  ideal; in char 2 (the branch is keyed on `characteristic()`, so it covers any
  nonsingular field-like char-2 metric, Nimber the main one) a separate no-half path
  takes blade idempotents like `e_i e_j` when they shrink the ideal and otherwise
  keeps the complete left-regular action. `spinor_rep`/`SpinorRep` build the explicit matrix up
  to `MAX_EXPLICIT_SPINOR_DIM`; `lazy_spinor_rep`/`LazySpinorRep` give the sparse,
  unbounded-dimension left-regular action beyond that cap. Clifford relations hold.
- **`spinor_norm.rs`** — the spinor norm `N : O(Q)→F*/F*²` (= norm2 mod squares) +
  the generic `versor_grade_parity` (Dickson; `forms::dickson_of_versor` delegates
  here) + `classify_versor` → `VersorClass<S>` (the spinor-norm + Dickson-parity
  record). Char-2 codomain is `F/℘(F)`.

## Hard rules (clifford-specific)

1. **The metric carries `q` and `b` independently — do not collapse them.**
   `q[i] = eᵢ²` (quadratic form); `b[(i,j)] = {eᵢ,eⱼ}` (polar/anticommutator, i<j).
   In char ≠ 2 they're linked; in char 2 they are NOT — `b` is alternating
   (`b(i,i)=0`) yet `q[i]` can be nonzero. Collapsing to one symmetric bilinear form
   silently makes every char-2 algebra commutative. There is a THIRD, *optional*
   field `a[(i,j)]` (i<j): the in-order/asymmetric contraction that lifts the engine
   to a general (non-symmetric) bilinear form `B` — `e_i e_j = e_i∧e_j + a_{ij}` for
   i<j; `b` stays the symmetric anticommutator. `a` empty ⇒ the ordinary Clifford
   algebra. Build with `Metric::new`, `::diagonal`, `::grassmann`, or
   `::general(q, b, a)`, never the bare struct literal (`a` is keyed i<j only).
2. **Signs go through the scalar's own `neg()`, never a literal `-1` or a
   `characteristic()` branch.** The product emits `S::one().neg()` from the wedge
   antisymmetry. For nimbers `neg` is identity, so `-1 = 1` and char-2 sign-vanishing
   falls out for free. Hardcoding signs breaks char 2.
3. **Verify, don't claim.** The `associativity_*` tests (incl.
   `associativity_general_bilinear_form`) catch product bugs, and
   `general_product_reproduces_reduce_word_when_a_empty` pins the general engine to
   the independent oracle. Add a test before trusting a new operation.
4. **Respect the width contract.** Non-index fixed-width integer payloads are
   `u128`/`i128`. Use `usize` for basis indices and dimensions only.

## Things that look like bugs but are not (clifford layer)

- **Char-2 Clifford over an orthogonal basis is commutative.** `e0*e1 == e1*e0` when
  `b` is empty and the scalar is a nimber. Correct: `{e0,e1}=2B=0` and `-1=1`. Set an
  off-diagonal `b[(i,j)]` to get non-commutativity.
- **`versor_inverse` succeeds iff the spinor norm `v ṽ` is a scalar *and* a monomial**
  (over surreals): `1/(ω+1)` is an infinite Hahn series, so a non-monomial norm has no
  representable inverse. Use `multivector_inverse` (the general linalg solve) when
  `1+B` isn't a versor — that's why the Cayley transform calls it.
- **The `neg_one` branch in `Multivector::display` never fires for nimbers.**
  `neg(one)=one` in char 2, so the `coeff==one` branch catches it first. Harmless.
- **`divided_power.rs` is a standalone parallel algebra, not a layer on the blade
  engine.** Γ(V) is the dual of Sym, not a view of the exterior algebra; Python
  therefore binds it as separate `DividedPowerAlgebra`/`DpVector` classes.
- **`Ordinal` is a Clifford scalar only inside its checked nim-product boundary.**
  `Scalar::mul` panics if a product escapes the source-verified Kummer tower; use
  `Ordinal::nim_mul` before constructing computations that need an explicit `Option`
  boundary.
