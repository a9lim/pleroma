# ROADMAP — cross-domain connections

This file is the *ambition* document: cross-pillar bridges worth building before
or shortly after the first public release. It is deliberately distinct from
`OPEN.md`:

- **`OPEN.md`** holds *genuine research problems* — things with no known answer
  (the natural Gold-quadric game rule, a game-native quadratic deformation of
  `GameExterior`, transfinite nim excesses past the verified table, and the
  transfinite Arf/Witt question for ordinal-nimber coefficients).
- **`ROADMAP.md`** (this file) holds *buildable bridges* — connections between the
  four mature pillars whose mathematics is largely standard. It now has a
  **built first wave** (Bridges A–D), a **partly built second wave** (Bridges
  E/H/I implemented, Bridge F still proposed), and the deferred Bridge G note.
  This document keeps the mathematical contract, the implemented or proposed
  surfaces, and the remaining honest boundaries in one place. Where a bridge
  brushes against an open question, it says so and points back to `OPEN.md`.

Use the project's claim-level discipline (`AGENTS.md` → "Claim levels and
non-claims") when these land: label each piece **standard math** / **implemented
and tested** / **interpretation** / **open**.

## Why these four

The four pillars currently connect like this:

```
            scalar ───coefficients──── clifford
              │  ╲                        │
        Hackenbush╲  trace_form/Gold      │ classifies
        Turning-  ╲      (forms)          │
         Corners   ╲        │             │
              │     ╲       │             │
            games ──Gold/Arf,──── forms ──┘
                    tropical       │
                    thermography   │
                                integral
```

Before this bridge pass, four edges were conspicuously **missing or partial**:

1. **`integral ↔ clifford` had no computational seam.** The lattice pillar and the
   Clifford engine now meet through `IntegralForm::clifford_metric*` and
   `integral::DiscriminantForm`. → **Bridge A.**
2. **The char-2 classifier spanned only one coefficient field.** It now classifies
   both `Nimber` and supported `Fpn<2,N>` metrics through the Arf façade. →
   **Bridge B.**
3. **`scalar` Galois theory and `clifford` outermorphisms were latent twins.** New
   Frobenius linear-map constructors feed the outermorphism spectral machinery. →
   **Bridge C.**
4. **The `No ↔ On₂` mirror was incomplete at the Clifford layer.** `Ordinal` now
   implements the checked/panic-on-escape `Scalar` surface, so
   `CliffordAlgebra<Ordinal>` builds and tests. → **Bridge D.**

Building the four closes the pillar graph: every pair of pillars that *can* talk
(modulo the game-group-isn't-a-ring constraint) then does.

---

## Bridge A — Lattice ↔ Clifford ↔ Brauer–Wall, via Milgram's Gauss sum

**Pillars:** `forms/integral/` ↔ `clifford/` ↔ `forms/witt/` ↔ `forms/char0`.
**Claim level:** standard math (Milgram/van der Blij; Conway–Sloane) made
computational. The headline bridge — it proves the project's spine crosses pillars.

### The mathematics

For an **even** integral lattice `L` (Gram `G`, so `G[i][i]` even), three objects
now meet in `integral/lattice.rs` and `integral/discriminant.rs`:

- the **signature** `σ = p − q`, computed by exact rational diagonalization,
- the **dual** `L# = G⁻¹L`, using the exact `Rational` inverse already used by `level`,
- the **discriminant group** `A_L = L#/L ≅ ⨁ ℤ/dᵢ`, `|A_L| = |det G|`, exposed
  through invariant factors and represented computationally as `Z^n / GZ^n`.

The bridge datum is the **discriminant quadratic form**

```text
q_L : A_L → ℚ/2ℤ,   q_L(x + L) = xᵀ G x   (mod 2ℤ),   x ∈ L#
b_L : A_L × A_L → ℚ/ℤ,   b_L(x,y) = xᵀ G y   (mod ℤ)
```

well-defined precisely because `L` is even. Its **Gauss sum**

```text
γ(q_L) = |A_L|^(−1/2) · Σ_{x ∈ A_L} exp(π i · q_L(x))
```

is a unit complex number, and **Milgram / van der Blij**:

```text
γ(q_L) = exp(2π i · σ / 8)
```

So the discriminant Gauss-sum **phase is the signature mod 8** — the *same* `ℤ/8`
that `witt/brauer_wall::bw_class_real` computes as the Bott index `(q−p) mod 8`,
that the char-0 8-fold table cycles through, and that makes `E₈` (signature 8 ≡ 0,
trivial `A_L`, `γ = 1`) the rank-8 even unimodular lattice. The bridge turns the
existing prose ("E₈ is where Bott and the lattice world coincide", `root_lattices.rs`)
into a theorem with a computation.

There is a **free internal oracle**: `genus.rs` already computes the `p=2` *oddity*
(trace mod 8), and the Conway–Sloane oddity formula `σ ≡ oddity − Σ_p p-excess
(mod 8)` must agree with the Milgram phase. Two independent routes to `σ mod 8`,
cross-checking each other.

### Implemented surface

- `integral/lattice.rs`
  - `IntegralForm::signature(&self) -> (usize, usize)` diagonalizes `G` over `ℚ`
    and counts signs of the rational pivots, so indefinite lattices are supported.
  - `IntegralForm::clifford_metric(&self) -> Metric<Rational>` — the warm-up rung:
    `q[i] = G[i][i]`, `b[(i,j)] = 2·G[i][j]`. Feeds `CliffordAlgebra<Rational>` and
    `classify_real`. `E₈ → Cl(8,0) → M₁₆(ℝ)`. Also a mod-2 reduction
    `clifford_metric_f2(&self) -> Option<Metric<Nimber>>` for even lattices,
    using `Q/2 mod 2` on the diagonal and `G_ij mod 2` off-diagonal.
- `integral/discriminant.rs`
  - `DiscriminantForm { group, reps, gram_inv }` is built from a nonsingular even
    `IntegralForm` using the standard `A_L ~= Z^n / GZ^n` presentation. The
    representative enumeration uses normalized integer relation rows rather than
    extending Smith normal form with transform matrices.
  - `quadratic_value_mod2`, `bilinear_value_mod1`, `GaussSum::phase_mod8`, and
    `milgram_signature_mod8() -> Option<i128>` make the finite quadratic module
    executable.
  - `verify_milgram(lattice) -> Option<bool>` compares the Gauss-sum phase to the
    exact signature and to the independent Conway-Sloane oddity route in `genus.rs`.

### Oracles / tests

Implemented tests cover `A_n`, `D_4`, `E₈`, `E₈ ⊕ E₈`, odd-lattice rejection, exact
signature on indefinite forms, and the rational / char-2 Clifford metric rungs.
The Milgram phase is checked against the exact signature and genus oddity route.

### Scope / caveats

- The clean Milgram statement is for **even** lattices. Odd (type-I) lattices need
  the oddity-corrected version; ship even-only first, document the boundary, and
  lean on the existing `genus.rs` oddity for the odd case rather than duplicating.
- The Gauss sum is an algebraic number; we compute it in `f64` and verify
  `|γ| = 1` + phase `= σ·45°`. An exact cyclotomic representation is a nice-to-have,
  not required for the check.

---

## Bridge B — the char-2 Arf classifier over the `Fpn<2,N>` fields

**Pillars:** `clifford/` (over `Fpn<2,N>`) ↔ `forms/char2/`.
**Claim level:** implemented-and-tested (standard Arf theory over finite char-2
fields); the *bridge* is new code, the math is classical.

### What landed

`CliffordAlgebra<Fpn<2,3>>` — a Clifford algebra over **F₈** (degree 3, which the
`u128` nimber backend cannot reach: it only holds subfields of 2-power degree) —
now builds **and** classifies. `Nimber` keeps its optimized `nim_trace` path, while
supported `Fpn<2,N>` fields use the same symplectic-reduction algorithm over
generic scalar operations plus the absolute trace.

### Implemented surface

- `char2/arf.rs`
  - `arf_char2<F: FiniteChar2Field>(metric) -> Option<ArfResult>` runs generic
    char-2 symplectic reduction over `Fp<2>` / `Fpn<2,N>`.
  - `arf_fpn_char2<const P, const N>(metric)` is the const-generic façade helper:
    it returns `None` unless `P = 2` and the extension polynomial is supported.
  - `ArfResult::arf` and the Artin-Schreier class are carried as `u128` bits, in
    line with the repo-wide integer-width policy.
- `classify.rs`
  - `Fpn<P,N>` now classifies to `FiniteFieldClass::{Odd, Char2}`, so the same
    monomorphized façade works for odd extensions and characteristic-2 extensions.
  - `WittClassify`, `IsometryClassify`, and `BrauerWallClassify` dispatch to the
    char-2 Arf invariant when `P = 2`.

### Oracles / tests

Implemented tests cross-check `arf_char2` against `arf_f2` when all entries are in
`F₂`, exercise genuine `F₈` coefficients through the absolute trace, verify
additivity over `⊥`, and brute-force the `F₈` zero-count bias for planes.

### Scope / caveats

Honest non-claim (`AGENTS.md`): this is *not* a new classification theorem for all
char-2 Clifford algebras — it computes Arf/BW for the finite `Fpn<2,N>` fields,
the same status the README states for the implemented finite char-2 legs.

---

## Bridge C — Frobenius as an outermorphism

**Pillars:** `scalar/finite_field` (Galois) ↔ `clifford/outermorphism` ↔
`forms/trace_form`.
**Claim level:** implemented-and-tested (the theorems are standard finite-field
theory); the bridge code and the cross-checks are new.

### The mathematics

The Frobenius `σ : F_{p^m} → F_{p^m}, x ↦ x^p` is `F_p`-**linear**. Pick an
`F_p`-basis (the project has them: `FiniteField` / `CyclicGaloisExtension::basis`),
form the matrix `M_σ`, and feed it as a `clifford::LinearMap<Fp<p>>` to the
outermorphism machinery. Then `outermorphism.rs` computes — char-faithfully, no
sign hardcoded — the full spectral suite of `σ`:

- **Characteristic polynomial.** By the normal basis theorem `F_{p^m}` is a free
  `F_p[σ]/(σ^m − 1)`-module of rank 1, so `char_poly(σ) = xᵐ − 1` (over `F₂`,
  `xᵐ + 1`). A clean, exact prediction `char_poly` must reproduce.
- **Vanishing intermediate exterior traces.** Since `xᵐ − 1` has no middle terms,
  the elementary symmetric functions `eₖ(σ) = tr Λᵏσ` satisfy `e₁ = … = e_{m−1} = 0`
  and `e_m = ±1`. Frobenius has a "flat" exterior spectrum — a striking,
  one-line-checkable consequence (`exterior_power_trace(alg, σ, k) == 0` for
  `0 < k < m`).
- **Determinant** `det(σ) = ∏ (m-th roots of unity) = ±1` — the constant term of
  the char poly; verifiable.

### The tie to `trace_form.rs`

`trace_form.rs` builds the **Frobenius-twisted** form `Tr_{E/F}(x · σᵏ(x))` (the
norm form over `Surcomplex`, the Gold form over the nim-fields). The trace itself
is `Tr = 1 + σ + σ² + … + σ^{m−1}` — a *polynomial in the very `σ` this bridge
realizes as a linear map*. So the bridge gives an outermorphism-level reading of
the trace-form construction: lift `σ` to the exterior algebra of `E`-as-`F`-space,
and the `Λᵏ` action organizes the twisted forms across grades. This is a genuine
conceptual link, not just a spectral cross-check.

### Implemented surface

- `clifford/frobenius.rs`
  - `CoordinateCyclicGaloisExtension` extends the cyclic Galois basis with a
    coordinate extractor.
  - `galois_linear_map::<E>(k)` and `frobenius_linear_map::<E>()` build
    `LinearMap<E::Base>` from the chosen basis.
  - `nimber_subfield_frobenius_linear_map(m, k)` gives small exact matrices for
    the represented nimber subfields, avoiding a 128-dimensional exterior-power
    computation when a four- or sixteen-dimensional one is the intended oracle.

Tests pin `char_poly = xᵐ ± 1`, the vanishing middle `Λᵏ`-traces, `det = ±1`, and
composition of Frobenius powers over `Fpn<2,m>`, odd-characteristic `Fpn`, and a
small nimber subfield.

### Scope / caveats

Pure cross-domain wiring + verification; no new theorem. Its value is that it makes
three pillars share one computation and gives `trace_form` a structural home.

---

## Bridge D — transfinite char-2 Clifford (`OrdinalAlgebra`)

**Pillars:** `scalar/big/ordinal` ↔ `clifford/`.
**Claim level:** implemented-and-tested for the checked engine/symmetry completion.
Classification of genuinely transfinite coefficients is still out of scope and
tracked in `OPEN.md`.

### The target and the totality boundary

`CliffordAlgebra<Ordinal>` would be the char-2 mirror of `SurrealAlgebra` (the
transfinite char-0 Clifford algebra), completing `No ↔ On₂` at the Clifford layer
exactly as `NimberGame` completed it at the games layer. A metric like
`q = [ω, ω+1]` would carry genuinely **infinite nimber squares**.

`Ordinal` now implements `Scalar`, but the totality issue remains explicit:
`Scalar::mul` is panic-on-escape and `Ordinal::nim_mul` is the non-panicking
mathematical surface. Products inside the source-verified Kummer tower are exact;
products past the verified table or outside the staged segment are rejected.

### The honest design

`Scalar for Ordinal` follows the **`Rational` precedent** (`Rational` is already an
overflow-prone `i128` engine-validation scalar, not the "real" char-0 home — that
is `Surreal`). The `mul` panic message names the verified-tower escape, while
`nim_mul` / `checked_inv` are available for callers that need an explicit `Option`
boundary.

### What it actually adds (be honest)

The finite odd-degree char-2 fields (`F₈`, `F₃₂`, …) are **already** reachable as
Clifford coefficients via `Fpn<2,N>` (and, with Bridge B, classifiable). So the
*genuine* novelty of `OrdinalAlgebra` is narrow but real: **transfinite**
coefficients — `ω`, `ω+1` as squares — the exact char-2 twin of `SurrealAlgebra`'s
`ω`/`ε`. It is a symmetry-completion and a demo of the `No ↔ On₂` mirror, not a new
computational capability over the finite case.

### Classification boundary

This bridge does not try to classify every `Metric<Ordinal>`.

- Purely finite ordinal entries delegate to the existing `Nimber` Arf route.
- Entries in the first transfinite finite window `F_4(ω) = F_64` use the same
  generic symplectic reduction and the six-term absolute trace.
- Larger staged finite fields and genuinely transfinite coefficients return `None`
  for Arf/Witt/Brauer-Wall. The general finite-subfield detector and the
  transfinite classifier are separate work; the latter remains an open problem.

### Implemented surface

- `scalar/big/ordinal/` — `impl Scalar for Ordinal` (panic-on-escape `mul`,
  `neg = id`, `characteristic() = 2`, `nim_mul`, and `checked_inv`).
- `clifford` tests build `CliffordAlgebra<Ordinal>` over `q = [ω, ω+1]`, check the
  Clifford relations, and exercise associativity over the transfinite metric.
- `forms/char2/arf.rs` and the classifier façade expose finite-window
  `Metric<Ordinal>` classification and deliberately return `None` outside it.

---

## Status Snapshot

All four bridges are independently implemented and tested in the Rust core:

- **A:** lattice signature, rational/char-2 Clifford metrics, discriminant forms,
  Milgram Gauss sums, and genus oddity cross-checks.
- **B:** generic finite characteristic-2 Arf classification over supported
  `Fpn<2,N>` fields, wired into classify/Witt/isometry/Brauer-Wall façades.
- **C:** Frobenius/Galois automorphisms as Clifford `LinearMap`s with
  outermorphism spectral tests.
- **D:** `Ordinal` as a checked/panic-on-escape `Scalar`, `CliffordAlgebra<Ordinal>`
  engine tests, and finite-window ordinal Arf classification.

The second-wave bridges **E, H, and I** are now implemented and tested in the Rust
core: theta/modular forms, code↔lattice Construction A, and the discriminant-form
Weil representation. Bridge F (the rational Brauer/Clifford invariant correction)
remains a proposed build target.

Remaining open edges are not implementation TODOs inside this roadmap: the natural
Gold-quadric game rule, game-native quadratic deformation of `GameExterior`, and
the genuinely transfinite Arf/Witt classifier all stay in `OPEN.md`.

---

# Second wave — E/H/I implemented, F proposed

The first wave (A–D) closed the *pillar graph*: every pair of pillars that can talk
now does. The second wave **deepens the spine** — it strengthens the mod-8 / `E₈` /
local↔global thread the project is already built around, rather than reaching for a
new pillar. Bridges **H, E, and I** below are now standard math made computational
in the core; Bridge **F** remains design-only.

Claim-level discipline still applies: each proposed bridge is **standard math made
computational**, the same status A–D shipped at — *not* a new theorem. Where the
naive statement is subtly wrong, the corrected statement is given inline (Bridge F
in particular: the Hasse invariant is **not** simply the Brauer class of the
Clifford algebra).

**Build order: H → E → I → F.** `codes.rs` (H) is the substrate and yields the
`D₁₆⁺` lattice that the Bridge E headline needs; E is the visible punchline; I
connects E back to the already-built Bridge A. Those three are built. F is the
most careful remaining work and is independent of the other three. Bridge **G**
(spinor genus) is noted at the end as a *deferred* bridge — classical but not
buildable from the current surface.

```
            (built A–E, H, I; F still proposed)
   codes ──Construction A── integral/lattice ──θ series── modular forms   (E, H)
     │  MacWilliams              │   │                          ▲
   weight enum ↔ theta          │   └── discriminant form ──Weil rep──┘   (I)
                                 │        (Bridge A)
   clifford even-subalgebra ──Clifford invariant── local_global Hilbert    (F)
                                              └── witt/Brauer (rational)
```

## Bridge E — theta series, modular forms, and the Milnor isospectral pair

**Pillars:** `forms/integral/` ↔ a small new modular-forms layer.
**Claim level:** IMPLEMENTED AND TESTED — standard math (Hecke; Milnor 1964; Conway–Sloane
Ch. 7) made computational. **The headline bridge of the second wave.**

### The mathematics

For a **positive-definite even** lattice `L` of rank `n` (Gram `G`), the theta
series is the generating function of representation numbers

```text
θ_L(τ) = Σ_{v ∈ L} q^{Q(v)/2} = Σ_{m ≥ 0} r_L(m) q^m,   q = e^{2πiτ},
r_L(m) = #{ v ∈ L : Q(v) = 2m }   (even ⇒ Q(v) ∈ 2ℤ, so the exponents are integers).
```

When `L` is even **unimodular** (so `n ≡ 0 (mod 8)`), `θ_L` is a modular form of
weight `n/2` for the **full** modular group:

```text
θ_L ∈ M_{n/2}(SL₂(ℤ)),    M_*(SL₂ℤ) = ℂ[E₄, E₆],
E₄ = 1 + 240 Σ σ₃(m) qᵐ,    E₆ = 1 − 504 Σ σ₅(m) qᵐ,    Δ = (E₄³ − E₆²)/1728.
```

The spaces are tiny: `dim M₄ = dim M₈ = 1`, `dim M₁₂ = 2`. Because `θ_L` has
constant term `1` (the zero vector), low-dimensionality forces *exact* identities:

- **n = 8:** `θ_{E₈} = E₄` (forced, `dim M₄ = 1`). The `q¹` coefficient is
  `r_{E₈}(1) = 240 = 240·σ₃(1)` — the 240 roots / kissing number already computed in
  `root_lattices.rs`.
- **n = 16 — the Milnor punchline.** `E₈ ⊕ E₈` and `D₁₆⁺` are the two even
  unimodular lattices of rank 16. Both `θ` are weight-8 with constant term 1, and
  `dim M₈ = 1`, so

  ```text
  θ_{E₈⊕E₈} = θ_{D₁₆⁺} = E₄² = 1 + 480 q + 61920 q² + …
  ```

  identically — yet the two lattices are **not isometric** (this is Milnor's
  example of isospectral non-isometric flat tori, "you can't hear the shape of a
  16-dimensional drum"). The shared `q¹` coefficient `480` is both root systems'
  count. The equality holds to **all** orders because `dim M₈ = 1` — the test
  checks finitely many coefficients; the mathematics supplies the rest.
- **n = 24 — Leech as a free oracle.** `Λ₂₄` is already built (`mass_formula::leech`)
  and has **no roots** (`r(1) = 0`). In `M₁₂ = ⟨E₄³, Δ⟩` the unique form with
  constant term 1 and zero `q¹` coefficient is `E₄³ − 720Δ`, so `θ_{Leech} = E₄³ −
  720Δ` is *pinned by the existing rootlessness check* — a strong internal oracle
  that needs no new lattice.

**Siegel–Weil (second rung, honest).** The mass-weighted average of `θ` over a
genus equals an Eisenstein series. At `n = 16` this is **consistent but degenerate**:
both class representatives have `θ = E₄²`, so the average is trivially `E₄²`. The
genuinely non-trivial check needs a genus whose classes have *different* theta
series (`n = 24`'s 24 Niemeier classes, or a small multi-class non-unimodular
genus). Ship the `n = 16` consistency check, document the degeneracy, and mark the
non-trivial Siegel–Weil as a further rung.

### Implemented surface

- `forms/integral/theta.rs`
  - `IntegralForm::theta_series(&self, terms: usize) -> Option<Vec<i128>>` — the
    first `terms` representation numbers, bucketing `short_vectors(2·(terms−1))` by
    `Q/2`. `None` for indefinite lattices (the same boundary `minimum`/`short_vectors`
    already draw). Exact integer counts.
- `forms/integral/modular.rs`
  - `eisenstein_e4(terms)`, `eisenstein_e6(terms) -> Vec<Rational>` — exact
    q-expansions via `σ₃`/`σ₅`.
  - `mk_basis(weight, terms) -> Vec<Vec<Rational>>` — the monomial basis
    `{ E₄ᵃ E₆ᵇ : 4a + 6b = weight }` of `M_{weight}(SL₂ℤ)`.
  - `as_modular_form(q_expansion, weight, terms) -> Option<Vec<Rational>>` — solve
    for the basis coordinates on the first `dim M_weight` coefficients, then assert
    the remaining computed coefficients match. This is the **rigorous** bridge:
    equality of two weight-`k` forms agreeing through `dim M_k` coefficients is
    exact, not numerical.
- `d16_plus()` via Bridge H's `construction_a` on the indecomposable Type II
  length-16 code.

### Oracles / implemented tests

- `θ_{E₈} = E₄`; `r(1) = 240`.
- `θ_{E₈⊕E₈} = θ_{D₁₆⁺} = E₄²` to many terms, while `Genus`/isometry confirm the two
  lattices are **in the same genus but not isometric** — the Milnor pair, executable.
- `θ_{Leech} = E₄³ − 720Δ`, pinned by `r(1) = 0`.
- `as_modular_form` round-trips each of the above into `mk_basis` coordinates.
- Siegel–Weil `n = 16` consistency (degenerate), with the closed-form `|Aut|`
  constants (`|W(E₈)|`, `|Aut(D₁₆⁺)| = 2¹⁵·16!`) recorded as constants — brute-force
  `automorphism_group_order` returns `None` past its node budget, so this follows the
  `LEECH_AUT_ORDER` convention.

### Scope / caveats

- Positive-definite only (indefinite theta is not a holomorphic modular form).
- Even lattices for the clean full-level statement; odd lattices and level-`N`
  lattices give `Γ₀(N)` forms — a documented boundary tied to the existing `level()`.
- All coefficients exact (integer counts; rational Eisenstein). No floating point —
  the identification is by finite-dimensionality, not numerical agreement.

---

## Bridge H — Construction A: codes ↔ lattices, MacWilliams ↔ theta transformation

**Pillars:** a new `forms/integral/codes.rs` ↔ `forms/integral/` (lattices, theta)
↔ `forms/char2/` and `clifford_metric_f2` (the F₂ refinement).
**Claim level:** IMPLEMENTED AND TESTED — standard math (Conway–Sloane Ch. 7; MacWilliams). The
**most on-spine** second-wave idea: it is "the same duality read three ways."

### The mathematics

A binary linear code `C ⊆ F₂ⁿ` of dimension `k`. **Construction A**:

```text
L_C = (1/√2) · { x ∈ ℤⁿ : (x mod 2) ∈ C }.
```

- `det L_C = 2^{n − 2k}`; `C` **self-dual** (`k = n/2`) ⇒ `L_C` **unimodular**.
- `C` **doubly-even** (every weight `≡ 0 mod 4`) and self-dual ⇒ `L_C` **even
  unimodular** ⇒ (Bridge E) `θ_{L_C} ∈ M_{n/2}(SL₂ℤ)`.
- The Hamming weight enumerator `W_C(x,y) = Σ_{c∈C} x^{n−wt(c)} y^{wt(c)}` determines
  the theta series through the Jacobi theta constants:

  ```text
  θ_{L_C}(τ) = W_C( θ₃(2τ), θ₂(2τ) ),
  θ₃(τ) = Σ_m q^{m²},   θ₂(τ) = Σ_m q^{(m+1/2)²}.
  ```

- **MacWilliams identity** `W_{C⊥}(x,y) = |C|⁻¹ · W_C(x+y, x−y)` is the *finite*
  shadow of the modular transformation `θ(−1/τ) ↔ τ^{n/2} θ(τ)`: code duality,
  lattice unimodularity, and modular invariance are **one** phenomenon. For a
  doubly-even self-dual code the enumerator is fixed by the order-8 Gleason group —
  the discrete reflection of `M_*(SL₂ℤ) = ℂ[E₄, E₆]`.

**Corrections (caught in review — do not ship the naive versions):**

1. The `1/√2` scaling is **required**: without it self-dual codes do not give
   unimodular lattices. Since `IntegralForm` wants an integer Gram, build an integer
   basis of the preimage `{x ∈ ℤⁿ : x mod 2 ∈ C}` and carry the `1/2` in the
   dot-product — exactly the trick `leech()` uses when it divides its Gram by 8.
2. **Golay Construction A is *not* Leech.** Bare Construction A on the extended
   Golay `[24,12,8]` code gives an even unimodular rank-24 lattice, but it **has
   roots** (the images of `2eᵢ` have norm 2). The Leech lattice is the *refined*
   glue/shift construction already in `mass_formula::leech`. Phrase H as the code↔
   lattice **interface**, with Leech as its known rootless refinement — never
   "Golay → Leech."

### Implemented surface

- `forms/integral/codes.rs`
  - `BinaryCode` (checked row-reduced F₂ row space).
  - `dual`, `is_self_dual`, `is_self_orthogonal`, `is_doubly_even`, `minimum_distance`,
    `weight_enumerator(&self) -> Vec<i128>`, `macwilliams_transform(&self) -> Vec<i128>`.
  - `construction_a(&self) -> Option<IntegralForm>` (integer Gram, `1/2`-scaled;
    `None` outside the integral-Gram boundary).
  - `theta_series_via_weight_enumerator(&self, terms) -> Option<Vec<i128>>`.
  - `golay_code()` (shared with `mass_formula::leech`), `hamming_code()`,
    `extended_hamming_code()`, the split `E₈⊕E₈` Type II length-16 code, and the
    indecomposable Type II length-16 code that yields `D₁₆⁺` for Bridge E.

### Oracles / implemented tests

- MacWilliams: `code.macwilliams_transform() == code.dual().weight_enumerator()` on
  Hamming `[7,4]` and Golay `[24,12]`.
- A doubly-even self-dual code ⇒ `construction_a(C).is_even() && .is_unimodular()`.
- `W_C(θ₃(2τ), θ₂(2τ)) == construction_a(C).theta_series(…)` on small codes — the
  bridge to E.
- The Type II length-16 code's `construction_a` is `D₁₆⁺`, feeding Bridge E's Milnor
  test; and Golay's `construction_a` is even unimodular rank 24 **with** roots
  (`short_vectors(2)` nonempty), pinned **distinct** from `leech()`.

### Scope / caveats

Binary codes and Construction A only (not B/D/E); the weight-enumerator↔theta
identity uses the Hamming enumerator and the exact `θ₂`/`θ₃` q-expansions.

---

## Bridge I — the Weil representation of the discriminant form

**Pillars:** `forms/integral/discriminant.rs` (Bridge A) ↔ `forms/integral/theta.rs`
(Bridge E) ↔ `forms/witt/brauer_wall` (the mod-8 phase).
**Claim level:** IMPLEMENTED AND TESTED — standard math (Weil; Nikulin; Borcherds). The elegant
connector: it makes the **already-built** Bridge A the local-global "bulk" whose
unimodular boundary is exactly Bridge E.

### The mathematics

The finite quadratic module `(A_L, q_L)` of Bridge A carries the **Weil
representation** `ρ_L` of (a metaplectic cover of) `SL₂(ℤ)` on `ℂ[A_L] = ⊕_{γ∈A_L}
ℂ·e_γ`, generated by the two standard generators `T = [[1,1],[0,1]]`,
`S = [[0,−1],[1,0]]`:

```text
ρ_L(T) e_γ = e^{ πi · q_L(γ) } · e_γ                                  (diagonal)
ρ_L(S) e_γ = (σ / √|A_L|) · Σ_{δ ∈ A_L} e^{ −2πi · b_L(γ,δ) } · e_δ   (finite Fourier)
σ = e^{ −2πi · sign(L) / 8 }   (the conjugate of the positive Milgram phase
                                  convention used by `GaussSum`).
```

The **vector-valued theta** `Θ_L = Σ_γ θ_{L+γ} e_γ` transforms under `ρ_L`. When `L`
is **unimodular**, `A_L = 0`, `ℂ[A_L] = ℂ`, `ρ_L` is the scalar weight-`(sign/2)`
multiplier, and `Θ_L` collapses to the scalar modular form of Bridge E. So Bridge I
is the bulk and Bridge E is its boundary.

The payoff is a **third independent route to `sign mod 8`** (after the rational
signature and the genus oddity that Bridge A already cross-checks): the `S`
prefactor is the conjugate phase, and `weil_s_recovers_milgram_phase_mod8` recovers
Bridge A's positive `phase_mod8`. The honest metaplectic relations are
`S² = σ²·(γ ↦ −γ)`, `S⁴ = σ⁴·I`, and `(ST)³ = S²`; for unimodular signature
`0 mod 8` they collapse to the familiar scalar relations.

### Implemented surface

- `forms/integral/discriminant.rs`
  - `Complex64` — dependency-free complex entries for Gauss sums and Weil matrices.
  - `DiscriminantForm::weil_t(&self)` — the diagonal `T`-multipliers `e^{πi q_L(γ)}`.
  - `DiscriminantForm::weil_s(&self) -> Option<Vec<Vec<Complex64>>>` — the `S`
    matrix (`f64`; exact cyclotomic storage remains unnecessary here).
  - `weil_s_prefactor_phase_mod8` and `weil_s_recovers_milgram_phase_mod8`.
  - `verify_weil_relations(&self) -> bool` — the corrected metaplectic relations
    above plus the Milgram phase recovery.

### Oracles / implemented tests

- The metaplectic relations on the `A_n`/`D_4`/`E_8` discriminant forms already
  exercised by Bridge A.
- `ρ(S)` prefactor recovers Bridge A's Milgram `phase_mod8` after conjugating back.
- Unimodular `E₈` ⇒ `|A_L| = 1`, a `1×1` scalar collapse whose weight matches Bridge
  E's `θ_{E₈} = E₄`.

### Scope / caveats

Even lattices (so `q_L` is well-defined), matching Bridge A's boundary; matrices in
`f64` with verified unit modulus, the same convention the Gauss sum uses.

---

## Bridge F — the rational Brauer class: Hasse invariant vs Clifford invariant

**Pillars:** `clifford/` (even subalgebra) ↔ `forms/local_global/` (Hilbert symbols)
↔ a new rational Brauer class in `forms/witt/`.
**Claim level:** PROPOSED — standard math (Lam, *Introduction to Quadratic Forms
over Fields*, Ch. V; Serre). The char-0/odd mirror of Bridge B (which classified
the **char-2** Clifford algebra by its Arf/Brauer–Wall bit). **Read the corrected
statement below** — the naive "Hasse invariant = Brauer class of the Clifford
algebra" is *false*, and the codebase already declines to claim it
(`forms/char0.rs` notes rational classification is not a full Brauer/BW class).

### The mathematics (corrected)

Over `ℚ`, the quadratic-form invariants live in `Br(ℚ)[2]`, which by
Hasse–Brauer–Noether injects into `⊕_v Br(ℚ_v)[2] = ⊕_v {±1}` — a finite set of
ramified places of even cardinality (`∏_v = +1`, Hilbert reciprocity, already an
oracle in `local_global/`). Two **distinct** invariants of `⟨a₁,…,aₙ⟩`:

```text
Hasse–Witt   s(q) = ∏_{i<j} (aᵢ, aⱼ)_v          (Serre; the per-place pieces are
                                                  already in hasse_at_place / hilbert_product)
Clifford     c(q) = [ Cl⁰(q) ] ∈ Br[2]          (the class of the even Clifford algebra)
```

They are **not equal**. They differ by an explicit factor built from `(−1,−1)`,
`(−1, d)`, `(d, d)` (`d = disc q`) determined by `n mod 8` — **Lam, Prop. V.3.20**
(table). The honest bridge therefore verifies the *correction*, not an identity:

1. forms side: `s(q)` from Hilbert products, then apply the `n mod 8`/`disc`
   correction to obtain `c(q)`;
2. clifford side: read the Brauer class of `Cl⁰(q)` directly for small forms (e.g.
   identify the quaternion factor `(a, b)` of a ternary/quaternary form) as an
   independent oracle.

This is precisely the char-0 analogue of Bridge B: the algebra the `clifford` pillar
builds, classified by the symbols the `forms` pillar computes — done correctly.

### Proposed surface

- `forms/witt/brauer_rational.rs`
  - `Brauer2Class { ramified: BTreeSet<Place> }` with XOR (symmetric-difference)
    addition — the rational 2-torsion Brauer class as its ramification set.
  - `hasse_brauer_class(entries: &[i128]) -> Brauer2Class` (Hilbert-symbol product
    over all places of ℚ).
  - `clifford_brauer_class(entries: &[i128]) -> Brauer2Class` (`hasse` + the
    `n mod 8`/`disc` correction table).
- A `clifford`-side reader for small forms (via `even_subalgebra` / quaternion
  identification) as the independent oracle.

### Oracles / proposed tests

- Reciprocity: every `Brauer2Class` has `|ramified|` even.
- Known algebras: `⟨1,−1⟩` split (∅ ramified); `⟨−1,−1,−1⟩` → Hamilton quaternions,
  ramified `{2, ∞}`; a spread of ternary/quaternary forms across each `n mod 8`.
- The correction table itself: `c(q)` vs `s(q)` per dimension class.
- Agreement with `bw_class_real` / Witt `e₂` where the surfaces overlap.

### Scope / caveats

`ℚ` (and `ℚ_v`) only; 2-torsion only (quadratic-form Brauer classes are 2-torsion).
**Do not** conflate `Brauer2Class` (ungraded Brauer) with the graded
`BrauerWallClass` until a rational Brauer–Wall story is separately modeled — keeping
them distinct is the whole reason `char0.rs` currently stops short, and F is what
would add the ungraded rational class correctly.

---

## G — spinor genus (deferred, noted for completeness)

Refining `genus → spinor genus → isometry class` via the spinor norm is classical
(Eichler; Cassels–Hall), and the `clifford/spinor_norm.rs` map is the right
primitive in spirit. But it is **not buildable from the current surface**:
`spinor_norm` computes one versor's norm, whereas the spinor genus needs the local
spinor-norm *images* `θ(O(L ⊗ ℤ_p))` at every prime plus adelic class-group
bookkeeping and the proper/improper class distinction. The one cheap, honest piece
is **Eichler's theorem** as a documented predicate — *indefinite, rank ≥ 3* ⇒ spinor
genus = isometry class — which would let `Genus` upgrade to a class statement in
exactly that regime. The full definite-lattice computation is a larger build; it
stays out of the second wave, adjacent to `OPEN.md` rather than scheduled here.
