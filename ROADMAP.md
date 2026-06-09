# ROADMAP — cross-domain connections

This file is the *ambition* document: cross-pillar bridges worth building before
or shortly after the first public release. It is deliberately distinct from
`OPEN.md`:

- **`OPEN.md`** holds *genuine research problems* — things with no known answer
  (the natural Gold-quadric game rule, a game-native quadratic deformation of
  `GameExterior`, transfinite nim excesses past the verified table).
- **`ROADMAP.md`** (this file) holds *buildable bridges* — connections between the
  four mature pillars whose mathematics is largely standard, but whose **code does
  not exist yet** and whose existence would make the project's symmetry spine
  *computational* rather than merely documented. Where a bridge brushes against an
  open question, it says so and points back to `OPEN.md`.

Use the project's claim-level discipline (`AGENTS.md` → "Claim levels and
non-claims") when these land: label each piece **standard math** / **implemented
and tested** / **interpretation** / **open**.

## Why these four

The five pillars currently connect like this:

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

Trace the edges and four are conspicuously **missing or partial**:

1. **`integral ↔ clifford` does not exist at all** — the lattice pillar and the
   Clifford engine never reference each other (verified: zero cross-refs). Yet the
   project's central **mod-8 spine** (`BW(ℝ) ≅ ℤ/8`, the 8-fold Clifford table,
   Bott periodicity, `E₈`) lives on *both* sides and is currently linked only by a
   prose comment. → **Bridge A.**
2. **The char-2 classifier spans only one of its coefficient fields.** The Clifford
   engine *builds* `CliffordAlgebra<Fpn<2,N>>` (Clifford algebras over `F₈`, `F₃₂`,
   …) but `arf_*` is hardcoded to `Metric<Nimber>`, so those algebras **cannot be
   classified**. The odd-char leg already classifies both `Fp` *and* `Fpn`; char-2
   should too. → **Bridge B.**
3. **`scalar` Galois theory and `clifford` outermorphisms are latent twins.** The
   Frobenius is an `F_p`-linear map; `clifford/outermorphism.rs` computes
   `char_poly`/`det`/`Λᵏ`-traces of any linear map char-faithfully. Nothing wires
   them together, though the result feeds straight back into `trace_form.rs`. →
   **Bridge C.**
4. **The `No ↔ On₂` mirror is incomplete at the Clifford layer.** There is a
   transfinite *char-0* Clifford algebra (`SurrealAlgebra`, exact ω/ε squares) but
   no transfinite *char-2* one — even though `NimberGame` already completed the
   same mirror at the *games* layer. → **Bridge D.**

Building the four closes the pillar graph: every pair of pillars that *can* talk
(modulo the game-group-isn't-a-ring constraint) then does.

> **Structural note (already done).** The `forms/` layout was regrouped so the
> bridges below have clean homes: the invariant groups now live in `forms/witt/`
> (`class.rs`/`ring.rs`/`brauer_wall.rs`) and the valuation decomposition in
> `forms/springer/` (`local`/`padic`/`laurent`/`surreal`/`char2`); the numeric
> field invariants were renamed `forms/field_invariants.rs`. Public paths are
> unchanged (each hub re-exports flat). Bridge A lands in `integral/` + a new seam
> to `clifford/` and reads `witt/brauer_wall`; Bridge B lands in `char2/`.

---

## Bridge A — Lattice ↔ Clifford ↔ Brauer–Wall, via Milgram's Gauss sum

**Pillars:** `forms/integral/` ↔ `clifford/` ↔ `forms/witt/` ↔ `forms/char0`.
**Claim level:** standard math (Milgram/van der Blij; Conway–Sloane) made
computational. The headline bridge — it proves the project's spine crosses pillars.

### The mathematics

For an **even** integral lattice `L` (Gram `G`, so `G[i][i]` even), three objects
already half-exist in `integral/lattice.rs`:

- the **signature** `σ = p − q` (we have `is_positive_definite`; need general `(p,q)`),
- the **dual** `L# = G⁻¹L` (we have the exact `Rational` inverse used by `level`),
- the **discriminant group** `A_L = L#/L ≅ ⨁ ℤ/dᵢ`, `|A_L| = |det G|` (we have
  `invariant_factors` via Smith normal form).

The missing datum is the **discriminant quadratic form**

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

### Implementation surface

- `integral/lattice.rs`
  - `IntegralForm::signature(&self) -> (usize, usize)` — diagonalize `G` over `ℚ`
    (reuse `forms/diagonalize.rs` via `Metric<Rational>`), count sign of each pivot.
    Generalizes the current pos-definite-only path; indefinite lattices welcome.
  - `IntegralForm::clifford_metric(&self) -> Metric<Rational>` — the warm-up rung:
    `q[i] = G[i][i]`, `b[(i,j)] = 2·G[i][j]`. Feeds `CliffordAlgebra<Rational>` and
    `classify_real`. `E₈ → Cl(8,0) → M₁₆(ℝ)`. Also a mod-2 reduction
    `clifford_metric_f2(&self) -> Metric<Nimber>` for the lattice's char-2 Arf
    (Bridge B classifies it; ties type-I/II to the Arf bit).
- new `integral/discriminant.rs`
  - `DiscriminantForm { group: Vec<i128> /* the dᵢ */, values: …, gram_inv: … }`
    built from `IntegralForm` (coset-rep generators come from the SNF transform
    matrices — extend `linalg::integer::smith_normal_form` to return the unimodular
    `U`, `V`).
  - `gauss_sum(&DiscriminantForm) -> Complex64` (and an exact 8th-root-of-unity
    phase extractor; the magnitude must be `1.0` to float tolerance — assert it).
  - `milgram_signature_mod8(&DiscriminantForm) -> i64` and a checked
    `verify_milgram(lattice) -> bool` comparing it to `signature()` **and** to the
    `genus.rs` oddity route.

### Oracles / tests

`E₈` (σ≡0, `γ=1`, `Cl(8,0)=M₁₆(ℝ)`); `A_n` (`A_L = ℤ/(n+1)`, σ=n); `D_4`
(`A_L=(ℤ/2)²`, σ=4, `γ = e^{iπ} = −1`); `E₈ ⊕ E₈` vs `D₁₆⁺` (both σ≡0, distinct
genus but same `γ` — shows `γ` is coarser than genus, as it should be). The
Milgram phase must match `bw_class_real` on the lattice's own Clifford algebra.

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

### The gap

`CliffordAlgebra<Fpn<2,3>>` — a Clifford algebra over **F₈** (degree 3, which the
`u128` nimber backend cannot reach: it only holds subfields of 2-power degree) —
*builds today* but **cannot be classified**: `char2/arf.rs::arf_nimber` is
hardcoded to `Metric<Nimber>` (`.0` accessors, `nim_mul`/`nim_trace` on `u128`).
The odd-char leg (`oddchar/`) already classifies both `Fp` and `Fpn` through one
trait-generic path; char-2 should mirror that.

### Implementation surface

- `char2/arf.rs`
  - `arf_char2<F: Char2ArfField>(metric: &Metric<F>) -> Option<ArfResult>` — the
    *same* symplectic-reduction algorithm as `arf_nimber`, but over generic field
    ops (`add`/`mul`/`inv`) plus the absolute trace `Tr_{F/F₂}`. The trace already
    exists: `forms/char2/field.rs::FiniteChar2Field::artin_schreier_class` is
    `Tr_{F_q/F₂}`, impl'd for `Fp<2>`/`Fpn<2,N>`.
  - Keep `arf_nimber` for `Nimber` (its `nim_trace` is the large-field-optimized
    route, and `Nimber` is deliberately **not** `FiniteChar2Field` — respect that
    boundary; introduce a small `Char2ArfField` interface that both satisfy, or
    parameterize over the two trace sources).
- `classify.rs`
  - extend `ClassifyForm` / `BrauerWallClassify` with an `Fpn<2,N> → ArfResult`
    impl, so `metric.classify()` / `.bw_class()` work for `F₈`/`F₃₂` Clifford
    algebras, completing the façade across the whole finite char-2 leg.

### Oracles / tests

Cross-check `arf_char2` against `arf_f2` (bitmask) when all entries are in `F₂`;
genuine `F₈` forms (entries with `*4`-degree elements) hand-checked via the trace;
additivity over `⊥`; the zero-count bias `#{Q=0} = 2^{2r−1} + (−1)^Arf 2^{r−1}`
brute-forced over `F₈^{2r}` for small `r`.

### Scope / caveats

Honest non-claim (`AGENTS.md`): this is *not* a new classification theorem for all
char-2 Clifford algebras — it computes Arf/BW for the finite `Fpn<2,N>` fields,
the same status the README already states for the nimber leg.

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

### Implementation surface

- new `clifford/frobenius.rs` (or a `scalar`-side helper feeding `clifford`)
  - `frobenius_linear_map<E: CyclicGaloisExtension>() -> LinearMap<E::Base>` —
    columns are `σ(eᵢ)` in the chosen basis.
  - tests pinning `char_poly = xᵐ ± 1`, the vanishing middle `Λᵏ`-traces, and
    `det = ±1`, over `Fpn<2,m>` and `Fpn<p,m>` (odd `p`), and over `Nimber`
    (`m | 128`) using its bit basis.

### Scope / caveats

Pure cross-domain wiring + verification; no new theorem. Its value is that it makes
three pillars share one computation and gives `trace_form` a structural home.

---

## Bridge D — transfinite char-2 Clifford (`OrdinalAlgebra`)

**Pillars:** `scalar/big/ordinal` ↔ `clifford/`.
**Claim level:** the *engine* is implementable-and-tested; the *classification of
genuinely transfinite coefficients* is **open** (cross-ref `OPEN.md`). The purest
symmetry-completion, but the most caveated — lowest priority.

### The target and the blocker

`CliffordAlgebra<Ordinal>` would be the char-2 mirror of `SurrealAlgebra` (the
transfinite char-0 Clifford algebra), completing `No ↔ On₂` at the Clifford layer
exactly as `NimberGame` completed it at the games layer. A metric like
`q = [ω, ω+1]` would carry genuinely **infinite nimber squares**.

**Blocker (verified):** `Ordinal` does *not* implement `Scalar`. `Scalar::mul`
returns `Self` (total), but `Ordinal::nim_mul` returns `Option` — it is partial,
returning `None` past the source-verified Kummer tower (`α₄₇`+ or exponent `≥
ω^ω`). So `Ordinal` cannot honestly carry a *total* `mul`.

### The honest design

Implement `Scalar for Ordinal` with a **panic-on-escape `mul`**, documented under
the **`Rational` precedent** (`Rational` is already an overflow-prone `i128`
engine-validation scalar, not the "real" char-0 home — that is `Surreal`). The
`mul` panics with a precise message when a product escapes the verified tower;
products that stay inside it (which includes *every* finite nim-subfield, e.g.
`F₈ = F₂(ω)`) are exact. Add a checked entry point for the eventual Python binding
so the runtime surface never panics silently.

### What it actually adds (be honest)

The finite odd-degree char-2 fields (`F₈`, `F₃₂`, …) are **already** reachable as
Clifford coefficients via `Fpn<2,N>` (and, with Bridge B, classifiable). So the
*genuine* novelty of `OrdinalAlgebra` is narrow but real: **transfinite**
coefficients — `ω`, `ω+1` as squares — the exact char-2 twin of `SurrealAlgebra`'s
`ω`/`ε`. It is a symmetry-completion and a demo of the `No ↔ On₂` mirror, not a new
computational capability over the finite case.

### Classification subtlety (where it touches `OPEN.md`)

The Arf invariant reduces to `F₂` via `Tr_{F/F₂}`, which needs a **finite** degree.
- If all metric entries lie in a common **finite** nim-subfield `F_{2^d} ⊂ On₂`
  (including odd-degree ones reached through `ω`), Arf is well-defined with degree
  `d` — reuse Bridge B's generic `arf_char2` after detecting the subfield.
- For **genuinely transfinite** entries the absolute trace is over an infinite
  extension and Arf-to-`F₂` is *not* defined as-is. That is a real open question
  (the transfinite Witt/Arf theory), and belongs in `OPEN.md` alongside the
  existing transfinite-nim thread — add a note there rather than forcing an answer.

### Implementation surface

- `scalar/big/ordinal/` — `impl Scalar for Ordinal` (panic-on-escape `mul`,
  `neg = id`, `characteristic() = 2`, `inv` via `nim_inv` analogue returning
  `None` when partial), plus axiom tests on the verified-tower subclass (the
  existing `ordinal_partial_field_axioms` test is the seed).
- engine tests: `CliffordAlgebra<Ordinal>` Clifford relations with `q = [ω, …]`,
  the associativity suite over a transfinite metric, `e_i² = ω`.
- defer the Python `OrdinalAlgebra` binding until after Codex's pass (and only via
  the checked, non-panicking entry point).

---

## Suggested sequencing

1. **B first** — small, self-contained, unblocks classification of the `Fpn`
   Clifford algebras the engine already builds, and is a prerequisite for the
   classification half of D.
2. **A** — the headline; the most on-theme and the most visible for a public
   release (the mod-8 spine made computational). Needs the `signature()`
   generalization and the SNF transform-matrix return.
3. **C** — quick, elegant, ties three pillars; good "small beautiful thing".
4. **D** — last; design-bless the panic-`mul` with a9 first, ship the engine,
   route classification through B, push the transfinite-Arf question to `OPEN.md`.

Each is independently shippable and independently testable; none requires the
others except D's classifier leaning on B.
