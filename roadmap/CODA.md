# Cross-pillar bridges — CODA (the built-bridge record & formal notes)

This is the **long-form record** of the cross-pillar bridges that were built during the
construction era — every connection between the four mature pillars that is
**implemented and tested in the Rust core**, with its mathematical contract, implemented
surface, oracles, honest boundaries, and (for Bridges J and K) the full
formalization-pass appendices with proofs. It is the "structural proofs and stuff that
doesn't belong in AGENTS.md": the terse working-notes summary of all this lives in the
`AGENTS.md` files (root + per-pillar); this is the reference behind it.

The bridges recorded here: **A–D** (first wave), **E/F/H/I** (second wave), **J** (third
wave), **M/N/O** (fourth wave), **K** (fifth wave). Companion documents:

- **`roadmap/CODA.md`** (this file) — the built-bridge record + formal appendices.
- **`roadmap/TODO.md`** — the game-valued ledger of **buildable** work (numbers,
  switches, ups) plus the **deferred** stars: `*1` (spinor genus), `*2` (the char-`p`
  Drinfeld/Carlitz mirror), `*4` (the wild local symbol) — real and on-thesis, not
  scheduled.
- **`roadmap/DONE.md`** — the **go-forward ledger** for newly completed work after the
  construction-era record migrated here.
- **`OPEN.md`** — *genuine research problems* with no known answer (the natural
  Gold-quadric game rule, a game-native quadratic deformation of `GameExterior`,
  transfinite nim excesses past the verified table, and the transfinite Arf/Witt
  question). Where a bridge brushes against one, it says so and points back to `OPEN.md`.

Use the project's claim-level discipline (`AGENTS.md` → "Claim levels and
non-claims"): every piece here is labeled **standard math** / **implemented and
tested** / **interpretation** / **open**.

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
now meet in `integral/lattice/` and `integral/discriminant/`:

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

- `integral/lattice/`
  - `IntegralForm::signature(&self) -> (usize, usize)` diagonalizes `G` over `ℚ`
    and counts signs of the rational pivots, so indefinite lattices are supported.
  - `IntegralForm::clifford_metric(&self) -> Metric<Rational>` — the warm-up rung:
    `q[i] = G[i][i]`, `b[(i,j)] = 2·G[i][j]`. Feeds `CliffordAlgebra<Rational>` and
    `classify_real`. `E₈ → Cl(8,0) → M₁₆(ℝ)`. Also a mod-2 reduction
    `clifford_metric_f2(&self) -> Option<Metric<Nimber>>` for even lattices,
    using `Q/2 mod 2` on the diagonal and `G_ij mod 2` off-diagonal.
- `integral/discriminant/`
  - `DiscriminantForm { group, reps, gram_inv }` is built from a nonsingular even
    `IntegralForm` using the standard `A_L ~= Z^n / GZ^n` presentation. The
    representative enumeration uses normalized integer relation rows rather than
    extending Smith normal form with transform matrices.
  - `quadratic_value_mod2`, `bilinear_value_mod1`, `GaussSum::phase_mod8`,
    `fqm_gauss_phase() -> Option<FqmGaussPhase>`, and
    `milgram_signature_mod8_fqm() -> Option<i128>` make the finite quadratic module's
    p-primary Milgram/Brown phase executable. The old `GaussSum` phase stays as the
    floating oracle.
  - `FiniteQuadraticModule` is the native cyclic-product presentation of a finite
    quadratic module, and `DiscriminantForm::fqm_witt_class()` /
    `is_fqm_witt_equivalent()` give the Wall/Nikulin Witt class by reducing
    p-primary modules to canonical anisotropic cores.
  - `verify_milgram(lattice) -> Option<bool>` compares the FQM phase to the legacy
    floating Gauss-sum route, the exact signature, and the independent Conway-Sloane
    oddity route in `genus.rs`.

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
  - `arf_char2<F: FiniteChar2Field>(metric) -> Option<ArfInvariants>` runs generic
    char-2 symplectic reduction over `Fp<2>` / `Fpn<2,N>`.
  - `arf_fpn_char2<const P, const N>(metric)` is the const-generic façade helper:
    it returns `None` unless `P = 2` and the extension polynomial is supported.
  - `ArfInvariants::arf` and the Artin-Schreier class are carried as `u128` bits, in
    line with the repo-wide integer-width policy.
- `classify.rs`
  - `Fpn<P,N>` now classifies to `FiniteFieldInvariants::{Odd, Char2}`, so the same
    monomorphized façade works for odd extensions and characteristic-2 extensions.
  - `ClassifyWitt`, `ClassifyIsometry`, and `ClassifyBrauerWall` dispatch to the
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
- Otherwise the finite-subfield detector computes the minimal represented
  `F_{2^m}` by generator support plus the Frobenius fixed-field test, and the
  ordinal Arf route uses the `m`-term absolute trace. This includes the old
  `F_4(ω) = F_64` window and later staged finite fields such as the `ω^ω`
  degree-20 cell.
- Coefficients outside the source-verified staged segment still return `None` for
  Arf/Witt/Brauer-Wall. The genuinely transfinite classifier remains an open
  problem.

### Implemented surface

- `scalar/big/ordinal/` — `impl Scalar for Ordinal` (panic-on-escape `mul`,
  `neg = id`, `characteristic() = 2`, `nim_mul`, finite-subfield detection, and
  `checked_inv`).
- `clifford` tests build `CliffordAlgebra<Ordinal>` over `q = [ω, ω+1]`, check the
  Clifford relations, and exercise associativity over the transfinite metric.
- `forms/char2/arf.rs` and the classifier façade expose finite-subfield
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

The second-wave bridges **E, F, H, and I** are now implemented and tested in the
Rust core: theta/modular forms, code↔lattice Construction A, the discriminant-form
Weil representation, and the rational Brauer/Clifford invariant correction
(`forms/witt/brauer_rational.rs`).

Remaining open edges are not implementation TODOs inside this roadmap: the natural
Gold-quadric game rule, a game-native source for quadratic data beyond the checked
integer `GameClifford` surface, and the genuinely transfinite Arf/Witt classifier
all stay in `OPEN.md`.

---

# Second wave — E/F/H/I implemented

The first wave (A–D) closed the *pillar graph*: every pair of pillars that can talk
now does. The second wave **deepens the spine** — it strengthens the mod-8 / `E₈` /
local↔global thread the project is already built around, rather than reaching for a
new pillar. Bridges **E, F, H, and I** below are now standard math made
computational in the core.

Claim-level discipline still applies: each proposed bridge is **standard math made
computational**, the same status A–D shipped at — *not* a new theorem. Where the
naive statement is subtly wrong, the corrected statement is given inline (Bridge F
in particular: the Hasse invariant is **not** simply the Brauer class of the
Clifford algebra).

**Build order: H → E → I → F.** `codes.rs` (H) is the substrate and yields the
`D₁₆⁺` lattice that the Bridge E headline needs; E is the visible punchline; I
connects E back to the already-built Bridge A. F is the most careful piece (the
`n mod 8`/disc correction) and is independent of the other three. All four are now
built. Bridge **G** (spinor genus) is noted at the end as a *deferred* bridge —
classical but not buildable from the current surface.

```
            (built A–I: A–D, then E, F, H, I)
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

**Pillars:** `forms/integral/discriminant/` (Bridge A) ↔ `forms/integral/theta.rs`
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

- `forms/integral/discriminant/`
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
↔ a rational Brauer class in `forms/witt/brauer_rational.rs`.
**Claim level:** IMPLEMENTED AND TESTED — standard math (Lam, *Introduction to
Quadratic Forms over Fields*, GSM 67, pp. 117–119; Serre). The char-0/odd mirror of
Bridge B (which classified the **char-2** Clifford algebra by its Arf/Brauer–Wall
bit). The naive "Hasse invariant = Brauer class of the Clifford algebra" is *false*,
and the codebase already declined to claim it (`forms/char0.rs` notes rational
classification is not a full Brauer/BW class); F adds the **corrected** ungraded
rational class.

### The mathematics (corrected)

Over `ℚ`, the quadratic-form invariants live in `Br(ℚ)[2]`, which by
Hasse–Brauer–Noether injects into `⊕_v Br(ℚ_v)[2] = ⊕_v {±1}` — a finite set of
ramified places of even cardinality (`∏_v = +1`, Hilbert reciprocity, already an
oracle in `local_global/`). Two **distinct** invariants of `⟨a₁,…,aₙ⟩`:

```text
Hasse–Witt   s(q) = Σ_{i<j} (aᵢ, aⱼ)            (Serre; the per-place pieces are
                                                  already in hasse_at_place / hilbert_product)
Clifford     c(q) = [ C(q) ]   (n even)         (the Brauer class of the Clifford algebra;
             c(q) = [ C₀(q) ]  (n odd)            the even part in odd rank)
```

They are **not equal**. They differ by an explicit `n mod 8` / discriminant term
`δ` built from `(−1,−1)` and `(−1, d)` (`d = a₁·…·aₙ`, the **unsigned** disc) —
Lam, GSM 67, pp. 117–119 (the same table SageMath's `clifford_invariant`
implements). Additively in `Br(ℚ)[2]`:

```text
c(q) = s(q) + δ(n mod 8, d),   δ =  0                  for n ≡ 1, 2
                                    (−1,−1) + (−1, d)   for n ≡ 3, 4
                                    (−1,−1)             for n ≡ 5, 6
                                    (−1, d)             for n ≡ 7, 0
```

The honest bridge verifies the *correction*, not an identity:

1. forms side: `s(q)` from Hilbert products, then apply the `n mod 8`/`disc`
   correction `δ` to obtain `c(q)`;
2. clifford side: read the Brauer class of the Clifford algebra directly for small
   forms — `C(⟨a,b⟩) ≅ (a,b)` (n=2) and `C₀(⟨a,b,c⟩) ≅ (−ab, −ac)` (n=3, the
   quaternion factor of the even subalgebra) — as the **independent** oracle.

This is precisely the char-0 analogue of Bridge B: the algebra the `clifford` pillar
builds, classified by the symbols the `forms` pillar computes — done correctly.

### Implemented surface

- `forms/witt/brauer_rational.rs`
  - `Brauer2Class { ramified: BTreeSet<Place> }` (private field) with `add` =
    symmetric difference (XOR), `split`/`is_split`, `local_invariant`,
    `satisfies_reciprocity`, and the `quaternion(a, b)` constructor (the class of
    `(a,b)` over ℚ). The rational 2-torsion Brauer class as its ramification set.
  - `hasse_brauer_class(entries: &[i128]) -> Option<Brauer2Class>` — the per-place
    Hasse invariant collected into a ramification set.
  - `clifford_brauer_class(entries: &[i128]) -> Option<Brauer2Class>` — `hasse` +
    the `n mod 8`/`disc` correction `δ`. `None` on a zero entry (radical) or
    bounded-arithmetic overflow.
- `Place` (in `local_global/padic.rs`) gained `Ord`/`PartialOrd` so the
  ramification set is a `BTreeSet` (ℝ before `Q_2`, `Q_3`, …).

### Oracles / implemented tests

- Reciprocity: every class has `|ramified|` even (`satisfies_reciprocity`), over a
  sweep of rank-2…6 forms.
- Known algebras: `⟨1,−1⟩` split (∅ ramified); `⟨−1,−1,−1⟩` and `⟨1,1,1⟩` →
  Hamilton quaternions, ramified `{ℝ, Q_2}` — with `⟨1,1,1⟩` showing `s = 0` while
  `c = (−1,−1)`, the sharpest demonstration that `c ≠ s`.
- The **independent** clifford-side oracle, over sweeps: `clifford(⟨a,b⟩) = (a,b)`
  (n=2) and `clifford(⟨a,b,c⟩) = (−ab,−ac)` (n=3); rank-1 always split.
- The correction table itself: `c(q) = s(q) + δ` checked across `n = 1…8`, with `δ`
  recomputed independently in the test from `Brauer2Class::quaternion`.

### Scope / caveats

`ℚ` (and `ℚ_v`) only; 2-torsion only (quadratic-form Brauer classes are 2-torsion).
**Do not** conflate `Brauer2Class` (ungraded Brauer) with the graded
`BrauerWallClass` until a rational Brauer–Wall story is separately modeled — keeping
them distinct is the whole reason `char0.rs` currently stops short, and F is what
would add the ungraded rational class correctly.


---

# Third wave — Bridge J (built)

The third wave came out of a deliberate "deepen, don't sprawl" review: the project is
near-saturated on the **place axis**, so the high-leverage moves are no longer *new
number systems* but closing threads already half-drawn. Bridge **J** below is the
built member of that wave — it names the valuation as the tropicalization the
`scalar/tropical` layer already defines, and adds Newton polygons over the valued
legs, with the slope ⟺ Springer-residue-layer cross-check. The other two third-wave
bridges — **K** (the full `ℚ/ℤ` cyclic-algebra Brauer invariant) and **L** (the
char-`p` Drinfeld/Carlitz mirror) — are now built (K, the fifth wave, recorded later in
this file) and deferred (L, now `*2` in `roadmap/TODO.md`) respectively.

Claim-level discipline still applies: J is **standard math made computational**, the
same status A–I shipped at — not a new theorem.

## Bridge J — the valuation as tropicalization; Newton polygons as tropical curves

**Pillars:** `scalar/tropical` ↔ `scalar/valued` ↔ `scalar/newton` ↔ the local-field
backends (`small/`, `functor/`, `global/`) ↔ `forms/springer`.
**Claim level:** IMPLEMENTED AND TESTED — standard math (tropical geometry;
Newton–Puiseux; valuation theory) made computational. The on-thesis **twin of the
already-shipped "thermography = tropical arithmetic" identity**, applied to the
*place axis* instead of the game axis.

### The mathematics

`scalar/tropical.rs` (the `Semiring`, min-plus / max-plus) is currently consumed
**only** by `games/tropical_thermography` — it is marooned on the games side. Yet the
valuation `v : K* → Γ` on every discretely-valued backend tropicalizes `K`: it is a
**homomorphism of multiplicative monoids** into `(Γ ∪ {∞}, min, +)`, **lax (subadditive)
for addition**, strict off the tropical vanishing locus:

```text
v(x·y)  = v(x) + v(y)                       (the tropical ⊗ — strict)
v(x + y) ≥ min(v(x), v(y))                  (the tropical ⊕ — lax)
v(x + y) = min(v(x), v(y))   if v(x) ≠ v(y) (strict off the vanishing locus)
```

So the whole `Valued` stack already **is** the tropicalization map; the project computes
it everywhere and names it as such nowhere. (**Honest correction from the formalization
pass:** "*is* the tropicalization" is meant **laxly** — no discretely-valued field admits
a *strict* additive homomorphism onto `ℤ_trop`; strictness is restored only by the
tropical **hyperfield** [Viro 2010], or by taking the three lines above as the
*definition* of a valuation [Maclagan–Sturmfels Ch. 2]. The slogan must not claim
strictness.) The payoff object is the **Newton
polygon**: for `f = Σ aᵢ xⁱ ∈ K[x]`, the lower convex hull of `(i, v(aᵢ))` is a
tropical curve whose **slopes are exactly the valuations of the roots** (horizontal
length = multiplicity), and whose break structure controls factorization into pieces
of distinct root-valuation — the discrete-valuation refinement `poly_factor` / Hensel
already half-use. The Springer decomposition's "two residue layers survive because the
value group is `ℤ`" is precisely the **graded pieces of the valuation/tropical
filtration**: each Newton slope *is* a residue layer. This closes a real asymmetry —
thermography names its option-fold `⊕` and cooling `⊗`; the valuation does the
identical algebra on the scalar side and currently says so nowhere.

### Implemented surface

- `scalar/valued.rs` — the `Valued` trait docs name `valuation` as the (lax)
  tropicalization morphism into `Tropical<MinPlus>`, plus the free adaptor
  `tropicalize<K: Valued>(x: &K) -> Tropical<MinPlus>` (no new math — it names the
  existing map; its tests are truncation-safe).
- `scalar/newton.rs` — `NewtonPolygon::of(coeffs: &[K]) -> Option<NewtonPolygon>`
  over any `K: Valued` (the lower convex hull of `(i, v(aᵢ))`; `None` for the zero
  polynomial). **Orientation trap (caught in the formalization pass):** with points
  `(i, v(aᵢ))`, a side of slope `−λ` carries roots of valuation `+λ`, so
  `root_valuations() -> Vec<(Rational, u128)>` returns the **negated** slopes (with
  horizontal lengths = multiplicities) while `slopes()` is the literal hull view;
  slopes are `Rational`, since root valuations can be fractional even though `Γ = ℤ`
  (the `Ramified` `xᴱ − ϖ` case). Also `zero_root_multiplicity()` (roots at `0`,
  valuation `+∞`) and `degree()`. Exact over `Qp`/`Qq`/`Laurent`/`Ramified`,
  exact-outright over the `F_q(t)` completion (the `Laurent` leg).
- a slope ↔ Springer-residue-layer cross-check (in `forms/springer/local.rs` tests):
  the Newton polygon **is** the Springer decomposition under tropicalization — it
  sees `(valuation, dim)` per layer and forgets the residue square class, the
  forgetful hierarchy `NP(f_q) ≺ {in_λ(f_q)} ≺ q`.

### Oracles / implemented tests

- The tropicalization laws (J.1): multiplicativity, the `⊕`-internal subadditivity,
  and equality off the vanishing locus — over `Qp`/`Qq`/`Laurent`, truncation-safe.
- Eisenstein `xᴱ − p`: a single slope, every root valuation `1/E`, cross-checked
  against the `Ramified` renormalization `Ramified::<…, E>::pi().valuation() = 1`.
- `x² − p` over `Q_p`: root valuation `1/2`, agreeing with `Qp::is_square = false`.
- Dumas additivity: a product of distinct-slope factors reconstructs the polygon.
- a monic integral polynomial has an all-flat polygon ⟺ `a₀` a unit ⟺ unit roots;
  zero roots (`+∞`) tracked separately; negative-valuation (pole) roots.
- `polygon_is_the_springer_shadow`: the side multiset `{(valuation, mult)}` equals
  the Springer buckets `{(valuation, dim)}` over `Q_5`/`Q_9`/`F_7((t))`, and the
  parity grouping reproduces `parity_layer`; `polygon_outlives_springer`: over
  residue char 2 the polygon succeeds while Springer returns `None`.

### Scope / caveats

- Discretely-valued legs only. The **divisible**-value-group surreal leg has no integer
  Newton lattice — the same boundary `springer/surreal.rs` already documents, and itself
  an instance of the local↔global symmetry, not a gap.
- The capped-precision models give Newton data valid to their precision horizon; flag the
  truncation as those backends already do.
- Tropical here is `MinPlus` (valuations); the `MaxPlus` dual is the thermography
  convention. Note the sign mirror rather than duplicating the semiring.

### Formalized

The full lemmas — J.1 (valuation↔tropical dictionary, with the lax/strict subtlety),
J.3 (graded ring `gr_v K ≅ k[u,u⁻¹]`), J.5 (slope theorem, with proof), J.6 (Dumas
additivity), J.7 (Eisenstein ↔ the `Ramified` renormalization), J.12 (each Newton slope
**is** a Springer residue layer) — with proofs, the witness tests, and references
(Springer; Lam; Koblitz; Neukirch; Dumas; Serre; Maclagan–Sturmfels; Viro; Stichtenoth)
are in the formal-proofs appendix below.


---

## Bridge J — formal statements and proofs (formalization-pass appendix)

> Moved here from the former `BRIDGES-DRAFT.md` (a parallel formalization front).
> Standard math made computational unless marked; the lemma/theorem numbering (J.1,
> J.5, …) is the one the Bridge J section above refers to.

**Status.** Everything below is **standard math** (no new theorems), per the third-wave discipline in `roadmap/CODA.md` (this file). Items marked ⟦implemented⟧ are witnessed by tests in this checkout; items marked ⟦proposed⟧ name the tests that would witness the proposed `NewtonPolygon` surface. Nothing here is at *interpretation* or *open* level.

## 0. Setup and notation

Throughout, $(K, v)$ is a field with a **normalized discrete valuation**: $v : K^\times \twoheadrightarrow \mathbb{Z}$ with $v(xy) = v(x) + v(y)$ and $v(x+y) \ge \min(v(x), v(y))$, extended by $v(0) = +\infty$. Write $\mathcal{O} = \{v \ge 0\}$, $\mathfrak{m} = \{v \ge 1\}$, residue field $k = \mathcal{O}/\mathfrak{m}$, and fix the uniformizer $\varpi$ (so $v(\varpi) = 1$). The **angular component** of $x \ne 0$ is $\mathrm{ac}(x) = \overline{x\varpi^{-v(x)}} \in k^\times$ (it depends on the choice of $\varpi$).

$\mathbb{T}$ denotes the min-plus tropical semiring $(\mathbb{Q} \cup \{+\infty\},\ \oplus = \min,\ \otimes = +)$, with $\oplus$-identity $\infty$ and $\otimes$-identity $0$.

Dictionary to the code (all in `/Users/a9lim/Work/ogdoad`):

| math | code |
|---|---|
| $v$, $\varpi$ | `Valued::valuation` (`None` = $\infty$), `Valued::uniformizer` — `src/scalar/valued.rs` |
| $\mathbb{T}$ | `Tropical<MinPlus>` — `src/scalar/tropical.rs` (`Semiring`; ⟦implemented⟧, fuzzed in `tests/tropical_axioms.rs`) |
| $k$, $\mathrm{ac}$ | `ResidueField::Residue`, `residue_unit` — `src/scalar/residue.rs` |
| discretely-valued legs | `Qp<P,K>` ($v(p){=}1$), `Qq<P,N,F>` (unramified, $v(p){=}1$), `Laurent<S,K>` ($v(t){=}1$), `Ramified<S,E>` (renormalized $v(\pi){=}1$, value group $\mathbb{Z}$), `Gauss<S>` ($v(t){=}0$) |
| $\mathbb{F}_q(t)$ per place | `try_valuation_at_ff`, `FFPlace::{Finite(π), Infinite}` — `src/forms/local_global/function_field.rs` |
| Springer buckets | `springer_decompose_local`, `LocalResidueForm { valuation, dim, disc_is_square }`, `parity_layer` — `src/forms/springer/local.rs` |
| Gauss valuation on $K[y]$ | `Poly::min_coeff_valuation` (`src/scalar/poly.rs`), coefficientwise reduction at the minimum (`reduce_poly_at_min` in `src/scalar/functor/gauss.rs`) |

---

## 1. (a) The valuation is the tropicalization

**Lemma J.1 (valuation–tropical dictionary).** ⟦standard math⟧ Define $\tau : K \to \mathbb{T}$ by $\tau(x) = v(x)$ (so $\tau(0) = \infty$). Then:

$$
\begin{aligned}
\text{(i)}\quad & \tau(xy) \;=\; \tau(x) \otimes \tau(y) \quad\text{for all } x, y \in K \text{ (including } 0\text{, by absorption)};\\
\text{(ii)}\quad & \tau(x+y) \,\oplus\, \bigl(\tau(x) \oplus \tau(y)\bigr) \;=\; \tau(x) \oplus \tau(y) \quad\text{i.e.}\quad v(x+y) \ge \min(v(x), v(y));\\
\text{(iii)}\quad & \tau(x+y) \;=\; \tau(x) \oplus \tau(y) \quad\text{whenever } \tau(x) \neq \tau(y);\\
\text{(iv)}\quad & \tau(1) = 0 = 1_{\mathbb{T}}, \qquad \tau(0) = \infty = 0_{\mathbb{T}}.
\end{aligned}
$$

*Proof.* (i), (ii), (iv) restate the valuation axioms in the $(\min,+)$ dictionary; the $\oplus$-internal phrasing of (ii) uses $a \ge b \iff a \oplus b = b$ in $(\mathbb{Q}\cup\{\infty\}, \min)$. For (iii): note first $v(-1) = 0$ (since $2\,v(-1) = v(1) = 0$ in $\mathbb{Z}$), so $v(-y) = v(y)$. Assume WLOG $v(x) < v(y)$, and suppose $v(x+y) > v(x)$. Then $v(x) = v\bigl((x+y) + (-y)\bigr) \ge \min(v(x+y), v(y)) > v(x)$, a contradiction. $\blacksquare$

**Remark J.2 (how "semiring homomorphism" is meant — a non-claim).** $\tau$ is a homomorphism of multiplicative monoids $(K, \cdot, 1, 0) \to (\mathbb{T}, \otimes, 1_\mathbb{T}, 0_\mathbb{T})$ and is **lax** for addition: (ii) with equality (iii) exactly off the *tropical vanishing locus* (the locus where the minimum is attained at least twice — e.g. $v(1 + (-1)) = \infty \ne 0$). No discretely-valued field admits a *strict* additive homomorphism onto $\mathbb{T}$; strict functoriality is restored by replacing $\mathbb{T}$ with the tropical **hyperfield** [Viro 2010], or by taking Lemma J.1(i)–(iii) as the *definition* of a valuation, as in [Maclagan–Sturmfels, Ch. 2]. the Bridge J section's slogan "the valuation **is** the tropicalization" has Lemma J.1 as its precise content; prose should not claim strictness.

**Lemma J.3 (graded ring of the valuation filtration).** ⟦standard math⟧ Let $\mathfrak{m}^\lambda = \{x : v(x) \ge \lambda\}$ for $\lambda \in \mathbb{Z}$ (fractional ideals). The associated graded ring of the filtration,
$$
\mathrm{gr}_v(K) \;=\; \bigoplus_{\lambda \in \mathbb{Z}} \mathfrak{m}^{\lambda}/\mathfrak{m}^{\lambda+1},
$$
is, after the choice of $\varpi$, isomorphic to $k[u, u^{-1}]$ ($u = $ class of $\varpi$), and the leading-form map $\sigma : K^\times \to \mathrm{gr}_v(K)$, $\sigma(x) = x \bmod \mathfrak{m}^{v(x)+1}$, is multiplicative, with
$$
\sigma(x) \;=\; \mathrm{ac}(x)\, u^{v(x)}.
$$

*Proof.* Write $x = \varpi^{v(x)} u_x$ with $u_x \in \mathcal{O}^\times$; then $\mathrm{ac}(x) = \bar{u}_x$, each graded piece is a one-dimensional $k$-vector space spanned by $u^\lambda$, and multiplicativity of $\sigma$ is multiplicativity of $v$ and of the residue map on units ($k$ is a field, so there is no cancellation of leading terms). $\blacksquare$

The two lemmas together say: **the valuation/tropical filtration of $K$ has tropical shadow $\tau$ and graded pieces $k \cdot u^\lambda$** — the "residue layers" of part (c).

**Witness tests (a).**
- ⟦implemented⟧ `src/scalar/valued.rs::tests::{uniformizers_have_valuation_one, zero_valuation_is_none}` (J.1(iv) and the $\infty$ convention); `src/scalar/functor/ramified.rs::tests::valuation_is_additive_under_multiplication` (J.1(i) on the ramified leg); `tests/tropical_axioms.rs` ($\mathbb{T}$ is a semiring, both conventions).
- ⟦proposed⟧ `tests/tropicalization.rs`, with the thin adaptor (the Bridge J surface):
  ```rust
  fn trop<K: Valued>(x: &K) -> Tropical<MinPlus> {
      match x.valuation() { Some(v) => Tropical::int(v), None => Tropical::infinity() }
  }
  ```
  proptest over `Qp<5,8>`, `Qq<3,4,2>`, `Laurent<Fp<7>,8>`, `Ramified<Qp<3,8>,2>`, `Gauss<Qp<5,6>>`:
  - `tropicalize_is_multiplicative`: `trop(x.mul(&y)) == trop(&x).mul(&trop(&y))` — exact, zero included;
  - `tropicalize_is_subadditive`: `let s = trop(&x).add(&trop(&y)); trop(&x.add(&y)).add(&s) == s` — the $\oplus$-internal J.1(ii), **truncation-safe**: if a deep cancellation renders the sum as the represented $0$, the left side is $\infty$ and the identity still holds;
  - `tropicalize_equality_off_vanishing_locus`: `if trop(&x) != trop(&y) { trop(&x.add(&y)) == trop(&x).add(&trop(&y)) }` — exact even in the capped models, since the leading term survives truncation.

---

## 2. (b) The Newton-polygon slope theorem

**Definition J.4 (Newton polygon).** For $f = \sum_{i=0}^{n} a_i x^i \in K[x]$ with $a_0 a_n \ne 0$, the **Newton polygon** $\mathrm{NP}(f)$ is the lower boundary of the convex hull of $\{(i, v(a_i)) : a_i \ne 0\} \subset \mathbb{R}^2$, a convex piecewise-linear chain from $(0, v(a_0))$ to $(n, v(a_n))$ with strictly increasing side slopes in $\mathbb{Q}$. (If $a_0 = 0$, factor out $x^m$ first; those $m$ roots are $0$, "valuation $\infty$".)

*Orientation convention — an implementation trap.* With points $(i, v(a_i))$, a side of slope $-\lambda$ corresponds to roots of valuation $+\lambda$. To keep the public surface matching the Bridge J section's "slopes are the valuations of the roots", the proposed type should expose `root_valuations() -> Vec<(Rational, u128)>` (negated slopes with horizontal lengths) rather than asking callers to negate; slopes are `Rational` (ratios of `i128`) since root valuations can be fractional even though $\Gamma = \mathbb{Z}$.

**Theorem J.5 (slope theorem).** ⟦standard math: Koblitz, GTM 58, Ch. IV; Neukirch, Ch. II⟧ Let $K$ be **complete** (henselian suffices) with respect to the discrete valuation $v$, let $f \in K[x]$ with $a_0 a_n \neq 0$, let $L$ be a splitting field of $f$, and let $w$ be the unique extension of $v$ to $L$. If $\mathrm{NP}(f)$ has a side of slope $-\lambda$ with horizontal length $\ell$, then $f$ has **exactly $\ell$ roots $r \in L$ (with multiplicity) with $w(r) = \lambda$**, and every root arises this way. In particular $\sum_{\text{sides}} \ell = n$ and the multiset of root valuations is determined by the coefficient valuations alone.

*Proof.* Existence/uniqueness of $w$ on the finite extension $L/K$ is the standard consequence of completeness, $w = \tfrac{1}{[L:K]}\, v \circ N_{L/K}$ [Neukirch, Ch. II]. Normalize $f$ monic (dividing by $a_n$ translates the polygon vertically; slopes and lengths are unchanged). Write $f = \prod_{j=1}^n (x - r_j)$ with $w(r_1) \le \cdots \le w(r_n)$. The coefficients are signed elementary symmetric functions: $a_{n-m} = \pm e_m(r_1, \dots, r_n)$, so by J.1(ii)–(iii) applied in $(L, w)$:
$$
v(a_{n-m}) \;=\; w(e_m) \;\ge\; \min_{|S| = m} \sum_{j \in S} w(r_j) \;=\; \sum_{j \le m} w(r_j),
$$
with **equality whenever the minimizing $m$-subset is unique**, i.e. whenever $w(r_m) < w(r_{m+1})$, and unconditionally at $m = 0$ and $m = n$ (a unique subset each). Let $h(i) := \sum_{j \le n-i} w(r_j)$ for $i = 0, \dots, n$ (height as a function of the point index $i = n - m$). Its successive slopes are $h(i+1) - h(i) = -w(r_{n-i})$, non-decreasing in $i$ because the $w(r_j)$ are sorted — so the graph of $h$ is convex; it lies on or below every point $(i, v(a_i))$; and it touches them at $i \in \{0, n\}$ and at every index where the sorted valuations jump — exactly the vertices of the graph of $h$. Hence the lower convex hull of the points **is** the graph of $h$, and the side of slope $-\lambda$ spans exactly the indices $i$ with $w(r_{n-i}) = \lambda$, of horizontal length $\#\{j : w(r_j) = \lambda\}$. $\blacksquare$

**Lemma J.6 (additivity; Dumas).** ⟦standard math: Dumas 1906⟧ For $f, g \in K[x]$ with nonzero constant terms, the sides of $\mathrm{NP}(fg)$ are obtained by concatenating the sides of $\mathrm{NP}(f)$ and $\mathrm{NP}(g)$ in increasing slope order; per-slope horizontal lengths add.

*Proof (complete case, which is all the project legs need).* Immediate from Theorem J.5: the root multiset of $fg$ in a common splitting field is the union of the two root multisets. (Dumas's original proof is a direct coefficient estimate and needs no completeness.) $\blacksquare$

**Corollary J.7 (Eisenstein).** ⟦standard math: Serre, *Local Fields*, Ch. I⟧ If $f$ is monic of degree $n$ with $v(a_i) \ge 1$ for $i < n$ and $v(a_0) = 1$, then $\mathrm{NP}(f)$ is the single side from $(0,1)$ to $(n,0)$, so every root has valuation $1/n$; $f$ is irreducible, and a root generates a totally ramified extension of degree $n$.

*Proof.* The polygon claim is immediate (all interior points lie on or above the segment). If $h \mid f$ is monic of degree $d$, then $v(h(0)) = \sum_{d \text{ roots}} w(r) = d/n \in \mathbb{Z}$ forces $d \in \{0, n\}$. The value group of $K(r)$ contains $\tfrac{1}{n}\mathbb{Z}$, so $e = n = [K(r):K]$. $\blacksquare$

This is exactly the project's `Ramified<S, E>` ($x^E - \varpi$): its *renormalized* valuation $\min_i\,(E \cdot v_S(a_i) + i)$ rescales the slope-$\tfrac{1}{E}$ root to $v(\pi) = 1$, restoring $\Gamma = \mathbb{Z}$ — which is why the Newton lattice stays integral on that leg.

**Corollary J.8 (unit roots ⟺ flat polygon).** For monic $f \in \mathcal{O}[x]$: all roots of $f$ are units of (the integral closure of $\mathcal{O}$ in) $L$ $\iff$ $\mathrm{NP}(f)$ is the single horizontal side at height $0$ $\iff$ $v(a_0) = 0$ $\iff$ the residue reduction $\bar{f} \in k[x]$ has $\bar{f}(0) \ne 0$.

*Proof.* $v(a_0) = \sum_j w(r_j)$ with every $w(r_j) \ge 0$ (monic, integral coefficients, J.5), so the sum vanishes iff every term does. $\blacksquare$

**Corollary J.9 (per-place polygons over the global $\mathbb{F}_q(t)$).** ⟦standard math: Stichtenoth, GTM 254, Ch. 1⟧ For $f \in \mathbb{F}_q(t)[x]$ and a place $P$ of $\mathbb{F}_q(t)$ (a monic irreducible $\pi$, or $\infty$ with $v_\infty = \deg \mathrm{den} - \deg \mathrm{num}$), the polygon $\mathrm{NP}_P(f)$ computed from the **exact** valuations $v_P(a_i)$ equals the Newton polygon of $f$ over the completion $\mathbb{F}_q(t)_P \cong \mathbb{F}_{q^{\deg P}}((\pi))$, and Theorem J.5 applies there. (The completion at a degree-1 finite place is literally the `Laurent` backend; coefficient valuations are insensitive to completion, so the global leg's polygon is exact with no precision model at all.)

**Witness tests (b)** — all ⟦proposed⟧, on `NewtonPolygon::of(coeffs: &[K]) -> NewtonPolygon` for `K: Valued`:
- `eisenstein_single_slope`: $\mathrm{NP}(x^E - p)$ over `Qp<5,8>` has one side, `root_valuations() == [(1/E, E)]`; cross-check `Ramified::<Qp<5,8>, E>::pi().valuation() == Some(1)` (J.7 ↔ the renormalization).
- `sqrt_p_slope_half`: $\mathrm{NP}(x^2 - p)$ over `Qp<5,8>` gives root valuation $\tfrac12 \notin \mathbb{Z}$; cross-check `Qp::<5,8>::from_i128(5).is_square() == Some(false)` (odd valuation ⇒ nonsquare; `src/scalar/small/analytic.rs`).
- `dumas_additivity`: for $f, g$ with distinct slopes over `Qp`/`Laurent`, per-slope lengths of $\mathrm{NP}(fg)$ are the sums (J.6).
- `flat_polygon_iff_unit_roots`: monic integral $f$; all-zero slopes $\iff$ `a₀.valuation() == Some(0)` $\iff$ the residue reduction has nonzero constant term (J.8, via `ResidueField::residue`).
- `ff_place_polygon_matches_completion`: $f$ over `RationalFunction<Fp<5>>` at the place $t$: polygon from `try_valuation_at_ff` equals the polygon of the coefficientwise image in `Laurent<Fp<5>, K>` (J.9 — the exact-global vs local-model agreement).

---

## 3. (c) Slopes are the Springer residue layers

**Theorem J.10 (Springer).** ⟦standard math: Springer, Indag. Math. 17 (1955); Lam, GSM 67, Ch. VI⟧ Let $K$ be complete discretely valued with $\operatorname{char} k \ne 2$, and fix $\varpi$. Every nondegenerate diagonal form over $K$ is isometric to $q_0 \perp \varpi\, q_1$ with $q_0, q_1$ having unit diagonal entries, and the two **residue homomorphisms** $\partial_0, \partial_1$ (sending $\langle u \rangle \mapsto \langle \bar{u} \rangle$ and $\langle \varpi u \rangle \mapsto \langle \bar{u} \rangle$ respectively) induce a group isomorphism
$$
(\partial_0, \partial_1) : W(K) \;\xrightarrow{\ \sim\ }\; W(k) \oplus W(k),
$$
where $\partial_1$ (not $\partial_0$) depends on the choice of $\varpi$. The two summands are indexed by $\Gamma/2\Gamma = \mathbb{Z}/2$ — they exist *because* the value group is not 2-divisible: $\langle \varpi^2 a \rangle \cong \langle a \rangle$, while $\langle \varpi a \rangle \not\cong \langle a \rangle$ in general.

This is the theorem behind `springer_decompose_local` + `parity_layer` ⟦implemented: `src/forms/springer/local.rs::tests::*`⟧; the code records, per valuation $\lambda$, the layer $(\lambda, \dim, \mathrm{disc\ square\text{-}class})$, and `parity_layer(ε)` is the data of $\partial_\varepsilon$.

**Definition J.11 ($\lambda$-initial form — the graded/tropical piece).** For $\lambda \in \mathbb{Z}$ and $f = \sum a_i x^i \in K[x]$, let
$$
m_\lambda(f) \;=\; \min_i \bigl(v(a_i) + i\lambda\bigr) \;=\; \bigoplus_i \tau(a_i) \otimes \lambda^{\otimes i} \quad(\text{the tropicalized } f \text{ evaluated at } \lambda),
$$
and define the **initial form** $\mathrm{in}_\lambda(f) \in k[y]$ as the coefficientwise reduction of $\varpi^{-m_\lambda(f)} f(\varpi^\lambda y)$ — i.e. substitute $x = \varpi^\lambda y$, then take the Gauss-valuation angular component (in the code: a $\varpi^\lambda$-shift, `Poly::min_coeff_valuation`, and the reduce-at-the-minimum step that `reduce_poly_at_min` in `src/scalar/functor/gauss.rs` already performs — `Gauss<S>` *is* the Gauss valuation this construction lives in). Two standard facts: $\lambda$ is the negative of a slope of $\mathrm{NP}(f)$ iff $\deg \mathrm{in}_\lambda(f) > \operatorname{ord}_y \mathrm{in}_\lambda(f)$ (the minimum is attained at two distinct $i$ — the **tropical-root** criterion [Maclagan–Sturmfels, Ch. 2–3]); and $\mathrm{in}_\lambda(fg) = \mathrm{in}_\lambda(f)\,\mathrm{in}_\lambda(g)$, since the Gauss valuation is a valuation on $K[y]$ and its angular component into the domain $k[y]$ is multiplicative (Lemma J.3 applied to $\mathrm{Gauss}$).

**Proposition J.12 (slope ⟺ residue layer, for diagonal forms).** ⟦standard math; elementary given J.5/J.6 + J.10⟧ Let $q = \langle a_1, \dots, a_n \rangle$ with all $a_i \in K^\times$ (zero entries are the radical, tracked separately as `radical_dim`), and let $f_q(x) = \prod_{i=1}^n (x - a_i)$. Then:

**(i) (the polygon is the bucket shadow).** $\mathrm{NP}(f_q)$ has a side of slope $-\lambda$ and horizontal length $\ell$ $\iff$ $\#\{i : v(a_i) = \lambda\} = \ell$. Hence the side multiset of $\mathrm{NP}(f_q)$ equals the multiset $\{(\texttt{g.valuation}, \texttt{g.dim})\}$ of the Springer decomposition — every Newton slope **is** a residue layer, and conversely.

**(ii) (the initial form is the residue layer's contents).** For each such $\lambda$,
$$
\mathrm{in}_\lambda(f_q) \;=\; c\, \cdot\, y^{\,\#\{i\,:\,v(a_i) > \lambda\}} \prod_{i\,:\,v(a_i) = \lambda} \bigl(y - \mathrm{ac}(a_i)\bigr), \qquad c = \prod_{i\,:\,v(a_i) < \lambda} \bigl(-\mathrm{ac}(a_i)\bigr) \in k^\times,
$$
so the nonzero roots of $\mathrm{in}_\lambda(f_q)$ in $\bar{k}$ are exactly the angular components of the layer, and the layer discriminant is recovered as $\prod_{v(a_i) = \lambda} \mathrm{ac}(a_i)$, whose $k$-square class is `disc_is_square`.

**(iii) (the Witt-level collapse).** If moreover $\operatorname{char} k \ne 2$, the Witt class of $q$ depends only on the layers grouped by $\lambda \bmod 2$: since $\langle a \rangle \cong \langle \varpi^{\,v(a) \bmod 2}\, u_a \rangle$, one gets $\partial_\varepsilon[q] = \bigl[\bigoplus_{v(a_i) \equiv \varepsilon (2)} \langle \mathrm{ac}(a_i) \rangle\bigr] \in W(k)$, and $(\partial_0, \partial_1)$ is Springer's isomorphism. `parity_layer(ε)` computes exactly the data of $\partial_\varepsilon$.

*Proof.* (i): each factor $(x - a_i)$ has the two-point polygon with the single side of slope $-v(a_i)$ and length 1 (using $v(-a_i) = v(a_i)$); apply Lemma J.6. (ii): $\mathrm{in}_\lambda(x - a) = y - \mathrm{ac}(a)$, $y$, or $-\mathrm{ac}(a)$ according as $v(a) = \lambda$, $> \lambda$, $< \lambda$ (compute $m_\lambda = \min(\lambda, v(a))$ directly); multiply, using multiplicativity of $\mathrm{in}_\lambda$ (Definition J.11). (iii): $a = \bigl(\varpi^{\lfloor v(a)/2 \rfloor}\bigr)^2\, \varpi^{\,v(a) \bmod 2}\, u_a$ and, for units, $\langle u \rangle \cong \langle u' \rangle$ over $K$ iff $\bar{u}/\bar{u}'$ is a square in $k$ (Hensel's lemma lifts residue squares when $\operatorname{char} k \ne 2$); then apply Theorem J.10. $\blacksquare$

**Remark J.13 (the forgetful hierarchy — what each level sees).** The data refine strictly:
$$
\underbrace{\mathrm{NP}(f_q)}_{\text{tropical shadow: } (\lambda, \dim) \text{ per layer}} \;\prec\; \underbrace{\{\mathrm{in}_\lambda(f_q)\}_\lambda}_{\text{graded pieces: } + \text{ angular components, hence } \texttt{disc\_is\_square}} \;\prec\; \underbrace{q \text{ itself}}_{\text{the form}}
$$
The polygon is precisely the image of the Springer decomposition under the tropicalization of Lemma J.1 — it sees valuations and dimensions and forgets the residue square classes. This is the exact sense of the Bridge J section's "the Springer layers are the graded pieces of the valuation/tropical filtration"; it is the place-axis twin of the games-side identity (thermography in $\mathbb{T}_{\max}$; the sign mirror `MinPlus`↔`MaxPlus` is a convention flip, not a second semiring — `src/scalar/tropical.rs` already enforces the two-type separation).

**Witness tests (c).**
- ⟦implemented⟧ `src/forms/springer/local.rs::tests::{one_engine_decomposes_every_discrete_leg, unramified_qq_reads_extension_residue, residue_char_two_is_rejected_uniformly}` — the bucket engine, the extension-residue square class, and the char-2 boundary.
- ⟦proposed⟧ `polygon_is_the_springer_shadow`: diagonal $\langle a_i \rangle$ over `Qp<5,8>`, `Qq<3,3,2>`, `Laurent<Fp<7>,8>`; build $f_q = \prod (x - a_i)$ via `Poly`; assert the side multiset `{(root_valuation, length)}` equals `{(g.valuation, g.dim)}` from `springer_decompose_local`, and that grouping sides by slope parity reproduces `parity_layer(0)`/`parity_layer(1)` cardinalities (J.12(i), (iii)).
- ⟦proposed⟧ `initial_form_recovers_layer_discriminant`: compute $\mathrm{in}_\lambda(f_q)$ by the shift + `min_coeff_valuation` + reduce-at-min recipe; assert the product of its nonzero roots (equivalently $\pm$ its lowest nonvanishing coefficient ratio) has `is_square_finite::<K::Residue>` equal to the layer's `disc_is_square` (J.12(ii)).
- ⟦proposed⟧ `polygon_outlives_springer`: over `Qp<2,8>` (residue char 2) and `Gauss<Qp<5,6>>` (infinite residue field), `NewtonPolygon::of` succeeds while `springer_decompose_local` returns `None` — J.12(i)–(ii) need no Witt theory; only (iii) does.

---

## 4. Scope boundaries and non-claims

- **Discretely-valued legs only.** The surreal leg has 2-divisible value group: the second Springer layer collapses ($W(\mathrm{No}) = W(\mathbb{R})$, `springer/surreal.rs`) and there is no integer Newton lattice. Polygons over divisible $\Gamma$ are definable but are *not claimed or scheduled* — the same boundary the Springer engine already documents, and itself an instance of the local↔global symmetry.
- **Char-2 residue fields.** J.5/J.6/J.12(i)–(ii) hold for any residue characteristic; J.10/J.12(iii) require $\operatorname{char} k \ne 2$. The char-2 local Witt theory is the separate Aravire–Jacob layer (`springer/char2/`) and is outside Bridge J.
- **Precision.** On the capped-relative models (`Qp`/`Qq`/`Laurent`/`Ramified`/`Gauss`), valuations of *represented nonzero* elements are exact, so polygons of represented coefficients are exact; a coefficient whose true valuation exceeds the precision horizon renders as $0$ (vertex absent). J.1(ii) is truncation-safe; equality claims hold off the vanishing locus. The $\mathbb{F}_q(t)$ leg (Corollary J.9) is exact outright.
- **Choice of $\varpi$.** $\mathrm{ac}$, $\mathrm{in}_\lambda$, and $\partial_1$ depend on it; the code pins it to `Valued::uniformizer` via `residue_unit`. $\partial_0$ and the polygon do not.
- **No strictness claim** for "$v$ is a semiring homomorphism" (Remark J.2). No new theorem anywhere in this bridge: J is standard math made computational, the same status as shipped bridges A–I.

## 5. References

- T. A. Springer, *Quadratic forms over fields with a discrete valuation I*, Indag. Math. **17** (1955).
- T. Y. Lam, *Introduction to Quadratic Forms over Fields*, GSM 67, AMS, 2005 — Ch. VI (residue homomorphisms, Springer's theorem).
- N. Koblitz, *p-adic Numbers, p-adic Analysis, and Zeta-Functions*, GTM 58, Springer, 2nd ed. 1984 — Ch. IV (Newton polygons).
- J. Neukirch, *Algebraic Number Theory*, Grundlehren 322, Springer, 1999 — Ch. II (complete/henselian valued fields, unique extension of valuations).
- G. Dumas, *Sur quelques cas d'irréductibilité des polynômes à coefficients rationnels*, J. Math. Pures Appl., 1906 (polygon additivity; the irreducibility criterion).
- J.-P. Serre, *Local Fields*, GTM 67, Springer, 1979 — Ch. I (Eisenstein polynomials, total ramification).
- D. Maclagan, B. Sturmfels, *Introduction to Tropical Geometry*, GSM 161, AMS, 2015 — Ch. 2–3 (valuations as tropicalization; tropical roots/Kapranov in rank 1).
- O. Viro, *Hyperfields for tropical geometry I. Hyperfields and dequantization*, arXiv:1006.3034, 2010 (strict functoriality via the tropical hyperfield).
- H. Stichtenoth, *Algebraic Function Fields and Codes*, GTM 254, Springer, 2009 — Ch. 1 (places of $\mathbb{F}_q(t)$).

---

---

# Fourth wave — M, N, O (built)

The fourth-wave review asked where the **symmetry table** itself (README → "The
symmetries") was still uneven, rather than where a new number system could go. It
proposed three bridges — **M** (the Brown `ℤ/8` invariant, the char-2 cell of the
mod-8 spine), **N** (the unification pass), and **O** (lexicodes) — and all three
are now built and tested.

Claim-level discipline still applies: every item is **standard math made
computational**, the same status A–J shipped at — not a new theorem.

## Bridge M — the Brown invariant: the char-2 cell of the mod-8 spine

**Pillars:** `forms/char2/` (Arf) ↔ `forms/integral/discriminant/` (Milgram,
Bridge A) ↔ `forms/witt/brauer_wall.rs` (the mod-8 cycle).
**Claim level:** IMPLEMENTED AND TESTED — standard math (E. H. Brown, *Generalizations of
the Kervaire invariant*, Ann. of Math. **95** (1972); C. T. C. Wall, *Quadratic forms
on finite groups*, Topology **2** (1963); Milgram/van der Blij) made computational.

### The asymmetry it repairs

The mod-8 spine otherwise lives entirely on the char-0 side: the exact rational
signature, the genus oddity (`genus_signature_mod8`), the Milgram Gauss-sum phase
(`milgram_signature_mod8`, Bridge A), and the Weil `S` prefactor (Bridge I) are four
routes to `σ mod 8`. The char-2 side carried only the `ℤ/2` Arf bit. The classical
object filling the char-2 mod-8 cell is the **Brown invariant** of `ℤ/4`-valued
quadratic refinements.

### The mathematics

A `ℤ/4`-quadratic form `q : V → ℤ/4` on an `F₂`-space satisfies
`q(x+y) = q(x) + q(y) + 2·b(x,y)` with `b : V×V → F₂` symmetric (and `b_ii = q_i mod 2`,
so **not** alternating). For `b` nondegenerate the Gauss sum is a `ℤ[i]`-integer of
absolute value `2^{n/2}`,

```text
Σ_{x ∈ V} i^{q(x)} = 2^{n/2} · ζ₈^β,    ζ₈ = e^{2πi/8},
```

and `β ∈ ℤ/8` is the **Brown invariant**: additive under `⊥`, a complete invariant up
to split planes, making the Witt group of the category cyclic of order 8 generated by
`⟨1⟩` (`q(x)=1`). Three identifications make this the missing cell, not a fifth pillar:

1. **Arf is the 2-torsion.** Doubling a classical nonsingular char-2 form `q′ : V → F₂`
   gives `2q′ : V → ℤ/4` with `Σ (−1)^{q′} = (−1)^{Arf}·2^{n/2}`, so `β(2q′) = 4·Arf(q′)`
   — the shipped Arf bit embeds as `{0,4} ⊂ ℤ/8`.
2. **Milgram on the 2-elementary slice is Brown.** For an even lattice `L` with
   2-elementary `A_L`, `t ↦ 2t` identifies `(A_L, 2q_L)` with a `ℤ/4`-quadratic form
   whose Brown sum *is* the Milgram Gauss sum, so `β(2q_L) ≡ sign(L) (mod 8)` — computed
   from the **integer value-counts** `(n₀−n₂)+i(n₁−n₃)`, a **fifth route to `σ mod 8`**
   and the first with no floating point (the `GaussSum` route is `f64`).
3. **The generators are shipped lattices.** `a_n(1)` (`A₁`): `β = 1 ≡ σ`; `e_7()`:
   `β = 7 ≡ σ`; `d_n(4)`: `β = 4 ≡ σ`; the unimodular `e_8()`: `β = 0`.

### Implemented surface

- `forms/char2/brown.rs`
  - `brown_f2(n, q4: &[u128], bmat: &[u128]) -> BrownInvariants` — the `arf_f2` idiom with
    `q4` (values mod 4) replacing the diagonal; `bmat` is the **off-diagonal** symmetric
    polar (the diagonal `b_ii = q4[i] mod 2` is derived). `BrownInvariants { beta, rank,
    radical_dim, radical_anisotropic }` mirrors `ArfInvariants` field-for-field.
  - **Reduction route** with exact-integer oracles: split off `rad(b)` (`q|rad` is
    linear into `{0,2}`, so `Σ_V` factors), then reduce the nonsingular core into odd
    lines (`β = 1`/`7`) and even planes (`β = 0`/`4`) and add the phases in `ℤ/8`.
    An anisotropic radical vanishes the full sum; `beta` still reports the core. The
    old direct Gauss-sum enumeration is retained as a test-only oracle through the
    former `rank ≤ 26` budget edge.
  - `double_f2(qd, bmat)` — the `q′ ↦ 2q′` embedding from `arf_f2` input data.
- `forms/integral/discriminant/`
  - `DiscriminantForm::brown_invariant(&self) -> Option<BrownInvariants>` — `Some` only for
    **2-elementary** `A_L` (read off the invariant factors), enumerating `A_L` directly
    via `quadratic_value_mod2`. `b_L` is nondegenerate on `A_L`, so this slice has no
    radical and `β ≡ sign(L) mod 8`.
  - `DiscriminantForm::fqm_gauss_phase(&self) -> Option<FqmGaussPhase>` — the
    p-primary Milgram/Brown phase projection over all represented discriminant groups,
    with `FqmPrimaryPhase { prime, order, exponent, phase_mod8 }` and total
    `phase_mod8`. This extends the phase computation past the 2-elementary Brown slice
    (`A_3`, `E_6`, mixed-primary sums, ...).
- `forms/integral/fqm_witt.rs`
  - `FiniteQuadraticModule::new` / `::cyclic` / `::direct_sum` — native finite
    quadratic modules in a cyclic product presentation, validated for nonsingularity
    and the quadratic law.
  - `FiniteQuadraticModule::witt_class` and `DiscriminantForm::fqm_witt_class` —
    p-primary Wall/Nikulin Witt classes, represented by canonical anisotropic cores
    after quotienting isotropic cyclic subgroups.
  - `DiscriminantForm::is_fqm_witt_equivalent` — equality in the FQM Witt group.

### Oracles / implemented tests

- `double_f2(q′).beta == 4 * arf_f2(q′).arf` across nonsingular metrics; doubled forms
  land in `{0,4}`.
- The generators `⟨1⟩ → β=1`, `⟨−1⟩ → β=7`, and the order-8 relation `⟨1⟩^{⊥8} → β=0`;
  the split objects (the even hyperbolic plane and `⟨1⟩ ⊥ ⟨−1⟩`) have `β=0`; additivity
  under `⊥` across a spread of components; anisotropic-radical detection.
- `brown_invariant` of `a_n(1)`/`e_7()`/`d_n(4)`/`d_n(8)`/`e_8()` gives `β ≡ sign mod 8`,
  cross-checked against `fqm_gauss_phase` and the shipped f64 `milgram_signature_mod8`;
  non-2-elementary forms (`a_n(2)`, `a_n(3)`, `e_6()`) return `None` for Brown but still
  have FQM phases.
- `fqm_gauss_phase` reports primary factors on `A_1 ⊕ A_2`, extends the 2-primary phase
  to `A_3` (`Z/4`), covers odd torsion such as `E_6` (`Z/3`), and matches the exact
  signature, genus oddity route, and legacy float oracle across the ADE zoo.
- `fqm_witt_class` agrees with the native cyclic presentation on `A_1`, separates
  `A_1` from `E_7`, reduces `⟨1/2⟩ ⊕ ⟨3/2⟩` and `A_2 ⊕ E_6` to the trivial Witt class,
  and keeps `fqm_gauss_phase` as its phase projection.

### Scope / caveats

- **Category trap (load-bearing):** Brown's `b` is symmetric-not-alternating with
  `b_ii = q_i mod 2`, **not** the engine's alternating char-2 polar; `double_f2` is the
  only bridge between the two categories. Kept distinct from the graded
  `BrauerWallClass`/Arf exactly as Bridge F insists for its Brauer class.
- The Brown lattice tie is **2-elementary discriminant groups only**; higher 2-power and
  odd torsion have the `FqmGaussPhase` Milgram/Brown **phase projection**, and the
  separate `FqmWittClass` surface carries the full finite-quadratic-module Witt class
  within its explicit finite enumeration budget.
- No new theorem: Brown 1972 is the source; the bridge is the wiring to Arf (shipped)
  and Milgram (Bridge A).

## Bridge N — the unification pass: four joins of already-shipped parts

**Pillars:** vary per item — each joins surfaces that already exist. **Claim level:**
IMPLEMENTED AND TESTED — standard math; each item is assembly + verification of
shipped machinery, deliberately smaller than a headline bridge.

### N.1 — Milnor's exact sequence: the Springer residues go global

**Pillars:** `forms/springer/` ↔ `forms/witt/` ↔ the integral pillar's signature.
The Witt-group-level statement of the local residue engine:

```text
0 → W(ℤ) → W(ℚ) →∂ ⊕_p W(F_p) → 0     (exact; Milnor–Husemoller Ch. IV; Lam GSM 67 Ch. IX)
```

`forms/witt/milnor.rs::global_residues(entries: &[i128]) -> Option<(i128,
BTreeMap<u128, WittClassG>)>` returns the **signature** (`W(ℤ) ≅ ℤ`, the kernel) and
the nonzero residues `∂_p`. For odd `p`, these are the second Springer residues,
computed exactly from the `i128` entries (`v_p` + Legendre + the signed-discriminant
square class, matching the `finite_odd_witt` convention) so `p` stays runtime while
`Fp<P>` is const-generic. For `p=2`, Milnor's hand-defined boundary contributes the
parity of diagonal lines with odd dyadic valuation, represented in the existing
`W(F_2) ≅ Z/2` carrier `WittClassG::Char2 { field_degree: 1, arf }`.

The equal-characteristic twin now ships too. For odd constant fields, the split
affine-line form

```text
W(F_q(t)) ≅ W(F_q) ⊕ ⊕_π W(F_q[t]/π)
```

is exposed as `forms/witt/milnor.rs::global_residues_ff(entries:
&[RationalFunction<S>]) -> Option<FunctionFieldMilnorResidues<S>>`.
The first component is the `W(F_q)` class selected by the even-valuation layer at
the degree place `∞`; the vector contains the nonzero second residues at finite
monic irreducible places. It reuses the exact `F_q(t)` place helpers
(`try_valuation_at_ff`, `try_residue_unit_at`, `try_chi_kappa`) rather than the
capped local-field models.

- **Oracles:** finite support (`∂_p = 0` for `p ∤ ∏aᵢ`, plus zero dyadic parity);
  square/hyperbolic invariance of `(signature, residues)`; residues distinguish
  `⟨1⟩` from `⟨3⟩` and `⟨1⟩` from `⟨2⟩`, cross-checked against the shipped
  Hasse–Minkowski `try_is_isotropic_q`; `∂₅` matches an independent computation
  through `springer_decompose_qp` on the capped `Q₅` model; `⟨2⟩`/`⟨1,2⟩`/`⟨−2⟩`
  pin the dyadic cell; and the function-field leg pins constants, the `t` place,
  nonsquare constants, a degree-2 finite place, square-multiple invariance,
  hyperbolic cancellation, and radical-entry rejection.
- **Boundary:** `global_residues_ff` is odd-characteristic only (`FiniteOddField`).
  Characteristic-2 function fields keep using the separate Artin-Schreier /
  Aravire-Jacob layer; tame and wild norm-residue symbols are Bridge K follow-ons,
  not part of this Witt-residue map.

### N.2 — the Scharlau transfer, named

**Pillars:** `scalar/extension` (`CyclicGaloisExtension`) ↔ `forms/trace_form`. The
existing `trace_twisted_form::<E>(0)` is `s_*(⟨1⟩)` for the transfer `s_* : W(E) →
W(F)` along `Tr_{E/F}` (Lam GSM 67 Ch. VII; Scharlau Ch. 2). New
`trace_form::transfer_diagonal<E: CyclicGaloisExtension>(entries: &[E]) ->
Metric<E::Base>` builds `s_*(⟨λ₁,…,λᵣ⟩) = ⟂ᵢ (x,y) ↦ Tr(λᵢ·x·y)` through the shipped
`assemble_twisted_form` core.

- **Oracles:** the `k=0` twisted form equals `transfer_diagonal(&[1])`; the transfer of
  a hyperbolic form splits; **Frobenius reciprocity** `s_*(r*(x)·y) = x·s_*(y)` (the
  form-level `Tr(c·λ·z) = c·Tr(λ·z)`); and **Springer's odd-degree theorem** —
  restriction `r*` is injective for odd `[E:F]`, witnessed by `⟨1,1⟩` staying
  anisotropic from `F₃` to `F₂₇`.
- **Boundary:** char ≠ 2 (the `Tr(x·σ(x)) = 2N = 0` trap the module documents); the
  char-2 transfer is the Artin–Schreier route in `function_field_char2.rs`.

### N.3 — Nikulin: genus ⟺ (signature, discriminant form)

**Pillars:** `forms/integral/genus` ↔ `forms/integral/discriminant`. Nikulin's
criterion (Izv. Akad. Nauk SSSR **43** (1979), Cor. 1.9.4) upgrades the mod-8 phase
comparison of Bridges A/I to a classification equivalence: two **even** lattices
share a genus iff they have equal signature pairs and isomorphic discriminant
quadratic forms. The missing piece — `DiscriminantForm::is_isomorphic(&self, other)
-> Option<bool>` — matches invariant factors, then runs a **budgeted** homomorphism-
extension search (minimal generators by maximal order → image assignment pruned by
order and `q`-value → BFS extension → `q`-preservation on every element), mirroring
`automorphism_group_order_bounded`'s `None`-past-budget pattern.

- **Oracles:** `are_in_same_genus(a,b) == (equal signatures ∧ q_a ≅ q_b)` across the
  zoo (`a_n`, `d_n`, `e_6/7/8`, sums), pinned by the **Milnor pair** (`E₈⊕E₈` vs
  `D₁₆⁺`: same genus, non-isometric, both trivial disc form) and easy separations
  (`A₂`: ℤ/3 vs `A₁⊕A₁`: (ℤ/2)²). `q`-sensitivity is pinned directly: `A₁` and `E₇`
  share the group ℤ/2 but have `q`-values `1/2` vs `3/2` and are **not** isomorphic.
- **Boundary:** even lattices only (the `from_lattice` boundary); the brute-force
  budget is honest (`None` past `ISO_GROUP_CAP`/node budget) — a cross-check of two
  shipped routes, not a p-adic-symbol reimplementation.

### N.4 — one Bernoulli source for Eisenstein and mass

**Pillars:** `forms/integral/mass_formula` ↔ `forms/integral/modular`. The mass
constants and the Eisenstein constants `240 = −8/B₄`, `−504 = −12/B₆` are the same
Bernoulli numbers. The Akiyama–Tanigawa helper in `mass_formula.rs` is now the shared
`pub(crate) bernoulli` source; `modular.rs::eisenstein_e4/e6` derive their constants
from it via `c_{2k} = −4k/B_{2k}`, with the literals kept as the pinned oracle
(TABLES.md discipline: derived value asserted equal to curated constant).

- **Oracles:** `eisenstein_constant(2) == 240`, `eisenstein_constant(3) == −504`; the
  von Staudt–Clausen denominators `B₂…B₈` as a free check; `mass_even_unimodular(8)`
  through the shared helper still `= 1/E8_WEYL_GROUP_ORDER`.

## Bridge O — lexicodes: greedy = mex, the games ↔ integral edge

**Pillars:** `games/` (mex) ↔ `forms/integral/codes` (Bridge H) → Construction A /
theta (Bridges H/E). **Claim level:** IMPLEMENTED AND TESTED — standard math
(Conway–Sloane, *Lexicographic codes…*, IEEE Trans. Inform. Theory **32** (1986)
337–348). Closes the one pillar edge the bridge graph still lacked: games ↔ integral.

The lexicode `L(n,d)` greedily keeps every vector at Hamming distance `≥ d` from those
kept so far; Conway–Sloane prove the result is **linear** by Sprague–Grundy theory.
`games/lexicode.rs` ships two routes:

- `lexicode_naive(n,d)` — the literal greedy scan for small `n`, **discover-don't-
  assert**: collect greedily, verify XOR-closure, `None` on a closure failure (which
  would *falsify* linearity rather than hide it).
- `lexicode(n,d)` — the production route, carrying the full distance array
  `dist[v] = d(v,C)` and updating it in one `O(2ⁿ)` pass per generator via the coset
  recurrence `d(v, C ∪ (g⊕C)) = min(d(v,C), d(v⊕g,C))` with a monotone cursor (so the
  `n=24` build is fast), budgeted by `LEXICODE_NODE_BUDGET`.
- `nim_lexicode_naive(2^k,n,d)` (spelled by exponent `k`) — the literal greedy scan
  over the nim alphabet `{0,...,2^k-1}`, returning `NimLexicode` after verifying
  coordinatewise nim-addition closure. `NimLexicode::is_closed_under_nim_scalars`
  asks the stronger field-linearity question by multiplying coordinates with finite
  nim multiplication.

The greedy step is shown to be `mex(Forbidden)` (the union of radius-`(d−1)` balls)
via [`grundy::mex`] and a toy-`n` witness; the deeper Conway–Sloane turning-game
realization is cited for transcription in a formalization pass, **subordinate to
`OPEN.md` §1** (the solved degree-1 shadow, not progress on the open question).

- **Oracles:** `lexicode_naive == lexicode` (n ≤ 12); `d=1 → F₂ⁿ`, `d=2 → even-weight`;
  `lexicode(7,3)`/`lexicode(8,4)` reproduce the Hamming weight enumerators;
  `lexicode(24,8)` is `[24,12,8]` doubly-even self-dual with the **Golay** weight
  enumerator (uniqueness of the Type II `[24,12,8]` code closes "is Golay"); and the
  chain rung `lexicode(24,8).construction_a()` is even unimodular rank 24 **with**
  roots — re-pinning Bridge H's Golay ≠ Leech boundary from the games side. The
  q-ary/nim route checks repetition lexicodes over bases `4`, `8`, and `16`: all are
  nim-additive, bases `4` and `16` are nim-scalar closed, and base `8` is not.
- **Scope:** the optimized production route remains binary. The base-`2^k` route is
  literal and budgeted, with lexicographic order = standard digit order (coordinate 0
  the most significant digit); a permuted coordinate order gives an equivalent code.
  The deeper Conway-Sloane turning-game realization is still cited for transcription
  in a formalization pass.

---

# Fifth wave — Bridge K (built)

Bridge K was the last unbuilt **non-deferred** bridge — the natural completion of the
Brauer thread. It lifts the shipped 2-torsion rational Brauer surface (`adelic.rs`,
Bridge F) to the **full local Brauer group** `Br(K_v) ≅ ℚ/ℤ`, via the cyclic-algebra
invariant of local class field theory, built from the Galois data Bridge C already
exposes (`CyclicGaloisExtension`). Standard math made computational — not a new theorem.

## Bridge K — cyclic algebras: the full `ℚ/ℤ` Brauer invariant

**Pillars:** `scalar/extension` (`CyclicGaloisExtension`) ↔ a new ungraded Brauer class
in `forms/witt/cyclic.rs` ↔ `forms/local_global/{adelic,function_field}` (the
reciprocity sequence) ↔ `forms/trace_form` (the degree-2 norm-form oracle).
**Claim level:** IMPLEMENTED AND TESTED — standard math (Serre, *Local Fields*, Ch. XII;
Gille–Szamuely §§6.3–6.4; Reiner §§31–32; Tate in Cassels–Fröhlich Ch. VII). Lifts the
**2-torsion** Brauer surface to the full **`Br(K_v) = ℚ/ℤ`** image; Bridge F's rational
Clifford invariant sits inside as the `½`-slice.

### The mathematics

A cyclic extension `E/K` of degree `n` with distinguished generator `σ` and `a ∈ K*`
defines the cyclic algebra `(χ_σ, a) = ⊕_{i<n} E·uⁱ`, `uⁿ = a`, `u·x = σ(x)·u`. Over a
local field with `E/K` **unramified** and `σ` the arithmetic Frobenius (the convention
every `CyclicGaloisExtension::sigma` uses), the class-field-theory invariant map gives

```text
inv_K[(χ_σ, a)] = v(a)/n  (mod ℤ) ∈ (1/n)ℤ/ℤ ⊂ ℚ/ℤ
```

— the **full** local Brauer group, not just its 2-torsion. The value reads only `v(a)`
and `n`; `σ` fixes the sign convention (`χ_σ(σ) = +1/n`). The quaternion case `n=2`
reproduces the shipped `brauer_local_invariants` place-by-place. Globally the
Albert–Brauer–Hasse–Noether sequence `0 → Br(K) → ⊕_v Br(K_v) → ℚ/ℤ → 0` gives
`∑_v inv_v ≡ 0`. Over `ℚ`, Minkowski forces every cyclic `E/ℚ` of degree `>1` to ramify,
so `n>2` reciprocity over `ℚ` would need ramified symbols — out of scope; the clean route
is `F_q(t)`, where the **constant extension** `F_{qⁿ}(t)` is unramified at *every* place
with `Frob_v = σ^{deg v}`, so `inv_v = deg(v)·v(a)/n` and `∑_v inv_v = deg(div a)/n = 0` —
full `ℚ/ℤ` reciprocity reduced to the product formula the function-field layer embodies.

### Implemented surface

- `forms/witt/cyclic.rs`
  - `BrauerClass` (private `local: BTreeMap<Place, Rational>`; values in `[0,1)`, zeros
    omitted) with `add` (entrywise mod ℤ), `invariant_sum`, `local`/`local_invariant`,
    `from_local`, `split`/`is_split`, and the Bridge F embedding
    `from_two_torsion(&Brauer2Class)` / `two_torsion() -> Option<BTreeSet<Place>>` (the
    `½`-slice and its inverse). `Place` already derives `Ord` (Bridge F shipped it).
  - `cyclic_algebra_invariant::<E: CyclicGaloisExtension>(a: &E::Base) -> Option<Rational>`
    where `E::Base: Valued` — `v(a)/n mod ℤ` for the unramified local class. Monomorphized
    at `E = Qq<P,N,F>` over `Q_p = Qq<P,N,1>` (the only `CyclicGaloisExtension` with a
    `Valued` base); exact even over the capped model (reads only the valuation), `None` on
    `a=0` / precision loss — never a wrong value.
- `forms/trace_form.rs`
  - `cyclic_algebra_trace_form::<E: CyclicGaloisExtension>(a: &E::Base) -> Metric<E::Base>`
    — the literal cyclic-algebra trace form `T_A(z)=Trd_A(z²)` for
    `A=(E/F,σ,a)=⊕E·u^i`, in the `E·u^i` line basis. The `u^0` and, for even degree,
    `u^{n/2}` self-blocks reuse `assemble_twisted_form`; the `i`/`n-i` line pairs are
    pure polar blocks.
- `forms/local_global/function_field.rs`
  - `constant_extension_invariants::<S: FiniteOddField>(n, a) -> Option<Vec<(FFPlace<S>, Rational)>>`
    — `inv_v = deg(v)·v(a)/n mod ℤ`, the exact full-`ℚ/ℤ` reciprocity oracle (everywhere
    unramified, no ramified symbols). A `Vec` since `FFPlace` is not `Ord`. Plus
    `constant_extension_invariant_sum` (`∑_v inv_v = 0`).

### Oracles / implemented tests

- **Degree-2 compatibility (the lift is a lift):** `cyclic_algebra_invariant::<Qq<5,4,2>>`
  matches the shipped `brauer_local_invariants(d, a)` at `Prime(5)` (`d=2`, a nonsquare
  unit) across `v_5(a) ∈ {0,1,2,3}`.
- **Splitting law:** `inv = 0 ⇔ n ∣ v(a)`; `(χ_σ, N_{E/K}(x))` splits (real norms over
  `Qq<3,3,2>` via the shipped `FieldExtension::norm`).
- **Image / additivity / convention:** over `n=3` the image is the full `(1/3)ℤ/ℤ`
  (`inv(a)=1/3`, `inv(a²)=2/3` — pinning `+v/n` against the negated convention), with
  additivity and `n·inv ≡ 0`.
- **Full-strength reciprocity over `F_q(t)`:** `n ∈ {2,3,4,5}`, `Σ inv_v ≡ 0`, with the
  independent `deg(div a) = 0` check; a degree-2 place carries `deg(v)=2` (e.g. `2/3` at
  `t²+2`), a value the 2-torsion surface cannot see.
- **Bridge F embedding:** `from_two_torsion`/`two_torsion` round-trip and additivity (XOR
  of ramification sets ↦ add of `½`-slices); the shipped quaternion reciprocity re-read
  through `BrauerClass::invariant_sum() = 0`.
- **Degree-2 norm-form oracle (§6 trace-form tie):** the cyclic class `(χ_σ,a) = (−1,a)_ℚ`
  over `E = ℚ(i)` splits at `v` ⇔ its reduced-norm form `⟨1,1,−a,−a⟩` (built from
  `trace_twisted_form::<Surcomplex<Rational>>(1) = ⟨2,2⟩`) is isotropic over `ℚ_v` ⇔
  `inv_v = 0` — tying the invariant to the shipped Hasse–Minkowski layer.
- **Cyclic trace-form oracle (§6(c)):** for `A=(ℚ(i)/ℚ, conjugation, a)`,
  `cyclic_algebra_trace_form` gives the literal `Trd(z²)` form
  `⟨2,-2,2a,2a⟩`, not the reduced norm. The test pins the honest degree-2 relation
  `Trd(z²)=Trd(z)^2-2·Nrd(z)` pointwise against `Nrd=⟨1,1,-a,-a⟩`, and checks over
  `F_27/F_3` that the `u`/`u²` cross-pair block is Witt-hyperbolic.

### Scope / caveats

- **Unramified-at-`v` only** for `v(a)/n`; ramified local symbols are out of scope (the
  `F_q(t)` route delivers full `ℚ/ℤ` strength without them). Reads only `v(a)`, `n`,
  `deg(v)`, so exact even over the capped-precision local models.
- **Ungraded** Brauer group — kept strictly distinct from the graded `BrauerWallClass`,
  exactly as Bridge F insists.
- **Finite legs carry no Brauer content** (Wedderburn): over `Nimber`/`Fpn` every central
  simple algebra splits, so the Gold forms have no `inv` (their classifier is
  Arf/Brauer–Wall, Bridge B). K lives only on the local/global legs (`Qq`, `F_q(t)`, and
  the real place via the 2-torsion embedding).
- `cyclic_algebra_trace_form` is **not** the reduced norm for general `n` (and is not
  equal to it for quaternions); it is the quadratic trace companion `Trd(z²)`. The
  degree-2 reduced-norm identity remains the separate high-value tie above.

---

## Bridge K — formal statements (formalization-pass appendix)

> Moved here from `roadmap/TODO.md` on building K. Standard math made computational; the
> theorems below are LCFT, the surface that realizes them shipped as in the section above.

**Status:** IMPLEMENTED AND TESTED. Every theorem is **standard math** (local/global class
field theory); the bridge made it computational on surfaces the crate already ships.

### 1. The cyclic algebra *(standard math)*

For a cyclic Galois `E/K` of degree `n` with generator `σ` and character
`χ_σ : Gal(E/K) → (1/n)ℤ/ℤ`, `χ_σ(σ) = 1/n`, and `a ∈ K*`, the **cyclic algebra**
`(χ_σ, a) = ⊕_{i<n} E·uⁱ`, `uⁿ = a`, `u·x = σ(x)·u` is central simple of degree `n`,
containing `E` as a maximal subfield (Gille–Szamuely, Ch. 2):

- `(χ_σ, a) ⊗ (χ_σ, b) ∼ (χ_σ, ab)` in `Br(K)`;
- `(χ_σ, a)` splits `⟺ a ∈ N_{E/K}(E*)`; in particular `(χ_σ, N_{E/K}(x))` splits;
- `a ↦ [(χ_σ, a)]` induces `K*/N_{E/K}(E*) ≅ Br(E/K)`;
- `n = 2`, char ≠ 2: `(χ_σ, a)` **is** the quaternion `(d, a)_K` for `E = K(√d)`; char 2:
  the Artin–Schreier symbol `[d, a)` already in `function_field_char2.rs`.

`CyclicGaloisExtension` carries exactly this data: `basis()`, `sigma()`/`sigma_power(k)`,
`FieldExtension::{trace, norm, extension_degree}`.

### 2. The local invariant *(standard math)*

For `K` nonarchimedean local, `E/K` **unramified** of degree `n`, `σ` the arithmetic
Frobenius, the invariant isomorphism `inv_K : Br(K) ≅ ℚ/ℤ` satisfies

```text
inv_K[(χ_σ, a)] = v(a)/n   (mod ℤ),
```

and every class arises this way (every CSA over a local field has an unramified splitting
field). References: Serre, *Local Fields* Ch. XII; Gille–Szamuely §6.3–6.4; Reiner §31.
Consequences: `(χ_σ, a)` splits at `K` iff `n ∣ v(a)`; the image is the full cyclic group
`(1/n)ℤ/ℤ`, not just its 2-torsion.

**Convention warning.** The sign of `inv` depends on the *arithmetic* Frobenius and
`χ_σ(σ) = +1/n`; the geometric-Frobenius convention negates it. Every `sigma()` impl
(`Fpn::frobenius`, the Witt–Frobenius on `Qq`, nim-squaring on `Nimber`) is arithmetic, so
`+v(a)/n` is the consistent choice. Reciprocity (§3) is convention-independent; degree-2
compatibility (§4) is not — pinned by the `n=3` asymmetric test (`inv(a²)=2/3 ≠ 1/3`).

**Archimedean place.** `Br(ℝ) = ½ℤ/ℤ`; for `E = ℂ`, `σ =` conjugation,
`inv_ℝ[(χ_σ, a)] = ½ iff a < 0`. No valuation to read — special-cased exactly as the
shipped quaternion route, and entered through the 2-torsion `from_two_torsion` embedding.
`Br(ℂ) = 0`.

**Ramified caveat (load-bearing).** If `E/K_v` is ramified, `v(a)/n` is **not** the
invariant; the general local symbol is needed. The shipped surface is scoped to
unramified-at-`v` data, which suffices for everything below.

### 3. Global reciprocity *(standard math)*

For a global field `K` the Albert–Brauer–Hasse–Noether sequence
`0 → Br(K) → ⊕_v Br(K_v) → ℚ/ℤ → 0` (Reiner §32; Tate in Cassels–Fröhlich Ch. VII) gives
`∑_v inv_v(A ⊗ K_v) ≡ 0 (mod ℤ)`, finitely supported. For a global cyclic `(χ_σ, a)` and
`v` unramified with `Frob_v = σ^{m_v}`, the local term is `inv_v = m_v·v(a)/n`.

**Scope fact, not a gap:** over `ℚ`, Minkowski ⇒ every cyclic `E/ℚ` of degree `>1`
ramifies, so a full-strength `n>2` reciprocity test over `ℚ` needs ramified symbols. The
crate uses the clean alternative `K = F_q(t)`: the **constant extension** `F_{qⁿ}(t)` is
unramified at *every* place (incl. `∞`), `Frob_v = σ^{deg v}`, so
`∑_v inv_v = (1/n)·∑_v deg(v)·v(a) = (1/n)·deg(div a) = 0` — full `ℚ/ℤ` reciprocity reduced
to "principal divisors have degree 0", the product formula already shipped
(`constant_extension_invariant_sum`). (`Br(F_q(t))` via residues: Faddeev, Gille–Szamuely
§6.4, using `Br(F_q) = 0`.)

### 4. The degree-2 lift of the shipped 2-torsion surface

Quaternions are the `n=2` cyclic algebras. For `p` odd and `d` a nonsquare unit at `p`,
`E = ℚ_p(√d)` is the unramified quadratic and
`inv_p[(χ_σ, a)] = v_p(a)/2 ≡ ½·[v_p(a) odd]`, while `(d,a)_p = (d/p)^{v_p(a)} =
(−1)^{v_p(a)}`, so the degree-2 cyclic invariant reproduces the shipped quaternion
`brauer_local_invariants` place-by-place (test 1). The new class type replaces "a set of
ramified places" by "a `ℚ/ℤ`-valued divisor of places", with the shipped surface as its
`{0, ½}` slice.

### 5. Bridge F as the 2-torsion part

`Brauer2Class` (a `BTreeSet<Place>`, symmetric-difference addition) embeds via
`from_two_torsion`: `v ↦ ½·[v ∈ ramified]`, a group monomorphism onto the 2-torsion of
`⊕_v ℚ/ℤ`. Quadratic-form Brauer classes are 2-torsion, so all of Bridge F (Hasse–Witt
`s(q)`, the even-Clifford class `c(q)`, the Lam `n mod 8`/disc correction) lands inside the
`BrauerClass` type; K supplies the full `ℚ/ℤ` ambient and the `n>2` classes F cannot see.
One ambient group, two constructors. Reciprocity restricted to the `½`-slice is
"`|ramified|` even". Kept **ungraded**, strictly distinct from `BrauerWallClass`.

### 6. The tie to `trace_form.rs` *(standard math)*

The honest statements behind the one-line gloss:

**(a) `n=2`, char ≠ 2.** `Nrd(x + yu) = N_{E/K}(x) − a·N_{E/K}(y)`. Since
`Q_1(x) := Tr(x·σ(x)) = 2·N_{E/K}(x)`, `Nrd ≅ ½Q_1 ⊥ (−a/2)Q_1`. With
`trace_twisted_form::<Surcomplex<Rational>>(1) = ⟨2,2⟩`, `Nrd[(−1,a)_ℚ] = ⟨1,1,−a,−a⟩`, and
`(χ_σ, a)` splits at `v` iff this form is isotropic over `K_v` iff `inv_v = 0` — the
**shipped degree-2 norm-form oracle** (test 6), tying `inv` to the Hasse–Minkowski layer.

**(b) `n=2`, char 2.** `Q_1 = Tr(xσ(x)) = 2N = 0` and `Tr(x²)` has vanishing polar — the
char-2 trap. The reduced-norm form of `[d, a)` is the 2-fold Pfister `[1,d] ⊥ a[1,d]`,
**already implemented** in `function_field_char2.rs` (Schmid's residue formula); that layer
*is* the char-2, `n=2` instance of Bridge K, shipped.

**(c) General `n`.** `Nrd` is degree `n`; the quadratic companion is `T_A(z) = Trd(z²)`.
Since reduced trace sees only the `u^0` coefficient, `T_A` decomposes over the line
pairs `Eu^i`/`Eu^{n-i}`: the `u^0` block, and the `u^{n/2}` block when `n` is even,
are `assemble_twisted_form` instances, while the remaining pairs are pure polar blocks.
This is now shipped as `cyclic_algebra_trace_form`; for `n=2` it satisfies
`Trd(z²)=Trd(z)^2-2·Nrd(z)`, so it is adjacent to but not equal to the reduced norm.

**(d) Non-tie, for honesty.** Over `Nimber`/`Fpn` every CSA splits (Wedderburn), so the
Gold forms carry **no** `inv`; their classifier is Arf/Brauer–Wall (Bridge B). K lives only
on the local/global legs (`Qq`, `F_q(t)`, `ℝ`).

### 7–9. Surface, tests, scope

As built — see "Implemented surface", "Oracles / implemented tests", and "Scope / caveats"
in the section above. References: §§1–3, 6 standard math (Serre, *Local Fields* Ch. XII,
XIV §5; Gille–Szamuely Ch. 2, §§6.3–6.4, §9.2; Reiner §§31–32; Tate in Cassels–Fröhlich
Ch. VII; Lam, GSM 67, Ch. III, V). No interpretation- or open-level claims are introduced.

---

## DONE — status snapshot

Implemented and tested in the Rust core:

- **First wave (A–D):** lattice/Clifford/Brauer–Wall via Milgram's Gauss sum (A);
  char-2 Arf over the `Fpn<2,N>` fields (B); Frobenius as an outermorphism (C);
  transfinite char-2 Clifford `OrdinalAlgebra` on the checked tower (D).
- **Second wave (E/F/H/I):** theta/modular forms and the Milnor isospectral pair (E);
  Construction A codes↔lattices with MacWilliams↔theta (H); the discriminant-form
  Weil representation (I); the rational Brauer/Clifford invariant correction (F).
- **Third wave (J):** the valuation as tropicalization plus Newton polygons, with the
  slope ⟺ Springer-residue-layer cross-check; formal proofs in the appendix above.
- **Fourth wave (M, N, O):** the Brown `ℤ/8` invariant — the char-2 cell of the
  mod-8 spine, `β = 4·Arf` and `β ≡ sign mod 8` on 2-elementary discriminant forms
  (M); the unification pass — Milnor's global residues (N.1), the Scharlau transfer
  (N.2), Nikulin's genus criterion (N.3), one Bernoulli source (N.4); and lexicodes
  (greedy = mex, the `[24,12,8]` lexicode is Golay; O).
- **Fifth wave (K):** the full `ℚ/ℤ` cyclic-algebra Brauer invariant — `BrauerClass`
  and `cyclic_algebra_invariant` (`v(a)/n`, the unramified local class) over the `Qq`
  leg, `constant_extension_invariants` (full-strength reciprocity over `F_q(t)`), and
  the degree-2 norm-form oracle; Bridge F embeds as the `½`-slice.

Buildable work and the deferred stars (`*1` spinor genus, `*2` Drinfeld/Carlitz,
`*4` the wild local symbol) live in `roadmap/TODO.md` — the game-valued ledger;
newly completed work goes in the `roadmap/DONE.md` ledger; the genuine open
problems stay in `OPEN.md`, loopy-valued: `tis`/`tisn`, `on`/`off`, `over`/`under`
(the old numerals §1–§4 survive as aliases).
