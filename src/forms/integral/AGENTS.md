# AGENTS.md — `src/forms/integral/`

The **arithmetic view** of the forms pillar. The classifiers in the parent work
over a *field* (square classes / Witt / Arf). An **integral lattice** is the
complementary object: a free ℤ-module with an integer Gram matrix. Its invariants
are arithmetic (det, level, minimum, kissing number, |Aut|), and the coarse
classification is the **genus** (local equivalence at every place), which reuses
the `local_global/padic.rs` + `local_global/adelic.rs` primitives.

`mod.rs` re-exports flat. This is the visible meeting point of the char-0 mod-8
spine (`BW(ℝ)=ℤ/8`, Bott, the 8-fold table) and the lattice world — `E₈` is the
unique rank-8 even unimodular lattice. Convention: **norm** `Q(x) = xᵀGx` (a
"norm-2 root" has `Q=2`).

- **`lattice.rs`** — `IntegralForm { gram: Vec<Vec<i128>> }` (private Gram, built via
  `new` (square+symmetric-checked) / `diagonal`, never a struct literal).
  `determinant` (fraction-free **Bareiss**, exact), `is_even`/`is_unimodular`,
  `is_positive_definite` (Sylvester leading-minors via Bareiss), `signature` (exact
  rational diagonalization), `invariant_factors` (SNF → discriminant group `L#/L`),
  `level` (smallest `N` with `N·G⁻¹` even-integral, via the exact `Rational` inverse),
  `clifford_metric` (rational Clifford metric), `clifford_metric_f2` (even-lattice
  mod-2 char-2 metric), `direct_sum`. The positive-definite geometry: `short_vectors`
  (unimodular size-reduction, then **Fincke–Pohst**: float LDLᵀ bounds the box, exact
  i128 norm filters the leaves, vectors mapped back to the original basis — float error
  can't add/drop a vector), `minimum`/`minimal_vectors`/`kissing_number`, and
  `automorphism_group_order` (closed-form diagonal/ADE/root-system fast paths, else
  backtracking over basis-vector images — every complete assignment is an
  automorphism, so the count is exact). **Looks like a bug, isn't:** (a) the geometry
  methods return `None` for indefinite lattices on purpose (infinitely many vectors of
  each norm); (b) |Aut| is bounded by an explicit node budget (`AUTO_NODE_BUDGET`) and
  returns `None` past it (`automorphism_group_order_bounded` exposes the budget) — an
  honest `None`, not silent truncation; (c) `level(⟨1⟩)=2`, not 1 — `ℤ` is odd. Oracles:
  `A_2`/`A_3`/`D_4`/`E_8` det, kissing (6/12/24/240), |Aut| (12/48/1152), level (3/·/·/1),
  `Z^n` (|Aut| `2ⁿ·n!`).
- **`diagonal.rs`** — `pub(crate)` exact-rational diagonalization helpers shared by
  `lattice`, `genus`, and `discriminant` (signature, Sylvester minors, p-adic
  Gram–Schmidt). Not a public surface.
- **`root_lattices.rs`** — the ADE catalogue: `a_n` (Cartan matrix), `d_n` (`B·Bᵀ`
  from the geometric basis `{eᵢ−e_{i+1}}∪{e_{n-2}+e_{n-1}}`), `e_6`/`e_7`/`e_8` (Dynkin
  edge lists). `coxeter_number = #roots/rank` (computed). `is_root_lattice` (min 2 +
  roots generate `L`, index off the HNF pivots). Det/kissing/Coxeter oracles protect
  every construction; |Aut| oracles include `A_n`→`2(n+1)!` (n≥2; `A_1`→2), `D_4`→1152,
  `D_5`→3840, and the named constant `E8_WEYL_GROUP_ORDER = 696729600`.
- **`discriminant.rs`** — the even-lattice discriminant form bridge: `DiscriminantForm
  { group, reps, gram_inv }` represents `A_L = L#/L` as `Z^n/GZ^n`;
  `quadratic_value_mod2`, `bilinear_value_mod1`, and `GaussSum::phase_mod8` compute the
  finite quadratic module; `verify_milgram` compares the Gauss-sum phase to the exact
  signature plus the genus oddity route. `Complex64`, `weil_t`, `weil_s`,
  `weil_s_prefactor_phase_mod8`, `weil_s_recovers_milgram_phase_mod8`, and
  `verify_weil_relations` implement the discriminant-form Weil representation.
  **Looks like a bug, isn't:** the standard Weil `S` prefactor is the conjugate of the
  positive Milgram phase stored by `GaussSum`; the verifier checks `S² = σ²·(γ↦−γ)`,
  `S⁴ = σ⁴·I`, and `(ST)³ = S²`, not the oversimplified `S⁴ = I`. The lattice ↔
  Clifford/Brauer-Wall mod-8 seam. Even-lattice only; odd type-I refinements stay a
  documented boundary.
- **`genus.rs`** — the **genus** = (signature, det, per-prime Conway–Sloane symbol).
  Engine: the p-adic Jordan decomposition (`jordan_blocks`, exact over `Rational`):
  odd `p` diagonalizes (valuation-ordered Gram–Schmidt); `p=2` peels 1-dim type-I lines
  and 2-dim even type-II planes by Schur complement. Per scale: `(dim, det mod 8, type,
  oddity = trace mod 8)` at `p=2`; odd `p` uses `(dim, det square class)`. `Genus::of` /
  `are_in_same_genus`. **Looks like a bug, isn't:** the comparison is **exact for odd
  `p`** (no sign-walking) and uses the full Conway–Sloane/Allcock fine-symbol reduction
  at `p=2` (normalize det residues, fuse compartment oddities, sign-walk left along
  trains adding `4` to crossed compartment oddities). The `Z⁸` (`1₀^{+8}`, type I) vs
  `E_8` (`1_{II}^{+8}`, type II), Sage canonical-symbol examples, and randomised `Uᵀ G U`
  isometry invariance pin the engine.
- **`mass_formula.rs`** — the **Minkowski–Siegel mass** of the even-unimodular genus,
  `mass(n) = |B_{n/2}|/n · ∏_{j<n/2} |B_{2j}|/(4j)` (Bernoulli by exact recurrence,
  `None` past the i128 ceiling). `mass(8) = 1/696729600 = 1/|W(E_8)|` — the formula
  *recovers* the `E_8` automorphism order the brute-force counter refuses; `n = 16, 24`
  match the published Niemeier values (i128 reaches exactly to 24). Plus the **Leech
  lattice** `leech()`: built from the Golay `[24,12,8]` code
  (`extended_golay_generator_rows`, `[I₁₂|A]`) → a `√8·Λ₂₄ ⊂ ℤ²⁴` spanning set → HNF
  basis `B` → `Gram = B·Bᵀ/8`. **Validated, not trusted:** rank-24 even unimodular with
  no roots *is* Leech (Niemeier), so the test checks `det=1`, even, `short_vectors(2)`
  empty (cheap; the full kissing 196560 is not enumerated). `|Aut(Λ₂₄)| = |Co₀|` is the
  factorized constant `LEECH_AUT_ORDER`.
- **`codes.rs`** — binary linear codes and Construction A: `BinaryCode` stores a checked
  row-reduced F₂ generator matrix; `dual`, `is_self_dual`, `is_self_orthogonal`,
  `is_doubly_even`, `minimum_distance`, `weight_enumerator`, `macwilliams_transform` are
  exact. `construction_a` uses the `1/sqrt(2)` scaling (HNF basis of `{x ∈ Z^n : x mod 2
  ∈ C}`, dot products /2); returns `None` when the scaled Gram is not integral. Shipped
  constructors: `hamming_code`, `extended_hamming_code`, `golay_code`,
  `type_ii_e8_sum_code`, `type_ii_len16_code`, `d16_plus`. **Looks like a bug, isn't:**
  bare Golay Construction A is even unimodular rank 24 **with roots**; it is not Leech.
- **`theta.rs` / `modular.rs`** — exact theta and modular-form bridge.
  `IntegralForm::theta_series(terms)` buckets short vectors by `Q/2`, `None` outside the
  positive-definite even-lattice boundary. `eisenstein_e4`, `eisenstein_e6`, `delta`,
  `mk_basis`, `as_modular_form` identify q-expansions exactly in `ℂ[E4,E6]`. Oracles pin
  `theta_E8 = E4`, `theta_{E8+E8} = theta_{D16+} = E4²`, Leech's rootless `q^1`
  coefficient in `E4³ - 720·Δ`, and the degenerate rank-16 Siegel–Weil consistency using
  `E8_WEYL_GROUP_ORDER` and `D16_PLUS_AUT_ORDER`.
