# AGENTS.md — `src/forms/`

The PILLAR of quadratic forms and their invariants. The organizing principle is
the **characteristic trichotomy**: the classification of a quadratic form
(equivalently, of the Clifford algebra it builds) is *one* theory split three ways
by `char F`. This axis cuts ACROSS the place table that organizes `scalar/`.

> Read root `OPEN.md` before touching `char2/`, `quadric_fit.rs`, `char0.rs`,
> `witt/`, or anything feeding the open play-semantics question.

`mod.rs` re-exports the legs + `classify` + diagonalize/equivalence + the `witt/`
invariant-group shelf (Witt group/ring + Brauer–Wall) + the `springer/`
valuation-graded (local↔global) decomposition + `local_global/` + `integral/` +
`field_invariants` + the symplectic/hermitian "form + involution" siblings, all
flat.

The directory groups the cross-cutting machinery into shelves that mirror how
`scalar/` is grouped by place: the trichotomy legs (`char0.rs`/`oddchar/`/`char2/`),
the invariant **groups** (`witt/`), the valuation-graded **decomposition**
(`springer/`), the **local↔global** layer (`local_global/`), and the **integral**
arithmetic view (`integral/`). Each multi-file cluster is its own subdir with a hub
`mod.rs` that re-exports flat, so public paths stay shallow (`forms::bw_class_real`,
`forms::springer_decompose_qp`, …) regardless of which shelf a symbol lives on.

Fixed-width mathematical payloads in this pillar are `u128`/`i128`: Arf bits,
Artin-Schreier classes, Dickson parities, Hilbert/Hasse signs, discriminant
residues, lattice entries, automorphism counts, and node budgets all follow that
contract. `usize` remains for dimensions and matrix indices.

## The façade

- **`classify.rs`** — the classifier FAÇADE: `ClassifyForm` + `WittClassify` +
  `IsometryClassify` + `WittDecompose` + `BrauerWallClassify`, keyed on the scalar
  so `metric.classify()` / `.witt_class()` / `.isometric_to()` / `.witt_decompose()`
  / `.bw_class()` pick the right leg **at compile time** (Surreal→CliffordType,
  Fp/Fpn→`FiniteFieldClass::{Odd, Char2}`, Nimber→ArfResult, finite-window
  Ordinal→ArfResult, …). `BrauerWallClassify` covers Surreal, Surcomplex, odd
  finite fields, nonsingular Nimber metrics, supported `Fpn<2,N>` metrics, and the
  documented finite ordinal windows. Rational &
  Surcomplex impl `ClassifyForm` but not `WittClassify` (their Witt data isn't a
  single `WittClassG` — honest, not a gap).
- **`diagonalize.rs`** — congruence diagonalization (char ≠ 2): `gram`,
  `diagonalize`, `as_diagonal`. Returns `None` in char 2 (nonsingular char-2 forms
  have an alternating polar form and are NOT diagonalizable — use the char-2
  symplectic Arf reduction). This is what lets char0/oddchar classify ARBITRARY
  (non-diagonal) metrics.
- **`equivalence.rs`** — isometry (per backend, via the complete invariant) + Witt
  decomposition (k·H ⊥ anisotropic kernel) over ℝ and F_q.

## The three legs

- **`char0.rs`** — the char-0 Clifford classifier: Cl(p,q) → matrix algebra over
  ℝ/ℂ/ℍ via the 8-fold table (real-closed surreal/rational) and the 2-fold table
  (surcomplex). `classify_real(p,q,r)` / `classify_complex(n,r)` are the bare-
  signature entry points (no metric needed); non-diagonal metrics are diagonalized
  first.
- **`oddchar/`** — odd-characteristic forms (re-exported flat): `field.rs`
  (`FiniteOddField` unifies Fp and Fpn square classes), `invariants.rs`
  (`classify_finite_odd`/`finite_odd_witt`/`discriminant`/`hasse` ≡ +1 over finite
  fields — ONE generic implementation keyed off the trait, Fp and Fpn share the
  path). dim + disc complete.
- **`char2/`** — characteristic-2 invariants (re-exported flat): `arf.rs` (the Arf
  invariant: `arf_f2` F₂ bitmask, `arf_nimber` for the represented nimber field,
  `arf_char2`/`arf_fpn_char2` for supported finite char-2 fields, and
  `arf_ordinal_finite` for the documented finite ordinal windows; all use
  symplectic reduction + trace and return `ArfResult { arf: u128, ... }`),
  `dickson.rs` (`dickson_matrix = rank(g−I) mod 2`, ker = SO;
  `dickson_of_versor` delegates to the generic versor grade parity), `field.rs`
  (`FiniteChar2Field` — the **additive** mirror of `FiniteOddField`: carries
  `artin_schreier_class = Tr_{F_q/F₂}` instead of `is_square_value`, since in char 2
  the multiplicative square class is trivial and the working datum is `F/℘(F) ≅ F₂`;
  impl for `Fp<2>`/`Fpn<2,N>`, NOT `Nimber` — same boundary as `FiniteOddField`).

The char0↔char2 classifier **symmetry** (the real 8-fold table mirrored by the
char-2 Arf/Brauer–Wall story) is one of the project's central threads.

## `witt/` — the invariant groups (Witt group, Witt ring, Brauer–Wall)

The three abelian groups the classifiers land in, grouped into one shelf (`mod.rs`
re-exports flat). This is the home of the **mod-8 spine**: `BW(ℝ) ≅ ℤ/8` is the
same periodicity as the char-0 8-fold table, Bott, and `E₈` in `integral/`.

- **`witt/class.rs`** — `WittClass`: the Witt group `W_q(F) ≅ ℤ/2` of a finite
  char-2 field/window, Arf-classified with `u128` bits. Plus `WittClassG`: the
  Char0/OddChar/Char2 trichotomy enum (odd-char is order-4) with checked group and
  ring operations; `try_mul` rejects Char2 because `W_q` is a module, not a ring.
- **`witt/ring.rs`** — the Witt RING: `tensor_form`, Pfister forms, fundamental
  ideal Iⁿ, the eₙ staircase (e0=dim, e1=disc, e2=Hasse). Stabilization per field
  (I²=0 over F_q; infinite ℝ tower via `e_real`). DON'T claim Arf=e2 (char-2
  indexing is Kato's, pinned).
- **`witt/brauer_wall.rs`** — the Brauer–Wall group BW(F): `bw_class_real` (Bott
  index (q−p) mod 8 ⇒ BW(ℝ)=ℤ/8), `bw_class_complex` (ℤ/2), `bw_class_oddchar`
  (order-4 ≅ W(F_q), DISCOVERED not asserted), `bw_class_nimber`, and façade
  dispatch for supported finite char-2 fields/windows (char-2 Arf/Witt class
  `ℤ/2`, nonsingular metrics only). Law = graded_tensor.

(The *numeric* field invariants the ring implies — level, u-invariant — live
separately in `field_invariants.rs`, not here; the name was disambiguated from the
old top-level `invariants.rs` for exactly that reason.)

## `springer/` — the discrete-valuation decomposition (a local–global symmetry)

One generic engine for the discretely-valued legs + the surreal odd-one-out +
the char-2 mirror, grouped into one shelf (`mod.rs` re-exports flat). The five
files were promoted from flat `springer_*.rs` at the `forms/` top level into this
subdir so the cluster reads as one local↔global story, matching the subdir
treatment of `oddchar/`/`char2/`/`local_global/`/`integral/`.

- **`springer/local.rs`** — the GENERIC engine `springer_decompose_local<K:
  ResidueField>` (+ `LocalResidueForm`/`LocalSpringerDecomp`/`parity_layer`), keyed
  off the `scalar::ResidueField` trait. ONE implementation; the residue field `k` is
  read through the trait (`residue_unit` = the angular component), the square-class
  via `is_square_finite`. Odd residue char only. (The siblings reach it via
  `super::local`, not the flat re-export, since they consume its types directly.)
- **`springer/padic.rs`** — the **mixed-characteristic** named entry points (thin
  wrappers returning `LocalSpringerDecomp`): `springer_decompose_qp` over `Q_p` (residue F_p) AND
  `springer_decompose_qq` over `Q_q` (residue F_q, the unramified extension — `F=1`
  recovers Q_p). Value group ℤ NOT 2-divisible ⇒ TWO residue layers survive
  (`parity_layer`) = W=W(k)². Adding Q_q makes this leg reach general F_q residues,
  matching the Laurent leg.
- **`springer/laurent.rs`** — the **equal-characteristic** named entry point (thin
  wrapper returning `LocalSpringerDecomp`): `springer_decompose_laurent` over `F_q((t))` (char p, residue
  F_q). Same two-layer story; residue char 2 REJECTED (the char-2 Witt boundary).
  Used by `local_global/function_field.rs` as an independent oracle.
- **`springer/surreal.rs`** — over the surreals (char 0, residue ℝ). The ONE that
  does NOT fit the generic engine: value group 2-divisible ⇒ W(No)=W(ℝ)=ℤ (second
  layer collapses), residue ℝ is a signature not a finite square-class. Keeps its
  own engine (owns the flat `ResidueForm`/`SpringerDecomp`/`springer_decompose`
  names) — that mismatch IS the symmetry, not a gap. So it stays out of
  `ResidueField`.
- **`springer/char2.rs`** — the **char-2 local Witt/Springer decomposition**
  (Aravire–Jacob `(φ₀, ψ, φ₁)`) + global isotropy over `F_q(t)`. Detailed below in
  the local–global section, because it is tightly coupled to
  `local_global/function_field_char2.rs` (it reuses that engine's `pub(crate)`
  helpers).

## Local–global (`local_global/`)

- **`local_global/global_field.rs`** — the `GlobalField` TRAIT: the local–global principle
  written ONCE over the two kinds of global field, `Rational` (ℚ, a number field)
  and `RationalFunction<S>` (`F_q(t)`, a function field). Five per-field primitives
  (`relevant_places`/`hilbert_symbol_at`/`is_local_square`/`is_global_square`/
  `is_isotropic_at_place`) + four DEFAULT theorem methods (`hasse_at_place`/
  `reciprocity_product`/`ramified_places`/`is_isotropic_global` = Hasse–Minkowski).
  The arithmetic primitives stay per-field (ℚ is i128 number theory with an
  archimedean place; `F_q(t)` is `F_q[t]` factorization with none — the missing
  real place IS the content), so `padic`/`adelic`/`function_field` keep their named
  public modules and functions, now thin wrappers over the trait. NOT a `Valued`
  abstraction (a global field carries all places at once, like `RationalFunction`/
  `Adele`). The mirror of what `ResidueField` did for the discrete Springer engine.
- **`local_global/padic.rs`** — the GENUINE Hilbert symbol over Q_p (odd-p + p=2 mod-8) — nontrivial
  unlike oddchar's +1 — + checked Hasse–Minkowski `try_is_isotropic_q` over ℚ.
  Oracle: Hilbert reciprocity `∏_v=+1`.
- **`local_global/adelic.rs`** — local–global rational helpers: `hilbert_product` over all places,
  rank≥3 adelic Hasse–Minkowski breakdown (`isotropy_over_adeles`/`AdelicIsotropy`),
  Brauer local invariant sums. Reuses `local_global/padic.rs`.
- **`local_global/function_field.rs`** — the **equal-characteristic mirror** of
  `local_global/padic.rs` + `local_global/adelic.rs` over the global function field
  `F_q(t)` (`scalar::RationalFunction`). Places
  `FFPlace{Infinite, Finite(π)}` (monic irreducibles + the degree place), the **tame**
  Hilbert symbol `hilbert_symbol_ff` (the odd-`p` `try_hilbert_symbol_qp` branch with the
  residue Legendre → `χ_κ`; **no `p=2` branch** since `q` is odd), reciprocity
  `hilbert_reciprocity_product_ff`, `is_isotropic_ff`/`is_isotropic_at_place`/
  `isotropy_over_ff_adeles` (Hasse–Minkowski, u-invariant 4 like `Q_p`, but **no
  archimedean place** ⇒ no definiteness condition), and `ramified_places_ff` (even
  count). Names carry `_ff` where `local_global/padic.rs` collides (e.g.
  `hasse_at_place_ff`). Exact (the product formula is `deg`-counting); odd residue
  char only — the `springer_laurent` boundary. Cross-checked against
  `springer_decompose_laurent`.
- **`local_global/function_field_char2.rs`** — the **equal-characteristic-2** mirror: the
  **asymmetric Artin–Schreier symbol** `[a,b)` over `F_{2^m}(t)` (`a` additive mod
  `℘`, `b` multiplicative), NOT the tame symbol. Local invariant = the **Schmid
  formula** `s_v(a,b) = Tr_{κ/F₂}(Res_v(a·dlog b))` (`as_symbol_at`), via a from-scratch
  **residue-of-differentials engine** (Hensel series `T(u)`, `P(T)=u`; the `∞` place by
  `u=1/t`). Reciprocity `∑_v s_v = 0` (`as_symbol_reciprocity_sum`, the gold oracle) +
  even ramification (`as_symbol_ramified_places`). Generic over `FiniteChar2Field`
  (so `F₂(t)`, `F₄(t)`, `F₈(t)` share one engine). Names carry `as_symbol_*` / `Char2Place`
  to avoid colliding with the odd `function_field` flat re-exports. The crate-private
  engine helpers (`strip_factor`/`inverse_mod`/`trace_kappa_to_f2`, and the factoriser
  `char2_monic_irreducible_factors` — a thin wrapper over the shared
  `poly_factor` finite-field factorizer, renamed off the odd-char
  `monic_irreducible_factors` so the flat `forms::*` glob stays unambiguous) are
  `pub(crate)` so `springer/char2.rs` reuses them.
- **`springer/char2.rs`** — the **char-2 local Witt/Springer decomposition**, the
  equal-char-2 mirror of `springer/local.rs` (but NOT the odd story at `p=2`: the wild
  `R_π` summand the `W=W(k)²` grading misses). `springer_decompose_local_char2(form,
  place)` gives the **Aravire–Jacob** `(φ₀, ψ, φ₁)` (`Char2LocalDecomp`): split each
  block coeff by Laurent-parity (`K=K²⊕πK²`), apply `[a,b]≅[1,a_ev·b]⊥⟨π⟩[1,a_odd·b]`,
  push each `[1,c]` to **Artin–Schreier normal form** (`asnf`: drop positive poles,
  clear even neg poles via `c_{n/2}+=√c_n`, keep the `κ`-constant Arf bit + odd neg
  poles `R_π`). Local isotropy `local_anisotropic_dim_char2`/`local_is_isotropic_char2`
  is invariant-driven (`ab∈℘(K_v)` for binary blocks, the AJ kernel for nonsingular
  parts, valuation parity for totally-singular tails, and the odd-dimensional
  Clifford invariant `Σ s_v(a_i b_i, c/a_i)` for one-class radical tails;
  `u(K_v)=4` ⇒ rank ≥ 5 isotropic). The form is `Char2QuadForm` (binary blocks + a
  totally-singular tail). **Read root AGENTS.md and OPEN.md** before touching: this is the corrected
  three-layer decomposition (the naive `W_q(k)²` was rightly avoided), pinned to
  source-derived oracles. **Global isotropy** `is_isotropic_global_char2(form) →
  Option<bool>` is Hasse–Minkowski over `F_q(t)` itself, on three ingredients past the
  symbol: `global_is_pe(f)` (`f ∈ ℘(F_q(t))`? — finite sweep of `f`'s poles + `∞`,
  settles rank 2: `[a,b]` iso ⟺ `ab ∈ ℘`), `ff_is_square(f)` (`f ∈ K²`? — all
  odd-degree coeffs of `num·den` vanish, settles the totally-singular part via
  `[K:K²]=2`), and a bad-place sweep over `relevant_places_char2(form)` for rank 3/4
  (good places isotropic by Chevalley–Warning). `u(F_q(t))=4` (`C₂`) ⇒ rank ≥ 5
  isotropic. **Looks like a bug, isn't:** rank 2 is NOT a finite bad-place scan — the constant-trace `℘`-obstruction
  (`[1,1]/F₂(t)`) lives at infinitely many odd-degree places, caught only by the global
  `℘` test.

## The "form + involution" siblings

- **`symplectic.rs`** — alternating forms: `SymplecticForm`, `hyperbolic`,
  `direct_sum`, `classify` (rank + radical_dim — the complete invariant, char-
  uniform). `classify_symplectic(gram)` convenience. The char-2 polar form of a
  nonsingular quadratic form lives here.
- **`hermitian.rs`** — Hermitian forms over Surcomplex (the involution `conj()` the
  symmetric leg never used): `HermitianForm` (conj-symmetric Gram), unitary
  congruence diagonalize → real diagonal, signature (Sylvester, the complete
  invariant = U(p,q)). `from_skew` handles the skew-Hermitian case via mult by i.

## Field invariants + the game bench

- **`field_invariants.rs`** — numeric FIELD invariants the Witt ring implies:
  level/Stufe s(F), pythagoras_number, u_invariant, is_sum_of_n_squares — computed
  over finite F_p (level≤2, u=2); ℝ/Q_p textbook constants documented. (Renamed from
  the old top-level `invariants.rs` to disambiguate from the `witt/` invariant-group
  shelf and `oddchar/invariants.rs`.)
- **`quadric_fit.rs`** — the "is this P-set a quadric?" research BENCH (split from the
  char2 classifier): `fit_f2_quadratic` (Boolean ANF/Mobius transform over the
  2^k membership table) + `QuadricFit` + `is_genuinely_quadratic`. The instrument the game
  probes / misère_quotient / octal_hunt / loopy_quadric feed P-positions into —
  distinct from the classifier.

## The trace-form bridge

- **`trace_form.rs`** — the seam from the `scalar::CyclicGaloisExtension` layer to
  the classifiers. `trace_twisted_form::<E>(k) -> Metric<E::Base>` builds the
  **Frobenius-twisted** trace form `Q_k(x) = Tr_{E/F}(x·σ^k(x))` (q on the diagonal,
  the alternating polar `Tr(eᵢσ^k eⱼ + eⱼσ^k eᵢ)` off it). NOT the naive `Tr(x²)`,
  whose polar form vanishes in char 2 (Frobenius is additive) — that's the trap the
  twist avoids. Instances: `Surcomplex` k=1 → the **norm form** `⟨2,2⟩`; unramified
  `Qq/Qp` via the Teichmuller-lifted residue basis; odd `Fpn` → a diagonalizable
  trace form. Two char-2 entry points to the **Gold form**
  `Tr(x^{1+2^a})`, classified → `ArfResult` (rank `= m − gcd(2a,m)`, Arf → the
  zero-count): `trace_form_arf::<E: …<Base=Fp<2>>>(k)` (the typed `Fpn<2,m>` path —
  build over `F_2`, lift `F_2 ↪ Nimber` via `Metric::map`), and `gold_form(m, a)`
  (the nim-native path over the subfield `F_{2^m} ⊂ Nimber`, m a power of two ≤ 128,
  reaching F_16/F_256/… that `Fpn` can't). This pulls the Python Gold thread
  (`experiments/trace_form_arf.py`, `gold_form_from_games.py`) into the typed core.
  The form has dim `[E:F]`, capped at `MAX_BASIS_DIM=128` (= the full nim-field's
  degree); the small power-of-two `m` keep the tests fast.
  The same `CyclicGaloisExtension` basis/generator data also feeds
  `clifford::frobenius::{galois_linear_map, frobenius_linear_map}`, giving the
  trace-form bridge a Clifford outermorphism oracle.

## Integral lattices (`integral/` — the arithmetic view / Arc 4)

The classifiers above work over a *field* (square classes / Witt / Arf). An
**integral lattice** is the complementary object: a free ℤ-module with an
integer Gram matrix. Its invariants are arithmetic (det, level, minimum, kissing
number, |Aut|), and the coarse classification is the **genus** (local
equivalence at every place), which reuses the `local_global/padic.rs` and
`local_global/adelic.rs` primitives. Staged M1→M4, all landed:
`integral/lattice.rs`, `integral/root_lattices.rs`, `integral/discriminant.rs`,
`integral/genus.rs`, `integral/mass_formula.rs` (+ the Leech lattice),
`integral/codes.rs`, `integral/theta.rs`, and `integral/modular.rs`.

- **`integral/lattice.rs`** (M1) — `IntegralForm { gram: Vec<Vec<i128>> }` (private Gram,
  built via `new` (square+symmetric-checked) / `diagonal`, never a struct literal).
  Convention: **norm** `Q(x) = xᵀGx` (a "norm-2 root" has `Q=2`; twice the value of
  `½Q`). `determinant` (fraction-free **Bareiss**, exact), `is_even`/`is_unimodular`,
  `is_positive_definite` (Sylvester leading-minors via Bareiss, exact),
  `signature` (exact rational diagonalization), `invariant_factors` (SNF →
  discriminant group `L#/L`), `level` (smallest `N` with `N·G⁻¹` even-integral, via
  the exact `Rational` inverse), `clifford_metric` (rational Clifford metric),
  `clifford_metric_f2` (even-lattice mod-2 char-2 metric), `direct_sum`. The
  positive-definite geometry: `short_vectors` (unimodular size-reduction first, then
  **Fincke–Pohst**: float LDLᵀ bounds the box, exact i128 norm filters the leaves,
  and vectors are mapped back to the original basis — float error can't add/drop a
  vector), `minimum`/`minimal_vectors`/`kissing_number`, and
  `automorphism_group_order` (closed-form diagonal/ADE/root-system fast paths first;
  otherwise backtracking over basis-vector images: `vᵢ` a lattice vector of norm `G[i][i]` with
  `⟨vᵢ,vⱼ⟩ = G[i][j]` — every complete assignment is an automorphism, so the count is
  exact). **Looks like a bug, isn't:** (a) the geometry methods return `None`
  for indefinite lattices on purpose (infinitely many vectors of each norm); (b) |Aut|
  is bounded by an explicit node budget (`AUTO_NODE_BUDGET`) and returns `None` past
  it (Leech is too big to brute-force; `E₈` is recognized by the ADE/root-system
  fast path) — an honest `None`, not a silent truncation
  (`automorphism_group_order_bounded` exposes the budget); (c) `level(⟨1⟩)=2`, not 1 —
  `ℤ` is odd, so the smallest `N` making `N·G⁻¹` even-diagonal is 2 (cf. `A_1=⟨2⟩→4`,
  `E₈→1`). Oracles: `A_2`/`A_3`/`D_4`/`E_8` det, kissing (6/12/24/240), |Aut|
  (12/48/1152), level (3/·/·/1), `Z^n` (|Aut| `2ⁿ·n!`).
- **`integral/root_lattices.rs`** (M2) — the ADE catalogue: `a_n` (Cartan matrix), `d_n`
  (`B·Bᵀ` from the geometric basis `{eᵢ−e_{i+1}}∪{e_{n-2}+e_{n-1}}`, sidestepping the
  fork-indexing), `e_6`/`e_7`/`e_8` (Dynkin edge lists). `coxeter_number = #roots/rank`
  (computed, not tabulated — an irreducible root system has `n·h` roots). `is_root_lattice`
  (min 2 + roots generate `L`, index read off the HNF pivots). `E_8` is the unique rank-8
  even unimodular lattice — the visible meeting point of the char-0 mod-8 Bott /
  `brauer_wall` BW(ℝ)=ℤ/8 story and the lattice world (root AGENTS.md). Det/kissing/Coxeter
  oracles protect every construction; |Aut| oracles include the small ones
  (`A_n`→`2(n+1)!`, `D_4`→1152, `D_5`→3840) and the named constant
  `E8_WEYL_GROUP_ORDER = 696729600`.
- **`integral/discriminant.rs`** — the even-lattice discriminant form bridge:
  `DiscriminantForm { group, reps, gram_inv }` represents `A_L = L#/L` as
  `Z^n/GZ^n`; `quadratic_value_mod2`, `bilinear_value_mod1`, and
  `GaussSum::phase_mod8` compute the finite quadratic module; and
  `verify_milgram` compares the Gauss-sum phase to the exact signature plus the
  genus oddity route. `Complex64`, `weil_t`, `weil_s`,
  `weil_s_prefactor_phase_mod8`, `weil_s_recovers_milgram_phase_mod8`, and
  `verify_weil_relations` implement the discriminant-form Weil representation.
  **Looks like a bug, isn't:** the standard Weil `S` prefactor is the conjugate of
  the positive Milgram phase stored by `GaussSum`; the verifier checks
  `S² = σ²·(γ↦−γ)`, `S⁴ = σ⁴·I`, and `(ST)³ = S²`, not the oversimplified
  `S⁴ = I`. This is the lattice ↔ Clifford/Brauer-Wall mod-8 seam. The
  implementation is even-lattice only; odd type-I refinements stay a documented
  boundary.
- **`integral/genus.rs`** (M3) — the **genus** = (signature, det, per-prime Conway–Sloane
  symbol). Engine: the **p-adic Jordan decomposition** (`jordan_blocks`, exact over
  `Rational`): odd `p` diagonalizes (valuation-ordered Gram–Schmidt, `e_i←e_i+e_j` to
  pull a diagonal pivot to the min valuation — `2` a unit); `p=2` peels 1-dim type-I
  lines and 2-dim even type-II planes by Schur complement. Per scale: `(dim,
  det mod 8, type, oddity = trace mod 8)` at `p=2`; odd `p` still uses `(dim,
  det square class)`. `Genus::of` / `are_in_same_genus`. **Looks like a bug,
  isn't:** (a) the comparison is **exact for odd `p`** (no sign-walking there) and
  uses the full Conway–Sloane/Allcock fine-symbol reduction at `p=2`: normalize
  determinant residues, fuse compartment oddities, then sign-walk left along trains
  while adding `4` to crossed compartment oddities (the giver/receiver bookkeeping);
  (b) signs/oddity are unused for odd `p`. The `Z⁸` (`1₀^{+8}`, type I) vs `E_8`
  (`1_{II}^{+8}`, type II), Sage canonical-symbol examples, and randomised
  `Uᵀ G U` isometry invariance pin the engine.
- **`integral/mass_formula.rs`** (M4) — the **Minkowski–Siegel mass** of the even-unimodular
  genus, `mass(n) = |B_{n/2}|/n · ∏_{j<n/2} |B_{2j}|/(4j)` (Bernoulli numbers generated
  by exact recurrence, checked cross-reduced rational mul → exact `(num, den)` or `None` past
  the i128 ceiling). `mass(8) = 1/696729600 = 1/|W(E_8)|` — the formula *recovers* the
  `E_8` automorphism order the brute-force counter refuses; `n = 16, 24` match the
  published Niemeier values (the i128 model reaches exactly to 24). Plus the **Leech
  lattice** `leech()`: built from the Golay `[24,12,8]` code
  (`extended_golay_generator_rows`, `[I₁₂|A]`) → a `√8·Λ₂₄ ⊂ ℤ²⁴` spanning set
  (`2·`Golay rows, `4(e₀+eᵢ)`, one odd
  `(−3,1²³)`) → HNF basis `B` → `Gram = B·Bᵀ/8`. **Validated, not trusted:** rank-24
  even unimodular with **no roots** *is* Leech (Niemeier), so the test checks `det=1`,
  even, `short_vectors(2)` empty (cheap — bound 2 < min 4; the full kissing 196560 is
  not enumerated). `|Aut(Λ₂₄)| = |Co₀|` is computed from the factorized constant `LEECH_AUT_ORDER` (far
  past brute force). Monster stays a thematic remark (moonshine, not a form computation).
- **`integral/codes.rs`** — binary linear codes and Construction A:
  `BinaryCode` stores a checked row-reduced F₂ generator matrix; `dual`,
  `is_self_dual`, `is_self_orthogonal`, `is_doubly_even`, `minimum_distance`,
  `weight_enumerator`, and `macwilliams_transform` are exact. `construction_a`
  uses the required `1/sqrt(2)` scaling by taking an HNF basis of
  `{x ∈ Z^n : x mod 2 ∈ C}` and dividing dot products by 2; it returns `None` when
  the scaled Gram is not integral. Shipped constructors: `hamming_code`,
  `extended_hamming_code`, `golay_code`, `type_ii_e8_sum_code`,
  `type_ii_len16_code`, and `d16_plus`. **Looks like a bug, isn't:** bare Golay
  Construction A is even unimodular rank 24 **with roots**; it is not Leech.
- **`integral/theta.rs` / `integral/modular.rs`** — exact theta and modular-form
  bridge. `IntegralForm::theta_series(terms)` buckets short vectors by `Q/2` and
  returns `None` outside the positive-definite even-lattice boundary. `eisenstein_e4`,
  `eisenstein_e6`, `delta`, `mk_basis`, and `as_modular_form` identify q-expansions
  exactly in `C[E4,E6]`. Oracles pin `theta_E8 = E4`,
  `theta_{E8+E8} = theta_{D16+} = E4^2`, Leech's rootless `q^1` coefficient in
  `E4^3 - 720 Delta`, and the degenerate rank-16 Siegel-Weil consistency using
  the recorded `E8_WEYL_GROUP_ORDER` and `D16_PLUS_AUT_ORDER` constants.

## Things that look like bugs but are not (forms layer)

- **`diagonalize`/`as_diagonal` return `None` in characteristic 2.** Not a bug: a
  nonsingular char-2 form has an alternating polar form and is not diagonalizable.
  The char-2 leg classifies via the symplectic Arf reduction (`char2/`) on the full
  (q, b) metric instead.
- **The odd-char Hasse invariant is ≡ +1** over a finite field — genuinely trivial
  there, unlike the p-adic Hilbert symbol in `local_global/padic.rs` (where Hasse does real work).
- **Rational & Surcomplex impl `ClassifyForm` but not `WittClassify`** — their Witt
  data isn't a single `WittClassG`. Honest, not a gap.
