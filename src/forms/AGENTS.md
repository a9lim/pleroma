# AGENTS.md — `src/forms/`

The PILLAR of quadratic forms and their invariants. The organizing principle is
the **characteristic trichotomy**: the classification of a quadratic form
(equivalently, of the Clifford algebra it builds) is *one* theory split three ways
by `char F`. This axis cuts ACROSS the place table that organizes `scalar/`.

> Read root `OPEN.md` before touching `char2/`, `quadric_fit.rs`, `char0.rs`,
> `witt/`, or anything feeding the open play-semantics question.

`mod.rs` re-exports the legs + `classify` + diagonalize/equivalence + the `witt/`
invariant-group shelf + the `springer/` valuation-graded decomposition +
`local_global/` + `integral/` + `field_invariants` + the
symplectic/hermitian "form + involution" siblings, all flat (`poly_factor` is
`pub(crate)` machinery, not re-exported). The cross-cutting
machinery is grouped into shelves mirroring how `scalar/` groups by place: the
trichotomy legs (`char0.rs`/`oddchar/`/`char2/`), the invariant **groups**
(`witt/`), the valuation-graded **decomposition** (`springer/`), the **local↔global**
layer (`local_global/`), and the **integral** arithmetic view (`integral/`). Each
multi-file cluster is a subdir with a hub `mod.rs` re-exporting flat, so public paths
stay shallow (`forms::bw_class_real`, `forms::springer_decompose_qp`, …).

`integral/` has its own [`AGENTS.md`](integral/AGENTS.md) (the lattice/genus/
mass/code/theta/Weil chain).

Fixed-width payloads here are `u128`/`i128`: Arf bits, Artin-Schreier classes,
Dickson parities, Hilbert/Hasse signs, discriminant residues, lattice entries,
automorphism counts, node budgets. `usize` is for dimensions and matrix indices.

## The façade

- **`classify.rs`** — the classifier FAÇADE: `ClassifyForm` + `WittClassify` +
  `IsometryClassify` + `WittDecompose` + `BrauerWallClassify`, keyed on the scalar so
  `metric.classify()` / `.witt_class()` / `.isometric_to()` / `.witt_decompose()` /
  `.bw_class()` pick the right leg **at compile time**:
  - `Surreal` → `CliffordType`
  - `Fp<P>` (odd primes only) → `OddCharType`. `Fp<2>` is OUTSIDE the façade:
    `classify_finite_odd` returns `None` for `P=2` (the char-2 Arf path requires
    the `(q,b)` metric, not the char-0 diagonalizer). Use `Fpn<2,N>` for char-2
    extension fields.
  - `Fpn<P,N>` → `FiniteFieldClass::{Odd, Char2}` (dispatched on `P==2` at runtime)
  - `Nimber` → `ArfResult`
  - finite-window `Ordinal` → `ArfResult`

  `WittDecompose` returns the leg-specific decomp record
  (`RealWittDecomp` / `OddWittDecomp` / `Char2WittDecomp`, with `Fpn` wrapping the last
  two in `FiniteFieldWittDecomp { Odd, Char2 }`). `BrauerWallClassify`
  covers Surreal, Surcomplex, odd finite fields, nonsingular Nimber metrics,
  supported `Fpn<2,N>` metrics, and the documented finite ordinal windows. Rational &
  Surcomplex impl `ClassifyForm` but not `WittClassify` (their Witt data isn't a
  single `WittClassG` — honest, not a gap).

  **Cross-leg asymmetry in `BrauerWallClassify` for singular metrics.** The
  char-2 legs (`bw_class_nimber`, `BrauerWallClassify for Fpn<2,N>`) return `None`
  for singular (non-nonsingular-polar) metrics.  The char-0 legs (`bw_class_real`,
  `bw_class_complex`) and the odd-char leg (`bw_class_finite_odd`) silently project
  the radical away and return the Brauer-Wall class of the nondegenerate core
  `Cl(Q/rad)`. The rustdocs on those three functions state this projection
  explicitly.
- **`diagonalize.rs`** — congruence diagonalization (char ≠ 2): `gram`,
  `diagonalize`, `as_diagonal`. Returns `None` in char 2 (nonsingular char-2 forms
  have an alternating polar form and are NOT diagonalizable — use the char-2
  symplectic Arf reduction). This is what lets char0/oddchar classify ARBITRARY
  (non-diagonal) metrics.
- **`equivalence.rs`** — isometry per backend (via the complete invariant:
  `isometric_finite_odd`, `isometric_finite_char2`/`isometric_fpn_char2`, real/
  rational/surcomplex/nimber) + Witt decomposition (`witt_decompose_real` →
  `RealWittDecomp`, `witt_decompose_finite_odd` → `OddWittDecomp`: k·H ⊥ anisotropic
  kernel).

## The three legs

- **`char0.rs`** — the char-0 Clifford classifier: Cl(p,q) → matrix algebra over
  ℝ/ℂ/ℍ via the 8-fold table (real-closed surreal/rational) and the 2-fold table
  (surcomplex). `classify_real(p,q,r)` / `classify_complex(n,r)` are the
  bare-signature entry points (no metric needed); non-diagonal metrics are
  diagonalized first.
- **`oddchar/`** — odd-characteristic forms (re-exported flat): `field.rs`
  (`FiniteOddField` unifies Fp and Fpn square classes), `invariants.rs`
  (`classify_finite_odd`/`finite_odd_witt`/`discriminant_finite_odd`/
  `hasse_invariant_finite_odd` ≡ +1 over finite fields — ONE generic implementation
  keyed off the trait; Fp and Fpn share the path). dim + disc complete.
- **`char2/`** — characteristic-2 invariants (re-exported flat): `arf.rs` (the Arf
  invariant: `arf_f2` F₂ bitmask, `arf_nimber` for the represented nimber field,
  `arf_char2`/`arf_fpn_char2` for supported finite char-2 fields, `arf_ordinal_finite`
  for the documented finite ordinal windows; all use symplectic reduction + trace and
  return `ArfResult { arf: u128, ... }`), `dickson.rs` (`dickson_matrix = rank(g−I)
  mod 2`, ker = SO; `dickson_of_versor` validates the input is a versor then delegates
  to the generic versor grade parity), `field.rs` (`FiniteChar2Field` — the
  **additive** mirror of `FiniteOddField`: carries `artin_schreier_class =
  Tr_{F_q/F₂}` instead of `is_square_value`, since in char 2 the multiplicative square
  class is trivial and the working datum is `F/℘(F) ≅ F₂`; impl for `Fp<2>`/`Fpn<2,N>`,
  NOT `Nimber` — same boundary as `FiniteOddField`).

The char0↔char2 classifier **symmetry** (the real 8-fold table mirrored by the
char-2 Arf/Brauer–Wall story) is one of the project's central threads.

## `witt/` — the invariant groups (Witt group, Witt ring, Brauer–Wall)

The three abelian groups the classifiers land in, one shelf (`mod.rs` re-exports
flat). Home of the **mod-8 spine**: `BW(ℝ) ≅ ℤ/8` is the same periodicity as the
char-0 8-fold table, Bott, and `E₈` in `integral/`.

- **`witt/class.rs`** — `WittClass`: the Witt group `W_q(F) ≅ ℤ/2` of the nimber
  field, Arf-classified with `u128` bits (`try_from_metric` takes `Metric<Nimber>`;
  ordinal-window metrics go through the `classify` façade's `WittClassG::Char2`). Plus
  `WittClassG`: the Char0/OddChar/Char2 trichotomy enum (odd-char is order-4) with
  checked group and ring operations; `try_mul` rejects Char2 because `W_q` is a
  module, not a ring.
- **`witt/ring.rs`** — the Witt RING: `tensor_form`, Pfister forms, fundamental ideal
  Iⁿ, the eₙ staircase (e0=dim, e1=disc, e2=Hasse). Stabilization per field (I²=0 over
  F_q; infinite ℝ tower via `e_real`). DON'T claim Arf=e2 (char-2 indexing is Kato's,
  pinned).
- **`witt/brauer_wall.rs`** — the Brauer–Wall group BW(F): `bw_class_real` (Bott index
  (q−p) mod 8 ⇒ BW(ℝ)=ℤ/8), `bw_class_complex` (ℤ/2), `bw_class_finite_odd` (order-4 ≅
  W(F_q)), `bw_class_nimber`, and façade dispatch for supported finite char-2
  fields/windows (char-2 Arf/Witt class `ℤ/2`, nonsingular metrics only). Law =
  graded_tensor.

(The *numeric* field invariants the ring implies — level, u-invariant — live in
`field_invariants.rs`, not here.)

## `springer/` — the discrete-valuation decomposition (a local–global symmetry)

One generic engine for the discretely-valued legs + the surreal odd-one-out + the
char-2 mirror, one shelf (`mod.rs` re-exports flat).

- **`springer/local.rs`** — the GENERIC engine `springer_decompose_local<K:
  ResidueField>` (+ `LocalResidueForm`/`LocalSpringerDecomp`/`parity_layer`), keyed
  off `scalar::ResidueField`. ONE implementation; the residue field `k` is read
  through the trait (`residue_unit` = the angular component), the square-class via
  `is_square_finite`. Odd residue char only.
- **`springer/padic.rs`** — the **mixed-characteristic** named entry points (thin
  wrappers returning `LocalSpringerDecomp`): `springer_decompose_qp` over `Q_p`
  (residue F_p) AND `springer_decompose_qq` over `Q_q` (residue F_q; `F=1` recovers
  Q_p). Value group ℤ NOT 2-divisible ⇒ TWO residue layers survive (`parity_layer`) =
  W=W(k)².
- **`springer/laurent.rs`** — the **equal-characteristic** named entry point:
  `springer_decompose_laurent` over `F_q((t))` (char p, residue F_q). Same two-layer
  story; residue char 2 REJECTED (the char-2 Witt boundary). Used by
  `local_global/function_field.rs` as an independent oracle.
- **`springer/surreal.rs`** — over the surreals (char 0, residue ℝ). The ONE that
  does NOT fit the generic engine: value group 2-divisible ⇒ W(No)=W(ℝ)=ℤ (second
  layer collapses), residue ℝ is a signature not a finite square-class. Keeps its own
  engine (the flat `ResidueForm`/`SpringerDecomp`/`springer_decompose` names) — that
  mismatch IS the symmetry, not a gap. So it stays out of `ResidueField`.
- **`springer/char2.rs`** — the **char-2 local Witt/Springer decomposition**
  (Aravire–Jacob `(φ₀, ψ, φ₁)`) + global isotropy over `F_q(t)`. Tightly coupled to
  `local_global/function_field_char2.rs` (it reuses that engine's `pub(crate)`
  helpers); detailed in the local–global section.

## Local–global (`local_global/`)

- **`local_global/global_field.rs`** — the `GlobalField` TRAIT: the local–global
  principle written ONCE over the two kinds of global field, `Rational` (ℚ) and
  `RationalFunction<S>` (`F_q(t)`). Five checked per-field primitives
  (`try_relevant_places`/`try_hilbert_symbol_at`/`try_is_local_square`/
  `try_is_global_square`/`try_is_isotropic_at_place`) + four checked DEFAULT theorem
  methods (`try_hasse_at_place`/`try_reciprocity_product`/`try_ramified_places`/
  `try_is_isotropic_global` = Hasse–Minkowski). The arithmetic primitives stay
  per-field (ℚ is i128 number theory with an archimedean place; `F_q(t)` is `F_q[t]`
  factorization with none — the missing real place IS the content), so
  `padic`/`adelic`/`function_field` keep their named public modules as thin wrappers
  over the trait. NOT a `Valued` abstraction (a global field carries all places at
  once, like `Adele`). The mirror of what `ResidueField` does for the discrete
  Springer engine.
- **`local_global/padic.rs`** — the GENUINE Hilbert symbol over Q_p (odd-p + p=2
  mod-8) — nontrivial unlike oddchar's +1 — + checked Hasse–Minkowski
  `try_is_isotropic_q` over ℚ. Oracle: Hilbert reciprocity `∏_v=+1`.
- **`local_global/adelic.rs`** — local–global rational helpers: `hilbert_product` over
  all places, rank≥3 adelic Hasse–Minkowski breakdown
  (`isotropy_over_adeles`/`AdelicIsotropy`), Brauer local invariant sums. Reuses
  `local_global/padic.rs`.
- **`local_global/function_field.rs`** — the **equal-characteristic mirror** of
  `padic.rs` + `adelic.rs` over `F_q(t)`. Places `FFPlace{Infinite, Finite(π)}`
  (monic irreducibles + the degree place), the **tame** Hilbert symbol
  `try_hilbert_symbol_ff` (the odd-`p` branch with the residue Legendre → `χ_κ`; no
  `p=2` branch since `q` is odd), reciprocity `try_hilbert_reciprocity_product_ff`,
  `try_is_isotropic_ff`/`try_is_isotropic_at_place_ff`/`try_isotropy_over_ff_adeles`
  (Hasse–Minkowski, u-invariant 4 like `Q_p`, but **no archimedean place** ⇒ no
  definiteness condition), `try_ramified_places_ff` (even count). Names carry `_ff`
  where `padic.rs` collides. Exact (the product formula is `deg`-counting); odd
  residue char only. Cross-checked against `springer_decompose_laurent`.
- **`local_global/function_field_char2.rs`** — the **equal-characteristic-2** mirror:
  the **asymmetric Artin–Schreier symbol** `[a,b)` over `F_{2^m}(t)` (`a` additive mod
  `℘`, `b` multiplicative), NOT the tame symbol. Local invariant = the **Schmid
  formula** `s_v(a,b) = Tr_{κ/F₂}(Res_v(a·dlog b))` (`as_symbol_at`), via a
  from-scratch residue-of-differentials engine (Hensel series `T(u)`, `P(T)=u`; the
  `∞` place by `u=1/t`). Reciprocity `∑_v s_v = 0` (`as_symbol_reciprocity_sum`, the
  gold oracle) + even ramification (`as_symbol_ramified_places`). Generic over
  `FiniteChar2Field` (so `F₂(t)`, `F₄(t)`, `F₈(t)` share one engine). Names carry
  `as_symbol_*` / `Char2Place`. The crate-private helpers (`strip_factor`/
  `inverse_mod`/`trace_kappa_to_f2`, and `char2_monic_irreducible_factors` — a thin
  wrapper over the shared `poly_factor` finite-field factorizer) are `pub(crate)` so
  `springer/char2.rs` reuses them.
- **`springer/char2.rs`** (detail) — the equal-char-2 mirror of `springer/local.rs`
  (but NOT the odd story at `p=2`: the wild `R_π` summand the `W=W(k)²` grading
  misses). `springer_decompose_local_char2(form, place)` gives the **Aravire–Jacob**
  `(φ₀, ψ, φ₁)` (`Char2LocalDecomp`): split each block coeff by Laurent-parity
  (`K=K²⊕πK²`), apply `[a,b]≅[1,a_ev·b]⊥⟨π⟩[1,a_odd·b]`, push each `[1,c]` to
  **Artin–Schreier normal form** (`asnf`: drop positive poles, clear even neg poles
  via `c_{n/2}+=√c_n`, keep the `κ`-constant Arf bit + odd neg poles `R_π`). Local
  isotropy `local_anisotropic_dim_char2`/`local_is_isotropic_char2` is
  invariant-driven (`ab∈℘(K_v)` for binary blocks, the AJ kernel for nonsingular
  parts, valuation parity for totally-singular tails, the odd-dimensional Clifford
  invariant `Σ s_v(a_i b_i, c/a_i)` for one-class radical tails; `u(K_v)=4` ⇒ rank ≥ 5
  isotropic). The form is `Char2QuadForm` (binary blocks + a totally-singular tail).
  **Global isotropy** `is_isotropic_global_char2(form) → Option<bool>` is
  Hasse–Minkowski over `F_q(t)` itself, on three ingredients past the symbol:
  `global_is_pe(f)` (`f ∈ ℘(F_q(t))`? — finite sweep of `f`'s poles + `∞`, settles
  rank 2: `[a,b]` iso ⟺ `ab ∈ ℘`), `ff_is_square(f)` (`f ∈ K²`? — all odd-degree
  coeffs of `num·den` vanish, settles the totally-singular part via `[K:K²]=2`), and a
  bad-place sweep over `relevant_places_char2(form)` for rank 3/4 (good places
  isotropic by Chevalley–Warning). `u(F_q(t))=4` (`C₂`) ⇒ rank ≥ 5 isotropic.
  **Looks like a bug, isn't:** rank 2 is NOT a finite bad-place scan — the
  constant-trace `℘`-obstruction (`[1,1]/F₂(t)`) lives at infinitely many odd-degree
  places, caught only by the global `℘` test.

## The "form + involution" siblings

- **`symplectic.rs`** — alternating forms: `SymplecticForm`, `hyperbolic`,
  `direct_sum`, `classify` (rank + radical_dim — the complete invariant,
  char-uniform). `classify_symplectic(gram)` convenience. The char-2 polar form of a
  nonsingular quadratic form lives here.
- **`hermitian.rs`** — Hermitian forms over Surcomplex (the involution `conj()` the
  symmetric leg never used): `HermitianForm` (conj-symmetric Gram), unitary congruence
  diagonalize → real diagonal, signature (Sylvester, the complete invariant = U(p,q)).
  `from_skew` handles the skew-Hermitian case via mult by i.

## Field invariants, the trace-form bridge, and the game bench

- **`field_invariants.rs`** — numeric FIELD invariants the Witt ring implies:
  level/Stufe s(F), pythagoras_number, u_invariant, is_sum_of_n_squares — computed
  over finite F_p (level≤2, u=2); ℝ/Q_p textbook constants documented.
- **`poly_factor.rs`** — the shared finite-field polynomial factorizer
  (`pub(crate) monic_irreducible_factor_support`), the place-finding primitive behind
  both function-field symbol layers; the public wrappers `monic_irreducible_factors`
  (`local_global/function_field.rs`, odd) and `char2_monic_irreducible_factors`
  (`local_global/function_field_char2.rs`, char-2) call into it.
- **`trace_form.rs`** — the seam from `scalar::CyclicGaloisExtension` to the
  classifiers. `trace_twisted_form::<E>(k) -> Metric<E::Base>` builds the
  **Frobenius-twisted** trace form `Q_k(x) = Tr_{E/F}(x·σ^k(x))` (q on the diagonal,
  the alternating polar `Tr(eᵢσ^k eⱼ + eⱼσ^k eᵢ)` off it). NOT the naive `Tr(x²)`,
  whose polar form vanishes in char 2 (Frobenius is additive) — the trap the twist
  avoids. Instances: `Surcomplex` k=1 → the **norm form** `⟨2,2⟩`; unramified `Qq/Qp`
  via the Teichmuller-lifted residue basis; odd `Fpn` → a diagonalizable trace form.
  Two char-2 entry points to the **Gold form** `Tr(x^{1+2^a})`, classified →
  `ArfResult` (rank `= m − gcd(2a,m)`, Arf → the zero-count): `trace_form_arf::<E:
  …<Base=Fp<2>>>(k)` (the typed `Fpn<2,m>` path — build over `F_2`, lift `F_2 ↪
  Nimber` via `Metric::map`), and `gold_form(m, a)` (the nim-native path over the
  subfield `F_{2^m} ⊂ Nimber`, m a power of two ≤ 128, reaching F_16/F_256/… that
  `Fpn` can't). The form has dim `[E:F]`, capped at `MAX_BASIS_DIM=128`. The same
  `CyclicGaloisExtension` basis/generator data also feeds
  `clifford::frobenius::{galois_linear_map, frobenius_linear_map}`, giving the bridge
  a Clifford outermorphism oracle.
- **`quadric_fit.rs`** — the "is this P-set a quadric?" research BENCH (split from the
  char2 classifier): `fit_f2_quadratic` (Boolean ANF/Möbius transform over the 2^k
  membership table) + `QuadricFit` + `is_genuinely_quadratic`. The instrument the game
  probes (misère_quotient / octal_hunt / loopy_quadric) feed P-positions into —
  distinct from the classifier.

## Things that look like bugs but are not (forms layer)

- **`diagonalize`/`as_diagonal` return `None` in characteristic 2.** Not a bug: a
  nonsingular char-2 form has an alternating polar form and is not diagonalizable. The
  char-2 leg classifies via the symplectic Arf reduction (`char2/`) on the full (q, b)
  metric instead.
- **The odd-char Hasse invariant is ≡ +1** over a finite field — genuinely trivial
  there, unlike the p-adic Hilbert symbol in `local_global/padic.rs` (where Hasse does
  real work).
- **Rational & Surcomplex impl `ClassifyForm` but not `WittClassify`** — their Witt
  data isn't a single `WittClassG`. Honest, not a gap.
