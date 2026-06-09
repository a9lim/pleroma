# AGENTS.md — `src/py/`

The PyO3 bindings, behind the `python` feature (`pyo3` is an optional dep;
`extension-module` only enabled here). This is the ONLY place `use pyo3` may
appear — keeping it out of the core is what stops `cargo test` linking libpython.

Python exposes Rust `u128`/`i128` payloads as normal Python integers. Keep that
width contract at the Rust boundary for form invariants, lattice entries, field
orders, game values, and budgets. `usize` is still required for dimensions, indices,
and PyO3 slots such as `__hash__`.

`mod.rs` is the `#[pymodule]`; it chains each submodule's `pub(crate) register()`.
Split per pillar:

- **`scalars.rs`** — the scalar pyclasses (`Nimber`, `NimberPoly`,
  `NimberRationalFunction`, fixed prime-field function-field rows
  `Fp2Poly`/`Fp2RationalFunction` through `Fp13Poly`/`Fp13RationalFunction`,
  fixed finite fields `Fp2`/`Fp3`/`Fp5`/`Fp7`/`Fp11`/
  `Fp13` and `F4`/`F8`/`F16`/`F9`/`F25`/`F27`, fixed p-adic slices
  `Zp2_4`/`Zp3_4`/`Zp5_4`/`Zp7_4`/`Zp11_4`/`Zp13_4` and
  `Qp2_4`/`Qp3_4`/`Qp5_4`/`Qp7_4`/`Qp11_4`/`Qp13_4`, fixed unramified slices
  `WittVec2_4_2`/`WittVec2_4_3`/`WittVec2_4_4`/`WittVec3_4_2`/
  `WittVec5_4_2`/`WittVec3_4_3` and
  `Qq2_4_2`/`Qq2_4_3`/`Qq2_4_4`/`Qq3_4_2`/`Qq5_4_2`/`Qq3_4_3`,
  fixed Laurent slices `LaurentRational_6`, `LaurentFp3_6`/`LaurentFp5_6`/
  `LaurentFp7_6`/`LaurentFp11_6`/`LaurentFp13_6`, `LaurentF9_6`/
  `LaurentF25_6`/`LaurentF27_6`, and fixed ramified quadratic/cubic p-adic slices
  `RamifiedQp2_4_E2`/`RamifiedQp3_4_E2`/`RamifiedQp5_4_E2`/
  `RamifiedQp7_4_E2`/`RamifiedQp11_4_E2`/`RamifiedQp13_4_E2` plus
  `RamifiedQp2_4_E3`/`RamifiedQp3_4_E3`/`RamifiedQp5_4_E3`/
  `RamifiedQp7_4_E3`/`RamifiedQp11_4_E3`/`RamifiedQp13_4_E3`,
  fixed Gauss slices `GaussQp2_4`/`GaussQp3_4`/`GaussQp5_4`/
  `GaussQp7_4`/`GaussQp11_4`/`GaussQp13_4`,
  `Rational`, `Surreal`, `Surcomplex`, `Integer`, `Omnific`, `Ordinal`,
  `SignExpansion`,
  `LocalQp`, `Adele`, `MaxPlusTropical`, `MinPlusTropical`) + constructors +
  nim-field free fns.
  Threads the `parse_*`/`wrap_*` hooks the `backend!` macro consumes. Surreal also
  exposes the simplicity bridge; Nimber exposes `pow`/`__pow__`/`frobenius`/`sqrt`
  alongside the `nim_*` Galois free fns. `LocalQp` is the runtime-prime p-adic cell
  (validated before calling the panic-guarded Rust constructors); `Adele` exposes
  the diagonal embedding, finite-place corrections, local components, idele
  predicates, the typed-only `AdelePlace` value used by Rust-name
  `absolute_value_at`, the product-formula surface, and the module helper
  `adele_prec`; `Rational` is the exact-ℚ scalar, including the Rust-name
  `try_new` constructor, `numer`/`denom` accessors, and Python value ordering;
  rational-valued surreals are constructed explicitly with
  `Surreal.from_rational`; `Surreal` exposes raw CNF terms, Python value ordering, sign/rational/dyadic projections,
  monomial construction, both finite and transfinite sign-expansion round-trips,
  the `SignExpansion` value type (`from_runs`, `from_finite`, `runs`, `length`,
  `as_finite`),
  plus the ordinal embedding/projection and lazy analytic helpers
  (`sqrt_to_terms`, `nth_root_to_terms`); `Ordinal` exposes CNF terms, Python ordinal ordering, the
  `<ω³` coefficient bridge, staged `nim_mul` and `checked_inv`, and feeds the
  `OrdinalAlgebra` backend; the
  finite-field classes expose the shipped `Fp`/`Fpn` scalar monomorphs, exact
  `Fpn::from_coeffs`/`into_coeffs` coefficient round-trips, exact
  square roots, Rust validator names (`assert_prime_modulus`, `is_supported_field`,
  `assert_supported_field`), reduction-polynomial metadata, group-order data, Frobenius
  iteration, Rust-name `FiniteField` methods (`ext_degree`, `min_poly_monic`,
  `relative_trace_over`, `relative_norm_over`, `multiplicative_order`), and their
  Galois toolkit; the fixed p-adic classes expose the
  `Scalar` surface plus exact `Qp::from_i128`, valuation/unit/square-root predicates
  and checked `is_square`/`sqrt` methods for `Z/p^4` and
  `Q_p` at `p ∈ {2,3,5,7,11,13}`, Rust validator names
  (`assert_supported_ring`, `assert_supported_field`), the local-field package
  (`uniformizer`, `residue`, `residue_unit`/`angular_component`, `teichmuller`,
  and, for `Qp`/`Qq`, `is_integral`/`to_integer`), and for the fixed unramified
  `W_4(F_q)` / `Q_q`, Laurent, and ramified quadratic/cubic cells listed above.
  The fixed `Q_q` classes also expose the Rust `FieldExtension`/`CyclicGaloisExtension`
  surface (`extension_degree`, `embed`, `trace`, `norm`, `basis`, `sigma`,
  `sigma_power`), with trace/norm returned in the matching `Qp*_4` Python base;
  arbitrary const-generic `Qq<P,N,F>` families remain Rust-only until there is an
  explicit runtime dispatch design. The fixed `GaussQp*_4` classes expose the
  valued/integrality surface plus residue/Teichmuller through the fixed
  `Fp*RationalFunction` residue classes;
  `NimberPoly`/`NimberRationalFunction` and the fixed
  `Fp*Poly`/`Fp*RationalFunction` pairs expose the exact
  `F_{2^128}[t]` / `F_{2^128}(t)` and prime-field `F_p[t]` / `F_p(t)` rows,
  with arithmetic through Python operators; the Rust-name `min_coeff_valuation` is exposed as a
  homogeneous typed coefficient-list dispatcher for the bound valued scalar
  classes (`Qp`/`Qq`/`Laurent`/`Ramified`/`Gauss`) because Python has no generic
  runtime `Poly<S: Valued>` type parameter; module-level nim arithmetic free
  functions expose the Rust
  names `nim_add`, `nim_mul`, `nim_pow`, `nim_square`, `nim_frobenius_iter`, and
  `nim_inv` alongside the existing trace/Galois toolkit; the tropical classes expose the dual
  semiring endpoints that thermography names, with semiring operations through `+`/`*`.
  Bound scalar classes expose the shared runtime `Scalar` surface
  (`zero`/`one`/`characteristic` where applicable, `is_zero`, partial
  inverses/division, and owned arithmetic operators), not only the operations
  each demo happened to need. Finite `Fpn`
  classes return the typed `ReductionPolynomialKind` metadata object from the
  Rust-name `reduction_polynomial_kind`.
- **`engine.rs`** — the `backend!` macro → `<World>Algebra` + `<World>MV` pairs
  (Nimber, the fixed finite fields, the fixed p-adic `Zp*_4`/`Qp*_4` slices,
  the fixed unramified `WittVec*_4_*`/`Qq*_4_*` slices, the fixed
  `Laurent*_6` slices, the fixed `RamifiedQp*_4_E2`/`RamifiedQp*_4_E3` slices,
  the fixed `GaussQp*_4` slices, NimberPoly/NimberRationalFunction,
  `Fp*Poly`/`Fp*RationalFunction`, Rational,
  Adele, Surreal/Surcomplex, Integer/Omnific, Ordinal) + conformal GA over every bound
  characteristic-zero scalar world with a matching MV carrier (`SurrealCga`,
  `RationalCga`, `AdeleCga`, `SurcomplexCga`,
  fixed `Qp*_4Cga`, fixed `Qq*_4_*Cga`, `LaurentRational6Cga`,
  `RamifiedQp*_4_E2Cga`/`RamifiedQp*_4_E3Cga`, and `GaussQp*_4Cga`). CGA is
  intentionally absent for positive-characteristic worlds because Rust `Cga::new`
  rejects them, for `Integer`/`Omnific` because `1/2` is not invertible, for
  `LaurentFp*`/`LaurentFq*` because their base characteristic is positive, and for
  runtime `LocalQp` because the Rust core deliberately does not implement
  `Scalar` for a value whose world `(p,k)` is only known at construction. MV methods
  cover the full Arc-C suite (clifford_conjugate,
  scalar_product, commutator, anticommutator, undual, meet, is_blade,
  blade_subspace, factor_blade, cayley, spinor_norm, versor_grade_parity,
  classify_versor, plus the named `VersorClass` record) and expose raw
  `(blade_mask, coeff)` terms; algebra methods add
  trace/char_poly/determinant/exterior_power_trace/outermorphism/apply_outermorphism/
  inverse_outermorphism, typed `<World>LinearMap` pyclasses for Rust
  `LinearMap<S>` (`from_columns`, `identity`, `n`, `cols`, `image`, `compose`)
  with algebra spectral/outhermorphism methods accepting only the matching typed
  map object, named `SpinorRep` records via `spinor_rep`, and lazy spinor action.
  The module-level `bits`/`grade` helpers expose the same mask utilities used by
  the Rust basis layer.
  Each algebra class also has Rust-name metric constructors/helpers (`general`,
  `grassmann`, `q`, `b_terms`, `a_terms`, `map`, `q_val`, `has_upper`, `is_orthogonal`),
  the named `pga(n)` constructor for `Cl(n,0,1)`, the
  symmetric-form `gram`/`diagonalize`/`as_diagonal` façade, and tensor-square plus
  explicit graded-tensor embedding helpers.
  They also expose the generic Witt-ring representative operations
  (`tensor_form`, `pfister1`, `pfister`, `in_fundamental_ideal`) for the bound scalar
  world, plus fixed-dispatch Frobenius/Galois matrix constructors
  Rust-name fixed-dispatch typed-map constructors
  (`frobenius_linear_map`, `galois_linear_map`,
  `nimber_subfield_frobenius_linear_map`) that feed the existing
  `trace`/`char_poly`/`determinant`/`apply_outermorphism` methods.
  `spinor_rep` reaches the nonsingular nimber char-2 path as well
  as the supported char-0 path; `lazy_spinor_rep` returns a lightweight
  `LazySpinorRep` Python façade, and `apply_generator` / `apply_vector` expose
  the Rust lazy left-regular action beyond the
	  explicit matrix cap. The same scalar set also gets `DividedPowerAlgebra`/`DpVector`
	  pyclasses for Γ(V), the symmetric Hopf mirror of the exterior algebra. Ordinal
	  products inherit the Rust core's checked Kummer boundary; use the scalar
	  `nim_mul` when you need `None` instead of an algebraic product boundary error.
	  MV objects expose the Rust-name helpers `grade_part`, `versor_inverse`, and
	  `multivector_inverse`.
	  The Python `LazySpinorRep` object delegates back to the owning algebra, so the
	  existing per-backend type checks still guard every applied vector.
- **`forms.rs`** — classify / witt / dickson / springer bindings, `FiniteFieldForm`
  (the runtime odd-char Fp/Fpn form wrapper) and `Char2FiniteFieldForm` (the
  runtime `Fpn<2,N>` Arf/Witt/Brauer-Wall/isometry wrapper for supported
  `N=1..4`), the Brauer–Wall classes (`bw_class_real`,
  `bw_class_complex`, `bw_class_nimber`, `bw_class_ordinal`, `bw_class_finite_odd`),
  `classify_real`/`classify_complex`, raw represented-subdomain char-0 helpers
  (`surreal_signature`, `surcomplex_rank`), Rust-name base-field isometry helpers
  (`isometric_real`, `isometric_rational`, `isometric_surcomplex`,
  `isometric_nimber`),
  `classify_rational` (with typed `BaseField`,
  `RationalPlace`, and rational local Hasse `RationalPlaceInvariant` rows) plus
  `arf_ordinal_finite`, the unified finite-field `FiniteFieldClass`
  wrapper alongside the leg-specific finite classifiers, constructible
  `WittClass`/`WittClassG`/`WittClassError`/`BrauerWallClass` invariant records,
  typed Witt-error inspection helpers, fixed-precision p-adic Springer
  decompositions `springer_decompose_qp`, `springer_decompose_qq`,
  and `springer_decompose_laurent`, the Rust-name generic local dispatcher
  `springer_decompose_local` over already-bound local algebra objects,
  `springer_decompose_ramified_qp4_e2`,
  and `springer_decompose_ramified_qp4_e3`
  with named `ResidueForm`/`LocalResidueForm` layers,
  `hilbert_symbol_at`/`hasse_at_place`/`hilbert_product`/
  `hilbert_reciprocity_product`, checked p-adic helpers (`is_square_qp`,
  `hilbert_symbol_qp`, `is_isotropic_at_p`, `is_isotropic_q`), plus rational Brauer
  local-invariant helpers,
  `classify_finite_algebra`/`witt_finite_algebra`/`bw_class_finite_algebra`/
  `isometric_finite_algebra`, raw `arf_f2`, finite odd-characteristic
  discriminants, finite char-2 Artin-Schreier classes, the odd-characteristic
  `F_q(t)` local-global helpers (`monic_irreducible_factors`,
  `relevant_places`, `valuation_at`, `is_local_square`, `is_global_square_ff`,
  `hilbert_symbol_ff`, `hasse_at_place_ff`, `hilbert_reciprocity_product_ff`,
  `ramified_places_ff`, `is_isotropic_at_place`, `is_isotropic_ff`, and
  `isotropy_over_ff_adeles`) with named `FunctionFieldLocalIsotropy` adelic rows, the
  characteristic-2 Artin-Schreier function-field symbol helpers
  (`char2_monic_irreducible_factors`, `as_symbol_at`, `as_symbol_places`,
  `as_symbol_reciprocity_sum`, `as_symbol_ramified_places`, `global_is_pe`) and
  `Char2FunctionFieldForm`/`Char2LocalDecomp` with named `Char2PsiTerm` wild
  summands for the Aravire-Jacob
  decomposition plus local/global char-2 isotropy, finite characteristic-2 Arf
  and isometry helpers (`arf_nimber`, `arf_fpn_char2`, `arf_ordinal_finite`,
  `isometric_fpn_char2`, `isometric_ordinal_finite`) and char-2 function-field helpers
  (`relevant_places_char2`,
  `springer_decompose_local_char2`, `local_anisotropic_dim_char2`,
	  `local_is_isotropic_char2`, `is_isotropic_global_char2`),
	  `Char2FunctionFieldForm.from_blocks`,
	  `isotropy_over_adeles`/`AdelicIsotropy`, finite-form `isometric_to`,
	  `SymplecticForm.from_gram`,
	  `HermitianForm.from_gram`, the
  finite-prime-field numeric invariants (`level`, `pythagoras_number`,
  `u_invariant`, finite-field `hilbert_symbol`, sum-of-squares), the
  odd-characteristic finite-field form helpers (`classify_finite_odd`, `finite_odd_witt`,
  `discriminant_finite_odd`, `hasse_invariant_finite_odd`,
  `witt_decompose_finite_odd`, `isometric_finite_odd`, prime-field `is_square`,
  extension-field `is_square_finite`, `bw_class_finite_odd`, `e_staircase_finite_odd`), the quadric bench (`fit_f2_quadratic`,
  `QuadricFit`), trace/Gold-form helpers (`trace_twisted_form`,
  `trace_form_arf`, `gold_form_arf`, `gold_form`), the
	  named Witt/staircase result records (`RealWittDecomp`, `OddWittDecomp`,
	  `EnStaircase`), `WittClass` /
	  `WittClassG` checked metric constructors (`try_from_metric`,
	  `try_char2_from_metric`) plus
	  checked arithmetic through Python operators, the integral-lattice layer
  (`IntegralForm`, ADE constructors under Rust `a_n`/`d_n`/
  `e_6`/`e_7`/`e_8` names, `d16_plus`, `Genus`/`ScaleSymbol`,
  even-unimodular mass, lattice/code automorphism constants, theta series,
  rational/F₂ Clifford metrics, root-lattice predicates/Coxeter numbers, genus
  signature mod 8), `BinaryCode`/Construction
  A, exact modular q-expansion helpers (`E4`, `E6`, `Delta`, `as_modular_form`), and
  `DiscriminantForm`/Milgram/Weil `S`/`T` surfaces, including the raw
  `extended_golay_generator_rows` generator matrix.
- **`games.rs`** — `Game` / `NumberGame` / `NimberGame` (the char-2 transfinite
  Nim-heap mirror) / `GameExterior` (`GameRelation`, quotient-aware add/scale/wedge,
  relation-search certificates, explicit integer relations, and the raw free Grassmann algebra) /
  `Hackenbush` + typed `Color` values +
  `nim_mul_mex` / concrete coin-turning and Tartan probes (`coin_companions`,
  exact `singleton_companions`/`turtles_companions`, `coin_turning_grundy`,
  `coin_turning_tartan_grundy`, plus callback-backed Rust-name
  `grundy_1d`/`tartan_grundy`) /
  `grundy_graph` / callback-backed Rust-name `grundy` / `mex`; the kernel outcome surface
  (`outcomes`/`p_positions`/`scoring_values`, typed
  `Outcome` values, scoring intervals as named `ScoreInterval` records); the misère/
  octal surface (`nim_canonical`, `misere_nim_p_predicted`,
  callback-backed Rust-name `try_misere_is_n`/`misere_is_n`/`misere_is_p`,
  `nim_moves`, `octal_moves`, `octal_misere_quotient`, `Quotient` with class
  product/signature helpers, `AbstractGame`, module-level `misere_quotient`); and the loopy graph engine
  (`LoopyGraph` including Rust-name callback constructor `from_rule`, typed
  `LoopyNimber` values,
  `loopy_nim_values_certified` / `LoopyNimCertificate`, typed
  `PartizanOutcome`, callback-backed Rust-name
	  `loopy_decision_sets`/`loopy_quadric_probe`) plus the
	  `LoopyValue` stopper catalogue (`on`/`off`/`over`/`under`/`dud`). NumberGame
  exposes the bidirectional transfinite birthday/sign-expansion bridge plus Python
  arithmetic/order operators; `NimberGame` exposes the matching operators and
  `turning_corners`; `Game` exposes Python arithmetic/equality operators,
  canonical and structural fingerprints,
  the ordinary thermograph plus the
  named tropical mirror and cooled stops, the all-small/atomic-weight calculus
	  both as methods and as module-level Rust-name functions, and the exact
	  `Pl`/`Thermograph` object API exposes the Rust wall data and tropical wall
	  operations directly. `GameExterior` exposes `with_relation_search` and
	  `relation_search_certificate`. `Hackenbush` and `LoopyGraph` expose their
	  underlying edge/adjacency payloads as well as value analyses.

## Binding Policy

The Python surface is **runtime-friendly parity**: everything that is a plain
runtime type is bound. What remains Rust-only is structural rather than an
ordinary binding backlog:

- **Open-ended const-generic backends are not instantiable from Python.** The
   shipped finite-field monomorphs are bound (`Fp2`, `Fp3`, `Fp5`, `Fp7`, `Fp11`,
   `Fp13`, `F4`, `F8`, `F16`, `F9`, `F25`, `F27`, with matching Clifford and
   divided-power backends). The p-adic integer/field rows now have a fixed Python
   dispatch slice at precision `4` for `p ∈ {2,3,5,7,11,13}` (`Zp*_4`/`Qp*_4`,
   with matching Clifford and divided-power backends) plus `springer_decompose_qp`
   for the odd-residue Springer layer. The unramified Witt/field rows also have a
   fixed precision-`4` dispatch slice for `(p, residue degree) = (2,2), (2,3),
   (2,4), (3,2), (5,2), and (3,3)`, again with matching Clifford and
   divided-power backends; `springer_decompose_qq` exposes the odd-residue
   theorem cases among those cells. The Laurent row has fixed precision-`6`
   dispatch over `Rational`, `F_3`, `F_5`, `F_7`, `F_11`, `F_13`, `F_9`, `F_25`,
   and `F_27`, with matching Clifford/divided-power backends and
   `springer_decompose_laurent` for the odd finite residue cases. The ramified
   row has fixed quadratic and cubic `Ramified<Qp<P,4>,{2,3}>` dispatch for
   `P ∈ {2,3,5,7,11,13}`, again with matching Clifford/divided-power backends;
   `springer_decompose_ramified_qp4_e2` and `springer_decompose_ramified_qp4_e3`
   expose the odd-residue theorem cases.
   The Gauss row has fixed `Gauss<Qp<P,4>>` dispatch for the same primes, with
   matching Clifford/divided-power backends; there is intentionally no Springer
   helper for it because its residue field is `F_p(tbar)`, not a finite field.
   Remaining families such as arbitrary
   `Zp<P,K>`, arbitrary `Qp<P,K>`, arbitrary `Qq<P,N,F>`, arbitrary
   `WittVec<P,N,F>`, arbitrary `Laurent<S,K>`, arbitrary `Ramified<S,E>`, and
   arbitrary `Gauss<S>` still take *compile-time* parameters; Python is runtime. There is no general
   `Qp(p=5, k=20)` or `Qq(p=3, k=20, f=7)` Clifford scalar without a dispatch macro
   enumerating instances or a runtime redesign. Where a runtime entry point was
   worth it, the project already built one — `FiniteFieldForm`/
   `Char2FiniteFieldForm` (runtime Fp/Fpn form façades), `LocalQp` (a scalar-only
   runtime-prime p-adic cell, because it is not a Rust `Scalar`), and `Adele`
   (an adelic `Scalar`, with matching Clifford, divided-power, and CGA backends).
   The rest stay Rust-only.
   Consequently const-generic field invariants beyond the fixed Python dispatch
   set are also Rust-only. The prime
   field numeric invariants for `F_2`, `F_3`, `F_5`, `F_7`, `F_11`, and `F_13`
   are exposed under their Rust names (`level`, `pythagoras_number`,
   `u_invariant`, `is_sum_of_n_squares`).

- **The exposed function-field rows are fixed-dispatch, not open-ended
  const-generics.** `NimberPoly`/`NimberRationalFunction` expose the exact
  `F_{2^128}[t]` / `F_{2^128}(t)` row, including Clifford and divided-power
  backends, while `Fp*Poly`/`Fp*RationalFunction` expose the prime-field rows used
  by the fixed `GaussQp*_4` residue fields, also with matching Clifford and
  divided-power backends. The odd-characteristic finite-field
  `F_q(t)` local-global form layer is
  exposed through a runtime dispatch over the same shipped odd finite fields
  (`F_3`, `F_5`, `F_7`, `F_11`, `F_13`, `F_9`, `F_25`, `F_27`), with coefficients
  encoded by finite-field element index. The characteristic-2 Artin-Schreier
  symbol layer gets the same fixed-dispatch treatment over `F_{2^N}(t)` for
  `N=1..4`, again with coefficients encoded by finite-field element index; the
  same dispatch also exposes the char-2 Springer `Char2QuadForm` API as
  `Char2FunctionFieldForm`.

- **Name-level Rust-only leftovers are structural, not user-facing math.** Python
  exposes `Metric::map` as a same-backend algebra method; cross-backend metric
  maps still require a compile-time target scalar. `springer_decompose_local<K>`
  is exposed as a runtime dispatcher over the fixed Python local algebra set.
  `Poly<S>::min_coeff_valuation` is exposed by the module-level
  `min_coeff_valuation` typed coefficient-list dispatcher. Crate-private and
  `pub(super)` kernels (`solve`, `inverse_matrix`, `smith_normal_form`,
  `geom_product_blades`, `wedge_sign`, `merge`, `reduce_word`, p-adic square-class
  reducers, thermograph wall combinators, and similar helpers) stay Rust-only
  because they are implementation kernels behind already-bound high-level APIs,
  not stable Python objects.

## Rules

- **Never `use pyo3` outside this module; never make it non-optional.** A green
  `cargo test` does NOT compile this feature — after touching `py/` *or any core API
  the bindings call* (e.g. renaming a `Scalar`/`FiniteField` method), run
  `cargo check --features python` AND `cargo clippy --features python --all-targets`.
  Display changes don't surface in `cargo test` either — rebuild + run `demo.py`.
- **Per-backend, no mixing.** Each Python backend monomorphises the generic engine
  to one concrete scalar type. Mixing scalar worlds in one algebra raises
  `TypeError` by construction — intended; do not add a runtime-tagged "any scalar"
  path.
- **Python operators:** `*` geometric, `^` wedge, `<<`/`>>` left/right contraction,
  `~` reverse, `/` divide (scalar or versor), `**` power, `+`/`-`, `==`.
- The smoke test is `demo.py` (rebuild via `maturin develop` first); add a section
  there when you bind something new.
