# AGENTS.md — `src/py/`

The PyO3 bindings, behind the `python` feature (`pyo3` is an optional dep;
`extension-module` only enabled here). This is the ONLY place `use pyo3` may
appear — keeping it out of the core is what stops `cargo test` linking libpython.

Python exposes Rust `u128`/`i128` payloads as normal Python integers. Keep that
width contract at the Rust boundary for form invariants, lattice entries, field
orders, game values, and budgets. `usize` is still required for dimensions,
indices, and PyO3 slots such as `__hash__`.

The exhaustive list of *which* concrete monomorphs are bound is data, not prose —
it lives in **`catalog.rs`**, not here. This file documents the structure and the
policy; consult `catalog.rs` for the actual instance set when you need it.

## Files

- **`mod.rs`** — the `#[pymodule]`; chains each submodule's `pub(crate) register()`
  (`scalars` → `engine` → `forms` → `games`) + the module-level `version()`. Declares
  `#[macro_use] mod catalog`.
- **`catalog.rs`** — the backend **manifest**: three `macro_rules!` lists
  (`py_engine_backends!`, `py_divided_power_backends!`, `py_cga_backends!`) that name
  every bound `<World>Algebra`/`MV`/`LinearMap`/`DpVector`/`Cga` monomorph and its
  `parse_*`/`wrap_*` hooks. `engine.rs` invokes these to stamp out the classes, so a
  new bound backend is added HERE, in one place. The single source of truth for the
  bound instance set; the most important structural file in `src/py/`.
- **`scalars.rs`** — the scalar pyclasses + constructors + nim-field free fns.
  Families bound: `Nimber`; the fixed finite fields (`Fp2`…`Fp13`, `F4`/`F8`/`F16`/
  `F9`/`F25`/`F27`); the fixed p-adic slices (`Zp*_4`/`Qp*_4` for p∈{2,3,5,7,11,13});
  the fixed unramified slices (`WittVec*_4_*`/`Qq*_4_*`); the fixed functor slices
  (`Laurent*_6`, `RamifiedQp*_4_E{2,3}`, `GaussQp*_4`); the function-field rows
  (`NimberPoly`/`NimberRationalFunction`, `Fp*Poly`/`Fp*RationalFunction`);
  `Rational`, `Surreal`, `Surcomplex`, `Integer`, `Omnific`, `Ordinal`,
  `SignExpansion`; the runtime cells `LocalQp`/`Adele`; and the tropical endpoints
  `MaxPlusTropical`/`MinPlusTropical`. Threads the `parse_*`/`wrap_*` hooks the
  `backend!` macro consumes. Bound scalar classes expose the shared runtime `Scalar`
  surface (`zero`/`one`/`characteristic` where applicable, `is_zero`, partial
  inverses/division, owned operators), plus per-world extras: Surreal's simplicity
  bridge + sign-expansion round-trips + lazy analytic helpers; Nimber's
  `pow`/`frobenius`/`sqrt` + the `nim_*` Galois toolkit; the Qp/Qq local-field package
  (`uniformizer`/`residue`/`residue_unit`/`teichmuller`, `is_integral`/`to_integer`)
  and the Qq `FieldExtension`/`CyclicGaloisExtension` surface; Ordinal's CNF terms +
  staged `nim_mul`/`checked_inv`; the `Fpn` Galois/reduction-poly metadata
  (`ReductionPolynomialKind`).
- **`engine.rs`** — the `backend!` macro → `<World>Algebra`/`<World>MV`/`<World>LinearMap`
  triplets (driven by `catalog.rs`), plus conformal GA (`<World>Cga`) over every bound
  char-0 scalar world with a matching MV carrier. MV methods cover the full GA suite
  (clifford_conjugate, scalar_product, commutator, anticommutator, undual, meet,
  is_blade, blade_subspace, factor_blade, cayley, cayley_inverse, spinor_norm, versor_grade_parity,
  classify_versor → `VersorClass`, plus raw `(blade_mask, coeff)` terms, `grade_part`,
  `versor_inverse`, `multivector_inverse`). Algebra methods add
  trace/char_poly/determinant/exterior_power_trace/apply_outermorphism/inverse_outermorphism, the typed
  `<World>LinearMap` pyclass, fixed-dispatch Frobenius/Galois map constructors
  (`frobenius_linear_map`, `galois_linear_map`,
  `nimber_subfield_frobenius_linear_map`), `spinor_rep`/`SpinorRep` (incl. the
  nonsingular nimber char-2 path), the lazy `lazy_spinor_rep`/`LazySpinorRep` (with
  `apply_generator`/`apply_vector` beyond the explicit matrix cap), the metric
  constructors/helpers (`general`/`grassmann`/`q`/`b_terms`/`a_terms`/`map`/`q_val`/…),
  the `pga(n)` constructor, the `gram`/`diagonalize`/`as_diagonal` façade, tensor-square
  + graded-tensor embeddings, and the Witt-ring representative ops (`tensor_form`,
  `pfister1`, `pfister`, `in_fundamental_ideal`). The same scalar set also gets
  `DividedPowerAlgebra`/`DpVector` (Γ(V), the symmetric Hopf mirror). Module-level
  `bits`/`grade` expose the mask utilities.
- **`forms.rs`** — the classifier / invariant / lattice bindings: `classify_real`/
  `classify_complex`/`classify_rational`, the Brauer–Wall classes, the runtime form
  façades `OddFiniteFieldForm` (Fp/Fpn) and `Char2FiniteFieldForm` (`Fpn<2,N>`,
  N=1..4), `FiniteFieldClass` + leg-specific classifiers, constructible
  `WittClass`/`WittClassG`/`BrauerWallClass` records, base-field isometry helpers, the
  Springer decompositions (the Surreal `springer_decompose`; `springer_decompose_qp`/
  `_qq`/`_laurent`/`_ramified_qp4_e{2,3}` + the generic `springer_decompose_local` and
  `springer_decompose_local_char2` dispatchers), the
  rational/p-adic local-global helpers (`hilbert_*`, `hasse_at_place`,
  `is_isotropic_q`, …), the odd `F_q(t)` layer (`try_*_ff`, `FunctionFieldLocalIsotropy`)
  and the char-2 Artin-Schreier layer (`as_symbol_*`, `Char2FunctionFieldForm`/
  `Char2LocalDecomp` with `Char2PsiTerm`, local/global char-2 isotropy), the
  symplectic/hermitian constructors, the field numeric invariants
  (`level`/`pythagoras_number`/`u_invariant`/sum-of-squares), the quadric bench
  (`fit_f2_quadratic`/`QuadricFit`), the trace/Gold-form helpers (`trace_twisted_form`,
  `trace_form_arf`, `gold_form_arf`, `gold_form`), and the integral-lattice layer
  (`IntegralForm`, the ADE constructors `a_n`/`d_n`/`e_6`/`e_7`/`e_8`/`d16_plus`,
  `Genus`/`ScaleSymbol`, mass/automorphism constants, `BinaryCode`/Construction A,
  theta + modular q-expansion helpers `eisenstein_e4`/`eisenstein_e6`/`delta`/`as_modular_form`,
  `DiscriminantForm`/Milgram/Weil `S`/`T`).
- **`games.rs`** — `Game`/`NumberGame`/`NimberGame`/`GameExterior`/`Hackenbush` +
  typed `Color`; the `GameExterior` relation surface (`GameRelation`,
  `GameRelationCertificate`, `RelationSearchCertificate`); `nim_mul_mex` + the
  coin-turning/Tartan probes; `grundy_graph`/`grundy`/`mex`; the kernel surface
  (`outcomes`/`p_positions`/`scoring_values`, typed `Outcome`, `ScoreInterval`); the
  misère/octal surface (`misere_quotient`, `Quotient`, `AbstractGame`, octal helpers);
  and the loopy engine (`LoopyGraph`, `LoopyNimber`,
  `loopy_nim_values_certified`/`LoopyNimCertificate`, `loopy_decision_sets`/
  `loopy_quadric_probe`, the `LoopyValue` stopper catalogue + typed `PartizanOutcome`).
  The games carry Python arithmetic/order operators, the thermograph + tropical-mirror
  + atomic-weight calculus, and the exact `Pl`/`Thermograph` wall API. Callback-backed
  Rust-name variants (`grundy`/`try_misere_is_n`/`loopy_quadric_probe`/…) accept a
  Python move-generator.

## Binding policy

The Python surface is **runtime-friendly parity**: everything that is a plain
runtime type is bound. What stays Rust-only is structural, not a backlog:

- **Open-ended const-generic backends are not instantiable from Python.** Python is
  runtime; `Qp<P,K>`, `Qq<P,N,F>`, `WittVec<P,N,F>`, `Laurent<S,K>`, `Ramified<S,E>`,
  `Gauss<S>`, `Fpn<P,N>`, … take *compile-time* parameters. The shipped Python set is
  a fixed dispatch slice (enumerated in `catalog.rs`): finite fields `Fp2`…`F27`;
  p-adic `Zp*_4`/`Qp*_4` (p∈{2,3,5,7,11,13}); unramified `WittVec*_4_*`/`Qq*_4_*`;
  Laurent `*_6`; ramified `RamifiedQp*_4_E{2,3}`; Gauss `GaussQp*_4`; each with matching
  Clifford/divided-power backends, and the odd-residue Springer helper where it
  applies. There is no `Qp(p=5, k=20)` Clifford scalar without extending the dispatch
  macro or a runtime redesign. Where a runtime entry point was worth it, the project
  built one: `OddFiniteFieldForm`/`Char2FiniteFieldForm` (runtime Fp/Fpn form façades),
  `LocalQp` (a scalar-only runtime-prime p-adic cell — it is not a Rust `Scalar`), and
  `Adele` (an adelic `Scalar`, with Clifford/divided-power/CGA backends).
- **The function-field rows are fixed-dispatch.** `NimberPoly`/`NimberRationalFunction`
  expose `F_{2^128}[t]`/`F_{2^128}(t)`; `Fp*Poly`/`Fp*RationalFunction` expose the
  prime-field rows the fixed `GaussQp*_4` residue fields use. The odd `F_q(t)`
  local-global layer dispatches over the shipped odd finite fields, and the char-2
  Artin-Schreier symbol layer over `F_{2^N}(t)` (N=1..4), coefficients encoded by
  finite-field element index.
- **CGA is absent where the math rejects it.** CGA needs `1/2`, so it is bound only for
  char-0 worlds with an MV carrier (Surreal/Rational/Adele/Surcomplex, the fixed
  `Qp*_4`/`Qq*_4_*`/`LaurentRational_6`/`RamifiedQp*`/`GaussQp*` cells). Absent for
  positive-characteristic worlds (incl. `LaurentFp*`/`LaurentFq*`), for
  `Integer`/`Omnific` (1/2 not invertible), and for runtime `LocalQp` (not a `Scalar`).
- **Name-level Rust-only leftovers are structural.** `Metric::map` is bound as a
  same-backend algebra method (cross-backend maps need a compile-time target);
  `springer_decompose_local<K>` is a runtime dispatcher over the fixed local set;
  `Poly<S>::min_coeff_valuation` is exposed via a typed coefficient-list dispatcher.
  Crate-private kernels (`solve`, `inverse_matrix`, `smith_normal_form`,
  `geom_product_blades`, `reduce_word`, the thermograph wall combinators, …) stay
  Rust-only — implementation kernels behind already-bound high-level APIs.

## Rules

- **Never `use pyo3` outside this module; never make it non-optional.** A green
  `cargo test` does NOT compile this feature — after touching `py/` *or any core API
  the bindings call* (e.g. renaming a `Scalar`/`FiniteField` method), run
  `cargo check --features python` AND `cargo clippy --features python --all-targets`.
  Display changes don't surface in `cargo test` either — rebuild + run `demo.py`.
- **Per-backend, no mixing.** Each Python backend monomorphises the generic engine to
  one concrete scalar type. Mixing scalar worlds in one algebra raises `TypeError` by
  construction — intended; do not add a runtime-tagged "any scalar" path.
- **Python operators:** `*` geometric, `^` wedge, `<<`/`>>` left/right contraction,
  `~` reverse, `/` divide (scalar or versor), `**` power, `+`/`-`, `==`.
- The smoke test is `demo.py` (rebuild via `maturin develop` first); add a section
  there when you bind something new, and a backend to `catalog.rs` when you add one.
