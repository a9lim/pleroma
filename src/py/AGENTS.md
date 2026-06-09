# AGENTS.md — `src/py/`

The PyO3 bindings, behind the `python` feature (`pyo3` is an optional dep;
`extension-module` only enabled here). This is the ONLY place `use pyo3` may
appear — keeping it out of the core is what stops `cargo test` linking libpython.

`mod.rs` is the `#[pymodule]`; it chains each submodule's `pub(crate) register()`.
Split per pillar:

- **`scalars.rs`** — the scalar pyclasses (`Nimber`, `Surreal`, `Surcomplex`,
  `Integer`, `Omnific`, `Ordinal`) + constructors + nim-field free fns. Threads the
  `parse_*`/`wrap_*` hooks the `backend!` macro consumes. Surreal also exposes the
  simplicity bridge; Nimber exposes `pow`/`__pow__`/`frobenius`/`sqrt` alongside the
  `nim_*` Galois free fns.
- **`engine.rs`** — the `backend!` macro → `<World>Algebra` + `<World>MV` pairs
  (Nimber/Surreal/Surcomplex/Integer/Omnific) + conformal GA (`Cga`). MV methods
  cover the full Arc-C suite (clifford_conjugate, scalar_product, commutator,
  anticommutator, undual, meet, is_blade, factor_blade, cayley); algebra methods add
  trace/char_poly/determinant/outermorphism/spinor_rep. `spinor_rep` reaches the
  nonsingular nimber char-2 path as well as the supported char-0 path.
- **`forms.rs`** — classify / witt / dickson / springer bindings, `FiniteFieldForm`
  (the RUNTIME Fp/Fpn form wrapper — the pattern that sidesteps const generics for
  the odd-char leg), the Brauer–Wall classes (`bw_class_real`, `bw_class_complex`,
  `bw_class_nimber`, `bw_class_oddchar`), `classify_real`/`classify_complex`,
  `hilbert_product`, and `isotropy_over_adeles`/`AdelicIsotropy`.
- **`games.rs`** — `Game` / `NumberGame` / `NimberGame` (the char-2 transfinite
  Nim-heap mirror) / `GameExterior` / `Hackenbush` +
  `nim_mul_mex` / `grundy_graph` / `mex`; the kernel outcome surface
  (`outcomes`/`p_positions`/`scoring_values`, Win/Loss/Draw as strings); the misère/
  octal surface (`nim_canonical`, `misere_nim_p_predicted`, `nim_moves`,
  `octal_moves`, `octal_misere_quotient`, `Quotient`, `AbstractGame`); and the loopy
  graph engine (`LoopyGraph`, `loopy_nim_values` with Side→None).

## Binding policy & the two structural walls

The Python surface is **runtime-friendly parity**: everything that is a plain
runtime type is bound. Two classes of thing are deliberately left Rust-only —
not gaps, structural limits documented here so nobody re-audits them:

1. **Const-generic backends are not instantiable from Python.** `Zp<P,K>`,
   `Qp<P,K>`, `Fp<P>`, `Fpn<P,N>`, `Qq<P,N,F>`, `WittVec<P,N,F>`, `Laurent<S,K>`,
   `Ramified<S,E>`, `Gauss<S>` all take *compile-time* parameters; Python is
   runtime. There is no `Qp(p=5, k=20)` without a dispatch macro enumerating
   instances or a runtime redesign. Where a runtime entry point was worth it, the
   project already built one — `FiniteFieldForm` (runtime Fp/Fpn for the odd-char
   forms leg) and `scalar::LocalQp` (runtime-prime p-adic). The rest stay Rust-only.
   Consequently `springer_padic`/`springer_laurent` and the const-generic field
   invariants (`level`/`pythagoras_number`/`u_invariant`) are also Rust-only.
2. **Closure-taking higher-order fns need callback adapters.** `grundy_1d`,
   `tartan_grundy`, the generic `grundy`, and the generic `misere_is_*`/
   `loopy_decision_sets`/`loopy_quadric_probe` take Rust closures. The concrete
   specializations ARE bound (`grundy_graph`, `nim_moves`, `octal_moves`,
   `LoopyGraph`, `octal_misere_quotient`); the closure-generic forms are not.

Other deliberate omissions:

- **No `Rational` pyclass.** The Rust `Rational` is an engine-validation-only scalar
  (overflow-prone i128); `Surreal` is the intended char-0 surface and subsumes ℚ
  (every rational embeds via `Surreal::from_rational`). The `rational(num, den)`
  helper returns a `Surreal`. So ℚ char-0 Clifford goes through `SurrealAlgebra`,
  not a separate `RationalAlgebra`.
- **`divided_power.rs`** (Γ(V), the symmetric mirror of the exterior Hopf algebra)
  is a standalone parallel algebra and stays Rust-only — binding it means a second
  `DividedPowerAlgebra` + `DpVector` class hierarchy, not worth it.
- **The `LoopyValue` stopper catalogue** (on/off/over/under/dud hand-arithmetic) is
  Rust-only; the computational loopy core (`LoopyGraph`, `loopy_nim_values`) is bound.
- **Full `HermitianForm`/`SymplecticForm`** need Gram-matrix-of-scalars construction
  in Python; only the signature/classify read-outs would be ergonomic, so they stay
  Rust-only for now.

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
