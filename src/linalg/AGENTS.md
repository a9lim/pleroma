# AGENTS.md — `src/linalg/`

Crate-private shared linear algebra, deliberately placed BELOW the mathematical
pillars rather than exposed as a public API. `mod.rs` is `pub(crate)` only.

Fixed-width arithmetic payloads here are `u128`/`i128`; `usize` is only for matrix
dimensions and indices. Keep relation rows, Smith/Hermite pivots, and integer
solver data on the repo-wide width contract.

- **`field.rs`** — Gaussian `solve` / `inverse_matrix` / `unit_pivot_nullspace` over
  any `Scalar` (a field, a local ring, or a precision model): the kernels pivot only on
  entries whose `Scalar::inv` exists, so over a field it is ordinary Gauss-Jordan and
  over a ring they return `None` when a required nonunit pivot appears. The crate's
  generic linear solver: used by the Clifford GA solves (`clifford::multivector_inverse`,
  blade analysis, `inverse_outermorphism`, the spinor builder) and the
  integral-lattice/symplectic/ramified layers.
- **`f2.rs`** — `nim_rank`: Gauss-Jordan row rank over `F_{2^128}` (concrete `u128`
  nimbers), the characteristic-2 row-kernel primitive.
- **`integer.rs`** — exact integer linear algebra over ℤ:
  - `normalize_relation_rows` (the crate's row **Hermite normal form**: increasing
    leading columns, positive pivots, zeros below each pivot, above-pivot entries
    reduced mod the pivot) + `reduce_integer_vector`. `normalize_relation_rows` is
    consumed by the integral-lattice layer (`forms/integral/`); `reduce_integer_vector`
    by the game exterior algebra's lattice quotient (`games/game_exterior/`).
  - `ext_gcd` (Bézout `a·x + b·y = gcd`) and `smith_normal_form` (invariant factors
    `d₀ | d₁ | …` via unimodular `ext_gcd`-based row/column combines; `∏ dᵢ = |det|`,
    cokernel `ℤⁿ/Mℤⁿ ≅ ⨁ ℤ/dᵢ`). Used by the integral-lattice layer:
    `forms/integral/lattice/` reads invariant factors off SNF;
    `forms/integral/discriminant/` enumerates `Z^n/GZ^n` representatives via the
    normalized relation rows.

This module is matrix-heavy and walks index-parallel arrays; the crate-level
`#![allow(clippy::needless_range_loop)]` in `lib.rs` covers it (explicit indices
read clearer than iterator adapters when the body indexes several arrays by the
same `i`).
