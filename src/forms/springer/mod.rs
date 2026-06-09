//! Springer decomposition — the valuation-graded (local↔global) decomposition of
//! a quadratic form across the complete valued fields, with one generic engine
//! and the surreal odd-one-out.
//!
//! Over a complete discretely-valued field `K` with residue field `k`, a form
//! splits by the valuation of its diagonal entries into residue layers; how many
//! layers survive is controlled by whether the value group is 2-divisible:
//!
//! ```text
//!  K                value group   residue k   layers (Springer)
//!  No (surreal)     2-divisible   ℝ           1   — W(No) = W(ℝ) = ℤ
//!  Q_p / Q_q        ℤ             F_p / F_q   2   — W = W(k)²
//!  F_q((t))         ℤ             F_q         2   — W = W(k)²
//! ```
//!
//! The discretely-valued legs share **one** engine, [`springer_decompose_local`]
//! (in [`local`]), keyed off the [`ResidueField`](crate::scalar::ResidueField)
//! trait — the residue field `k` is read through the trait, never hardcoded:
//!
//! * [`padic`] — the mixed-characteristic entry points
//!   [`springer_decompose_qp`] (`Q_p`, residue `F_p`) and
//!   [`springer_decompose_qq`] (`Q_q`, residue `F_q`).
//! * [`laurent`] — the equal-characteristic entry point
//!   [`springer_decompose_laurent`] (`F_q((t))`, residue `F_q`).
//! * [`char2`] — the equal-characteristic-**2** mirror,
//!   [`springer_decompose_local_char2`]: the Aravire–Jacob `(φ₀, ψ, φ₁)`
//!   three-layer decomposition (the wild `R_π` summand the naive `W = W(k)²`
//!   grading misses), plus global isotropy over `F_q(t)` itself.
//! * [`surreal`] — [`springer_decompose`] over the surreals (char 0, residue ℝ),
//!   the ONE that does *not* fit the generic engine: its value group is
//!   2-divisible, so the second residue map collapses and `W(No) = W(ℝ) = ℤ`. It
//!   keeps its own engine; that mismatch *is* the local–global symmetry, not a gap.
//!
//! Children are private modules re-exported flat, so the public API stays shallow
//! (`forms::springer_decompose_qp`, `forms::springer_decompose_local`, …).

mod char2;
mod laurent;
mod local;
mod padic;
mod surreal;

pub use char2::*;
pub use laurent::*;
pub use local::*;
pub use padic::*;
pub use surreal::*;
