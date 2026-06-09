//! The invariant *groups* of quadratic forms: the Witt group, the Witt ring, and
//! the Brauer–Wall group — the three abelian groups the classifiers land in.
//!
//! * [`class`] — the Witt **group** `W_q(F)`: [`WittClass`] (the order-2 group of
//!   a finite nim-field, Arf-classified) and [`WittClassG`], the Char0/OddChar/Char2
//!   trichotomy enum that the classifier façade returns. (The Char2 leg is a
//!   *module*, not a ring — its `mul` panics; see [`ring`] for why.)
//! * [`ring`] — the Witt **ring**: [`tensor_form`], Pfister forms, the fundamental
//!   ideal `Iⁿ`, and the `eₙ` staircase (`e0 = dim`, `e1 = disc`, `e2 = Hasse`),
//!   with per-field stabilisation (`I² = 0` over `F_q`; the infinite ℝ tower).
//! * [`brauer_wall`] — the Brauer–Wall group `BW(F)`: [`bw_class_real`] (the Bott
//!   index `(q−p) mod 8`, so `BW(ℝ) ≅ ℤ/8`), [`bw_class_complex`] (`ℤ/2`),
//!   [`bw_class_oddchar`] (order-4, `≅ W(F_q)`), and [`bw_class_nimber`] (the
//!   char-2 Arf/Witt class `ℤ/2`, nonsingular metrics only). The law is the
//!   graded tensor product.
//!
//! The mod-8 spine lives here: `BW(ℝ) ≅ ℤ/8` is the same periodicity as the char-0
//! 8-fold Clifford table, Bott periodicity, and `E₈` as the rank-8 even unimodular
//! lattice (see [`integral`](crate::forms::integral)).
//!
//! Children are private modules re-exported flat, so the public API stays shallow
//! (`forms::WittClassG`, `forms::tensor_form`, `forms::bw_class_real`, …). The
//! numeric field invariants (level, u-invariant) the ring *implies* live separately
//! in [`field_invariants`](crate::forms::field_invariants).

mod brauer_wall;
mod class;
mod ring;

pub use brauer_wall::*;
pub use class::*;
pub use ring::*;
