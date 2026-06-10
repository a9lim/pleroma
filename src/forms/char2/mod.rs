//! Characteristic-2 quadratic-form invariants.
//!
//! Characteristic 2 has two different but adjacent invariants:
//!
//! * `arf` classifies the quadratic form / Clifford algebra through the Arf
//!   invariant.
//! * `dickson` classifies orthogonal transformations by the Dickson invariant,
//!   the determinant replacement in characteristic 2.
//!
//! plus `field`, the [`FiniteChar2Field`] capability trait — the additive
//! (Artin–Schreier) mirror of [`FiniteOddField`](crate::forms::FiniteOddField)
//! that the char-2 local–global layer is generic over.
//!
//! The public exports stay flat (`forms::arf_invariant`,
//! `forms::dickson_matrix`, `forms::FiniteChar2Field`, …), matching the rest of the
//! forms pillar.

mod arf;
mod dickson;
mod field;

pub use arf::*;
pub use dickson::*;
pub use field::*;
