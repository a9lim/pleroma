//! Characteristic-2 quadratic-form invariants.
//!
//! Characteristic 2 has two different but adjacent invariants:
//!
//! * [`arf`] classifies the quadratic form / Clifford algebra through the Arf
//!   invariant.
//! * [`dickson`] classifies orthogonal transformations by the Dickson invariant,
//!   the determinant replacement in characteristic 2.
//!
//! The public exports stay flat (`forms::arf_invariant`,
//! `forms::dickson_matrix`, …), matching the rest of the forms pillar.

mod arf;
mod dickson;

pub use arf::*;
pub use dickson::*;
