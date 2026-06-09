//! Integral lattices: the arithmetic form world over `Z`.
//!
//! This submodule is the forms pillar's integral complement to field-valued
//! quadratic-form classification. It keeps the lattice object, ADE catalogue,
//! genus computation, and mass/Leech layer together while the parent
//! `forms` module re-exports both the modules and their public items flat.

pub mod genus;
pub mod lattice;
pub mod mass_formula;
pub mod root_lattices;

pub use genus::*;
pub use lattice::*;
pub use mass_formula::*;
pub use root_lattices::*;
