//! The p-primary Milgram/Brown Gauss-sum phase projection of a finite quadratic module.
//!
//! `FqmPrimaryPhase` and `FqmGaussPhase` are the public types that carry the phase
//! decomposition. The computation lives in `form.rs` (which has access to the full
//! group tables); these types are separated here so that modules importing only the
//! type records do not need the full cyclotomic arithmetic.

/// One p-primary Milgram/Brown phase of a finite quadratic module.
///
/// This is the **Gauss-sum phase projection** of the finite-quadratic-module Witt
/// class, not Wall's full generator-and-relation normal form.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FqmPrimaryPhase {
    /// The prime `p` of the primary subgroup.
    pub prime: u128,
    /// The cardinality of the p-primary subgroup.
    pub order: usize,
    /// The largest order of an element in this p-primary subgroup.
    pub exponent: u128,
    /// The normalized Gauss-sum phase `ζ_8^phase_mod8`.
    pub phase_mod8: i128,
}

/// The Milgram/Brown `Z/8` phase projection of a finite quadratic module.
///
/// The full Witt group of finite quadratic modules has finer Wall/Nikulin/
/// Kawauchi-Kojima generator data. This record intentionally exposes only the
/// p-local normalized Gauss-sum phases and their total.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FqmGaussPhase {
    /// The cardinality of the full finite quadratic module.
    pub order: usize,
    /// The total phase, i.e. the value congruent to the lattice signature mod 8.
    pub phase_mod8: i128,
    /// The p-primary phase factors whose sum is `phase_mod8` in `Z/8`.
    pub primary: Vec<FqmPrimaryPhase>,
}
