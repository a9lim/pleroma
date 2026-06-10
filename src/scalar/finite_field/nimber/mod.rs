//! Nimbers: the ordinals under nim-addition and nim-multiplication, Conway's
//! Field On_2 of characteristic 2. Restricted here to `u128`, i.e. nimbers
//! below 2^128 — which *is* exactly the finite nim-field F_{2^128}, and contains
//! every smaller F_{2^{2^k}} (k <= 7: F_4, F_16, F_256, ... F_{2^128}).
//!
//! The implementation is split along the actual layers:
//!
//! * `arithmetic` — XOR addition, nim multiplication, Frobenius/sqrt, inverse.
//! * `artin_schreier` — trace and the `y² + y = c` solver.
//! * `galois` — degree, conjugates, minimal polynomial, relative trace/norm,
//!   multiplicative order, primitive elements, and discrete log.
//!
//! The public `nim_*` functions stay re-exported from this module, so callers can
//! keep using `scalar::nim_mul`, `scalar::nim_trace`, etc.

mod arithmetic;
mod artin_schreier;
mod galois;

pub use arithmetic::*;
pub use artin_schreier::*;
pub use galois::*;

use crate::scalar::Scalar;

/// A nimber, i.e. an element of On_2 truncated to F_{2^128}.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct Nimber(pub u128);

impl std::fmt::Debug for Nimber {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "*{}", self.0)
    }
}

impl Scalar for Nimber {
    fn zero() -> Self {
        Nimber(0)
    }
    fn one() -> Self {
        Nimber(1)
    }
    fn add(&self, rhs: &Self) -> Self {
        Nimber(nim_add(self.0, rhs.0))
    }
    fn neg(&self) -> Self {
        // characteristic 2: every element is its own additive inverse
        *self
    }
    fn mul(&self, rhs: &Self) -> Self {
        Nimber(nim_mul(self.0, rhs.0))
    }
    fn characteristic() -> u128 {
        2
    }
    fn inv(&self) -> Option<Self> {
        nim_inv(self.0).map(Nimber)
    }
}

#[cfg(test)]
mod tests;
