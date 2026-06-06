//! The integers ℤ as a `Scalar`: the coefficient ring for the exterior algebra
//! of the *game group* (`games/partizan.rs`). Games form an abelian group — a
//! ℤ-module — but not a ring, so an exterior algebra (which needs only a
//! commutative ring of coefficients and a module of generators) is exactly the
//! Clifford-adjacent structure that lives on *all* of game-world, not only the
//! field-like cores. Only `±1` are invertible, which is fine: the Grassmann
//! product never calls `inv`.

use crate::scalar::Scalar;
use std::fmt;

#[derive(Clone, Copy, PartialEq, Eq)]
pub struct Integer(pub i128);

impl fmt::Debug for Integer {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl Scalar for Integer {
    fn zero() -> Self {
        Integer(0)
    }
    fn one() -> Self {
        Integer(1)
    }
    fn add(&self, rhs: &Self) -> Self {
        Integer(self.0 + rhs.0)
    }
    fn neg(&self) -> Self {
        Integer(-self.0)
    }
    fn mul(&self, rhs: &Self) -> Self {
        Integer(self.0 * rhs.0)
    }
    fn characteristic() -> u128 {
        0
    }
    fn inv(&self) -> Option<Self> {
        match self.0 {
            1 | -1 => Some(*self),
            _ => None, // ℤ has only the units ±1
        }
    }
}
