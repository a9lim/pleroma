//! The integers ℤ as a `Scalar`: the coefficient ring for the exterior algebra
//! of the *game group* (`games/partizan.rs`). Games form an abelian group — a
//! ℤ-module — but not a ring, so an exterior algebra (which needs only a
//! commutative ring of coefficients and a module of generators) is exactly the
//! Clifford-adjacent structure that lives on *all* of game-world, not only the
//! field-like cores. Only `±1` are invertible, which is fine: the Grassmann
//! product never calls `inv`.

use crate::scalar::Scalar;
use std::fmt;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn integer_basic_ops() {
        let a = Integer(5);
        let b = Integer(3);
        assert_eq!(a.add(&b), Integer(8));
        assert_eq!(a.mul(&b), Integer(15));
        assert_eq!(a.neg(), Integer(-5));
        assert_eq!(Integer(i128::MAX).neg(), Integer(i128::MIN + 1));
    }

    #[test]
    #[should_panic(expected = "Integer addition overflowed i128")]
    fn integer_add_overflows_loudly() {
        let _ = Integer(i128::MAX).add(&Integer(1));
    }

    #[test]
    #[should_panic(expected = "Integer negation overflowed i128")]
    fn integer_neg_overflows_loudly() {
        // i128::MIN has no positive counterpart in i128
        let _ = Integer(i128::MIN).neg();
    }

    #[test]
    #[should_panic(expected = "Integer multiplication overflowed i128")]
    fn integer_mul_overflows_loudly() {
        let _ = Integer(i128::MAX).mul(&Integer(2));
    }
}

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
        Integer(
            self.0
                .checked_add(rhs.0)
                .expect("Integer addition overflowed i128"),
        )
    }
    fn neg(&self) -> Self {
        Integer(
            self.0
                .checked_neg()
                .expect("Integer negation overflowed i128"),
        )
    }
    fn mul(&self, rhs: &Self) -> Self {
        Integer(
            self.0
                .checked_mul(rhs.0)
                .expect("Integer multiplication overflowed i128"),
        )
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
