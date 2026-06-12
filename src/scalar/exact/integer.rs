//! The integers ℤ as a `Scalar`: the coefficient ring for the exterior algebra
//! of the *game group* (`games/partizan.rs`). Games form an abelian group — a
//! ℤ-module — but not a ring, so an exterior algebra (which needs only a
//! commutative ring of coefficients and a module of generators) is exactly the
//! Clifford-adjacent structure that lives on *all* of game-world, not only the
//! field-like cores. Only `±1` are invertible, which is fine: the Grassmann
//! product never calls `inv`.

use crate::scalar::Scalar;
use std::cmp::Ordering;
use std::fmt;

/// Failure mode for exact Euclidean division in [`Integer`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum IntegerDivExactError {
    /// Division by zero.
    DivisionByZero,
    /// The divisor was nonzero but did not divide exactly; carries the
    /// Euclidean remainder `0 <= r < |divisor|`.
    Remainder(Integer),
}

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

    #[test]
    fn euclidean_division_uses_nonnegative_remainders() {
        assert_eq!(
            Integer(7).divrem(&Integer(3)),
            Some((Integer(2), Integer(1)))
        );
        assert_eq!(
            Integer(-7).divrem(&Integer(3)),
            Some((Integer(-3), Integer(2)))
        );
        assert_eq!(
            Integer(7).divrem(&Integer(-3)),
            Some((Integer(-2), Integer(1)))
        );
        assert_eq!(Integer(7).rem(&Integer(0)), None);
    }

    #[test]
    fn exact_division_reports_the_remainder() {
        assert_eq!(Integer(6).div_exact(&Integer(3)), Ok(Integer(2)));
        assert_eq!(
            Integer(7).div_exact(&Integer(3)),
            Err(IntegerDivExactError::Remainder(Integer(1)))
        );
        assert_eq!(
            Integer(7).div_exact(&Integer(0)),
            Err(IntegerDivExactError::DivisionByZero)
        );
    }

    #[test]
    fn standard_order_is_the_integer_order() {
        assert!(Integer(-2) < Integer(5));
        assert_eq!(
            std::cmp::Ord::cmp(&Integer(4), &Integer(4)),
            Ordering::Equal
        );
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub struct Integer(pub i128);

impl Integer {
    /// Euclidean division `self = q * divisor + r`, with `0 <= r < |divisor|`.
    ///
    /// Returns `None` for division by zero. Quotient overflow (`i128::MIN / -1`)
    /// panics like the rest of this fixed-width backend's arithmetic.
    pub fn divrem(&self, divisor: &Self) -> Option<(Self, Self)> {
        if divisor.0 == 0 {
            return None;
        }
        let q = self
            .0
            .checked_div_euclid(divisor.0)
            .expect("Integer Euclidean quotient overflowed i128");
        let r = self
            .0
            .checked_rem_euclid(divisor.0)
            .expect("Integer Euclidean remainder overflowed i128");
        Some((Integer(q), Integer(r)))
    }

    /// Euclidean remainder `self mod divisor`, with `0 <= r < |divisor|`.
    ///
    /// Returns `None` for division by zero.
    pub fn rem(&self, divisor: &Self) -> Option<Self> {
        self.divrem(divisor).map(|(_, r)| r)
    }

    /// Exact Euclidean division, returning the quotient iff the remainder is
    /// zero. Non-exact division carries the remainder for caller diagnostics.
    pub fn div_exact(&self, divisor: &Self) -> Result<Self, IntegerDivExactError> {
        let (q, r) = self
            .divrem(divisor)
            .ok_or(IntegerDivExactError::DivisionByZero)?;
        if r.is_zero() {
            Ok(q)
        } else {
            Err(IntegerDivExactError::Remainder(r))
        }
    }
}

impl From<i128> for Integer {
    /// The ℤ-embedding: the identity homomorphism ℤ → ℤ.
    fn from(n: i128) -> Self {
        Integer(n)
    }
}

impl fmt::Display for Integer {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl fmt::Debug for Integer {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(self, f)
    }
}

impl PartialOrd for Integer {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(std::cmp::Ord::cmp(self, other))
    }
}

impl Ord for Integer {
    fn cmp(&self, other: &Self) -> Ordering {
        self.0.cmp(&other.0)
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
    /// Faster direct construction; semantically identical to the default double-and-add.
    fn from_int(n: i128) -> Self {
        Integer(n)
    }
}
