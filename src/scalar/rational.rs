//! Exact rational ℚ — *not* a game backend, just the char-0 scalar used to
//! validate the geometric-product engine against the known Cl(p,q)
//! classification before trusting the exotic backends. (The surreal backend is
//! the real char-0 home.)

use crate::scalar::Scalar;
use std::cmp::Ordering;
use std::fmt;

/// Exact rational over i128, used only for engine validation. Overflow is a
/// known limitation — fine for the small forms in the test suite, not meant
/// for serious arithmetic. (The surreal backend is the real char-0 home.)
#[derive(Clone)]
pub struct Rational {
    num: i128,
    den: i128, // always > 0, gcd(num, den) == 1
}

fn gcd(a: i128, b: i128) -> i128 {
    let (mut a, mut b) = (a.abs(), b.abs());
    while b != 0 {
        let t = b;
        b = a % b;
        a = t;
    }
    a
}

impl Rational {
    pub fn new(num: i128, den: i128) -> Self {
        assert!(den != 0, "zero denominator");
        let sign = if den < 0 { -1 } else { 1 };
        let (num, den) = (num * sign, den * sign);
        let g = gcd(num, den).max(1);
        Rational {
            num: num / g,
            den: den / g,
        }
    }

    pub fn int(n: i128) -> Self {
        Rational { num: n, den: 1 }
    }

    /// Sign as an Ordering relative to zero (den is always > 0).
    pub fn sign(&self) -> Ordering {
        self.num.cmp(&0)
    }

    /// True iff this rational is a (rational) integer, i.e. its denominator is 1.
    /// Used by the omnific-integer backend to test the constant CNF term.
    pub fn is_integer(&self) -> bool {
        self.den == 1
    }

    /// The numerator (in lowest terms; carries the sign).
    pub fn numer(&self) -> i128 {
        self.num
    }

    /// The denominator (in lowest terms; always > 0).
    pub fn denom(&self) -> i128 {
        self.den
    }

    /// Total order on values (denominator is always positive).
    pub fn cmp(&self, other: &Self) -> Ordering {
        self.sub(other).sign()
    }

    /// The greatest integer ≤ this rational.
    pub fn floor(&self) -> i128 {
        self.num.div_euclid(self.den)
    }
}

impl fmt::Debug for Rational {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.den == 1 {
            write!(f, "{}", self.num)
        } else {
            write!(f, "{}/{}", self.num, self.den)
        }
    }
}

impl PartialEq for Rational {
    fn eq(&self, other: &Self) -> bool {
        // both are in lowest terms with positive denominator
        self.num == other.num && self.den == other.den
    }
}

impl Scalar for Rational {
    fn zero() -> Self {
        Rational { num: 0, den: 1 }
    }
    fn one() -> Self {
        Rational { num: 1, den: 1 }
    }
    fn add(&self, rhs: &Self) -> Self {
        Rational::new(self.num * rhs.den + rhs.num * self.den, self.den * rhs.den)
    }
    fn neg(&self) -> Self {
        Rational {
            num: -self.num,
            den: self.den,
        }
    }
    fn mul(&self, rhs: &Self) -> Self {
        Rational::new(self.num * rhs.num, self.den * rhs.den)
    }
    fn characteristic() -> u128 {
        0
    }
    fn inv(&self) -> Option<Self> {
        if self.num == 0 {
            None
        } else {
            Some(Rational::new(self.den, self.num))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scalar::Scalar;

    #[test]
    fn rational_arithmetic() {
        let half = Rational::new(1, 2);
        let third = Rational::new(1, 3);
        assert_eq!(half.add(&third), Rational::new(5, 6));
        assert_eq!(half.mul(&third), Rational::new(1, 6));
        assert_eq!(half.sub(&half), Rational::zero());
        assert_eq!(half.add(&half), Rational::one());
        assert_eq!(Rational::new(2, 4), Rational::new(1, 2)); // reduction
    }
}
