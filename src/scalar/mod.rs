//! The scalar interface every Clifford backend implements.
//!
//! A Clifford algebra needs a *commutative ring* of scalars. The whole point
//! of this project is that combinatorial games only supply such a ring on
//! their field-like subclasses — nimbers, surreals, surcomplex — so each of
//! those is a `Scalar` impl, and the multivector engine in `clifford.rs` is
//! written once, generic over this trait.
//!
//! `Rational` lives here too: it is *not* a headline game-backend, just an
//! exact characteristic-0 scalar used to validate the geometric-product
//! engine against the known Cl(p,q) classification before we trust the
//! exotic backends.

// The coefficient worlds. Each is a commutative-ring `Scalar` backend; this
// module is the shared trait plus the exact `Rational`/`Integer` used to
// validate the engine. Re-exported flat so call sites read `scalar::Nimber`.
pub mod fp;
pub mod nimber;
pub mod omnific;
pub mod onag;
pub mod surcomplex;
pub mod surreal;

pub use fp::*;
pub use nimber::*;
pub use omnific::*;
pub use onag::*;
pub use surcomplex::*;
pub use surreal::*;

use std::cmp::Ordering;
use std::fmt;
use std::fmt::Debug;

pub trait Scalar: Clone + PartialEq + Debug {
    fn zero() -> Self;
    fn one() -> Self;
    fn add(&self, rhs: &Self) -> Self;
    fn neg(&self) -> Self;
    fn mul(&self, rhs: &Self) -> Self;

    /// Field characteristic: 0 for surreal / surcomplex / rational, 2 for
    /// nimbers. The engine reads this because char 2 is genuinely different:
    /// `-1 == 1` so blade-reordering signs vanish, and the quadratic form Q
    /// (the squares) must be carried independently of the alternating
    /// off-diagonal bilinear form B — see the notes in `clifford.rs`.
    fn characteristic() -> u32;

    /// Multiplicative inverse, or `None` if not invertible (zero) or not
    /// finitely representable in this backend (e.g. a non-monomial surreal,
    /// whose inverse is an infinite Hahn series).
    fn inv(&self) -> Option<Self>;

    fn is_zero(&self) -> bool {
        *self == Self::zero()
    }

    fn sub(&self, rhs: &Self) -> Self {
        self.add(&rhs.neg())
    }
}

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
    fn characteristic() -> u32 {
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

/// The integers ℤ as a `Scalar`. Used as the coefficient ring for the exterior
/// algebra of the *game group* (`partizan.rs`): games form an abelian group — a
/// ℤ-module — but not a ring, so an exterior algebra (which needs only a
/// commutative ring of coefficients and a module of generators) is exactly the
/// Clifford-adjacent structure that lives on *all* of game-world, not only the
/// field-like cores. Only `±1` are invertible, which is fine: the Grassmann
/// product never calls `inv`.
#[derive(Clone, Copy, PartialEq, Eq)]
pub struct Integer(pub i64);

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
    fn characteristic() -> u32 {
        0
    }
    fn inv(&self) -> Option<Self> {
        match self.0 {
            1 | -1 => Some(*self),
            _ => None, // ℤ has only the units ±1
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
