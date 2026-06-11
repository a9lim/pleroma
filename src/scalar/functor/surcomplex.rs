//! Surcomplex numbers: adjoin i (with i² = −1) to any scalar backend.
//!
//! Generic over `S`, so `Surcomplex<Surreal>` is the complexification of the
//! implemented finite-support surreal backend. `Surcomplex<Rational>` is the
//! Gaussian rationals `ℚ[i]`, handy for tests.
//!
//! Over a *characteristic-2* backend this construction is **degenerate**, and
//! the tool demonstrates exactly why: i² = −1 = 1, so (1+i)² = 1 + 2i + i²
//! = 1 + 0 + 1 = 0. `1+i` is a nonzero nilpotent ⇒ zero divisors ⇒ not a
//! field. This is the concrete reason surcomplex only does something useful
//! over characteristic-0 scalar worlds here. Full `On₂` is algebraically closed,
//! while the fixed-width `Nimber` backend is `F_{2^128}`; neither makes this
//! char-2 adjunction a field. (See `tests::nimber_surcomplex_is_degenerate`.)

use crate::scalar::Scalar;

#[derive(Clone, PartialEq)]
pub struct Surcomplex<S: Scalar> {
    pub re: S,
    pub im: S,
}

impl<S: Scalar> std::fmt::Display for Surcomplex<S> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.im.is_zero() {
            write!(f, "{}", self.re)
        } else if self.re.is_zero() {
            write!(f, "{}·i", self.im)
        } else {
            write!(f, "{} + {}·i", self.re, self.im)
        }
    }
}

impl<S: Scalar> std::fmt::Debug for Surcomplex<S> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        std::fmt::Display::fmt(self, f)
    }
}

impl<S: Scalar> Surcomplex<S> {
    pub fn new(re: S, im: S) -> Self {
        Surcomplex { re, im }
    }

    /// The imaginary unit i.
    pub fn i() -> Self {
        Surcomplex {
            re: S::zero(),
            im: S::one(),
        }
    }

    /// Complex conjugate a − b·i.
    pub fn conj(&self) -> Self {
        Surcomplex {
            re: self.re.clone(),
            im: self.im.neg(),
        }
    }
}

impl<S: Scalar> Scalar for Surcomplex<S> {
    fn zero() -> Self {
        Surcomplex {
            re: S::zero(),
            im: S::zero(),
        }
    }
    fn one() -> Self {
        Surcomplex {
            re: S::one(),
            im: S::zero(),
        }
    }
    fn add(&self, rhs: &Self) -> Self {
        Surcomplex {
            re: self.re.add(&rhs.re),
            im: self.im.add(&rhs.im),
        }
    }
    fn neg(&self) -> Self {
        Surcomplex {
            re: self.re.neg(),
            im: self.im.neg(),
        }
    }
    fn mul(&self, rhs: &Self) -> Self {
        // (a+bi)(c+di) = (ac − bd) + (ad + bc) i
        let ac = self.re.mul(&rhs.re);
        let bd = self.im.mul(&rhs.im);
        let ad = self.re.mul(&rhs.im);
        let bc = self.im.mul(&rhs.re);
        Surcomplex {
            re: ac.sub(&bd),
            im: ad.add(&bc),
        }
    }
    fn characteristic() -> u128 {
        // adjoining i does not change the characteristic
        S::characteristic()
    }
    fn inv(&self) -> Option<Self> {
        // (a+bi)^{-1} = (a − bi)/(a²+b²), valid when the norm a²+b² inverts in S.
        let n = self.re.mul(&self.re).add(&self.im.mul(&self.im));
        let ninv = n.inv()?;
        Some(Surcomplex {
            re: self.re.mul(&ninv),
            im: self.im.neg().mul(&ninv),
        })
    }
    fn is_zero(&self) -> bool {
        self.re.is_zero() && self.im.is_zero()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scalar::Nimber;
    use crate::scalar::Rational;

    type Gauss = Surcomplex<Rational>;

    fn g(re: i128, im: i128) -> Gauss {
        Surcomplex::new(Rational::int(re), Rational::int(im))
    }

    #[test]
    fn gaussian_arithmetic() {
        let i = Gauss::i();
        assert_eq!(i.mul(&i), g(-1, 0)); // i^2 = -1
        let one_plus_i = g(1, 1);
        assert_eq!(one_plus_i.mul(&one_plus_i), g(0, 2)); // (1+i)^2 = 2i
        assert_eq!(one_plus_i.mul(&one_plus_i.conj()), g(2, 0)); // |1+i|^2 = 2
    }

    #[test]
    fn gaussian_inverse() {
        let z = g(1, 1); // 1 + i
        assert_eq!(z.mul(&z.inv().unwrap()), Gauss::one());
        let z2 = g(3, -2);
        assert_eq!(z2.mul(&z2.inv().unwrap()), Gauss::one());
        assert!(Gauss::zero().inv().is_none());
    }

    #[test]
    fn nimber_surcomplex_is_degenerate() {
        // Over char 2: i^2 = -1 = 1, so (1+i)^2 = 0. Nonzero nilpotent ⇒
        // not a field. This is the theorem made executable.
        type NC = Surcomplex<Nimber>;
        let i = NC::i();
        assert_eq!(i.mul(&i), NC::one()); // i^2 = 1, not -1≠1
        let one_plus_i = NC::new(Nimber(1), Nimber(1));
        assert!(one_plus_i.mul(&one_plus_i).is_zero()); // (1+i)^2 = 0
        assert!(!one_plus_i.is_zero()); // but 1+i itself ≠ 0
    }
}
