//! Prime fields `F_p` of odd characteristic — a comparison backend.
//!
//! Like `Rational` and `Integer`, this is not a game-world core; it is here to
//! complete the **characteristic trichotomy** the rest of the library realizes:
//! char 0 (surreal/surcomplex, classified by signature → matrix algebra), char 2
//! (nimbers, classified by the Arf invariant), and now **odd characteristic**
//! (classified by dimension + discriminant; see `forms::oddchar`). Putting `F_p` in the
//! same generic `Scalar` engine lets the odd-char classifier run on the very same
//! `Metric`/`CliffordAlgebra` machinery.
//!
//! ## The const-generic modulus
//!
//! `Scalar::zero()`/`one()` take no `self`, so the modulus cannot live in the
//! value alone. We carry it in the **type**: `Fp<P>` is the field of `P`
//! elements. A different prime is a different type — exactly the per-backend,
//! no-mixing discipline the rest of the crate already uses (you cannot
//! accidentally add an `Fp<3>` to an `Fp<5>`). `P` must be prime; scalar
//! operations assert this instead of silently turning field-theory APIs into
//! arithmetic over `Z/PZ`.
//!
//! Unlike the nimbers, `neg` here is a *genuine* negation (`P − a ≠ a` for
//! `a ≠ 0`), so the Clifford antisymmetry signs are real — a useful contrast to
//! the char-2 backend where `−1 = 1`.

use crate::scalar::{is_prime_u128, Scalar};
use std::fmt;

/// An element of the prime field `F_P` (invariant: `0 ≤ value < P`).
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct Fp<const P: u128>(u128);

pub(crate) fn add_mod<const P: u128>(a: u128, b: u128) -> u128 {
    debug_assert!(P > 0 && a < P && b < P);
    if a >= P - b {
        a - (P - b)
    } else {
        a + b
    }
}

pub(crate) fn mul_mod<const P: u128>(mut a: u128, mut b: u128) -> u128 {
    debug_assert!(P > 0 && a < P && b < P);
    let mut acc = 0u128;
    while b > 0 {
        if b & 1 == 1 {
            acc = add_mod::<P>(acc, a);
        }
        b >>= 1;
        if b > 0 {
            a = add_mod::<P>(a, a);
        }
    }
    acc
}

impl<const P: u128> Fp<P> {
    pub fn modulus_is_prime() -> bool {
        is_prime_u128(P)
    }

    pub fn assert_prime_modulus() {
        assert!(Self::modulus_is_prime(), "Fp<P> needs prime P, got {P}");
    }

    /// Reduce an integer (possibly negative) into `F_P`.
    pub fn new(n: i128) -> Self {
        Self::assert_prime_modulus();
        let v = if n >= 0 {
            (n as u128) % P
        } else {
            let r = n.unsigned_abs() % P;
            if r == 0 {
                0
            } else {
                P - r
            }
        };
        Fp(v)
    }

    /// Reduce an unsigned integer into `F_P`.
    pub fn from_u128(n: u128) -> Self {
        Self::assert_prime_modulus();
        Fp(n % P)
    }

    /// The canonical representative in `[0, P)`.
    pub fn value(self) -> u128 {
        self.0
    }
}

impl<const P: u128> fmt::Debug for Fp<P> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl<const P: u128> Scalar for Fp<P> {
    fn zero() -> Self {
        Self::assert_prime_modulus();
        Fp(0)
    }
    fn one() -> Self {
        Self::assert_prime_modulus();
        Fp(1 % P)
    }
    fn add(&self, rhs: &Self) -> Self {
        Self::assert_prime_modulus();
        Fp(add_mod::<P>(self.0, rhs.0))
    }
    fn neg(&self) -> Self {
        Self::assert_prime_modulus();
        if self.0 == 0 {
            Fp(0)
        } else {
            Fp(P - self.0)
        }
    }
    fn mul(&self, rhs: &Self) -> Self {
        Self::assert_prime_modulus();
        Fp(mul_mod::<P>(self.0, rhs.0))
    }
    fn characteristic() -> u128 {
        Self::assert_prime_modulus();
        P
    }
    fn inv(&self) -> Option<Self> {
        Self::assert_prime_modulus();
        if self.0 == 0 {
            return None;
        }
        Some(self.pow(P - 2))
    }
}

impl<const P: u128> Fp<P> {
    pub fn pow(&self, mut e: u128) -> Self {
        Self::assert_prime_modulus();
        let mut base = *self;
        let mut acc = Self::one();
        while e > 0 {
            if e & 1 == 1 {
                acc = acc.mul(&base);
            }
            e >>= 1;
            if e > 0 {
                base = base.mul(&base);
            }
        }
        acc
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::clifford::{CliffordAlgebra, Metric};

    fn elems<const P: u128>() -> Vec<Fp<P>> {
        (0..P).map(Fp::<P>::from_u128).collect()
    }

    fn check_field_axioms<const P: u128>() {
        let es = elems::<P>();
        for &a in &es {
            for &b in &es {
                // commutativity
                assert_eq!(a.add(&b), b.add(&a));
                assert_eq!(a.mul(&b), b.mul(&a));
                for &c in &es {
                    // associativity + distributivity
                    assert_eq!(a.add(&b).add(&c), a.add(&b.add(&c)));
                    assert_eq!(a.mul(&b).mul(&c), a.mul(&b.mul(&c)));
                    assert_eq!(a.mul(&b.add(&c)), a.mul(&b).add(&a.mul(&c)));
                }
            }
            // additive identity/inverse
            assert_eq!(a.add(&Fp::<P>::zero()), a);
            assert_eq!(a.add(&a.neg()), Fp::<P>::zero());
            // multiplicative inverse for nonzero
            if !a.is_zero() {
                let ai = a.inv().expect("nonzero is invertible in a field");
                assert_eq!(a.mul(&ai), Fp::<P>::one());
            } else {
                assert!(a.inv().is_none());
            }
        }
    }

    #[test]
    fn field_axioms_f5_f7() {
        check_field_axioms::<5>();
        check_field_axioms::<7>();
        check_field_axioms::<13>();
    }

    #[test]
    fn inverse_matches_brute_force() {
        for a in elems::<11>() {
            let brute = elems::<11>()
                .into_iter()
                .find(|b| a.mul(b) == Fp::<11>::one());
            assert_eq!(a.inv(), brute);
        }
    }

    #[test]
    fn negation_is_genuine() {
        // unlike nimbers, neg is a real negation: −1 = P−1 ≠ 1 for odd P.
        let one = Fp::<5>::one();
        assert_eq!(one.neg(), Fp::<5>::from_u128(4));
        assert_ne!(one.neg(), one);
        assert_eq!(Fp::<5>::new(-1), Fp::<5>::from_u128(4));
        assert_eq!(Fp::<5>::characteristic(), 5);
    }

    #[test]
    fn clifford_over_f3_monomorphises() {
        // Cl over F_3 with q = [1, 2]: real antisymmetry (−1 = 2), and
        // (e0e1)² = −(q0 q1) = −2 = 1 (mod 3).
        let alg = CliffordAlgebra::new(
            2,
            Metric::diagonal(vec![Fp::<3>::from_u128(1), Fp::<3>::from_u128(2)]),
        );
        let (e0, e1) = (alg.gen(0), alg.gen(1));
        assert_eq!(alg.mul(&e0, &e0), alg.scalar(Fp::<3>::from_u128(1)));
        assert_eq!(alg.mul(&e1, &e1), alg.scalar(Fp::<3>::from_u128(2)));
        // e0 e1 = −(e1 e0), and −1 = 2 in F_3
        assert_eq!(
            alg.mul(&e0, &e1),
            alg.scalar_mul(&Fp::<3>::new(-1), &alg.mul(&e1, &e0))
        );
        // (e0e1)² = 1
        let e0e1 = alg.mul(&e0, &e1);
        assert_eq!(alg.mul(&e0e1, &e0e1), alg.scalar(Fp::<3>::from_u128(1)));
    }

    #[test]
    fn composite_modulus_is_rejected() {
        assert!(std::panic::catch_unwind(Fp::<4>::one).is_err());
        assert!(std::panic::catch_unwind(|| Fp::<9>::new(2)).is_err());
    }
}
