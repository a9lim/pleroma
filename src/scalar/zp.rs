//! The `p`-adic integers `Z_p`, truncated to precision `k` — the ring `Z/p^k`.
//!
//! This is the **ring of integers** of the field `Q_p`, not `Q_p` itself: `p` is a
//! **non-unit** here (it generates the maximal ideal), so `Z_p` is a *local ring*,
//! not a field — exactly the [`Omnific`](crate::scalar::Omnific)/`Integer` posture,
//! and the reason it is named `Zp`, not `Qp`. A Clifford algebra needs only a
//! commutative *ring* of scalars, so `Z/p^k` supports the
//! Clifford-with-nilpotents / exterior structure; because `p` is a non-unit, a
//! Clifford algebra over `Z/p^k` is a genuine zero-divisor / non-semisimple object —
//! the engine's nilpotent path exercised at the *scalar* level.
//!
//! Where the companion [`forms::padic`](crate::forms::padic) module supplies the
//! quadratic-form payoff over `Q_p` (the Hilbert symbol, Hasse–Minkowski) at the
//! forms layer, this is "the integers underneath" as an actual `Scalar` backend.
//!
//! ## The const-generic modulus, two parameters
//!
//! Like `Fp`/`Fpn`, both the prime `p` and the precision `k` live in the **type**
//! (`Scalar::zero()/one()` take no `self`): `Zp<const P: u64, const K: u32>` is
//! `Z/p^k`, carried as the residue in `[0, p^k)`. `characteristic()` returns **0** —
//! `Z/p^k` is a length-`k` truncation of the characteristic-0 ring `Z_p`, and is not
//! a finite field of characteristic `p`; reporting 0 keeps it out of the
//! finite-field classifier, which would be wrong here.

use crate::scalar::Scalar;
use std::fmt;

/// An element of `Z/p^k` (the `p`-adic integers to precision `k`): the residue in
/// `[0, p^k)`.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct Zp<const P: u64, const K: u32>(pub u64);

impl<const P: u64, const K: u32> Zp<P, K> {
    /// The modulus `p^k`.
    pub fn modulus() -> u128 {
        (P as u128).pow(K)
    }

    /// Reduce an integer (possibly negative) into `Z/p^k`.
    pub fn new(n: i128) -> Self {
        let m = Self::modulus() as i128;
        Zp((((n % m) + m) % m) as u64)
    }

    /// The `p`-adic valuation of this element, capped at the precision `k`
    /// (`v_p(0)` reads as `k`, the precision floor).
    pub fn valuation(&self) -> u32 {
        if self.0 == 0 {
            return K;
        }
        let mut n = self.0;
        let mut v = 0;
        while n % P == 0 {
            n /= P;
            v += 1;
        }
        v
    }

    /// Whether this element is a unit (invertible) in `Z/p^k`: iff `p ∤ a`.
    pub fn is_unit(&self) -> bool {
        self.0 % P != 0
    }
}

impl<const P: u64, const K: u32> fmt::Debug for Zp<P, K> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} (mod {}^{})", self.0, P, K)
    }
}

impl<const P: u64, const K: u32> Scalar for Zp<P, K> {
    fn zero() -> Self {
        Zp(0)
    }

    fn one() -> Self {
        Zp((1 % Self::modulus()) as u64)
    }

    fn add(&self, rhs: &Self) -> Self {
        let m = Self::modulus();
        Zp(((self.0 as u128 + rhs.0 as u128) % m) as u64)
    }

    fn neg(&self) -> Self {
        if self.0 == 0 {
            Zp(0)
        } else {
            Zp((Self::modulus() - self.0 as u128) as u64)
        }
    }

    fn mul(&self, rhs: &Self) -> Self {
        let m = Self::modulus();
        Zp(((self.0 as u128 * rhs.0 as u128) % m) as u64)
    }

    fn characteristic() -> u32 {
        // Z/p^k is a truncation of the characteristic-0 ring Z_p, and is a local
        // *ring*, not a field — so 0, matching Integer/Omnific (NOT the prime p).
        0
    }

    fn inv(&self) -> Option<Self> {
        // Local ring: a unit iff p ∤ a. Invert units by extended Euclid mod p^k;
        // return None for non-units (p | a, including 0) — the Omnific discipline,
        // never leaving the ring with a spurious 1/p.
        if !self.is_unit() {
            return None;
        }
        let m = Self::modulus() as i128;
        let (mut t, mut newt) = (0i128, 1i128);
        let (mut r, mut newr) = (m, self.0 as i128);
        while newr != 0 {
            let quot = r / newr;
            t -= quot * newt;
            std::mem::swap(&mut t, &mut newt);
            r -= quot * newr;
            std::mem::swap(&mut r, &mut newr);
        }
        // r = gcd = 1 for a unit; t is the inverse.
        Some(Zp((((t % m) + m) % m) as u64))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::clifford::{CliffordAlgebra, Metric};

    fn elems<const P: u64, const K: u32>() -> Vec<Zp<P, K>> {
        (0..Zp::<P, K>::modulus() as u64).map(Zp::<P, K>).collect()
    }

    fn check_ring_axioms<const P: u64, const K: u32>() {
        let es = elems::<P, K>();
        let zero = Zp::<P, K>::zero();
        let one = Zp::<P, K>::one();
        for &a in &es {
            assert_eq!(a.add(&zero), a);
            assert_eq!(a.add(&a.neg()), zero);
            assert_eq!(a.mul(&one), a);
            // unit iff p ∤ a; units invert, non-units do not (ring, not field).
            if a.is_unit() {
                let ai = a.inv().expect("a unit must be invertible");
                assert_eq!(a.mul(&ai), one);
            } else {
                assert!(a.inv().is_none(), "a non-unit must not invert");
            }
            for &b in &es {
                assert_eq!(a.add(&b), b.add(&a));
                assert_eq!(a.mul(&b), b.mul(&a));
                for &c in &es {
                    assert_eq!(a.add(&b).add(&c), a.add(&b.add(&c)));
                    assert_eq!(a.mul(&b).mul(&c), a.mul(&b.mul(&c)));
                    assert_eq!(a.mul(&b.add(&c)), a.mul(&b).add(&a.mul(&c)));
                }
            }
        }
    }

    #[test]
    fn ring_axioms_z8_z9_z27_z16() {
        check_ring_axioms::<2, 3>(); // Z/8
        check_ring_axioms::<3, 2>(); // Z/9
        check_ring_axioms::<3, 3>(); // Z/27
        check_ring_axioms::<2, 4>(); // Z/16
    }

    #[test]
    fn p_is_a_non_unit_the_defining_property() {
        // p (and its multiples) are non-units — this is what makes it a ring, not a
        // field, and distinguishes Z_p from Q_p.
        assert_eq!(Zp::<2, 3>(2).inv(), None); // 2 ∤-invertible in Z/8
        assert_eq!(Zp::<2, 3>(4).inv(), None);
        assert_eq!(Zp::<2, 3>(0).inv(), None);
        assert_eq!(Zp::<3, 3>(3).inv(), None);
        // odd residues ARE units in Z/8.
        assert_eq!(Zp::<2, 3>(3).inv(), Some(Zp::<2, 3>(3))); // 3·3 = 9 ≡ 1 mod 8
        assert_eq!(Zp::<2, 3>(7).inv(), Some(Zp::<2, 3>(7))); // 7·7 = 49 ≡ 1 mod 8
    }

    #[test]
    fn characteristic_is_zero_not_p() {
        // It is the integers underneath, not a char-p field.
        assert_eq!(Zp::<2, 3>::characteristic(), 0);
        assert_eq!(Zp::<3, 3>::characteristic(), 0);
    }

    #[test]
    fn inverse_reduces_to_the_mod_p_inverse() {
        // Hensel consistency: the unit inverse in Z/p^k reduces mod p to the F_p
        // inverse. In Z/27, 2⁻¹ = 14 (2·14 = 28 ≡ 1); 14 ≡ 2 mod 3, and 2⁻¹ = 2 in F_3.
        let two = Zp::<3, 3>(2);
        let inv = two.inv().unwrap();
        assert_eq!(two.mul(&inv), Zp::<3, 3>::one());
        assert_eq!(inv.0 % 3, 2);
    }

    #[test]
    fn clifford_over_z4_runs_with_a_nonunit_scalar() {
        // Cl over Z/4 with q = [1, 2]: the engine runs over a non-field ring; 2 is a
        // zero divisor (2·2 = 0 in Z/4), so this is a genuinely non-semisimple
        // Clifford algebra — the nilpotent path exercised at the scalar level.
        let alg = CliffordAlgebra::new(2, Metric::diagonal(vec![Zp::<2, 2>(1), Zp::<2, 2>(2)]));
        let (e0, e1) = (alg.gen(0), alg.gen(1));
        assert_eq!(alg.mul(&e0, &e0), alg.scalar(Zp::<2, 2>(1)));
        assert_eq!(alg.mul(&e1, &e1), alg.scalar(Zp::<2, 2>(2)));
        // (e1)⁴ = q1² = 2² = 0 in Z/4: a nilpotent generator.
        let e1sq = alg.mul(&e1, &e1);
        assert_eq!(alg.mul(&e1sq, &e1sq), alg.scalar(Zp::<2, 2>::zero()));
    }
}
