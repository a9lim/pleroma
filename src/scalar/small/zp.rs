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
//! (`Scalar::zero()/one()` take no `self`): `Zp<const P: u128, const K: u128>` is
//! `Z/p^k`, carried as the residue in `[0, p^k)`. Its characteristic is the
//! modulus `p^k` (the additive order of `1`), even though it is a truncation of
//! the characteristic-0 ring `Z_p` and not a field of characteristic `p`.

use crate::scalar::{mod_inverse_u128, Fp, Scalar};
use std::fmt;

/// An element of `Z/p^k` (the `p`-adic integers to precision `k`): the residue in
/// `[0, p^k)`.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct Zp<const P: u128, const K: u128>(pub u128);

impl<const P: u128, const K: u128> Zp<P, K> {
    pub fn assert_supported_ring() {
        assert!(
            Fp::<P>::modulus_is_prime() && K > 0,
            "Zp<P,K> needs prime P and positive precision K, got P={P}, K={K}"
        );
    }

    /// The modulus `p^k`.
    pub fn modulus() -> u128 {
        Self::assert_supported_ring();
        let mut acc = 1u128;
        for _ in 0..K {
            acc = acc.checked_mul(P).expect("Zp modulus exceeds u128");
        }
        acc
    }

    /// Reduce an integer (possibly negative) into `Z/p^k`.
    pub fn new(n: i128) -> Self {
        Self::assert_supported_ring();
        let m = Self::modulus() as i128;
        Zp((((n % m) + m) % m) as u128)
    }

    /// The `p`-adic valuation of this element, capped at the precision `k`
    /// (`v_p(0)` reads as `k`, the precision floor).
    pub fn valuation(&self) -> u128 {
        Self::assert_supported_ring();
        if self.0 == 0 {
            return K;
        }
        let mut n = self.0;
        let mut v = 0;
        while n.is_multiple_of(P) {
            n /= P;
            v += 1;
        }
        v
    }

    /// Whether this element is a unit (invertible) in `Z/p^k`: iff `p ∤ a`.
    pub fn is_unit(&self) -> bool {
        Self::assert_supported_ring();
        !self.0.is_multiple_of(P)
    }
}

impl<const P: u128, const K: u128> fmt::Debug for Zp<P, K> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} (mod {}^{})", self.0, P, K)
    }
}

impl<const P: u128, const K: u128> Scalar for Zp<P, K> {
    fn zero() -> Self {
        Self::assert_supported_ring();
        Zp(0)
    }

    fn one() -> Self {
        Self::assert_supported_ring();
        Zp(1 % Self::modulus())
    }

    fn add(&self, rhs: &Self) -> Self {
        Self::assert_supported_ring();
        let m = Self::modulus();
        Zp(self.0.checked_add(rhs.0).expect("Zp addition exceeds u128") % m)
    }

    fn neg(&self) -> Self {
        Self::assert_supported_ring();
        if self.0 == 0 {
            Zp(0)
        } else {
            Zp(Self::modulus() - self.0)
        }
    }

    fn mul(&self, rhs: &Self) -> Self {
        Self::assert_supported_ring();
        let m = Self::modulus();
        Zp(self
            .0
            .checked_mul(rhs.0)
            .expect("Zp multiplication exceeds u128")
            % m)
    }

    fn characteristic() -> u128 {
        Self::assert_supported_ring();
        // The finite quotient Z/p^k has characteristic p^k: p^k · 1 = 0, and no
        // smaller positive multiple of 1 vanishes.
        Self::modulus()
    }

    fn inv(&self) -> Option<Self> {
        Self::assert_supported_ring();
        // Local ring: a unit iff p ∤ a. Invert units by extended Euclid mod p^k;
        // return None for non-units (p | a, including 0) — the Omnific discipline,
        // never leaving the ring with a spurious 1/p.
        if !self.is_unit() {
            return None;
        }
        Some(Zp(mod_inverse_u128(self.0, Self::modulus())?))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::clifford::{CliffordAlgebra, Metric};

    fn elems<const P: u128, const K: u128>() -> Vec<Zp<P, K>> {
        (0..Zp::<P, K>::modulus()).map(Zp::<P, K>).collect()
    }

    fn check_ring_axioms<const P: u128, const K: u128>() {
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
    fn characteristic_is_the_modulus_not_the_prime() {
        // It is not the char-p field F_p, but the finite quotient still has a
        // literal ring characteristic: the additive order of 1.
        assert_eq!(Zp::<2, 3>::characteristic(), 8);
        assert_eq!(Zp::<3, 3>::characteristic(), 27);
    }

    #[test]
    fn invalid_parameters_are_rejected() {
        assert!(std::panic::catch_unwind(Zp::<4, 3>::one).is_err());
        assert!(std::panic::catch_unwind(Zp::<5, 0>::one).is_err());
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
