//! The `p`-adic numbers `Q_p`, truncated to precision `k` — the **field** of
//! fractions of [`Zp`](crate::scalar::Zp).
//!
//! Where [`Zp`](crate::scalar::Zp) is the *ring* `Z/p^k` (with `p` a non-unit),
//! this is the *field* `Q_p`: every nonzero element is invertible because we
//! carry an explicit **valuation** alongside a unit mantissa, so `1/p` exists
//! (it is `p^{-1}`). This is the exact p-adic mirror of the
//! [`Omnific`](crate::scalar::Omnific) ⊂ [`Surreal`](crate::scalar::Surreal)
//! relationship on the ω-side: a ring sitting inside its field of fractions,
//! the field obtained by allowing *negative* leading exponents.
//!
//! Representation: a nonzero `x = p^{val} · unit`, with `unit` a p-adic **unit**
//! (`p ∤ unit`) carried mod `p^k` (so `k` significant p-adic digits), and `val`
//! a *signed* integer. Zero is the sentinel `{ unit: 0, val: 0 }`.
//!
//! Characteristic is **0** (it is a genuine char-0 field), distinguishing it
//! from `Zp`, whose `characteristic()` is the modulus `p^k`. A Clifford algebra
//! over `Q_p` is therefore semisimple — the companion
//! [`forms::padic`](crate::forms::padic) / `springer::padic` modules read their
//! Hilbert-symbol / residue-form payoff off this backend.
//!
//! ## Precision contract (capped-relative — read this)
//!
//! There is no finite-memory *exact* model of `Q_p`: a unit can have an inverse
//! of infinite p-adic support (`1/(p+1) = 1 − p + p² − …`), so any concrete
//! backend truncates. This one uses the standard **capped-relative** model
//! (as in SageMath's default p-adics): `k` significant mantissa digits relative
//! to the valuation. Multiplication and inversion are **exact** (valuations add;
//! the mantissa is a genuine unit of `Z/p^k`); **addition is not associative
//! across precision boundaries** — additive cancellation below the retained
//! window reads as `0`, exactly like floating point. `Qp` is therefore a
//! *precision model*, not an exact commutative ring: it is used at the forms
//! layer (valuation + residue square class, both robust to relative precision)
//! and is deliberately **excluded from the exact-ring fuzz suite**.

use crate::scalar::{mod_inverse_u128, Fp, Rational, Scalar};
use std::fmt;

/// An element of `Q_p` to precision `k`: `p^{val} · unit` with `p ∤ unit` carried
/// mod `p^k`, or the sentinel `{ 0, 0 }` for the field zero.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct Qp<const P: u128, const K: u128> {
    /// The unit mantissa in `[0, p^k)` (`0` only for the field zero; otherwise
    /// `p ∤ unit`).
    unit: u128,
    /// The (signed) p-adic valuation of the element (`0` for the field zero).
    val: i128,
}

/// `P^e`, checked against `u128` overflow.
fn p_pow<const P: u128>(e: u128) -> u128 {
    let mut acc = 1u128;
    for _ in 0..e {
        acc = acc.checked_mul(P).expect("Qp: p-power exceeds u128");
    }
    acc
}

impl<const P: u128, const K: u128> Qp<P, K> {
    pub fn assert_supported_field() {
        assert!(
            Fp::<P>::modulus_is_prime() && K > 0,
            "Qp<P,K> needs prime P and positive precision K, got P={P}, K={K}"
        );
        let mut acc = 1u128;
        for _ in 0..K {
            acc = acc.checked_mul(P).expect("Qp modulus exceeds u128");
            assert!(
                acc <= i128::MAX as u128,
                "Qp<P,K> modulus must fit i128-backed embeddings, got P={P}, K={K}"
            );
        }
    }

    /// The mantissa modulus `p^k`.
    pub fn modulus() -> u128 {
        Self::assert_supported_field();
        p_pow::<P>(K)
    }

    /// Build `p^{val} · unit`, normalizing: factor any `p` out of `unit` into the
    /// valuation, reduce the mantissa mod `p^k`. A `unit` reducing to `0 mod p^k`
    /// yields the field zero.
    fn normalized(unit_raw: u128, val: i128) -> Self {
        let m = Self::modulus();
        let mut u = unit_raw % m;
        if u == 0 {
            return Qp { unit: 0, val: 0 };
        }
        let mut v = val;
        while u.is_multiple_of(P) {
            u /= P;
            v += 1;
        }
        Qp { unit: u, val: v }
    }

    /// Embed a (signed) integer, extracting its p-adic valuation.
    pub fn from_i128(n: i128) -> Self {
        Self::assert_supported_field();
        if n == 0 {
            return Qp { unit: 0, val: 0 };
        }
        let pp = P as i128;
        let mut w = 0i128;
        let mut nn = n;
        while nn % pp == 0 {
            nn /= pp;
            w += 1;
        }
        let m = Self::modulus() as i128;
        let unit = (((nn % m) + m) % m) as u128;
        Qp { unit, val: w }
    }

    /// `p^v` — the pure power, mantissa `1`. `from_p_power(-1)` is `1/p`.
    pub fn from_p_power(v: i128) -> Self {
        Self::assert_supported_field();
        Qp {
            unit: 1 % Self::modulus(),
            val: v,
        }
    }

    /// Embed a rational number into `Q_p`: `from_i128(num) · from_i128(den)^{-1}`.
    pub fn from_rational(q: &Rational) -> Self {
        let num = Self::from_i128(q.numer());
        let den = Self::from_i128(q.denom());
        num.mul(&den.inv().expect("Qp::from_rational: nonzero denominator"))
    }

    /// The p-adic valuation, or `None` for zero (whose valuation is `+∞`).
    pub fn valuation(&self) -> Option<i128> {
        Self::assert_supported_field();
        if self.unit == 0 {
            None
        } else {
            Some(self.val)
        }
    }

    /// The unit mantissa in `[0, p^k)`.
    pub fn unit(&self) -> u128 {
        Self::assert_supported_field();
        self.unit
    }
}

impl<const P: u128, const K: u128> fmt::Debug for Qp<P, K> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.unit == 0 {
            return write!(f, "0 (Q_{})", P);
        }
        if self.val == 0 {
            write!(f, "{} (mod {}^{})", self.unit, P, K)
        } else {
            write!(f, "{}·{}^{} (mod {}^{})", self.unit, P, self.val, P, K)
        }
    }
}

impl<const P: u128, const K: u128> Scalar for Qp<P, K> {
    fn zero() -> Self {
        Self::assert_supported_field();
        Qp { unit: 0, val: 0 }
    }

    fn one() -> Self {
        Self::assert_supported_field();
        Qp {
            unit: 1 % Self::modulus(),
            val: 0,
        }
    }

    fn add(&self, rhs: &Self) -> Self {
        Self::assert_supported_field();
        if self.unit == 0 {
            return *rhs;
        }
        if rhs.unit == 0 {
            return *self;
        }
        let m = Self::modulus();
        // Align on the smaller valuation: x + y = p^{vlo}·(ulo + p^{d}·uhi).
        let (lo, hi) = if self.val <= rhs.val {
            (self, rhs)
        } else {
            (rhs, self)
        };
        let d = (hi.val - lo.val) as u128;
        let shifted = if d >= K {
            0 // below precision — the higher-valuation term vanishes
        } else {
            p_pow::<P>(d)
                .checked_mul(hi.unit)
                .expect("Qp addition mantissa product exceeds u128")
                % m
        };
        let b = lo
            .unit
            .checked_add(shifted)
            .expect("Qp addition mantissa sum exceeds u128")
            % m;
        if b == 0 {
            return Qp { unit: 0, val: 0 }; // cancelled below precision
        }
        Self::normalized(b, lo.val)
    }

    fn neg(&self) -> Self {
        Self::assert_supported_field();
        if self.unit == 0 {
            return *self;
        }
        Qp {
            unit: Self::modulus() - self.unit,
            val: self.val,
        }
    }

    fn mul(&self, rhs: &Self) -> Self {
        Self::assert_supported_field();
        if self.unit == 0 || rhs.unit == 0 {
            return Qp { unit: 0, val: 0 };
        }
        // Product of units is a unit: no renormalization needed.
        let m = Self::modulus();
        Qp {
            unit: self
                .unit
                .checked_mul(rhs.unit)
                .expect("Qp multiplication mantissa product exceeds u128")
                % m,
            val: self
                .val
                .checked_add(rhs.val)
                .expect("Qp multiplication valuation exceeds i128"),
        }
    }

    fn characteristic() -> u128 {
        Self::assert_supported_field();
        0 // a genuine field of characteristic 0 — unlike Zp's modulus p^k
    }

    fn inv(&self) -> Option<Self> {
        Self::assert_supported_field();
        // Total on nonzero: (p^v·u)^{-1} = p^{-v}·u^{-1}. THE field property,
        // versus Zp::inv which is None for any p-divisible element.
        if self.unit == 0 {
            return None;
        }
        let uinv = mod_inverse_u128(self.unit, Self::modulus())?;
        Some(Qp {
            unit: uinv,
            val: -self.val,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    type Q5 = Qp<5, 4>; // Q_5 to 4 digits
    type Q2 = Qp<2, 6>; // Q_2 to 6 digits
    type Q3 = Qp<3, 3>; // Q_3 to 3 digits

    #[test]
    fn one_over_p_exists_and_is_a_field() {
        // The defining win over Zp: p is a unit here, 1/p = p^{-1}.
        let p = Q5::from_i128(5);
        let pinv = p.inv().unwrap();
        assert_eq!(pinv, Q5::from_p_power(-1));
        assert_eq!(p.mul(&pinv), Q5::one());
        // And 1/p genuinely has negative valuation.
        assert_eq!(pinv.valuation(), Some(-1));
        // Zero is the only non-invertible element.
        assert_eq!(Q5::zero().inv(), None);
    }

    #[test]
    fn from_rational_matches_integer_embedding_and_denominator_inverse() {
        let x = Q5::from_rational(&Rational::new(50, 3));
        let expected = Q5::from_i128(50).mul(&Q5::from_i128(3).inv().unwrap());
        assert_eq!(x, expected);
        assert_eq!(x.valuation(), Some(2));

        let y = Q5::from_rational(&Rational::new(3, 50));
        assert_eq!(y.valuation(), Some(-2));
        assert_eq!(x.mul(&y), Q5::one());
    }

    #[test]
    fn every_nonzero_inverts() {
        // Sample units × a spread of valuations; all invert (field), unlike Zp.
        for u in 1..Q3::modulus() {
            if u % 3 == 0 {
                continue; // not a unit mantissa
            }
            for v in -2i128..=2 {
                let x = Qp::<3, 3>::normalized(u, v);
                let xi = x.inv().expect("Q_p: every nonzero inverts");
                assert_eq!(x.mul(&xi), Q3::one(), "x·x⁻¹ ≠ 1 for {x:?}");
            }
        }
    }

    #[test]
    fn valuation_is_additive_under_multiplication() {
        let a = Q2::from_i128(12); // 2^2 · 3 ⇒ val 2
        let b = Q2::from_i128(20); // 2^2 · 5 ⇒ val 2
        assert_eq!(a.valuation(), Some(2));
        assert_eq!(b.valuation(), Some(2));
        assert_eq!(a.mul(&b).valuation(), Some(4));
        // 1/p drops the valuation.
        assert_eq!(a.mul(&Q2::from_p_power(-1)).valuation(), Some(1));
    }

    #[test]
    fn characteristic_is_zero_not_the_modulus() {
        assert_eq!(Q5::characteristic(), 0);
        assert_eq!(Q2::characteristic(), 0);
    }

    #[test]
    fn invalid_parameters_are_rejected() {
        assert!(std::panic::catch_unwind(Qp::<4, 3>::one).is_err());
        assert!(std::panic::catch_unwind(Qp::<5, 0>::one).is_err());
        assert!(std::panic::catch_unwind(Qp::<2, 127>::one).is_err());
    }

    #[test]
    fn multiplication_is_an_exact_abelian_group_on_nonzero() {
        // Multiplication is exact (valuations add, mantissa is a genuine unit of
        // Z/p^k): associative, commutative, with one and inverses — a real field
        // multiplicative group, unaffected by the relative-precision caveat.
        let es: Vec<Q3> = {
            let mut v = Vec::new();
            for u in (1..Q3::modulus()).filter(|u| u % 3 != 0) {
                for val in -2i128..=2 {
                    v.push(Qp::<3, 3>::normalized(u, val));
                }
            }
            v
        };
        let one = Q3::one();
        for a in &es {
            assert_eq!(a.mul(&one), *a);
            assert_eq!(a.mul(&a.inv().unwrap()), one);
            for b in &es {
                assert_eq!(a.mul(b), b.mul(a)); // commutative
                for c in &es {
                    assert_eq!(a.mul(b).mul(c), a.mul(&b.mul(c))); // associative
                }
            }
        }
    }

    #[test]
    fn addition_exact_facts_hold() {
        // The precision-safe additive facts (identity, negation, commutativity)
        // hold exactly for every element; only cross-precision associativity is
        // sacrificed (documented below).
        let es: Vec<Q3> = (0..Q3::modulus() as i128)
            .map(Qp::<3, 3>::from_i128)
            .collect();
        let zero = Q3::zero();
        for a in &es {
            assert_eq!(a.add(&zero), *a);
            assert_eq!(a.add(&a.neg()), zero);
            for b in &es {
                assert_eq!(a.add(b), b.add(a)); // commutative
            }
        }
    }

    #[test]
    fn relative_precision_cancellation_is_intended() {
        // INTENDED, not a bug: 1 + (p^k − 1) at valuation 0 cancels in the top
        // k digits, lands at valuation k = outside the retained window, and so
        // reads as 0 — the capped-relative contract (cf. floating point).
        let one = Q3::one();
        let almost = Qp::<3, 3>::from_i128(Q3::modulus() as i128 - 1); // p^k − 1, a unit
        assert_eq!(one.add(&almost), Q3::zero());
        // Whereas 1 + p (a precision-safe sum) is exact and nonzero.
        assert_eq!(one.add(&Qp::<3, 3>::from_i128(3)), Qp::<3, 3>::from_i128(4));
    }

    #[test]
    fn p_times_one_over_p_is_one_each_prime() {
        assert_eq!(Q2::from_i128(2).mul(&Q2::from_p_power(-1)), Q2::one());
        assert_eq!(Q3::from_i128(3).mul(&Q3::from_p_power(-1)), Q3::one());
        assert_eq!(Q5::from_i128(5).mul(&Q5::from_p_power(-1)), Q5::one());
    }
}
