//! The unramified extension `Q_q = Frac(W_N(F_q))` of `Q_p` — the **field of
//! fractions** of the Witt vectors, truncated to precision `N`.
//!
//! This completes the (field, ring of integers) pairing on the unramified leg.
//! Where [`WittVec`] is the *ring* `W_N(F_q)` (with `p` a
//! non-unit, residue field `F_q = F_{p^F}`), this is the *field* `Q_q`: the
//! unique unramified extension of `Q_p` of residue degree `F`. It is to `WittVec`
//! exactly what [`Qp`](crate::scalar::Qp) is to [`Zp`](crate::scalar::Zp) —
//! adjoin an explicit **valuation** so `1/p` exists — and `Q_q` for `F = 1` *is*
//! `Q_p`.
//!
//! Representation: a nonzero `x = p^{val} · unit`, with `unit` a **Witt unit**
//! (residue `≠ 0` in `F_q`) carried in `W_N(F_q)`, and `val` a *signed* integer.
//! Zero is the sentinel `{ unit: 0, val: 0 }`. Multiplication uses the genuine
//! unramified-ring product of `WittVec` (the residue of a product of units is the
//! product of residues, still nonzero — so units stay units with no carry); the
//! valuations add.
//!
//! Characteristic is **0** (a genuine char-0 field), distinguishing it from
//! `WittVec`, whose `characteristic()` is the precision modulus `p^N`. A Clifford
//! algebra over `Q_q` is therefore semisimple.
//!
//! ## Precision contract (capped-relative — read this)
//!
//! Like [`Qp`](crate::scalar::Qp) (and for the same reason — `1/(p+1)` has
//! infinite support) this is the standard **capped-relative** model: `N`
//! significant `p`-adic digits relative to the valuation. Multiplication and
//! inversion are **exact** (valuations add; the mantissa is a genuine unit of
//! `W_N(F_q)`); **addition is not associative across precision boundaries** —
//! cancellation below the retained window reads as `0`. `Qq` is therefore a
//! *precision model*, not an exact ring, and is **excluded from the exact-ring
//! fuzz suite**.

use crate::scalar::{Fpn, Scalar, WittVec};
use std::fmt;

/// An element of `Q_q = Frac(W_N(F_q))`: `p^{val} · unit` with `unit` a Witt unit
/// (residue `≠ 0`), or the sentinel `{ 0, 0 }` for the field zero.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct Qq<const P: u128, const N: usize, const F: usize> {
    /// The unit mantissa in `W_N(F_q)` (residue `≠ 0`, except the field zero).
    unit: WittVec<P, N, F>,
    /// The (signed) `p`-adic valuation (`0` for the field zero).
    val: i128,
}

impl<const P: u128, const N: usize, const F: usize> Qq<P, N, F> {
    /// Canonicalize `p^{val} · w`: peel every factor of `p` from `w` into the
    /// valuation (capped-relative — the retained window slides up), landing on a
    /// Witt unit or the field zero.
    fn normalized(mut unit: WittVec<P, N, F>, mut val: i128) -> Self {
        loop {
            if unit.is_zero() {
                return Qq {
                    unit: WittVec::zero(),
                    val: 0,
                };
            }
            match unit.try_divide_by_p() {
                Some(d) => {
                    unit = d;
                    val += 1;
                }
                None => return Qq { unit, val }, // residue ≠ 0 ⇒ a unit
            }
        }
    }

    /// Embed a (signed) integer, extracting its `p`-adic valuation.
    pub fn from_int(n: i128) -> Self {
        Self::normalized(WittVec::from_int(n), 0)
    }

    /// Embed a Witt vector (a `W_N(F_q)` element) into its field of fractions.
    pub fn from_witt(w: WittVec<P, N, F>) -> Self {
        Self::normalized(w, 0)
    }

    /// `p^v` — the pure power, unit mantissa `1`. `from_p_power(-1)` is `1/p`.
    pub fn from_p_power(v: i128) -> Self {
        Qq {
            unit: WittVec::one(),
            val: v,
        }
    }

    /// The `p`-adic valuation, or `None` for zero (whose valuation is `+∞`).
    pub fn valuation(&self) -> Option<i128> {
        if self.unit.is_zero() {
            None
        } else {
            Some(self.val)
        }
    }

    /// The residue of the unit mantissa in `F_q` (the residue square-class carrier),
    /// or `None` for zero.
    pub fn unit_residue(&self) -> Option<Fpn<P, F>> {
        if self.unit.is_zero() {
            None
        } else {
            Some(self.unit.residue())
        }
    }

    /// The unit mantissa (a Witt unit, or the zero vector for the field zero).
    pub fn unit(&self) -> WittVec<P, N, F> {
        self.unit
    }
}

impl<const P: u128, const N: usize, const F: usize> fmt::Debug for Qq<P, N, F> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.unit.is_zero() {
            return write!(f, "0 (Q_{}^{})", P, F);
        }
        if self.val == 0 {
            write!(f, "{:?}", self.unit)
        } else {
            write!(f, "{:?}·{}^{}", self.unit, P, self.val)
        }
    }
}

impl<const P: u128, const N: usize, const F: usize> Scalar for Qq<P, N, F> {
    fn zero() -> Self {
        Qq {
            unit: WittVec::zero(),
            val: 0,
        }
    }

    fn one() -> Self {
        Qq {
            unit: WittVec::one(),
            val: 0,
        }
    }

    fn add(&self, rhs: &Self) -> Self {
        if self.unit.is_zero() {
            return *rhs;
        }
        if rhs.unit.is_zero() {
            return *self;
        }
        // Align on the smaller valuation: x + y = p^{vlo}·(ulo + p^{d}·uhi). In
        // W_N(F_q) the factor p^d is automatically 0 once d ≥ N (p^N = 0), so the
        // higher-valuation term vanishes below precision exactly like Qp.
        let (lo, hi) = if self.val <= rhs.val {
            (self, rhs)
        } else {
            (rhs, self)
        };
        let d = (hi.val - lo.val) as usize;
        let shifted = if d >= N {
            WittVec::zero()
        } else {
            let p = WittVec::<P, N, F>::from_int(P as i128);
            let mut s = hi.unit;
            for _ in 0..d {
                s = s.mul(&p);
            }
            s
        };
        Self::normalized(lo.unit.add(&shifted), lo.val)
    }

    fn neg(&self) -> Self {
        Qq {
            unit: self.unit.neg(),
            val: self.val,
        }
    }

    fn mul(&self, rhs: &Self) -> Self {
        if self.unit.is_zero() || rhs.unit.is_zero() {
            return Self::zero();
        }
        // The product of two Witt units is a Witt unit (residue = product of
        // residues ≠ 0), so no renormalization is needed; valuations add.
        Qq {
            unit: self.unit.mul(&rhs.unit),
            val: self.val + rhs.val,
        }
    }

    fn characteristic() -> u128 {
        0 // a genuine char-0 field, unlike WittVec's precision modulus p^N
    }

    fn inv(&self) -> Option<Self> {
        // Total on nonzero: (p^v·u)^{-1} = p^{-v}·u^{-1}, and the Witt unit u
        // always inverts in W_N(F_q) (residue ≠ 0). THE field property.
        if self.unit.is_zero() {
            return None;
        }
        Some(Qq {
            unit: self.unit.inv().expect("a Witt unit must invert"),
            val: -self.val,
        })
    }

    fn is_zero(&self) -> bool {
        self.unit.is_zero()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scalar::Qp;

    type Q2 = Qq<2, 4, 1>; // F = 1 ⇒ Q_2 itself, 4 digits
    type Q4 = Qq<2, 4, 2>; // unramified quadratic ext of Q_2 (residue F_4)
    type Q9 = Qq<3, 3, 2>; // unramified quadratic ext of Q_3 (residue F_9)

    #[test]
    fn reduces_to_qp_when_residue_degree_is_one() {
        // ORACLE: Q_q with F = 1 is Q_p — match the Qp backend digit-for-digit.
        for u in (1..16u128).filter(|u| u % 2 != 0) {
            for v in -2i128..=2 {
                let a = Qq::<2, 4, 1>::normalized(WittVec::from_int(u as i128), v);
                let b = Qp::<2, 4>::from_i128(u as i128).mul(&Qp::<2, 4>::from_p_power(v));
                // compare via the field invariants both expose
                assert_eq!(a.valuation(), b.valuation());
                let ai = a.inv().unwrap();
                assert_eq!(ai.valuation(), b.inv().unwrap().valuation());
                assert_eq!(a.mul(&ai), Q2::one());
            }
        }
    }

    #[test]
    fn one_over_p_exists_and_is_a_field() {
        let p = Q4::from_int(2);
        let pinv = p.inv().unwrap();
        assert_eq!(pinv, Q4::from_p_power(-1));
        assert_eq!(p.mul(&pinv), Q4::one());
        assert_eq!(pinv.valuation(), Some(-1));
        assert_eq!(Q4::zero().inv(), None);
    }

    #[test]
    fn valuation_is_additive_under_multiplication() {
        let a = Q4::from_int(12); // 2^2 · 3 ⇒ val 2
        let b = Q4::from_int(20); // 2^2 · 5 ⇒ val 2
        assert_eq!(a.valuation(), Some(2));
        assert_eq!(b.valuation(), Some(2));
        assert_eq!(a.mul(&b).valuation(), Some(4));
        assert_eq!(a.mul(&Q4::from_p_power(-1)).valuation(), Some(1));
    }

    #[test]
    fn residue_degree_two_has_genuine_fq_residues() {
        // A unit whose residue is the F_4 generator (not in the prime field):
        // its inverse's residue is the F_4 inverse — genuinely unramified content.
        let g = WittVec::<2, 4, 2>([0, 1]); // residue = t ∈ F_4 (Fpn generator)
        let x = Q4::from_witt(g);
        assert_eq!(x.valuation(), Some(0));
        let xi = x.inv().unwrap();
        assert_eq!(x.mul(&xi), Q4::one());
        // residue of the inverse = (residue)^{-1} in F_4
        assert_eq!(
            xi.unit_residue().unwrap(),
            g.residue().inv().unwrap(),
            "residue map commutes with inversion"
        );
    }

    #[test]
    fn characteristic_is_zero_not_the_modulus() {
        assert_eq!(Q4::characteristic(), 0);
        assert_eq!(Q9::characteristic(), 0);
        // contrast: the ring of integers has characteristic p^N
        assert_eq!(WittVec::<2, 4, 2>::characteristic(), 16);
    }

    #[test]
    fn multiplicative_group_is_exact() {
        // Multiplication is exact (the relative-precision caveat is additive only):
        // a small unramified multiplicative group over Q_9 is a genuine abelian group.
        let units: Vec<Q9> = [[1u128, 0], [0, 1], [1, 1], [2, 1], [1, 2]]
            .iter()
            .flat_map(|c| (-1i128..=1).map(move |v| Qq::<3, 3, 2>::normalized(WittVec(*c), v)))
            .collect();
        let one = Q9::one();
        for a in &units {
            assert_eq!(a.mul(&one), *a);
            assert_eq!(a.mul(&a.inv().unwrap()), one);
            for b in &units {
                assert_eq!(a.mul(b), b.mul(a));
                for c in &units {
                    assert_eq!(a.mul(b).mul(c), a.mul(&b.mul(c)));
                }
            }
        }
    }
}
