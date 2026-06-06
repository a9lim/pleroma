//! The omnific integers `Oz` — Conway's integers of the surreal world.
//!
//! `Oz` is to the surreals as ℤ is to ℝ: the ring of surreals with no
//! fractional part. Concretely, a surreal in Conway normal form
//! `Σ ω^{y_i} r_i` is an **omnific integer** iff
//!
//!   * it has **no infinitesimal terms** (every exponent `y_i ≥ 0`), and
//!   * its **constant term** (the `y = 0` coefficient, if present) is an
//!     ordinary integer.
//!
//! So `ω`, `ω²+3`, `½ω`, `ω^ω − ω + 7` are omnific integers; `ε = ω⁻¹`, `ω + ½`,
//! `5/3` are not. (The coefficients on *infinite* terms may be any rational —
//! only the finite part is constrained, exactly as `⌊x⌋` is unconstrained above
//! the units digit.)
//!
//! This is the surreal mirror of the `Integer` backend: a **transfinite
//! commutative ring**, not a field. A Clifford algebra needs only a commutative
//! ring of scalars, so `Oz` supports the Clifford-with-nilpotents / exterior
//! structure (the headline being the exterior algebra with *transfinite*
//! coefficients) — while the clean Cl(p,q) classification, which needs a field,
//! does not apply. Only `±1` are units, so `inv` returns `Some` only there
//! (delegating to `Surreal::inv` would leave the ring: `ω ↦ ε`).

use crate::scalar::Scalar;
use crate::scalar::Surreal;

/// An omnific integer: a surreal with no infinitesimal part and an integer
/// constant term. The inner surreal is private so every value is validated at
/// construction; arithmetic preserves the property automatically.
#[derive(Clone, PartialEq)]
pub struct Omnific(Surreal);

impl std::fmt::Debug for Omnific {
    // delegate to the inner surreal so multivector displays read `ω·e0e1`, not
    // `Omnific(ω)·e0e1`.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self.0)
    }
}

/// True iff the surreal `s` is an omnific integer (see the module docs).
pub fn is_omnific_integer(s: &Surreal) -> bool {
    for (exp, coeff) in s.terms() {
        use std::cmp::Ordering::*;
        match exp.sign() {
            Less => return false, // an infinitesimal term ⇒ not an integer
            Equal => {
                // the constant term must be an integer coefficient
                if !coeff.is_integer() {
                    return false;
                }
            }
            Greater => {} // infinite term: any rational coefficient is fine
        }
    }
    true
}

impl Omnific {
    /// The omnific integer `n`.
    pub fn from_int(n: i128) -> Self {
        Omnific(Surreal::from_int(n))
    }

    /// `ω`, the simplest infinite omnific integer.
    pub fn omega() -> Self {
        Omnific(Surreal::omega())
    }

    /// Wrap a surreal, returning `None` unless it is an omnific integer.
    pub fn from_surreal(s: Surreal) -> Option<Self> {
        if is_omnific_integer(&s) {
            Some(Omnific(s))
        } else {
            None
        }
    }

    /// The underlying surreal.
    pub fn inner(&self) -> &Surreal {
        &self.0
    }
}

impl Scalar for Omnific {
    fn zero() -> Self {
        Omnific(Surreal::zero())
    }
    fn one() -> Self {
        Omnific(Surreal::one())
    }
    fn add(&self, rhs: &Self) -> Self {
        // sums of omnific integers are omnific integers — no re-validation needed
        Omnific(self.0.add(&rhs.0))
    }
    fn neg(&self) -> Self {
        Omnific(self.0.neg())
    }
    fn mul(&self, rhs: &Self) -> Self {
        Omnific(self.0.mul(&rhs.0))
    }
    fn characteristic() -> u32 {
        0
    }
    fn inv(&self) -> Option<Self> {
        // a ring, not a field: only ±1 are units. (Do NOT delegate to
        // Surreal::inv, which would send ω ↦ ε and leave the ring.)
        if self.0 == Surreal::one() || self.0 == Surreal::one().neg() {
            Some(self.clone())
        } else {
            None
        }
    }
    fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::clifford::{CliffordAlgebra, Metric};
    use crate::scalar::{Integer, Rational};

    fn surr_int(n: i128) -> Surreal {
        Surreal::from_int(n)
    }

    #[test]
    fn validator_accepts_omnific_integers() {
        assert!(is_omnific_integer(&Surreal::omega())); // ω
        assert!(is_omnific_integer(
            &Surreal::omega_pow(surr_int(2)).add(&surr_int(3))
        )); // ω²+3
            // ½ω: an infinite term with a fractional coefficient is fine
        assert!(is_omnific_integer(&Surreal::monomial(
            surr_int(1),
            Rational::new(1, 2)
        )));
        assert!(is_omnific_integer(&surr_int(7)));
        assert!(is_omnific_integer(&Surreal::zero()));
    }

    #[test]
    fn validator_rejects_non_integers() {
        assert!(!is_omnific_integer(&Surreal::epsilon())); // ε = ω⁻¹
                                                           // ω + ½: fractional constant term
        assert!(!is_omnific_integer(
            &Surreal::omega().add(&Surreal::from_rational(Rational::new(1, 2)))
        ));
        // 5/3: fractional constant
        assert!(!is_omnific_integer(&Surreal::from_rational(Rational::new(
            5, 3
        ))));
        assert!(Omnific::from_surreal(Surreal::epsilon()).is_none());
    }

    #[test]
    fn ring_axioms_on_a_sample() {
        let a = Omnific::omega(); // ω
        let b = Omnific::from_int(3);
        let c = Omnific::from_surreal(Surreal::omega_pow(surr_int(2))).unwrap(); // ω²
                                                                                 // associativity + distributivity
        assert_eq!(a.add(&b).add(&c), a.add(&b.add(&c)));
        assert_eq!(a.mul(&b).mul(&c), a.mul(&b.mul(&c)));
        assert_eq!(a.mul(&b.add(&c)), a.mul(&b).add(&a.mul(&c)));
        // commutativity
        assert_eq!(a.add(&b), b.add(&a));
        assert_eq!(a.mul(&b), b.mul(&a));
        // identity + inverse-of-add
        assert_eq!(a.add(&Omnific::zero()), a);
        assert_eq!(a.mul(&Omnific::one()), a);
        assert_eq!(a.sub(&a), Omnific::zero());
        // closure: every result is still an omnific integer
        assert!(is_omnific_integer(a.mul(&b).add(&c).inner()));
    }

    #[test]
    fn only_plus_minus_one_are_units() {
        assert!(Omnific::one().inv().is_some());
        assert!(Omnific::one().neg().inv().is_some());
        assert!(Omnific::from_int(2).inv().is_none());
        assert!(Omnific::omega().inv().is_none()); // ω is not a unit (1/ω = ε ∉ Oz)
        assert!(Omnific::zero().inv().is_none());
    }

    #[test]
    fn exterior_algebra_over_oz_with_transfinite_coefficients() {
        // Λ over the transfinite ring Oz: nilpotent generators, antisymmetry,
        // and ω-scale coefficients flow through the wedge.
        let alg = CliffordAlgebra::new(3, Metric::<Omnific>::grassmann(3));
        let (e0, e1, e2) = (alg.gen(0), alg.gen(1), alg.gen(2));
        // e_i² = 0
        assert!(alg.mul(&e0, &e0).is_zero());
        // antisymmetry: e0 e1 + e1 e0 = 0 (b = 0)
        assert_eq!(alg.add(&alg.mul(&e0, &e1), &alg.mul(&e1, &e0)), alg.zero());
        // a transfinite coefficient survives: (ω·e0) ∧ e1 = ω·e0e1
        let we0 = alg.scalar_mul(&Omnific::omega(), &e0);
        let prod = alg.wedge(&we0, &e1);
        assert_eq!(
            prod,
            alg.scalar_mul(&Omnific::omega(), &alg.wedge(&e0, &e1))
        );
        // and a triple wedge is the top blade
        let triple = alg.wedge(&alg.wedge(&e0, &e1), &e2);
        assert_eq!(triple, alg.blade(&[0, 1, 2]));
    }

    #[test]
    fn matches_integer_backend_on_integer_inputs() {
        // The exterior structure over Oz agrees with the ℤ backend when the
        // coefficients are ordinary integers: 2 e0 ∧ 3 e1 = 6 e0e1.
        let oz = CliffordAlgebra::new(2, Metric::<Omnific>::grassmann(2));
        let oz_prod = oz.wedge(
            &oz.scalar_mul(&Omnific::from_int(2), &oz.gen(0)),
            &oz.scalar_mul(&Omnific::from_int(3), &oz.gen(1)),
        );
        assert_eq!(
            oz_prod,
            oz.scalar_mul(&Omnific::from_int(6), &oz.wedge(&oz.gen(0), &oz.gen(1)))
        );

        let zz = CliffordAlgebra::new(2, Metric::<Integer>::grassmann(2));
        let zz_prod = zz.wedge(
            &zz.scalar_mul(&Integer(2), &zz.gen(0)),
            &zz.scalar_mul(&Integer(3), &zz.gen(1)),
        );
        // same blade, coefficient 6 in each backend's own ring
        assert_eq!(*zz_prod.terms.get(&0b11).unwrap(), Integer(6));
        assert_eq!(*oz_prod.terms.get(&0b11).unwrap(), Omnific::from_int(6));
    }
}
