//! The (field, ring of integers) pairing, made **structural**.
//!
//! The "any number" table in [`scalar`](crate::scalar) is organised around one
//! recurring relationship: almost every field ships beside its **ring of
//! integers** — `ℚ`/`ℤ`, `No`/`Oz`, `Q_p`/`Z_p`, `Q_q`/`W_N(F_q)`. Until now
//! that pairing lived only in doc comments. These two traits promote it to the
//! type system, so the relationship is checkable rather than merely described:
//!
//!   * [`HasFractionField`] — a ring `R` knows its field of fractions and the
//!     canonical embedding `R ↪ Frac(R)`.
//!   * [`HasRingOfIntegers`] — a field `K` knows its ring of integers (the
//!     valuation / integrality subring) and the integrality test `K → R ∪ {⊥}`.
//!
//! They are linked: `HasRingOfIntegers::Int` is bounded by
//! `HasFractionField<Frac = Self>`, so the ring of integers of `K` is a ring
//! whose fraction field is `K` again. That closes the loop at the type level, and
//! the generic round-trip law `frac ∘ int = id` (an embedded ring element is
//! integral and recovers itself) is exercised in [`tests`] for every pair.
//!
//! ## What is and isn't paired
//!
//! The four pairs above are exactly the table rows where the field and its ring
//! of integers are **distinct backends**. The [`Laurent`](crate::scalar::Laurent)
//! functor is the one case where they share a type: the ring of integers `F_q[[t]]`
//! of `F_q((t))` is the valuation subring (`Laurent::is_integral`, valuation `≥ 0`)
//! *inside* the same `Laurent<S, K>`, not a separate world — so it stays outside
//! this trait pairing, honestly, rather than pointing `Int` at itself.
//!
//! ## Functors split on this exact line
//!
//! The two kinds of root-level functor differ precisely in how they treat this
//! pairing, which is a clean structural dichotomy:
//!
//!   * **Algebraic functors transport the distinct-type pairing.**
//!     [`Surcomplex`](crate::scalar::Surcomplex) is `i`-adjunction; if the base `R`
//!     has a fraction field, so does `R[i]` (componentwise), and if `K` has a ring
//!     of integers `O_K`, then `K[i]` has the order `O_K[i]`. The blanket impls
//!     below make this functorial: the `S = Rational` instance is the **Gaussian**
//!     row `ℤ[i] ⊂ ℚ[i]`, the `S = Surreal` instance is `Omnific[i] ⊂ Surcomplex`,
//!     all from one pair of impls. (For `ℚ(i)` the order `ℤ[i]` is the maximal
//!     order; over a general base it is `O_K[i]`, which we do not claim is maximal.)
//!   * **Valuation functors keep a same-type subring.**
//!     [`Laurent`](crate::scalar::Laurent) and
//!     [`Ramified`](crate::scalar::Ramified) adjoin (a transcendental / a
//!     ramified root with) a *valuation*; their ring of integers is the
//!     valuation-`≥ 0` subring of the *same* type (`is_integral`), so they stay
//!     out of the pairing — the same honesty as `Laurent` above.

use crate::scalar::{
    mul_mod_u128, Integer, Omnific, Poly, Qp, Qq, Rational, RationalFunction, Scalar, Surcomplex,
    Surreal, WittVec, Zp,
};

/// A (commutative) ring that knows its field of fractions.
pub trait HasFractionField: Scalar {
    /// The field of fractions `Frac(R)`.
    type Frac: Scalar;
    /// The canonical ring embedding `R ↪ Frac(R)`.
    fn to_fraction(&self) -> Self::Frac;
}

/// A field that knows its ring of integers — the valuation / integrality subring.
pub trait HasRingOfIntegers: Scalar {
    /// The ring of integers, itself a ring whose fraction field is `Self`.
    type Int: HasFractionField<Frac = Self>;
    /// Whether this element lies in the ring of integers.
    fn is_integral(&self) -> bool;
    /// This element as a ring-of-integers element, or `None` if it is not integral.
    fn to_integer(&self) -> Option<Self::Int>;
}

// ───────────────────────── ℤ ⊂ ℚ ─────────────────────────

impl HasFractionField for Integer {
    type Frac = Rational;
    fn to_fraction(&self) -> Rational {
        Rational::int(self.0)
    }
}

impl HasRingOfIntegers for Rational {
    type Int = Integer;
    fn is_integral(&self) -> bool {
        self.is_integer()
    }
    fn to_integer(&self) -> Option<Integer> {
        if self.is_integer() {
            Some(Integer(self.numer()))
        } else {
            None
        }
    }
}

// ───────────────────────── Oz ⊂ No ─────────────────────────

impl HasFractionField for Omnific {
    type Frac = Surreal;
    fn to_fraction(&self) -> Surreal {
        self.inner().clone()
    }
}

impl HasRingOfIntegers for Surreal {
    type Int = Omnific;
    fn is_integral(&self) -> bool {
        Omnific::from_surreal(self.clone()).is_some()
    }
    fn to_integer(&self) -> Option<Omnific> {
        Omnific::from_surreal(self.clone())
    }
}

// ───────────────────────── Z_p ⊂ Q_p ─────────────────────────

impl<const P: u128, const K: u128> HasFractionField for Zp<P, K> {
    type Frac = Qp<P, K>;
    fn to_fraction(&self) -> Qp<P, K> {
        Qp::from_i128((self.0 % Zp::<P, K>::modulus()) as i128)
    }
}

impl<const P: u128, const K: u128> HasRingOfIntegers for Qp<P, K> {
    type Int = Zp<P, K>;
    fn is_integral(&self) -> bool {
        // valuation ≥ 0 (zero has valuation +∞, hence integral).
        self.valuation().is_none_or(|v| v >= 0)
    }
    fn to_integer(&self) -> Option<Zp<P, K>> {
        let Some(v) = self.valuation() else {
            return Some(Zp(0)); // zero
        };
        if v < 0 {
            return None;
        }
        // residue = unit · p^v  (mod p^k)
        let m = Qp::<P, K>::modulus();
        let mut acc = self.unit() % m;
        for _ in 0..v {
            acc = mul_mod_u128(acc, P % m, m);
        }
        Some(Zp(acc))
    }
}

// ───────────────────────── W_N(F_q) ⊂ Q_q ─────────────────────────

impl<const P: u128, const N: usize, const F: usize> HasFractionField for WittVec<P, N, F> {
    type Frac = Qq<P, N, F>;
    fn to_fraction(&self) -> Qq<P, N, F> {
        Qq::from_witt(*self)
    }
}

impl<const P: u128, const N: usize, const F: usize> HasRingOfIntegers for Qq<P, N, F> {
    type Int = WittVec<P, N, F>;
    fn is_integral(&self) -> bool {
        self.valuation().is_none_or(|v| v >= 0)
    }
    fn to_integer(&self) -> Option<WittVec<P, N, F>> {
        let Some(v) = self.valuation() else {
            return Some(WittVec::zero()); // zero
        };
        if v < 0 {
            return None;
        }
        // ring element = unit · p^v  in W_N(F_q)
        let p = WittVec::<P, N, F>::from_int(P as i128);
        let mut acc = self.unit();
        for _ in 0..v {
            acc = acc.mul(&p);
        }
        Some(acc)
    }
}

// ───────────────────────── F_q[t] ⊂ F_q(t) ─────────────────────────
//
// The char-p mirror of `ℤ ⊂ ℚ`, and the one **distinct-type** pairing on the
// finite/function-field row (the local `F_q[[t]] ⊂ F_q((t))` of `Laurent` is the
// same-type valuation subring, so it stays out — see the module note above). Here
// the ring of integers `F_q[t] = Poly<S>` is a genuinely separate backend from the
// global field `RationalFunction<S> = F_q(t)`, so the pairing is structural.

impl<S: Scalar> HasFractionField for Poly<S> {
    type Frac = RationalFunction<S>;
    fn to_fraction(&self) -> RationalFunction<S> {
        RationalFunction::from_poly(self.clone())
    }
}

impl<S: Scalar> HasRingOfIntegers for RationalFunction<S> {
    type Int = Poly<S>;
    fn is_integral(&self) -> bool {
        // integral ⟺ a polynomial ⟺ the (monic) denominator divides the numerator.
        self.num().rem(self.den()).is_zero()
    }
    fn to_integer(&self) -> Option<Poly<S>> {
        let (quot, rem) = self.num().divrem(self.den());
        rem.is_zero().then_some(quot)
    }
}

// ───────────── functorial: Surcomplex transports the pairing ─────────────
//
// The algebraic `i`-adjunction functor preserves the (field, ring of integers)
// relationship: `Frac(R[i]) = Frac(R)[i]` and `O_{K[i]} = O_K[i]`. The round-trip
// law closes because `(K::Int)::Frac = K` already holds for the base pair.

impl<R: HasFractionField> HasFractionField for Surcomplex<R> {
    type Frac = Surcomplex<R::Frac>;
    fn to_fraction(&self) -> Surcomplex<R::Frac> {
        Surcomplex::new(self.re.to_fraction(), self.im.to_fraction())
    }
}

impl<K: HasRingOfIntegers> HasRingOfIntegers for Surcomplex<K> {
    type Int = Surcomplex<K::Int>;
    fn is_integral(&self) -> bool {
        self.re.is_integral() && self.im.is_integral()
    }
    fn to_integer(&self) -> Option<Surcomplex<K::Int>> {
        Some(Surcomplex::new(
            self.re.to_integer()?,
            self.im.to_integer()?,
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fmt::Debug;

    /// The pairing law in one shot: a ring element embeds to an integral field
    /// element that recovers itself (`frac ∘ int = id`), and a non-integral field
    /// element is correctly rejected.
    fn assert_pairs<R>(r: &R)
    where
        R: HasFractionField + PartialEq + Debug,
        R::Frac: HasRingOfIntegers<Int = R>,
    {
        let x = r.to_fraction();
        assert!(
            x.is_integral(),
            "embedded ring element must be integral: {x:?}"
        );
        assert_eq!(
            x.to_integer().as_ref(),
            Some(r),
            "frac∘int round-trip failed"
        );
    }

    #[test]
    fn integer_rational_pairing() {
        for n in -6i128..=6 {
            assert_pairs(&Integer(n));
        }
        // a genuine fraction is not integral
        let half = Rational::new(1, 2);
        assert!(!half.is_integral());
        assert_eq!(half.to_integer(), None);
        // an integer-valued rational is
        assert!(Rational::int(4).is_integral());
        assert_eq!(Rational::int(4).to_integer(), Some(Integer(4)));
    }

    #[test]
    fn omnific_surreal_pairing() {
        assert_pairs(&Omnific::from_int(3));
        assert_pairs(&Omnific::omega()); // ω is an omnific integer
        assert_pairs(&Omnific::from_surreal(Surreal::omega_pow(Surreal::from_int(2))).unwrap());
        // ε and a fractional number are not integral surreals
        assert!(!Surreal::epsilon().is_integral());
        assert_eq!(Surreal::epsilon().to_integer(), None);
        assert!(!Surreal::from_rational(Rational::new(1, 2)).is_integral());
    }

    #[test]
    fn zp_qp_pairing() {
        for r in 0..27u128 {
            assert_pairs(&Zp::<3, 3>(r));
        }
        // 1/p is a genuine fraction: valuation -1, not integral.
        let inv_p = Qp::<3, 3>::from_p_power(-1);
        assert!(!inv_p.is_integral());
        assert_eq!(inv_p.to_integer(), None);
        // p itself IS integral and lands on Zp(p).
        let p = Qp::<3, 3>::from_i128(3);
        assert!(p.is_integral());
        assert_eq!(p.to_integer(), Some(Zp::<3, 3>(3)));
    }

    #[test]
    fn qp_to_integer_uses_modular_multiplication_at_the_boundary() {
        type Q = Qp<3, 80>;
        let x = Q::from_i128(-1).mul(&Q::from_i128(3));
        assert_eq!(x.valuation(), Some(1));
        assert_eq!(x.to_integer(), Some(Zp::<3, 80>(Q::modulus() - 3)));
    }

    #[test]
    fn wittvec_qq_pairing() {
        // W_2(F_4) ⊂ Q_4: every ring element round-trips through the fraction field.
        for code in 0..16u128 {
            assert_pairs(&WittVec::<2, 2, 2>([code & 3, (code >> 2) & 3]));
        }
        // 1/p is not integral in Q_4.
        let inv_p = Qq::<2, 4, 2>::from_p_power(-1);
        assert!(!inv_p.is_integral());
        assert_eq!(inv_p.to_integer(), None);
        // a Witt unit with genuine F_4 residue is integral and recovers itself.
        let u = WittVec::<2, 4, 2>([1, 1]);
        assert!(u.to_fraction().is_integral());
        assert_eq!(u.to_fraction().to_integer(), Some(u));
    }

    #[test]
    fn poly_rational_function_pairing() {
        use crate::scalar::Fp;
        type P = Poly<Fp<5>>;
        // every polynomial round-trips through F_5(t) = Frac(F_5[t]).
        let samples = [
            P::constant(Fp::<5>::new(3)),
            P::x(),
            Poly::new(vec![Fp::<5>::new(1), Fp::<5>::new(0), Fp::<5>::new(2)]), // 2t² + 1
        ];
        for p in &samples {
            assert_pairs(p);
        }
        // a genuine rational function 1/t is not a polynomial.
        let inv_t = RationalFunction::<Fp<5>>::t().inv().unwrap();
        assert!(!inv_t.is_integral());
        assert_eq!(inv_t.to_integer(), None);
        // t²/t IS integral and recovers the polynomial t (the stored form is unreduced).
        let t2_over_t = RationalFunction::new(
            vec![Fp::<5>::new(0), Fp::<5>::new(0), Fp::<5>::new(1)],
            vec![Fp::<5>::new(0), Fp::<5>::new(1)],
        );
        assert!(t2_over_t.is_integral());
        assert_eq!(t2_over_t.to_integer(), Some(P::x()));
    }

    #[test]
    fn surcomplex_transports_the_gaussian_pairing() {
        // ℤ[i] ⊂ ℚ[i], the S = Rational instance of the blanket impls — every
        // Gaussian integer round-trips through the Gaussian rationals.
        for re in -3i128..=3 {
            for im in -3i128..=3 {
                assert_pairs(&Surcomplex::new(Integer(re), Integer(im)));
            }
        }
        // a genuine Gaussian fraction (½i) is not integral.
        let half_i = Surcomplex::new(Rational::int(0), Rational::new(1, 2));
        assert!(!half_i.is_integral());
        assert_eq!(half_i.to_integer(), None);
        // an integer-valued surcomplex recovers itself.
        let g = Surcomplex::new(Rational::int(3), Rational::int(-2));
        assert!(g.is_integral());
        assert_eq!(
            g.to_integer(),
            Some(Surcomplex::new(Integer(3), Integer(-2)))
        );
    }
}
