//! The **global rational function field** `F_q(t)` — the equal-characteristic
//! mirror of `ℚ` as a global field.
//!
//! This is the char-`p` analogue of the (field, ring-of-integers) pairing `ℚ ⊃ ℤ`:
//! the place-organized table's finite/function-field row gets a genuine global
//! **field** `F_q(t)` whose ring of integers is the polynomial ring `F_q[t]`. Its
//! completions are the local fields the rest of the crate already carries —
//! [`Laurent`](crate::scalar::Laurent) `= F_q((t))` at each monic irreducible place,
//! and `F_q((1/t))` at the degree place `∞` — and its arithmetic feeds the
//! local–global form layer [`forms::function_field`](crate::forms) (Hilbert
//! reciprocity `∏_v (a,b)_v = +1` and Hasse–Minkowski over `F_q(t)`), the exact
//! char-`p` mirror of [`forms::padic`](crate::forms)/[`forms::adelic`](crate::forms).
//!
//! ## Exact, unlike the local precision models
//!
//! Every other function-field-adjacent backend (`Laurent`, `Ramified`, `Gauss`,
//! `Qp`, `Adele`) is a **capped-relative precision model** and is excluded from the
//! exact-ring fuzz. `RationalFunction` is **exact**: a genuine commutative field
//! over an exact finite base, so it *joins* the fuzz suite. The product formula it
//! ultimately witnesses — `deg(zeros) = deg(poles)` — is combinatorial and exact,
//! the cleaner mirror of the `ℚ`-adele's archimedean absolute value.
//!
//! ## Representation
//!
//! `num(t) / den(t)` as a pair of [`Poly`]s, with numerator and denominator
//! gcd-reduced and the denominator normalized monic. This differs deliberately from
//! [`Gauss`](crate::scalar::Gauss), whose capped-precision valuation model keeps
//! unreduced fractions to avoid precision-unstable cancellation. `RationalFunction`
//! is exact, so canonical reduction is safe and keeps global-place arithmetic from
//! growing unnecessary common factors. Like [`Adele`](crate::scalar::Adele), it is
//! deliberately **not** `Valued` (no single canonical uniformizer); the forms layer
//! computes per-place valuations from [`num`](RationalFunction::num) and
//! [`den`](RationalFunction::den).

use crate::scalar::{ExactFieldScalar, Poly, Scalar};
use std::fmt;

/// An element of `F_q(t)` (more generally `S(t)` over any field `S`): `num / den`
/// with `den` monic.
#[derive(Clone)]
pub struct RationalFunction<S: ExactFieldScalar> {
    num: Poly<S>,
    den: Poly<S>,
}

impl<S: ExactFieldScalar> RationalFunction<S> {
    /// Assemble `num / den` (already-`Poly`), gcd-reducing and normalizing the
    /// denominator to monic.
    fn from_polys(num: Poly<S>, den: Poly<S>) -> Self {
        assert!(!den.is_zero(), "RationalFunction: zero denominator");
        if num.is_zero() {
            return RationalFunction {
                num: Poly::zero(),
                den: Poly::one(),
            };
        }
        let gcd = num.gcd(&den);
        let (num, den) = if gcd == Poly::one() {
            (num, den)
        } else {
            let (nq, nr) = num.divrem(&gcd);
            let (dq, dr) = den.divrem(&gcd);
            debug_assert!(nr.is_zero() && dr.is_zero(), "gcd must divide both");
            (nq, dq)
        };
        let lead_inv = den
            .leading()
            .unwrap()
            .inv()
            .expect("a field's nonzero leading coefficient inverts");
        RationalFunction {
            num: num.scale(&lead_inv),
            den: den.scale(&lead_inv),
        }
    }

    /// Build `num / den` from low-degree-first coefficient vectors over `S`.
    pub fn new(num: Vec<S>, den: Vec<S>) -> Self {
        RationalFunction::from_polys(Poly::new(num), Poly::new(den))
    }

    /// A polynomial as a rational function `p / 1`.
    pub fn from_poly(p: Poly<S>) -> Self {
        RationalFunction::from_polys(p, Poly::one())
    }

    /// Embed a base scalar as the constant `s / 1`.
    pub fn from_base(s: S) -> Self {
        RationalFunction::from_poly(Poly::constant(s))
    }

    /// The indeterminate `t`.
    pub fn t() -> Self {
        RationalFunction::from_poly(Poly::x())
    }

    /// The numerator polynomial.
    pub fn num(&self) -> &Poly<S> {
        &self.num
    }

    /// The (monic) denominator polynomial.
    pub fn den(&self) -> &Poly<S> {
        &self.den
    }
}

impl<S: ExactFieldScalar> PartialEq for RationalFunction<S> {
    /// Cross-multiplication: `a/b = c/d ⇔ a·d = c·b`.
    fn eq(&self, other: &Self) -> bool {
        self.num.mul(&other.den) == other.num.mul(&self.den)
    }
}

impl<S: ExactFieldScalar> fmt::Debug for RationalFunction<S> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.den == Poly::one() {
            write!(f, "{:?}", self.num)
        } else {
            write!(f, "[{:?}] / [{:?}]", self.num, self.den)
        }
    }
}

impl<S: ExactFieldScalar> Scalar for RationalFunction<S> {
    fn zero() -> Self {
        RationalFunction {
            num: Poly::zero(),
            den: Poly::one(),
        }
    }

    fn one() -> Self {
        RationalFunction {
            num: Poly::one(),
            den: Poly::one(),
        }
    }

    fn add(&self, rhs: &Self) -> Self {
        // a/b + c/d = (a·d + c·b) / (b·d)
        let num = self.num.mul(&rhs.den).add(&rhs.num.mul(&self.den));
        let den = self.den.mul(&rhs.den);
        RationalFunction::from_polys(num, den)
    }

    fn neg(&self) -> Self {
        RationalFunction {
            num: self.num.neg(),
            den: self.den.clone(),
        }
    }

    fn mul(&self, rhs: &Self) -> Self {
        RationalFunction::from_polys(self.num.mul(&rhs.num), self.den.mul(&rhs.den))
    }

    fn characteristic() -> u128 {
        S::characteristic()
    }

    fn inv(&self) -> Option<Self> {
        if self.num.is_zero() {
            return None;
        }
        // (num/den)⁻¹ = den/num — total on nonzero, no gcd needed.
        Some(RationalFunction::from_polys(
            self.den.clone(),
            self.num.clone(),
        ))
    }

    fn is_zero(&self) -> bool {
        self.num.is_zero()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scalar::Fp;

    type F = RationalFunction<Fp<5>>;

    fn rf(num: &[i128], den: &[i128]) -> F {
        RationalFunction::new(
            num.iter().map(|&n| Fp::<5>::new(n)).collect(),
            den.iter().map(|&n| Fp::<5>::new(n)).collect(),
        )
    }

    #[test]
    fn is_an_exact_field() {
        let samples = [
            F::t(),
            F::from_base(Fp::<5>::new(2)),
            rf(&[1, 1], &[1]),       // 1 + t
            rf(&[1], &[0, 1]),       // 1/t
            rf(&[2, 0, 1], &[1, 1]), // (2 + t²)/(1 + t)
        ];
        for x in &samples {
            let xi = x.inv().expect("nonzero inverts in a field");
            assert_eq!(x.mul(&xi), F::one(), "x·x⁻¹ ≠ 1 for {x:?}");
        }
        assert_eq!(F::zero().inv(), None);
        assert_eq!(F::characteristic(), 5);
    }

    #[test]
    fn cross_multiplication_equality() {
        // t/t = 1; (2t)/2 = t; common factors are removed on construction.
        assert_eq!(rf(&[0, 1], &[0, 1]), F::one());
        assert_eq!(rf(&[0, 2], &[2]), F::t());
        assert_ne!(F::t(), F::one());
    }

    #[test]
    fn fractions_are_gcd_reduced_and_denominator_monic() {
        // (t + 1)(t + 2) / (2(t + 1)) = (t + 2) / 2 = 1 + 3t over F_5.
        let x = rf(&[2, 3, 1], &[2, 2]);
        assert_eq!(x.den(), &Poly::one());
        assert_eq!(x.num(), &Poly::new(vec![Fp::<5>::new(1), Fp::<5>::new(3)]));
    }

    #[test]
    fn ring_axioms_on_a_sample() {
        let es = [
            F::zero(),
            F::one(),
            F::t(),
            F::from_base(Fp::<5>::new(3)),
            rf(&[1, 1], &[1]), // 1 + t
            rf(&[1], &[0, 1]), // 1/t
        ];
        for a in &es {
            assert_eq!(a.add(&F::zero()), *a);
            assert_eq!(a.add(&a.neg()), F::zero());
            assert_eq!(a.mul(&F::one()), *a);
            for b in &es {
                assert_eq!(a.add(b), b.add(a));
                assert_eq!(a.mul(b), b.mul(a));
                for d in &es {
                    assert_eq!(a.add(b).add(d), a.add(&b.add(d)));
                    assert_eq!(a.mul(b).mul(d), a.mul(&b.mul(d)));
                    assert_eq!(a.mul(&b.add(d)), a.mul(b).add(&a.mul(d)));
                }
            }
        }
    }

    #[test]
    fn num_den_accessors_expose_polys_for_the_forms_layer() {
        let x = rf(&[0, 1], &[1, 1]); // t / (1 + t)
        assert_eq!(x.num(), &Poly::new(vec![Fp::<5>::new(0), Fp::<5>::new(1)]));
        assert_eq!(x.den(), &Poly::new(vec![Fp::<5>::new(1), Fp::<5>::new(1)]));
    }
}
