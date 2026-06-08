//! Surreal numbers in Conway normal form, with *surreal* exponents.
//!
//! A nonzero surreal is uniquely
//!     x = Σ_{i} ω^{y_i} · r_i ,   y_0 > y_1 > ... (surreal exponents),  r_i ≠ 0
//! Here we keep **finite support**, take coefficients `r_i` in ℚ (the exact
//! finite stand-in for ℝ — the honest truncation; true CNF allows any real),
//! and let the exponents `y_i` be fully recursive `Surreal`s. So ω, ε = ω⁻¹,
//! √ω = ω^{1/2}, and even ω^ω are all representable.
//!
//! ## Why it terminates
//!
//! Arithmetic realises the Hahn-series structure ℝ((ω^No)): the ω-map is a
//! group homomorphism ω^a · ω^b = ω^{a+b}, so multiplication adds exponents
//! and convolves coefficients. Every operation (add, mul, compare) recurses
//! only on the **exponents**, which are strictly simpler than the number
//! itself — so for any finite-depth surreal the recursion bottoms out.
//!
//! ## Order
//!
//! The leading (largest-exponent) term dominates everything below it, so
//! sign(x) = sign of the leading coefficient and x < y ⇔ sign(x − y) < 0.
//!
//! ## Layout
//!
//! This module is the CNF **core** — the representation, the canonical form, and
//! the ring/field arithmetic ([`Scalar`]). Three companion files carry the
//! theory built on top, all as further `impl Surreal` blocks:
//!
//!   * `simplicity` — the `{L|R}` / simplicity bridge (dyadic recognition,
//!     birthdays, `simplest_*`) and `floor`/`frac` (the bridge to `Oz`).
//!   * `sign_expansion` — the sign-expansion encoding, finite and (Gonshor)
//!     transfinite, plus the [`SignExpansion`] type and `as_ordinal`.
//!   * `analytic` — the lazy/truncated field layer: Neumann-series inverse and
//!     real `k`-th roots of a non-monomial Hahn series.

mod analytic;
mod sign_expansion;
mod simplicity;

pub use sign_expansion::SignExpansion;

use crate::scalar::{Rational, Scalar};
use std::cmp::Ordering;
use std::fmt;

#[derive(Clone)]
pub struct Surreal {
    /// (exponent, coefficient), strictly descending by exponent, coeffs ≠ 0.
    terms: Vec<(Surreal, Rational)>,
}

/// Sort raw (exponent, coeff) terms into canonical form: descending by
/// exponent **value** (the surreal order — a field operation, since `ω−1 < ω`
/// despite being structurally longer), like exponents merged by ℚ-addition, zero
/// coefficients dropped. The descending-merge recipe is shared with the ordinal
/// backend via [`cnf::merge_descending`](super::cnf::merge_descending); only
/// these three primitives are surreal-specific.
fn canonicalize(raw: Vec<(Surreal, Rational)>) -> Vec<(Surreal, Rational)> {
    super::cnf::merge_descending(raw, |a, b| a.cmp(b), |x, y| x.add(y), |c| c.is_zero())
}

impl Surreal {
    /// The constant ω^0 · q. Zero if q == 0.
    pub fn from_rational(q: Rational) -> Self {
        if q.is_zero() {
            Surreal { terms: Vec::new() }
        } else {
            Surreal {
                terms: vec![(Surreal::zero(), q)],
            }
        }
    }

    pub fn from_int(n: i128) -> Self {
        Surreal::from_rational(Rational::int(n))
    }

    /// A single monomial coeff · ω^exp.
    pub fn monomial(exp: Surreal, coeff: Rational) -> Self {
        if coeff.is_zero() {
            Surreal { terms: Vec::new() }
        } else {
            Surreal {
                terms: vec![(exp, coeff)],
            }
        }
    }

    /// ω^exp (coefficient 1).
    pub fn omega_pow(exp: Surreal) -> Self {
        Surreal::monomial(exp, Rational::one())
    }

    /// ω, the simplest infinite surreal.
    pub fn omega() -> Self {
        Surreal::omega_pow(Surreal::one())
    }

    /// ε = ω⁻¹, an infinitesimal.
    pub fn epsilon() -> Self {
        Surreal::omega_pow(Surreal::from_int(-1))
    }

    /// Total order on surreal *values*: sign of the difference.
    pub fn cmp(&self, other: &Surreal) -> Ordering {
        self.sub(other).sign()
    }

    /// Sign of this number: the sign of its dominant (leading) coefficient.
    pub fn sign(&self) -> Ordering {
        match self.terms.first() {
            None => Ordering::Equal,
            Some((_, c)) => c.sign(),
        }
    }

    pub fn terms(&self) -> &[(Surreal, Rational)] {
        &self.terms
    }

    /// Keep the `n` leading (largest-exponent) terms. Terms are stored strictly
    /// descending, so this is the top-`n` of the Hahn series. Used by the
    /// [`analytic`] layer (and its tests) to bound working precision.
    fn truncate(&self, n: usize) -> Surreal {
        if self.terms.len() <= n {
            self.clone()
        } else {
            Surreal {
                terms: self.terms[..n].to_vec(),
            }
        }
    }
}

impl PartialEq for Surreal {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl Scalar for Surreal {
    fn zero() -> Self {
        Surreal { terms: Vec::new() }
    }

    fn one() -> Self {
        Surreal {
            terms: vec![(Surreal::zero(), Rational::one())],
        }
    }

    fn add(&self, rhs: &Self) -> Self {
        let mut raw = self.terms.clone();
        raw.extend(rhs.terms.iter().cloned());
        Surreal {
            terms: canonicalize(raw),
        }
    }

    fn neg(&self) -> Self {
        Surreal {
            terms: self
                .terms
                .iter()
                .map(|(e, c)| (e.clone(), c.neg()))
                .collect(),
        }
    }

    fn mul(&self, rhs: &Self) -> Self {
        // (Σ ω^{a} r)(Σ ω^{b} s) = Σ ω^{a+b} (r·s)   — exponents add, coeffs multiply
        let mut raw = Vec::with_capacity(self.terms.len() * rhs.terms.len());
        for (a, r) in &self.terms {
            for (b, s) in &rhs.terms {
                raw.push((a.add(b), r.mul(s)));
            }
        }
        Surreal {
            terms: canonicalize(raw),
        }
    }

    fn characteristic() -> u128 {
        0
    }

    fn inv(&self) -> Option<Self> {
        // A monomial coeff·ω^e inverts exactly to (1/coeff)·ω^{-e}. A genuine
        // sum has an inverse of infinite Hahn support (e.g. 1/(ω+1) =
        // ω⁻¹ − ω⁻² + ω⁻³ − …), which this finite representation can't hold.
        // [`Surreal::inv_to_terms`] gives that inverse to a chosen precision.
        if self.terms.len() == 1 {
            let (e, c) = &self.terms[0];
            let cinv = c.inv()?;
            Some(Surreal {
                terms: vec![(e.neg(), cinv)],
            })
        } else {
            None
        }
    }

    fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }
}

/// Format coeff·ω^exp for a *non-negative* magnitude coefficient.
fn fmt_term_mag(e: &Surreal, mag: &Rational) -> String {
    if e.is_zero() {
        return format!("{:?}", mag); // a plain constant
    }
    let base = if *e == Surreal::one() {
        "ω".to_string()
    } else if e.terms.len() == 1 && e.terms[0].0.is_zero() {
        // exponent is a bare rational: ω^2, ω^-1, ω^(1/2) — no parens needed
        format!("ω^{:?}", e.terms[0].1)
    } else {
        format!("ω^({:?})", e)
    };
    if *mag == Rational::one() {
        base
    } else {
        format!("{:?}{}", mag, base)
    }
}

impl fmt::Debug for Surreal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.terms.is_empty() {
            return write!(f, "0");
        }
        let mut s = String::new();
        for (idx, (e, c)) in self.terms.iter().enumerate() {
            let neg = c.sign() == Ordering::Less;
            let mag = if neg { c.neg() } else { c.clone() };
            let term = fmt_term_mag(e, &mag);
            if idx == 0 {
                if neg {
                    s.push('-');
                }
                s.push_str(&term);
            } else {
                s.push_str(if neg { " - " } else { " + " });
                s.push_str(&term);
            }
        }
        write!(f, "{}", s)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::clifford::{CliffordAlgebra, Metric};
    use crate::scalar::Ordinal;

    fn int(n: i128) -> Surreal {
        Surreal::from_int(n)
    }

    #[test]
    fn rational_constants_behave() {
        assert_eq!(int(1).add(&int(1)), int(2));
        assert_eq!(int(3).mul(&int(4)), int(12));
        assert_eq!(int(5).sub(&int(5)), Surreal::zero());
        assert!(int(0).is_zero());
    }

    #[test]
    fn omega_is_bigger_than_every_integer() {
        let w = Surreal::omega();
        assert_eq!(w.cmp(&int(1_000_000)), Ordering::Greater);
        // ω − 1 is still infinite and still below ω
        let w_minus_1 = w.sub(&int(1));
        assert_eq!(w_minus_1.cmp(&int(1_000_000)), Ordering::Greater);
        assert_eq!(w_minus_1.cmp(&w), Ordering::Less);
    }

    #[test]
    fn epsilon_is_a_positive_infinitesimal() {
        let eps = Surreal::epsilon();
        assert_eq!(eps.sign(), Ordering::Greater); // ε > 0
                                                   // ε < any positive rational
        let tiny = Surreal::from_rational(Rational::new(1, 1_000_000));
        assert_eq!(eps.cmp(&tiny), Ordering::Less);
    }

    #[test]
    fn omega_times_epsilon_is_one() {
        // ω^1 · ω^{-1} = ω^0 = 1
        assert_eq!(Surreal::omega().mul(&Surreal::epsilon()), Surreal::one());
    }

    #[test]
    fn omega_squared_and_sqrt_omega() {
        let w = Surreal::omega();
        assert_eq!(w.mul(&w), Surreal::omega_pow(int(2))); // ω·ω = ω^2
        let root = Surreal::omega_pow(Surreal::from_rational(Rational::new(1, 2)));
        assert_eq!(root.mul(&root), w); // (√ω)^2 = ω
    }

    #[test]
    fn difference_of_two_infinites() {
        // (ω + 1)(ω − 1) = ω^2 − 1
        let w = Surreal::omega();
        let lhs = w.add(&int(1)).mul(&w.sub(&int(1)));
        let rhs = Surreal::omega_pow(int(2)).sub(&int(1));
        assert_eq!(lhs, rhs);
    }

    #[test]
    fn recursive_exponent_omega_to_the_omega() {
        // ω^ω: the exponent is itself ω. It dominates ω^n for every finite n.
        let w_to_w = Surreal::omega_pow(Surreal::omega());
        let w_to_100 = Surreal::omega_pow(int(100));
        assert_eq!(w_to_w.cmp(&w_to_100), Ordering::Greater);
        // ω^ω · ω^ω = ω^(2ω)
        assert_eq!(
            w_to_w.mul(&w_to_w),
            Surreal::omega_pow(Surreal::omega().mul(&int(2)))
        );
    }

    #[test]
    fn monomial_inverse() {
        assert_eq!(Surreal::omega().inv().unwrap(), Surreal::epsilon()); // ω⁻¹ = ε
        assert_eq!(Surreal::epsilon().inv().unwrap(), Surreal::omega()); // ε⁻¹ = ω
        let three = int(3);
        assert_eq!(three.mul(&three.inv().unwrap()), Surreal::one()); // 3·⅓ = 1
        let w2 = Surreal::omega_pow(int(2));
        assert_eq!(w2.mul(&w2.inv().unwrap()), Surreal::one()); // ω²·ω⁻² = 1
                                                                // a genuine sum has no finite-support inverse
        assert!(Surreal::omega().add(&int(1)).inv().is_none());
        assert!(Surreal::zero().inv().is_none());
    }

    #[test]
    fn distributivity() {
        let a = Surreal::omega().add(&int(2));
        let b = Surreal::epsilon().add(&int(3));
        let c = Surreal::omega_pow(int(2));
        let lhs = a.mul(&b.add(&c));
        let rhs = a.mul(&b).add(&a.mul(&c));
        assert_eq!(lhs, rhs);
    }

    #[test]
    fn clifford_with_infinite_and_infinitesimal_squares() {
        // e0^2 = ω (infinite), e1^2 = ε (infinitesimal), orthogonal.
        // Then (e0 e1)^2 = -(e0^2)(e1^2) = -(ω·ε) = -1. A unit bivector from a
        // metric with no finite entries at all.
        let alg = CliffordAlgebra::new(
            2,
            Metric::diagonal(vec![Surreal::omega(), Surreal::epsilon()]),
        );
        let e0e1 = alg.mul(&alg.gen(0), &alg.gen(1));
        let sq = alg.mul(&e0e1, &e0e1);
        assert_eq!(sq, alg.scalar(int(-1)));
    }

    fn dyadic(num: i128, den: i128) -> Surreal {
        Surreal::from_rational(Rational::new(num, den))
    }

    #[test]
    fn dyadic_recognition() {
        assert_eq!(int(5).as_dyadic(), Some((5, 0)));
        assert_eq!(dyadic(3, 4).as_dyadic(), Some((3, 2)));
        assert_eq!(dyadic(2, 4).as_dyadic(), Some((1, 1))); // reduces to ½
        assert_eq!(dyadic(1, 3).as_dyadic(), None); // not dyadic
        assert_eq!(Surreal::omega().as_dyadic(), None); // infinite
        assert_eq!(Surreal::epsilon().as_dyadic(), None); // infinitesimal
        assert!(int(0).is_dyadic());
    }

    #[test]
    fn dyadic_birthdays_match_construction() {
        // 0 born day 0; ±n day n; ½ day 2; ¼,¾ day 3; ⅜ day 4.
        assert_eq!(int(0).dyadic_birthday(), Some(0));
        assert_eq!(int(3).dyadic_birthday(), Some(3));
        assert_eq!(int(-2).dyadic_birthday(), Some(2));
        assert_eq!(dyadic(1, 2).dyadic_birthday(), Some(2));
        assert_eq!(dyadic(1, 4).dyadic_birthday(), Some(3));
        assert_eq!(dyadic(3, 4).dyadic_birthday(), Some(3));
        assert_eq!(dyadic(3, 8).dyadic_birthday(), Some(4));
    }

    #[test]
    fn simplest_between_picks_the_root() {
        let s = |a: Surreal, b: Surreal| Surreal::simplest_between(&a, &b).unwrap();
        assert_eq!(s(int(0), int(2)), int(1)); // integer 1
        assert_eq!(s(int(-1), int(1)), int(0)); // 0 straddled ⇒ 0
        assert_eq!(s(int(0), int(1)), dyadic(1, 2)); // ½
        assert_eq!(s(int(1), int(2)), dyadic(3, 2)); // 3/2
        assert_eq!(s(dyadic(1, 3), dyadic(2, 3)), dyadic(1, 2)); // ½ between ⅓,⅔
        assert_eq!(s(dyadic(1, 4), dyadic(1, 2)), dyadic(3, 8)); // 3/8
        assert_eq!(s(int(-2), int(-1)), dyadic(-3, 2)); // negatives
                                                        // disordered / non-finite endpoints ⇒ None
        assert!(Surreal::simplest_between(&int(2), &int(1)).is_none());
        assert!(Surreal::simplest_between(&int(0), &Surreal::omega()).is_none());
    }

    #[test]
    fn simplest_above_and_below() {
        assert_eq!(int(2).simplest_above().unwrap(), int(3)); // {2|} = 3
        assert_eq!(dyadic(1, 2).simplest_above().unwrap(), int(1)); // {½|} = 1
        assert_eq!(int(-1).simplest_above().unwrap(), int(0)); // {-1|} = 0
        assert_eq!(int(-2).simplest_below().unwrap(), int(-3)); // {|-2} = -3
        assert_eq!(int(1).simplest_below().unwrap(), int(0)); // {|1} = 0
    }

    #[test]
    fn floor_and_frac() {
        use crate::scalar::is_omnific_integer;
        let w = Surreal::omega();
        let eps = Surreal::epsilon();
        let half = dyadic(1, 2);
        let cases = [
            (w.add(&half), w.clone()),      // ⌊ω+½⌋ = ω
            (w.sub(&half), w.sub(&int(1))), // ⌊ω−½⌋ = ω−1
            (dyadic(3, 2), int(1)),         // ⌊3/2⌋ = 1
            (dyadic(-3, 2), int(-2)),       // ⌊−3/2⌋ = −2
            (eps.clone(), int(0)),          // ⌊ε⌋ = 0
            (eps.neg(), int(-1)),           // ⌊−ε⌋ = −1
            (int(1).sub(&eps), int(0)),     // ⌊1−ε⌋ = 0
            (w.sub(&eps), w.sub(&int(1))),  // ⌊ω−ε⌋ = ω−1
            (int(5), int(5)),               // ⌊5⌋ = 5
            (w.clone(), w.clone()),         // ⌊ω⌋ = ω
            (
                Surreal::monomial(int(1), Rational::new(1, 2)),
                Surreal::monomial(int(1), Rational::new(1, 2)),
            ), // ⌊½ω⌋ = ½ω
        ];
        for (x, expected) in cases {
            let f = x.floor();
            assert_eq!(f, expected, "floor of {:?}", x);
            assert!(is_omnific_integer(&f), "floor of {:?} not omnific", x);
            // ⌊x⌋ ≤ x < ⌊x⌋ + 1
            assert!(f.cmp(&x) != Ordering::Greater);
            assert!(x.cmp(&f.add(&int(1))) == Ordering::Less);
            // fractional part in [0,1)
            let fr = x.frac();
            assert!(fr.sign() != Ordering::Less);
            assert!(fr.cmp(&int(1)) == Ordering::Less);
            assert_eq!(f.add(&fr), x); // x = ⌊x⌋ + {x}
        }
    }

    #[test]
    fn sign_expansion_round_trips() {
        let cases: [(Surreal, Vec<bool>); 6] = [
            (int(0), vec![]),
            (int(1), vec![true]),
            (int(2), vec![true, true]),
            (int(-1), vec![false]),
            (dyadic(1, 2), vec![true, false]),
            (dyadic(3, 4), vec![true, false, true]),
        ];
        for (s, signs) in &cases {
            assert_eq!(
                s.sign_expansion().as_ref(),
                Some(signs),
                "sign exp of {:?}",
                s
            );
            assert_eq!(&Surreal::from_sign_expansion(signs), s);
            // length is the birthday
            assert_eq!(signs.len() as u128, s.dyadic_birthday().unwrap());
        }
        // a longer sweep of dyadics round-trips
        for num in -8i128..=8 {
            for k in 0..4u32 {
                let s = dyadic(num, 1i128 << k);
                let signs = s.sign_expansion().unwrap();
                assert_eq!(Surreal::from_sign_expansion(&signs), s);
                assert_eq!(signs.len() as u128, s.dyadic_birthday().unwrap());
            }
        }
        // non-dyadics have no finite sign expansion
        assert!(Surreal::from_rational(Rational::new(1, 3))
            .sign_expansion()
            .is_none());
        assert!(Surreal::omega().sign_expansion().is_none());
        assert!(Surreal::epsilon().sign_expansion().is_none());
    }

    #[test]
    fn birthday_ordinal_is_finite_for_dyadics() {
        assert_eq!(
            dyadic(3, 4).birthday_ordinal().unwrap().as_finite(),
            Some(3)
        );
        assert_eq!(int(0).birthday_ordinal().unwrap().as_finite(), Some(0));
        // ω now has a transfinite birthday (was None before P3): birthday(ω) = ω.
        assert_eq!(Surreal::omega().birthday_ordinal(), Some(Ordinal::omega()));
    }

    #[test]
    fn transfinite_sign_expansions_and_birthdays() {
        let w = Surreal::omega();
        // Every ordinal is α-many pluses; its birthday is α.
        assert_eq!(
            w.transfinite_sign_expansion().unwrap().runs(),
            &[(true, Ordinal::omega())]
        );
        assert_eq!(w.birthday_ordinal(), Some(Ordinal::omega()));
        // ω + 1 is a (longer) ordinal, distinct from ε's length.
        let w1 = w.add(&int(1));
        assert_eq!(
            w1.birthday_ordinal(),
            Some(Ordinal::omega().ord_add(&Ordinal::from_u128(1)))
        );
        // ω^ω is an ordinal too — handled, not None.
        let wtw = Surreal::omega_pow(Surreal::omega());
        assert_eq!(
            wtw.birthday_ordinal(),
            Some(Ordinal::omega_pow(Ordinal::omega()))
        );
        // ε = ω⁻¹ : +(−)^ω, length 1 + ω = ω.
        let eps = Surreal::epsilon();
        assert_eq!(
            eps.transfinite_sign_expansion().unwrap().runs(),
            &[(true, Ordinal::from_u128(1)), (false, Ordinal::omega())]
        );
        assert_eq!(eps.birthday_ordinal(), Some(Ordinal::omega())); // 1 + ω = ω
                                                                    // The honest boundary: √ω, ω−1, ½ω are outside the subclass ⇒ None.
        assert!(w.sqrt_to_terms(4).unwrap().birthday_ordinal().is_none()); // √ω = ω^{1/2}
        assert!(w.sub(&int(1)).birthday_ordinal().is_none()); // ω − 1
        assert!(Surreal::monomial(int(1), Rational::new(1, 2))
            .birthday_ordinal()
            .is_none()); // ½ω
                         // Finite dyadics still agree with the flat sign expansion.
        let s = dyadic(3, 4);
        assert_eq!(
            s.transfinite_sign_expansion().unwrap().as_finite(),
            s.sign_expansion()
        );
        // as_ordinal round-trips a few ordinals.
        assert_eq!(int(5).as_ordinal(), Some(Ordinal::from_u128(5)));
        assert_eq!(w.as_ordinal(), Some(Ordinal::omega()));
        assert_eq!(w.sub(&int(1)).as_ordinal(), None);
    }

    #[test]
    fn truncated_inverse_neumann_series() {
        let w = Surreal::omega();
        let x = w.add(&int(1)); // ω + 1
                                // 1/(ω+1) = ω⁻¹ − ω⁻² + ω⁻³ − … : the three leading terms.
        let inv3 = x.inv_to_terms(3).unwrap();
        let expected = Surreal::monomial(int(-1), Rational::one())
            .add(&Surreal::monomial(int(-2), Rational::int(-1)))
            .add(&Surreal::monomial(int(-3), Rational::one()));
        assert_eq!(inv3, expected);
        // x · (1/x) = 1 in the leading term (truncation error below it).
        let prod = x.inv_to_terms(10).unwrap().mul(&x);
        assert_eq!(prod.truncate(1), Surreal::one());
        // a monomial inverts exactly and matches Scalar::inv.
        assert_eq!(w.inv_to_terms(5).unwrap(), Surreal::epsilon());
        // zero has no inverse, at any precision.
        assert!(Surreal::zero().inv_to_terms(3).is_none());
    }

    #[test]
    fn surreal_square_roots() {
        let w = Surreal::omega();
        // √ω = ω^{1/2}, exact (monomial radicand), and squares back to ω.
        let root_w = w.sqrt_to_terms(4).unwrap();
        assert_eq!(
            root_w,
            Surreal::omega_pow(Surreal::from_rational(Rational::new(1, 2)))
        );
        assert_eq!(root_w.mul(&root_w), w);
        // √4 = 2.
        assert_eq!(int(4).sqrt_to_terms(4).unwrap(), int(2));
        // √(ω²+2ω+1) = ω+1 in the two leading terms (perfect square, square lead).
        let perfect = Surreal::omega_pow(int(2))
            .add(&Surreal::monomial(int(1), Rational::int(2)))
            .add(&int(1)); // ω² + 2ω + 1
        assert_eq!(perfect.sqrt_to_terms(2).unwrap(), w.add(&int(1)));
        // The honest ℚ-coefficient boundary: leading coeff not a perfect square
        // ⇒ None (√2, √(2ω)); negative ⇒ None (√−ω). √0 = 0.
        assert!(int(2).sqrt_to_terms(4).is_none());
        assert!(Surreal::monomial(int(1), Rational::int(2))
            .sqrt_to_terms(4)
            .is_none());
        assert!(w.neg().sqrt_to_terms(4).is_none());
        assert_eq!(Surreal::zero().sqrt_to_terms(4).unwrap(), Surreal::zero());
    }

    #[test]
    fn surreal_nth_roots() {
        let w = Surreal::omega();
        // ∛(ω³) = ω (monomial radicand, exact).
        assert_eq!(
            Surreal::omega_pow(int(3)).nth_root_to_terms(3, 4).unwrap(),
            w
        );
        // ∛(−8) = −2 (odd root of a negative is allowed).
        assert_eq!(int(-8).nth_root_to_terms(3, 4).unwrap(), int(-2));
        // ∛2 is irrational ⇒ None.
        assert!(int(2).nth_root_to_terms(3, 4).is_none());
        // (1+ε)³ = 1 + 3ε + 3ε² + ε³ ; the cube root recovers 1+ε in 2 terms.
        let base = Surreal::one().add(&Surreal::epsilon());
        let cubed = base.mul(&base).mul(&base);
        assert_eq!(cubed.nth_root_to_terms(3, 2).unwrap(), base);
    }
}
