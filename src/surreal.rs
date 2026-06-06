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

use crate::scalar::{Rational, Scalar};
use std::cmp::Ordering;
use std::fmt;

#[derive(Clone)]
pub struct Surreal {
    /// (exponent, coefficient), strictly descending by exponent, coeffs ≠ 0.
    terms: Vec<(Surreal, Rational)>,
}

/// Sort raw (exponent, coeff) terms into canonical form: descending by
/// exponent value, like exponents merged, zero coefficients dropped.
fn canonicalize(mut raw: Vec<(Surreal, Rational)>) -> Vec<(Surreal, Rational)> {
    raw.sort_by(|a, b| b.0.cmp(&a.0)); // descending by exponent
    let mut out: Vec<(Surreal, Rational)> = Vec::new();
    for (exp, coeff) in raw {
        if let Some(last) = out.last_mut() {
            if last.0.cmp(&exp) == Ordering::Equal {
                last.1 = last.1.add(&coeff);
                continue;
            }
        }
        out.push((exp, coeff));
    }
    out.retain(|(_, c)| !c.is_zero());
    out
}

impl Surreal {
    /// The constant ω^0 · q. Zero if q == 0.
    pub fn from_rational(q: Rational) -> Self {
        if q.is_zero() {
            Surreal { terms: Vec::new() }
        } else {
            Surreal { terms: vec![(Surreal::zero(), q)] }
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
            Surreal { terms: vec![(exp, coeff)] }
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
        Surreal { terms: vec![(Surreal::zero(), Rational::one())] }
    }

    fn add(&self, rhs: &Self) -> Self {
        let mut raw = self.terms.clone();
        raw.extend(rhs.terms.iter().cloned());
        Surreal { terms: canonicalize(raw) }
    }

    fn neg(&self) -> Self {
        Surreal {
            terms: self.terms.iter().map(|(e, c)| (e.clone(), c.neg())).collect(),
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
        Surreal { terms: canonicalize(raw) }
    }

    fn characteristic() -> u32 {
        0
    }

    fn inv(&self) -> Option<Self> {
        // A monomial coeff·ω^e inverts exactly to (1/coeff)·ω^{-e}. A genuine
        // sum has an inverse of infinite Hahn support (e.g. 1/(ω+1) =
        // ω⁻¹ − ω⁻² + ω⁻³ − …), which this finite representation can't hold.
        if self.terms.len() == 1 {
            let (e, c) = &self.terms[0];
            let cinv = c.inv()?;
            Some(Surreal { terms: vec![(e.neg(), cinv)] })
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
}
