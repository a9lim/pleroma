//! The **lazy / truncated field** layer for surreals: the Hahn-series operations
//! whose exact result has infinite support, returned to a chosen precision `n`
//! (the surreal analogue of the precision-`k` truncation in `Zp`/`Qp`).
//!
//!   * [`Surreal::inv_to_terms`] — the multiplicative inverse of a *non-monomial*
//!     as a Neumann series (where [`crate::scalar::Scalar::inv`] returns `None`).
//!   * [`Surreal::sqrt`] / [`Surreal::nth_root`] — real roots via the binomial
//!     series, `Some` when the requested leading window can be certified inside
//!     the finite-support/i128-coefficient representation.

use super::Surreal;
use crate::scalar::{Rational, Scalar};
use std::cmp::Ordering;

const SERIES_POWER_LIMIT: usize = 4096;
const SERIES_TERM_BUDGET: usize = 100_000;

impl Surreal {
    /// The **truncated multiplicative inverse**: the `n` leading terms of `1/x`,
    /// summed as the Neumann series of its infinite Hahn expansion. Where
    /// [`Scalar::inv`] returns `None` for any non-monomial (the exact inverse has
    /// infinite support), this returns that inverse to a chosen precision `n` —
    /// the surreal analogue of the precision-`k` truncation in
    /// [`Zp`](crate::scalar::Zp)/[`Qp`](crate::scalar::Qp). `None` for `0`, or
    /// when cancellations do not expose a stable `n`-term window inside the
    /// finite search budget.
    ///
    /// Method: factor `x = m·(1+r)` with `m` the leading monomial and `r` an
    /// infinitesimal (leading exponent `< 0`); then `1/x = m⁻¹·Σ_{k≥0}(−r)^k`,
    /// which converges in the Hahn (valuation) sense because `(−r)^k` leads at
    /// `k·deg(r) → −∞`. Example: `1/(ω+1) = ω⁻¹ − ω⁻² + ω⁻³ − …`.
    pub fn inv_to_terms(&self, n: usize) -> Option<Surreal> {
        if self.is_zero() {
            return None;
        }
        if n == 0 {
            return Some(Surreal::zero());
        }
        let (e0, c0) = self.terms[0].clone();
        let m_inv = Surreal::monomial(e0.neg(), c0.inv()?); // ℚ unit: always Some
        let r = m_inv.mul(self).sub(&Surreal::one()); // x = m·(1+r)
        if r.is_zero() {
            return Some(m_inv); // x was a monomial — exact inverse
        }
        let neg_r = r.neg();
        let mut series = Surreal::one();
        let mut power = Surreal::one();
        for _ in 0..SERIES_POWER_LIMIT {
            power = power.mul(&neg_r);
            if power.is_zero() {
                return Some(m_inv.mul(&series).truncate(n));
            }
            if power.terms.len() > SERIES_TERM_BUDGET {
                return None;
            }
            if leading_below_known_window(&series, n, &power) {
                return Some(m_inv.mul(&series).truncate(n));
            }
            series = checked_surreal_add(&series, &power)?;
            if series.terms.len() > SERIES_TERM_BUDGET {
                return None;
            }
        }
        None
    }

    /// The **truncated real square root** to `n` leading terms, or `None`. `None`
    /// when `self < 0`, or its leading coefficient is not a perfect square in ℚ
    /// (the deliberate ℚ-coefficient boundary: `√2` and `√(2ω)` are `None`),
    /// or when deep cancellation, i128 coefficient overflow, or series-budget
    /// exhaustion prevents constructing the requested window. `Some` is guaranteed
    /// when the leading coefficient is an exact ℚ square and the binomial series
    /// converges within the budget — for example `√ω = ω^{1/2}` (monomial) and
    /// `√(ω²+2ω+1) = ω+1` are always exact in their leading terms.
    ///
    /// This is the lazy ([`SeriesRoots`](crate::scalar::SeriesRoots)) primitive;
    /// for the *exact* value (no precision argument) see the
    /// [`ExactRoots::sqrt`](crate::scalar::ExactRoots::sqrt) impl, which squares
    /// these truncations back until one matches.
    pub fn sqrt_to_terms(&self, n: usize) -> Option<Surreal> {
        self.nth_root_to_terms(2, n)
    }

    /// The **truncated real `k`-th root** to `n` leading terms (`k ≥ 1`), or
    /// `None`. `None` when `k = 0`, the leading coefficient is not a perfect ℚ
    /// `k`-th power, for even `k` when `self ≤ 0`, or when deep cancellation,
    /// i128 coefficient overflow, or series-budget exhaustion prevents constructing
    /// the requested window. See [`sqrt_to_terms`](Self::sqrt_to_terms) for scope.
    pub fn nth_root_to_terms(&self, k: u128, n: usize) -> Option<Surreal> {
        if k == 0 {
            return None;
        }
        if self.is_zero() {
            return Some(Surreal::zero());
        }
        if k.is_multiple_of(2) && self.sign() == Ordering::Less {
            return None; // no real even root of a negative
        }
        let (e0, c0) = self.terms[0].clone();
        // leading root: ω^{e0/k} · c0^{1/k}, the latter exact-in-ℚ or None.
        let root_c0 = c0.nth_root(k)?;
        let e0_over_k = e0.mul(&Surreal::from_rational(Rational::new(1, k as i128)));
        let root_m = Surreal::monomial(e0_over_k, root_c0);
        // (1+r)^{1/k} via the binomial series; r infinitesimal.
        let m_inv = Surreal::monomial(e0.neg(), c0.inv()?);
        let r = m_inv.mul(self).sub(&Surreal::one());
        if r.is_zero() {
            return Some(root_m); // exact (monomial radicand)
        }
        let alpha = Rational::new(1, k as i128);
        let series = binomial_series(&r, alpha, n)?;
        Some(root_m.mul(&series).truncate(n))
    }
}

/// `Σ_j binom(α, j) · r^j` truncated to (about) `n` leading terms, with `r` an
/// infinitesimal (leading exponent `< 0`) so the series converges in the Hahn
/// sense. `binom(α,j) = binom(α,j−1)·(α−(j−1))/j` is accumulated with checked
/// i128 arithmetic; overflow means the requested window is outside this backend's
/// represented coefficient range.
fn binomial_series(r: &Surreal, alpha: Rational, n: usize) -> Option<Surreal> {
    if n == 0 {
        return Some(Surreal::zero());
    }
    let mut series = Surreal::one();
    let mut power = Surreal::one(); // r^j
    let mut coeff = Rational::one(); // binom(α, j)
    for j in 1..=SERIES_POWER_LIMIT {
        power = power.mul(r);
        if power.is_zero() {
            return Some(series.truncate(n));
        }
        if power.terms.len() > SERIES_TERM_BUDGET {
            return None;
        }
        if leading_below_known_window(&series, n, &power) {
            return Some(series.truncate(n));
        }
        coeff = next_binomial_coeff(&coeff, &alpha, j)?;
        if coeff.is_zero() {
            continue; // α a nonneg integer ⇒ later binomials vanish, but keep marching
        }
        let contrib = checked_surreal_scale(&coeff, &power)?;
        series = checked_surreal_add(&series, &contrib)?;
        if series.terms.len() > SERIES_TERM_BUDGET {
            return None;
        }
    }
    None
}

fn leading_below_known_window(series: &Surreal, n: usize, next_power: &Surreal) -> bool {
    n == 0
        || series
            .terms
            .get(n - 1)
            .is_some_and(|(nth_exp, _)| next_power.terms[0].0.cmp(nth_exp) == Ordering::Less)
}

fn gcd_u128(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        let t = b;
        b = a % b;
        a = t;
    }
    a
}

fn checked_rational_mul(a: &Rational, b: &Rational) -> Option<Rational> {
    let mut an = a.numer();
    let mut ad = a.denom();
    let mut bn = b.numer();
    let mut bd = b.denom();

    let g1 = gcd_u128(an.unsigned_abs(), bd as u128);
    if g1 > 1 {
        let g1 = i128::try_from(g1).ok()?;
        an /= g1;
        bd /= g1;
    }
    let g2 = gcd_u128(bn.unsigned_abs(), ad as u128);
    if g2 > 1 {
        let g2 = i128::try_from(g2).ok()?;
        bn /= g2;
        ad /= g2;
    }

    Rational::try_new(an.checked_mul(bn)?, ad.checked_mul(bd)?)
}

fn checked_rational_add(a: &Rational, b: &Rational) -> Option<Rational> {
    let g = gcd_u128(a.denom() as u128, b.denom() as u128).max(1);
    let g = i128::try_from(g).ok()?;
    let lhs_scale = b.denom() / g;
    let rhs_scale = a.denom() / g;
    let num = a
        .numer()
        .checked_mul(lhs_scale)?
        .checked_add(b.numer().checked_mul(rhs_scale)?)?;
    let den = a.denom().checked_mul(lhs_scale)?;
    Rational::try_new(num, den)
}

fn checked_rational_div_usize(a: &Rational, d: usize) -> Option<Rational> {
    let d = i128::try_from(d).ok()?;
    let g = gcd_u128(a.numer().unsigned_abs(), d as u128);
    let g = i128::try_from(g).ok()?;
    let num = a.numer() / g;
    let den_factor = d / g;
    Rational::try_new(num, a.denom().checked_mul(den_factor)?)
}

fn rational_sub_usize(a: &Rational, rhs: usize) -> Option<Rational> {
    let rhs = i128::try_from(rhs).ok()?;
    let scaled_rhs = rhs.checked_mul(a.denom())?;
    Rational::try_new(a.numer().checked_sub(scaled_rhs)?, a.denom())
}

fn next_binomial_coeff(prev: &Rational, alpha: &Rational, j: usize) -> Option<Rational> {
    let shifted = rational_sub_usize(alpha, j - 1)?;
    let num = checked_rational_mul(prev, &shifted)?;
    checked_rational_div_usize(&num, j)
}

fn checked_surreal_scale(coeff: &Rational, x: &Surreal) -> Option<Surreal> {
    let mut terms = Vec::with_capacity(x.terms.len());
    for (exp, c) in &x.terms {
        let scaled = checked_rational_mul(coeff, c)?;
        if !scaled.is_zero() {
            terms.push((exp.clone(), scaled));
        }
    }
    Some(Surreal { terms })
}

fn checked_surreal_add(a: &Surreal, b: &Surreal) -> Option<Surreal> {
    let mut terms = Vec::with_capacity(a.terms.len() + b.terms.len());
    let mut i = 0;
    let mut j = 0;
    while i < a.terms.len() && j < b.terms.len() {
        match a.terms[i].0.cmp(&b.terms[j].0) {
            Ordering::Greater => {
                terms.push(a.terms[i].clone());
                i += 1;
            }
            Ordering::Less => {
                terms.push(b.terms[j].clone());
                j += 1;
            }
            Ordering::Equal => {
                let coeff = checked_rational_add(&a.terms[i].1, &b.terms[j].1)?;
                if !coeff.is_zero() {
                    terms.push((a.terms[i].0.clone(), coeff));
                }
                i += 1;
                j += 1;
            }
        }
    }
    terms.extend_from_slice(&a.terms[i..]);
    terms.extend_from_slice(&b.terms[j..]);
    Some(Surreal { terms })
}
