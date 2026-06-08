//! The **lazy / truncated field** layer for surreals: the Hahn-series operations
//! whose exact result has infinite support, returned to a chosen precision `n`
//! (the surreal analogue of the precision-`k` truncation in `Zp`/`Qp`).
//!
//!   * [`Surreal::inv_to_terms`] — the multiplicative inverse of a *non-monomial*
//!     as a Neumann series (where [`crate::scalar::Scalar::inv`] returns `None`).
//!   * [`Surreal::sqrt`] / [`Surreal::nth_root`] — real roots via the binomial
//!     series, `Some` exactly when the leading coefficient is a perfect ℚ-power.

use super::Surreal;
use crate::scalar::{Rational, Scalar};
use std::cmp::Ordering;

impl Surreal {
    /// The **truncated multiplicative inverse**: the `n` leading terms of `1/x`,
    /// summed as the Neumann series of its infinite Hahn expansion. Where
    /// [`Scalar::inv`] returns `None` for any non-monomial (the exact inverse has
    /// infinite support), this returns that inverse to a chosen precision `n` —
    /// the surreal analogue of the precision-`k` truncation in
    /// [`Zp`](crate::scalar::Zp)/[`Qp`](crate::scalar::Qp). `None` only for `0`.
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
        let w = 2 * n + 8; // internal working width, final trimmed to n
        let mut series = Surreal::one();
        let mut power = Surreal::one();
        for _ in 0..(4 * w + 16) {
            power = power.mul(&neg_r).truncate(w);
            if power.is_zero() {
                break;
            }
            if series.terms.len() >= w
                && power.terms[0].0.cmp(&series.terms[w - 1].0) == Ordering::Less
            {
                break; // this (and all smaller) powers no longer reach the window
            }
            series = series.add(&power).truncate(w);
        }
        Some(m_inv.mul(&series).truncate(n))
    }

    /// The **truncated real square root** to `n` leading terms, or `None`. `Some`
    /// iff `self ≥ 0` **and** its leading coefficient is a perfect square in ℚ —
    /// the deliberate ℚ-coefficient boundary: `√2` and `√(2ω)` are `None`
    /// (`√2` is not a finite-CNF-with-ℚ-coeffs surreal), while `√ω = ω^{1/2}`
    /// and `√(ω²+2ω+1) = ω+1` are exact in their leading terms.
    ///
    /// This is the lazy ([`SeriesRoots`](crate::scalar::SeriesRoots)) primitive;
    /// for the *exact* value (no precision argument) see the
    /// [`ExactRoots::sqrt`](crate::scalar::ExactRoots::sqrt) impl, which squares
    /// these truncations back until one matches.
    pub fn sqrt_to_terms(&self, n: usize) -> Option<Surreal> {
        self.nth_root_to_terms(2, n)
    }

    /// The **truncated real `k`-th root** to `n` leading terms (`k ≥ 1`), or
    /// `None`. `Some` iff the leading coefficient is a perfect ℚ `k`-th power
    /// (and, for even `k`, `self > 0`). See [`sqrt_to_terms`](Self::sqrt_to_terms)
    /// for the scope.
    pub fn nth_root_to_terms(&self, k: u32, n: usize) -> Option<Surreal> {
        if k == 0 {
            return None;
        }
        if self.is_zero() {
            return Some(Surreal::zero());
        }
        if k % 2 == 0 && self.sign() == Ordering::Less {
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
        let series = binomial_series(&r, alpha, n);
        Some(root_m.mul(&series).truncate(n))
    }
}

/// `Σ_j binom(α, j) · r^j` truncated to (about) `n` leading terms, with `r` an
/// infinitesimal (leading exponent `< 0`) so the series converges in the Hahn
/// sense. `binom(α,j) = binom(α,j−1)·(α−(j−1))/j` accumulated over ℚ.
fn binomial_series(r: &Surreal, alpha: Rational, n: usize) -> Surreal {
    let w = 2 * n + 8;
    let mut series = Surreal::one();
    let mut power = Surreal::one(); // r^j
    let mut coeff = Rational::one(); // binom(α, j)
    for j in 1..=(4 * w + 16) {
        let jm1 = Rational::int((j - 1) as i128);
        let jr = Rational::int(j as i128);
        coeff = coeff.mul(&alpha.sub(&jm1)).mul(&jr.inv().unwrap());
        power = power.mul(r).truncate(w);
        if power.is_zero() {
            break;
        }
        if coeff.is_zero() {
            continue; // α a nonneg integer ⇒ later binomials vanish, but keep marching
        }
        let contrib = Surreal::monomial(Surreal::zero(), coeff.clone())
            .mul(&power)
            .truncate(w);
        if series.terms.len() >= w
            && contrib.terms[0].0.cmp(&series.terms[w - 1].0) == Ordering::Less
        {
            break;
        }
        series = series.add(&contrib).truncate(w);
    }
    series
}
