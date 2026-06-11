//! Formal Laurent series `S((t))`: adjoin a transcendental `t` (with a valuation)
//! to any scalar backend.
//!
//! This is the **second root-level functor**, the transcendental twin of
//! [`Surcomplex`](crate::scalar::Surcomplex). The two are the two ways to grow a
//! field:
//!
//!   * `Surcomplex<S>` adjoins an **algebraic** element — a root of `x²+1`.
//!   * `Laurent<S, K>` adjoins a **transcendental** element `t` together with a
//!     `t`-adic **valuation**.
//!
//! Both sit *orthogonal* to the "any number" table (they are functors, not
//! concrete worlds). What `Laurent` does that nothing else in the table does is
//! fill the **equal-characteristic local** cell: over a finite field
//! `Laurent<Fpn<P, N>, K>` is `F_{p^N}((t))`, the characteristic-`p` mirror of
//! [`Qp`](crate::scalar::Qp) (the characteristic-0 / mixed-characteristic local
//! field). Its ring of integers is the power-series ring `F_q[[t]]` (the
//! non-negative-valuation elements). `Laurent<Rational, K>` is `ℚ((t))`;
//! `Laurent<Surreal, K>` is an exotic transfinite-coefficient local field.
//!
//! `Laurent<S, K>` is a **field iff `S` is** — `inv` succeeds exactly when the
//! leading coefficient inverts in `S`.
//!
//! ## Precision contract (capped-relative — read this)
//!
//! There is no finite-memory *exact* model of `S((t))`: an inverse can have
//! infinite `t`-adic support (`1/(1+t) = 1 − t + t² − …`), so any concrete
//! backend truncates. Like [`Qp`](crate::scalar::Qp) this uses the standard
//! **capped-relative** model — `K` significant coefficients relative to the
//! valuation. It is in fact *cleaner* than the p-adic case: coefficients live in
//! `S` independently, so there are no inter-digit carries. Multiplication and
//! inversion are **exact** (valuations add; the unit series stays a unit);
//! **addition is not associative across precision boundaries** — additive
//! cancellation below the retained window reads as `0`, exactly like floating
//! point and exactly like `Qp`. `Laurent` is therefore a *precision model*, not
//! an exact commutative ring: it is used at the forms layer (valuation + residue,
//! both robust to relative precision) and is deliberately **excluded from the
//! exact-ring fuzz suite**.

use crate::scalar::Scalar;
use std::fmt;

/// An element of `S((t))` to relative precision `K`: `t^{val} · unit`, where
/// `unit = u₀ + u₁t + … + u_{m−1}t^{m−1}` (`m ≤ K`, `u₀ ≠ 0`) is a unit of the
/// power-series ring. The field zero is the canonical sentinel `{ unit: [], val: 0 }`.
#[derive(Clone, PartialEq, Eq, Hash)]
pub struct Laurent<S: Scalar, const K: usize> {
    /// Coefficients of the unit series, lowest exponent first; leading coeff
    /// `unit[0] ≠ 0`, no trailing zeros, length `≤ K`. Empty ⇒ the zero element.
    unit: Vec<S>,
    /// The (signed) `t`-adic valuation; `0` for zero.
    val: i128,
}

impl<S: Scalar, const K: usize> Laurent<S, K> {
    pub fn assert_supported_precision() {
        assert!(K > 0, "Laurent<S,K> needs positive precision K, got K={K}");
    }

    /// The relative precision (number of retained significant coefficients).
    pub fn precision() -> usize {
        Self::assert_supported_precision();
        K
    }

    /// Canonicalize a raw `t^{val} · coeffs` into the reduced form: cap to `K`
    /// significant coefficients, strip trailing zeros, then strip leading zeros
    /// (folding them into the valuation). All-zero ⇒ the zero sentinel.
    fn normalized(coeffs: Vec<S>, val: i128) -> Self {
        Self::assert_supported_precision();
        // Leading zeros raise the valuation (the relative-precision window slides
        // up; we keep at most K coefficients from the new leading term).
        let lead = coeffs.iter().position(|c| !c.is_zero());
        let Some(lead) = lead else {
            return Laurent {
                unit: Vec::new(),
                val: 0,
            };
        };
        let mut unit: Vec<S> = coeffs[lead..].to_vec();
        unit.truncate(K);
        // Trailing zeros carry no information and would break canonical equality.
        while unit.last().map(|c| c.is_zero()).unwrap_or(false) {
            unit.pop();
        }
        Laurent {
            unit,
            val: val + lead as i128,
        }
    }

    /// Build `t^{val} · (c₀ + c₁t + …)` from raw coefficients (lowest first).
    pub fn from_coeffs(coeffs: Vec<S>, val: i128) -> Self {
        Self::normalized(coeffs, val)
    }

    /// Embed a scalar as the constant series `s` (valuation `0`).
    pub fn from_scalar(s: S) -> Self {
        Self::normalized(vec![s], 0)
    }

    /// The uniformizer `t = t¹`.
    pub fn t() -> Self {
        Self::assert_supported_precision();
        Laurent {
            unit: vec![S::one()],
            val: 1,
        }
    }

    /// The pure power `t^v` (unit series `1`). `from_t_power(-1)` is `1/t`.
    pub fn from_t_power(v: i128) -> Self {
        Self::assert_supported_precision();
        Laurent {
            unit: vec![S::one()],
            val: v,
        }
    }

    /// The `t`-adic valuation, or `None` for zero (whose valuation is `+∞`).
    pub fn valuation(&self) -> Option<i128> {
        if self.unit.is_empty() {
            None
        } else {
            Some(self.val)
        }
    }

    /// The leading coefficient `u₀` (the residue square-class carrier the forms
    /// layer reads), or `None` for zero.
    pub fn leading_coeff(&self) -> Option<S> {
        self.unit.first().cloned()
    }

    /// Whether this series lies in the power-series ring `S[[t]]` — the **ring of
    /// integers** of `S((t))` — i.e. has non-negative valuation (zero included).
    /// This is the one row of the "any number" table where the ring of integers is
    /// a valuation subring of the *same* type rather than a separate backend, so it
    /// is exposed here as an inherent predicate instead of via the
    /// [`HasRingOfIntegers`](crate::scalar::HasRingOfIntegers) pairing.
    pub fn is_integral(&self) -> bool {
        self.valuation().is_none_or(|v| v >= 0)
    }

    /// The retained coefficients of the unit series (lowest exponent first).
    pub fn unit_coeffs(&self) -> &[S] {
        &self.unit
    }

    /// The coefficient of `t^exp`, or `S::zero()` if outside the retained window.
    pub fn coeff(&self, exp: i128) -> S {
        if self.unit.is_empty() {
            return S::zero();
        }
        let i = exp - self.val;
        // Guard both ends as i128 before any usize cast: a huge positive i
        // would truncate-wrap on the cast and alias a small index.
        if i < 0 || i >= self.unit.len() as i128 {
            S::zero()
        } else {
            self.unit[i as usize].clone()
        }
    }
}

impl<S: Scalar, const K: usize> fmt::Display for Laurent<S, K> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.unit.is_empty() {
            return write!(f, "0 (S((t)))");
        }
        let mut first = true;
        for (i, c) in self.unit.iter().enumerate() {
            if c.is_zero() {
                continue;
            }
            if !first {
                write!(f, " + ")?;
            }
            first = false;
            let e = self.val + i as i128;
            match e {
                0 => write!(f, "{c}")?,
                1 => write!(f, "{c}·t")?,
                _ => write!(f, "{c}·t^{e}")?,
            }
        }
        write!(f, " + O(t^{})", self.val + self.unit.len() as i128)
    }
}

impl<S: Scalar, const K: usize> fmt::Debug for Laurent<S, K> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(self, f)
    }
}

impl<S: Scalar, const K: usize> Scalar for Laurent<S, K> {
    fn zero() -> Self {
        Self::assert_supported_precision();
        Laurent {
            unit: Vec::new(),
            val: 0,
        }
    }

    fn one() -> Self {
        Self::assert_supported_precision();
        Laurent {
            unit: vec![S::one()],
            val: 0,
        }
    }

    fn add(&self, rhs: &Self) -> Self {
        Self::assert_supported_precision();
        if self.unit.is_empty() {
            return rhs.clone();
        }
        if rhs.unit.is_empty() {
            return self.clone();
        }
        // Align on the smaller valuation and add within that relative window:
        // x + y = t^{vlo}·(lo + t^{d}·hi), keeping K coefficients from vlo. Any
        // hi-term shifted past the window vanishes (capped-relative, like Qp).
        let (lo, hi) = if self.val <= rhs.val {
            (self, rhs)
        } else {
            (rhs, self)
        };
        // always ≥ 0 by the sort above; compare against K as i128 BEFORE casting:
        // a gap ≥ K means every hi-coefficient falls outside the retained window (the
        // same precision guard Qp::add applies). Casting a huge i128 to usize first
        // truncates-wraps on 64-bit, turning a correct "hi vanishes" into a silently
        // wrong index offset (and panicking in debug mode).
        let d_i128 = hi.val - lo.val;
        if d_i128 >= K as i128 {
            // hi is entirely outside the relative window — result is lo only.
            let mut coeffs = vec![S::zero(); lo.unit.len().min(K)];
            for (i, c) in lo.unit.iter().take(K).enumerate() {
                coeffs[i] = c.clone();
            }
            return Self::normalized(coeffs, lo.val);
        }
        let d = d_i128 as usize; // safe: 0 ≤ d_i128 < K ≤ usize::MAX
        let len = lo.unit.len().max(d + hi.unit.len()).min(K);
        let mut coeffs = vec![S::zero(); len];
        for (i, c) in lo.unit.iter().enumerate() {
            if i < len {
                coeffs[i] = coeffs[i].add(c);
            }
        }
        for (j, c) in hi.unit.iter().enumerate() {
            let i = j + d;
            if i < len {
                coeffs[i] = coeffs[i].add(c);
            }
        }
        Self::normalized(coeffs, lo.val)
    }

    fn neg(&self) -> Self {
        Self::assert_supported_precision();
        Laurent {
            unit: self.unit.iter().map(|c| c.neg()).collect(),
            val: self.val,
        }
    }

    fn mul(&self, rhs: &Self) -> Self {
        Self::assert_supported_precision();
        if self.unit.is_empty() || rhs.unit.is_empty() {
            return Self::zero();
        }
        // t^{a+b} · (U·V), convolution truncated to K coefficients. The product
        // of unit series is a unit series over a domain; over a base ring with
        // zero divisors the leading coeff may cancel, so we renormalize.
        let len = (self.unit.len() + rhs.unit.len() - 1).min(K);
        let mut coeffs = vec![S::zero(); len];
        for (i, a) in self.unit.iter().enumerate() {
            if i >= len {
                break;
            }
            for (j, b) in rhs.unit.iter().enumerate() {
                let k = i + j;
                if k >= len {
                    break;
                }
                coeffs[k] = coeffs[k].add(&a.mul(b));
            }
        }
        Self::normalized(coeffs, self.val + rhs.val)
    }

    fn characteristic() -> u128 {
        Self::assert_supported_precision();
        // Adjoining a transcendental t does not change the characteristic:
        // F_q((t)) has characteristic p, ℚ((t)) characteristic 0.
        S::characteristic()
    }

    fn inv(&self) -> Option<Self> {
        Self::assert_supported_precision();
        // (t^a·U)^{-1} = t^{-a}·U^{-1}. The unit-series inverse is the standard
        // recurrence w₀ = u₀⁻¹, wₙ = −u₀⁻¹·Σ_{i=1}^{n} uᵢ·w_{n−i}, carried to K
        // terms. Total on nonzero iff the leading coeff inverts in S — THE field
        // property (when S is a field), exactly mirroring Qp::inv.
        let u0 = self.unit.first()?;
        let u0inv = u0.inv()?;
        let mut w = vec![S::zero(); K];
        w[0] = u0inv.clone();
        for n in 1..K {
            let mut acc = S::zero();
            for i in 1..=n {
                if i < self.unit.len() {
                    acc = acc.add(&self.unit[i].mul(&w[n - i]));
                }
            }
            w[n] = u0inv.mul(&acc).neg();
        }
        // w₀ ≠ 0, so normalization only trims trailing zeros / caps length.
        Some(Self::normalized(w, -self.val))
    }

    fn is_zero(&self) -> bool {
        Self::assert_supported_precision();
        self.unit.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scalar::{Fpn, Nimber, Rational};

    type L = Laurent<Rational, 6>; // ℚ((t)) to 6 significant terms

    fn r(n: i128) -> Rational {
        Rational::int(n)
    }

    fn lc(coeffs: &[i128], val: i128) -> L {
        Laurent::from_coeffs(coeffs.iter().map(|&n| r(n)).collect(), val)
    }

    #[test]
    fn ring_basics_hold() {
        let x = lc(&[1, 2, 3], 0); // 1 + 2t + 3t²
        assert_eq!(x.add(&L::zero()), x);
        assert_eq!(x.mul(&L::one()), x);
        assert_eq!(x.add(&x.neg()), L::zero());
        // (1 + t)(1 − t) = 1 − t²
        let a = lc(&[1, 1], 0);
        let b = lc(&[1, -1], 0);
        assert_eq!(a.mul(&b), lc(&[1, 0, -1], 0));
    }

    #[test]
    fn valuations_add_under_multiplication() {
        let a = lc(&[2, 1], 1); // t·(2 + t),  val 1
        let b = lc(&[3], -2); // 3·t^{-2},     val -2
        assert_eq!(a.valuation(), Some(1));
        assert_eq!(b.valuation(), Some(-2));
        assert_eq!(a.mul(&b).valuation(), Some(-1));
        // 1/t drops the valuation by one.
        assert_eq!(a.mul(&L::from_t_power(-1)).valuation(), Some(0));
    }

    #[test]
    fn inverse_is_the_neumann_series() {
        // 1/(1 − t) = 1 + t + t² + t³ + … to K terms.
        let one_minus_t = lc(&[1, -1], 0);
        let inv = one_minus_t.inv().unwrap();
        assert_eq!(inv, lc(&[1, 1, 1, 1, 1, 1], 0));
        assert_eq!(one_minus_t.mul(&inv), L::one());
    }

    #[test]
    fn every_nonzero_inverts_when_base_is_a_field() {
        // 1/t exists (negative valuation) — the defining win of the *field*
        // S((t)) over the *ring* S[[t]].
        let t = L::t();
        let ti = t.inv().unwrap();
        assert_eq!(ti, L::from_t_power(-1));
        assert_eq!(t.mul(&ti), L::one());
        assert_eq!(L::zero().inv(), None);
        // a spread of units × valuations all invert
        for lead in 1..4 {
            for v in -2i128..=2 {
                let x = lc(&[lead, 1, 2], v);
                let xi = x.inv().expect("nonzero must invert over a field");
                assert_eq!(x.mul(&xi), L::one(), "x·x⁻¹ ≠ 1 for {x:?}");
            }
        }
    }

    #[test]
    fn equal_characteristic_local_field_over_finite_base() {
        // The headline cell: F_{p^N}((t)) has characteristic p (not 0 like Qp).
        type F8t = Laurent<Fpn<2, 3>, 4>; // F_8((t))
        assert_eq!(F8t::characteristic(), 2);
        type F9t = Laurent<Fpn<3, 2>, 4>; // F_9((t))
        assert_eq!(F9t::characteristic(), 3);
        // and it is a field: 1/t and unit inverses exist
        let x = F8t::t();
        assert_eq!(x.mul(&x.inv().unwrap()), F8t::one());
    }

    #[test]
    fn relative_precision_cancellation_is_intended() {
        // INTENDED, not a bug (the capped-relative contract, cf. Qp): a unit term
        // beyond the K-window of the smaller valuation vanishes on addition.
        // t^0·1 + t^6·1 at precision 6: the t^6 term is outside [0,6) ⇒ reads as 1.
        let one = L::one();
        let far = L::from_t_power(6);
        assert_eq!(one.add(&far), one);
        // whereas a within-window sum is exact.
        assert_eq!(one.add(&L::from_t_power(2)), lc(&[1, 0, 1], 0));
    }

    #[test]
    fn characteristic_two_base_threads_neg_through_scalar() {
        // Over a char-2 base, neg is identity, so x + x = 0 coefficient-wise.
        type Ln = Laurent<Nimber, 4>;
        let x = Laurent::<Nimber, 4>::from_coeffs(vec![Nimber(3), Nimber(5)], 0);
        assert_eq!(x.add(&x), Ln::zero());
        assert_eq!(x.neg(), x);
    }

    #[test]
    fn zero_precision_is_rejected() {
        type L0 = Laurent<Rational, 0>;
        assert!(std::panic::catch_unwind(L0::one).is_err());
        assert!(std::panic::catch_unwind(L0::t).is_err());
        assert!(std::panic::catch_unwind(|| L0::from_coeffs(vec![Rational::one()], 0)).is_err());
    }

    #[test]
    fn m2_huge_valuation_gap_does_not_panic_or_corrupt() {
        // Regression (audit M-2): `(hi.val - lo.val) as usize` for a gap ≥ K (e.g.
        // i128::MAX) previously overflowed the intermediate `d + hi.unit.len()`
        // computation and could panic in debug mode or silently corrupt in release.
        // The fix compares the gap as i128 against K before any usize cast.
        //
        // Mathematically: t^0·1 + t^{i128::MAX}·1 = 1 (hi is outside the K=6 window).
        let one = L::one();
        let far = L::from_t_power(i128::MAX);
        assert_eq!(one.add(&far), one, "huge gap: hi term must vanish");
        assert_eq!(far.add(&one), one, "symmetric");
        // coeff: the coefficient at t^{i128::MAX} should be 0 for `one`.
        assert_eq!(one.coeff(i128::MAX), Rational::zero());
        // A gap equal to K (= 6 here) is the boundary: result is lo only.
        let at_k = L::from_t_power(6);
        assert_eq!(one.add(&at_k), one, "gap == K: hi vanishes too");
    }
}
