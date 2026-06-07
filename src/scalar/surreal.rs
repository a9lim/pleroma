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

use crate::scalar::{Ordinal, Rational, Scalar};
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

    // -- The {L|R} / simplicity bridge: finite rationals, dyadics, birthday --

    /// This surreal as a finite rational, if it is one — a single constant
    /// (`ω⁰`) term, or zero. `None` for anything carrying an `ω`-term
    /// (infinite/infinitesimal), which no short game can reach.
    pub fn as_rational(&self) -> Option<Rational> {
        match self.terms.as_slice() {
            [] => Some(Rational::zero()),
            [(e, c)] if e.is_zero() => Some(c.clone()),
            _ => None,
        }
    }

    /// This surreal as a dyadic rational `num / 2^k` — exactly the values a short
    /// partizan game can take ([`crate::games::Game::number_value`]). Returns
    /// `(num, k)` with `num` odd whenever `k > 0`. `None` for non-dyadics.
    pub fn as_dyadic(&self) -> Option<(i128, u32)> {
        let q = self.as_rational()?;
        let den = q.denom();
        if den & (den - 1) != 0 {
            return None; // denominator is not a power of two
        }
        Some((q.numer(), den.trailing_zeros()))
    }

    /// True iff this surreal is a dyadic rational.
    pub fn is_dyadic(&self) -> bool {
        self.as_dyadic().is_some()
    }

    /// The birthday of a dyadic rational — the day it is born in the surreal
    /// construction (`0`↦0, `±n`↦n, `½`↦2, `¼`/`¾`↦3, …), equal to the
    /// [birthday](crate::games::Game::birthday) of its canonical game. `None`
    /// for non-dyadics, whose birthday is an infinite ordinal outside this
    /// finite-support representation.
    pub fn dyadic_birthday(&self) -> Option<u128> {
        let (num, k) = self.as_dyadic()?;
        Some(birthday_dyadic(num, k))
    }

    /// The simplest surreal strictly greater than `self` — the value of `{self|}`
    /// — when `self` is a finite rational. `None` if `self` carries an `ω`-term.
    pub fn simplest_above(&self) -> Option<Surreal> {
        let q = self.as_rational()?;
        let v = if q.sign() == Ordering::Less {
            Rational::zero() // 0 is the simplest number above any negative
        } else {
            Rational::int(q.floor() + 1) // the least integer strictly above q ≥ 0
        };
        Some(Surreal::from_rational(v))
    }

    /// The simplest surreal strictly less than `self` — the value of `{|self}` —
    /// when `self` is a finite rational. `None` if `self` carries an `ω`-term.
    pub fn simplest_below(&self) -> Option<Surreal> {
        Some(self.neg().simplest_above()?.neg())
    }

    /// The unique simplest surreal strictly between `a` and `b` (Conway's
    /// simplicity theorem), realised when that value is dyadic — i.e. when `a`
    /// and `b` are finite rationals with `a < b`. The result is always dyadic.
    /// `None` if either endpoint carries an `ω`-term or `a ≥ b`.
    pub fn simplest_between(a: &Surreal, b: &Surreal) -> Option<Surreal> {
        let (qa, qb) = (a.as_rational()?, b.as_rational()?);
        if qa.cmp(&qb) != Ordering::Less {
            return None;
        }
        Some(Surreal::from_rational(simplest_rational_between(qa, qb)))
    }

    // -- floor / fractional part : the bridge to the omnific integers Oz --

    /// The **floor** ⌊x⌋ — the greatest omnific integer ≤ `x`, as a `Surreal`.
    /// Concretely: keep every infinite term (`ω`-exponent `> 0`, any rational
    /// coefficient), floor the finite constant, and drop every infinitesimal
    /// term (`ω`-exponent `< 0`). If the finite constant is already an integer,
    /// a negative infinitesimal tail borrows one from that integer part. The
    /// result is always an omnific integer ([`crate::scalar::Omnific`]);
    /// `Omnific::floor` wraps it as one. Satisfies `⌊x⌋ ≤ x < ⌊x⌋ + 1`.
    pub fn floor(&self) -> Surreal {
        let mut terms: Vec<(Surreal, Rational)> = Vec::new();
        let mut constant = Rational::zero();
        let mut saw_constant = false;
        let mut infinitesimal_sign = Ordering::Equal;
        for (e, c) in &self.terms {
            match e.sign() {
                Ordering::Greater => terms.push((e.clone(), c.clone())), // infinite term kept
                Ordering::Equal => {
                    constant = c.clone();
                    saw_constant = true;
                }
                Ordering::Less if infinitesimal_sign == Ordering::Equal => {
                    infinitesimal_sign = c.sign();
                }
                Ordering::Less => {} // lower infinitesimals are dominated
            }
        }
        let mut f = constant.floor();
        if (!saw_constant || constant.is_integer()) && infinitesimal_sign == Ordering::Less {
            f -= 1;
        }
        if f != 0 {
            terms.push((Surreal::zero(), Rational::int(f)));
        }
        // terms stay strictly descending (a subset of self's, same order)
        Surreal { terms }
    }

    /// The **fractional part** `x − ⌊x⌋`, always in `[0, 1)` (it may be an
    /// infinitesimal-carrying value such as `½ + ε`).
    pub fn frac(&self) -> Surreal {
        self.sub(&self.floor())
    }

    // -- sign expansion : the canonical surreal-tree encoding (finite case) --

    /// The **sign expansion** of a *dyadic* surreal: the sequence of left/right
    /// turns (`true = +`, `false = −`) on the path from the root `0` to `x` in
    /// the surreal tree. Its length is exactly the
    /// [birthday](Self::dyadic_birthday). `None` for non-dyadics (`1/3`,
    /// `ω`, `ε`, …), whose sign expansions are transfinite and so not finitely
    /// listable here. Inverse of [`from_sign_expansion`](Self::from_sign_expansion).
    ///
    /// Examples: `0 ↦ []`, `1 ↦ [+]`, `2 ↦ [+,+]`, `½ ↦ [+,−]`, `¾ ↦ [+,−,+]`.
    pub fn sign_expansion(&self) -> Option<Vec<bool>> {
        if !self.is_dyadic() {
            return None;
        }
        let x = self.as_rational().unwrap();
        let (mut lo, mut hi): (Option<Rational>, Option<Rational>) = (None, None);
        let mut signs = Vec::new();
        loop {
            let v = simplest_in_cut(&lo, &hi);
            match x.cmp(&v) {
                Ordering::Equal => break,
                Ordering::Greater => {
                    signs.push(true);
                    lo = Some(v);
                }
                Ordering::Less => {
                    signs.push(false);
                    hi = Some(v);
                }
            }
        }
        Some(signs)
    }

    /// The dyadic surreal with the given finite sign expansion (`true = +`), by
    /// walking the surreal tree. The empty sequence is `0`. Inverse of
    /// [`sign_expansion`](Self::sign_expansion).
    pub fn from_sign_expansion(signs: &[bool]) -> Surreal {
        let (mut lo, mut hi): (Option<Rational>, Option<Rational>) = (None, None);
        for &s in signs {
            let v = simplest_in_cut(&lo, &hi);
            if s {
                lo = Some(v);
            } else {
                hi = Some(v);
            }
        }
        Surreal::from_rational(simplest_in_cut(&lo, &hi))
    }

    /// This surreal as a (non-negative) **ordinal**, if it is one: an ordinal is
    /// exactly a surreal whose CNF has all non-negative ordinal exponents and
    /// positive *integer* coefficients (so the surreal value equals the Cantor
    /// normal form). Covers `0`, every natural, `ω`, `ω·n`, `ω^k`, and the
    /// transfinite `ω^ω`, `ω^{ω^ω}`, …. `None` for anything with a negative or
    /// fractional coefficient (`ω−1`, `½ω`) or a non-ordinal exponent (`√ω =
    /// ω^{1/2}`). Recurses only on the strictly-simpler exponents.
    pub fn as_ordinal(&self) -> Option<Ordinal> {
        let mut result = Ordinal::zero();
        for (e, c) in &self.terms {
            if !c.is_integer() || c.sign() != Ordering::Greater {
                return None; // coefficient must be a positive natural
            }
            if e.sign() == Ordering::Less {
                return None; // exponent must be ≥ 0 to be an ordinal power
            }
            let eord = e.as_ordinal()?; // recursion: exponent is strictly simpler
                                        // terms are descending, so ord_add appends in CNF order.
            result = result.ord_add(&Ordinal::monomial(eord, c.numer() as u128));
        }
        Some(result)
    }

    /// The **(possibly transfinite) sign expansion** over the *representable
    /// subclass* — the run-length-encoded ±-sequence whose length is the
    /// birthday. Confident Gonshor cases: `0` (empty); dyadics (the exact finite
    /// path); every non-negative ordinal `α` ↦ `α` pluses, and its negative ↦
    /// `α` minuses (covers `ω`, `ω·n`, `ω^ω`, …); and `ε = ω⁻¹ ↦ +(−)^ω`.
    /// Returns `None` outside that subclass — the honest boundary: `√ω`,
    /// `ω−1`, `½ω`, mixed ordinal+infinitesimal — rather than emitting an
    /// unverified interleaving.
    pub fn transfinite_sign_expansion(&self) -> Option<SignExpansion> {
        if self.is_zero() {
            return Some(SignExpansion { runs: Vec::new() });
        }
        // Dyadic / finite: the exact tree walk, run-length encoded.
        if let Some(signs) = self.sign_expansion() {
            return Some(SignExpansion::from_finite(&signs));
        }
        // A non-negative ordinal is α pluses; its negation, α minuses.
        if let Some(alpha) = self.as_ordinal() {
            if !alpha.is_zero() {
                return Some(SignExpansion {
                    runs: vec![(true, alpha)],
                });
            }
        }
        if let Some(alpha) = self.neg().as_ordinal() {
            if !alpha.is_zero() {
                return Some(SignExpansion {
                    runs: vec![(false, alpha)],
                });
            }
        }
        // ε = ω⁻¹ : one plus, then ω minuses (Gonshor). The one confident
        // infinitesimal; ω^{-k} for k ≥ 2 and rational multiples are out of scope.
        if *self == Surreal::epsilon() {
            return Some(SignExpansion {
                runs: vec![(true, Ordinal::from_u128(1)), (false, Ordinal::omega())],
            });
        }
        None
    }

    /// The **birthday** as an [`Ordinal`]. Dyadics use the fast finite path;
    /// otherwise the birthday is the ordinal *length* of the
    /// [transfinite sign expansion](Self::transfinite_sign_expansion) — so
    /// `ω ↦ ω`, `ω+1 ↦ ω+1`, `ε ↦ ω`, `ω^ω ↦ ω^ω`. `None` outside the
    /// representable subclass (`√ω`, …).
    pub fn birthday_ordinal(&self) -> Option<Ordinal> {
        if let Some(b) = self.dyadic_birthday() {
            return Some(Ordinal::from_u128(b));
        }
        Some(self.transfinite_sign_expansion()?.length())
    }

    // -- lazy/truncated field arithmetic : Hahn-series inversion and roots --

    /// Keep the `n` leading (largest-exponent) terms. Terms are stored strictly
    /// descending, so this is the top-`n` of the Hahn series.
    fn truncate(&self, n: usize) -> Surreal {
        if self.terms.len() <= n {
            self.clone()
        } else {
            Surreal {
                terms: self.terms[..n].to_vec(),
            }
        }
    }

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
    pub fn sqrt(&self, n: usize) -> Option<Surreal> {
        self.nth_root(2, n)
    }

    /// The **truncated real `k`-th root** to `n` leading terms (`k ≥ 1`), or
    /// `None`. `Some` iff the leading coefficient is a perfect ℚ `k`-th power
    /// (and, for even `k`, `self > 0`). See [`sqrt`](Self::sqrt) for the scope.
    pub fn nth_root(&self, k: u32, n: usize) -> Option<Surreal> {
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

/// The simplest dyadic strictly **below** `h` (the value of the cut `{|h}`).
fn simplest_below_rat(h: &Rational) -> Rational {
    if h.sign() == Ordering::Greater {
        Rational::zero() // 0 is the simplest number below any positive
    } else {
        let f = h.floor();
        if Rational::int(f).cmp(h) == Ordering::Less {
            Rational::int(f) // h non-integer: ⌊h⌋ is the closest-to-0 integer below it
        } else {
            Rational::int(f - 1) // h an integer: the next integer down
        }
    }
}

/// The simplest dyadic strictly **above** `l` (the value of the cut `{l|}`).
fn simplest_above_rat(l: &Rational) -> Rational {
    simplest_below_rat(&l.neg()).neg()
}

/// The simplest dyadic strictly inside the open cut `(lo, hi)`; `None` bounds are
/// `∓∞`. This is the surreal-tree node selected at each step of a sign-expansion
/// walk.
fn simplest_in_cut(lo: &Option<Rational>, hi: &Option<Rational>) -> Rational {
    match (lo, hi) {
        (None, None) => Rational::zero(),
        (None, Some(h)) => simplest_below_rat(h),
        (Some(l), None) => simplest_above_rat(l),
        (Some(l), Some(h)) => simplest_rational_between(l.clone(), h.clone()),
    }
}

/// Strip factors of two from a dyadic `num / 2^k` to put it in lowest terms.
fn reduce_dyadic(mut num: i128, mut k: u32) -> (i128, u32) {
    while k > 0 && num % 2 == 0 {
        num /= 2;
        k -= 1;
    }
    (num, k)
}

/// Birthday of the dyadic `num / 2^k` via the canonical `{L|R}` recursion: an
/// integer `n` is born on day `|n|`; a non-integer dyadic on `1 +` the later of
/// its two nearest-dyadic options at `±1/2^k`.
fn birthday_dyadic(num: i128, k: u32) -> u128 {
    if k == 0 {
        return num.unsigned_abs();
    }
    let (ln, lk) = reduce_dyadic(num - 1, k);
    let (rn, rk) = reduce_dyadic(num + 1, k);
    1 + birthday_dyadic(ln, lk).max(birthday_dyadic(rn, rk))
}

/// The simplest dyadic strictly between two rationals `a < b` (the shallowest
/// node of the surreal tree inside the interval).
fn simplest_rational_between(a: Rational, b: Rational) -> Rational {
    // Reflect so the descent only handles the non-negative spine.
    if b.sign() != Ordering::Greater {
        return simplest_rational_between(b.neg(), a.neg()).neg();
    }
    if a.sign() == Ordering::Less {
        return Rational::zero(); // a < 0 < b: 0 is the root, simplest of all
    }
    // 0 ≤ a < b. The least integer strictly above a:
    let c = a.floor() + 1;
    if Rational::int(c).cmp(&b) == Ordering::Less {
        return Rational::int(c); // an integer lives in (a,b); c is closest to 0
    }
    // a and b lie inside one open unit interval (m, m+1).
    let m = a.floor();
    let off = Rational::int(m);
    off.add(&simplest_in_unit(a.sub(&off), b.sub(&off)))
}

/// The shallowest dyadic in `(a, b)` with `0 ≤ a < b ≤ 1`, by binary
/// subdivision of the unit interval.
fn simplest_in_unit(a: Rational, b: Rational) -> Rational {
    let half = Rational::new(1, 2);
    let mut lo = Rational::zero();
    let mut hi = Rational::one();
    loop {
        let mid = lo.add(&hi).mul(&half);
        let below_b = mid.cmp(&b) == Ordering::Less;
        let above_a = a.cmp(&mid) == Ordering::Less;
        if above_a && below_b {
            return mid;
        } else if !above_a {
            lo = mid; // mid ≤ a: search the upper half
        } else {
            hi = mid; // mid ≥ b: search the lower half
        }
    }
}

/// A (possibly transfinite) sign expansion as **runs**: `(sign, length)` pairs,
/// `true = +`, lengths ordinals. A finite expansion is just runs with finite
/// lengths; `ω`-many pluses is a single run `(true, ω)`. The total length (the
/// ordinary ordinal sum of the run lengths) is the surreal's birthday.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SignExpansion {
    runs: Vec<(bool, Ordinal)>,
}

impl SignExpansion {
    /// The runs `(sign, length)`, left to right.
    pub fn runs(&self) -> &[(bool, Ordinal)] {
        &self.runs
    }

    /// The total ordinal length = the birthday (ordinary ordinal sum of runs).
    pub fn length(&self) -> Ordinal {
        let mut len = Ordinal::zero();
        for (_, l) in &self.runs {
            len = len.ord_add(l);
        }
        len
    }

    /// Run-length-encode a finite ±-sequence (`true = +`).
    pub fn from_finite(signs: &[bool]) -> Self {
        let mut runs: Vec<(bool, Ordinal)> = Vec::new();
        for &s in signs {
            if let Some(last) = runs.last_mut() {
                if last.0 == s {
                    last.1 = last.1.ord_add(&Ordinal::from_u128(1));
                    continue;
                }
            }
            runs.push((s, Ordinal::from_u128(1)));
        }
        SignExpansion { runs }
    }

    /// The flat ±-sequence, when every run length is finite; `None` if any run
    /// is transfinite (e.g. `ω`-many signs).
    pub fn as_finite(&self) -> Option<Vec<bool>> {
        let mut out = Vec::new();
        for (s, l) in &self.runs {
            let n = l.as_finite()?;
            for _ in 0..n {
                out.push(*s);
            }
        }
        Some(out)
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
        assert!(w.sqrt(4).unwrap().birthday_ordinal().is_none()); // √ω = ω^{1/2}
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
        let root_w = w.sqrt(4).unwrap();
        assert_eq!(
            root_w,
            Surreal::omega_pow(Surreal::from_rational(Rational::new(1, 2)))
        );
        assert_eq!(root_w.mul(&root_w), w);
        // √4 = 2.
        assert_eq!(int(4).sqrt(4).unwrap(), int(2));
        // √(ω²+2ω+1) = ω+1 in the two leading terms (perfect square, square lead).
        let perfect = Surreal::omega_pow(int(2))
            .add(&Surreal::monomial(int(1), Rational::int(2)))
            .add(&int(1)); // ω² + 2ω + 1
        assert_eq!(perfect.sqrt(2).unwrap(), w.add(&int(1)));
        // The honest ℚ-coefficient boundary: leading coeff not a perfect square
        // ⇒ None (√2, √(2ω)); negative ⇒ None (√−ω). √0 = 0.
        assert!(int(2).sqrt(4).is_none());
        assert!(Surreal::monomial(int(1), Rational::int(2))
            .sqrt(4)
            .is_none());
        assert!(w.neg().sqrt(4).is_none());
        assert_eq!(Surreal::zero().sqrt(4).unwrap(), Surreal::zero());
    }

    #[test]
    fn surreal_nth_roots() {
        let w = Surreal::omega();
        // ∛(ω³) = ω (monomial radicand, exact).
        assert_eq!(Surreal::omega_pow(int(3)).nth_root(3, 4).unwrap(), w);
        // ∛(−8) = −2 (odd root of a negative is allowed).
        assert_eq!(int(-8).nth_root(3, 4).unwrap(), int(-2));
        // ∛2 is irrational ⇒ None.
        assert!(int(2).nth_root(3, 4).is_none());
        // (1+ε)³ = 1 + 3ε + 3ε² + ε³ ; the cube root recovers 1+ε in 2 terms.
        let base = Surreal::one().add(&Surreal::epsilon());
        let cubed = base.mul(&base).mul(&base);
        assert_eq!(cubed.nth_root(3, 2).unwrap(), base);
    }
}
