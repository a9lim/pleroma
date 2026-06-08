//! The **analytic layer**, unified: root-taking and lazy inversion across every
//! coefficient world, behind two traits that split on the one honest difference
//! between the backends — *where the precision lives*.
//!
//! Until now "the analytic layer" was four disconnected sets of inherent methods
//! with two incompatible signatures: the exact `Option<Self>` roots of
//! [`Rational`](crate::scalar::Rational) / the p-adics / the finite fields, and
//! the precision-argument `(n) -> Option<Self>` series of
//! [`Surreal`](crate::scalar::Surreal). This module promotes that split to the
//! type system, the same way [`valued`](crate::scalar::valued) and
//! [`integrality`](crate::scalar::integrality) promoted valuation and the
//! (field, ring-of-integers) pairing.
//!
//! # The two traits
//!
//!   * [`ExactRoots`] — `is_square` + `sqrt`, returning `Option<Self>` with **no
//!     precision argument**. The result is either exact (ℚ, finite fields) or
//!     exact *to the type's own precision* (the const-`K` p-adics, `Laurent`). A
//!     non-square is an honest `None`. Implemented for `Rational`, `Nimber`,
//!     `Fp`, `Fpn`, `Zp`, `Qp`, `Qq`, `WittVec`, `Surreal`, `Laurent`, and —
//!     functorially — `Surcomplex` over an ordered base.
//!   * [`SeriesRoots`] — `sqrt_to_terms` / `nth_root_to_terms` / `inv_to_terms`,
//!     taking a **caller-chosen** term count `n`. This is the lazy-field interface
//!     of [`Surreal`] (the one world with unbounded, not type-fixed, precision),
//!     whose exact inverse / root would have infinite Hahn support.
//!
//! `Surreal` implements **both**: its `SeriesRoots` methods are the primitive
//! (truncated) operations, and its [`ExactRoots::sqrt`] is the exact value
//! recovered by squaring back the truncations until one matches — the
//! `SeriesRoots → ExactRoots` bridge on the representable subdomain.
//!
//! # Surcomplex roots fall out functorially
//!
//! The headline payoff: the algebraic-closure square root of a surcomplex number
//! used to be a private helper buried in [`forms::char0`](crate::forms::char0).
//! Here it is one blanket impl —
//! `impl<R: ExactRoots + Ordered> ExactRoots for Surcomplex<R>` — exactly the
//! shape [`integrality`](crate::scalar::integrality) uses to transport the
//! Gaussian pairing. The complex `√(a+bi)` formula needs the base to be *ordered*
//! (to pick the branch) with exact non-negative roots, so the bound is
//! [`Ordered`], satisfied by `Rational` (→ the Gaussian rationals ℚ[i]) and
//! `Surreal` (→ the surreal-complex field). `Surcomplex<Surreal>` additionally
//! gets a lazy [`inv_to_terms`](Surcomplex::inv_to_terms), making division
//! first-class even when the norm `a²+b²` is a non-monomial.
//!
//! # Honest boundaries
//!
//!   * The const-`K` p-adic **fields/rings assert odd residue characteristic** in
//!     `is_square`/`sqrt` (the dyadic case is the forms mod-8 story, not a Newton
//!     lift); the trait impls inherit that assertion. The finite fields and
//!     `Laurent` handle characteristic 2 natively (inverse Frobenius).
//!   * [`Gauss`](crate::scalar::Gauss) (rational functions) and
//!     [`Ramified`](crate::scalar::Ramified) get **no** `ExactRoots` here: a
//!     rational-function square root needs polynomial perfect-square detection, a
//!     different machine, and the ramified case its own uniformizer bookkeeping.
//!     Left out honestly rather than stubbed.

use crate::scalar::{
    Fp, Fpn, Laurent, Nimber, Qp, Qq, Rational, Scalar, Surcomplex, Surreal, WittVec, Zp,
};
use std::cmp::Ordering;

// ─────────────────────────────── the traits ───────────────────────────────

/// Exact root-taking: a square root that is either exact, exact to the type's
/// fixed precision, or honestly absent. No caller-chosen precision (contrast
/// [`SeriesRoots`]).
pub trait ExactRoots: Scalar {
    /// Whether this element is a square in its world.
    fn is_square(&self) -> bool;

    /// A square root, or `None` if this is not a square (in the represented
    /// subdomain). When there are two roots, the canonical / residue-lifted one.
    fn sqrt(&self) -> Option<Self>;
}

/// Lazy / truncated analytic operations carrying a caller-chosen term count `n`,
/// for worlds whose exact result has infinite support (the Hahn-series surreals).
pub trait SeriesRoots: Scalar {
    /// The `n` leading terms of a real square root, or `None`. See
    /// [`Surreal::sqrt_to_terms`](crate::scalar::Surreal::sqrt_to_terms).
    fn sqrt_to_terms(&self, n: usize) -> Option<Self>;

    /// The `n` leading terms of a real `k`-th root, or `None`. See
    /// [`Surreal::nth_root_to_terms`](crate::scalar::Surreal::nth_root_to_terms).
    fn nth_root_to_terms(&self, k: u32, n: usize) -> Option<Self>;

    /// The `n` leading terms of the multiplicative inverse (a Neumann series for a
    /// non-monomial). See [`Surreal::inv_to_terms`](crate::scalar::Surreal::inv_to_terms).
    fn inv_to_terms(&self, n: usize) -> Option<Self>;
}

/// A scalar carrying a sign — an ordered world. The datum the
/// [`Surcomplex`] blanket [`ExactRoots`] needs to pick the right branch of the
/// complex square root. Deliberately *not* a [`Scalar`] supertrait (the finite
/// and p-adic worlds are unordered), same discipline as
/// [`Valued`](crate::scalar::Valued).
pub trait Ordered: Scalar {
    /// `Greater` / `Less` / `Equal` against zero (the sign of the dominant term).
    fn sign(&self) -> Ordering;
}

// ─────────────────── residue-field square roots (shared) ───────────────────
//
// The Tonelli–Shanks residue-field roots live here, at the analytic root, because
// they are field primitives with no p-adic dependency: `ExactRoots for Fp/Fpn`
// uses them directly, and the Hensel lift in `small/analytic.rs` imports them as
// the lift's seed. (They were previously private to `small/analytic.rs`.)

/// `base^e mod m` (the residue fields here are tiny, so `u128` products suffice).
fn mod_pow(mut base: u128, mut e: u128, m: u128) -> u128 {
    let mut acc = 1u128 % m;
    base %= m;
    while e > 0 {
        if e & 1 == 1 {
            acc = (acc * base) % m;
        }
        base = (base * base) % m;
        e >>= 1;
    }
    acc
}

/// Whether `a` is a square in `F_p` (odd `p`): `a = 0` or `a^{(p−1)/2} = 1`.
pub(crate) fn fp_is_square(a: u128, p: u128) -> bool {
    let a = a % p;
    a == 0 || mod_pow(a, (p - 1) / 2, p) == 1
}

/// A square root of `a` in `F_p` (odd `p`) via Tonelli–Shanks, or `None` if `a`
/// is a non-residue. The returned root is the seed for the Hensel lift.
pub(crate) fn fp_sqrt(a: u128, p: u128) -> Option<u128> {
    let a = a % p;
    if a == 0 {
        return Some(0);
    }
    if mod_pow(a, (p - 1) / 2, p) != 1 {
        return None; // non-residue
    }
    if p % 4 == 3 {
        return Some(mod_pow(a, (p + 1) / 4, p)); // the fast branch
    }
    // p ≡ 1 (mod 4): full Tonelli–Shanks. Write p−1 = q·2^s with q odd.
    let (mut q, mut s) = (p - 1, 0u32);
    while q % 2 == 0 {
        q /= 2;
        s += 1;
    }
    // Find a quadratic non-residue z.
    let mut z = 2u128;
    while mod_pow(z, (p - 1) / 2, p) != p - 1 {
        z += 1;
    }
    let mut m = s;
    let mut c = mod_pow(z, q, p);
    let mut t = mod_pow(a, q, p);
    let mut r = mod_pow(a, q.div_ceil(2), p);
    loop {
        if t == 1 {
            return Some(r);
        }
        // least i in 1..m with t^{2^i} = 1
        let mut i = 1u32;
        let mut t2 = (t * t) % p;
        while t2 != 1 {
            t2 = (t2 * t2) % p;
            i += 1;
        }
        let b = mod_pow(c, 1u128 << (m - i - 1), p);
        m = i;
        c = (b * b) % p;
        t = (t * c) % p;
        r = (r * b) % p;
    }
}

/// A square root of `a` in `F_q = F_{p^N}` via Tonelli–Shanks over the field, or
/// `None` for a non-square. Uses a primitive element as the guaranteed quadratic
/// non-residue. Works in **either characteristic**: in char 2 every element is a
/// square and `q − 1` is odd (`s = 0`), so the loop returns `a^{q/2}` (the inverse
/// Frobenius) on the first step.
pub(crate) fn fq_sqrt<const P: u128, const N: usize>(a: Fpn<P, N>) -> Option<Fpn<P, N>> {
    use crate::scalar::FiniteField;
    if a.is_zero() {
        return Some(Fpn::zero());
    }
    if !a.is_square() {
        return None;
    }
    let one = Fpn::<P, N>::one();
    // q−1 = q'·2^s with q' odd.
    let (mut qodd, mut s) = (Fpn::<P, N>::order() - 1, 0u32);
    while qodd % 2 == 0 {
        qodd /= 2;
        s += 1;
    }
    let z = Fpn::<P, N>::primitive_element(); // a generator ⇒ a non-residue
    let mut m = s;
    let mut c = z.pow(qodd);
    let mut t = a.pow(qodd);
    let mut r = a.pow(qodd.div_ceil(2));
    loop {
        if t == one {
            return Some(r);
        }
        let mut i = 1u32;
        let mut t2 = t.mul(&t);
        while t2 != one {
            t2 = t2.mul(&t2);
            i += 1;
        }
        let b = c.pow(1u128 << (m - i - 1));
        m = i;
        c = b.mul(&b);
        t = t.mul(&c);
        r = r.mul(&b);
    }
}

// ─────────────────────────── ExactRoots impls ───────────────────────────

impl ExactRoots for Rational {
    fn is_square(&self) -> bool {
        // inherent `Rational::sqrt` (a perfect ℚ-square ⇒ Some).
        self.sqrt().is_some()
    }
    fn sqrt(&self) -> Option<Self> {
        // Inherent shadows the trait method here, so this delegates, not recurses.
        Rational::sqrt(self)
    }
}

impl ExactRoots for Nimber {
    fn is_square(&self) -> bool {
        true // char 2: the Frobenius x ↦ x² is a bijection ⇒ every element is a square
    }
    fn sqrt(&self) -> Option<Self> {
        Some(Nimber(crate::scalar::nim_sqrt(self.0)))
    }
}

impl<const P: u128> ExactRoots for Fp<P> {
    fn is_square(&self) -> bool {
        if P == 2 {
            return true; // Frobenius onto in char 2
        }
        fp_is_square(self.0, P)
    }
    fn sqrt(&self) -> Option<Self> {
        if self.0 == 0 {
            return Some(Fp(0));
        }
        if P == 2 {
            return Some(*self); // x² = x for x ∈ {0,1}
        }
        fp_sqrt(self.0, P).map(Fp)
    }
}

impl<const P: u128, const N: usize> ExactRoots for Fpn<P, N> {
    fn is_square(&self) -> bool {
        // inherent `Fpn::is_square` (Euler in odd char, `true` in char 2).
        Fpn::is_square(self)
    }
    fn sqrt(&self) -> Option<Self> {
        fq_sqrt(*self)
    }
}

impl<const P: u128, const K: u128> ExactRoots for Zp<P, K> {
    fn is_square(&self) -> bool {
        // inherent (asserts odd p — the dyadic case is the forms mod-8 story).
        Zp::is_square(self)
    }
    fn sqrt(&self) -> Option<Self> {
        Zp::sqrt(self)
    }
}

impl<const P: u128, const K: u128> ExactRoots for Qp<P, K> {
    fn is_square(&self) -> bool {
        Qp::is_square(self)
    }
    fn sqrt(&self) -> Option<Self> {
        Qp::sqrt(self)
    }
}

impl<const P: u128, const N: usize, const F: usize> ExactRoots for WittVec<P, N, F> {
    fn is_square(&self) -> bool {
        WittVec::is_square(self)
    }
    fn sqrt(&self) -> Option<Self> {
        WittVec::sqrt(self)
    }
}

impl<const P: u128, const N: usize, const F: usize> ExactRoots for Qq<P, N, F> {
    fn is_square(&self) -> bool {
        Qq::is_square(self)
    }
    fn sqrt(&self) -> Option<Self> {
        Qq::sqrt(self)
    }
}

impl ExactRoots for Surreal {
    fn is_square(&self) -> bool {
        ExactRoots::sqrt(self).is_some()
    }
    /// The **exact** real square root on the represented subdomain: square back the
    /// truncated [`SeriesRoots`] roots at growing precision until one squares to
    /// `self`. `None` for negatives, and for radicands outside the
    /// finite-CNF-with-ℚ-coefficients subclass (e.g. `√2`).
    fn sqrt(&self) -> Option<Self> {
        if self.is_zero() {
            return Some(Surreal::zero());
        }
        if self.sign() != Ordering::Greater {
            return None;
        }
        let base = self.terms().len().max(1);
        for n in 1..=(8 * base + 32) {
            let root = self.sqrt_to_terms(n)?;
            if root.mul(&root) == *self {
                return Some(root);
            }
        }
        None
    }
}

/// `Laurent<S, K>` over an [`ExactRoots`] base is itself [`ExactRoots`] to relative
/// precision `K` — the equal-characteristic local-field mirror of the p-adic
/// `Qp::sqrt`. A series is a square iff its valuation is even and its unit series
/// is a square; the unit-series root is Newton's iteration in odd/zero
/// characteristic and the even-exponent inverse-Frobenius in characteristic 2.
impl<S: ExactRoots, const K: usize> ExactRoots for Laurent<S, K> {
    fn is_square(&self) -> bool {
        self.sqrt().is_some()
    }
    fn sqrt(&self) -> Option<Self> {
        if self.is_zero() {
            return Some(Self::zero());
        }
        let v = self.valuation().expect("nonzero has a valuation");
        if v % 2 != 0 {
            return None; // odd valuation ⇒ never a square
        }
        let unit = self.unit_coeffs();
        let root_unit = if S::characteristic() == 2 {
            // (Σ wⱼ tʲ)² = Σ wⱼ² t^{2j}: a square iff every odd-exponent coeff
            // vanishes; then wⱼ = √(u_{2j}) (a square in a finite char-2 field).
            for (i, c) in unit.iter().enumerate() {
                if i % 2 == 1 && !c.is_zero() {
                    return None;
                }
            }
            let mut w = Vec::with_capacity(unit.len().div_ceil(2));
            let mut j = 0;
            while 2 * j < unit.len() {
                w.push(unit[2 * j].sqrt()?);
                j += 1;
            }
            Laurent::<S, K>::from_coeffs(w, 0)
        } else {
            // Newton over the unit series: y ← (y + U·y⁻¹)·½, doubling correct
            // terms each step, seeded by the leading-coefficient root.
            let two_inv = S::one().add(&S::one()).inv()?; // 1/2 (a unit in odd/0 char)
            let half = Laurent::<S, K>::from_scalar(two_inv);
            let u0 = unit[0].clone();
            let seed = u0.sqrt()?; // leading coeff must be a square in S
            let unit_series = Laurent::<S, K>::from_coeffs(unit.to_vec(), 0);
            let mut y = Laurent::<S, K>::from_scalar(seed);
            for _ in 0..64 {
                let yi = y.inv()?;
                let next = unit_series.mul(&yi).add(&y).mul(&half);
                if next == y {
                    break;
                }
                y = next;
            }
            y
        };
        // reattach t^{v/2}
        Some(Laurent::<S, K>::from_t_power(v / 2).mul(&root_unit))
    }
}

/// Functorial: the **algebraic closure** square root. `Surcomplex<R>` over an
/// ordered [`ExactRoots`] base is [`ExactRoots`], via the classical
/// `√(a+bi) = ±(√((|z|+a)/2) + sgn(b)·√((|z|−a)/2) i)`. This is the operation that
/// used to be a private helper in `forms::char0`; it now lives where it belongs,
/// and the classifier calls the trait.
impl<R: ExactRoots + Ordered> ExactRoots for Surcomplex<R> {
    fn is_square(&self) -> bool {
        self.sqrt().is_some()
    }
    fn sqrt(&self) -> Option<Self> {
        if self.is_zero() {
            return Some(Surcomplex::zero());
        }
        let root = if self.im.is_zero() {
            match self.re.sign() {
                Ordering::Greater => Surcomplex::new(self.re.sqrt()?, R::zero()),
                Ordering::Less => Surcomplex::new(R::zero(), self.re.neg().sqrt()?),
                Ordering::Equal => Surcomplex::zero(),
            }
        } else {
            let half = R::one().add(&R::one()).inv()?; // 1/2 (ordered ⇒ char 0)
            let norm_sq = self.re.mul(&self.re).add(&self.im.mul(&self.im));
            let norm = norm_sq.sqrt()?;
            let a2 = norm.add(&self.re).mul(&half);
            let b2 = norm.sub(&self.re).mul(&half);
            let a = a2.sqrt()?;
            let mut b = b2.sqrt()?;
            if self.im.sign() == Ordering::Less {
                b = b.neg();
            }
            Surcomplex::new(a, b)
        };
        // Guard the represented subdomain: only accept an exactly-verifying root.
        if root.mul(&root) == *self {
            Some(root)
        } else {
            None
        }
    }
}

// ─────────────────────────── SeriesRoots impls ───────────────────────────

impl SeriesRoots for Surreal {
    fn sqrt_to_terms(&self, n: usize) -> Option<Self> {
        // Inherent shadows the trait method, so these delegate, not recurse.
        Surreal::sqrt_to_terms(self, n)
    }
    fn nth_root_to_terms(&self, k: u32, n: usize) -> Option<Self> {
        Surreal::nth_root_to_terms(self, k, n)
    }
    fn inv_to_terms(&self, n: usize) -> Option<Self> {
        Surreal::inv_to_terms(self, n)
    }
}

// ─────────────────────────────── Ordered impls ───────────────────────────────

impl Ordered for Rational {
    fn sign(&self) -> Ordering {
        Rational::sign(self)
    }
}

impl Ordered for Surreal {
    fn sign(&self) -> Ordering {
        Surreal::sign(self)
    }
}

// ───────────── lazy surcomplex inversion (the SeriesRoots companion) ─────────
//
// `Surcomplex::inv` (the `Scalar` method) needs the norm `a²+b²` to invert
// exactly — i.e. to be a monomial surreal. When it is a non-monomial, the inverse
// has infinite Hahn support; this gives it to `n` leading terms, reusing the
// base's lazy inverse. The full `SeriesRoots` trait is left Surreal-only: lazy
// complex *roots* would add little over the exact `ExactRoots` `sqrt` above.

impl<S: SeriesRoots> Surcomplex<S> {
    /// The **truncated inverse** `1/(a+bi)` to `n` leading terms: `(a−bi)/(a²+b²)`
    /// with the base norm inverted lazily. Works where [`Scalar::inv`] returns
    /// `None` (a non-monomial norm). `None` only for `0`.
    pub fn inv_to_terms(&self, n: usize) -> Option<Surcomplex<S>> {
        let norm = self.re.mul(&self.re).add(&self.im.mul(&self.im));
        let ninv = norm.inv_to_terms(n)?;
        Some(Surcomplex::new(
            self.re.mul(&ninv),
            self.im.neg().mul(&ninv),
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scalar::Scalar;

    // ---------- the exact-root family ----------

    #[test]
    fn rational_exact_roots() {
        assert!(ExactRoots::is_square(&Rational::int(4)));
        assert_eq!(Rational::int(4).sqrt(), Some(Rational::int(2)));
        assert!(!ExactRoots::is_square(&Rational::int(2)));
        assert_eq!(ExactRoots::sqrt(&Rational::int(2)), None);
        assert_eq!(Rational::new(9, 16).sqrt(), Some(Rational::new(3, 4)));
    }

    #[test]
    fn nimber_is_perfect_field() {
        // every element of F_{2^128} is a square; the root squares back.
        for x in [1u128, 2, 3, 255, 1 << 40] {
            let n = Nimber(x);
            assert!(ExactRoots::is_square(&n));
            let r = ExactRoots::sqrt(&n).unwrap();
            assert_eq!(r.mul(&r), n);
        }
    }

    #[test]
    fn fp_fpn_exact_roots() {
        // F_7: 2 is a square (3²=2), 3 is not.
        assert!(ExactRoots::is_square(&Fp::<7>::new(2)));
        let r = ExactRoots::sqrt(&Fp::<7>::new(2)).unwrap();
        assert_eq!(r.mul(&r), Fp::<7>::new(2));
        assert!(!ExactRoots::is_square(&Fp::<7>::new(3)));
        // F_2 (char 2): every element is a square, sqrt(x) = x.
        assert_eq!(ExactRoots::sqrt(&Fp::<2>::new(1)), Some(Fp::<2>::new(1)));
        // F_8 (char 2): a generator is a square (inverse Frobenius), roots back.
        let g = Fpn::<2, 3>::generator();
        assert!(ExactRoots::is_square(&g));
        let rg = ExactRoots::sqrt(&g).unwrap();
        assert_eq!(rg.mul(&rg), g);
        // F_9 (odd char): squares round-trip, non-squares decline.
        let mut squares = 0;
        for code in 0..9u128 {
            let x = {
                let (c0, c1) = (code % 3, code / 3);
                Fpn::<3, 2>::from_coeffs(&[c0, c1])
            };
            match ExactRoots::sqrt(&x) {
                Some(r) => {
                    assert_eq!(r.mul(&r), x);
                    squares += 1;
                }
                None => assert!(!ExactRoots::is_square(&x)),
            }
        }
        assert_eq!(squares, 5); // 0 and the four nonzero QRs of F_9
    }

    #[test]
    fn padic_delegation_matches_inherent() {
        let x = Qp::<5, 5>::from_i128(4);
        assert_eq!(ExactRoots::sqrt(&x), Qp::sqrt(&x));
        assert!(ExactRoots::is_square(&x));
        let z = Zp::<7, 3>(2);
        assert_eq!(ExactRoots::is_square(&z), Zp::is_square(&z));
    }

    #[test]
    fn surreal_exact_sqrt() {
        // perfect square (ω+1)² round-trips; ω is a monomial square; 2 declines.
        let w = Surreal::omega();
        let perfect = w.add(&Surreal::one()).mul(&w.add(&Surreal::one()));
        assert_eq!(ExactRoots::sqrt(&perfect), Some(w.add(&Surreal::one())));
        assert!(ExactRoots::is_square(&w)); // √ω = ω^{1/2}
        assert_eq!(ExactRoots::sqrt(&Surreal::from_int(2)), None); // √2 outside ℚ-CNF
        assert_eq!(ExactRoots::sqrt(&w.neg()), None); // negative
        assert_eq!(
            ExactRoots::sqrt(&Surreal::from_int(4)),
            Some(Surreal::from_int(2))
        );
    }

    // ---------- surcomplex roots (the functorial payoff) ----------

    #[test]
    fn gaussian_sqrt() {
        type G = Surcomplex<Rational>;
        let g = |re: i128, im: i128| Surcomplex::new(Rational::int(re), Rational::int(im));
        // (2+i)² = 3+4i, so √(3+4i) = 2+i (im > 0 branch).
        let r = ExactRoots::sqrt(&g(3, 4)).unwrap();
        assert_eq!(r.mul(&r), g(3, 4));
        assert!(ExactRoots::is_square(&g(3, 4)));
        // √(-1) = i.
        assert_eq!(ExactRoots::sqrt(&g(-1, 0)), Some(G::i()));
        // a non-Gaussian-square declines.
        assert_eq!(ExactRoots::sqrt(&g(2, 0)), None); // √2 not in ℚ[i]
    }

    #[test]
    fn surcomplex_surreal_sqrt() {
        // √ω over the surreal-complex field is exact (ω^{1/2}); a negative real
        // gives a pure imaginary; an imaginary radicand round-trips.
        let w = Surcomplex::<Surreal>::new(Surreal::omega(), Surreal::zero());
        let rw = ExactRoots::sqrt(&w).unwrap();
        assert_eq!(rw.mul(&rw), w);
        let z = Surcomplex::<Surreal>::new(Surreal::from_int(3), Surreal::from_int(4));
        let rz = ExactRoots::sqrt(&z).unwrap();
        assert_eq!(rz.mul(&rz), z);
    }

    #[test]
    fn surcomplex_lazy_inverse() {
        // 1/(1+i) where the norm 2 is a monomial — but exercise the lazy path:
        // 1/(ω+1 + i) has a non-monomial norm (ω+1)²+1, so Scalar::inv is None
        // while inv_to_terms succeeds and (z · z⁻¹) ≈ 1 to the kept precision.
        let z = Surcomplex::<Surreal>::new(Surreal::omega().add(&Surreal::one()), Surreal::one());
        assert!(z.inv().is_none()); // non-monomial norm
        let zi = z.inv_to_terms(8).unwrap();
        let prod = z.mul(&zi);
        // real part leads with 1 (exponent 0); imaginary part cancels exactly.
        assert_eq!(prod.re.terms()[0].1, Rational::one());
        assert_eq!(prod.re.terms()[0].0, Surreal::zero());
        assert!(prod.im.is_zero());
    }

    // ---------- Laurent roots (the equal-characteristic local field) ----------

    #[test]
    fn laurent_sqrt_char0() {
        type L = Laurent<Rational, 8>;
        let r = |n: i128| Rational::int(n);
        // (1 + t)² = 1 + 2t + t²; its sqrt recovers 1 + t to precision.
        let base = Laurent::<Rational, 8>::from_coeffs(vec![r(1), r(1)], 0);
        let sq = base.mul(&base);
        let root = ExactRoots::sqrt(&sq).unwrap();
        assert_eq!(root.mul(&root), sq);
        // t² · (square) ⇒ square at half valuation.
        let shifted = Laurent::<Rational, 8>::from_t_power(2).mul(&sq);
        let rs = ExactRoots::sqrt(&shifted).unwrap();
        assert_eq!(rs.valuation(), Some(1));
        assert_eq!(rs.mul(&rs), shifted);
        // odd valuation ⇒ not a square.
        assert_eq!(ExactRoots::sqrt(&L::t()), None);
        // a non-square leading coefficient (2) declines.
        assert_eq!(
            ExactRoots::sqrt(&Laurent::<Rational, 8>::from_scalar(r(2))),
            None
        );
    }

    #[test]
    fn laurent_sqrt_char2() {
        // F_8((t)): a square has only even-exponent terms; odd terms ⇒ not a square.
        type L = Laurent<Fpn<2, 3>, 6>;
        let g = Fpn::<2, 3>::generator();
        let one = Fpn::<2, 3>::one();
        let base = Laurent::<Fpn<2, 3>, 6>::from_coeffs(vec![one, g], 0); // 1 + g·t
        let sq = base.mul(&base); // = 1 + g²·t² (cross term 2·g = 0)
        let root = ExactRoots::sqrt(&sq).unwrap();
        assert_eq!(root.mul(&root), sq);
        // an odd-exponent term blocks square-ness.
        let odd = Laurent::<Fpn<2, 3>, 6>::from_coeffs(vec![one, one], 0); // 1 + t
        assert_eq!(ExactRoots::sqrt(&odd), None);
        assert!(!ExactRoots::is_square(&odd));
        let _ = L::one();
    }

    // ---------- generic use ----------

    #[test]
    fn exact_roots_is_generic() {
        fn round_trips<S: ExactRoots>(squares: &[S]) {
            for s in squares {
                if let Some(r) = ExactRoots::sqrt(s) {
                    assert_eq!(r.mul(&r), *s);
                    assert!(ExactRoots::is_square(s));
                }
            }
        }
        round_trips(&[Rational::int(9), Rational::int(2)]);
        round_trips(&[Nimber(7), Nimber(255)]);
        round_trips(&[Fp::<11>::new(3), Fp::<11>::new(5)]);
    }
}
