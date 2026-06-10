//! **Finite fields** — the residue worlds, where the field is finite. The whole
//! char trichotomy's finite leg, plus the unramified ring of integers that mirrors
//! `Z_p`:
//!
//!   * [`fp`] — `F_p`, the prime fields (odd characteristic): the residue field of
//!     `Z_p`, and the base of every extension here.
//!   * [`fpn`] — `F_{p^n}`, finite extension fields via shipped `(p,n)`-keyed
//!     reduction polynomials. Completes the odd-char tower *and* the supported
//!     char-2 odd-degree fields the nimbers cannot reach (currently `F_8`).
//!   * [`nimber`] — `On₂` truncated to `F_{2^128}`: the char-2 nim-field where
//!     `add = XOR` and `mul` is the coin-turning game product. The main char-2
//!     backend; the only finite field that is also a game-value field.
//!   * [`wittvec`] — `W_N(F_q)`, the truncated Witt vectors `(Z/p^N)[t]/(f̃)`: the
//!     unramified ring of integers over the residue field `F_q`. The char-p mirror
//!     of `Z_p` (which is `W(F_p)`) — completing the (field, ring of integers)
//!     pattern in positive characteristic.
//!
//! [`nimber`] and [`fpn`] share a finite-field analysis toolkit (Frobenius orbit,
//! degree, minimal polynomial, relative trace/norm, multiplicative order, discrete
//! log). That shared algorithm is the [`FiniteField`] trait below; the per-backend
//! impls supply only the Frobenius map, the field shape, and the group order.

pub mod fp;
pub mod fpn;
pub mod nimber;
pub mod wittvec;

pub use fp::*;
pub use fpn::*;
pub use nimber::*;
pub use wittvec::*;

use crate::scalar::Scalar;

/// The **finite-field analysis toolkit**, shared by [`Nimber`] (`F_{2^128}`) and
/// [`Fpn`] (`F_{p^n}`). The two backends differ only in how a field element is
/// stored and how the Frobenius `x ↦ x^p` acts; *every* derived notion — Galois
/// degree, conjugates, minimal polynomial, relative trace/norm, multiplicative
/// order, primitivity, discrete log — is one algorithm over that data, written
/// once here as default methods.
///
/// An impl supplies five things: the Frobenius map, integer exponentiation, the
/// extension degree `[F : F_p]`, the order of `F*`, and the prime factors of that
/// order. Backends with a sharper algorithm for a derived method (nimber's
/// Pohlig–Hellman [`discrete_log`](FiniteField::discrete_log)) override it.
pub trait FiniteField: Scalar + Copy {
    /// The Frobenius endomorphism `x ↦ x^p` — the generator of `Gal(F / F_p)`.
    fn frobenius(&self) -> Self;

    /// Exponentiation `self^e` by an ordinary integer exponent.
    fn pow(&self, e: u128) -> Self;

    /// The extension degree `[F : F_p]`, so `|F| = p^{ext_degree}`.
    fn ext_degree() -> usize;

    /// The order of the multiplicative group `F* = |F| − 1`.
    fn group_order() -> u128;

    /// The distinct prime factors of [`group_order`](FiniteField::group_order).
    fn group_order_factors() -> Vec<u128>;

    /// The Frobenius applied `k` times: `x ↦ x^{p^k}`.
    fn frobenius_iter(&self, k: usize) -> Self {
        let mut x = *self;
        for _ in 0..k {
            x = x.frobenius();
        }
        x
    }

    /// The **degree** of `self` over `F_p`: the least `d | ext_degree` with
    /// `x^{p^d} = x`, i.e. the dimension of the smallest subfield containing it.
    fn degree(&self) -> usize {
        for d in divisors(Self::ext_degree()) {
            if self.frobenius_iter(d) == *self {
                return d;
            }
        }
        Self::ext_degree()
    }

    /// The distinct **Galois conjugates** `x, x^p, …, x^{p^{d-1}}` (`d = degree`) —
    /// the roots of the minimal polynomial, each once.
    fn conjugates(&self) -> Vec<Self> {
        let d = self.degree();
        let mut out = Vec::with_capacity(d);
        let mut c = *self;
        for _ in 0..d {
            out.push(c);
            c = c.frobenius();
        }
        out
    }

    /// The monic **minimal polynomial** over `F_p`, as field elements (each
    /// landing in the prime subfield), constant term first. Formed as
    /// `∏ (X − xᵢ)` over the conjugates; the Galois-closed orbit guarantees the
    /// coefficients lie in `F_p`. Backends project this to their native
    /// coefficient type (`nim_min_poly`, `Fpn::min_poly`).
    fn min_poly_monic(&self) -> Vec<Self> {
        let mut poly = vec![Self::one()]; // the constant polynomial 1
        for c in self.conjugates() {
            let neg_c = c.neg(); // in char 2, −c = c (Nimber::neg is identity)
            let mut next = vec![Self::zero(); poly.len() + 1];
            for (i, a) in poly.iter().enumerate() {
                next[i + 1] = next[i + 1].add(a); // X · a
                next[i] = next[i].add(&neg_c.mul(a)); // (−c) · a
            }
            poly = next;
        }
        poly
    }

    /// The **relative trace** `Tr_{F_{p^m}/F_{p^e}}(x) = Σ_i x^{p^{e·i}}`, viewing
    /// `x` inside the degree-`m` subfield. Requires `e | m`.
    fn relative_trace_over(&self, m: usize, e: usize) -> Self {
        assert!(e > 0 && m.is_multiple_of(e), "relative trace needs e | m");
        let mut acc = Self::zero();
        let mut t = *self;
        for _ in 0..(m / e) {
            acc = acc.add(&t);
            t = t.frobenius_iter(e);
        }
        acc
    }

    /// The **relative norm** `N_{F_{p^m}/F_{p^e}}(x) = ∏_i x^{p^{e·i}}` — the
    /// multiplicative companion of [`relative_trace_over`](Self::relative_trace_over). Requires `e | m`.
    fn relative_norm_over(&self, m: usize, e: usize) -> Self {
        assert!(e > 0 && m.is_multiple_of(e), "relative norm needs e | m");
        let mut acc = Self::one();
        let mut t = *self;
        for _ in 0..(m / e) {
            acc = acc.mul(&t);
            t = t.frobenius_iter(e);
        }
        acc
    }

    /// The trace from the **full** field to `F_{p^e}` (`e = 1` is the absolute
    /// trace to the prime field).
    fn relative_trace(&self, e: usize) -> Self {
        self.relative_trace_over(Self::ext_degree(), e)
    }

    /// The norm from the **full** field to `F_{p^e}`.
    fn relative_norm(&self, e: usize) -> Self {
        self.relative_norm_over(Self::ext_degree(), e)
    }

    /// The **multiplicative order** of `self` in `F*` (`None` for `0`); divides
    /// [`group_order`](FiniteField::group_order).
    fn multiplicative_order(&self) -> Option<u128> {
        if self.is_zero() {
            return None;
        }
        let mut ord = Self::group_order();
        for p in Self::group_order_factors() {
            while ord % p == 0 && self.pow(ord / p) == Self::one() {
                ord /= p;
            }
        }
        Some(ord)
    }

    /// Whether `self` generates the whole multiplicative group `F*`.
    fn is_primitive(&self) -> bool {
        self.multiplicative_order() == Some(Self::group_order())
    }

    /// **Discrete logarithm** to base `self`: the least `e ≥ 0` with `self^e = x`,
    /// or `None` if `x ∉ ⟨self⟩`. The default brute-forces the cyclic subgroup
    /// `⟨self⟩`; backends with a faster route (nimber's Pohlig–Hellman) override.
    fn discrete_log(&self, x: Self) -> Option<u128> {
        if self.is_zero() {
            return None;
        }
        let n = self.multiplicative_order()?;
        let mut cur = Self::one();
        for e in 0..n {
            if cur == x {
                return Some(e);
            }
            cur = cur.mul(self);
        }
        None
    }
}

/// The divisors of `n` in ascending order (small `n` = a field's extension degree).
fn divisors(n: usize) -> Vec<usize> {
    (1..=n).filter(|d| n.is_multiple_of(*d)).collect()
}
