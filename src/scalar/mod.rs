//! The scalar interface every Clifford backend implements.
//!
//! A Clifford algebra needs a *commutative ring* of scalars. The whole point of
//! this project is that combinatorial games only supply such a ring on their
//! field-like subclasses — nimbers, surreals, surcomplex — so each of those is a
//! `Scalar` impl, and the multivector engine in `clifford/` is written once,
//! generic over this trait.
//!
//! This module is the trait; every coefficient world is a descendant module,
//! re-exported flat (`scalar::Nimber`, `scalar::Surreal`, …) so public paths stay
//! shallow regardless of how deep the family tree goes.
//!
//! # The "any number" table
//!
//! The backends are grouped by *place* — the kind of number — and almost every
//! field ships with its **ring of integers**, the same (field, ring) pattern four
//! times over:
//!
//! | place | field | ring of integers | residue |
//! |---|---|---|---|
//! | [`exact`]        — Archimedean        | `Rational` ℚ       | `Integer` ℤ   | — |
//! | [`big`]          — transfinite        | `Surreal` No       | `Omnific` Oz  | ≈ℝ |
//! | [`big`]          — transfinite char-2 | `Ordinal` On₂      | (itself)      | — |
//! | [`small`]        — p-adic             | `Qp` Q_p           | `Zp` Z_p      | F_p |
//! | [`small`]        — p-adic, unramified | `Qq` Q_q           | `WittVec` W_N | F_q |
//! | [`finite_field`] — finite             | `Fp`/`Fpn` F_{p^n} | (itself)      | — |
//! | [`finite_field`] — char-2 nim         | `Nimber` F_2¹²⁸    | (itself)      | — |
//! | [`global`]       — all places at once | `Adele` A_Q model  | integral predicate | — |
//!
//! The **residue** column is itself structural, via [`residue`] ([`ResidueField`]):
//! the discretely-valued local fields know their residue field `k = 𝒪/𝔪` and the
//! reduction `𝒪 → k` (`Qp → F_p`, `Qq → F_q`, `Laurent → F_q`), the last piece of
//! the local-field package `(K, 𝒪, 𝔪, k, Γ, ϖ)` to leave the doc comments — joining
//! [`integrality`] (the `𝒪`/`K` pairing), [`valued`] (`Γ`, `ϖ`), and [`analytic`]
//! (roots). It is what lets the discrete Springer decomposition be written once.
//!
//! The [`global`] family is the place-organized table's local-global row: every
//! other row picks *one* place, while `Adele` is a finite-precision model of the
//! restricted product over all rational places (product formula, Hilbert
//! reciprocity, adelic Hasse–Minkowski; see [`forms::adelic`](crate::forms::adelic)).
//! Its runtime-prime cell [`LocalQp`] fills the const-generic gap the table
//! otherwise cannot represent.
//!
//! `Ordinal` On₂ (the transfinite nimbers, [`big::ordinal`]) is algebraically
//! closed of characteristic 2, not a local field — so its ring-of-integers cell is
//! "(itself)", honestly vacuous, exactly as for the finite fields.
//!
//! The **equal-characteristic local** cell — `F_q((t))` over `F_q[[t]]`, the
//! char-`p` mirror of the `Qp`/`Zp` row — is filled by the [`Laurent`] functor
//! (below), not a row of its own.
//!
//! The [`functor`] module sits *orthogonal* to the table — the ways to grow a
//! field, by an algebraic root or a transcendental, residue- or value-extending
//! (see [`functor`] for the full 2×2 square):
//!   * [`Surcomplex`] is `Surcomplex<S>` — a generic *i-adjunction* functor
//!     (adjoin a root of `x²+1`) over any backend, not a concrete world.
//!   * [`Laurent`] is `Laurent<S, K>` — a generic *t-adjunction* functor (adjoin a
//!     transcendental `t` with a valuation), the formal Laurent field `S((t))`.
//!     Applied to a finite field it fills the **equal-characteristic local** cell
//!     (`F_q((t))`, the char-`p` mirror of `Qp`); its ring of integers is `F_q[[t]]`.
//!   * [`Ramified`] is `Ramified<S, E>` — a generic *ramified* `π`-adjunction
//!     functor (adjoin a root of the Eisenstein polynomial `xᴱ − ϖ`) over a
//!     [`Valued`] base. It fills the **ramified** local cell: `Q_p(p^{1/E})` over
//!     `Qp`, the ramified twin of the unramified `Qq`. The valuation datum it
//!     needs from the base is abstracted by the [`Valued`] trait.
//!   * [`Gauss`] is `Gauss<S>` — a generic *t-adjunction* with the **Gauss
//!     valuation** over a [`Valued`] base, the rational function field `S(t)` with
//!     `v(t) = 0`. The residue-extending twin of `Laurent` (residue field `k(t̄)`,
//!     value group unchanged); the fourth, last corner of the functor square.
//!
//! And [`ordinal`](big::ordinal)'s nimbers are the **char-2 mirror of the
//! surreals** — the transfinite "big" number in characteristic 2 — so they sit
//! in [`big`] alongside `Surreal`/`Omnific`, not with the finite nim-field.
//!
//! The characteristic trichotomy that organises [`crate::forms`] cuts *across*
//! this table (char 0 in `exact`/`big`/`small`, char 2 in `nimber`/`ordinal`, odd
//! and even in `finite_field`); the two pillars are complementary views of the
//! same backends.

pub mod analytic;
pub mod big;
pub mod exact;
pub mod exactness;
pub mod extension;
pub mod finite_field;
pub mod functor;
pub mod global;
pub mod integrality;
pub mod poly;
pub mod residue;
pub mod small;
pub mod valued;

pub use analytic::*;
pub use big::*;
pub use exact::*;
pub use exactness::*;
pub use extension::*;
pub use finite_field::*;
pub use functor::*;
pub use global::*;
pub use integrality::*;
pub use poly::*;
pub use residue::*;
pub use small::*;
pub use valued::*;

use std::fmt::Debug;
use std::ops::{Add, Mul, Neg, Sub};

/// Generate the owned-value operators `+`, `-` (binary and unary), and `*` for a
/// [`Scalar`] backend by forwarding to its trait methods, so downstream code can
/// write `a + b`, `a * b`, `-a` instead of `a.add(&b)`, `a.mul(&b)`, `a.neg()`.
///
/// Deliberately *not* a [`Scalar`] supertrait bound: these are concrete-type
/// conveniences for callers (`Surreal + Surreal`, `-nimber`), so generic engine
/// code over `S: Scalar` keeps resolving `.add(&x)` / `.mul(&x)` to the `&self`
/// trait methods — operators-on-`S` would shadow them at owned-receiver sites and
/// force clones the borrow-based engine avoids. Division stays a method
/// ([`Scalar::inv`] is partial — `Div` would have to panic), and the
/// by-reference operator forms are omitted for the same reason.
///
/// The generic-backend form takes its generic clause in brackets, e.g.
/// `impl_scalar_ops!([const P: u128] Fp<P>)` or `impl_scalar_ops!([S: Scalar]
/// Surcomplex<S>)`; the bare form `impl_scalar_ops!(Rational)` is for the
/// concrete ones. (Brackets, not `<…>`, so the matcher stays unambiguous.)
macro_rules! impl_scalar_ops {
    ([$($gen:tt)*] $ty:ty) => {
        impl<$($gen)*> Add for $ty {
            type Output = $ty;
            #[inline]
            fn add(self, rhs: $ty) -> $ty { <$ty as $crate::scalar::Scalar>::add(&self, &rhs) }
        }
        impl<$($gen)*> Sub for $ty {
            type Output = $ty;
            #[inline]
            fn sub(self, rhs: $ty) -> $ty { <$ty as $crate::scalar::Scalar>::sub(&self, &rhs) }
        }
        impl<$($gen)*> Mul for $ty {
            type Output = $ty;
            #[inline]
            fn mul(self, rhs: $ty) -> $ty { <$ty as $crate::scalar::Scalar>::mul(&self, &rhs) }
        }
        impl<$($gen)*> Neg for $ty {
            type Output = $ty;
            #[inline]
            fn neg(self) -> $ty { <$ty as $crate::scalar::Scalar>::neg(&self) }
        }
    };
    ($ty:ty) => {
        impl Add for $ty {
            type Output = $ty;
            #[inline]
            fn add(self, rhs: $ty) -> $ty { <$ty as $crate::scalar::Scalar>::add(&self, &rhs) }
        }
        impl Sub for $ty {
            type Output = $ty;
            #[inline]
            fn sub(self, rhs: $ty) -> $ty { <$ty as $crate::scalar::Scalar>::sub(&self, &rhs) }
        }
        impl Mul for $ty {
            type Output = $ty;
            #[inline]
            fn mul(self, rhs: $ty) -> $ty { <$ty as $crate::scalar::Scalar>::mul(&self, &rhs) }
        }
        impl Neg for $ty {
            type Output = $ty;
            #[inline]
            fn neg(self) -> $ty { <$ty as $crate::scalar::Scalar>::neg(&self) }
        }
    };
}

pub trait Scalar: Clone + PartialEq + Debug {
    fn zero() -> Self;
    fn one() -> Self;
    fn add(&self, rhs: &Self) -> Self;
    fn neg(&self) -> Self;
    fn mul(&self, rhs: &Self) -> Self;

    /// Ring characteristic: 0 for characteristic-0 domains, a positive additive
    /// order of `1` for finite fields and finite quotient rings (`Z/p^k`,
    /// truncated Witt vectors, etc.). The engine itself gets signs from
    /// [`Scalar::neg`]; callers that care about characteristic must distinguish
    /// fields from local rings separately.
    fn characteristic() -> u128;

    /// Multiplicative inverse, or `None` if not invertible (zero) or not
    /// finitely representable in this backend (e.g. a non-monomial surreal,
    /// whose inverse is an infinite Hahn series).
    fn inv(&self) -> Option<Self>;

    fn is_zero(&self) -> bool {
        *self == Self::zero()
    }

    fn sub(&self, rhs: &Self) -> Self {
        self.add(&rhs.neg())
    }
}

// The operator manifest: every backend gets `+ - *` and unary `-` forwarded to
// its `Scalar` methods, so the whole table reads uniformly (`a + b`, `-a`). The
// const-/type-generic backends carry their generic clause; the concrete ones
// don't. (See [`impl_scalar_ops`].)
impl_scalar_ops!(Rational);
impl_scalar_ops!(Integer);
impl_scalar_ops!(Surreal);
impl_scalar_ops!(Omnific);
impl_scalar_ops!(Nimber);
impl_scalar_ops!([const P: u128] Fp<P>);
impl_scalar_ops!([const P: u128, const N: usize] Fpn<P, N>);
impl_scalar_ops!([const P: u128, const N: usize, const F: usize] WittVec<P, N, F>);
impl_scalar_ops!([const P: u128, const K: u128] Qp<P, K>);
impl_scalar_ops!([const P: u128, const K: u128] Zp<P, K>);
impl_scalar_ops!([const P: u128, const N: usize, const F: usize] Qq<P, N, F>);
impl_scalar_ops!([S: Scalar] Surcomplex<S>);
impl_scalar_ops!([S: Scalar, const K: usize] Laurent<S, K>);
impl_scalar_ops!([S: Valued, const E: usize] Ramified<S, E>);
impl_scalar_ops!([S: Valued] Gauss<S>);
impl_scalar_ops!(Adele);
impl_scalar_ops!([S: Scalar] RationalFunction<S>);
impl_scalar_ops!([S: Scalar] Poly<S>);

#[cfg(test)]
mod ops_tests {
    use super::*;

    /// Operators must agree with the `Scalar` trait methods they forward to —
    /// over a char-0 field and a char-2 one (where `-a = a`).
    #[test]
    fn operators_match_trait_methods() {
        let (a, b) = (Rational::new(2, 3), Rational::new(1, 6));
        assert_eq!(a.clone() + b.clone(), Scalar::add(&a, &b));
        assert_eq!(a.clone() - b.clone(), Scalar::sub(&a, &b));
        assert_eq!(a.clone() * b.clone(), Scalar::mul(&a, &b));
        assert_eq!(-a.clone(), Scalar::neg(&a));
        assert_eq!(a.clone() - a.clone(), Rational::zero());

        // char 2: `+` is XOR, `*` the nim product, and unary `-` is identity.
        let (x, y) = (Nimber(6), Nimber(3));
        assert_eq!(x + y, Scalar::add(&x, &y));
        assert_eq!(x * y, Scalar::mul(&x, &y));
        assert_eq!(-x, x);
    }
}
