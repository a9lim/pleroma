//! The tropical (min-plus / max-plus) semiring — the algebraic structure that
//! combinatorial-game **thermography** already computes unnamed.
//!
//! A *semiring* keeps a ring's two operations and distributivity but **drops the
//! additive inverse**: tropical ⊕ is idempotent (`a ⊕ a = a`), so there is no
//! `−a`. That is exactly why [`Tropical`] is **not** a [`Scalar`] (and never
//! enters `clifford/`): a Clifford algebra needs a commutative *ring* of scalars,
//! and an idempotent ⊕ has no inverse to give one. The boundary mirrors the one
//! the games pillar already draws — games under disjunctive sum are an abelian
//! *group*, not a ring, so the Clifford story lives only on the field-like cores.
//! [`Semiring`] is therefore a **sibling trait** (like [`Valued`](crate::scalar::Valued)
//! /[`ResidueField`](crate::scalar::ResidueField)), not a `Scalar` supertrait; its
//! impl delegates to inherent methods of the same name.
//!
//! # Two dual conventions
//!
//! Tropical arithmetic comes in two mirror-image flavours, and thermography uses
//! **both at once** (the two scaffold walls live in dual semirings — see
//! [`crate::games::tropical_thermography`]):
//!
//! | convention | ⊕ (`add`) | ⊕-identity `zero` | ⊗ (`mul`) | ⊗-identity `one` |
//! |---|---|---|---|---|
//! | [`MaxPlus`] | `max`     | `−∞`              | `+`       | `0`              |
//! | [`MinPlus`] | `min`     | `+∞`              | `+`       | `0`              |
//!
//! The convention is a compile-time marker, so [`Tropical<MaxPlus>`] and
//! [`Tropical<MinPlus>`] are **distinct types** — the type system forbids mixing
//! the dual walls — that share one implementation body (the `Surcomplex<S>` /
//! `Laurent<S, K>` move).

use crate::scalar::{Rational, Scalar};
use std::cmp::Ordering;
use std::fmt::{self, Debug};
use std::marker::PhantomData;

/// A commutative **semiring**: two associative/commutative operations `⊕`
/// ([`add`](Semiring::add)) and `⊗` ([`mul`](Semiring::mul)) with identities,
/// `⊗` distributing over `⊕`, and `0` (`⊕`-identity) absorbing under `⊗`.
///
/// Deliberately **not** a [`Scalar`] supertrait: a
/// semiring need not have additive inverses (the tropical one does not), so it
/// cannot be a ring. Implementors should make the inherent methods of the same
/// name shadow these (inherent-shadows-trait ⇒ the trait bodies delegate rather
/// than recurse), exactly as the other sibling traits do.
pub trait Semiring: Clone + PartialEq + Debug {
    /// The `⊕`-identity (and `⊗`-absorbing element).
    fn zero() -> Self;
    /// The `⊗`-identity.
    fn one() -> Self;
    /// `⊕` — the additive operation (idempotent in the tropical case).
    fn add(&self, rhs: &Self) -> Self;
    /// `⊗` — the multiplicative operation.
    fn mul(&self, rhs: &Self) -> Self;
    /// Whether this is the `⊕`-identity.
    fn is_zero(&self) -> bool {
        *self == Self::zero()
    }
}

mod sealed {
    pub trait Sealed {}
}

/// Which way a tropical semiring rounds. `MAX = true` is `(max, +)` (⊕-identity
/// `−∞`); `MAX = false` is `(min, +)` (⊕-identity `+∞`). Sealed: the only two
/// inhabitants are [`MaxPlus`] and [`MinPlus`].
pub trait TropicalConvention: sealed::Sealed + Clone + PartialEq + Debug + 'static {
    /// `true` for `(max, +)`, `false` for `(min, +)`.
    const MAX: bool;
    /// How the ⊕-identity (the infinity that anchors this convention) prints.
    const ZERO_GLYPH: &'static str;
}

/// The `(max, +)` convention: ⊕ = `max`, ⊕-identity `−∞`.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Hash)]
pub struct MaxPlus;

/// The `(min, +)` convention: ⊕ = `min`, ⊕-identity `+∞`.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Hash)]
pub struct MinPlus;

impl sealed::Sealed for MaxPlus {}
impl sealed::Sealed for MinPlus {}

impl TropicalConvention for MaxPlus {
    const MAX: bool = true;
    const ZERO_GLYPH: &'static str = "−∞";
}
impl TropicalConvention for MinPlus {
    const MAX: bool = false;
    const ZERO_GLYPH: &'static str = "∞";
}

/// The underlying value: a finite rational, or the convention's infinity (the
/// ⊕-identity). Private — callers use [`Tropical::value`]/[`Tropical::is_infinity`].
#[derive(Clone, PartialEq)]
enum TropVal {
    Finite(Rational),
    Infinity,
}

/// A tropical number in convention `C`. `Tropical<MaxPlus>` lives in `(max, +)`,
/// `Tropical<MinPlus>` in `(min, +)`; they are distinct, non-interoperating types.
///
/// The unused `C` is carried by a [`PhantomData`] (a bare enum with an unused
/// type parameter would not compile).
#[derive(Clone, PartialEq)]
pub struct Tropical<C: TropicalConvention = MaxPlus> {
    inner: TropVal,
    _c: PhantomData<C>,
}

impl<C: TropicalConvention> Tropical<C> {
    /// The finite tropical number `r`.
    pub fn finite(r: Rational) -> Self {
        Tropical {
            inner: TropVal::Finite(r),
            _c: PhantomData,
        }
    }

    /// The convention's infinity — the ⊕-identity (`−∞` for max-plus, `+∞` for
    /// min-plus) and ⊗-absorbing element.
    pub fn infinity() -> Self {
        Tropical {
            inner: TropVal::Infinity,
            _c: PhantomData,
        }
    }

    /// The finite tropical integer `n`.
    pub fn int(n: i128) -> Self {
        Self::finite(Rational::int(n))
    }

    /// The finite value as a rational, or `None` at infinity.
    pub fn value(&self) -> Option<Rational> {
        match &self.inner {
            TropVal::Finite(r) => Some(r.clone()),
            TropVal::Infinity => None,
        }
    }

    /// Whether this is the convention's infinity (the ⊕-identity).
    pub fn is_infinity(&self) -> bool {
        matches!(self.inner, TropVal::Infinity)
    }

    // --- the semiring operations (inherent; the `Semiring` impl delegates here) ---

    /// The ⊕-identity: `−∞` (max-plus) or `+∞` (min-plus).
    pub fn zero() -> Self {
        Self::infinity()
    }

    /// The ⊗-identity: the finite `0`.
    pub fn one() -> Self {
        Self::int(0)
    }

    /// Tropical ⊕ — `max` (max-plus) or `min` (min-plus). Infinity is the
    /// identity on both sides.
    pub fn add(&self, rhs: &Self) -> Self {
        match (&self.inner, &rhs.inner) {
            // Infinity is the ⊕-identity: it loses to any finite value.
            (TropVal::Infinity, _) => rhs.clone(),
            (_, TropVal::Infinity) => self.clone(),
            (TropVal::Finite(a), TropVal::Finite(b)) => {
                let keep_self = if C::MAX {
                    a.cmp(b) != Ordering::Less // a ≥ b ⇒ keep a (the max)
                } else {
                    a.cmp(b) != Ordering::Greater // a ≤ b ⇒ keep a (the min)
                };
                if keep_self {
                    self.clone()
                } else {
                    rhs.clone()
                }
            }
        }
    }

    /// Tropical ⊗ — ordinary `+` of the values; infinity absorbs.
    pub fn mul(&self, rhs: &Self) -> Self {
        match (&self.inner, &rhs.inner) {
            (TropVal::Infinity, _) | (_, TropVal::Infinity) => Self::infinity(),
            (TropVal::Finite(a), TropVal::Finite(b)) => Self::finite(a.add(b)),
        }
    }

    /// Whether this is the ⊕-identity (infinity).
    pub fn is_zero(&self) -> bool {
        self.is_infinity()
    }

    fn fmt_inner(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self.inner {
            TropVal::Infinity => f.write_str(C::ZERO_GLYPH),
            // `Rational` is `Debug`-only (no `Display`).
            TropVal::Finite(r) => write!(f, "{r:?}"),
        }
    }
}

impl<C: TropicalConvention> fmt::Display for Tropical<C> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.fmt_inner(f)
    }
}

impl<C: TropicalConvention> fmt::Debug for Tropical<C> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.fmt_inner(f)
    }
}

impl<C: TropicalConvention> Semiring for Tropical<C> {
    // Each delegates to the inherent method of the same name (inherent shadows
    // trait in path position, so this is not a recursive call).
    fn zero() -> Self {
        Tropical::zero()
    }
    fn one() -> Self {
        Tropical::one()
    }
    fn add(&self, rhs: &Self) -> Self {
        Tropical::add(self, rhs)
    }
    fn mul(&self, rhs: &Self) -> Self {
        Tropical::mul(self, rhs)
    }
    fn is_zero(&self) -> bool {
        Tropical::is_zero(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rat(n: i128, d: i128) -> Rational {
        Rational::new(n, d)
    }

    #[test]
    fn identities() {
        // 0 (∞) is the ⊕-identity; 1 (finite 0) is the ⊗-identity.
        let five = Tropical::<MaxPlus>::int(5);
        assert_eq!(five.add(&Tropical::zero()), five);
        assert_eq!(Tropical::<MaxPlus>::zero().add(&five), five);
        assert_eq!(five.mul(&Tropical::one()), five);
        assert_eq!(Tropical::<MaxPlus>::one().mul(&five), five);

        let m = Tropical::<MinPlus>::int(5);
        assert_eq!(m.add(&Tropical::zero()), m);
        assert_eq!(m.mul(&Tropical::one()), m);
    }

    #[test]
    fn zero_absorbs_under_otimes() {
        // ∞ ⊗ a = ∞ in both conventions.
        let a = Tropical::<MaxPlus>::int(7);
        assert!(Tropical::<MaxPlus>::zero().mul(&a).is_zero());
        assert!(a.mul(&Tropical::<MaxPlus>::zero()).is_zero());

        let b = Tropical::<MinPlus>::int(7);
        assert!(Tropical::<MinPlus>::zero().mul(&b).is_zero());
    }

    #[test]
    fn oplus_is_idempotent() {
        let a = Tropical::<MaxPlus>::finite(rat(3, 2));
        assert_eq!(a.add(&a), a);
        let b = Tropical::<MinPlus>::finite(rat(3, 2));
        assert_eq!(b.add(&b), b);
        // ∞ ⊕ ∞ = ∞
        assert!(Tropical::<MaxPlus>::zero().add(&Tropical::zero()).is_zero());
    }

    #[test]
    fn max_vs_min_hand_checked() {
        let (two_x, five_x) = (Tropical::<MaxPlus>::int(2), Tropical::<MaxPlus>::int(5));
        let (two_n, five_n) = (Tropical::<MinPlus>::int(2), Tropical::<MinPlus>::int(5));
        // 2 ⊕ 5 = 5 (max) vs 2 (min)
        assert_eq!(two_x.add(&five_x), five_x);
        assert_eq!(two_n.add(&five_n), two_n);
        // 2 ⊗ 5 = 7 in both
        assert_eq!(two_x.mul(&five_x), Tropical::<MaxPlus>::int(7));
        assert_eq!(two_n.mul(&five_n), Tropical::<MinPlus>::int(7));
    }

    #[test]
    fn distributivity() {
        let a = Tropical::<MaxPlus>::int(2);
        let b = Tropical::<MaxPlus>::int(5);
        let c = Tropical::<MaxPlus>::int(3);
        // a ⊗ (b ⊕ c) = (a ⊗ b) ⊕ (a ⊗ c)
        assert_eq!(a.mul(&b.add(&c)), a.mul(&b).add(&a.mul(&c)));
        // with an infinity in the mix
        let inf = Tropical::<MaxPlus>::zero();
        assert_eq!(a.mul(&b.add(&inf)), a.mul(&b).add(&a.mul(&inf)));
    }

    #[test]
    fn display_smoke() {
        assert_eq!(format!("{}", Tropical::<MaxPlus>::int(3)), "3");
        assert_eq!(format!("{}", Tropical::<MaxPlus>::finite(rat(3, 2))), "3/2");
        assert_eq!(format!("{}", Tropical::<MaxPlus>::infinity()), "−∞");
        assert_eq!(format!("{}", Tropical::<MinPlus>::infinity()), "∞");
        assert_eq!(format!("{:?}", Tropical::<MinPlus>::int(-4)), "-4");
    }

    #[test]
    fn value_and_infinity_accessors() {
        assert_eq!(Tropical::<MaxPlus>::int(9).value(), Some(Rational::int(9)));
        assert_eq!(Tropical::<MaxPlus>::infinity().value(), None);
        assert!(Tropical::<MinPlus>::infinity().is_infinity());
        assert!(!Tropical::<MinPlus>::int(0).is_infinity());
    }
}
