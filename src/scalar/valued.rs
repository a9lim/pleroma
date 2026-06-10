//! The [`Valued`] trait: a scalar carrying a discrete valuation and a canonical
//! uniformizer.
//!
//! Every backend in the **non-Archimedean local** part of the "any number" table
//! already exposes an *inherent* `valuation()` and a way to name its prime element
//! ([`Qp`]/[`Qq`] via `from_p_power(1)`,
//! [`Laurent`] via `t()`). This trait promotes that shared
//! shape to the type system so the [`Ramified`](crate::scalar::Ramified)
//! ramified-extension functor can fold a *generic* base valuation — it adjoins a
//! uniformizer `π` with `πᴱ = ϖ`, and `ϖ = S::uniformizer()` is exactly the datum
//! it needs from the base field.
//!
//! Deliberately **not** a [`Scalar`] supertrait (same reasoning as the operator
//! manifest): only the discretely-valued local fields are `Valued`. The exact
//! Archimedean worlds (`Rational`, `Surreal`) carry no canonical uniformizer and
//! are intentionally left out. The rings of integers (`Zp`, `WittVec`) are also
//! left out: a `Ramified` base must be a *field* so its `inv` is total on
//! nonzero.

use crate::scalar::{Laurent, Qp, Qq, Scalar};

/// A scalar with a discrete valuation `v : K → ℤ ∪ {∞}` and a canonical
/// uniformizer `ϖ` (the valuation-`1` element). The valuation here is the same
/// one each backend exposes inherently; this trait just makes it generic.
pub trait Valued: Scalar {
    /// The valuation of this element, or `None` for zero (valuation `+∞`).
    fn valuation(&self) -> Option<i128>;

    /// The canonical uniformizer `ϖ` — the prime element of valuation `1`
    /// (`p` for `Qp`/`Qq`, `t` for `Laurent`).
    fn uniformizer() -> Self;
}

impl<const P: u128, const K: u128> Valued for Qp<P, K> {
    fn valuation(&self) -> Option<i128> {
        // Inherent `Qp::valuation` shadows this trait method in method-call
        // position, so this delegates rather than recursing.
        Qp::valuation(self)
    }
    fn uniformizer() -> Self {
        Qp::from_p_power(1)
    }
}

impl<const P: u128, const N: usize, const F: usize> Valued for Qq<P, N, F> {
    fn valuation(&self) -> Option<i128> {
        Qq::valuation(self)
    }
    fn uniformizer() -> Self {
        Qq::from_p_power(1)
    }
}

impl<S: Scalar, const K: usize> Valued for Laurent<S, K> {
    fn valuation(&self) -> Option<i128> {
        Laurent::valuation(self)
    }
    fn uniformizer() -> Self {
        Laurent::t()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scalar::{Fp, Rational};

    #[test]
    fn uniformizers_have_valuation_one() {
        assert_eq!(Qp::<5, 4>::uniformizer().valuation(), Some(1));
        assert_eq!(Qq::<3, 4, 2>::uniformizer().valuation(), Some(1));
        assert_eq!(Laurent::<Rational, 6>::uniformizer().valuation(), Some(1));
        assert_eq!(Laurent::<Fp<7>, 6>::uniformizer().valuation(), Some(1));
    }

    #[test]
    fn zero_valuation_is_none() {
        assert_eq!(<Qp<5, 4> as Valued>::valuation(&Qp::zero()), None);
        assert_eq!(
            <Laurent<Rational, 6> as Valued>::valuation(&Laurent::zero()),
            None
        );
    }

    #[test]
    fn trait_valuation_matches_inherent() {
        let x = Qp::<5, 4>::from_i128(50); // 2·5²  ⇒ valuation 2
        assert_eq!(<Qp<5, 4> as Valued>::valuation(&x), x.valuation());
        assert_eq!(<Qp<5, 4> as Valued>::valuation(&x), Some(2));
    }
}
