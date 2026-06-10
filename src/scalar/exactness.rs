//! Exactness markers for scalar backends.
//!
//! [`Scalar`] is the operational contract the Clifford engine needs: a
//! commutative scalar world with `+`, `-`, `*`, zero, one, characteristic, and a
//! partial inverse. It deliberately does **not** say whether the represented
//! arithmetic is exact. This module names that second axis explicitly:
//!
//! * [`ExactScalar`] — arithmetic is exact in the represented scalar ring. Finite
//!   quotient rings such as [`Zp`] and
//!   [`WittVec`] count: they are exact finite rings, even
//!   though they model a truncated infinite object.
//! * [`PrecisionScalar`] — the backend is a capped-relative precision model.
//!   These worlds are useful at the valuation/forms layer, but they are omitted
//!   from the exact-ring fuzz because addition can lose information across the
//!   retained window.
//!
//! Both traits are markers, not [`Scalar`] supertraits. The engine keeps accepting
//! any `Scalar`; tests and generic helpers that require exact ring laws can opt in
//! to [`ExactScalar`] when that stronger promise matters.

use crate::scalar::{
    Adele, Fp, Fpn, Gauss, Integer, Laurent, Nimber, Omnific, Poly, Qp, Qq, Ramified, Rational,
    RationalFunction, Scalar, Surcomplex, Surreal, Valued, WittVec, Zp,
};

/// A scalar backend whose represented arithmetic is exact.
pub trait ExactScalar: Scalar {}

/// An exact scalar backend that is intended to behave as a field in the represented
/// domain. This is a marker for generic constructions such as `S(t)` that need
/// nonzero leading coefficients to invert; it is intentionally narrower than
/// [`ExactScalar`].
pub trait ExactFieldScalar: ExactScalar {}

/// A scalar backend with capped-relative precision. These are still valid
/// [`Scalar`] implementations for the geometric engine, but their addition is a
/// finite-precision model rather than an exact ring operation on the ideal object.
pub trait PrecisionScalar: Scalar {}

// Exact rings / represented fields.
impl ExactScalar for Integer {}
impl ExactScalar for Rational {}
impl ExactFieldScalar for Rational {}
impl ExactScalar for Surreal {}
impl ExactScalar for Omnific {}
impl ExactScalar for Nimber {}
impl ExactFieldScalar for Nimber {}
impl<const P: u128> ExactScalar for Fp<P> {}
impl<const P: u128> ExactFieldScalar for Fp<P> {}
impl<const P: u128, const N: usize> ExactScalar for Fpn<P, N> {}
impl<const P: u128, const N: usize> ExactFieldScalar for Fpn<P, N> {}
impl<const P: u128, const K: u128> ExactScalar for Zp<P, K> {}
impl<const P: u128, const N: usize, const F: usize> ExactScalar for WittVec<P, N, F> {}
impl<S: ExactScalar> ExactScalar for Surcomplex<S> {}
impl<S: ExactScalar> ExactScalar for Poly<S> {}
impl<S: ExactFieldScalar> ExactScalar for RationalFunction<S> {}
impl<S: ExactFieldScalar> ExactFieldScalar for RationalFunction<S> {}

// Capped-relative precision models.
impl<const P: u128, const K: u128> PrecisionScalar for Qp<P, K> {}
impl<const P: u128, const N: usize, const F: usize> PrecisionScalar for Qq<P, N, F> {}
impl<S: Scalar, const K: usize> PrecisionScalar for Laurent<S, K> {}
impl<S: PrecisionScalar + Valued, const E: usize> PrecisionScalar for Ramified<S, E> {}
impl<S: PrecisionScalar + Valued> PrecisionScalar for Gauss<S> {}
impl PrecisionScalar for Adele {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scalar::{Fp, Laurent, Qp, Rational};

    fn assert_exact<S: ExactScalar>() {}
    fn assert_exact_field<S: ExactFieldScalar>() {}
    fn assert_precision<S: PrecisionScalar>() {}

    #[test]
    fn marker_traits_cover_the_table() {
        assert_exact::<Integer>();
        assert_exact::<Rational>();
        assert_exact_field::<Rational>();
        assert_exact::<Surreal>();
        assert_exact::<Omnific>();
        assert_exact::<Nimber>();
        assert_exact_field::<Nimber>();
        assert_exact::<Fp<7>>();
        assert_exact_field::<Fp<7>>();
        assert_exact::<Fpn<3, 2>>();
        assert_exact_field::<Fpn<3, 2>>();
        assert_exact::<Zp<5, 4>>();
        assert_exact::<WittVec<3, 3, 2>>();
        assert_exact::<Surcomplex<Surreal>>();
        assert_exact::<Poly<Rational>>();
        assert_exact::<RationalFunction<Fp<7>>>();
        assert_exact_field::<RationalFunction<Fp<7>>>();

        assert_precision::<Qp<5, 4>>();
        assert_precision::<Qq<3, 3, 2>>();
        assert_precision::<Laurent<Fp<5>, 6>>();
        assert_precision::<Ramified<Qp<3, 6>, 2>>();
        assert_precision::<Gauss<Qp<3, 6>>>();
        assert_precision::<Gauss<Laurent<Fp<5>, 6>>>();
        assert_precision::<Adele>();
    }
}
