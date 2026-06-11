//! Frobenius/Galois automorphisms as grade-1 linear maps.
//!
//! [`CyclicGaloisExtension`] already gives
//! a basis and a generator `sigma` of the cyclic Galois group. This module is the
//! Clifford-side bridge: express `sigma` in that basis, feed the resulting matrix
//! to [`LinearMap`], and then the existing outermorphism machinery computes
//! determinant, characteristic polynomial, and exterior-power traces.

use crate::clifford::LinearMap;
use crate::scalar::{CyclicGaloisExtension, FieldExtension, Fp, Fpn, Nimber};

/// A cyclic Galois extension whose distinguished basis has a concrete coordinate
/// extractor. This is deliberately narrower than `CyclicGaloisExtension`: trace
/// forms only need a basis and `sigma`, while a matrix needs coordinates.
pub trait CoordinateCyclicGaloisExtension: CyclicGaloisExtension {
    /// Coordinates of `x` in [`CyclicGaloisExtension::basis`], as base-field
    /// scalars.
    fn coordinates(x: &Self) -> Vec<Self::Base>;
}

impl<const P: u128, const N: usize> CoordinateCyclicGaloisExtension for Fpn<P, N> {
    fn coordinates(x: &Self) -> Vec<Self::Base> {
        x.coeffs().iter().map(|&c| Fp::<P>::from_u128(c)).collect()
    }
}

impl CoordinateCyclicGaloisExtension for Nimber {
    fn coordinates(x: &Self) -> Vec<Self::Base> {
        (0..Self::extension_degree())
            .map(|i| Fp::<2>::from_u128((x.0 >> i) & 1))
            .collect()
    }
}

/// The base-field linear map of `sigma^k` on `E` in the distinguished basis.
pub fn galois_linear_map<E>(k: usize) -> LinearMap<E::Base>
where
    E: CoordinateCyclicGaloisExtension,
{
    let cols = E::basis()
        .into_iter()
        .map(|e| E::coordinates(&e.sigma_power(k)))
        .collect();
    LinearMap::from_columns(cols)
}

/// The base-field linear map of the Frobenius/Galois generator `sigma`.
pub fn frobenius_linear_map<E>() -> LinearMap<E::Base>
where
    E: CoordinateCyclicGaloisExtension,
{
    galois_linear_map::<E>(1)
}

/// Frobenius on the nimber subfield `F_{2^m} ⊂ F_{2^128}`, using the bit basis
/// `{1,2,...,2^(m-1)}`. This keeps the outermorphism spectral checks small; the
/// full `Nimber` extension has degree 128 and is intentionally too large for
/// exterior-power characteristic-polynomial computation.
pub fn nimber_subfield_frobenius_linear_map(m: usize, k: usize) -> LinearMap<Fp<2>> {
    assert!(
        m.is_power_of_two() && m <= 128,
        "nimber subfields represented by low bits require m a power of two <= 128"
    );
    let cols = (0..m)
        .map(|i| {
            let mut x = Nimber(1u128 << i);
            for _ in 0..k {
                x = x.sigma();
            }
            (0..m).map(|j| Fp::<2>::from_u128((x.0 >> j) & 1)).collect()
        })
        .collect();
    LinearMap::from_columns(cols)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::clifford::{char_poly, determinant, exterior_power_trace, CliffordAlgebra, Metric};
    use crate::scalar::{Fpn, Scalar};

    fn exterior_alg<S: Scalar>(n: usize) -> CliffordAlgebra<S> {
        CliffordAlgebra::new(n, Metric::grassmann(n))
    }

    fn expected_xn_minus_one<S: Scalar>(n: usize) -> Vec<S> {
        let mut out = vec![S::zero(); n + 1];
        out[0] = S::one();
        out[n] = S::one().neg();
        out
    }

    fn check_frobenius_spectrum<S>(alg: &CliffordAlgebra<S>, f: &LinearMap<S>)
    where
        S: Scalar + std::fmt::Debug + PartialEq,
    {
        let n = alg.dim();
        assert_eq!(char_poly(alg, f), expected_xn_minus_one::<S>(n));
        for k in 1..n {
            assert_eq!(exterior_power_trace(alg, f, k), S::zero(), "grade {k}");
        }
        let det = determinant(alg, f);
        let expected_det = if (n + 1) % 2 == 1 {
            S::one().neg()
        } else {
            S::one()
        };
        assert_eq!(det, expected_det);
    }

    #[test]
    fn fpn_frobenius_has_xn_minus_one_char_poly() {
        type F8 = Fpn<2, 3>;
        let f8 = frobenius_linear_map::<F8>();
        assert_eq!(f8.n, 3);
        check_frobenius_spectrum(&exterior_alg::<Fp<2>>(3), &f8);

        type F9 = Fpn<3, 2>;
        let f9 = frobenius_linear_map::<F9>();
        assert_eq!(f9.n, 2);
        check_frobenius_spectrum(&exterior_alg::<Fp<3>>(2), &f9);

        type F27 = Fpn<3, 3>;
        let f27 = frobenius_linear_map::<F27>();
        assert_eq!(f27.n, 3);
        check_frobenius_spectrum(&exterior_alg::<Fp<3>>(3), &f27);
    }

    #[test]
    fn nimber_subfield_frobenius_uses_the_same_outermorphism_oracle() {
        let f16 = nimber_subfield_frobenius_linear_map(4, 1);
        assert_eq!(f16.n, 4);
        check_frobenius_spectrum(&exterior_alg::<Fp<2>>(4), &f16);
    }

    #[test]
    fn frobenius_power_composes_as_expected() {
        type F8 = Fpn<2, 3>;
        let sigma = frobenius_linear_map::<F8>();
        let sigma2 = galois_linear_map::<F8>(2);
        let sigma3 = galois_linear_map::<F8>(3);
        assert_eq!(sigma.compose(&sigma), sigma2);
        assert_eq!(sigma.compose(&sigma2), sigma3);
        assert_eq!(sigma3, LinearMap::identity(3));
    }
}
