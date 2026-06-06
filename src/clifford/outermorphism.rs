//! Outermorphisms: the canonical lift of a grade-1 linear map to the whole
//! algebra, and the determinant that falls out of it.
//!
//! A linear map `f: V → V` extends uniquely to an algebra endomorphism of the
//! exterior structure by `f(a ∧ b) = f(a) ∧ f(b)` — an *outermorphism*. It acts
//! blade by blade: send each generator `e_i` to `f(e_i)` and wedge the images in
//! order. Because the engine's `wedge` already carries the reordering sign
//! through `S::neg()`, the lift is character-faithful for free — in particular
//! the **determinant** computed here is the ordinary determinant in char 0 and
//! the char-2 determinant (= permanent) over the nimbers, with no sign hardcoded.
//!
//! The determinant is read off the top grade: `f(I) = det(f)·I` for the unit
//! pseudoscalar `I` (Grassmann's original definition of the determinant). This
//! is a structurally independent computation from cofactor expansion, so it
//! doubles as a check on the engine's `wedge`.

use crate::clifford::{bits, CliffordAlgebra, Multivector};
use crate::scalar::Scalar;

/// A linear map `V → V` on grade 1, stored column-major: `cols[i]` is the image
/// `f(e_i)` as a length-`n` coefficient vector over `e_0..e_{n-1}` (so
/// `cols[i][j]` is the coefficient of `e_j` in `f(e_i)`, i.e. the matrix entry
/// `M[j][i]`).
#[derive(Clone, Debug, PartialEq)]
pub struct LinearMap<S: Scalar> {
    pub n: usize,
    pub cols: Vec<Vec<S>>,
}

impl<S: Scalar> LinearMap<S> {
    /// Build from columns `cols[i] = f(e_i)`; panics if not square `n×n`.
    pub fn from_columns(cols: Vec<Vec<S>>) -> Self {
        let n = cols.len();
        assert!(
            cols.iter().all(|c| c.len() == n),
            "LinearMap must be square: each column has length n"
        );
        LinearMap { n, cols }
    }

    /// The identity map on `n` generators.
    pub fn identity(n: usize) -> Self {
        let cols = (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| if i == j { S::one() } else { S::zero() })
                    .collect()
            })
            .collect();
        LinearMap { n, cols }
    }

    /// `f(e_i)` as a grade-1 multivector in `alg`.
    pub fn image(&self, alg: &CliffordAlgebra<S>, i: usize) -> Multivector<S> {
        let mut out = alg.zero();
        for (j, c) in self.cols[i].iter().enumerate() {
            if !c.is_zero() {
                out = alg.add(&out, &alg.scalar_mul(c, &alg.gen(j)));
            }
        }
        out
    }

    /// The composite `self ∘ inner` (apply `inner`, then `self`): the ordinary
    /// matrix product `M_self · M_inner`.
    pub fn compose(&self, inner: &LinearMap<S>) -> LinearMap<S> {
        assert_eq!(self.n, inner.n, "dimension mismatch in compose");
        let n = self.n;
        // cols_{f∘g}[i][j] = Σ_k cols_f[k][j] · cols_g[i][k]
        let cols = (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| {
                        let mut acc = S::zero();
                        for k in 0..n {
                            acc = acc.add(&self.cols[k][j].mul(&inner.cols[i][k]));
                        }
                        acc
                    })
                    .collect()
            })
            .collect();
        LinearMap { n, cols }
    }
}

/// Apply the outermorphism of `f` to a multivector: every blade `e_{i_0..i_k}`
/// maps to `f(e_{i_0}) ∧ … ∧ f(e_{i_k})` (grade 0 passes through, since the empty
/// wedge is the unit scalar).
pub fn apply_outermorphism<S: Scalar>(
    alg: &CliffordAlgebra<S>,
    f: &LinearMap<S>,
    mv: &Multivector<S>,
) -> Multivector<S> {
    debug_assert_eq!(f.n, alg.dim, "LinearMap dimension must match the algebra");
    let mut out = alg.zero();
    for (&mask, coeff) in &mv.terms {
        // Fold f(e_i) over the set bits in ascending order, starting at 1.
        let mut acc = alg.scalar(S::one());
        for i in bits(mask) {
            acc = alg.wedge(&acc, &f.image(alg, i));
        }
        out = alg.add(&out, &alg.scalar_mul(coeff, &acc));
    }
    out
}

/// The determinant of `f`: the scalar by which its outermorphism scales the unit
/// pseudoscalar, `f(I) = det(f)·I`.
pub fn determinant<S: Scalar>(alg: &CliffordAlgebra<S>, f: &LinearMap<S>) -> S {
    let pseudo = alg.pseudoscalar();
    let image = apply_outermorphism(alg, f, &pseudo);
    // Pseudoscalar mask = the single key of `pseudo`.
    let mask = *pseudo.terms.keys().next().expect("pseudoscalar is nonzero");
    image.terms.get(&mask).cloned().unwrap_or_else(S::zero)
}

/// The inverse outermorphism, if `f` is invertible over `S`: returns the
/// `LinearMap` of `f⁻¹` (Gauss–Jordan over the scalar). `None` if any pivot is
/// not invertible in the backend (e.g. a non-monomial surreal), which includes
/// the singular case.
pub fn inverse_outermorphism<S: Scalar>(f: &LinearMap<S>) -> Option<LinearMap<S>> {
    let n = f.n;
    // Row-major working matrix `m[r][c] = M[r][c] = cols[c][r]`, augmented with I.
    let mut m: Vec<Vec<S>> = (0..n)
        .map(|r| (0..n).map(|c| f.cols[c][r].clone()).collect())
        .collect();
    let mut inv = LinearMap::<S>::identity(n).matrix_rows();
    for col in 0..n {
        // Find a pivot row at/after `col` with an invertible entry in `col`.
        let piv = (col..n).find(|&r| m[r][col].inv().is_some())?;
        m.swap(col, piv);
        inv.swap(col, piv);
        let pinv = m[col][col].inv()?;
        // Scale the pivot row to make the pivot 1.
        for c in 0..n {
            m[col][c] = m[col][c].mul(&pinv);
            inv[col][c] = inv[col][c].mul(&pinv);
        }
        // Eliminate the pivot column from every other row.
        for r in 0..n {
            if r == col {
                continue;
            }
            let factor = m[r][col].clone();
            if factor.is_zero() {
                continue;
            }
            for c in 0..n {
                m[r][c] = m[r][c].sub(&factor.mul(&m[col][c]));
                inv[r][c] = inv[r][c].sub(&factor.mul(&inv[col][c]));
            }
        }
    }
    // inv is now M⁻¹ in row-major form; convert back to columns.
    let cols = (0..n)
        .map(|i| (0..n).map(|j| inv[j][i].clone()).collect())
        .collect();
    Some(LinearMap { n, cols })
}

impl<S: Scalar> LinearMap<S> {
    /// Row-major matrix `rows[r][c] = M[r][c] = cols[c][r]`.
    fn matrix_rows(&self) -> Vec<Vec<S>> {
        (0..self.n)
            .map(|r| (0..self.n).map(|c| self.cols[c][r].clone()).collect())
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::clifford::Metric;
    use crate::scalar::Nimber;
    use crate::scalar::Rational;

    fn r(n: i128) -> Rational {
        Rational::int(n)
    }

    fn euclid(n: usize) -> CliffordAlgebra<Rational> {
        CliffordAlgebra::new(n, Metric::diagonal(vec![r(1); n]))
    }

    #[test]
    fn determinant_of_identity_is_one() {
        let alg = euclid(3);
        let id = LinearMap::identity(3);
        assert_eq!(determinant(&alg, &id), r(1));
    }

    #[test]
    fn determinant_2x2_matches_hand() {
        // f(e0)=2e0+3e1, f(e1)=e0+4e1 ⇒ M=[[2,1],[3,4]], det = 8−3 = 5.
        let alg = euclid(2);
        let f = LinearMap::from_columns(vec![vec![r(2), r(3)], vec![r(1), r(4)]]);
        assert_eq!(determinant(&alg, &f), r(5));
    }

    #[test]
    fn determinant_3x3_matches_hand() {
        // cols: f(e0)=e0+4e2, f(e1)=3e1, f(e2)=2e0+5e2 ⇒ det = −9.
        let alg = euclid(3);
        let f = LinearMap::from_columns(vec![
            vec![r(1), r(0), r(4)],
            vec![r(0), r(3), r(0)],
            vec![r(2), r(0), r(5)],
        ]);
        assert_eq!(determinant(&alg, &f), r(-9));
    }

    #[test]
    fn determinant_is_multiplicative_rational() {
        let alg = euclid(3);
        let f = LinearMap::from_columns(vec![
            vec![r(1), r(2), r(0)],
            vec![r(0), r(1), r(3)],
            vec![r(2), r(0), r(1)],
        ]);
        let g = LinearMap::from_columns(vec![
            vec![r(2), r(0), r(1)],
            vec![r(1), r(3), r(0)],
            vec![r(0), r(1), r(2)],
        ]);
        let fg = f.compose(&g);
        let lhs = determinant(&alg, &fg);
        let rhs = determinant(&alg, &f).mul(&determinant(&alg, &g));
        assert_eq!(lhs, rhs);
    }

    #[test]
    fn determinant_is_multiplicative_nimber() {
        // The char-2 determinant: neg = id, so this is the permanent — still
        // multiplicative, and computed with no special-casing.
        let alg = CliffordAlgebra::new(3, Metric::diagonal(vec![Nimber(1); 3]));
        let f = LinearMap::from_columns(vec![
            vec![Nimber(1), Nimber(2), Nimber(0)],
            vec![Nimber(3), Nimber(1), Nimber(2)],
            vec![Nimber(0), Nimber(1), Nimber(5)],
        ]);
        let g = LinearMap::from_columns(vec![
            vec![Nimber(2), Nimber(0), Nimber(1)],
            vec![Nimber(1), Nimber(4), Nimber(0)],
            vec![Nimber(0), Nimber(3), Nimber(2)],
        ]);
        let fg = f.compose(&g);
        let lhs = determinant(&alg, &fg);
        let rhs = determinant(&alg, &f).mul(&determinant(&alg, &g));
        assert_eq!(lhs, rhs);
    }

    #[test]
    fn outermorphism_is_an_algebra_map_on_wedge() {
        let alg = euclid(3);
        let f = LinearMap::from_columns(vec![
            vec![r(1), r(2), r(0)],
            vec![r(0), r(1), r(3)],
            vec![r(2), r(0), r(1)],
        ]);
        let e0e1 = alg.wedge(&alg.gen(0), &alg.gen(1));
        let lhs = apply_outermorphism(&alg, &f, &e0e1);
        let rhs = alg.wedge(&f.image(&alg, 0), &f.image(&alg, 1));
        assert_eq!(lhs, rhs);
    }

    #[test]
    fn outermorphism_respects_composition() {
        let alg = euclid(3);
        let f = LinearMap::from_columns(vec![
            vec![r(1), r(2), r(0)],
            vec![r(0), r(1), r(3)],
            vec![r(2), r(0), r(1)],
        ]);
        let g = LinearMap::from_columns(vec![
            vec![r(2), r(0), r(1)],
            vec![r(1), r(3), r(0)],
            vec![r(0), r(1), r(2)],
        ]);
        let fg = f.compose(&g);
        let mv = alg.add(&alg.gen(0), &alg.wedge(&alg.gen(1), &alg.gen(2)));
        let lhs = apply_outermorphism(&alg, &fg, &mv);
        let rhs = apply_outermorphism(&alg, &f, &apply_outermorphism(&alg, &g, &mv));
        assert_eq!(lhs, rhs);
    }

    #[test]
    fn inverse_outermorphism_inverts() {
        let f = LinearMap::from_columns(vec![
            vec![r(1), r(2), r(0)],
            vec![r(0), r(1), r(3)],
            vec![r(2), r(0), r(1)],
        ]);
        let finv = inverse_outermorphism(&f).unwrap();
        let prod = f.compose(&finv);
        assert_eq!(prod, LinearMap::identity(3));
        // det(f⁻¹) = 1/det(f)
        let alg = euclid(3);
        let d = determinant(&alg, &f);
        let dinv = determinant(&alg, &finv);
        assert_eq!(d.mul(&dinv), r(1));
    }

    #[test]
    fn singular_map_has_no_inverse() {
        // Two equal columns ⇒ rank-deficient ⇒ no inverse.
        let f = LinearMap::from_columns(vec![
            vec![r(1), r(1), r(0)],
            vec![r(1), r(1), r(0)],
            vec![r(0), r(0), r(1)],
        ]);
        assert!(inverse_outermorphism(&f).is_none());
        let alg = euclid(3);
        assert_eq!(determinant(&alg, &f), r(0));
    }
}
