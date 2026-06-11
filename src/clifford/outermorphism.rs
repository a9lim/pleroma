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
use crate::linalg::field;
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
    debug_assert_eq!(f.n, alg.dim(), "LinearMap dimension must match the algebra");
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

/// The grade-`k` basis blade masks over `n` generators (the `C(n,k)` subsets),
/// enumerated by Gosper's hack. Exponential in `n` summed over all grades, so
/// the spectral routines below are for modest dimensions.
fn grade_k_masks(n: usize, k: usize) -> Vec<u128> {
    if k == 0 {
        return vec![0];
    }
    if k > n {
        return vec![];
    }
    assert!(n <= u128::BITS as usize, "basis masks fit in u128");
    if k == u128::BITS as usize {
        return vec![u128::MAX];
    }
    let mut out = Vec::new();
    let mut c: u128 = (1u128 << k) - 1;
    let limit = (n < u128::BITS as usize).then(|| 1u128 << n);
    loop {
        out.push(c);
        let u = c & c.wrapping_neg();
        let v = c.checked_add(u);
        match v {
            Some(v) if v != 0 => {
                let next = v + (((v ^ c) / u) >> 2);
                if limit.is_some_and(|lim| next >= lim) {
                    break;
                }
                c = next;
            }
            _ => break,
        }
    }
    out
}

/// The trace of the `k`-th exterior power `Λᵏf` — the `k`-th elementary
/// symmetric function of the eigenvalues, equivalently the sum of the `k×k`
/// principal minors. `Λ⁰f` has trace `1`, `Λ¹f` is the ordinary trace, and
/// `Λⁿf` is the [`determinant`]. Computed straight from the outermorphism:
/// `tr Λᵏf = Σ_{|S|=k} ⟨e_S , f(e_S)⟩`, so it is character-faithful for free.
pub fn exterior_power_trace<S: Scalar>(alg: &CliffordAlgebra<S>, f: &LinearMap<S>, k: usize) -> S {
    debug_assert_eq!(f.n, alg.dim(), "LinearMap dimension must match the algebra");
    let mut acc = S::zero();
    for mask in grade_k_masks(alg.dim(), k) {
        let blade = alg.blade(&bits(mask));
        let img = apply_outermorphism(alg, f, &blade);
        // ⟨e_S , f(e_S)⟩ — the diagonal entry of Λᵏf at this blade.
        if let Some(c) = img.terms.get(&mask) {
            // `blade` may carry a ±1 from ordering; normalise by it.
            let sign = blade.terms.get(&mask).cloned().unwrap_or_else(S::one);
            acc = acc.add(&c.mul(&sign));
        }
    }
    acc
}

/// The ordinary trace of `f` (`= tr Λ¹f = Σᵢ Mᵢᵢ`).
pub fn trace<S: Scalar>(alg: &CliffordAlgebra<S>, f: &LinearMap<S>) -> S {
    exterior_power_trace(alg, f, 1)
}

/// The characteristic polynomial `det(t·I − f)`, returned as coefficients in
/// **descending** degree: `[1, −c₁, c₂, …, (−1)ⁿcₙ]`, where `cₖ = tr Λᵏf`. The
/// leading coefficient is `1` (monic) and the constant term is `(−1)ⁿ det(f)`.
/// Char-faithful — over the nimbers every sign collapses, giving the char-2
/// characteristic polynomial with no special-casing.
pub fn char_poly<S: Scalar>(alg: &CliffordAlgebra<S>, f: &LinearMap<S>) -> Vec<S> {
    let n = alg.dim();
    (0..=n)
        .map(|k| {
            let ck = exterior_power_trace(alg, f, k);
            if k % 2 == 1 {
                ck.neg()
            } else {
                ck
            }
        })
        .collect()
}

/// The inverse outermorphism, if `f` is invertible over `S`: returns the
/// `LinearMap` of `f⁻¹` (Gauss–Jordan over the scalar). `None` if any pivot is
/// not invertible in the backend (e.g. a non-monomial surreal), which includes
/// the singular case.
pub fn inverse_outermorphism<S: Scalar>(f: &LinearMap<S>) -> Option<LinearMap<S>> {
    let n = f.n;
    // Row-major working matrix `m[r][c] = M[r][c] = cols[c][r]`.
    let m: Vec<Vec<S>> = (0..n)
        .map(|r| (0..n).map(|c| f.cols[c][r].clone()).collect())
        .collect();
    let inv = field::inverse_matrix(m)?;
    // inv is now M⁻¹ in row-major form; convert back to columns.
    let cols = (0..n)
        .map(|i| (0..n).map(|j| inv[j][i].clone()).collect())
        .collect();
    Some(LinearMap { n, cols })
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
    fn char_poly_of_identity_is_binomial() {
        // char poly of I_3 is (t−1)³ = t³ − 3t² + 3t − 1.
        let alg = euclid(3);
        let id = LinearMap::identity(3);
        assert_eq!(char_poly(&alg, &id), vec![r(1), r(-3), r(3), r(-1)]);
        assert_eq!(trace(&alg, &id), r(3));
    }

    #[test]
    fn grade_masks_cover_the_full_u128_basis_window() {
        let one_blades = grade_k_masks(128, 1);
        assert_eq!(one_blades.len(), 128);
        assert_eq!(one_blades[0], 1);
        assert_eq!(one_blades[127], 1u128 << 127);
        assert_eq!(grade_k_masks(128, 128), vec![u128::MAX]);
    }

    #[test]
    fn trace_of_identity_at_dim_128_is_128() {
        let alg = euclid(128);
        let id = LinearMap::identity(128);
        assert_eq!(trace(&alg, &id), r(128));
    }

    #[test]
    fn char_poly_matches_trace_and_determinant() {
        // M = [[2,1],[3,4]]: trace 6, det 5, char poly t² − 6t + 5.
        let alg = euclid(2);
        let f = LinearMap::from_columns(vec![vec![r(2), r(3)], vec![r(1), r(4)]]);
        let p = char_poly(&alg, &f);
        assert_eq!(p, vec![r(1), r(-6), r(5)]);
        assert_eq!(trace(&alg, &f), r(6));
        // constant term = (−1)² det = det.
        assert_eq!(*p.last().unwrap(), determinant(&alg, &f));
        // exterior-power traces are the elementary symmetric functions.
        assert_eq!(exterior_power_trace(&alg, &f, 0), r(1));
        assert_eq!(exterior_power_trace(&alg, &f, 1), r(6));
        assert_eq!(exterior_power_trace(&alg, &f, 2), r(5));
    }

    #[test]
    fn char_poly_constant_term_is_signed_determinant_3x3() {
        let alg = euclid(3);
        let f = LinearMap::from_columns(vec![
            vec![r(1), r(0), r(4)],
            vec![r(0), r(3), r(0)],
            vec![r(2), r(0), r(5)],
        ]);
        let p = char_poly(&alg, &f);
        // n=3 ⇒ constant term = (−1)³ det = −det = −(−9) = 9.
        assert_eq!(*p.last().unwrap(), r(9));
        assert_eq!(*p.last().unwrap(), determinant(&alg, &f).neg());
    }

    #[test]
    fn char_poly_is_char_faithful_over_nimbers() {
        // Over the nimbers neg = id, so every coefficient is the bare Λᵏ-trace
        // (the char-2 characteristic polynomial), and trace is the XOR diagonal.
        let alg = CliffordAlgebra::new(2, Metric::diagonal(vec![Nimber(1); 2]));
        let f =
            LinearMap::from_columns(vec![vec![Nimber(2), Nimber(3)], vec![Nimber(1), Nimber(4)]]);
        assert_eq!(trace(&alg, &f), Nimber(2 ^ 4));
        let p = char_poly(&alg, &f);
        assert_eq!(p[0], Nimber(1));
        assert_eq!(*p.last().unwrap(), determinant(&alg, &f)); // (−1)²=1 anyway
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
