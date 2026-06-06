//! Blade analysis: deciding whether a multivector is *decomposable* (a wedge of
//! vectors) and, when it is, factoring it back into them.
//!
//! A grade-`k` multivector `A` is a **blade** iff it equals `v₁ ∧ … ∧ vₖ` for
//! some vectors. The clean intrinsic test is via its **blade subspace**
//!
//!   `⟨A⟩ = { x ∈ V : x ∧ A = 0 }`,
//!
//! the set of vectors "inside" `A`. For a genuine `k`-blade this subspace is
//! exactly `k`-dimensional (and `A` is, up to scale, the wedge of any basis of
//! it); for a non-blade like `e₀∧e₁ + e₂∧e₃` it collapses to `{0}`. So `A` is a
//! blade iff `dim⟨A⟩ = grade(A)`. [`factor_blade`] then returns a basis of `⟨A⟩`
//! rescaled so its wedge is `A` on the nose.
//!
//! Everything is exact and char-faithful (the wedge carries its signs through
//! `S::neg`); the only field operation used is `inv`, for the nullspace solve.

use crate::clifford::{grade, CliffordAlgebra, Multivector};
use crate::scalar::Scalar;
use std::collections::BTreeSet;

/// The common grade of a homogeneous multivector, or `None` if it is zero or of
/// mixed grade.
fn homogeneous_grade<S: Scalar>(a: &Multivector<S>) -> Option<usize> {
    let mut g: Option<usize> = None;
    for &mask in a.terms.keys() {
        let gk = grade(mask);
        match g {
            None => g = Some(gk),
            Some(x) if x == gk => {}
            _ => return None, // mixed grade
        }
    }
    g // None ⇔ zero (no terms)
}

/// A basis of the right nullspace `{ x : M x = 0 }` of a row-major matrix with
/// `ncols` columns, by reduction to RREF over the field.
fn nullspace<S: Scalar>(mut m: Vec<Vec<S>>, ncols: usize) -> Vec<Vec<S>> {
    let nrows = m.len();
    let mut pivot_cols: Vec<usize> = Vec::new();
    let mut row = 0;
    for col in 0..ncols {
        let Some(piv) = (row..nrows).find(|&r| m[r][col].inv().is_some()) else {
            continue;
        };
        m.swap(row, piv);
        let pinv = m[row][col].inv().expect("pivot is invertible");
        for c in 0..ncols {
            m[row][c] = m[row][c].mul(&pinv);
        }
        for r in 0..nrows {
            if r == row {
                continue;
            }
            let f = m[r][col].clone();
            if f.is_zero() {
                continue;
            }
            for c in 0..ncols {
                let sub = f.mul(&m[row][c]);
                m[r][c] = m[r][c].sub(&sub);
            }
        }
        pivot_cols.push(col);
        row += 1;
        if row == nrows {
            break;
        }
    }
    let mut basis = Vec::new();
    for fc in (0..ncols).filter(|c| !pivot_cols.contains(c)) {
        let mut x = vec![S::zero(); ncols];
        x[fc] = S::one();
        for (ri, &pc) in pivot_cols.iter().enumerate() {
            x[pc] = m[ri][fc].neg();
        }
        basis.push(x);
    }
    basis
}

/// A grade-1 multivector from a coefficient vector over `e_0..e_{n-1}`.
fn vector_from_coeffs<S: Scalar>(alg: &CliffordAlgebra<S>, x: &[S]) -> Multivector<S> {
    let mut out = alg.zero();
    for (i, c) in x.iter().enumerate() {
        if !c.is_zero() {
            out = alg.add(&out, &alg.scalar_mul(c, &alg.gen(i)));
        }
    }
    out
}

/// A basis of the **blade subspace** `⟨A⟩ = { x : x ∧ A = 0 }` of a homogeneous
/// multivector `A`, as coefficient vectors over the generators. `None` if `A` is
/// zero or of mixed grade. (A scalar returns the empty basis — its blade
/// subspace is `{0}`.)
pub fn blade_subspace<S: Scalar>(
    alg: &CliffordAlgebra<S>,
    a: &Multivector<S>,
) -> Option<Vec<Vec<S>>> {
    let k = homogeneous_grade(a)?;
    let n = alg.dim;
    if k == 0 {
        return Some(vec![]);
    }
    // Columns of the linear map x ↦ x ∧ A are e_i ∧ A (grade k+1).
    let cols: Vec<Multivector<S>> = (0..n).map(|i| alg.wedge(&alg.gen(i), a)).collect();
    let mut maskset: BTreeSet<u128> = BTreeSet::new();
    for c in &cols {
        maskset.extend(c.terms.keys().copied());
    }
    let masks: Vec<u128> = maskset.into_iter().collect();
    let mat: Vec<Vec<S>> = masks
        .iter()
        .map(|&mask| {
            (0..n)
                .map(|i| cols[i].terms.get(&mask).cloned().unwrap_or_else(S::zero))
                .collect()
        })
        .collect();
    Some(nullspace(mat, n))
}

/// Whether `A` is a **blade** (a decomposable homogeneous multivector). Scalars
/// and vectors always are; `A` of grade `k ≥ 2` is a blade iff `dim⟨A⟩ = k`.
/// Zero and mixed-grade multivectors are not blades.
pub fn is_blade<S: Scalar>(alg: &CliffordAlgebra<S>, a: &Multivector<S>) -> bool {
    match homogeneous_grade(a) {
        None => false,
        Some(0) => true,
        Some(k) => blade_subspace(alg, a).map_or(false, |s| s.len() == k),
    }
}

/// Factor a blade into grade-1 vectors whose wedge is exactly `A`. `None` if `A`
/// is not a blade. A scalar returns itself as the single factor; a grade-`k`
/// blade returns `k` vectors.
pub fn factor_blade<S: Scalar>(
    alg: &CliffordAlgebra<S>,
    a: &Multivector<S>,
) -> Option<Vec<Multivector<S>>> {
    let k = homogeneous_grade(a)?;
    if k == 0 {
        return Some(vec![a.clone()]);
    }
    let basis = blade_subspace(alg, a)?;
    if basis.len() != k {
        return None; // not decomposable
    }
    let mut vecs: Vec<Multivector<S>> = basis.iter().map(|x| vector_from_coeffs(alg, x)).collect();
    // The wedge of the basis is some nonzero scalar multiple λ·A; rescale one
    // factor by 1/λ so the product is A exactly.
    let mut w = alg.scalar(S::one());
    for v in &vecs {
        w = alg.wedge(&w, v);
    }
    let mask = *a.terms.keys().next().expect("nonzero blade");
    let wa = w.terms.get(&mask).cloned().unwrap_or_else(S::zero);
    let aa = a.terms.get(&mask).cloned().expect("mask present in A");
    let lambda = wa.mul(&aa.inv()?);
    let linv = lambda.inv()?;
    vecs[0] = alg.scalar_mul(&linv, &vecs[0]);
    Some(vecs)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::clifford::Metric;
    use crate::scalar::Rational;

    fn r(n: i128) -> Rational {
        Rational::int(n)
    }
    fn euclid(n: usize) -> CliffordAlgebra<Rational> {
        CliffordAlgebra::new(n, Metric::diagonal(vec![r(1); n]))
    }

    #[test]
    fn simple_wedges_are_blades() {
        let alg = euclid(4);
        assert!(is_blade(&alg, &alg.scalar(r(3)))); // scalar
        assert!(is_blade(&alg, &alg.gen(1))); // vector
        let e01 = alg.wedge(&alg.gen(0), &alg.gen(1));
        assert!(is_blade(&alg, &e01));
        assert_eq!(blade_subspace(&alg, &e01).unwrap().len(), 2);
        let e012 = alg.wedge(&e01, &alg.gen(2));
        assert!(is_blade(&alg, &e012));
        assert_eq!(blade_subspace(&alg, &e012).unwrap().len(), 3);
    }

    #[test]
    fn non_simple_bivector_is_not_a_blade() {
        // e0∧e1 + e2∧e3 in R⁴ is the canonical non-decomposable 2-vector.
        let alg = euclid(4);
        let a = alg.add(
            &alg.wedge(&alg.gen(0), &alg.gen(1)),
            &alg.wedge(&alg.gen(2), &alg.gen(3)),
        );
        assert!(!is_blade(&alg, &a));
        assert_eq!(blade_subspace(&alg, &a).unwrap().len(), 0);
        assert!(factor_blade(&alg, &a).is_none());
    }

    #[test]
    fn factor_reconstructs_the_blade() {
        let alg = euclid(4);
        // a "skew" 2-blade: (e0+e1) ∧ (e2 + 2e3).
        let v = alg.add(&alg.gen(0), &alg.gen(1));
        let w = alg.add(&alg.gen(2), &alg.scalar_mul(&r(2), &alg.gen(3)));
        let blade = alg.wedge(&v, &w);
        let factors = factor_blade(&alg, &blade).unwrap();
        assert_eq!(factors.len(), 2);
        // wedging the factors back together reproduces the blade exactly.
        let mut prod = alg.scalar(r(1));
        for f in &factors {
            prod = alg.wedge(&prod, f);
        }
        assert_eq!(prod, blade);
    }

    #[test]
    fn mixed_grade_and_zero_are_not_blades() {
        let alg = euclid(3);
        let mixed = alg.add(&alg.scalar(r(1)), &alg.gen(0));
        assert!(!is_blade(&alg, &mixed));
        assert!(factor_blade(&alg, &mixed).is_none());
        assert!(!is_blade(&alg, &alg.zero()));
    }

    #[test]
    fn pseudoscalar_is_a_top_blade() {
        let alg = euclid(3);
        let i = alg.pseudoscalar();
        assert!(is_blade(&alg, &i));
        let factors = factor_blade(&alg, &i).unwrap();
        assert_eq!(factors.len(), 3);
        let mut prod = alg.scalar(r(1));
        for f in &factors {
            prod = alg.wedge(&prod, f);
        }
        assert_eq!(prod, i);
    }
}
