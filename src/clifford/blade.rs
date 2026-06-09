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
//! `is_blade` is ring-native: it checks the quadratic Pluecker relations, so it
//! does not rely on field division and correctly handles nonunit scalar
//! multiples such as `2e0` over `Z`. `blade_subspace` and the general
//! factorization path still return a free basis only when unit-pivot linear
//! algebra suffices; monomial blades and vectors are factored without division.

use crate::clifford::{bits, grade, CliffordAlgebra, Multivector};
use crate::linalg::field;
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

fn combinations(n: usize, k: usize) -> Vec<u128> {
    fn rec(out: &mut Vec<u128>, n: usize, k: usize, start: usize, mask: u128) {
        if k == 0 {
            out.push(mask);
            return;
        }
        for i in start..=n - k {
            rec(out, n, k - 1, i + 1, mask | (1u128 << i));
        }
    }
    if k > n {
        return Vec::new();
    }
    let mut out = Vec::new();
    rec(&mut out, n, k, 0, 0);
    out
}

fn coeff<S: Scalar>(a: &Multivector<S>, mask: u128) -> S {
    a.terms.get(&mask).cloned().unwrap_or_else(S::zero)
}

fn higher_bits(mask: u128, i: usize) -> usize {
    if i + 1 >= u128::BITS as usize {
        0
    } else {
        (mask >> (i + 1)).count_ones() as usize
    }
}

/// The Pluecker equations for the affine cone over `Gr(k,n)`. They are
/// polynomial, hence valid over any commutative scalar ring.
fn plucker_relations_hold<S: Scalar>(
    alg: &CliffordAlgebra<S>,
    a: &Multivector<S>,
    k: usize,
) -> bool {
    let n = alg.dim;
    if k == 0 || k == 1 || k == n {
        return true;
    }
    for i_mask in combinations(n, k - 1) {
        for j_mask in combinations(n, k + 1) {
            let mut acc = S::zero();
            let mut jj = j_mask;
            let mut pos = 0usize;
            while jj != 0 {
                let j = jj.trailing_zeros() as usize;
                let bit = 1u128 << j;
                jj &= jj - 1;
                if i_mask & bit != 0 {
                    pos += 1;
                    continue;
                }
                let mut term = coeff(a, i_mask | bit).mul(&coeff(a, j_mask ^ bit));
                if (pos + higher_bits(i_mask, j)) & 1 == 1 {
                    term = term.neg();
                }
                acc = acc.add(&term);
                pos += 1;
            }
            if !acc.is_zero() {
                return false;
            }
        }
    }
    true
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
    if k == 1 {
        let mut x = vec![S::zero(); n];
        for (&mask, c) in &a.terms {
            let i = mask.trailing_zeros() as usize;
            x[i] = c.clone();
        }
        return Some(vec![x]);
    }
    if a.terms.len() == 1 {
        let (&mask, _) = a.terms.iter().next().expect("single term");
        let gens = bits(mask);
        if gens.len() == k {
            let mut basis = Vec::with_capacity(k);
            for g in gens {
                let mut x = vec![S::zero(); n];
                x[g] = S::one();
                basis.push(x);
            }
            return Some(basis);
        }
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
    field::unit_pivot_nullspace(mat, n)
}

/// Whether `A` is a **blade** (a decomposable homogeneous multivector). Scalars,
/// vectors, and top-grade homogeneous multivectors always are; intermediate
/// grades are checked by the Pluecker relations over the scalar ring. Zero and
/// mixed-grade multivectors are not blades.
pub fn is_blade<S: Scalar>(alg: &CliffordAlgebra<S>, a: &Multivector<S>) -> bool {
    match homogeneous_grade(a) {
        None => false,
        Some(0) => true,
        Some(k) if k <= alg.dim => plucker_relations_hold(alg, a, k),
        Some(_) => false,
    }
}

fn monomial_factor<S: Scalar>(
    alg: &CliffordAlgebra<S>,
    a: &Multivector<S>,
    k: usize,
) -> Option<Vec<Multivector<S>>> {
    if a.terms.len() != 1 {
        return None;
    }
    let (&mask, coeff) = a.terms.iter().next()?;
    let gens = bits(mask);
    if gens.len() != k || gens.is_empty() {
        return None;
    }
    let mut out = Vec::with_capacity(k);
    out.push(alg.scalar_mul(coeff, &alg.gen(gens[0])));
    for &g in &gens[1..] {
        out.push(alg.gen(g));
    }
    Some(out)
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
    if k == 1 {
        return Some(vec![a.clone()]);
    }
    if !is_blade(alg, a) {
        return None;
    }
    if let Some(factors) = monomial_factor(alg, a, k) {
        return Some(factors);
    }
    let basis = blade_subspace(alg, a)?;
    if basis.len() != k {
        return None;
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
    use crate::scalar::{Integer, Rational};

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

    #[test]
    fn integer_nonunit_multiples_are_blades() {
        let alg = CliffordAlgebra::new(3, Metric::<Integer>::grassmann(3));
        let two = Integer(2);
        let v = alg.scalar_mul(&two, &alg.gen(0));
        assert!(is_blade(&alg, &v));
        assert_eq!(factor_blade(&alg, &v).unwrap(), vec![v.clone()]);

        let e01 = alg.wedge(&alg.gen(0), &alg.gen(1));
        let two_e01 = alg.scalar_mul(&two, &e01);
        assert!(is_blade(&alg, &two_e01));
        assert_eq!(blade_subspace(&alg, &two_e01).unwrap().len(), 2);
        let factors = factor_blade(&alg, &two_e01).unwrap();
        let mut prod = alg.scalar(Integer(1));
        for f in &factors {
            prod = alg.wedge(&prod, f);
        }
        assert_eq!(prod, two_e01);
    }

    #[test]
    fn pluecker_rejects_integer_non_simple_bivector() {
        let alg = CliffordAlgebra::new(4, Metric::<Integer>::grassmann(4));
        let a = alg.add(
            &alg.wedge(&alg.gen(0), &alg.gen(1)),
            &alg.wedge(&alg.gen(2), &alg.gen(3)),
        );
        assert!(!is_blade(&alg, &a));
        assert!(factor_blade(&alg, &a).is_none());
    }

    #[test]
    fn integer_blade_subspace_refuses_nonunit_kernel_pivot() {
        let alg = CliffordAlgebra::new(3, Metric::<Integer>::grassmann(3));
        let minus_two = Integer(-2);
        let e01 = alg.wedge(&alg.gen(0), &alg.gen(1));
        let e02 = alg.wedge(&alg.gen(0), &alg.gen(2));
        let a = alg.add(
            &alg.scalar_mul(&minus_two, &e01),
            &alg.scalar_mul(&minus_two, &e02),
        );

        assert!(is_blade(&alg, &a));
        assert!(blade_subspace(&alg, &a).is_none());
        assert!(factor_blade(&alg, &a).is_none());
    }
}
