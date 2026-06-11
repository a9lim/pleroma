//! Integral lattices: the ℤ-Gram-matrix view of a quadratic form.
//!
//! The forms pillar elsewhere classifies a quadratic form *over a field* (by its
//! square classes / Witt class / Arf invariant). An **integral lattice** is the
//! complementary object: a free ℤ-module `L ≅ ℤⁿ` with an integer-valued
//! symmetric bilinear form, recorded by its Gram matrix `G = (⟨eᵢ, eⱼ⟩)`. Its
//! invariants are arithmetic, not just field-theoretic — the determinant, the
//! level, the minimum and kissing number, the automorphism group order — and the
//! coarse classification is the **genus** (local equivalence at every place),
//! built on the same p-adic primitives `local_global/padic.rs` and
//! `local_global/adelic.rs` already carry. This module is the M1 core (the
//! geometry of one lattice); `integral/root_lattices.rs`, `integral/genus.rs`,
//! and `integral/mass_formula.rs` build the A/D/E catalogue, the genus
//! equivalence, and the Conway–Sloane mass formula on top.
//!
//! Conventions. The **norm** of `x ∈ L` is `Q(x) = xᵀ G x` (so a "norm-2 vector"
//! has `Q = 2`, matching the root-lattice literature; this is twice the value of
//! the associated quadratic form `½Q` when the lattice is even). The geometric
//! routines — [`IntegralForm::minimum`], [`minimal_vectors`](IntegralForm::minimal_vectors),
//! [`kissing_number`](IntegralForm::kissing_number),
//! [`automorphism_group_order`](IntegralForm::automorphism_group_order) — assume the
//! lattice is **positive definite** and return `None` otherwise (an indefinite
//! lattice has infinitely many vectors of every norm and an infinite
//! automorphism group). Vectors are reported in lattice (basis) coordinates as
//! integer vectors, both signs included.
//!
//! Honest cutoff. Short-vector enumeration first tries an exact rational ellipsoid
//! box from `G⁻¹` when the box is small enough; larger boxes apply a conservative
//! unimodular size-reduction pass (integral shears/swaps, so the lattice is
//! unchanged), then run Fincke–Pohst (an LDL-bounded box search with exact norm
//! filtering) and map the vectors back to the original coordinates. Automorphism
//! counting first checks closed-form families: diagonal signed-permutation
//! lattices, literal `A`/`D`/`E` Cartan bases, and then basis-independent root
//! systems recovered from the norm-2 roots. Everything else falls back to a
//! backtracking search over basis images, which is **exponential** in general.
//! The fallback is bounded by an explicit node budget ([`AUTO_NODE_BUDGET`]);
//! when the search exceeds it the count is reported as `None` rather than
//! silently truncated. Use
//! [`automorphism_group_order_bounded`](IntegralForm::automorphism_group_order_bounded)
//! to choose the budget explicitly.
//!
//! # Module layout
//!
//! - `core` — [`IntegralForm`] struct + basic arithmetic (constructors,
//!   det/signature/level/Clifford metrics/direct_sum).
//! - `geometry` — short-vector enumeration (Fincke–Pohst), automorphism
//!   counting, and [`AUTO_NODE_BUDGET`].

mod core;
mod geometry;

pub use core::IntegralForm;
pub use geometry::AUTO_NODE_BUDGET;

// Re-export the test-visible exact-bounded helper used by lattice tests.
#[cfg(test)]
use geometry::SHORT_VECTOR_EXACT_ENUM_LIMIT;

#[cfg(test)]
mod tests {
    use super::*;

    fn a_n(n: usize) -> IntegralForm {
        // A_n Cartan matrix: 2 on the diagonal, -1 on the off-diagonals.
        let mut g = vec![vec![0i128; n]; n];
        for i in 0..n {
            g[i][i] = 2;
            if i + 1 < n {
                g[i][i + 1] = -1;
                g[i + 1][i] = -1;
            }
        }
        IntegralForm::new(g).unwrap()
    }

    fn d4() -> IntegralForm {
        IntegralForm::new(vec![
            vec![2, -1, 0, 0],
            vec![-1, 2, -1, -1],
            vec![0, -1, 2, 0],
            vec![0, -1, 0, 2],
        ])
        .unwrap()
    }

    fn e8() -> IntegralForm {
        // E_8 Cartan matrix (Bourbaki labelling): even unimodular, det 1.
        IntegralForm::new(vec![
            vec![2, -1, 0, 0, 0, 0, 0, 0],
            vec![-1, 2, -1, 0, 0, 0, 0, 0],
            vec![0, -1, 2, -1, 0, 0, 0, 0],
            vec![0, 0, -1, 2, -1, 0, 0, 0],
            vec![0, 0, 0, -1, 2, -1, 0, -1],
            vec![0, 0, 0, 0, -1, 2, -1, 0],
            vec![0, 0, 0, 0, 0, -1, 2, 0],
            vec![0, 0, 0, 0, -1, 0, 0, 2],
        ])
        .unwrap()
    }

    fn permute_basis(l: &IntegralForm, perm: &[usize]) -> IntegralForm {
        let n = l.dim();
        assert_eq!(perm.len(), n);
        let mut g = vec![vec![0i128; n]; n];
        for i in 0..n {
            for j in 0..n {
                g[i][j] = l.gram()[perm[i]][perm[j]];
            }
        }
        IntegralForm::new(g).unwrap()
    }

    #[test]
    fn rejects_non_symmetric() {
        assert!(IntegralForm::new(vec![vec![1, 2], vec![3, 4]]).is_none());
        assert!(IntegralForm::new(vec![vec![1, 2, 3], vec![2, 4]]).is_none());
        assert!(IntegralForm::new(vec![vec![2, -1], vec![-1, 2]]).is_some());
    }

    #[test]
    fn determinants_and_evenness() {
        assert_eq!(a_n(2).determinant(), 3);
        assert_eq!(a_n(3).determinant(), 4);
        assert_eq!(d4().determinant(), 4);
        assert_eq!(e8().determinant(), 1);
        assert!(e8().is_unimodular());
        assert!(e8().is_even());
        assert!(a_n(2).is_even());
        // Z^3 is odd unimodular.
        let z3 = IntegralForm::diagonal(&[1, 1, 1]);
        assert_eq!(z3.determinant(), 1);
        assert!(z3.is_unimodular());
        assert!(!z3.is_even());
    }

    #[test]
    fn invariant_factors_track_discriminant_group() {
        assert_eq!(a_n(2).invariant_factors(), vec![1, 3]); // ℤ/3
        assert_eq!(d4().invariant_factors(), vec![1, 1, 2, 2]); // (ℤ/2)²
        assert_eq!(e8().invariant_factors(), vec![1, 1, 1, 1, 1, 1, 1, 1]);
        // product of nonzero factors = |det|
        let prod: i128 = d4().invariant_factors().iter().product();
        assert_eq!(prod, d4().determinant().abs());
    }

    #[test]
    fn levels_match_known_values() {
        assert_eq!(IntegralForm::diagonal(&[2]).level(), Some(4)); // A_1 = ⟨2⟩
        assert_eq!(a_n(2).level(), Some(3)); // hexagonal lattice, level 3
        assert_eq!(e8().level(), Some(1)); // even unimodular
                                           // ℤ = ⟨1⟩ is odd: G⁻¹ = [1] has odd diagonal, so the smallest N making
                                           // N·G⁻¹ even-integral is 2 (cf. A_1 = ⟨2⟩ → 4).
        assert_eq!(IntegralForm::diagonal(&[1]).level(), Some(2));
    }

    #[test]
    fn signature_handles_indefinite_and_skew_bases() {
        assert_eq!(IntegralForm::diagonal(&[1, 1, -1]).signature(), (2, 1));
        let hyp = IntegralForm::new(vec![vec![0, 1], vec![1, 0]]).unwrap();
        assert_eq!(hyp.signature(), (1, 1));
        assert_eq!(
            IntegralForm::new(vec![vec![0, 0], vec![0, 0]])
                .unwrap()
                .signature(),
            (0, 0)
        );
    }

    #[test]
    fn lattice_clifford_metrics_preserve_q_and_polar_data() {
        use crate::scalar::{Nimber, Rational};
        let a2 = a_n(2);
        let rat = a2.clifford_metric();
        assert_eq!(rat.q, vec![Rational::int(2), Rational::int(2)]);
        assert_eq!(rat.b[&(0, 1)], Rational::int(-2));

        let f2 = a2.clifford_metric_f2().unwrap();
        assert_eq!(f2.q, vec![Nimber(1), Nimber(1)]);
        assert_eq!(f2.b[&(0, 1)], Nimber(1));
        assert!(IntegralForm::diagonal(&[1]).clifford_metric_f2().is_none());
    }

    #[test]
    fn minimum_and_kissing_numbers() {
        // Root lattices: minimum 2, kissing = number of roots.
        assert_eq!(a_n(2).minimum(), Some(2));
        assert_eq!(a_n(2).kissing_number(), Some(6)); // n(n+1) = 6
        assert_eq!(a_n(3).kissing_number(), Some(12)); // 3·4
        assert_eq!(d4().minimum(), Some(2));
        assert_eq!(d4().kissing_number(), Some(24)); // 2n(n-1) = 24
        assert_eq!(e8().minimum(), Some(2));
        assert_eq!(e8().kissing_number(), Some(240));
        // ℤ²: minimum 1, the four ±eᵢ.
        let z2 = IntegralForm::diagonal(&[1, 1]);
        assert_eq!(z2.minimum(), Some(1));
        assert_eq!(z2.kissing_number(), Some(4));
    }

    #[test]
    fn short_vectors_return_original_coordinates_after_basis_reduction() {
        // Uᵀ I U for U = [[1, 10], [0, 1]] is a badly skewed basis of Z².
        // The norm-1 vectors in this basis are ±(1,0) and ±(-10,1).
        let g = IntegralForm::new(vec![vec![1, 10], vec![10, 101]]).unwrap();
        let mut exact = g
            .short_vectors_exact_bounded(1, SHORT_VECTOR_EXACT_ENUM_LIMIT)
            .expect("small rational ellipsoid box is enumerated exactly");
        exact.sort();
        let mut vecs = g.short_vectors(1).unwrap();
        vecs.sort();
        assert_eq!(exact, vecs);
        assert_eq!(
            vecs,
            vec![vec![-10, 1], vec![-1, 0], vec![1, 0], vec![10, -1]]
        );
        assert!(vecs.iter().all(|v| g.norm(v) == 1));
    }

    #[test]
    fn short_vectors_are_indefinite_safe() {
        // An indefinite form has no finite short-vector set.
        let hyp = IntegralForm::new(vec![vec![0, 1], vec![1, 0]]).unwrap();
        assert!(!hyp.is_positive_definite());
        assert_eq!(hyp.short_vectors(4), None);
        assert_eq!(hyp.minimum(), None);
        assert_eq!(hyp.automorphism_group_order(), None);
    }

    #[test]
    fn automorphism_orders_match_known() {
        // Aut(Z^n) = signed permutations = 2^n · n!.
        assert_eq!(
            IntegralForm::diagonal(&[1, 1]).automorphism_group_order(),
            Some(8)
        );
        assert_eq!(
            IntegralForm::diagonal(&[1, 1, 1]).automorphism_group_order(),
            Some(48)
        );
        // Aut(A_2) = dihedral of order 12 (W(A_2)=S_3 times ±1).
        assert_eq!(a_n(2).automorphism_group_order(), Some(12));
        // Aut(A_3) = W(A_3) × {±1} = 24 · 2 = 48.
        assert_eq!(a_n(3).automorphism_group_order(), Some(48));
        // |Aut(D_4)| = 1152.
        assert_eq!(d4().automorphism_group_order(), Some(1152));
        // E_8 is recognized by its standard Cartan basis instead of brute-forced.
        assert_eq!(e8().automorphism_group_order_bounded(1), Some(696_729_600));
    }

    #[test]
    fn automorphism_budget_cutoff_reports_none() {
        // Permuted root bases are now recognized by the root-system fast path,
        // independent of the standard Cartan syntax.
        let d4_permuted = permute_basis(&d4(), &[2, 0, 1, 3]);
        assert_eq!(d4_permuted.automorphism_group_order_bounded(1), Some(1152));

        // A tiny budget still forces the fallback search to give up rather than
        // silently truncating on a non-root lattice: an honest None, not a wrong count.
        let generic = IntegralForm::new(vec![vec![2, 1], vec![1, 3]]).unwrap();
        assert_eq!(generic.automorphism_group_order_bounded(0), None);
    }

    #[test]
    fn direct_sum_is_block_diagonal() {
        let sum = a_n(2).direct_sum(&IntegralForm::diagonal(&[1]));
        assert_eq!(sum.dim(), 3);
        assert_eq!(sum.determinant(), 3); // det(A_2) · det(⟨1⟩)
                                          // E_8 ⟂ E_8 is rank-16 even unimodular.
        let e8e8 = e8().direct_sum(&e8());
        assert_eq!(e8e8.dim(), 16);
        assert_eq!(e8e8.determinant(), 1);
        assert!(e8e8.is_even());
        for i in 0..8 {
            for j in 8..16 {
                assert_eq!(e8e8.gram()[i][j], 0);
            }
        }
    }

    #[test]
    fn ldl_returns_none_on_indefinite_gram() {
        // The internal ldl() helper must return None rather than producing a
        // non-positive pivot when called on an indefinite Gram matrix. This
        // guards against the search silently dropping short vectors due to a
        // divide-by-zero or negative-sqrt in the float bound.
        let hyp = IntegralForm::new(vec![vec![0, 1], vec![1, 0]]).unwrap();
        assert!(hyp.ldl().is_none());
        // A positive-definite lattice must produce a valid decomposition.
        assert!(a_n(2).ldl().is_some());
        let (d, _) = a_n(2).ldl().unwrap();
        assert!(d.iter().all(|&di| di > 0.0));
    }
}
