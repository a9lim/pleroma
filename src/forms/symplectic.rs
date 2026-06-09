//! **Symplectic (alternating) forms** — the skew member of the "form + involution"
//! family, completing it beside the symmetric bilinear forms (the rest of this
//! pillar) and the [`HermitianForm`](crate::forms::HermitianForm).
//!
//! An alternating bilinear form `B(x, x) = 0` (equivalently `Bᵀ = −B` *and* zero
//! diagonal — the diagonal condition is the genuine constraint in characteristic
//! 2, where `−B = B`) has the simplest classification in all of form theory: over
//! **any** field it is congruent to an orthogonal sum of hyperbolic planes and a
//! zero radical,
//!
//! ```text
//! B  ≅  (rank/2) · H  ⟂  0^{radical}
//! ```
//!
//! so the complete invariant is just `(rank, radical_dim)` with `rank` always even
//! — there is no characteristic trichotomy to dispatch (unlike the symmetric and
//! Hermitian cases), so this is a single generic routine. The radical is the
//! kernel `{x : Bx = 0}` (left and right kernels coincide for an alternating form).
//! Classification returns `None` over ring backends when the shared unit-pivot
//! solver encounters a nonunit pivot.

use crate::scalar::Scalar;

/// A symplectic (alternating) form, carried by its alternating Gram matrix.
#[derive(Debug, Clone, PartialEq)]
pub struct SymplecticForm<S: Scalar> {
    gram: Vec<Vec<S>>,
}

/// The complete invariant of an alternating form: its rank (always even, twice the
/// number of hyperbolic planes) and the dimension of its radical.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SymplecticClass {
    /// `2 × (number of hyperbolic planes)` — always even.
    pub rank: usize,
    /// Dimension of the radical (the kernel of the form).
    pub radical_dim: usize,
}

impl SymplecticClass {
    /// The number of hyperbolic planes in the canonical decomposition.
    pub fn planes(&self) -> usize {
        self.rank / 2
    }
}

impl<S: Scalar> SymplecticForm<S> {
    /// Build from a Gram matrix, checking it is square and **alternating**: zero
    /// diagonal and `A[i][j] = −A[j][i]`. Returns `None` otherwise. (In char 2 the
    /// off-diagonal condition reads as symmetry, so the explicit zero-diagonal
    /// check is what distinguishes alternating from merely symmetric there.)
    pub fn from_gram(gram: Vec<Vec<S>>) -> Option<Self> {
        let n = gram.len();
        for row in &gram {
            if row.len() != n {
                return None;
            }
        }
        for i in 0..n {
            if !gram[i][i].is_zero() {
                return None;
            }
            for j in (i + 1)..n {
                if gram[i][j] != gram[j][i].neg() {
                    return None;
                }
            }
        }
        Some(SymplecticForm { gram })
    }

    /// The standard symplectic form `r · H` on `2r` generators: the block-diagonal
    /// sum of `r` hyperbolic planes `[[0, 1], [−1, 0]]`.
    pub fn hyperbolic(r: usize) -> Self {
        let n = 2 * r;
        let mut gram = vec![vec![S::zero(); n]; n];
        for k in 0..r {
            gram[2 * k][2 * k + 1] = S::one();
            gram[2 * k + 1][2 * k] = S::one().neg();
        }
        SymplecticForm { gram }
    }

    pub fn dim(&self) -> usize {
        self.gram.len()
    }

    /// The orthogonal direct sum (block-diagonal Gram).
    pub fn direct_sum(&self, other: &SymplecticForm<S>) -> SymplecticForm<S> {
        let (n, m) = (self.dim(), other.dim());
        let mut gram = vec![vec![S::zero(); n + m]; n + m];
        for i in 0..n {
            for j in 0..n {
                gram[i][j] = self.gram[i][j].clone();
            }
        }
        for i in 0..m {
            for j in 0..m {
                gram[n + i][n + j] = other.gram[i][j].clone();
            }
        }
        SymplecticForm { gram }
    }

    /// Classify the form: `(rank, radical_dim)`, the complete invariant over
    /// fields. The radical is the nullspace of the Gram; the rank is
    /// `dim − radical_dim` and is always even. Returns `None` when unit-pivot
    /// elimination cannot decide the kernel over a non-field scalar ring.
    pub fn classify(&self) -> Option<SymplecticClass> {
        let n = self.dim();
        let radical_dim = crate::linalg::field::unit_pivot_nullspace(self.gram.clone(), n)?.len();
        Some(SymplecticClass {
            rank: n - radical_dim,
            radical_dim,
        })
    }
}

/// Classify an alternating Gram matrix directly, or `None` if it is not square and
/// alternating. Convenience over [`SymplecticForm::from_gram`] + `classify`.
pub fn classify_symplectic<S: Scalar>(gram: Vec<Vec<S>>) -> Option<SymplecticClass> {
    SymplecticForm::from_gram(gram)?.classify()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scalar::{Nimber, Rational};

    fn r(n: i128) -> Rational {
        Rational::int(n)
    }

    #[test]
    fn hyperbolic_plane_has_rank_two() {
        let h = SymplecticForm::<Rational>::hyperbolic(1);
        assert_eq!(
            h.classify().unwrap(),
            SymplecticClass {
                rank: 2,
                radical_dim: 0
            }
        );
        assert_eq!(h.classify().unwrap().planes(), 1);
    }

    #[test]
    fn rank_is_always_even_and_radical_splits_off() {
        // 2 planes ⟂ a 1-dim radical: rank 4, radical 1.
        let f = SymplecticForm::<Rational>::hyperbolic(2).direct_sum(
            &SymplecticForm::from_gram(vec![vec![r(0)]]).unwrap(), // the zero form on 1 gen
        );
        let c = f.classify().unwrap();
        assert_eq!((c.rank, c.radical_dim), (4, 1));
        assert_eq!(c.rank % 2, 0);
    }

    #[test]
    fn non_alternating_is_rejected() {
        // nonzero diagonal: not alternating.
        assert!(SymplecticForm::from_gram(vec![vec![r(1), r(0)], vec![r(0), r(0)]]).is_none());
        // symmetric off-diagonal over char 0: not alternating (A[0][1] ≠ −A[1][0]).
        assert!(SymplecticForm::from_gram(vec![vec![r(0), r(1)], vec![r(1), r(0)]]).is_none());
    }

    #[test]
    fn char_two_alternating_is_symmetric_with_zero_diagonal() {
        // Over a nim-field, −1 = 1, so an alternating form is a symmetric matrix
        // with zero diagonal. [[0,1],[1,0]] is a hyperbolic plane.
        let h =
            SymplecticForm::from_gram(vec![vec![Nimber(0), Nimber(1)], vec![Nimber(1), Nimber(0)]])
                .unwrap();
        assert_eq!(
            h.classify().unwrap(),
            SymplecticClass {
                rank: 2,
                radical_dim: 0
            }
        );
        // but a nonzero diagonal is still rejected (alternating ⊋ symmetric).
        assert!(SymplecticForm::from_gram(vec![
            vec![Nimber(1), Nimber(1)],
            vec![Nimber(1), Nimber(0)],
        ])
        .is_none());
    }

    #[test]
    fn degenerate_form_is_all_radical() {
        // the zero form on 3 generators: rank 0, radical 3.
        let z = SymplecticForm::<Rational>::from_gram(vec![vec![r(0); 3]; 3]).unwrap();
        assert_eq!(
            z.classify().unwrap(),
            SymplecticClass {
                rank: 0,
                radical_dim: 3
            }
        );
    }

    #[test]
    fn free_function_matches_method() {
        let g = SymplecticForm::<Rational>::hyperbolic(3);
        assert_eq!(classify_symplectic(vec_gram(&g)), g.classify());
    }

    #[test]
    fn nonfield_nonunit_pivot_is_refused() {
        use crate::scalar::Integer;

        let gram = vec![vec![Integer(0), Integer(2)], vec![Integer(-2), Integer(0)]];
        let f = SymplecticForm::from_gram(gram).unwrap();
        assert_eq!(f.classify(), None);
    }

    fn vec_gram(f: &SymplecticForm<Rational>) -> Vec<Vec<Rational>> {
        (0..f.dim())
            .map(|i| (0..f.dim()).map(|j| f.gram[i][j].clone()).collect())
            .collect()
    }
}
