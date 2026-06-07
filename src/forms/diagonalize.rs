//! Congruence diagonalization of a symmetric bilinear form (characteristic ≠ 2).
//!
//! The classifiers (`char0`, `oddchar`) read a form off the *diagonal* of a
//! metric. A metric with a nonzero off-diagonal polar form `b` is the same
//! quadratic form in a skew basis; in characteristic ≠ 2 it can always be
//! brought to an orthogonal basis by a congruence `P → PᵀMP` (symmetric Gaussian
//! elimination). This module performs that reduction, so every classifier in
//! `forms/` accepts an *arbitrary* (symmetric) metric, not only a diagonal one.
//!
//! Characteristic 2 is genuinely different — a nonsingular char-2 form is *not*
//! diagonalizable (its polar form is alternating), which is exactly why that leg
//! uses the symplectic Arf reduction (`forms::char2`) instead. [`diagonalize`]
//! returns `None` there, on the nose.
//!
//! The antisymmetric (`a`) part of a general bilinear form is a gauge — it does
//! not change the quadratic form — so it is ignored here: only `q` (the squares)
//! and `b` (the symmetric polar form) determine the diagonalization.

use crate::clifford::Metric;
use crate::scalar::Scalar;

/// The Gram matrix of the symmetric bilinear form `B` with `B(eᵢ,eᵢ) = q[i]` and
/// `B(eᵢ,eⱼ) = b[(i,j)]/2` for `i<j`. `None` in characteristic 2 (no `½`).
pub fn gram<S: Scalar>(m: &Metric<S>) -> Option<Vec<Vec<S>>> {
    let two = S::one().add(&S::one());
    let half = two.inv()?; // None ⇔ characteristic 2
    let n = m.q.len();
    let mut g = vec![vec![S::zero(); n]; n];
    for (i, qi) in m.q.iter().enumerate() {
        g[i][i] = qi.clone();
    }
    for (&(i, j), bij) in &m.b {
        let off = bij.mul(&half);
        g[j][i] = off.clone();
        g[i][j] = off;
    }
    Some(g)
}

/// Swap generators `k` and `m`: simultaneous row and column swap (a congruence).
fn swap_sym<S: Scalar>(g: &mut [Vec<S>], k: usize, m: usize) {
    g.swap(k, m);
    for row in g.iter_mut() {
        row.swap(k, m);
    }
}

/// Replace generator `i` by `eᵢ + eⱼ`: add row/column `j` into row/column `i`.
/// Used to manufacture a nonzero diagonal pivot when the diagonal vanishes.
fn add_sym<S: Scalar>(g: &mut [Vec<S>], i: usize, j: usize) {
    let n = g.len();
    for t in 0..n {
        g[i][t] = g[i][t].add(&g[j][t].clone());
    }
    for t in 0..n {
        g[t][i] = g[t][i].add(&g[t][j].clone());
    }
}

/// Ensure `g[k][k]` is nonzero by a congruence, if the trailing block `[k..]` is
/// not entirely zero. Returns `false` when that block is all zero (the radical).
fn ensure_pivot<S: Scalar>(g: &mut [Vec<S>], k: usize) -> bool {
    let n = g.len();
    if !g[k][k].is_zero() {
        return true;
    }
    // A nonzero diagonal entry further down: swap it up.
    for m in (k + 1)..n {
        if !g[m][m].is_zero() {
            swap_sym(g, k, m);
            return true;
        }
    }
    // All diagonals zero: find a nonzero off-diagonal, add to make a diagonal.
    for i in k..n {
        for j in (i + 1)..n {
            if !g[i][j].is_zero() {
                add_sym(g, i, j); // g[i][i] becomes 2·g[i][j] ≠ 0 (char ≠ 2)
                if i != k {
                    swap_sym(g, k, i);
                }
                return true;
            }
        }
    }
    false // the whole trailing block is zero ⇒ radical
}

/// An isometric **diagonal** metric for a symmetric form, by congruence
/// (Gram–Schmidt over the field). `None` in characteristic 2, where nonsingular
/// forms are not diagonalizable (use the Arf reduction in `forms::char2`), or when
/// a generic scalar-ring call encounters a nonzero nonunit pivot.
///
/// The diagonal entries are the squares of an orthogonal basis; null entries are
/// the radical. Equal as a quadratic form to the input — the classifiers may be
/// run on the result.
pub fn diagonalize<S: Scalar>(m: &Metric<S>) -> Option<Metric<S>> {
    let mut g = gram(m)?;
    let n = g.len();
    for k in 0..n {
        if !ensure_pivot(&mut g, k) {
            break; // remaining block is the radical (all zero)
        }
        let pivot_inv = g[k][k].inv()?;
        for r in (k + 1)..n {
            if g[r][k].is_zero() {
                continue;
            }
            let factor = g[r][k].mul(&pivot_inv);
            // row r -= factor · row k, then col r -= factor · col k (congruence)
            let row_k: Vec<S> = g[k].clone();
            for t in 0..n {
                let sub = factor.mul(&row_k[t]);
                g[r][t] = g[r][t].sub(&sub);
            }
            let col_k: Vec<S> = (0..n).map(|t| g[t][k].clone()).collect();
            for t in 0..n {
                let sub = factor.mul(&col_k[t]);
                g[t][r] = g[t][r].sub(&sub);
            }
        }
    }
    let diag: Vec<S> = (0..n).map(|i| g[i][i].clone()).collect();
    Some(Metric::diagonal(diag))
}

/// The input unchanged if already diagonal (`b`/`a` empty), else its congruence
/// [`diagonalize`]. `None` only when diagonalization is needed but impossible
/// (characteristic 2).
pub fn as_diagonal<S: Scalar>(m: &Metric<S>) -> Option<Metric<S>> {
    if m.b.is_empty() && m.a.is_empty() {
        Some(m.clone())
    } else {
        diagonalize(m)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forms::{classify_oddchar, classify_rational, classify_surreal};
    use crate::scalar::{Fp, Rational, Surreal, Zp};
    use std::collections::BTreeMap;

    fn rat(n: i128) -> Rational {
        Rational::int(n)
    }

    #[test]
    fn hyperbolic_plane_diagonalizes_to_disc_minus_one() {
        // q = [0,0], {e0,e1} = 2 ⇒ B(e0,e1)=1: the hyperbolic plane H.
        let mut b = BTreeMap::new();
        b.insert((0, 1), rat(2));
        let m = Metric::new(vec![rat(0), rat(0)], b);
        let d = diagonalize(&m).unwrap();
        // det of the diagonalization is −1 mod squares (the disc of H).
        let det = d.q.iter().fold(rat(1), |acc, x| acc.mul(x));
        // det = a · b with the two diagonal entries; over ℚ it should be −(square).
        assert_eq!(det.sign(), std::cmp::Ordering::Less);
        // and it classifies as the same algebra as the diagonal ⟨1,−1⟩.
        assert_eq!(
            classify_rational(&m).unwrap(),
            classify_rational(&Metric::diagonal(vec![rat(1), rat(-1)])).unwrap()
        );
    }

    #[test]
    fn off_diagonal_real_form_keeps_its_signature() {
        // ⟨2,2⟩ skewed by an off-diagonal: q=[2,2], B(e0,e1)=1 (b=2). Gram
        // [[2,1],[1,2]] is positive-definite (eigenvalues 1,3) ⇒ signature (2,0).
        let mut b = BTreeMap::new();
        b.insert((0, 1), rat(2));
        let m = Metric::new(vec![rat(2), rat(2)], b);
        let d = diagonalize(&m).unwrap();
        assert!(d.q.iter().all(|x| x.sign() == std::cmp::Ordering::Greater));
        assert_eq!(
            classify_surreal(&Metric::diagonal(
                d.q.iter()
                    .map(|x| Surreal::from_rational(x.clone()))
                    .collect()
            )),
            classify_surreal(&Metric::diagonal(vec![
                Surreal::from_int(1),
                Surreal::from_int(1)
            ]))
        );
    }

    #[test]
    fn classifiers_now_accept_nondiagonal_metrics() {
        // The same hyperbolic plane, over F_5: classify it without pre-diagonalizing.
        const P: u128 = 5;
        let mut b = BTreeMap::new();
        b.insert((0, 1), Fp::<P>::new(2));
        let m = Metric::new(vec![Fp::<P>::new(0), Fp::<P>::new(0)], b);
        let got = classify_oddchar(&m).unwrap();
        let want =
            classify_oddchar(&Metric::diagonal(vec![Fp::<P>::new(1), Fp::<P>::new(-1)])).unwrap();
        assert_eq!(got.dim, want.dim);
        assert_eq!(got.disc_is_square, want.disc_is_square);
    }

    #[test]
    fn characteristic_two_is_not_diagonalizable() {
        use crate::scalar::Nimber;
        let mut b = BTreeMap::new();
        b.insert((0, 1), Nimber(1));
        let m = Metric::new(vec![Nimber(1), Nimber(1)], b);
        assert!(diagonalize(&m).is_none());
        assert!(as_diagonal(&m).is_none());
    }

    #[test]
    fn nonfield_nonunit_pivot_returns_none() {
        let m = Metric::diagonal(vec![Zp::<3, 2>::new(3)]);
        assert!(diagonalize(&m).is_none());
    }
}
