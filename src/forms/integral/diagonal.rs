//! Shared exact diagonalization routines for integral Gram matrices.

use crate::scalar::{Rational, Scalar};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum DegenerateBehavior {
    /// Stop when the remaining active block is the radical.
    StopAtRadical,
    /// Treat a radical as an internal invariant violation.
    RequireNonsingular,
}

fn rdiv(a: &Rational, b: &Rational) -> Rational {
    a.mul(
        &b.inv()
            .expect("division by zero rational in congruence diagonalization"),
    )
}

/// Diagonal entries reached by exact rational congruence diagonalization.
///
/// When every active diagonal entry is zero but an off-diagonal entry remains,
/// the basis shear `e_r <- e_r + e_s` creates a nonzero diagonal pivot. This is
/// the common Jacobi/Sylvester step used by the signature, genus, and
/// discriminant layers.
pub(crate) fn rational_congruence_diagonal(
    gram: &[Vec<i128>],
    degenerate: DegenerateBehavior,
) -> Vec<Rational> {
    let n = gram.len();
    let mut a: Vec<Vec<Rational>> = gram
        .iter()
        .map(|row| row.iter().map(|&x| Rational::int(x)).collect())
        .collect();
    let mut active: Vec<usize> = (0..n).collect();
    let mut out = Vec::with_capacity(n);
    while !active.is_empty() {
        if !active.iter().any(|&r| !a[r][r].is_zero()) {
            let mut pair = None;
            'find: for (ai, &r) in active.iter().enumerate() {
                for &s in &active[ai + 1..] {
                    if !a[r][s].is_zero() {
                        pair = Some((r, s));
                        break 'find;
                    }
                }
            }
            let Some((r, s)) = pair else {
                match degenerate {
                    DegenerateBehavior::StopAtRadical => break,
                    DegenerateBehavior::RequireNonsingular => {
                        panic!("nondegenerate form has a nonzero entry")
                    }
                }
            };
            for &c in &active {
                a[r][c] = a[r][c].add(&a[s][c].clone());
            }
            for &rr in &active {
                a[rr][r] = a[rr][r].add(&a[rr][s].clone());
            }
        }
        let Some(i) = active.iter().copied().find(|&r| !a[r][r].is_zero()) else {
            match degenerate {
                DegenerateBehavior::StopAtRadical => break,
                DegenerateBehavior::RequireNonsingular => panic!("a diagonal pivot now exists"),
            }
        };
        let pivot = a[i][i].clone();
        out.push(pivot.clone());
        let rest: Vec<usize> = active.iter().copied().filter(|&r| r != i).collect();
        for &r in &rest {
            for &s in &rest {
                let corr = rdiv(&a[r][i].mul(&a[i][s]), &pivot);
                a[r][s] = a[r][s].sub(&corr);
            }
        }
        active = rest;
    }
    out
}

pub(crate) fn signature_from_diagonal(diag: &[Rational]) -> (usize, usize) {
    let (mut pos, mut neg) = (0usize, 0usize);
    for d in diag {
        match d.sign() {
            std::cmp::Ordering::Greater => pos += 1,
            std::cmp::Ordering::Less => neg += 1,
            std::cmp::Ordering::Equal => {}
        }
    }
    (pos, neg)
}
