//! **Hermitian forms** over the surcomplex field — the natural quadratic-form
//! structure over a field carrying an involution, which the rest of the forms
//! pillar (symmetric/bilinear) never used.
//!
//! [`Surcomplex`] carries the conjugation `i ↦ −i` ([`Surcomplex::conj`]); a
//! Hermitian form has a conjugate-symmetric Gram matrix `H* = H` (so the diagonal
//! is real). Over the ideal algebraically-closed complexification of a real-closed
//! base, a nondegenerate Hermitian form is classified completely by its
//! **signature** `(p, q)` — Sylvester's law of inertia, the unitary-group
//! `U(p,q)` analogue of the orthogonal signature. The implemented finite-support
//! backend should be read with the same exact-representability caveat as
//! `forms::char0`. We reduce by **unitary (conjugate) congruence**
//! `H ↦ M* H M`, which keeps the form Hermitian and drives it to a real diagonal,
//! then read the signs.

use crate::scalar::{Scalar, Surcomplex};
use std::cmp::Ordering;

/// A Hermitian form, carried by its conjugate-symmetric Gram matrix over
/// `Surcomplex<S>`.
#[derive(Debug, Clone, PartialEq)]
pub struct HermitianForm<S: Scalar> {
    gram: Vec<Vec<Surcomplex<S>>>,
}

/// The signature of a Hermitian form: `(#positive, #negative, #radical)` real
/// diagonal entries after unitary diagonalization. The complete invariant over
/// the surcomplex field.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct HermitianSignature {
    pub pos: usize,
    pub neg: usize,
    pub radical: usize,
}

/// Congruence by the elementary unit `E = I + λ·E_{source,target}`: `H ↦ E* H E`,
/// i.e. `col_target += λ·col_source` then `row_target += conj(λ)·row_source`.
/// Preserves Hermitian-ness.
fn combine<S: Scalar>(
    h: &mut [Vec<Surcomplex<S>>],
    target: usize,
    source: usize,
    lambda: &Surcomplex<S>,
) {
    let n = h.len();
    for r in 0..n {
        let add = lambda.mul(&h[r][source]);
        h[r][target] = h[r][target].add(&add);
    }
    let cl = lambda.conj();
    for c in 0..n {
        let add = cl.mul(&h[source][c]);
        h[target][c] = h[target][c].add(&add);
    }
}

/// Congruence permutation: swap rows `k,i` and columns `k,i`.
fn swap_rows_cols<S: Scalar>(h: &mut [Vec<Surcomplex<S>>], k: usize, i: usize) {
    h.swap(k, i);
    for row in h.iter_mut() {
        row.swap(k, i);
    }
}

/// Make `h[k][k]` a nonzero (real) pivot by congruence, or report that the whole
/// trailing block `[k..]` is zero (radical).
fn ensure_pivot<S: Scalar>(h: &mut [Vec<Surcomplex<S>>], k: usize) -> bool {
    let n = h.len();
    if !h[k][k].is_zero() {
        return true;
    }
    // a nonzero diagonal entry further down → swap it up.
    for i in (k + 1)..n {
        if !h[i][i].is_zero() {
            swap_rows_cols(h, k, i);
            return true;
        }
    }
    // all trailing diagonals zero: combine in an off-diagonal partner. With
    // λ = conj(H[k][j]), the new H[k][k] = H[k][j]·conj(H[k][j]) +
    // conj(H[k][j])·H[k][j] = 2|H[k][j]|² ≠ 0 (real).
    for j in (k + 1)..n {
        if !h[k][j].is_zero() {
            let lambda = h[k][j].conj();
            combine(h, k, j, &lambda);
            return true;
        }
    }
    false // the trailing block is entirely zero
}

impl<S: Scalar> HermitianForm<S> {
    /// Build from a Gram matrix, checking it is square, conjugate-symmetric, and
    /// real on the diagonal. `None` otherwise.
    pub fn from_gram(gram: Vec<Vec<Surcomplex<S>>>) -> Option<Self> {
        let n = gram.len();
        for row in &gram {
            if row.len() != n {
                return None;
            }
        }
        for i in 0..n {
            if !gram[i][i].im.is_zero() {
                return None; // Hermitian diagonal must be real
            }
            for j in 0..n {
                if gram[i][j] != gram[j][i].conj() {
                    return None; // H* = H
                }
            }
        }
        Some(HermitianForm { gram })
    }

    /// Build from a **skew-Hermitian** Gram matrix (`H* = −H`, so the diagonal is
    /// purely imaginary), returning the *Hermitian* form `iH` that classifies it.
    ///
    /// Over a field carrying the conjugation `i ↦ −i`, multiplication by `i` is a
    /// bijection `{skew-Hermitian} → {Hermitian}` (`conj(i) = −i` makes `(iH)* =
    /// iH` exactly when `H* = −H`), so the entire signature machinery transports —
    /// the skew-Hermitian invariant is the signature of `iH`. `None` if the input
    /// is not square and skew-Hermitian.
    pub fn from_skew(gram: Vec<Vec<Surcomplex<S>>>) -> Option<Self> {
        let n = gram.len();
        for row in &gram {
            if row.len() != n {
                return None;
            }
        }
        for i in 0..n {
            if !gram[i][i].re.is_zero() {
                return None; // skew-Hermitian diagonal is purely imaginary
            }
            for j in 0..n {
                if gram[i][j] != gram[j][i].conj().neg() {
                    return None; // H* = −H
                }
            }
        }
        let i_unit = Surcomplex::i();
        let h: Vec<Vec<Surcomplex<S>>> = gram
            .iter()
            .map(|row| row.iter().map(|x| i_unit.mul(x)).collect())
            .collect();
        Self::from_gram(h)
    }

    /// A diagonal Hermitian form from real entries.
    pub fn diagonal(reals: Vec<S>) -> Self {
        let n = reals.len();
        let mut gram = vec![vec![Surcomplex::zero(); n]; n];
        for (i, r) in reals.into_iter().enumerate() {
            gram[i][i] = Surcomplex::new(r, S::zero());
        }
        HermitianForm { gram }
    }

    pub fn dim(&self) -> usize {
        self.gram.len()
    }

    /// The orthogonal direct sum (block-diagonal Gram).
    pub fn direct_sum(&self, other: &HermitianForm<S>) -> HermitianForm<S> {
        let (n, m) = (self.dim(), other.dim());
        let mut gram = vec![vec![Surcomplex::zero(); n + m]; n + m];
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
        HermitianForm { gram }
    }

    /// Unitary (conjugate) congruence to a **real diagonal** — the diagonal
    /// entries (`re`-parts; the `im` parts vanish) whose signs are the signature.
    pub fn diagonalize(&self) -> Vec<S> {
        let n = self.dim();
        let mut h = self.gram.clone();
        for k in 0..n {
            if !ensure_pivot(&mut h, k) {
                continue; // h[k][k] stays 0: a radical direction
            }
            let pinv = h[k][k]
                .inv()
                .expect("nonzero real pivot inverts in a field");
            for i in (k + 1)..n {
                if !h[i][k].is_zero() {
                    let mu = h[k][i].neg().mul(&pinv); // −H[k][i]/H[k][k]
                    combine(&mut h, i, k, &mu);
                }
            }
        }
        (0..n).map(|k| h[k][k].re.clone()).collect()
    }

    /// The Hermitian signature `(pos, neg, radical)`, reading the sign of each
    /// real diagonal entry through `sign` (e.g. `|x| x.sign()` over the surreals
    /// or rationals — the ordered base field). The complete isometry invariant.
    pub fn signature(&self, sign: impl Fn(&S) -> Ordering) -> HermitianSignature {
        let mut sig = HermitianSignature {
            pos: 0,
            neg: 0,
            radical: 0,
        };
        for d in self.diagonalize() {
            match sign(&d) {
                Ordering::Greater => sig.pos += 1,
                Ordering::Less => sig.neg += 1,
                Ordering::Equal => sig.radical += 1,
            }
        }
        sig
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scalar::{Rational, Surreal};

    type GC = Surcomplex<Rational>;

    fn gc(re: i128, im: i128) -> GC {
        Surcomplex::new(Rational::int(re), Rational::int(im))
    }
    fn rsign(x: &Rational) -> Ordering {
        x.sign()
    }

    #[test]
    fn diagonal_real_form_has_sylvester_signature() {
        // ⟨1,1,−1⟩ → (2,1,0); a real-entry Hermitian form is just the symmetric one.
        let h = HermitianForm::<Rational>::diagonal(vec![
            Rational::int(1),
            Rational::int(1),
            Rational::int(-1),
        ]);
        assert_eq!(
            h.signature(rsign),
            HermitianSignature {
                pos: 2,
                neg: 1,
                radical: 0
            }
        );
    }

    #[test]
    fn off_diagonal_hermitian_diagonalizes() {
        // H = [[2, i], [−i, 2]] is Hermitian (H[1][0] = conj(i) = −i), positive
        // definite (det = 4 − 1 = 3 > 0, trace 4 > 0) ⇒ signature (2,0).
        let h = HermitianForm::from_gram(vec![vec![gc(2, 0), gc(0, 1)], vec![gc(0, -1), gc(2, 0)]])
            .unwrap();
        // diagonalizes to [2, 3/2]; both positive.
        assert_eq!(h.diagonalize(), vec![Rational::int(2), Rational::new(3, 2)]);
        assert_eq!(
            h.signature(rsign),
            HermitianSignature {
                pos: 2,
                neg: 0,
                radical: 0
            }
        );
        // a non-Hermitian matrix is rejected.
        assert!(HermitianForm::from_gram(vec![
            vec![gc(2, 0), gc(0, 1)],
            vec![gc(0, 1), gc(2, 0)], // should be −i to be Hermitian
        ])
        .is_none());
    }

    #[test]
    fn off_diagonal_pivot_uses_conjugate_partner() {
        let h = HermitianForm::from_gram(vec![vec![gc(0, 0), gc(1, 1)], vec![gc(1, -1), gc(0, 0)]])
            .unwrap();
        assert_eq!(
            h.diagonalize(),
            vec![Rational::int(4), Rational::new(-1, 2)]
        );
        assert_eq!(
            h.signature(rsign),
            HermitianSignature {
                pos: 1,
                neg: 1,
                radical: 0
            }
        );
    }

    #[test]
    fn indefinite_and_radical() {
        // [[1,0],[0,−1]] → (1,1,0); a zero diagonal entry is radical.
        let h = HermitianForm::from_gram(vec![vec![gc(1, 0), gc(0, 0)], vec![gc(0, 0), gc(-1, 0)]])
            .unwrap();
        assert_eq!(h.signature(rsign).pos, 1);
        assert_eq!(h.signature(rsign).neg, 1);
        let rad = HermitianForm::<Rational>::diagonal(vec![Rational::int(0), Rational::int(5)]);
        assert_eq!(h.direct_sum(&h).signature(rsign).pos, 2); // additive
        assert_eq!(rad.signature(rsign).radical, 1);
    }

    #[test]
    fn skew_hermitian_signature_via_multiplication_by_i() {
        // The real skew-symmetric form [[0,1],[−1,0]] is skew-Hermitian; iH =
        // [[0,i],[−i,0]] is the standard Hermitian form with eigenvalues ±1, so
        // the skew-Hermitian signature is (1,1).
        let h = HermitianForm::<Rational>::from_skew(vec![
            vec![gc(0, 0), gc(1, 0)],
            vec![gc(-1, 0), gc(0, 0)],
        ])
        .unwrap();
        // det(iH) = −1 < 0 ⇒ indefinite ⇒ signature (1,1) (the exact diagonal
        // values depend on the congruence path, only the signs are invariant).
        let sig = h.signature(rsign);
        assert_eq!((sig.pos, sig.neg), (1, 1));
        // a purely-imaginary diagonal is allowed (skew-Hermitian); a real one is not.
        assert!(HermitianForm::<Rational>::from_skew(vec![
            vec![gc(0, 2), gc(0, 0)],
            vec![gc(0, 0), gc(0, -3)],
        ])
        .is_some());
        assert!(HermitianForm::<Rational>::from_skew(vec![
            vec![gc(1, 0), gc(0, 0)],
            vec![gc(0, 0), gc(0, 0)],
        ])
        .is_none());
    }

    #[test]
    fn signature_over_surreal_base() {
        // Hermitian forms over the surreal-complex field, with exact infinite
        // entries: ⟨ω, −ε⟩ Hermitian signature (1,1).
        let h =
            HermitianForm::<Surreal>::diagonal(vec![Surreal::omega(), Surreal::epsilon().neg()]);
        let sig = h.signature(|x| x.sign());
        assert_eq!(sig.pos, 1);
        assert_eq!(sig.neg, 1);
    }
}
