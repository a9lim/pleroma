//! The Dickson invariant: the characteristic-2 replacement for determinant on
//! orthogonal groups.

use crate::clifford::{CliffordAlgebra, Multivector};
use crate::linalg::f2;
use crate::scalar::{nim_add, Nimber};

/// The **Dickson invariant** `D(g) ∈ F₂` of an orthogonal transformation `g`,
/// given as an n×n matrix over a nim-field: `D(g) = dim Im(g − I) mod 2`
/// (`= rank(g + I) mod 2`, since `−1 = 1`).
///
/// In characteristic 2 the determinant of any `g ∈ O(Q)` is forced to `1`, so it
/// cannot separate rotations from reflections — the Dickson invariant is the
/// replacement, with `SO(Q) = ker D`. A single reflection has `D = 1`; a product
/// of `k` reflections has `D = k mod 2`. It is the companion to the Arf
/// invariant: **Arf classifies the form, Dickson classifies `O(Q)`.**
pub fn dickson_matrix(g: &[Vec<u128>]) -> u8 {
    let n = g.len();
    let mut m: Vec<Vec<u128>> = g.to_vec();
    for i in 0..n {
        m[i][i] = nim_add(m[i][i], 1); // g − I  (= g + I in char 2)
    }
    (f2::nim_rank(m) % 2) as u8
}

/// The Dickson invariant of a Clifford **versor** (a product of vectors) acting
/// by the twisted adjoint: it is the ℤ₂-grade parity of the versor — an even
/// versor (rotor) lies in `SO` with `D = 0`, an odd versor (e.g. a single vector,
/// a reflection) has `D = 1`. Returns `None` if the multivector is not of
/// homogeneous grade parity (hence not a versor) or is zero.
pub fn dickson_of_versor(alg: &CliffordAlgebra<Nimber>, v: &Multivector<Nimber>) -> Option<u8> {
    // The Dickson invariant of a versor is its grade parity, a fact independent of
    // the scalar field — so this is the char-2 specialisation of the generic
    // `clifford::versor_grade_parity`.
    let dickson = crate::clifford::versor_grade_parity(v)?;
    alg.versor_inverse(v)?;
    Some(dickson)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::clifford::{CliffordAlgebra, Metric};
    use crate::scalar::nim_mul;

    #[test]
    fn dickson_separates_rotations_from_reflections() {
        // identity is a rotation: D = 0.
        assert_eq!(dickson_matrix(&[vec![1, 0], vec![0, 1]]), 0);
        // the swap (0 1; 1 0) preserves the hyperbolic form x0 x1 and is a
        // reflection (odd): D = 1.
        assert_eq!(dickson_matrix(&[vec![0, 1], vec![1, 0]]), 1);
        // a hyperbolic "rotation" diag(t, t⁻¹) preserves x0 x1; for t=*2 in F₄,
        // t⁻¹ = *3, so g = diag(*2,*3): D = 0 (in SO).
        assert_eq!(dickson_matrix(&[vec![2, 0], vec![0, 3]]), 0);
        // composing two reflections (here swap∘swap = identity) gives D = 0.
        let swap = [[0u128, 1], [1, 0]];
        let mut comp = vec![vec![0u128; 2]; 2];
        for i in 0..2 {
            for j in 0..2 {
                let mut acc = 0u128;
                for k in 0..2 {
                    acc ^= nim_mul(swap[i][k], swap[k][j]);
                }
                comp[i][j] = acc;
            }
        }
        assert_eq!(dickson_matrix(&comp), 0);
    }

    #[test]
    fn dickson_of_versor_is_grade_parity() {
        let alg = CliffordAlgebra::new(2, Metric::diagonal(vec![Nimber(1), Nimber(1)]));
        let scalar_one = alg.scalar(Nimber(1));
        let e0 = alg.gen(0);
        let e0e1 = alg.mul(&alg.gen(0), &alg.gen(1));
        assert_eq!(dickson_of_versor(&alg, &scalar_one), Some(0)); // identity rotor
        assert_eq!(dickson_of_versor(&alg, &e0), Some(1)); // a vector = a reflection
        assert_eq!(dickson_of_versor(&alg, &e0e1), Some(0)); // a bivector = a rotor
                                                             // mixed parity ⇒ not a versor
        let mixed = alg.add(&e0, &e0e1);
        assert_eq!(dickson_of_versor(&alg, &mixed), None);

        let null_alg = CliffordAlgebra::new(1, Metric::grassmann(1));
        assert_eq!(dickson_of_versor(&null_alg, &null_alg.gen(0)), None);
    }
}
