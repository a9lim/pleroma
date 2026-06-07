//! Concrete spinor modules: left ideals and the matrices that realize Clifford
//! generators as operators on column spinors.
//!
//! The char-0 classifier (`forms::char0`) *names* the real-closed matrix algebra
//! `Cl(p,q) ≅ M_d(K)`; this module builds concrete operator matrices. It searches
//! for an idempotent `f` as a product of commuting "halves" `½(1 + w)` with
//! `w² = +1`, takes the left ideal `S = Cl·f`, picks a basis, and reads off the
//! matrix of left multiplication by each generator on that basis. Those matrices
//! satisfy the Clifford relations `Mᵢ² = qᵢ·I`, `MᵢMⱼ + MⱼMᵢ = bᵢⱼ·I`
//! automatically. When the idempotent search reaches a minimal ideal in the
//! standard real-closed table, its dimension matches `matrix_dim · dim_ℝ(K)`;
//! otherwise the constructor keeps the complete left-regular representation
//! rather than returning an incomplete guess.
//!
//! ## Scope
//!
//! Nondegenerate **orthogonal** metrics (`b`,`a` empty, no null `qᵢ`) in
//! characteristic 0. The constructor first searches for commuting idempotent
//! cuts `½(1+w)` and uses the resulting left ideal when that shrinks the module.
//! If no further explicit cut is found, it still returns the complete
//! left-regular representation (`f = 1`). Degenerate, nonorthogonal,
//! positive-characteristic, or non-enumerable dimensions return `None`.

use crate::clifford::MAX_BASIS_DIM;
use crate::clifford::{bits, CliffordAlgebra, Multivector};
use crate::scalar::Scalar;

/// Explicit spinor matrices grow exponentially (`basis_dim²` entries per
/// generator), so this constructor is intentionally capped instead of pretending
/// a 128-generator representation is materializable.
const MAX_EXPLICIT_SPINOR_DIM: usize = 10;

/// A concrete spinor representation of a Clifford algebra.
pub struct SpinorRep<S: Scalar> {
    /// The idempotent `f` (`f² = f`) generating the represented left ideal.
    pub idempotent: Multivector<S>,
    /// A basis of the left ideal `Cl·f` (in reduced echelon form). If
    /// `is_left_regular` is true, this is the whole algebra.
    pub basis: Vec<Multivector<S>>,
    /// `gen_matrices[i]` is the matrix of left multiplication by `eᵢ` on `basis`
    /// (indexed `[row][col]`; column `j` is the action on `basis[j]`).
    pub gen_matrices: Vec<Vec<Vec<S>>>,
    /// True when the constructor fell back to `f = 1`, i.e. the complete
    /// left-regular representation.
    pub is_left_regular: bool,
}

fn is_idempotent<S: Scalar>(alg: &CliffordAlgebra<S>, f: &Multivector<S>) -> bool {
    &alg.mul(f, f) == f
}

fn commutes<S: Scalar>(alg: &CliffordAlgebra<S>, x: &Multivector<S>, y: &Multivector<S>) -> bool {
    alg.mul(x, y) == alg.mul(y, x)
}

/// Reduced row-echelon basis of the span of `vectors` (each a multivector),
/// keyed by pivot blade-mask, normalized so each pivot coefficient is 1 and is 0
/// in every other basis vector. `None` if a needed pivot is not invertible.
fn rref<S: Scalar>(
    alg: &CliffordAlgebra<S>,
    vectors: &[Multivector<S>],
) -> Option<Vec<(u128, Multivector<S>)>> {
    let mut basis: Vec<(u128, Multivector<S>)> = Vec::new();
    for v in vectors {
        let mut v = v.clone();
        // reduce by existing pivots
        for (p, bvec) in &basis {
            if let Some(c) = v.terms.get(p).cloned() {
                v = alg.add(&v, &alg.scalar_mul(&c.neg(), bvec));
            }
        }
        if v.is_zero() {
            continue;
        }
        let pivot = *v.terms.keys().next().unwrap(); // smallest mask (BTreeMap)
        let lead = v.terms.get(&pivot).cloned().unwrap();
        let linv = lead.inv()?;
        v = alg.scalar_mul(&linv, &v);
        // eliminate this pivot from the existing basis vectors (full reduction)
        for (_, bvec) in &mut basis {
            if let Some(c) = bvec.terms.get(&pivot).cloned() {
                *bvec = alg.add(bvec, &alg.scalar_mul(&c.neg(), &v));
            }
        }
        basis.push((pivot, v));
    }
    basis.sort_by_key(|(p, _)| *p);
    Some(basis)
}

fn blade_count(dim: usize) -> Option<u128> {
    if dim >= MAX_BASIS_DIM || dim > MAX_EXPLICIT_SPINOR_DIM {
        None
    } else {
        Some(1u128 << dim)
    }
}

/// All `blade · f` for blades of the algebra — a spanning set for the left ideal.
fn ideal_spanning_set<S: Scalar>(
    alg: &CliffordAlgebra<S>,
    f: &Multivector<S>,
) -> Option<Vec<Multivector<S>>> {
    let count = blade_count(alg.dim)?;
    Some(
        (0..count)
            .map(|mask| alg.mul(&alg.blade(&bits(mask)), f))
            .collect(),
    )
}

fn ideal_dim<S: Scalar>(alg: &CliffordAlgebra<S>, f: &Multivector<S>) -> usize {
    let Some(spanning) = ideal_spanning_set(alg, f) else {
        return 0;
    };
    rref(alg, &spanning).map(|b| b.len()).unwrap_or(0)
}

/// Coordinates of `target` in a reduced-echelon `basis` (pivot coefficients).
fn coords<S: Scalar>(
    alg: &CliffordAlgebra<S>,
    basis: &[(u128, Multivector<S>)],
    target: &Multivector<S>,
) -> Option<Vec<S>> {
    let coords: Vec<S> = basis
        .iter()
        .map(|(p, _)| target.terms.get(p).cloned().unwrap_or_else(S::zero))
        .collect();
    let mut recon = alg.zero();
    for (c, (_, b)) in coords.iter().zip(basis.iter()) {
        recon = alg.add(&recon, &alg.scalar_mul(c, b));
    }
    if recon == *target {
        Some(coords)
    } else {
        None
    }
}

/// Build a concrete spinor representation. `None` for non-orthogonal,
/// degenerate, positive-characteristic, or non-enumerable metrics (see the
/// module docs).
pub fn spinor_rep<S: Scalar>(alg: &CliffordAlgebra<S>) -> Option<SpinorRep<S>> {
    if !alg.metric.b.is_empty() || !alg.metric.a.is_empty() {
        return None; // orthogonal metrics only
    }
    if S::characteristic() != 0 {
        return None;
    }
    blade_count(alg.dim)?;
    if (0..alg.dim).any(|i| alg.metric.q.get(i).map(|x| x.is_zero()).unwrap_or(true)) {
        return None; // nondegenerate only (no null generators)
    }
    let half = S::one().add(&S::one()).inv()?; // needs ½ (char 0)
    let one = alg.scalar(S::one());

    // Greedily multiply in commuting ½(1+w) factors (w² = +1) while they shrink
    // the represented left ideal. If no cut applies, f=1 gives the full regular
    // representation.
    let mut f = one.clone();
    let mut chosen: Vec<Multivector<S>> = Vec::new();
    let mut cur = ideal_dim(alg, &f);
    loop {
        let mut progressed = false;
        for mask in 1..blade_count(alg.dim)? {
            let w = alg.blade(&bits(mask));
            if alg.mul(&w, &w) != one {
                continue; // need w² = +1
            }
            if !chosen.iter().all(|c| commutes(alg, c, &w)) {
                continue;
            }
            let half_factor = alg.scalar_mul(&half, &alg.add(&one, &w));
            let f2 = alg.mul(&f, &half_factor);
            if !is_idempotent(alg, &f2) {
                continue;
            }
            let d2 = ideal_dim(alg, &f2);
            if d2 < cur {
                f = f2;
                chosen.push(w);
                cur = d2;
                progressed = true;
                break;
            }
        }
        if !progressed {
            break;
        }
    }

    let is_left_regular = f == one;
    let basis = rref(alg, &ideal_spanning_set(alg, &f)?)?;
    let k = basis.len();

    // gen_matrices[i][row][col]: left multiplication by e_i on the basis.
    let mut gen_matrices = vec![vec![vec![S::zero(); k]; k]; alg.dim];
    for i in 0..alg.dim {
        for (col, (_, bvec)) in basis.iter().enumerate() {
            let target = alg.mul(&alg.gen(i), bvec);
            let cs = coords(alg, &basis, &target)?;
            for (row, c) in cs.into_iter().enumerate() {
                gen_matrices[i][row][col] = c;
            }
        }
    }

    let basis_vectors = basis.into_iter().map(|(_, v)| v).collect();
    Some(SpinorRep {
        idempotent: f,
        basis: basis_vectors,
        gen_matrices,
        is_left_regular,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::clifford::Metric;
    use crate::forms::{classify_rational, BaseField};
    use crate::scalar::Rational;

    fn r(n: i128) -> Rational {
        Rational::int(n)
    }

    fn cl(qs: &[i128]) -> CliffordAlgebra<Rational> {
        CliffordAlgebra::new(
            qs.len(),
            Metric::diagonal(qs.iter().map(|&x| r(x)).collect()),
        )
    }

    fn mat_mul(a: &[Vec<Rational>], b: &[Vec<Rational>]) -> Vec<Vec<Rational>> {
        let n = a.len();
        let m = b[0].len();
        let k = b.len();
        let mut out = vec![vec![r(0); m]; n];
        for (i, row) in out.iter_mut().enumerate() {
            for (j, cell) in row.iter_mut().enumerate() {
                let mut acc = r(0);
                for t in 0..k {
                    acc = acc.add(&a[i][t].mul(&b[t][j]));
                }
                *cell = acc;
            }
        }
        out
    }

    fn mat_add(a: &[Vec<Rational>], b: &[Vec<Rational>]) -> Vec<Vec<Rational>> {
        a.iter()
            .zip(b)
            .map(|(ra, rb)| ra.iter().zip(rb).map(|(x, y)| x.add(y)).collect())
            .collect()
    }

    fn scalar_id(s: Rational, n: usize) -> Vec<Vec<Rational>> {
        (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| if i == j { s.clone() } else { r(0) })
                    .collect()
            })
            .collect()
    }

    /// Expected minimal-ideal real dimension = matrix_dim · dim_ℝ(base).
    fn expected_ideal_dim(qs: &[i128]) -> usize {
        let t = classify_rational(&cl(qs).metric).unwrap().real_closure;
        let base = match t.base {
            BaseField::R => 1,
            BaseField::C => 2,
            BaseField::H => 4,
        };
        t.matrix_dim * base
    }

    fn check_clifford_relations(qs: &[i128]) {
        let alg = cl(qs);
        let rep = spinor_rep(&alg).unwrap();
        let k = rep.basis.len();
        assert!(is_idempotent(&alg, &rep.idempotent), "f² ≠ f");
        assert_eq!(
            k,
            expected_ideal_dim(qs),
            "ideal dim ≠ classifier, q={qs:?}"
        );
        // Mᵢ² = qᵢ·I
        for (i, &qi) in qs.iter().enumerate() {
            let mi = &rep.gen_matrices[i];
            assert_eq!(mat_mul(mi, mi), scalar_id(r(qi), k), "M{i}² ≠ q{i}·I");
        }
        // MᵢMⱼ + MⱼMᵢ = 0 (orthogonal, i≠j)
        for i in 0..qs.len() {
            for j in (i + 1)..qs.len() {
                let mi = &rep.gen_matrices[i];
                let mj = &rep.gen_matrices[j];
                let anti = mat_add(&mat_mul(mi, mj), &mat_mul(mj, mi));
                assert_eq!(anti, scalar_id(r(0), k), "{{M{i},M{j}}} ≠ 0");
            }
        }
    }

    #[test]
    fn cl20_spinors_are_two_by_two_real() {
        // Cl(2,0) ≅ M₂(ℝ): minimal left ideal is the 2-dim column space.
        check_clifford_relations(&[1, 1]);
    }

    #[test]
    fn cl30_pauli_spinors() {
        // Cl(3,0) ≅ M₂(ℂ): the Pauli representation, real ideal dim 4.
        check_clifford_relations(&[1, 1, 1]);
    }

    #[test]
    fn cl02_quaternion_spinors() {
        // Cl(0,2) ≅ ℍ: no +1-square blade, so f = 1 and the ideal is all of ℍ;
        // M₀² = M₁² = −I and they anticommute (the quaternion relations).
        check_clifford_relations(&[-1, -1]);
        let alg = cl(&[-1, -1]);
        let rep = spinor_rep(&alg).unwrap();
        assert_eq!(rep.basis.len(), 4);
        assert!(rep.is_left_regular);
        // f = 1 (no idempotent factor was found)
        assert_eq!(rep.idempotent, alg.scalar(r(1)));
    }

    #[test]
    fn cl11_split_spinors() {
        // Cl(1,1) ≅ M₂(ℝ): mixed signature, ideal dim 2.
        check_clifford_relations(&[1, -1]);
    }

    #[test]
    fn cl40_spinors() {
        // Cl(4,0) ≅ M₂(ℍ): ideal real dim = 2·4 = 8.
        check_clifford_relations(&[1, 1, 1, 1]);
    }

    #[test]
    fn degenerate_and_nonorthogonal_are_rejected() {
        // null generator
        assert!(spinor_rep(&cl(&[1, 0])).is_none());
        // non-orthogonal
        let mut b = std::collections::BTreeMap::new();
        b.insert((0usize, 1usize), r(1));
        let alg = CliffordAlgebra::new(2, Metric::new(vec![r(1), r(1)], b));
        assert!(spinor_rep(&alg).is_none());
    }

    #[test]
    fn nonsquare_rational_metrics_get_complete_regular_representation() {
        let alg = cl(&[2]);
        let rep = spinor_rep(&alg).unwrap();
        assert!(rep.is_left_regular);
        assert_eq!(rep.basis.len(), 2);
        let m0 = &rep.gen_matrices[0];
        assert_eq!(mat_mul(m0, m0), scalar_id(r(2), rep.basis.len()));
    }

    #[test]
    fn positive_characteristic_and_non_enumerable_dims_are_rejected() {
        use crate::scalar::Fp;
        let fp_alg = CliffordAlgebra::new(1, Metric::diagonal(vec![Fp::<3>::one()]));
        assert!(spinor_rep(&fp_alg).is_none());

        let large = CliffordAlgebra::new(
            MAX_EXPLICIT_SPINOR_DIM + 1,
            Metric::diagonal(vec![r(1); MAX_EXPLICIT_SPINOR_DIM + 1]),
        );
        assert!(spinor_rep(&large).is_none());
    }
}
