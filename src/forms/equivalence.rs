//! Isometry testing and Witt (hyperbolic) decomposition — the constructive
//! companions to the classifiers.
//!
//! Two quadratic forms over the same field are **isometric** iff they share the
//! complete invariant of that field's leg of the trichotomy: the full signature
//! over the exact-square `Surreal` subdomain, the rank over the exact-square
//! `Surcomplex` subdomain, `(dim, discriminant)` over a finite odd-characteristic
//! field, and the characteristic-2 Arf/radical data over a nim-field. Each
//! `isometric_*` here is exactly that comparison, run on the diagonalized form so
//! it accepts arbitrary (non-diagonal) metrics.
//!
//! **Witt decomposition** writes a form as `k · H ⊥ (anisotropic)`: `k`
//! hyperbolic planes (the Witt index) plus an anisotropic kernel unique up to
//! isometry (Witt's theorem). The anisotropic kernel is the form's class in the
//! Witt group made concrete.

use crate::clifford::Metric;
use crate::forms::char2::{
    arf_nimber_at_degree, min_field_degree, nimber_metric_max_val, ordinal_to_nimber_metric,
};
use crate::forms::{
    arf_char2, arf_fpn_char2, arf_ordinal_finite, as_diagonal, classify_finite_odd,
};
use crate::forms::{FiniteChar2Field, FiniteOddField};
use crate::scalar::{Fpn, Nimber, Ordinal, Rational, Surcomplex, Surreal};

// ----------------------------------------------------------------------------
// Isometry
// ----------------------------------------------------------------------------

/// Are two real (surreal-scalar) forms isometric? `Some(true/false)`, or `None`
/// if either fails to diagonalize.
pub fn isometric_real(m1: &Metric<Surreal>, m2: &Metric<Surreal>) -> Option<bool> {
    let s1 = crate::forms::char0::surreal_signature(m1)?;
    let s2 = crate::forms::char0::surreal_signature(m2)?;
    Some(s1 == s2)
}

/// Are two rational forms isometric by the Hasse--Minkowski invariant package:
/// dimension, discriminant square-class, and Hasse invariant at every place.
pub fn isometric_rational(m1: &Metric<Rational>, m2: &Metric<Rational>) -> Option<bool> {
    Some(crate::forms::classify_rational(m1)? == crate::forms::classify_rational(m2)?)
}

/// Are two surcomplex forms isometric on the exact-square subdomain? Over an
/// algebraically closed field the only invariants are rank and radical dimension.
pub fn isometric_surcomplex(
    m1: &Metric<Surcomplex<Surreal>>,
    m2: &Metric<Surcomplex<Surreal>>,
) -> Option<bool> {
    Some(crate::forms::char0::surcomplex_rank(m1)? == crate::forms::char0::surcomplex_rank(m2)?)
}

/// Are two forms over the same finite odd field isometric? Over a finite field
/// `(dim, discriminant square-class)` is a complete invariant.
pub fn isometric_finite_odd<F: FiniteOddField>(m1: &Metric<F>, m2: &Metric<F>) -> Option<bool> {
    Some(classify_finite_odd(m1)? == classify_finite_odd(m2)?)
}

/// Are two nim-field (characteristic 2) forms isometric? In the nondefective
/// case the Arf invariant of the symplectic complement is part of the invariant.
/// If the polar radical is defective (`Q` nonzero on the radical), adding that
/// radical direction to a symplectic pair toggles the complement's Arf value, so
/// the complement Arf is not an isometry invariant and is deliberately ignored.
///
/// Both Arf invariants are computed using the **same** field degree — the
/// smallest nim-subfield containing all entries of *both* metrics — so that the
/// trace `F_{2^m} → F₂` is consistent.  Computing each independently with its
/// own minimal field degree can yield different Arf bits for isometric forms
/// whose entries span different subfields (the trace of a fixed element differs
/// depending on which extension it is traced from).
pub fn isometric_nimber(m1: &Metric<Nimber>, m2: &Metric<Nimber>) -> Option<bool> {
    let maxv = nimber_metric_max_val(m1).max(nimber_metric_max_val(m2));
    let m = min_field_degree(maxv);
    let a1 = arf_nimber_at_degree(m1, m)?;
    let a2 = arf_nimber_at_degree(m2, m)?;
    Some(same_char2_isometry_invariant(&a1, &a2))
}

/// Are two forms over a supported finite field of characteristic 2 isometric?
/// Same invariant as the nimber path: rank, radical data, and Arf unless the
/// radical is defective.
pub fn isometric_finite_char2<F: FiniteChar2Field>(m1: &Metric<F>, m2: &Metric<F>) -> Option<bool> {
    let a1 = arf_char2(m1)?;
    let a2 = arf_char2(m2)?;
    Some(same_char2_isometry_invariant(&a1, &a2))
}

/// The `Fpn<P,N>` façade helper; returns `None` unless `P = 2`.
pub fn isometric_fpn_char2<const P: u128, const N: usize>(
    m1: &Metric<Fpn<P, N>>,
    m2: &Metric<Fpn<P, N>>,
) -> Option<bool> {
    let a1 = arf_fpn_char2(m1)?;
    let a2 = arf_fpn_char2(m2)?;
    Some(same_char2_isometry_invariant(&a1, &a2))
}

/// Are two supported finite-window ordinal-nimber forms isometric? Returns
/// `None` for ordinal coefficients outside the detected finite subfields.
///
/// For forms whose entries are purely finite ordinals (the inner nimber case),
/// both Arf invariants are computed using the same field degree — the smallest
/// nim-subfield containing entries of *both* metrics — mirroring the
/// `isometric_nimber` consistency guarantee.  The F_64 ordinal-transfinite
/// case uses a fixed six-term trace for both and is unaffected.
pub fn isometric_ordinal_finite(m1: &Metric<Ordinal>, m2: &Metric<Ordinal>) -> Option<bool> {
    // Try the pure-finite (inner Nimber) path: compute both using a common m.
    if let (Some(n1), Some(n2)) = (ordinal_to_nimber_metric(m1), ordinal_to_nimber_metric(m2)) {
        let maxv = nimber_metric_max_val(&n1).max(nimber_metric_max_val(&n2));
        let m = min_field_degree(maxv);
        let a1 = arf_nimber_at_degree(&n1, m)?;
        let a2 = arf_nimber_at_degree(&n2, m)?;
        return Some(same_char2_isometry_invariant(&a1, &a2));
    }
    // Fall back to the independent-trace path (F_64 or cross-case).
    let a1 = arf_ordinal_finite(m1)?;
    let a2 = arf_ordinal_finite(m2)?;
    Some(same_char2_isometry_invariant(&a1, &a2))
}

fn same_char2_isometry_invariant(
    a1: &crate::forms::ArfResult,
    a2: &crate::forms::ArfResult,
) -> bool {
    a1.rank == a2.rank
        && a1.radical_dim == a2.radical_dim
        && a1.radical_anisotropic == a2.radical_anisotropic
        && (a1.radical_anisotropic || a1.arf == a2.arf)
}

// ----------------------------------------------------------------------------
// Witt decomposition
// ----------------------------------------------------------------------------

/// Witt decomposition of a real form: `witt_index` hyperbolic planes plus an
/// anisotropic kernel that is definite (all `+` or all `−`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RealWittDecomp {
    /// Number of hyperbolic planes split off (`min(p, q)`).
    pub witt_index: usize,
    /// `+1` directions in the anisotropic kernel.
    pub anisotropic_pos: usize,
    /// `−1` directions in the anisotropic kernel.
    pub anisotropic_neg: usize,
    /// Dimension of the radical (null directions).
    pub radical_dim: usize,
}

/// Witt decomposition over the exact-square surreal subdomain:
/// `form ≅ k·H ⊥ ⟨±1⟩^{|p−q|}` plus the radical. `k = min(p, q)`.
pub fn witt_decompose_real(m: &Metric<Surreal>) -> Option<RealWittDecomp> {
    let (p, q, r) = crate::forms::char0::surreal_signature(m)?;
    let k = p.min(q);
    Some(RealWittDecomp {
        witt_index: k,
        anisotropic_pos: p - k,
        anisotropic_neg: q - k,
        radical_dim: r,
    })
}

/// Witt decomposition of an odd-characteristic form: `witt_index` hyperbolic
/// planes plus an anisotropic kernel of dimension 0, 1, or 2.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct OddWittDecomp {
    /// Characteristic prime.
    pub p: u128,
    /// Field order `q`; equal to `p` for prime fields and `p^n` for extensions.
    pub field_order: u128,
    /// Number of hyperbolic planes split off.
    pub witt_index: usize,
    /// Dimension of the anisotropic kernel: `0` (hyperbolic), `1` (odd dim), or
    /// `2` (an anisotropic plane).
    pub anisotropic_dim: usize,
    /// Whether the anisotropic kernel's discriminant is a square.
    pub anisotropic_disc_is_square: bool,
    /// Dimension of the radical.
    pub radical_dim: usize,
}

/// Witt decomposition over any finite field `F_q` of odd characteristic: every
/// form of odd dimension has anisotropic kernel `⟨d⟩` (dim 1); an even-dimensional
/// form is hyperbolic (dim 0) iff its discriminant matches the hyperbolic one,
/// else its anisotropic kernel is the unique anisotropic plane (dim 2).
pub fn witt_decompose_finite_odd<F: FiniteOddField>(m: &Metric<F>) -> Option<OddWittDecomp> {
    F::ensure_supported()?;
    let d = as_diagonal(m)?;
    let nonzero: Vec<F> = d.q.into_iter().filter(|x| !x.is_zero()).collect();
    let dim = nonzero.len();
    let radical_dim = m.q.len() - dim;
    let det = nonzero.iter().fold(F::one(), |acc, x| acc.mul(x));

    let anisotropic_dim = if dim % 2 == 1 {
        1
    } else {
        // even dim 2k: hyperbolic iff (−1)^k · det is a square.
        let k = dim / 2;
        let sign = if k % 2 == 1 {
            F::from_i128(-1)
        } else {
            F::one()
        };
        if F::is_square_value(sign.mul(&det)) {
            0
        } else {
            2
        }
    };
    let witt_index = (dim - anisotropic_dim) / 2;
    // anisotropic kernel disc = det · (−1)^{witt_index} (mod squares).
    let twist = if witt_index % 2 == 1 {
        F::from_i128(-1)
    } else {
        F::one()
    };
    Some(OddWittDecomp {
        p: F::characteristic_prime(),
        field_order: F::field_order(),
        witt_index,
        anisotropic_dim,
        anisotropic_disc_is_square: F::is_square_value(det.mul(&twist)),
        radical_dim,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forms::arf_invariant;
    use crate::scalar::Fp;
    use std::collections::BTreeMap;

    fn rsur(xs: &[i128]) -> Metric<Surreal> {
        Metric::diagonal(xs.iter().map(|&x| Surreal::from_int(x)).collect())
    }
    fn ofp<const P: u128>(xs: &[i128]) -> Metric<Fp<P>> {
        Metric::diagonal(xs.iter().map(|&x| Fp::<P>::new(x)).collect())
    }

    #[test]
    fn real_isometry_is_signature_equality() {
        // ⟨1,−1⟩ ≅ ⟨1,−1⟩ but ≇ ⟨1,1⟩.
        assert_eq!(isometric_real(&rsur(&[1, -1]), &rsur(&[-1, 1])), Some(true));
        assert_eq!(isometric_real(&rsur(&[1, -1]), &rsur(&[1, 1])), Some(false));
        // The implemented Surreal backend cannot rescale ⟨2⟩ to ⟨1⟩ exactly.
        assert_eq!(isometric_real(&rsur(&[1]), &rsur(&[2])), None);
        // the skewed hyperbolic plane is isometric to ⟨1,−1⟩.
        let mut b = BTreeMap::new();
        b.insert((0, 1), Surreal::from_int(1));
        let h = Metric::new(vec![Surreal::from_int(0), Surreal::from_int(0)], b);
        assert_eq!(isometric_real(&h, &rsur(&[1, -1])), Some(true));
    }

    #[test]
    fn rational_isometry_sees_square_classes() {
        let q1 = Metric::diagonal(vec![Rational::int(1)]);
        let q2 = Metric::diagonal(vec![Rational::int(2)]);
        assert_eq!(isometric_rational(&q1, &q1), Some(true));
        assert_eq!(isometric_rational(&q1, &q2), Some(false));
    }

    #[test]
    fn real_witt_decomposition_splits_hyperbolics() {
        // ⟨1,1,1,−1,−1⟩ : witt index 2, anisotropic ⟨1⟩.
        let d = witt_decompose_real(&rsur(&[1, 1, 1, -1, -1])).unwrap();
        assert_eq!(
            d,
            RealWittDecomp {
                witt_index: 2,
                anisotropic_pos: 1,
                anisotropic_neg: 0,
                radical_dim: 0,
            }
        );
        // ⟨1,−1⟩ ⊥ radical: pure hyperbolic plane + a null line.
        let d = witt_decompose_real(&rsur(&[1, -1, 0])).unwrap();
        assert_eq!(d.witt_index, 1);
        assert_eq!((d.anisotropic_pos, d.anisotropic_neg), (0, 0));
        assert_eq!(d.radical_dim, 1);
    }

    #[test]
    fn oddchar_isometry_and_witt() {
        const P: u128 = 5;
        // ⟨1,1⟩ ≅ H over F_5 (−1 is a square), so it is hyperbolic, witt index 1.
        let d = witt_decompose_finite_odd(&ofp::<P>(&[1, 1])).unwrap();
        assert_eq!(d.anisotropic_dim, 0);
        assert_eq!(d.witt_index, 1);
        // odd dim ⟨1,1,1⟩: anisotropic kernel dim 1.
        let d = witt_decompose_finite_odd(&ofp::<P>(&[1, 1, 1])).unwrap();
        assert_eq!(d.anisotropic_dim, 1);
        assert_eq!(d.witt_index, 1);
        // isometry by (dim, disc): ⟨1,1⟩ ≅ ⟨2,3⟩? det 1 vs 6=1, both square-class
        assert_eq!(
            isometric_finite_odd(&ofp::<P>(&[1, 1]), &ofp::<P>(&[2, 3])),
            Some(true)
        );
        assert_eq!(
            isometric_finite_odd(&ofp::<P>(&[1, 1]), &ofp::<P>(&[1, 2])),
            Some(false)
        );
    }

    #[test]
    fn oddchar_anisotropic_plane_over_f3() {
        const P: u128 = 3;
        // ⟨1,1⟩ over F_3 is anisotropic (−1 nonsquare): dim-2 kernel, witt index 0.
        let d = witt_decompose_finite_odd(&ofp::<P>(&[1, 1])).unwrap();
        assert_eq!(d.anisotropic_dim, 2);
        assert_eq!(d.witt_index, 0);
    }

    #[test]
    fn nimber_isometry_by_arf() {
        // Over F_2: ⟨1,1⟩ with polar 1 is the anisotropic plane (Arf 1); the
        // hyperbolic plane ⟨0,0⟩ with polar 1 is Arf 0 — not isometric.
        let plane = |q0, q1| {
            let mut b = BTreeMap::new();
            b.insert((0, 1), Nimber(1));
            Metric::new(vec![Nimber(q0), Nimber(q1)], b)
        };
        assert_eq!(isometric_nimber(&plane(1, 1), &plane(1, 1)), Some(true));
        assert_eq!(isometric_nimber(&plane(1, 1), &plane(0, 0)), Some(false));
    }

    // Witness test for M-4: forms over different nim-subfields must be compared
    // using the same field degree for the trace.  Key insight: a form that is
    // anisotropic over F_4 can become isotropic (and hence isometric to the
    // hyperbolic plane) over F_16, because the Artin-Schreier obstruction
    // Tr_{F_16/F_2}(q0·q1) can vanish even when Tr_{F_4/F_2}(q0·q1) = 1.
    //
    // Before the fix, isometric_nimber used each form's own minimal trace,
    // causing the same pair to compare unequal depending on how the form's
    // entries were written — a basis-change invariance failure.
    #[test]
    fn nimber_cross_subfield_isometry_witness() {
        use crate::scalar::nim_mul;

        // Form A (F_4 entries): q=[2,2], b01=1.
        // Tr_{F_4/F_2}(2*2) = Tr_{F_4/F_2}(3) = 3 XOR 2 = 1 → Arf 1 (anisotropic over F_4).
        // But Tr_{F_16/F_2}(3) = 3 XOR 2 XOR 3 XOR 2 = 0 → Arf 0 over F_16.
        // So this form becomes isotropic when viewed over F_16.
        let plane_f4 = {
            let mut b = BTreeMap::new();
            b.insert((0usize, 1usize), Nimber(1));
            Metric::new(vec![Nimber(2), Nimber(2)], b)
        };

        // Form B: apply the basis change diag(α, 1) with α = 4 ∈ F_{16} \ F_4.
        //   q_B[0] = α² * 2 = nim_mul(6, 2),  q_B[1] = 2,  b_B[0,1] = 4.
        // Form B is isometric to A by construction (change of basis over F_16).
        let alpha: u128 = 4;
        let alpha_sq = nim_mul(alpha, alpha); // = 6
        let q_b0 = nim_mul(alpha_sq, 2); // in F_16
        let b_b01 = alpha; // = 4, in F_{16} \ F_4

        assert!(q_b0 >= 4 || b_b01 >= 4, "expected F_16 entries");

        let plane_f16 = {
            let mut b = BTreeMap::new();
            b.insert((0usize, 1usize), Nimber(b_b01));
            Metric::new(vec![Nimber(q_b0), Nimber(2)], b)
        };

        // Standalone arf_invariant uses each form's own minimal field, so the
        // raw Arf bits can differ (Arf=1 for F_4 minimal, Arf=0 for F_16 minimal).
        // This is by-design for the standalone classifier; isometric_nimber must
        // compensate by using the joint field degree.
        let a_f4_standalone = arf_invariant(&plane_f4).unwrap();
        let a_f16_standalone = arf_invariant(&plane_f16).unwrap();
        // They will disagree; record without asserting to document the contrast.
        let _ = (a_f4_standalone.arf, a_f16_standalone.arf);

        // isometric_nimber uses the joint field degree → correctly reports isometric.
        assert_eq!(
            isometric_nimber(&plane_f4, &plane_f16),
            Some(true),
            "isometric forms (related by a basis change) must compare equal"
        );

        // Pure F_2 forms: ⟨1,1⟩ vs ⟨0,0⟩ use joint m=1; should still distinguish.
        let aniso_f2 = {
            let mut b = BTreeMap::new();
            b.insert((0usize, 1usize), Nimber(1));
            Metric::new(vec![Nimber(1), Nimber(1)], b)
        };
        let hyp_f2 = {
            let mut b = BTreeMap::new();
            b.insert((0usize, 1usize), Nimber(1));
            Metric::new(vec![Nimber(0), Nimber(0)], b)
        };
        assert_eq!(
            isometric_nimber(&aniso_f2, &hyp_f2),
            Some(false),
            "same-field anisotropic vs hyperbolic must remain distinguished"
        );

        // plane_f4 viewed jointly with hyp_f2: joint m = max(m_f4, m_f2) = 2.
        // Over F_4, Tr_{F_4/F_2}(3) = 1 → anisotropic.  Hyp has Arf 0 → not isometric.
        assert_eq!(
            isometric_nimber(&plane_f4, &hyp_f2),
            Some(false),
            "F_4 anisotropic plane must not be isometric to F_2 hyperbolic (joint m=2)"
        );

        // plane_f4 vs a hyperbolic plane written with F_16 entries:
        // joint m=4. Tr_{F_16/F_2}(3)=0, so plane_f4 looks hyperbolic at m=4.
        // The F_16 hyperbolic plane also has Arf=0 at m=4. → isometric over F_16.
        let hyp_f16 = {
            let mut b = BTreeMap::new();
            b.insert((0usize, 1usize), Nimber(b_b01)); // b01 = 4 ∈ F_16
            Metric::new(vec![Nimber(0), Nimber(0)], b)
        };
        assert_eq!(
            isometric_nimber(&plane_f4, &hyp_f16),
            Some(true),
            "F_4 anisotropic plane is isometric to the F_16 hyperbolic plane (joint m=4 \
             makes the obstruction vanish)"
        );
    }

    #[test]
    fn defective_radical_ignores_complement_arf() {
        // With a defective radical r (Q(r)=1), replacing both symplectic vectors
        // by a+r and b+r toggles the complement Arf but preserves the whole form.
        let mut b = BTreeMap::new();
        b.insert((0, 1), Nimber(1));
        let split_complement = Metric::new(vec![Nimber(0), Nimber(0), Nimber(1)], b.clone());
        let anisotropic_complement = Metric::new(vec![Nimber(1), Nimber(1), Nimber(1)], b);
        let a1 = arf_invariant(&split_complement).unwrap();
        let a2 = arf_invariant(&anisotropic_complement).unwrap();
        assert_ne!(a1.arf, a2.arf);
        assert!(a1.radical_anisotropic && a2.radical_anisotropic);
        assert_eq!(
            isometric_nimber(&split_complement, &anisotropic_complement),
            Some(true)
        );
    }
}
