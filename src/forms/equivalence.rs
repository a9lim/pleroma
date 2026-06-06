//! Isometry testing and Witt (hyperbolic) decomposition — the constructive
//! companions to the classifiers.
//!
//! Two quadratic forms over the same field are **isometric** iff they share the
//! complete invariant of that field's leg of the trichotomy: the full signature
//! over a real-closed field, the rank over an algebraically closed one, `(dim,
//! discriminant)` over a finite odd-characteristic field, and the Arf data over
//! a nim-field. Each `isometric_*` here is exactly that comparison, run on the
//! diagonalized form so it accepts arbitrary (non-diagonal) metrics.
//!
//! **Witt decomposition** writes a form as `k · H ⊥ (anisotropic)`: `k`
//! hyperbolic planes (the Witt index) plus an anisotropic kernel unique up to
//! isometry (Witt's theorem). The anisotropic kernel is the form's class in the
//! Witt group made concrete.

use crate::clifford::Metric;
use crate::forms::{arf_invariant, as_diagonal, classify_oddchar};
use crate::scalar::{Fp, Nimber, Rational, Scalar, Surcomplex, Surreal};

// ----------------------------------------------------------------------------
// Isometry
// ----------------------------------------------------------------------------

/// Are two real (surreal-scalar) forms isometric? `Some(true/false)`, or `None`
/// if either fails to diagonalize.
pub fn isometric_real(m1: &Metric<Surreal>, m2: &Metric<Surreal>) -> Option<bool> {
    let s1 = crate::forms::char0::signature(m1, |x| x.sign())?;
    let s2 = crate::forms::char0::signature(m2, |x| x.sign())?;
    Some(s1 == s2)
}

/// Are two rational forms isometric *as real forms* (same signature)? Rational
/// isometry is finer (it would track square classes too); this is the
/// real-closed comparison the surreal backend realises.
pub fn isometric_rational(m1: &Metric<Rational>, m2: &Metric<Rational>) -> Option<bool> {
    let s1 = crate::forms::char0::signature(m1, |x| x.sign())?;
    let s2 = crate::forms::char0::signature(m2, |x| x.sign())?;
    Some(s1 == s2)
}

/// Are two surcomplex forms isometric? Over an algebraically closed field the
/// only invariants are rank and radical dimension.
pub fn isometric_surcomplex(
    m1: &Metric<Surcomplex<Surreal>>,
    m2: &Metric<Surcomplex<Surreal>>,
) -> Option<bool> {
    let rank = |m: &Metric<Surcomplex<Surreal>>| -> Option<(usize, usize)> {
        let d = as_diagonal(m)?;
        let nz = d.q.iter().filter(|z| !z.is_zero()).count();
        Some((nz, d.q.len() - nz))
    };
    Some(rank(m1)? == rank(m2)?)
}

/// Are two odd-characteristic forms isometric? Over a finite field `(dim,
/// discriminant square-class)` is a complete invariant.
pub fn isometric_oddchar<const P: u128>(m1: &Metric<Fp<P>>, m2: &Metric<Fp<P>>) -> Option<bool> {
    Some(classify_oddchar(m1)? == classify_oddchar(m2)?)
}

/// Are two nim-field (characteristic 2) forms isometric? The Arf data
/// `(arf, rank, radical_dim, radical_anisotropic)` is the complete invariant.
pub fn isometric_nimber(m1: &Metric<Nimber>, m2: &Metric<Nimber>) -> Option<bool> {
    let a1 = arf_invariant(m1)?;
    let a2 = arf_invariant(m2)?;
    Some(a1 == a2)
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

/// Witt decomposition over the real-closed surreals: `form ≅ k·H ⊥ ⟨±1⟩^{|p−q|}`
/// plus the radical. `k = min(p, q)`.
pub fn witt_decompose_real(m: &Metric<Surreal>) -> Option<RealWittDecomp> {
    let (p, q, r) = crate::forms::char0::signature(m, |x| x.sign())?;
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
    pub p: u128,
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

/// Witt decomposition over a finite field `F_q` (odd `q`): every form of odd
/// dimension has anisotropic kernel `⟨d⟩` (dim 1); an even-dimensional form is
/// hyperbolic (dim 0) iff its discriminant matches the hyperbolic one, else its
/// anisotropic kernel is the unique anisotropic plane (dim 2).
pub fn witt_decompose_oddchar<const P: u128>(m: &Metric<Fp<P>>) -> Option<OddWittDecomp> {
    let d = as_diagonal(m)?;
    let nonzero: Vec<Fp<P>> = d.q.into_iter().filter(|x| !x.is_zero()).collect();
    let dim = nonzero.len();
    let radical_dim = m.q.len() - dim;
    let det = nonzero.iter().fold(Fp::<P>::one(), |acc, x| acc.mul(x));

    let anisotropic_dim = if dim % 2 == 1 {
        1
    } else {
        // even dim 2k: hyperbolic iff (−1)^k · det is a square.
        let k = dim / 2;
        let sign = if k % 2 == 1 {
            Fp::<P>::new(-1)
        } else {
            Fp::<P>::one()
        };
        if super::is_square(sign.mul(&det)) {
            0
        } else {
            2
        }
    };
    let witt_index = (dim - anisotropic_dim) / 2;
    // anisotropic kernel disc = det · (−1)^{witt_index} (mod squares).
    let twist = if witt_index % 2 == 1 {
        Fp::<P>::new(-1)
    } else {
        Fp::<P>::one()
    };
    Some(OddWittDecomp {
        p: P,
        witt_index,
        anisotropic_dim,
        anisotropic_disc_is_square: super::is_square(det.mul(&twist)),
        radical_dim,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
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
        // the skewed hyperbolic plane is isometric to ⟨1,−1⟩.
        let mut b = BTreeMap::new();
        b.insert((0, 1), Surreal::from_int(2));
        let h = Metric::new(vec![Surreal::from_int(0), Surreal::from_int(0)], b);
        assert_eq!(isometric_real(&h, &rsur(&[1, -1])), Some(true));
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
        let d = witt_decompose_oddchar(&ofp::<P>(&[1, 1])).unwrap();
        assert_eq!(d.anisotropic_dim, 0);
        assert_eq!(d.witt_index, 1);
        // odd dim ⟨1,1,1⟩: anisotropic kernel dim 1.
        let d = witt_decompose_oddchar(&ofp::<P>(&[1, 1, 1])).unwrap();
        assert_eq!(d.anisotropic_dim, 1);
        assert_eq!(d.witt_index, 1);
        // isometry by (dim, disc): ⟨1,1⟩ ≅ ⟨2,3⟩? det 1 vs 6=1, both square-class
        assert_eq!(
            isometric_oddchar(&ofp::<P>(&[1, 1]), &ofp::<P>(&[2, 3])),
            Some(true)
        );
        assert_eq!(
            isometric_oddchar(&ofp::<P>(&[1, 1]), &ofp::<P>(&[1, 2])),
            Some(false)
        );
    }

    #[test]
    fn oddchar_anisotropic_plane_over_f3() {
        const P: u128 = 3;
        // ⟨1,1⟩ over F_3 is anisotropic (−1 nonsquare): dim-2 kernel, witt index 0.
        let d = witt_decompose_oddchar(&ofp::<P>(&[1, 1])).unwrap();
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
}
