//! Springer decomposition over the **equal-characteristic** local field `F_q((t))`
//! — the third discretely-valued sibling, beside the surreal leg
//! ([`springer_decompose`](crate::forms::springer_decompose)) and the `Q_p`/`Q_q`
//! leg ([`springer_decompose_qp`](crate::forms::springer_decompose_qp)), and a named
//! entry point into the generic engine
//! [`springer_decompose_local`](crate::forms::springer_decompose_local). See that
//! engine for the full discrete-valuation trichotomy table.
//!
//! The novelty against the `Q_p` twin is twofold: this is the **equal
//! characteristic** case (the field itself has characteristic `p`, not `0`), and the
//! residue field is a general **`F_q = F_{p^n}`** — the same general-residue reach
//! the `Q_q` sibling gives the mixed-characteristic leg. The per-layer discriminant
//! square-class lives in `F_q*/(F_q*)²` and genuinely exercises the extension-field
//! square-class (`Fpn::is_square`), not just `F_p`.
//!
//! Like the `Q_p`/`Q_q` siblings: the value group `ℤ` is **not** 2-divisible, so
//! scaling an entry by `t²` is a square but scaling by `t` is not — only the
//! valuation *parity* matters for the Witt class, and the two parities give the two
//! `W(F_q)` summands ([`LocalSpringerDecomp::parity_layer`]). Requires an **odd**
//! residue characteristic and an already-diagonal metric (`Laurent` is a precision
//! model).
//!
//! ## The residue-characteristic-2 boundary (honest scope)
//!
//! Residue characteristic 2 — `F_{2^n}((t))` — is **rejected** (returns `None`),
//! exactly as the `Q_p`/`Q_q` siblings reject residue char 2. This is not laziness:
//! Springer's second residue map requires residue characteristic `≠ 2`, and a
//! *diagonal* char-2 form is totally singular (the polar form vanishes), so the
//! clean `W(F_q((t))) = W(F_q) ⊕ W(F_q)` grading genuinely does not hold there. The
//! char-2 Witt/Arf theory lives in [`char2`](crate::forms::char2), over the full
//! `(q, b)` metric, not through this valuation filtration.

use crate::clifford::Metric;
use crate::forms::FiniteOddField;
use crate::scalar::Laurent;

use super::local::springer_decompose_local;
pub use super::local::{
    LocalResidueForm as LaurentResidueForm, LocalSpringerDecomp as LaurentSpringerDecomp,
};

/// Decompose an `F_q((t))` diagonal quadratic form by `t`-adic valuation. `None` if
/// the residue field is not a supported finite field of odd characteristic, or the
/// metric is non-diagonal. A thin wrapper over [`springer_decompose_local`] (residue
/// field `S = F_q`).
pub fn springer_decompose_laurent<S: FiniteOddField, const K: usize>(
    metric: &Metric<Laurent<S, K>>,
) -> Option<LaurentSpringerDecomp> {
    springer_decompose_local(metric)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scalar::{Fp, Fpn, Scalar};

    type L5 = Laurent<Fp<5>, 4>; // F_5((t))
    type L9 = Laurent<Fpn<3, 2>, 4>; // F_9((t)) — genuine extension residue

    fn l5(coeffs: &[i128], val: i128) -> L5 {
        Laurent::from_coeffs(coeffs.iter().map(|&n| Fp::<5>::new(n)).collect(), val)
    }

    #[test]
    fn two_residue_layers_survive() {
        // ⟨1, t⟩ over F_5((t)): valuation 0 and valuation 1 — the two layers do NOT
        // collapse (value group ℤ is not 2-divisible), mirroring the Q_p sibling.
        let m = Metric::diagonal(vec![l5(&[1], 0), l5(&[1], 1)]);
        let d = springer_decompose_laurent(&m).unwrap();
        assert_eq!(d.graded.len(), 2);
        assert_eq!(d.graded[0].valuation, 1); // descending
        assert_eq!(d.graded[1].valuation, 0);
        assert_eq!(
            d.graded.iter().map(|g| g.dim).sum::<usize>() + d.radical_dim,
            2
        );
        assert_eq!(d.parity_layer(0).len(), 1); // valuation 0
        assert_eq!(d.parity_layer(1).len(), 1); // valuation 1
        assert!(d.graded[0].disc_is_square); // residue 1 is a square
        assert!(d.graded[1].disc_is_square);
    }

    #[test]
    fn residue_square_class_tracks_nonsquares_in_f5() {
        // ⟨2, 2t⟩: residues 2 (val 0) and 2 (val 1), both nonsquares mod 5.
        let m = Metric::diagonal(vec![l5(&[2], 0), l5(&[2], 1)]);
        let d = springer_decompose_laurent(&m).unwrap();
        assert!(!d.graded[0].disc_is_square);
        assert!(!d.graded[1].disc_is_square);
        // ⟨2, 3⟩ at valuation 0: disc = 2·3 = 6 ≡ 1, a square ⇒ the layer is square.
        let m2 = Metric::diagonal(vec![l5(&[2], 0), l5(&[3], 0)]);
        let d2 = springer_decompose_laurent(&m2).unwrap();
        assert_eq!(d2.graded.len(), 1);
        assert_eq!(d2.graded[0].dim, 2);
        assert!(d2.graded[0].disc_is_square);
    }

    #[test]
    fn extension_residue_field_f9_square_class() {
        // The genuine generalization over the p-adic sibling: the residue field is
        // F_9, so the square-class is computed in F_9*/(F_9*)², not F_3. Take an
        // actual F_9 nonsquare `ns`: `ns²` is a square, `ns` is not — content
        // invisible to an F_p-only decomposition.
        let ns = (0..9u128)
            .map(|c| Fpn::<3, 2>::from_coeffs(&[c % 3, c / 3]))
            .find(|x| !x.is_zero() && !x.is_square())
            .expect("F_9 has nonsquares");
        let sq = ns.mul(&ns); // a square by construction
        let m = Metric::diagonal(vec![
            Laurent::<Fpn<3, 2>, 4>::from_coeffs(vec![sq], 1),
            Laurent::<Fpn<3, 2>, 4>::from_coeffs(vec![ns], 0),
        ]);
        let d = springer_decompose_laurent(&m).unwrap();
        assert_eq!(d.graded.len(), 2);
        assert_eq!(d.graded[0].valuation, 1); // descending
        assert!(d.graded[0].disc_is_square, "ns² is a square in F_9");
        assert!(!d.graded[1].disc_is_square, "ns is a nonsquare in F_9");
        let _ = L9::one(); // (type alias exercised)
    }

    #[test]
    fn radical_and_char2_rejection() {
        // a genuine zero entry is radical, not a residue layer.
        let m = Metric::diagonal(vec![L5::zero(), l5(&[1], 1)]);
        let d = springer_decompose_laurent(&m).unwrap();
        assert_eq!(d.radical_dim, 1);
        assert_eq!(d.graded.len(), 1);
        // residue characteristic 2 (F_8((t))) is rejected — the Springer boundary.
        let m2 = Metric::diagonal(vec![Laurent::<Fpn<2, 3>, 4>::one()]);
        assert!(springer_decompose_laurent(&m2).is_none());
    }
}
