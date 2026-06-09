//! Springer decomposition over the **mixed-characteristic** local fields `Q_p`
//! (residue `F_p`) and its unramified extensions `Q_q` (residue `F_q`) — two named
//! entry points into the generic engine
//! [`springer_decompose_local`](crate::forms::springer_decompose_local), and the
//! discretely-valued mirror of the surreal leg
//! ([`springer_decompose`](crate::forms::springer_decompose)).
//!
//! `Q_p` is Henselian with residue field `F_p` and value group `ℤ`. Springer's
//! theorem gives `W(Q_p) ≅ W(F_p) ⊕ (W(F_p) ⊗ ℤ/2ℤ)`. The genuine novelty against
//! the surreal twin: there the value group `No` is **2-divisible**, so the second
//! summand vanishes (`W(No) = W(ℝ)`); here `ℤ/2ℤ ≠ 0`, so **two residue layers
//! survive** — the valuation-even and valuation-odd parts are independent residue
//! summands. Scaling an entry by `p²` is a square in `Q_p`, so only the valuation
//! *parity* matters for the Witt class, and the two parities give the two `W(F_p)`
//! copies ([`LocalSpringerDecomp::parity_layer`]).
//!
//! `Q_q = Frac(W_N(F_q))`, the unramified extension of residue degree `F`, is the
//! same story with residue field `F_q = F_{p^F}` in place of `F_p` — and `Q_q` with
//! `F = 1` *is* `Q_p`. Adding it is what makes the mixed-characteristic leg reach
//! general `F_q` residues, matching the equal-characteristic Laurent leg
//! ([`springer_decompose_laurent`](crate::forms::springer_decompose_laurent)) which
//! already did: the
//! per-layer discriminant square-class then lives in `F_q*/(F_q*)²` and genuinely
//! exercises the extension-field square-class, not just `F_p`.
//!
//! Both read the filtration off a diagonal metric, bucketing entries by valuation;
//! both require an **odd** residue characteristic (`p = 2`, resp. residue char 2,
//! returns `None`) and an already-diagonal metric (`Qp`/`Qq` are precision models,
//! so we do not congruence-diagonalize over them). The residue-characteristic-2
//! boundary is the same one documented on the generic engine and the Laurent
//! sibling.

use crate::clifford::Metric;
use crate::scalar::{Qp, Qq};

use super::local::springer_decompose_local;
pub use super::local::{
    LocalResidueForm as PadicResidueForm, LocalSpringerDecomp as PadicSpringerDecomp,
};

/// Decompose a `Q_p` diagonal quadratic form by p-adic valuation. `None` if `p` is
/// not an odd prime, or the metric is non-diagonal. A thin wrapper over
/// [`springer_decompose_local`] (residue field `F_p`).
pub fn springer_decompose_qp<const P: u128, const K: u128>(
    metric: &Metric<Qp<P, K>>,
) -> Option<PadicSpringerDecomp> {
    springer_decompose_local(metric)
}

/// Decompose a `Q_q` diagonal quadratic form by p-adic valuation, reading the
/// per-layer square-class in the residue field `F_q = F_{p^F}`. `None` if the
/// residue characteristic is `2` or the residue field is unsupported, or the metric
/// is non-diagonal. The unramified generalization of [`springer_decompose_qp`]
/// (`F = 1` recovers it).
pub fn springer_decompose_qq<const P: u128, const N: usize, const F: usize>(
    metric: &Metric<Qq<P, N, F>>,
) -> Option<PadicSpringerDecomp> {
    springer_decompose_local(metric)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scalar::Scalar;

    type Q5 = Qp<5, 4>;

    #[test]
    fn two_residue_layers_survive() {
        // ⟨1, 5⟩ over Q_5: valuation 0 (unit 1) and valuation 1 (unit 1) — the two
        // layers do NOT collapse (unlike the 2-divisible surreal value group).
        let m = Metric::diagonal(vec![Q5::from_i128(1), Q5::from_i128(5)]);
        let d = springer_decompose_qp(&m).unwrap();
        assert_eq!(d.graded.len(), 2);
        assert_eq!(d.graded[0].valuation, 1); // sorted descending
        assert_eq!(d.graded[1].valuation, 0);
        // dims + radical recover the total dimension.
        assert_eq!(
            d.graded.iter().map(|g| g.dim).sum::<usize>() + d.radical_dim,
            2
        );
        // the even-valuation and odd-valuation layers are the two W(F_5) summands.
        assert_eq!(d.parity_layer(0).len(), 1); // valuation 0
        assert_eq!(d.parity_layer(1).len(), 1); // valuation 1
                                                // residue unit 1 is a square in F_5.
        assert!(d.graded[0].disc_is_square);
        assert!(d.graded[1].disc_is_square);
    }

    #[test]
    fn residue_square_class_tracks_nonsquares() {
        // ⟨2, 10⟩ over Q_5: residues 2 (val 0) and 2 (val 1), both nonsquares.
        let m = Metric::diagonal(vec![Q5::from_i128(2), Q5::from_i128(10)]);
        let d = springer_decompose_qp(&m).unwrap();
        assert_eq!(d.graded.len(), 2);
        assert!(!d.graded[0].disc_is_square); // residue 2 is a nonsquare mod 5
        assert!(!d.graded[1].disc_is_square);
        // ⟨2, 3⟩ at valuation 0: disc = 2·3 = 6 ≡ 1, a square ⇒ the layer is square.
        let m2 = Metric::diagonal(vec![Q5::from_i128(2), Q5::from_i128(3)]);
        let d2 = springer_decompose_qp(&m2).unwrap();
        assert_eq!(d2.graded.len(), 1);
        assert_eq!(d2.graded[0].dim, 2);
        assert!(d2.graded[0].disc_is_square); // 2·3 = 6 ≡ 1 mod 5
    }

    #[test]
    fn radical_and_rejections() {
        // a genuine zero entry is radical, not a residue layer.
        let m = Metric::diagonal(vec![Q5::zero(), Q5::from_i128(5)]);
        let d = springer_decompose_qp(&m).unwrap();
        assert_eq!(d.radical_dim, 1);
        assert_eq!(d.graded.len(), 1);
        // p = 2 is rejected (residue square-class theory needs odd p).
        assert!(springer_decompose_qp(&Metric::diagonal(vec![Qp::<2, 4>::from_i128(1)])).is_none());
    }

    #[test]
    fn unramified_qq_recovers_qp_when_residue_degree_is_one() {
        // Q_q with F = 1 IS Q_p — the decomposition must agree with the Q_p path.
        type Q5Unram = Qq<5, 4, 1>;
        let mqq = Metric::diagonal(vec![Q5Unram::from_int(1), Q5Unram::from_int(5)]);
        let dqq = springer_decompose_qq(&mqq).unwrap();
        let mqp = Metric::diagonal(vec![Q5::from_i128(1), Q5::from_i128(5)]);
        let dqp = springer_decompose_qp(&mqp).unwrap();
        assert_eq!(dqq, dqp);
    }

    #[test]
    fn unramified_qq_reads_f9_square_class() {
        // The genuine extension-residue content: residue field F_9, square-class in
        // F_9*/(F_9*)², invisible to a bare-Q_p (F_3) decomposition.
        use crate::scalar::{Fpn, WittVec};
        type Q9 = Qq<3, 3, 2>;
        let ns = (0..9u128)
            .map(|c| Fpn::<3, 2>::from_coeffs(&[c % 3, c / 3]))
            .find(|x| !x.is_zero() && !x.is_square())
            .expect("F_9 has nonsquares");
        // ⟨ns², ns·p⟩: val 0 carries ns² (square), val 1 carries ns (nonsquare).
        let m = Metric::diagonal(vec![
            Q9::from_witt(WittVec::<3, 3, 2>(ns.mul(&ns).into_coeffs())),
            Q9::from_witt(WittVec::<3, 3, 2>(ns.into_coeffs())).mul(&Q9::from_int(3)),
        ]);
        let d = springer_decompose_qq(&m).unwrap();
        assert_eq!(d.graded.len(), 2);
        assert_eq!(d.graded[0].valuation, 1); // descending
        assert!(!d.graded[0].disc_is_square, "ns is a nonsquare in F_9");
        assert!(d.graded[1].disc_is_square, "ns² is a square in F_9");
        // residue characteristic 2 (Q_q over F_{2^…}) is rejected.
        assert!(
            springer_decompose_qq(&Metric::diagonal(vec![Qq::<2, 4, 2>::from_int(1)])).is_none()
        );
    }
}
