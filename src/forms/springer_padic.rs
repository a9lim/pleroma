//! Springer decomposition of a quadratic form over the **p-adic field** `Q_p` —
//! the discretely-valued mirror of [`springer`](crate::forms::springer) (which
//! lives over the surreals).
//!
//! `Q_p` is Henselian with residue field `F_p` and value group `ℤ`. Springer's
//! theorem gives `W(Q_p) ≅ W(F_p) ⊕ (W(F_p) ⊗ ℤ/2ℤ)`. The genuine novelty
//! against the surreal twin: there the value group `No` is **2-divisible**, so
//! the second summand vanishes (`W(No) = W(ℝ)`); here `ℤ/2ℤ ≠ 0`, so **two
//! residue layers survive** — the valuation-even and valuation-odd parts are
//! independent `F_p`-form summands. Scaling an entry by `p²` is a square in
//! `Q_p`, so only the valuation *parity* matters for the Witt class, and the two
//! parities give the two `W(F_p)` copies.
//!
//! This reads the filtration off a [`Qp`](crate::scalar::Qp)-valued diagonal
//! metric: bucket the entries by p-adic valuation, each bucket a residue
//! `F_p`-form recorded by dimension and discriminant square-class. Requires an
//! odd prime `p` (the residue square-class theory) and an already-diagonal metric
//! (`Qp` is a precision model, so we don't congruence-diagonalize over it).

use crate::clifford::Metric;
use crate::forms::is_square;
use crate::scalar::{Fp, Qp};

/// One graded piece of a p-adic Springer decomposition: a residue `F_p`-form at a
/// fixed p-adic valuation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PadicResidueForm {
    /// The p-adic valuation of this graded piece.
    pub valuation: i128,
    /// The residue form's dimension (number of entries at this valuation).
    pub dim: usize,
    /// Whether the residue form's discriminant (product of the residue units mod
    /// `p`) is a square in `F_p` — the `H¹` datum of this layer.
    pub disc_is_square: bool,
}

/// A p-adic Springer decomposition: the valuation-graded residue forms (sorted
/// most-negative-valuation... actually most *infinite* first, i.e. descending),
/// and the radical (genuinely zero entries).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PadicSpringerDecomp {
    pub graded: Vec<PadicResidueForm>,
    pub radical_dim: usize,
}

impl PadicSpringerDecomp {
    /// The residue layers whose valuation has the given parity (`0` = even,
    /// `1` = odd) — the two summands `W(F_p) ⊕ W(F_p)` of `W(Q_p)`.
    pub fn parity_layer(&self, parity: u8) -> Vec<&PadicResidueForm> {
        self.graded
            .iter()
            .filter(|g| (g.valuation.rem_euclid(2) as u8) == parity)
            .collect()
    }
}

/// Decompose a `Q_p` diagonal quadratic form by p-adic valuation. `None` if `p`
/// is not an odd prime, or the metric is non-diagonal.
pub fn springer_decompose_qp<const P: u128, const K: u128>(
    metric: &Metric<Qp<P, K>>,
) -> Option<PadicSpringerDecomp> {
    if P == 2 || !Fp::<P>::modulus_is_prime() {
        return None;
    }
    if !metric.b.is_empty() || metric.has_upper() {
        return None; // already-diagonal only (Qp is a precision model)
    }
    let mut buckets: Vec<(i128, usize, bool)> = Vec::new(); // (valuation, dim, disc_is_square)
    let mut radical_dim = 0usize;
    for x in &metric.q {
        match x.valuation() {
            None => radical_dim += 1, // a genuine zero
            Some(v) => {
                let residue = x.unit() % P; // leading p-adic digit ∈ F_p*
                let sq = is_square::<P>(Fp::<P>(residue));
                match buckets.iter_mut().find(|(bv, _, _)| *bv == v) {
                    Some((_, dim, disc)) => {
                        *dim += 1;
                        *disc = *disc == sq; // square-class is multiplicative (XNOR)
                    }
                    None => buckets.push((v, 1, sq)),
                }
            }
        }
    }
    buckets.sort_by(|a, b| b.0.cmp(&a.0)); // descending valuation
    let graded = buckets
        .into_iter()
        .map(|(valuation, dim, disc_is_square)| PadicResidueForm {
            valuation,
            dim,
            disc_is_square,
        })
        .collect();
    Some(PadicSpringerDecomp {
        graded,
        radical_dim,
    })
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
}
