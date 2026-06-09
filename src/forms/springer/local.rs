//! The **generic** discrete-valuation Springer decomposition — one engine, keyed
//! off the [`ResidueField`](crate::scalar::ResidueField) trait, shared by all three
//! discretely-valued legs. The named entry points
//! [`springer_decompose_qp`](crate::forms::springer_decompose_qp),
//! [`springer_decompose_qq`](crate::forms::springer_decompose_qq), and
//! [`springer_decompose_laurent`](crate::forms::springer_decompose_laurent) are thin
//! wrappers over [`springer_decompose_local`].
//!
//! Springer's theorem decomposes a form over a complete discretely-valued field by
//! its valuation filtration. The discrete-valuation legs share **one** structure —
//! a per-valuation residue `k`-form recorded by dimension and discriminant
//! square-class — differing only in the residue field `k`, which the
//! [`ResidueField`](crate::scalar::ResidueField) trait abstracts:
//!
//! | field | char | value group | residue `k` | second layer |
//! |---|---|---|---|---|
//! | `No`        (surreal) | 0 | `No` (2-divisible) | ℝ | collapses — `W(No)=W(ℝ)`, see [`springer_decompose`](crate::forms::springer_decompose) |
//! | `Q_p`       (p-adic)  | 0 | ℤ | `F_p` | survives — `W(Q_p)=W(F_p)²` |
//! | `Q_q`       (unramified p-adic) | 0 | ℤ | `F_q` | survives — `W(Q_q)=W(F_q)²` |
//! | `F_q((t))`  (Laurent) | `p` | ℤ | `F_q` | survives — `W(F_q((t)))=W(F_q)²` |
//!
//! The surreal leg is the one that does **not** fit here, and that is the content
//! of the symmetry, not a gap: its value group `No` is 2-divisible, so the second
//! residue map vanishes and the residue is ℝ (a signature, not a finite
//! square-class). It is `Valued` only in the loose ω-adic sense and is deliberately
//! not a [`ResidueField`], so it keeps its own engine — the surreal
//! [`springer_decompose`](crate::forms::springer_decompose).
//!
//! Among the three that *do* fit, the residue field is `F_p` only for the bare
//! `Q_p`; for `Q_q` (unramified, residue degree `F`) and `F_q((t))` it is a general
//! `F_q`, so the per-layer discriminant square-class genuinely exercises the
//! extension-field square-class. Adding `Q_q` is what makes the mixed-characteristic
//! leg reach general `F_q` residues, matching what the equal-characteristic Laurent
//! leg already did — the two legs are now symmetric in their residue reach.
//!
//! ## The residue-characteristic-2 boundary (honest scope)
//!
//! All three legs require **odd** residue characteristic (the generic engine
//! returns `None` when `K::Residue` is not a supported odd finite field). Springer's
//! second residue map needs residue characteristic `≠ 2`, and a *diagonal* char-2
//! form is totally singular (its polar form vanishes), so the clean
//! `W = W(k) ⊕ W(k)` grading genuinely does not hold there. The char-2 Witt/Arf
//! theory lives in [`char2`](crate::forms::char2), over the full `(q, b)` metric.

use crate::clifford::Metric;
use crate::forms::{is_square_finite, FiniteOddField};
use crate::scalar::ResidueField;

/// One graded piece of a discrete Springer decomposition: a residue `k`-form at a
/// fixed valuation, recorded by dimension and discriminant square-class.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LocalResidueForm {
    /// The valuation of this graded piece.
    pub valuation: i128,
    /// The residue form's dimension (number of entries at this valuation).
    pub dim: usize,
    /// Whether the residue form's discriminant (product of the residue units in
    /// `k`) is a square in `k` — the `H¹` datum of this layer.
    pub disc_is_square: bool,
}

/// A discrete Springer decomposition: the valuation-graded residue forms (sorted
/// descending by valuation), and the radical (genuinely zero entries).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LocalSpringerDecomp {
    pub graded: Vec<LocalResidueForm>,
    pub radical_dim: usize,
}

impl LocalSpringerDecomp {
    /// The residue layers whose valuation has the given parity (`0` = even,
    /// `1` = odd) — the two summands `W(k) ⊕ W(k)` (the value group `ℤ` is not
    /// 2-divisible, so scaling an entry by `ϖ²` is a square but by `ϖ` is not, and
    /// only the valuation parity matters for the Witt class).
    pub fn parity_layer(&self, parity: u8) -> Vec<&LocalResidueForm> {
        self.graded
            .iter()
            .filter(|g| (g.valuation.rem_euclid(2) as u8) == parity)
            .collect()
    }
}

/// Decompose a diagonal quadratic form over any [`ResidueField`] by its valuation
/// filtration: bucket the entries by valuation, each bucket a residue `k`-form
/// recorded by dimension and discriminant square-class. `None` if the residue field
/// is not a supported finite field of odd characteristic, or the metric is
/// non-diagonal (the local backends are precision models, so we do not
/// congruence-diagonalize over them).
pub fn springer_decompose_local<K>(metric: &Metric<K>) -> Option<LocalSpringerDecomp>
where
    K: ResidueField,
    K::Residue: FiniteOddField,
{
    if !K::Residue::is_supported_odd_field() {
        return None; // odd residue characteristic only (see the char-2 boundary)
    }
    if !metric.b.is_empty() || metric.has_upper() {
        return None; // already-diagonal only
    }
    let mut buckets: Vec<(i128, usize, bool)> = Vec::new(); // (valuation, dim, disc_is_square)
    let mut radical_dim = 0usize;
    for x in &metric.q {
        match x.valuation() {
            None => radical_dim += 1, // a genuine zero
            Some(v) => {
                let unit = x
                    .residue_unit()
                    .expect("a nonzero element has an angular component");
                let sq = is_square_finite::<K::Residue>(unit);
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
    buckets.sort_by_key(|x| std::cmp::Reverse(x.0)); // descending valuation
    let graded = buckets
        .into_iter()
        .map(|(valuation, dim, disc_is_square)| LocalResidueForm {
            valuation,
            dim,
            disc_is_square,
        })
        .collect();
    Some(LocalSpringerDecomp {
        graded,
        radical_dim,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scalar::{Fp, Laurent, Qp, Qq, Scalar};

    /// The engine is genuinely generic: the same call decomposes a `Q_p` form and
    /// an `F_q((t))` form, reading each one's residue field through the trait.
    #[test]
    fn one_engine_decomposes_every_discrete_leg() {
        // ⟨1, 5⟩ over Q_5: two valuation layers, both residue-square.
        let qp = Metric::diagonal(vec![Qp::<5, 4>::from_i128(1), Qp::<5, 4>::from_i128(5)]);
        let dp = springer_decompose_local(&qp).unwrap();
        assert_eq!(dp.graded.len(), 2);
        assert_eq!(dp.parity_layer(0).len(), 1);
        assert_eq!(dp.parity_layer(1).len(), 1);

        // ⟨1, t⟩ over F_5((t)): the mirror, equal characteristic.
        let lt = Metric::diagonal(vec![
            Laurent::<Fp<5>, 4>::from_coeffs(vec![Fp::<5>::new(1)], 0),
            Laurent::<Fp<5>, 4>::from_coeffs(vec![Fp::<5>::new(1)], 1),
        ]);
        let dl = springer_decompose_local(&lt).unwrap();
        assert_eq!(dl.graded.len(), 2);
        // structurally identical decompositions from one engine.
        assert_eq!(dp.graded, dl.graded);
    }

    /// Residue characteristic 2 is rejected on every leg (the Springer boundary).
    #[test]
    fn residue_char_two_is_rejected_uniformly() {
        assert!(springer_decompose_local(&Metric::diagonal(vec![Qp::<2, 4>::one()])).is_none());
        assert!(springer_decompose_local(&Metric::diagonal(vec![Laurent::<
            crate::scalar::Fpn<2, 3>,
            4,
        >::one()]))
        .is_none());
    }

    /// The new sibling: `Q_q` (unramified, residue degree 2) reads its square-class
    /// in `F_9`, not `F_3` — content invisible to a bare-`Q_p` decomposition.
    #[test]
    fn unramified_qq_reads_extension_residue() {
        use crate::scalar::{Fpn, WittVec};
        type Q9 = Qq<3, 3, 2>;
        let ns = (0..9u128)
            .map(|c| Fpn::<3, 2>::from_coeffs(&[c % 3, c / 3]))
            .find(|x| !x.is_zero() && !x.is_square())
            .expect("F_9 has nonsquares");
        // ⟨ns·p, ns²⟩: valuation 1 carries ns (a nonsquare), valuation 0 carries ns² (a
        // square). The naive lift `WittVec(x.0)` is a Witt unit with residue x.
        let m = Metric::diagonal(vec![
            Q9::from_witt(WittVec::<3, 3, 2>(ns.into_coeffs())).mul(&Q9::from_int(3)),
            Q9::from_witt(WittVec::<3, 3, 2>(ns.mul(&ns).into_coeffs())),
        ]);
        let d = springer_decompose_local(&m).unwrap();
        assert_eq!(d.graded.len(), 2);
        assert_eq!(d.graded[0].valuation, 1);
        assert!(!d.graded[0].disc_is_square, "ns is a nonsquare in F_9");
        assert!(d.graded[1].disc_is_square, "ns² is a square in F_9");
    }
}
