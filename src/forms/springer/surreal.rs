//! Springer-style decomposition of a diagonal quadratic form over the implemented
//! **surreal** Hahn/CNF backend — the non-Archimedean valuation structure this
//! scalar world uniquely exposes.
//!
//! A surreal Hahn series carries the ω-adic valuation
//! `v(Σ ω^{y_i} r_i) = y_0` (the leading exponent). In this finite-support model
//! the residue coefficient is rational, while the ideal mathematical comparison
//! point has residue field ℝ and value group `No` under addition. Springer's
//! theorem decomposes a form over such a field into valuation-graded residue
//! forms.
//!
//! ## The honest headline: no bigger Witt group
//!
//! Springer's theorem gives `W(F) ≅ W(k) ⊕ (W(k) ⊗ Γ/2Γ)` for a Henselian valued
//! field with residue field `k` and value group `Γ`. For the full surreal value
//! group, `Γ` is **2-divisible** (every surreal has a half: `a/2` exists), so
//! `Γ/2Γ = 0` and the second summand vanishes. In this crate, treat the
//! filtration as implemented valuation data, not as proof that the finite-support
//! scalar model is itself a full real-closed field.
//!
//! What *is* new, and what this module exposes, is the **valuation filtration**
//! itself: the form's entries grouped by ω-adic valuation, each graded piece a
//! residue ℝ-form read off by sign. No Archimedean Clifford library sees this
//! structure, because over ℝ every nonzero entry has valuation 0. The built-in
//! cross-check is that the residue signatures sum to the ordinary signature
//! (`classify::classify_surreal`).

use crate::clifford::Metric;
use crate::scalar::Scalar;
use crate::scalar::Surreal;
use std::cmp::Ordering;

/// One graded piece of a Springer decomposition: a residue ℝ-form at a fixed
/// ω-adic valuation, recorded by its signature `(p, q)` = (#positive, #negative)
/// residue squares.
#[derive(Debug, Clone, PartialEq)]
pub struct ResidueForm {
    /// The ω-adic valuation (leading exponent) of this graded piece.
    pub valuation: Surreal,
    /// The residue ℝ-form's signature `(p, q)`.
    pub signature: (usize, usize),
}

/// A Springer decomposition: the valuation-graded residue forms (sorted by
/// valuation, most-infinite first), the radical dimension (genuinely zero
/// entries), and the total ordinary signature (which must equal the form's
/// `classify_surreal` signature).
#[derive(Debug, Clone, PartialEq)]
pub struct SpringerDecomp {
    pub graded: Vec<ResidueForm>,
    pub radical_dim: usize,
    pub total_signature: (usize, usize),
}

/// Decompose a surreal quadratic form by ω-adic valuation. Non-orthogonal
/// symmetric metrics are first congruence-diagonalized, matching
/// [`classify_surreal`](crate::forms::classify_surreal). Returns `None` only if
/// that diagonalization is impossible in the finite surreal representation.
pub fn springer_decompose(metric: &Metric<Surreal>) -> Option<SpringerDecomp> {
    let metric = crate::forms::as_diagonal(metric)?;
    // buckets: (valuation, (p, q)) — Surreal is not Ord/Hash, so we bucket with
    // an O(n²) cmp-based scan (fine for the small metrics this is used on).
    let mut buckets: Vec<(Surreal, (usize, usize))> = Vec::new();
    let mut radical_dim = 0usize;

    for x in &metric.q {
        if x.is_zero() {
            radical_dim += 1;
            continue;
        }
        // leading term = (ω-adic valuation, residue coefficient)
        let (exp, coeff) = x.terms().first().expect("nonzero surreal has a term");
        let sign = coeff.sign();
        let slot = buckets
            .iter_mut()
            .find(|(v, _)| v.cmp(exp) == Ordering::Equal);
        let (p, q) = match slot {
            Some((_, pq)) => pq,
            None => {
                buckets.push((exp.clone(), (0, 0)));
                &mut buckets.last_mut().unwrap().1
            }
        };
        match sign {
            Ordering::Greater => *p += 1,
            Ordering::Less => *q += 1,
            Ordering::Equal => unreachable!("nonzero surreal has nonzero leading coeff"),
        }
    }

    // sort by valuation, descending (most infinite first)
    buckets.sort_by(|a, b| b.0.cmp(&a.0));

    let mut total = (0usize, 0usize);
    let graded: Vec<ResidueForm> = buckets
        .into_iter()
        .map(|(valuation, signature)| {
            total.0 += signature.0;
            total.1 += signature.1;
            ResidueForm {
                valuation,
                signature,
            }
        })
        .collect();

    Some(SpringerDecomp {
        graded,
        radical_dim,
        total_signature: total,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forms::classify_surreal;
    use crate::scalar::Rational;

    fn w(n: i128) -> Surreal {
        Surreal::from_int(n)
    }

    #[test]
    fn three_valuation_levels() {
        // [ω, ε, 1, −1]: v=1 → (1,0), v=0 → (1,1), v=−1 → (1,0); total (3,1).
        let m = Metric::diagonal(vec![Surreal::omega(), Surreal::epsilon(), w(1), w(-1)]);
        let d = springer_decompose(&m).unwrap();
        assert_eq!(d.graded.len(), 3);
        // sorted most-infinite first: valuations 1, 0, −1
        assert_eq!(d.graded[0].valuation, w(1));
        assert_eq!(d.graded[0].signature, (1, 0));
        assert_eq!(d.graded[1].valuation, w(0));
        assert_eq!(d.graded[1].signature, (1, 1));
        assert_eq!(d.graded[2].valuation, w(-1));
        assert_eq!(d.graded[2].signature, (1, 0));
        assert_eq!(d.total_signature, (3, 1));
        // the built-in cross-check: residue signatures sum to the ordinary one
        assert_eq!(d.total_signature, classify_surreal(&m).unwrap().signature);
    }

    #[test]
    fn single_valuation_bucket() {
        // [ω, 2ω, −ω]: all valuation 1, residue signs +,+,− ⇒ one bucket (2,1).
        let two_omega = Surreal::monomial(Surreal::one(), Rational::int(2));
        let m = Metric::diagonal(vec![Surreal::omega(), two_omega, Surreal::omega().neg()]);
        let d = springer_decompose(&m).unwrap();
        assert_eq!(d.graded.len(), 1);
        assert_eq!(d.graded[0].valuation, w(1));
        assert_eq!(d.graded[0].signature, (2, 1));
        assert_eq!(d.total_signature, (2, 1));
    }

    #[test]
    fn reads_only_the_leading_term() {
        // [ω+1]: valuation 1, residue + (the dominated +1 is invisible).
        let m = Metric::diagonal(vec![Surreal::omega().add(&w(1))]);
        let d = springer_decompose(&m).unwrap();
        assert_eq!(d.graded.len(), 1);
        assert_eq!(d.graded[0].valuation, w(1));
        assert_eq!(d.graded[0].signature, (1, 0));
    }

    #[test]
    fn radical_is_counted_separately() {
        let m = Metric::diagonal(vec![w(0), Surreal::omega()]);
        let d = springer_decompose(&m).unwrap();
        assert_eq!(d.radical_dim, 1);
        assert_eq!(d.graded.len(), 1);
        assert_eq!(d.total_signature, (1, 0));
    }

    #[test]
    fn witt_class_is_just_the_signature() {
        // The honesty test: because the value group is 2-divisible, the only
        // Witt invariant is the total signature p−q (W(No)=W(ℝ)=ℤ). The graded
        // filtration is richer than ℝ, but it does not yield a bigger group.
        let m = Metric::diagonal(vec![Surreal::omega(), w(1), Surreal::epsilon().neg()]);
        let d = springer_decompose(&m).unwrap();
        let witt = d.total_signature.0 as isize - d.total_signature.1 as isize;
        let (p, q) = classify_surreal(&m).unwrap().signature;
        assert_eq!(witt, p as isize - q as isize);
    }

    #[test]
    fn nonorthogonal_metric_is_diagonalized_first() {
        let mut b = std::collections::BTreeMap::new();
        b.insert((0usize, 1usize), Surreal::from_int(1));
        let m = Metric::new(vec![w(0), w(0)], b);
        let d = springer_decompose(&m).unwrap();
        assert_eq!(d.total_signature, (1, 1));
        assert_eq!(d.total_signature, classify_surreal(&m).unwrap().signature);
    }
}
