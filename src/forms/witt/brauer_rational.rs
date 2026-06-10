//! The **ungraded rational Brauer class** of a quadratic form over `ℚ` — the
//! char-0 / odd mirror of Bridge B (which classified the *char-2* Clifford algebra
//! by its Arf/Brauer–Wall bit), done **correctly**: the Hasse–Witt invariant is
//! *not* the Brauer class of the Clifford algebra, and this module computes the
//! exact correction between them.
//!
//! ## Two distinct invariants
//!
//! Over `ℚ`, quadratic-form Brauer invariants live in `Br(ℚ)[2]`, which by
//! Hasse–Brauer–Noether injects into `⊕_v Br(ℚ_v)[2] = ⊕_v {±1}` — a finite set of
//! ramified places of **even** cardinality (`∏_v = +1`, Hilbert reciprocity, the
//! oracle already in [`local_global`](crate::forms::local_global)). A
//! [`Brauer2Class`] *is* that ramification set. For `q = ⟨a₁,…,aₙ⟩` two **distinct**
//! 2-torsion classes:
//!
//! ```text
//! Hasse–Witt   s(q) = Σ_{i<j} (aᵢ, aⱼ)         (the per-place pieces are the Hilbert
//!                                                products already in `try_hasse_at_place`)
//! Clifford     c(q) = [ C(q) ]   (n even)       (the Brauer class of the Clifford algebra;
//!              c(q) = [ C₀(q) ]  (n odd)         even part in odd rank)
//! ```
//!
//! They are **not equal**. They differ by an explicit `n mod 8` / discriminant
//! correction (Lam, *Introduction to Quadratic Forms over Fields*, GSM 67, pp.
//! 117–119; the same table SageMath's `clifford_invariant` implements). Writing
//! `d = a₁·…·aₙ ∈ ℚ*/ℚ*²` (the **unsigned** discriminant) and additively in
//! `Br(ℚ)[2]`:
//!
//! ```text
//! c(q) = s(q) + δ(n mod 8, d),     δ =  0                  for n ≡ 1, 2
//!                                       (−1,−1) + (−1, d)   for n ≡ 3, 4
//!                                       (−1,−1)             for n ≡ 5, 6
//!                                       (−1, d)             for n ≡ 7, 0
//! ```
//!
//! So [`hasse_brauer_class`] reads `s(q)` off the Hilbert products, and
//! [`clifford_brauer_class`] applies the correction to obtain `c(q)`. The honest
//! bridge verifies the **correction**, not an identity. The independent oracle for
//! small forms is the direct Clifford structure: `C(⟨a,b⟩) ≅ (a,b)` and
//! `C₀(⟨a,b,c⟩) ≅ (−ab, −ac)` (the quaternion factor of the even subalgebra) — the
//! "clifford-side reader" the bridge proposes, exercised in the tests.
//!
//! ## Scope (honest boundaries)
//!
//! `ℚ` (and `ℚ_v`) only; **2-torsion only** (quadratic-form Brauer classes are
//! 2-torsion). This is the **ungraded** Brauer class — kept strictly distinct from
//! the graded [`BrauerWallClass`](crate::forms::bw_class_real); conflating them is
//! exactly what `char0.rs` declines to do, and the rational ungraded class is what
//! Bridge F adds. The full `ℚ/ℤ` lift via cyclic algebras (Bridge K) would embed
//! this as its 2-torsion `½`-slice (one shared class type, two constructors); that
//! generalization is tracked in `roadmap/TBD.md`, not built here.

use std::collections::BTreeSet;

use crate::forms::{
    relevant_primes, try_disc_class, try_hasse_at_place, try_hilbert_symbol_at, Place,
};
use crate::scalar::{Rational, Scalar};

/// The ungraded rational **2-torsion Brauer class**, represented by its set of
/// **ramified places** (the places `v` where `inv_v = ½`). The split class is the
/// empty set; the group law is symmetric difference (XOR of indicator sets =
/// addition of `½`'s mod `1`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Brauer2Class {
    ramified: BTreeSet<Place>,
}

impl Brauer2Class {
    /// The split (trivial) class: ramified nowhere.
    pub fn split() -> Self {
        Brauer2Class {
            ramified: BTreeSet::new(),
        }
    }

    /// Whether this is the split class.
    pub fn is_split(&self) -> bool {
        self.ramified.is_empty()
    }

    /// The set of ramified places (where `inv_v = ½`).
    pub fn ramified_places(&self) -> &BTreeSet<Place> {
        &self.ramified
    }

    /// The Brauer-group sum (tensor product of algebras): symmetric difference of
    /// the ramification sets — the 2-torsion addition law (XOR).
    pub fn add(&self, other: &Self) -> Self {
        Brauer2Class {
            ramified: self
                .ramified
                .symmetric_difference(&other.ramified)
                .copied()
                .collect(),
        }
    }

    /// The local invariant `inv_v ∈ {0, ½} ⊂ ℚ/ℤ` at a place.
    pub fn local_invariant(&self, place: Place) -> Rational {
        if self.ramified.contains(&place) {
            Rational::try_new(1, 2).expect("1/2 is a valid rational")
        } else {
            Rational::zero()
        }
    }

    /// **Hilbert reciprocity**, additively: a global Brauer class ramifies at an
    /// **even** number of places (`∑_v inv_v ≡ 0 mod ℤ`). True for every class this
    /// module builds.
    pub fn satisfies_reciprocity(&self) -> bool {
        self.ramified.len().is_multiple_of(2)
    }

    /// The class of the **quaternion algebra** `(a, b)` over `ℚ`: ramified exactly
    /// at the places `v` where the Hilbert symbol `(a, b)_v = −1`. `None` if a
    /// Hilbert symbol is undefined (an argument is `0` or square-class arithmetic
    /// overflows `i128`).
    pub fn quaternion(a: i128, b: i128) -> Option<Self> {
        if a == 0 || b == 0 {
            return None;
        }
        let mut ramified = BTreeSet::new();
        if try_hilbert_symbol_at(a, b, Place::Real)? == -1 {
            ramified.insert(Place::Real);
        }
        for p in relevant_primes(&[a, b]) {
            if try_hilbert_symbol_at(a, b, Place::Prime(p))? == -1 {
                ramified.insert(Place::Prime(p));
            }
        }
        Some(Brauer2Class { ramified })
    }
}

/// The **Hasse–Witt invariant** `s(q) = Σ_{i<j} (aᵢ, aⱼ)` as a Brauer class: the
/// places where the per-place Hasse invariant `∏_{i<j}(aᵢ,aⱼ)_v` is `−1`. The
/// entries are nonzero integer square-class representatives of a nondegenerate
/// diagonal form (a rational form scales to integer entries without changing the
/// class). `None` on a zero entry (the radical is not part of the nondegenerate
/// Brauer class) or on bounded-arithmetic overflow.
pub fn hasse_brauer_class(entries: &[i128]) -> Option<Brauer2Class> {
    if entries.contains(&0) {
        return None;
    }
    let mut ramified = BTreeSet::new();
    if try_hasse_at_place(entries, Place::Real)? == -1 {
        ramified.insert(Place::Real);
    }
    for p in relevant_primes(entries) {
        if try_hasse_at_place(entries, Place::Prime(p))? == -1 {
            ramified.insert(Place::Prime(p));
        }
    }
    Some(Brauer2Class { ramified })
}

/// The `n mod 8` / discriminant correction `δ` between the Hasse–Witt and Clifford
/// invariants (Lam, GSM 67, pp. 117–119): `(−1, d)` for `n ≡ 0,3,4,7`, plus
/// `(−1,−1)` for `n ≡ 3,4,5,6`.
fn clifford_correction(n: usize, d: i128) -> Option<Brauer2Class> {
    let r = n % 8;
    let mut delta = Brauer2Class::split();
    if matches!(r, 0 | 3 | 4 | 7) {
        delta = delta.add(&Brauer2Class::quaternion(-1, d)?);
    }
    if matches!(r, 3..=6) {
        delta = delta.add(&Brauer2Class::quaternion(-1, -1)?);
    }
    Some(delta)
}

/// The **Clifford invariant** `c(q) = [C(q)]` (`n` even) / `[C₀(q)]` (`n` odd) as a
/// rational Brauer class: the Hasse–Witt class corrected by the `n mod 8` /
/// discriminant term `δ`. The honest char-0 analogue of Bridge B — the algebra the
/// `clifford` pillar builds, classified by the symbols the `forms` pillar computes.
///
/// Entries as in [`hasse_brauer_class`]. `None` on a zero entry or overflow.
pub fn clifford_brauer_class(entries: &[i128]) -> Option<Brauer2Class> {
    if entries.contains(&0) {
        return None;
    }
    let s = hasse_brauer_class(entries)?;
    let d = if entries.is_empty() {
        1
    } else {
        try_disc_class(entries)?
    };
    Some(s.add(&clifford_correction(entries.len(), d)?))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn places(ps: &[Place]) -> BTreeSet<Place> {
        ps.iter().copied().collect()
    }

    fn clifford(entries: &[i128]) -> Brauer2Class {
        clifford_brauer_class(entries).expect("test square classes fit i128")
    }

    fn hasse(entries: &[i128]) -> Brauer2Class {
        hasse_brauer_class(entries).expect("test square classes fit i128")
    }

    fn quat(a: i128, b: i128) -> Brauer2Class {
        Brauer2Class::quaternion(a, b).expect("test quaternion is defined")
    }

    // --- the XOR group law ---

    #[test]
    fn add_is_xor_of_ramification_sets() {
        let h = quat(-1, -1); // {ℝ, Q_2}
        assert!(h.add(&h).is_split(), "x + x = 0 (2-torsion)");
        assert_eq!(h.add(&Brauer2Class::split()), h, "0 is the identity");
        // commutativity is free from BTreeSet symmetric difference; spot-check a sum.
        let k = quat(2, 5);
        assert_eq!(h.add(&k), k.add(&h));
    }

    // --- anchors: known algebras ---

    #[test]
    fn split_form_is_split() {
        // ⟨1,−1⟩: C(q) = M₂(ℚ), split. n=2 ⇒ no correction, and (1,−1) is split.
        assert!(clifford(&[1, -1]).is_split());
        assert!(hasse(&[1, -1]).is_split());
    }

    #[test]
    fn hamilton_quaternions_ramify_at_2_and_infinity() {
        // ⟨−1,−1,−1⟩: C₀(q) ≅ Hamilton (−1,−1), ramified {ℝ, Q_2}.
        assert_eq!(
            *clifford(&[-1, -1, -1]).ramified_places(),
            places(&[Place::Real, Place::Prime(2)])
        );
        // ⟨1,1,1⟩: also Hamilton — but its *Hasse* class is split (s = 0); only the
        // n≡3 correction (−1,−1)+(−1,1) = (−1,−1) supplies the Clifford class. The
        // sharpest demonstration that c ≠ s.
        assert!(hasse(&[1, 1, 1]).is_split());
        assert_eq!(
            *clifford(&[1, 1, 1]).ramified_places(),
            places(&[Place::Real, Place::Prime(2)])
        );
    }

    // --- the independent clifford-side oracle (the quaternion factor) ---

    #[test]
    fn rank_two_clifford_is_the_quaternion_algebra() {
        // C(⟨a,b⟩) ≅ (a,b): the n=2 Clifford invariant is literally the quaternion
        // class — an oracle that never goes through the s+correction route.
        for a in [-3i128, -2, -1, 1, 2, 3, 5, 6, 7] {
            for b in [-5i128, -3, -1, 1, 2, 3, 5, 7, 10] {
                assert_eq!(clifford(&[a, b]), quat(a, b), "C(⟨{a},{b}⟩) ≠ (a,b)");
            }
        }
    }

    #[test]
    fn rank_three_even_clifford_is_minus_ab_minus_ac() {
        // C₀(⟨a,b,c⟩) ≅ (−ab, −ac): the independent ternary oracle (the even
        // subalgebra's quaternion factor), validating the n≡3 correction.
        for a in [-3i128, -1, 1, 2, 3, 5] {
            for b in [-3i128, -1, 1, 2, 5, 7] {
                for c in [-5i128, -1, 1, 3, 6] {
                    assert_eq!(
                        clifford(&[a, b, c]),
                        quat(-a * b, -a * c),
                        "C₀(⟨{a},{b},{c}⟩) ≠ (−ab,−ac)"
                    );
                }
            }
        }
    }

    #[test]
    fn rank_one_is_always_split() {
        // C₀(⟨a⟩) = ℚ, c = 0, for every a.
        for a in [-7i128, -2, -1, 1, 2, 3, 5] {
            assert!(clifford(&[a]).is_split(), "⟨{a}⟩ should be split");
        }
    }

    // --- reciprocity: even ramification across a sweep of forms ---

    #[test]
    fn every_class_satisfies_reciprocity() {
        let forms: &[&[i128]] = &[
            &[1, -1],
            &[2, 3],
            &[-1, -1, -1],
            &[1, 1, 1],
            &[2, 3, 5],
            &[1, -2, -5],
            &[1, 1, 1, 1],
            &[1, 1, 1, -1],
            &[2, 3, 5, 7],
            &[1, 1, 1, 1, 1],
            &[-1, -2, -3, -5, -7],
            &[2, 3, 5, 7, 11, 13],
        ];
        for f in forms {
            assert!(
                clifford(f).satisfies_reciprocity(),
                "c({f:?}) ramifies oddly"
            );
            assert!(hasse(f).satisfies_reciprocity(), "s({f:?}) ramifies oddly");
        }
    }

    // --- the correction table itself: c(q) vs s(q) per n mod 8 ---

    /// The independently-tabulated correction (Lam GSM 67 pp. 117–119), recomputed
    /// in the test from `Brauer2Class::quaternion` — so a mis-encoded match arm in
    /// `clifford_correction` is caught.
    fn expected_correction(n: usize, d: i128) -> Brauer2Class {
        let mut delta = Brauer2Class::split();
        if matches!(n % 8, 3..=6) {
            delta = delta.add(&quat(-1, -1));
        }
        if matches!(n % 8, 0 | 3 | 4 | 7) {
            delta = delta.add(&quat(-1, d));
        }
        delta
    }

    #[test]
    fn clifford_is_hasse_plus_the_documented_correction() {
        let forms: &[&[i128]] = &[
            &[1],                        // n=1
            &[2, 3],                     // n=2
            &[1, 2, 3],                  // n=3
            &[1, 2, 3, 5],               // n=4
            &[1, 2, 3, 5, 7],            // n=5
            &[1, 2, 3, 5, 7, 11],        // n=6
            &[2, 3, 5, 7, 11, 13, 1],    // n=7
            &[1, 2, 3, 5, 7, 11, 13, 1], // n=8 ≡ 0
        ];
        for f in forms {
            let d = try_disc_class(f).expect("disc fits i128");
            let expected = hasse(f).add(&expected_correction(f.len(), d));
            assert_eq!(clifford(f), expected, "correction mismatch for {f:?}");
        }
    }

    #[test]
    fn n_one_and_two_have_no_correction() {
        // n ≡ 1, 2: c(q) = s(q) exactly (δ = 0).
        for f in [&[3i128] as &[i128], &[5], &[2, 7], &[-3, 5]] {
            assert_eq!(clifford(f), hasse(f), "δ should vanish for {f:?}");
        }
    }

    #[test]
    fn rejects_degenerate_and_overflow() {
        assert_eq!(clifford_brauer_class(&[1, 0, 1]), None); // radical present
        assert_eq!(hasse_brauer_class(&[0]), None);
    }
}
