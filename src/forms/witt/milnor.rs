//! Bridge N.1 — Milnor's exact sequence: the Springer residues assembled globally.
//!
//! The shipped Springer engine (`springer/`) computes per-place residue buckets and
//! the local–global layer decides per-form isotropy; this module assembles the
//! Witt-**group**-level global statement. Milnor's exact sequence supplies it
//! (Milnor–Husemoller, *Symmetric Bilinear Forms*, Ch. IV; Lam, GSM 67, Ch. IX):
//!
//! ```text
//! 0 → W(ℤ) → W(ℚ) →∂ ⊕_p W(F_p) → 0        (exact)
//! ```
//!
//! The kernel `W(ℤ) ≅ ℤ` is detected by the **signature**; for odd `p`, the boundary
//! `∂_p` is the **second Springer residue** lifted from `LocalResidueForm` buckets to
//! Witt classes. For `p = 2`, Milnor's hand-defined boundary lands in
//! `W(F₂) ≅ ℤ/2`: a diagonal line contributes exactly when its `2`-adic valuation is
//! odd (the residue unit is then the unique nonzero element of `F₂`). So
//! `(signature, (∂_p)_p)` is a *complete* invariant of `W(ℚ)`: two rational diagonal
//! forms are Witt-equivalent over `ℚ` iff they share a signature and all residues —
//! the sequence ties three pillar surfaces together (the Springer residues, the
//! global field layer, and the integral pillar's signature).
//!
//! The equal-characteristic odd leg uses the split form of the same idea:
//!
//! ```text
//! W(F_q(t)) ≅ W(F_q) ⊕ ⊕_π W(F_q[t]/π).
//! ```
//!
//! [`global_residues_ff`] returns the `W(F_q)` summand from the even-valuation layer
//! at the degree place `∞`, plus the nonzero second residues at finite monic
//! irreducible places. This is exact on the shipped `RationalFunction`/`Poly`
//! backend and uses the same `FFPlace` arithmetic as the function-field Hilbert and
//! Hasse–Minkowski layers.
//!
//! **Claim level:** standard math (Milnor; Lam GSM 67, Ch. IX) made computational.
//! The residue is computed directly from the `i128` entries (`v_p`, the Legendre
//! symbol, and the signed-discriminant square class), matching the
//! [`finite_odd_witt`](crate::forms::finite_odd_witt) convention, so it is **exact**;
//! `springer_decompose_qp` on the capped `Q_p` model is the cross-check oracle.
//!
//! **The `∂₂` boundary (load-bearing).** `∂₂` (residue characteristic 2) is **not**
//! Springer's second residue — Milnor defines it by hand in Ch. IV. This module uses
//! the crate's existing char-2 [`WittClassG`] carrier as the `W(F₂) ≅ ℤ/2` target:
//! `Char2 { field_degree: 1, arf }`, with `arf` the parity of odd dyadic valuation
//! lines. The char-2 constant fields of `F_q(t)` are a separate matter (the
//! Aravire–Jacob layer in `springer/char2.rs`), and tame/wild norm-residue symbols
//! stay with the cyclic-Brauer follow-ons rather than this Witt-residue map.

use crate::forms::local_global::padic::{legendre, relevant_primes, unit_part, val_p};
use crate::forms::{
    try_chi_kappa, try_kappa_order, try_relevant_places_ff, try_residue_unit_at,
    try_valuation_at_ff, FFPlace, FiniteOddField, WittClassG,
};
use crate::scalar::{Poly, RationalFunction, Scalar};
use std::collections::BTreeMap;

/// The split Milnor invariant of a diagonal form over odd `F_q(t)`.
///
/// The first component is the constant-field class selected at `∞`; the vector is
/// the finite-place support of nonzero second residues.
pub type FunctionFieldMilnorResidues<S> = (WittClassG, Vec<(FFPlace<S>, WittClassG)>);

/// The second residue `∂_p⟨a_1,…,a_n⟩` at an **odd** prime `p`, as a Witt class over
/// `F_p`. It collects the residue units of the entries of **odd** `p`-valuation and
/// returns the Witt class of `⟂ ⟨ū_i⟩` over `F_p`, using the multiplicativity of the
/// Legendre symbol (so no product overflows): `∏ (u_i | p)` times the
/// `(−1)^{m(m−1)/2}` signed-discriminant correction gives the square class.
fn second_residue_at(entries: &[i128], p: u128) -> WittClassG {
    let pi = p as i128;
    let mut leg_prod: i128 = 1; // ∏ (u_i | p) over odd-valuation entries
    let mut m: i128 = 0; // dimension of the residue form
    for &a in entries {
        if val_p(a, pi) % 2 == 1 {
            leg_prod *= legendre(unit_part(a, pi), pi);
            m += 1;
        }
    }
    let leg_neg1 = legendre(-1, pi); // (−1 | p): +1 iff p ≡ 1 (mod 4)
    let signed_leg = if ((m * (m - 1) / 2) & 1) == 1 {
        leg_prod * leg_neg1
    } else {
        leg_prod
    };
    WittClassG::OddChar {
        field_order: p,
        kappa: if leg_neg1 == 1 { 0 } else { 1 },
        e0: (m & 1) as u128,
        sclass: if signed_leg == 1 { 0 } else { 1 },
    }
}

/// Milnor's hand-defined dyadic residue `∂₂ : W(ℚ) → W(F₂) ≅ ℤ/2`.
/// Since every odd unit reduces to `1 ∈ F₂`, only the parity of entries with odd
/// `2`-adic valuation survives.
fn dyadic_residue_at(entries: &[i128]) -> WittClassG {
    let arf = entries.iter().filter(|&&a| val_p(a, 2) % 2 == 1).count() as u128 & 1;
    WittClassG::Char2 {
        field_degree: 1,
        arf,
    }
}

/// Whether a Witt class over `F_p` is the zero class (even dimension and square signed
/// discriminant ⇒ hyperbolic).
fn is_zero_residue(w: &WittClassG) -> bool {
    matches!(
        w,
        WittClassG::OddChar {
            e0: 0,
            sclass: 0,
            ..
        } | WittClassG::Char2 { arf: 0, .. }
    )
}

/// The image of the rational diagonal form `⟨a_1,…,a_n⟩` (nonzero `i128` entries)
/// under the Milnor map `W(ℚ) → ℤ ⊕ ⊕_p W(F_p)`: the **signature** `(#positive −
/// #negative)` and the nonzero residues `∂_p`, keyed by prime. Zero residues are
/// omitted, so the map of an everywhere-good integral form is empty.
///
/// `None` if any entry is zero (a radical — the form is degenerate). Two forms with
/// equal `global_residues` are Witt-equivalent over `ℚ`; a difference at any prime,
/// or in the signature, witnesses inequivalence.
pub fn global_residues(entries: &[i128]) -> Option<(i128, BTreeMap<u128, WittClassG>)> {
    if entries.contains(&0) {
        return None;
    }
    let signature: i128 = entries.iter().map(|&a| a.signum()).sum();
    let mut residues = BTreeMap::new();
    for p in relevant_primes(entries) {
        let w = if p == 2 {
            dyadic_residue_at(entries)
        } else {
            second_residue_at(entries, p)
        };
        if !is_zero_residue(&w) {
            residues.insert(p, w);
        }
    }
    Some((signature, residues))
}

fn oddchar_witt_from_residue_units<S: FiniteOddField>(
    units: &[Poly<S>],
    place: &FFPlace<S>,
) -> Option<WittClassG> {
    let mut chi_prod: i128 = 1;
    for unit in units {
        chi_prod *= try_chi_kappa(unit, place)?;
    }
    let m = i128::try_from(units.len()).ok()?;
    let field_order = try_kappa_order(place)?;
    let chi_neg1 = if field_order % 4 == 1 { 1 } else { -1 };
    let signed_chi = if ((m * (m - 1) / 2) & 1) == 1 {
        chi_prod * chi_neg1
    } else {
        chi_prod
    };
    Some(WittClassG::OddChar {
        field_order,
        kappa: if chi_neg1 == 1 { 0 } else { 1 },
        e0: (m & 1) as u128,
        sclass: if signed_chi == 1 { 0 } else { 1 },
    })
}

fn second_residue_at_ff<S: FiniteOddField>(
    entries: &[RationalFunction<S>],
    place: &FFPlace<S>,
) -> Option<WittClassG> {
    let mut units = Vec::new();
    for entry in entries {
        if try_valuation_at_ff(entry, place)?.rem_euclid(2) != 0 {
            units.push(try_residue_unit_at(entry, place)?);
        }
    }
    oddchar_witt_from_residue_units(&units, place)
}

fn constant_class_at_infinity_ff<S: FiniteOddField>(
    entries: &[RationalFunction<S>],
) -> Option<WittClassG> {
    let place = FFPlace::Infinite;
    let mut units = Vec::new();
    for entry in entries {
        if try_valuation_at_ff(entry, &place)?.rem_euclid(2) == 0 {
            units.push(try_residue_unit_at(entry, &place)?);
        }
    }
    oddchar_witt_from_residue_units(&units, &place)
}

/// The split Milnor map for a diagonal form over `F_q(t)` with odd `q`:
/// `W(F_q(t)) ≅ W(F_q) ⊕ ⊕_π W(F_q[t]/π)`.
///
/// The first component is the `W(F_q)` class obtained by the even-valuation
/// layer at the degree place `∞`; the vector contains the nonzero second
/// residues at finite monic irreducible places. Zero residues are omitted.
///
/// `None` if any entry is zero. Characteristic-2 function fields use the
/// separate Artin-Schreier/Aravire-Jacob layer, not this tame odd-residue
/// sequence.
pub fn global_residues_ff<S: FiniteOddField>(
    entries: &[RationalFunction<S>],
) -> Option<FunctionFieldMilnorResidues<S>> {
    if entries.iter().any(|entry| entry.is_zero()) {
        return None;
    }
    let constant = constant_class_at_infinity_ff(entries)?;
    let mut residues = Vec::new();
    for place in try_relevant_places_ff(entries)? {
        if matches!(place, FFPlace::Infinite) {
            continue;
        }
        let w = second_residue_at_ff(entries, &place)?;
        if !is_zero_residue(&w) {
            residues.push((place, w));
        }
    }
    Some((constant, residues))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::clifford::Metric;
    use crate::forms::{springer_decompose_qp, try_is_isotropic_q};
    use crate::scalar::{Fp, Qp, RationalFunction};

    /// `∂₅` via the capped `Q₅` Springer engine: the Witt class of the odd-valuation
    /// (parity-1) residue layer, built independently of the `i128` route.
    fn springer_residue_q5(entries: &[i128]) -> WittClassG {
        type Q5 = Qp<5, 6>;
        let metric = Metric::diagonal(entries.iter().map(|&a| Q5::from_i128(a)).collect());
        let decomp = springer_decompose_qp(&metric).unwrap();
        let mut dim = 0usize;
        let mut disc_sq = true; // running square class of the residue discriminant
        for form in decomp.parity_layer(1) {
            dim += form.dim;
            disc_sq = disc_sq == form.disc_is_square; // XNOR of square classes
        }
        let m = dim as i128;
        let leg_neg1 = legendre(-1, 5); // +1 (5 ≡ 1 mod 4)
        let signed_sq = if ((m * (m - 1) / 2) & 1) == 1 && leg_neg1 != 1 {
            !disc_sq
        } else {
            disc_sq
        };
        WittClassG::OddChar {
            field_order: 5,
            kappa: if leg_neg1 == 1 { 0 } else { 1 },
            e0: (dim & 1) as u128,
            sclass: if signed_sq { 0 } else { 1 },
        }
    }

    fn f2_class(arf: u128) -> WittClassG {
        WittClassG::Char2 {
            field_degree: 1,
            arf,
        }
    }

    type F5 = RationalFunction<Fp<5>>;
    type Poly5 = Poly<Fp<5>>;

    fn rf(num: &[i128], den: &[i128]) -> F5 {
        RationalFunction::new(
            num.iter().map(|&n| Fp::<5>::new(n)).collect(),
            den.iter().map(|&n| Fp::<5>::new(n)).collect(),
        )
    }

    fn poly(c: &[i128]) -> Poly5 {
        Poly::new(c.iter().map(|&n| Fp::<5>::new(n)).collect())
    }

    fn odd_class(field_order: u128, e0: u128, sclass: u128) -> WittClassG {
        WittClassG::OddChar {
            field_order,
            kappa: if field_order % 4 == 1 { 0 } else { 1 },
            e0,
            sclass,
        }
    }

    fn residue_at<'a>(
        residues: &'a [(FFPlace<Fp<5>>, WittClassG)],
        place: &FFPlace<Fp<5>>,
    ) -> Option<&'a WittClassG> {
        residues.iter().find(|(pl, _)| pl == place).map(|(_, w)| w)
    }

    #[test]
    fn second_residue_matches_springer_over_q5() {
        // The exact i128 residue and the capped-Q₅ Springer residue agree on forms
        // exercising even/odd valuations and square/nonsquare units at 5.
        for entries in [
            vec![1, 5],
            vec![2, 10],
            vec![3, 15, 5],
            vec![1, 1],
            vec![7, 5, 25, 2],
        ] {
            assert_eq!(
                second_residue_at(&entries, 5),
                springer_residue_q5(&entries),
                "∂₅ mismatch on {entries:?}"
            );
        }
    }

    #[test]
    fn dyadic_residue_is_milnors_hand_boundary() {
        // Over F_2 every odd unit reduces to 1, so ∂_2 only sees the parity of
        // odd 2-adic valuation lines.
        assert_eq!(dyadic_residue_at(&[1]), f2_class(0));
        assert_eq!(dyadic_residue_at(&[2]), f2_class(1));
        assert_eq!(dyadic_residue_at(&[-2]), f2_class(1));
        assert_eq!(dyadic_residue_at(&[1, 2]), f2_class(1));
        assert_eq!(dyadic_residue_at(&[2, -2]), f2_class(0));
    }

    #[test]
    fn global_residues_include_the_dyadic_cell() {
        for (entries, signature) in [(&[2i128][..], 1), (&[1, 2], 2), (&[-2], -1)] {
            let (sig, res) = global_residues(entries).unwrap();
            assert_eq!(sig, signature);
            assert_eq!(res.get(&2), Some(&f2_class(1)), "entries={entries:?}");
        }

        let (sig, res) = global_residues(&[2, -2]).unwrap();
        assert_eq!(sig, 0);
        assert!(
            res.is_empty(),
            "the hyperbolic pair <2,-2> has zero residues"
        );

        let (_, mixed) = global_residues(&[6]).unwrap();
        assert_eq!(
            mixed.keys().copied().collect::<Vec<_>>(),
            vec![2, 3],
            "<6> has both dyadic and odd-prime residues"
        );
    }

    #[test]
    fn residues_have_finite_support_at_dividing_primes() {
        // ∂_p = 0 for p ∤ ∏ a_i: ⟨1,1,1⟩ has no odd residues.
        let (sig, res) = global_residues(&[1, 1, 1]).unwrap();
        assert_eq!(sig, 3);
        assert!(res.is_empty());
        // ⟨3, 5⟩: residues exactly at 3 and 5 (each an odd-valuation unit line).
        let (sig, res) = global_residues(&[3, 5]).unwrap();
        assert_eq!(sig, 2);
        assert_eq!(res.keys().copied().collect::<Vec<_>>(), vec![3, 5]);
    }

    #[test]
    fn radical_entry_is_rejected() {
        assert_eq!(global_residues(&[1, 0, 2]), None);
    }

    #[test]
    fn function_field_residues_split_at_infinity() {
        let (constant, residues) = global_residues_ff(&[rf(&[1], &[1])]).unwrap();
        assert_eq!(constant, odd_class(5, 1, 0));
        assert!(
            residues.is_empty(),
            "constant forms have no finite residues"
        );

        let (constant, residues) = global_residues_ff(&[rf(&[0, 1], &[1])]).unwrap();
        assert_eq!(constant, odd_class(5, 0, 0));
        assert_eq!(
            residue_at(&residues, &FFPlace::Finite(poly(&[0, 1]))),
            Some(&odd_class(5, 1, 0))
        );

        let (constant, residues) = global_residues_ff(&[rf(&[1], &[0, 1])]).unwrap();
        assert_eq!(constant, odd_class(5, 0, 0));
        assert_eq!(
            residue_at(&residues, &FFPlace::Finite(poly(&[0, 1]))),
            Some(&odd_class(5, 1, 0))
        );

        let (constant, residues) = global_residues_ff(&[rf(&[2], &[1])]).unwrap();
        assert_eq!(constant, odd_class(5, 1, 1), "2 is nonsquare in F_5");
        assert!(residues.is_empty());
    }

    #[test]
    fn function_field_residues_see_degree_two_places() {
        let place = FFPlace::Finite(poly(&[2, 0, 1])); // t^2 + 2 irreducible over F_5
        let (constant, residues) = global_residues_ff(&[rf(&[2, 0, 1], &[1])]).unwrap();
        assert_eq!(constant, odd_class(5, 1, 0));
        assert_eq!(residue_at(&residues, &place), Some(&odd_class(25, 1, 0)));
    }

    #[test]
    fn function_field_residues_are_square_and_hyperbolic_stable() {
        let base = global_residues_ff(&[rf(&[0, 1], &[1])]).unwrap();
        let square = rf(&[1, 1], &[1]).mul(&rf(&[1, 1], &[1]));
        let square_multiple = global_residues_ff(&[rf(&[0, 1], &[1]).mul(&square)]).unwrap();
        assert_eq!(square_multiple, base);

        let hyperbolic = global_residues_ff(&[rf(&[0, 1], &[1]), rf(&[0, 4], &[1])]).unwrap();
        assert_eq!(hyperbolic.0, odd_class(5, 0, 0));
        assert!(hyperbolic.1.is_empty());
    }

    #[test]
    fn function_field_residues_reject_radical_entries() {
        assert_eq!(global_residues_ff(&[rf(&[1], &[1]), rf(&[0], &[1])]), None);
    }

    #[test]
    fn witt_invariants_are_square_and_hyperbolic_stable() {
        // ⟨3⟩ ≅ ⟨12⟩ (12 = 3·4, a square multiple) and adding a hyperbolic plane
        // ⟨1,−1⟩ changes nothing — all three share signature and residues.
        let base = global_residues(&[3]).unwrap();
        assert_eq!(global_residues(&[12]).unwrap(), base);
        assert_eq!(global_residues(&[3, 1, -1]).unwrap(), base);
        // Same at the dyadic prime: ⟨2⟩ ≅ ⟨8⟩, and ⟨1,-1⟩ is still hyperbolic.
        let dyadic = global_residues(&[2]).unwrap();
        assert_eq!(global_residues(&[8]).unwrap(), dyadic);
        assert_eq!(global_residues(&[2, 1, -1]).unwrap(), dyadic);
    }

    #[test]
    fn residues_distinguish_inequivalent_forms() {
        // ⟨1⟩ and ⟨3⟩ have equal signature but differ at p = 3 ⇒ not Witt-equivalent.
        let one = global_residues(&[1]).unwrap();
        let three = global_residues(&[3]).unwrap();
        assert_eq!(one.0, three.0, "same signature");
        assert_ne!(one.1, three.1, "different residue at 3");
        // Cross-check with Hasse–Minkowski: ⟨1,−3⟩ is anisotropic over ℚ (3 is not a
        // square), so ⟨1⟩ ⊥ ⟨−3⟩ is not hyperbolic — they are genuinely inequivalent.
        assert_eq!(try_is_isotropic_q(&[1, -3]), Some(false));

        // Same signature, dyadic residue differs: ⟨1⟩ and ⟨2⟩ are not equivalent.
        let two = global_residues(&[2]).unwrap();
        assert_eq!(one.0, two.0, "same signature");
        assert_ne!(one.1, two.1, "different dyadic residue");
        assert_eq!(try_is_isotropic_q(&[1, -2]), Some(false));
    }

    #[test]
    fn reconstruction_agrees_with_hasse_minkowski() {
        // Equal residues + equal signature ⇒ Witt-equivalent ⇒ a ⊥ (−b) hyperbolic,
        // hence isotropic. ⟨3⟩ vs ⟨12⟩: ⟨3,−12⟩ is isotropic (x = 2y).
        assert_eq!(
            global_residues(&[3]).unwrap(),
            global_residues(&[12]).unwrap()
        );
        assert_eq!(try_is_isotropic_q(&[3, -12]), Some(true));

        // ⟨3,5⟩ vs ⟨12,45⟩ (entrywise square multiples): same residues at 3 and 5,
        // and ⟨3,5,−12,−45⟩ is isotropic ((x,z) = (2,1): 3·4 − 12 = 0).
        assert_eq!(
            global_residues(&[3, 5]).unwrap(),
            global_residues(&[12, 45]).unwrap()
        );
        assert_eq!(try_is_isotropic_q(&[3, 5, -12, -45]), Some(true));

        // Dyadic reconstruction: ⟨2⟩ vs ⟨8⟩ differ by a square multiple, so the
        // difference form is isotropic; ⟨2⟩ vs ⟨1⟩ has a dyadic-residue mismatch.
        assert_eq!(
            global_residues(&[2]).unwrap(),
            global_residues(&[8]).unwrap()
        );
        assert_eq!(try_is_isotropic_q(&[2, -8]), Some(true));
        assert_ne!(
            global_residues(&[2]).unwrap(),
            global_residues(&[1]).unwrap()
        );
        assert_eq!(try_is_isotropic_q(&[2, -1]), Some(false));
    }
}
