//! Discriminant quadratic forms of integral lattices and Milgram's Gauss sums.
//!
//! For a nonsingular even lattice `L` with Gram matrix `G`, this module uses the
//! standard presentation
//!
//! ```text
//! A_L = L#/L ~= Z^n / G Z^n,    y |-> G^{-1} y
//! q_L(y) = y^T G^{-1} y mod 2Z.
//! ```
//!
//! The normalized Gauss sum of `q_L` has phase `exp(2*pi*i*signature/8)`.
//! Odd lattices use the parallel [`OddDiscriminantForm`] surface, where
//! `q_L(y) = y^T G^{-1} y mod Z`; their signature congruence needs the
//! Conway-Sloane oddity correction recorded by [`odd_milgram_report`].
//!
//! # Module layout
//!
//! - `complex` — the hand-rolled [`Complex64`] (dependency-free; deliberately
//!   shadows `num_complex::Complex64`).
//! - `gauss_sum` — [`GaussSum`] and the matrix helpers for the Weil representation.
//! - `form` — [`DiscriminantForm`], [`OddDiscriminantForm`],
//!   [`genus_signature_mod8`], [`verify_milgram`], [`odd_milgram_report`].
//! - `phases` — [`FqmPrimaryPhase`], [`FqmGaussPhase`]: the p-primary
//!   Milgram/Brown phase projection of a finite quadratic module.

mod complex;
mod form;
mod gauss_sum;
mod phases;

pub use complex::Complex64;
pub use form::{
    genus_signature_mod8, odd_milgram_report, verify_milgram, verify_odd_milgram, DiscriminantForm,
    OddDiscriminantForm, OddMilgramReport,
};
pub(crate) use form::{phase_mod8_from_q_values, IsoTables};
pub use gauss_sum::GaussSum;
pub use phases::{FqmGaussPhase, FqmPrimaryPhase};

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forms::{
        a_n, are_in_same_genus, d16_plus, d_n, e_6, e_7, e_8, repetition_code,
        type_i_z2_plus_e8_code, IntegralForm,
    };
    use crate::scalar::{Rational, Scalar};

    /// Nikulin's right-hand side: equal signature pairs and isomorphic discriminant
    /// quadratic forms. Both lattices must be even (the `from_lattice` boundary).
    fn nikulin_rhs(a: &IntegralForm, b: &IntegralForm) -> bool {
        if a.signature() != b.signature() {
            return false;
        }
        let qa = DiscriminantForm::from_lattice(a).expect("even lattice a");
        let qb = DiscriminantForm::from_lattice(b).expect("even lattice b");
        qa.is_isomorphic(&qb) == Some(true)
    }

    #[test]
    fn discriminant_iso_is_reflexive_and_q_sensitive() {
        for l in [a_n(1), a_n(3), d_n(4), e_6(), e_7(), e_8()] {
            let q = DiscriminantForm::from_lattice(&l).unwrap();
            assert_eq!(q.is_isomorphic(&q), Some(true), "reflexive");
        }
        // A_1 and E_7 share the group ℤ/2 but have q-values 1/2 vs 3/2 — *not*
        // isomorphic forms. The search must see q, not just the group.
        let a1 = DiscriminantForm::from_lattice(&a_n(1)).unwrap();
        let e7 = DiscriminantForm::from_lattice(&e_7()).unwrap();
        assert_eq!(a1.group, e7.group, "same invariant factors ℤ/2");
        assert_eq!(a1.is_isomorphic(&e7), Some(false), "q distinguishes them");
        // Different groups: ℤ/3 (A_2) vs (ℤ/2)² (A_1 ⊕ A_1).
        let a2 = DiscriminantForm::from_lattice(&a_n(2)).unwrap();
        let a1a1 = DiscriminantForm::from_lattice(&a_n(1).direct_sum(&a_n(1))).unwrap();
        assert_eq!(a2.is_isomorphic(&a1a1), Some(false));
    }

    #[test]
    fn nikulin_genus_iff_signature_and_discriminant_form() {
        // The Milnor pair: even unimodular rank 16, same genus, non-isometric, both
        // with trivial discriminant form — Nikulin says same genus, and it is.
        let e8e8 = e_8().direct_sum(&e_8());
        let d16 = d16_plus();
        assert!(nikulin_rhs(&e8e8, &d16));
        assert!(are_in_same_genus(&e8e8, &d16));

        // are_in_same_genus ⟺ (equal signatures ∧ isomorphic discriminant forms)
        // across the even-lattice zoo.
        let zoo = [
            a_n(1),
            a_n(2),
            a_n(3),
            a_n(1).direct_sum(&a_n(1)),
            d_n(4),
            e_6(),
            e_7(),
            e_8(),
        ];
        for (i, a) in zoo.iter().enumerate() {
            for b in &zoo[i..] {
                assert_eq!(
                    are_in_same_genus(a, b),
                    nikulin_rhs(a, b),
                    "Nikulin equivalence failed for a pair"
                );
            }
        }
    }

    #[test]
    fn a1_discriminant_form_has_quarter_turn_phase() {
        let a1 = a_n(1);
        let disc = DiscriminantForm::from_lattice(&a1).unwrap();
        assert_eq!(disc.group, vec![2]);
        assert_eq!(disc.reps.len(), 2);
        assert_eq!(disc.quadratic_value_mod2(&[1]), Rational::new(1, 2));
        assert_eq!(disc.milgram_signature_mod8(), Some(1));
        assert_eq!(disc.weil_s_prefactor_phase_mod8(), Some(7));
        assert_eq!(disc.weil_s_recovers_milgram_phase_mod8(), Some(1));
        assert!(disc.verify_weil_relations());
        assert_eq!(verify_milgram(&a1), Some(true));
    }

    #[test]
    fn ade_root_lattices_match_milgram_phase() {
        for n in 1..=5 {
            let a = a_n(n);
            let disc = DiscriminantForm::from_lattice(&a).unwrap();
            assert_eq!(disc.group, vec![n as i128 + 1]);
            assert_eq!(disc.milgram_signature_mod8_fqm(), Some(n as i128 % 8));
            assert_eq!(disc.milgram_signature_mod8(), Some(n as i128 % 8));
            assert!(disc.verify_weil_relations(), "Weil relations A_{n}");
            assert_eq!(verify_milgram(&a), Some(true), "A_{n}");
        }

        let d4 = d_n(4);
        let disc = DiscriminantForm::from_lattice(&d4).unwrap();
        assert_eq!(disc.group, vec![2, 2]);
        assert_eq!(disc.milgram_signature_mod8_fqm(), Some(4));
        assert_eq!(disc.milgram_signature_mod8(), Some(4));
        let gs = disc.gauss_sum();
        assert!((gs.re + 1.0).abs() < 1e-8 && gs.im.abs() < 1e-8);
        assert_eq!(disc.weil_s_recovers_milgram_phase_mod8(), Some(4));
        assert!(disc.verify_weil_relations());
        assert_eq!(verify_milgram(&d4), Some(true));
    }

    #[test]
    fn e8_is_unimodular_and_milgram_trivial() {
        let e8 = e_8();
        let disc = DiscriminantForm::from_lattice(&e8).unwrap();
        assert!(disc.group.is_empty());
        assert_eq!(disc.reps, vec![vec![0; 8]]);
        assert_eq!(disc.milgram_signature_mod8(), Some(0));
        assert_eq!(disc.weil_t(), vec![Complex64::one()]);
        assert_eq!(disc.weil_s().unwrap(), vec![vec![Complex64::one()]]);
        assert!(disc.verify_weil_relations());
        assert_eq!(verify_milgram(&e8), Some(true));

        let e8e8 = e8.direct_sum(&e8);
        assert_eq!(
            DiscriminantForm::from_lattice(&e8e8)
                .unwrap()
                .milgram_signature_mod8_fqm(),
            Some(0)
        );
        assert_eq!(
            DiscriminantForm::from_lattice(&e8e8)
                .unwrap()
                .milgram_signature_mod8(),
            Some(0)
        );
        assert_eq!(verify_milgram(&e8e8), Some(true));
    }

    #[test]
    fn fqm_gauss_phase_reports_primary_factors() {
        let a1a2 = a_n(1).direct_sum(&a_n(2));
        let disc = DiscriminantForm::from_lattice(&a1a2).unwrap();
        let phase = disc.fqm_gauss_phase().unwrap();
        assert_eq!(phase.order, 6);
        assert_eq!(phase.phase_mod8, 3);
        assert_eq!(
            phase.primary,
            vec![
                FqmPrimaryPhase {
                    prime: 2,
                    order: 2,
                    exponent: 2,
                    phase_mod8: 1,
                },
                FqmPrimaryPhase {
                    prime: 3,
                    order: 3,
                    exponent: 3,
                    phase_mod8: 2,
                },
            ]
        );
    }

    #[test]
    fn fqm_phase_extends_past_2_elementary_brown_slice() {
        // A_3 has discriminant group Z/4, so the old 2-elementary Brown bridge
        // declines. The p-primary FQM phase still sees the Milgram signature.
        let a3 = DiscriminantForm::from_lattice(&a_n(3)).unwrap();
        assert_eq!(a3.group, vec![4]);
        assert_eq!(a3.brown_invariant(), None);
        assert_eq!(a3.milgram_signature_mod8_fqm(), Some(3));
        assert_eq!(a3.fqm_gauss_phase().unwrap().primary[0].prime, 2);

        // E_6 is odd torsion (Z/3): outside Brown's char-2 slice, inside the FQM
        // Gauss phase projection.
        let e6 = DiscriminantForm::from_lattice(&e_6()).unwrap();
        assert_eq!(e6.group, vec![3]);
        assert_eq!(e6.brown_invariant(), None);
        assert_eq!(e6.milgram_signature_mod8_fqm(), Some(6));
        assert_eq!(e6.fqm_gauss_phase().unwrap().primary[0].prime, 3);
    }

    #[test]
    fn fqm_phase_matches_signature_genus_and_float_oracle_on_zoo() {
        let zoo = [
            a_n(1),
            a_n(2),
            a_n(3),
            a_n(4),
            a_n(5),
            d_n(4),
            d_n(5),
            d_n(8),
            e_6(),
            e_7(),
            e_8(),
        ];
        for l in zoo {
            let disc = DiscriminantForm::from_lattice(&l).unwrap();
            let fqm = disc.milgram_signature_mod8_fqm().unwrap();
            let float = disc.milgram_signature_mod8().unwrap();
            let (pos, neg) = l.signature();
            let sig = (pos as i128 - neg as i128).rem_euclid(8);
            assert_eq!(fqm, sig, "FQM phase mismatch for group {:?}", disc.group);
            assert_eq!(
                float, sig,
                "float phase mismatch for group {:?}",
                disc.group
            );
            assert_eq!(genus_signature_mod8(&l), Some(sig), "genus route mismatch");
            assert_eq!(verify_milgram(&l), Some(true), "Milgram verifier mismatch");
        }
    }

    #[test]
    fn brown_invariant_recovers_signature_mod8_on_2_elementary_forms() {
        // β ≡ sign(L) mod 8 — the fifth route to σ mod 8, exact-integer (Bridge M).
        // 2-elementary generators: A_1 (ℤ/2, β=1), E_7 (ℤ/2, β=7), D_4 ((ℤ/2)², β=4),
        // D_8 ((ℤ/2)², β=0), and the unimodular E_8 (β=0).
        for (l, want) in [
            (a_n(1), 1u128),
            (e_7(), 7),
            (d_n(4), 4),
            (d_n(8), 0),
            (e_8(), 0),
        ] {
            let disc = DiscriminantForm::from_lattice(&l).unwrap();
            let brown = disc.brown_invariant().expect("2-elementary");
            assert_eq!(brown.beta, want, "β mismatch");
            assert_eq!(brown.radical_dim, 0, "discriminant b is nondegenerate");
            // cross-check against the shipped f64 Milgram phase.
            let milgram = disc.milgram_signature_mod8().unwrap().rem_euclid(8) as u128;
            assert_eq!(brown.beta, milgram, "β ≢ Milgram phase");
        }
    }

    #[test]
    fn brown_invariant_is_none_off_the_2_elementary_slice() {
        // A_2 has discriminant group ℤ/3 (odd torsion); A_3 has ℤ/4 (exponent 4).
        // Neither is 2-elementary — the Brown slice declines, honestly.
        assert_eq!(
            DiscriminantForm::from_lattice(&a_n(2))
                .unwrap()
                .brown_invariant(),
            None
        );
        assert_eq!(
            DiscriminantForm::from_lattice(&a_n(3))
                .unwrap()
                .brown_invariant(),
            None
        );
        // E_6 has discriminant group ℤ/3 as well.
        assert_eq!(
            DiscriminantForm::from_lattice(&e_6())
                .unwrap()
                .brown_invariant(),
            None
        );
    }

    #[test]
    fn odd_lattices_have_no_even_discriminant_quadratic_form() {
        assert!(DiscriminantForm::from_lattice(&IntegralForm::diagonal(&[1])).is_none());
    }

    #[test]
    fn odd_discriminant_form_uses_q_mod_one() {
        let z = IntegralForm::diagonal(&[1]);
        let zd = OddDiscriminantForm::from_lattice(&z).unwrap();
        assert!(zd.group.is_empty());
        assert_eq!(zd.reps, vec![vec![0]]);
        assert_eq!(zd.quadratic_value_mod1(&[0]), Rational::zero());
        assert_eq!(zd.gauss_phase_mod8(), Some(0));
        assert!(DiscriminantForm::from_lattice(&z).is_none());

        let three = IntegralForm::diagonal(&[3]);
        let od = OddDiscriminantForm::from_lattice(&three).unwrap();
        assert_eq!(od.group, vec![3]);
        assert_eq!(od.quadratic_value_mod1(&[1]), Rational::new(1, 3));
        assert_eq!(od.quadratic_value_mod1(&[2]), Rational::new(1, 3));
        assert_eq!(od.bilinear_value_mod1(&[1], &[1]), Rational::new(1, 3));
        assert_eq!(od.gauss_phase_mod8(), Some(2));
    }

    #[test]
    fn odd_milgram_report_matches_signature_with_oddity_correction() {
        let cases = [
            IntegralForm::diagonal(&[1]),
            IntegralForm::diagonal(&[3]),
            IntegralForm::diagonal(&[1, 2]),
            IntegralForm::diagonal(&[1]).direct_sum(&e_8()),
        ];
        for l in cases {
            let report = odd_milgram_report(&l).unwrap();
            assert!(report.verified(), "odd Milgram report failed: {report:?}");
            assert_eq!(verify_odd_milgram(&l), Some(true));
        }

        let z = odd_milgram_report(&IntegralForm::diagonal(&[1])).unwrap();
        assert_eq!(z.oddity_mod8, 1);
        assert_eq!(z.p_excess_mod8, 0);
        assert_eq!(z.corrected_signature_mod8, 1);

        let three = odd_milgram_report(&IntegralForm::diagonal(&[3])).unwrap();
        assert_eq!(three.oddity_mod8, 3);
        assert_eq!(three.p_excess_mod8, 2);
        assert_eq!(three.corrected_signature_mod8, 1);
        assert_eq!(verify_odd_milgram(&e_8()), None);
    }

    #[test]
    fn odd_construction_a_witnesses_type_i_unimodular_lattices() {
        let z2_code = repetition_code(2).unwrap();
        assert!(z2_code.is_self_dual());
        assert!(!z2_code.is_doubly_even());
        let z2 = z2_code.construction_a().unwrap();
        assert_eq!(z2.dim(), 2);
        assert!(z2.is_unimodular());
        assert!(!z2.is_even());
        assert_eq!(z2.minimum(), Some(1));
        assert_eq!(z2.kissing_number(), Some(4));
        assert_eq!(verify_odd_milgram(&z2), Some(true));

        let z2_e8 = type_i_z2_plus_e8_code().construction_a().unwrap();
        assert_eq!(z2_e8.dim(), 10);
        assert!(z2_e8.is_unimodular());
        assert!(!z2_e8.is_even());
        assert_eq!(z2_e8.minimum(), Some(1));
        assert_eq!(z2_e8.kissing_number(), Some(4));
        assert_eq!(verify_odd_milgram(&z2_e8), Some(true));
    }
}
