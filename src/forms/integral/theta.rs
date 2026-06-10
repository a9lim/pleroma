//! Theta series of positive-definite even integral lattices.
//!
//! For an even lattice `L`, the scalar theta series is
//! `theta_L = sum_v q^{Q(v)/2}`. The implementation enumerates short vectors up
//! to the requested norm bound and buckets them exactly by `Q/2`.

use super::lattice::IntegralForm;

impl IntegralForm {
    /// The first `terms` coefficients of `theta_L(q) = sum_m r_L(m) q^m`, where
    /// `r_L(m) = #{v in L : Q(v) = 2m}`.
    ///
    /// Returns `None` unless the lattice is positive definite and even. The
    /// constant term is the zero vector; [`IntegralForm::short_vectors`] supplies
    /// the nonzero terms.
    pub fn theta_series(&self, terms: usize) -> Option<Vec<i128>> {
        if terms == 0 {
            return Some(Vec::new());
        }
        if !self.is_even() || !self.is_positive_definite() {
            return None;
        }
        let mut out = vec![0i128; terms];
        out[0] = 1;
        let bound = 2i128
            .checked_mul((terms - 1) as i128)
            .expect("theta bound exceeds i128");
        for v in self.short_vectors(bound)? {
            let q = self.norm(&v);
            debug_assert_eq!(q % 2, 0);
            let idx = usize::try_from(q / 2).ok()?;
            if idx < terms {
                out[idx] += 1;
            }
        }
        Some(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forms::{
        as_modular_form, delta, e_8, eisenstein_e4, leech, mass_even_unimodular, modular_qexp_mul,
        modular_qexp_scale, modular_qexp_sub, qexp_from_i128, type_ii_e8_sum_code,
        type_ii_len16_code, D16_PLUS_AUT_ORDER, E8_WEYL_GROUP_ORDER,
    };
    use crate::scalar::{Rational, Scalar};

    #[test]
    fn theta_series_respects_even_positive_boundary() {
        assert_eq!(
            IntegralForm::diagonal(&[2]).theta_series(5),
            Some(vec![1, 2, 0, 0, 2])
        );
        assert!(IntegralForm::diagonal(&[1]).theta_series(3).is_none());
        assert!(IntegralForm::diagonal(&[2, -2]).theta_series(3).is_none());
    }

    #[test]
    fn theta_series_counts_the_headline_lattices() {
        assert_eq!(e_8().theta_series(3), Some(vec![1, 240, 2160]));
        assert_eq!(
            type_ii_e8_sum_code().theta_series_via_weight_enumerator(2),
            Some(vec![1, 480])
        );
        assert_eq!(
            type_ii_len16_code().theta_series_via_weight_enumerator(2),
            Some(vec![1, 480])
        );
        assert_eq!(leech().theta_series(2), Some(vec![1, 0]));
    }

    #[test]
    fn theta_series_identifies_the_full_modular_forms() {
        let terms = 3;
        let e4 = eisenstein_e4(terms);
        let e8_theta = qexp_from_i128(&e_8().theta_series(terms).unwrap());
        assert_eq!(e8_theta, e4);
        assert_eq!(
            as_modular_form(&e8_theta, 4, terms),
            Some(vec![Rational::one()])
        );

        let e4_squared = modular_qexp_mul(&e4, &e4, 3);
        let split = qexp_from_i128(
            &type_ii_e8_sum_code()
                .theta_series_via_weight_enumerator(3)
                .unwrap(),
        );
        let d16 = qexp_from_i128(
            &type_ii_len16_code()
                .theta_series_via_weight_enumerator(3)
                .unwrap(),
        );
        assert_eq!(split, e4_squared);
        assert_eq!(d16, e4_squared);
        assert_eq!(as_modular_form(&d16, 8, 3), Some(vec![Rational::one()]));

        // Rank 16 Siegel-Weil is degenerate but consistent: the two class
        // representatives already have the same theta series, so the
        // mass-weighted average is the same `E4^2`.
        //
        // Pin |Aut(E8⊕E8)| = 2·|W(E8)|² (factor 2 from the swap automorphism).
        let w2 = E8_WEYL_GROUP_ORDER
            .checked_mul(E8_WEYL_GROUP_ORDER)
            .expect("|W(E8)|^2 exceeds u128");
        assert_eq!(w2, 485_432_135_516_160_000);
        let aut_e8_e8 = 2u128
            .checked_mul(w2)
            .expect("|Aut(E8+E8)| = 2·|W(E8)|^2 exceeds u128");
        assert_eq!(aut_e8_e8, 970_864_271_032_320_000);
        assert_eq!(D16_PLUS_AUT_ORDER, 685_597_979_049_984_000);
    }

    #[test]
    fn siegel_weil_rank16_mass_identity_is_exact() {
        // For the rank-16 even-unimodular genus with two classes (E8⊕E8 and D16+),
        // the Siegel-Weil identity requires:
        //   1/|Aut(E8⊕E8)| + 1/|Aut(D16+)| = mass_even_unimodular(16).
        // |Aut(E8⊕E8)| = 2·|W(E8)|² because the two summands can be swapped.
        let w = E8_WEYL_GROUP_ORDER;
        let aut_e8_e8 = 2u128.checked_mul(w).unwrap().checked_mul(w).unwrap();
        let d = D16_PLUS_AUT_ORDER;
        let (mass_num, mass_den) = mass_even_unimodular(16).unwrap();
        // Clear denominators: mass_num * aut_e8_e8 * d == (d + aut_e8_e8) * mass_den
        // (all values fit in i128 after the lcm reduction in mass_even_unimodular).
        let lhs = mass_num
            .checked_mul(aut_e8_e8 as i128)
            .and_then(|x| x.checked_mul(d as i128));
        let rhs = (d as i128)
            .checked_add(aut_e8_e8 as i128)
            .and_then(|s| s.checked_mul(mass_den));
        assert_eq!(
            lhs, rhs,
            "Siegel-Weil: 1/|Aut(E8+E8)| + 1/|Aut(D16+)| != mass(16)"
        );
    }

    #[test]
    fn leech_theta_is_pinned_by_rootlessness_in_weight_12() {
        let terms = 2;
        let e4 = eisenstein_e4(terms);
        let e4_cubed = modular_qexp_mul(&modular_qexp_mul(&e4, &e4, terms), &e4, terms);
        let leech_form = modular_qexp_sub(
            &e4_cubed,
            &modular_qexp_scale(&delta(terms), Rational::int(720), terms),
            terms,
        );
        let leech_theta = qexp_from_i128(&leech().theta_series(terms).unwrap());
        assert_eq!(leech_theta, leech_form);
        assert_eq!(
            as_modular_form(&leech_theta, 12, terms),
            Some(vec![Rational::new(7, 12), Rational::new(5, 12)])
        );
    }
}
