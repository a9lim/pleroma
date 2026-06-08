//! The multivector engine, generic over any `Scalar` backend.
//!
//! ## Metric data — characteristic-faithful by design
//!
//! A blade is a `u128` bitmask over basis generators e_0..e_127. The algebra is
//! defined by two independent pieces of data, *not* a single bilinear form:
//!
//!   * `q[i]`      = e_i²                      (the quadratic form / squares)
//!   * `b[(i,j)]`  = e_i e_j + e_j e_i  (i<j)  (the polar / anticommutator form)
//!
//! In characteristic ≠ 2 these are linked (`b = 2·offdiag`, `q = diag`), so an
//! orthogonal basis just sets `b = 0`. In characteristic 2 they are genuinely
//! independent: the polar form is *alternating* (`b(i,i)=0`) yet `q[i]` can be
//! nonzero, and a nonzero off-diagonal `b[(i,j)]` is exactly what makes the
//! nim-Clifford algebra *non-commutative*. Carrying both is the faithful thing.
//!
//! "With nilpotents": set `q[i] = 0` and you get a null generator, e_i² = 0.
//! All `q = 0`, all `b = 0` ⇒ the exterior/Grassmann algebra.
//!
//! ## Product
//!
//! Two canonical blades multiply by concatenating their (ascending) generator
//! lists into a word and reducing to canonical form with the rules
//!   e_i e_i  → q[i]                            (equal adjacent: contract)
//!   e_i e_j  → b[(j,i)] − e_j e_i   (i>j)      (out of order: swap, emit polar)
//! The `−` goes through the scalar's own `neg()`, so in characteristic 2 it is
//! `+` automatically and signs vanish — no special-casing. Termination: each
//! step lowers (word length, inversion count) lexicographically.

mod algebra;
mod basis;
mod inverse;
mod metric;
mod multivector;
mod product;
mod terms;

pub use algebra::CliffordAlgebra;
pub use basis::{bits, grade, MAX_BASIS_DIM};
pub use metric::Metric;
pub use multivector::Multivector;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scalar::{Nimber, Rational, Scalar, Surreal};
    use std::collections::BTreeMap;

    fn r(n: i128) -> Rational {
        Rational::int(n)
    }

    #[test]
    #[should_panic(expected = "at most 128 generators")]
    fn algebra_dimension_must_fit_blade_mask() {
        let _ = CliffordAlgebra::new(129, Metric::<Rational>::grassmann(129));
    }

    #[test]
    #[should_panic(expected = "b-keys must satisfy i < j")]
    fn metric_rejects_reversed_or_diagonal_polar_keys() {
        let mut b = BTreeMap::new();
        b.insert((1usize, 0usize), r(1));
        let _ = Metric::new(vec![r(1), r(1)], b);
    }

    #[test]
    #[should_panic(expected = "generator index 2 out of range")]
    fn generator_index_must_be_in_the_algebra() {
        let alg = CliffordAlgebra::new(2, Metric::diagonal(vec![r(1), r(1)]));
        let _ = alg.gen(2);
    }

    #[test]
    fn complex_numbers_cl01() {
        let alg = CliffordAlgebra::new(1, Metric::diagonal(vec![r(-1)]));
        assert_eq!(alg.mul(&alg.gen(0), &alg.gen(0)), alg.scalar(r(-1)));
    }

    #[test]
    fn cl20_bivector_squares_to_minus_one() {
        let alg = CliffordAlgebra::new(2, Metric::diagonal(vec![r(1), r(1)]));
        let e0 = alg.gen(0);
        let e1 = alg.gen(1);
        let e0e1 = alg.mul(&e0, &e1);
        let e1e0 = alg.mul(&e1, &e0);
        assert_eq!(e0e1, alg.scalar_mul(&r(-1), &e1e0));
        assert_eq!(alg.mul(&e0e1, &e0e1), alg.scalar(r(-1)));
    }

    #[test]
    fn orthogonal_blade_product_handles_repeated_indices_directly() {
        let metric = Metric::diagonal(vec![r(2), r(3), r(5)]);
        let mut expect = BTreeMap::new();
        // e_0e_1 · e_1e_2 = q_1 e_0e_2.
        expect.insert(0b101, r(3));
        assert_eq!(metric.geom_product_blades(0b011, 0b110), expect);

        let mut expect_scalar = BTreeMap::new();
        // e_0e_1 · e_0e_1 = -q_0q_1.
        expect_scalar.insert(0, r(-6));
        assert_eq!(metric.geom_product_blades(0b011, 0b011), expect_scalar);
    }

    #[test]
    fn grassmann_generators_are_nilpotent() {
        let alg = CliffordAlgebra::new(3, Metric::grassmann(3));
        for i in 0..3 {
            let ei = alg.gen(i);
            assert!(alg.mul(&ei, &ei).is_zero(), "e{i}^2 should be 0");
        }
        let (e0, e1) = (alg.gen(0), alg.gen(1));
        assert_eq!(alg.mul(&e0, &e1), alg.wedge(&e0, &e1));
        assert_eq!(
            alg.mul(&e0, &e1),
            alg.scalar_mul(&r(-1), &alg.mul(&e1, &e0))
        );
    }

    #[test]
    fn nimber_orthogonal_is_commutative() {
        let alg = CliffordAlgebra::new(2, Metric::diagonal(vec![Nimber(2), Nimber(3)]));
        let e0 = alg.gen(0);
        let e1 = alg.gen(1);
        assert_eq!(alg.mul(&e0, &e1), alg.mul(&e1, &e0));
        assert_eq!(alg.mul(&e0, &e0), alg.scalar(Nimber(2)));
    }

    #[test]
    fn nimber_offdiagonal_is_noncommutative() {
        let mut b = BTreeMap::new();
        b.insert((0usize, 1usize), Nimber(1));
        let alg = CliffordAlgebra::new(2, Metric::new(vec![Nimber(0), Nimber(0)], b));
        let e0 = alg.gen(0);
        let e1 = alg.gen(1);
        let anti = alg.add(&alg.mul(&e0, &e1), &alg.mul(&e1, &e0));
        assert_eq!(anti, alg.scalar(Nimber(1)));
        assert_ne!(alg.mul(&e0, &e1), alg.mul(&e1, &e0));
    }

    fn assert_associative<S: Scalar>(alg: &CliffordAlgebra<S>, gens: &[Multivector<S>]) {
        for a in gens {
            for b in gens {
                for c in gens {
                    assert_eq!(alg.mul(&alg.mul(a, b), c), alg.mul(a, &alg.mul(b, c)));
                }
            }
        }
    }

    #[test]
    fn associativity_rational_nonorthogonal() {
        let mut b = BTreeMap::new();
        b.insert((0usize, 1usize), r(1));
        b.insert((1usize, 2usize), r(-1));
        let alg = CliffordAlgebra::new(3, Metric::new(vec![r(1), r(-1), r(2)], b));
        let gens = [
            alg.gen(0),
            alg.gen(1),
            alg.gen(2),
            alg.mul(&alg.gen(0), &alg.gen(1)),
            alg.add(&alg.gen(0), &alg.scalar(r(3))),
        ];
        assert_associative(&alg, &gens);
    }

    #[test]
    fn associativity_nimber_nonorthogonal() {
        let mut b = BTreeMap::new();
        b.insert((0usize, 1usize), Nimber(1));
        b.insert((0usize, 2usize), Nimber(3));
        let alg = CliffordAlgebra::new(3, Metric::new(vec![Nimber(2), Nimber(1), Nimber(0)], b));
        let gens = [
            alg.gen(0),
            alg.gen(1),
            alg.gen(2),
            alg.mul(&alg.gen(0), &alg.gen(1)),
            alg.add(&alg.gen(2), &alg.scalar(Nimber(5))),
        ];
        assert_associative(&alg, &gens);
    }

    #[test]
    fn vector_inverse() {
        let alg = CliffordAlgebra::new(3, Metric::diagonal(vec![r(1), r(1), r(1)]));
        let v = alg.gen(0);
        let vi = alg.versor_inverse(&v).unwrap();
        assert_eq!(alg.mul(&v, &vi), alg.scalar(r(1)));
        assert_eq!(vi, v);

        let alg2 = CliffordAlgebra::new(1, Metric::diagonal(vec![r(2)]));
        let e0 = alg2.gen(0);
        assert_eq!(
            alg2.mul(&e0, &alg2.versor_inverse(&e0).unwrap()),
            alg2.scalar(r(1))
        );

        let alg0 = CliffordAlgebra::new(1, Metric::<Rational>::grassmann(1));
        assert!(alg0.versor_inverse(&alg0.gen(0)).is_none());
    }

    #[test]
    fn reflection_fixes_and_negates() {
        let alg = CliffordAlgebra::new(2, Metric::diagonal(vec![r(1), r(1)]));
        let (e0, e1) = (alg.gen(0), alg.gen(1));
        assert_eq!(alg.reflect(&e1, &e0).unwrap(), e0);
        assert_eq!(alg.reflect(&e1, &e1).unwrap(), alg.scalar_mul(&r(-1), &e1));
    }

    #[test]
    fn rotor_preserves_norm() {
        let alg = CliffordAlgebra::new(2, Metric::diagonal(vec![r(1), r(1)]));
        let rotor = alg.mul(&alg.gen(0), &alg.gen(1));
        let x = alg.add(
            &alg.scalar_mul(&r(3), &alg.gen(0)),
            &alg.scalar_mul(&r(4), &alg.gen(1)),
        );
        let rx = alg.sandwich(&rotor, &x).unwrap();
        assert_eq!(alg.norm2(&rx), alg.norm2(&x));
    }

    #[test]
    fn twisted_adjoint_matches_reflect_and_sandwich() {
        let alg = CliffordAlgebra::new(2, Metric::diagonal(vec![r(1), r(1)]));
        let (e0, e1) = (alg.gen(0), alg.gen(1));
        let x = alg.add(&alg.scalar_mul(&r(3), &e0), &alg.scalar_mul(&r(4), &e1));
        assert_eq!(
            alg.twisted_sandwich(&e1, &x).unwrap(),
            alg.reflect(&e1, &x).unwrap()
        );
        let rotor = alg.mul(&e0, &e1);
        assert_eq!(
            alg.twisted_sandwich(&rotor, &x).unwrap(),
            alg.sandwich(&rotor, &x).unwrap()
        );
    }

    #[test]
    fn left_contraction_lowers_grade() {
        let alg = CliffordAlgebra::new(3, Metric::diagonal(vec![r(1), r(1), r(1)]));
        let e0 = alg.gen(0);
        let e0e1 = alg.mul(&alg.gen(0), &alg.gen(1));
        assert_eq!(alg.left_contract(&e0, &e0e1), alg.gen(1));
        let three = alg.scalar(r(3));
        assert_eq!(
            alg.left_contract(&three, &e0e1),
            alg.scalar_mul(&r(3), &e0e1)
        );
    }

    #[test]
    fn dual_of_vector_is_bivector_in_3d() {
        let alg = CliffordAlgebra::new(3, Metric::diagonal(vec![r(1), r(1), r(1)]));
        let d = alg.dual(&alg.gen(0)).unwrap();
        assert!(!d.is_zero());
        assert_eq!(alg.grade_part(&d, 2), d);
    }

    #[test]
    fn grade_involution_signs() {
        let alg = CliffordAlgebra::new(2, Metric::diagonal(vec![r(1), r(1)]));
        let v = alg.add(
            &alg.scalar(r(5)),
            &alg.add(&alg.gen(0), &alg.mul(&alg.gen(0), &alg.gen(1))),
        );
        let expect = alg.add(
            &alg.scalar(r(5)),
            &alg.add(
                &alg.scalar_mul(&r(-1), &alg.gen(0)),
                &alg.mul(&alg.gen(0), &alg.gen(1)),
            ),
        );
        assert_eq!(alg.grade_involution(&v), expect);
    }

    #[test]
    fn versor_over_surreal_metric() {
        let alg = CliffordAlgebra::new(
            2,
            Metric::diagonal(vec![Surreal::omega(), Surreal::epsilon()]),
        );
        let e0 = alg.gen(0);
        let inv = alg.versor_inverse(&e0).unwrap();
        assert_eq!(alg.mul(&e0, &inv), alg.scalar(Surreal::one()));
    }

    #[test]
    fn even_subalgebra_of_cl30_is_quaternions() {
        let alg = CliffordAlgebra::new(3, Metric::diagonal(vec![r(1), r(1), r(1)]));
        let even = alg.even_subalgebra().unwrap();
        assert_eq!(even.dim, 2);
        let (f0, f1) = (even.gen(0), even.gen(1));
        assert_eq!(even.mul(&f0, &f0), even.scalar(r(-1)));
        assert_eq!(even.mul(&f1, &f1), even.scalar(r(-1)));
        assert_eq!(
            even.mul(&f0, &f1),
            even.scalar_mul(&r(-1), &even.mul(&f1, &f0))
        );
    }

    #[test]
    fn even_part_projection() {
        let alg = CliffordAlgebra::new(2, Metric::diagonal(vec![r(1), r(1)]));
        let v = alg.add(
            &alg.scalar(r(5)),
            &alg.add(
                &alg.scalar_mul(&r(2), &alg.gen(0)),
                &alg.mul(&alg.gen(0), &alg.gen(1)),
            ),
        );
        let expect = alg.add(&alg.scalar(r(5)), &alg.mul(&alg.gen(0), &alg.gen(1)));
        assert_eq!(alg.even_part(&v), expect);
    }

    #[test]
    fn graded_tensor_blocks_are_orthogonal() {
        let left = CliffordAlgebra::new(1, Metric::diagonal(vec![r(1)]));
        let right = CliffordAlgebra::new(1, Metric::diagonal(vec![r(-1)]));
        let alg = left.graded_tensor(&right);
        assert_eq!(alg.dim, 2);
        let e0 = alg.gen(0);
        let e1 = alg.gen(1);
        assert_eq!(alg.mul(&e0, &e0), alg.scalar(r(1)));
        assert_eq!(alg.mul(&e1, &e1), alg.scalar(r(-1)));
        assert_eq!(
            alg.mul(&e0, &e1),
            alg.scalar_mul(&r(-1), &alg.mul(&e1, &e0))
        );
        assert_eq!(alg.embed_first(&left.gen(0)), e0);
        assert_eq!(alg.embed_second(&right.gen(0), left.dim), e1);
    }

    #[test]
    fn general_product_reproduces_reduce_word_when_a_empty() {
        let mut b = BTreeMap::new();
        b.insert((0usize, 1usize), r(1));
        b.insert((1usize, 2usize), r(-1));
        b.insert((0usize, 2usize), r(2));
        let m = Metric::new(vec![r(1), r(-1), r(2)], b);
        for ba in 0u128..8 {
            for bb in 0u128..8 {
                let word: Vec<usize> = bits(ba).into_iter().chain(bits(bb)).collect();
                assert_eq!(
                    m.geom_product_blades(ba, bb),
                    m.reduce_word(&word),
                    "mismatch on blades {ba:#b}·{bb:#b}"
                );
            }
        }
    }

    #[test]
    fn general_bilinear_in_order_contraction() {
        let mut a = BTreeMap::new();
        a.insert((0usize, 1usize), r(5));
        let alg = CliffordAlgebra::new(2, Metric::general(vec![r(1), r(1)], BTreeMap::new(), a));
        let (e0, e1) = (alg.gen(0), alg.gen(1));
        let blade = alg.wedge(&e0, &e1);
        assert_eq!(alg.mul(&e0, &e1), alg.add(&blade, &alg.scalar(r(5))));
        assert_eq!(alg.add(&alg.mul(&e0, &e1), &alg.mul(&e1, &e0)), alg.zero());
    }

    #[test]
    fn associativity_general_bilinear_form() {
        let mut b = BTreeMap::new();
        b.insert((0usize, 1usize), r(1));
        b.insert((1usize, 2usize), r(2));
        let mut a = BTreeMap::new();
        a.insert((0usize, 1usize), r(3));
        a.insert((0usize, 2usize), r(-1));
        let alg = CliffordAlgebra::new(3, Metric::general(vec![r(2), r(-1), r(1)], b, a));
        let gens = [
            alg.gen(0),
            alg.gen(1),
            alg.gen(2),
            alg.mul(&alg.gen(0), &alg.gen(1)),
            alg.add(&alg.gen(2), &alg.scalar(r(3))),
        ];
        assert_associative(&alg, &gens);

        let mut bn = BTreeMap::new();
        bn.insert((0usize, 1usize), Nimber(1));
        let mut an = BTreeMap::new();
        an.insert((0usize, 1usize), Nimber(2));
        an.insert((1usize, 2usize), Nimber(3));
        let algn = CliffordAlgebra::new(
            3,
            Metric::general(vec![Nimber(2), Nimber(1), Nimber(0)], bn, an),
        );
        let gensn = [
            algn.gen(0),
            algn.gen(1),
            algn.gen(2),
            algn.mul(&algn.gen(0), &algn.gen(1)),
        ];
        assert_associative(&algn, &gensn);
    }
}
