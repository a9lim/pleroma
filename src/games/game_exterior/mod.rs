//! The **exterior algebra of the game group**: `Λ` over `ℤ` on a chosen tuple of
//! games. This is the Clifford-adjacent structure that lives on *all* of
//! game-world — not just the field-like numbers — because the partizan games form
//! an abelian group (a `ℤ`-module), and the Grassmann algebra is the exterior
//! algebra of that module.
//!
//! Three layers, re-exported flat so every public path is unchanged:
//!
//!   * [`relations`] — [`GameRelation`], [`GameRelationCertificate`],
//!     [`RelationSearchCertificate`]: the relation and certificate record types.
//!   * [`lambda`] — [`GameExterior`]: the free Grassmann engine quotiented by
//!     integer game relations such as `2⋆=0`.
//!   * [`clifford`] — [`GameCliffordError`] and [`GameClifford`]: the checked
//!     integer-valued Clifford deformation surface; constructors verify that every
//!     game relation is null and polar-radical before accepting the metric.
//!
//! Generators may be non-numbers (`⋆`, `↑`, switches) — exactly where the
//! Clifford/scalar story cannot go — which is the point: the
//! [`Game`](crate::games::Game) group is not a ring, but it *is* a `ℤ`-module,
//! and that is enough for `Λ`. The stronger question of a natural game-native
//! source for the quadratic data remains open in `docs/OPEN.md`.

pub mod clifford;
pub mod lambda;
pub mod relations;

pub use clifford::*;
pub use lambda::*;
pub use relations::*;
#[cfg(test)]
mod tests {
    use super::*;

    use crate::games::Game;
    use crate::scalar::Integer;
    use std::collections::BTreeMap;

    #[test]
    fn exterior_algebra_lives_on_non_numbers() {
        // Generators that are NOT numbers — exactly where the Clifford/scalar story
        // cannot go — yet the quotient exterior algebra is well defined on them.
        let ext = GameExterior::new(vec![Game::star(), Game::up(), Game::switch(1, -1)]);
        assert!(!ext.game(0).is_number()); // ⋆
        assert!(!ext.game(1).is_number()); // ↑
        let (e0, e1) = (ext.generator(0), ext.generator(1));
        let alg = ext.algebra();
        // the wedge is antisymmetric and nonzero, but quotient-aware operations
        // still remember that it may carry torsion inherited from ⋆.
        let e01 = ext.wedge(&e0, &e1);
        assert!(!e01.is_zero());
        assert_eq!(e01, alg.scalar_mul(&Integer(-1), &alg.wedge(&e1, &e0)));
        assert!(alg.wedge(&e0, &e0).is_zero()); // e_i ∧ e_i = 0
    }

    #[test]
    fn grade1_is_the_game_group() {
        // Λ¹ → game group is a group homomorphism, recovering disjunctive sum.
        let ext = GameExterior::new(vec![Game::star(), Game::up()]);
        let (e0, e1) = (ext.generator(0), ext.generator(1));
        let alg = ext.algebra();
        // value(e0 + e1) = ⋆ + ↑
        let sum = alg.add(&e0, &e1);
        assert!(ext.value_of_grade1(&sum).eq(&Game::star().add(&Game::up())));
        // value(2·e0) = ⋆ + ⋆ = 0  (the 2-torsion of ⋆ shows up as a relation)
        let two_e0 = alg.scalar_mul(&Integer(2), &e0);
        assert!(ext.value_of_grade1(&two_e0).eq(&Game::zero()));
        // value(e0 − e1) = ⋆ − ↑
        let diff = alg.add(&e0, &alg.scalar_mul(&Integer(-1), &e1));
        assert!(ext
            .value_of_grade1(&diff)
            .eq(&Game::star().add(&Game::up().neg())));
    }

    #[test]
    fn game_relations_propagate_through_the_exterior_ideal() {
        let ext = GameExterior::new(vec![Game::star(), Game::up()]);
        assert!(ext.relations().iter().any(|r| r.coeffs == vec![2, 0]));
        let (star, up) = (ext.generator(0), ext.generator(1));
        let star_wedge_up = ext.wedge(&star, &up);
        assert!(!ext.is_zero(&star_wedge_up));
        assert!(ext.is_zero(&ext.scalar_mul(2, &star_wedge_up)));
    }

    #[test]
    fn duplicate_game_generators_are_quotiented_before_wedging() {
        let ext = GameExterior::new(vec![Game::star(), Game::star()]);
        assert!(ext
            .relations()
            .iter()
            .any(|r| r.coeffs == vec![1, -1] || r.coeffs == vec![-1, 1]));
        let e0 = ext.generator(0);
        let e1 = ext.generator(1);
        assert_eq!(ext.reduce(&e0), ext.reduce(&e1));
        assert!(ext.is_zero(&ext.wedge(&e0, &e1)));
    }

    #[test]
    fn relation_search_finds_three_generator_cross_relations() {
        let star = Game::star();
        let up = Game::up();
        let sum = star.add(&up);
        let ext = GameExterior::with_relation_search(vec![star, up, sum], 1);
        assert!(ext.relation_search_complete());
        assert!(ext
            .relations()
            .iter()
            .any(|r| r.coeffs == vec![1, 1, -1] || r.coeffs == vec![-1, -1, 1]));
        let e0 = ext.generator(0);
        let e1 = ext.generator(1);
        let e2 = ext.generator(2);
        assert_eq!(ext.add(&e0, &e1), e2);
    }

    #[test]
    fn relation_search_certificate_records_the_zero_rows() {
        let star = Game::star();
        let ext = GameExterior::with_relation_search(vec![star.clone(), star], 1);
        let cert = ext.relation_search_certificate();
        let zero_key = Game::zero().canonical_string();
        assert_eq!(cert.bound, 1);
        assert!(cert.exhaustive);
        assert_eq!(cert.candidate_count, Some(8)); // 3^2 - 1
        assert!(cert.relations.iter().all(|r| r.value_key == zero_key));
        assert!(cert.relations.iter().all(|r| r.independent));
        assert!(cert
            .relations
            .iter()
            .any(|r| r.coeffs == vec![1, -1] || r.coeffs == vec![-1, 1]));
    }

    #[test]
    fn explicit_relation_certificate_marks_dependent_rows() {
        let star = Game::star();
        let up = Game::up();
        let ext = GameExterior::with_relations(
            vec![star, up],
            vec![GameRelation::new(vec![2, 0]), GameRelation::new(vec![4, 0])],
        );
        let cert = ext.relation_search_certificate();
        assert_eq!(cert.relations.len(), 2);
        assert!(cert.relations[0].independent);
        assert!(!cert.relations[1].independent);
    }

    #[test]
    fn checked_game_clifford_accepts_free_quadratic_data() {
        let mut b = BTreeMap::new();
        b.insert((0, 1), 3);
        let cl = GameClifford::free(vec![Game::up(), Game::switch(1, -1)], vec![1, 0], b).unwrap();
        let alg = cl.algebra();
        let e0 = cl.generator(0);
        let e1 = cl.generator(1);

        assert_eq!(cl.mul(&e0, &e0), alg.scalar(Integer(1)));
        let anticommutator = cl.add(&cl.mul(&e0, &e1), &cl.mul(&e1, &e0));
        assert_eq!(anticommutator, alg.scalar(Integer(3)));
    }

    #[test]
    fn checked_game_clifford_rejects_torsion_quadratic_data() {
        let rel = GameRelation::new(vec![2, 0]);
        let err = GameClifford::with_quadratic_data(
            vec![Game::star(), Game::up()],
            vec![rel.clone()],
            vec![1, 0],
            BTreeMap::new(),
        )
        .err()
        .unwrap();
        assert!(matches!(
            err,
            GameCliffordError::RelationPolarNonzero {
                relation_index: 0,
                generator: 0,
                value: 4
            }
        ));

        let mut b = BTreeMap::new();
        b.insert((0, 1), 1);
        let err = GameClifford::with_quadratic_data(
            vec![Game::star(), Game::up()],
            vec![rel],
            vec![0, 0],
            b,
        )
        .err()
        .unwrap();
        assert!(matches!(
            err,
            GameCliffordError::RelationPolarNonzero {
                relation_index: 0,
                generator: 1,
                value: 2
            }
        ));
    }

    #[test]
    fn checked_game_clifford_accepts_torsion_vanishings() {
        let cl = GameClifford::with_quadratic_data(
            vec![Game::star(), Game::up()],
            vec![GameRelation::new(vec![2, 0])],
            vec![0, 5],
            BTreeMap::new(),
        )
        .unwrap();
        let star = cl.generator(0);
        let up = cl.generator(1);

        assert!(cl.is_zero(&cl.scalar_mul(2, &star)));
        assert_eq!(cl.mul(&up, &up), cl.algebra().scalar(Integer(5)));
        let star_times_up = cl.mul(&star, &up);
        assert!(!cl.is_zero(&star_times_up));
        assert!(cl.is_zero(&cl.scalar_mul(2, &star_times_up)));
    }

    #[test]
    fn checked_game_clifford_handles_duplicate_generators() {
        let mut incompatible_b = BTreeMap::new();
        incompatible_b.insert((0, 1), 2);
        let err = GameClifford::with_quadratic_data(
            vec![Game::star(), Game::star()],
            vec![GameRelation::new(vec![1, -1])],
            vec![1, 2],
            incompatible_b,
        )
        .err()
        .unwrap();
        assert!(matches!(
            err,
            GameCliffordError::RelationPolarNonzero {
                relation_index: 0,
                generator: 1,
                value: -2
            }
        ));

        let mut compatible_b = BTreeMap::new();
        compatible_b.insert((0, 1), 2);
        let cl = GameClifford::with_quadratic_data(
            vec![Game::star(), Game::star()],
            vec![GameRelation::new(vec![1, -1])],
            vec![1, 1],
            compatible_b,
        )
        .unwrap();
        let e0 = cl.generator(0);
        let e1 = cl.generator(1);
        assert_eq!(cl.reduce(&e0), cl.reduce(&e1));

        let e0e1 = cl.mul(&e0, &e1);
        let one = cl.algebra().scalar(Integer(1));
        assert!(cl.is_zero(&cl.add(&e0e1, &cl.scalar_mul(-1, &one))));
    }

    #[test]
    fn checked_game_clifford_relation_search_finds_torsion() {
        let cl = GameClifford::with_relation_search(
            vec![Game::star(), Game::up()],
            2,
            vec![0, 0],
            BTreeMap::new(),
        )
        .unwrap();
        assert!(cl.relation_search_complete());
        assert!(cl.relations().iter().any(|r| r.coeffs == vec![2, 0]));
        assert!(cl.is_zero(&cl.scalar_mul(2, &cl.generator(0))));
    }
}
