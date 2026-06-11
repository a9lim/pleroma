//! Relation and certificate record types: [`GameRelation`],
//! [`GameRelationCertificate`], and [`RelationSearchCertificate`].

use crate::games::partizan::Game;
use crate::linalg::integer::reduce_integer_vector;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GameRelation {
    pub coeffs: Vec<i128>,
}

impl GameRelation {
    pub fn new(coeffs: Vec<i128>) -> Self {
        assert!(
            coeffs.iter().any(|&c| c != 0),
            "game relation must be nonzero"
        );
        GameRelation { coeffs }
    }
}

/// A stored witness for an accepted game-group relation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GameRelationCertificate {
    /// The relation vector `Σ coeffs[i]·g_i = 0`.
    pub coeffs: Vec<i128>,
    /// Canonical value key of the evaluated relation. Accepted relations have
    /// the same value key as [`Game::zero`].
    pub value_key: String,
    /// Whether this row added new information modulo earlier accepted rows.
    pub independent: bool,
}

/// Audit trail for bounded automatic relation discovery.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RelationSearchCertificate {
    /// Coefficient box `[-bound, bound]` used for automatic discovery (`0` for
    /// explicitly supplied relations).
    pub bound: i128,
    /// `true` iff the whole coefficient box was searched.
    pub exhaustive: bool,
    /// Number of nonzero candidates in the coefficient box, if it fit in `usize`.
    pub candidate_count: Option<usize>,
    /// Accepted relation rows, in the order they were imposed.
    pub relations: Vec<GameRelationCertificate>,
}

// ---------------------------------------------------------------------------
// Certificate helpers (pub(super) — only the lambda/clifford constructors call them)
// ---------------------------------------------------------------------------

pub(super) fn relation_search_certificate(
    gens: &[Game],
    bound: i128,
    exhaustive: bool,
    candidate_count: Option<usize>,
    relations: &[GameRelation],
    independent: bool,
) -> RelationSearchCertificate {
    RelationSearchCertificate {
        bound,
        exhaustive,
        candidate_count,
        relations: relation_certificates(gens, relations, independent),
    }
}

pub(super) fn relation_certificates(
    gens: &[Game],
    relations: &[GameRelation],
    trust_independent: bool,
) -> Vec<GameRelationCertificate> {
    let mut previous = Vec::new();
    relations
        .iter()
        .map(|rel| {
            let mut reduced = rel.coeffs.clone();
            reduce_integer_vector(&mut reduced, previous.clone());
            let independent = trust_independent || reduced.iter().any(|&c| c != 0);
            previous.push(rel.coeffs.clone());
            GameRelationCertificate {
                coeffs: rel.coeffs.clone(),
                value_key: eval_relation(gens, &rel.coeffs).canonical_string(),
                independent,
            }
        })
        .collect()
}

pub(super) fn eval_relation(gens: &[Game], coeffs: &[i128]) -> Game {
    let mut acc = Game::zero();
    for (g, &c) in gens.iter().zip(coeffs) {
        acc = acc.add(&g.times_int(c));
    }
    acc
}
