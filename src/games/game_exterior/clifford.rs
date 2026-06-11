//! The checked integer-valued Clifford deformation: [`GameCliffordError`] and
//! [`GameClifford`].

use crate::clifford::{bits, CliffordAlgebra, Metric, Multivector};
use crate::games::partizan::Game;
use crate::linalg::integer::reduce_integer_vector;
use crate::scalar::Integer;
use std::collections::BTreeMap;
use std::fmt;

use super::lambda::{discover_relations, grade_masks, relation_multivector};
use super::relations::{
    eval_relation, relation_search_certificate, GameRelation, RelationSearchCertificate,
};

/// Why a checked game-Clifford deformation was rejected.
///
/// The target here is an integer-valued Clifford deformation on the chosen game
/// subgroup. Relations in the game group are imposed as Clifford-ideal relations,
/// so each relation vector must be both null for `Q` and radical for the polar
/// pairing. Over the torsion-free target `Z`, this is what forces documented
/// vanishings such as `2* = 0` killing every pairing involving `*`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum GameCliffordError {
    QuadraticLength {
        expected: usize,
        got: usize,
    },
    BilinearKeyInvalid {
        i: usize,
        j: usize,
        dim: usize,
    },
    RelationLength {
        relation_index: usize,
        expected: usize,
        got: usize,
    },
    RelationNotZero {
        relation_index: usize,
        value_key: String,
    },
    RelationPolarNonzero {
        relation_index: usize,
        generator: usize,
        value: i128,
    },
    RelationQuadraticNonzero {
        relation_index: usize,
        value: i128,
    },
    ArithmeticOverflow {
        relation_index: usize,
        context: &'static str,
    },
}

impl fmt::Display for GameCliffordError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            GameCliffordError::QuadraticLength { expected, got } => write!(
                f,
                "quadratic diagonal length must match generator count: expected {expected}, got {got}"
            ),
            GameCliffordError::BilinearKeyInvalid { i, j, dim } => write!(
                f,
                "bilinear key ({i},{j}) must satisfy i < j < {dim}"
            ),
            GameCliffordError::RelationLength {
                relation_index,
                expected,
                got,
            } => write!(
                f,
                "game relation #{relation_index} length must match generator count: expected {expected}, got {got}"
            ),
            GameCliffordError::RelationNotZero {
                relation_index,
                value_key,
            } => write!(
                f,
                "game relation #{relation_index} does not evaluate to zero (value {value_key})"
            ),
            GameCliffordError::RelationPolarNonzero {
                relation_index,
                generator,
                value,
            } => write!(
                f,
                "game relation #{relation_index} has nonzero polar pairing with generator {generator}: {value}"
            ),
            GameCliffordError::RelationQuadraticNonzero {
                relation_index,
                value,
            } => write!(
                f,
                "game relation #{relation_index} has nonzero quadratic value: {value}"
            ),
            GameCliffordError::ArithmeticOverflow {
                relation_index,
                context,
            } => write!(
                f,
                "integer overflow while checking game relation #{relation_index} ({context})"
            ),
        }
    }
}

impl std::error::Error for GameCliffordError {}

/// A checked integer-valued Clifford deformation of a game subgroup.
///
/// This is deliberately an engineering object, not a claim that arbitrary games
/// form Clifford scalars. The caller supplies integer quadratic data: diagonal
/// values `Q(e_i)` and off-diagonal polar values `{e_i,e_j}` for `i < j`. The
/// constructor verifies that every imposed game-group relation is null and polar
/// radical for that data before quotienting the ordinary integer Clifford algebra
/// by the generated relation ideal.
#[derive(Clone)]
pub struct GameClifford {
    alg: CliffordAlgebra<Integer>,
    gens: Vec<Game>,
    relations: Vec<GameRelation>,
    relation_search_complete: bool,
    relation_certificate: RelationSearchCertificate,
}

impl GameClifford {
    /// Build using bounded automatic relation discovery, matching
    /// [`GameExterior::new`](crate::games::GameExterior::new).
    pub fn new(
        gens: Vec<Game>,
        q: Vec<i128>,
        b: BTreeMap<(usize, usize), i128>,
    ) -> Result<GameClifford, GameCliffordError> {
        GameClifford::with_relation_search(gens, 3, q, b)
    }

    /// The free integer Clifford algebra on the chosen generators, with no
    /// game-group relations imposed. This is useful as an ambient object, but it
    /// does not check torsion or duplicate-generator constraints.
    pub fn free(
        gens: Vec<Game>,
        q: Vec<i128>,
        b: BTreeMap<(usize, usize), i128>,
    ) -> Result<GameClifford, GameCliffordError> {
        GameClifford::with_quadratic_data(gens, vec![], q, b)
    }

    /// Build from all bounded discovered relations `Σ c_i g_i = 0`, then verify
    /// those relations against the supplied quadratic data.
    pub fn with_relation_search(
        gens: Vec<Game>,
        bound: i128,
        q: Vec<i128>,
        b: BTreeMap<(usize, usize), i128>,
    ) -> Result<GameClifford, GameCliffordError> {
        let (relations, complete, candidate_count) = discover_relations(&gens, bound);
        let relation_certificate =
            relation_search_certificate(&gens, bound, complete, candidate_count, &relations, true);
        let mut out = GameClifford::with_quadratic_data(gens, relations, q, b)?;
        out.relation_search_complete = complete;
        out.relation_certificate = relation_certificate;
        Ok(out)
    }

    /// Build from explicit game-group relations and hand-supplied integer
    /// quadratic data. Every relation is checked both in the game group and in
    /// the Clifford data:
    ///
    /// * `Σ c_i g_i = 0` in the game group;
    /// * `Q(Σ c_i e_i) = 0`;
    /// * the polar pairing of `Σ c_i e_i` with every basis generator is zero.
    pub fn with_quadratic_data(
        gens: Vec<Game>,
        relations: Vec<GameRelation>,
        q: Vec<i128>,
        b: BTreeMap<(usize, usize), i128>,
    ) -> Result<GameClifford, GameCliffordError> {
        let n = gens.len();
        validate_quadratic_shape(n, &q, &b)?;
        for (relation_index, rel) in relations.iter().enumerate() {
            validate_game_relation(relation_index, &gens, rel)?;
            validate_quadratic_relation(relation_index, rel, &q, &b)?;
        }
        let relation_certificate =
            relation_search_certificate(&gens, 0, true, None, &relations, false);
        let metric = Metric::new(
            q.into_iter().map(Integer).collect(),
            b.into_iter().map(|(key, value)| (key, Integer(value))),
        );
        Ok(GameClifford {
            alg: CliffordAlgebra::new(n, metric),
            gens,
            relation_certificate,
            relations,
            relation_search_complete: true,
        })
    }

    /// The underlying free integer Clifford algebra before quotienting by
    /// game-group relations.
    pub fn algebra(&self) -> &CliffordAlgebra<Integer> {
        &self.alg
    }

    pub fn relations(&self) -> &[GameRelation] {
        &self.relations
    }

    pub fn relation_search_complete(&self) -> bool {
        self.relation_search_complete
    }

    pub fn relation_search_certificate(&self) -> &RelationSearchCertificate {
        &self.relation_certificate
    }

    /// The grade-1 generator `e_i` (corresponding to the game `g_i`), reduced in
    /// the checked Clifford quotient.
    pub fn generator(&self, i: usize) -> Multivector<Integer> {
        self.reduce(&self.alg.gen(i))
    }

    /// The game `g_i` a generator stands for.
    pub fn game(&self, i: usize) -> &Game {
        &self.gens[i]
    }

    /// The module map from grade-1 elements to the game group. Panics if the
    /// reduced multivector is not purely grade 1.
    pub fn value_of_grade1(&self, mv: &Multivector<Integer>) -> Game {
        let mut acc = Game::zero();
        let mv = self.reduce(mv);
        for (&blade, coeff) in &mv.terms {
            assert_eq!(
                blade.count_ones(),
                1,
                "value_of_grade1 expects a grade-1 element"
            );
            let i = blade.trailing_zeros() as usize;
            acc = acc.add(&self.gens[i].times_int(coeff.0));
        }
        acc
    }

    pub fn add(&self, a: &Multivector<Integer>, b: &Multivector<Integer>) -> Multivector<Integer> {
        self.reduce(&self.alg.add(a, b))
    }

    pub fn scalar_mul(&self, s: i128, a: &Multivector<Integer>) -> Multivector<Integer> {
        self.reduce(&self.alg.scalar_mul(&Integer(s), a))
    }

    /// Quotient-aware Clifford product.
    pub fn mul(&self, a: &Multivector<Integer>, b: &Multivector<Integer>) -> Multivector<Integer> {
        self.reduce(&self.alg.mul(a, b))
    }

    /// Metric-independent exterior product, reduced in the same checked quotient.
    pub fn wedge(
        &self,
        a: &Multivector<Integer>,
        b: &Multivector<Integer>,
    ) -> Multivector<Integer> {
        self.reduce(&self.alg.wedge(a, b))
    }

    pub fn is_zero(&self, mv: &Multivector<Integer>) -> bool {
        self.reduce(mv).is_zero()
    }

    /// Reduce a free Clifford multivector by the two-sided ideal generated by the
    /// stored grade-1 game relations. Constructor compatibility checks ensure
    /// these relation vectors are null and polar-radical, so this quotient is the
    /// intended integer Clifford deformation of the game subgroup.
    pub fn reduce(&self, mv: &Multivector<Integer>) -> Multivector<Integer> {
        reduce_by_clifford_relation_ideal(&self.alg, self.gens.len(), &self.relations, mv)
    }
}

fn validate_quadratic_shape(
    n: usize,
    q: &[i128],
    b: &BTreeMap<(usize, usize), i128>,
) -> Result<(), GameCliffordError> {
    if q.len() != n {
        return Err(GameCliffordError::QuadraticLength {
            expected: n,
            got: q.len(),
        });
    }
    for &(i, j) in b.keys() {
        if i >= j || j >= n {
            return Err(GameCliffordError::BilinearKeyInvalid { i, j, dim: n });
        }
    }
    Ok(())
}

fn validate_game_relation(
    relation_index: usize,
    gens: &[Game],
    rel: &GameRelation,
) -> Result<(), GameCliffordError> {
    if rel.coeffs.len() != gens.len() {
        return Err(GameCliffordError::RelationLength {
            relation_index,
            expected: gens.len(),
            got: rel.coeffs.len(),
        });
    }
    let value = eval_relation(gens, &rel.coeffs);
    if !value.eq(&Game::zero()) {
        return Err(GameCliffordError::RelationNotZero {
            relation_index,
            value_key: value.canonical_string(),
        });
    }
    Ok(())
}

fn validate_quadratic_relation(
    relation_index: usize,
    rel: &GameRelation,
    q: &[i128],
    b: &BTreeMap<(usize, usize), i128>,
) -> Result<(), GameCliffordError> {
    for j in 0..q.len() {
        let value = relation_polar_value(relation_index, &rel.coeffs, q, b, j)?;
        if value != 0 {
            return Err(GameCliffordError::RelationPolarNonzero {
                relation_index,
                generator: j,
                value,
            });
        }
    }
    let value = relation_quadratic_value(relation_index, &rel.coeffs, q, b)?;
    if value != 0 {
        return Err(GameCliffordError::RelationQuadraticNonzero {
            relation_index,
            value,
        });
    }
    Ok(())
}

fn relation_polar_value(
    relation_index: usize,
    coeffs: &[i128],
    q: &[i128],
    b: &BTreeMap<(usize, usize), i128>,
    j: usize,
) -> Result<i128, GameCliffordError> {
    let mut acc = 0i128;
    for (i, &c) in coeffs.iter().enumerate() {
        if c == 0 {
            continue;
        }
        let polar_entry = if i == j {
            checked_mul_i128(relation_index, q[i], 2, "diagonal polar entry")?
        } else {
            let key = if i < j { (i, j) } else { (j, i) };
            *b.get(&key).unwrap_or(&0)
        };
        let term = checked_mul_i128(relation_index, c, polar_entry, "polar term")?;
        acc = checked_add_i128(relation_index, acc, term, "polar sum")?;
    }
    Ok(acc)
}

fn relation_quadratic_value(
    relation_index: usize,
    coeffs: &[i128],
    q: &[i128],
    b: &BTreeMap<(usize, usize), i128>,
) -> Result<i128, GameCliffordError> {
    let mut acc = 0i128;
    for (i, &c) in coeffs.iter().enumerate() {
        if c == 0 || q[i] == 0 {
            continue;
        }
        let square = checked_mul_i128(relation_index, c, c, "quadratic square")?;
        let term = checked_mul_i128(relation_index, square, q[i], "diagonal quadratic term")?;
        acc = checked_add_i128(relation_index, acc, term, "quadratic sum")?;
    }
    for i in 0..coeffs.len() {
        for j in i + 1..coeffs.len() {
            let bij = *b.get(&(i, j)).unwrap_or(&0);
            if coeffs[i] == 0 || coeffs[j] == 0 || bij == 0 {
                continue;
            }
            let coeff_product =
                checked_mul_i128(relation_index, coeffs[i], coeffs[j], "cross coefficient")?;
            let term =
                checked_mul_i128(relation_index, coeff_product, bij, "cross quadratic term")?;
            acc = checked_add_i128(relation_index, acc, term, "quadratic sum")?;
        }
    }
    Ok(acc)
}

fn checked_add_i128(
    relation_index: usize,
    a: i128,
    b: i128,
    context: &'static str,
) -> Result<i128, GameCliffordError> {
    a.checked_add(b)
        .ok_or(GameCliffordError::ArithmeticOverflow {
            relation_index,
            context,
        })
}

fn checked_mul_i128(
    relation_index: usize,
    a: i128,
    b: i128,
    context: &'static str,
) -> Result<i128, GameCliffordError> {
    a.checked_mul(b)
        .ok_or(GameCliffordError::ArithmeticOverflow {
            relation_index,
            context,
        })
}

fn reduce_by_clifford_relation_ideal(
    alg: &CliffordAlgebra<Integer>,
    dim: usize,
    relations: &[GameRelation],
    mv: &Multivector<Integer>,
) -> Multivector<Integer> {
    if relations.is_empty() || mv.is_zero() {
        return mv.clone();
    }
    let basis = all_blade_masks(dim);
    let index: BTreeMap<u128, usize> = basis.iter().enumerate().map(|(i, &m)| (m, i)).collect();
    let mut v = vec![0i128; basis.len()];
    for (&blade, coeff) in &mv.terms {
        if let Some(&i) = index.get(&blade) {
            v[i] += coeff.0;
        }
    }
    let rows = relation_rows_for_clifford_ideal(alg, relations, &basis, &index);
    reduce_integer_vector(&mut v, rows);
    let terms = basis
        .into_iter()
        .zip(v)
        .filter(|&(_, coeff)| coeff != 0)
        .map(|(blade, coeff)| (blade, Integer(coeff)))
        .collect();
    Multivector { terms }
}

fn relation_rows_for_clifford_ideal(
    alg: &CliffordAlgebra<Integer>,
    relations: &[GameRelation],
    basis: &[u128],
    index: &BTreeMap<u128, usize>,
) -> Vec<Vec<i128>> {
    let mut rows = Vec::new();
    for rel in relations {
        let rel_mv = relation_multivector(rel);
        for &mask in basis {
            let blade = alg.blade(&bits(mask));
            push_clifford_relation_row(alg.mul(&rel_mv, &blade), index, &mut rows);
            push_clifford_relation_row(alg.mul(&blade, &rel_mv), index, &mut rows);
        }
    }
    rows
}

fn push_clifford_relation_row(
    mv: Multivector<Integer>,
    index: &BTreeMap<u128, usize>,
    rows: &mut Vec<Vec<i128>>,
) {
    let mut row = vec![0i128; index.len()];
    for (blade, coeff) in mv.terms {
        if let Some(&i) = index.get(&blade) {
            row[i] += coeff.0;
        }
    }
    if row.iter().any(|&x| x != 0) {
        rows.push(row);
    }
}

fn all_blade_masks(n: usize) -> Vec<u128> {
    let mut out = Vec::new();
    for grade in 0..=n {
        out.extend(grade_masks(n, grade));
    }
    out
}
