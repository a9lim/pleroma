//! The **exterior algebra of the game group**: `Λ` over `ℤ` on a chosen tuple of
//! games. This is the Clifford-adjacent structure that lives on *all* of
//! game-world — not just the field-like numbers — because the partizan games form
//! an abelian group (a `ℤ`-module), and the Grassmann algebra is the exterior
//! algebra of that module.
//!
//! [`GameExterior`] wraps the free Grassmann engine ([`CliffordAlgebra`] over
//! [`Integer`] with the all-zero metric) and quotients it by the integer
//! relations that actually hold among the chosen generators (e.g. `2⋆ = 0`), so a
//! relation propagates through the exterior ideal: `2⋆ = 0 ⟹ 2(⋆∧↑) = 0`. The
//! relations are either supplied explicitly or discovered by a small bounded
//! search; the integer-lattice reduction that imposes them is the row machinery
//! at the bottom of this file.
//!
//! Generators may be non-numbers (`⋆`, `↑`, switches) — exactly where the
//! Clifford/scalar story cannot go — which is the point: the [`Game`] group is not
//! a ring, but it *is* a `ℤ`-module, and that is enough for `Λ`.
//!
//! [`GameClifford`] is the checked deformation wrapper: the caller supplies an
//! integer-valued Clifford metric on the same generator tuple, and the constructor
//! verifies that every game relation is null and polar-radical for that data before
//! imposing the Clifford ideal. This is a quotient-compatible engineering surface,
//! not a claim that the metric was game-native.

use super::Game;
use crate::clifford::{bits, CliffordAlgebra, Metric, Multivector};
use crate::linalg::integer::reduce_integer_vector;
use crate::scalar::Integer;
use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

const DEFAULT_RELATION_BOUND: i128 = 3;
const MAX_AUTO_RELATION_CANDIDATES: usize = 100;
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

/// The exterior algebra generated by a chosen tuple of games, quotienting the
/// free Grassmann algebra by known integer relations among those games.
///
/// The raw [`algebra`](Self::algebra) is still the free Grassmann engine. The
/// quotient-aware operations ([`reduce`](Self::reduce), [`wedge`](Self::wedge),
/// [`add`](Self::add), [`is_zero`](Self::is_zero)) impose the exterior ideal
/// generated by the stored grade-1 relations, so a relation such as `2⋆ = 0`
/// propagates to `2(⋆∧↑) = 0`.
#[derive(Clone)]
pub struct GameExterior {
    alg: CliffordAlgebra<Integer>,
    gens: Vec<Game>,
    relations: Vec<GameRelation>,
    relation_search_complete: bool,
    relation_certificate: RelationSearchCertificate,
}

impl GameExterior {
    pub fn new(gens: Vec<Game>) -> GameExterior {
        GameExterior::with_relation_search(gens, DEFAULT_RELATION_BOUND)
    }

    /// The free Grassmann algebra on the chosen generators, with no game-group
    /// relations imposed. Useful as the ambient object when explicit quotienting
    /// is not desired.
    pub fn free(gens: Vec<Game>) -> GameExterior {
        GameExterior::with_relations(gens, vec![])
    }

    /// Build the quotient using all bounded discovered relations `Σ c_i g_i = 0`
    /// with coefficients in `[-bound, bound]`, when that finite search is small
    /// enough to run exhaustively. If the coefficient box is too large, automatic
    /// discovery falls back to singleton torsion and
    /// [`relation_search_complete`](Self::relation_search_complete) reports
    /// `false`; use [`with_relations`](Self::with_relations) for known larger
    /// cross-generator relations.
    pub fn with_relation_search(gens: Vec<Game>, bound: i128) -> GameExterior {
        let (relations, complete, candidate_count) = discover_relations(&gens, bound);
        let relation_certificate =
            relation_search_certificate(&gens, bound, complete, candidate_count, &relations, true);
        let mut ext = GameExterior::with_relations(gens, relations);
        ext.relation_search_complete = complete;
        ext.relation_certificate = relation_certificate;
        ext
    }

    /// Build the quotient from explicit integer relations among the supplied
    /// generators. Each relation is verified against the game group before it is
    /// accepted.
    pub fn with_relations(gens: Vec<Game>, relations: Vec<GameRelation>) -> GameExterior {
        let n = gens.len();
        for rel in &relations {
            assert_eq!(
                rel.coeffs.len(),
                n,
                "game relation length must match generator count"
            );
            assert!(
                eval_relation(&gens, &rel.coeffs).eq(&Game::zero()),
                "declared game relation does not evaluate to zero"
            );
        }
        let relation_certificate =
            relation_search_certificate(&gens, 0, true, None, &relations, false);
        GameExterior {
            alg: CliffordAlgebra::new(n, Metric::grassmann(n)),
            gens,
            relation_certificate,
            relations,
            relation_search_complete: true,
        }
    }

    /// The underlying free Grassmann algebra. Use the quotient-aware methods on
    /// `GameExterior` when game-group relations should be imposed.
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

    /// The grade-1 generator `e_i` (corresponding to the game `g_i`).
    pub fn generator(&self, i: usize) -> Multivector<Integer> {
        self.reduce(&self.alg.gen(i))
    }

    /// The game `g_i` a generator stands for.
    pub fn game(&self, i: usize) -> &Game {
        &self.gens[i]
    }

    /// The module map `Λ¹ → (game group)`: send a grade-1 element `Σ c_i e_i` to
    /// the game `Σ c_i · g_i`. Defined entirely with the game *group* operations
    /// (sum, negation, integer multiple) and no game product — so it is valid for
    /// non-number generators. Panics if `mv` is not purely grade 1.
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

    /// Reduce a free Grassmann multivector modulo the exterior ideal generated by
    /// the stored game relations.
    pub fn reduce(&self, mv: &Multivector<Integer>) -> Multivector<Integer> {
        if self.relations.is_empty() || mv.is_zero() {
            return mv.clone();
        }
        let mut out = self.alg.zero();
        let mut by_grade: BTreeMap<usize, BTreeMap<u128, i128>> = BTreeMap::new();
        for (&blade, coeff) in &mv.terms {
            by_grade
                .entry(blade.count_ones() as usize)
                .or_default()
                .insert(blade, coeff.0);
        }
        for (grade, terms) in by_grade {
            let reduced = self.reduce_grade(grade, &terms);
            for (blade, coeff) in reduced {
                if coeff != 0 {
                    out.terms.insert(blade, Integer(coeff));
                }
            }
        }
        out
    }

    fn reduce_grade(&self, grade: usize, terms: &BTreeMap<u128, i128>) -> BTreeMap<u128, i128> {
        if grade == 0 {
            return terms.clone();
        }
        let basis = grade_masks(self.gens.len(), grade);
        if basis.is_empty() {
            return BTreeMap::new();
        }
        let index: BTreeMap<u128, usize> = basis.iter().enumerate().map(|(i, &m)| (m, i)).collect();
        let mut v = vec![0i128; basis.len()];
        for (&blade, &coeff) in terms {
            if let Some(&i) = index.get(&blade) {
                v[i] += coeff;
            }
        }
        let rows = self.relation_rows_for_grade(grade, &basis, &index);
        reduce_integer_vector(&mut v, rows);
        basis.into_iter().zip(v).filter(|&(_, c)| c != 0).collect()
    }

    fn relation_rows_for_grade(
        &self,
        grade: usize,
        basis: &[u128],
        index: &BTreeMap<u128, usize>,
    ) -> Vec<Vec<i128>> {
        let mut rows = Vec::new();
        let lower_basis = grade_masks(self.gens.len(), grade - 1);
        for rel in &self.relations {
            let rel_mv = relation_multivector(rel);
            for mask in &lower_basis {
                let blade = self.alg.blade(&bits(*mask));
                let wedged = self.alg.wedge(&rel_mv, &blade);
                let mut row = vec![0i128; basis.len()];
                for (&b, coeff) in &wedged.terms {
                    if let Some(&i) = index.get(&b) {
                        row[i] += coeff.0;
                    }
                }
                if row.iter().any(|&x| x != 0) {
                    rows.push(row);
                }
            }
        }
        rows
    }
}

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
    /// [`GameExterior::new`].
    pub fn new(
        gens: Vec<Game>,
        q: Vec<i128>,
        b: BTreeMap<(usize, usize), i128>,
    ) -> Result<GameClifford, GameCliffordError> {
        GameClifford::with_relation_search(gens, DEFAULT_RELATION_BOUND, q, b)
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

fn relation_multivector(rel: &GameRelation) -> Multivector<Integer> {
    let mut terms = BTreeMap::new();
    for (i, &coeff) in rel.coeffs.iter().enumerate() {
        if coeff != 0 {
            terms.insert(1u128 << i, Integer(coeff));
        }
    }
    Multivector { terms }
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

fn eval_relation(gens: &[Game], coeffs: &[i128]) -> Game {
    let mut acc = Game::zero();
    for (g, &c) in gens.iter().zip(coeffs) {
        acc = acc.add(&g.times_int(c));
    }
    acc
}

fn canonical_relation(mut coeffs: Vec<i128>) -> Option<Vec<i128>> {
    let first = coeffs.iter().position(|&c| c != 0)?;
    if coeffs[first] < 0 {
        for c in &mut coeffs {
            *c = -*c;
        }
    }
    Some(coeffs)
}

fn discover_relations(gens: &[Game], bound: i128) -> (Vec<GameRelation>, bool, Option<usize>) {
    if gens.is_empty() || bound <= 0 {
        return (Vec::new(), true, Some(0));
    }
    let n = gens.len();
    let mut seen = BTreeSet::new();
    let mut out = Vec::new();

    for i in 0..n {
        for c in 1..=bound {
            let mut coeffs = vec![0i128; n];
            coeffs[i] = c;
            if push_relation_if_independent(gens, coeffs, &mut seen, &mut out) {
                break;
            }
        }
    }

    let Some(count) = bounded_relation_candidate_count(n, bound) else {
        return (out, false, None);
    };
    if count > MAX_AUTO_RELATION_CANDIDATES {
        return (out, false, Some(count));
    }

    let mut candidates = Vec::new();
    enumerate_bounded_relations(n, bound, &mut |coeffs| {
        if let Some(key) = canonical_relation(coeffs) {
            candidates.push(key);
        }
    });
    candidates.sort_by_key(|v| (v.iter().map(|c| c.abs()).sum::<i128>(), v.clone()));
    for coeffs in candidates {
        push_relation_if_independent(gens, coeffs, &mut seen, &mut out);
    }
    (out, true, Some(count))
}

fn relation_search_certificate(
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

fn relation_certificates(
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

fn bounded_relation_candidate_count(n: usize, bound: i128) -> Option<usize> {
    let width = usize::try_from(bound.checked_mul(2)?.checked_add(1)?).ok()?;
    let mut count = 1usize;
    for _ in 0..n {
        count = count.checked_mul(width)?;
    }
    count.checked_sub(1)
}

fn enumerate_bounded_relations(n: usize, bound: i128, f: &mut impl FnMut(Vec<i128>)) {
    fn rec(i: usize, n: usize, bound: i128, coeffs: &mut [i128], f: &mut impl FnMut(Vec<i128>)) {
        if i == n {
            if coeffs.iter().any(|&c| c != 0) {
                f(coeffs.to_vec());
            }
            return;
        }
        for c in -bound..=bound {
            coeffs[i] = c;
            rec(i + 1, n, bound, coeffs, f);
        }
    }
    let mut coeffs = vec![0i128; n];
    rec(0, n, bound, &mut coeffs, f);
}

fn push_relation_if_independent(
    gens: &[Game],
    coeffs: Vec<i128>,
    seen: &mut BTreeSet<Vec<i128>>,
    out: &mut Vec<GameRelation>,
) -> bool {
    let Some(key) = canonical_relation(coeffs) else {
        return false;
    };
    if !seen.insert(key.clone()) {
        return false;
    }
    if !eval_relation(gens, &key).eq(&Game::zero()) {
        return false;
    }
    let mut reduced = key.clone();
    let rows: Vec<Vec<i128>> = out.iter().map(|r| r.coeffs.clone()).collect();
    reduce_integer_vector(&mut reduced, rows);
    if reduced.iter().all(|&c| c == 0) {
        return false;
    }
    out.push(GameRelation::new(key));
    true
}

fn grade_masks(n: usize, grade: usize) -> Vec<u128> {
    if grade > n {
        return Vec::new();
    }
    fn rec(n: usize, grade: usize, start: usize, mask: u128, out: &mut Vec<u128>) {
        if grade == 0 {
            out.push(mask);
            return;
        }
        for i in start..=n - grade {
            rec(n, grade - 1, i + 1, mask | (1u128 << i), out);
        }
    }
    let mut out = Vec::new();
    rec(n, grade, 0, 0, &mut out);
    out
}

fn all_blade_masks(n: usize) -> Vec<u128> {
    let mut out = Vec::new();
    for grade in 0..=n {
        out.extend(grade_masks(n, grade));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

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
