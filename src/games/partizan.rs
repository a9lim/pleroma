//! Short partizan combinatorial games, and the exterior algebra of the **game
//! group**.
//!
//! Conway's games `G = { G^L | G^R }` form, under disjunctive sum, a partially
//! ordered abelian group — but *not a ring* (the product is only a congruence on
//! the numbers). That is the obstruction the whole project lives around: a
//! Clifford algebra needs a commutative scalar *ring*, so it only reaches the
//! field-like cores (nimbers, surreals, surcomplex).
//!
//! The **exterior algebra**, by contrast, needs only a commutative ring of
//! *coefficients* (here ℤ) and a *module* of generators — and the game group is a
//! ℤ-module. So `Λ(game group)` is well defined on ALL of game-world, the one
//! Clifford-adjacent structure that does not require the game product to exist.
//! Concretely: pick generator games `g_1, …, g_n`, build the free Grassmann
//! algebra over the formal basis `e_i`, detect or accept integer relations among
//! the `g_i`, and quotient by the exterior ideal those linear relations generate.
//! The grade-1 map `e_i ↦ g_i` uses only disjunctive sum, so it works for
//! generators that are not even numbers (`⋆`, `↑`), where Conway multiplication
//! is undefined.
//!
//! This module ships a small short-game engine (sum, negation, the recursive
//! order, birthday, number test) good enough to construct example games and prove
//! the generators are genuine games, plus the `GameExterior` bridge.

use crate::clifford::{bits, CliffordAlgebra, Metric, Multivector};
use crate::scalar::{Integer, Scalar, Surreal};
use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

/// A short partizan game `{ left | right }`. Reference-counted (atomically, so the
/// PyO3 wrapper is `Send + Sync`) — options are shared cheaply across the
/// recursive sum/negation.
#[derive(Clone)]
pub struct Game(Arc<GameData>);

struct GameData {
    left: Vec<Game>,
    right: Vec<Game>,
}

impl Game {
    pub fn new(left: Vec<Game>, right: Vec<Game>) -> Game {
        Game(Arc::new(GameData { left, right }))
    }

    pub fn left(&self) -> &[Game] {
        &self.0.left
    }
    pub fn right(&self) -> &[Game] {
        &self.0.right
    }

    /// `0 = { | }` — the empty game (second player wins).
    pub fn zero() -> Game {
        Game::new(vec![], vec![])
    }

    /// `⋆ = { 0 | 0 }` — a single Nim-heap of size 1. Fuzzy with 0; *not* a number.
    pub fn star() -> Game {
        let z = Game::zero();
        Game::new(vec![z.clone()], vec![z])
    }

    /// The integer game `n`: `{ n−1 | }` for n>0, `{ | n+1 }` for n<0, `0` for 0.
    pub fn integer(n: i128) -> Game {
        if n == 0 {
            Game::zero()
        } else if n > 0 {
            Game::new(vec![Game::integer(n - 1)], vec![])
        } else {
            Game::new(vec![], vec![Game::integer(n + 1)])
        }
    }

    /// `↑ = { 0 | ⋆ }` — "up", a positive infinitesimal; `0 < ↑` but `↑ < x` for
    /// every positive number x. A non-number.
    pub fn up() -> Game {
        Game::new(vec![Game::zero()], vec![Game::star()])
    }

    /// The switch `{ a | b }` (e.g. `{1 | -1}` is `±1`). A non-number when a ≥ b.
    pub fn switch(a: i128, b: i128) -> Game {
        Game::new(vec![Game::integer(a)], vec![Game::integer(b)])
    }

    /// Negation `−G = { −G^R | −G^L }` (the additive inverse in the game group).
    pub fn neg(&self) -> Game {
        Game::new(
            self.right().iter().map(|g| g.neg()).collect(),
            self.left().iter().map(|g| g.neg()).collect(),
        )
    }

    /// Disjunctive sum `G + H` — the group operation.
    pub fn add(&self, other: &Game) -> Game {
        let mut left = Vec::new();
        for gl in self.left() {
            left.push(gl.add(other));
        }
        for hl in other.left() {
            left.push(self.add(hl));
        }
        let mut right = Vec::new();
        for gr in self.right() {
            right.push(gr.add(other));
        }
        for hr in other.right() {
            right.push(self.add(hr));
        }
        Game::new(left, right)
    }

    /// The order: `G ≤ H ⟺ (∄ G^L ≥ H) ∧ (∄ H^R ≤ G)`. Recurses on options
    /// (strictly simpler games), so it terminates.
    pub fn le(&self, other: &Game) -> bool {
        self.left().iter().all(|gl| !other.le(gl)) && other.right().iter().all(|hr| !hr.le(self))
    }

    /// Value equality: `G = H ⟺ G ≤ H ≤ G`.
    pub fn eq(&self, other: &Game) -> bool {
        self.le(other) && other.le(self)
    }

    /// Confused/incomparable: `G ‖ H` (neither `≤` holds) — the hallmark of a
    /// non-number relative to its options.
    pub fn fuzzy(&self, other: &Game) -> bool {
        !self.le(other) && !other.le(self)
    }

    /// The birthday (formation day): `0` for `{|}`, else `1 + max` over options.
    pub fn birthday(&self) -> u128 {
        self.left()
            .iter()
            .chain(self.right())
            .map(|g| g.birthday())
            .max()
            .map_or(0, |m| m + 1)
    }

    /// The integer multiple `n · G` in the game group (repeated sum / negation).
    pub fn times_int(&self, n: i128) -> Game {
        if n == 0 {
            Game::zero()
        } else if n > 0 {
            let mut acc = self.clone();
            for _ in 1..n {
                acc = acc.add(self);
            }
            acc
        } else {
            self.neg().times_int(-n)
        }
    }

    /// A readable structural form: `0` for `{|}`, else `{L|R}` recursively.
    pub fn display(&self) -> String {
        if self.left().is_empty() && self.right().is_empty() {
            return "0".to_string();
        }
        let l: Vec<String> = self.left().iter().map(|g| g.display()).collect();
        let r: Vec<String> = self.right().iter().map(|g| g.display()).collect();
        format!("{{{}|{}}}", l.join(","), r.join(","))
    }

    /// Whether `G` is a (surreal) *number*: all options are numbers and every left
    /// option is strictly below every right option. Numbers are exactly the games
    /// the Conway product (and hence the Clifford story) can reach.
    pub fn is_number(&self) -> bool {
        self.left().iter().all(|g| g.is_number())
            && self.right().iter().all(|g| g.is_number())
            && self
                .left()
                .iter()
                .all(|gl| self.right().iter().all(|gr| gl.le(gr) && !gr.le(gl)))
    }

    // ---- Canonical form (Conway's simplicity theorem) ----

    /// The **canonical form**: the unique simplest game equal in value to `self`.
    /// Options are first put in canonical form, then dominated options are
    /// removed and reversible options bypassed, repeatedly, until stable. Two
    /// short games are equal iff their canonical forms are
    /// [structurally identical](Self::structural_eq), so this is the normal form
    /// that makes equality a syntactic check and `birthday` the true (least)
    /// formation day.
    pub fn canonical(&self) -> Game {
        let left: Vec<Game> = self.left().iter().map(Game::canonical).collect();
        let right: Vec<Game> = self.right().iter().map(Game::canonical).collect();
        let mut cur = Game::new(left, right);
        loop {
            let (bypassed, bypassed_any) = cur.bypass_reversible_once();
            let reduced = bypassed.remove_dominated();
            let removed_any = reduced.left().len() != bypassed.left().len()
                || reduced.right().len() != bypassed.right().len();
            cur = reduced;
            if !bypassed_any && !removed_any {
                return cur;
            }
        }
    }

    /// One pass of reversibility bypass: a Left option `G^L` with some Right
    /// option `G^LR ≤ G` is replaced by all Left options of that `G^LR`
    /// (symmetrically on the Right). Returns the new game and whether anything
    /// was bypassed.
    fn bypass_reversible_once(&self) -> (Game, bool) {
        let mut changed = false;
        let mut new_left = Vec::new();
        for l in self.left() {
            if let Some(lr) = l.right().iter().find(|lr| lr.le(self)) {
                changed = true;
                new_left.extend(lr.left().iter().cloned());
            } else {
                new_left.push(l.clone());
            }
        }
        let mut new_right = Vec::new();
        for r in self.right() {
            if let Some(rl) = r.left().iter().find(|rl| self.le(rl)) {
                changed = true;
                new_right.extend(rl.right().iter().cloned());
            } else {
                new_right.push(r.clone());
            }
        }
        (Game::new(new_left, new_right), changed)
    }

    /// Drop dominated options: keep only the order-maximal Left options and the
    /// order-minimal Right options (one representative per equal value).
    fn remove_dominated(&self) -> Game {
        Game::new(maximal_games(self.left()), minimal_games(self.right()))
    }

    /// An order-independent string `{L|R}` of the game *as given* (options sorted
    /// recursively) — a structural fingerprint, **not** reduced to canonical
    /// form. Use [`canonical_string`](Self::canonical_string) for a value key.
    pub fn structural_string(&self) -> String {
        let mut l: Vec<String> = self.left().iter().map(Game::structural_string).collect();
        let mut r: Vec<String> = self.right().iter().map(Game::structural_string).collect();
        l.sort();
        r.sort();
        format!("{{{}|{}}}", l.join(","), r.join(","))
    }

    /// The canonical-form string — a true fingerprint of the game's *value*: it
    /// canonicalizes first, then renders order-independently. Two games are equal
    /// in value iff their `canonical_string`s match, so this is a hashable key for
    /// game values.
    pub fn canonical_string(&self) -> String {
        self.canonical().structural_string()
    }

    /// Structural identity up to option ordering, of the games *as given* (no
    /// canonicalization). Most useful as `a.canonical().structural_eq(&b.canonical())`;
    /// for value equality prefer [`eq`](Self::eq) or matching
    /// [`canonical_string`](Self::canonical_string)s.
    pub fn structural_eq(&self, other: &Game) -> bool {
        self.structural_string() == other.structural_string()
    }

    /// Whether `self` is already in canonical form.
    pub fn is_canonical(&self) -> bool {
        self.structural_eq(&self.canonical())
    }

    // ---- The game ↔ surreal bridge (numbers only) ----

    /// The surreal value of a number-valued game, by the simplicity theorem:
    /// `value({G^L | G^R})` is the simplest surreal strictly between the largest
    /// Left value and the smallest Right value. `None` if `self` is not a number
    /// (`⋆`, `↑`, switches, …) or its value is not dyadic. Inverse of
    /// [`from_surreal`](Self::from_surreal) on dyadics.
    pub fn number_value(&self) -> Option<Surreal> {
        if !self.is_number() {
            return None;
        }
        let lvals: Vec<Surreal> = self
            .left()
            .iter()
            .map(Game::number_value)
            .collect::<Option<_>>()?;
        let rvals: Vec<Surreal> = self
            .right()
            .iter()
            .map(Game::number_value)
            .collect::<Option<_>>()?;
        let lmax = lvals
            .into_iter()
            .reduce(|a, b| if a.cmp(&b) == Ordering::Less { b } else { a });
        let rmin = rvals
            .into_iter()
            .reduce(|a, b| if a.cmp(&b) == Ordering::Greater { b } else { a });
        match (lmax, rmin) {
            (None, None) => Some(Surreal::zero()),
            (Some(l), None) => l.simplest_above(),
            (None, Some(r)) => r.simplest_below(),
            (Some(l), Some(r)) => Surreal::simplest_between(&l, &r),
        }
    }

    /// The canonical game of a dyadic-rational surreal — the `{L|R}` form Conway's
    /// construction gives that number. `None` if `s` is not dyadic (infinite,
    /// infinitesimal, or a non-dyadic rational, none of which is a short game).
    /// Inverse of [`number_value`](Self::number_value).
    pub fn from_surreal(s: &Surreal) -> Option<Game> {
        let (num, k) = s.as_dyadic()?;
        Some(game_of_dyadic(num, k))
    }
}

/// Keep only the order-maximal games (Left options of a canonical form): drop any
/// option dominated by — or equal to — a kept one.
fn maximal_games(opts: &[Game]) -> Vec<Game> {
    let mut kept: Vec<Game> = Vec::new();
    for cand in opts {
        if kept.iter().any(|k| cand.le(k)) {
            continue; // dominated by (or equal to) a kept option
        }
        kept.retain(|k| !k.le(cand)); // drop kept options strictly below cand
        kept.push(cand.clone());
    }
    kept
}

/// Keep only the order-minimal games (Right options of a canonical form).
fn minimal_games(opts: &[Game]) -> Vec<Game> {
    let mut kept: Vec<Game> = Vec::new();
    for cand in opts {
        if kept.iter().any(|k| k.le(cand)) {
            continue; // dominated by (or equal to) a kept option
        }
        kept.retain(|k| !cand.le(k)); // drop kept options strictly above cand
        kept.push(cand.clone());
    }
    kept
}

/// Strip factors of two from a dyadic `num / 2^k` to put it in lowest terms.
fn reduce_dyadic_pair(mut num: i128, mut k: u32) -> (i128, u32) {
    while k > 0 && num % 2 == 0 {
        num /= 2;
        k -= 1;
    }
    (num, k)
}

/// The canonical game of the dyadic `num / 2^k`: an integer for `k = 0`, else
/// `{ (num-1)/2^k | (num+1)/2^k }` with the options reduced to lowest terms.
fn game_of_dyadic(num: i128, k: u32) -> Game {
    if k == 0 {
        return Game::integer(num);
    }
    let (ln, lk) = reduce_dyadic_pair(num - 1, k);
    let (rn, rk) = reduce_dyadic_pair(num + 1, k);
    Game::new(vec![game_of_dyadic(ln, lk)], vec![game_of_dyadic(rn, rk)])
}

const DEFAULT_RELATION_BOUND: i128 = 3;
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

    /// Build the quotient using small discovered relations `Σ c_i g_i = 0`.
    /// Automatic discovery always checks singleton torsion and checks pair/full
    /// coefficient searches only for two-generator presentations; larger
    /// presentations should use [`with_relations`](Self::with_relations) for
    /// known cross-generator relations.
    pub fn with_relation_search(gens: Vec<Game>, bound: i128) -> GameExterior {
        let relations = discover_relations(&gens, bound);
        GameExterior::with_relations(gens, relations)
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
        GameExterior {
            alg: CliffordAlgebra::new(n, Metric::grassmann(n)),
            gens,
            relations,
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

fn relation_multivector(rel: &GameRelation) -> Multivector<Integer> {
    let mut terms = BTreeMap::new();
    for (i, &coeff) in rel.coeffs.iter().enumerate() {
        if coeff != 0 {
            terms.insert(1u128 << i, Integer(coeff));
        }
    }
    Multivector { terms }
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

fn discover_relations(gens: &[Game], bound: i128) -> Vec<GameRelation> {
    if gens.is_empty() || bound <= 0 {
        return Vec::new();
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

    if n > 2 {
        return out;
    }

    let mut candidates = Vec::new();
    for i in 0..n {
        for j in (i + 1)..n {
            for a in -bound..=bound {
                for b in -bound..=bound {
                    if a == 0 && b == 0 {
                        continue;
                    }
                    let mut coeffs = vec![0i128; n];
                    coeffs[i] = a;
                    coeffs[j] = b;
                    let Some(key) = canonical_relation(coeffs) else {
                        continue;
                    };
                    candidates.push(key);
                }
            }
        }
    }
    candidates.sort_by_key(|v| (v.iter().map(|c| c.abs()).sum::<i128>(), v.clone()));
    for coeffs in candidates {
        push_relation_if_independent(gens, coeffs, &mut seen, &mut out);
    }
    out
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

fn leading(row: &[i128]) -> Option<usize> {
    row.iter().position(|&x| x != 0)
}

fn row_is_zero(row: &[i128]) -> bool {
    row.iter().all(|&x| x == 0)
}

fn checked_abs(x: i128) -> i128 {
    x.checked_abs()
        .expect("integer relation coefficient magnitude exceeds i128")
}

fn negate_row(row: &mut [i128]) {
    for x in row {
        *x = x
            .checked_neg()
            .expect("integer relation coefficient magnitude exceeds i128");
    }
}

fn sub_row_multiple(target: &mut [i128], source: &[i128], q: i128) {
    for (t, &s) in target.iter_mut().zip(source) {
        let delta = q
            .checked_mul(s)
            .expect("integer relation row operation exceeds i128");
        *t = t
            .checked_sub(delta)
            .expect("integer relation row operation exceeds i128");
    }
}

/// Row Hermite normal form for an integer row lattice.
///
/// The returned rows generate exactly the same submodule as the input rows, have
/// increasing leading columns, positive pivots, zeros below each pivot, and
/// residues above pivots reduced modulo the pivot. This gives
/// [`reduce_integer_vector`] a canonical quotient representative for
/// `Z^n / <rows>`.
fn normalize_relation_rows(mut rows: Vec<Vec<i128>>) -> Vec<Vec<i128>> {
    let width = rows.first().map_or(0, Vec::len);
    assert!(
        rows.iter().all(|r| r.len() == width),
        "integer relation rows must have equal width"
    );
    rows.retain(|r| !row_is_zero(r));
    let mut rank = 0usize;
    for col in 0..width {
        let Some(pivot) = (rank..rows.len()).find(|&r| rows[r][col] != 0) else {
            continue;
        };
        rows.swap(rank, pivot);
        if rows[rank][col] < 0 {
            negate_row(&mut rows[rank]);
        }

        loop {
            let Some(r) = ((rank + 1)..rows.len()).find(|&r| rows[r][col] != 0) else {
                break;
            };
            let pivot_val = rows[rank][col];
            let q = rows[r][col].div_euclid(pivot_val);
            let source = rows[rank].clone();
            sub_row_multiple(&mut rows[r], &source, q);
            if rows[r][col] != 0 && checked_abs(rows[r][col]) < checked_abs(rows[rank][col]) {
                rows.swap(rank, r);
                if rows[rank][col] < 0 {
                    negate_row(&mut rows[rank]);
                }
            }
        }

        if rows[rank][col] < 0 {
            negate_row(&mut rows[rank]);
        }
        let pivot_val = rows[rank][col];
        let source = rows[rank].clone();
        for r in 0..rows.len() {
            if r == rank || rows[r][col] == 0 {
                continue;
            }
            let q = rows[r][col].div_euclid(pivot_val);
            sub_row_multiple(&mut rows[r], &source, q);
        }
        rank += 1;
    }
    rows.retain(|r| !row_is_zero(r));
    rows.sort_by_key(|r| leading(r).unwrap_or(usize::MAX));
    rows
}

fn reduce_integer_vector(v: &mut [i128], rows: Vec<Vec<i128>>) {
    for row in normalize_relation_rows(rows) {
        let Some(lead) = leading(&row) else {
            continue;
        };
        let pivot = row[lead];
        debug_assert!(pivot > 0);
        let q = v[lead].div_euclid(pivot);
        if q != 0 {
            for i in 0..v.len() {
                v[i] -= q * row[i];
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn game_group_basics() {
        // 1 + (−1) = 0; the integers embed and add correctly.
        assert!(Game::integer(1).add(&Game::integer(-1)).eq(&Game::zero()));
        assert!(Game::integer(2)
            .add(&Game::integer(3))
            .eq(&Game::integer(5)));
        assert!(Game::integer(1).le(&Game::integer(2)));
        // ⋆ is fuzzy with 0 (a non-number), and ⋆ + ⋆ = 0 (it is 2-torsion).
        assert!(Game::star().fuzzy(&Game::zero()));
        assert!(!Game::star().is_number());
        assert!(Game::star().add(&Game::star()).eq(&Game::zero()));
        // ↑ is a positive infinitesimal: 0 < ↑, but ↑ is not a number.
        assert!(Game::zero().le(&Game::up()) && !Game::up().le(&Game::zero()));
        assert!(!Game::up().is_number());
        // birthdays
        assert_eq!(Game::zero().birthday(), 0);
        assert_eq!(Game::star().birthday(), 1);
        assert_eq!(Game::integer(3).birthday(), 3);
    }

    #[test]
    fn canonical_removes_dominated_options() {
        // {0, −1 | } : the Left option −1 is dominated by 0, so it drops out and
        // the game collapses to {0 | } = 1.
        let g = Game::new(vec![Game::integer(0), Game::integer(-1)], vec![]);
        assert!(g.canonical().structural_eq(&Game::integer(1)));
        // {0 | 2} is the number 1, whose canonical form is {0 | } (the 2 reverses
        // out through its Left option 1 ≥ G).
        let g = Game::new(vec![Game::integer(0)], vec![Game::integer(2)]);
        assert!(g.canonical().structural_eq(&Game::integer(1)));
    }

    #[test]
    fn canonical_fixes_the_already_simple_games() {
        // ⋆, ↑ and switches are already canonical.
        for g in [
            Game::star(),
            Game::up(),
            Game::switch(1, -1),
            Game::integer(4),
        ] {
            assert!(g.is_canonical(), "{} should be canonical", g.display());
            assert!(g.canonical().structural_eq(&g));
        }
    }

    #[test]
    fn canonical_of_g_minus_g_is_zero() {
        // G − G = 0 for every game ⇒ its canonical form is the empty game {|}.
        for g in [
            Game::up(),
            Game::switch(3, -2),
            Game::star(),
            Game::integer(2),
        ] {
            let z = g.add(&g.neg());
            assert!(z.eq(&Game::zero()));
            assert!(z.canonical().structural_eq(&Game::zero()));
        }
    }

    #[test]
    fn canonical_is_idempotent_and_value_preserving() {
        let g = Game::new(
            vec![Game::integer(0), Game::integer(-1), Game::switch(2, 0)],
            vec![Game::integer(3)],
        );
        let c = g.canonical();
        assert!(c.eq(&g)); // value preserved
        assert!(c.canonical().structural_eq(&c)); // idempotent
    }

    #[test]
    fn number_value_round_trips_through_games() {
        use crate::scalar::{Rational, Surreal};
        let dy = |n: i128, d: i128| Surreal::from_rational(Rational::new(n, d));
        // surreal → canonical game → surreal is the identity on dyadics
        for s in [dy(0, 1), dy(1, 1), dy(-3, 1), dy(1, 2), dy(3, 4), dy(-5, 8)] {
            let g = Game::from_surreal(&s).unwrap();
            assert_eq!(g.number_value(), Some(s.clone()));
            assert!(g.is_canonical()); // the dyadic game is born canonical
            assert_eq!(s.dyadic_birthday(), Some(g.birthday()));
        }
        // a number game reduces to the canonical game of its value
        let g = Game::new(vec![Game::integer(0)], vec![Game::integer(1)]); // ½
        assert_eq!(g.number_value(), Some(dy(1, 2)));
        // non-numbers have no surreal value
        assert_eq!(Game::star().number_value(), None);
        assert_eq!(Game::up().number_value(), None);
        assert_eq!(Game::switch(1, -1).number_value(), None);
    }

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
    fn integer_relation_reduction_uses_the_full_row_lattice() {
        let rows = vec![vec![2, 0], vec![3, 0]];
        assert_eq!(normalize_relation_rows(rows.clone()), vec![vec![1, 0]]);

        let mut v = vec![5, 7];
        reduce_integer_vector(&mut v, rows);
        assert_eq!(v, vec![0, 7]);
    }

    #[test]
    fn integer_relation_reduction_handles_coupled_relations() {
        let rows = vec![vec![2, 4], vec![6, 10]];
        let mut v = vec![8, 14]; // one copy of each relation
        reduce_integer_vector(&mut v, rows.clone());
        assert_eq!(v, vec![0, 0]);

        let mut shifted = vec![9, 14];
        reduce_integer_vector(&mut shifted, rows);
        assert_ne!(shifted, vec![0, 0]);
    }
}
