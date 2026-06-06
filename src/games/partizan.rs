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
//! Concretely: pick generator games `g_1, …, g_n`, take the Grassmann algebra
//! `Λ` over ℤ on `n` generators (the shipped engine with the null metric), and
//! the module map `Λ¹ → (game group)`, `e_i ↦ g_i`, is built with no game
//! product anywhere — so it works for generators that are not even numbers
//! (`⋆`, `↑`), where Conway multiplication is undefined.
//!
//! This module ships a small short-game engine (sum, negation, the recursive
//! order, birthday, number test) good enough to construct example games and prove
//! the generators are genuine games, plus the `GameExterior` bridge.

use crate::clifford::{CliffordAlgebra, Metric, Multivector};
use crate::scalar::Integer;
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
    pub fn integer(n: i64) -> Game {
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
    pub fn switch(a: i64, b: i64) -> Game {
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
    pub fn birthday(&self) -> u32 {
        self.left()
            .iter()
            .chain(self.right())
            .map(|g| g.birthday())
            .max()
            .map_or(0, |m| m + 1)
    }

    /// The integer multiple `n · G` in the game group (repeated sum / negation).
    pub fn times_int(&self, n: i64) -> Game {
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
}

/// The exterior algebra `Λ` of the game group, on a chosen tuple of generator
/// games. Built on the shipped Grassmann engine (`CliffordAlgebra<Integer>` with
/// the all-null metric), it carries the full wedge structure (antisymmetry,
/// grading) over generators that need *no* product of their own.
pub struct GameExterior {
    alg: CliffordAlgebra<Integer>,
    gens: Vec<Game>,
}

impl GameExterior {
    pub fn new(gens: Vec<Game>) -> GameExterior {
        let n = gens.len();
        GameExterior {
            alg: CliffordAlgebra::new(n, Metric::grassmann(n)),
            gens,
        }
    }

    /// The underlying Grassmann algebra (for wedge, grading, etc.).
    pub fn algebra(&self) -> &CliffordAlgebra<Integer> {
        &self.alg
    }

    /// The grade-1 generator `e_i` (corresponding to the game `g_i`).
    pub fn generator(&self, i: usize) -> Multivector<Integer> {
        self.alg.gen(i)
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
    fn exterior_algebra_lives_on_non_numbers() {
        // Generators that are NOT numbers — exactly where the Clifford/scalar story
        // cannot go — yet the exterior algebra is perfectly well defined on them.
        let ext = GameExterior::new(vec![Game::star(), Game::up(), Game::switch(1, -1)]);
        assert!(!ext.game(0).is_number()); // ⋆
        assert!(!ext.game(1).is_number()); // ↑
        let (e0, e1) = (ext.generator(0), ext.generator(1));
        let alg = ext.algebra();
        // the wedge is antisymmetric and nonzero — genuine grade-2 structure.
        let e01 = alg.wedge(&e0, &e1);
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
}
