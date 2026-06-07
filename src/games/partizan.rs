//! Short partizan combinatorial games.
//!
//! Conway's games `G = { G^L | G^R }` form, under disjunctive sum, a partially
//! ordered abelian group — but *not a ring* (the product is only a congruence on
//! the numbers). That is the obstruction the whole project lives around: a
//! Clifford algebra needs a commutative scalar *ring*, so it only reaches the
//! field-like cores (nimbers, surreals, surcomplex).
//!
//! This module ships the short-game engine — sum, negation, the recursive order,
//! birthday, the number test, and the **canonical form** (dominated/reversible
//! reduction, the value key) — together with the game↔surreal bridge
//! ([`Game::number_value`] / [`Game::from_surreal`], numbers only) and
//! [`NumberGame`], the numbers-only view that carries *transfinite* number games
//! (`ω`, `ε`) by their [`Surreal`] value with no infinite option tree.
//!
//! The exterior algebra of the game group — the one Clifford-adjacent structure
//! that lives on *all* of game-world because it needs only the ℤ-module structure,
//! not the game product — is the sibling module
//! [`game_exterior`](crate::games::game_exterior).

use crate::scalar::{Ordinal, Scalar, Surreal};
use std::cmp::Ordering;
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

    /// The Nim-heap `⋆n = { ⋆0, …, ⋆(n−1) | ⋆0, …, ⋆(n−1) }` (impartial, so the
    /// Left and Right options coincide). `⋆0 = 0`, `⋆1 = ⋆`. Used as the **remote
    /// (far) star** in the atomic-weight calculus.
    pub fn nim_heap(n: u128) -> Game {
        let opts: Vec<Game> = (0..n).map(Game::nim_heap).collect();
        Game::new(opts.clone(), opts)
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

    /// The **ordinal sum** `G : H` ("`G` then `H`"): play in the subordinate `H`
    /// freely, but a move into the base `G` discards `H` entirely. Recursively
    /// `G : H = { G^L, G:H^L | G^R, G:H^R }` — the `G`-moves go to the bare
    /// `G`-options (no `:H`), the `H`-moves keep the base. Not commutative, and
    /// distinct from the disjunctive [`add`](Self::add). (Berlekamp's Hackenbush
    /// strings are ordinal sums of single edges.)
    pub fn ordinal_sum(&self, h: &Game) -> Game {
        let mut left: Vec<Game> = self.left().to_vec();
        for hl in h.left() {
            left.push(self.ordinal_sum(hl));
        }
        let mut right: Vec<Game> = self.right().to_vec();
        for hr in h.right() {
            right.push(self.ordinal_sum(hr));
        }
        Game::new(left, right)
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

    /// Whether `G` is **all-small**: at every position, there is a Left option iff
    /// there is a Right option. The all-small games are the infinitesimally-small
    /// ones (built from `0`, `⋆`, `↑`, …) on which the atomic weight is defined;
    /// numbers and switches are *not* all-small.
    pub fn is_all_small(&self) -> bool {
        if self.left().is_empty() != self.right().is_empty() {
            return false;
        }
        self.left()
            .iter()
            .chain(self.right())
            .all(|g| g.is_all_small())
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

/// A transfinite **number-valued** game, carried by its surreal value rather than
/// a (necessarily infinite) option tree. Numbers are the one transfinite class
/// needing no materialized options: value, birthday, and the group/order
/// operations all come from [`Surreal`]. The finite [`Game`] engine is untouched
/// — `NumberGame` is a parallel *view*, not a `Game`, the numbers-only honoring
/// of "games of transfinite birthday" (`ω = {0,1,2,…|}` is a number).
#[derive(Clone, Debug, PartialEq)]
pub struct NumberGame {
    value: Surreal,
}

impl NumberGame {
    /// The number-game of a surreal value (always succeeds — no options built).
    pub fn from_surreal(s: &Surreal) -> NumberGame {
        NumberGame { value: s.clone() }
    }

    /// The exact surreal value.
    pub fn value(&self) -> &Surreal {
        &self.value
    }

    /// The birthday as an [`Ordinal`], via [`Surreal::birthday_ordinal`]. `None`
    /// when the value is outside the representable sign-expansion subclass (e.g.
    /// `√ω`).
    pub fn birthday(&self) -> Option<Ordinal> {
        self.value.birthday_ordinal()
    }

    /// Negation (additive inverse) — surreal negation.
    pub fn neg(&self) -> NumberGame {
        NumberGame {
            value: self.value.neg(),
        }
    }

    /// Disjunctive sum: for numbers this is exactly surreal addition (no options
    /// materialized).
    pub fn add(&self, other: &NumberGame) -> NumberGame {
        NumberGame {
            value: self.value.add(&other.value),
        }
    }

    /// The game order = the surreal order on values.
    pub fn cmp(&self, other: &NumberGame) -> Ordering {
        self.value.cmp(&other.value)
    }

    /// Bridge to the finite engine: `Some(short Game)` iff the value is dyadic;
    /// `None` for genuinely transfinite numbers (`ω`, `ε`, …), which have no
    /// finite option tree. On dyadics this agrees with
    /// [`Game::from_surreal`]/[`Game::number_value`].
    pub fn to_finite_game(&self) -> Option<Game> {
        Game::from_surreal(&self.value)
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
    fn ordinal_sum_basics() {
        // 0 is a left identity, and `G:0 = G` structurally.
        assert!(Game::zero().ordinal_sum(&Game::up()).eq(&Game::up()));
        assert!(Game::switch(2, -1)
            .ordinal_sum(&Game::zero())
            .structural_eq(&Game::switch(2, -1)));
        // ⋆ : ⋆ = ⋆2 (a 2-edge green Hackenbush path).
        let star2 = Game::new(
            vec![Game::integer(0), Game::star()],
            vec![Game::integer(0), Game::star()],
        );
        assert!(Game::star().ordinal_sum(&Game::star()).eq(&star2));
        // ordinal sum of positive integers is ordinary addition: 1:1 = 2.
        assert!(Game::integer(1)
            .ordinal_sum(&Game::integer(1))
            .eq(&Game::integer(2)));
        // not commutative in general: 1:⋆ ≠ ⋆:1.
        assert!(!Game::integer(1)
            .ordinal_sum(&Game::star())
            .eq(&Game::star().ordinal_sum(&Game::integer(1))));
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
    fn number_game_transfinite_bridge() {
        use crate::scalar::{Ordinal, Rational, Surreal};
        let w = Surreal::omega();
        let ng = NumberGame::from_surreal(&w);
        assert_eq!(ng.value(), &w);
        assert_eq!(ng.birthday(), Some(Ordinal::omega())); // birthday(ω) = ω
        assert!(ng.to_finite_game().is_none()); // ω is not a short game
                                                // ω + 1 by pure surreal delegation; order against a big finite number.
        let one = NumberGame::from_surreal(&Surreal::from_int(1));
        assert_eq!(ng.add(&one).value(), &w.add(&Surreal::from_int(1)));
        assert_eq!(
            ng.cmp(&NumberGame::from_surreal(&Surreal::from_int(1_000_000))),
            Ordering::Greater
        );
        assert_eq!(ng.neg().value(), &w.neg());
        // On the dyadic overlap the transfinite birthday matches the finite game's,
        // and the downcast to a short game succeeds.
        let d = Surreal::from_rational(Rational::new(3, 4));
        let ngd = NumberGame::from_surreal(&d);
        let fin = Game::from_surreal(&d).unwrap();
        assert_eq!(ngd.birthday().unwrap().as_finite(), Some(fin.birthday()));
        assert!(ngd.to_finite_game().is_some());
    }

}
