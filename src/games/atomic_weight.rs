//! **Atomic weight** of all-small games — the piece of temperature theory that
//! [`thermography`](crate::games::thermography) explicitly deferred.
//!
//! For an *all-small* game `G` (every position has a Left move iff a Right move),
//! the atomic weight measures its size in copies of `↑`. It is computed by
//! Siegel's **Constructive Atomic Weight** theorem (the "two-ahead rule"; Conway
//! *ONAG* ch. 16 / Winning Ways; the implementable restatement is Larsson &
//! Nowakowski, *arXiv:2007.03949*, Thm 10).
//!
//! The calculus. Build the candidate from the *raw* option weights:
//! ```text
//!   A = { aw(G^L) − 2 | aw(G^R) + 2 }            (as a game)
//!   if A is not an integer:  aw(G) = A
//!   else (A ∈ ℤ), with the far star ⋆N (N > every nimber in G):
//!       G ‖ ⋆N (or G = ⋆N):   aw(G) = 0
//!       G ⊳ ⋆N (G > ⋆N):      aw(G) = max{ n ∈ ℤ : n ◁| A^R for every A^R }
//!       G ⊲ ⋆N (G < ⋆N):      aw(G) = min{ n ∈ ℤ : n ▷| A^L for every A^L }
//! ```
//! where `n ◁| x` ("less than or confused with") is `!x.le(n)` and `n ▷| x` is
//! `!n.le(x)`. The `±2` shift is the "two-ahead". The integer-branch magnitude
//! is a predicate over `A`'s **own option games** — crucially comparing the
//! *integer* `n` against the (possibly non-integer) option games `A^R = aw(G^R)+2`,
//! so it stays correct when an option's atomic weight is a fraction (e.g. `½`):
//! a naive `1 + max_R aw(G^R)` is wrong there. The positive branch is bounded by
//! the *tightest* (smallest) right option.
//!
//! **Additive.** On all-small games atomic weight is a homomorphism:
//! `aw(G+H) = aw(G) + aw(H)` and `aw(−G) = −aw(G)` (Larsson–Nowakowski Thm 1,
//! restating Siegel) — see `atomic_weight_is_additive`.

use crate::games::Game;

/// `G` as an integer value, if it is an integer-valued game; else `None`.
fn game_as_int(g: &Game) -> Option<i128> {
    let (num, k) = g.number_value()?.as_dyadic()?;
    (k == 0).then_some(num)
}

/// The **atomic weight** of an all-small game, as a `Game` value (usually an
/// integer, occasionally a non-integer game). `None` if `G` is not all-small
/// (the calculus is undefined there).
pub fn atomic_weight(g: &Game) -> Option<Game> {
    if !g.is_all_small() {
        return None;
    }
    let g = g.canonical();
    if g.left().is_empty() && g.right().is_empty() {
        return Some(Game::integer(0)); // aw(0) = 0
    }
    // Option atomic weights (all-small ⇒ each is Some) — may be non-integer games.
    let awl: Vec<Game> = g.left().iter().map(atomic_weight).collect::<Option<_>>()?;
    let awr: Vec<Game> = g.right().iter().map(atomic_weight).collect::<Option<_>>()?;

    // The candidate A = { aw(G^L) − 2 | aw(G^R) + 2 }, kept as RAW option games:
    // the integer-branch predicate compares against these, not their values.
    let a_left: Vec<Game> = awl.iter().map(|a| a.add(&Game::integer(-2))).collect();
    let a_right: Vec<Game> = awr.iter().map(|a| a.add(&Game::integer(2))).collect();
    let a_canon = Game::new(a_left.clone(), a_right.clone()).canonical();

    // If A is not an integer, the candidate value stands.
    let a_int = match game_as_int(&a_canon) {
        None => return Some(a_canon),
        Some(k) => k,
    };

    // Integer case: resolve by the far star ⋆N, N > every nimber in G.
    let far = Game::nim_heap(g.birthday() + 1);
    let le_gf = g.le(&far);
    let le_fg = far.le(&g);
    if le_fg && !le_gf {
        // G ⊳ ⋆N : max{ n : n ◁| A^R for every right option } (n ◁| x ⇔ !x.le(n)).
        let pred = |n: i128| {
            let gn = Game::integer(n);
            a_right.iter().all(|r| !r.le(&gn))
        };
        let mut n = a_int;
        while !pred(n) {
            n -= 1;
        }
        while pred(n + 1) {
            n += 1;
        }
        Some(Game::integer(n))
    } else if le_gf && !le_fg {
        // G ⊲ ⋆N : min{ n : n ▷| A^L for every left option } (n ▷| x ⇔ !n.le(x)).
        let pred = |n: i128| {
            let gn = Game::integer(n);
            a_left.iter().all(|l| !gn.le(l))
        };
        let mut n = a_int;
        while !pred(n) {
            n += 1;
        }
        while pred(n - 1) {
            n -= 1;
        }
        Some(Game::integer(n))
    } else {
        // G ‖ ⋆N (or G = ⋆N): atomic weight 0.
        Some(Game::integer(0))
    }
}

/// The atomic weight as an integer, when it is one — `None` if `G` is not
/// all-small, or its atomic weight is a genuine non-integer game.
pub fn atomic_weight_int(g: &Game) -> Option<i128> {
    game_as_int(&atomic_weight(g)?)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn up_n(n: i128) -> Game {
        Game::up().times_int(n)
    }
    fn aw(g: &Game) -> Option<i128> {
        atomic_weight_int(g)
    }

    #[test]
    fn is_all_small_predicate() {
        assert!(Game::zero().is_all_small());
        assert!(Game::star().is_all_small());
        assert!(Game::nim_heap(2).is_all_small());
        assert!(Game::up().is_all_small());
        assert!(Game::up().add(&Game::star()).is_all_small());
        // numbers and switches are NOT all-small.
        assert!(!Game::integer(3).is_all_small());
        assert!(!Game::switch(1, -1).is_all_small());
        assert!(!Game::from_surreal(&crate::scalar::Surreal::from_int(1))
            .unwrap()
            .is_all_small());
    }

    #[test]
    fn atomic_weight_oracle_table() {
        // The known values (Winning Ways / ONAG).
        assert_eq!(aw(&Game::zero()), Some(0));
        assert_eq!(aw(&Game::star()), Some(0)); // ⋆1
        assert_eq!(aw(&Game::nim_heap(2)), Some(0)); // ⋆2
        assert_eq!(aw(&Game::nim_heap(3)), Some(0)); // ⋆3
        assert_eq!(aw(&Game::up()), Some(1)); // ↑
        assert_eq!(aw(&Game::up().add(&Game::star())), Some(1)); // ↑*
        assert_eq!(aw(&up_n(2)), Some(2)); // ⇑ — the case that breaks a naive ±1 rule
        assert_eq!(aw(&up_n(3)), Some(3)); // ↑↑↑
        assert_eq!(aw(&Game::up().neg()), Some(-1)); // ↓
        assert_eq!(aw(&up_n(2).neg()), Some(-2)); // ⇓
        assert_eq!(aw(&Game::up().neg().add(&Game::star())), Some(-1)); // ↓*
        assert_eq!(aw(&Game::up().add(&Game::nim_heap(2))), Some(1)); // ↑+⋆2
    }

    #[test]
    fn atomic_weight_negation_symmetry() {
        // aw(−G) = −aw(G) on the oracle set.
        for g in [
            Game::up(),
            up_n(2),
            Game::up().add(&Game::star()),
            Game::nim_heap(2),
            Game::zero(),
        ] {
            let a = aw(&g).unwrap();
            assert_eq!(aw(&g.neg()), Some(-a));
        }
    }

    #[test]
    fn non_all_small_has_no_atomic_weight() {
        assert!(atomic_weight(&Game::integer(3)).is_none());
        assert!(atomic_weight(&Game::switch(2, 0)).is_none());
        assert_eq!(atomic_weight_int(&Game::integer(3)), None);
    }

    #[test]
    fn atomic_weight_is_additive() {
        // Atomic weight IS a homomorphism on all-small games (Larsson–Nowakowski
        // Thm 1, restating Siegel): aw(G+H) = aw(G)+aw(H).
        let parts = [
            Game::up(),
            Game::up().times_int(2),
            Game::star(),
            Game::nim_heap(2),
            Game::up().add(&Game::star()),
            Game::up().neg(),
        ];
        for g in &parts {
            for h in &parts {
                assert_eq!(
                    aw(&g.add(h)).unwrap(),
                    aw(g).unwrap() + aw(h).unwrap(),
                    "aw not additive on a pair"
                );
            }
        }
    }

    #[test]
    fn integer_branch_handles_fractional_option_weights() {
        // Codex's counterexample to a naive `1 + max_R aw(G^R)` rule. The right
        // option {⇑ | ↓} has the NON-integer atomic weight ½, so the candidate is
        // A = {0 | 5/2} = 1 (integer), and the far-star branch must give
        // max{ n : n < 5/2 } = 2 — comparing the integer n against the fractional
        // option game, not assuming option weights are integers.
        let h = Game::new(vec![up_n(2)], vec![Game::up().neg()]); // {⇑ | ↓}
        assert!(atomic_weight_int(&h).is_none()); // aw(h) = ½, a non-integer game
        let g = Game::new(vec![up_n(2)], vec![h]); // {⇑ | {⇑|↓}}
        assert_eq!(aw(&g), Some(2));
    }
}
