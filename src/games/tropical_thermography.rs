//! Thermography **is** tropical arithmetic — named and machine-checked.
//!
//! The thermograph recursion in [`crate::games::thermography`] already *computes*
//! a tropical structure without naming it. This module names it and pins the
//! naming faithful against the existing (golden-tested) [`thermograph`](crate::games::thermograph).
//!
//! ## The correspondence (standard math, not a new theorem)
//!
//! For a short partizan game `G` the two scaffold walls are built by folding the
//! options' walls:
//!
//! - the **left** raw wall is the pointwise `max` over Left options of each
//!   option's *right* wall — that is tropical `⊕` in the **(max, +)** semiring
//!   ([`MaxPlus`](crate::scalar::MaxPlus));
//! - the **right** raw wall is the pointwise `min` over Right options of each
//!   option's *left* wall — tropical `⊕` in the **dual (min, +)** semiring
//!   ([`MinPlus`](crate::scalar::MinPlus));
//! - **cooling** by a temperature shifts a wall's value by `±t` — tropical `⊗`
//!   (tropical multiplication is ordinary `+`), the named [`Pl::otimes`].
//!
//! The two walls genuinely live in **dual** semirings, which is why the scalar
//! layer makes `Tropical<MaxPlus>` and `Tropical<MinPlus>` distinct types. This
//! is the classical content of temperature theory (Berlekamp–Conway–Guy *Winning
//! Ways*; Conway *ONAG*; Siegel *Combinatorial Game Theory*) made explicit at the
//! type level — **claim level: standard math, implemented-and-tested**.
//!
//! [`thermograph_via_tropical`] runs the shared thermograph recursion with the
//! two option folds routed through the **named** [`Pl::oplus_max`]/[`Pl::oplus_min`]
//! wrappers. Its sole job is to prove the `⊕`-naming is faithful, so it is pinned
//! **equal** to [`thermograph`](crate::games::thermograph) — it is not a second implementation of cooling.

use crate::games::piecewise::{add_pl, combine, Pl};
use crate::games::thermography::{walls_with, Thermograph};
use crate::games::Game;

impl Pl {
    /// Tropical `⊕` in the **(max, +)** convention — the exact pointwise maximum
    /// of two walls. This is the fold the thermograph applies to the Left
    /// options' right walls when building `left_raw`.
    pub fn oplus_max(&self, other: &Pl) -> Pl {
        combine(self, other, true)
    }

    /// Tropical `⊕` in the dual **(min, +)** convention — the exact pointwise
    /// minimum. The fold over the Right options' left walls (`right_raw`).
    pub fn oplus_min(&self, other: &Pl) -> Pl {
        combine(self, other, false)
    }

    /// Tropical `⊗` — the pointwise sum of two walls (tropical multiplication is
    /// ordinary `+` on values). Below the lower of two games' temperatures the
    /// wall of a disjunctive sum is the `⊗` of the component walls, which is why
    /// the mean value is additive.
    pub fn otimes(&self, other: &Pl) -> Pl {
        add_pl(self, other)
    }
}

/// The thermograph of `g`, computed with the option folds named as tropical `⊕`
/// in the dual `(max, +)`/`(min, +)` semirings. Pinned **equal** to
/// [`thermograph`](crate::games::thermograph) (the inline tests are the proof); `None` on the same
/// degenerate positions.
pub fn thermograph_via_tropical(g: &Game) -> Option<Thermograph> {
    let (left_wall, right_wall, mast, temperature) = walls_with(g, |acc, wall, is_max| {
        if is_max {
            acc.oplus_max(wall)
        } else {
            acc.oplus_min(wall)
        }
    })?;
    Some(Thermograph {
        mast,
        temperature,
        left_wall,
        right_wall,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::games::piecewise::req;
    use crate::games::thermography::thermograph;
    use crate::scalar::Rational;

    /// Two thermographs are equal iff their masts, temperatures, and both walls
    /// (as breakpoint lists) agree. The walls are byte-identical here because the
    /// bridge runs the *same* `combine`/`freeze`, so structural `==` is exact.
    fn same(a: &Thermograph, b: &Thermograph) -> bool {
        req(&a.mast, &b.mast)
            && req(&a.temperature, &b.temperature)
            && a.left_wall.points() == b.left_wall.points()
            && a.right_wall.points() == b.right_wall.points()
    }

    fn star2() -> Game {
        // *2 = {0,* | 0,*}
        Game::new(
            vec![Game::integer(0), Game::star()],
            vec![Game::integer(0), Game::star()],
        )
    }

    #[test]
    fn via_tropical_matches_thermograph() {
        let mut games = vec![
            Game::integer(-3),
            Game::integer(0),
            Game::integer(5),
            // ½ = {0|1}
            Game::new(vec![Game::integer(0)], vec![Game::integer(1)]),
            Game::star(),
            Game::up(),
            Game::up().neg(),
            star2(),
            // nested hot game {3 | {1|−1}}
            Game::new(vec![Game::integer(3)], vec![Game::switch(1, -1)]),
        ];
        for (a, b) in [(1i128, -1i128), (2, -2), (3, -1), (0, -4), (5, 1)] {
            games.push(Game::switch(a, b));
        }
        // a sum of two hot games
        games.push(Game::switch(4, -4).add(&Game::switch(2, 0)));

        for g in &games {
            let golden = thermograph(g);
            let named = thermograph_via_tropical(g);
            match (&golden, &named) {
                (Some(x), Some(y)) => assert!(same(x, y), "mismatch on {}", g.display()),
                (None, None) => {}
                _ => panic!("Some/None disagreement on {}", g.display()),
            }
        }
    }

    #[test]
    fn oplus_wrappers_pin_to_combine() {
        // Use genuine (non-constant) walls from a switch's thermograph.
        let th = thermograph(&Game::switch(3, -1)).unwrap();
        let (l, r) = (&th.left_wall, &th.right_wall);
        assert_eq!(l.oplus_max(r).points(), combine(l, r, true).points());
        assert_eq!(l.oplus_min(r).points(), combine(l, r, false).points());
    }

    #[test]
    fn otimes_pins_to_add_pl() {
        let th = thermograph(&Game::switch(3, -1)).unwrap();
        let (l, r) = (&th.left_wall, &th.right_wall);
        assert_eq!(l.otimes(r).points(), add_pl(l, r).points());
    }

    #[test]
    fn dual_semiring_recovers_classic_stops() {
        // For {a|b} (a > b): the left wall (a (max,+)/⊕ fold) freezes to LS = a,
        // the right wall (a dual (min,+)/⊕ fold) to RS = b. The two stops thus
        // come out of the two *dual* tropical semirings — the reason the scalar
        // layer keeps `Tropical<MaxPlus>` and `Tropical<MinPlus>` distinct types.
        for (a, b) in [(1i128, -1i128), (5, 1), (3, -1)] {
            let th = thermograph_via_tropical(&Game::switch(a, b)).unwrap();
            assert!(req(&th.left_stop(), &Rational::int(a)), "LS {{{a}|{b}}}");
            assert!(req(&th.right_stop(), &Rational::int(b)), "RS {{{a}|{b}}}");
        }
    }
}
