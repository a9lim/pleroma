//! Transfinite number-valued games carried by their surreal value.

use crate::games::Game;
use crate::scalar::{Ordinal, Scalar, SignExpansion, Surreal};
use std::cmp::Ordering;

/// A transfinite **number-valued** game, carried by its surreal value rather than
/// a (necessarily infinite) option tree. Numbers are a transfinite class needing no
/// materialized options (see also [`NimberGame`](crate::games::NimberGame) for the
/// characteristic-2 impartial mirror): value, birthday, and the group/order
/// operations all come from [`Surreal`]. The finite [`Game`] engine is untouched
/// — `NumberGame` is a parallel *view*, not a `Game`, the numbers-only honoring
/// of "games of transfinite birthday" (`ω = {0,1,2,...|}` is a number).
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
    /// `sqrt(omega)`).
    pub fn birthday(&self) -> Option<Ordinal> {
        self.value.birthday_ordinal()
    }

    /// The **sign expansion** — the canonical ±-path from `0` to this number in
    /// the surreal tree, run-length-encoded (its length is the birthday). This is
    /// the finite encoding of the game's (transfinitely deep) `{Left | Right}`
    /// tree: a transfinite number like `ω = {0,1,2,…|}` has an *infinite* option
    /// set that cannot be listed, but its sign expansion `+^ω` is finite data.
    /// `None` outside the representable subclass.
    pub fn sign_expansion(&self) -> Option<SignExpansion> {
        self.value.transfinite_sign_expansion()
    }

    /// Reconstruct a number-game from a sign expansion — the inverse of
    /// [`sign_expansion`](Self::sign_expansion), closing the transfinite
    /// surreal↔game round trip *through* the canonical birthday path rather than a
    /// stored value: `from_sign_expansion(g.sign_expansion()?) == Some(g)`. This is
    /// the transfinite analogue of the dyadic [`Game::from_surreal`] /
    /// [`Game::number_value`] bridge. `None` outside the representable subclass.
    pub fn from_sign_expansion(se: &SignExpansion) -> Option<NumberGame> {
        Surreal::from_transfinite_sign_expansion(se).map(|value| NumberGame { value })
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
    // Inherent value-order, deliberately kept off `std::cmp::Ord` to mirror
    // `Surreal::cmp` and the partial `Game` order (see AGENTS.md).
    #[allow(clippy::should_implement_trait)]
    pub fn cmp(&self, other: &NumberGame) -> Ordering {
        self.value.cmp(&other.value)
    }

    /// Bridge to the finite engine: `Some(short Game)` iff the value is dyadic;
    /// `None` for genuinely transfinite numbers (`omega`, `epsilon`, ...), which
    /// have no finite option tree. On dyadics this agrees with
    /// [`Game::from_surreal`]/[`Game::number_value`].
    pub fn to_finite_game(&self) -> Option<Game> {
        Game::from_surreal(&self.value)
    }
}

impl std::ops::Add for NumberGame {
    type Output = NumberGame;

    fn add(self, rhs: NumberGame) -> NumberGame {
        NumberGame::add(&self, &rhs)
    }
}

impl std::ops::Neg for NumberGame {
    type Output = NumberGame;

    fn neg(self) -> NumberGame {
        NumberGame::neg(&self)
    }
}

impl std::ops::Mul for NumberGame {
    type Output = NumberGame;

    fn mul(self, rhs: NumberGame) -> NumberGame {
        NumberGame {
            value: <Surreal as Scalar>::mul(&self.value, &rhs.value),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scalar::{Rational, Surreal};

    #[test]
    fn transfinite_bridge() {
        let w = Surreal::omega();
        let ng = NumberGame::from_surreal(&w);
        assert_eq!(ng.value(), &w);
        assert_eq!(ng.birthday(), Some(Ordinal::omega()));
        assert!(ng.to_finite_game().is_none());

        let one = NumberGame::from_surreal(&Surreal::from_int(1));
        assert_eq!(ng.add(&one).value(), &w.add(&Surreal::from_int(1)));
        assert_eq!(
            ng.cmp(&NumberGame::from_surreal(&Surreal::from_int(1_000_000))),
            Ordering::Greater
        );
        assert_eq!(ng.neg().value(), &w.neg());

        let d = Surreal::from_rational(Rational::new(3, 4));
        let ngd = NumberGame::from_surreal(&d);
        let fin = Game::from_surreal(&d).unwrap();
        assert_eq!(ngd.birthday().unwrap().as_finite(), Some(fin.birthday()));
        assert!(ngd.to_finite_game().is_some());
    }

    #[test]
    fn operator_traits_delegate_to_surreal_arithmetic() {
        let two = NumberGame::from_surreal(&Surreal::from_int(2));
        let three = NumberGame::from_surreal(&Surreal::from_int(3));

        assert_eq!((two.clone() + three.clone()).value(), &Surreal::from_int(5));
        assert_eq!((-two.clone()).value(), &Surreal::from_int(-2));
        assert_eq!((two * three).value(), &Surreal::from_int(6));
    }

    #[test]
    fn sign_expansion_round_trip_through_the_game_tree() {
        // The full transfinite round trip: serialize each number-game to its
        // canonical birthday ±-path and reconstruct it, value-for-value — going
        // *through* the sign expansion, not a stored value. Spans dyadic, ordinal,
        // negative-ordinal, and the infinitesimal ε.
        let cases = [
            Surreal::from_int(0),
            Surreal::from_rational(Rational::new(3, 4)),
            Surreal::from_rational(Rational::new(-5, 8)),
            Surreal::omega(),                            // ω = {0,1,2,…|}
            Surreal::omega().add(&Surreal::from_int(1)), // ω+1
            Surreal::omega_pow(Surreal::omega()),        // ω^ω
            Surreal::omega().neg(),                      // −ω
            Surreal::epsilon(),                          // ε
        ];
        for v in &cases {
            let g = NumberGame::from_surreal(v);
            let se = g.sign_expansion().expect("representable");
            let back = NumberGame::from_sign_expansion(&se).expect("reconstructible");
            assert_eq!(back, g, "game sign-expansion round trip failed: {v:?}");
            assert_eq!(back.value(), v);
        }
        // sqrt(ω) is outside the representable subclass — honestly None, no crash.
        let root_omega = Surreal::omega_pow(Surreal::from_rational(Rational::new(1, 2)));
        assert!(NumberGame::from_surreal(&root_omega)
            .sign_expansion()
            .is_none());
    }
}
