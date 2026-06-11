//! Transfinite **nimber-valued** (impartial) games carried by their ordinal
//! Grundy value — the characteristic-2 mirror of
//! [`number_game`](crate::games::number_game).
//!
//! This closes the project's central **`No ↔ On₂`** symmetry at the *games* layer.
//! That mirror already lives at the scalar layer — [`Surreal`](crate::scalar::Surreal)
//! (the char-0 field `No`) and [`Ordinal`] (the char-2
//! algebraically-closed field `On₂`) share one Cantor-normal-form core
//! ([`big::cnf`](crate::scalar::big)) — but the games column was lopsided: surreal
//! numbers had [`NumberGame`](crate::games::NumberGame) (transfinite games carried
//! by a surreal value), while the impartial side had only the *finite*
//! [`Game::nim_heap`](crate::games::Game::nim_heap). `NimberGame` is the missing
//! half: a transfinite Nim heap `⋆α`, carried by its ordinal Grundy value.
//!
//! The mirror is exact. The Grundy value of the Nim heap of ordinal size `α` is `α`
//! itself, so — exactly as a number game needs no materialized option tree because
//! [`Surreal`](crate::scalar::Surreal) carries everything — a single [`Ordinal`] (the `On₂` backend) carries
//! the Grundy value, the disjunctive sum (Sprague–Grundy XOR = nim-addition), and
//! the Turning-Corners product (nim-multiplication). Where `NumberGame` delegates to
//! `Surreal`, `NimberGame` delegates to `Ordinal`; the two differ exactly where `No`
//! and `On₂` differ — the order vs. the coefficient merge (`+` vs. `XOR`) — which is
//! why this is a parallel *view*, not a shared type. The finite [`Game`] engine is
//! untouched.
//!
//! The one genuine asymmetry is honest and structural: in characteristic 2 every
//! impartial game is its own additive inverse (`G + G = 0`), so [`neg`](NimberGame::neg)
//! is the identity, and the product is a *separate* game (Turning-Corners), not the
//! disjunctive sum — whereas surreal addition and multiplication are both there in
//! the field. That mismatch is the `+`-vs-`XOR` / field-vs-nonfield content of the
//! `No ↔ On₂` mirror, not a gap.

use crate::games::Game;
use crate::scalar::Ordinal;
use std::cmp::Ordering;

/// A transfinite **nimber-valued** (impartial) game — the Nim heap `⋆α` — carried
/// by its ordinal Grundy value rather than a (necessarily infinite) option set. The
/// char-2 mirror of [`NumberGame`](crate::games::NumberGame): that one carries a
/// [`Surreal`](crate::scalar::Surreal) (`No`, char 0), this carries an
/// [`Ordinal`] (`On₂`, char 2), the two sharing one CNF core.
#[derive(Clone, Debug, PartialEq)]
pub struct NimberGame {
    grundy: Ordinal,
}

impl NimberGame {
    /// The transfinite Nim heap `⋆α` of a given ordinal Grundy value (always
    /// succeeds — no options built). The char-2 mirror of
    /// [`NumberGame::from_surreal`](crate::games::NumberGame::from_surreal).
    pub fn from_ordinal(o: &Ordinal) -> NimberGame {
        NimberGame { grundy: o.clone() }
    }

    /// The finite Nim heap `⋆n`.
    pub fn nim_heap(n: u128) -> NimberGame {
        NimberGame {
            grundy: Ordinal::from_u128(n),
        }
    }

    /// The exact Grundy value (a transfinite nimber). The char-2 mirror of
    /// [`NumberGame::value`](crate::games::NumberGame::value); the game is a
    /// P-position (previous-player win) iff this is `0`.
    pub fn grundy(&self) -> &Ordinal {
        &self.grundy
    }

    /// Disjunctive sum: the Sprague–Grundy XOR of the Grundy values
    /// (nim-addition on `On₂`). Always defined — nim-addition is complete on the
    /// represented CNF (`⋆α + ⋆α = 0`, `⋆ω + ⋆1 = ⋆(ω+1)`). The mirror of
    /// [`NumberGame::add`](crate::games::NumberGame::add) (surreal addition).
    pub fn add(&self, other: &NimberGame) -> NimberGame {
        NimberGame {
            grundy: self.grundy.nim_add(&other.grundy),
        }
    }

    /// Negation — the **identity**. Every impartial game is its own additive inverse
    /// (`G + G = 0`; in char 2, `−α = α`). The mirror of
    /// [`NumberGame::neg`](crate::games::NumberGame::neg), trivial on this leg by the
    /// char-2 structure.
    pub fn neg(&self) -> NimberGame {
        self.clone()
    }

    /// The **Turning-Corners product**: Conway's coin game realizing
    /// nim-multiplication (the transfinite extension of
    /// [`coin_turning::nim_mul_mex`](crate::games::nim_mul_mex)). Defined across the
    /// `On₂` prime-power tower, including the non-scalar Kummer branching (`α_7 = ω+1`,
    /// …); `None` only when a Kummer carry needs a prime `> 47` (past the verified
    /// excess table) or at `≥ ⋆ω^(ω^ω)` (see [`big::ordinal`](crate::scalar::big)).
    /// Unlike the surreal leg — where the product is field multiplication — for nimbers
    /// the product is a *separate* game from the disjunctive sum; this is the seam where
    /// the game pillar meets the nimber field (`⋆ω ⊗ ⋆ω ⊗ ⋆ω = ⋆2`, Conway's `ω³ = 2`).
    pub fn turning_corners(&self, other: &NimberGame) -> Option<NimberGame> {
        self.grundy
            .nim_mul(&other.grundy)
            .map(|grundy| NimberGame { grundy })
    }

    /// The **heap-size order** (the ordinal order on Grundy values). The mirror of
    /// [`NumberGame::cmp`](crate::games::NumberGame::cmp). Note this is the *ordinal*
    /// order, deliberately distinct from the nim-value structure (`On₂` as a field is
    /// unordered) — see [`Ordinal::cmp`](crate::scalar::Ordinal::cmp).
    // Inherent value-order, kept off `std::cmp::Ord` to mirror `Ordinal::cmp` and
    // the partial `Game` order (see AGENTS.md).
    #[allow(clippy::should_implement_trait)]
    pub fn cmp(&self, other: &NimberGame) -> Ordering {
        self.grundy.cmp(&other.grundy)
    }

    /// Bridge to the finite engine: `Some(Game::nim_heap(n))` iff the Grundy value
    /// is finite (`< ω`); `None` for genuinely transfinite heaps, which have no
    /// finite option tree. The char-2 mirror of
    /// [`NumberGame::to_finite_game`](crate::games::NumberGame::to_finite_game)
    /// (`Some` iff dyadic). On finite heaps the finite game's value agrees.
    pub fn to_finite_game(&self) -> Option<Game> {
        self.grundy.as_finite().map(Game::nim_heap)
    }
}

impl std::ops::Add for NimberGame {
    type Output = NimberGame;

    fn add(self, rhs: NimberGame) -> NimberGame {
        NimberGame::add(&self, &rhs)
    }
}

impl std::ops::Neg for NimberGame {
    type Output = NimberGame;

    fn neg(self) -> NimberGame {
        NimberGame::neg(&self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn transfinite_bridge() {
        // ⋆ω: Grundy value ω, no finite option tree (mirror of the ω number game).
        let w = Ordinal::omega();
        let g = NimberGame::from_ordinal(&w);
        assert_eq!(g.grundy(), &w);
        assert!(g.to_finite_game().is_none());

        // disjunctive sum = nim-add; ⋆ω + ⋆ω = 0 (a P-position); neg is identity.
        assert!(g.add(&g).grundy().is_zero());
        assert_eq!(g.neg(), g);
        let one = NimberGame::nim_heap(1);
        assert_eq!(g.add(&one).grundy(), &w.nim_add(&Ordinal::from_u128(1))); // ω+1

        // the heap-size (ordinal) order: ⋆ω dominates every finite heap.
        assert_eq!(g.cmp(&NimberGame::nim_heap(1_000_000)), Ordering::Greater);

        // a finite heap bridges to the finite engine, value-for-value.
        let fin = NimberGame::nim_heap(2);
        assert!(fin.to_finite_game().unwrap().eq(&Game::nim_heap(2)));
    }

    #[test]
    fn turning_corners_is_nim_multiplication() {
        // Below ω^ω the product is the genuine On₂ nim-product. ⋆2 ⊗ ⋆3:
        let two = NimberGame::nim_heap(2);
        let three = NimberGame::nim_heap(3);
        let prod = two.turning_corners(&three).unwrap();
        assert_eq!(
            prod.grundy().as_finite(),
            Ordinal::from_u128(2)
                .nim_mul(&Ordinal::from_u128(3))
                .unwrap()
                .as_finite()
        );

        // Conway's ω³ = 2: ⋆ω ⊗ ⋆ω ⊗ ⋆ω = ⋆2, built from coin games.
        let w = NimberGame::from_ordinal(&Ordinal::omega());
        let w3 = w
            .turning_corners(&w)
            .and_then(|w2| w2.turning_corners(&w))
            .unwrap();
        assert_eq!(w3.grundy(), &Ordinal::from_u128(2));

        // ⋆ω^ω is now the degree-5 generator χ_5: ⋆ω^ω ⊗ ⋆ω^ω = ⋆ω^(ω·2) (was staged
        // under the old ω^ω boundary).
        let ww = NimberGame::from_ordinal(&Ordinal::omega_pow(Ordinal::omega()));
        assert_eq!(
            ww.turning_corners(&ww).unwrap().grundy(),
            &Ordinal::omega_pow(Ordinal::monomial(Ordinal::from_u128(1), 2))
        );
        // the staged boundary is now the non-scalar Kummer carry: ⋆ω^(ω^ω) is None.
        let www =
            NimberGame::from_ordinal(&Ordinal::omega_pow(Ordinal::omega_pow(Ordinal::omega())));
        assert!(www.turning_corners(&www).is_none());
    }

    #[test]
    fn additive_operator_traits_delegate_to_nim_arithmetic() {
        let two = NimberGame::nim_heap(2);
        let three = NimberGame::nim_heap(3);

        assert_eq!(
            (two.clone() + three.clone()).grundy(),
            &Ordinal::from_u128(2).nim_add(&Ordinal::from_u128(3))
        );
        assert_eq!(-two.clone(), two);
        assert_eq!(two.turning_corners(&three), Some(NimberGame::nim_heap(1)));
    }

    #[test]
    fn mirrors_the_number_game_on_the_shared_cnf() {
        // The structural point: NimberGame is to On₂ what NumberGame is to No, and
        // both read the same CNF tower. ω+1 here (XOR-merge) vs. there (+-merge):
        // the heaps add by nim-addition, never collapsing equal ω-powers the way
        // ordinary ordinal addition would (ω + ω = ω·2, but ⋆ω + ⋆ω = 0).
        let w = NimberGame::from_ordinal(&Ordinal::omega());
        assert!(w.add(&w).grundy().is_zero(), "⋆ω + ⋆ω = 0 (XOR, not ω·2)");
        let wp1 = w.add(&NimberGame::nim_heap(1));
        assert_eq!(format!("{:?}", wp1.grundy()), "*(ω + 1)");
    }
}
