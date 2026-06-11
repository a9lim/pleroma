//! The canonical-stopper catalogue: [`LoopyWinner`], [`LoopyPartizanOutcome`],
//! [`PartizanOutcome`], and [`LoopyValue`].

use std::cmp::Ordering;

/// The winner of one of the two starter questions in a finite loopy partizan
/// graph.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum LoopyWinner {
    /// Left can force a win from this starter state.
    Left,
    /// Right can force a win from this starter state.
    Right,
    /// Neither player can force a win; optimal play can be drawn forever.
    Draw,
}

/// The exact two-sided outcome of a partizan loopy position: one result when Left
/// is to move, and one result when Right is to move.
///
/// The classical five outcome classes embed as the cases where the pair is
/// `(Right, Left)` (`P`), `(Left, Right)` (`N`), `(Left, Left)` (`L`),
/// `(Right, Right)` (`R`), or `(Draw, Draw)` (`Draw`). Mixed cases such as
/// `tis = (Left, Draw)` are real loopy-partizan values and deliberately do not
/// collapse to a [`PartizanOutcome`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct LoopyPartizanOutcome {
    pub left_to_move: LoopyWinner,
    pub right_to_move: LoopyWinner,
}

impl LoopyPartizanOutcome {
    pub fn new(left_to_move: LoopyWinner, right_to_move: LoopyWinner) -> Self {
        Self {
            left_to_move,
            right_to_move,
        }
    }

    /// The classical partizan outcome class, when this two-sided result lies in
    /// the classical five-class image.
    pub fn partizan_class(&self) -> Option<PartizanOutcome> {
        use LoopyWinner::*;
        match (self.left_to_move, self.right_to_move) {
            (Right, Left) => Some(PartizanOutcome::P),
            (Left, Right) => Some(PartizanOutcome::N),
            (Left, Left) => Some(PartizanOutcome::L),
            (Right, Right) => Some(PartizanOutcome::R),
            (Draw, Draw) => Some(PartizanOutcome::Draw),
            _ => None,
        }
    }

    pub fn has_draw(&self) -> bool {
        self.left_to_move == LoopyWinner::Draw || self.right_to_move == LoopyWinner::Draw
    }
}

/// The outcome class of a (partizan, possibly loopy) game value: who wins under
/// optimal play. Unlike the impartial [`Outcome`](crate::games::Outcome) (which is keyed on the player to
/// move), this names the partizan class directly, and adds [`Draw`](Self::Draw)
/// for loopy values like `dud`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PartizanOutcome {
    /// Previous player wins (the player who *just* moved) — i.e. the player to move
    /// loses. The class of `0`.
    P,
    /// Next player wins (the player to move). The class of `∗`.
    N,
    /// Left wins regardless of who moves first.
    L,
    /// Right wins regardless of who moves first.
    R,
    /// Neither player can force a win — a draw under best play. The class of `dud`.
    Draw,
}

/// A catalogue of named loopy values, plus integer onside/offside (`s&t`) tags.
/// This is not a complete equality theory for loopy games; arithmetic returns
/// `None` whenever a sum leaves the represented catalogue or would require a
/// non-local sidling/equality proof.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum LoopyValue {
    /// `0 = {|}` — the second player (previous mover) wins.
    Zero,
    /// `∗ = {0|0}` — the first player (next mover) wins.
    Star,
    /// `on = {on|}` — Right has no move and loses; Left wins regardless. Larger
    /// than every stopper.
    On,
    /// `off = {|off} = −on` — Left has no move and loses; Right wins regardless.
    Off,
    /// `over = {0|over}` — a positive infinitesimal: `0 < over < x` for every
    /// positive number `x`. Left wins.
    Over,
    /// `under = {under|0} = −over` — a negative infinitesimal. Right wins.
    Under,
    /// `± = {on|off}` — the hot next-player loopy switch.
    PlusMinus,
    /// `tis` (`this is`), the left-swinging non-stopper. In this repo's finite
    /// tag convention it records the two-sided result `(Left, Draw)` and the
    /// sidled sides `1&0`.
    Tis,
    /// `tisn` (`this isn't`), the right-swinging dual of [`Tis`](Self::Tis).
    /// It records `(Draw, Right)` and sidled sides `0&-1`.
    Tisn,
    /// A finite onside/offside tag `s&t`. Addition and negation are carried on the
    /// pair itself; equality with arbitrary loopy games is not decided here.
    OnsideOffside { onside: i128, offside: i128 },
    /// `dud = {dud|dud}` — the "deathless universal draw": both players loop
    /// forever, neither wins. Absorbing under sum; confused with every value.
    Dud,
}

impl LoopyValue {
    /// Build an onside/offside value `s&t`, normalizing `0&0` to `0`.
    pub fn onside_offside(onside: i128, offside: i128) -> LoopyValue {
        if onside == 0 && offside == 0 {
            LoopyValue::Zero
        } else {
            LoopyValue::OnsideOffside { onside, offside }
        }
    }

    /// The `{Left | Right}` form, for display.
    pub fn form(&self) -> String {
        match self {
            LoopyValue::Zero => "{|}".to_string(),
            LoopyValue::Star => "{0|0}".to_string(),
            LoopyValue::On => "{on|}".to_string(),
            LoopyValue::Off => "{|off}".to_string(),
            LoopyValue::Over => "{0|over}".to_string(),
            LoopyValue::Under => "{under|0}".to_string(),
            LoopyValue::PlusMinus => "{on|off}".to_string(),
            LoopyValue::Tis => "{0|tisn}".to_string(),
            LoopyValue::Tisn => "{tis|0}".to_string(),
            LoopyValue::OnsideOffside { onside, offside } => {
                format!("{onside}&{offside}")
            }
            LoopyValue::Dud => "{dud|dud}".to_string(),
        }
    }

    /// The conventional name.
    pub fn name(&self) -> String {
        match self {
            LoopyValue::Zero => "0".to_string(),
            LoopyValue::Star => "*".to_string(),
            LoopyValue::On => "on".to_string(),
            LoopyValue::Off => "off".to_string(),
            LoopyValue::Over => "over".to_string(),
            LoopyValue::Under => "under".to_string(),
            LoopyValue::PlusMinus => "±".to_string(),
            LoopyValue::Tis => "tis".to_string(),
            LoopyValue::Tisn => "tisn".to_string(),
            LoopyValue::OnsideOffside { onside, offside } => {
                format!("{onside}&{offside}")
            }
            LoopyValue::Dud => "dud".to_string(),
        }
    }

    /// Who wins under optimal play for each starter. Use
    /// [`partizan_outcome`](Self::partizan_outcome) when you need the classical
    /// five-class projection.
    pub fn outcome(&self) -> LoopyPartizanOutcome {
        use LoopyWinner::*;
        match self {
            LoopyValue::Zero => LoopyPartizanOutcome::new(Right, Left),
            LoopyValue::Star | LoopyValue::PlusMinus => LoopyPartizanOutcome::new(Left, Right),
            LoopyValue::On | LoopyValue::Over => LoopyPartizanOutcome::new(Left, Left),
            LoopyValue::Off | LoopyValue::Under => LoopyPartizanOutcome::new(Right, Right),
            LoopyValue::Tis => LoopyPartizanOutcome::new(Left, Draw),
            LoopyValue::Tisn => LoopyPartizanOutcome::new(Draw, Right),
            LoopyValue::OnsideOffside { onside, offside } => {
                LoopyPartizanOutcome::new(winner_from_sign(*onside), winner_from_sign(*offside))
            }
            LoopyValue::Dud => LoopyPartizanOutcome::new(Draw, Draw),
        }
    }

    /// The classical partizan outcome class, when this value has one. Values such
    /// as `tis` and `tisn` have a mixed draw/win starter pair, so they return
    /// `None` rather than being flattened into a false five-class answer.
    pub fn partizan_outcome(&self) -> Option<PartizanOutcome> {
        self.outcome().partizan_class()
    }

    /// The sidled onside/offside pair when this finite tag carries one.
    pub fn sides(&self) -> Option<(i128, i128)> {
        match *self {
            LoopyValue::Tis => Some((1, 0)),
            LoopyValue::Tisn => Some((0, -1)),
            LoopyValue::OnsideOffside { onside, offside } => Some((onside, offside)),
            _ => None,
        }
    }

    /// Negation (swap the Left/Right roles): `−on = off`, `−over = under`, and the
    /// self-negating `0`, `∗`, `±`, `dud`.
    pub fn neg(&self) -> LoopyValue {
        match self {
            LoopyValue::Zero => LoopyValue::Zero,
            LoopyValue::Star => LoopyValue::Star,
            LoopyValue::On => LoopyValue::Off,
            LoopyValue::Off => LoopyValue::On,
            LoopyValue::Over => LoopyValue::Under,
            LoopyValue::Under => LoopyValue::Over,
            LoopyValue::PlusMinus => LoopyValue::PlusMinus,
            LoopyValue::Tis => LoopyValue::Tisn,
            LoopyValue::Tisn => LoopyValue::Tis,
            LoopyValue::OnsideOffside { onside, offside } => {
                LoopyValue::onside_offside(-*offside, -*onside)
            }
            LoopyValue::Dud => LoopyValue::Dud,
        }
    }

    /// Whether this value is a **stopper** (guaranteed to end when played in
    /// isolation). The named non-stoppers here are `dud`, `tis`, and `tisn`.
    pub fn is_stopper(&self) -> bool {
        !matches!(self, LoopyValue::Dud | LoopyValue::Tis | LoopyValue::Tisn)
    }

    /// The disjunctive sum, where it is defined on this catalogue. Returns `None`
    /// when the sum leaves the catalogue or when this small catalogue deliberately
    /// refuses a drawn value not represented by its named tags.
    ///
    /// The closed cases: `dud` absorbs everything (`dud + G = dud`); `on + off =
    /// dud`; `on`/`off` absorb every other represented stopper (`on` is `>` every
    /// stopper); `∗ + ∗ = 0`; `over + over = over`, `under + under = under`,
    /// `∗ + over = over`, `∗ + under = under`; `s&t + u&v = (s+u)&(t+v)`;
    /// and `0` is the identity.
    pub fn add(&self, other: &LoopyValue) -> Option<LoopyValue> {
        use LoopyValue::*;
        let r = match (*self, *other) {
            (Dud, _) | (_, Dud) => Dud,
            (Zero, x) | (x, Zero) => x,
            (On, On) => On,
            (Off, Off) => Off,
            (On, Off) | (Off, On) => Dud,
            (On, Star)
            | (Star, On)
            | (On, Over)
            | (Over, On)
            | (On, Under)
            | (Under, On)
            | (On, PlusMinus)
            | (PlusMinus, On) => On,
            (Off, Star)
            | (Star, Off)
            | (Off, Over)
            | (Over, Off)
            | (Off, Under)
            | (Under, Off)
            | (Off, PlusMinus)
            | (PlusMinus, Off) => Off,
            (Star, Star) => Zero,
            (Over, Over) | (Star, Over) | (Over, Star) => Over,
            (Under, Under) | (Star, Under) | (Under, Star) => Under,
            (
                OnsideOffside {
                    onside: a,
                    offside: b,
                },
                OnsideOffside {
                    onside: c,
                    offside: d,
                },
            ) => LoopyValue::onside_offside(a + c, b + d),
            (Over, Under) | (Under, Over) => return None,
            (PlusMinus, PlusMinus)
            | (PlusMinus, Star)
            | (Star, PlusMinus)
            | (PlusMinus, Over)
            | (Over, PlusMinus)
            | (PlusMinus, Under)
            | (Under, PlusMinus)
            | (Tis, _)
            | (_, Tis)
            | (Tisn, _)
            | (_, Tisn)
            | (OnsideOffside { .. }, _)
            | (_, OnsideOffside { .. }) => return None,
        };
        Some(r)
    }
}

pub(super) fn winner_from_sign(x: i128) -> LoopyWinner {
    if x > 0 {
        LoopyWinner::Left
    } else if x < 0 {
        LoopyWinner::Right
    } else {
        LoopyWinner::Draw
    }
}

impl PartialOrd for LoopyValue {
    /// The conservative partial order on the catalogue. The comparable core is the
    /// chain `off < under < ∗ < over < on`, with `0` confused with `∗` and between
    /// `under` and `over`. `on` sits above and `off` below every other non-`dud`
    /// value. `dud` is confused with
    /// everything (comparable only to itself). Incomparable ⇒ `None`.
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        use LoopyValue::*;
        if self == other {
            return Some(Ordering::Equal);
        }
        match (*self, *other) {
            // dud is confused with every other value.
            (Dud, _) | (_, Dud) => None,
            // The extended tags need a genuine comparison proof; equality was
            // handled above, so keep the catalogue order conservative.
            (PlusMinus, _)
            | (_, PlusMinus)
            | (Tis, _)
            | (_, Tis)
            | (Tisn, _)
            | (_, Tisn)
            | (OnsideOffside { .. }, _)
            | (_, OnsideOffside { .. }) => None,
            // on is the top, off the bottom (over all non-dud values).
            (On, _) => Some(Ordering::Greater),
            (_, On) => Some(Ordering::Less),
            (Off, _) => Some(Ordering::Less),
            (_, Off) => Some(Ordering::Greater),
            // star is confused with 0, but sits between under and over.
            (Star, Zero) | (Zero, Star) => None,
            (Star, Over) | (Under, Star) => Some(Ordering::Less),
            (Over, Star) | (Star, Under) => Some(Ordering::Greater),
            // the remaining comparable chain under < 0 < over.
            (a, b) => {
                let rank = |v: LoopyValue| match v {
                    Under => -1i128,
                    Zero => 0,
                    Over => 1,
                    _ => unreachable!("on/off/star/dud handled above"),
                };
                Some(rank(a).cmp(&rank(b)))
            }
        }
    }
}
