//! Temperature theory for short partizan games: stops, cooling, and the
//! **thermograph** (mean value + temperature).
//!
//! This is the half of combinatorial game theory that complements canonical
//! form and Grundy theory: where canonical form answers "*what game is this?*",
//! thermography answers "*how big, and how urgent, is the move?*".
//!
//! ## The construction
//!
//! Cool a game by a temperature `t`:
//! ```text
//! G_t = { G^L_t − t | G^R_t + t }     (until it becomes a number, then it freezes)
//! ```
//! As `t` rises the Left and Right *scaffolds* march toward each other; the
//! temperature `t(G)` is where they meet and the game freezes to its **mean
//! value** (the mast). We compute the whole **thermograph** — the two scaffolds
//! as piecewise-linear functions of `t` — exactly, in dyadic rationals.
//!
//! Standard facts the tests pin:
//! - a switch `{a|b}` (`a > b`) has temperature `(a−b)/2` and mean `(a+b)/2`;
//! - a *number* has temperature `−1` (the coldest) and is its own mean;
//! - the all-small games `⋆`, `↑` have temperature `0`, mean `0`;
//! - the **mean value is additive**: `mean(G+H) = mean(G) + mean(H)`.
//!
//! The base case feeds option *walls* (not temperatures) into the parent
//! recursion, so the `−1` number convention is for reporting only and never
//! contaminates a hot computation.

pub use crate::games::piecewise::Pl;
use crate::games::piecewise::{combine, rdiv, sub_pl};
use crate::games::Game;
use crate::scalar::{Rational, Scalar};
use std::cmp::Ordering;

/// Least `t ≥ 0` where `D(t) = E(t) − 2t = 0`, given `D(0) ≥ 0` and `D → −∞`.
/// Here `E = left_raw − right_raw`, so this is the temperature at which the
/// scaffolds meet.
pub(crate) fn meeting_temperature(e: &Pl) -> Rational {
    let two = Rational::int(2);
    let d_at = |t: &Rational| e.value_at(t).sub(&two.mul(t));

    let t0 = Rational::zero();
    if d_at(&t0).sign() != Ordering::Greater {
        return t0; // D(0) ≤ 0 ⇒ already meeting (infinitesimal/tepid game)
    }
    let ts: Vec<Rational> = e.pts.iter().map(|(t, _)| t.clone()).collect();
    for w in ts.windows(2) {
        let (a, b) = (&w[0], &w[1]);
        let (da, db) = (d_at(a), d_at(b));
        if da.sign() == Ordering::Equal {
            return a.clone();
        }
        if da.sign() == Ordering::Greater && db.sign() != Ordering::Greater {
            // zero in (a,b]: t = a + (b−a)·da/(da−db)
            let denom = da.sub(&db);
            return a.add(&b.sub(a).mul(&rdiv(&da, &denom)));
        }
    }
    // final ray: E constant = last value, D = last − 2t, zero at last/2
    let last = e.pts.last().unwrap().1.clone();
    rdiv(&last, &two)
}

/// Freeze a Left scaffold: `left_raw(t) − t` below `τ`, then the mast.
pub(crate) fn freeze(raw: &Pl, tau: &Rational, mast: &Rational, left: bool) -> Pl {
    let mut pts = Vec::new();
    for (t, v) in &raw.pts {
        if t.cmp(tau) == Ordering::Less {
            let shifted = if left { v.sub(t) } else { v.add(t) };
            pts.push((t.clone(), shifted));
        }
    }
    pts.push((tau.clone(), mast.clone()));
    Pl { pts }.cleaned()
}

/// The thermograph of a short partizan game: its mean value, temperature, and
/// the two scaffold walls as functions of temperature.
#[derive(Clone, Debug)]
pub struct Thermograph {
    /// The mean (mast) value — the number the game freezes to.
    pub mast: Rational,
    /// The temperature: where the scaffolds meet. `−1` for a number.
    pub temperature: Rational,
    /// The Left scaffold (left wall), constant `= mast` above the temperature.
    pub left_wall: Pl,
    /// The Right scaffold (right wall).
    pub right_wall: Pl,
}

impl Thermograph {
    /// The mean value (mast).
    pub fn mean(&self) -> Rational {
        self.mast.clone()
    }
    /// The Left stop `LS(G)` — the left wall at temperature `0`.
    pub fn left_stop(&self) -> Rational {
        self.left_wall.value_at(&Rational::zero())
    }
    /// The Right stop `RS(G)` — the right wall at temperature `0`.
    pub fn right_stop(&self) -> Rational {
        self.right_wall.value_at(&Rational::zero())
    }
    /// The cooled stops `(LS, RS)` at temperature `t`. Both equal the mast once
    /// `t` reaches the temperature (the game has frozen).
    pub fn cooled_stops(&self, t: &Rational) -> (Rational, Rational) {
        (self.left_wall.value_at(t), self.right_wall.value_at(t))
    }
}

/// Shared recursion: returns `(left_wall, right_wall, mast, temperature)`, or
/// `None` if a non-number canonical position has an empty option side (outside
/// the domain of ordinary temperature theory). The caller supplies the fold used
/// to combine option walls; ordinary thermography uses `combine`, while the
/// tropical naming layer routes the same recursion through `oplus_max/min`.
pub(crate) fn walls_with<F>(g: &Game, fold: F) -> Option<(Pl, Pl, Rational, Rational)>
where
    F: Fn(&Pl, &Pl, bool) -> Pl + Copy,
{
    let g = g.canonical();
    if g.is_number() {
        let v = g.number_value()?.as_rational()?; // dyadic ⇒ rational
        let c = Pl::constant(v.clone());
        return Some((c.clone(), c, v, Rational::int(-1)));
    }
    if g.left().is_empty() || g.right().is_empty() {
        return None;
    }
    // left_raw = max over Left options of the option's RIGHT wall
    let mut left_raw: Option<Pl> = None;
    for l in g.left() {
        let rw = walls_with(l, fold)?.1;
        left_raw = Some(match left_raw {
            None => rw,
            Some(acc) => fold(&acc, &rw, true),
        });
    }
    // right_raw = min over Right options of the option's LEFT wall
    let mut right_raw: Option<Pl> = None;
    for r in g.right() {
        let lw = walls_with(r, fold)?.0;
        right_raw = Some(match right_raw {
            None => lw,
            Some(acc) => fold(&acc, &lw, false),
        });
    }
    let (left_raw, right_raw) = (left_raw.unwrap(), right_raw.unwrap());
    let e = sub_pl(&left_raw, &right_raw);
    let tau = meeting_temperature(&e);
    let mast = left_raw.value_at(&tau).sub(&tau);
    let left_wall = freeze(&left_raw, &tau, &mast, true);
    let right_wall = freeze(&right_raw, &tau, &mast, false);
    Some((left_wall, right_wall, mast, tau))
}

fn walls(g: &Game) -> Option<(Pl, Pl, Rational, Rational)> {
    walls_with(g, combine)
}

/// The thermograph of `g`, or `None` for the (rare, post-canonicalization)
/// degenerate positions outside temperature theory's domain.
pub fn thermograph(g: &Game) -> Option<Thermograph> {
    let (left_wall, right_wall, mast, temperature) = walls(g)?;
    Some(Thermograph {
        mast,
        temperature,
        left_wall,
        right_wall,
    })
}

/// The temperature `t(G)` (`−1` for a number).
pub fn temperature(g: &Game) -> Option<Rational> {
    thermograph(g).map(|t| t.temperature)
}

/// The mean (mast) value of `g`.
pub fn mean_value(g: &Game) -> Option<Rational> {
    thermograph(g).map(|t| t.mast)
}

/// The Left stop `LS(G)` at temperature 0.
pub fn left_stop(g: &Game) -> Option<Rational> {
    thermograph(g).map(|t| t.left_stop())
}

/// The Right stop `RS(G)` at temperature 0.
pub fn right_stop(g: &Game) -> Option<Rational> {
    thermograph(g).map(|t| t.right_stop())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::games::piecewise::req;

    fn rat(n: i128, d: i128) -> Rational {
        Rational::new(n, d)
    }
    fn int(n: i128) -> Rational {
        Rational::int(n)
    }

    #[test]
    fn numbers_are_cold() {
        for n in [-3i128, 0, 5] {
            let th = thermograph(&Game::integer(n)).unwrap();
            assert!(req(&th.mast, &int(n)));
            assert!(req(&th.temperature, &int(-1)));
            assert!(req(&th.left_stop(), &int(n)) && req(&th.right_stop(), &int(n)));
        }
        // a dyadic number ½ = {0|1}
        let half = Game::new(vec![Game::integer(0)], vec![Game::integer(1)]);
        let th = thermograph(&half).unwrap();
        assert!(req(&th.mast, &rat(1, 2)));
        assert!(req(&th.temperature, &int(-1)));
    }

    #[test]
    fn switches_have_classic_thermographs() {
        // {a|b}, a > b integers: temperature (a−b)/2, mean (a+b)/2, stops (a,b).
        for (a, b) in [(1i128, -1i128), (2, -2), (3, -1), (0, -4), (5, 1)] {
            let g = Game::switch(a, b);
            let th = thermograph(&g).unwrap();
            assert!(req(&th.temperature, &rat(a - b, 2)), "temp {{{a}|{b}}}");
            assert!(req(&th.mast, &rat(a + b, 2)), "mean {{{a}|{b}}}");
            assert!(req(&th.left_stop(), &int(a)), "LS {{{a}|{b}}}");
            assert!(req(&th.right_stop(), &int(b)), "RS {{{a}|{b}}}");
        }
    }

    #[test]
    fn infinitesimals_have_temperature_zero() {
        for g in [Game::star(), Game::up(), Game::up().neg()] {
            let th = thermograph(&g).unwrap();
            assert!(req(&th.temperature, &int(0)), "temp of {}", g.display());
            assert!(req(&th.mast, &int(0)), "mean of {}", g.display());
        }
        // *2 = {0,*|0,*}
        let star2 = Game::new(
            vec![Game::integer(0), Game::star()],
            vec![Game::integer(0), Game::star()],
        );
        let th = thermograph(&star2).unwrap();
        assert!(req(&th.temperature, &int(0)));
        assert!(req(&th.mast, &int(0)));
    }

    #[test]
    fn nested_hot_game() {
        // {3 | {1|−1}}: right option is the switch ±1 (temp 1, mean 0).
        let g = Game::new(vec![Game::integer(3)], vec![Game::switch(1, -1)]);
        let th = thermograph(&g).unwrap();
        assert!(req(&th.mast, &rat(3, 2)), "mast {:?}", th.mast);
        assert!(
            req(&th.temperature, &rat(3, 2)),
            "temp {:?}",
            th.temperature
        );
        assert!(req(&th.left_stop(), &int(3)));
        assert!(req(&th.right_stop(), &int(1))); // RS = LS(switch) = 1
    }

    #[test]
    fn mean_value_is_additive() {
        let g = Game::switch(1, -1); // mean 0
        let h = Game::switch(2, 0); // mean 1
        let mg = mean_value(&g).unwrap();
        let mh = mean_value(&h).unwrap();
        let mgh = mean_value(&g.add(&h)).unwrap();
        assert!(req(&mgh, &mg.add(&mh)), "mean(G+H) = {:?}", mgh);

        // a colder + hotter pair
        let a = Game::switch(4, -4); // mean 0
        let b = Game::integer(3); // mean 3
        assert!(req(&mean_value(&a.add(&b)).unwrap(), &int(3)));
    }

    #[test]
    fn cooled_stops_freeze_at_temperature() {
        let g = Game::switch(3, -1); // temp 2, mean 1
        let th = thermograph(&g).unwrap();
        // below temperature the stops are still apart
        let (l1, r1) = th.cooled_stops(&int(1));
        assert!(l1.cmp(&r1) == Ordering::Greater);
        // at/above temperature both equal the mast
        let (l2, r2) = th.cooled_stops(&int(2));
        assert!(req(&l2, &th.mast) && req(&r2, &th.mast));
        let (l3, r3) = th.cooled_stops(&int(5));
        assert!(req(&l3, &th.mast) && req(&r3, &th.mast));
    }
}
