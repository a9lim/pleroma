//! Loopy combinatorial games — games whose move graph may contain cycles, so
//! play need not terminate. This is the third escape (beside the interactive
//! [`kernel`](crate::games::kernel) route and the [`misere`](crate::games::misere)
//! route) from the XOR-linear P-sets of normal-play disjunctive sums: a cyclic
//! rule admits a **Draw** outcome — a position from which neither player can force
//! a win — and the Draw-set is a genuinely new degree of freedom to test against
//! the Gold quadric `{Q=0}` (see `NOTES.md`, the Tier-2 open question).
//!
//! Three layers, in weight order:
//!
//!   1. [`LoopyGraph`] — the graph-level engine. A thin, fully-computable wrapper
//!      over [`kernel::outcomes`](crate::games::outcomes), which already performs
//!      retrograde Win/Loss/**Draw** analysis on cyclic graphs. This is the
//!      load-bearing part.
//!   2. [`loopy_nim_values`] — the impartial loopy nim-values: Draw positions are
//!      `Side` (the loopy `∞`), the rest carry an ordinary nimber. Acyclic
//!      non-Draw regions use the usual DAG recursion; small cyclic non-Draw
//!      regions use a bounded sidling solver for the finite mex equations.
//!   3. [`LoopyValue`] — a small catalogue of the canonical stoppers (`on`, `off`,
//!      `over`, `under`, `dud`, `∗`) with their outcomes, negation, partial order,
//!      and the partial sum-monoid. A finite tag carrying an infinite object — the
//!      same discipline as [`NumberGame`](crate::games::NumberGame).
//!
//! And the payoff for this project, [`loopy_decision_sets`] / [`loopy_quadric_probe`]:
//! take an arbitrary cyclic move rule on positions `F₂^k` and read off **both** its
//! Loss-set and its Draw-set, fitting each with
//! [`fit_f2_quadratic`](crate::forms::fit_f2_quadratic). A B-coupled cyclic rule
//! whose *Draw-set* is `{Q=0}` would be a Tier-2 witness even if its Loss-set is
//! not — structurally impossible for the acyclic `interactive_kernel` probe.
//!
//! Deliberately **out of scope** here: [`Game`](crate::games::Game) stays an acyclic
//! `Arc` tree (it cannot represent cycles, by construction), and
//! [`thermography`](crate::games::thermography) stays finite-game-only — loopy games
//! never freeze to a number, so classical temperature does not apply. Partizan
//! loopy outcomes (a two-sided Left/Right retrograde solver), unbounded sidling,
//! and the `±`/`tis`/`tisn` stopper arithmetic are honestly deferred.

use std::cmp::Ordering;

use crate::forms::{fit_f2_quadratic, QuadricFit};
use crate::games::grundy::mex;
use crate::games::kernel::{self, Outcome};

const MAX_SIDLING_ASSIGNMENTS: usize = 200_000;

// ---------------------------------------------------------------------------
// 1. The canonical-stopper catalogue.
// ---------------------------------------------------------------------------

/// The outcome class of a (partizan, possibly loopy) game value: who wins under
/// optimal play. Unlike the impartial [`Outcome`] (which is keyed on the player to
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

/// A catalogue of the canonical loopy stoppers (plus `dud`, the canonical
/// non-stopper draw). These are the values with finite names; general loopy values
/// need the onside/offside (`s&t`) machinery, which is out of scope here.
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
    /// `dud = {dud|dud}` — the "deathless universal draw": both players loop
    /// forever, neither wins. Absorbing under sum; confused with every value.
    Dud,
}

impl LoopyValue {
    /// The `{Left | Right}` form, for display.
    pub fn form(&self) -> &'static str {
        match self {
            LoopyValue::Zero => "{|}",
            LoopyValue::Star => "{0|0}",
            LoopyValue::On => "{on|}",
            LoopyValue::Off => "{|off}",
            LoopyValue::Over => "{0|over}",
            LoopyValue::Under => "{under|0}",
            LoopyValue::Dud => "{dud|dud}",
        }
    }

    /// The conventional name.
    pub fn name(&self) -> &'static str {
        match self {
            LoopyValue::Zero => "0",
            LoopyValue::Star => "*",
            LoopyValue::On => "on",
            LoopyValue::Off => "off",
            LoopyValue::Over => "over",
            LoopyValue::Under => "under",
            LoopyValue::Dud => "dud",
        }
    }

    /// Who wins under optimal play.
    pub fn outcome(&self) -> PartizanOutcome {
        match self {
            LoopyValue::Zero => PartizanOutcome::P,
            LoopyValue::Star => PartizanOutcome::N,
            LoopyValue::On | LoopyValue::Over => PartizanOutcome::L,
            LoopyValue::Off | LoopyValue::Under => PartizanOutcome::R,
            LoopyValue::Dud => PartizanOutcome::Draw,
        }
    }

    /// Negation (swap the Left/Right roles): `−on = off`, `−over = under`, and the
    /// self-negating `0`, `∗`, `dud`.
    pub fn neg(&self) -> LoopyValue {
        match self {
            LoopyValue::Zero => LoopyValue::Zero,
            LoopyValue::Star => LoopyValue::Star,
            LoopyValue::On => LoopyValue::Off,
            LoopyValue::Off => LoopyValue::On,
            LoopyValue::Over => LoopyValue::Under,
            LoopyValue::Under => LoopyValue::Over,
            LoopyValue::Dud => LoopyValue::Dud,
        }
    }

    /// Whether this value is a **stopper** (guaranteed to end when played in
    /// isolation). Everything here is a stopper except `dud`.
    pub fn is_stopper(&self) -> bool {
        !matches!(self, LoopyValue::Dud)
    }

    /// The disjunctive sum, where it is defined on this catalogue. Returns `None`
    /// when the sum leaves the catalogue (e.g. `over + over`, a distinct
    /// infinitesimal not named here) — honestly partial, not wrong.
    ///
    /// The closed cases: `dud` absorbs everything (`dud + G = dud`); `on + off =
    /// dud`; `on`/`off` absorb every other stopper (`on` is `>` every stopper);
    /// `∗ + ∗ = 0`; `over + under = 0`; and `0` is the identity.
    pub fn add(&self, other: &LoopyValue) -> Option<LoopyValue> {
        use LoopyValue::*;
        let r = match (*self, *other) {
            (Dud, _) | (_, Dud) => Dud,
            (Zero, x) | (x, Zero) => x,
            (On, On) => On,
            (Off, Off) => Off,
            (On, Off) | (Off, On) => Dud,
            (On, Star) | (Star, On) | (On, Over) | (Over, On) | (On, Under) | (Under, On) => On,
            (Off, Star) | (Star, Off) | (Off, Over) | (Over, Off) | (Off, Under) | (Under, Off) => {
                Off
            }
            (Star, Star) => Zero,
            (Over, Under) | (Under, Over) => Zero,
            // over+over, under+under, star+over, star+under leave the catalogue.
            _ => return None,
        };
        Some(r)
    }
}

impl PartialOrd for LoopyValue {
    /// The conservative partial order on the catalogue. The comparable core is the
    /// chain `off < under < 0 < over < on`, with `on` above and `off` below every
    /// other (non-`dud`) value. `∗` is confused with `0`, `over`, `under` (only
    /// comparable to the extremes `on`/`off`), and `dud` is confused with
    /// everything (comparable only to itself). Incomparable ⇒ `None`.
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        use LoopyValue::*;
        if self == other {
            return Some(Ordering::Equal);
        }
        match (*self, *other) {
            // dud is confused with every other value.
            (Dud, _) | (_, Dud) => None,
            // on is the top, off the bottom (over all non-dud values).
            (On, _) => Some(Ordering::Greater),
            (_, On) => Some(Ordering::Less),
            (Off, _) => Some(Ordering::Less),
            (_, Off) => Some(Ordering::Greater),
            // star is confused with 0/over/under.
            (Star, _) | (_, Star) => None,
            // the remaining comparable chain under < 0 < over.
            (a, b) => {
                let rank = |v: LoopyValue| match v {
                    Under => -1i8,
                    Zero => 0,
                    Over => 1,
                    _ => unreachable!("on/off/star/dud handled above"),
                };
                Some(rank(a).cmp(&rank(b)))
            }
        }
    }
}

// ---------------------------------------------------------------------------
// 2. The graph-level engine.
// ---------------------------------------------------------------------------

/// A loopy game as a finite move graph (`succ[v]` = the positions reachable from
/// `v` in one move). The move graph may be cyclic; outcomes are computed by the
/// retrograde [`kernel::outcomes`](crate::games::outcomes) (Win / Loss / Draw,
/// where **Loss = P-position** and **Draw = the loopy escape**).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LoopyGraph {
    succ: Vec<Vec<usize>>,
}

impl LoopyGraph {
    /// Build from explicit adjacency lists.
    pub fn new(succ: Vec<Vec<usize>>) -> LoopyGraph {
        LoopyGraph { succ }
    }

    /// Build from a move rule on positions `0..n` (the rule may produce cycles).
    pub fn from_rule<F: Fn(usize) -> Vec<usize>>(n: usize, moves: F) -> LoopyGraph {
        LoopyGraph {
            succ: (0..n).map(moves).collect(),
        }
    }

    /// The adjacency lists.
    pub fn succ(&self) -> &[Vec<usize>] {
        &self.succ
    }

    /// Win / Loss / Draw of every position (retrograde analysis).
    pub fn outcomes(&self) -> Vec<Outcome> {
        kernel::outcomes(&self.succ)
    }

    /// The Loss positions = **P-positions** (the player to move loses).
    pub fn loss_set(&self) -> Vec<usize> {
        self.indices_with(Outcome::Loss)
    }

    /// The Win positions = N-positions (the player to move wins).
    pub fn win_set(&self) -> Vec<usize> {
        self.indices_with(Outcome::Win)
    }

    /// The Draw positions — the loopy degree of freedom (neither player can force a
    /// win). Empty iff the game is effectively non-loopy.
    pub fn draw_set(&self) -> Vec<usize> {
        self.indices_with(Outcome::Draw)
    }

    fn indices_with(&self, want: Outcome) -> Vec<usize> {
        self.outcomes()
            .into_iter()
            .enumerate()
            .filter(|(_, o)| *o == want)
            .map(|(i, _)| i)
            .collect()
    }

    /// A coarse reading of a position as a catalogue [`LoopyValue`], via its
    /// impartial outcome only: a **Loss** is `0`, a **Draw** is `dud`. A **Win** is
    /// `None` — its value is a nonzero loopy nimber (use [`loopy_nim_values`]), not
    /// a named catalogue stopper. This is deliberately partial: an impartial move
    /// graph cannot express the Left/Right asymmetry of `on`/`off`/`over`/`under`.
    pub fn classify(&self, v: usize) -> Option<LoopyValue> {
        match self.outcomes().get(v)? {
            Outcome::Loss => Some(LoopyValue::Zero),
            Outcome::Draw => Some(LoopyValue::Dud),
            Outcome::Win => None,
        }
    }
}

// ---------------------------------------------------------------------------
// 3. Impartial loopy nim-values (partial sidling).
// ---------------------------------------------------------------------------

/// A loopy nim-value: an ordinary nimber, or `Side` (the loopy `∞`) for a drawn
/// position.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LoopyNimber {
    /// A genuine nimber (the position terminates under optimal impartial play).
    Value(u128),
    /// The "side" value `∞`: a Draw position, from which play can be sustained
    /// forever.
    Side,
}

/// Certificate for [`loopy_nim_values_certified`]: the outcome split, the positions
/// promoted to `Side`, and whether the bounded sidling solver was needed.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LoopyNimCertificate {
    pub outcomes: Vec<Outcome>,
    pub side_positions: Vec<usize>,
    pub used_sidling_solver: bool,
    pub sidling_assignments_examined: usize,
}

/// Loopy nim-values of an impartial game graph. Draw positions (per
/// [`kernel::outcomes`](crate::games::outcomes)) are `Side`; the rest carry an
/// ordinary nimber `mex`-computed over their non-`Side` options.
///
/// **Exact** when the non-Draw subgraph is acyclic — there `Value(0) ⟺ Loss` and
/// the values agree with [`grundy_graph`](crate::games::grundy_graph). Returns
/// `None` when the non-Draw subgraph still contains a cycle (e.g. a cycle whose
/// members are resolved to Win/Loss by an exit): assigning finite nim-values there
/// requires the full *sidling* fixpoint, which is deferred. This mirrors
/// `grundy_graph` returning `None` on a cycle, refined by first pulling out the
/// draws as `Side`.
pub fn loopy_nim_values(succ: &[Vec<usize>]) -> Option<Vec<LoopyNimber>> {
    loopy_nim_values_certified(succ).map(|(values, _)| values)
}

/// [`loopy_nim_values`] plus a small certificate explaining the outcome split and
/// whether cyclic non-Draw sidling was solved by the bounded mex-equation search.
pub fn loopy_nim_values_certified(
    succ: &[Vec<usize>],
) -> Option<(Vec<LoopyNimber>, LoopyNimCertificate)> {
    let n = succ.len();
    let out = kernel::outcomes(succ);
    let is_side: Vec<bool> = out.iter().map(|o| *o == Outcome::Draw).collect();
    let mut val = vec![0u128; n];
    let mut state = vec![0u8; n]; // 0 unvisited, 1 visiting, 2 done
    let mut needs_sidling = false;

    fn dfs(
        succ: &[Vec<usize>],
        is_side: &[bool],
        v: usize,
        state: &mut [u8],
        val: &mut [u128],
    ) -> Option<()> {
        match state[v] {
            2 => return Some(()),
            1 => return None, // back-edge among non-Side nodes ⇒ defer to full sidling
            _ => {}
        }
        state[v] = 1;
        let mut opts = Vec::new();
        for &w in &succ[v] {
            if is_side[w] {
                continue; // a Side option neither blocks a mex value nor forces a loss
            }
            dfs(succ, is_side, w, state, val)?;
            opts.push(val[w]);
        }
        val[v] = mex(opts);
        state[v] = 2;
        Some(())
    }

    for v in 0..n {
        if !is_side[v] && dfs(succ, &is_side, v, &mut state, &mut val).is_none() {
            needs_sidling = true;
            break;
        }
    }

    let mut assignments = 0usize;
    if needs_sidling {
        let (sidled, count) = solve_mex_sidling(succ, &is_side)?;
        val = sidled;
        assignments = count;
    }

    let values: Vec<LoopyNimber> = (0..n)
        .map(|v| {
            if is_side[v] {
                LoopyNimber::Side
            } else {
                LoopyNimber::Value(val[v])
            }
        })
        .collect();
    let cert = LoopyNimCertificate {
        outcomes: out,
        side_positions: is_side
            .iter()
            .enumerate()
            .filter_map(|(i, &side)| side.then_some(i))
            .collect(),
        used_sidling_solver: needs_sidling,
        sidling_assignments_examined: assignments,
    };
    Some((values, cert))
}

fn solve_mex_sidling(succ: &[Vec<usize>], is_side: &[bool]) -> Option<(Vec<u128>, usize)> {
    let n = succ.len();
    let finite: Vec<usize> = (0..n).filter(|&v| !is_side[v]).collect();
    let mut order = finite.clone();
    order.sort_by_key(|&v| succ[v].iter().filter(|&&w| !is_side[w]).count());
    let mut assigned = vec![false; n];
    for (v, &side) in is_side.iter().enumerate() {
        if side {
            assigned[v] = true;
        }
    }
    let values = vec![0u128; n];
    let max_for: Vec<u128> = (0..n)
        .map(|v| succ[v].iter().filter(|&&w| !is_side[w]).count() as u128)
        .collect();
    let examined = 0usize;

    struct Solver<'a> {
        order: Vec<usize>,
        succ: &'a [Vec<usize>],
        is_side: &'a [bool],
        max_for: Vec<u128>,
        assigned: Vec<bool>,
        values: Vec<u128>,
        examined: usize,
    }

    impl Solver<'_> {
        fn rec(&mut self, idx: usize) -> Option<bool> {
            if self.examined > MAX_SIDLING_ASSIGNMENTS {
                return None;
            }
            if idx == self.order.len() {
                return Some(all_mex_equations_hold(
                    self.succ,
                    self.is_side,
                    &self.values,
                ));
            }
            let v = self.order[idx];
            for candidate in 0..=self.max_for[v] {
                self.examined += 1;
                if self.examined > MAX_SIDLING_ASSIGNMENTS {
                    return None;
                }
                self.values[v] = candidate;
                self.assigned[v] = true;
                if partial_mex_equations_hold(self.succ, self.is_side, &self.assigned, &self.values)
                {
                    match self.rec(idx + 1) {
                        Some(true) => return Some(true),
                        Some(false) => {}
                        None => return None,
                    }
                }
                self.assigned[v] = false;
            }
            Some(false)
        }
    }

    let mut solver = Solver {
        order,
        succ,
        is_side,
        max_for,
        assigned,
        values,
        examined,
    };
    match solver.rec(0) {
        Some(true) => Some((solver.values, solver.examined)),
        Some(false) | None => None,
    }
}

fn partial_mex_equations_hold(
    succ: &[Vec<usize>],
    is_side: &[bool],
    assigned: &[bool],
    values: &[u128],
) -> bool {
    for v in 0..succ.len() {
        if is_side[v] || !assigned[v] {
            continue;
        }
        if succ[v].iter().any(|&w| !is_side[w] && !assigned[w]) {
            continue;
        }
        if values[v] != mex_value(succ, is_side, values, v) {
            return false;
        }
    }
    true
}

fn all_mex_equations_hold(succ: &[Vec<usize>], is_side: &[bool], values: &[u128]) -> bool {
    (0..succ.len())
        .filter(|&v| !is_side[v])
        .all(|v| values[v] == mex_value(succ, is_side, values, v))
}

fn mex_value(succ: &[Vec<usize>], is_side: &[bool], values: &[u128], v: usize) -> u128 {
    mex(succ[v]
        .iter()
        .filter_map(|&w| (!is_side[w]).then_some(values[w])))
}

// ---------------------------------------------------------------------------
// 4. The research instrument: Loss-set AND Draw-set of a cyclic rule.
// ---------------------------------------------------------------------------

/// Given a move rule on positions `0..n` (cycles allowed), return its
/// `(loss_set, draw_set)` — the P-positions and the loopy Draw positions. The
/// acyclic analogue (`examples/interactive_kernel.rs`) discards the Draw count;
/// here both sets are first-class, which is the point: a cyclic rule can carve a
/// non-XOR-linear Draw-set.
pub fn loopy_decision_sets<F: Fn(usize) -> Vec<usize>>(
    n: usize,
    moves: F,
) -> (Vec<usize>, Vec<usize>) {
    let g = LoopyGraph::from_rule(n, moves);
    (g.loss_set(), g.draw_set())
}

/// Probe a cyclic move rule on `F₂^k` (positions `0..2^k`) for a quadric P-set or
/// Draw-set: returns `(loss_fit, draw_fit)`, each the
/// [`fit_f2_quadratic`](crate::forms::fit_f2_quadratic) of the corresponding set
/// (or `None` if that set is not the zero-set of any `F₂` quadratic form). A
/// genuinely-quadratic Draw-set ([`QuadricFit::is_genuinely_quadratic`]) is the
/// Tier-2 target.
pub fn loopy_quadric_probe<F: Fn(usize) -> Vec<usize>>(
    k: usize,
    moves: F,
) -> (Option<QuadricFit>, Option<QuadricFit>) {
    assert!(k <= 20, "loopy_quadric_probe is exponential in k");
    let n = 1usize << k;
    let (loss, draw) = loopy_decision_sets(n, moves);
    let loss_u: Vec<u128> = loss.iter().map(|&v| v as u128).collect();
    let draw_u: Vec<u128> = draw.iter().map(|&v| v as u128).collect();
    (fit_f2_quadratic(&loss_u, k), fit_f2_quadratic(&draw_u, k))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::games::grundy_graph;

    // --- the catalogue ---

    #[test]
    fn negation_is_an_involution_and_swaps_sides() {
        use LoopyValue::*;
        for v in [Zero, Star, On, Off, Over, Under, Dud] {
            assert_eq!(v.neg().neg(), v);
        }
        assert_eq!(On.neg(), Off);
        assert_eq!(Over.neg(), Under);
        assert_eq!(Dud.neg(), Dud);
    }

    #[test]
    fn outcomes_of_the_stoppers() {
        use LoopyValue::*;
        assert_eq!(Zero.outcome(), PartizanOutcome::P);
        assert_eq!(Star.outcome(), PartizanOutcome::N);
        assert_eq!(On.outcome(), PartizanOutcome::L);
        assert_eq!(Off.outcome(), PartizanOutcome::R);
        assert_eq!(Over.outcome(), PartizanOutcome::L);
        assert_eq!(Under.outcome(), PartizanOutcome::R);
        assert_eq!(Dud.outcome(), PartizanOutcome::Draw);
        assert!(!Dud.is_stopper());
        assert!(On.is_stopper());
    }

    #[test]
    fn the_closed_sums() {
        use LoopyValue::*;
        // 0 is the identity.
        for v in [Zero, Star, On, Off, Over, Under, Dud] {
            assert_eq!(Zero.add(&v), Some(v));
        }
        // dud absorbs everything.
        for v in [Zero, Star, On, Off, Over, Under, Dud] {
            assert_eq!(Dud.add(&v), Some(Dud));
            assert_eq!(v.add(&Dud), Some(Dud));
        }
        assert_eq!(On.add(&Off), Some(Dud)); // on + off = dud
        assert_eq!(On.add(&On), Some(On));
        assert_eq!(Off.add(&Off), Some(Off));
        assert_eq!(On.add(&Star), Some(On)); // on absorbs stoppers
        assert_eq!(On.add(&Over), Some(On));
        assert_eq!(Star.add(&Star), Some(Zero));
        assert_eq!(Over.add(&Under), Some(Zero));
        // honestly partial outside the catalogue.
        assert_eq!(Over.add(&Over), None);
        assert_eq!(Star.add(&Over), None);
    }

    #[test]
    fn the_partial_order() {
        use LoopyValue::*;
        // the comparable chain off < under < 0 < over < on.
        assert!(Off < Under && Under < Zero && Zero < Over && Over < On);
        assert!(Off < On);
        // on/off are the extremes (over every non-dud value).
        assert!(On > Star && Off < Star);
        // star is confused with 0/over/under; dud with everything.
        assert_eq!(Star.partial_cmp(&Zero), None);
        assert_eq!(Star.partial_cmp(&Over), None);
        assert_eq!(Dud.partial_cmp(&Zero), None);
        assert_eq!(Dud.partial_cmp(&On), None);
        assert_eq!(Dud.partial_cmp(&Dud), Some(Ordering::Equal));
    }

    // --- the graph engine ---

    #[test]
    fn two_cycle_is_all_draws() {
        let g = LoopyGraph::new(vec![vec![1], vec![0]]);
        assert_eq!(g.outcomes(), vec![Outcome::Draw, Outcome::Draw]);
        assert_eq!(g.draw_set(), vec![0, 1]);
        assert_eq!(g.classify(0), Some(LoopyValue::Dud));
    }

    #[test]
    fn nim_heap_path_has_no_draws() {
        // The Nim heap of size n is the path n → {n-1, …, 0}: only 0 is a Loss.
        let n = 6usize;
        let succ: Vec<Vec<usize>> = (0..=n).map(|h| (0..h).collect()).collect();
        let g = LoopyGraph::new(succ);
        assert_eq!(g.loss_set(), vec![0]);
        assert!(g.draw_set().is_empty());
        assert_eq!(g.classify(0), Some(LoopyValue::Zero));
    }

    // --- loopy nim-values ---

    #[test]
    fn loopy_nim_values_match_grundy_on_acyclic_graphs() {
        // No draws ⇒ the non-Side subgraph is the whole (acyclic) graph.
        let succ = vec![vec![1, 2], vec![3], vec![3], vec![]];
        let lv = loopy_nim_values(&succ).unwrap();
        let g = grundy_graph(&succ).unwrap();
        for v in 0..succ.len() {
            assert_eq!(lv[v], LoopyNimber::Value(g[v]));
        }
    }

    #[test]
    fn draws_are_side_and_value_zero_is_loss() {
        // 0↔1 a drawn 2-cycle; 2→3, 3 terminal (Loss). 2 is a Win (→ Loss 3).
        let succ = vec![vec![1], vec![0], vec![3], vec![]];
        let lv = loopy_nim_values(&succ).unwrap();
        assert_eq!(lv[0], LoopyNimber::Side);
        assert_eq!(lv[1], LoopyNimber::Side);
        assert_eq!(lv[3], LoopyNimber::Value(0)); // terminal ⇒ Loss ⇒ 0
        assert_eq!(lv[2], LoopyNimber::Value(1)); // mex{0} = 1
    }

    #[test]
    fn cyclic_non_draw_subgraph_uses_bounded_sidling() {
        // cycle-with-exit: 0→1, 1→{2,0}, 2 terminal. kernel resolves 0,1 to
        // Loss/Win (non-Draw), and the bounded sidling solver finds the finite mex
        // fixed point g = [0, 1, 0].
        let succ = vec![vec![1], vec![2, 0], vec![]];
        let (values, cert) = loopy_nim_values_certified(&succ).unwrap();
        assert_eq!(
            values,
            vec![
                LoopyNimber::Value(0),
                LoopyNimber::Value(1),
                LoopyNimber::Value(0)
            ]
        );
        assert!(cert.used_sidling_solver);
        assert!(cert.sidling_assignments_examined > 0);
        // but the outcome analysis is still exact.
        let g = LoopyGraph::new(succ);
        assert_eq!(
            g.outcomes(),
            vec![Outcome::Loss, Outcome::Win, Outcome::Loss]
        );
    }

    // --- the research instrument ---

    #[test]
    fn decision_sets_recover_an_acyclic_loss_set_with_no_draws() {
        // A downward (terminating) rule: move v → any w < v. Then 0 is the only
        // Loss and there are no Draws — matching the acyclic interactive probe.
        let n = 8;
        let (loss, draw) = loopy_decision_sets(n, |v| (0..v).collect());
        assert_eq!(loss, vec![0]);
        assert!(draw.is_empty());
    }

    #[test]
    fn quadric_probe_reads_both_sets() {
        // A cyclic rule on F₂² that makes {0} a Loss and pairs the rest into a draw
        // cycle — exercising both fit slots. Here we just check the plumbing: the
        // loss-set fits (a point) and the call returns without panicking.
        let (loss_fit, _draw_fit) = loopy_quadric_probe(2, |v| {
            if v == 0 {
                vec![] // terminal ⇒ Loss
            } else {
                vec![0] // everyone moves to 0 ⇒ Win, no draws
            }
        });
        // {0} as a P-set over F₂² is the anisotropic quadric (Arf 1).
        let f = loss_fit.expect("{0} is a quadric");
        assert!(f.is_genuinely_quadratic());
        assert_eq!(f.arf.arf, 1);
    }
}
