//! Impartial loopy nim-values: [`LoopyNimber`], [`LoopyNimCertificate`],
//! [`loopy_nim_values`], and [`loopy_nim_values_certified`].

use crate::games::grundy::mex;
use crate::games::kernel::{self, Outcome};

const MAX_SIDLING_ASSIGNMENTS: usize = 200_000;

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
/// promoted to `Side`, whether the bounded sidling solver was needed, and the
/// checked recovery condition for additive finite-nimber claims.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LoopyNimCertificate {
    pub outcomes: Vec<Outcome>,
    pub side_positions: Vec<usize>,
    pub used_sidling_solver: bool,
    pub sidling_assignments_examined: usize,
    /// True when every finite-valued position has only finite-valued options.
    /// Under this checked condition the emitted finite nimbers are ordinary
    /// Sprague-Grundy labels on a closed subgame, so additivity claims are local
    /// checked facts instead of prose caveats.
    pub recovery_condition_holds: bool,
    /// Finite-valued positions with at least one `Side` option. These are exactly
    /// the blockers for the checked recovery condition above.
    pub recovery_blockers: Vec<usize>,
}

/// Loopy nim-values of an impartial game graph. Draw positions (per
/// [`kernel::outcomes`](crate::games::outcomes)) are `Side`; the rest carry an
/// ordinary nimber `mex`-computed over their non-`Side` options.
///
/// **Exact** when the non-Draw subgraph is acyclic — there `Value(0) ⟺ Loss` and
/// the values agree with [`grundy_graph`](crate::games::grundy_graph). If the
/// non-Draw subgraph is cyclic, a bounded sidling search is accepted only when the
/// finite mex equations have a **unique** solution; ambiguous or over-budget
/// cyclic systems return `None` rather than choosing an order-dependent value.
///
/// **Recovery check**: when a position has Draw (Side) options the emitted
/// `Value(k)` is the Grundy value of the Draw-deleted subgraph at that vertex.
/// The certificate records a checked finite recovery condition:
/// `recovery_condition_holds` iff all finite-valued positions have only
/// finite-valued successors. Only under that condition should additivity-over-sums
/// be cited for the finite nimbers. The `Side` values themselves have no additive
/// nimber arithmetic.
pub fn loopy_nim_values(succ: &[Vec<usize>]) -> Option<Vec<LoopyNimber>> {
    loopy_nim_values_certified(succ).map(|(values, _)| values)
}

/// [`loopy_nim_values`] plus a small certificate explaining the outcome split and
/// whether cyclic non-Draw sidling was solved uniquely by the bounded mex-equation
/// search.
pub fn loopy_nim_values_certified(
    succ: &[Vec<usize>],
) -> Option<(Vec<LoopyNimber>, LoopyNimCertificate)> {
    let n = succ.len();
    let out = kernel::outcomes(succ);
    let is_side: Vec<bool> = out.iter().map(|o| *o == Outcome::Draw).collect();
    let mut val = vec![0u128; n];
    let mut state = vec![0u128; n]; // 0 unvisited, 1 visiting, 2 done
    let mut needs_sidling = false;

    fn dfs(
        succ: &[Vec<usize>],
        is_side: &[bool],
        v: usize,
        state: &mut [u128],
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
    let recovery_blockers: Vec<usize> = (0..n)
        .filter(|&v| !is_side[v] && succ[v].iter().any(|&w| is_side[w]))
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
        recovery_condition_holds: recovery_blockers.is_empty(),
        recovery_blockers,
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
        fn rec(&mut self, idx: usize, solution: &mut Option<Vec<u128>>) -> Option<bool> {
            if self.examined > MAX_SIDLING_ASSIGNMENTS {
                return None;
            }
            if idx == self.order.len() {
                if all_mex_equations_hold(self.succ, self.is_side, &self.values) {
                    if solution.is_some() {
                        return Some(false); // multiple fixed points: not canonical
                    }
                    *solution = Some(self.values.clone());
                }
                return Some(true);
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
                    match self.rec(idx + 1, solution) {
                        Some(true) => {}
                        Some(false) => return Some(false),
                        None => return None,
                    }
                }
                self.assigned[v] = false;
            }
            Some(true)
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
    let mut solution = None;
    match solver.rec(0, &mut solution) {
        Some(true) => solution.map(|values| (values, solver.examined)),
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
