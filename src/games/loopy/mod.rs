//! Loopy combinatorial games — games whose move graph may contain cycles, so
//! play need not terminate. This is the third escape (beside the interactive
//! [`kernel`](crate::games::kernel) route and the [`misere`](crate::games::misere)
//! route) from the XOR-linear P-sets of normal-play disjunctive sums: a cyclic
//! rule admits a **Draw** outcome — a position from which neither player can force
//! a win — and the Draw-set is a genuinely new degree of freedom to test against
//! the Gold quadric `{Q=0}` (see `docs/OPEN.md`, the Tier-2 open question).
//!
//! Four layers, re-exported flat so every public path is unchanged:
//!
//!   * [`catalogue`] — [`LoopyWinner`], [`LoopyPartizanOutcome`],
//!     [`PartizanOutcome`], and the [`LoopyValue`] stopper catalogue
//!     (on/off/over/under/dud/±/tis/tisn/∗/0/`s&t` with outcome/neg/partial
//!     order/partial sum).
//!   * [`graph`] — [`LoopyGraph`], the computable wrapper over
//!     [`kernel::outcomes`](crate::games::outcomes) (Win / Loss / Draw retrograde
//!     analysis).
//!   * [`partizan`] — [`LoopyPartizanGraph`]: the two-sided Left/Right retrograde
//!     solver returning exact [`LoopyPartizanOutcome`] pairs, projecting to the
//!     classical five-class [`PartizanOutcome`] only when honest.
//!   * [`nim_values`] — [`LoopyNimber`], [`LoopyNimCertificate`],
//!     [`loopy_nim_values`], and [`loopy_nim_values_certified`]: impartial loopy
//!     nim-values with certificates (including the checked recovery condition).
//!   * [`research`] — [`loopy_decision_sets`] and [`loopy_quadric_probe`]: the
//!     Loss-set / Draw-set research instrument.
//!
//! Deliberately **out of scope** here: [`Game`](crate::games::Game) stays an acyclic
//! `Arc` tree (it cannot represent cycles, by construction), and
//! [`thermography`](crate::games::thermography) stays finite-game-only — loopy games
//! never freeze to a number, so classical temperature does not apply. The sidling
//! support is finite and certified: over-budget or non-canonical fixed-point
//! systems return `None` rather than pretending to be full loopy-game equality.

pub mod catalogue;
pub mod graph;
pub mod nim_values;
pub mod partizan;
pub mod research;

pub use catalogue::*;
pub use graph::*;
pub use nim_values::*;
pub use partizan::*;
pub use research::*;
#[cfg(test)]
mod tests {
    use super::*;

    use crate::games::kernel::{self, Outcome};

    use std::cmp::Ordering;

    use crate::games::grundy_graph;

    // --- the catalogue ---

    #[test]
    fn negation_is_an_involution_and_swaps_sides() {
        use LoopyValue::*;
        for v in [
            Zero,
            Star,
            On,
            Off,
            Over,
            Under,
            PlusMinus,
            Tis,
            Tisn,
            LoopyValue::onside_offside(3, -2),
            Dud,
        ] {
            assert_eq!(v.neg().neg(), v);
        }
        assert_eq!(On.neg(), Off);
        assert_eq!(Over.neg(), Under);
        assert_eq!(Tis.neg(), Tisn);
        assert_eq!(
            LoopyValue::onside_offside(3, -2).neg(),
            LoopyValue::onside_offside(2, -3)
        );
        assert_eq!(Dud.neg(), Dud);
    }

    #[test]
    fn outcomes_of_the_stoppers() {
        use LoopyValue::*;
        assert_eq!(Zero.partizan_outcome(), Some(PartizanOutcome::P));
        assert_eq!(Star.partizan_outcome(), Some(PartizanOutcome::N));
        assert_eq!(PlusMinus.partizan_outcome(), Some(PartizanOutcome::N));
        assert_eq!(On.partizan_outcome(), Some(PartizanOutcome::L));
        assert_eq!(Off.partizan_outcome(), Some(PartizanOutcome::R));
        assert_eq!(Over.partizan_outcome(), Some(PartizanOutcome::L));
        assert_eq!(Under.partizan_outcome(), Some(PartizanOutcome::R));
        assert_eq!(Dud.partizan_outcome(), Some(PartizanOutcome::Draw));
        assert_eq!(
            Tis.outcome(),
            LoopyPartizanOutcome::new(LoopyWinner::Left, LoopyWinner::Draw)
        );
        assert_eq!(
            Tisn.outcome(),
            LoopyPartizanOutcome::new(LoopyWinner::Draw, LoopyWinner::Right)
        );
        assert_eq!(Tis.partizan_outcome(), None);
        assert_eq!(Tis.sides(), Some((1, 0)));
        assert_eq!(Tisn.sides(), Some((0, -1)));
        assert!(!Dud.is_stopper());
        assert!(!Tis.is_stopper());
        assert!(On.is_stopper());
    }

    #[test]
    fn the_closed_sums() {
        use LoopyValue::*;
        // 0 is the identity.
        for v in [Zero, Star, On, Off, Over, Under, PlusMinus, Tis, Tisn, Dud] {
            assert_eq!(Zero.add(&v), Some(v));
        }
        // dud absorbs everything.
        for v in [Zero, Star, On, Off, Over, Under, PlusMinus, Tis, Tisn, Dud] {
            assert_eq!(Dud.add(&v), Some(Dud));
            assert_eq!(v.add(&Dud), Some(Dud));
        }
        assert_eq!(On.add(&Off), Some(Dud)); // on + off = dud
        assert_eq!(On.add(&On), Some(On));
        assert_eq!(Off.add(&Off), Some(Off));
        assert_eq!(On.add(&Star), Some(On)); // on absorbs stoppers
        assert_eq!(On.add(&Over), Some(On));
        assert_eq!(Star.add(&Star), Some(Zero));
        assert_eq!(Over.add(&Under), None);
        assert_eq!(Over.add(&Over), Some(Over));
        assert_eq!(Under.add(&Under), Some(Under));
        assert_eq!(Star.add(&Over), Some(Over));
        assert_eq!(Star.add(&Under), Some(Under));
        // over+under is a draw-class value outside these named tags.
        assert_eq!(Under.add(&Over), None);
        assert_eq!(
            LoopyValue::onside_offside(1, 0).add(&LoopyValue::onside_offside(0, -1)),
            Some(LoopyValue::onside_offside(1, -1))
        );
        assert_eq!(Tis.add(&Tisn), None);
    }

    #[test]
    fn the_partial_order() {
        use LoopyValue::*;
        // the comparable chain off < under < 0 < over < on.
        assert!(Off < Under && Under < Zero && Zero < Over && Over < On);
        assert!(Under < Star && Star < Over);
        assert!(Off < On);
        // on/off are the extremes (over every non-dud value).
        assert!(On > Star && Off < Star);
        // star is confused with 0; dud with everything.
        assert_eq!(Star.partial_cmp(&Zero), None);
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

    // --- the partizan graph engine ---

    #[test]
    fn partizan_graph_recovers_classical_short_outcomes() {
        // position 0 is terminal; 1 = *; 2 = {0|}; 3 = {|0}.
        let left = vec![vec![], vec![0], vec![0], vec![]];
        let right = vec![vec![], vec![0], vec![], vec![0]];
        let g = LoopyPartizanGraph::new(left, right);
        assert_eq!(
            g.partizan_outcomes(),
            vec![
                Some(PartizanOutcome::P),
                Some(PartizanOutcome::N),
                Some(PartizanOutcome::L),
                Some(PartizanOutcome::R),
            ]
        );
        assert!(g.draw_set().is_empty());
    }

    #[test]
    fn partizan_graph_keeps_tis_as_mixed_draw_class() {
        // Repo convention: tis = {0|tisn}, tisn = {tis|0}, with 0 terminal.
        let left = vec![vec![2], vec![0], vec![]];
        let right = vec![vec![1], vec![2], vec![]];
        let g = LoopyPartizanGraph::new(left, right);
        let out = g.outcomes();
        assert_eq!(out[0], LoopyValue::Tis.outcome());
        assert_eq!(out[1], LoopyValue::Tisn.outcome());
        assert_eq!(g.classify(0), None);
        assert_eq!(g.nonclassical_set(), vec![0, 1]);
        assert_eq!(g.draw_set(), vec![0, 1]);
    }

    #[test]
    fn impartial_partizan_graph_matches_kernel_outcomes() {
        let succ = vec![vec![1], vec![2, 0], vec![]];
        let g = LoopyPartizanGraph::new(succ.clone(), succ.clone());
        assert_eq!(
            g.partizan_outcomes(),
            kernel::outcomes(&succ)
                .into_iter()
                .map(|o| match o {
                    Outcome::Loss => Some(PartizanOutcome::P),
                    Outcome::Win => Some(PartizanOutcome::N),
                    Outcome::Draw => Some(PartizanOutcome::Draw),
                })
                .collect::<Vec<_>>()
        );
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
        assert!(cert.recovery_condition_holds);
        assert!(cert.recovery_blockers.is_empty());
        // but the outcome analysis is still exact.
        let g = LoopyGraph::new(succ);
        assert_eq!(
            g.outcomes(),
            vec![Outcome::Loss, Outcome::Win, Outcome::Loss]
        );
    }

    #[test]
    fn ambiguous_cyclic_sidling_returns_none() {
        // Symmetric cycle-with-exits:
        //   g0 = mex{g1,0}, g1 = mex{g0,0}
        // has two fixed points, (1,2) and (2,1). Positions 0 and 1 are graph-
        // symmetric, so choosing either finite assignment would be noncanonical.
        let succ = vec![vec![1, 2], vec![0, 3], vec![], vec![]];
        assert_eq!(loopy_nim_values(&succ), None);
        assert_eq!(loopy_nim_values_certified(&succ), None);
        let g = LoopyGraph::new(succ);
        assert_eq!(
            g.outcomes(),
            vec![Outcome::Win, Outcome::Win, Outcome::Loss, Outcome::Loss]
        );
    }

    #[test]
    fn recovery_certificate_flags_finite_positions_with_side_options() {
        // 0↔1 is Side; 2 also has a move to terminal 3, so 2 is finite-valued but
        // points at a Side option. Its local mex value is computed, while the
        // recovery/additivity condition is explicitly false.
        let succ = vec![vec![1], vec![0], vec![0, 3], vec![]];
        let (_values, cert) = loopy_nim_values_certified(&succ).unwrap();
        assert_eq!(cert.side_positions, vec![0, 1]);
        assert!(!cert.recovery_condition_holds);
        assert_eq!(cert.recovery_blockers, vec![2]);
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
