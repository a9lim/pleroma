//! Outcomes of a finite impartial game graph — the instrument for the
//! "interactive" route to the open question.
//!
//! Normal-play disjunctive sums give XOR-linear P-sets (subspaces); the escape is
//! an *interactive* game whose move graph is not a disjunctive sum. For any such
//! game on a finite position set, this computes the normal-play outcome of every
//! position by retrograde analysis (Win / Loss / Draw), where the **Loss
//! positions are the P-positions** (the player to move loses). With that we can
//! take any candidate move rule — e.g. one coupled through a quadratic form's
//! polar — and ask whether its P-set is the quadric `{Q=0}`.
//!
//! Convention (normal play): a position with no moves is a **Loss** (the player to
//! move cannot move and loses). A position is a **Win** if some move leads to a
//! Loss; a **Loss** if every move leads to a Win; otherwise (only reachable in a
//! cyclic graph) a **Draw** — the player can avoid losing forever.

use std::collections::VecDeque;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Outcome {
    /// The player to move loses — a **P-position**.
    Loss,
    /// The player to move wins — an **N-position**.
    Win,
    /// Neither player can force a win (loopy / cyclic game).
    Draw,
}

/// Normal-play outcomes of a finite game graph given as adjacency lists
/// (`succ[v]` = the positions reachable from `v` in one move). Retrograde
/// analysis: O(V + E).
pub fn outcomes(succ: &[Vec<usize>]) -> Vec<Outcome> {
    let n = succ.len();
    let mut pred = vec![Vec::new(); n];
    for (u, outs) in succ.iter().enumerate() {
        for &v in outs {
            pred[v].push(u);
        }
    }
    let mut remaining: Vec<usize> = succ.iter().map(|o| o.len()).collect();
    let mut label: Vec<Option<Outcome>> = vec![None; n];
    let mut queue: VecDeque<usize> = VecDeque::new();

    // terminal positions are losses (can't move ⇒ lose)
    for v in 0..n {
        if succ[v].is_empty() {
            label[v] = Some(Outcome::Loss);
            queue.push_back(v);
        }
    }
    while let Some(v) = queue.pop_front() {
        let lv = label[v].unwrap();
        for &u in &pred[v] {
            if label[u].is_some() {
                continue;
            }
            match lv {
                // u can move to a Loss ⇒ u is a Win
                Outcome::Loss => {
                    label[u] = Some(Outcome::Win);
                    queue.push_back(u);
                }
                // every move of u now accounted Win; if all are, u is a Loss
                Outcome::Win => {
                    remaining[u] -= 1;
                    if remaining[u] == 0 {
                        label[u] = Some(Outcome::Loss);
                        queue.push_back(u);
                    }
                }
                Outcome::Draw => {}
            }
        }
    }
    label
        .into_iter()
        .map(|l| l.unwrap_or(Outcome::Draw))
        .collect()
}

/// The P-positions (Loss positions) of a game graph, as node indices.
pub fn p_positions(succ: &[Vec<usize>]) -> Vec<usize> {
    outcomes(succ)
        .into_iter()
        .enumerate()
        .filter(|(_, o)| *o == Outcome::Loss)
        .map(|(i, _)| i)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn single_edge_and_terminal() {
        // 0 -> 1, 1 terminal: 1 is Loss (can't move), 0 is Win.
        let succ = vec![vec![1], vec![]];
        assert_eq!(outcomes(&succ), vec![Outcome::Win, Outcome::Loss]);
    }

    #[test]
    fn two_cycle_is_a_draw() {
        // 0 <-> 1 with no exit: neither can force a win ⇒ both Draw.
        let succ = vec![vec![1], vec![0]];
        assert_eq!(outcomes(&succ), vec![Outcome::Draw, Outcome::Draw]);
    }

    #[test]
    fn nim_heap_as_a_path_matches_normal_play() {
        // A single Nim heap of size n is the path n -> n-1 -> … -> 0. Normal play:
        // only the empty heap (0) is a Loss (P); every n≥1 is a Win.
        let n = 6;
        let succ: Vec<Vec<usize>> = (0..=n).map(|h| (0..h).collect()).collect();
        let out = outcomes(&succ);
        assert_eq!(out[0], Outcome::Loss);
        assert!((1..=n).all(|h| out[h] == Outcome::Win));
    }

    #[test]
    fn cycle_with_exit_resolves() {
        // 0 -> 1 -> 2 (terminal), and 1 -> 0 (back edge). 2 Loss; 1 Win (→2);
        // 0 Win (→1? no, 1 is Win) — 0's only move is to 1 (Win) ⇒ 0 is Loss.
        let succ = vec![vec![1], vec![2, 0], vec![]];
        assert_eq!(
            outcomes(&succ),
            vec![Outcome::Loss, Outcome::Win, Outcome::Loss]
        );
    }
}
