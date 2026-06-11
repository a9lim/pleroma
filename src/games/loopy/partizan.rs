//! The two-sided partizan loopy graph engine: [`LoopyPartizanGraph`] and its
//! exact-outcome solver.

use std::collections::VecDeque;

use super::catalogue::{LoopyPartizanOutcome, LoopyWinner, PartizanOutcome};

/// A finite loopy partizan game graph. `left[v]` are Left's legal moves from
/// position `v`; `right[v]` are Right's legal moves. Cycles are allowed.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LoopyPartizanGraph {
    left: Vec<Vec<usize>>,
    right: Vec<Vec<usize>>,
}

impl LoopyPartizanGraph {
    /// Build from explicit Left and Right adjacency lists.
    pub fn new(left: Vec<Vec<usize>>, right: Vec<Vec<usize>>) -> LoopyPartizanGraph {
        assert_eq!(
            left.len(),
            right.len(),
            "left/right move tables must have the same number of positions"
        );
        LoopyPartizanGraph { left, right }
    }

    /// Build from move rules on positions `0..n`.
    pub fn from_rules<L, R>(n: usize, left_moves: L, right_moves: R) -> LoopyPartizanGraph
    where
        L: Fn(usize) -> Vec<usize>,
        R: Fn(usize) -> Vec<usize>,
    {
        LoopyPartizanGraph {
            left: (0..n).map(left_moves).collect(),
            right: (0..n).map(right_moves).collect(),
        }
    }

    /// Left's adjacency lists.
    pub fn left(&self) -> &[Vec<usize>] {
        &self.left
    }

    /// Right's adjacency lists.
    pub fn right(&self) -> &[Vec<usize>] {
        &self.right
    }

    /// Exact two-sided loopy-partizan outcome of every position.
    pub fn outcomes(&self) -> Vec<LoopyPartizanOutcome> {
        solve_partizan_outcomes(&self.left, &self.right)
    }

    /// Classical partizan outcome classes where the exact two-sided outcome lies
    /// in the five-class image. Mixed loopy starter pairs (`tis`, `tisn`, …)
    /// return `None`.
    pub fn partizan_outcomes(&self) -> Vec<Option<PartizanOutcome>> {
        self.outcomes()
            .into_iter()
            .map(|o| o.partizan_class())
            .collect()
    }

    /// The classical class of position `v`, if it has one.
    pub fn classify(&self, v: usize) -> Option<PartizanOutcome> {
        self.outcomes().get(v).and_then(|o| o.partizan_class())
    }

    /// Positions whose exact starter pair contains a draw for at least one player
    /// to move.
    pub fn draw_set(&self) -> Vec<usize> {
        self.outcomes()
            .into_iter()
            .enumerate()
            .filter_map(|(i, o)| o.has_draw().then_some(i))
            .collect()
    }

    /// Positions whose exact outcome is outside the classical five classes.
    pub fn nonclassical_set(&self) -> Vec<usize> {
        self.outcomes()
            .into_iter()
            .enumerate()
            .filter_map(|(i, o)| o.partizan_class().is_none().then_some(i))
            .collect()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Turn {
    Left,
    Right,
}

fn state(v: usize, turn: Turn) -> usize {
    2 * v
        + match turn {
            Turn::Left => 0,
            Turn::Right => 1,
        }
}

fn state_parts(s: usize) -> (usize, Turn) {
    (s / 2, if s & 1 == 0 { Turn::Left } else { Turn::Right })
}

fn owner_winner(turn: Turn) -> LoopyWinner {
    match turn {
        Turn::Left => LoopyWinner::Left,
        Turn::Right => LoopyWinner::Right,
    }
}

fn opponent_winner(turn: Turn) -> LoopyWinner {
    match turn {
        Turn::Left => LoopyWinner::Right,
        Turn::Right => LoopyWinner::Left,
    }
}

fn solve_partizan_outcomes(left: &[Vec<usize>], right: &[Vec<usize>]) -> Vec<LoopyPartizanOutcome> {
    assert_eq!(
        left.len(),
        right.len(),
        "left/right move tables must have the same number of positions"
    );
    let n = left.len();
    let states = 2 * n;
    let mut succ = vec![Vec::new(); states];
    let mut pred = vec![Vec::new(); states];
    for v in 0..n {
        for &w in &left[v] {
            let s = state(v, Turn::Left);
            let t = state(w, Turn::Right);
            succ[s].push(t);
            pred[t].push(s);
        }
        for &w in &right[v] {
            let s = state(v, Turn::Right);
            let t = state(w, Turn::Left);
            succ[s].push(t);
            pred[t].push(s);
        }
    }

    let mut remaining: Vec<usize> = succ.iter().map(Vec::len).collect();
    let mut label: Vec<Option<LoopyWinner>> = vec![None; states];
    let mut queue = VecDeque::new();

    for s in 0..states {
        if succ[s].is_empty() {
            let (_, turn) = state_parts(s);
            label[s] = Some(opponent_winner(turn));
            queue.push_back(s);
        }
    }

    while let Some(s) = queue.pop_front() {
        let winner = label[s].unwrap();
        for &p in &pred[s] {
            if label[p].is_some() {
                continue;
            }
            let (_, turn) = state_parts(p);
            if winner == owner_winner(turn) {
                label[p] = Some(winner);
                queue.push_back(p);
            } else {
                remaining[p] -= 1;
                if remaining[p] == 0 {
                    label[p] = Some(winner);
                    queue.push_back(p);
                }
            }
        }
    }

    (0..n)
        .map(|v| {
            LoopyPartizanOutcome::new(
                label[state(v, Turn::Left)].unwrap_or(LoopyWinner::Draw),
                label[state(v, Turn::Right)].unwrap_or(LoopyWinner::Draw),
            )
        })
        .collect()
}
