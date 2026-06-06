//! Sprague–Grundy theory: the normal-play impartial center.
//!
//! Every finite impartial game (normal play) has a **Grundy value** — a nimber
//! — given by the *mex* (minimal excludant) of its options' Grundy values:
//!
//!   g(G) = mex { g(G') : G' a move from G }.
//!
//! A position is a P-position (the player to move loses) iff `g(G) = 0`, so this
//! is the impartial refinement of [`kernel::outcomes`](crate::games::outcomes):
//! `outcomes` returns a Win/Loss/Draw *bit*, the Grundy value the full nimber.
//!
//! The **Sprague–Grundy theorem** is then `g(G + H) = g(G) ⊕ g(H)`: under
//! disjunctive sum every impartial game behaves as a single Nim heap of size
//! `g(G)`. This module computes Grundy values two ways — over an explicit finite
//! game graph ([`grundy_graph`], the mirror of `kernel::outcomes`) and over a
//! move-generating closure ([`grundy`], the normal-play mirror of
//! [`misere::misere_is_n`](crate::games::misere_is_n)) — and the tests pin the
//! theorem against Bouton's XOR for Nim.

use std::collections::{HashMap, HashSet};
use std::hash::Hash;

/// The minimal excludant of a set of nimbers: the least non-negative integer not
/// present.
pub fn mex<I: IntoIterator<Item = u128>>(values: I) -> u128 {
    let seen: HashSet<u128> = values.into_iter().collect();
    let mut m = 0u128;
    while seen.contains(&m) {
        m += 1;
    }
    m
}

fn grundy_dfs(succ: &[Vec<usize>], v: usize, state: &mut [u8], g: &mut [u128]) -> Option<()> {
    match state[v] {
        2 => return Some(()),
        1 => return None, // back-edge ⇒ a cycle ⇒ Grundy value undefined
        _ => {}
    }
    state[v] = 1;
    for &w in &succ[v] {
        grundy_dfs(succ, w, state, g)?;
    }
    g[v] = mex(succ[v].iter().map(|&w| g[w]));
    state[v] = 2;
    Some(())
}

/// Grundy values of a finite **acyclic** impartial game graph, given as
/// adjacency lists (`succ[v]` = the positions reachable from `v` in one move).
/// Returns `None` if the graph has a cycle — Grundy values are only defined on
/// terminating (loopy-free) games; use [`outcomes`](crate::games::outcomes) for
/// the Win/Loss/Draw analysis of loopy games.
///
/// Position `v` is a P-position (Loss) iff `result[v] == 0`.
pub fn grundy_graph(succ: &[Vec<usize>]) -> Option<Vec<u128>> {
    let n = succ.len();
    let mut state = vec![0u8; n];
    let mut g = vec![0u128; n];
    for v in 0..n {
        grundy_dfs(succ, v, &mut state, &mut g)?;
    }
    Some(g)
}

/// The Grundy value of a position given a move generator, memoised. `moves(pos)`
/// returns the positions reachable in one move; the game must be finite and
/// acyclic (no position reachable from itself), exactly as for
/// [`misere::misere_is_n`](crate::games::misere_is_n). This is the way to compute
/// Grundy values of games defined by a rule (Nim, octal games, …) without
/// materialising the whole graph.
pub fn grundy<P, F>(pos: &P, moves: &F, memo: &mut HashMap<P, u128>) -> u128
where
    P: Eq + Hash + Clone,
    F: Fn(&P) -> Vec<P>,
{
    if let Some(&v) = memo.get(pos) {
        return v;
    }
    let options = moves(pos);
    let values: Vec<u128> = options.iter().map(|p| grundy(p, moves, memo)).collect();
    let g = mex(values);
    memo.insert(pos.clone(), g);
    g
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::games::{outcomes, Outcome};

    #[test]
    fn mex_basics() {
        assert_eq!(mex([]), 0);
        assert_eq!(mex([0, 1, 2]), 3);
        assert_eq!(mex([1, 2, 3]), 0);
        assert_eq!(mex([0, 2, 3]), 1);
        assert_eq!(mex([0, 0, 1, 1]), 2); // duplicates ignored
    }

    #[test]
    fn single_nim_heap_grundy_is_its_size() {
        // The heap of size h is the path h → {h-1, …, 0}; g(h) = mex{0..h} = h.
        let n = 8;
        let succ: Vec<Vec<usize>> = (0..=n).map(|h| (0..h).collect()).collect();
        let g = grundy_graph(&succ).unwrap();
        assert_eq!(g, (0..=n as u128).collect::<Vec<_>>());
    }

    #[test]
    fn grundy_zero_iff_loss() {
        // Cross-check against kernel::outcomes on an acyclic graph: Grundy 0 ⟺ Loss.
        let succ = vec![vec![1, 2], vec![3], vec![3], vec![]];
        let g = grundy_graph(&succ).unwrap();
        let o = outcomes(&succ);
        for v in 0..succ.len() {
            assert_eq!(g[v] == 0, o[v] == Outcome::Loss);
        }
    }

    #[test]
    fn cyclic_graph_has_no_grundy_value() {
        assert_eq!(grundy_graph(&[vec![1], vec![0]]), None);
    }

    #[test]
    fn sprague_grundy_theorem_is_boutons_xor() {
        // Nim: a position is a multiset of heap sizes; a move shrinks one heap.
        // The Sprague–Grundy theorem says g(position) = XOR of the heap sizes.
        fn nim_moves(heaps: &Vec<u128>) -> Vec<Vec<u128>> {
            let mut out = Vec::new();
            for (i, &h) in heaps.iter().enumerate() {
                for new in 0..h {
                    let mut next = heaps.clone();
                    next[i] = new;
                    next.retain(|&x| x != 0);
                    next.sort_unstable();
                    out.push(next);
                }
            }
            out
        }
        let mut memo = HashMap::new();
        for a in 0..=4u128 {
            for b in 0..=4u128 {
                for c in 0..=4u128 {
                    let mut heaps: Vec<u128> = vec![a, b, c];
                    heaps.retain(|&x| x != 0);
                    heaps.sort_unstable();
                    let g = grundy(&heaps, &nim_moves, &mut memo);
                    assert_eq!(g, a ^ b ^ c, "heaps {a},{b},{c}");
                }
            }
        }
    }
}
