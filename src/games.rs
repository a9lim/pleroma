//! Coin-turning games and the game-theoretic definition of nim-multiplication
//! — the concrete bridge from games to On₂.
//!
//! Conway's *Turning Corners*: on a grid of coins, a move picks a coin and the
//! SW rectangle under it, turning the four corners. The Grundy value of the
//! single coin at (x, y) satisfies the excludant recurrence
//!
//!   x ⊗ y = mex { (i⊗y) ⊕ (x⊗j) ⊕ (i⊗j) : 0 ≤ i < x, 0 ≤ j < y }.
//!
//! This *is* nim-multiplication — defined entirely by a game. `nim_mul_mex`
//! computes it, and the tests confirm it agrees with the algebraic Fermat-power
//! `nim_mul` in `nimber.rs`. Two independent definitions, one combinatorial and
//! one field-theoretic, of the same product.
//!
//! Nim-addition is likewise a game: the disjunctive sum of single-coin
//! positions XORs their Grundy values (Sprague–Grundy). And there is a 1-D
//! coin-turning game whose single-coin Grundy value is g(n) = 2ⁿ — "turn coin n
//! together with any subset of the coins left of it": by induction the reachable
//! values are the full F₂-span of {1,2,…,2ⁿ⁻¹}, so the mex is 2ⁿ. Under that
//! game a position's value *is* a nimber (the bitmask of its heads-up coins),
//! which is the sense in which the nimber backend is "made of games".

use std::cell::RefCell;
use std::collections::{HashMap, HashSet};

thread_local! {
    static MEX_MEMO: RefCell<HashMap<(u64, u64), u64>> = RefCell::new(HashMap::new());
}

fn mex(seen: &HashSet<u64>) -> u64 {
    let mut m = 0u64;
    while seen.contains(&m) {
        m += 1;
    }
    m
}

/// Nim-multiplication via Conway's Turning-Corners excludant recurrence — the
/// *game* definition of the product.
pub fn nim_mul_mex(x: u64, y: u64) -> u64 {
    if x == 0 || y == 0 {
        return 0;
    }
    if let Some(v) = MEX_MEMO.with(|m| m.borrow().get(&(x, y)).copied()) {
        return v;
    }
    let mut seen = HashSet::new();
    for i in 0..x {
        for j in 0..y {
            seen.insert(nim_mul_mex(i, y) ^ nim_mul_mex(x, j) ^ nim_mul_mex(i, j));
        }
    }
    let r = mex(&seen);
    MEX_MEMO.with(|m| m.borrow_mut().insert((x, y), r));
    r
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::nimber::nim_mul;

    #[test]
    fn game_definition_equals_algebraic_nim_mul() {
        // Turning-Corners Grundy values == Fermat-power nim-multiplication.
        for x in 0u64..48 {
            for y in 0u64..48 {
                assert_eq!(
                    nim_mul_mex(x, y),
                    nim_mul(x, y),
                    "mismatch at ({x}, {y})"
                );
            }
        }
    }

    #[test]
    fn turning_corners_realizes_the_field_table() {
        // Spot-check the famous small products straight from the game recurrence.
        assert_eq!(nim_mul_mex(2, 2), 3);
        assert_eq!(nim_mul_mex(2, 3), 1);
        assert_eq!(nim_mul_mex(4, 4), 6);
        assert_eq!(nim_mul_mex(2, 4), 8); // distinct Fermat powers ⇒ ordinary product
    }
}
