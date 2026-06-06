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
    static MEX_MEMO: RefCell<HashMap<(u128, u128), u128>> = RefCell::new(HashMap::new());
}

fn mex(seen: &HashSet<u128>) -> u128 {
    let mut m = 0u128;
    while seen.contains(&m) {
        m += 1;
    }
    m
}

fn lower_mask(n: u128) -> u128 {
    assert!(n < 128, "coin bitmask positions must be < 128");
    if n == 0 {
        0
    } else {
        (1u128 << n) - 1
    }
}

/// Nim-multiplication via Conway's Turning-Corners excludant recurrence — the
/// *game* definition of the product.
pub fn nim_mul_mex(x: u128, y: u128) -> u128 {
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

// ---------------------------------------------------------------------------
// General coin-turning games and the Tartan (2-D) product
// ---------------------------------------------------------------------------
//
// A 1-D coin-turning game is specified by, for a single heads coin at position
// `n`, which subsets `S ⊆ {0,…,n−1}` of strictly-lower coins may be co-turned
// along with `n` (which is the rightmost coin of the move, the one going H→T).
// Its single-coin Grundy value is
//
//     g(n) = mex { ⊕_{i∈S} g(i) : S a legal companion set of n }.
//
// The **Tartan product** `G ⊗ H` is the 2-D coin-turning game whose move on a
// coin at (x,y) is a product `T_A × T_B` of a move `T_A` of G (rightmost x) and
// a move `T_B` of H (rightmost y): it turns the rectangle of cells (a,b) with
// a ∈ T_A, b ∈ T_B. The **Tartan/Product theorem** (Berlekamp–Conway–Guy,
// *Winning Ways*) says its single-coin Grundy value factors as the *nim-product*
// of the component values:
//
//     (G ⊗ H) grundy at (x,y) = g_G(x) ⊗ g_H(y).
//
// Turning Corners is the special case `G = H =` "turn n and exactly one lower
// coin" (which has g(n) = n), recovering `nim_mul_mex(x,y) = x ⊗ y`. The
// functions below compute the tartan Grundy from the 2-D excludant *directly*,
// so the theorem can be checked rather than assumed — see the tests.

/// "Turn coin `n` and exactly one strictly-lower coin." Single-coin Grundy
/// `g(n) = n`; its tartan square is Turning Corners.
pub fn singleton_companions(n: u128) -> Vec<u128> {
    assert!(n < 128, "coin bitmask positions must be < 128");
    (0..n).map(|i| 1u128 << i).collect()
}

/// Turning Turtles: "turn coin `n` and optionally one strictly-lower coin."
/// Single-coin Grundy `g(n) = n + 1`.
pub fn turtles_companions(n: u128) -> Vec<u128> {
    assert!(n < 128, "coin bitmask positions must be < 128");
    let mut v = vec![0u128]; // the empty companion set (turn n alone)
    v.extend((0..n).map(|i| 1u128 << i));
    v
}

/// Single-coin Grundy value of a 1-D coin-turning game (memoised). `companions`
/// returns the legal companion sets (bitmasks ⊆ {0..n-1}) for a coin at `n`.
pub fn grundy_1d<F: Fn(u128) -> Vec<u128>>(
    companions: &F,
    n: u128,
    memo: &mut HashMap<u128, u128>,
) -> u128 {
    let allowed = lower_mask(n);
    if let Some(&v) = memo.get(&n) {
        return v;
    }
    let mut seen = HashSet::new();
    for s in companions(n) {
        assert!(
            s & !allowed == 0,
            "companion mask contains a coin that is not strictly lower than n"
        );
        let mut acc = 0u128;
        let mut ss = s;
        while ss != 0 {
            let i = ss.trailing_zeros() as u128;
            ss &= ss - 1;
            acc ^= grundy_1d(companions, i, memo);
        }
        seen.insert(acc);
    }
    let g = mex(&seen);
    memo.insert(n, g);
    g
}

/// Single-coin Grundy value of the Tartan product of two 1-D coin-turning games
/// at cell (x,y), computed from the 2-D excludant directly (memoised). Compare
/// against `nim_mul` of the component 1-D Grundy values to check the theorem.
pub fn tartan_grundy<A, B>(
    comp_a: &A,
    comp_b: &B,
    x: u128,
    y: u128,
    memo: &mut HashMap<(u128, u128), u128>,
) -> u128
where
    A: Fn(u128) -> Vec<u128>,
    B: Fn(u128) -> Vec<u128>,
{
    assert!(x < 128 && y < 128, "tartan bitmask positions must be < 128");
    let allowed_x = lower_mask(x);
    let allowed_y = lower_mask(y);
    if let Some(&v) = memo.get(&(x, y)) {
        return v;
    }
    let mut seen = HashSet::new();
    for ta in comp_a(x) {
        assert!(
            ta & !allowed_x == 0,
            "row companion mask contains a coin that is not strictly lower than x"
        );
        let acells = ta | (1u128 << x); // rows turned (rightmost = x)
        for tb in comp_b(y) {
            assert!(
                tb & !allowed_y == 0,
                "column companion mask contains a coin that is not strictly lower than y"
            );
            let bcells = tb | (1u128 << y); // cols turned (rightmost = y)
                                            // XOR the values of every turned cell except (x,y) itself.
            let mut acc = 0u128;
            let mut aa = acells;
            while aa != 0 {
                let a = aa.trailing_zeros() as u128;
                aa &= aa - 1;
                let mut bb = bcells;
                while bb != 0 {
                    let b = bb.trailing_zeros() as u128;
                    bb &= bb - 1;
                    if a == x && b == y {
                        continue;
                    }
                    acc ^= tartan_grundy(comp_a, comp_b, a, b, memo);
                }
            }
            seen.insert(acc);
        }
    }
    let g = mex(&seen);
    memo.insert((x, y), g);
    g
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scalar::nim_mul;

    #[test]
    fn game_definition_equals_algebraic_nim_mul() {
        // Turning-Corners Grundy values == Fermat-power nim-multiplication.
        for x in 0u128..48 {
            for y in 0u128..48 {
                assert_eq!(nim_mul_mex(x, y), nim_mul(x, y), "mismatch at ({x}, {y})");
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

    #[test]
    fn one_d_grundy_sequences() {
        let mut m = HashMap::new();
        for n in 0..10 {
            assert_eq!(grundy_1d(&singleton_companions, n, &mut m), n); // g(n) = n
        }
        let mut m2 = HashMap::new();
        for n in 0..10 {
            assert_eq!(grundy_1d(&turtles_companions, n, &mut m2), n + 1); // g(n) = n+1
        }
    }

    #[test]
    #[should_panic(expected = "positions must be < 128")]
    fn companion_bitmasks_are_bounded_by_u128_width() {
        let _ = singleton_companions(128);
    }

    #[test]
    #[should_panic(expected = "not strictly lower than n")]
    fn custom_companion_masks_must_use_lower_coins() {
        let bad = |n| vec![1u128 << n];
        let mut memo = HashMap::new();
        let _ = grundy_1d(&bad, 3, &mut memo);
    }

    #[test]
    fn tartan_square_of_singleton_game_is_turning_corners() {
        // tartan(g(n)=n, itself) at (x,y) = x ⊗ y = Turning Corners = nim_mul_mex.
        let mut tm = HashMap::new();
        for x in 0u128..6 {
            for y in 0u128..6 {
                assert_eq!(
                    tartan_grundy(&singleton_companions, &singleton_companions, x, y, &mut tm),
                    nim_mul_mex(x, y),
                    "tartan ≠ Turning Corners at ({x},{y})"
                );
            }
        }
    }

    #[test]
    fn tartan_product_theorem() {
        // The Tartan/Product theorem on mixed component games: the 2-D Grundy
        // (from the excludant) equals the nim-product of the 1-D Grundy values.
        fn check<A: Fn(u128) -> Vec<u128>, B: Fn(u128) -> Vec<u128>>(ga: &A, gb: &B) {
            let (mut ma, mut mb, mut tm) = (HashMap::new(), HashMap::new(), HashMap::new());
            for x in 0u128..5 {
                for y in 0u128..5 {
                    let direct = tartan_grundy(ga, gb, x, y, &mut tm);
                    let factored = nim_mul(grundy_1d(ga, x, &mut ma), grundy_1d(gb, y, &mut mb));
                    assert_eq!(direct, factored, "Tartan theorem failed at ({x},{y})");
                }
            }
        }
        check(&singleton_companions, &singleton_companions);
        check(&singleton_companions, &turtles_companions);
        check(&turtles_companions, &singleton_companions);
        check(&turtles_companions, &turtles_companions);
    }
}
