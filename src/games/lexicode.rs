//! Bridge O — binary **lexicodes**: greedy = mex, the games ↔ integral edge.
//!
//! The lexicode `L(n, d)`: scan `F₂ⁿ` in lexicographic order and greedily keep every
//! vector at Hamming distance `≥ d` from all vectors kept so far. Conway–Sloane
//! (*Lexicographic codes: error-correcting codes from game theory*, IEEE Trans.
//! Inform. Theory **32** (1986) 337–348) prove the resulting set is **linear**, and
//! the proof is game theory: the greedy scan is the **mex** rule, the codewords are
//! the Grundy-value-0 positions of a turning game, and XOR-closure is Sprague–Grundy
//! theory rather than coding theory. The celebrated instances are shipped codes —
//! the `[7,4,3]` Hamming code, the `[8,4,4]` extended Hamming code, and the
//! `[24,12,8]` extended binary **Golay** code are all lexicodes — so each acquires a
//! second, game-theoretic construction. That makes a full chain executable:
//!
//! ```text
//! mex → lexicode → Golay → Construction A → even unimodular rank 24 (with roots) → theta
//! ```
//!
//! Every link past the first is already shipped (Bridges H/E); this file supplies the
//! first, closing the one pillar edge the bridge graph still lacked.
//!
//! **The mex step (executable, self-contained).** After keeping a code `C`, let
//! `Forbidden = ⋃_{c∈C} { m : d(m,c) < d }` be the union of radius-`(d−1)` Hamming
//! balls. The next greedy codeword is exactly `mex(Forbidden)` — the least vector not
//! excluded ([`crate::games::grundy::mex`]). The deeper Conway–Sloane turning-game
//! realization (the Grundy-value theorem) is cited, not reconstructed here: it is to
//! be transcribed from the 1986 paper in the formalization pass.
//!
//! **Relation to `OPEN.md` §1 (interpretation level).** `OPEN.md` §1 records that
//! normal-play P-sets are *linear* in Grundy coordinates. Lexicodes are the classical
//! demonstration of the **solved** (degree-1) side of that line: a fixed, natural,
//! non-tautological greedy rule whose P-set is a rich linear code. This bridge
//! supplies that degree-1 case as executable context; it does **not** touch the open
//! Gold-quadric question and must not be cited as progress on it.
//!
//! **Convention.** The lexicographic order is the standard bit order with
//! **coordinate 0 the most significant bit**, so scanning integers `0,1,2,…` upward
//! *is* the lexicographic scan. A permuted coordinate order gives a different
//! (equivalent) code. Binary only; the nim-field `2^{2^k}` linearity statement is
//! documented context, not a shipped surface.

use crate::forms::BinaryCode;

/// Backstop on the incremental search (codeword comparisons), mirroring
/// [`crate::forms::AUTO_NODE_BUDGET`]: past it, [`lexicode`] reports `None` rather
/// than running unbounded — an honest budget, not a silent cap.
pub const LEXICODE_NODE_BUDGET: u128 = 50_000_000_000;

/// Decode the integer `mask` into a length-`n` codeword row (coordinate 0 = MSB), the
/// `Vec<u8>` form [`BinaryCode::new`] expects.
fn mask_to_row(mask: u32, n: usize) -> Vec<u8> {
    (0..n).map(|i| ((mask >> (n - 1 - i)) & 1) as u8).collect()
}

/// A GF(2) basis of the span of `vectors` (integers as bit-vectors), by XOR
/// elimination keyed on the highest set bit.
fn bitmask_basis(vectors: &[u32]) -> Vec<u32> {
    let mut basis: Vec<u32> = Vec::new();
    for &v in vectors {
        let mut x = v;
        for &b in &basis {
            let hb = 31 - b.leading_zeros();
            if (x >> hb) & 1 == 1 {
                x ^= b;
            }
        }
        if x != 0 {
            basis.push(x);
        }
    }
    basis
}

/// The **literal** greedy lexicode `L(n, d)` for small `n` (`≤ 14`): scan every
/// vector of `F₂ⁿ` in lexicographic order, keep those at Hamming distance `≥ d` from
/// all kept so far, then **discover-don't-assert** — verify the kept set is closed
/// under XOR (the linearity theorem) and return `None` on a closure failure (which
/// would *falsify* linearity rather than hide it). The reference implementation that
/// pins the production [`lexicode`]. `None` for `n = 0` or `n > 14` (the `2ⁿ`
/// budget) and on a closure failure.
pub fn lexicode_naive(n: usize, d: usize) -> Option<BinaryCode> {
    if n == 0 || n > 14 {
        return None;
    }
    let size: u32 = 1u32 << n;
    let mut kept: Vec<u32> = Vec::new();
    for m in 0..size {
        if kept.iter().all(|&c| (m ^ c).count_ones() as usize >= d) {
            kept.push(m);
        }
    }
    // discover-don't-assert: the kept set must be linear (Conway–Sloane).
    let set: std::collections::HashSet<u32> = kept.iter().copied().collect();
    for &a in &kept {
        for &b in &kept {
            if !set.contains(&(a ^ b)) {
                return None;
            }
        }
    }
    let basis = bitmask_basis(&kept);
    BinaryCode::new(n, basis.iter().map(|&v| mask_to_row(v, n)).collect())
}

/// The production lexicode `L(n, d)`, built incrementally on linearity: the next
/// generator is the lex-least vector `v` whose coset `v + C` has minimum weight `≥ d`
/// (i.e. `d(v, C) ≥ d`).
///
/// Rather than recompute `d(v, C)` by scanning codewords, it carries the whole
/// distance array `dist[v] = d(v, C)` and updates it in one `O(2ⁿ)` pass per
/// generator via the coset recurrence
///
/// ```text
/// d(v, C ∪ (g ⊕ C)) = min( d(v, C), d(v ⊕ g, C) )
/// ```
///
/// (safe to apply in place: the `g`-twin's stale read only ever re-supplies a term
/// already in the `min`). The generator cursor is monotone — a vector killed to
/// `dist < d` never revives, since `dist` only decreases — so the total scan is a
/// single sweep. Budgeted by [`LEXICODE_NODE_BUDGET`] (`None` past it). `None` for
/// `n = 0` or `n > 26` (the `2ⁿ`-byte distance array).
pub fn lexicode(n: usize, d: usize) -> Option<BinaryCode> {
    lexicode_bounded(n, d, LEXICODE_NODE_BUDGET)
}

/// [`lexicode`] with an explicit operation budget (distance-array reads/writes).
pub fn lexicode_bounded(n: usize, d: usize, node_budget: u128) -> Option<BinaryCode> {
    if n == 0 || n > 26 {
        return None;
    }
    let size: usize = 1usize << n;
    // dist[v] = d(v, C); initially C = {0}, so dist[v] = weight(v).
    let mut dist: Vec<u8> = (0..size).map(|v| (v as u32).count_ones() as u8).collect();
    let mut basis: Vec<u32> = Vec::new();
    let mut budget = node_budget;
    let mut cursor: usize = 1; // v = 0 is already in the code (dist 0)
    loop {
        // lex-least v ≥ cursor still at distance ≥ d from the code is the next generator.
        while cursor < size && (dist[cursor] as usize) < d {
            cursor += 1;
        }
        if cursor >= size {
            break;
        }
        let g = cursor;
        basis.push(g as u32);
        // dist[v] ← min(dist[v], dist[v ⊕ g]) over the whole array (one pass).
        for v in 0..size {
            if budget == 0 {
                return None;
            }
            budget -= 1;
            let alt = dist[v ^ g];
            if alt < dist[v] {
                dist[v] = alt;
            }
        }
        cursor += 1;
    }
    BinaryCode::new(n, basis.iter().map(|&g| mask_to_row(g, n)).collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forms::{extended_hamming_code, golay_code, hamming_code};
    use crate::games::grundy::mex;

    /// Brute-force greedy kept-set as raw masks (for the mex witness).
    fn greedy_masks(n: usize, d: usize) -> Vec<u32> {
        let mut kept: Vec<u32> = Vec::new();
        for m in 0..(1u32 << n) {
            if kept.iter().all(|&c| (m ^ c).count_ones() as usize >= d) {
                kept.push(m);
            }
        }
        kept
    }

    #[test]
    fn greedy_step_is_mex_of_the_forbidden_balls() {
        // At every step the next kept vector is mex(Forbidden), Forbidden = the union
        // of radius-(d−1) balls around the kept codewords. Reconstruct the greedy scan
        // purely through `mex` and check it reproduces the direct scan.
        let (n, d) = (6usize, 3usize);
        let direct = greedy_masks(n, d);

        let mut kept: Vec<u32> = Vec::new();
        loop {
            // Forbidden = vectors within distance < d of some kept codeword.
            let forbidden: Vec<u128> = (0..(1u32 << n))
                .filter(|&m| kept.iter().any(|&c| ((m ^ c).count_ones() as usize) < d))
                .map(u128::from)
                .collect();
            let next = mex(forbidden);
            if next >= (1u128 << n) {
                break; // every remaining vector is forbidden — scan complete
            }
            kept.push(next as u32);
        }
        assert_eq!(
            kept, direct,
            "mex reconstruction must equal the greedy scan"
        );
    }

    #[test]
    fn naive_and_production_agree_for_small_n() {
        // The literal greedy scan pins the production (distance-array) route.
        for n in 1..=12 {
            for d in 1..=4 {
                let a = lexicode_naive(n, d);
                let b = lexicode(n, d);
                assert_eq!(a, b, "lexicode_naive vs lexicode at (n={n}, d={d})");
            }
        }
    }

    #[test]
    fn distance_one_is_the_full_space_and_two_is_even_weight() {
        // d = 1: no constraint ⇒ all of F₂ⁿ.
        let full = lexicode(5, 1).unwrap();
        assert_eq!(full.len(), 5);
        assert_eq!(full.dim(), 5);
        // d = 2: the even-weight code [n, n−1, 2].
        let even = lexicode(5, 2).unwrap();
        assert_eq!((even.len(), even.dim()), (5, 4));
        assert_eq!(even.minimum_distance(), Some(2));
    }

    #[test]
    fn lexicode_reproduces_hamming_codes() {
        // [7,4,3] Hamming and [8,4,4] extended Hamming as lexicodes.
        let h = lexicode(7, 3).unwrap();
        assert_eq!((h.len(), h.dim(), h.minimum_distance()), (7, 4, Some(3)));
        let eh = lexicode(8, 4).unwrap();
        assert_eq!((eh.len(), eh.dim(), eh.minimum_distance()), (8, 4, Some(4)));
        // Permutation-invariant identity bundle (the bit order may differ from the
        // shipped constructors' — assert weight enumerator, the equivalence invariant).
        assert_eq!(h.weight_enumerator(), hamming_code().weight_enumerator());
        assert_eq!(
            eh.weight_enumerator(),
            extended_hamming_code().weight_enumerator()
        );
    }

    #[test]
    fn lexicode_24_8_is_golay_and_chains_to_a_lattice_with_roots() {
        let g = lexicode(24, 8).expect("lexicode(24,8) within budget");
        // The [24,12,8] doubly-even self-dual predicate bundle.
        assert_eq!(g.len(), 24);
        assert_eq!(g.dim(), 12);
        assert_eq!(g.minimum_distance(), Some(8));
        assert!(g.is_doubly_even());
        assert!(g.is_self_dual());
        // Uniqueness of the [24,12,8] Type II code (MacWilliams–Sloane; Pless) upgrades
        // the bundle to "is the Golay code": equal weight enumerators ⇒ equivalent.
        assert_eq!(g.weight_enumerator(), golay_code().weight_enumerator());

        // The chain rung: Construction A of the length-24 lexicode is even unimodular
        // rank 24 *with* roots (≠ Leech) — re-pinning Bridge H's boundary from games.
        let lattice = g
            .construction_a()
            .expect("doubly-even self-dual ⇒ integral Gram");
        assert!(lattice.is_even());
        assert!(lattice.is_unimodular());
        let roots = lattice.short_vectors(2).expect("definite ⇒ enumerable");
        assert!(
            !roots.is_empty(),
            "Golay Construction A has roots, unlike Leech"
        );
    }
}
