//! Misère play: where disjunctive sums stop being linear.
//!
//! NOTES.md's open question needs a game whose P-positions are the *quadric*
//! `{Q=0}` of a Gold form. Normal-play disjunctive sums can't supply one: their
//! outcomes are XOR-linear (P ⟺ ⊕ of Grundy values = 0), so the P-set is always
//! a *subspace*. The two escape routes are interactive games and **misère** play
//! (last player to move loses), where sums are genuinely non-linear — Grundy
//! values no longer determine the outcome of a sum.
//!
//! This module is the instrument for the misère route: a memoised misère-outcome
//! evaluator for any finite impartial game (given a `moves` function), plus the
//! canonical witness that misère is non-linear — misère Nim, whose P-set is
//! provably *not* `{⊕ = 0}` and not even a coset. That clears the bar normal play
//! fails. Whether a misère game's P-set is an actual Gold quadric is the part
//! that stays open; this gives the tooling to test candidates.

use std::collections::HashMap;
use std::hash::Hash;

/// Misère outcome of a finite impartial game: `true` = **N-position** (the player
/// to move wins under misère, last-to-move-loses), `false` = **P-position** (the
/// previous player wins). `moves(p)` lists the positions reachable in one move; a
/// position with no moves is terminal, and under misère the player who *cannot*
/// move **wins**, so a terminal position is an N-position. Memoised on positions.
pub fn misere_is_n<P, F>(pos: &P, moves: &F, memo: &mut HashMap<P, bool>) -> bool
where
    P: Clone + Eq + Hash,
    F: Fn(&P) -> Vec<P>,
{
    if let Some(&v) = memo.get(pos) {
        return v;
    }
    let nexts = moves(pos);
    // terminal ⇒ N (can't-move wins); otherwise N ⟺ some move reaches a P.
    let result = nexts.is_empty() || nexts.iter().any(|q| !misere_is_n(q, moves, memo));
    memo.insert(pos.clone(), result);
    result
}

/// Convenience: `true` iff `pos` is a misère P-position (second player wins).
pub fn misere_is_p<P, F>(pos: &P, moves: &F, memo: &mut HashMap<P, bool>) -> bool
where
    P: Clone + Eq + Hash,
    F: Fn(&P) -> Vec<P>,
{
    !misere_is_n(pos, moves, memo)
}

/// A Nim position: heap sizes, kept sorted ascending with empty heaps dropped so
/// equal positions share a memo key.
pub fn nim_canonical(mut heaps: Vec<u128>) -> Vec<u128> {
    heaps.retain(|&h| h != 0);
    heaps.sort_unstable();
    heaps
}

/// The moves of Nim: reduce any one heap to any strictly smaller size.
pub fn nim_moves(pos: &Vec<u128>) -> Vec<Vec<u128>> {
    let mut out = Vec::new();
    for i in 0..pos.len() {
        for v in 0..pos[i] {
            let mut q = pos.clone();
            q[i] = v;
            out.push(nim_canonical(q));
        }
    }
    out
}

/// The misère-Nim theorem (Bouton): a position is a misère P-position iff either
/// every heap is ≤ 1 and there is an *odd* number of heaps, or some heap is ≥ 2
/// and the nim-sum (XOR) of the heaps is 0. (The empty position is N.)
pub fn misere_nim_p_predicted(heaps: &[u128]) -> bool {
    let xor = heaps.iter().fold(0u128, |a, &h| a ^ h);
    let max = heaps.iter().copied().max().unwrap_or(0);
    if max <= 1 {
        heaps.len() % 2 == 1
    } else {
        xor == 0
    }
}

// ---------------------------------------------------------------------------
// The misère indistinguishability quotient (Plambeck–Siegel), bounded
// ---------------------------------------------------------------------------
//
// Two positions G, H are *indistinguishable* if outcome(G+X) = outcome(H+X) for
// every test X; the equivalence classes form a commutative monoid (the misère
// quotient) carrying a distinguished P-set. We compute it *bounded*: positions
// are sums of atoms up to `elem_bound`, tested against sums up to `test_bound`.
// For a game with a finite quotient this is exact once the bounds exceed its
// pretension; otherwise it is a finite observational approximation of the
// congruence (bounded tests may merge more positions than the true quotient,
// i.e. return a coarser quotient). The point of the
// instrument is to ask, of the resulting P-set, the question the project cares
// about: is it a quadric, and what is its Arf (win-bias)?

/// An abstract finite impartial game: position 0 is the empty game (the identity
/// under disjunctive sum, with no moves); positions `1..moves.len()` carry option
/// index-lists `moves[p]` (each option is a position index; 0 = move to empty).
pub struct AbstractGame {
    pub moves: Vec<Vec<usize>>,
}

impl AbstractGame {
    /// Moves of a disjunctive sum (a multiset of nonzero component positions): in
    /// any one component, replace it by one of its options (dropping the empty).
    fn sum_moves(&self, pos: &Vec<usize>) -> Vec<Vec<usize>> {
        let mut out = Vec::new();
        for idx in 0..pos.len() {
            for &q in &self.moves[pos[idx]] {
                let mut np = pos.clone();
                if q == 0 {
                    np.remove(idx);
                } else {
                    np[idx] = q;
                }
                np.sort_unstable();
                out.push(np);
            }
        }
        out
    }

    fn canon(pos: &[usize]) -> Vec<usize> {
        let mut v: Vec<usize> = pos.iter().copied().filter(|&p| p != 0).collect();
        v.sort_unstable();
        v
    }

    /// Misère outcome of a sum (multiset of component positions): `true` = N.
    pub fn misere_outcome(&self, pos: &[usize], memo: &mut HashMap<Vec<usize>, bool>) -> bool {
        let canon = Self::canon(pos);
        misere_is_n(&canon, &|p| self.sum_moves(p), memo)
    }
}

/// All sorted multisets of `atoms` (assumed sorted) with total length `0..=max`.
fn multisets(atoms: &[usize], max: usize) -> Vec<Vec<usize>> {
    let mut result = vec![vec![]];
    let mut frontier = vec![vec![]];
    for _ in 0..max {
        let mut next = Vec::new();
        for m in &frontier {
            let last = m.last().copied().unwrap_or(0);
            for &a in atoms.iter().filter(|&&a| a >= last) {
                let mut nm = m.clone();
                nm.push(a);
                next.push(nm);
            }
        }
        result.extend(next.iter().cloned());
        frontier = next;
    }
    result
}

/// A bounded misère indistinguishability quotient.
#[derive(Debug, Clone)]
pub struct Quotient {
    /// The enumerated elements (sorted multisets of atoms, up to `elem_bound`).
    pub elements: Vec<Vec<usize>>,
    /// Class id of each element (parallel to `elements`).
    pub class_of: Vec<usize>,
    /// Number of distinct classes found.
    pub num_classes: usize,
    /// A representative multiset for each class.
    pub class_rep: Vec<Vec<usize>>,
    /// P-status of each class (`true` = a misère P-position / second-player win).
    pub class_is_p: Vec<bool>,
}

/// Build a quotient from `elements` and a `tests` set, given an `outcome`
/// function (`true` = N) on atom-multisets. `outcome` carries its own memo. Two
/// elements share a class iff `outcome(G+T)` agrees for every test `T`.
fn build_quotient(
    elements: Vec<Vec<usize>>,
    tests: &[Vec<usize>],
    mut outcome: impl FnMut(&[usize]) -> bool,
) -> Quotient {
    let signatures: Vec<Vec<bool>> = elements
        .iter()
        .map(|g| {
            tests
                .iter()
                .map(|t| {
                    let mut gt = g.clone();
                    gt.extend_from_slice(t);
                    outcome(&gt)
                })
                .collect()
        })
        .collect();

    let mut class_of = vec![0usize; elements.len()];
    let mut uniq: Vec<Vec<bool>> = Vec::new();
    let mut class_rep: Vec<Vec<usize>> = Vec::new();
    for (i, sig) in signatures.iter().enumerate() {
        match uniq.iter().position(|s| s == sig) {
            Some(c) => class_of[i] = c,
            None => {
                class_of[i] = uniq.len();
                uniq.push(sig.clone());
                class_rep.push(elements[i].clone());
            }
        }
    }
    let class_is_p: Vec<bool> = class_rep.iter().map(|r| !outcome(r)).collect();

    Quotient {
        num_classes: uniq.len(),
        elements,
        class_of,
        class_rep,
        class_is_p,
    }
}

/// Compute the bounded misère quotient of `game` over the generating `atoms`,
/// distinguishing elements (sums up to `elem_bound`) by their outcomes against
/// tests (sums up to `test_bound`).
pub fn misere_quotient(
    game: &AbstractGame,
    atoms: &[usize],
    elem_bound: usize,
    test_bound: usize,
) -> Quotient {
    let mut atoms_sorted = atoms.to_vec();
    atoms_sorted.sort_unstable();
    let elements = multisets(&atoms_sorted, elem_bound);
    let tests = multisets(&atoms_sorted, test_bound);
    let mut memo: HashMap<Vec<usize>, bool> = HashMap::new();
    build_quotient(elements, &tests, |g| game.misere_outcome(g, &mut memo))
}

// ---------------------------------------------------------------------------
// Octal games — the wild hunting ground for a quadric P-set
// ---------------------------------------------------------------------------

/// Moves of an octal game `0.d₁d₂…` (`code[k-1] = dₖ`) on a heap-multiset. From a
/// heap of size n, remove k tokens (1 ≤ k ≤ n): leaving the heap empty needs
/// `dₖ & 1`; leaving one nonempty heap `n−k` needs `dₖ & 2`; splitting `n−k` into
/// two nonempty heaps needs `dₖ & 4`. (Nim is `0.333…`, Dawson's chess `0.137`.)
pub fn octal_moves(code: &[u8], pos: &[u128]) -> Vec<Vec<u128>> {
    let mut out = Vec::new();
    for idx in 0..pos.len() {
        let n = pos[idx];
        let base: Vec<u128> = pos
            .iter()
            .enumerate()
            .filter(|&(i, _)| i != idx)
            .map(|(_, &h)| h)
            .collect();
        for k in 1..=n {
            let d = *code.get((k - 1) as usize).unwrap_or(&0);
            let rem = n - k;
            if rem == 0 {
                if d & 1 != 0 {
                    let mut p = base.clone();
                    p.sort_unstable();
                    out.push(p);
                }
            } else {
                if d & 2 != 0 {
                    let mut p = base.clone();
                    p.push(rem);
                    p.sort_unstable();
                    out.push(p);
                }
                if d & 4 != 0 {
                    for a in 1..=rem / 2 {
                        let mut p = base.clone();
                        p.push(a);
                        p.push(rem - a);
                        p.sort_unstable();
                        out.push(p);
                    }
                }
            }
        }
    }
    out
}

/// The bounded misère quotient of an octal game, over single heaps of size
/// `1..=max_heap` as atoms (a heap-multiset is a sum). Splitting moves are handled
/// (a heap can become two), so the position type is the heap-multiset itself.
pub fn octal_misere_quotient(
    code: &[u8],
    max_heap: usize,
    elem_bound: usize,
    test_bound: usize,
) -> Quotient {
    let atoms: Vec<usize> = (1..=max_heap).collect();
    let elements = multisets(&atoms, elem_bound);
    let tests = multisets(&atoms, test_bound);
    let mut memo: HashMap<Vec<u128>, bool> = HashMap::new();
    let moves = |p: &Vec<u128>| octal_moves(code, p);
    build_quotient(elements, &tests, |g| {
        let mut pos: Vec<u128> = g.iter().map(|&x| x as u128).collect();
        pos.sort_unstable();
        misere_is_n(&pos, &moves, &mut memo)
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn misere_nim_matches_boutons_theorem() {
        // Verify the tree evaluator against the closed-form theorem over all Nim
        // positions with up to 4 heaps of size ≤ 4.
        let mut memo: HashMap<Vec<u128>, bool> = HashMap::new();
        fn rec(prefix: &mut Vec<u128>, depth: usize, memo: &mut HashMap<Vec<u128>, bool>) {
            if depth == 0 {
                let pos = nim_canonical(prefix.clone());
                let is_p = misere_is_p(&pos, &nim_moves, memo);
                assert_eq!(
                    is_p,
                    misere_nim_p_predicted(&pos),
                    "misère Nim mismatch at {pos:?}"
                );
                return;
            }
            for h in 0..=4u128 {
                prefix.push(h);
                rec(prefix, depth - 1, memo);
                prefix.pop();
            }
        }
        rec(&mut Vec::new(), 4, &mut memo);
    }

    #[test]
    fn star_misere_quotient_is_z2() {
        // ⋆ = position 1, moving only to 0 (empty). Its misère quotient is the
        // group ℤ/2 = {1, a | a²=1}, with P-set {a} (an odd number of ⋆'s).
        let star = AbstractGame {
            moves: vec![vec![], vec![0]],
        };
        let q = misere_quotient(&star, &[1], 5, 3);
        assert_eq!(q.num_classes, 2, "⋆ quotient should be order 2 (ℤ/2)");
        // the empty position is N (not P); a single ⋆ is P.
        let empty_class = q.class_of[q.elements.iter().position(|e| e.is_empty()).unwrap()];
        let star_class = q.class_of[q.elements.iter().position(|e| e == &vec![1]).unwrap()];
        assert!(!q.class_is_p[empty_class]);
        assert!(q.class_is_p[star_class]);
        // a²=1: two ⋆'s fall in the identity (empty) class.
        let two = q.class_of[q.elements.iter().position(|e| e == &vec![1, 1]).unwrap()];
        assert_eq!(two, empty_class);
        // exactly one P-class (the win-bias is a single coset)
        assert_eq!(q.class_is_p.iter().filter(|&&p| p).count(), 1);
    }

    #[test]
    fn octal_nim_matches_misere_nim() {
        // 0.333… is Nim: octal moves' misère outcomes match Bouton's theorem.
        let code = [3u8, 3, 3, 3];
        let mut memo: HashMap<Vec<u128>, bool> = HashMap::new();
        for heaps in [
            vec![1u128],
            vec![1, 1],
            vec![2],
            vec![2, 1],
            vec![3, 2, 1],
            vec![2, 2],
            vec![3, 3],
        ] {
            let mut h = heaps.clone();
            h.sort_unstable();
            let is_n = misere_is_n(&h, &|p| octal_moves(&code, p), &mut memo);
            assert_eq!(
                is_n,
                !misere_nim_p_predicted(&heaps),
                "octal Nim ≠ Bouton at {heaps:?}"
            );
        }
    }

    #[test]
    fn octal_star_quotient_is_z2() {
        // Nim restricted to heaps of size 1 (just ⋆) ⇒ the ℤ/2 quotient again.
        let q = octal_misere_quotient(&[3, 3, 3], 1, 5, 3);
        assert_eq!(q.num_classes, 2);
    }

    #[test]
    fn misere_is_genuinely_nonlinear() {
        // The normal-play P-set is exactly {XOR = 0} — a subspace. The misère
        // P-set is NOT: it contains a XOR≠0 point and excludes a XOR=0 point, so
        // it is neither {XOR=0} nor a coset of any subspace. This is precisely the
        // non-linearity normal-play sums lack (and that a quadratic P-set needs).
        let mut memo: HashMap<Vec<u128>, bool> = HashMap::new();
        let one = nim_canonical(vec![1]); // XOR = 1, but misère-P (you must take the last coin)
        let oneone = nim_canonical(vec![1, 1]); // XOR = 0, but misère-N
        assert!(misere_is_p(&one, &nim_moves, &mut memo));
        assert!(!misere_is_p(&oneone, &nim_moves, &mut memo));
        // 0 ∈ P-set?  empty position is terminal ⇒ N, so 0 ∉ misère-P. A subspace
        // (or its outcome set) would contain 0; a coset structure is impossible
        // because [1] (xor 1) is P while [1,1]+[1,1]-style xor-0 combos are N.
        let empty = nim_canonical(vec![]);
        assert!(!misere_is_p(&empty, &nim_moves, &mut memo));

        // Concrete subspace-failure witness: u=[1], v=[1] are both in the P-set
        // under the all-ones regime, but their nim-sum (xor) leaves the regime.
        // (Here the point is structural: outcome is not an XOR-linear function.)
        let three_ones = nim_canonical(vec![1, 1, 1]); // XOR = 1, misère-P (odd count)
        assert!(misere_is_p(&three_ones, &nim_moves, &mut memo));
    }
}
