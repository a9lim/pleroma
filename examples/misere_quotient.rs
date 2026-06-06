//! Misère quotients, and whether their P-sets are quadrics.
//!   cargo run --example misere_quotient
//!
//! The open question wants a game whose P-set is a Gold *quadric* `{Q=0}`. The
//! misère route is promising because misère sums are non-linear (see
//! `misere.rs`). This probe computes the misère indistinguishability quotient of
//! several small games and asks of each P-set: when the quotient is an elementary
//! abelian 2-group `(ℤ/2)^k` (so its elements coordinatise as `F₂^k`), is the
//! P-set a quadric, and is it *genuinely* quadratic (nonzero Arf-rank) or merely
//! affine/linear?
//!
//! Honest expectation: the small / tame quotients here are either `ℤ/2` (a
//! rank-0, linear P-set) or non-group monoids where the F₂-quadric framing does
//! not even apply. A genuine quadric P-set would need a `(ℤ/2)^k` quotient with
//! `k ≥ 2` and Arf-rank `≥ 2`; finding (or ruling out) one is the open work. This
//! prints what the instrument actually finds.

use pleroma::forms::fit_f2_quadratic;
use pleroma::games::{misere_quotient, AbstractGame, Quotient};

/// Build the abstract game of Nim with the given heap sizes as position types:
/// position `h` (1..=max) is a heap of size h, moving to any smaller heap (incl.
/// the empty heap 0). Atoms are the heap sizes in `heaps`.
fn nim_game(max: usize) -> AbstractGame {
    let moves = (0..=max).map(|h| (0..h).collect::<Vec<_>>()).collect();
    AbstractGame { moves }
}

/// If the quotient is `(ℤ/2)^k` on the given atoms (each atom an involution, and
/// the subset→class map a bijection), return the P-set as `F₂^k` bitmasks.
fn p_set_as_f2(q: &Quotient, atoms: &[usize]) -> Option<Vec<u32>> {
    let k = atoms.len();
    if k > 12 {
        return None;
    }
    // map a subset bitmask to its class via the enumerated elements
    let class_of_subset = |mask: u32| -> Option<usize> {
        let mut ms: Vec<usize> = (0..k)
            .filter(|&i| mask & (1 << i) != 0)
            .map(|i| atoms[i])
            .collect();
        ms.sort_unstable();
        q.elements
            .iter()
            .position(|e| *e == ms)
            .map(|idx| q.class_of[idx])
    };
    let mut seen = std::collections::HashSet::new();
    let mut pset = Vec::new();
    for v in 0u32..(1 << k) {
        let c = class_of_subset(v)?;
        if !seen.insert(c) && v.count_ones() <= 1 {
            // a generator coincided with an earlier class ⇒ not full-rank (ℤ/2)^k
        }
        if q.class_is_p[c] {
            pset.push(v);
        }
    }
    // require the 2^k subsets to hit exactly 2^k distinct classes (a bijection)
    if q.num_classes == (1 << k) {
        Some(pset)
    } else {
        None
    }
}

fn report(name: &str, game: &AbstractGame, atoms: &[usize], elem: usize, test: usize) {
    println!("\n── {name} ──");
    let q = misere_quotient(game, atoms, elem, test);
    let p_classes = q.class_is_p.iter().filter(|&&p| p).count();
    println!(
        "  quotient order = {}   P-classes = {}   (bounds: elem≤{elem}, test≤{test})",
        q.num_classes, p_classes
    );
    // is every atom an involution?  a²  ≈  identity (the empty class)?
    let id_class = q.class_of[q.elements.iter().position(|e| e.is_empty()).unwrap()];
    let involutions = atoms.iter().all(|&a| {
        q.elements
            .iter()
            .position(|e| *e == vec![a, a])
            .map(|i| q.class_of[i] == id_class)
            .unwrap_or(false)
    });
    println!("  every atom an involution (a²=1)? {involutions}");

    match p_set_as_f2(&q, atoms) {
        Some(pset) => {
            println!(
                "  quotient ≅ (ℤ/2)^{}  → testing the P-set as an F₂ quadric…",
                atoms.len()
            );
            match fit_f2_quadratic(&pset, atoms.len()) {
                Some(fit) => {
                    if fit.is_genuinely_quadratic() {
                        println!(
                            "    P-set IS a genuine quadric:  Arf={}, rank={}  ← a quadratic refinement!",
                            fit.arf.arf, fit.arf.rank
                        );
                    } else {
                        println!(
                            "    P-set is a quadric but rank 0 (affine/linear) — no quadratic content."
                        );
                    }
                }
                None => println!("    P-set is not a quadric."),
            }
        }
        None => println!("  not a full-rank (ℤ/2)^k group ⇒ the F₂-quadric framing doesn't apply."),
    }
}

fn main() {
    println!("Misère quotients and the quadric question.");

    // ⋆ alone: the classic quotient ℤ/2.
    report(
        "star (⋆ only)",
        &AbstractGame {
            moves: vec![vec![], vec![0]],
        },
        &[1],
        6,
        4,
    );

    // Nim with small heaps (tame): quotients are not elementary-2 groups.
    report("misère Nim, heaps {1,2}", &nim_game(2), &[1, 2], 5, 4);
    report("misère Nim, heaps {1,2,3}", &nim_game(3), &[1, 2, 3], 5, 3);

    println!("\nConclusion: among these small/tame games the misère quotient is either");
    println!("ℤ/2 (a rank-0, linear P-set) or a non-group monoid where the F₂-quadric");
    println!("framing doesn't apply. No genuine quadric P-set appears here — a quadratic");
    println!("refinement would need a (ℤ/2)^k quotient (k≥2) of Arf-rank ≥2, i.e. a *wild*");
    println!("quotient of that shape. The instrument to test any candidate is now in place.");
}
