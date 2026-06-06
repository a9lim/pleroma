//! Hunt octal games for a misère quotient that hosts a Gold-style quadric.
//!   cargo run --release --example octal_hunt
//!
//! The target shape (the one that would *close* the open question): a misère
//! quotient that is an elementary abelian 2-group `(ℤ/2)^k` — so its elements
//! coordinatise as `F₂^k` — whose P-set is a *genuine* quadric (Arf-rank ≥ 2).
//! Among the tame games tried earlier none had it. Octal games include the wild
//! ones, so we sweep a range of octal codes, compute each bounded misère
//! quotient over single-heap atoms, and (when the quotient is `(ℤ/2)^k`) run the
//! P-set through `fit_f2_quadratic`. Any hit is printed loudly; otherwise the
//! summary records how close anything got.

use pleroma::forms::fit_f2_quadratic;
use pleroma::games::{octal_misere_quotient, Quotient};

/// If the quotient on atoms `1..=k` is a full-rank `(ℤ/2)^k` (the `2^k` squarefree
/// subsets hit all `2^k` classes bijectively), return its P-set as `F₂^k` masks.
fn p_set_as_f2(q: &Quotient, k: usize) -> Option<Vec<u32>> {
    if q.num_classes != (1 << k) {
        return None;
    }
    let class_of_subset = |mask: u32| -> Option<usize> {
        let mut ms: Vec<usize> = (0..k)
            .filter(|&i| mask & (1 << i) != 0)
            .map(|i| i + 1)
            .collect();
        ms.sort_unstable();
        q.elements
            .iter()
            .position(|e| *e == ms)
            .map(|idx| q.class_of[idx])
    };
    let mut hit = std::collections::HashSet::new();
    let mut pset = Vec::new();
    for v in 0u32..(1 << k) {
        let c = class_of_subset(v)?;
        hit.insert(c);
        if q.class_is_p[c] {
            pset.push(v);
        }
    }
    if hit.len() == (1 << k) {
        Some(pset)
    } else {
        None
    }
}

fn code_str(code: &[u8]) -> String {
    format!(
        "0.{}",
        code.iter().map(|d| d.to_string()).collect::<String>()
    )
}

fn main() {
    let max_heap = 4usize;
    let (elem, test) = (5usize, 4usize);

    // Sweep: all codes of length 1–3 whose first digit allows taking a whole heap
    // of 1 (d₁ odd) — otherwise heaps of size 1 are inert and the game degenerates.
    let mut codes: Vec<Vec<u8>> = Vec::new();
    for d1 in [1u8, 3, 5, 7] {
        codes.push(vec![d1]);
        for d2 in 0u8..8 {
            codes.push(vec![d1, d2]);
            for d3 in 0u8..8 {
                codes.push(vec![d1, d2, d3]);
            }
        }
    }
    println!(
        "Hunting {} octal codes (max_heap={max_heap}, bounds elem≤{elem}/test≤{test})…\n",
        codes.len()
    );

    let mut hits = 0;
    let mut two_groups = 0; // quotients that came out (ℤ/2)^k with k≥2
    let mut order_hist: std::collections::BTreeMap<usize, usize> =
        std::collections::BTreeMap::new();

    for code in &codes {
        // scan the heap cutoff: a clean group may appear only among small heaps.
        for k in 2..=max_heap {
            let q = octal_misere_quotient(code, k, elem, test);
            *order_hist.entry(q.num_classes).or_insert(0) += 1;
            if let Some(pset) = p_set_as_f2(&q, k) {
                two_groups += 1;
                if let Some(fit) = fit_f2_quadratic(&pset, k) {
                    if fit.is_genuinely_quadratic() {
                        hits += 1;
                        println!(
                            "  ★ HIT: {} on heaps 1..={k}: (ℤ/2)^{k}, P-set quadric Arf={} rank={}",
                            code_str(code),
                            fit.arf.arf,
                            fit.arf.rank
                        );
                    } else {
                        println!(
                            "  · {} heaps 1..={k}: (ℤ/2)^{k} but P-set affine/linear (rank 0)",
                            code_str(code)
                        );
                    }
                }
            }
        }
    }

    println!("\n── summary ──");
    println!(
        "codes swept (×heap-cutoffs): {}",
        codes.len() * (max_heap - 1)
    );
    println!("(ℤ/2)^k quotients found (k≥2): {two_groups}");
    println!("genuine-quadric P-set HITS:   {hits}");
    println!("quotient-order histogram (order: count):");
    for (ord, cnt) in &order_hist {
        println!("  {ord:>3}: {cnt}");
    }
    if hits == 0 {
        println!("\nNo octal game in this range has a (ℤ/2)^k misère quotient whose P-set is a");
        println!("genuine quadric. The open question survives the hunt: the quadric P-set, if it");
        println!("exists, lives outside the elementary-2-abelian octal quotients reachable here.");
    }
}
