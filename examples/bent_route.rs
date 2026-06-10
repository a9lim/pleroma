//! The bent route to the open question.
//!   cargo run --example bent_route
//!
//! Bent (nondegenerate) game-realizable forms are a clean Tier-2 target. The
//! polar form B then has trivial radical R(B) = {0}, so:
//!   (i)  the symmetric-B loopy rule whose Loss-set is R(B) (loopy_quadric.rs)
//!        collapses to {0} — the radical route is empty, no coincidence can fake a
//!        hit (cf. the (m,a)=(4,1) artifact where R(B) happened to equal {Q=0});
//!   (ii) the frame-blind Sp(B) no-go applies without a degenerate radical layer.
//!        So any candidate rule in this setting must use more than B alone.
//!
//! Read the form as an ISING ENERGY over F_2:
//!     Q(v) = Σ_{i<j} B_ij v_i v_j  +  Σ_i q_i v_i,
//! pairwise couplings B (game-realizable: the coin-turning polar form) plus a
//! per-coin local field q_i = Q(e_i) (the single-coin self-Gold value, also
//! game-realizable). {Q=0} is the even-energy set; bent ⇔ B nondegenerate.
//!
//! The decisive new probe is a LOCAL SPIN-FLIP rule: flip a single coin i when the
//! local energy change ΔQ_i(v) = q_i ⊕ B(v,e_i) is 1. It reads ONLY the couplings
//! B and the per-coin field q_i — Tier-2 data (B + a diagonal frame), never the
//! global Q. If its P-set is {Q=0}, the open question reduces to whether "pairwise
//! coupling + per-coin field, played as spin flips" counts as natural. We use a
//! bent Gold COMPONENT Tr(λ x^{1+2^a}) (bent for 2/3 of λ; see gold_family_survey).

use ogdoad::forms::fit_f2_quadratic;
use ogdoad::games::loopy_decision_sets;

mod common;
use common::{bent_gold as gold, bent_polar as polar, p_set};

fn describe(label: &str, p: &[u128], zero: &[u128], draws: usize, m: u128) {
    let n = 1usize << m;
    let ps: std::collections::HashSet<u128> = p.iter().copied().collect();
    let zs: std::collections::HashSet<u128> = zero.iter().copied().collect();
    let agree = (0..n as u128)
        .filter(|v| ps.contains(v) == zs.contains(v))
        .count();
    print!(
        "  {label:<28} |P|={:<3} draws={draws:<3} agree {agree}/{n} with {{Q=0}}",
        p.len()
    );
    if p == zero {
        println!("   ← EXACTLY {{Q=0}} !!");
    } else {
        match fit_f2_quadratic(p, m as usize) {
            Some(f) if f.is_genuinely_quadratic() => {
                let bent = if f.arf.rank == m as usize {
                    ", BENT"
                } else {
                    ""
                };
                println!("   quadric (Arf={}, rank={}{bent})", f.arf.arf, f.arf.rank)
            }
            Some(_) => println!("   affine/linear (a subspace coset)"),
            None => println!("   not a quadric"),
        }
    }
}

fn main() {
    let (m, a) = (8u128, 1u128);
    let n = 1usize << m;

    // Find a bent witness: |{Q=0}| = 2^{m-1} ± 2^{m/2-1} ⇔ nondegenerate.
    let half = 1u128 << (m - 1);
    let off = 1u128 << (m / 2 - 1);
    let lam = (1..1u128 << m)
        .find(|&l| {
            let z = (0..1u128 << m).filter(|&v| gold(v, l, a, m) == 0).count() as u128;
            z == half + off || z == half - off
        })
        .expect("a bent component exists");
    let zero: Vec<u128> = (0..n as u128)
        .filter(|&v| gold(v, lam, a, m) == 0)
        .collect();
    let z = zero.len() as u128;
    let arf = if z == half + off { 0 } else { 1 };
    println!(
        "Bent Gold component  Q(v) = Tr({lam}·v^{{1+2^{a}}})  on F_2^{m}:  bent, Arf={arf}, |{{Q=0}}|={z}\n"
    );

    // (i) the radical route is dead: symmetric-B loopy Loss-set = R(B) = {0}.
    let (loss, _draw) = loopy_decision_sets(n, |v| {
        (1..n)
            .filter(|&d| polar(v as u128, d as u128, lam, a, m) == 1)
            .map(|d| v ^ d)
            .collect()
    });
    println!(
        "(i)  symmetric-B loopy Loss-set = R(B) = {:?}  ⇒ radical route empty (bent).\n",
        loss
    );

    // Sanity: the local energy change ΔQ_i(v) = q_i ⊕ B(v,e_i) = Q(v⊕e_i) ⊕ Q(v).
    let q_diag: Vec<u128> = (0..m).map(|i| gold(1 << i, lam, a, m)).collect();
    let delta = |v: u128, i: u128| q_diag[i as usize] ^ polar(v, 1 << i, lam, a, m);
    assert!((0..n as u128).all(
        |v| (0..m).all(|i| delta(v, i) == (gold(v ^ (1 << i), lam, a, m) ^ gold(v, lam, a, m)))
    ));
    println!("(ii) candidate rules (downward: turn OFF set bits only, so play terminates):\n");

    // Rule A — LOCAL SPIN-FLIP / Ising: turn off bit i iff the local energy changes,
    //          ΔQ_i(v) = q_i ⊕ B(v,e_i). Reads couplings B + per-coin field q_i only.
    let ra: Vec<Vec<usize>> = (0..n)
        .map(|v| {
            (0..m)
                .filter(|&i| v as u128 & (1 << i) != 0 && delta(v as u128, i) == 1)
                .map(|i| (v as u128 ^ (1 << i)) as usize)
                .collect()
        })
        .collect();
    let (pa, da) = p_set(&ra);
    describe("A local spin-flip (B+field)", &pa, &zero, da, m);

    // Rule B — B-only single-bit (the old Tier-1.5 baseline, no diagonal field).
    let rb: Vec<Vec<usize>> = (0..n)
        .map(|v| {
            (0..m)
                .filter(|&i| v as u128 & (1 << i) != 0 && polar(v as u128, 1 << i, lam, a, m) == 1)
                .map(|i| (v as u128 ^ (1 << i)) as usize)
                .collect()
        })
        .collect();
    let (pb, db) = p_set(&rb);
    describe("B single-bit B-only", &pb, &zero, db, m);

    // Rule C — B-coupled descent to any smaller w with B(v, v⊕w)=1.
    let rc: Vec<Vec<usize>> = (0..n)
        .map(|v| {
            (0..v)
                .filter(|&w| polar(v as u128, (v ^ w) as u128, lam, a, m) == 1)
                .collect()
        })
        .collect();
    let (pc, dc) = p_set(&rc);
    describe("C B-coupled descent", &pc, &zero, dc, m);

    // Rule D — Q-changing descent (tautological baseline: references global Q).
    let rd: Vec<Vec<usize>> = (0..n)
        .map(|v| {
            (0..v)
                .filter(|&w| gold(w as u128, lam, a, m) != gold(v as u128, lam, a, m))
                .collect()
        })
        .collect();
    let (pd, dd) = p_set(&rd);
    describe("D Q-changing (tautological)", &pd, &zero, dd, m);

    println!("\nReading (what the bent form reveals). Rule D (global Q) is the tautological");
    println!("baseline and hits {{Q=0}} exactly. The genuine results are B and A:");
    println!(" • Rule B reads ONLY the couplings B in the bit frame (no diagonal, no Q) and");
    println!("   already produces a genuine BENT quadric of the CORRECT Arf — the right");
    println!("   isometry class — but a DIFFERENT member of it (agreement at chance, 128/256).");
    println!("   So B+frame reaches the right kind of quadric; the residual gap to the specific");
    println!("   Gold {{Q=0}} is alignment within the O(Q)-orbit, i.e. the diagonal framing.");
    println!(" • Rule A is the natural Ising completion — the local spin-flip that ADDS the");
    println!("   per-coin field q_i to B. It does NOT align B's quadric to {{Q=0}}; it leaves");
    println!("   the quadric variety entirely. So the naive local-field assembly fails: the");
    println!("   diagonal framing must enter some other way than a per-coin spin-flip gate.");
    println!("This sharpens the open question on this bent case: B+frame can reach a");
    println!("right-Arf quadric class, but aligning to the specific Gold quadric remains open.");
}
