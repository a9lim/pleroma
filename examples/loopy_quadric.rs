//! The loopy route to the open question — and what it reveals.
//!   cargo run --example loopy_quadric
//!
//! `interactive_kernel` orients moves strictly downward so the game terminates —
//! an acyclic graph, hence no Draws, a Win/Loss bit per position. A *loopy* rule
//! keeps both flip directions, so positions sit on cycles and a third outcome
//! appears: **Draw**. OPEN.md's Tier-2 obstruction is that normal-play P-sets are
//! XOR-linear; the Draw-set of a cyclic rule is a new degree of freedom, not bound
//! by that linearity. `loopy_decision_sets` exposes both the Loss-set and the
//! Draw-set; `fit_f2_quadratic` names each.
//!
//! What the sweep actually finds is instructive, and honest: a *symmetric* B-only
//! rule (B is symmetric, so the move graph is undirected) collapses to detecting
//! the **radical of B**. In an undirected loopy graph the only Losses are the
//! isolated vertices — and `v` is isolated exactly when `B(v,·) ≡ 0`, i.e.
//! `v ∈ R(B)`. So Loss-set = R(B) and Draw-set = its complement, regardless of Q.
//! At `(m,a)=(4,1)` it *happens* that `R(B) = {Q=0}` (both 4 points), which looks
//! like a Tier-2 hit but is a coincidence of the radical — it breaks at `m=8`. And
//! R(B) is exactly the degenerate part where the frame-blind no-go is *silent*. So
//! the loopy B-only rule confirms the obstruction from a new angle rather than
//! breaking it; a genuine witness must hit `{Q=0}` where it is NOT the radical.

use ogdoad::forms::{fit_f2_quadratic, QuadricFit};
use ogdoad::games::loopy_decision_sets;

mod common;
use common::{gold, polar};

/// The radical R(B) = { v : B(v,d) = 0 for every direction d }.
fn radical(a: u128, m: u128) -> Vec<u128> {
    let n = 1u128 << m;
    (0..n)
        .filter(|&v| (0..n).all(|d| polar(v, d, a, m) == 0))
        .collect()
}

fn name_fit(fit: &Option<QuadricFit>) -> String {
    match fit {
        Some(f) if f.is_genuinely_quadratic() => {
            format!("quadric (Arf={}, rank={})", f.arf.arf, f.arf.rank)
        }
        Some(_) => "affine/linear (subspace coset)".to_string(),
        None => "not a quadric".to_string(),
    }
}

fn run(m: u128, a: u128) {
    let k = m as usize;
    let n = 1usize << k;
    let zero: Vec<u128> = (0..n as u128).filter(|&v| gold(v, a, m) == 0).collect();
    let rad = radical(a, m);

    // Symmetric B-coupling: move v → v⊕d for every direction d with B(v,d)=1. B is
    // symmetric so the graph is undirected (hence loopy).
    let (loss, draw) = loopy_decision_sets(n, |v| {
        (1..n)
            .filter(|&d| polar(v as u128, d as u128, a, m) == 1)
            .map(|d| v ^ d)
            .collect()
    });
    let loss_u: Vec<u128> = loss.iter().map(|&v| v as u128).collect();
    let draw_u: Vec<u128> = draw.iter().map(|&v| v as u128).collect();

    println!(
        "F_2^{m}, Gold Q_{a}:  |{{Q=0}}|={}  |R(B)|={}",
        zero.len(),
        rad.len()
    );
    println!(
        "  symmetric B-coupling:  |Loss|={}  |Draw|={}",
        loss.len(),
        draw.len()
    );
    println!(
        "    Loss-set = R(B)? {}     Loss-set = {{Q=0}}? {}",
        loss_u == rad,
        loss_u == zero
    );
    println!(
        "    Loss-set is a {}",
        name_fit(&fit_f2_quadratic(&loss_u, k))
    );
    println!(
        "    Draw-set is a {}",
        name_fit(&fit_f2_quadratic(&draw_u, k))
    );
}

fn main() {
    println!("Loopy (cyclic) B-only rules — both flip directions, so Draws appear.\n");
    println!("(m,a)=(4,1): the radical coincidentally equals {{Q=0}} —");
    run(4, 1);
    println!("\n(m,a)=(8,1): the coincidence breaks — Loss is still R(B), not {{Q=0}} —");
    run(8, 1);

    println!("\nConclusion. The symmetric B-only loopy rule detects R(B), the radical — at");
    println!("(4,1) that equals {{Q=0}} (4 points each), at (8,1) it does not (|R(B)|=4 vs");
    println!("|{{Q=0}}|=112). R(B) is precisely where the Sp(B) frame-blind no-go is silent, so");
    println!("the loopy instrument reproduces the obstruction rather than escaping it. The open");
    println!("question stands; the Draw-set is now a first-class target the probes can sweep.");
}
