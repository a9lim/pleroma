//! The interactive route: can a natural B-coupled game have P-set {Q=0}?
//!   cargo run --example interactive_kernel
//!
//! Normal-play disjunctive sums give XOR-linear (subspace) P-sets, so they can't
//! produce the Gold quadric {Q=0}. The escape is an *interactive* game whose move
//! graph is not a disjunctive sum. Two facts frame the search:
//!
//!  (a) ANY set S is the P-set of *some* acyclic game — trivially: send every
//!      v∉S to a fixed loss in S, and every v∈S only to wins. So existence is
//!      free; the open question is whether a *natural*, uniform, B-coupled rule
//!      does it.
//!
//!  (b) A natural candidate couples moves through the polar form B of Q (the
//!      exact obstruction identified earlier). We orient moves downward (to
//!      strictly smaller integers) so the game terminates, gate them by B or Q,
//!      compute the P-set with the kernel solver, and compare it to {Q=0} — using
//!      `fit_f2_quadratic` to say what the P-set actually *is*.
//!
//! This prints (a) as a sanity check and then the result of (b): the natural
//! rules tried, how close their P-sets come to {Q=0}, and what they are instead.

use pleroma::forms::fit_f2_quadratic;
use pleroma::games::{outcomes, Outcome};
use pleroma::scalar::{nim_add, nim_mul, nim_square, nim_trace};

/// Gold form Q_a(v) = Tr(v^{1+2^a}) over F_{2^m}, valued in {0,1}.
fn gold(v: u64, a: u32, m: u32) -> u64 {
    let mut g = v;
    for _ in 0..a {
        g = nim_square(g);
    }
    nim_trace(nim_mul(v, g), m)
}

/// Polar form B(u,v) = Q(u⊕v) ⊕ Q(u) ⊕ Q(v) ∈ {0,1}.
fn polar(u: u64, v: u64, a: u32, m: u32) -> u64 {
    gold(nim_add(u, v), a, m) ^ gold(u, a, m) ^ gold(v, a, m)
}

fn p_set(succ: &[Vec<usize>]) -> (Vec<u32>, usize) {
    let out = outcomes(succ);
    let draws = out.iter().filter(|&&o| o == Outcome::Draw).count();
    let p = out
        .iter()
        .enumerate()
        .filter(|(_, o)| **o == Outcome::Loss)
        .map(|(i, _)| i as u32)
        .collect();
    (p, draws)
}

fn agreement(p: &[u32], zero: &[u32], n: usize) -> usize {
    let ps: std::collections::HashSet<u32> = p.iter().copied().collect();
    let zs: std::collections::HashSet<u32> = zero.iter().copied().collect();
    (0..n as u32)
        .filter(|v| ps.contains(v) == zs.contains(v))
        .count()
}

fn describe_pset(label: &str, p: &[u32], zero: &[u32], draws: usize, m: u32) {
    let n = 1usize << m;
    let agree = agreement(p, zero, n);
    let plen = p.len();
    let bias = plen as i64 - (n as i64 / 2);
    print!(
        "  {label:<26} |P|={plen:<3} draws={draws:<3} agree {agree}/{n} with {{Q=0}}  bias={bias:+}"
    );
    if p == zero {
        println!("   ← EXACTLY {{Q=0}}!");
    } else {
        match fit_f2_quadratic(p, m as usize) {
            Some(f) if f.is_genuinely_quadratic() => {
                println!(
                    "   P-set is a quadric (Arf={}, rank={})",
                    f.arf.arf, f.arf.rank
                )
            }
            Some(_) => println!("   P-set is affine/linear (a subspace coset)"),
            None => println!("   P-set is not a quadric"),
        }
    }
}

fn main() {
    let (m, a) = (4u32, 1u32);
    let n = 1usize << m;
    let zero: Vec<u32> = (0..n as u64)
        .filter(|&v| gold(v, a, m) == 0)
        .map(|v| v as u32)
        .collect();
    println!(
        "Gold form Q_{a} on F_2^{m}:  |{{Q=0}}| = {} of {n}",
        zero.len()
    );

    // (a) existence is trivial: hand-build an acyclic game with P-set = {Q=0}.
    let zset: std::collections::HashSet<u32> = zero.iter().copied().collect();
    let adhoc: Vec<Vec<usize>> = (0..n)
        .map(|v| {
            let vv = v as u32;
            if zset.contains(&vv) {
                // a target loss: move only to non-{Q=0} (wins)
                (0..v).filter(|&w| !zset.contains(&(w as u32))).collect()
            } else {
                vec![0] // a win: move to 0 ∈ {Q=0} (a loss)
            }
        })
        .collect();
    let (adhoc_p, _) = p_set(&adhoc);
    println!(
        "\n(a) hand-built acyclic game reproduces {{Q=0}} exactly: {}",
        adhoc_p == zero
    );
    println!("    ⇒ existence is free; the open question is a NATURAL uniform rule.");

    // (b) natural uniform B-/Q-coupled descent rules (downward ⇒ terminating).
    println!("\n(b) natural uniform rules (move v→w only for w<v, so the game terminates):");

    // Rule 1: move legal iff the polar form pairs v with the flip direction.
    let r1: Vec<Vec<usize>> = (0..n)
        .map(|v| {
            (0..v)
                .filter(|&w| polar(v as u64, nim_add(v as u64, w as u64), a, m) == 1)
                .collect()
        })
        .collect();
    let (p1, d1) = p_set(&r1);
    describe_pset("B-coupled descent", &p1, &zero, d1, m);

    // Rule 2: move legal iff it changes the value of Q.
    let r2: Vec<Vec<usize>> = (0..n)
        .map(|v| {
            (0..v)
                .filter(|&w| gold(w as u64, a, m) != gold(v as u64, a, m))
                .collect()
        })
        .collect();
    let (p2, d2) = p_set(&r2);
    describe_pset("Q-changing descent", &p2, &zero, d2, m);

    // Rule 3: turn off one set bit, gated by B with the turned-off direction.
    let r3: Vec<Vec<usize>> = (0..n)
        .map(|v| {
            (0..m)
                .filter(|&i| v as u64 & (1 << i) != 0)
                .map(|i| (v as u64 ^ (1 << i)) as usize)
                .filter(|&w| polar(v as u64, nim_add(v as u64, w as u64), a, m) == 1)
                .collect()
        })
        .collect();
    let (p3, d3) = p_set(&r3);
    describe_pset("single-bit B-gated turn", &p3, &zero, d3, m);

    println!("\nConclusion. One uniform rule — 'Q-changing descent' — reproduces {{Q=0}} exactly,");
    println!("but ONLY because it references Q directly in the move legality (move iff you flip");
    println!("Q): it bakes the form into the rules, so it is tautological, not a discovery.");
    println!("The rules coupled through B — the polar form, which is the legitimately");
    println!("game-realizable (coin-turning) ingredient — do NOT give {{Q=0}}: B-coupled descent");
    println!("yields an affine subspace, and the single-bit B-gated turn yields a *different*");
    println!("quadric (wrong Arf). So the sharp open question stands: a game whose moves are");
    println!("built from the combinatorial ingredients (B / coin-turning) ALONE — not from Q");
    println!("itself — with P-set {{Q=0}}. The kernel solver + fit_f2_quadratic are the test");
    println!("bench; the gap is now precisely a B-only rule that integrates up to the Q-quadric.");
}
