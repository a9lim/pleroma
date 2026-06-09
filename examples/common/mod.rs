#![allow(dead_code)] // Each example binary imports a different subset of these helpers.

use pleroma::games::{outcomes, Outcome, Quotient};
use pleroma::scalar::{nim_add, nim_mul, nim_square, nim_trace};

fn frob(mut x: u128, a: u128) -> u128 {
    for _ in 0..a {
        x = nim_square(x);
    }
    x
}

/// Gold form `Q_a(v) = Tr(v^{1+2^a})` over `F_{2^m}`, valued in `{0,1}`.
pub fn gold(v: u128, a: u128, m: u128) -> u128 {
    nim_trace(nim_mul(v, frob(v, a)), m)
}

/// Bent Gold component `Tr(lambda * v^{1+2^a})`, valued in `{0,1}`.
pub fn bent_gold(v: u128, lambda: u128, a: u128, m: u128) -> u128 {
    nim_trace(nim_mul(lambda, nim_mul(v, frob(v, a))), m)
}

/// Polar form `B(u,v) = Q(u+v) + Q(u) + Q(v)`.
pub fn polar(u: u128, v: u128, a: u128, m: u128) -> u128 {
    gold(nim_add(u, v), a, m) ^ gold(u, a, m) ^ gold(v, a, m)
}

/// Polar form for a bent Gold component.
pub fn bent_polar(u: u128, v: u128, lambda: u128, a: u128, m: u128) -> u128 {
    bent_gold(nim_add(u, v), lambda, a, m) ^ bent_gold(u, lambda, a, m) ^ bent_gold(v, lambda, a, m)
}

pub fn p_set(succ: &[Vec<usize>]) -> (Vec<u128>, usize) {
    let out = outcomes(succ);
    let draws = out.iter().filter(|&&o| o == Outcome::Draw).count();
    let p = out
        .iter()
        .enumerate()
        .filter(|(_, o)| **o == Outcome::Loss)
        .map(|(i, _)| i as u128)
        .collect();
    (p, draws)
}

/// If `q` is an elementary abelian 2-group on `atoms`, return its P-set as
/// `F_2^k` bitmasks.
pub fn p_set_as_f2(q: &Quotient, atoms: &[usize]) -> Option<Vec<u128>> {
    let k = atoms.len();
    if k > 12 || q.num_classes != (1 << k) {
        return None;
    }
    let class_of_subset = |mask: u128| -> Option<usize> {
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
    let mut hit = std::collections::HashSet::new();
    let mut pset = Vec::new();
    for v in 0u128..(1 << k) {
        let c = class_of_subset(v)?;
        hit.insert(c);
        if q.class_is_p[c] {
            pset.push(v);
        }
    }
    (hit.len() == (1 << k)).then_some(pset)
}
