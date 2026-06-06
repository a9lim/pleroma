//! The Arf invariant — the complete classifying invariant of a nonsingular
//! quadratic form over F₂, and (Bertram et al., "Real Clifford algebras and
//! quadratic forms over F₂", arXiv:1601.07664) the invariant that classifies
//! the characteristic-2 Clifford algebra the form defines: two such algebras
//! are isomorphic iff their F₂ forms share an Arf invariant.
//!
//! A nim-Clifford metric (q, b) restricted to F₂ entries *is* a quadratic form
//! over F₂ on the generator space: for x ∈ F₂ⁿ,
//!     Q(x) = Σ_i x_i q_i  +  Σ_{i<j} x_i x_j b_{ij}        (x_i² = x_i)
//! with polar form B(e_i,e_j) = b_{ij} (alternating, B(e_i,e_i)=0). We compute
//! a symplectic basis {a_k,b_k} for B (peeling hyperbolic pairs, leaving the
//! radical) and return Arf = Σ_k Q(a_k) Q(b_k) ∈ F₂.
//!
//! `arf_f2` is the F₂ case (u32 bitmask vectors over ≤32 generators).
//! `arf_nimber` handles a form over any nim-subfield F_{2^{2^k}}: symplectic
//! reduction over the field (normalising pairs with `nim_inv`), then the Arf
//! sum is pushed to F₂ by the field trace. `arf_invariant` uses the latter.

use crate::clifford::Metric;
use crate::nimber::{nim_add, nim_inv, nim_mul, Nimber};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ArfResult {
    /// Arf invariant of the nonsingular core (0 or 1).
    pub arf: u8,
    /// Rank of the polar form B = 2 × (number of hyperbolic pairs).
    pub rank: usize,
    /// Dimension of the polar-form radical (where B vanishes).
    pub radical_dim: usize,
    /// Whether Q is nonzero somewhere on the radical (a "defective" direction).
    pub radical_anisotropic: bool,
    /// Orthogonal type of the nonsingular core: "O+" (split) iff Arf=0.
    pub o_type: &'static str,
}

/// Bits of `mask` strictly above position `i`.
fn above(i: usize) -> u32 {
    if i >= 31 {
        0
    } else {
        (!0u32) << (i + 1)
    }
}

/// Q(v) for a bitmask vector v.
fn q_of(v: u32, qd: &[bool], bmat: &[u32]) -> bool {
    let mut acc = false;
    let mut vv = v;
    while vv != 0 {
        let i = vv.trailing_zeros() as usize;
        vv &= vv - 1;
        if qd[i] {
            acc ^= true;
        }
        // pairs (i, j), j>i, both in v, with b_{ij}=1
        let inter = bmat[i] & v & above(i);
        if inter.count_ones() & 1 == 1 {
            acc ^= true;
        }
    }
    acc
}

/// Polar form B(u, v) = Σ_{i≠j} u_i v_j b_{ij}.
fn b_of(u: u32, v: u32, bmat: &[u32]) -> bool {
    let mut acc = false;
    let mut uu = u;
    while uu != 0 {
        let i = uu.trailing_zeros() as usize;
        uu &= uu - 1;
        if (bmat[i] & v).count_ones() & 1 == 1 {
            acc ^= true;
        }
    }
    acc
}

/// Arf invariant of an F₂ quadratic form given by diagonal `qd` (the squares)
/// and symmetric adjacency `bmat` (the polar form; bmat[i] bit j ⇔ b_{ij}=1).
pub fn arf_f2(n: usize, qd: &[bool], bmat: &[u32]) -> ArfResult {
    let mut vectors: Vec<u32> = (0..n).map(|i| 1u32 << i).collect();
    let mut arf = false;
    let mut pairs = 0usize;
    let mut radical: Vec<u32> = Vec::new();

    while let Some(a) = vectors.pop() {
        if let Some(pos) = vectors.iter().position(|&w| b_of(a, w, bmat)) {
            let bb = vectors.swap_remove(pos);
            // make every remaining vector orthogonal to both a and bb
            for w in vectors.iter_mut() {
                let mut nw = *w;
                if b_of(*w, bb, bmat) {
                    nw ^= a;
                }
                if b_of(*w, a, bmat) {
                    nw ^= bb;
                }
                *w = nw;
            }
            if q_of(a, qd, bmat) && q_of(bb, qd, bmat) {
                arf ^= true;
            }
            pairs += 1;
        } else {
            radical.push(a); // orthogonal to everything ⇒ radical
        }
    }

    let radical_anisotropic = radical.iter().any(|&v| q_of(v, qd, bmat));
    ArfResult {
        arf: arf as u8,
        rank: 2 * pairs,
        radical_dim: radical.len(),
        radical_anisotropic,
        o_type: if arf { "O-" } else { "O+" },
    }
}

// ---------------------------------------------------------------------------
// General nim-field version (any On₂ entries, reduced to F₂ via the trace)
// ---------------------------------------------------------------------------

/// Smallest extension degree m = 2^k over F₂ such that the nim-subfield
/// F_{2^m} (the nimbers below 2^m) contains `max_val`.
fn min_field_degree(max_val: u64) -> u32 {
    let mut m = 1u32; // 2^k, starting k = 0  (F_2)
    loop {
        if m >= 64 {
            return 64;
        }
        if max_val < (1u64 << m) {
            return m;
        }
        m <<= 1;
    }
}

/// Field trace F_{2^m} → F₂:  Tr(x) = x + x² + x⁴ + … + x^{2^{m-1}} ∈ {0,1}.
/// The Arf invariant lives in k/℘(k); for a finite field this trace realises
/// the canonical isomorphism k/℘(k) ≅ F₂.
fn nim_trace(x: u64, m: u32) -> u64 {
    let mut acc = x;
    let mut t = x;
    for _ in 1..m {
        t = nim_mul(t, t);
        acc ^= t;
    }
    acc
}

fn vscale(c: u64, v: &[u64]) -> Vec<u64> {
    v.iter().map(|&x| nim_mul(c, x)).collect()
}
fn vadd(u: &[u64], v: &[u64]) -> Vec<u64> {
    u.iter().zip(v).map(|(&a, &b)| nim_add(a, b)).collect()
}

/// Q(v) = Σ_i v_i² q_i + Σ_{i<j} v_i v_j b_{ij}, over the nim-field.
fn qf(v: &[u64], q: &[u64], bmat: &[Vec<u64>]) -> u64 {
    let n = v.len();
    let mut acc = 0u64;
    for i in 0..n {
        acc ^= nim_mul(nim_mul(v[i], v[i]), q[i]);
        for j in (i + 1)..n {
            acc ^= nim_mul(nim_mul(v[i], v[j]), bmat[i][j]);
        }
    }
    acc
}

/// Polar form B(u,v) = Σ_{i<j} (u_i v_j + u_j v_i) b_{ij}, over the nim-field.
fn bf(u: &[u64], v: &[u64], bmat: &[Vec<u64>]) -> u64 {
    let n = u.len();
    let mut acc = 0u64;
    for i in 0..n {
        for j in (i + 1)..n {
            let cross = nim_add(nim_mul(u[i], v[j]), nim_mul(u[j], v[i]));
            acc ^= nim_mul(cross, bmat[i][j]);
        }
    }
    acc
}

/// Arf invariant of a nimber Clifford metric over its field of definition (the
/// smallest nim-subfield containing all entries), reduced to F₂ via the trace.
/// Works for any nimber metric — F₂ is the special case where the trace is the
/// identity. Symplectic reduction normalises each pair with `nim_inv`.
pub fn arf_nimber(metric: &Metric<Nimber>) -> ArfResult {
    let n = metric.q.len();
    let q: Vec<u64> = metric.q.iter().map(|x| x.0).collect();
    let mut bmat = vec![vec![0u64; n]; n];
    for (&(i, j), v) in &metric.b {
        bmat[i][j] = v.0;
        bmat[j][i] = v.0;
    }

    let mut maxv = q.iter().copied().max().unwrap_or(0);
    for row in &bmat {
        maxv = maxv.max(row.iter().copied().max().unwrap_or(0));
    }
    let m = min_field_degree(maxv);

    let mut vectors: Vec<Vec<u64>> = (0..n)
        .map(|i| {
            let mut e = vec![0u64; n];
            e[i] = 1;
            e
        })
        .collect();

    let mut s = 0u64; // Σ Q(a_k) Q(b_k), a field element
    let mut pairs = 0usize;
    let mut radical_dim = 0usize;
    let mut radical_anisotropic = false;

    while let Some(a) = vectors.pop() {
        if let Some(pos) = vectors.iter().position(|w| bf(&a, w, &bmat) != 0) {
            let braw = vectors.swap_remove(pos);
            let c = bf(&a, &braw, &bmat);
            let b = vscale(nim_inv(c).unwrap(), &braw); // rescale so B(a,b) = 1
            for w in vectors.iter_mut() {
                let wb = bf(w, &b, &bmat);
                let wa = bf(w, &a, &bmat);
                let mut nw = w.clone();
                if wb != 0 {
                    nw = vadd(&nw, &vscale(wb, &a));
                }
                if wa != 0 {
                    nw = vadd(&nw, &vscale(wa, &b));
                }
                *w = nw;
            }
            s ^= nim_mul(qf(&a, &q, &bmat), qf(&b, &q, &bmat));
            pairs += 1;
        } else {
            radical_dim += 1;
            if qf(&a, &q, &bmat) != 0 {
                radical_anisotropic = true;
            }
        }
    }

    let arf = nim_trace(s, m) as u8;
    ArfResult {
        arf,
        rank: 2 * pairs,
        radical_dim,
        radical_anisotropic,
        o_type: if arf == 1 { "O-" } else { "O+" },
    }
}

/// Arf invariant of a nimber Clifford metric (the char-2 Clifford classifier).
pub fn arf_invariant(metric: &Metric<Nimber>) -> ArfResult {
    arf_nimber(metric)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeMap;

    fn metric(qs: &[u64], bs: &[((usize, usize), u64)]) -> Metric<Nimber> {
        let q = qs.iter().map(|&x| Nimber(x)).collect();
        let mut b = BTreeMap::new();
        for &((i, j), v) in bs {
            b.insert((i, j), Nimber(v));
        }
        Metric { q, b }
    }
    fn b1(pairs: &[(usize, usize)]) -> Vec<((usize, usize), u64)> {
        pairs.iter().map(|&p| (p, 1)).collect()
    }

    #[test]
    fn hyperbolic_plane_is_o_plus() {
        // Q = x0 x1: a single hyperbolic pair, Arf 0.
        let r = arf_invariant(&metric(&[0, 0], &b1(&[(0, 1)])));
        assert_eq!((r.arf, r.rank, r.radical_dim, r.o_type), (0, 2, 0, "O+"));
    }

    #[test]
    fn anisotropic_plane_is_o_minus() {
        // Q = x0² + x0 x1 + x1²: Arf 1.
        let r = arf_invariant(&metric(&[1, 1], &b1(&[(0, 1)])));
        assert_eq!((r.arf, r.rank, r.o_type), (1, 2, "O-"));
    }

    #[test]
    fn the_two_planes_are_distinguished() {
        let h = arf_invariant(&metric(&[0, 0], &b1(&[(0, 1)])));
        let a = arf_invariant(&metric(&[1, 1], &b1(&[(0, 1)])));
        assert_ne!(h.arf, a.arf); // exactly what classifies them
    }

    #[test]
    fn arf_is_additive_over_orthogonal_sum() {
        // H⊕H = O+,  H⊕A = O-,  A⊕A = O+  (two anisotropic planes ≅ two hyperbolic)
        let hh = arf_invariant(&metric(&[0, 0, 0, 0], &b1(&[(0, 1), (2, 3)])));
        let ha = arf_invariant(&metric(&[0, 0, 1, 1], &b1(&[(0, 1), (2, 3)])));
        let aa = arf_invariant(&metric(&[1, 1, 1, 1], &b1(&[(0, 1), (2, 3)])));
        assert_eq!((hh.arf, hh.rank), (0, 4));
        assert_eq!((ha.arf, ha.rank), (1, 4));
        assert_eq!((aa.arf, aa.rank), (0, 4)); // A⊕A ≅ H⊕H
    }

    #[test]
    fn radical_is_detected() {
        // Q = x0 x1 + x2²: rank-2 core ⊕ a defective radical direction.
        let r = arf_invariant(&metric(&[0, 0, 1], &b1(&[(0, 1)])));
        assert_eq!((r.rank, r.radical_dim, r.radical_anisotropic, r.arf), (2, 1, true, 0));
    }

    #[test]
    fn f4_forms_via_trace() {
        // Genuine F₄ forms (entries up to *3), hand-computed via the trace:
        //   q=[*2,*3], b01=*1:  S = *2⊗*3 = *1,  Tr_{F₄/F₂}(*1) = *1+*1 = 0  ⇒ O+
        let r1 = arf_invariant(&metric(&[2, 3], &b1(&[(0, 1)])));
        assert_eq!((r1.arf, r1.o_type, r1.rank), (0, "O+", 2));
        //   q=[*2,*2], b01=*1:  S = *2⊗*2 = *3,  Tr(*3) = *3+*2 = *1       ⇒ O-
        let r2 = arf_invariant(&metric(&[2, 2], &b1(&[(0, 1)])));
        assert_eq!((r2.arf, r2.o_type, r2.rank), (1, "O-", 2));
    }

    #[test]
    fn general_agrees_with_f2_bitmask() {
        // The general nim-field path must match the F₂ bitmask version on every
        // F₂ form (arf, rank, radical_dim, anisotropy, type all invariant).
        let cases: &[(&[u64], &[(usize, usize)])] = &[
            (&[0, 0], &[(0, 1)]),
            (&[1, 1], &[(0, 1)]),
            (&[0, 0, 1], &[(0, 1)]),
            (&[1, 0, 1, 1], &[(0, 1), (2, 3)]),
            (&[1, 1, 1, 1, 0], &[(0, 1), (2, 3)]),
        ];
        for (qs, ps) in cases {
            let general = arf_nimber(&metric(qs, &b1(ps)));
            let n = qs.len();
            let qd: Vec<bool> = qs.iter().map(|&x| x == 1).collect();
            let mut bmat = vec![0u32; n];
            for &(i, j) in *ps {
                bmat[i] |= 1 << j;
                bmat[j] |= 1 << i;
            }
            assert_eq!(general, arf_f2(n, &qd, &bmat), "mismatch on q={:?}", qs);
        }
    }
}
