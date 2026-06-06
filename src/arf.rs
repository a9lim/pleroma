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
//! Vectors are u32 bitmasks over the (≤32) generators. This is the F₂ case;
//! a general nim-field form reduces to it via the field trace (not yet wired).

use crate::clifford::Metric;
use crate::nimber::Nimber;

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

/// Arf invariant of a nimber Clifford metric, provided every q/b entry is in
/// F₂ = {*0, *1}. Returns `None` if any entry is a higher nimber (that form
/// lives over a larger nim-field and reduces to F₂ only via the trace).
pub fn arf_invariant(metric: &Metric<Nimber>) -> Option<ArfResult> {
    let n = metric.q.len();
    if n > 32 {
        return None;
    }
    let mut qd = vec![false; n];
    for (i, slot) in qd.iter_mut().enumerate() {
        *slot = match metric.q[i].0 {
            0 => false,
            1 => true,
            _ => return None,
        };
    }
    let mut bmat = vec![0u32; n];
    for (&(i, j), v) in &metric.b {
        match v.0 {
            0 => {}
            1 => {
                bmat[i] |= 1 << j;
                bmat[j] |= 1 << i;
            }
            _ => return None,
        }
    }
    Some(arf_f2(n, &qd, &bmat))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeMap;

    fn metric(qs: &[u64], bs: &[(usize, usize)]) -> Metric<Nimber> {
        let q = qs.iter().map(|&x| Nimber(x)).collect();
        let mut b = BTreeMap::new();
        for &(i, j) in bs {
            b.insert((i, j), Nimber(1));
        }
        Metric { q, b }
    }

    #[test]
    fn hyperbolic_plane_is_o_plus() {
        // Q = x0 x1: a single hyperbolic pair, Arf 0.
        let r = arf_invariant(&metric(&[0, 0], &[(0, 1)])).unwrap();
        assert_eq!(r.arf, 0);
        assert_eq!(r.rank, 2);
        assert_eq!(r.radical_dim, 0);
        assert_eq!(r.o_type, "O+");
    }

    #[test]
    fn anisotropic_plane_is_o_minus() {
        // Q = x0² + x0 x1 + x1²: Arf 1.
        let r = arf_invariant(&metric(&[1, 1], &[(0, 1)])).unwrap();
        assert_eq!(r.arf, 1);
        assert_eq!(r.rank, 2);
        assert_eq!(r.o_type, "O-");
    }

    #[test]
    fn the_two_planes_are_distinguished() {
        let h = arf_invariant(&metric(&[0, 0], &[(0, 1)])).unwrap();
        let a = arf_invariant(&metric(&[1, 1], &[(0, 1)])).unwrap();
        assert_ne!(h.arf, a.arf); // exactly what classifies them
    }

    #[test]
    fn arf_is_additive_over_orthogonal_sum() {
        // H ⊕ H = O+,  H ⊕ A = O-,  A ⊕ A = O+  (two anisotropic planes ≅ two hyperbolic)
        let hh = arf_invariant(&metric(&[0, 0, 0, 0], &[(0, 1), (2, 3)])).unwrap();
        let ha = arf_invariant(&metric(&[0, 0, 1, 1], &[(0, 1), (2, 3)])).unwrap();
        let aa = arf_invariant(&metric(&[1, 1, 1, 1], &[(0, 1), (2, 3)])).unwrap();
        assert_eq!((hh.arf, hh.rank), (0, 4));
        assert_eq!((ha.arf, ha.rank), (1, 4));
        assert_eq!((aa.arf, aa.rank), (0, 4)); // A ⊕ A ≅ H ⊕ H
    }

    #[test]
    fn radical_is_detected() {
        // Q = x0 x1 + x2²: rank-2 core ⊕ a defective radical direction.
        let r = arf_invariant(&metric(&[0, 0, 1], &[(0, 1)])).unwrap();
        assert_eq!(r.rank, 2);
        assert_eq!(r.radical_dim, 1);
        assert!(r.radical_anisotropic);
        assert_eq!(r.arf, 0);
    }

    #[test]
    fn non_f2_entries_rejected() {
        assert!(arf_invariant(&metric(&[2], &[])).is_none()); // q0 = *2 ∉ F₂
    }
}
