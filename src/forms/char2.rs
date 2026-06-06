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

use crate::clifford::{Metric, Multivector};
use crate::scalar::{nim_add, nim_inv, nim_mul, nim_trace, Nimber};

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

// ---------------------------------------------------------------------------
// The Dickson invariant — the characteristic-2 determinant replacement
// ---------------------------------------------------------------------------

/// Rank of a matrix over the nim-field F_{2^64}, by Gaussian elimination with nim
/// arithmetic. Rows are dense u64 vectors (all the same length).
fn nim_matrix_rank(mut rows: Vec<Vec<u64>>) -> usize {
    let nrows = rows.len();
    if nrows == 0 {
        return 0;
    }
    let ncols = rows[0].len();
    let mut pr = 0usize; // current pivot row
    for col in 0..ncols {
        let Some(p) = (pr..nrows).find(|&r| rows[r][col] != 0) else {
            continue;
        };
        rows.swap(pr, p);
        let inv = nim_inv(rows[pr][col]).unwrap();
        for c in col..ncols {
            rows[pr][c] = nim_mul(rows[pr][c], inv);
        }
        for r in 0..nrows {
            if r != pr && rows[r][col] != 0 {
                let f = rows[r][col];
                for c in col..ncols {
                    rows[r][c] = nim_add(rows[r][c], nim_mul(f, rows[pr][c]));
                }
            }
        }
        pr += 1;
        if pr == nrows {
            break;
        }
    }
    pr
}

/// The **Dickson invariant** `D(g) ∈ F₂` of an orthogonal transformation `g`,
/// given as an n×n matrix over a nim-field: `D(g) = dim Im(g − I) mod 2`
/// (`= rank(g + I) mod 2`, since `−1 = 1`).
///
/// In characteristic 2 the determinant of any `g ∈ O(Q)` is forced to `1`, so it
/// cannot separate rotations from reflections — the Dickson invariant is the
/// replacement, with `SO(Q) = ker D`. A single reflection has `D = 1`; a product
/// of `k` reflections has `D = k mod 2`. It is the companion to the Arf
/// invariant: **Arf classifies the form, Dickson classifies `O(Q)`.**
pub fn dickson_matrix(g: &[Vec<u64>]) -> u8 {
    let n = g.len();
    let mut m: Vec<Vec<u64>> = g.to_vec();
    for i in 0..n {
        m[i][i] = nim_add(m[i][i], 1); // g − I  (= g + I in char 2)
    }
    (nim_matrix_rank(m) % 2) as u8
}

/// The Dickson invariant of a Clifford **versor** (a product of vectors) acting
/// by the twisted adjoint: it is the ℤ₂-grade parity of the versor — an even
/// versor (rotor) lies in `SO` with `D = 0`, an odd versor (e.g. a single vector,
/// a reflection) has `D = 1`. Returns `None` if the multivector is not of
/// homogeneous grade parity (hence not a versor) or is zero.
pub fn dickson_of_versor(v: &Multivector<Nimber>) -> Option<u8> {
    // The Dickson invariant of a versor is its grade parity, a fact independent of
    // the scalar field — so this is the char-2 specialisation of the generic
    // `clifford::versor_grade_parity`.
    crate::clifford::versor_grade_parity(v)
}

// ---------------------------------------------------------------------------
// Fitting an F₂ quadratic form to a set — the "is this P-set a quadric?" test
// ---------------------------------------------------------------------------

/// The result of fitting a quadratic form to a subset of F₂^k.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct QuadricFit {
    /// Constant term: false ⇒ `0 ∈ set` (form through the origin); true ⇒ affine
    /// offset (`set = {Q = 1}` for the homogeneous part below).
    pub constant: bool,
    /// Diagonal q_i (the linear/`x_i` coefficients = squares over F₂).
    pub qd: Vec<bool>,
    /// Polar form bmat (the `x_i x_j` coefficients), as adjacency rows.
    pub bmat: Vec<u32>,
    /// Arf classification of the homogeneous quadratic part.
    pub arf: ArfResult,
}

impl QuadricFit {
    /// Whether the fitted form has genuine quadratic content (nonzero polar form
    /// rank). `false` ⇒ the set is an affine flat / linear condition, no quadratic
    /// refinement.
    pub fn is_genuinely_quadratic(&self) -> bool {
        self.arf.rank > 0
    }
}

/// F₂ scalar product of two coefficient vectors stored as `u64` bitmasks.
fn f2_dot(a: u64, b: u64) -> bool {
    (a & b).count_ones() & 1 == 1
}

/// Try to fit a quadratic form `Q(x) = c ⊕ Σ q_i x_i ⊕ Σ_{i<j} b_ij x_i x_j` over
/// F₂ on `k` variables whose zero set is exactly `set` (a list of bitmask points
/// of F₂^k). Returns `None` if no quadratic form has that zero set. Solved by
/// Gaussian elimination over F₂ on the `2^k` membership equations.
///
/// This is the instrument both game probes feed their P-positions into: it answers
/// "is this P-set a quadric, and if so what is its Arf (win-bias)?", and
/// distinguishes a genuine quadric (`is_genuinely_quadratic`) from a mere affine
/// subspace (the XOR-linear case normal play already produces).
pub fn fit_f2_quadratic(set: &[u32], k: usize) -> Option<QuadricFit> {
    assert!(k <= 12, "fit_f2_quadratic is exponential in k");
    // Coefficient layout: bit 0 = constant; bits 1..=k = linear x_i;
    // then one bit per pair (i<j) for the quadratic terms.
    let mut pair_index = vec![vec![0usize; k]; k];
    let mut nbits = 1 + k;
    for i in 0..k {
        for j in (i + 1)..k {
            pair_index[i][j] = nbits;
            nbits += 1;
        }
    }
    // Feature vector φ(v) over the coefficient layout (as a u64 bitmask).
    let phi = |v: u32| -> u64 {
        let mut f = 1u64; // constant
        for i in 0..k {
            if v & (1 << i) != 0 {
                f |= 1u64 << (1 + i);
            }
        }
        for i in 0..k {
            for j in (i + 1)..k {
                if v & (1 << i) != 0 && v & (1 << j) != 0 {
                    f |= 1u64 << pair_index[i][j];
                }
            }
        }
        f
    };
    let in_set: std::collections::HashSet<u32> = set.iter().copied().collect();

    // Build the augmented system: rows (φ(v) | target), target = 0 iff v ∈ set.
    let mut rows: Vec<(u64, bool)> = (0..(1u32 << k))
        .map(|v| (phi(v), !in_set.contains(&v)))
        .collect();

    // Gaussian elimination over F₂ (coefficient bits 0..nbits as pivots).
    let mut pivots: Vec<(usize, usize)> = Vec::new(); // (bit, row index)
    let mut r = 0usize;
    for bit in 0..nbits {
        if let Some(p) = (r..rows.len()).find(|&i| rows[i].0 & (1u64 << bit) != 0) {
            rows.swap(r, p);
            let (prow, ptgt) = rows[r];
            for i in 0..rows.len() {
                if i != r && rows[i].0 & (1u64 << bit) != 0 {
                    rows[i].0 ^= prow;
                    rows[i].1 ^= ptgt;
                }
            }
            pivots.push((bit, r));
            r += 1;
        }
    }
    // Consistency: any all-zero coefficient row must have target 0.
    if rows.iter().any(|&(coef, tgt)| coef == 0 && tgt) {
        return None; // not a quadric
    }
    // Read off one solution (free variables = 0).
    let mut sol = 0u64;
    for &(bit, row) in &pivots {
        if rows[row].1 {
            sol |= 1u64 << bit;
        }
    }
    // (Sanity is guaranteed by construction; the form below reproduces `set`.)
    let _ = f2_dot; // (kept for clarity of the dot-product convention)

    let constant = sol & 1 != 0;
    let qd: Vec<bool> = (0..k).map(|i| sol & (1u64 << (1 + i)) != 0).collect();
    let mut bmat = vec![0u32; k];
    for i in 0..k {
        for j in (i + 1)..k {
            if sol & (1u64 << pair_index[i][j]) != 0 {
                bmat[i] |= 1 << j;
                bmat[j] |= 1 << i;
            }
        }
    }
    let arf = arf_f2(k, &qd, &bmat);
    Some(QuadricFit {
        constant,
        qd,
        bmat,
        arf,
    })
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
        Metric::new(q, b)
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
    fn arf_additive_over_graded_tensor() {
        // The same A⊕A ≅ H⊕H fact, but built with the `direct_sum` *operation*
        // rather than a hand-written 4-generator metric: arf is additive over ⟂.
        let a = metric(&[1, 1], &b1(&[(0, 1)])); // anisotropic plane, Arf 1
        let h = metric(&[0, 0], &b1(&[(0, 1)])); // hyperbolic plane,  Arf 0
        let aa = arf_invariant(&a.direct_sum(&a));
        let hh = arf_invariant(&h.direct_sum(&h));
        let ah = arf_invariant(&a.direct_sum(&h));
        assert_eq!(aa.arf, 0); // 1 + 1 = 0
        assert_eq!(hh.arf, 0); // 0 + 0 = 0  ⇒  A⊕A ≅ H⊕H
        assert_eq!(ah.arf, 1); // 1 + 0 = 1
        assert_eq!((aa.rank, hh.rank, ah.rank), (4, 4, 4));
    }

    #[test]
    fn radical_is_detected() {
        // Q = x0 x1 + x2²: rank-2 core ⊕ a defective radical direction.
        let r = arf_invariant(&metric(&[0, 0, 1], &b1(&[(0, 1)])));
        assert_eq!(
            (r.rank, r.radical_dim, r.radical_anisotropic, r.arf),
            (2, 1, true, 0)
        );
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
    fn dickson_separates_rotations_from_reflections() {
        // identity is a rotation: D = 0.
        assert_eq!(dickson_matrix(&[vec![1, 0], vec![0, 1]]), 0);
        // the swap (0 1; 1 0) preserves the hyperbolic form x0 x1 and is a
        // reflection (odd): D = 1.
        assert_eq!(dickson_matrix(&[vec![0, 1], vec![1, 0]]), 1);
        // a hyperbolic "rotation" diag(t, t⁻¹) preserves x0 x1; for t=*2 in F₄,
        // t⁻¹ = *3, so g = diag(*2,*3): D = 0 (in SO).
        assert_eq!(dickson_matrix(&[vec![2, 0], vec![0, 3]]), 0);
        // composing two reflections (here swap∘swap = identity) gives D = 0.
        let swap = [[0u64, 1], [1, 0]];
        let mut comp = vec![vec![0u64; 2]; 2];
        for i in 0..2 {
            for j in 0..2 {
                let mut acc = 0u64;
                for k in 0..2 {
                    acc ^= nim_mul(swap[i][k], swap[k][j]);
                }
                comp[i][j] = acc;
            }
        }
        assert_eq!(dickson_matrix(&comp), 0);
    }

    #[test]
    fn dickson_of_versor_is_grade_parity() {
        use crate::clifford::{CliffordAlgebra, Metric};
        let alg = CliffordAlgebra::new(
            2,
            Metric::new(vec![Nimber(1), Nimber(1)], {
                let mut b = BTreeMap::new();
                b.insert((0, 1), Nimber(1));
                b
            }),
        );
        let scalar_one = alg.scalar(Nimber(1));
        let e0 = alg.gen(0);
        let e0e1 = alg.mul(&alg.gen(0), &alg.gen(1));
        assert_eq!(dickson_of_versor(&scalar_one), Some(0)); // identity rotor
        assert_eq!(dickson_of_versor(&e0), Some(1)); // a vector = a reflection
        assert_eq!(dickson_of_versor(&e0e1), Some(0)); // a bivector = a rotor
                                                       // mixed parity ⇒ not a versor
        let mixed = alg.add(&e0, &e0e1);
        assert_eq!(dickson_of_versor(&mixed), None);
    }

    // Evaluate a fitted form Q at v and return Q(v) ∈ {false,true}.
    fn eval_fit(fit: &QuadricFit, v: u32) -> bool {
        let mut acc = fit.constant;
        for i in 0..fit.qd.len() {
            if fit.qd[i] && v & (1 << i) != 0 {
                acc ^= true;
            }
        }
        for i in 0..fit.qd.len() {
            for j in (i + 1)..fit.qd.len() {
                if fit.bmat[i] & (1 << j) != 0 && v & (1 << i) != 0 && v & (1 << j) != 0 {
                    acc ^= true;
                }
            }
        }
        acc
    }

    #[test]
    fn fit_recovers_known_quadrics() {
        // hyperbolic Q = x0 x1: zero set {00,01,10}; genuine quadric, Arf 0.
        let h = fit_f2_quadratic(&[0, 1, 2], 2).unwrap();
        assert!(h.is_genuinely_quadratic());
        assert_eq!(h.arf.arf, 0);
        assert!(!h.constant);
        // anisotropic Q = x0²+x0x1+x1²: zero set {00}; Arf 1.
        let a = fit_f2_quadratic(&[0], 2).unwrap();
        assert!(a.is_genuinely_quadratic());
        assert_eq!(a.arf.arf, 1);
        // a LINEAR condition x0⊕x1=0: zero set {00,11}; a quadric but rank 0
        // (affine flat), i.e. NOT genuinely quadratic.
        let lin = fit_f2_quadratic(&[0, 3], 2).unwrap();
        assert!(!lin.is_genuinely_quadratic());
        assert_eq!(lin.arf.rank, 0);
    }

    #[test]
    fn quadric_count_and_roundtrip_over_f2_cubed() {
        // Over F₂³ there are 2^(1+3+3) = 128 quadratic forms but 2^8 = 256 subsets,
        // so exactly 128 subsets are quadrics — and each fit must reproduce its set.
        let mut count = 0;
        for s in 0u32..(1 << 8) {
            let set: Vec<u32> = (0..8u32).filter(|&v| s & (1 << v) != 0).collect();
            if let Some(fit) = fit_f2_quadratic(&set, 3) {
                count += 1;
                let recovered: Vec<u32> = (0..8u32).filter(|&v| !eval_fit(&fit, v)).collect();
                assert_eq!(recovered, set, "fit did not reproduce its own set");
            }
        }
        assert_eq!(count, 128, "expected exactly 2^7 quadrics over F₂³");
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
