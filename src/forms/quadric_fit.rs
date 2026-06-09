//! Fitting an F₂ quadratic form to a point set — the **"is this P-set a
//! quadric?"** test bench.
//!
//! This is not a classifier (that is [`char2`](crate::forms::char2)'s Arf): it is
//! the research instrument the game probes feed their P-positions into. Given a
//! subset `S ⊆ F₂^k`, [`fit_f2_quadratic`] decides whether `S` is the zero set of
//! *some* quadratic form and, if so, returns that form together with its
//! [Arf](crate::forms::ArfResult) — distinguishing a genuine quadric (nonzero
//! polar rank) from a mere affine flat (the XOR-linear case normal play already
//! produces). It is the bench behind the `misere_quotient` and `octal_hunt`
//! examples and the open-question probes; see root `OPEN.md`.

use crate::forms::{arf_f2, ArfResult};

/// The result of fitting a quadratic form to a subset of F₂^k.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct QuadricFit {
    /// Constant term: false ⇒ `0 ∈ set` (form through the origin); true ⇒ affine
    /// offset (`set = {Q = 1}` for the homogeneous part below).
    pub constant: bool,
    /// Diagonal q_i (the linear/`x_i` coefficients = squares over F₂).
    pub qd: Vec<bool>,
    /// Polar form bmat (the `x_i x_j` coefficients), as adjacency rows.
    pub bmat: Vec<u128>,
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

/// Try to fit a quadratic form `Q(x) = c ⊕ Σ q_i x_i ⊕ Σ_{i<j} b_ij x_i x_j` over
/// F₂ on `k` variables whose zero set is exactly `set` (a list of bitmask points
/// of F₂^k). Returns `None` if no quadratic form has that zero set. The unique
/// Boolean algebraic normal form is computed by a fast Mobius transform on the
/// truth table; fitting succeeds exactly when every coefficient of degree `> 2`
/// vanishes.
///
/// This is the instrument both game probes feed their P-positions into: it answers
/// "is this P-set a quadric, and if so what is its Arf (win-bias)?", and
/// distinguishes a genuine quadric ([`QuadricFit::is_genuinely_quadratic`]) from a
/// mere affine subspace (the XOR-linear case normal play already produces).
pub fn fit_f2_quadratic(set: &[u128], k: usize) -> Option<QuadricFit> {
    const MAX_ANF_DIM: usize = 20;
    assert!(
        k <= MAX_ANF_DIM,
        "fit_f2_quadratic is exponential in k; max supported k is {MAX_ANF_DIM}"
    );
    let n = 1usize << k;
    let domain_mask = if k == 0 { 0 } else { (1u128 << k) - 1 };
    assert!(
        set.iter().all(|&v| v & !domain_mask == 0),
        "fit_f2_quadratic received a point outside F_2^{k}"
    );

    // Truth table for the target function Q(v): zero on `set`, one off it.
    let mut coeffs = vec![true; n];
    for &v in set {
        coeffs[v as usize] = false;
    }

    // Mobius transform: truth table -> algebraic normal form coefficients.
    for i in 0..k {
        let bit = 1usize << i;
        for mask in 0..n {
            if mask & bit != 0 {
                coeffs[mask] ^= coeffs[mask ^ bit];
            }
        }
    }

    if coeffs
        .iter()
        .enumerate()
        .any(|(mask, &c)| c && mask.count_ones() > 2)
    {
        return None;
    }

    let constant = coeffs[0];
    let qd: Vec<bool> = (0..k).map(|i| coeffs[1usize << i]).collect();
    let mut bmat = vec![0u128; k];
    for i in 0..k {
        for j in (i + 1)..k {
            if coeffs[(1usize << i) | (1usize << j)] {
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

    // Evaluate a fitted form Q at v and return Q(v) ∈ {false,true}.
    fn eval_fit(fit: &QuadricFit, v: u128) -> bool {
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
    fn fit_supports_beyond_the_old_coefficient_layout_ceiling() {
        let set: Vec<u128> = (0..(1u128 << 16)).collect();
        let fit = fit_f2_quadratic(&set, 16).unwrap();
        assert_eq!(fit.qd.len(), 16);
        assert_eq!(fit.arf.rank, 0);
        assert!(!fit.constant);
    }

    #[test]
    fn quadric_count_and_roundtrip_over_f2_cubed() {
        // Over F₂³ there are 2^(1+3+3) = 128 quadratic forms but 2^8 = 256 subsets,
        // so exactly 128 subsets are quadrics — and each fit must reproduce its set.
        let mut count = 0;
        for s in 0u128..(1 << 8) {
            let set: Vec<u128> = (0..8u128).filter(|&v| s & (1 << v) != 0).collect();
            if let Some(fit) = fit_f2_quadratic(&set, 3) {
                count += 1;
                let recovered: Vec<u128> = (0..8u128).filter(|&v| !eval_fit(&fit, v)).collect();
                assert_eq!(recovered, set, "fit did not reproduce its own set");
            }
        }
        assert_eq!(count, 128, "expected exactly 2^7 quadrics over F₂³");
    }

    #[test]
    fn cubic_truth_table_is_rejected() {
        // Q(v) = x0 x1 x2 has zero set all points except 111; it is not quadratic.
        let set: Vec<u128> = (0..8u128).filter(|&v| v != 7).collect();
        assert!(fit_f2_quadratic(&set, 3).is_none());
    }
}
