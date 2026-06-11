//! The `IntegralForm` type and its basic arithmetic invariants.
//!
//! Covers: constructors, determinant (Bareiss), evenness, unimodularity,
//! positive-definiteness (Sylvester), signature, invariant factors (Smith),
//! level, Clifford metrics, and direct sum. The geometry
//! (short-vector enumeration, Fincke–Pohst, automorphism counting) lives in
//! [`super::geometry`].

use crate::forms::integral::diagonal::{
    rational_congruence_diagonal, signature_from_diagonal, DegenerateBehavior,
};
use crate::linalg::field::inverse_matrix;
use crate::linalg::integer::smith_normal_form;
use crate::scalar::{Nimber, Rational};
use std::collections::BTreeMap;

// ── small arithmetic helpers ──

pub(super) fn gcd_i128(a: i128, b: i128) -> i128 {
    let (mut a, mut b) = (a.abs(), b.abs());
    while b != 0 {
        let t = b;
        b = a % b;
        a = t;
    }
    a
}

pub(super) fn lcm_i128(a: i128, b: i128) -> i128 {
    if a == 0 || b == 0 {
        return 0;
    }
    let g = gcd_i128(a, b);
    (a / g)
        .checked_mul(b)
        .expect("lattice level exceeds i128")
        .abs()
}

// ── Bareiss determinant ──

/// Fraction-free (Bareiss) determinant of a square integer matrix — exact, no
/// rational intermediates. Overflow on the integer intermediates is the same
/// i128 limit the rest of the crate carries.
pub(super) fn bareiss_det(mut a: Vec<Vec<i128>>) -> i128 {
    let n = a.len();
    if n == 0 {
        return 1;
    }
    let mut sign = 1i128;
    let mut prev = 1i128;
    for k in 0..n - 1 {
        if a[k][k] == 0 {
            match (k + 1..n).find(|&r| a[r][k] != 0) {
                Some(r) => {
                    a.swap(k, r);
                    sign = -sign;
                }
                None => return 0,
            }
        }
        for i in k + 1..n {
            for j in k + 1..n {
                let p1 = a[i][j]
                    .checked_mul(a[k][k])
                    .expect("Bareiss determinant exceeds i128");
                let p2 = a[i][k]
                    .checked_mul(a[k][j])
                    .expect("Bareiss determinant exceeds i128");
                a[i][j] = (p1 - p2) / prev; // exact by the Bareiss identity
            }
        }
        prev = a[k][k];
    }
    sign * a[n - 1][n - 1]
}

// ── IntegralForm ──

/// A positive-definite or indefinite integral lattice, recorded by its symmetric integer
/// Gram matrix `G`. Construct with [`IntegralForm::new`] (validates square +
/// symmetric) or [`IntegralForm::diagonal`]; the Gram is kept private so the
/// symmetry invariant cannot be broken by a bare struct literal.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IntegralForm {
    pub(super) gram: Vec<Vec<i128>>,
}

impl IntegralForm {
    /// Build a lattice from a symmetric integer Gram matrix. Returns `None` if
    /// the matrix is not square or not symmetric.
    pub fn new(gram: Vec<Vec<i128>>) -> Option<Self> {
        let n = gram.len();
        if gram.iter().any(|row| row.len() != n) {
            return None;
        }
        for i in 0..n {
            for j in 0..n {
                if gram[i][j] != gram[j][i] {
                    return None;
                }
            }
        }
        Some(IntegralForm { gram })
    }

    /// The diagonal lattice `⟨d₀, d₁, …⟩` (an orthogonal sum of rank-1 forms).
    pub fn diagonal(diag: &[i128]) -> Self {
        let n = diag.len();
        let mut gram = vec![vec![0i128; n]; n];
        for (i, &d) in diag.iter().enumerate() {
            gram[i][i] = d;
        }
        IntegralForm { gram }
    }

    /// The rank `n` of the lattice.
    pub fn dim(&self) -> usize {
        self.gram.len()
    }

    /// The Gram matrix `G = (⟨eᵢ, eⱼ⟩)`.
    pub fn gram(&self) -> &[Vec<i128>] {
        &self.gram
    }

    /// The bilinear pairing `⟨x, y⟩ = xᵀ G y`.
    pub fn inner(&self, x: &[i128], y: &[i128]) -> i128 {
        let n = self.dim();
        debug_assert_eq!(x.len(), n);
        debug_assert_eq!(y.len(), n);
        let mut acc = 0i128;
        for i in 0..n {
            if x[i] == 0 {
                continue;
            }
            let mut row = 0i128;
            for j in 0..n {
                row = row
                    .checked_add(
                        self.gram[i][j]
                            .checked_mul(y[j])
                            .expect("lattice inner product exceeds i128"),
                    )
                    .expect("lattice inner product exceeds i128");
            }
            acc = acc
                .checked_add(x[i].checked_mul(row).expect("lattice norm exceeds i128"))
                .expect("lattice norm exceeds i128");
        }
        acc
    }

    /// The norm `Q(x) = xᵀ G x`.
    pub fn norm(&self, x: &[i128]) -> i128 {
        self.inner(x, x)
    }

    /// The determinant `det G` (Bareiss; exact). For a positive-definite lattice
    /// this is the squared covolume and is positive; `|det G|` is the order of
    /// the discriminant group `L# / L`.
    pub fn determinant(&self) -> i128 {
        bareiss_det(self.gram.clone())
    }

    /// `|det G| = 1`: the lattice is unimodular (`L# = L`, self-dual).
    pub fn is_unimodular(&self) -> bool {
        self.determinant().abs() == 1
    }

    /// Every diagonal Gram entry is even, i.e. `Q(x)` is even for all `x` (an
    /// *even* lattice). Off-diagonal symmetry already makes the cross terms
    /// `2⟨eᵢ, eⱼ⟩` even.
    pub fn is_even(&self) -> bool {
        (0..self.dim()).all(|i| self.gram[i][i] % 2 == 0)
    }

    /// Positive definiteness, via Sylvester's criterion: every leading principal
    /// minor is `> 0` (computed exactly with Bareiss).
    pub fn is_positive_definite(&self) -> bool {
        let n = self.dim();
        for k in 1..=n {
            let minor: Vec<Vec<i128>> = (0..k).map(|i| self.gram[i][..k].to_vec()).collect();
            if bareiss_det(minor) <= 0 {
                return false;
            }
        }
        true
    }

    /// The real signature `(p, q)`: positive and negative dimensions after exact
    /// rational congruence diagonalization. Degenerate directions, if any, are
    /// omitted from the pair.
    pub fn signature(&self) -> (usize, usize) {
        let diag = rational_congruence_diagonal(&self.gram, DegenerateBehavior::StopAtRadical);
        signature_from_diagonal(&diag)
    }

    /// The invariant factors `d₀ | d₁ | …` of the discriminant group (Smith
    /// normal form of `G`): `L# / L ≅ ⨁ ℤ/dᵢ`. For a nonsingular lattice the
    /// nonzero factors multiply to `|det G|`.
    pub fn invariant_factors(&self) -> Vec<i128> {
        smith_normal_form(self.gram.clone())
    }

    /// The **level** `N`: the smallest positive integer with `N·G⁻¹` an even
    /// integral matrix (integral, with even diagonal). Returns `None` if `G` is
    /// singular. An even unimodular lattice has level 1. For even lattices this
    /// equals the level of the modular form the theta series of `L` belongs to;
    /// for odd lattices this is the lattice-theoretic level only — the theta
    /// series of an odd lattice is not a standard modular form for this level.
    pub fn level(&self) -> Option<i128> {
        let n = self.dim();
        if n == 0 {
            return Some(1);
        }
        let mat: Vec<Vec<Rational>> = self
            .gram
            .iter()
            .map(|row| row.iter().map(|&x| Rational::int(x)).collect())
            .collect();
        let inv = inverse_matrix(mat)?;
        let mut level = 1i128;
        for i in 0..n {
            for j in 0..n {
                let e = &inv[i][j];
                let den = e.denom(); // > 0, coprime to numerator
                                     // N·(num/den) ∈ ℤ ⟺ den | N. On the diagonal also need it even:
                                     // (N/den)·num even, which forces a further factor of 2 when num is odd.
                let modulus = if i == j && e.numer() % 2 != 0 {
                    den.checked_mul(2).expect("lattice level exceeds i128")
                } else {
                    den
                };
                level = lcm_i128(level, modulus);
            }
        }
        Some(level)
    }

    /// The rational Clifford metric attached to the lattice bilinear form:
    /// `e_i^2 = G_ii` and `{e_i,e_j} = 2G_ij`.
    pub fn clifford_metric(&self) -> crate::clifford::Metric<Rational> {
        let n = self.dim();
        let q = (0..n).map(|i| Rational::int(self.gram[i][i])).collect();
        let mut b = BTreeMap::new();
        for i in 0..n {
            for j in (i + 1)..n {
                let v = self.gram[i][j]
                    .checked_mul(2)
                    .expect("lattice Clifford metric exceeds i128");
                if v != 0 {
                    b.insert((i, j), Rational::int(v));
                }
            }
        }
        crate::clifford::Metric::new(q, b)
    }

    /// The characteristic-2 quadratic refinement of an even lattice, reduced
    /// modulo 2 from `Q/2`: `q_i = G_ii/2 (mod 2)` and `b_ij = G_ij (mod 2)`.
    /// Returns `None` for odd lattices, where `Q/2` is not integral.
    pub fn clifford_metric_f2(&self) -> Option<crate::clifford::Metric<Nimber>> {
        if !self.is_even() {
            return None;
        }
        let n = self.dim();
        let q = (0..n)
            .map(|i| Nimber((self.gram[i][i] / 2).rem_euclid(2) as u128))
            .collect();
        let mut b = BTreeMap::new();
        for i in 0..n {
            for j in (i + 1)..n {
                let v = self.gram[i][j].rem_euclid(2) as u128;
                if v != 0 {
                    b.insert((i, j), Nimber(v));
                }
            }
        }
        Some(crate::clifford::Metric::new(q, b))
    }

    /// The orthogonal direct sum `L ⟂ M` (block-diagonal Gram).
    pub fn direct_sum(&self, other: &IntegralForm) -> IntegralForm {
        let n = self.dim();
        let m = other.dim();
        let mut gram = vec![vec![0i128; n + m]; n + m];
        for i in 0..n {
            for j in 0..n {
                gram[i][j] = self.gram[i][j];
            }
        }
        for i in 0..m {
            for j in 0..m {
                gram[n + i][n + j] = other.gram[i][j];
            }
        }
        IntegralForm { gram }
    }

    /// `G·x` as an integer vector.
    pub(super) fn matvec(&self, x: &[i128]) -> Vec<i128> {
        let n = self.dim();
        (0..n)
            .map(|i| {
                let mut acc = 0i128;
                for j in 0..n {
                    acc = acc
                        .checked_add(
                            self.gram[i][j]
                                .checked_mul(x[j])
                                .expect("lattice matvec exceeds i128"),
                        )
                        .expect("lattice matvec exceeds i128");
                }
                acc
            })
            .collect()
    }
}

pub(super) fn dot(a: &[i128], b: &[i128]) -> i128 {
    let mut acc = 0i128;
    for (&x, &y) in a.iter().zip(b) {
        acc = acc
            .checked_add(x.checked_mul(y).expect("lattice dot exceeds i128"))
            .expect("lattice dot exceeds i128");
    }
    acc
}
