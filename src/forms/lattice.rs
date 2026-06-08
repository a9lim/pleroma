//! Integral lattices: the ℤ-Gram-matrix view of a quadratic form.
//!
//! The forms pillar elsewhere classifies a quadratic form *over a field* (by its
//! square classes / Witt class / Arf invariant). An **integral lattice** is the
//! complementary object: a free ℤ-module `L ≅ ℤⁿ` with an integer-valued
//! symmetric bilinear form, recorded by its Gram matrix `G = (⟨eᵢ, eⱼ⟩)`. Its
//! invariants are arithmetic, not just field-theoretic — the determinant, the
//! level, the minimum and kissing number, the automorphism group order — and the
//! coarse classification is the **genus** (local equivalence at every place),
//! built on the same p-adic primitives `padic.rs`/`adelic.rs` already carry. This
//! module is the M1 core (the geometry of one lattice); `root_lattices.rs`,
//! `genus.rs`, and `mass_formula.rs` build the A/D/E catalogue, the genus
//! equivalence, and the Conway–Sloane mass formula on top.
//!
//! Conventions. The **norm** of `x ∈ L` is `Q(x) = xᵀ G x` (so a "norm-2 vector"
//! has `Q = 2`, matching the root-lattice literature; this is twice the value of
//! the associated quadratic form `½Q` when the lattice is even). The geometric
//! routines — [`IntegralForm::minimum`], [`minimal_vectors`](IntegralForm::minimal_vectors),
//! [`kissing_number`](IntegralForm::kissing_number),
//! [`automorphism_group_order`](IntegralForm::automorphism_group_order) — assume the
//! lattice is **positive definite** and return `None` otherwise (an indefinite
//! lattice has infinitely many vectors of every norm and an infinite
//! automorphism group). Vectors are reported in lattice (basis) coordinates as
//! integer vectors, both signs included.
//!
//! Honest cutoff. Short-vector enumeration first applies a conservative unimodular
//! size-reduction pass (integral shears/swaps, so the lattice is unchanged), then
//! runs Fincke–Pohst (an LDL-bounded box search with exact norm filtering) and maps
//! the vectors back to the original coordinates. Automorphism counting first checks
//! closed-form families: diagonal signed-permutation lattices, literal `A`/`D`/`E`
//! Cartan bases, and then basis-independent root systems recovered from the norm-2
//! roots. Everything else falls back to a backtracking search over basis images,
//! which is **exponential** in general. The fallback is bounded by an explicit node
//! budget ([`AUTO_NODE_BUDGET`]); when the search exceeds it the count is reported as
//! `None` rather than silently truncated. Use
//! [`automorphism_group_order_bounded`](IntegralForm::automorphism_group_order_bounded)
//! to choose the budget explicitly.

use crate::linalg::field::inverse_matrix;
use crate::linalg::integer::smith_normal_form;
use crate::scalar::Rational;
use std::collections::{BTreeMap, VecDeque};

/// The default node budget for [`IntegralForm::automorphism_group_order`]. Beyond
/// this many backtracking nodes the search reports `None` (the lattice is too
/// large for brute-force automorphism enumeration — e.g. `E₈`, whose Weyl group
/// has order ~7·10⁸, or the Leech lattice). The bound is explicit, not silent.
pub const AUTO_NODE_BUDGET: u64 = 100_000_000;

/// A positive-definite or indefinite integral lattice, recorded by its symmetric integer
/// Gram matrix `G`. Construct with [`IntegralForm::new`] (validates square +
/// symmetric) or [`IntegralForm::diagonal`]; the Gram is kept private so the
/// symmetry invariant cannot be broken by a bare struct literal.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IntegralForm {
    gram: Vec<Vec<i128>>,
}

fn gcd_i128(a: i128, b: i128) -> i128 {
    let (mut a, mut b) = (a.abs(), b.abs());
    while b != 0 {
        let t = b;
        b = a % b;
        a = t;
    }
    a
}

fn lcm_i128(a: i128, b: i128) -> i128 {
    if a == 0 || b == 0 {
        return 0;
    }
    let g = gcd_i128(a, b);
    (a / g)
        .checked_mul(b)
        .expect("lattice level exceeds i128")
        .abs()
}

fn checked_factorial(n: usize) -> Option<u128> {
    let mut acc = 1u128;
    for k in 2..=n {
        acc = acc.checked_mul(k as u128)?;
    }
    Some(acc)
}

fn checked_pow2(n: usize) -> Option<u128> {
    1u128.checked_shl(u32::try_from(n).ok()?)
}

fn signed_permutation_order(n: usize) -> Option<u128> {
    checked_pow2(n)?.checked_mul(checked_factorial(n)?)
}

fn round_div_nearest(num: i128, den: i128) -> i128 {
    debug_assert!(den > 0);
    let q = num.div_euclid(den);
    let r = num.rem_euclid(den);
    if r.checked_mul(2).expect("rounding residue exceeds i128") >= den {
        q + 1
    } else {
        q
    }
}

fn identity_i128(n: usize) -> Vec<Vec<i128>> {
    let mut out = vec![vec![0i128; n]; n];
    for (i, row) in out.iter_mut().enumerate() {
        row[i] = 1;
    }
    out
}

fn map_coords(u: &[Vec<i128>], y: &[i128]) -> Vec<i128> {
    let n = y.len();
    let mut out = vec![0i128; n];
    for i in 0..n {
        let mut acc = 0i128;
        for (j, &yj) in y.iter().enumerate() {
            acc = acc
                .checked_add(u[i][j].checked_mul(yj).expect("basis map exceeds i128"))
                .expect("basis map exceeds i128");
        }
        out[i] = acc;
    }
    out
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum RootComponentKind {
    A(usize),
    D(usize),
    E(usize),
}

impl RootComponentKind {
    fn from_rank_and_roots(rank: usize, roots: usize) -> Option<Self> {
        if rank >= 1 && roots == rank.checked_mul(rank + 1)? {
            return Some(RootComponentKind::A(rank));
        }
        if rank >= 4 && roots == 2usize.checked_mul(rank)?.checked_mul(rank - 1)? {
            return Some(RootComponentKind::D(rank));
        }
        match (rank, roots) {
            (6, 72) => Some(RootComponentKind::E(6)),
            (7, 126) => Some(RootComponentKind::E(7)),
            (8, 240) => Some(RootComponentKind::E(8)),
            _ => None,
        }
    }

    fn automorphism_order(self) -> Option<u128> {
        match self {
            RootComponentKind::A(n) => {
                if n == 1 {
                    Some(2)
                } else {
                    checked_factorial(n + 1)?.checked_mul(2)
                }
            }
            RootComponentKind::D(n) => match n {
                4 => Some(1152),
                _ if n >= 5 => checked_pow2(n)?.checked_mul(checked_factorial(n)?),
                _ => None,
            },
            RootComponentKind::E(6) => Some(103_680),
            RootComponentKind::E(7) => Some(2_903_040),
            RootComponentKind::E(8) => Some(696_729_600),
            RootComponentKind::E(_) => None,
        }
    }
}

fn canonical_root(mut v: Vec<i128>) -> Vec<i128> {
    if let Some(&first) = v.iter().find(|&&x| x != 0) {
        if first < 0 {
            for x in &mut v {
                *x = -*x;
            }
        }
    }
    v
}

fn rows_generate_full_lattice(rows: &[Vec<i128>], n: usize) -> bool {
    let hnf = crate::linalg::integer::normalize_relation_rows(rows.to_vec());
    if hnf.len() != n {
        return false;
    }
    let mut index = 1i128;
    for (i, row) in hnf.iter().enumerate() {
        index = index
            .checked_mul(row[i].abs())
            .expect("root-lattice index exceeds i128");
    }
    index == 1
}

fn simple_laced_cartan_matches(gram: &[Vec<i128>], edges: &[(usize, usize)]) -> bool {
    let n = gram.len();
    if gram.iter().any(|row| row.len() != n) {
        return false;
    }
    let mut adjacent = vec![vec![false; n]; n];
    for &(a, b) in edges {
        if a >= n || b >= n || a == b {
            return false;
        }
        adjacent[a][b] = true;
        adjacent[b][a] = true;
    }
    for i in 0..n {
        for j in 0..n {
            let expected = if i == j {
                2
            } else if adjacent[i][j] {
                -1
            } else {
                0
            };
            if gram[i][j] != expected {
                return false;
            }
        }
    }
    true
}

/// Fraction-free (Bareiss) determinant of a square integer matrix — exact, no
/// rational intermediates. Overflow on the integer intermediates is the same
/// i128 limit the rest of the crate carries.
fn bareiss_det(mut a: Vec<Vec<i128>>) -> i128 {
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

    /// The invariant factors `d₀ | d₁ | …` of the discriminant group (Smith
    /// normal form of `G`): `L# / L ≅ ⨁ ℤ/dᵢ`. For a nonsingular lattice the
    /// nonzero factors multiply to `|det G|`.
    pub fn invariant_factors(&self) -> Vec<i128> {
        smith_normal_form(self.gram.clone())
    }

    /// The **level** `N`: the smallest positive integer with `N·G⁻¹` an even
    /// integral matrix (integral, with even diagonal). This is the level of the
    /// modular form the theta series of `L` belongs to. Returns `None` if `G` is
    /// singular. An even unimodular lattice has level 1.
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

    // ---- positive-definite geometry (Fincke–Pohst + backtracking) ----

    /// The LDLᵀ decomposition in floating point: returns `(d, u)` with `d[i] > 0`
    /// and an upper unit factor `u[i][j]` (`j > i`) such that
    /// `Q(x) = Σᵢ d[i]·(x[i] + Σ_{j>i} u[i][j]·x[j])²`. The geometric search uses
    /// this only to *bound* the integer box; every reported vector is then
    /// checked with the exact integer norm, so float error cannot admit or omit a
    /// vector.
    fn ldl(&self) -> (Vec<f64>, Vec<Vec<f64>>) {
        let n = self.dim();
        let mut d = vec![0.0f64; n];
        let mut l = vec![vec![0.0f64; n]; n]; // unit lower triangular
        for j in 0..n {
            let mut dj = self.gram[j][j] as f64;
            for k in 0..j {
                dj -= l[j][k] * l[j][k] * d[k];
            }
            d[j] = dj;
            l[j][j] = 1.0;
            for i in j + 1..n {
                let mut s = self.gram[i][j] as f64;
                for k in 0..j {
                    s -= l[i][k] * l[j][k] * d[k];
                }
                l[i][j] = if dj != 0.0 { s / dj } else { 0.0 };
            }
        }
        // u[i][j] = L[j][i] for j > i
        let mut u = vec![vec![0.0f64; n]; n];
        for i in 0..n {
            for j in i + 1..n {
                u[i][j] = l[j][i];
            }
        }
        (d, u)
    }

    fn shear_basis(
        gram: &mut [Vec<i128>],
        transform: &mut [Vec<i128>],
        i: usize,
        j: usize,
        k: i128,
    ) {
        if k == 0 {
            return;
        }
        let n = gram.len();
        let gij = gram[i][j];
        let gii = gram[i][i];
        let new_jj = gram[j][j]
            .checked_add(
                k.checked_mul(2)
                    .and_then(|x| x.checked_mul(gij))
                    .expect("basis reduction exceeds i128"),
            )
            .and_then(|x| {
                k.checked_mul(k)
                    .and_then(|kk| kk.checked_mul(gii))
                    .and_then(|term| x.checked_add(term))
            })
            .expect("basis reduction exceeds i128");
        let mut new_col = vec![0i128; n];
        for l in 0..n {
            new_col[l] = gram[l][j]
                .checked_add(
                    k.checked_mul(gram[l][i])
                        .expect("basis reduction exceeds i128"),
                )
                .expect("basis reduction exceeds i128");
        }
        for l in 0..n {
            gram[l][j] = new_col[l];
            gram[j][l] = new_col[l];
        }
        gram[j][j] = new_jj;
        for row in transform {
            row[j] = row[j]
                .checked_add(k.checked_mul(row[i]).expect("basis map exceeds i128"))
                .expect("basis map exceeds i128");
        }
    }

    fn swap_basis(gram: &mut [Vec<i128>], transform: &mut [Vec<i128>], i: usize, j: usize) {
        if i == j {
            return;
        }
        gram.swap(i, j);
        for row in gram.iter_mut() {
            row.swap(i, j);
        }
        for row in transform {
            row.swap(i, j);
        }
    }

    /// A conservative integral size reduction of a positive-definite basis. The
    /// returned transform `U` maps reduced-basis coordinates back to `self`'s
    /// coordinates and has determinant `±1`.
    fn size_reduced_basis(&self) -> (IntegralForm, Vec<Vec<i128>>) {
        let n = self.dim();
        let mut gram = self.gram.clone();
        let mut transform = identity_i128(n);
        let max_passes = 8 * n.saturating_mul(n).saturating_add(1);
        for _ in 0..max_passes {
            let mut changed = false;
            for i in 0..n {
                if gram[i][i] <= 0 {
                    continue;
                }
                for j in i + 1..n {
                    let k = -round_div_nearest(gram[i][j], gram[i][i]);
                    if k != 0 {
                        Self::shear_basis(&mut gram, &mut transform, i, j, k);
                        changed = true;
                    }
                }
            }
            for i in 0..n.saturating_sub(1) {
                if gram[i + 1][i + 1] < gram[i][i] {
                    Self::swap_basis(&mut gram, &mut transform, i, i + 1);
                    changed = true;
                }
            }
            if !changed {
                break;
            }
        }
        (IntegralForm { gram }, transform)
    }

    /// All nonzero lattice vectors `x` with `0 < Q(x) ≤ bound`, in lattice
    /// coordinates and including both signs. `None` if the lattice is not
    /// positive definite (the count would be infinite).
    pub fn short_vectors(&self, bound: i128) -> Option<Vec<Vec<i128>>> {
        if !self.is_positive_definite() {
            return None;
        }
        if self.dim() == 0 || bound <= 0 {
            return Some(Vec::new());
        }
        let (reduced, transform) = self.size_reduced_basis();
        let vecs = reduced.short_vectors_raw(bound)?;
        Some(
            vecs.into_iter()
                .map(|v| map_coords(&transform, &v))
                .collect(),
        )
    }

    fn short_vectors_raw(&self, bound: i128) -> Option<Vec<Vec<i128>>> {
        if !self.is_positive_definite() {
            return None;
        }
        let n = self.dim();
        if n == 0 || bound <= 0 {
            return Some(Vec::new());
        }
        let (d, u) = self.ldl();
        let mut out = Vec::new();
        let mut x = vec![0i128; n];
        // Pad the float radius outward so rounding can only over-collect; the
        // exact integer filter at the leaf removes any spurious vectors.
        let eps = 1e-9 * (bound as f64).max(1.0) + 1e-9;
        self.fp_search(n, bound, &d, &u, eps, 0.0, &mut x, &mut out);
        Some(out)
    }

    #[allow(clippy::too_many_arguments)]
    fn fp_search(
        &self,
        i: usize,
        bound: i128,
        d: &[f64],
        u: &[Vec<f64>],
        eps: f64,
        tail: f64,
        x: &mut [i128],
        out: &mut Vec<Vec<i128>>,
    ) {
        if i == 0 {
            let q = self.norm(x);
            if q > 0 && q <= bound {
                out.push(x.to_vec());
            }
            return;
        }
        let idx = i - 1;
        let mut center = 0.0f64;
        for j in i..d.len() {
            center += u[idx][j] * x[j] as f64;
        }
        let remaining = bound as f64 - tail;
        if remaining < -eps {
            return;
        }
        let radius = (remaining.max(0.0) / d[idx]).sqrt() + eps;
        let lo = (-center - radius).ceil() as i128;
        let hi = (-center + radius).floor() as i128;
        for xi in lo..=hi {
            x[idx] = xi;
            let coord = xi as f64 + center;
            self.fp_search(idx, bound, d, u, eps, tail + d[idx] * coord * coord, x, out);
        }
        x[idx] = 0;
    }

    /// The minimum `min { Q(x) : x ∈ L, x ≠ 0 }`, or `None` if the lattice is
    /// empty or not positive definite.
    pub fn minimum(&self) -> Option<i128> {
        if self.dim() == 0 {
            return None;
        }
        let min_diag = (0..self.dim()).map(|i| self.gram[i][i]).min()?;
        let vecs = self.short_vectors(min_diag)?;
        vecs.iter().map(|v| self.norm(v)).min()
    }

    /// All minimal vectors (norm equal to [`minimum`](IntegralForm::minimum)),
    /// both signs included. `None` if not positive definite.
    pub fn minimal_vectors(&self) -> Option<Vec<Vec<i128>>> {
        let m = self.minimum()?;
        let vecs = self.short_vectors(m)?;
        Some(vecs.into_iter().filter(|v| self.norm(v) == m).collect())
    }

    /// The kissing number: the count of minimal vectors. `None` if not positive
    /// definite.
    pub fn kissing_number(&self) -> Option<usize> {
        self.minimal_vectors().map(|v| v.len())
    }

    /// The order of the automorphism group `Aut(L) = { U ∈ GLₙ(ℤ) : UᵀGU = G }`,
    /// or `None` if the lattice is not positive definite or the unrecognized
    /// fallback search exceeds [`AUTO_NODE_BUDGET`]. Recognized diagonal and
    /// `A`/`D`/`E` root-lattice bases use closed-form Weyl/root-system orders. See
    /// [`automorphism_group_order_bounded`](IntegralForm::automorphism_group_order_bounded).
    pub fn automorphism_group_order(&self) -> Option<u128> {
        self.automorphism_group_order_bounded(AUTO_NODE_BUDGET)
    }

    /// [`automorphism_group_order`](IntegralForm::automorphism_group_order) with
    /// an explicit node budget. An automorphism is determined by where it sends a
    /// basis `e₀, …, e_{n-1}`; the images must be lattice vectors `vᵢ` with
    /// `⟨vᵢ, vⱼ⟩ = ⟨eᵢ, eⱼ⟩`. The search enumerates candidate images (lattice
    /// vectors of each basis norm) and backtracks on the Gram constraints; each
    /// complete assignment is an automorphism, so the count is exact. Returns
    /// `None` if more than `node_budget` candidate-nodes are visited.
    pub fn automorphism_group_order_bounded(&self, node_budget: u64) -> Option<u128> {
        let n = self.dim();
        if n == 0 {
            return Some(1);
        }
        if !self.is_positive_definite() {
            return None;
        }
        if let Some(order) = self.automorphism_group_order_fast() {
            return Some(order);
        }
        let max_diag = (0..n).map(|i| self.gram[i][i]).max().unwrap();
        let cands = self.short_vectors(max_diag)?;
        // Precompute G·v for each candidate so inner products are plain dot
        // products: ⟨v_a, v_b⟩ = v_aᵀ G v_b = v_a · (G v_b).
        let gv: Vec<Vec<i128>> = cands.iter().map(|v| self.matvec(v)).collect();
        let norms: Vec<i128> = cands.iter().zip(&gv).map(|(v, gvb)| dot(v, gvb)).collect();
        let diag: Vec<i128> = (0..n).map(|i| self.gram[i][i]).collect();
        let per_level: Vec<Vec<usize>> = (0..n)
            .map(|lvl| {
                (0..cands.len())
                    .filter(|&c| norms[c] == diag[lvl])
                    .collect()
            })
            .collect();
        let mut count: u128 = 0;
        let mut nodes: u64 = 0;
        let mut chosen: Vec<usize> = Vec::with_capacity(n);
        let ok = self.aut_backtrack(
            0,
            &per_level,
            &cands,
            &gv,
            &mut chosen,
            &mut count,
            &mut nodes,
            node_budget,
        );
        if ok {
            Some(count)
        } else {
            None
        }
    }

    /// Closed-form automorphism orders. Literal standard bases are checked first
    /// because they are cheap; then norm-2 roots are used for a basis-independent
    /// simply-laced root-system classifier before falling back to exact search.
    fn automorphism_group_order_fast(&self) -> Option<u128> {
        let n = self.dim();
        if n == 0 {
            return Some(1);
        }

        // d·I_n has signed permutations, including the root lattice A_1^n at d=2.
        let d = self.gram[0][0];
        if d > 0
            && (0..n).all(|i| {
                (0..n).all(|j| {
                    let expected = if i == j { d } else { 0 };
                    self.gram[i][j] == expected
                })
            })
        {
            return signed_permutation_order(n);
        }

        if self.matches_a_cartan() {
            return if n == 1 {
                Some(2)
            } else {
                checked_factorial(n + 1)?.checked_mul(2)
            };
        }
        if self.matches_d_cartan() {
            return match n {
                0 | 1 => None,
                2 => signed_permutation_order(2),
                3 => Some(48),
                4 => Some(1152),
                _ => checked_pow2(n)?.checked_mul(checked_factorial(n)?),
            };
        }
        if self.matches_e6_cartan() {
            return Some(103_680);
        }
        if self.matches_e7_cartan() {
            return Some(2_903_040);
        }
        if self.matches_e8_cartan() {
            return Some(696_729_600);
        }
        self.root_system_automorphism_order()
    }

    fn root_system_automorphism_order(&self) -> Option<u128> {
        if !self.is_even() || self.minimum()? != 2 {
            return None;
        }
        let n = self.dim();
        let mut roots: Vec<Vec<i128>> = Vec::new();
        for root in self.minimal_vectors()? {
            let root = canonical_root(root);
            if !roots.contains(&root) {
                roots.push(root);
            }
        }
        if !rows_generate_full_lattice(&roots, n) {
            return None;
        }

        let mut seen = vec![false; roots.len()];
        let mut kinds = Vec::new();
        for start in 0..roots.len() {
            if seen[start] {
                continue;
            }
            let mut queue = VecDeque::from([start]);
            seen[start] = true;
            let mut component = Vec::new();
            while let Some(i) = queue.pop_front() {
                component.push(i);
                for j in 0..roots.len() {
                    if !seen[j] && self.inner(&roots[i], &roots[j]) != 0 {
                        seen[j] = true;
                        queue.push_back(j);
                    }
                }
            }
            let component_roots: Vec<Vec<i128>> =
                component.into_iter().map(|i| roots[i].clone()).collect();
            let rank =
                crate::linalg::integer::normalize_relation_rows(component_roots.clone()).len();
            let root_count = component_roots.len().checked_mul(2)?;
            kinds.push(RootComponentKind::from_rank_and_roots(rank, root_count)?);
        }

        let mut order = 1u128;
        let mut multiplicities: BTreeMap<RootComponentKind, usize> = BTreeMap::new();
        for kind in kinds {
            order = order.checked_mul(kind.automorphism_order()?)?;
            *multiplicities.entry(kind).or_insert(0) += 1;
        }
        for mult in multiplicities.values() {
            order = order.checked_mul(checked_factorial(*mult)?)?;
        }
        Some(order)
    }

    fn matches_a_cartan(&self) -> bool {
        let n = self.dim();
        if n == 0 {
            return false;
        }
        let edges: Vec<(usize, usize)> = (0..n.saturating_sub(1)).map(|i| (i, i + 1)).collect();
        simple_laced_cartan_matches(&self.gram, &edges)
    }

    fn matches_d_cartan(&self) -> bool {
        let n = self.dim();
        if n < 2 {
            return false;
        }
        if n == 2 {
            return simple_laced_cartan_matches(&self.gram, &[]);
        }
        if n == 3 {
            return simple_laced_cartan_matches(&self.gram, &[(0, 1), (0, 2)]);
        }
        let mut edges: Vec<(usize, usize)> = (0..n - 3).map(|i| (i, i + 1)).collect();
        edges.push((n - 3, n - 2));
        edges.push((n - 3, n - 1));
        simple_laced_cartan_matches(&self.gram, &edges)
    }

    fn matches_e6_cartan(&self) -> bool {
        simple_laced_cartan_matches(&self.gram, &[(0, 1), (1, 2), (2, 3), (3, 4), (2, 5)])
    }

    fn matches_e7_cartan(&self) -> bool {
        simple_laced_cartan_matches(
            &self.gram,
            &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (2, 6)],
        )
    }

    fn matches_e8_cartan(&self) -> bool {
        simple_laced_cartan_matches(
            &self.gram,
            &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (4, 7)],
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn aut_backtrack(
        &self,
        level: usize,
        per_level: &[Vec<usize>],
        cands: &[Vec<i128>],
        gv: &[Vec<i128>],
        chosen: &mut Vec<usize>,
        count: &mut u128,
        nodes: &mut u64,
        budget: u64,
    ) -> bool {
        if level == self.dim() {
            *count += 1;
            return true;
        }
        for &c in &per_level[level] {
            *nodes += 1;
            if *nodes > budget {
                return false;
            }
            let mut ok = true;
            for (b, &cb) in chosen.iter().enumerate() {
                // ⟨v_level, v_b⟩ must equal G[level][b].
                if dot(&cands[c], &gv[cb]) != self.gram[level][b] {
                    ok = false;
                    break;
                }
            }
            if ok {
                chosen.push(c);
                if !self.aut_backtrack(
                    level + 1,
                    per_level,
                    cands,
                    gv,
                    chosen,
                    count,
                    nodes,
                    budget,
                ) {
                    return false;
                }
                chosen.pop();
            }
        }
        true
    }

    /// `G·x` as an integer vector.
    fn matvec(&self, x: &[i128]) -> Vec<i128> {
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

fn dot(a: &[i128], b: &[i128]) -> i128 {
    let mut acc = 0i128;
    for (&x, &y) in a.iter().zip(b) {
        acc = acc
            .checked_add(x.checked_mul(y).expect("lattice dot exceeds i128"))
            .expect("lattice dot exceeds i128");
    }
    acc
}

#[cfg(test)]
mod tests {
    use super::*;

    fn a_n(n: usize) -> IntegralForm {
        // A_n Cartan matrix: 2 on the diagonal, -1 on the off-diagonals.
        let mut g = vec![vec![0i128; n]; n];
        for i in 0..n {
            g[i][i] = 2;
            if i + 1 < n {
                g[i][i + 1] = -1;
                g[i + 1][i] = -1;
            }
        }
        IntegralForm::new(g).unwrap()
    }

    fn d4() -> IntegralForm {
        IntegralForm::new(vec![
            vec![2, -1, 0, 0],
            vec![-1, 2, -1, -1],
            vec![0, -1, 2, 0],
            vec![0, -1, 0, 2],
        ])
        .unwrap()
    }

    fn e8() -> IntegralForm {
        // E_8 Cartan matrix (Bourbaki labelling): even unimodular, det 1.
        IntegralForm::new(vec![
            vec![2, -1, 0, 0, 0, 0, 0, 0],
            vec![-1, 2, -1, 0, 0, 0, 0, 0],
            vec![0, -1, 2, -1, 0, 0, 0, 0],
            vec![0, 0, -1, 2, -1, 0, 0, 0],
            vec![0, 0, 0, -1, 2, -1, 0, -1],
            vec![0, 0, 0, 0, -1, 2, -1, 0],
            vec![0, 0, 0, 0, 0, -1, 2, 0],
            vec![0, 0, 0, 0, -1, 0, 0, 2],
        ])
        .unwrap()
    }

    fn permute_basis(l: &IntegralForm, perm: &[usize]) -> IntegralForm {
        let n = l.dim();
        assert_eq!(perm.len(), n);
        let mut g = vec![vec![0i128; n]; n];
        for i in 0..n {
            for j in 0..n {
                g[i][j] = l.gram()[perm[i]][perm[j]];
            }
        }
        IntegralForm::new(g).unwrap()
    }

    #[test]
    fn rejects_non_symmetric() {
        assert!(IntegralForm::new(vec![vec![1, 2], vec![3, 4]]).is_none());
        assert!(IntegralForm::new(vec![vec![1, 2, 3], vec![2, 4]]).is_none());
        assert!(IntegralForm::new(vec![vec![2, -1], vec![-1, 2]]).is_some());
    }

    #[test]
    fn determinants_and_evenness() {
        assert_eq!(a_n(2).determinant(), 3);
        assert_eq!(a_n(3).determinant(), 4);
        assert_eq!(d4().determinant(), 4);
        assert_eq!(e8().determinant(), 1);
        assert!(e8().is_unimodular());
        assert!(e8().is_even());
        assert!(a_n(2).is_even());
        // Z^3 is odd unimodular.
        let z3 = IntegralForm::diagonal(&[1, 1, 1]);
        assert_eq!(z3.determinant(), 1);
        assert!(z3.is_unimodular());
        assert!(!z3.is_even());
    }

    #[test]
    fn invariant_factors_track_discriminant_group() {
        assert_eq!(a_n(2).invariant_factors(), vec![1, 3]); // ℤ/3
        assert_eq!(d4().invariant_factors(), vec![1, 1, 2, 2]); // (ℤ/2)²
        assert_eq!(e8().invariant_factors(), vec![1, 1, 1, 1, 1, 1, 1, 1]);
        // product of nonzero factors = |det|
        let prod: i128 = d4().invariant_factors().iter().product();
        assert_eq!(prod, d4().determinant().abs());
    }

    #[test]
    fn levels_match_known_values() {
        assert_eq!(IntegralForm::diagonal(&[2]).level(), Some(4)); // A_1 = ⟨2⟩
        assert_eq!(a_n(2).level(), Some(3)); // hexagonal lattice, level 3
        assert_eq!(e8().level(), Some(1)); // even unimodular
                                           // ℤ = ⟨1⟩ is odd: G⁻¹ = [1] has odd diagonal, so the smallest N making
                                           // N·G⁻¹ even-integral is 2 (cf. A_1 = ⟨2⟩ → 4).
        assert_eq!(IntegralForm::diagonal(&[1]).level(), Some(2));
    }

    #[test]
    fn minimum_and_kissing_numbers() {
        // Root lattices: minimum 2, kissing = number of roots.
        assert_eq!(a_n(2).minimum(), Some(2));
        assert_eq!(a_n(2).kissing_number(), Some(6)); // n(n+1) = 6
        assert_eq!(a_n(3).kissing_number(), Some(12)); // 3·4
        assert_eq!(d4().minimum(), Some(2));
        assert_eq!(d4().kissing_number(), Some(24)); // 2n(n-1) = 24
        assert_eq!(e8().minimum(), Some(2));
        assert_eq!(e8().kissing_number(), Some(240));
        // ℤ²: minimum 1, the four ±eᵢ.
        let z2 = IntegralForm::diagonal(&[1, 1]);
        assert_eq!(z2.minimum(), Some(1));
        assert_eq!(z2.kissing_number(), Some(4));
    }

    #[test]
    fn short_vectors_return_original_coordinates_after_basis_reduction() {
        // Uᵀ I U for U = [[1, 10], [0, 1]] is a badly skewed basis of Z².
        // The norm-1 vectors in this basis are ±(1,0) and ±(-10,1).
        let g = IntegralForm::new(vec![vec![1, 10], vec![10, 101]]).unwrap();
        let mut vecs = g.short_vectors(1).unwrap();
        vecs.sort();
        assert_eq!(
            vecs,
            vec![vec![-10, 1], vec![-1, 0], vec![1, 0], vec![10, -1]]
        );
        assert!(vecs.iter().all(|v| g.norm(v) == 1));
    }

    #[test]
    fn short_vectors_are_indefinite_safe() {
        // An indefinite form has no finite short-vector set.
        let hyp = IntegralForm::new(vec![vec![0, 1], vec![1, 0]]).unwrap();
        assert!(!hyp.is_positive_definite());
        assert_eq!(hyp.short_vectors(4), None);
        assert_eq!(hyp.minimum(), None);
        assert_eq!(hyp.automorphism_group_order(), None);
    }

    #[test]
    fn automorphism_orders_match_known() {
        // Aut(Z^n) = signed permutations = 2^n · n!.
        assert_eq!(
            IntegralForm::diagonal(&[1, 1]).automorphism_group_order(),
            Some(8)
        );
        assert_eq!(
            IntegralForm::diagonal(&[1, 1, 1]).automorphism_group_order(),
            Some(48)
        );
        // Aut(A_2) = dihedral of order 12 (W(A_2)=S_3 times ±1).
        assert_eq!(a_n(2).automorphism_group_order(), Some(12));
        // Aut(A_3) = W(A_3) × {±1} = 24 · 2 = 48.
        assert_eq!(a_n(3).automorphism_group_order(), Some(48));
        // |Aut(D_4)| = 1152.
        assert_eq!(d4().automorphism_group_order(), Some(1152));
        // E_8 is recognized by its standard Cartan basis instead of brute-forced.
        assert_eq!(e8().automorphism_group_order_bounded(1), Some(696_729_600));
    }

    #[test]
    fn automorphism_budget_cutoff_reports_none() {
        // Permuted root bases are now recognized by the root-system fast path,
        // independent of the standard Cartan syntax.
        let d4_permuted = permute_basis(&d4(), &[2, 0, 1, 3]);
        assert_eq!(d4_permuted.automorphism_group_order_bounded(1), Some(1152));

        // A tiny budget still forces the fallback search to give up rather than
        // silently truncating on a non-root lattice: an honest None, not a wrong count.
        let generic = IntegralForm::new(vec![vec![2, 1], vec![1, 3]]).unwrap();
        assert_eq!(generic.automorphism_group_order_bounded(0), None);
    }

    #[test]
    fn direct_sum_is_block_diagonal() {
        let sum = a_n(2).direct_sum(&IntegralForm::diagonal(&[1]));
        assert_eq!(sum.dim(), 3);
        assert_eq!(sum.determinant(), 3); // det(A_2) · det(⟨1⟩)
                                          // E_8 ⟂ E_8 is rank-16 even unimodular.
        let e8e8 = e8().direct_sum(&e8());
        assert_eq!(e8e8.dim(), 16);
        assert_eq!(e8e8.determinant(), 1);
        assert!(e8e8.is_even());
        assert_eq!(e8e8.kissing_number(), Some(480)); // 240 + 240
    }
}
