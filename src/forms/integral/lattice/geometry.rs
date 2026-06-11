//! Short-vector enumeration (Fincke–Pohst) and automorphism counting for
//! positive-definite integral lattices. These are the expensive geometric
//! routines that build on top of the basic [`IntegralForm`](super::core::IntegralForm)
//! arithmetic in `core.rs`.

use super::core::{dot, IntegralForm};
use crate::linalg::field::inverse_matrix;
use crate::scalar::{Rational, Scalar};
use std::collections::{BTreeMap, VecDeque};

/// The default node budget for [`IntegralForm::automorphism_group_order`]. Beyond
/// this many backtracking nodes the search reports `None` (the lattice is too
/// large for brute-force automorphism enumeration — e.g. `E₈`, whose Weyl group
/// has order ~7·10⁸, or the Leech lattice). The bound is explicit, not silent.
pub const AUTO_NODE_BUDGET: u128 = 100_000_000;
pub(super) const SHORT_VECTOR_EXACT_ENUM_LIMIT: u128 = 2_000_000;

// ── small combinatorial helpers ──

pub(super) fn checked_factorial(n: usize) -> Option<u128> {
    let mut acc = 1u128;
    for k in 2..=n {
        acc = acc.checked_mul(k as u128)?;
    }
    Some(acc)
}

pub(super) fn checked_pow2(n: usize) -> Option<u128> {
    if n >= 128 {
        None
    } else {
        Some(1u128 << n)
    }
}

pub(super) fn signed_permutation_order(n: usize) -> Option<u128> {
    checked_pow2(n)?.checked_mul(checked_factorial(n)?)
}

pub(super) fn a_root_automorphism_order(n: usize) -> Option<u128> {
    if n == 0 {
        None
    } else if n == 1 {
        Some(2)
    } else {
        checked_factorial(n + 1)?.checked_mul(2)
    }
}

pub(super) fn d_root_automorphism_order(n: usize) -> Option<u128> {
    match n {
        0 | 1 => None,
        2 => signed_permutation_order(2),  // D_2 = A_1 x A_1.
        3 => a_root_automorphism_order(3), // D_3 = A_3.
        4 => checked_pow2(3)?
            .checked_mul(checked_factorial(4)?)?
            .checked_mul(6), // triality: Aut(D_4 diagram) = S_3.
        _ => checked_pow2(n)?.checked_mul(checked_factorial(n)?),
    }
}

fn square_ge_scaled(r: u128, num: u128, den: u128) -> bool {
    match r.checked_mul(r).and_then(|rr| rr.checked_mul(den)) {
        Some(lhs) => lhs >= num,
        None => true,
    }
}

fn ceil_sqrt_rational(x: &Rational) -> Option<i128> {
    if x.sign() != std::cmp::Ordering::Greater {
        return Some(0);
    }
    let num = u128::try_from(x.numer()).ok()?;
    let den = u128::try_from(x.denom()).ok()?;
    let approx = ((num as f64) / (den as f64)).sqrt().ceil();
    let mut hi = if approx.is_finite() && approx >= 0.0 {
        (approx as u128).saturating_add(2).max(1)
    } else {
        1
    };
    while !square_ge_scaled(hi, num, den) {
        hi = hi.checked_mul(2)?;
    }
    let mut lo = 0u128;
    while lo < hi {
        let mid = lo + (hi - lo) / 2;
        if square_ge_scaled(mid, num, den) {
            hi = mid;
        } else {
            lo = mid + 1;
        }
    }
    i128::try_from(lo).ok()
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
            RootComponentKind::A(n) => a_root_automorphism_order(n),
            RootComponentKind::D(n) => d_root_automorphism_order(n),
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

// ── geometry methods on IntegralForm ──

impl IntegralForm {
    /// The LDLᵀ decomposition in floating point: returns `Some((d, u))` where
    /// `d[i]` is the `i`-th pivot and `u[i][j]` (`j > i`) is the upper unit factor,
    /// giving `Q(x) ≈ Σᵢ d[i]·(x[i] + Σ_{j>i} u[i][j]·x[j])²`. Returns `None` if
    /// any pivot is non-positive, which signals unexpected loss of definiteness due to
    /// floating-point error (the caller already guards against indefinite lattices, so
    /// this is a safety fallback rather than an expected code path). Every candidate
    /// vector produced by the float bound is rechecked with the exact integer norm;
    /// false positives are removed, but when a pivot rounds to zero or below the
    /// corresponding branch may be skipped — use `short_vectors_exact_bounded` for
    /// small lattices where the float bound is not needed.
    pub(super) fn ldl(&self) -> Option<(Vec<f64>, Vec<Vec<f64>>)> {
        let n = self.dim();
        let mut d = vec![0.0f64; n];
        let mut l = vec![vec![0.0f64; n]; n]; // unit lower triangular
        for j in 0..n {
            let mut dj = self.gram[j][j] as f64;
            for k in 0..j {
                dj -= l[j][k] * l[j][k] * d[k];
            }
            if dj <= 0.0 {
                return None;
            }
            d[j] = dj;
            l[j][j] = 1.0;
            for i in j + 1..n {
                let mut s = self.gram[i][j] as f64;
                for k in 0..j {
                    s -= l[i][k] * l[j][k] * d[k];
                }
                l[i][j] = s / dj;
            }
        }
        // u[i][j] = L[j][i] for j > i
        let mut u = vec![vec![0.0f64; n]; n];
        for i in 0..n {
            for j in i + 1..n {
                u[i][j] = l[j][i];
            }
        }
        Some((d, u))
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
        if let Some(vecs) = self.short_vectors_exact_bounded(bound, SHORT_VECTOR_EXACT_ENUM_LIMIT) {
            return Some(vecs);
        }
        let (reduced, transform) = self.size_reduced_basis();
        let vecs = reduced.short_vectors_raw(bound)?;
        Some(
            vecs.into_iter()
                .map(|v| map_coords(&transform, &v))
                .collect(),
        )
    }

    pub(super) fn short_vectors_exact_bounded(
        &self,
        bound: i128,
        limit: u128,
    ) -> Option<Vec<Vec<i128>>> {
        let n = self.dim();
        let mat: Vec<Vec<Rational>> = self
            .gram
            .iter()
            .map(|row| row.iter().map(|&x| Rational::int(x)).collect())
            .collect();
        let inv = inverse_matrix(mat)?;
        let mut ranges = Vec::with_capacity(n);
        let mut count = 1u128;
        for i in 0..n {
            let radius2 = Rational::int(bound).mul(&inv[i][i]);
            let r = ceil_sqrt_rational(&radius2)?;
            let ru = u128::try_from(r).ok()?;
            let width = ru.checked_mul(2)?.checked_add(1)?;
            count = count.checked_mul(width)?;
            if count > limit {
                return None;
            }
            ranges.push(r);
        }
        let mut out = Vec::new();
        let mut x = vec![0i128; n];
        self.enumerate_exact_box(&ranges, 0, bound, &mut x, &mut out);
        Some(out)
    }

    fn enumerate_exact_box(
        &self,
        ranges: &[i128],
        idx: usize,
        bound: i128,
        x: &mut [i128],
        out: &mut Vec<Vec<i128>>,
    ) {
        if idx == ranges.len() {
            let q = self.norm(x);
            if q > 0 && q <= bound {
                out.push(x.to_vec());
            }
            return;
        }
        for xi in -ranges[idx]..=ranges[idx] {
            x[idx] = xi;
            self.enumerate_exact_box(ranges, idx + 1, bound, x, out);
        }
        x[idx] = 0;
    }

    fn short_vectors_raw(&self, bound: i128) -> Option<Vec<Vec<i128>>> {
        if !self.is_positive_definite() {
            return None;
        }
        let n = self.dim();
        if n == 0 || bound <= 0 {
            return Some(Vec::new());
        }
        // `ldl` returns None if any pivot rounds to <= 0 (unexpected loss of
        // definiteness under floating-point). Fall back to None so the caller
        // can return an error rather than silently omitting vectors.
        let (d, u) = self.ldl()?;
        let mut out = Vec::new();
        let mut x = vec![0i128; n];
        // Pad the float radius outward; the exact integer filter at the leaf
        // removes any spurious vectors admitted by the float bound. Small boxes
        // are handled above by exact rational bounds; this path is for larger
        // enumerations where the float bound is the practical cutoff.
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
    pub fn automorphism_group_order_bounded(&self, node_budget: u128) -> Option<u128> {
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
        let mut nodes: u128 = 0;
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
            return a_root_automorphism_order(n);
        }
        if self.matches_d_cartan() {
            return d_root_automorphism_order(n);
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
        nodes: &mut u128,
        budget: u128,
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
}
