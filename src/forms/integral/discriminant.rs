//! Discriminant quadratic forms of even integral lattices and Milgram's Gauss sum.
//!
//! For a nonsingular even lattice `L` with Gram matrix `G`, this module uses the
//! standard presentation
//!
//! ```text
//! A_L = L#/L ~= Z^n / G Z^n,    y |-> G^{-1} y
//! q_L(y) = y^T G^{-1} y mod 2Z.
//! ```
//!
//! The normalized Gauss sum of `q_L` has phase `exp(2*pi*i*signature/8)`.

use crate::forms::integral::diagonal::{rational_congruence_diagonal, DegenerateBehavior};
use crate::forms::integral::{Genus, IntegralForm};
use crate::linalg::field::inverse_matrix;
use crate::linalg::integer::{normalize_relation_rows, reduce_integer_vector};
use crate::scalar::{Rational, Scalar};
use std::collections::{BTreeSet, HashSet};

/// A normalized complex Gauss sum, kept dependency-free.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GaussSum {
    pub re: f64,
    pub im: f64,
}

impl GaussSum {
    pub fn abs(&self) -> f64 {
        self.re.hypot(self.im)
    }

    /// Phase as an eighth-root index: `0` for `1`, `1` for `exp(pi*i/4)`, ... .
    /// Returns `None` if the magnitude or angle is not close to an eighth root.
    pub fn phase_mod8(&self, tol: f64) -> Option<i128> {
        if (self.abs() - 1.0).abs() > tol {
            return None;
        }
        let step = std::f64::consts::FRAC_PI_4;
        let raw = (self.im.atan2(self.re) / step).round() as i128;
        let k = raw.rem_euclid(8);
        let target = (k as f64) * step;
        if (self.re - target.cos()).abs() <= tol && (self.im - target.sin()).abs() <= tol {
            Some(k)
        } else {
            None
        }
    }
}

/// A tiny dependency-free complex number for Gauss sums and Weil matrices.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Complex64 {
    pub re: f64,
    pub im: f64,
}

impl Complex64 {
    pub fn zero() -> Self {
        Complex64 { re: 0.0, im: 0.0 }
    }

    pub fn one() -> Self {
        Complex64 { re: 1.0, im: 0.0 }
    }

    pub fn cis(theta: f64) -> Self {
        Complex64 {
            re: theta.cos(),
            im: theta.sin(),
        }
    }

    /// `exp(pi*i*k/4)`.
    pub fn eighth_root(k: i128) -> Self {
        Complex64::cis((k.rem_euclid(8) as f64) * std::f64::consts::FRAC_PI_4)
    }

    pub fn abs(&self) -> f64 {
        self.re.hypot(self.im)
    }

    pub fn add(&self, rhs: &Self) -> Self {
        Complex64 {
            re: self.re + rhs.re,
            im: self.im + rhs.im,
        }
    }

    pub fn sub(&self, rhs: &Self) -> Self {
        Complex64 {
            re: self.re - rhs.re,
            im: self.im - rhs.im,
        }
    }

    pub fn mul(&self, rhs: &Self) -> Self {
        Complex64 {
            re: self.re * rhs.re - self.im * rhs.im,
            im: self.re * rhs.im + self.im * rhs.re,
        }
    }

    pub fn scale(&self, c: f64) -> Self {
        Complex64 {
            re: self.re * c,
            im: self.im * c,
        }
    }

    pub fn approx_eq(&self, rhs: &Self, tol: f64) -> bool {
        self.sub(rhs).abs() <= tol
    }
}

/// The finite discriminant quadratic module of an even lattice.
#[derive(Clone, Debug, PartialEq)]
pub struct DiscriminantForm {
    /// Nontrivial invariant factors of `A_L`.
    pub group: Vec<i128>,
    /// Canonical representatives `y` for `Z^n / GZ^n`.
    pub reps: Vec<Vec<i128>>,
    /// The exact inverse Gram matrix.
    pub gram_inv: Vec<Vec<Rational>>,
}

fn mat_identity(n: usize) -> Vec<Vec<Complex64>> {
    let mut out = vec![vec![Complex64::zero(); n]; n];
    for (i, row) in out.iter_mut().enumerate() {
        row[i] = Complex64::one();
    }
    out
}

fn mat_mul(a: &[Vec<Complex64>], b: &[Vec<Complex64>]) -> Vec<Vec<Complex64>> {
    let n = a.len();
    let m = b.first().map_or(0, Vec::len);
    let inner = b.len();
    let mut out = vec![vec![Complex64::zero(); m]; n];
    for i in 0..n {
        for k in 0..inner {
            for j in 0..m {
                out[i][j] = out[i][j].add(&a[i][k].mul(&b[k][j]));
            }
        }
    }
    out
}

fn mat_pow(a: &[Vec<Complex64>], exp: usize) -> Vec<Vec<Complex64>> {
    let mut out = mat_identity(a.len());
    for _ in 0..exp {
        out = mat_mul(a, &out);
    }
    out
}

fn mat_scale(a: &[Vec<Complex64>], c: Complex64) -> Vec<Vec<Complex64>> {
    a.iter()
        .map(|row| row.iter().map(|x| x.mul(&c)).collect())
        .collect()
}

fn mat_approx_eq(a: &[Vec<Complex64>], b: &[Vec<Complex64>], tol: f64) -> bool {
    a.len() == b.len()
        && a.iter().zip(b).all(|(ra, rb)| {
            ra.len() == rb.len() && ra.iter().zip(rb).all(|(x, y)| x.approx_eq(y, tol))
        })
}

fn rational_mod_int(x: Rational, modulus: i128) -> Rational {
    debug_assert!(modulus > 0);
    let den = x.denom();
    let mden = den
        .checked_mul(modulus)
        .expect("rational modulus exceeds i128");
    Rational::new(x.numer().rem_euclid(mden), den)
}

fn rational_to_f64(x: &Rational) -> f64 {
    (x.numer() as f64) / (x.denom() as f64)
}

fn dot_inv(v: &[i128], inv: &[Vec<Rational>], w: &[i128]) -> Rational {
    let n = v.len();
    let mut acc = Rational::zero();
    for i in 0..n {
        if v[i] == 0 {
            continue;
        }
        for (j, &wj) in w.iter().enumerate() {
            if wj == 0 {
                continue;
            }
            acc = acc.add(&Rational::int(v[i]).mul(&inv[i][j]).mul(&Rational::int(wj)));
        }
    }
    acc
}

fn enumerate_hnf_reps(rows: &[Vec<i128>]) -> Option<Vec<Vec<i128>>> {
    let n = rows.len();
    if n == 0 {
        return Some(vec![Vec::new()]);
    }
    if rows.iter().any(|r| r.len() != n) {
        return None;
    }
    let mut pivots = Vec::with_capacity(n);
    for (i, row) in rows.iter().enumerate() {
        let lead = row.iter().position(|&x| x != 0)?;
        if lead != i || row[i] <= 0 {
            return None;
        }
        pivots.push(row[i]);
    }

    let mut reps = BTreeSet::new();
    let mut cur = vec![0i128; n];
    fn rec(
        idx: usize,
        pivots: &[i128],
        cur: &mut [i128],
        rows: &[Vec<i128>],
        reps: &mut BTreeSet<Vec<i128>>,
    ) {
        if idx == pivots.len() {
            let mut v = cur.to_vec();
            reduce_integer_vector(&mut v, rows.to_vec());
            reps.insert(v);
            return;
        }
        for x in 0..pivots[idx] {
            cur[idx] = x;
            rec(idx + 1, pivots, cur, rows, reps);
        }
        cur[idx] = 0;
    }
    rec(0, &pivots, &mut cur, rows, &mut reps);
    Some(reps.into_iter().collect())
}

/// Largest discriminant group [`DiscriminantForm::is_isomorphic`] will tabulate; past
/// it the Cayley-table build is refused with `None` (an honest budget, like
/// [`crate::forms::AUTO_NODE_BUDGET`]), never a wrong answer.
const ISO_GROUP_CAP: usize = 256;

/// Default node budget for the isomorphism search (candidate generator-images tried).
const ISO_NODE_BUDGET: u128 = 50_000_000;

/// The finite-abelian-group data of a discriminant form needed to compare two of
/// them: the identity index, each element's `q_L` value and additive order, and the
/// full Cayley addition table (indices into `reps`).
struct IsoTables {
    zero: usize,
    q: Vec<Rational>,
    order: Vec<usize>,
    add: Vec<Vec<usize>>,
}

/// The subgroup generated by `gens`, as the set of element indices.
fn subgroup_closure(t: &IsoTables, gens: &[usize]) -> HashSet<usize> {
    let mut set: HashSet<usize> = HashSet::new();
    set.insert(t.zero);
    let mut frontier = vec![t.zero];
    while let Some(x) = frontier.pop() {
        for &g in gens {
            let nx = t.add[x][g];
            if set.insert(nx) {
                frontier.push(nx);
            }
        }
    }
    set
}

/// A minimal generating set, chosen greedily by maximal order (which realizes the
/// invariant-factor count for a finite abelian group).
fn min_generators(t: &IsoTables) -> Vec<usize> {
    let n = t.order.len();
    let mut gens: Vec<usize> = Vec::new();
    let mut covered = subgroup_closure(t, &gens);
    while covered.len() < n {
        let g = (0..n)
            .filter(|i| !covered.contains(i))
            .max_by_key(|&i| t.order[i])
            .expect("a non-covered element exists while |covered| < |A|");
        gens.push(g);
        covered = subgroup_closure(t, &gens);
    }
    gens
}

/// Given images for the generators of `lt`, extend by the homomorphism property
/// (BFS over `lt`'s generator steps) and check the result is a `q`-preserving
/// bijection `lt → mt`. Returns `false` on any inconsistency.
fn verify_iso(lt: &IsoTables, mt: &IsoTables, gens: &[usize], img: &[usize]) -> bool {
    let n = lt.order.len();
    let mut phi = vec![usize::MAX; n];
    phi[lt.zero] = mt.zero;
    let mut frontier = vec![lt.zero];
    while let Some(x) = frontier.pop() {
        for (t, &g) in gens.iter().enumerate() {
            let nx = lt.add[x][g];
            let nimg = mt.add[phi[x]][img[t]];
            if phi[nx] == usize::MAX {
                phi[nx] = nimg;
                frontier.push(nx);
            } else if phi[nx] != nimg {
                return false; // not a well-defined homomorphism
            }
        }
    }
    if phi.contains(&usize::MAX) {
        return false; // gens did not generate (should not happen)
    }
    // Injective ⇒ bijective (equal finite cardinality).
    let mut seen: HashSet<usize> = HashSet::new();
    if !phi.iter().all(|&p| seen.insert(p)) {
        return false;
    }
    // q preserved on *every* element (the complete quadratic-form check; homomorphism
    // + matching q on generators alone does not force it).
    (0..n).all(|i| mt.q[phi[i]] == lt.q[i])
}

/// DFS over generator-image assignments (pruned by equal order and equal `q`),
/// returning `Some(true)` on the first valid isomorphism, `Some(false)` if exhausted,
/// `None` if the node budget runs out first.
fn search_iso(
    lt: &IsoTables,
    mt: &IsoTables,
    gens: &[usize],
    img: &mut Vec<usize>,
    budget: &mut u128,
) -> Option<bool> {
    let depth = img.len();
    if depth == gens.len() {
        return Some(verify_iso(lt, mt, gens, img));
    }
    let g = gens[depth];
    for cand in 0..mt.order.len() {
        if mt.order[cand] != lt.order[g] || mt.q[cand] != lt.q[g] {
            continue;
        }
        if *budget == 0 {
            return None;
        }
        *budget -= 1;
        img.push(cand);
        match search_iso(lt, mt, gens, img, budget) {
            Some(true) => return Some(true),
            Some(false) => {}
            None => return None,
        }
        img.pop();
    }
    Some(false)
}

impl DiscriminantForm {
    /// Build `q_L` for a nonsingular even lattice. Odd lattices return `None`
    /// because `x^T G x mod 2Z` is not well-defined on `L#/L`.
    pub fn from_lattice(lattice: &IntegralForm) -> Option<Self> {
        if !lattice.is_even() || lattice.determinant() == 0 {
            return None;
        }
        let mat: Vec<Vec<Rational>> = lattice
            .gram()
            .iter()
            .map(|row| row.iter().map(|&x| Rational::int(x)).collect())
            .collect();
        let gram_inv = inverse_matrix(mat)?;
        let hnf = normalize_relation_rows(lattice.gram().to_vec());
        let reps = enumerate_hnf_reps(&hnf)?;
        let det = lattice.determinant().unsigned_abs() as usize;
        if reps.len() != det {
            return None;
        }
        let group = lattice
            .invariant_factors()
            .into_iter()
            .filter(|&d| d > 1)
            .collect();
        Some(DiscriminantForm {
            group,
            reps,
            gram_inv,
        })
    }

    /// `q_L(y) = y^T G^{-1} y mod 2Z`, represented in `[0, 2)`.
    pub fn quadratic_value_mod2(&self, y: &[i128]) -> Rational {
        rational_mod_int(dot_inv(y, &self.gram_inv, y), 2)
    }

    /// `b_L(y,z) = y^T G^{-1} z mod Z`, represented in `[0, 1)`.
    pub fn bilinear_value_mod1(&self, y: &[i128], z: &[i128]) -> Rational {
        rational_mod_int(dot_inv(y, &self.gram_inv, z), 1)
    }

    /// The normalized Gauss sum
    /// `|A_L|^{-1/2} * sum_x exp(pi*i*q_L(x))`.
    pub fn gauss_sum(&self) -> GaussSum {
        let mut re = 0.0f64;
        let mut im = 0.0f64;
        for y in &self.reps {
            let theta = std::f64::consts::PI * rational_to_f64(&self.quadratic_value_mod2(y));
            re += theta.cos();
            im += theta.sin();
        }
        let scale = 1.0 / (self.reps.len() as f64).sqrt();
        GaussSum {
            re: re * scale,
            im: im * scale,
        }
    }

    /// The Milgram phase as `signature mod 8`, extracted from the Gauss sum.
    pub fn milgram_signature_mod8(&self) -> Option<i128> {
        self.gauss_sum().phase_mod8(1e-8)
    }

    /// The `reps` index of the coset containing the raw integer vector `v`.
    fn element_index(&self, v: &[i128]) -> Option<usize> {
        self.reps
            .iter()
            .position(|r| self.equivalent_mod_lattice(r, v))
    }

    /// Tabulate the finite abelian group `(A_L, +)` with each element's `q_L` value
    /// and order, plus the full addition table. `None` past `ISO_GROUP_CAP`.
    fn iso_tables(&self) -> Option<IsoTables> {
        let n = self.reps.len();
        if n > ISO_GROUP_CAP {
            return None;
        }
        let zero = self.reps.iter().position(|r| r.iter().all(|&x| x == 0))?;
        let q: Vec<Rational> = self
            .reps
            .iter()
            .map(|r| self.quadratic_value_mod2(r))
            .collect();
        let mut add = vec![vec![0usize; n]; n];
        for i in 0..n {
            for j in 0..n {
                let s: Vec<i128> = self.reps[i]
                    .iter()
                    .zip(&self.reps[j])
                    .map(|(&a, &b)| a + b)
                    .collect();
                add[i][j] = self.element_index(&s)?;
            }
        }
        let mut order = vec![1usize; n];
        for i in 0..n {
            let mut cur = i;
            let mut k = 1usize;
            while cur != zero {
                cur = add[cur][i];
                k += 1;
            }
            order[i] = k;
        }
        Some(IsoTables {
            zero,
            q,
            order,
            add,
        })
    }

    /// Whether two discriminant quadratic forms `(A_L, q_L)` and `(A_M, q_M)` are
    /// **isomorphic** — equal invariant factors plus a `q`-preserving group
    /// isomorphism. This is the finite-quadratic-module half of **Nikulin's
    /// criterion** (Nikulin, *Integral symmetric bilinear forms…*, Izv. Akad. Nauk
    /// SSSR **43** (1979), Cor. 1.9.4): two **even** lattices share a genus iff their
    /// signature pairs agree and their discriminant forms are isomorphic. Both inputs
    /// are even-lattice discriminant forms (the [`from_lattice`](Self::from_lattice)
    /// boundary); the signature half is checked separately by the caller.
    ///
    /// `Some(true)`/`Some(false)` is a decided answer; `None` only past the budget
    /// (group larger than `ISO_GROUP_CAP`, or the search exceeding the node budget)
    /// — an honest unknown, never a wrong value. A cross-check of two shipped routes
    /// (this and `are_in_same_genus`), not a p-adic-symbol reimplementation.
    pub fn is_isomorphic(&self, other: &Self) -> Option<bool> {
        self.is_isomorphic_bounded(other, ISO_NODE_BUDGET)
    }

    /// [`is_isomorphic`](Self::is_isomorphic) with an explicit node budget.
    pub fn is_isomorphic_bounded(&self, other: &Self, node_budget: u128) -> Option<bool> {
        if self.reps.len() != other.reps.len() {
            return Some(false);
        }
        let mut g1 = self.group.clone();
        let mut g2 = other.group.clone();
        g1.sort_unstable();
        g2.sort_unstable();
        if g1 != g2 {
            return Some(false);
        }
        let lt = self.iso_tables()?;
        let mt = other.iso_tables()?;
        // Necessary: the q-value multisets must agree (canonical reps ⇒ exact keys).
        let mut ql: Vec<(i128, i128)> = lt.q.iter().map(|x| (x.numer(), x.denom())).collect();
        let mut qm: Vec<(i128, i128)> = mt.q.iter().map(|x| (x.numer(), x.denom())).collect();
        ql.sort_unstable();
        qm.sort_unstable();
        if ql != qm {
            return Some(false);
        }
        let gens = min_generators(&lt);
        let mut budget = node_budget;
        let mut img: Vec<usize> = Vec::with_capacity(gens.len());
        search_iso(&lt, &mt, &gens, &mut img, &mut budget)
    }

    fn equivalent_mod_lattice(&self, a: &[i128], b: &[i128]) -> bool {
        let n = self.gram_inv.len();
        if a.len() != n || b.len() != n {
            return false;
        }
        let diff: Vec<i128> = a.iter().zip(b).map(|(&x, &y)| x - y).collect();
        for row in &self.gram_inv {
            let mut coord = Rational::zero();
            for (r, &d) in row.iter().zip(&diff) {
                if d != 0 {
                    coord = coord.add(&r.mul(&Rational::int(d)));
                }
            }
            if !coord.is_integer() {
                return false;
            }
        }
        true
    }

    fn negation_matrix(&self) -> Option<Vec<Vec<Complex64>>> {
        let n = self.reps.len();
        let mut out = vec![vec![Complex64::zero(); n]; n];
        for (col, gamma) in self.reps.iter().enumerate() {
            let neg_gamma: Vec<i128> = gamma.iter().map(|&x| -x).collect();
            let row = self
                .reps
                .iter()
                .position(|delta| self.equivalent_mod_lattice(delta, &neg_gamma))?;
            out[row][col] = Complex64::one();
        }
        Some(out)
    }

    fn weil_t_matrix(&self) -> Vec<Vec<Complex64>> {
        let t = self.weil_t();
        let mut out = vec![vec![Complex64::zero(); t.len()]; t.len()];
        for (i, z) in t.into_iter().enumerate() {
            out[i][i] = z;
        }
        out
    }

    /// The diagonal Weil `T` multipliers `exp(pi*i*q_L(gamma))`.
    pub fn weil_t(&self) -> Vec<Complex64> {
        self.reps
            .iter()
            .map(|gamma| {
                let theta =
                    std::f64::consts::PI * rational_to_f64(&self.quadratic_value_mod2(gamma));
                Complex64::cis(theta)
            })
            .collect()
    }

    /// The phase index of the `S` prefactor in the standard Weil convention:
    /// `exp(-2*pi*i*sign/8)`. The existing Milgram Gauss sum stores the conjugate
    /// phase `exp(+2*pi*i*sign/8)`, so this returns `-sign mod 8`.
    pub fn weil_s_prefactor_phase_mod8(&self) -> Option<i128> {
        Some((-self.milgram_signature_mod8()?).rem_euclid(8))
    }

    /// Recover the positive Milgram signature phase from the Weil `S` prefactor.
    pub fn weil_s_recovers_milgram_phase_mod8(&self) -> Option<i128> {
        Some((-self.weil_s_prefactor_phase_mod8()?).rem_euclid(8))
    }

    /// The Weil `S` matrix in the basis of discriminant representatives:
    /// `(sigma/sqrt(|A|)) * exp(-2*pi*i*b_L(gamma,delta))`.
    pub fn weil_s(&self) -> Option<Vec<Vec<Complex64>>> {
        let n = self.reps.len();
        if n == 0 {
            return None;
        }
        let sigma = Complex64::eighth_root(self.weil_s_prefactor_phase_mod8()?);
        let scale = 1.0 / (n as f64).sqrt();
        let mut out = vec![vec![Complex64::zero(); n]; n];
        for (col, gamma) in self.reps.iter().enumerate() {
            for (row, delta) in self.reps.iter().enumerate() {
                let theta = -2.0
                    * std::f64::consts::PI
                    * rational_to_f64(&self.bilinear_value_mod1(gamma, delta));
                out[row][col] = sigma.mul(&Complex64::cis(theta)).scale(scale);
            }
        }
        Some(out)
    }

    /// Verify the finite Weil representation bookkeeping for this discriminant
    /// form. With the standard `S` prefactor, the honest metaplectic relations are
    /// `S^2 = sigma^2 * (gamma -> -gamma)`, `S^4 = sigma^4 * I`, and
    /// `(ST)^3 = S^2`; for unimodular signature `0 mod 8` these collapse to the
    /// familiar scalar relations.
    pub fn verify_weil_relations(&self) -> bool {
        let Some(s_phase) = self.weil_s_prefactor_phase_mod8() else {
            return false;
        };
        if self.weil_s_recovers_milgram_phase_mod8() != self.milgram_signature_mod8() {
            return false;
        }
        let Some(s) = self.weil_s() else {
            return false;
        };
        let t = self.weil_t_matrix();
        let Some(neg) = self.negation_matrix() else {
            return false;
        };
        let tol = 1e-8;
        if self.weil_t().iter().any(|z| (z.abs() - 1.0).abs() > tol) {
            return false;
        }
        let s2 = mat_pow(&s, 2);
        let s4 = mat_pow(&s, 4);
        let st3 = mat_pow(&mat_mul(&s, &t), 3);
        let s2_target = mat_scale(&neg, Complex64::eighth_root(2 * s_phase));
        let s4_target = mat_scale(
            &mat_identity(self.reps.len()),
            Complex64::eighth_root(4 * s_phase),
        );
        mat_approx_eq(&s2, &s2_target, tol)
            && mat_approx_eq(&s4, &s4_target, tol)
            && mat_approx_eq(&st3, &s2, tol)
    }
}

fn pow_mod8(mut base: i128, mut exp: u128) -> i128 {
    base = base.rem_euclid(8);
    let mut acc = 1i128;
    while exp > 0 {
        if exp & 1 == 1 {
            acc = (acc * base).rem_euclid(8);
        }
        base = (base * base).rem_euclid(8);
        exp >>= 1;
    }
    acc
}

fn v_p_i128(mut x: i128, p: i128) -> i128 {
    debug_assert!(x != 0);
    let mut k = 0i128;
    while x % p == 0 {
        x /= p;
        k += 1;
    }
    k
}

fn unit_part_i128(mut x: i128, p: i128) -> i128 {
    while x % p == 0 {
        x /= p;
    }
    x
}

fn rat_val(r: &Rational, p: i128) -> i128 {
    v_p_i128(r.numer(), p) - v_p_i128(r.denom(), p)
}

fn unit_mod8(r: &Rational) -> i128 {
    let a = unit_part_i128(r.numer(), 2).rem_euclid(8);
    let b = unit_part_i128(r.denom(), 2).rem_euclid(8);
    (a * b).rem_euclid(8)
}

fn is_antisquare_2(u: i128) -> bool {
    matches!(u.rem_euclid(8), 3 | 5)
}

fn diagonal_entries(lattice: &IntegralForm) -> Option<Vec<Rational>> {
    if lattice.determinant() == 0 {
        return None;
    }
    Some(rational_congruence_diagonal(
        lattice.gram(),
        DegenerateBehavior::RequireNonsingular,
    ))
}

fn two_adic_oddity(diag: &[Rational]) -> i128 {
    diag.iter()
        .map(|d| {
            let u = unit_mod8(d);
            let antisquare = rat_val(d, 2).rem_euclid(2) != 0 && is_antisquare_2(u);
            (u + if antisquare { 4 } else { 0 }).rem_euclid(8)
        })
        .sum::<i128>()
        .rem_euclid(8)
}

fn symbol_p_excess_mod8(p: u128, scale: u128, dim: usize, sign: i128) -> i128 {
    let q = pow_mod8(p as i128, scale);
    let antisquare = if scale % 2 == 1 && sign < 0 { 4 } else { 0 };
    ((dim as i128) * (q - 1) + antisquare).rem_euclid(8)
}

/// Signature mod 8 from the Conway-Sloane oddity formula, using exact rational
/// diagonalization as an independent check on Milgram's Gauss sum.
pub fn genus_signature_mod8(lattice: &IntegralForm) -> Option<i128> {
    let genus = Genus::of(lattice)?;
    let diag = diagonal_entries(lattice)?;
    let oddity = two_adic_oddity(&diag);
    let p_excess = genus
        .primes()
        .into_iter()
        .filter(|&p| p != 2)
        .flat_map(|p| {
            genus
                .symbol_at(p)
                .iter()
                .map(move |s| symbol_p_excess_mod8(p, s.scale, s.dim, s.sign))
        })
        .sum::<i128>()
        .rem_euclid(8);
    Some((oddity - p_excess).rem_euclid(8))
}

/// Verify Milgram/van der Blij for an even lattice, comparing the discriminant
/// Gauss-sum phase against both exact real signature and the genus oddity route.
pub fn verify_milgram(lattice: &IntegralForm) -> Option<bool> {
    let disc = DiscriminantForm::from_lattice(lattice)?;
    let phase = disc.milgram_signature_mod8()?;
    let (pos, neg) = lattice.signature();
    let sig = (pos as i128 - neg as i128).rem_euclid(8);
    let genus_sig = genus_signature_mod8(lattice)?;
    Some(phase == sig && genus_sig == sig)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forms::{a_n, are_in_same_genus, d16_plus, d_n, e_6, e_7, e_8, IntegralForm};

    /// Nikulin's right-hand side: equal signature pairs and isomorphic discriminant
    /// quadratic forms. Both lattices must be even (the `from_lattice` boundary).
    fn nikulin_rhs(a: &IntegralForm, b: &IntegralForm) -> bool {
        if a.signature() != b.signature() {
            return false;
        }
        let qa = DiscriminantForm::from_lattice(a).expect("even lattice a");
        let qb = DiscriminantForm::from_lattice(b).expect("even lattice b");
        qa.is_isomorphic(&qb) == Some(true)
    }

    #[test]
    fn discriminant_iso_is_reflexive_and_q_sensitive() {
        for l in [a_n(1), a_n(3), d_n(4), e_6(), e_7(), e_8()] {
            let q = DiscriminantForm::from_lattice(&l).unwrap();
            assert_eq!(q.is_isomorphic(&q), Some(true), "reflexive");
        }
        // A_1 and E_7 share the group ℤ/2 but have q-values 1/2 vs 3/2 — *not*
        // isomorphic forms. The search must see q, not just the group.
        let a1 = DiscriminantForm::from_lattice(&a_n(1)).unwrap();
        let e7 = DiscriminantForm::from_lattice(&e_7()).unwrap();
        assert_eq!(a1.group, e7.group, "same invariant factors ℤ/2");
        assert_eq!(a1.is_isomorphic(&e7), Some(false), "q distinguishes them");
        // Different groups: ℤ/3 (A_2) vs (ℤ/2)² (A_1 ⊕ A_1).
        let a2 = DiscriminantForm::from_lattice(&a_n(2)).unwrap();
        let a1a1 = DiscriminantForm::from_lattice(&a_n(1).direct_sum(&a_n(1))).unwrap();
        assert_eq!(a2.is_isomorphic(&a1a1), Some(false));
    }

    #[test]
    fn nikulin_genus_iff_signature_and_discriminant_form() {
        // The Milnor pair: even unimodular rank 16, same genus, non-isometric, both
        // with trivial discriminant form — Nikulin says same genus, and it is.
        let e8e8 = e_8().direct_sum(&e_8());
        let d16 = d16_plus();
        assert!(nikulin_rhs(&e8e8, &d16));
        assert!(are_in_same_genus(&e8e8, &d16));

        // are_in_same_genus ⟺ (equal signatures ∧ isomorphic discriminant forms)
        // across the even-lattice zoo.
        let zoo = [
            a_n(1),
            a_n(2),
            a_n(3),
            a_n(1).direct_sum(&a_n(1)),
            d_n(4),
            e_6(),
            e_7(),
            e_8(),
        ];
        for (i, a) in zoo.iter().enumerate() {
            for b in &zoo[i..] {
                assert_eq!(
                    are_in_same_genus(a, b),
                    nikulin_rhs(a, b),
                    "Nikulin equivalence failed for a pair"
                );
            }
        }
    }

    #[test]
    fn a1_discriminant_form_has_quarter_turn_phase() {
        let a1 = a_n(1);
        let disc = DiscriminantForm::from_lattice(&a1).unwrap();
        assert_eq!(disc.group, vec![2]);
        assert_eq!(disc.reps.len(), 2);
        assert_eq!(disc.quadratic_value_mod2(&[1]), Rational::new(1, 2));
        assert_eq!(disc.milgram_signature_mod8(), Some(1));
        assert_eq!(disc.weil_s_prefactor_phase_mod8(), Some(7));
        assert_eq!(disc.weil_s_recovers_milgram_phase_mod8(), Some(1));
        assert!(disc.verify_weil_relations());
        assert_eq!(verify_milgram(&a1), Some(true));
    }

    #[test]
    fn ade_root_lattices_match_milgram_phase() {
        for n in 1..=5 {
            let a = a_n(n);
            let disc = DiscriminantForm::from_lattice(&a).unwrap();
            assert_eq!(disc.group, vec![n as i128 + 1]);
            assert_eq!(disc.milgram_signature_mod8(), Some(n as i128 % 8));
            assert!(disc.verify_weil_relations(), "Weil relations A_{n}");
            assert_eq!(verify_milgram(&a), Some(true), "A_{n}");
        }

        let d4 = d_n(4);
        let disc = DiscriminantForm::from_lattice(&d4).unwrap();
        assert_eq!(disc.group, vec![2, 2]);
        assert_eq!(disc.milgram_signature_mod8(), Some(4));
        let gs = disc.gauss_sum();
        assert!((gs.re + 1.0).abs() < 1e-8 && gs.im.abs() < 1e-8);
        assert_eq!(disc.weil_s_recovers_milgram_phase_mod8(), Some(4));
        assert!(disc.verify_weil_relations());
        assert_eq!(verify_milgram(&d4), Some(true));
    }

    #[test]
    fn e8_is_unimodular_and_milgram_trivial() {
        let e8 = e_8();
        let disc = DiscriminantForm::from_lattice(&e8).unwrap();
        assert!(disc.group.is_empty());
        assert_eq!(disc.reps, vec![vec![0; 8]]);
        assert_eq!(disc.milgram_signature_mod8(), Some(0));
        assert_eq!(disc.weil_t(), vec![Complex64::one()]);
        assert_eq!(disc.weil_s().unwrap(), vec![vec![Complex64::one()]]);
        assert!(disc.verify_weil_relations());
        assert_eq!(verify_milgram(&e8), Some(true));

        let e8e8 = e8.direct_sum(&e8);
        assert_eq!(
            DiscriminantForm::from_lattice(&e8e8)
                .unwrap()
                .milgram_signature_mod8(),
            Some(0)
        );
        assert_eq!(verify_milgram(&e8e8), Some(true));
    }

    #[test]
    fn odd_lattices_have_no_even_discriminant_quadratic_form() {
        assert!(DiscriminantForm::from_lattice(&IntegralForm::diagonal(&[1])).is_none());
    }
}
