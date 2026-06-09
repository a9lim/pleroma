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
use std::collections::BTreeSet;

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

fn mod_pow_u128(mut base: u128, mut exp: u128, modulus: u128) -> u128 {
    let mut acc = 1u128;
    base %= modulus;
    while exp > 0 {
        if exp & 1 == 1 {
            acc = (acc * base) % modulus;
        }
        base = (base * base) % modulus;
        exp >>= 1;
    }
    acc
}

fn is_square_mod_odd_p(unit: i128, p: i128) -> bool {
    let u = unit.rem_euclid(p) as u128;
    if u == 0 {
        return true;
    }
    mod_pow_u128(u, ((p as u128) - 1) / 2, p as u128) == 1
}

fn unit_is_antisquare_odd(r: &Rational, p: i128) -> bool {
    let a = unit_part_i128(r.numer(), p);
    let b = unit_part_i128(r.denom(), p);
    let unit = (a.rem_euclid(p) * b.rem_euclid(p)).rem_euclid(p);
    !is_square_mod_odd_p(unit, p)
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

fn relevant_odd_primes(det: i128) -> Vec<u128> {
    let mut n = det.unsigned_abs();
    while n.is_multiple_of(2) {
        n /= 2;
    }
    let mut ps = Vec::new();
    let mut d = 3u128;
    while d <= n / d {
        if n.is_multiple_of(d) {
            ps.push(d);
            while n.is_multiple_of(d) {
                n /= d;
            }
        }
        d += 2;
    }
    if n > 1 {
        ps.push(n);
    }
    ps
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

fn odd_p_excess(diag: &[Rational], p: u128) -> i128 {
    let p_i = p as i128;
    let p_sig = diag
        .iter()
        .map(|d| {
            let v = rat_val(d, p_i).unsigned_abs();
            let pow = pow_mod8(p as i128, v);
            (pow + if unit_is_antisquare_odd(d, p_i) { 4 } else { 0 }).rem_euclid(8)
        })
        .sum::<i128>()
        .rem_euclid(8);
    ((diag.len() as i128) - p_sig).rem_euclid(8)
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

#[allow(dead_code)]
fn rational_p_excess_mod8(lattice: &IntegralForm, diag: &[Rational]) -> i128 {
    relevant_odd_primes(lattice.determinant())
        .into_iter()
        .map(|p| odd_p_excess(diag, p))
        .sum::<i128>()
        .rem_euclid(8)
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
    use crate::forms::{a_n, d_n, e_8};

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
