//! The normalized Gauss sum `GaussSum` and its phase extraction, plus the matrix
//! helpers (`mat_identity`, `mat_mul`, `mat_pow`, `mat_scale`, `mat_approx_eq`) that
//! the Weil-representation builder in the parent module needs.

use super::complex::Complex64;

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

// ── matrix helpers for the Weil representation ──

pub(super) fn mat_identity(n: usize) -> Vec<Vec<Complex64>> {
    let mut out = vec![vec![Complex64::zero(); n]; n];
    for (i, row) in out.iter_mut().enumerate() {
        row[i] = Complex64::one();
    }
    out
}

pub(super) fn mat_mul(a: &[Vec<Complex64>], b: &[Vec<Complex64>]) -> Vec<Vec<Complex64>> {
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

pub(super) fn mat_pow(a: &[Vec<Complex64>], exp: usize) -> Vec<Vec<Complex64>> {
    let mut out = mat_identity(a.len());
    for _ in 0..exp {
        out = mat_mul(a, &out);
    }
    out
}

pub(super) fn mat_scale(a: &[Vec<Complex64>], c: Complex64) -> Vec<Vec<Complex64>> {
    a.iter()
        .map(|row| row.iter().map(|x| x.mul(&c)).collect())
        .collect()
}

pub(super) fn mat_approx_eq(a: &[Vec<Complex64>], b: &[Vec<Complex64>], tol: f64) -> bool {
    a.len() == b.len()
        && a.iter().zip(b).all(|(ra, rb)| {
            ra.len() == rb.len() && ra.iter().zip(rb).all(|(x, y)| x.approx_eq(y, tol))
        })
}
