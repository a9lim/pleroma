//! A tiny dependency-free complex number for Gauss sums and Weil matrices.
//!
//! This is a deliberate name-shadow of `num_complex::Complex<f64>` — the crate
//! carries no `num-complex` dependency, and adding one solely for `(f64, f64)`
//! arithmetic over a small finite discriminant group is not worth the transitive
//! cost. The type is `pub` inside the crate and re-exported at the `integral` level.

/// A tiny dependency-free complex number for Gauss sums and Weil matrices.
///
/// Deliberately shadows `num_complex::Complex64`; the discriminant-form Weil
/// representation needs only basic `f64` arithmetic over small finite groups, and
/// this keeps the crate dependency-free.
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
