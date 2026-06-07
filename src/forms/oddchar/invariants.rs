//! Quadratic-form invariants over finite fields of odd characteristic.

use super::{ensure_odd_prime, FiniteOddField};
use crate::clifford::Metric;
use crate::forms::{as_diagonal, WittClassG};
use crate::scalar::Fp;

/// The classification of a nondegenerate-plus-radical diagonal form over `F_P`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct OddCharType {
    /// Characteristic prime.
    pub p: u128,
    /// Field order `q`; equal to `p` for prime fields and `p^n` for extensions.
    pub field_order: u128,
    /// Nondegenerate dimension (number of nonzero diagonal entries).
    pub dim: usize,
    /// Radical (null) dimension.
    pub radical_dim: usize,
    /// Discriminant square-class: `true` if `det` of the nondegenerate part is a
    /// square. With `dim`, a complete isometry invariant over a finite field.
    pub disc_is_square: bool,
    /// The Hasse–Witt invariant — always `+1` over a finite field.
    pub hasse: i8,
}

impl OddCharType {
    pub fn display(&self) -> String {
        let d = if self.disc_is_square { "□" } else { "✶" };
        let field = format!("F_{}", self.field_order);
        let rad = if self.radical_dim > 0 {
            format!(" ⊗ Λ({}^{})", field, self.radical_dim)
        } else {
            String::new()
        };
        format!(
            "{}: dim {} disc {} hasse {:+}{}",
            field, self.dim, d, self.hasse, rad
        )
    }
}

/// The Hasse invariant `∏_{i<j} (q_i, q_j)` over a finite odd field. Finite
/// fields have trivial Brauer group, so every nonzero Hilbert symbol is `+1`;
/// the prime-field [`super::hilbert_symbol`] wrapper still keeps the brute-force
/// witness for tests and pedagogy.
pub fn hasse_invariant_finite_odd<F: FiniteOddField>(metric: &Metric<F>) -> Option<i8> {
    F::ensure_supported()?;
    let metric = as_diagonal(metric)?;
    let qs: Vec<F> = metric.q.iter().copied().filter(|x| !x.is_zero()).collect();
    let mut h = 1i8;
    for i in 0..qs.len() {
        for _j in (i + 1)..qs.len() {
            h *= 1;
        }
    }
    Some(h)
}

/// The discriminant of the nondegenerate diagonal part over any finite odd field.
pub fn discriminant_finite_odd<F: FiniteOddField>(metric: &Metric<F>) -> Option<F> {
    F::ensure_supported()?;
    let metric = as_diagonal(metric)?;
    let mut d = F::one();
    for x in &metric.q {
        if !x.is_zero() {
            d = d.mul(x);
        }
    }
    Some(d)
}

/// Classify a form over any finite field of odd characteristic.
pub fn classify_finite_odd<F: FiniteOddField>(metric: &Metric<F>) -> Option<OddCharType> {
    F::ensure_supported()?;
    let metric = as_diagonal(metric)?;
    let dim = metric.q.iter().filter(|x| !x.is_zero()).count();
    let radical_dim = metric.q.len() - dim;
    let disc = discriminant_finite_odd(&metric)?;
    Some(OddCharType {
        p: F::characteristic_prime(),
        field_order: F::field_order(),
        dim,
        radical_dim,
        disc_is_square: F::is_square_value(disc),
        hasse: hasse_invariant_finite_odd(&metric)?,
    })
}

/// The finite odd-field Witt class `(dim mod 2, signed discriminant class)`.
pub fn finite_odd_witt<F: FiniteOddField>(metric: &Metric<F>) -> Option<WittClassG> {
    F::ensure_supported()?;
    let metric = as_diagonal(metric)?;
    let mut det = F::one();
    let mut m = 0usize;
    for x in &metric.q {
        if !x.is_zero() {
            det = det.mul(x);
            m += 1;
        }
    }
    let signed = if (m * (m.wrapping_sub(1)) / 2) & 1 == 1 {
        det.neg()
    } else {
        det
    };
    let kappa = if F::is_square_value(F::from_i128(-1)) {
        0
    } else {
        1
    };
    Some(WittClassG::OddChar {
        kappa,
        e0: (m & 1) as u8,
        sclass: if F::is_square_value(signed) { 0 } else { 1 },
    })
}

/// The Hasse invariant `∏_{i<j} (q_i, q_j)` of a form (nonzero diagonal entries
/// only; non-diagonal metrics are congruence-diagonalized first). Always `+1`
/// over a finite field.
pub fn hasse_invariant<const P: u128>(metric: &Metric<Fp<P>>) -> Option<i8> {
    ensure_odd_prime::<P>()?;
    hasse_invariant_finite_odd(metric)
}

/// The discriminant (product of the nonzero diagonal entries = det of the
/// nondegenerate part). Non-diagonal metrics are congruence-diagonalized first.
pub fn discriminant<const P: u128>(metric: &Metric<Fp<P>>) -> Option<Fp<P>> {
    ensure_odd_prime::<P>()?;
    discriminant_finite_odd(metric)
}

/// Classify an odd-characteristic form. Non-diagonal metrics are
/// congruence-diagonalized first (char ≠ 2 always has `½`), so any symmetric
/// metric is accepted.
pub fn classify_oddchar<const P: u128>(metric: &Metric<Fp<P>>) -> Option<OddCharType> {
    ensure_odd_prime::<P>()?;
    classify_finite_odd(metric)
}

/// The odd-characteristic Witt class: `(dim mod 2, signed discriminant class)`,
/// with `kappa` = nonsquareness of `−1`. Non-diagonal metrics are
/// congruence-diagonalized first. The signed discriminant `(−1)^{m(m−1)/2}·det`
/// is the genuine Witt invariant; see `witt::WittClassG`.
pub fn oddchar_witt<const P: u128>(metric: &Metric<Fp<P>>) -> Option<WittClassG> {
    ensure_odd_prime::<P>()?;
    finite_odd_witt(metric)
}
