//! The characteristic-0 Clifford classifier — the symmetry-completing companion
//! to `arf.rs`. Where the Arf invariant returns the isomorphism class of a
//! char-2 (nimber) Clifford algebra, this returns the isomorphism class of a
//! char-0 one as a concrete matrix algebra over ℝ, ℂ, or ℍ.
//!
//! ## The two tables
//!
//! Over a **real-closed** field (the surreals, or ℚ as a stand-in) every
//! nonzero square can be rescaled to ±1 (the field is real-closed, so √ exists
//! for positive elements — even √ω = ω^{1/2}). So a *diagonal* metric is
//! classified by its signature `(p, q, r)` = (#positive, #negative, #null)
//! squares, and the nondegenerate `Cl(p,q)` follows the 8-fold Bott table
//! indexed by `s = (q − p) mod 8` (with `n = p+q`):
//!
//! | s | algebra            |   | s | algebra            |
//! |---|--------------------|---|---|--------------------|
//! | 0 | ℝ(2^{n/2})         |   | 4 | ℍ(2^{(n−2)/2})     |
//! | 1 | ℂ(2^{(n−1)/2})     |   | 5 | ℂ(2^{(n−1)/2})     |
//! | 2 | ℍ(2^{(n−2)/2})     |   | 6 | ℝ(2^{n/2})         |
//! | 3 | ℍ(2^{(n−3)/2})²    |   | 7 | ℝ(2^{(n−1)/2})²    |
//!
//! Over an **algebraically closed** field (surcomplex) all nonzero squares are
//! equivalent, so only `(n, r)` matter and the classification is 2-fold:
//! `Cl(n,ℂ) ≅ ℂ(2^{n/2})` for n even, `ℂ(2^{(n−1)/2})²` for n odd.
//!
//! The null directions (radical of dim `r`) contribute an exterior factor:
//! `Cl(p,q,r) ≅ Cl(p,q) ⊗ Λ(F^r)` over the ground field `F ∈ {ℝ, ℂ}`.

use crate::clifford::Metric;
use crate::scalar::{Rational, Scalar};
use crate::scalar::Surcomplex;
use crate::scalar::Surreal;
use std::cmp::Ordering;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BaseField {
    R,
    C,
    H,
}

impl BaseField {
    fn symbol(self) -> &'static str {
        match self {
            BaseField::R => "R",
            BaseField::C => "C",
            BaseField::H => "H",
        }
    }
}

/// The isomorphism class of a char-0 Clifford algebra: a matrix algebra (or a
/// direct sum of two of them) over ℝ/ℂ/ℍ, optionally tensored with the exterior
/// algebra of the metric's radical.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CliffordType {
    /// The division ring underlying the matrix algebra.
    pub base: BaseField,
    /// `m` such that the (semisimple) core is `M_m(base)` (or two copies of it).
    pub matrix_dim: usize,
    /// Whether the core is a direct sum of two equal matrix algebras (`⊕`).
    pub doubled: bool,
    /// Dimension of the metric radical (null directions): an `Λ(ground^r)` factor.
    pub radical_dim: usize,
    /// The ground field of the classification (ℝ for real, ℂ for surcomplex);
    /// the field over which the radical's exterior factor is taken.
    pub ground: BaseField,
    /// The nondegenerate signature `(p, q)` (positive, negative squares). For the
    /// complex case `q` is 0 and `p` is the nondegenerate dimension.
    pub signature: (usize, usize),
}

impl CliffordType {
    /// Human-readable name, e.g. `M_2(H)`, `M_4(R) ⊕ M_4(R)`, `C ⊗ Λ(R^1)`.
    pub fn display(&self) -> String {
        let unit = if self.matrix_dim == 1 {
            self.base.symbol().to_string()
        } else {
            format!("M_{}({})", self.matrix_dim, self.base.symbol())
        };
        let core = if self.doubled {
            format!("{unit} ⊕ {unit}")
        } else {
            unit
        };
        if self.radical_dim > 0 {
            format!("{core} ⊗ Λ({}^{})", self.ground.symbol(), self.radical_dim)
        } else {
            core
        }
    }
}

/// `2^k`.
fn p2(k: usize) -> usize {
    1usize << k
}

/// Classify the nondegenerate real Clifford algebra `Cl(p,q)` (no radical) by
/// the 8-fold Bott table. `radical_dim`/`ground` are filled in by the callers.
fn real_core(p: usize, q: usize) -> (BaseField, usize, bool) {
    let n = p + q;
    let s = (q as isize - p as isize).rem_euclid(8) as usize;
    match s {
        0 | 6 => (BaseField::R, p2(n / 2), false),
        1 | 5 => (BaseField::C, p2((n - 1) / 2), false),
        2 | 4 => (BaseField::H, p2((n - 2) / 2), false),
        3 => (BaseField::H, p2((n - 3) / 2), true),
        7 => (BaseField::R, p2((n - 1) / 2), true),
        _ => unreachable!(),
    }
}

/// Classify a real Clifford algebra from its signature `(p, q, r)`.
pub fn classify_real(p: usize, q: usize, r: usize) -> CliffordType {
    let (base, matrix_dim, doubled) = real_core(p, q);
    CliffordType {
        base,
        matrix_dim,
        doubled,
        radical_dim: r,
        ground: BaseField::R,
        signature: (p, q),
    }
}

/// Classify a complex Clifford algebra from `(n, r)` (nondegenerate dim, radical).
pub fn classify_complex(n: usize, r: usize) -> CliffordType {
    let (matrix_dim, doubled) = if n % 2 == 0 {
        (p2(n / 2), false)
    } else {
        (p2((n - 1) / 2), true)
    };
    CliffordType {
        base: BaseField::C,
        matrix_dim,
        doubled,
        radical_dim: r,
        ground: BaseField::C,
        signature: (n, 0),
    }
}

/// Extract `(p, q, r)` from a *diagonal* metric using a sign function. Returns
/// `None` if the metric is non-orthogonal (`b`/`a` nonempty) — diagonalization
/// is not attempted here.
fn signature<S: crate::scalar::Scalar>(
    metric: &Metric<S>,
    sign: impl Fn(&S) -> Ordering,
) -> Option<(usize, usize, usize)> {
    if !metric.b.is_empty() || !metric.a.is_empty() {
        return None;
    }
    let (mut p, mut q, mut r) = (0, 0, 0);
    for x in &metric.q {
        match sign(x) {
            Ordering::Greater => p += 1,
            Ordering::Less => q += 1,
            Ordering::Equal => r += 1,
        }
    }
    Some((p, q, r))
}

/// Classify a rational-scalar Clifford algebra (the validation backend). Diagonal
/// metrics only; `None` otherwise.
pub fn classify_rational(metric: &Metric<Rational>) -> Option<CliffordType> {
    let (p, q, r) = signature(metric, |x| x.sign())?;
    Some(classify_real(p, q, r))
}

/// Classify a surreal-scalar Clifford algebra. The surreals are real-closed, so
/// this is the genuine ℝ-Clifford classification — with metric entries allowed to
/// be infinite (ω) or infinitesimal (ε); only their *sign* matters. Diagonal
/// metrics only.
pub fn classify_surreal(metric: &Metric<Surreal>) -> Option<CliffordType> {
    let (p, q, r) = signature(metric, |x| x.sign())?;
    Some(classify_real(p, q, r))
}

/// Classify a surcomplex-scalar Clifford algebra. The field is algebraically
/// closed, so only nondegenerate-dimension and radical matter (2-fold). Diagonal
/// metrics only.
pub fn classify_surcomplex(metric: &Metric<Surcomplex<Surreal>>) -> Option<CliffordType> {
    if !metric.b.is_empty() || !metric.a.is_empty() {
        return None;
    }
    let nonzero = metric.q.iter().filter(|z| !z.is_zero()).count();
    let r = metric.q.len() - nonzero;
    Some(classify_complex(nonzero, r))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::clifford::{CliffordAlgebra, Metric};
    use crate::scalar::Scalar;

    fn rat(n: i128) -> Rational {
        Rational::int(n)
    }
    fn cl_real(qs: &[i128]) -> Option<CliffordType> {
        classify_rational(&Metric::diagonal(qs.iter().map(|&x| rat(x)).collect()))
    }
    fn name(qs: &[i128]) -> String {
        cl_real(qs).unwrap().display()
    }

    #[test]
    fn low_dimensional_real_clifford_table() {
        assert_eq!(name(&[]), "R"); // Cl(0,0) = ℝ
        assert_eq!(name(&[1]), "R ⊕ R"); // Cl(1,0) = ℝ⊕ℝ
        assert_eq!(name(&[-1]), "C"); // Cl(0,1) = ℂ
        assert_eq!(name(&[1, 1]), "M_2(R)"); // Cl(2,0) = M₂(ℝ)
        assert_eq!(name(&[1, -1]), "M_2(R)"); // Cl(1,1) = M₂(ℝ)
        assert_eq!(name(&[-1, -1]), "H"); // Cl(0,2) = ℍ
        assert_eq!(name(&[1, 1, 1]), "M_2(C)"); // Cl(3,0) = M₂(ℂ)
        assert_eq!(name(&[-1, -1, -1]), "H ⊕ H"); // Cl(0,3) = ℍ⊕ℍ
        assert_eq!(name(&[-1, -1, -1, -1]), "M_2(H)"); // Cl(0,4) = M₂(ℍ)
    }

    #[test]
    fn physics_signatures() {
        // Spacetime algebra Cl(1,3) ≅ M₂(ℍ); Cl(3,1) ≅ M₄(ℝ) (the two conventions
        // are genuinely different algebras — a classic subtlety the table shows).
        assert_eq!(name(&[1, -1, -1, -1]), "M_2(H)"); // Cl(1,3)
        assert_eq!(name(&[1, 1, 1, -1]), "M_4(R)"); // Cl(3,1)
                                                    // Conformal geometric algebra Cl(4,1) ≅ M₄(ℂ).
        assert_eq!(name(&[1, 1, 1, 1, -1]), "M_4(C)"); // Cl(4,1)
    }

    #[test]
    fn dimension_is_consistent() {
        // real-dim of the algebra must equal 2^n for every nondegenerate signature.
        for p in 0..=5usize {
            for q in 0..=5usize {
                let t = classify_real(p, q, 0);
                let unit = match t.base {
                    BaseField::R => 1,
                    BaseField::C => 2,
                    BaseField::H => 4,
                };
                let copies = if t.doubled { 2 } else { 1 };
                let real_dim = copies * unit * t.matrix_dim * t.matrix_dim;
                assert_eq!(real_dim, 1usize << (p + q), "Cl({p},{q})");
            }
        }
    }

    #[test]
    fn radical_gives_exterior_factor() {
        // Cl(0,1,2): ℂ tensor an exterior algebra on the 2 null directions.
        assert_eq!(name(&[-1, 0, 0]), "C ⊗ Λ(R^2)");
        // pure Grassmann Λ(R^3) = Cl(0,0,3): trivial core ⊗ Λ.
        assert_eq!(name(&[0, 0, 0]), "R ⊗ Λ(R^3)");
    }

    #[test]
    fn surreal_signs_classify_by_sign_only() {
        // Infinite/infinitesimal squares, but the signature is (1,1): Cl(1,1)=M₂(ℝ).
        let m = Metric::diagonal(vec![Surreal::omega(), Surreal::epsilon().neg()]);
        assert_eq!(classify_surreal(&m).unwrap().display(), "M_2(R)");
    }

    #[test]
    fn surcomplex_is_two_fold() {
        let even =
            Metric::<Surcomplex<Surreal>>::diagonal(vec![Surcomplex::one(), Surcomplex::one()]);
        assert_eq!(classify_surcomplex(&even).unwrap().display(), "M_2(C)"); // n=2
        let odd = Metric::<Surcomplex<Surreal>>::diagonal(vec![Surcomplex::one()]);
        assert_eq!(classify_surcomplex(&odd).unwrap().display(), "C ⊕ C"); // n=1
    }

    #[test]
    fn even_subalgebra_classification_drops_one_dimension() {
        // Cl(3,0)⁰ ≅ Cl(0,2) = ℍ — ties the classifier to even_subalgebra.
        let alg = CliffordAlgebra::new(3, Metric::diagonal(vec![rat(1), rat(1), rat(1)]));
        let even = alg.even_subalgebra().unwrap();
        assert_eq!(classify_rational(&even.metric).unwrap().display(), "H");
        // Cl(1,3)⁰ ≅ Cl(1,2) ... check it matches a direct signature classification.
        let st = CliffordAlgebra::new(4, Metric::diagonal(vec![rat(1), rat(-1), rat(-1), rat(-1)]));
        let st_even = st.even_subalgebra().unwrap();
        // pivot is the last non-null (a −1 direction): f_i² = −q_i·(−1) = q_i.
        // signature of the even part here is (1,2) ⇒ same class as Cl(1,2).
        assert_eq!(
            classify_rational(&st_even.metric).unwrap().display(),
            classify_real(1, 2, 0).display()
        );
    }
}
