//! Python bindings for the form classifiers across the characteristic
//! trichotomy: Arf (char 2), the char-0 Clifford type, the Witt classes,
//! Dickson, the odd-characteristic classifier, and the Springer decomposition.
//! These consume the `pub(crate)` algebra types stamped by [`super::engine`].

use super::engine::{NimberAlgebra, NimberMV, SurcomplexAlgebra, SurrealAlgebra};
use crate::clifford::Metric;
use crate::forms::{WittClass, WittClassG};
use crate::scalar::Fp;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

#[pyclass(name = "ArfResult", module = "pleroma")]
struct PyArfResult {
    inner: crate::forms::ArfResult,
}

#[pymethods]
impl PyArfResult {
    #[getter]
    fn arf(&self) -> u8 {
        self.inner.arf
    }
    #[getter]
    fn rank(&self) -> usize {
        self.inner.rank
    }
    #[getter]
    fn radical_dim(&self) -> usize {
        self.inner.radical_dim
    }
    #[getter]
    fn radical_anisotropic(&self) -> bool {
        self.inner.radical_anisotropic
    }
    #[getter]
    fn o_type(&self) -> &'static str {
        self.inner.o_type
    }
    fn __repr__(&self) -> String {
        format!(
            "ArfResult(arf={}, type={}, rank={}, radical_dim={}, radical_anisotropic={})",
            self.inner.arf,
            self.inner.o_type,
            self.inner.rank,
            self.inner.radical_dim,
            self.inner.radical_anisotropic,
        )
    }
}

/// Arf invariant (the char-2 Clifford classifier) of a nimber algebra whose
/// metric has F₂ entries.
#[pyfunction]
fn arf_invariant(alg: &NimberAlgebra) -> PyResult<PyArfResult> {
    let inner = crate::forms::arf_invariant(&alg.inner.metric).ok_or_else(|| {
        PyValueError::new_err("Arf invariant is undefined for general-bilinear metrics")
    })?;
    Ok(PyArfResult { inner })
}
// ---------------------------------------------------------------------------
// Char-0 classifier
// ---------------------------------------------------------------------------

#[pyclass(name = "CliffordType", module = "pleroma")]
struct PyCliffordType {
    inner: crate::forms::CliffordType,
}

#[pymethods]
impl PyCliffordType {
    #[getter]
    fn base(&self) -> String {
        format!("{:?}", self.inner.base)
    }
    #[getter]
    fn matrix_dim(&self) -> usize {
        self.inner.matrix_dim
    }
    #[getter]
    fn doubled(&self) -> bool {
        self.inner.doubled
    }
    #[getter]
    fn radical_dim(&self) -> usize {
        self.inner.radical_dim
    }
    #[getter]
    fn signature(&self) -> (usize, usize) {
        self.inner.signature
    }
    fn __repr__(&self) -> String {
        self.inner.display()
    }
}

/// Classify a surreal Clifford algebra (the genuine real classification) as a
/// matrix algebra over ℝ/ℂ/ℍ. Diagonal metrics only.
#[pyfunction]
fn classify_surreal(alg: &SurrealAlgebra) -> PyResult<PyCliffordType> {
    crate::forms::classify_surreal(&alg.inner.metric)
        .map(|t| PyCliffordType { inner: t })
        .ok_or_else(|| PyValueError::new_err("classifier needs a diagonal (orthogonal) metric"))
}

/// Classify a surcomplex Clifford algebra (the 2-fold complex classification).
/// Diagonal metrics only.
#[pyfunction]
fn classify_surcomplex(alg: &SurcomplexAlgebra) -> PyResult<PyCliffordType> {
    crate::forms::classify_surcomplex(&alg.inner.metric)
        .map(|t| PyCliffordType { inner: t })
        .ok_or_else(|| PyValueError::new_err("classifier needs a diagonal (orthogonal) metric"))
}

// ---------------------------------------------------------------------------
// Witt group + Dickson invariant
// ---------------------------------------------------------------------------

#[pyclass(name = "WittClass", module = "pleroma", from_py_object)]
#[derive(Clone)]
struct PyWittClass {
    inner: WittClass,
}

#[pymethods]
impl PyWittClass {
    #[getter]
    fn arf(&self) -> u8 {
        self.inner.arf
    }
    fn add(&self, other: &PyWittClass) -> PyWittClass {
        PyWittClass {
            inner: self.inner.add(&other.inner),
        }
    }
    fn __add__(&self, other: &PyWittClass) -> PyWittClass {
        self.add(other)
    }
    fn is_hyperbolic(&self) -> bool {
        self.inner.is_hyperbolic()
    }
    fn anisotropic_dim(&self) -> usize {
        self.inner.anisotropic_dim()
    }
    fn __eq__(&self, other: &PyWittClass) -> bool {
        self.inner == other.inner
    }
    fn __repr__(&self) -> String {
        self.inner.display()
    }
}

/// The Witt class (in `W_q ≅ ℤ/2`) of a nimber Clifford metric.
#[pyfunction]
fn witt_class(alg: &NimberAlgebra) -> PyWittClass {
    PyWittClass {
        inner: WittClass::from_metric(&alg.inner.metric),
    }
}

/// The Dickson invariant of an orthogonal matrix over the nim-field (the char-2
/// determinant replacement; `0` ⇒ rotation/SO, `1` ⇒ reflection).
#[pyfunction]
fn dickson_matrix(g: Vec<Vec<u128>>) -> u8 {
    crate::forms::dickson_matrix(&g)
}

/// The Dickson invariant of a nimber Clifford versor (= its grade parity).
#[pyfunction]
fn dickson_of_versor(v: &NimberMV) -> PyResult<u8> {
    crate::forms::dickson_of_versor(&v.alg, &v.mv)
        .ok_or_else(|| PyValueError::new_err("not an invertible homogeneous versor"))
}
// ---------------------------------------------------------------------------
// Odd-characteristic classifier (the trichotomy's third leg)
// ---------------------------------------------------------------------------

fn fp_diag<const P: u128>(q: &[i128]) -> Metric<Fp<P>> {
    Metric::diagonal(q.iter().map(|&x| Fp::<P>::new(x)).collect())
}

#[pyclass(name = "OddCharType", module = "pleroma")]
struct PyOddCharType {
    inner: crate::forms::OddCharType,
}

#[pymethods]
impl PyOddCharType {
    #[getter]
    fn p(&self) -> u128 {
        self.inner.p
    }
    #[getter]
    fn dim(&self) -> usize {
        self.inner.dim
    }
    #[getter]
    fn radical_dim(&self) -> usize {
        self.inner.radical_dim
    }
    #[getter]
    fn disc_is_square(&self) -> bool {
        self.inner.disc_is_square
    }
    #[getter]
    fn hasse(&self) -> i8 {
        self.inner.hasse
    }
    fn __repr__(&self) -> String {
        self.inner.display()
    }
}

#[pyclass(name = "WittClassG", module = "pleroma", from_py_object)]
#[derive(Clone)]
struct PyWittClassG {
    inner: WittClassG,
}

#[pymethods]
impl PyWittClassG {
    fn add(&self, other: &PyWittClassG) -> PyWittClassG {
        PyWittClassG {
            inner: self.inner.add(&other.inner),
        }
    }
    fn __add__(&self, other: &PyWittClassG) -> PyWittClassG {
        self.add(other)
    }
    fn __eq__(&self, other: &PyWittClassG) -> bool {
        self.inner == other.inner
    }
    fn __repr__(&self) -> String {
        match self.inner {
            WittClassG::Char0 { signature } => format!("WittClassG::Char0(signature={signature})"),
            WittClassG::OddChar { kappa, e0, sclass } => {
                format!("WittClassG::OddChar(kappa={kappa}, e0={e0}, sclass={sclass})")
            }
            WittClassG::Char2 { arf } => format!("WittClassG::Char2(arf={arf})"),
        }
    }
}

/// Classify a diagonal odd-characteristic form `q` over `F_p` (dimension +
/// discriminant + Hasse). Supported primes: 3, 5, 7, 11, 13.
#[pyfunction]
fn classify_oddchar(p: u128, q: Vec<i128>) -> PyResult<PyOddCharType> {
    let res = match p {
        3 => crate::forms::classify_oddchar(&fp_diag::<3>(&q)),
        5 => crate::forms::classify_oddchar(&fp_diag::<5>(&q)),
        7 => crate::forms::classify_oddchar(&fp_diag::<7>(&q)),
        11 => crate::forms::classify_oddchar(&fp_diag::<11>(&q)),
        13 => crate::forms::classify_oddchar(&fp_diag::<13>(&q)),
        _ => return Err(PyValueError::new_err("supported primes: 3, 5, 7, 11, 13")),
    };
    res.map(|t| PyOddCharType { inner: t })
        .ok_or_else(|| PyValueError::new_err("non-diagonal metric"))
}

/// The odd-characteristic Witt class of a diagonal form `q` over `F_p`.
#[pyfunction]
fn oddchar_witt(p: u128, q: Vec<i128>) -> PyResult<PyWittClassG> {
    let res = match p {
        3 => crate::forms::oddchar_witt(&fp_diag::<3>(&q)),
        5 => crate::forms::oddchar_witt(&fp_diag::<5>(&q)),
        7 => crate::forms::oddchar_witt(&fp_diag::<7>(&q)),
        11 => crate::forms::oddchar_witt(&fp_diag::<11>(&q)),
        13 => crate::forms::oddchar_witt(&fp_diag::<13>(&q)),
        _ => return Err(PyValueError::new_err("supported primes: 3, 5, 7, 11, 13")),
    };
    res.map(|w| PyWittClassG { inner: w })
        .ok_or_else(|| PyValueError::new_err("non-diagonal metric"))
}

/// Is `x` a square mod `p`? (Euler's criterion.) Supported primes: 3, 5, 7, 11, 13.
#[pyfunction]
fn is_square_mod(p: u128, x: i128) -> PyResult<bool> {
    Ok(match p {
        3 => crate::forms::is_square(Fp::<3>::new(x)),
        5 => crate::forms::is_square(Fp::<5>::new(x)),
        7 => crate::forms::is_square(Fp::<7>::new(x)),
        11 => crate::forms::is_square(Fp::<11>::new(x)),
        13 => crate::forms::is_square(Fp::<13>::new(x)),
        _ => return Err(PyValueError::new_err("supported primes: 3, 5, 7, 11, 13")),
    })
}

/// The Hasse–Witt invariant of a diagonal form `q` over `F_p` (always +1 over a
/// finite field). Supported primes: 3, 5, 7, 11, 13.
#[pyfunction]
fn hasse_invariant(p: u128, q: Vec<i128>) -> PyResult<i8> {
    let res = match p {
        3 => crate::forms::hasse_invariant(&fp_diag::<3>(&q)),
        5 => crate::forms::hasse_invariant(&fp_diag::<5>(&q)),
        7 => crate::forms::hasse_invariant(&fp_diag::<7>(&q)),
        11 => crate::forms::hasse_invariant(&fp_diag::<11>(&q)),
        13 => crate::forms::hasse_invariant(&fp_diag::<13>(&q)),
        _ => return Err(PyValueError::new_err("supported primes: 3, 5, 7, 11, 13")),
    };
    res.ok_or_else(|| PyValueError::new_err("non-diagonal metric"))
}

/// Witt decomposition of a diagonal odd-char form `q` over `F_p`: returns
/// `(witt_index, anisotropic_dim, anisotropic_disc_is_square, radical_dim)`.
/// Supported primes: 3, 5, 7, 11, 13.
#[pyfunction]
fn witt_decompose_oddchar(p: u128, q: Vec<i128>) -> PyResult<(usize, usize, bool, usize)> {
    let d = match p {
        3 => crate::forms::witt_decompose_oddchar(&fp_diag::<3>(&q)),
        5 => crate::forms::witt_decompose_oddchar(&fp_diag::<5>(&q)),
        7 => crate::forms::witt_decompose_oddchar(&fp_diag::<7>(&q)),
        11 => crate::forms::witt_decompose_oddchar(&fp_diag::<11>(&q)),
        13 => crate::forms::witt_decompose_oddchar(&fp_diag::<13>(&q)),
        _ => return Err(PyValueError::new_err("supported primes: 3, 5, 7, 11, 13")),
    }
    .ok_or_else(|| PyValueError::new_err("non-diagonal metric"))?;
    Ok((
        d.witt_index,
        d.anisotropic_dim,
        d.anisotropic_disc_is_square,
        d.radical_dim,
    ))
}

/// Are two diagonal odd-char forms over `F_p` isometric? `(dim, discriminant)`
/// is a complete invariant. Supported primes: 3, 5, 7, 11, 13.
#[pyfunction]
fn is_isometric_oddchar(p: u128, q1: Vec<i128>, q2: Vec<i128>) -> PyResult<bool> {
    let res = match p {
        3 => crate::forms::isometric_oddchar(&fp_diag::<3>(&q1), &fp_diag::<3>(&q2)),
        5 => crate::forms::isometric_oddchar(&fp_diag::<5>(&q1), &fp_diag::<5>(&q2)),
        7 => crate::forms::isometric_oddchar(&fp_diag::<7>(&q1), &fp_diag::<7>(&q2)),
        11 => crate::forms::isometric_oddchar(&fp_diag::<11>(&q1), &fp_diag::<11>(&q2)),
        13 => crate::forms::isometric_oddchar(&fp_diag::<13>(&q1), &fp_diag::<13>(&q2)),
        _ => return Err(PyValueError::new_err("supported primes: 3, 5, 7, 11, 13")),
    };
    res.ok_or_else(|| PyValueError::new_err("non-diagonal metric"))
}

/// Witt decomposition of a real (surreal) form: returns `(witt_index,
/// anisotropic_pos, anisotropic_neg, radical_dim)`.
#[pyfunction]
fn witt_decompose_real(alg: &SurrealAlgebra) -> PyResult<(usize, usize, usize, usize)> {
    let d = crate::forms::witt_decompose_real(&alg.inner.metric)
        .ok_or_else(|| PyValueError::new_err("could not diagonalize the metric"))?;
    Ok((
        d.witt_index,
        d.anisotropic_pos,
        d.anisotropic_neg,
        d.radical_dim,
    ))
}

// ---------------------------------------------------------------------------
// Non-Archimedean Springer decomposition (surreal)
// ---------------------------------------------------------------------------

#[pyclass(name = "SpringerDecomp", module = "pleroma")]
struct PySpringerDecomp {
    #[pyo3(get)]
    graded: Vec<(String, (usize, usize))>,
    #[pyo3(get)]
    radical_dim: usize,
    #[pyo3(get)]
    total_signature: (usize, usize),
}

#[pymethods]
impl PySpringerDecomp {
    fn __repr__(&self) -> String {
        format!(
            "SpringerDecomp(graded={:?}, radical_dim={}, total_signature={:?})",
            self.graded, self.radical_dim, self.total_signature
        )
    }
}

/// The non-Archimedean Springer decomposition of a diagonal surreal form: its
/// ω-adic valuation filtration into residue ℝ-signatures.
#[pyfunction]
fn springer_decompose(alg: &SurrealAlgebra) -> PyResult<PySpringerDecomp> {
    let d = crate::forms::springer_decompose(&alg.inner.metric)
        .ok_or_else(|| PyValueError::new_err("Springer decomposition needs a diagonal metric"))?;
    let graded = d
        .graded
        .iter()
        .map(|rf| (format!("{:?}", rf.valuation), rf.signature))
        .collect();
    Ok(PySpringerDecomp {
        graded,
        radical_dim: d.radical_dim,
        total_signature: d.total_signature,
    })
}

// ---------------------------------------------------------------------------
// Witt ring + cohomological invariant staircase (eₙ)
// ---------------------------------------------------------------------------

/// The cohomological invariant staircase `(e₀, e₁, e₂)` of a diagonal odd-char
/// form `q` over `F_p`, with the field's stabilization (`I² = 0` over a finite
/// field ⇒ `stabilizes_at = 2`, `e₂ = +1`). Returns `(e0, e1, e2, stabilizes_at)`.
/// Supported primes: 3, 5, 7, 11, 13.
#[pyfunction]
fn e_staircase_oddchar(p: u128, q: Vec<i128>) -> PyResult<(u8, u8, i8, usize)> {
    let s = match p {
        3 => crate::forms::e_staircase_oddchar(&fp_diag::<3>(&q)),
        5 => crate::forms::e_staircase_oddchar(&fp_diag::<5>(&q)),
        7 => crate::forms::e_staircase_oddchar(&fp_diag::<7>(&q)),
        11 => crate::forms::e_staircase_oddchar(&fp_diag::<11>(&q)),
        13 => crate::forms::e_staircase_oddchar(&fp_diag::<13>(&q)),
        _ => return Err(PyValueError::new_err("supported primes: 3, 5, 7, 11, 13")),
    }
    .ok_or_else(|| PyValueError::new_err("non-diagonal metric"))?;
    Ok((s.e0, s.e1, s.e2, s.stabilizes_at))
}

/// The real cohomological invariant `eₙ` of a form of signature `σ` over `ℝ`:
/// `Some((σ/2ⁿ) mod 2)` if the form is in `Iⁿ` (i.e. `2ⁿ | σ`), else `None`. The
/// staircase reads the 2-adic expansion of the signature (the infinite ℝ tower).
#[pyfunction]
fn e_real(signature: i128, n: usize) -> Option<u8> {
    crate::forms::e_real(signature, n)
}

// ---------------------------------------------------------------------------
// p-adic Hilbert symbol + Hasse–Minkowski over Q
// ---------------------------------------------------------------------------

/// The Hilbert symbol `(a, b)_p` over `Q_p` (`p`-adic). Unlike the finite-field
/// Hilbert symbol (always `+1`), this is genuinely nontrivial — e.g. `(−1,−1)_2 = −1`.
#[pyfunction]
fn hilbert_symbol_qp(a: i128, b: i128, p: u128) -> i8 {
    crate::forms::hilbert_symbol_qp(a, b, p)
}

/// The Hilbert symbol `(a, b)_∞` over `ℝ` (`−1` iff both are negative).
#[pyfunction]
fn hilbert_symbol_real(a: i128, b: i128) -> i8 {
    crate::forms::hilbert_symbol_real(a, b)
}

/// Is the integer `n` a square in `Q_p`?
#[pyfunction]
fn is_square_qp(n: i128, p: u128) -> bool {
    crate::forms::is_square_qp(n, p)
}

/// Is the diagonal rational/integer form `⟨a₁,…,aₙ⟩` isotropic over `Q`? By the
/// **Hasse–Minkowski** principle (isotropic over `ℝ` and every `Q_p`). E.g.
/// `⟨1,1,1⟩` is anisotropic, `⟨1,1,-1⟩` isotropic, `⟨1,1,-3⟩` anisotropic.
#[pyfunction]
fn is_isotropic_q(entries: Vec<i128>) -> bool {
    crate::forms::is_isotropic_q(&entries)
}

// ---------------------------------------------------------------------------
// Brauer–Wall group
// ---------------------------------------------------------------------------

#[pyclass(name = "BrauerWallClass", module = "pleroma", from_py_object)]
#[derive(Clone)]
struct PyBrauerWallClass {
    inner: crate::forms::BrauerWallClass,
}

#[pymethods]
impl PyBrauerWallClass {
    fn add(&self, other: &PyBrauerWallClass) -> PyBrauerWallClass {
        PyBrauerWallClass {
            inner: self.inner.add(&other.inner),
        }
    }
    fn __add__(&self, other: &PyBrauerWallClass) -> PyBrauerWallClass {
        self.add(other)
    }
    fn __eq__(&self, other: &PyBrauerWallClass) -> bool {
        self.inner == other.inner
    }
    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}

/// The Brauer–Wall class of a surreal (real) Clifford algebra: the Bott index
/// `s = (q − p) mod 8` in `BW(ℝ) ≅ ℤ/8` — the periodicity table as a group element.
#[pyfunction]
fn bw_class_real(alg: &SurrealAlgebra) -> PyResult<PyBrauerWallClass> {
    crate::forms::bw_class_real(&alg.inner.metric)
        .map(|c| PyBrauerWallClass { inner: c })
        .ok_or_else(|| PyValueError::new_err("Brauer–Wall class needs a diagonal metric"))
}

/// The Brauer–Wall class of a surcomplex (complex) Clifford algebra in
/// `BW(ℂ) ≅ ℤ/2` (the dimension parity).
#[pyfunction]
fn bw_class_complex(alg: &SurcomplexAlgebra) -> PyResult<PyBrauerWallClass> {
    crate::forms::bw_class_complex(&alg.inner.metric)
        .map(|c| PyBrauerWallClass { inner: c })
        .ok_or_else(|| PyValueError::new_err("Brauer–Wall class needs a diagonal metric"))
}

/// The Brauer–Wall class of a diagonal odd-char form `q` over `F_p` (the order-4
/// graded part, `BW(F_q) ≅ W(F_q)`). Supported primes: 3, 5, 7, 11, 13.
#[pyfunction]
fn bw_class_oddchar(p: u128, q: Vec<i128>) -> PyResult<PyBrauerWallClass> {
    let res = match p {
        3 => crate::forms::bw_class_oddchar(&fp_diag::<3>(&q)),
        5 => crate::forms::bw_class_oddchar(&fp_diag::<5>(&q)),
        7 => crate::forms::bw_class_oddchar(&fp_diag::<7>(&q)),
        11 => crate::forms::bw_class_oddchar(&fp_diag::<11>(&q)),
        13 => crate::forms::bw_class_oddchar(&fp_diag::<13>(&q)),
        _ => return Err(PyValueError::new_err("supported primes: 3, 5, 7, 11, 13")),
    };
    res.map(|c| PyBrauerWallClass { inner: c })
        .ok_or_else(|| PyValueError::new_err("non-diagonal metric"))
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyArfResult>()?;
    m.add_class::<PyCliffordType>()?;
    m.add_class::<PyWittClass>()?;
    m.add_class::<PyOddCharType>()?;
    m.add_class::<PyWittClassG>()?;
    m.add_class::<PySpringerDecomp>()?;
    m.add_class::<PyBrauerWallClass>()?;
    m.add_function(wrap_pyfunction!(arf_invariant, m)?)?;
    m.add_function(wrap_pyfunction!(classify_surreal, m)?)?;
    m.add_function(wrap_pyfunction!(classify_surcomplex, m)?)?;
    m.add_function(wrap_pyfunction!(witt_class, m)?)?;
    m.add_function(wrap_pyfunction!(dickson_matrix, m)?)?;
    m.add_function(wrap_pyfunction!(dickson_of_versor, m)?)?;
    m.add_function(wrap_pyfunction!(classify_oddchar, m)?)?;
    m.add_function(wrap_pyfunction!(oddchar_witt, m)?)?;
    m.add_function(wrap_pyfunction!(witt_decompose_oddchar, m)?)?;
    m.add_function(wrap_pyfunction!(witt_decompose_real, m)?)?;
    m.add_function(wrap_pyfunction!(is_isometric_oddchar, m)?)?;
    m.add_function(wrap_pyfunction!(is_square_mod, m)?)?;
    m.add_function(wrap_pyfunction!(hasse_invariant, m)?)?;
    m.add_function(wrap_pyfunction!(springer_decompose, m)?)?;
    m.add_function(wrap_pyfunction!(e_staircase_oddchar, m)?)?;
    m.add_function(wrap_pyfunction!(e_real, m)?)?;
    m.add_function(wrap_pyfunction!(hilbert_symbol_qp, m)?)?;
    m.add_function(wrap_pyfunction!(hilbert_symbol_real, m)?)?;
    m.add_function(wrap_pyfunction!(is_square_qp, m)?)?;
    m.add_function(wrap_pyfunction!(is_isotropic_q, m)?)?;
    m.add_function(wrap_pyfunction!(bw_class_real, m)?)?;
    m.add_function(wrap_pyfunction!(bw_class_complex, m)?)?;
    m.add_function(wrap_pyfunction!(bw_class_oddchar, m)?)?;
    Ok(())
}
