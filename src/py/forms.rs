//! Python bindings for the form classifiers across the characteristic
//! trichotomy: Arf (char 2), the char-0 Clifford type, the Witt classes,
//! Dickson, the odd-characteristic classifier, and the Springer decomposition.
//! These consume the `pub(crate)` algebra types stamped by [`super::engine`].

use super::engine::{NimberAlgebra, NimberMV, SurcomplexAlgebra, SurrealAlgebra};
use crate::clifford::Metric;
use crate::forms::{FiniteOddField, WittClass, WittClassG};
use crate::scalar::{Fp, Fpn, Rational};
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

/// Classify a surreal Clifford algebra on the exact-square real-table subdomain
/// as a matrix algebra over ℝ/ℂ/ℍ. Symmetric metrics are diagonalized when possible.
#[pyfunction]
fn classify_surreal(alg: &SurrealAlgebra) -> PyResult<PyCliffordType> {
    crate::forms::classify_surreal(&alg.inner.metric)
        .map(|t| PyCliffordType { inner: t })
        .ok_or_else(|| {
            PyValueError::new_err(
                "classifier could not diagonalize this metric or needs an unrepresented square root",
            )
        })
}

/// Classify a surcomplex Clifford algebra on the exact-square complex-table
/// subdomain. Symmetric metrics are diagonalized when possible.
#[pyfunction]
fn classify_surcomplex(alg: &SurcomplexAlgebra) -> PyResult<PyCliffordType> {
    crate::forms::classify_surcomplex(&alg.inner.metric)
        .map(|t| PyCliffordType { inner: t })
        .ok_or_else(|| {
            PyValueError::new_err(
                "classifier could not diagonalize this metric or needs an unrepresented square root",
            )
        })
}

/// Classify a real Clifford algebra directly from its signature `(p, q, r)`
/// (`p` plus-squares, `q` minus-squares, `r` null/radical dimensions) — the
/// 8-fold table, no metric needed. Complement to `classify_surreal`.
#[pyfunction]
#[pyo3(signature = (p, q, r=0))]
fn classify_real(p: usize, q: usize, r: usize) -> PyCliffordType {
    PyCliffordType {
        inner: crate::forms::classify_real(p, q, r),
    }
}

/// Classify a complex Clifford algebra directly from `(n, r)` (`n` nondegenerate
/// dimensions, `r` null/radical) — the 2-fold table. Complement to
/// `classify_surcomplex`.
#[pyfunction]
#[pyo3(signature = (n, r=0))]
fn classify_complex(n: usize, r: usize) -> PyCliffordType {
    PyCliffordType {
        inner: crate::forms::classify_complex(n, r),
    }
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
fn witt_class(alg: &NimberAlgebra) -> PyResult<PyWittClass> {
    WittClass::try_from_metric(&alg.inner.metric)
        .map(|inner| PyWittClass { inner })
        .map_err(|err| PyValueError::new_err(format!("Witt class is undefined: {err:?}")))
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

fn finite_diag<F: FiniteOddField>(q: &[i128]) -> Metric<F> {
    Metric::diagonal(q.iter().map(|&x| F::from_i128(x)).collect())
}

fn unsupported_finite_field_err() -> PyErr {
    PyValueError::new_err("supported odd finite fields: F_3, F_5, F_7, F_11, F_13, F_9, F_25, F_27")
}

fn finite_field_order(p: u128, degree: usize) -> PyResult<u128> {
    match (p, degree) {
        (3, 1) => Ok(3),
        (5, 1) => Ok(5),
        (7, 1) => Ok(7),
        (11, 1) => Ok(11),
        (13, 1) => Ok(13),
        (3, 2) => Ok(9),
        (5, 2) => Ok(25),
        (3, 3) => Ok(27),
        _ => Err(unsupported_finite_field_err()),
    }
}

macro_rules! with_finite_odd_metric {
    ($p:expr, $degree:expr, $q:expr, |$metric:ident| $body:expr) => {{
        match ($p, $degree) {
            (3, 1) => {
                let $metric = finite_diag::<Fp<3>>($q);
                $body
            }
            (5, 1) => {
                let $metric = finite_diag::<Fp<5>>($q);
                $body
            }
            (7, 1) => {
                let $metric = finite_diag::<Fp<7>>($q);
                $body
            }
            (11, 1) => {
                let $metric = finite_diag::<Fp<11>>($q);
                $body
            }
            (13, 1) => {
                let $metric = finite_diag::<Fp<13>>($q);
                $body
            }
            (3, 2) => {
                let $metric = finite_diag::<Fpn<3, 2>>($q);
                $body
            }
            (5, 2) => {
                let $metric = finite_diag::<Fpn<5, 2>>($q);
                $body
            }
            (3, 3) => {
                let $metric = finite_diag::<Fpn<3, 3>>($q);
                $body
            }
            _ => return Err(unsupported_finite_field_err()),
        }
    }};
}

macro_rules! with_finite_odd_metrics {
    ($p:expr, $degree:expr, $q1:expr, $q2:expr, |$m1:ident, $m2:ident| $body:expr) => {{
        match ($p, $degree) {
            (3, 1) => {
                let $m1 = finite_diag::<Fp<3>>($q1);
                let $m2 = finite_diag::<Fp<3>>($q2);
                $body
            }
            (5, 1) => {
                let $m1 = finite_diag::<Fp<5>>($q1);
                let $m2 = finite_diag::<Fp<5>>($q2);
                $body
            }
            (7, 1) => {
                let $m1 = finite_diag::<Fp<7>>($q1);
                let $m2 = finite_diag::<Fp<7>>($q2);
                $body
            }
            (11, 1) => {
                let $m1 = finite_diag::<Fp<11>>($q1);
                let $m2 = finite_diag::<Fp<11>>($q2);
                $body
            }
            (13, 1) => {
                let $m1 = finite_diag::<Fp<13>>($q1);
                let $m2 = finite_diag::<Fp<13>>($q2);
                $body
            }
            (3, 2) => {
                let $m1 = finite_diag::<Fpn<3, 2>>($q1);
                let $m2 = finite_diag::<Fpn<3, 2>>($q2);
                $body
            }
            (5, 2) => {
                let $m1 = finite_diag::<Fpn<5, 2>>($q1);
                let $m2 = finite_diag::<Fpn<5, 2>>($q2);
                $body
            }
            (3, 3) => {
                let $m1 = finite_diag::<Fpn<3, 3>>($q1);
                let $m2 = finite_diag::<Fpn<3, 3>>($q2);
                $body
            }
            _ => return Err(unsupported_finite_field_err()),
        }
    }};
}

macro_rules! with_finite_odd_value {
    ($p:expr, $degree:expr, $x:expr, |$value:ident| $body:expr) => {{
        match ($p, $degree) {
            (3, 1) => {
                let $value = <Fp<3> as FiniteOddField>::from_i128($x);
                $body
            }
            (5, 1) => {
                let $value = <Fp<5> as FiniteOddField>::from_i128($x);
                $body
            }
            (7, 1) => {
                let $value = <Fp<7> as FiniteOddField>::from_i128($x);
                $body
            }
            (11, 1) => {
                let $value = <Fp<11> as FiniteOddField>::from_i128($x);
                $body
            }
            (13, 1) => {
                let $value = <Fp<13> as FiniteOddField>::from_i128($x);
                $body
            }
            (3, 2) => {
                let $value = <Fpn<3, 2> as FiniteOddField>::from_i128($x);
                $body
            }
            (5, 2) => {
                let $value = <Fpn<5, 2> as FiniteOddField>::from_i128($x);
                $body
            }
            (3, 3) => {
                let $value = <Fpn<3, 3> as FiniteOddField>::from_i128($x);
                $body
            }
            _ => return Err(unsupported_finite_field_err()),
        }
    }};
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
    fn field_order(&self) -> u128 {
        self.inner.field_order
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

#[pyclass(name = "FiniteFieldForm", module = "pleroma", from_py_object)]
#[derive(Clone)]
struct PyFiniteFieldForm {
    p: u128,
    degree: usize,
    q: Vec<i128>,
}

#[pymethods]
impl PyFiniteFieldForm {
    #[new]
    #[pyo3(signature = (p, q, degree=1))]
    fn new(p: u128, q: Vec<i128>, degree: usize) -> PyResult<Self> {
        finite_field_order(p, degree)?;
        Ok(PyFiniteFieldForm { p, degree, q })
    }

    #[getter]
    fn p(&self) -> u128 {
        self.p
    }

    #[getter]
    fn degree(&self) -> usize {
        self.degree
    }

    #[getter]
    fn field_order(&self) -> PyResult<u128> {
        finite_field_order(self.p, self.degree)
    }

    #[getter]
    fn diagonal(&self) -> Vec<i128> {
        self.q.clone()
    }

    fn classify(&self) -> PyResult<PyOddCharType> {
        let res = with_finite_odd_metric!(self.p, self.degree, &self.q, |m| {
            crate::forms::classify_finite_odd(&m)
        });
        res.map(|inner| PyOddCharType { inner })
            .ok_or_else(|| PyValueError::new_err("non-diagonal metric"))
    }

    fn witt_class(&self) -> PyResult<PyWittClassG> {
        let res = with_finite_odd_metric!(self.p, self.degree, &self.q, |m| {
            crate::forms::finite_odd_witt(&m)
        });
        res.map(|inner| PyWittClassG { inner })
            .ok_or_else(|| PyValueError::new_err("non-diagonal metric"))
    }

    fn witt_decompose(&self) -> PyResult<(usize, usize, bool, usize)> {
        let d = with_finite_odd_metric!(self.p, self.degree, &self.q, |m| {
            crate::forms::witt_decompose_finite_odd(&m)
        })
        .ok_or_else(|| PyValueError::new_err("non-diagonal metric"))?;
        Ok((
            d.witt_index,
            d.anisotropic_dim,
            d.anisotropic_disc_is_square,
            d.radical_dim,
        ))
    }

    fn is_isometric(&self, other: &PyFiniteFieldForm) -> PyResult<bool> {
        if self.p != other.p || self.degree != other.degree {
            return Err(PyValueError::new_err(
                "isometry needs both forms over the same finite field",
            ));
        }
        with_finite_odd_metrics!(self.p, self.degree, &self.q, &other.q, |m1, m2| {
            crate::forms::isometric_finite_odd(&m1, &m2)
                .ok_or_else(|| PyValueError::new_err("non-diagonal metric"))
        })
    }

    fn is_square(&self, x: i128) -> PyResult<bool> {
        Ok(with_finite_odd_value!(self.p, self.degree, x, |value| {
            crate::forms::is_square_finite(value)
        }))
    }

    fn hasse_invariant(&self) -> PyResult<i8> {
        with_finite_odd_metric!(self.p, self.degree, &self.q, |m| {
            crate::forms::hasse_invariant_finite_odd(&m)
                .ok_or_else(|| PyValueError::new_err("non-diagonal metric"))
        })
    }

    fn e_staircase(&self) -> PyResult<(u8, u8, i8, usize)> {
        let s = with_finite_odd_metric!(self.p, self.degree, &self.q, |m| {
            crate::forms::e_staircase_finite_odd(&m)
        })
        .ok_or_else(|| PyValueError::new_err("non-diagonal metric"))?;
        Ok((s.e0, s.e1, s.e2, s.stabilizes_at))
    }

    fn bw_class(&self) -> PyResult<PyBrauerWallClass> {
        let res = with_finite_odd_metric!(self.p, self.degree, &self.q, |m| {
            crate::forms::bw_class_finite_odd(&m)
        });
        res.map(|inner| PyBrauerWallClass { inner })
            .ok_or_else(|| PyValueError::new_err("non-diagonal metric"))
    }

    fn __repr__(&self) -> PyResult<String> {
        Ok(format!(
            "FiniteFieldForm(F_{}, diagonal={:?})",
            self.field_order()?,
            self.q
        ))
    }
}

#[pyclass(name = "WittClassG", module = "pleroma", from_py_object)]
#[derive(Clone)]
struct PyWittClassG {
    inner: WittClassG,
}

#[pymethods]
impl PyWittClassG {
    fn add(&self, other: &PyWittClassG) -> PyResult<PyWittClassG> {
        self.inner
            .try_add(&other.inner)
            .map(|inner| PyWittClassG { inner })
            .map_err(PyValueError::new_err)
    }
    fn __add__(&self, other: &PyWittClassG) -> PyResult<PyWittClassG> {
        self.add(other)
    }
    /// The Witt-**ring** product (tensor of forms). Defined on the char-0 and
    /// odd-char legs; panics on a char-2 operand (`W_q` is a module, not a ring).
    fn mul(&self, other: &PyWittClassG) -> PyResult<PyWittClassG> {
        self.inner
            .try_mul(&other.inner)
            .map(|inner| PyWittClassG { inner })
            .map_err(PyValueError::new_err)
    }
    fn __mul__(&self, other: &PyWittClassG) -> PyResult<PyWittClassG> {
        self.mul(other)
    }
    fn __eq__(&self, other: &PyWittClassG) -> bool {
        self.inner == other.inner
    }
    fn __repr__(&self) -> String {
        match self.inner {
            WittClassG::Char0 { signature } => format!("WittClassG::Char0(signature={signature})"),
            WittClassG::OddChar {
                field_order,
                kappa,
                e0,
                sclass,
            } => {
                format!(
                    "WittClassG::OddChar(field_order={field_order}, kappa={kappa}, e0={e0}, sclass={sclass})"
                )
            }
            WittClassG::Char2 { arf } => format!("WittClassG::Char2(arf={arf})"),
        }
    }
}

/// Classify a diagonal odd-characteristic form `q` over `F_{p^degree}` (dimension
/// + discriminant + Hasse). Supported fields: F_3, F_5, F_7, F_11, F_13, F_9,
/// F_25, F_27.
#[pyfunction]
#[pyo3(signature = (p, q, degree=1))]
fn classify_oddchar(p: u128, q: Vec<i128>, degree: usize) -> PyResult<PyOddCharType> {
    PyFiniteFieldForm::new(p, q, degree)?.classify()
}

/// The odd-characteristic Witt class of a diagonal form `q` over `F_{p^degree}`.
#[pyfunction]
#[pyo3(signature = (p, q, degree=1))]
fn oddchar_witt(p: u128, q: Vec<i128>, degree: usize) -> PyResult<PyWittClassG> {
    PyFiniteFieldForm::new(p, q, degree)?.witt_class()
}

/// Is `x` a square in `F_{p^degree}`?
#[pyfunction]
#[pyo3(signature = (p, x, degree=1))]
fn is_square_mod(p: u128, x: i128, degree: usize) -> PyResult<bool> {
    PyFiniteFieldForm::new(p, Vec::new(), degree)?.is_square(x)
}

/// The Hasse–Witt invariant of a diagonal form `q` over `F_p` (always +1 over a
/// finite field).
#[pyfunction]
#[pyo3(signature = (p, q, degree=1))]
fn hasse_invariant(p: u128, q: Vec<i128>, degree: usize) -> PyResult<i8> {
    PyFiniteFieldForm::new(p, q, degree)?.hasse_invariant()
}

/// Witt decomposition of a diagonal odd-char form `q` over `F_{p^degree}`: returns
/// `(witt_index, anisotropic_dim, anisotropic_disc_is_square, radical_dim)`.
#[pyfunction]
#[pyo3(signature = (p, q, degree=1))]
fn witt_decompose_oddchar(
    p: u128,
    q: Vec<i128>,
    degree: usize,
) -> PyResult<(usize, usize, bool, usize)> {
    PyFiniteFieldForm::new(p, q, degree)?.witt_decompose()
}

/// Are two diagonal odd-char forms over `F_{p^degree}` isometric? `(dim,
/// discriminant)` is a complete invariant.
#[pyfunction]
#[pyo3(signature = (p, q1, q2, degree=1))]
fn is_isometric_oddchar(p: u128, q1: Vec<i128>, q2: Vec<i128>, degree: usize) -> PyResult<bool> {
    let f1 = PyFiniteFieldForm::new(p, q1, degree)?;
    let f2 = PyFiniteFieldForm::new(p, q2, degree)?;
    f1.is_isometric(&f2)
}

/// Witt decomposition of a surreal form on the exact-square real-table subdomain:
/// returns `(witt_index, anisotropic_pos, anisotropic_neg, radical_dim)`.
#[pyfunction]
fn witt_decompose_real(alg: &SurrealAlgebra) -> PyResult<(usize, usize, usize, usize)> {
    let d = crate::forms::witt_decompose_real(&alg.inner.metric).ok_or_else(|| {
        PyValueError::new_err("metric is outside the exact-square real-table subdomain")
    })?;
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

/// The non-Archimedean Springer decomposition of a surreal form: diagonalizes
/// first, then reads the ω-adic valuation filtration into residue ℝ-signatures.
#[pyfunction]
fn springer_decompose(alg: &SurrealAlgebra) -> PyResult<PySpringerDecomp> {
    let d = crate::forms::springer_decompose(&alg.inner.metric).ok_or_else(|| {
        PyValueError::new_err("Springer decomposition could not diagonalize this metric")
    })?;
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
    PyFiniteFieldForm::new(p, q, 1)?.e_staircase()
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
fn hilbert_symbol_qp(a: i128, b: i128, p: u128) -> PyResult<i8> {
    crate::forms::try_hilbert_symbol_qp(a, b, p).ok_or_else(|| {
        PyValueError::new_err(
            "Hilbert symbol needs prime p <= i128::MAX, nonzero arguments, and bounded square classes",
        )
    })
}

/// The Hilbert symbol `(a, b)_∞` over `ℝ` (`−1` iff both are negative).
#[pyfunction]
fn hilbert_symbol_real(a: i128, b: i128) -> i8 {
    crate::forms::hilbert_symbol_real(a, b)
}

/// Is the integer `n` a square in `Q_p`?
#[pyfunction]
fn is_square_qp(n: i128, p: u128) -> PyResult<bool> {
    crate::forms::try_is_square_qp(n, p)
        .ok_or_else(|| PyValueError::new_err("Q_p square test needs prime p <= i128::MAX"))
}

/// Is the diagonal rational/integer form `⟨a₁,…,aₙ⟩` isotropic over `Q`? By the
/// **Hasse–Minkowski** principle (isotropic over `ℝ` and every `Q_p`). E.g.
/// `⟨1,1,1⟩` is anisotropic, `⟨1,1,-1⟩` isotropic, `⟨1,1,-3⟩` anisotropic.
#[pyfunction]
fn is_isotropic_q(entries: Vec<i128>) -> PyResult<bool> {
    crate::forms::try_is_isotropic_q(&entries).ok_or_else(|| {
        PyValueError::new_err("rational isotropy overflowed bounded i128 arithmetic")
    })
}

/// The Hilbert-symbol product `∏_v (a, b)_v` over all places of `ℚ`, for `a, b ∈
/// ℚ^*` passed as `(num, den)` pairs. Equal to `+1` for all `a, b` — Hilbert
/// reciprocity, the multiplicative analogue of the adelic product formula.
#[pyfunction]
fn hilbert_product(a: (i128, i128), b: (i128, i128)) -> PyResult<i8> {
    let a = Rational::try_new(a.0, a.1).ok_or_else(|| {
        PyValueError::new_err("first rational has zero denominator or overflowed bounded i128")
    })?;
    let b = Rational::try_new(b.0, b.1).ok_or_else(|| {
        PyValueError::new_err("second rational has zero denominator or overflowed bounded i128")
    })?;
    Ok(crate::forms::hilbert_product(&a, &b))
}

/// The per-place isotropy breakdown of a `ℚ`-form (rank ≥ 3): isotropy at `ℝ` and
/// at each relevant prime. `is_global()` (isotropic everywhere) equals
/// `is_isotropic_q` (Hasse–Minkowski).
#[pyclass(name = "AdelicIsotropy", module = "pleroma")]
struct PyAdelicIsotropy {
    inner: crate::forms::AdelicIsotropy,
}

#[pymethods]
impl PyAdelicIsotropy {
    /// Isotropy over the Archimedean completion `ℝ`.
    #[getter]
    fn real(&self) -> bool {
        self.inner.real
    }
    /// Isotropy over `Q_p` at each relevant prime, as a `{p: bool}` dict.
    #[getter]
    fn local(&self) -> std::collections::BTreeMap<u128, bool> {
        self.inner.local.clone()
    }
    /// Isotropic over `ℚ` iff isotropic at every place (the local–global principle).
    fn is_global(&self) -> bool {
        self.inner.is_global()
    }
    fn __repr__(&self) -> String {
        format!(
            "AdelicIsotropy(real={}, local={:?}, is_global={})",
            self.inner.real,
            self.inner.local,
            self.inner.is_global()
        )
    }
}

/// The adelic Hasse–Minkowski decomposition of a diagonal integer form of **rank
/// ≥ 3**: isotropy at `ℝ` and each relevant prime. Errors on rank ≤ 2 (there
/// isotropy is a global-square condition — use `is_isotropic_q`).
#[pyfunction]
fn isotropy_over_adeles(entries: Vec<i128>) -> PyResult<PyAdelicIsotropy> {
    if entries.len() < 3 {
        return Err(PyValueError::new_err(
            "adelic isotropy decomposition needs rank ≥ 3 (use is_isotropic_q for rank ≤ 2)",
        ));
    }
    Ok(PyAdelicIsotropy {
        inner: crate::forms::isotropy_over_adeles(&entries),
    })
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
    fn add(&self, other: &PyBrauerWallClass) -> PyResult<PyBrauerWallClass> {
        self.inner
            .try_add(&other.inner)
            .map(|inner| PyBrauerWallClass { inner })
            .map_err(PyValueError::new_err)
    }
    fn __add__(&self, other: &PyBrauerWallClass) -> PyResult<PyBrauerWallClass> {
        self.add(other)
    }
    fn __eq__(&self, other: &PyBrauerWallClass) -> bool {
        self.inner == other.inner
    }
    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}

/// The Brauer–Wall class of a surreal Clifford algebra on the exact-square
/// real-table subdomain: the Bott index `s = (q − p) mod 8`.
#[pyfunction]
fn bw_class_real(alg: &SurrealAlgebra) -> PyResult<PyBrauerWallClass> {
    crate::forms::bw_class_real(&alg.inner.metric)
        .map(|c| PyBrauerWallClass { inner: c })
        .ok_or_else(|| {
            PyValueError::new_err("metric is outside the exact-square real-table subdomain")
        })
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
    PyFiniteFieldForm::new(p, q, 1)?.bw_class()
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyArfResult>()?;
    m.add_class::<PyCliffordType>()?;
    m.add_class::<PyWittClass>()?;
    m.add_class::<PyOddCharType>()?;
    m.add_class::<PyFiniteFieldForm>()?;
    m.add_class::<PyWittClassG>()?;
    m.add_class::<PySpringerDecomp>()?;
    m.add_class::<PyBrauerWallClass>()?;
    m.add_class::<PyAdelicIsotropy>()?;
    m.add_function(wrap_pyfunction!(arf_invariant, m)?)?;
    m.add_function(wrap_pyfunction!(classify_surreal, m)?)?;
    m.add_function(wrap_pyfunction!(classify_surcomplex, m)?)?;
    m.add_function(wrap_pyfunction!(classify_real, m)?)?;
    m.add_function(wrap_pyfunction!(classify_complex, m)?)?;
    m.add_function(wrap_pyfunction!(hilbert_product, m)?)?;
    m.add_function(wrap_pyfunction!(isotropy_over_adeles, m)?)?;
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
