//! Python bindings for the form classifiers across the characteristic
//! trichotomy: Arf (char 2), the char-0 Clifford type, the Witt classes,
//! Dickson, the odd-characteristic classifier, and the Springer decomposition.
//! These consume the `pub(crate)` algebra types stamped by [`super::engine`].

use super::engine::{NimberAlgebra, NimberMV, SurcomplexAlgebra, SurrealAlgebra};
use super::scalars::{parse_surcomplex, parse_surreal, wrap_surreal, PySurreal};
use crate::clifford::{CliffordAlgebra, Metric};
use crate::forms::{
    FiniteOddField, HermitianForm, IntegralForm, SymplecticForm, WittClass, WittClassG,
};
use crate::scalar::{Fp, Fpn, Nimber, Rational, Surcomplex, Surreal};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::sync::Arc;

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

#[pyclass(name = "QuadricFit", module = "pleroma")]
struct PyQuadricFit {
    inner: crate::forms::QuadricFit,
}

#[pymethods]
impl PyQuadricFit {
    #[getter]
    fn constant(&self) -> bool {
        self.inner.constant
    }
    #[getter]
    fn diagonal(&self) -> Vec<bool> {
        self.inner.qd.clone()
    }
    #[getter]
    fn polar_rows(&self) -> Vec<u128> {
        self.inner.bmat.clone()
    }
    fn arf(&self) -> PyArfResult {
        PyArfResult {
            inner: self.inner.arf.clone(),
        }
    }
    fn is_genuinely_quadratic(&self) -> bool {
        self.inner.is_genuinely_quadratic()
    }
    fn __repr__(&self) -> String {
        format!(
            "QuadricFit(constant={}, diagonal={:?}, polar_rows={:?}, arf={:?})",
            self.inner.constant, self.inner.qd, self.inner.bmat, self.inner.arf
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

/// Fit an F₂ quadratic form to a subset of `F_2^k`, returning the recovered
/// coefficients and Arf data if the set is a quadric.
#[pyfunction]
fn fit_f2_quadratic(set: Vec<u128>, k: usize) -> PyResult<Option<PyQuadricFit>> {
    const MAX_ANF_DIM: usize = 20;
    if k > MAX_ANF_DIM {
        return Err(PyValueError::new_err(format!(
            "fit_f2_quadratic is exponential in k; max supported k is {MAX_ANF_DIM}"
        )));
    }
    let domain_mask = if k == 0 { 0 } else { (1u128 << k) - 1 };
    if set.iter().any(|&v| v & !domain_mask != 0) {
        return Err(PyValueError::new_err(format!(
            "point outside F_2^{k} in fit_f2_quadratic input"
        )));
    }
    Ok(crate::forms::fit_f2_quadratic(&set, k).map(|inner| PyQuadricFit { inner }))
}

fn validate_gold_args(m: usize) -> PyResult<()> {
    if !m.is_power_of_two() || m > 128 {
        return Err(PyValueError::new_err(
            "Gold form needs m a positive power of two <= 128",
        ));
    }
    Ok(())
}

/// The Arf data of the Gold form `Q_a(x)=Tr(x^(1+2^a))` on the nim subfield
/// `F_{2^m}`.
#[pyfunction]
fn gold_form_arf(m: usize, a: usize) -> PyResult<PyArfResult> {
    validate_gold_args(m)?;
    let metric = crate::forms::gold_form(m, a);
    crate::forms::arf_invariant(&metric)
        .map(|inner| PyArfResult { inner })
        .ok_or_else(|| PyValueError::new_err("Gold form unexpectedly failed Arf classification"))
}

/// The Gold form as a `NimberAlgebra`, so Python can inspect the underlying
/// Clifford product as well as its Arf invariant.
#[pyfunction]
fn gold_form_algebra(m: usize, a: usize) -> PyResult<NimberAlgebra> {
    validate_gold_args(m)?;
    let metric = crate::forms::gold_form(m, a);
    Ok(NimberAlgebra {
        inner: Arc::new(CliffordAlgebra::new(metric.q.len(), metric)),
    })
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
// Alternating and Hermitian forms (the "form + involution" siblings)
// ---------------------------------------------------------------------------

#[pyclass(name = "SymplecticClass", module = "pleroma")]
struct PySymplecticClass {
    inner: crate::forms::SymplecticClass,
}

#[pymethods]
impl PySymplecticClass {
    #[getter]
    fn rank(&self) -> usize {
        self.inner.rank
    }
    #[getter]
    fn radical_dim(&self) -> usize {
        self.inner.radical_dim
    }
    fn planes(&self) -> usize {
        self.inner.planes()
    }
    fn __repr__(&self) -> String {
        format!(
            "SymplecticClass(rank={}, radical_dim={}, planes={})",
            self.inner.rank,
            self.inner.radical_dim,
            self.inner.planes()
        )
    }
}

fn rational_gram(gram: Vec<Vec<i128>>) -> Vec<Vec<Rational>> {
    gram.into_iter()
        .map(|row| row.into_iter().map(Rational::int).collect())
        .collect()
}

#[pyclass(name = "SymplecticForm", module = "pleroma", from_py_object)]
#[derive(Clone)]
struct PySymplecticForm {
    inner: SymplecticForm<Rational>,
}

#[pymethods]
impl PySymplecticForm {
    #[new]
    fn new(gram: Vec<Vec<i128>>) -> PyResult<Self> {
        SymplecticForm::from_gram(rational_gram(gram))
            .map(|inner| PySymplecticForm { inner })
            .ok_or_else(|| PyValueError::new_err("Gram matrix must be square and alternating"))
    }
    #[staticmethod]
    fn hyperbolic(r: usize) -> PySymplecticForm {
        PySymplecticForm {
            inner: SymplecticForm::hyperbolic(r),
        }
    }
    #[getter]
    fn dim(&self) -> usize {
        self.inner.dim()
    }
    fn direct_sum(&self, other: &PySymplecticForm) -> PySymplecticForm {
        PySymplecticForm {
            inner: self.inner.direct_sum(&other.inner),
        }
    }
    fn classify(&self) -> PySymplecticClass {
        PySymplecticClass {
            inner: self.inner.classify(),
        }
    }
    fn __repr__(&self) -> String {
        format!("SymplecticForm(dim={})", self.inner.dim())
    }
}

/// Classify an integer/rational alternating Gram matrix: complete invariant
/// `(rank, radical_dim)`.
#[pyfunction]
fn classify_symplectic(gram: Vec<Vec<i128>>) -> PyResult<PySymplecticClass> {
    crate::forms::classify_symplectic(rational_gram(gram))
        .map(|inner| PySymplecticClass { inner })
        .ok_or_else(|| PyValueError::new_err("Gram matrix must be square and alternating"))
}

/// The same alternating-form classifier over the nimber backend, where
/// alternating means symmetric with zero diagonal because `-1 = 1`.
#[pyfunction]
fn classify_symplectic_nimber(gram: Vec<Vec<u128>>) -> PyResult<PySymplecticClass> {
    let gram: Vec<Vec<Nimber>> = gram
        .into_iter()
        .map(|row| row.into_iter().map(Nimber).collect())
        .collect();
    crate::forms::classify_symplectic(gram)
        .map(|inner| PySymplecticClass { inner })
        .ok_or_else(|| PyValueError::new_err("Nimber Gram matrix must be square and alternating"))
}

#[pyclass(name = "HermitianSignature", module = "pleroma")]
struct PyHermitianSignature {
    inner: crate::forms::HermitianSignature,
}

#[pymethods]
impl PyHermitianSignature {
    #[getter]
    fn pos(&self) -> usize {
        self.inner.pos
    }
    #[getter]
    fn neg(&self) -> usize {
        self.inner.neg
    }
    #[getter]
    fn radical(&self) -> usize {
        self.inner.radical
    }
    fn as_tuple(&self) -> (usize, usize, usize) {
        (self.inner.pos, self.inner.neg, self.inner.radical)
    }
    fn __repr__(&self) -> String {
        format!(
            "HermitianSignature(pos={}, neg={}, radical={})",
            self.inner.pos, self.inner.neg, self.inner.radical
        )
    }
}

#[pyclass(name = "HermitianForm", module = "pleroma", from_py_object)]
#[derive(Clone)]
struct PyHermitianForm {
    inner: HermitianForm<Surreal>,
}

fn parse_surcomplex_gram(
    gram: Vec<Vec<Bound<'_, PyAny>>>,
) -> PyResult<Vec<Vec<Surcomplex<Surreal>>>> {
    let n = gram.len();
    let mut out = Vec::with_capacity(n);
    for row in &gram {
        if row.len() != n {
            return Err(PyValueError::new_err("Gram matrix must be square"));
        }
        let mut r = Vec::with_capacity(n);
        for x in row {
            r.push(parse_surcomplex(x)?);
        }
        out.push(r);
    }
    Ok(out)
}

#[pymethods]
impl PyHermitianForm {
    #[new]
    fn new(gram: Vec<Vec<Bound<'_, PyAny>>>) -> PyResult<Self> {
        HermitianForm::from_gram(parse_surcomplex_gram(gram)?)
            .map(|inner| PyHermitianForm { inner })
            .ok_or_else(|| PyValueError::new_err("Gram matrix must be Hermitian"))
    }
    #[staticmethod]
    fn from_skew(gram: Vec<Vec<Bound<'_, PyAny>>>) -> PyResult<PyHermitianForm> {
        HermitianForm::from_skew(parse_surcomplex_gram(gram)?)
            .map(|inner| PyHermitianForm { inner })
            .ok_or_else(|| PyValueError::new_err("Gram matrix must be skew-Hermitian"))
    }
    #[staticmethod]
    fn diagonal(reals: Vec<Bound<'_, PyAny>>) -> PyResult<PyHermitianForm> {
        let mut ds = Vec::with_capacity(reals.len());
        for x in &reals {
            ds.push(parse_surreal(x)?);
        }
        Ok(PyHermitianForm {
            inner: HermitianForm::diagonal(ds),
        })
    }
    #[getter]
    fn dim(&self) -> usize {
        self.inner.dim()
    }
    fn direct_sum(&self, other: &PyHermitianForm) -> PyHermitianForm {
        PyHermitianForm {
            inner: self.inner.direct_sum(&other.inner),
        }
    }
    fn diagonalize(&self) -> Vec<PySurreal> {
        self.inner
            .diagonalize()
            .into_iter()
            .map(wrap_surreal)
            .collect()
    }
    fn signature(&self) -> PyHermitianSignature {
        PyHermitianSignature {
            inner: self.inner.signature(|x| x.sign()),
        }
    }
    fn __repr__(&self) -> String {
        format!("HermitianForm(dim={})", self.inner.dim())
    }
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

fn unsupported_prime_field_err() -> PyErr {
    PyValueError::new_err("supported prime fields: F_2, F_3, F_5, F_7, F_11, F_13")
}

macro_rules! with_prime_field {
    ($p:expr, $body:ident) => {{
        match $p {
            2 => $body::<2>(),
            3 => $body::<3>(),
            5 => $body::<5>(),
            7 => $body::<7>(),
            11 => $body::<11>(),
            13 => $body::<13>(),
            _ => return Err(unsupported_prime_field_err()),
        }
    }};
}

fn level_for_prime<const P: u128>() -> PyResult<Option<usize>> {
    Ok(crate::forms::level::<P>())
}

fn pythagoras_for_prime<const P: u128>() -> PyResult<Option<usize>> {
    Ok(crate::forms::pythagoras_number::<P>())
}

fn u_invariant_for_prime<const P: u128>() -> PyResult<Option<usize>> {
    Ok(crate::forms::u_invariant::<P>())
}

fn sum_of_squares_for_prime<const P: u128>(x: i128, n: usize) -> PyResult<bool> {
    Ok(crate::forms::is_sum_of_n_squares::<P>(Fp::<P>::new(x), n))
}

/// The level/Stufe of the prime field `F_p`: least `n` with `-1` a sum of `n`
/// squares. Returns `None` only for the char-2/degenerate cases where the Rust
/// invariant deliberately declines; supported dispatch primes are finite.
#[pyfunction]
fn finite_field_level(p: u128) -> PyResult<Option<usize>> {
    with_prime_field!(p, level_for_prime)
}

/// The Pythagoras number of the prime field `F_p`: least `n` such that every sum
/// of squares is already a sum of `n` squares.
#[pyfunction]
fn finite_field_pythagoras_number(p: u128) -> PyResult<Option<usize>> {
    with_prime_field!(p, pythagoras_for_prime)
}

/// The u-invariant of the prime field `F_p`: largest dimension of an anisotropic
/// quadratic form. In characteristic 2 this returns `None` because the diagonal
/// odd-characteristic model is not the right form theory.
#[pyfunction]
fn finite_field_u_invariant(p: u128) -> PyResult<Option<usize>> {
    with_prime_field!(p, u_invariant_for_prime)
}

/// Is `x` a sum of exactly `n` squares in the prime field `F_p`?
#[pyfunction]
fn is_sum_of_n_squares(p: u128, x: i128, n: usize) -> PyResult<bool> {
    match p {
        2 => sum_of_squares_for_prime::<2>(x, n),
        3 => sum_of_squares_for_prime::<3>(x, n),
        5 => sum_of_squares_for_prime::<5>(x, n),
        7 => sum_of_squares_for_prime::<7>(x, n),
        11 => sum_of_squares_for_prime::<11>(x, n),
        13 => sum_of_squares_for_prime::<13>(x, n),
        _ => Err(unsupported_prime_field_err()),
    }
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
// Integral lattices, genus, ADE catalogue, mass / Leech
// ---------------------------------------------------------------------------

#[pyclass(name = "ScaleSymbol", module = "pleroma", skip_from_py_object)]
#[derive(Clone)]
struct PyScaleSymbol {
    inner: crate::forms::ScaleSymbol,
}

#[pymethods]
impl PyScaleSymbol {
    #[getter]
    fn scale(&self) -> u32 {
        self.inner.scale
    }
    #[getter]
    fn dim(&self) -> usize {
        self.inner.dim
    }
    #[getter]
    fn sign(&self) -> i8 {
        self.inner.sign
    }
    #[getter]
    fn det_mod8(&self) -> i64 {
        self.inner.det_mod8
    }
    #[getter]
    fn type_ii(&self) -> bool {
        self.inner.type_ii
    }
    #[getter]
    fn oddity(&self) -> i64 {
        self.inner.oddity
    }
    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}

#[pyclass(name = "Genus", module = "pleroma", skip_from_py_object)]
#[derive(Clone)]
struct PyGenus {
    inner: crate::forms::Genus,
}

#[pymethods]
impl PyGenus {
    #[getter]
    fn dim(&self) -> usize {
        self.inner.dim
    }
    #[getter]
    fn signature(&self) -> (usize, usize) {
        self.inner.signature
    }
    #[getter]
    fn det(&self) -> i128 {
        self.inner.det
    }
    fn primes(&self) -> Vec<u128> {
        self.inner.primes()
    }
    fn symbol_at(&self, p: u128) -> Vec<PyScaleSymbol> {
        self.inner
            .symbol_at(p)
            .iter()
            .cloned()
            .map(|inner| PyScaleSymbol { inner })
            .collect()
    }
    fn __repr__(&self) -> String {
        format!(
            "Genus(dim={}, signature={:?}, det={}, primes={:?})",
            self.inner.dim,
            self.inner.signature,
            self.inner.det,
            self.inner.primes()
        )
    }
}

#[pyclass(name = "IntegralForm", module = "pleroma", from_py_object)]
#[derive(Clone)]
struct PyIntegralForm {
    inner: IntegralForm,
}

fn check_lattice_vec(l: &IntegralForm, v: &[i128], name: &str) -> PyResult<()> {
    if v.len() != l.dim() {
        Err(PyValueError::new_err(format!(
            "{name} has length {}, expected {}",
            v.len(),
            l.dim()
        )))
    } else {
        Ok(())
    }
}

#[pymethods]
impl PyIntegralForm {
    #[new]
    fn new(gram: Vec<Vec<i128>>) -> PyResult<Self> {
        IntegralForm::new(gram)
            .map(|inner| PyIntegralForm { inner })
            .ok_or_else(|| PyValueError::new_err("Gram matrix must be square and symmetric"))
    }
    #[staticmethod]
    fn diagonal(diag: Vec<i128>) -> PyIntegralForm {
        PyIntegralForm {
            inner: IntegralForm::diagonal(&diag),
        }
    }
    #[staticmethod]
    fn a(n: usize) -> PyResult<PyIntegralForm> {
        if n < 1 {
            return Err(PyValueError::new_err("A_n requires n >= 1"));
        }
        Ok(PyIntegralForm {
            inner: crate::forms::a_n(n),
        })
    }
    #[staticmethod]
    fn d(n: usize) -> PyResult<PyIntegralForm> {
        if n < 2 {
            return Err(PyValueError::new_err("D_n requires n >= 2"));
        }
        Ok(PyIntegralForm {
            inner: crate::forms::d_n(n),
        })
    }
    #[staticmethod]
    fn e6() -> PyIntegralForm {
        PyIntegralForm {
            inner: crate::forms::e_6(),
        }
    }
    #[staticmethod]
    fn e7() -> PyIntegralForm {
        PyIntegralForm {
            inner: crate::forms::e_7(),
        }
    }
    #[staticmethod]
    fn e8() -> PyIntegralForm {
        PyIntegralForm {
            inner: crate::forms::e_8(),
        }
    }
    #[staticmethod]
    fn leech() -> PyIntegralForm {
        PyIntegralForm {
            inner: crate::forms::leech(),
        }
    }
    #[getter]
    fn dim(&self) -> usize {
        self.inner.dim()
    }
    #[getter]
    fn gram(&self) -> Vec<Vec<i128>> {
        self.inner.gram().to_vec()
    }
    fn inner(&self, x: Vec<i128>, y: Vec<i128>) -> PyResult<i128> {
        check_lattice_vec(&self.inner, &x, "x")?;
        check_lattice_vec(&self.inner, &y, "y")?;
        Ok(self.inner.inner(&x, &y))
    }
    fn norm(&self, x: Vec<i128>) -> PyResult<i128> {
        check_lattice_vec(&self.inner, &x, "x")?;
        Ok(self.inner.norm(&x))
    }
    fn determinant(&self) -> i128 {
        self.inner.determinant()
    }
    fn is_unimodular(&self) -> bool {
        self.inner.is_unimodular()
    }
    fn is_even(&self) -> bool {
        self.inner.is_even()
    }
    fn is_positive_definite(&self) -> bool {
        self.inner.is_positive_definite()
    }
    fn invariant_factors(&self) -> Vec<i128> {
        self.inner.invariant_factors()
    }
    fn level(&self) -> Option<i128> {
        self.inner.level()
    }
    fn direct_sum(&self, other: &PyIntegralForm) -> PyIntegralForm {
        PyIntegralForm {
            inner: self.inner.direct_sum(&other.inner),
        }
    }
    fn short_vectors(&self, bound: i128) -> Option<Vec<Vec<i128>>> {
        self.inner.short_vectors(bound)
    }
    fn minimum(&self) -> Option<i128> {
        self.inner.minimum()
    }
    fn minimal_vectors(&self) -> Option<Vec<Vec<i128>>> {
        self.inner.minimal_vectors()
    }
    fn kissing_number(&self) -> Option<usize> {
        self.inner.kissing_number()
    }
    fn automorphism_group_order(&self) -> Option<u128> {
        self.inner.automorphism_group_order()
    }
    fn automorphism_group_order_bounded(&self, node_budget: u64) -> Option<u128> {
        self.inner.automorphism_group_order_bounded(node_budget)
    }
    fn coxeter_number(&self) -> Option<i128> {
        crate::forms::coxeter_number(&self.inner)
    }
    fn is_root_lattice(&self) -> bool {
        crate::forms::is_root_lattice(&self.inner)
    }
    fn genus(&self) -> Option<PyGenus> {
        crate::forms::Genus::of(&self.inner).map(|inner| PyGenus { inner })
    }
    fn same_genus(&self, other: &PyIntegralForm) -> bool {
        crate::forms::are_in_same_genus(&self.inner, &other.inner)
    }
    fn __repr__(&self) -> String {
        format!("IntegralForm(gram={:?})", self.inner.gram())
    }
}

#[pyfunction]
fn root_lattice_a(n: usize) -> PyResult<PyIntegralForm> {
    PyIntegralForm::a(n)
}

#[pyfunction]
fn root_lattice_d(n: usize) -> PyResult<PyIntegralForm> {
    PyIntegralForm::d(n)
}

#[pyfunction]
fn root_lattice_e6() -> PyIntegralForm {
    PyIntegralForm::e6()
}

#[pyfunction]
fn root_lattice_e7() -> PyIntegralForm {
    PyIntegralForm::e7()
}

#[pyfunction]
fn root_lattice_e8() -> PyIntegralForm {
    PyIntegralForm::e8()
}

#[pyfunction]
fn leech_lattice() -> PyIntegralForm {
    PyIntegralForm::leech()
}

#[pyfunction]
fn are_in_same_genus(a: &PyIntegralForm, b: &PyIntegralForm) -> bool {
    crate::forms::are_in_same_genus(&a.inner, &b.inner)
}

#[pyfunction]
fn mass_even_unimodular(n: u32) -> Option<(i128, i128)> {
    crate::forms::mass_even_unimodular(n)
}

#[pyfunction]
fn leech_aut_order() -> u128 {
    crate::forms::LEECH_AUT_ORDER
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
        match self.inner {
            crate::forms::BrauerWallClass::Real(s) => {
                format!("BrauerWallClass::Real({s})")
            }
            crate::forms::BrauerWallClass::Complex(p) => {
                format!("BrauerWallClass::Complex({p})")
            }
            crate::forms::BrauerWallClass::OddChar {
                field_order,
                kappa,
                e0,
                sclass,
            } => {
                format!(
                    "BrauerWallClass::OddChar(field_order={field_order}, kappa={kappa}, e0={e0}, sclass={sclass})"
                )
            }
            crate::forms::BrauerWallClass::Char2 { arf } => {
                format!("BrauerWallClass::Char2(arf={arf})")
            }
        }
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

/// The Brauer-Wall class of a nonsingular nimber Clifford algebra in
/// `BW(F_{2^m}) ≅ W_q(F_{2^m}) ≅ Z/2` (the Arf/Witt class).
#[pyfunction]
fn bw_class_nimber(alg: &NimberAlgebra) -> PyResult<PyBrauerWallClass> {
    crate::forms::bw_class_nimber(&alg.inner.metric)
        .map(|c| PyBrauerWallClass { inner: c })
        .ok_or_else(|| {
            PyValueError::new_err(
                "char-2 Brauer-Wall class needs a nonsingular non-general-bilinear nimber metric",
            )
        })
}

/// The Brauer–Wall class of a diagonal odd-char form `q` over `F_p` (the order-4
/// graded part, `BW(F_q) ≅ W(F_q)`). Supported primes: 3, 5, 7, 11, 13.
#[pyfunction]
fn bw_class_oddchar(p: u128, q: Vec<i128>) -> PyResult<PyBrauerWallClass> {
    PyFiniteFieldForm::new(p, q, 1)?.bw_class()
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyArfResult>()?;
    m.add_class::<PyQuadricFit>()?;
    m.add_class::<PyCliffordType>()?;
    m.add_class::<PyWittClass>()?;
    m.add_class::<PyOddCharType>()?;
    m.add_class::<PyFiniteFieldForm>()?;
    m.add_class::<PyWittClassG>()?;
    m.add_class::<PySymplecticClass>()?;
    m.add_class::<PySymplecticForm>()?;
    m.add_class::<PyHermitianSignature>()?;
    m.add_class::<PyHermitianForm>()?;
    m.add_class::<PySpringerDecomp>()?;
    m.add_class::<PyBrauerWallClass>()?;
    m.add_class::<PyAdelicIsotropy>()?;
    m.add_class::<PyIntegralForm>()?;
    m.add_class::<PyScaleSymbol>()?;
    m.add_class::<PyGenus>()?;
    m.add_function(wrap_pyfunction!(arf_invariant, m)?)?;
    m.add_function(wrap_pyfunction!(fit_f2_quadratic, m)?)?;
    m.add_function(wrap_pyfunction!(gold_form_arf, m)?)?;
    m.add_function(wrap_pyfunction!(gold_form_algebra, m)?)?;
    m.add_function(wrap_pyfunction!(classify_surreal, m)?)?;
    m.add_function(wrap_pyfunction!(classify_surcomplex, m)?)?;
    m.add_function(wrap_pyfunction!(classify_real, m)?)?;
    m.add_function(wrap_pyfunction!(classify_complex, m)?)?;
    m.add_function(wrap_pyfunction!(hilbert_product, m)?)?;
    m.add_function(wrap_pyfunction!(isotropy_over_adeles, m)?)?;
    m.add_function(wrap_pyfunction!(root_lattice_a, m)?)?;
    m.add_function(wrap_pyfunction!(root_lattice_d, m)?)?;
    m.add_function(wrap_pyfunction!(root_lattice_e6, m)?)?;
    m.add_function(wrap_pyfunction!(root_lattice_e7, m)?)?;
    m.add_function(wrap_pyfunction!(root_lattice_e8, m)?)?;
    m.add_function(wrap_pyfunction!(leech_lattice, m)?)?;
    m.add_function(wrap_pyfunction!(are_in_same_genus, m)?)?;
    m.add_function(wrap_pyfunction!(mass_even_unimodular, m)?)?;
    m.add_function(wrap_pyfunction!(leech_aut_order, m)?)?;
    m.add_function(wrap_pyfunction!(witt_class, m)?)?;
    m.add_function(wrap_pyfunction!(dickson_matrix, m)?)?;
    m.add_function(wrap_pyfunction!(dickson_of_versor, m)?)?;
    m.add_function(wrap_pyfunction!(classify_symplectic, m)?)?;
    m.add_function(wrap_pyfunction!(classify_symplectic_nimber, m)?)?;
    m.add_function(wrap_pyfunction!(classify_oddchar, m)?)?;
    m.add_function(wrap_pyfunction!(oddchar_witt, m)?)?;
    m.add_function(wrap_pyfunction!(witt_decompose_oddchar, m)?)?;
    m.add_function(wrap_pyfunction!(witt_decompose_real, m)?)?;
    m.add_function(wrap_pyfunction!(is_isometric_oddchar, m)?)?;
    m.add_function(wrap_pyfunction!(is_square_mod, m)?)?;
    m.add_function(wrap_pyfunction!(finite_field_level, m)?)?;
    m.add_function(wrap_pyfunction!(finite_field_pythagoras_number, m)?)?;
    m.add_function(wrap_pyfunction!(finite_field_u_invariant, m)?)?;
    m.add_function(wrap_pyfunction!(is_sum_of_n_squares, m)?)?;
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
    m.add_function(wrap_pyfunction!(bw_class_nimber, m)?)?;
    m.add_function(wrap_pyfunction!(bw_class_oddchar, m)?)?;
    Ok(())
}
