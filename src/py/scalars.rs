//! Python bindings for the scalar worlds: the per-backend scalar types
//! (`Nimber`, `Surreal`, `Surcomplex`, `Integer`, `Omnific`, `Ordinal`), their
//! constructors, and the nim-field operations. `parse_*` / `wrap_*` are
//! `pub(crate)` because the `backend!` macro in [`super::engine`] threads them in
//! as the per-backend parse/wrap hooks.

use crate::scalar::{Integer, Nimber, Omnific, Ordinal, Rational, Scalar, Surcomplex, Surreal};
use pyo3::basic::CompareOp;
use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::IntoPyObjectExt;

// ---------------------------------------------------------------------------
// Scalar pyclasses + parsers
// ---------------------------------------------------------------------------

#[pyclass(name = "Nimber", module = "pleroma", from_py_object)]
#[derive(Clone)]
pub(crate) struct PyNimber {
    inner: Nimber,
}

#[pymethods]
impl PyNimber {
    #[new]
    fn new(value: u128) -> Self {
        PyNimber {
            inner: Nimber(value),
        }
    }
    #[getter]
    fn value(&self) -> u128 {
        self.inner.0
    }
    fn __add__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyNimber> {
        Ok(PyNimber {
            inner: self.inner.add(&parse_nimber(other)?),
        })
    }
    fn __radd__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyNimber> {
        self.__add__(other)
    }
    fn __mul__(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        // defer to the other operand (e.g. a multivector's __rmul__) if it isn't a scalar
        match parse_nimber(other) {
            Ok(o) => PyNimber {
                inner: self.inner.mul(&o),
            }
            .into_py_any(py),
            Err(_) => Ok(py.NotImplemented()),
        }
    }
    fn __rmul__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyNimber> {
        Ok(PyNimber {
            inner: self.inner.mul(&parse_nimber(other)?),
        })
    }
    fn inv(&self) -> PyResult<PyNimber> {
        self.inner
            .inv()
            .map(|n| PyNimber { inner: n })
            .ok_or_else(|| PyValueError::new_err("*0 has no inverse"))
    }
    fn __truediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyNimber> {
        let o = parse_nimber(other)?;
        let oi = o
            .inv()
            .ok_or_else(|| PyValueError::new_err("division by *0"))?;
        Ok(PyNimber {
            inner: self.inner.mul(&oi),
        })
    }
    fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
        matches!(parse_nimber(other), Ok(n) if n == self.inner)
    }
    fn __hash__(&self) -> isize {
        self.inner.0 as isize
    }
    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}

pub(crate) fn parse_nimber(obj: &Bound<'_, PyAny>) -> PyResult<Nimber> {
    if let Ok(n) = obj.cast::<PyNimber>() {
        return Ok(n.borrow().inner);
    }
    if let Ok(v) = obj.extract::<u128>() {
        return Ok(Nimber(v));
    }
    Err(PyTypeError::new_err("expected Nimber or non-negative int"))
}

#[pyclass(name = "Surreal", module = "pleroma", from_py_object)]
#[derive(Clone)]
pub(crate) struct PySurreal {
    inner: Surreal,
}

#[pymethods]
impl PySurreal {
    fn __add__(&self, other: &Bound<'_, PyAny>) -> PyResult<PySurreal> {
        Ok(PySurreal {
            inner: self.inner.add(&parse_surreal(other)?),
        })
    }
    fn __radd__(&self, other: &Bound<'_, PyAny>) -> PyResult<PySurreal> {
        self.__add__(other)
    }
    fn __sub__(&self, other: &Bound<'_, PyAny>) -> PyResult<PySurreal> {
        Ok(PySurreal {
            inner: self.inner.sub(&parse_surreal(other)?),
        })
    }
    fn __rsub__(&self, other: &Bound<'_, PyAny>) -> PyResult<PySurreal> {
        Ok(PySurreal {
            inner: parse_surreal(other)?.sub(&self.inner),
        })
    }
    fn __mul__(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        match parse_surreal(other) {
            Ok(o) => PySurreal {
                inner: self.inner.mul(&o),
            }
            .into_py_any(py),
            Err(_) => Ok(py.NotImplemented()),
        }
    }
    fn __rmul__(&self, other: &Bound<'_, PyAny>) -> PyResult<PySurreal> {
        Ok(PySurreal {
            inner: self.inner.mul(&parse_surreal(other)?),
        })
    }
    fn __neg__(&self) -> PySurreal {
        PySurreal {
            inner: self.inner.neg(),
        }
    }
    fn inv(&self) -> PyResult<PySurreal> {
        self.inner
            .inv()
            .map(|s| PySurreal { inner: s })
            .ok_or_else(|| {
                PyValueError::new_err(
                    "only monomials (coeff·ω^e) have a finite-support surreal inverse",
                )
            })
    }
    fn __truediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<PySurreal> {
        let o = parse_surreal(other)?;
        let oi = o
            .inv()
            .ok_or_else(|| PyValueError::new_err("divisor has no finite-support inverse"))?;
        Ok(PySurreal {
            inner: self.inner.mul(&oi),
        })
    }
    fn __pow__(&self, n: u128, _modulo: Option<&Bound<'_, PyAny>>) -> PySurreal {
        let mut acc = Surreal::one();
        for _ in 0..n {
            acc = acc.mul(&self.inner);
        }
        PySurreal { inner: acc }
    }
    fn __richcmp__(&self, other: &Bound<'_, PyAny>, op: CompareOp) -> PyResult<bool> {
        Ok(op.matches(self.inner.cmp(&parse_surreal(other)?)))
    }
    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}

pub(crate) fn parse_surreal(obj: &Bound<'_, PyAny>) -> PyResult<Surreal> {
    if let Ok(s) = obj.cast::<PySurreal>() {
        return Ok(s.borrow().inner.clone());
    }
    if let Ok(v) = obj.extract::<i128>() {
        return Ok(Surreal::from_int(v));
    }
    Err(PyTypeError::new_err("expected Surreal or int"))
}

#[pyclass(name = "Surcomplex", module = "pleroma", from_py_object)]
#[derive(Clone)]
pub(crate) struct PySurcomplex {
    inner: Surcomplex<Surreal>,
}

#[pymethods]
impl PySurcomplex {
    #[new]
    #[pyo3(signature = (re, im=None))]
    fn new(re: &Bound<'_, PyAny>, im: Option<&Bound<'_, PyAny>>) -> PyResult<Self> {
        let r = parse_surreal(re)?;
        let i = match im {
            Some(x) => parse_surreal(x)?,
            None => Surreal::zero(),
        };
        Ok(PySurcomplex {
            inner: Surcomplex::new(r, i),
        })
    }
    #[staticmethod]
    fn i() -> PySurcomplex {
        PySurcomplex {
            inner: Surcomplex::i(),
        }
    }
    #[getter]
    fn re(&self) -> PySurreal {
        PySurreal {
            inner: self.inner.re.clone(),
        }
    }
    #[getter]
    fn im(&self) -> PySurreal {
        PySurreal {
            inner: self.inner.im.clone(),
        }
    }
    fn conj(&self) -> PySurcomplex {
        PySurcomplex {
            inner: self.inner.conj(),
        }
    }
    fn __add__(&self, other: &Bound<'_, PyAny>) -> PyResult<PySurcomplex> {
        Ok(PySurcomplex {
            inner: self.inner.add(&parse_surcomplex(other)?),
        })
    }
    fn __radd__(&self, other: &Bound<'_, PyAny>) -> PyResult<PySurcomplex> {
        self.__add__(other)
    }
    fn __sub__(&self, other: &Bound<'_, PyAny>) -> PyResult<PySurcomplex> {
        Ok(PySurcomplex {
            inner: self.inner.sub(&parse_surcomplex(other)?),
        })
    }
    fn __mul__(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        match parse_surcomplex(other) {
            Ok(o) => PySurcomplex {
                inner: self.inner.mul(&o),
            }
            .into_py_any(py),
            Err(_) => Ok(py.NotImplemented()),
        }
    }
    fn __rmul__(&self, other: &Bound<'_, PyAny>) -> PyResult<PySurcomplex> {
        Ok(PySurcomplex {
            inner: self.inner.mul(&parse_surcomplex(other)?),
        })
    }
    fn __neg__(&self) -> PySurcomplex {
        PySurcomplex {
            inner: self.inner.neg(),
        }
    }
    fn inv(&self) -> PyResult<PySurcomplex> {
        self.inner
            .inv()
            .map(|s| PySurcomplex { inner: s })
            .ok_or_else(|| {
                PyValueError::new_err("inverse needs an invertible norm a²+b² (a monomial surreal)")
            })
    }
    fn __truediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<PySurcomplex> {
        let o = parse_surcomplex(other)?;
        let oi = o
            .inv()
            .ok_or_else(|| PyValueError::new_err("divisor has no representable inverse"))?;
        Ok(PySurcomplex {
            inner: self.inner.mul(&oi),
        })
    }
    fn __pow__(&self, n: u128, _modulo: Option<&Bound<'_, PyAny>>) -> PySurcomplex {
        let mut acc = Surcomplex::<Surreal>::one();
        for _ in 0..n {
            acc = acc.mul(&self.inner);
        }
        PySurcomplex { inner: acc }
    }
    fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
        matches!(parse_surcomplex(other), Ok(x) if x == self.inner)
    }
    fn __repr__(&self) -> String {
        if self.inner.im.is_zero() {
            format!("{:?}", self.inner.re)
        } else {
            format!("{:?} + ({:?})i", self.inner.re, self.inner.im)
        }
    }
}

pub(crate) fn parse_surcomplex(obj: &Bound<'_, PyAny>) -> PyResult<Surcomplex<Surreal>> {
    if let Ok(s) = obj.cast::<PySurcomplex>() {
        return Ok(s.borrow().inner.clone());
    }
    if let Ok(s) = obj.cast::<PySurreal>() {
        return Ok(Surcomplex::new(s.borrow().inner.clone(), Surreal::zero()));
    }
    if let Ok(v) = obj.extract::<i128>() {
        return Ok(Surcomplex::new(Surreal::from_int(v), Surreal::zero()));
    }
    Err(PyTypeError::new_err("expected Surcomplex, Surreal, or int"))
}

#[pyclass(name = "Integer", module = "pleroma", from_py_object)]
#[derive(Clone)]
pub(crate) struct PyInteger {
    inner: Integer,
}

#[pymethods]
impl PyInteger {
    #[new]
    fn new(value: i128) -> Self {
        PyInteger {
            inner: Integer(value),
        }
    }
    #[getter]
    fn value(&self) -> i128 {
        self.inner.0
    }
    fn __add__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyInteger> {
        Ok(PyInteger {
            inner: self.inner.add(&parse_integer(other)?),
        })
    }
    fn __radd__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyInteger> {
        self.__add__(other)
    }
    fn __sub__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyInteger> {
        Ok(PyInteger {
            inner: self.inner.sub(&parse_integer(other)?),
        })
    }
    fn __neg__(&self) -> PyInteger {
        PyInteger {
            inner: self.inner.neg(),
        }
    }
    fn __mul__(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        match parse_integer(other) {
            Ok(o) => PyInteger {
                inner: self.inner.mul(&o),
            }
            .into_py_any(py),
            Err(_) => Ok(py.NotImplemented()),
        }
    }
    fn __rmul__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyInteger> {
        Ok(PyInteger {
            inner: self.inner.mul(&parse_integer(other)?),
        })
    }
    fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
        matches!(parse_integer(other), Ok(n) if n == self.inner)
    }
    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}

pub(crate) fn parse_integer(obj: &Bound<'_, PyAny>) -> PyResult<Integer> {
    if let Ok(n) = obj.cast::<PyInteger>() {
        return Ok(n.borrow().inner);
    }
    if let Ok(v) = obj.extract::<i128>() {
        return Ok(Integer(v));
    }
    Err(PyTypeError::new_err("expected Integer or int"))
}

pub(crate) fn wrap_integer(i: Integer) -> PyInteger {
    PyInteger { inner: i }
}

pub(crate) fn wrap_nimber(n: Nimber) -> PyNimber {
    PyNimber { inner: n }
}
pub(crate) fn wrap_surreal(s: Surreal) -> PySurreal {
    PySurreal { inner: s }
}
pub(crate) fn wrap_surcomplex(s: Surcomplex<Surreal>) -> PySurcomplex {
    PySurcomplex { inner: s }
}

// --- Omnific integers Oz: the surreal integers, a transfinite ring ----------

#[pyclass(name = "Omnific", module = "pleroma", from_py_object)]
#[derive(Clone)]
pub(crate) struct PyOmnific {
    inner: Omnific,
}

#[pymethods]
impl PyOmnific {
    #[new]
    fn new(value: i128) -> Self {
        PyOmnific {
            inner: Omnific::from_int(value),
        }
    }
    /// The underlying surreal value.
    fn surreal(&self) -> PySurreal {
        PySurreal {
            inner: self.inner.inner().clone(),
        }
    }
    fn __add__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyOmnific> {
        Ok(PyOmnific {
            inner: self.inner.add(&parse_omnific(other)?),
        })
    }
    fn __radd__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyOmnific> {
        self.__add__(other)
    }
    fn __sub__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyOmnific> {
        Ok(PyOmnific {
            inner: self.inner.sub(&parse_omnific(other)?),
        })
    }
    fn __neg__(&self) -> PyOmnific {
        PyOmnific {
            inner: self.inner.neg(),
        }
    }
    fn __mul__(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        match parse_omnific(other) {
            Ok(o) => PyOmnific {
                inner: self.inner.mul(&o),
            }
            .into_py_any(py),
            Err(_) => Ok(py.NotImplemented()),
        }
    }
    fn __rmul__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyOmnific> {
        Ok(PyOmnific {
            inner: self.inner.mul(&parse_omnific(other)?),
        })
    }
    fn inv(&self) -> PyResult<PyOmnific> {
        self.inner
            .inv()
            .map(|o| PyOmnific { inner: o })
            .ok_or_else(|| PyValueError::new_err("Oz is a ring: only ±1 are invertible"))
    }
    fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
        matches!(parse_omnific(other), Ok(o) if o == self.inner)
    }
    fn __repr__(&self) -> String {
        format!("{:?}", self.inner.inner())
    }
}

pub(crate) fn parse_omnific(obj: &Bound<'_, PyAny>) -> PyResult<Omnific> {
    if let Ok(o) = obj.cast::<PyOmnific>() {
        return Ok(o.borrow().inner.clone());
    }
    if let Ok(s) = obj.cast::<PySurreal>() {
        return Omnific::from_surreal(s.borrow().inner.clone())
            .ok_or_else(|| PyValueError::new_err("surreal is not an omnific integer"));
    }
    if let Ok(v) = obj.extract::<i128>() {
        return Ok(Omnific::from_int(v));
    }
    Err(PyTypeError::new_err(
        "expected Omnific, omnific Surreal, or int",
    ))
}

pub(crate) fn wrap_omnific(o: Omnific) -> PyOmnific {
    PyOmnific { inner: o }
}

/// The omnific integer `n`.
#[pyfunction]
fn omnific(n: i128) -> PyOmnific {
    PyOmnific {
        inner: Omnific::from_int(n),
    }
}

/// `ω` as an omnific integer.
#[pyfunction]
fn omnific_omega() -> PyOmnific {
    PyOmnific {
        inner: Omnific::omega(),
    }
}

// ---------------------------------------------------------------------------
// Surreal builders
// ---------------------------------------------------------------------------

#[pyfunction]
fn omega() -> PySurreal {
    PySurreal {
        inner: Surreal::omega(),
    }
}

#[pyfunction]
fn epsilon() -> PySurreal {
    PySurreal {
        inner: Surreal::epsilon(),
    }
}

#[pyfunction]
fn omega_pow(exp: &Bound<'_, PyAny>) -> PyResult<PySurreal> {
    Ok(PySurreal {
        inner: Surreal::omega_pow(parse_surreal(exp)?),
    })
}

#[pyfunction]
fn rational(num: i128, den: i128) -> PyResult<PySurreal> {
    if den == 0 {
        return Err(PyValueError::new_err("zero denominator"));
    }
    Ok(PySurreal {
        inner: Surreal::from_rational(Rational::new(num, den)),
    })
}

#[pyfunction]
fn surreal(n: i128) -> PySurreal {
    PySurreal {
        inner: Surreal::from_int(n),
    }
}
// ---------------------------------------------------------------------------
// Nim field operations (the Artin–Schreier ↔ Arf bridge)
// ---------------------------------------------------------------------------

/// Nim square root (inverse Frobenius); always defined in char 2.
#[pyfunction]
fn nim_sqrt(x: u128) -> u128 {
    crate::scalar::nim_sqrt(x)
}

/// Field trace F_{2^m} → F₂ — the map the Arf invariant is read through and the
/// obstruction to solving `y²+y=c`.
#[pyfunction]
fn nim_trace(x: u128, m: u128) -> u128 {
    crate::scalar::nim_trace(x, m)
}

/// Solve the Artin–Schreier equation `y²+y=c` in F_{2^m} (`None` iff Tr(c)≠0).
#[pyfunction]
fn nim_solve_artin_schreier(c: u128, m: u128) -> Option<u128> {
    crate::scalar::nim_solve_artin_schreier(c, m)
}

/// Whether `y²+y=c` is solvable in F_{2^m} — i.e. `Tr(c)=0`.
#[pyfunction]
fn nim_is_artin_schreier_solvable(c: u128, m: u128) -> bool {
    crate::scalar::nim_is_artin_schreier_solvable(c, m)
}
// ---------------------------------------------------------------------------
// Transfinite (ordinal) nimbers
// ---------------------------------------------------------------------------

#[pyclass(name = "Ordinal", module = "pleroma", from_py_object)]
#[derive(Clone)]
struct PyOrdinal {
    inner: Ordinal,
}

#[pymethods]
impl PyOrdinal {
    #[new]
    fn new(n: u128) -> Self {
        PyOrdinal {
            inner: Ordinal::from_u128(n),
        }
    }
    /// `ω`, the first infinite ordinal nimber.
    #[staticmethod]
    fn omega() -> PyOrdinal {
        PyOrdinal {
            inner: Ordinal::omega(),
        }
    }
    /// `ω^exp` (coefficient 1).
    #[staticmethod]
    fn omega_pow(exp: &PyOrdinal) -> PyOrdinal {
        PyOrdinal {
            inner: Ordinal::omega_pow(exp.inner.clone()),
        }
    }
    /// `ω^exp · coeff`.
    #[staticmethod]
    fn monomial(exp: &PyOrdinal, coeff: u128) -> PyOrdinal {
        PyOrdinal {
            inner: Ordinal::monomial(exp.inner.clone(), coeff),
        }
    }
    /// Nim-addition (complete and exact): XOR of like-`ω`-power coefficients.
    fn nim_add(&self, other: &PyOrdinal) -> PyOrdinal {
        PyOrdinal {
            inner: self.inner.nim_add(&other.inner),
        }
    }
    /// Nim-multiplication (partial): exact for finite × finite; `None` when either
    /// operand is infinite (the general ordinal product is staged).
    fn nim_mul(&self, other: &PyOrdinal) -> Option<PyOrdinal> {
        self.inner
            .nim_mul(&other.inner)
            .map(|o| PyOrdinal { inner: o })
    }
    fn is_zero(&self) -> bool {
        self.inner.is_zero()
    }
    /// The finite nimber value, if this ordinal is finite.
    fn as_finite(&self) -> Option<u128> {
        self.inner.as_finite()
    }
    fn __richcmp__(&self, other: &PyOrdinal, op: CompareOp) -> bool {
        op.matches(self.inner.cmp(&other.inner))
    }
    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyNimber>()?;
    m.add_class::<PySurreal>()?;
    m.add_class::<PySurcomplex>()?;
    m.add_class::<PyInteger>()?;
    m.add_class::<PyOmnific>()?;
    m.add_class::<PyOrdinal>()?;
    m.add_function(wrap_pyfunction!(omnific, m)?)?;
    m.add_function(wrap_pyfunction!(omnific_omega, m)?)?;
    m.add_function(wrap_pyfunction!(omega, m)?)?;
    m.add_function(wrap_pyfunction!(epsilon, m)?)?;
    m.add_function(wrap_pyfunction!(omega_pow, m)?)?;
    m.add_function(wrap_pyfunction!(rational, m)?)?;
    m.add_function(wrap_pyfunction!(surreal, m)?)?;
    m.add_function(wrap_pyfunction!(nim_sqrt, m)?)?;
    m.add_function(wrap_pyfunction!(nim_trace, m)?)?;
    m.add_function(wrap_pyfunction!(nim_solve_artin_schreier, m)?)?;
    m.add_function(wrap_pyfunction!(nim_is_artin_schreier_solvable, m)?)?;
    Ok(())
}
