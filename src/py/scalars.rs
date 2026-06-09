//! Python bindings for the scalar worlds: the per-backend scalar types
//! (`Nimber`, `Surreal`, `Surcomplex`, `Integer`, `Omnific`, `Ordinal`), their
//! constructors, and the nim-field operations. `parse_*` / `wrap_*` are
//! `pub(crate)` because the `backend!` macro in [`super::engine`] threads them in
//! as the per-backend parse/wrap hooks.

use crate::scalar::{
    Adele, ExactRoots, FiniteField, Integer, LocalQp, MaxPlus, MinPlus, Nimber, Omnific, Ordinal,
    Rational, Scalar, Surcomplex, Surreal, Tropical,
};
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
    /// Degree over F₂ (dimension of the smallest nim-subfield containing it).
    fn degree(&self) -> u32 {
        crate::scalar::nim_degree(self.inner.0)
    }
    /// The Galois conjugates `x, x², x⁴, …` over F₂.
    fn conjugates(&self) -> Vec<PyNimber> {
        crate::scalar::nim_conjugates(self.inner.0)
            .into_iter()
            .map(|v| PyNimber { inner: Nimber(v) })
            .collect()
    }
    /// Minimal polynomial over F₂: coefficients `{0,1}` from the constant term up.
    fn min_poly(&self) -> Vec<u32> {
        crate::scalar::nim_min_poly(self.inner.0)
            .into_iter()
            .map(u32::from)
            .collect()
    }
    /// Multiplicative order in F_{2^128}* (`None` for `*0`).
    fn order(&self) -> Option<u128> {
        crate::scalar::nim_order(self.inner.0)
    }
    /// Whether this generates the full multiplicative group F_{2^128}*.
    fn is_primitive(&self) -> bool {
        crate::scalar::nim_is_primitive(self.inner.0)
    }
    /// Discrete log to base `self`: least `e` with `self**e == x`, else `None`.
    fn discrete_log(&self, x: &Bound<'_, PyAny>) -> PyResult<Option<u128>> {
        Ok(crate::scalar::nim_discrete_log(
            self.inner.0,
            parse_nimber(x)?.0,
        ))
    }
    /// `self` raised to the power `e` in `F_{2^128}` (fast exponentiation).
    fn pow(&self, e: u128) -> PyNimber {
        PyNimber {
            inner: self.inner.pow(e),
        }
    }
    fn __pow__(&self, e: u128, modulo: Option<u128>) -> PyResult<PyNimber> {
        if modulo.is_some() {
            return Err(PyValueError::new_err("nimber ** does not take a modulus"));
        }
        Ok(self.pow(e))
    }
    /// The Frobenius image `x²` — the generator of `Gal(F_{2^128}/F₂)`.
    fn frobenius(&self) -> PyNimber {
        PyNimber {
            inner: self.inner.frobenius(),
        }
    }
    /// The nim square root (unique in characteristic 2 — Frobenius is a bijection).
    fn sqrt(&self) -> PyNimber {
        PyNimber {
            inner: Nimber(crate::scalar::nim_sqrt(self.inner.0)),
        }
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
    /// True iff this surreal is a dyadic rational (a short-game number).
    fn is_dyadic(&self) -> bool {
        self.inner.is_dyadic()
    }
    /// The birthday of a dyadic rational; errors for non-dyadics (infinite
    /// birthday, outside this finite-support representation).
    fn dyadic_birthday(&self) -> PyResult<u128> {
        self.inner
            .dyadic_birthday()
            .ok_or_else(|| PyValueError::new_err("birthday is only finite for dyadic rationals"))
    }
    /// The simplest surreal strictly greater than this one (`{self|}`), when
    /// finite.
    fn simplest_above(&self) -> PyResult<PySurreal> {
        self.inner
            .simplest_above()
            .map(|inner| PySurreal { inner })
            .ok_or_else(|| PyValueError::new_err("simplest_above needs a finite rational"))
    }
    /// The simplest surreal strictly less than this one (`{|self}`), when finite.
    fn simplest_below(&self) -> PyResult<PySurreal> {
        self.inner
            .simplest_below()
            .map(|inner| PySurreal { inner })
            .ok_or_else(|| PyValueError::new_err("simplest_below needs a finite rational"))
    }
    /// The unique simplest surreal strictly between `a` and `b` (Conway's
    /// simplicity theorem), when it is dyadic. Errors if the endpoints are not
    /// finite rationals with `a < b`.
    #[staticmethod]
    fn simplest_between(a: &Bound<'_, PyAny>, b: &Bound<'_, PyAny>) -> PyResult<PySurreal> {
        let (a, b) = (parse_surreal(a)?, parse_surreal(b)?);
        Surreal::simplest_between(&a, &b)
            .map(|inner| PySurreal { inner })
            .ok_or_else(|| {
                PyValueError::new_err("no dyadic between (need finite rationals with a < b)")
            })
    }
    /// The floor ⌊x⌋ as a surreal — the greatest omnific integer ≤ x.
    fn floor(&self) -> PySurreal {
        PySurreal {
            inner: self.inner.floor(),
        }
    }
    /// The floor ⌊x⌋ as an `Omnific` integer.
    fn omnific_floor(&self) -> PyOmnific {
        wrap_omnific(Omnific::floor(&self.inner))
    }
    /// The fractional part `x − ⌊x⌋`, in `[0, 1)`.
    fn frac(&self) -> PySurreal {
        PySurreal {
            inner: self.inner.frac(),
        }
    }
    /// The sign expansion (`True = +`) of a dyadic surreal; `None` for
    /// non-dyadics (transfinite expansion). Length equals the birthday.
    fn sign_expansion(&self) -> Option<Vec<bool>> {
        self.inner.sign_expansion()
    }
    /// The dyadic surreal with the given finite sign expansion (`True = +`).
    #[staticmethod]
    fn from_sign_expansion(signs: Vec<bool>) -> PySurreal {
        PySurreal {
            inner: Surreal::from_sign_expansion(&signs),
        }
    }
    /// The **truncated inverse** `1/x` to `n` leading terms (Neumann series) —
    /// works for non-monomials too, unlike [`inv`](Self::inv). Errors on `0`.
    fn inv_to_terms(&self, n: usize) -> PyResult<PySurreal> {
        self.inner
            .inv_to_terms(n)
            .map(|inner| PySurreal { inner })
            .ok_or_else(|| PyValueError::new_err("0 has no inverse"))
    }
    /// The **truncated real square root** to `n` leading terms (the lazy
    /// `SeriesRoots` primitive); `None` unless the leading coefficient is a perfect
    /// ℚ-square and the value is ≥ 0 (so `√2` and `√(2ω)` are `None`, while
    /// `√ω = ω^{1/2}` is exact). For the precision-free exact value see
    /// [`exact_sqrt`](Self::exact_sqrt).
    fn sqrt(&self, n: usize) -> Option<PySurreal> {
        self.inner.sqrt_to_terms(n).map(|inner| PySurreal { inner })
    }
    /// The **truncated real `k`-th root** to `n` leading terms (same ℚ-power scope).
    fn nth_root(&self, k: u32, n: usize) -> Option<PySurreal> {
        self.inner
            .nth_root_to_terms(k, n)
            .map(|inner| PySurreal { inner })
    }
    /// Whether this is a square in the represented surreal subdomain (`ExactRoots`).
    fn is_square(&self) -> bool {
        ExactRoots::is_square(&self.inner)
    }
    /// The **exact** real square root (no precision argument): `Some` iff a finite
    /// represented surreal squares back to this — `√ω = ω^{1/2}`, `√4 = 2`, but
    /// `√2` is `None`. The exact companion to [`sqrt`](Self::sqrt).
    fn exact_sqrt(&self) -> Option<PySurreal> {
        ExactRoots::sqrt(&self.inner).map(|inner| PySurreal { inner })
    }
    /// The **birthday** as an `Ordinal` (transfinite-aware): `ω ↦ ω`, `ε ↦ ω`,
    /// `ω^ω ↦ ω^ω`. `None` outside the representable subclass (`√ω`, …).
    fn birthday_ordinal(&self) -> Option<PyOrdinal> {
        self.inner
            .birthday_ordinal()
            .map(|inner| PyOrdinal { inner })
    }
    /// The (possibly transfinite) **sign expansion** as runs `(sign, length)`
    /// (`True = +`, length an `Ordinal`); `None` outside the representable
    /// subclass.
    fn transfinite_sign_expansion(&self) -> Option<Vec<(bool, PyOrdinal)>> {
        self.inner.transfinite_sign_expansion().map(|se| {
            se.runs()
                .iter()
                .map(|(s, l)| (*s, PyOrdinal { inner: l.clone() }))
                .collect()
        })
    }
    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}

impl PySurreal {
    /// Wrap a core `Surreal` (used by the games↔surreal bridge).
    pub(crate) fn from_inner(inner: Surreal) -> PySurreal {
        PySurreal { inner }
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
    /// Whether this is a square in the surcomplex field (`ExactRoots`).
    fn is_square(&self) -> bool {
        ExactRoots::is_square(&self.inner)
    }
    /// The **exact algebraic-closure square root** `√(a+bi)`: the surcomplex field
    /// is algebraically closed over its real-closed base, so a represented number
    /// has a represented root. `None` outside the represented square subdomain
    /// (e.g. `√2`). The functorial companion of `Surreal.exact_sqrt`.
    fn sqrt(&self) -> Option<PySurcomplex> {
        ExactRoots::sqrt(&self.inner).map(|inner| PySurcomplex { inner })
    }
    /// The **truncated inverse** `1/(a+bi)` to `n` leading terms — succeeds where
    /// [`inv`](Self::inv) returns `None` because the norm `a²+b²` is a non-monomial
    /// surreal. Errors only on `0`.
    fn inv_to_terms(&self, n: usize) -> PyResult<PySurcomplex> {
        self.inner
            .inv_to_terms(n)
            .map(|inner| PySurcomplex { inner })
            .ok_or_else(|| PyValueError::new_err("0 has no inverse"))
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

// --- Runtime local/global precision models: LocalQp and Adele ---------------

fn is_prime_u128(p: u128) -> bool {
    if p < 2 {
        return false;
    }
    if p.is_multiple_of(2) {
        return p == 2;
    }
    let mut d = 3u128;
    while d <= p / d {
        if p.is_multiple_of(d) {
            return false;
        }
        d += 2;
    }
    true
}

fn checked_pow_u128(base: u128, exp: u128) -> Option<u128> {
    if exp > 128 {
        return None;
    }
    let mut acc = 1u128;
    for _ in 0..exp {
        acc = acc.checked_mul(base)?;
    }
    Some(acc)
}

fn validate_local_qp_world(p: u128, k: u128) -> PyResult<()> {
    if !is_prime_u128(p) || k == 0 {
        return Err(PyValueError::new_err(
            "LocalQp needs a prime p and positive precision k",
        ));
    }
    if p > i128::MAX as u128 {
        return Err(PyValueError::new_err(
            "LocalQp prime must fit the bounded i128 arithmetic model",
        ));
    }
    if checked_pow_u128(p, k).is_none() {
        return Err(PyValueError::new_err("LocalQp modulus p^k exceeds u128"));
    }
    let Some(double_k) = k.checked_mul(2) else {
        return Err(PyValueError::new_err(
            "LocalQp precision is too large for checked arithmetic",
        ));
    };
    if checked_pow_u128(p, double_k).is_none() {
        return Err(PyValueError::new_err(
            "LocalQp needs (p^k)^2 to fit u128 for safe mantissa arithmetic",
        ));
    }
    Ok(())
}

fn rational_from_pair(num: i128, den: i128) -> PyResult<Rational> {
    Rational::try_new(num, den)
        .ok_or_else(|| PyValueError::new_err("zero denominator or bounded i128 overflow"))
}

fn rational_pair(q: &Rational) -> (i128, i128) {
    (q.numer(), q.denom())
}

fn adele_precision_for_prime(p: u128) -> PyResult<u128> {
    if !is_prime_u128(p) || p > i128::MAX as u128 {
        return Err(PyValueError::new_err(
            "adele finite-place precision needs a prime p within i128",
        ));
    }
    Ok(crate::scalar::global::adele::adele_prec(p))
}

#[pyclass(name = "LocalQp", module = "pleroma", from_py_object)]
#[derive(Clone)]
pub(crate) struct PyLocalQp {
    inner: LocalQp,
}

fn parse_local_qp_in_world(obj: &Bound<'_, PyAny>, p: u128, k: u128) -> PyResult<LocalQp> {
    validate_local_qp_world(p, k)?;
    if let Ok(x) = obj.cast::<PyLocalQp>() {
        let x = x.borrow().inner;
        if x.prime() != p || x.precision() != k {
            return Err(PyValueError::new_err(format!(
                "cannot mix LocalQp worlds Q_{}@{} and Q_{}@{}",
                p,
                k,
                x.prime(),
                x.precision()
            )));
        }
        return Ok(x);
    }
    if let Ok(v) = obj.extract::<i128>() {
        return Ok(LocalQp::from_i128(p, k, v));
    }
    Err(PyTypeError::new_err(
        "expected LocalQp from the same (p,k) world or int",
    ))
}

#[pymethods]
impl PyLocalQp {
    #[new]
    #[pyo3(signature = (p, k, value=0))]
    fn new(p: u128, k: u128, value: i128) -> PyResult<Self> {
        validate_local_qp_world(p, k)?;
        Ok(PyLocalQp {
            inner: LocalQp::from_i128(p, k, value),
        })
    }
    #[staticmethod]
    fn zero(p: u128, k: u128) -> PyResult<PyLocalQp> {
        validate_local_qp_world(p, k)?;
        Ok(PyLocalQp {
            inner: LocalQp::zero(p, k),
        })
    }
    #[staticmethod]
    fn one(p: u128, k: u128) -> PyResult<PyLocalQp> {
        validate_local_qp_world(p, k)?;
        Ok(PyLocalQp {
            inner: LocalQp::one(p, k),
        })
    }
    #[staticmethod]
    fn from_p_power(p: u128, k: u128, v: i128) -> PyResult<PyLocalQp> {
        validate_local_qp_world(p, k)?;
        Ok(PyLocalQp {
            inner: LocalQp::from_p_power(p, k, v),
        })
    }
    #[staticmethod]
    fn from_rational(p: u128, k: u128, num: i128, den: i128) -> PyResult<PyLocalQp> {
        validate_local_qp_world(p, k)?;
        let q = rational_from_pair(num, den)?;
        Ok(PyLocalQp {
            inner: LocalQp::from_rational(p, k, &q),
        })
    }
    #[getter]
    fn prime(&self) -> u128 {
        self.inner.prime()
    }
    #[getter]
    fn precision(&self) -> u128 {
        self.inner.precision()
    }
    #[getter]
    fn unit(&self) -> u128 {
        self.inner.unit()
    }
    fn modulus(&self) -> u128 {
        self.inner.modulus()
    }
    fn valuation(&self) -> Option<i128> {
        self.inner.valuation()
    }
    fn is_zero(&self) -> bool {
        self.inner.is_zero()
    }
    fn inv(&self) -> PyResult<PyLocalQp> {
        self.inner
            .inv()
            .map(|inner| PyLocalQp { inner })
            .ok_or_else(|| PyValueError::new_err("0 has no p-adic inverse"))
    }
    fn __add__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyLocalQp> {
        let rhs = parse_local_qp_in_world(other, self.inner.prime(), self.inner.precision())?;
        Ok(PyLocalQp {
            inner: self.inner.add(&rhs),
        })
    }
    fn __radd__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyLocalQp> {
        self.__add__(other)
    }
    fn __sub__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyLocalQp> {
        let rhs = parse_local_qp_in_world(other, self.inner.prime(), self.inner.precision())?;
        Ok(PyLocalQp {
            inner: self.inner.add(&rhs.neg()),
        })
    }
    fn __rsub__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyLocalQp> {
        let lhs = parse_local_qp_in_world(other, self.inner.prime(), self.inner.precision())?;
        Ok(PyLocalQp {
            inner: lhs.add(&self.inner.neg()),
        })
    }
    fn __neg__(&self) -> PyLocalQp {
        PyLocalQp {
            inner: self.inner.neg(),
        }
    }
    fn __mul__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyLocalQp> {
        let rhs = parse_local_qp_in_world(other, self.inner.prime(), self.inner.precision())?;
        Ok(PyLocalQp {
            inner: self.inner.mul(&rhs),
        })
    }
    fn __rmul__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyLocalQp> {
        self.__mul__(other)
    }
    fn __truediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyLocalQp> {
        let rhs = parse_local_qp_in_world(other, self.inner.prime(), self.inner.precision())?;
        let rinv = rhs
            .inv()
            .ok_or_else(|| PyValueError::new_err("division by 0 in LocalQp"))?;
        Ok(PyLocalQp {
            inner: self.inner.mul(&rinv),
        })
    }
    fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
        matches!(other.cast::<PyLocalQp>(), Ok(x) if x.borrow().inner == self.inner)
    }
    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}

fn parse_adele(obj: &Bound<'_, PyAny>) -> PyResult<Adele> {
    if let Ok(a) = obj.cast::<PyAdele>() {
        return Ok(a.borrow().inner.clone());
    }
    if let Ok(v) = obj.extract::<i128>() {
        return Ok(Adele::from_rational(&Rational::int(v)));
    }
    Err(PyTypeError::new_err("expected Adele or int"))
}

#[pyclass(name = "Adele", module = "pleroma", from_py_object)]
#[derive(Clone)]
pub(crate) struct PyAdele {
    inner: Adele,
}

#[pymethods]
impl PyAdele {
    #[new]
    #[pyo3(signature = (num=0, den=1))]
    fn new(num: i128, den: i128) -> PyResult<Self> {
        let q = rational_from_pair(num, den)?;
        Ok(PyAdele {
            inner: Adele::from_rational(&q),
        })
    }
    #[staticmethod]
    fn from_rational(num: i128, den: i128) -> PyResult<PyAdele> {
        PyAdele::new(num, den)
    }
    #[staticmethod]
    fn finite_precision(p: u128) -> PyResult<u128> {
        adele_precision_for_prime(p)
    }
    fn with_correction(&self, p: u128, dev: &PyLocalQp) -> PyResult<PyAdele> {
        let expected = adele_precision_for_prime(p)?;
        if dev.inner.prime() != p || dev.inner.precision() != expected {
            return Err(PyValueError::new_err(format!(
                "Adele correction at p={p} must be LocalQp(p={p}, k={expected})"
            )));
        }
        Ok(PyAdele {
            inner: self.inner.clone().with_correction(p, dev.inner),
        })
    }
    fn with_archimedean(&self, num: i128, den: i128) -> PyResult<PyAdele> {
        Ok(PyAdele {
            inner: self
                .inner
                .clone()
                .with_archimedean(rational_from_pair(num, den)?),
        })
    }
    #[getter]
    fn principal(&self) -> (i128, i128) {
        rational_pair(self.inner.principal())
    }
    #[getter]
    fn archimedean(&self) -> (i128, i128) {
        rational_pair(self.inner.archimedean())
    }
    fn local_at(&self, p: u128) -> PyResult<PyLocalQp> {
        adele_precision_for_prime(p)?;
        Ok(PyLocalQp {
            inner: self.inner.local_at(p),
        })
    }
    fn is_principal(&self) -> bool {
        self.inner.is_principal()
    }
    fn is_idele(&self) -> bool {
        self.inner.is_idele()
    }
    fn is_integral(&self) -> bool {
        self.inner.is_integral()
    }
    fn absolute_value_real(&self) -> (i128, i128) {
        rational_pair(
            &self
                .inner
                .absolute_value_at(crate::scalar::AdelePlace::Real),
        )
    }
    fn absolute_value_at(&self, p: u128) -> PyResult<(i128, i128)> {
        adele_precision_for_prime(p)?;
        Ok(rational_pair(
            &self
                .inner
                .absolute_value_at(crate::scalar::AdelePlace::Prime(p)),
        ))
    }
    fn idele_norm(&self) -> (i128, i128) {
        rational_pair(&self.inner.idele_norm())
    }
    fn satisfies_product_formula(&self) -> bool {
        self.inner.satisfies_product_formula()
    }
    fn inv(&self) -> PyResult<PyAdele> {
        self.inner
            .inv()
            .map(|inner| PyAdele { inner })
            .ok_or_else(|| PyValueError::new_err("adele is not an idele"))
    }
    fn __add__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyAdele> {
        Ok(PyAdele {
            inner: self.inner.add(&parse_adele(other)?),
        })
    }
    fn __radd__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyAdele> {
        self.__add__(other)
    }
    fn __sub__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyAdele> {
        Ok(PyAdele {
            inner: self.inner.sub(&parse_adele(other)?),
        })
    }
    fn __rsub__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyAdele> {
        Ok(PyAdele {
            inner: parse_adele(other)?.sub(&self.inner),
        })
    }
    fn __neg__(&self) -> PyAdele {
        PyAdele {
            inner: self.inner.neg(),
        }
    }
    fn __mul__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyAdele> {
        Ok(PyAdele {
            inner: self.inner.mul(&parse_adele(other)?),
        })
    }
    fn __rmul__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyAdele> {
        self.__mul__(other)
    }
    fn __truediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyAdele> {
        let rhs = parse_adele(other)?;
        let rinv = rhs
            .inv()
            .ok_or_else(|| PyValueError::new_err("Adele divisor is not an idele"))?;
        Ok(PyAdele {
            inner: self.inner.mul(&rinv),
        })
    }
    fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
        matches!(other.cast::<PyAdele>(), Ok(a) if a.borrow().inner == self.inner)
    }
    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}

macro_rules! tropical_pyclass {
    ($py:ident, $name:literal, $conv:ty) => {
        #[pyclass(name = $name, module = "pleroma", from_py_object)]
        #[derive(Clone)]
        pub(crate) struct $py {
            inner: Tropical<$conv>,
        }

        #[pymethods]
        impl $py {
            #[new]
            #[pyo3(signature = (num=None, den=1))]
            fn new(num: Option<i128>, den: i128) -> PyResult<Self> {
                match num {
                    Some(num) => Ok($py {
                        inner: Tropical::<$conv>::finite(rational_from_pair(num, den)?),
                    }),
                    None => Ok($py {
                        inner: Tropical::<$conv>::infinity(),
                    }),
                }
            }
            #[staticmethod]
            fn finite(num: i128, den: i128) -> PyResult<$py> {
                Ok($py {
                    inner: Tropical::<$conv>::finite(rational_from_pair(num, den)?),
                })
            }
            #[staticmethod]
            fn infinity() -> $py {
                $py {
                    inner: Tropical::<$conv>::infinity(),
                }
            }
            #[staticmethod]
            fn zero() -> $py {
                $py {
                    inner: Tropical::<$conv>::zero(),
                }
            }
            #[staticmethod]
            fn one() -> $py {
                $py {
                    inner: Tropical::<$conv>::one(),
                }
            }
            fn value(&self) -> Option<(i128, i128)> {
                self.inner.value().map(|q| rational_pair(&q))
            }
            fn is_infinity(&self) -> bool {
                self.inner.is_infinity()
            }
            fn __add__(&self, other: &$py) -> $py {
                $py {
                    inner: self.inner.add(&other.inner),
                }
            }
            fn __mul__(&self, other: &$py) -> $py {
                $py {
                    inner: self.inner.mul(&other.inner),
                }
            }
            fn __eq__(&self, other: &$py) -> bool {
                self.inner == other.inner
            }
            fn __repr__(&self) -> String {
                format!("{:?}", self.inner)
            }
        }
    };
}

tropical_pyclass!(PyMaxPlusTropical, "MaxPlusTropical", MaxPlus);
tropical_pyclass!(PyMinPlusTropical, "MinPlusTropical", MinPlus);

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

// --- finite-field analysis toolkit ---

/// Degree of `x` over F₂ (the smallest nim-subfield F_{2^d} containing it).
#[pyfunction]
fn nim_degree(x: u128) -> u32 {
    crate::scalar::nim_degree(x)
}

/// The Galois conjugates `x, x², x⁴, …` over F₂.
#[pyfunction]
fn nim_conjugates(x: u128) -> Vec<u128> {
    crate::scalar::nim_conjugates(x)
}

/// Minimal polynomial of `x` over F₂: coefficients `{0,1}` from the constant up.
#[pyfunction]
fn nim_min_poly(x: u128) -> Vec<u32> {
    crate::scalar::nim_min_poly(x)
        .into_iter()
        .map(u32::from)
        .collect()
}

/// Relative trace `Tr_{F_{2^m}/F_{2^e}}(x)` (the `e=1` case is `nim_trace`).
#[pyfunction]
fn nim_relative_trace(x: u128, m: u32, e: u32) -> u128 {
    crate::scalar::nim_relative_trace(x, m, e)
}

/// Relative norm `N_{F_{2^m}/F_{2^e}}(x)` (norm to the prime field is trivial).
#[pyfunction]
fn nim_relative_norm(x: u128, m: u32, e: u32) -> u128 {
    crate::scalar::nim_relative_norm(x, m, e)
}

/// Multiplicative order of `x` in F_{2^128}* (`None` for `*0`).
#[pyfunction]
fn nim_order(x: u128) -> Option<u128> {
    crate::scalar::nim_order(x)
}

/// Whether `x` generates the full group F_{2^128}*.
#[pyfunction]
fn nim_is_primitive(x: u128) -> bool {
    crate::scalar::nim_is_primitive(x)
}

/// A primitive element (generator) of F_{2^128}*.
#[pyfunction]
fn nim_primitive_element() -> u128 {
    crate::scalar::nim_primitive_element()
}

/// Discrete log in F_{2^128}*: least `e` with `base**e == x`, else `None`.
#[pyfunction]
fn nim_discrete_log(base: u128, x: u128) -> Option<u128> {
    crate::scalar::nim_discrete_log(base, x)
}
// ---------------------------------------------------------------------------
// Transfinite (ordinal) nimbers
// ---------------------------------------------------------------------------

#[pyclass(name = "Ordinal", module = "pleroma", from_py_object)]
#[derive(Clone)]
pub(crate) struct PyOrdinal {
    inner: Ordinal,
}

impl PyOrdinal {
    pub(crate) fn from_inner(inner: Ordinal) -> PyOrdinal {
        PyOrdinal { inner }
    }
    pub(crate) fn as_ordinal(&self) -> &Ordinal {
        &self.inner
    }
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
    /// **Ordinary** (Cantor) ordinal addition — NOT nim: `1 + ω = ω` but
    /// `ω + ω = ω·2` (coefficients add as naturals, not XOR).
    fn ord_add(&self, other: &PyOrdinal) -> PyOrdinal {
        PyOrdinal {
            inner: self.inner.ord_add(&other.inner),
        }
    }
    /// **Ordinary** (Cantor) ordinal multiplication — NOT nim (`2·ω = ω`).
    fn ord_mul(&self, other: &PyOrdinal) -> PyOrdinal {
        PyOrdinal {
            inner: self.inner.ord_mul(&other.inner),
        }
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
    m.add_class::<PyLocalQp>()?;
    m.add_class::<PyAdele>()?;
    m.add_class::<PyMaxPlusTropical>()?;
    m.add_class::<PyMinPlusTropical>()?;
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
    m.add_function(wrap_pyfunction!(nim_degree, m)?)?;
    m.add_function(wrap_pyfunction!(nim_conjugates, m)?)?;
    m.add_function(wrap_pyfunction!(nim_min_poly, m)?)?;
    m.add_function(wrap_pyfunction!(nim_relative_trace, m)?)?;
    m.add_function(wrap_pyfunction!(nim_relative_norm, m)?)?;
    m.add_function(wrap_pyfunction!(nim_order, m)?)?;
    m.add_function(wrap_pyfunction!(nim_is_primitive, m)?)?;
    m.add_function(wrap_pyfunction!(nim_primitive_element, m)?)?;
    m.add_function(wrap_pyfunction!(nim_discrete_log, m)?)?;
    Ok(())
}
