//! Python bindings for the scalar worlds: the per-backend scalar types
//! (`Nimber`, finite fields, fixed p-adic slices, `Rational`, `Surreal`,
//! `Surcomplex`, `Integer`, `Omnific`, `Ordinal`), their constructors, and the
//! nim-field operations. `parse_*` / `wrap_*` are `pub(crate)` because the
//! `backend!` macro in [`super::engine`] threads them in as the per-backend
//! parse/wrap hooks.

use crate::scalar::{
    is_prime_u128, Adele, AdelePlace, CyclicGaloisExtension, ExactRoots, FieldExtension,
    FiniteField, Fp, Fpn, Gauss, HasFractionField, HasRingOfIntegers, Integer, Laurent, LocalQp,
    MaxPlus, MinPlus, Nimber, Omnific, Ordinal, Poly, Qp, Qq, Ramified, Rational, RationalFunction,
    ReductionPolynomialKind, ResidueField, Scalar, SignExpansion, Surcomplex, Surreal, Tropical,
    Valued, WittVec, Zp,
};
use pyo3::basic::CompareOp;
use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::IntoPyObjectExt;
use std::cmp::Ordering;

fn ordering_to_i8(ordering: Ordering) -> i8 {
    match ordering {
        Ordering::Less => -1,
        Ordering::Equal => 0,
        Ordering::Greater => 1,
    }
}

fn validate_relative_degrees<F: FiniteField>(x: &F, m: usize, e: usize) -> PyResult<()> {
    if e == 0 || !m.is_multiple_of(e) {
        return Err(PyValueError::new_err(
            "relative trace/norm needs e > 0 and e | m",
        ));
    }
    let ext = F::ext_degree();
    if m == 0 || !ext.is_multiple_of(m) {
        return Err(PyValueError::new_err(
            "relative trace/norm needs m a positive divisor of the represented field degree",
        ));
    }
    if !m.is_multiple_of(x.degree()) {
        return Err(PyValueError::new_err(
            "relative trace/norm input is not contained in the requested degree-m subfield",
        ));
    }
    Ok(())
}

fn qp_to_qq_base<const P: u128, const N: usize, const K: u128>(x: Qp<P, K>) -> Qq<P, N, 1> {
    debug_assert_eq!(N as u128, K, "Python fixed Qq/Qp precisions must match");
    match x.valuation() {
        None => Qq::<P, N, 1>::zero(),
        Some(v) => {
            Qq::<P, N, 1>::from_p_power(v).mul(&Qq::from_witt(WittVec::<P, N, 1>([x.unit()])))
        }
    }
}

fn qq_base_to_qp<const P: u128, const N: usize, const K: u128>(x: Qq<P, N, 1>) -> Qp<P, K> {
    debug_assert_eq!(N as u128, K, "Python fixed Qq/Qp precisions must match");
    match x.valuation() {
        None => Qp::<P, K>::zero(),
        Some(v) => {
            let unit = i128::try_from(x.unit().0[0]).expect("Python fixed Qq unit fits i128");
            Qp::<P, K>::from_i128(unit).mul(&Qp::<P, K>::from_p_power(v))
        }
    }
}

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
    #[staticmethod]
    fn zero() -> Self {
        wrap_nimber(Nimber::zero())
    }
    #[staticmethod]
    fn one() -> Self {
        wrap_nimber(Nimber::one())
    }
    #[staticmethod]
    fn characteristic() -> u128 {
        Nimber::characteristic()
    }
    #[getter]
    fn value(&self) -> u128 {
        self.inner.0
    }
    fn is_zero(&self) -> bool {
        self.inner.is_zero()
    }
    fn __add__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyNimber> {
        Ok(PyNimber {
            inner: self.inner.add(&parse_nimber(other)?),
        })
    }
    fn __radd__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyNimber> {
        self.__add__(other)
    }
    fn __sub__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyNimber> {
        Ok(PyNimber {
            inner: self.inner.sub(&parse_nimber(other)?),
        })
    }
    fn __rsub__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyNimber> {
        Ok(PyNimber {
            inner: parse_nimber(other)?.sub(&self.inner),
        })
    }
    fn __neg__(&self) -> PyNimber {
        PyNimber { inner: self.inner }
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
    fn __rtruediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyNimber> {
        let si = self
            .inner
            .inv()
            .ok_or_else(|| PyValueError::new_err("division by *0"))?;
        Ok(PyNimber {
            inner: parse_nimber(other)?.mul(&si),
        })
    }
    fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
        matches!(other.cast::<PyNimber>(), Ok(n) if n.borrow().inner == self.inner)
    }
    fn __hash__(&self) -> usize {
        self.inner.0 as usize
    }
    /// Degree over F₂ (dimension of the smallest nim-subfield containing it).
    fn degree(&self) -> u128 {
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
    fn min_poly(&self) -> Vec<u128> {
        crate::scalar::nim_min_poly(self.inner.0)
    }
    /// Monic minimal polynomial over F₂ as nimber coefficients.
    fn min_poly_monic(&self) -> Vec<PyNimber> {
        self.inner
            .min_poly_monic()
            .into_iter()
            .map(wrap_nimber)
            .collect()
    }
    /// Multiplicative order in F_{2^128}* (`None` for `*0`).
    fn multiplicative_order(&self) -> Option<u128> {
        self.inner.multiplicative_order()
    }
    /// Whether this generates the full multiplicative group F_{2^128}*.
    fn is_primitive(&self) -> bool {
        crate::scalar::nim_is_primitive(self.inner.0)
    }
    /// The extension degree `[F_{2^128}:F_2]`.
    #[staticmethod]
    fn ext_degree() -> usize {
        <Nimber as FiniteField>::ext_degree()
    }
    /// The multiplicative group order `|F_{2^128}*|`.
    #[staticmethod]
    fn group_order() -> u128 {
        <Nimber as FiniteField>::group_order()
    }
    /// The distinct prime factors of `group_order()`.
    #[staticmethod]
    fn group_order_factors() -> Vec<u128> {
        <Nimber as FiniteField>::group_order_factors()
    }
    /// Discrete log to base `self`: least `e` with `self**e == x`, else `None`.
    fn discrete_log(&self, x: &Bound<'_, PyAny>) -> PyResult<Option<u128>> {
        Ok(crate::scalar::nim_discrete_log(
            self.inner.0,
            parse_nimber(x)?.0,
        ))
    }
    /// Relative trace `Tr_{F_{2^m}/F_{2^e}}(self)`, returned as a nimber.
    fn relative_trace_over(&self, m: usize, e: usize) -> PyResult<PyNimber> {
        validate_relative_degrees(&self.inner, m, e)?;
        Ok(wrap_nimber(self.inner.relative_trace_over(m, e)))
    }
    /// Relative norm `N_{F_{2^m}/F_{2^e}}(self)`, returned as a nimber.
    fn relative_norm_over(&self, m: usize, e: usize) -> PyResult<PyNimber> {
        validate_relative_degrees(&self.inner, m, e)?;
        Ok(wrap_nimber(self.inner.relative_norm_over(m, e)))
    }
    /// Trace from the full field `F_{2^128}` to `F_{2^e}`.
    fn relative_trace(&self, e: usize) -> PyResult<PyNimber> {
        validate_relative_degrees(&self.inner, <Nimber as FiniteField>::ext_degree(), e)?;
        Ok(wrap_nimber(self.inner.relative_trace(e)))
    }
    /// Norm from the full field `F_{2^128}` to `F_{2^e}`.
    fn relative_norm(&self, e: usize) -> PyResult<PyNimber> {
        validate_relative_degrees(&self.inner, <Nimber as FiniteField>::ext_degree(), e)?;
        Ok(wrap_nimber(self.inner.relative_norm(e)))
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
    /// The Frobenius applied `k` times: `x ↦ x^(2^k)`.
    fn frobenius_iter(&self, k: usize) -> PyNimber {
        PyNimber {
            inner: self.inner.frobenius_iter(k),
        }
    }
    /// The nim square root (unique in characteristic 2 — Frobenius is a bijection).
    fn sqrt(&self) -> PyNimber {
        PyNimber {
            inner: Nimber(crate::scalar::nim_sqrt(self.inner.0)),
        }
    }
    /// Whether this is a square in `F_{2^128}` (always true in characteristic 2).
    fn is_square(&self) -> bool {
        ExactRoots::is_square(&self.inner)
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

#[pyclass(name = "NimberPoly", module = "pleroma", from_py_object)]
#[derive(Clone)]
pub(crate) struct PyNimberPoly {
    inner: Poly<Nimber>,
}

fn parse_nimber_coeffs(items: Vec<Bound<'_, PyAny>>) -> PyResult<Vec<Nimber>> {
    items.iter().map(parse_nimber).collect()
}

#[pymethods]
impl PyNimberPoly {
    #[new]
    fn new(coeffs: Vec<Bound<'_, PyAny>>) -> PyResult<Self> {
        Ok(PyNimberPoly {
            inner: Poly::new(parse_nimber_coeffs(coeffs)?),
        })
    }
    #[staticmethod]
    fn zero() -> Self {
        wrap_nimber_poly(Poly::zero())
    }
    #[staticmethod]
    fn one() -> Self {
        wrap_nimber_poly(Poly::one())
    }
    #[staticmethod]
    fn characteristic() -> u128 {
        <Poly<Nimber> as Scalar>::characteristic()
    }
    #[staticmethod]
    fn x() -> Self {
        wrap_nimber_poly(Poly::x())
    }
    #[staticmethod]
    fn constant(s: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(wrap_nimber_poly(Poly::constant(parse_nimber(s)?)))
    }
    #[staticmethod]
    fn monomial(deg: usize, coeff: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(wrap_nimber_poly(Poly::monomial(deg, parse_nimber(coeff)?)))
    }
    #[getter]
    fn coeffs(&self) -> Vec<PyNimber> {
        self.inner
            .coeffs()
            .iter()
            .copied()
            .map(wrap_nimber)
            .collect()
    }
    #[getter]
    fn degree(&self) -> Option<usize> {
        self.inner.degree()
    }
    fn leading(&self) -> Option<PyNimber> {
        self.inner.leading().copied().map(wrap_nimber)
    }
    fn coeff(&self, i: usize) -> PyNimber {
        wrap_nimber(self.inner.coeff(i))
    }
    fn is_zero(&self) -> bool {
        self.inner.is_zero()
    }
    fn eval(&self, x: &Bound<'_, PyAny>) -> PyResult<PyNimber> {
        Ok(wrap_nimber(self.inner.eval(&parse_nimber(x)?)))
    }
    fn scale(&self, s: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(wrap_nimber_poly(self.inner.scale(&parse_nimber(s)?)))
    }
    fn make_monic(&self) -> PyResult<Self> {
        if self.inner.is_zero() {
            return Err(PyValueError::new_err("make_monic of the zero polynomial"));
        }
        Ok(wrap_nimber_poly(self.inner.make_monic()))
    }
    fn divrem(&self, divisor: &PyNimberPoly) -> PyResult<(PyNimberPoly, PyNimberPoly)> {
        if divisor.inner.is_zero() {
            return Err(PyValueError::new_err("polynomial division by zero"));
        }
        let (q, r) = self.inner.divrem(&divisor.inner);
        Ok((wrap_nimber_poly(q), wrap_nimber_poly(r)))
    }
    fn rem(&self, divisor: &PyNimberPoly) -> PyResult<PyNimberPoly> {
        if divisor.inner.is_zero() {
            return Err(PyValueError::new_err("polynomial remainder by zero"));
        }
        Ok(wrap_nimber_poly(self.inner.rem(&divisor.inner)))
    }
    fn divides(&self, multiple: &PyNimberPoly) -> bool {
        self.inner.divides(&multiple.inner)
    }
    fn gcd(&self, other: &PyNimberPoly) -> PyNimberPoly {
        wrap_nimber_poly(self.inner.gcd(&other.inner))
    }
    fn mul_mod(&self, other: &PyNimberPoly, modulus: &PyNimberPoly) -> PyResult<PyNimberPoly> {
        if modulus.inner.is_zero() {
            return Err(PyValueError::new_err("polynomial modulus is zero"));
        }
        Ok(wrap_nimber_poly(
            self.inner.mul_mod(&other.inner, &modulus.inner),
        ))
    }
    fn pow_mod(&self, e: u128, modulus: &PyNimberPoly) -> PyResult<PyNimberPoly> {
        if modulus.inner.is_zero() {
            return Err(PyValueError::new_err("polynomial modulus is zero"));
        }
        Ok(wrap_nimber_poly(self.inner.pow_mod(e, &modulus.inner)))
    }
    fn inv(&self) -> PyResult<PyNimberPoly> {
        self.inner
            .inv()
            .map(wrap_nimber_poly)
            .ok_or_else(|| PyValueError::new_err("only nonzero constant polynomials invert"))
    }
    fn to_fraction(&self) -> PyNimberRationalFunction {
        wrap_nimber_rational_function(<Poly<Nimber> as HasFractionField>::to_fraction(&self.inner))
    }
    fn __add__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyNimberPoly> {
        Ok(wrap_nimber_poly(self.inner.add(&parse_nimber_poly(other)?)))
    }
    fn __radd__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyNimberPoly> {
        self.__add__(other)
    }
    fn __sub__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyNimberPoly> {
        Ok(wrap_nimber_poly(self.inner.sub(&parse_nimber_poly(other)?)))
    }
    fn __rsub__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyNimberPoly> {
        Ok(wrap_nimber_poly(parse_nimber_poly(other)?.sub(&self.inner)))
    }
    fn __neg__(&self) -> PyNimberPoly {
        wrap_nimber_poly(self.inner.neg())
    }
    fn __mul__(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        match parse_nimber_poly(other) {
            Ok(o) => wrap_nimber_poly(self.inner.mul(&o)).into_py_any(py),
            Err(_) => Ok(py.NotImplemented()),
        }
    }
    fn __rmul__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyNimberPoly> {
        Ok(wrap_nimber_poly(self.inner.mul(&parse_nimber_poly(other)?)))
    }
    fn __truediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyNimberPoly> {
        let o = parse_nimber_poly(other)?;
        let oi = o
            .inv()
            .ok_or_else(|| PyValueError::new_err("polynomial divisor is not a unit"))?;
        Ok(wrap_nimber_poly(self.inner.mul(&oi)))
    }
    fn __rtruediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyNimberPoly> {
        let si = self
            .inner
            .inv()
            .ok_or_else(|| PyValueError::new_err("polynomial divisor is not a unit"))?;
        Ok(wrap_nimber_poly(parse_nimber_poly(other)?.mul(&si)))
    }
    fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
        matches!(parse_nimber_poly(other), Ok(p) if p == self.inner)
    }
    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}

pub(crate) fn parse_nimber_poly(obj: &Bound<'_, PyAny>) -> PyResult<Poly<Nimber>> {
    if let Ok(p) = obj.cast::<PyNimberPoly>() {
        return Ok(p.borrow().inner.clone());
    }
    if let Ok(n) = obj.cast::<PyNimber>() {
        return Ok(Poly::constant(n.borrow().inner));
    }
    if let Ok(v) = obj.extract::<u128>() {
        return Ok(Poly::constant(Nimber(v)));
    }
    if let Ok(items) = obj.extract::<Vec<Bound<'_, PyAny>>>() {
        return Ok(Poly::new(parse_nimber_coeffs(items)?));
    }
    Err(PyTypeError::new_err(
        "expected NimberPoly, Nimber, non-negative int, or coefficient list",
    ))
}

pub(crate) fn wrap_nimber_poly(p: Poly<Nimber>) -> PyNimberPoly {
    PyNimberPoly { inner: p }
}

#[pyclass(name = "NimberRationalFunction", module = "pleroma", from_py_object)]
#[derive(Clone)]
pub(crate) struct PyNimberRationalFunction {
    inner: RationalFunction<Nimber>,
}

fn rational_function_from_polys(
    num: Poly<Nimber>,
    den: Poly<Nimber>,
) -> PyResult<RationalFunction<Nimber>> {
    if den.is_zero() {
        return Err(PyValueError::new_err(
            "NimberRationalFunction denominator is zero",
        ));
    }
    Ok(RationalFunction::new(
        num.coeffs().to_vec(),
        den.coeffs().to_vec(),
    ))
}

#[pymethods]
impl PyNimberRationalFunction {
    #[new]
    fn new(num: Vec<Bound<'_, PyAny>>, den: Vec<Bound<'_, PyAny>>) -> PyResult<Self> {
        Ok(wrap_nimber_rational_function(rational_function_from_polys(
            Poly::new(parse_nimber_coeffs(num)?),
            Poly::new(parse_nimber_coeffs(den)?),
        )?))
    }
    #[staticmethod]
    fn zero() -> Self {
        wrap_nimber_rational_function(RationalFunction::zero())
    }
    #[staticmethod]
    fn one() -> Self {
        wrap_nimber_rational_function(RationalFunction::one())
    }
    #[staticmethod]
    fn characteristic() -> u128 {
        <RationalFunction<Nimber> as Scalar>::characteristic()
    }
    #[staticmethod]
    fn t() -> Self {
        wrap_nimber_rational_function(RationalFunction::t())
    }
    #[staticmethod]
    fn from_base(s: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(wrap_nimber_rational_function(RationalFunction::from_base(
            parse_nimber(s)?,
        )))
    }
    #[staticmethod]
    fn from_poly(p: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(wrap_nimber_rational_function(RationalFunction::from_poly(
            parse_nimber_poly(p)?,
        )))
    }
    #[getter]
    fn num(&self) -> PyNimberPoly {
        wrap_nimber_poly(self.inner.num().clone())
    }
    #[getter]
    fn den(&self) -> PyNimberPoly {
        wrap_nimber_poly(self.inner.den().clone())
    }
    fn is_zero(&self) -> bool {
        self.inner.is_zero()
    }
    fn is_integral(&self) -> bool {
        <RationalFunction<Nimber> as HasRingOfIntegers>::is_integral(&self.inner)
    }
    fn to_integer(&self) -> Option<PyNimberPoly> {
        <RationalFunction<Nimber> as HasRingOfIntegers>::to_integer(&self.inner)
            .map(wrap_nimber_poly)
    }
    fn inv(&self) -> PyResult<PyNimberRationalFunction> {
        self.inner
            .inv()
            .map(wrap_nimber_rational_function)
            .ok_or_else(|| PyValueError::new_err("0 has no inverse in F_{2^128}(t)"))
    }
    fn __add__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyNimberRationalFunction> {
        Ok(wrap_nimber_rational_function(
            self.inner.add(&parse_nimber_rational_function(other)?),
        ))
    }
    fn __radd__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyNimberRationalFunction> {
        self.__add__(other)
    }
    fn __sub__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyNimberRationalFunction> {
        Ok(wrap_nimber_rational_function(
            self.inner.sub(&parse_nimber_rational_function(other)?),
        ))
    }
    fn __rsub__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyNimberRationalFunction> {
        Ok(wrap_nimber_rational_function(
            parse_nimber_rational_function(other)?.sub(&self.inner),
        ))
    }
    fn __neg__(&self) -> PyNimberRationalFunction {
        wrap_nimber_rational_function(self.inner.neg())
    }
    fn __mul__(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        match parse_nimber_rational_function(other) {
            Ok(o) => wrap_nimber_rational_function(self.inner.mul(&o)).into_py_any(py),
            Err(_) => Ok(py.NotImplemented()),
        }
    }
    fn __rmul__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyNimberRationalFunction> {
        Ok(wrap_nimber_rational_function(
            self.inner.mul(&parse_nimber_rational_function(other)?),
        ))
    }
    fn __truediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyNimberRationalFunction> {
        let o = parse_nimber_rational_function(other)?;
        let oi = o
            .inv()
            .ok_or_else(|| PyValueError::new_err("division by zero in F_{2^128}(t)"))?;
        Ok(wrap_nimber_rational_function(self.inner.mul(&oi)))
    }
    fn __rtruediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyNimberRationalFunction> {
        let si = self
            .inner
            .inv()
            .ok_or_else(|| PyValueError::new_err("division by zero in F_{2^128}(t)"))?;
        Ok(wrap_nimber_rational_function(
            parse_nimber_rational_function(other)?.mul(&si),
        ))
    }
    fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
        matches!(parse_nimber_rational_function(other), Ok(f) if f == self.inner)
    }
    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}

pub(crate) fn parse_nimber_rational_function(
    obj: &Bound<'_, PyAny>,
) -> PyResult<RationalFunction<Nimber>> {
    if let Ok(f) = obj.cast::<PyNimberRationalFunction>() {
        return Ok(f.borrow().inner.clone());
    }
    if let Ok(p) = obj.cast::<PyNimberPoly>() {
        return Ok(RationalFunction::from_poly(p.borrow().inner.clone()));
    }
    if let Ok(n) = obj.cast::<PyNimber>() {
        return Ok(RationalFunction::from_base(n.borrow().inner));
    }
    if let Ok(v) = obj.extract::<u128>() {
        return Ok(RationalFunction::from_base(Nimber(v)));
    }
    if let Ok(items) = obj.extract::<Vec<Bound<'_, PyAny>>>() {
        return Ok(RationalFunction::from_poly(Poly::new(parse_nimber_coeffs(
            items,
        )?)));
    }
    Err(PyTypeError::new_err(
        "expected NimberRationalFunction, NimberPoly, Nimber, non-negative int, or numerator coefficient list",
    ))
}

pub(crate) fn wrap_nimber_rational_function(
    f: RationalFunction<Nimber>,
) -> PyNimberRationalFunction {
    PyNimberRationalFunction { inner: f }
}

macro_rules! prime_field_pyclass {
    ($py:ident, $name:literal, $parse:ident, $wrap:ident, $p:literal) => {
        #[pyclass(name = $name, module = "pleroma", from_py_object)]
        #[derive(Clone)]
        pub(crate) struct $py {
            inner: Fp<$p>,
        }

        #[pymethods]
        impl $py {
            #[new]
            fn new(value: i128) -> Self {
                $wrap(Fp::<$p>::new(value))
            }
            #[staticmethod]
            fn modulus() -> u128 {
                $p
            }
            #[staticmethod]
            fn zero() -> Self {
                $wrap(Fp::<$p>::zero())
            }
            #[staticmethod]
            fn one() -> Self {
                $wrap(Fp::<$p>::one())
            }
            #[staticmethod]
            fn characteristic() -> u128 {
                Fp::<$p>::characteristic()
            }
            #[staticmethod]
            fn modulus_is_prime() -> bool {
                Fp::<$p>::modulus_is_prime()
            }
            #[staticmethod]
            fn assert_prime_modulus() {
                Fp::<$p>::assert_prime_modulus()
            }
            #[staticmethod]
            fn from_u128(value: u128) -> Self {
                $wrap(Fp::<$p>::from_u128(value))
            }
            #[getter]
            fn value(&self) -> u128 {
                self.inner.value()
            }
            fn is_zero(&self) -> bool {
                self.inner.is_zero()
            }
            fn is_square(&self) -> bool {
                ExactRoots::is_square(&self.inner)
            }
            fn sqrt(&self) -> Option<Self> {
                ExactRoots::sqrt(&self.inner).map($wrap)
            }
            fn pow(&self, e: u128) -> Self {
                $wrap(self.inner.pow(e))
            }
            fn __pow__(&self, e: u128, modulo: Option<u128>) -> PyResult<Self> {
                if modulo.is_some() {
                    return Err(PyValueError::new_err(
                        "finite-field exponentiation does not take a modulus",
                    ));
                }
                Ok(self.pow(e))
            }
            fn inv(&self) -> PyResult<Self> {
                self.inner
                    .inv()
                    .map($wrap)
                    .ok_or_else(|| PyValueError::new_err("0 has no inverse in this field"))
            }
            fn __add__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                Ok($wrap(self.inner.add(&$parse(other)?)))
            }
            fn __radd__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                self.__add__(other)
            }
            fn __sub__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                Ok($wrap(self.inner.sub(&$parse(other)?)))
            }
            fn __rsub__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                Ok($wrap($parse(other)?.sub(&self.inner)))
            }
            fn __neg__(&self) -> Self {
                $wrap(self.inner.neg())
            }
            fn __mul__(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
                match $parse(other) {
                    Ok(o) => $wrap(self.inner.mul(&o)).into_py_any(py),
                    Err(_) => Ok(py.NotImplemented()),
                }
            }
            fn __rmul__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                Ok($wrap(self.inner.mul(&$parse(other)?)))
            }
            fn __truediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                let o = $parse(other)?;
                let oi = o
                    .inv()
                    .ok_or_else(|| PyValueError::new_err("division by 0 in this field"))?;
                Ok($wrap(self.inner.mul(&oi)))
            }
            fn __rtruediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                let si = self
                    .inner
                    .inv()
                    .ok_or_else(|| PyValueError::new_err("division by 0 in this field"))?;
                Ok($wrap($parse(other)?.mul(&si)))
            }
            fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
                matches!(other.cast::<$py>(), Ok(x) if x.borrow().inner == self.inner)
            }
            fn __hash__(&self) -> usize {
                self.inner.value() as usize
            }
            fn __repr__(&self) -> String {
                format!("{:?}", self.inner)
            }
        }

        pub(crate) fn $parse(obj: &Bound<'_, PyAny>) -> PyResult<Fp<$p>> {
            if let Ok(x) = obj.cast::<$py>() {
                return Ok(x.borrow().inner);
            }
            if let Ok(v) = obj.extract::<i128>() {
                return Ok(Fp::<$p>::new(v));
            }
            Err(PyTypeError::new_err(concat!(
                "expected ",
                $name,
                " or int"
            )))
        }

        pub(crate) fn $wrap(inner: Fp<$p>) -> $py {
            $py { inner }
        }
    };
}

macro_rules! extension_field_pyclass {
    ($py:ident, $name:literal, $parse:ident, $wrap:ident, $p:literal, $n:literal) => {
        #[pyclass(name = $name, module = "pleroma", from_py_object)]
        #[derive(Clone)]
        pub(crate) struct $py {
            inner: Fpn<$p, $n>,
        }

        pub(crate) fn $parse(obj: &Bound<'_, PyAny>) -> PyResult<Fpn<$p, $n>> {
            if let Ok(x) = obj.cast::<$py>() {
                return Ok(x.borrow().inner);
            }
            if let Ok(cs) = obj.extract::<Vec<u128>>() {
                return parse_coeffs::<$p, $n>(&cs).map(|cs| Fpn::<$p, $n>::from_coeffs(&cs));
            }
            if let Ok(v) = obj.extract::<u128>() {
                return Ok(Fpn::<$p, $n>::constant(v));
            }
            Err(PyTypeError::new_err(concat!(
                "expected ",
                $name,
                ", non-negative int, or coefficient list"
            )))
        }

        pub(crate) fn $wrap(inner: Fpn<$p, $n>) -> $py {
            $py { inner }
        }

        #[pymethods]
        impl $py {
            #[staticmethod]
            fn from_coeffs(coeffs: Vec<u128>) -> PyResult<Self> {
                parse_coeffs::<$p, $n>(&coeffs)
                    .map(|cs| Fpn::<$p, $n>::from_coeffs(&cs))
                    .map($wrap)
            }
            #[staticmethod]
            fn constant(c: u128) -> Self {
                $wrap(Fpn::<$p, $n>::constant(c))
            }
            #[staticmethod]
            fn generator() -> Self {
                $wrap(Fpn::<$p, $n>::generator())
            }
            #[staticmethod]
            fn primitive_element() -> Self {
                $wrap(Fpn::<$p, $n>::primitive_element())
            }
            #[staticmethod]
            fn zero() -> Self {
                $wrap(Fpn::<$p, $n>::zero())
            }
            #[staticmethod]
            fn one() -> Self {
                $wrap(Fpn::<$p, $n>::one())
            }
            #[staticmethod]
            fn characteristic() -> u128 {
                Fpn::<$p, $n>::characteristic()
            }
            #[staticmethod]
            fn is_supported_field() -> bool {
                Fpn::<$p, $n>::is_supported_field()
            }
            #[staticmethod]
            fn assert_supported_field() {
                Fpn::<$p, $n>::assert_supported_field()
            }
            #[staticmethod]
            fn from_index(code: u128) -> PyResult<Self> {
                if code >= Fpn::<$p, $n>::field_order() {
                    return Err(PyValueError::new_err(format!(
                        "field element index {code} is outside F_{}",
                        Fpn::<$p, $n>::field_order()
                    )));
                }
                let mut coeffs = Vec::with_capacity($n);
                let mut x = code;
                for _ in 0..$n {
                    coeffs.push(x % $p);
                    x /= $p;
                }
                Ok($wrap(Fpn::<$p, $n>::from_coeffs(&coeffs)))
            }
            #[staticmethod]
            fn field_order() -> u128 {
                Fpn::<$p, $n>::field_order()
            }
            #[staticmethod]
            fn ext_degree() -> usize {
                <Fpn<$p, $n> as FiniteField>::ext_degree()
            }
            #[staticmethod]
            fn group_order() -> u128 {
                <Fpn<$p, $n> as FiniteField>::group_order()
            }
            #[staticmethod]
            fn group_order_factors() -> Vec<u128> {
                <Fpn<$p, $n> as FiniteField>::group_order_factors()
            }
            #[staticmethod]
            fn reduction_rule() -> Vec<u128> {
                Fpn::<$p, $n>::reduction_rule().to_vec()
            }
            #[staticmethod]
            fn reduction_polynomial_kind() -> PyReductionPolynomialKind {
                wrap_reduction_polynomial_kind(Fpn::<$p, $n>::reduction_polynomial_kind())
            }
            #[staticmethod]
            fn is_conway_polynomial() -> bool {
                Fpn::<$p, $n>::is_conway_polynomial()
            }
            #[getter]
            fn coeffs(&self) -> Vec<u128> {
                self.inner.coeffs().to_vec()
            }
            // Mirrors Rust's consuming `Fpn::into_coeffs` under Python's borrowed method model.
            #[allow(clippy::wrong_self_convention)]
            fn into_coeffs(&self) -> Vec<u128> {
                self.inner.into_coeffs().to_vec()
            }
            fn coeff(&self, i: usize) -> u128 {
                self.inner.coeff(i)
            }
            fn is_zero(&self) -> bool {
                self.inner.is_zero()
            }
            fn is_square(&self) -> bool {
                self.inner.is_square()
            }
            fn sqrt(&self) -> Option<Self> {
                ExactRoots::sqrt(&self.inner).map($wrap)
            }
            fn degree(&self) -> usize {
                self.inner.degree()
            }
            fn conjugates(&self) -> Vec<Self> {
                self.inner.conjugates().into_iter().map($wrap).collect()
            }
            fn min_poly(&self) -> Vec<u128> {
                self.inner.min_poly()
            }
            fn min_poly_monic(&self) -> Vec<Self> {
                self.inner.min_poly_monic().into_iter().map($wrap).collect()
            }
            fn frobenius(&self) -> Self {
                $wrap(self.inner.frobenius())
            }
            fn frobenius_iter(&self, k: usize) -> Self {
                $wrap(self.inner.frobenius_iter(k))
            }
            fn relative_trace_over(&self, m: usize, e: usize) -> PyResult<Self> {
                validate_relative_degrees(&self.inner, m, e)?;
                Ok($wrap(self.inner.relative_trace_over(m, e)))
            }
            fn relative_norm_over(&self, m: usize, e: usize) -> PyResult<Self> {
                validate_relative_degrees(&self.inner, m, e)?;
                Ok($wrap(self.inner.relative_norm_over(m, e)))
            }
            fn relative_trace(&self, e: usize) -> PyResult<Self> {
                validate_relative_degrees(
                    &self.inner,
                    <Fpn<$p, $n> as FiniteField>::ext_degree(),
                    e,
                )?;
                Ok($wrap(self.inner.relative_trace(e)))
            }
            fn relative_norm(&self, e: usize) -> PyResult<Self> {
                validate_relative_degrees(
                    &self.inner,
                    <Fpn<$p, $n> as FiniteField>::ext_degree(),
                    e,
                )?;
                Ok($wrap(self.inner.relative_norm(e)))
            }
            fn absolute_trace(&self) -> u128 {
                self.inner.relative_trace(1).coeff(0)
            }
            fn absolute_norm(&self) -> u128 {
                self.inner.relative_norm(1).coeff(0)
            }
            fn multiplicative_order(&self) -> Option<u128> {
                self.inner.multiplicative_order()
            }
            fn is_primitive(&self) -> bool {
                self.inner.is_primitive()
            }
            fn discrete_log(&self, x: &Bound<'_, PyAny>) -> PyResult<Option<u128>> {
                Ok(self.inner.discrete_log($parse(x)?))
            }
            fn pow(&self, e: u128) -> Self {
                $wrap(self.inner.pow(e))
            }
            fn __pow__(&self, e: u128, modulo: Option<u128>) -> PyResult<Self> {
                if modulo.is_some() {
                    return Err(PyValueError::new_err(
                        "finite-field exponentiation does not take a modulus",
                    ));
                }
                Ok(self.pow(e))
            }
            fn inv(&self) -> PyResult<Self> {
                self.inner
                    .inv()
                    .map($wrap)
                    .ok_or_else(|| PyValueError::new_err("0 has no inverse in this field"))
            }
            fn __add__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                Ok($wrap(self.inner.add(&$parse(other)?)))
            }
            fn __radd__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                self.__add__(other)
            }
            fn __sub__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                Ok($wrap(self.inner.sub(&$parse(other)?)))
            }
            fn __rsub__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                Ok($wrap($parse(other)?.sub(&self.inner)))
            }
            fn __neg__(&self) -> Self {
                $wrap(self.inner.neg())
            }
            fn __mul__(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
                match $parse(other) {
                    Ok(o) => $wrap(self.inner.mul(&o)).into_py_any(py),
                    Err(_) => Ok(py.NotImplemented()),
                }
            }
            fn __rmul__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                Ok($wrap(self.inner.mul(&$parse(other)?)))
            }
            fn __truediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                let o = $parse(other)?;
                let oi = o
                    .inv()
                    .ok_or_else(|| PyValueError::new_err("division by 0 in this field"))?;
                Ok($wrap(self.inner.mul(&oi)))
            }
            fn __rtruediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                let si = self
                    .inner
                    .inv()
                    .ok_or_else(|| PyValueError::new_err("division by 0 in this field"))?;
                Ok($wrap($parse(other)?.mul(&si)))
            }
            fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
                matches!(other.cast::<$py>(), Ok(x) if x.borrow().inner == self.inner)
            }
            fn __repr__(&self) -> String {
                format!("{:?}", self.inner)
            }
        }
    };
}

#[pyclass(name = "ReductionPolynomialKind", module = "pleroma", from_py_object)]
#[derive(Clone, Copy)]
struct PyReductionPolynomialKind {
    inner: ReductionPolynomialKind,
}

fn reduction_polynomial_kind_name(kind: ReductionPolynomialKind) -> &'static str {
    match kind {
        ReductionPolynomialKind::PrimeField => "PrimeField",
        ReductionPolynomialKind::Conway => "Conway",
        ReductionPolynomialKind::Irreducible => "Irreducible",
    }
}

fn wrap_reduction_polynomial_kind(inner: ReductionPolynomialKind) -> PyReductionPolynomialKind {
    PyReductionPolynomialKind { inner }
}

fn parse_reduction_polynomial_kind(obj: &Bound<'_, PyAny>) -> PyResult<ReductionPolynomialKind> {
    if let Ok(kind) = obj.cast::<PyReductionPolynomialKind>() {
        return Ok(kind.borrow().inner);
    }
    Err(PyTypeError::new_err("expected ReductionPolynomialKind"))
}

#[pymethods]
impl PyReductionPolynomialKind {
    #[staticmethod]
    fn prime_field() -> Self {
        wrap_reduction_polynomial_kind(ReductionPolynomialKind::PrimeField)
    }

    #[staticmethod]
    fn conway() -> Self {
        wrap_reduction_polynomial_kind(ReductionPolynomialKind::Conway)
    }

    #[staticmethod]
    fn irreducible() -> Self {
        wrap_reduction_polynomial_kind(ReductionPolynomialKind::Irreducible)
    }

    #[getter]
    fn name(&self) -> &'static str {
        reduction_polynomial_kind_name(self.inner)
    }

    #[getter]
    fn is_prime_field(&self) -> bool {
        self.inner == ReductionPolynomialKind::PrimeField
    }

    #[getter]
    fn is_conway(&self) -> bool {
        self.inner == ReductionPolynomialKind::Conway
    }

    #[getter]
    fn is_irreducible(&self) -> bool {
        self.inner == ReductionPolynomialKind::Irreducible
    }

    fn __str__(&self) -> &'static str {
        reduction_polynomial_kind_name(self.inner)
    }

    fn __repr__(&self) -> String {
        format!(
            "ReductionPolynomialKind.{}",
            reduction_polynomial_kind_name(self.inner)
        )
    }

    fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
        matches!(parse_reduction_polynomial_kind(other), Ok(kind) if kind == self.inner)
    }
}

fn parse_coeffs<const P: u128, const N: usize>(coeffs: &[u128]) -> PyResult<Vec<u128>> {
    if coeffs.iter().skip(N).any(|&c| c % P != 0) {
        return Err(PyValueError::new_err(format!(
            "coefficient list has nonzero terms beyond degree {}",
            N - 1
        )));
    }
    Ok(coeffs.to_vec())
}

prime_field_pyclass!(PyFp2, "Fp2", parse_fp2, wrap_fp2, 2);
prime_field_pyclass!(PyFp3, "Fp3", parse_fp3, wrap_fp3, 3);
prime_field_pyclass!(PyFp5, "Fp5", parse_fp5, wrap_fp5, 5);
prime_field_pyclass!(PyFp7, "Fp7", parse_fp7, wrap_fp7, 7);
prime_field_pyclass!(PyFp11, "Fp11", parse_fp11, wrap_fp11, 11);
prime_field_pyclass!(PyFp13, "Fp13", parse_fp13, wrap_fp13, 13);

extension_field_pyclass!(PyF4, "F4", parse_f4, wrap_f4, 2, 2);
extension_field_pyclass!(PyF8, "F8", parse_f8, wrap_f8, 2, 3);
extension_field_pyclass!(PyF16, "F16", parse_f16, wrap_f16, 2, 4);
extension_field_pyclass!(PyF9, "F9", parse_f9, wrap_f9, 3, 2);
extension_field_pyclass!(PyF25, "F25", parse_f25, wrap_f25, 5, 2);
extension_field_pyclass!(PyF27, "F27", parse_f27, wrap_f27, 3, 3);

macro_rules! function_field_pyclasses {
    (
        $poly_py:ident, $poly_name:literal, $parse_poly:ident, $wrap_poly:ident,
        $rf_py:ident, $rf_name:literal, $parse_rf:ident, $wrap_rf:ident,
        $base:ty, $base_py:ty, $base_parse:path, $base_wrap:path
    ) => {
        #[pyclass(name = $poly_name, module = "pleroma", from_py_object)]
        #[derive(Clone)]
        pub(crate) struct $poly_py {
            inner: Poly<$base>,
        }

        pub(crate) fn $parse_poly(obj: &Bound<'_, PyAny>) -> PyResult<Poly<$base>> {
            if let Ok(p) = obj.cast::<$poly_py>() {
                return Ok(p.borrow().inner.clone());
            }
            if let Ok(s) = $base_parse(obj) {
                return Ok(Poly::constant(s));
            }
            if let Ok(items) = obj.extract::<Vec<Bound<'_, PyAny>>>() {
                return items
                    .iter()
                    .map($base_parse)
                    .collect::<PyResult<Vec<_>>>()
                    .map(Poly::new);
            }
            Err(PyTypeError::new_err(concat!(
                "expected ",
                $poly_name,
                ", base scalar, int, or coefficient list"
            )))
        }

        pub(crate) fn $wrap_poly(inner: Poly<$base>) -> $poly_py {
            $poly_py { inner }
        }

        #[pymethods]
        impl $poly_py {
            #[new]
            fn new(coeffs: Vec<Bound<'_, PyAny>>) -> PyResult<Self> {
                coeffs
                    .iter()
                    .map($base_parse)
                    .collect::<PyResult<Vec<_>>>()
                    .map(Poly::new)
                    .map($wrap_poly)
            }
            #[staticmethod]
            fn zero() -> Self {
                $wrap_poly(Poly::zero())
            }
            #[staticmethod]
            fn one() -> Self {
                $wrap_poly(Poly::one())
            }
            #[staticmethod]
            fn characteristic() -> u128 {
                <Poly<$base> as Scalar>::characteristic()
            }
            #[staticmethod]
            fn x() -> Self {
                $wrap_poly(Poly::x())
            }
            #[staticmethod]
            fn constant(s: &Bound<'_, PyAny>) -> PyResult<Self> {
                Ok($wrap_poly(Poly::constant($base_parse(s)?)))
            }
            #[staticmethod]
            fn monomial(deg: usize, coeff: &Bound<'_, PyAny>) -> PyResult<Self> {
                Ok($wrap_poly(Poly::monomial(deg, $base_parse(coeff)?)))
            }
            #[getter]
            fn coeffs(&self) -> Vec<$base_py> {
                self.inner
                    .coeffs()
                    .iter()
                    .cloned()
                    .map($base_wrap)
                    .collect()
            }
            #[getter]
            fn degree(&self) -> Option<usize> {
                self.inner.degree()
            }
            fn leading(&self) -> Option<$base_py> {
                self.inner.leading().cloned().map($base_wrap)
            }
            fn coeff(&self, i: usize) -> $base_py {
                $base_wrap(self.inner.coeff(i))
            }
            fn is_zero(&self) -> bool {
                self.inner.is_zero()
            }
            fn eval(&self, x: &Bound<'_, PyAny>) -> PyResult<$base_py> {
                Ok($base_wrap(self.inner.eval(&$base_parse(x)?)))
            }
            fn scale(&self, s: &Bound<'_, PyAny>) -> PyResult<Self> {
                Ok($wrap_poly(self.inner.scale(&$base_parse(s)?)))
            }
            fn make_monic(&self) -> PyResult<Self> {
                if self.inner.is_zero() {
                    return Err(PyValueError::new_err("make_monic of the zero polynomial"));
                }
                Ok($wrap_poly(self.inner.make_monic()))
            }
            fn divrem(&self, divisor: &$poly_py) -> PyResult<($poly_py, $poly_py)> {
                if divisor.inner.is_zero() {
                    return Err(PyValueError::new_err("polynomial division by zero"));
                }
                let (q, r) = self.inner.divrem(&divisor.inner);
                Ok(($wrap_poly(q), $wrap_poly(r)))
            }
            fn rem(&self, divisor: &$poly_py) -> PyResult<$poly_py> {
                if divisor.inner.is_zero() {
                    return Err(PyValueError::new_err("polynomial remainder by zero"));
                }
                Ok($wrap_poly(self.inner.rem(&divisor.inner)))
            }
            fn divides(&self, multiple: &$poly_py) -> bool {
                self.inner.divides(&multiple.inner)
            }
            fn gcd(&self, other: &$poly_py) -> $poly_py {
                $wrap_poly(self.inner.gcd(&other.inner))
            }
            fn mul_mod(&self, other: &$poly_py, modulus: &$poly_py) -> PyResult<$poly_py> {
                if modulus.inner.is_zero() {
                    return Err(PyValueError::new_err("polynomial modulus is zero"));
                }
                Ok($wrap_poly(self.inner.mul_mod(&other.inner, &modulus.inner)))
            }
            fn pow_mod(&self, e: u128, modulus: &$poly_py) -> PyResult<$poly_py> {
                if modulus.inner.is_zero() {
                    return Err(PyValueError::new_err("polynomial modulus is zero"));
                }
                Ok($wrap_poly(self.inner.pow_mod(e, &modulus.inner)))
            }
            fn inv(&self) -> PyResult<$poly_py> {
                self.inner
                    .inv()
                    .map($wrap_poly)
                    .ok_or_else(|| PyValueError::new_err("only nonzero constant polynomials invert"))
            }
            fn to_fraction(&self) -> $rf_py {
                $wrap_rf(<Poly<$base> as HasFractionField>::to_fraction(&self.inner))
            }
            fn __add__(&self, other: &Bound<'_, PyAny>) -> PyResult<$poly_py> {
                Ok($wrap_poly(self.inner.add(&$parse_poly(other)?)))
            }
            fn __radd__(&self, other: &Bound<'_, PyAny>) -> PyResult<$poly_py> {
                self.__add__(other)
            }
            fn __sub__(&self, other: &Bound<'_, PyAny>) -> PyResult<$poly_py> {
                Ok($wrap_poly(self.inner.sub(&$parse_poly(other)?)))
            }
            fn __rsub__(&self, other: &Bound<'_, PyAny>) -> PyResult<$poly_py> {
                Ok($wrap_poly($parse_poly(other)?.sub(&self.inner)))
            }
            fn __neg__(&self) -> $poly_py {
                $wrap_poly(self.inner.neg())
            }
            fn __mul__(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
                match $parse_poly(other) {
                    Ok(o) => $wrap_poly(self.inner.mul(&o)).into_py_any(py),
                    Err(_) => Ok(py.NotImplemented()),
                }
            }
            fn __rmul__(&self, other: &Bound<'_, PyAny>) -> PyResult<$poly_py> {
                Ok($wrap_poly(self.inner.mul(&$parse_poly(other)?)))
            }
            fn __truediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<$poly_py> {
                let o = $parse_poly(other)?;
                let oi = o
                    .inv()
                    .ok_or_else(|| PyValueError::new_err("polynomial divisor is not a unit"))?;
                Ok($wrap_poly(self.inner.mul(&oi)))
            }
            fn __rtruediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<$poly_py> {
                let si = self
                    .inner
                    .inv()
                    .ok_or_else(|| PyValueError::new_err("polynomial divisor is not a unit"))?;
                Ok($wrap_poly($parse_poly(other)?.mul(&si)))
            }
            fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
                matches!($parse_poly(other), Ok(p) if p == self.inner)
            }
            fn __repr__(&self) -> String {
                format!("{:?}", self.inner)
            }
        }

        #[pyclass(name = $rf_name, module = "pleroma", from_py_object)]
        #[derive(Clone)]
        pub(crate) struct $rf_py {
            inner: RationalFunction<$base>,
        }

        pub(crate) fn $parse_rf(obj: &Bound<'_, PyAny>) -> PyResult<RationalFunction<$base>> {
            if let Ok(f) = obj.cast::<$rf_py>() {
                return Ok(f.borrow().inner.clone());
            }
            if let Ok(p) = obj.cast::<$poly_py>() {
                return Ok(RationalFunction::from_poly(p.borrow().inner.clone()));
            }
            if let Ok(s) = $base_parse(obj) {
                return Ok(RationalFunction::from_base(s));
            }
            if let Ok(items) = obj.extract::<Vec<Bound<'_, PyAny>>>() {
                return items
                    .iter()
                    .map($base_parse)
                    .collect::<PyResult<Vec<_>>>()
                    .map(Poly::new)
                    .map(RationalFunction::from_poly);
            }
            Err(PyTypeError::new_err(concat!(
                "expected ",
                $rf_name,
                ", ",
                $poly_name,
                ", base scalar, int, or numerator coefficient list"
            )))
        }

        pub(crate) fn $wrap_rf(inner: RationalFunction<$base>) -> $rf_py {
            $rf_py { inner }
        }

        #[pymethods]
        impl $rf_py {
            #[new]
            fn new(num: Vec<Bound<'_, PyAny>>, den: Vec<Bound<'_, PyAny>>) -> PyResult<Self> {
                let num = num
                    .iter()
                    .map($base_parse)
                    .collect::<PyResult<Vec<_>>>()
                    .map(Poly::new)?;
                let den = den
                    .iter()
                    .map($base_parse)
                    .collect::<PyResult<Vec<_>>>()
                    .map(Poly::new)?;
                if den.is_zero() {
                    return Err(PyValueError::new_err(concat!(
                        $rf_name,
                        " denominator is zero"
                    )));
                }
                Ok($wrap_rf(RationalFunction::new(
                    num.coeffs().to_vec(),
                    den.coeffs().to_vec(),
                )))
            }
            #[staticmethod]
            fn zero() -> Self {
                $wrap_rf(RationalFunction::zero())
            }
            #[staticmethod]
            fn one() -> Self {
                $wrap_rf(RationalFunction::one())
            }
            #[staticmethod]
            fn characteristic() -> u128 {
                <RationalFunction<$base> as Scalar>::characteristic()
            }
            #[staticmethod]
            fn t() -> Self {
                $wrap_rf(RationalFunction::t())
            }
            #[staticmethod]
            fn from_base(s: &Bound<'_, PyAny>) -> PyResult<Self> {
                Ok($wrap_rf(RationalFunction::from_base($base_parse(s)?)))
            }
            #[staticmethod]
            fn from_poly(p: &Bound<'_, PyAny>) -> PyResult<Self> {
                Ok($wrap_rf(RationalFunction::from_poly($parse_poly(p)?)))
            }
            #[getter]
            fn num(&self) -> $poly_py {
                $wrap_poly(self.inner.num().clone())
            }
            #[getter]
            fn den(&self) -> $poly_py {
                $wrap_poly(self.inner.den().clone())
            }
            fn is_zero(&self) -> bool {
                self.inner.is_zero()
            }
            fn is_integral(&self) -> bool {
                <RationalFunction<$base> as HasRingOfIntegers>::is_integral(&self.inner)
            }
            fn to_integer(&self) -> Option<$poly_py> {
                <RationalFunction<$base> as HasRingOfIntegers>::to_integer(&self.inner)
                    .map($wrap_poly)
            }
            fn inv(&self) -> PyResult<$rf_py> {
                self.inner
                    .inv()
                    .map($wrap_rf)
                    .ok_or_else(|| PyValueError::new_err(concat!("0 has no inverse in ", $rf_name)))
            }
            fn __add__(&self, other: &Bound<'_, PyAny>) -> PyResult<$rf_py> {
                Ok($wrap_rf(self.inner.add(&$parse_rf(other)?)))
            }
            fn __radd__(&self, other: &Bound<'_, PyAny>) -> PyResult<$rf_py> {
                self.__add__(other)
            }
            fn __sub__(&self, other: &Bound<'_, PyAny>) -> PyResult<$rf_py> {
                Ok($wrap_rf(self.inner.sub(&$parse_rf(other)?)))
            }
            fn __rsub__(&self, other: &Bound<'_, PyAny>) -> PyResult<$rf_py> {
                Ok($wrap_rf($parse_rf(other)?.sub(&self.inner)))
            }
            fn __neg__(&self) -> $rf_py {
                $wrap_rf(self.inner.neg())
            }
            fn __mul__(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
                match $parse_rf(other) {
                    Ok(o) => $wrap_rf(self.inner.mul(&o)).into_py_any(py),
                    Err(_) => Ok(py.NotImplemented()),
                }
            }
            fn __rmul__(&self, other: &Bound<'_, PyAny>) -> PyResult<$rf_py> {
                Ok($wrap_rf(self.inner.mul(&$parse_rf(other)?)))
            }
            fn __truediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<$rf_py> {
                let o = $parse_rf(other)?;
                let oi = o
                    .inv()
                    .ok_or_else(|| PyValueError::new_err(concat!("division by zero in ", $rf_name)))?;
                Ok($wrap_rf(self.inner.mul(&oi)))
            }
            fn __rtruediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<$rf_py> {
                let si = self
                    .inner
                    .inv()
                    .ok_or_else(|| PyValueError::new_err(concat!("division by zero in ", $rf_name)))?;
                Ok($wrap_rf($parse_rf(other)?.mul(&si)))
            }
            fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
                matches!($parse_rf(other), Ok(f) if f == self.inner)
            }
            fn __repr__(&self) -> String {
                format!("{:?}", self.inner)
            }
        }
    };
}

function_field_pyclasses!(
    PyFp2Poly,
    "Fp2Poly",
    parse_fp2_poly,
    wrap_fp2_poly,
    PyFp2RationalFunction,
    "Fp2RationalFunction",
    parse_fp2_rational_function,
    wrap_fp2_rational_function,
    Fp<2>,
    PyFp2,
    parse_fp2,
    wrap_fp2
);
function_field_pyclasses!(
    PyFp3Poly,
    "Fp3Poly",
    parse_fp3_poly,
    wrap_fp3_poly,
    PyFp3RationalFunction,
    "Fp3RationalFunction",
    parse_fp3_rational_function,
    wrap_fp3_rational_function,
    Fp<3>,
    PyFp3,
    parse_fp3,
    wrap_fp3
);
function_field_pyclasses!(
    PyFp5Poly,
    "Fp5Poly",
    parse_fp5_poly,
    wrap_fp5_poly,
    PyFp5RationalFunction,
    "Fp5RationalFunction",
    parse_fp5_rational_function,
    wrap_fp5_rational_function,
    Fp<5>,
    PyFp5,
    parse_fp5,
    wrap_fp5
);
function_field_pyclasses!(
    PyFp7Poly,
    "Fp7Poly",
    parse_fp7_poly,
    wrap_fp7_poly,
    PyFp7RationalFunction,
    "Fp7RationalFunction",
    parse_fp7_rational_function,
    wrap_fp7_rational_function,
    Fp<7>,
    PyFp7,
    parse_fp7,
    wrap_fp7
);
function_field_pyclasses!(
    PyFp11Poly,
    "Fp11Poly",
    parse_fp11_poly,
    wrap_fp11_poly,
    PyFp11RationalFunction,
    "Fp11RationalFunction",
    parse_fp11_rational_function,
    wrap_fp11_rational_function,
    Fp<11>,
    PyFp11,
    parse_fp11,
    wrap_fp11
);
function_field_pyclasses!(
    PyFp13Poly,
    "Fp13Poly",
    parse_fp13_poly,
    wrap_fp13_poly,
    PyFp13RationalFunction,
    "Fp13RationalFunction",
    parse_fp13_rational_function,
    wrap_fp13_rational_function,
    Fp<13>,
    PyFp13,
    parse_fp13,
    wrap_fp13
);

macro_rules! zp_pyclass {
    ($py:ident, $name:literal, $parse:ident, $wrap:ident, $p:literal, $k:literal) => {
        #[pyclass(name = $name, module = "pleroma", from_py_object)]
        #[derive(Clone)]
        pub(crate) struct $py {
            inner: Zp<$p, $k>,
        }

        #[pymethods]
        impl $py {
            #[new]
            fn new(value: i128) -> Self {
                $wrap(Zp::<$p, $k>::new(value))
            }
            #[staticmethod]
            fn zero() -> Self {
                $wrap(Zp::<$p, $k>::zero())
            }
            #[staticmethod]
            fn one() -> Self {
                $wrap(Zp::<$p, $k>::one())
            }
            #[staticmethod]
            fn prime() -> u128 {
                $p
            }
            #[staticmethod]
            fn precision() -> u128 {
                $k
            }
            #[staticmethod]
            fn modulus() -> u128 {
                Zp::<$p, $k>::modulus()
            }
            #[staticmethod]
            fn assert_supported_ring() {
                Zp::<$p, $k>::assert_supported_ring()
            }
            #[staticmethod]
            fn characteristic() -> u128 {
                Zp::<$p, $k>::characteristic()
            }
            #[getter]
            fn value(&self) -> u128 {
                self.inner.0
            }
            fn valuation(&self) -> u128 {
                self.inner.valuation()
            }
            fn is_unit(&self) -> bool {
                self.inner.is_unit()
            }
            fn is_zero(&self) -> bool {
                self.inner.is_zero()
            }
            fn is_square(&self) -> PyResult<bool> {
                self.inner.is_square().ok_or_else(|| {
                    PyValueError::new_err("squarehood is undecidable at this p-adic precision")
                })
            }
            fn sqrt(&self) -> PyResult<Option<Self>> {
                self.inner
                    .sqrt()
                    .map(|root| root.map($wrap))
                    .ok_or_else(|| {
                        PyValueError::new_err(
                            "square-root construction is not implemented for this p-adic case",
                        )
                    })
            }
            fn inv(&self) -> PyResult<Self> {
                self.inner
                    .inv()
                    .map($wrap)
                    .ok_or_else(|| PyValueError::new_err("non-unit has no inverse in Z/p^k"))
            }
            fn __add__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                Ok($wrap(self.inner.add(&$parse(other)?)))
            }
            fn __radd__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                self.__add__(other)
            }
            fn __sub__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                Ok($wrap(self.inner.sub(&$parse(other)?)))
            }
            fn __rsub__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                Ok($wrap($parse(other)?.sub(&self.inner)))
            }
            fn __neg__(&self) -> Self {
                $wrap(self.inner.neg())
            }
            fn __mul__(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
                match $parse(other) {
                    Ok(o) => $wrap(self.inner.mul(&o)).into_py_any(py),
                    Err(_) => Ok(py.NotImplemented()),
                }
            }
            fn __rmul__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                Ok($wrap(self.inner.mul(&$parse(other)?)))
            }
            fn __truediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                let o = $parse(other)?;
                let oi = o
                    .inv()
                    .ok_or_else(|| PyValueError::new_err("divisor is not a unit in Z/p^k"))?;
                Ok($wrap(self.inner.mul(&oi)))
            }
            fn __rtruediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                let si = self
                    .inner
                    .inv()
                    .ok_or_else(|| PyValueError::new_err("divisor is not a unit in Z/p^k"))?;
                Ok($wrap($parse(other)?.mul(&si)))
            }
            fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
                matches!(other.cast::<$py>(), Ok(x) if x.borrow().inner == self.inner)
            }
            fn __hash__(&self) -> usize {
                self.inner.0 as usize
            }
            fn __repr__(&self) -> String {
                format!("{:?}", self.inner)
            }
        }

        pub(crate) fn $parse(obj: &Bound<'_, PyAny>) -> PyResult<Zp<$p, $k>> {
            if let Ok(x) = obj.cast::<$py>() {
                return Ok(x.borrow().inner);
            }
            if let Ok(v) = obj.extract::<i128>() {
                return Ok(Zp::<$p, $k>::new(v));
            }
            Err(PyTypeError::new_err(concat!("expected ", $name, " or int")))
        }

        pub(crate) fn $wrap(inner: Zp<$p, $k>) -> $py {
            $py { inner }
        }
    };
}

macro_rules! qp_pyclass {
    (
        $py:ident, $name:literal, $parse:ident, $wrap:ident, $p:literal, $k:literal,
        $int_py:ty, $int_wrap:path,
        $res_py:ty, $res_parse:path, $res_wrap:path
    ) => {
        #[pyclass(name = $name, module = "pleroma", from_py_object)]
        #[derive(Clone)]
        pub(crate) struct $py {
            inner: Qp<$p, $k>,
        }

        #[pymethods]
        impl $py {
            #[staticmethod]
            fn from_i128(value: i128) -> Self {
                $wrap(Qp::<$p, $k>::from_i128(value))
            }
            #[staticmethod]
            fn zero() -> Self {
                $wrap(Qp::<$p, $k>::zero())
            }
            #[staticmethod]
            fn one() -> Self {
                $wrap(Qp::<$p, $k>::one())
            }
            #[staticmethod]
            fn from_p_power(v: i128) -> Self {
                $wrap(Qp::<$p, $k>::from_p_power(v))
            }
            #[staticmethod]
            fn uniformizer() -> Self {
                $wrap(<Qp<$p, $k> as Valued>::uniformizer())
            }
            #[staticmethod]
            fn teichmuller(residue: &Bound<'_, PyAny>) -> PyResult<Self> {
                Ok($wrap(<Qp<$p, $k> as ResidueField>::teichmuller(
                    $res_parse(residue)?,
                )))
            }
            #[staticmethod]
            fn from_rational(num: i128, den: i128) -> PyResult<Self> {
                let q = rational_from_pair(num, den)?;
                Ok($wrap(Qp::<$p, $k>::from_rational(&q)))
            }
            #[staticmethod]
            fn prime() -> u128 {
                $p
            }
            #[staticmethod]
            fn precision() -> u128 {
                $k
            }
            #[staticmethod]
            fn modulus() -> u128 {
                Qp::<$p, $k>::modulus()
            }
            #[staticmethod]
            fn assert_supported_field() {
                Qp::<$p, $k>::assert_supported_field()
            }
            #[staticmethod]
            fn characteristic() -> u128 {
                Qp::<$p, $k>::characteristic()
            }
            #[getter]
            fn unit(&self) -> u128 {
                self.inner.unit()
            }
            fn valuation(&self) -> Option<i128> {
                self.inner.valuation()
            }
            fn residue(&self) -> Option<$res_py> {
                <Qp<$p, $k> as ResidueField>::residue(&self.inner).map($res_wrap)
            }
            fn residue_unit(&self) -> Option<$res_py> {
                <Qp<$p, $k> as ResidueField>::residue_unit(&self.inner).map($res_wrap)
            }
            fn angular_component(&self) -> Option<$res_py> {
                self.residue_unit()
            }
            fn is_integral(&self) -> bool {
                <Qp<$p, $k> as HasRingOfIntegers>::is_integral(&self.inner)
            }
            fn to_integer(&self) -> Option<$int_py> {
                <Qp<$p, $k> as HasRingOfIntegers>::to_integer(&self.inner).map($int_wrap)
            }
            fn is_zero(&self) -> bool {
                self.inner.is_zero()
            }
            fn is_square(&self) -> PyResult<bool> {
                self.inner.is_square().ok_or_else(|| {
                    PyValueError::new_err("squarehood is undecidable at this p-adic precision")
                })
            }
            fn sqrt(&self) -> PyResult<Option<Self>> {
                self.inner
                    .sqrt()
                    .map(|root| root.map($wrap))
                    .ok_or_else(|| {
                        PyValueError::new_err(
                            "square-root construction is not implemented for this p-adic case",
                        )
                    })
            }
            fn inv(&self) -> PyResult<Self> {
                self.inner
                    .inv()
                    .map($wrap)
                    .ok_or_else(|| PyValueError::new_err("0 has no inverse in Q_p"))
            }
            fn __add__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                Ok($wrap(self.inner.add(&$parse(other)?)))
            }
            fn __radd__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                self.__add__(other)
            }
            fn __sub__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                Ok($wrap(self.inner.sub(&$parse(other)?)))
            }
            fn __rsub__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                Ok($wrap($parse(other)?.sub(&self.inner)))
            }
            fn __neg__(&self) -> Self {
                $wrap(self.inner.neg())
            }
            fn __mul__(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
                match $parse(other) {
                    Ok(o) => $wrap(self.inner.mul(&o)).into_py_any(py),
                    Err(_) => Ok(py.NotImplemented()),
                }
            }
            fn __rmul__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                Ok($wrap(self.inner.mul(&$parse(other)?)))
            }
            fn __truediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                let o = $parse(other)?;
                let oi = o
                    .inv()
                    .ok_or_else(|| PyValueError::new_err("division by 0 in Q_p"))?;
                Ok($wrap(self.inner.mul(&oi)))
            }
            fn __rtruediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                let si = self
                    .inner
                    .inv()
                    .ok_or_else(|| PyValueError::new_err("division by 0 in Q_p"))?;
                Ok($wrap($parse(other)?.mul(&si)))
            }
            fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
                matches!(other.cast::<$py>(), Ok(x) if x.borrow().inner == self.inner)
            }
            fn __repr__(&self) -> String {
                format!("{:?}", self.inner)
            }
        }

        pub(crate) fn $parse(obj: &Bound<'_, PyAny>) -> PyResult<Qp<$p, $k>> {
            if let Ok(x) = obj.cast::<$py>() {
                return Ok(x.borrow().inner);
            }
            if let Ok(v) = obj.extract::<i128>() {
                return Ok(Qp::<$p, $k>::from_i128(v));
            }
            Err(PyTypeError::new_err(concat!("expected ", $name, " or int")))
        }

        pub(crate) fn $wrap(inner: Qp<$p, $k>) -> $py {
            $py { inner }
        }
    };
}

zp_pyclass!(PyZp2_4, "Zp2_4", parse_zp2_4, wrap_zp2_4, 2, 4);
zp_pyclass!(PyZp3_4, "Zp3_4", parse_zp3_4, wrap_zp3_4, 3, 4);
zp_pyclass!(PyZp5_4, "Zp5_4", parse_zp5_4, wrap_zp5_4, 5, 4);
zp_pyclass!(PyZp7_4, "Zp7_4", parse_zp7_4, wrap_zp7_4, 7, 4);
zp_pyclass!(PyZp11_4, "Zp11_4", parse_zp11_4, wrap_zp11_4, 11, 4);
zp_pyclass!(PyZp13_4, "Zp13_4", parse_zp13_4, wrap_zp13_4, 13, 4);

qp_pyclass!(
    PyQp2_4,
    "Qp2_4",
    parse_qp2_4,
    wrap_qp2_4,
    2,
    4,
    PyZp2_4,
    wrap_zp2_4,
    PyFp2,
    parse_fp2,
    wrap_fp2
);
qp_pyclass!(
    PyQp3_4,
    "Qp3_4",
    parse_qp3_4,
    wrap_qp3_4,
    3,
    4,
    PyZp3_4,
    wrap_zp3_4,
    PyFp3,
    parse_fp3,
    wrap_fp3
);
qp_pyclass!(
    PyQp5_4,
    "Qp5_4",
    parse_qp5_4,
    wrap_qp5_4,
    5,
    4,
    PyZp5_4,
    wrap_zp5_4,
    PyFp5,
    parse_fp5,
    wrap_fp5
);
qp_pyclass!(
    PyQp7_4,
    "Qp7_4",
    parse_qp7_4,
    wrap_qp7_4,
    7,
    4,
    PyZp7_4,
    wrap_zp7_4,
    PyFp7,
    parse_fp7,
    wrap_fp7
);
qp_pyclass!(
    PyQp11_4,
    "Qp11_4",
    parse_qp11_4,
    wrap_qp11_4,
    11,
    4,
    PyZp11_4,
    wrap_zp11_4,
    PyFp11,
    parse_fp11,
    wrap_fp11
);
qp_pyclass!(
    PyQp13_4,
    "Qp13_4",
    parse_qp13_4,
    wrap_qp13_4,
    13,
    4,
    PyZp13_4,
    wrap_zp13_4,
    PyFp13,
    parse_fp13,
    wrap_fp13
);

fn parse_witt_coords<const P: u128, const N: usize, const F: usize>(
    coords: &[u128],
) -> PyResult<[u128; F]> {
    if coords.len() != F {
        return Err(PyValueError::new_err(format!(
            "expected exactly {F} unramified-ring coordinates"
        )));
    }
    let modulus = WittVec::<P, N, F>::modulus();
    let mut out = [0u128; F];
    for (i, &coord) in coords.iter().enumerate() {
        out[i] = coord % modulus;
    }
    Ok(out)
}

macro_rules! witt_vec_pyclass {
    (
        $py:ident, $name:literal, $parse:ident, $wrap:ident,
        $p:literal, $n:literal, $f:literal,
        $res_py:ty, $res_parse:path, $res_wrap:path
    ) => {
        #[pyclass(name = $name, module = "pleroma", from_py_object)]
        #[derive(Clone)]
        pub(crate) struct $py {
            inner: WittVec<$p, $n, $f>,
        }

        #[pymethods]
        impl $py {
            #[new]
            fn new(coords: Vec<u128>) -> PyResult<Self> {
                Ok($wrap(WittVec::<$p, $n, $f>(parse_witt_coords::<
                    $p,
                    $n,
                    $f,
                >(&coords)?)))
            }
            #[staticmethod]
            fn zero() -> Self {
                $wrap(WittVec::<$p, $n, $f>::zero())
            }
            #[staticmethod]
            fn one() -> Self {
                $wrap(WittVec::<$p, $n, $f>::one())
            }
            #[staticmethod]
            fn from_int(value: i128) -> Self {
                $wrap(WittVec::<$p, $n, $f>::from_int(value))
            }
            #[staticmethod]
            fn teichmuller(x: &Bound<'_, PyAny>) -> PyResult<Self> {
                Ok($wrap(WittVec::<$p, $n, $f>::teichmuller($res_parse(
                    x,
                )?)))
            }
            #[staticmethod]
            fn from_witt_components(xs: Vec<Bound<'_, PyAny>>) -> PyResult<Self> {
                if xs.len() != $n {
                    return Err(PyValueError::new_err(concat!(
                        "expected exactly ",
                        stringify!($n),
                        " Teichmuller digits"
                    )));
                }
                let mut parsed = Vec::with_capacity(xs.len());
                for x in &xs {
                    parsed.push($res_parse(x)?);
                }
                Ok($wrap(WittVec::<$p, $n, $f>::from_witt_components(
                    &parsed,
                )))
            }
            #[staticmethod]
            fn prime() -> u128 {
                $p
            }
            #[staticmethod]
            fn precision() -> usize {
                $n
            }
            #[staticmethod]
            fn residue_degree() -> usize {
                $f
            }
            #[staticmethod]
            fn modulus() -> u128 {
                WittVec::<$p, $n, $f>::modulus()
            }
            #[staticmethod]
            fn residue_order() -> u128 {
                WittVec::<$p, $n, $f>::residue_order()
            }
            #[staticmethod]
            fn characteristic() -> u128 {
                WittVec::<$p, $n, $f>::characteristic()
            }
            #[getter]
            fn coords(&self) -> Vec<u128> {
                self.inner.0.to_vec()
            }
            fn residue(&self) -> $res_py {
                $res_wrap(self.inner.residue())
            }
            fn witt_components(&self) -> Vec<$res_py> {
                self.inner
                    .witt_components()
                    .into_iter()
                    .map($res_wrap)
                    .collect()
            }
            fn p_valuation(&self) -> usize {
                self.inner.p_valuation()
            }
            fn try_divide_by_p(&self) -> Option<Self> {
                self.inner.try_divide_by_p().map($wrap)
            }
            fn is_zero(&self) -> bool {
                self.inner.is_zero()
            }
            fn is_square(&self) -> PyResult<bool> {
                self.inner.is_square().ok_or_else(|| {
                    PyValueError::new_err("squarehood is undecidable at this Witt precision")
                })
            }
            fn sqrt(&self) -> PyResult<Option<Self>> {
                self.inner
                    .sqrt()
                    .map(|root| root.map($wrap))
                    .ok_or_else(|| {
                        PyValueError::new_err(
                            "square-root construction is not implemented for this Witt case",
                        )
                    })
            }
            fn inv(&self) -> PyResult<Self> {
                self.inner
                    .inv()
                    .map($wrap)
                    .ok_or_else(|| PyValueError::new_err("non-unit has no inverse in W_N(F_q)"))
            }
            fn __add__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                Ok($wrap(self.inner.add(&$parse(other)?)))
            }
            fn __radd__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                self.__add__(other)
            }
            fn __sub__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                Ok($wrap(self.inner.sub(&$parse(other)?)))
            }
            fn __rsub__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                Ok($wrap($parse(other)?.sub(&self.inner)))
            }
            fn __neg__(&self) -> Self {
                $wrap(self.inner.neg())
            }
            fn __mul__(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
                match $parse(other) {
                    Ok(o) => $wrap(self.inner.mul(&o)).into_py_any(py),
                    Err(_) => Ok(py.NotImplemented()),
                }
            }
            fn __rmul__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                Ok($wrap(self.inner.mul(&$parse(other)?)))
            }
            fn __truediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                let o = $parse(other)?;
                let oi = o
                    .inv()
                    .ok_or_else(|| PyValueError::new_err("divisor is not a Witt unit"))?;
                Ok($wrap(self.inner.mul(&oi)))
            }
            fn __rtruediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                let si = self
                    .inner
                    .inv()
                    .ok_or_else(|| PyValueError::new_err("divisor is not a Witt unit"))?;
                Ok($wrap($parse(other)?.mul(&si)))
            }
            fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
                matches!(other.cast::<$py>(), Ok(x) if x.borrow().inner == self.inner)
            }
            fn __hash__(&self) -> usize {
                self.inner.0.iter().fold(0usize, |acc, &x| {
                    acc.wrapping_mul(257).wrapping_add(x as usize)
                })
            }
            fn __repr__(&self) -> String {
                format!("{:?}", self.inner)
            }
        }

        pub(crate) fn $parse(obj: &Bound<'_, PyAny>) -> PyResult<WittVec<$p, $n, $f>> {
            if let Ok(x) = obj.cast::<$py>() {
                return Ok(x.borrow().inner);
            }
            if let Ok(coords) = obj.extract::<Vec<u128>>() {
                return Ok(WittVec::<$p, $n, $f>(parse_witt_coords::<
                    $p,
                    $n,
                    $f,
                >(&coords)?));
            }
            if let Ok(v) = obj.extract::<i128>() {
                return Ok(WittVec::<$p, $n, $f>::from_int(v));
            }
            Err(PyTypeError::new_err(concat!(
                "expected ",
                $name,
                ", coordinate list, or int"
            )))
        }

        pub(crate) fn $wrap(inner: WittVec<$p, $n, $f>) -> $py {
            $py { inner }
        }
    };
}

macro_rules! qq_pyclass {
    (
        $py:ident, $name:literal, $parse:ident, $wrap:ident,
        $p:literal, $n:literal, $f:literal,
        $base_py:ty, $base_parse:path, $base_wrap:path,
        $witt_py:ty, $witt_parse:path, $witt_wrap:path,
        $res_py:ty, $res_parse:path, $res_wrap:path
    ) => {
        #[pyclass(name = $name, module = "pleroma", from_py_object)]
        #[derive(Clone)]
        pub(crate) struct $py {
            inner: Qq<$p, $n, $f>,
        }

        #[pymethods]
        impl $py {
            #[new]
            fn new(value: i128) -> Self {
                $wrap(Qq::<$p, $n, $f>::from_int(value))
            }
            #[staticmethod]
            fn zero() -> Self {
                $wrap(Qq::<$p, $n, $f>::zero())
            }
            #[staticmethod]
            fn one() -> Self {
                $wrap(Qq::<$p, $n, $f>::one())
            }
            #[staticmethod]
            fn from_int(value: i128) -> Self {
                $wrap(Qq::<$p, $n, $f>::from_int(value))
            }
            #[staticmethod]
            fn from_witt(w: &Bound<'_, PyAny>) -> PyResult<Self> {
                Ok($wrap(Qq::<$p, $n, $f>::from_witt($witt_parse(w)?)))
            }
            #[staticmethod]
            fn from_p_power(v: i128) -> Self {
                $wrap(Qq::<$p, $n, $f>::from_p_power(v))
            }
            #[staticmethod]
            fn uniformizer() -> Self {
                $wrap(<Qq<$p, $n, $f> as Valued>::uniformizer())
            }
            #[staticmethod]
            fn teichmuller(x: &Bound<'_, PyAny>) -> PyResult<Self> {
                Ok($wrap(Qq::<$p, $n, $f>::teichmuller($res_parse(x)?)))
            }
            #[staticmethod]
            fn prime() -> u128 {
                $p
            }
            #[staticmethod]
            fn precision() -> usize {
                $n
            }
            #[staticmethod]
            fn residue_degree() -> usize {
                $f
            }
            #[staticmethod]
            fn characteristic() -> u128 {
                Qq::<$p, $n, $f>::characteristic()
            }
            #[staticmethod]
            fn extension_degree() -> usize {
                <Qq<$p, $n, $f> as FieldExtension>::extension_degree()
            }
            #[staticmethod]
            fn embed(base: &Bound<'_, PyAny>) -> PyResult<Self> {
                let base = qp_to_qq_base::<$p, $n, $n>($base_parse(base)?);
                Ok($wrap(<Qq<$p, $n, $f> as FieldExtension>::embed(&base)))
            }
            fn valuation(&self) -> Option<i128> {
                self.inner.valuation()
            }
            fn unit(&self) -> $witt_py {
                $witt_wrap(self.inner.unit())
            }
            fn unit_residue(&self) -> Option<$res_py> {
                self.inner.unit_residue().map($res_wrap)
            }
            fn residue(&self) -> Option<$res_py> {
                <Qq<$p, $n, $f> as ResidueField>::residue(&self.inner).map($res_wrap)
            }
            fn residue_unit(&self) -> Option<$res_py> {
                <Qq<$p, $n, $f> as ResidueField>::residue_unit(&self.inner).map($res_wrap)
            }
            fn angular_component(&self) -> Option<$res_py> {
                self.residue_unit()
            }
            fn is_integral(&self) -> bool {
                <Qq<$p, $n, $f> as HasRingOfIntegers>::is_integral(&self.inner)
            }
            fn to_integer(&self) -> Option<$witt_py> {
                <Qq<$p, $n, $f> as HasRingOfIntegers>::to_integer(&self.inner).map($witt_wrap)
            }
            fn trace(&self) -> $base_py {
                $base_wrap(qq_base_to_qp::<$p, $n, $n>(
                    <Qq<$p, $n, $f> as FieldExtension>::trace(&self.inner),
                ))
            }
            fn norm(&self) -> $base_py {
                $base_wrap(qq_base_to_qp::<$p, $n, $n>(
                    <Qq<$p, $n, $f> as FieldExtension>::norm(&self.inner),
                ))
            }
            #[staticmethod]
            fn basis() -> Vec<Self> {
                <Qq<$p, $n, $f> as CyclicGaloisExtension>::basis()
                    .into_iter()
                    .map($wrap)
                    .collect()
            }
            fn sigma(&self) -> Self {
                $wrap(<Qq<$p, $n, $f> as CyclicGaloisExtension>::sigma(
                    &self.inner,
                ))
            }
            fn sigma_power(&self, k: usize) -> Self {
                $wrap(<Qq<$p, $n, $f> as CyclicGaloisExtension>::sigma_power(
                    &self.inner,
                    k,
                ))
            }
            fn is_zero(&self) -> bool {
                self.inner.is_zero()
            }
            fn is_square(&self) -> PyResult<bool> {
                self.inner.is_square().ok_or_else(|| {
                    PyValueError::new_err("squarehood is undecidable at this Qq precision")
                })
            }
            fn sqrt(&self) -> PyResult<Option<Self>> {
                self.inner
                    .sqrt()
                    .map(|root| root.map($wrap))
                    .ok_or_else(|| {
                        PyValueError::new_err(
                            "square-root construction is not implemented for this Qq case",
                        )
                    })
            }
            fn inv(&self) -> PyResult<Self> {
                self.inner
                    .inv()
                    .map($wrap)
                    .ok_or_else(|| PyValueError::new_err("0 has no inverse in Q_q"))
            }
            fn __add__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                Ok($wrap(self.inner.add(&$parse(other)?)))
            }
            fn __radd__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                self.__add__(other)
            }
            fn __sub__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                Ok($wrap(self.inner.sub(&$parse(other)?)))
            }
            fn __rsub__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                Ok($wrap($parse(other)?.sub(&self.inner)))
            }
            fn __neg__(&self) -> Self {
                $wrap(self.inner.neg())
            }
            fn __mul__(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
                match $parse(other) {
                    Ok(o) => $wrap(self.inner.mul(&o)).into_py_any(py),
                    Err(_) => Ok(py.NotImplemented()),
                }
            }
            fn __rmul__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                Ok($wrap(self.inner.mul(&$parse(other)?)))
            }
            fn __truediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                let o = $parse(other)?;
                let oi = o
                    .inv()
                    .ok_or_else(|| PyValueError::new_err("division by 0 in Q_q"))?;
                Ok($wrap(self.inner.mul(&oi)))
            }
            fn __rtruediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                let si = self
                    .inner
                    .inv()
                    .ok_or_else(|| PyValueError::new_err("division by 0 in Q_q"))?;
                Ok($wrap($parse(other)?.mul(&si)))
            }
            fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
                matches!($parse(other), Ok(x) if x == self.inner)
            }
            fn __repr__(&self) -> String {
                format!("{:?}", self.inner)
            }
        }

        pub(crate) fn $parse(obj: &Bound<'_, PyAny>) -> PyResult<Qq<$p, $n, $f>> {
            if let Ok(x) = obj.cast::<$py>() {
                return Ok(x.borrow().inner);
            }
            if let Ok(w) = $witt_parse(obj) {
                return Ok(Qq::<$p, $n, $f>::from_witt(w));
            }
            if let Ok(v) = obj.extract::<i128>() {
                return Ok(Qq::<$p, $n, $f>::from_int(v));
            }
            Err(PyTypeError::new_err(concat!(
                "expected ",
                $name,
                ", matching Witt vector, coordinate list, or int"
            )))
        }

        pub(crate) fn $wrap(inner: Qq<$p, $n, $f>) -> $py {
            $py { inner }
        }
    };
}

witt_vec_pyclass!(
    PyWittVec2_4_2,
    "WittVec2_4_2",
    parse_witt_vec2_4_2,
    wrap_witt_vec2_4_2,
    2,
    4,
    2,
    PyF4,
    parse_f4,
    wrap_f4
);
witt_vec_pyclass!(
    PyWittVec2_4_3,
    "WittVec2_4_3",
    parse_witt_vec2_4_3,
    wrap_witt_vec2_4_3,
    2,
    4,
    3,
    PyF8,
    parse_f8,
    wrap_f8
);
witt_vec_pyclass!(
    PyWittVec2_4_4,
    "WittVec2_4_4",
    parse_witt_vec2_4_4,
    wrap_witt_vec2_4_4,
    2,
    4,
    4,
    PyF16,
    parse_f16,
    wrap_f16
);
witt_vec_pyclass!(
    PyWittVec3_4_2,
    "WittVec3_4_2",
    parse_witt_vec3_4_2,
    wrap_witt_vec3_4_2,
    3,
    4,
    2,
    PyF9,
    parse_f9,
    wrap_f9
);
witt_vec_pyclass!(
    PyWittVec5_4_2,
    "WittVec5_4_2",
    parse_witt_vec5_4_2,
    wrap_witt_vec5_4_2,
    5,
    4,
    2,
    PyF25,
    parse_f25,
    wrap_f25
);
witt_vec_pyclass!(
    PyWittVec3_4_3,
    "WittVec3_4_3",
    parse_witt_vec3_4_3,
    wrap_witt_vec3_4_3,
    3,
    4,
    3,
    PyF27,
    parse_f27,
    wrap_f27
);

qq_pyclass!(
    PyQq2_4_2,
    "Qq2_4_2",
    parse_qq2_4_2,
    wrap_qq2_4_2,
    2,
    4,
    2,
    PyQp2_4,
    parse_qp2_4,
    wrap_qp2_4,
    PyWittVec2_4_2,
    parse_witt_vec2_4_2,
    wrap_witt_vec2_4_2,
    PyF4,
    parse_f4,
    wrap_f4
);
qq_pyclass!(
    PyQq2_4_3,
    "Qq2_4_3",
    parse_qq2_4_3,
    wrap_qq2_4_3,
    2,
    4,
    3,
    PyQp2_4,
    parse_qp2_4,
    wrap_qp2_4,
    PyWittVec2_4_3,
    parse_witt_vec2_4_3,
    wrap_witt_vec2_4_3,
    PyF8,
    parse_f8,
    wrap_f8
);
qq_pyclass!(
    PyQq2_4_4,
    "Qq2_4_4",
    parse_qq2_4_4,
    wrap_qq2_4_4,
    2,
    4,
    4,
    PyQp2_4,
    parse_qp2_4,
    wrap_qp2_4,
    PyWittVec2_4_4,
    parse_witt_vec2_4_4,
    wrap_witt_vec2_4_4,
    PyF16,
    parse_f16,
    wrap_f16
);
qq_pyclass!(
    PyQq3_4_2,
    "Qq3_4_2",
    parse_qq3_4_2,
    wrap_qq3_4_2,
    3,
    4,
    2,
    PyQp3_4,
    parse_qp3_4,
    wrap_qp3_4,
    PyWittVec3_4_2,
    parse_witt_vec3_4_2,
    wrap_witt_vec3_4_2,
    PyF9,
    parse_f9,
    wrap_f9
);
qq_pyclass!(
    PyQq5_4_2,
    "Qq5_4_2",
    parse_qq5_4_2,
    wrap_qq5_4_2,
    5,
    4,
    2,
    PyQp5_4,
    parse_qp5_4,
    wrap_qp5_4,
    PyWittVec5_4_2,
    parse_witt_vec5_4_2,
    wrap_witt_vec5_4_2,
    PyF25,
    parse_f25,
    wrap_f25
);
qq_pyclass!(
    PyQq3_4_3,
    "Qq3_4_3",
    parse_qq3_4_3,
    wrap_qq3_4_3,
    3,
    4,
    3,
    PyQp3_4,
    parse_qp3_4,
    wrap_qp3_4,
    PyWittVec3_4_3,
    parse_witt_vec3_4_3,
    wrap_witt_vec3_4_3,
    PyF27,
    parse_f27,
    wrap_f27
);

macro_rules! laurent_pyclass {
    (
        $py:ident, $name:literal, $parse:ident, $wrap:ident,
        $base:ty, $k:literal, $base_py:ty, $base_parse:path, $base_wrap:path
    ) => {
        #[pyclass(name = $name, module = "pleroma", from_py_object)]
        #[derive(Clone)]
        pub(crate) struct $py {
            inner: Laurent<$base, $k>,
        }

        pub(crate) fn $parse(obj: &Bound<'_, PyAny>) -> PyResult<Laurent<$base, $k>> {
            if let Ok(x) = obj.cast::<$py>() {
                return Ok(x.borrow().inner.clone());
            }
            if let Ok(s) = $base_parse(obj) {
                return Ok(Laurent::<$base, $k>::from_scalar(s));
            }
            if let Ok(items) = obj.extract::<Vec<Bound<'_, PyAny>>>() {
                let mut coeffs = Vec::with_capacity(items.len());
                for item in &items {
                    coeffs.push($base_parse(item)?);
                }
                return Ok(Laurent::<$base, $k>::from_coeffs(coeffs, 0));
            }
            Err(PyTypeError::new_err(concat!(
                "expected ",
                $name,
                ", base scalar, or coefficient list"
            )))
        }

        pub(crate) fn $wrap(inner: Laurent<$base, $k>) -> $py {
            $py { inner }
        }

        #[pymethods]
        impl $py {
            #[new]
            #[pyo3(signature = (coeffs, valuation=0))]
            fn new(coeffs: Vec<Bound<'_, PyAny>>, valuation: i128) -> PyResult<Self> {
                Self::from_coeffs(coeffs, valuation)
            }
            #[staticmethod]
            fn zero() -> Self {
                $wrap(Laurent::<$base, $k>::zero())
            }
            #[staticmethod]
            fn one() -> Self {
                $wrap(Laurent::<$base, $k>::one())
            }
            #[staticmethod]
            fn from_scalar(s: &Bound<'_, PyAny>) -> PyResult<Self> {
                Ok($wrap(Laurent::<$base, $k>::from_scalar($base_parse(s)?)))
            }
            #[staticmethod]
            fn teichmuller(residue: &Bound<'_, PyAny>) -> PyResult<Self> {
                Ok($wrap(<Laurent<$base, $k> as ResidueField>::teichmuller(
                    $base_parse(residue)?,
                )))
            }
            #[staticmethod]
            fn t() -> Self {
                $wrap(Laurent::<$base, $k>::t())
            }
            #[staticmethod]
            fn uniformizer() -> Self {
                $wrap(<Laurent<$base, $k> as Valued>::uniformizer())
            }
            #[staticmethod]
            fn from_t_power(v: i128) -> Self {
                $wrap(Laurent::<$base, $k>::from_t_power(v))
            }
            #[staticmethod]
            #[pyo3(signature = (coeffs, valuation=0))]
            fn from_coeffs(coeffs: Vec<Bound<'_, PyAny>>, valuation: i128) -> PyResult<Self> {
                let mut parsed = Vec::with_capacity(coeffs.len());
                for coeff in &coeffs {
                    parsed.push($base_parse(coeff)?);
                }
                Ok($wrap(Laurent::<$base, $k>::from_coeffs(parsed, valuation)))
            }
            #[staticmethod]
            fn precision() -> usize {
                Laurent::<$base, $k>::precision()
            }
            #[staticmethod]
            fn characteristic() -> u128 {
                Laurent::<$base, $k>::characteristic()
            }
            fn valuation(&self) -> Option<i128> {
                self.inner.valuation()
            }
            fn leading_coeff(&self) -> Option<$base_py> {
                self.inner.leading_coeff().map($base_wrap)
            }
            fn residue(&self) -> Option<$base_py> {
                <Laurent<$base, $k> as ResidueField>::residue(&self.inner).map($base_wrap)
            }
            fn residue_unit(&self) -> Option<$base_py> {
                <Laurent<$base, $k> as ResidueField>::residue_unit(&self.inner).map($base_wrap)
            }
            fn angular_component(&self) -> Option<$base_py> {
                self.residue_unit()
            }
            fn is_integral(&self) -> bool {
                self.inner.is_integral()
            }
            fn unit_coeffs(&self) -> Vec<$base_py> {
                self.inner
                    .unit_coeffs()
                    .iter()
                    .cloned()
                    .map($base_wrap)
                    .collect()
            }
            fn coeff(&self, exp: i128) -> $base_py {
                $base_wrap(self.inner.coeff(exp))
            }
            fn is_zero(&self) -> bool {
                self.inner.is_zero()
            }
            fn is_square(&self) -> bool {
                ExactRoots::is_square(&self.inner)
            }
            fn sqrt(&self) -> Option<Self> {
                ExactRoots::sqrt(&self.inner).map($wrap)
            }
            fn inv(&self) -> PyResult<Self> {
                self.inner
                    .inv()
                    .map($wrap)
                    .ok_or_else(|| PyValueError::new_err("0 has no inverse in this Laurent field"))
            }
            fn __add__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                Ok($wrap(self.inner.add(&$parse(other)?)))
            }
            fn __radd__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                self.__add__(other)
            }
            fn __sub__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                Ok($wrap(self.inner.sub(&$parse(other)?)))
            }
            fn __rsub__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                Ok($wrap($parse(other)?.sub(&self.inner)))
            }
            fn __neg__(&self) -> Self {
                $wrap(self.inner.neg())
            }
            fn __mul__(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
                match $parse(other) {
                    Ok(o) => $wrap(self.inner.mul(&o)).into_py_any(py),
                    Err(_) => Ok(py.NotImplemented()),
                }
            }
            fn __rmul__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                Ok($wrap(self.inner.mul(&$parse(other)?)))
            }
            fn __truediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                let o = $parse(other)?;
                let oi = o
                    .inv()
                    .ok_or_else(|| PyValueError::new_err("division by 0 in this Laurent field"))?;
                Ok($wrap(self.inner.mul(&oi)))
            }
            fn __rtruediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                let si = self
                    .inner
                    .inv()
                    .ok_or_else(|| PyValueError::new_err("division by 0 in this Laurent field"))?;
                Ok($wrap($parse(other)?.mul(&si)))
            }
            fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
                matches!($parse(other), Ok(x) if x == self.inner)
            }
            fn __repr__(&self) -> String {
                format!("{:?}", self.inner)
            }
        }
    };
}

laurent_pyclass!(
    PyLaurentRational6,
    "LaurentRational_6",
    parse_laurent_rational_6,
    wrap_laurent_rational_6,
    Rational,
    6,
    PyRational,
    parse_rational,
    wrap_rational
);
laurent_pyclass!(
    PyLaurentFp3_6,
    "LaurentFp3_6",
    parse_laurent_fp3_6,
    wrap_laurent_fp3_6,
    Fp<3>,
    6,
    PyFp3,
    parse_fp3,
    wrap_fp3
);
laurent_pyclass!(
    PyLaurentFp5_6,
    "LaurentFp5_6",
    parse_laurent_fp5_6,
    wrap_laurent_fp5_6,
    Fp<5>,
    6,
    PyFp5,
    parse_fp5,
    wrap_fp5
);
laurent_pyclass!(
    PyLaurentFp7_6,
    "LaurentFp7_6",
    parse_laurent_fp7_6,
    wrap_laurent_fp7_6,
    Fp<7>,
    6,
    PyFp7,
    parse_fp7,
    wrap_fp7
);
laurent_pyclass!(
    PyLaurentFp11_6,
    "LaurentFp11_6",
    parse_laurent_fp11_6,
    wrap_laurent_fp11_6,
    Fp<11>,
    6,
    PyFp11,
    parse_fp11,
    wrap_fp11
);
laurent_pyclass!(
    PyLaurentFp13_6,
    "LaurentFp13_6",
    parse_laurent_fp13_6,
    wrap_laurent_fp13_6,
    Fp<13>,
    6,
    PyFp13,
    parse_fp13,
    wrap_fp13
);
laurent_pyclass!(
    PyLaurentF9_6,
    "LaurentF9_6",
    parse_laurent_f9_6,
    wrap_laurent_f9_6,
    Fpn<3, 2>,
    6,
    PyF9,
    parse_f9,
    wrap_f9
);
laurent_pyclass!(
    PyLaurentF25_6,
    "LaurentF25_6",
    parse_laurent_f25_6,
    wrap_laurent_f25_6,
    Fpn<5, 2>,
    6,
    PyF25,
    parse_f25,
    wrap_f25
);
laurent_pyclass!(
    PyLaurentF27_6,
    "LaurentF27_6",
    parse_laurent_f27_6,
    wrap_laurent_f27_6,
    Fpn<3, 3>,
    6,
    PyF27,
    parse_f27,
    wrap_f27
);

macro_rules! ramified_pyclass {
    (
        $py:ident, $name:literal, $parse:ident, $wrap:ident,
        $base:ty, $e:literal, $base_py:ty, $base_parse:path, $base_wrap:path,
        $res_py:ty, $res_parse:path, $res_wrap:path
    ) => {
        #[pyclass(name = $name, module = "pleroma", from_py_object)]
        #[derive(Clone)]
        pub(crate) struct $py {
            inner: Ramified<$base, $e>,
        }

        pub(crate) fn $parse(obj: &Bound<'_, PyAny>) -> PyResult<Ramified<$base, $e>> {
            if let Ok(x) = obj.cast::<$py>() {
                return Ok(x.borrow().inner.clone());
            }
            if let Ok(s) = $base_parse(obj) {
                return Ok(Ramified::<$base, $e>::from_base(s));
            }
            if let Ok(items) = obj.extract::<Vec<Bound<'_, PyAny>>>() {
                let mut coeffs = Vec::with_capacity(items.len());
                for item in &items {
                    coeffs.push($base_parse(item)?);
                }
                return Ok(Ramified::<$base, $e>::new(coeffs));
            }
            Err(PyTypeError::new_err(concat!(
                "expected ",
                $name,
                ", base scalar, or component list"
            )))
        }

        pub(crate) fn $wrap(inner: Ramified<$base, $e>) -> $py {
            $py { inner }
        }

        #[pymethods]
        impl $py {
            #[new]
            fn new(components: Vec<Bound<'_, PyAny>>) -> PyResult<Self> {
                Self::from_components(components)
            }
            #[staticmethod]
            fn zero() -> Self {
                $wrap(Ramified::<$base, $e>::zero())
            }
            #[staticmethod]
            fn one() -> Self {
                $wrap(Ramified::<$base, $e>::one())
            }
            #[staticmethod]
            fn from_base(s: &Bound<'_, PyAny>) -> PyResult<Self> {
                Ok($wrap(Ramified::<$base, $e>::from_base($base_parse(s)?)))
            }
            #[staticmethod]
            fn teichmuller(residue: &Bound<'_, PyAny>) -> PyResult<Self> {
                Ok($wrap(<Ramified<$base, $e> as ResidueField>::teichmuller(
                    $res_parse(residue)?,
                )))
            }
            #[staticmethod]
            fn pi() -> Self {
                $wrap(Ramified::<$base, $e>::pi())
            }
            #[staticmethod]
            fn uniformizer() -> Self {
                $wrap(<Ramified<$base, $e> as Valued>::uniformizer())
            }
            #[staticmethod]
            fn degree() -> usize {
                $e
            }
            #[staticmethod]
            fn characteristic() -> u128 {
                Ramified::<$base, $e>::characteristic()
            }
            #[staticmethod]
            fn from_components(components: Vec<Bound<'_, PyAny>>) -> PyResult<Self> {
                let mut parsed = Vec::with_capacity(components.len());
                for component in &components {
                    parsed.push($base_parse(component)?);
                }
                Ok($wrap(Ramified::<$base, $e>::new(parsed)))
            }
            fn components(&self) -> Vec<$base_py> {
                self.inner
                    .components()
                    .iter()
                    .cloned()
                    .map($base_wrap)
                    .collect()
            }
            fn valuation(&self) -> Option<i128> {
                self.inner.valuation()
            }
            fn residue(&self) -> Option<$res_py> {
                <Ramified<$base, $e> as ResidueField>::residue(&self.inner).map($res_wrap)
            }
            fn residue_unit(&self) -> Option<$res_py> {
                <Ramified<$base, $e> as ResidueField>::residue_unit(&self.inner).map($res_wrap)
            }
            fn angular_component(&self) -> Option<$res_py> {
                self.residue_unit()
            }
            fn is_integral(&self) -> bool {
                self.inner.is_integral()
            }
            fn is_zero(&self) -> bool {
                self.inner.is_zero()
            }
            fn inv(&self) -> PyResult<Self> {
                self.inner
                    .inv()
                    .map($wrap)
                    .ok_or_else(|| PyValueError::new_err("0 has no inverse in this ramified field"))
            }
            fn __add__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                Ok($wrap(self.inner.add(&$parse(other)?)))
            }
            fn __radd__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                self.__add__(other)
            }
            fn __sub__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                Ok($wrap(self.inner.sub(&$parse(other)?)))
            }
            fn __rsub__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                Ok($wrap($parse(other)?.sub(&self.inner)))
            }
            fn __neg__(&self) -> Self {
                $wrap(self.inner.neg())
            }
            fn __mul__(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
                match $parse(other) {
                    Ok(o) => $wrap(self.inner.mul(&o)).into_py_any(py),
                    Err(_) => Ok(py.NotImplemented()),
                }
            }
            fn __rmul__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                Ok($wrap(self.inner.mul(&$parse(other)?)))
            }
            fn __truediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                let o = $parse(other)?;
                let oi = o.inv().ok_or_else(|| {
                    PyValueError::new_err("division by 0 in this ramified field")
                })?;
                Ok($wrap(self.inner.mul(&oi)))
            }
            fn __rtruediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                let si = self
                    .inner
                    .inv()
                    .ok_or_else(|| PyValueError::new_err("division by 0 in this ramified field"))?;
                Ok($wrap($parse(other)?.mul(&si)))
            }
            fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
                matches!($parse(other), Ok(x) if x == self.inner)
            }
            fn __repr__(&self) -> String {
                format!("{:?}", self.inner)
            }
        }
    };
}

ramified_pyclass!(
    PyRamifiedQp2_4E2,
    "RamifiedQp2_4_E2",
    parse_ramified_qp2_4_e2,
    wrap_ramified_qp2_4_e2,
    Qp<2, 4>,
    2,
    PyQp2_4,
    parse_qp2_4,
    wrap_qp2_4,
    PyFp2,
    parse_fp2,
    wrap_fp2
);
ramified_pyclass!(
    PyRamifiedQp3_4E2,
    "RamifiedQp3_4_E2",
    parse_ramified_qp3_4_e2,
    wrap_ramified_qp3_4_e2,
    Qp<3, 4>,
    2,
    PyQp3_4,
    parse_qp3_4,
    wrap_qp3_4,
    PyFp3,
    parse_fp3,
    wrap_fp3
);
ramified_pyclass!(
    PyRamifiedQp5_4E2,
    "RamifiedQp5_4_E2",
    parse_ramified_qp5_4_e2,
    wrap_ramified_qp5_4_e2,
    Qp<5, 4>,
    2,
    PyQp5_4,
    parse_qp5_4,
    wrap_qp5_4,
    PyFp5,
    parse_fp5,
    wrap_fp5
);
ramified_pyclass!(
    PyRamifiedQp7_4E2,
    "RamifiedQp7_4_E2",
    parse_ramified_qp7_4_e2,
    wrap_ramified_qp7_4_e2,
    Qp<7, 4>,
    2,
    PyQp7_4,
    parse_qp7_4,
    wrap_qp7_4,
    PyFp7,
    parse_fp7,
    wrap_fp7
);
ramified_pyclass!(
    PyRamifiedQp11_4E2,
    "RamifiedQp11_4_E2",
    parse_ramified_qp11_4_e2,
    wrap_ramified_qp11_4_e2,
    Qp<11, 4>,
    2,
    PyQp11_4,
    parse_qp11_4,
    wrap_qp11_4,
    PyFp11,
    parse_fp11,
    wrap_fp11
);
ramified_pyclass!(
    PyRamifiedQp13_4E2,
    "RamifiedQp13_4_E2",
    parse_ramified_qp13_4_e2,
    wrap_ramified_qp13_4_e2,
    Qp<13, 4>,
    2,
    PyQp13_4,
    parse_qp13_4,
    wrap_qp13_4,
    PyFp13,
    parse_fp13,
    wrap_fp13
);
ramified_pyclass!(
    PyRamifiedQp2_4E3,
    "RamifiedQp2_4_E3",
    parse_ramified_qp2_4_e3,
    wrap_ramified_qp2_4_e3,
    Qp<2, 4>,
    3,
    PyQp2_4,
    parse_qp2_4,
    wrap_qp2_4,
    PyFp2,
    parse_fp2,
    wrap_fp2
);
ramified_pyclass!(
    PyRamifiedQp3_4E3,
    "RamifiedQp3_4_E3",
    parse_ramified_qp3_4_e3,
    wrap_ramified_qp3_4_e3,
    Qp<3, 4>,
    3,
    PyQp3_4,
    parse_qp3_4,
    wrap_qp3_4,
    PyFp3,
    parse_fp3,
    wrap_fp3
);
ramified_pyclass!(
    PyRamifiedQp5_4E3,
    "RamifiedQp5_4_E3",
    parse_ramified_qp5_4_e3,
    wrap_ramified_qp5_4_e3,
    Qp<5, 4>,
    3,
    PyQp5_4,
    parse_qp5_4,
    wrap_qp5_4,
    PyFp5,
    parse_fp5,
    wrap_fp5
);
ramified_pyclass!(
    PyRamifiedQp7_4E3,
    "RamifiedQp7_4_E3",
    parse_ramified_qp7_4_e3,
    wrap_ramified_qp7_4_e3,
    Qp<7, 4>,
    3,
    PyQp7_4,
    parse_qp7_4,
    wrap_qp7_4,
    PyFp7,
    parse_fp7,
    wrap_fp7
);
ramified_pyclass!(
    PyRamifiedQp11_4E3,
    "RamifiedQp11_4_E3",
    parse_ramified_qp11_4_e3,
    wrap_ramified_qp11_4_e3,
    Qp<11, 4>,
    3,
    PyQp11_4,
    parse_qp11_4,
    wrap_qp11_4,
    PyFp11,
    parse_fp11,
    wrap_fp11
);
ramified_pyclass!(
    PyRamifiedQp13_4E3,
    "RamifiedQp13_4_E3",
    parse_ramified_qp13_4_e3,
    wrap_ramified_qp13_4_e3,
    Qp<13, 4>,
    3,
    PyQp13_4,
    parse_qp13_4,
    wrap_qp13_4,
    PyFp13,
    parse_fp13,
    wrap_fp13
);

fn checked_gauss<S: Valued>(name: &str, num: Vec<S>, den: Vec<S>) -> PyResult<Gauss<S>> {
    if den.iter().all(Scalar::is_zero) {
        Err(PyValueError::new_err(format!("{name}: zero denominator")))
    } else {
        Ok(Gauss::new(num, den))
    }
}

macro_rules! gauss_pyclass {
    (
        $py:ident, $name:literal, $parse:ident, $wrap:ident,
        $base:ty, $base_py:ty, $base_parse:path, $base_wrap:path,
        $res_py:ty, $res_parse:path, $res_wrap:path
    ) => {
        #[pyclass(name = $name, module = "pleroma", from_py_object)]
        #[derive(Clone)]
        pub(crate) struct $py {
            inner: Gauss<$base>,
        }

        pub(crate) fn $parse(obj: &Bound<'_, PyAny>) -> PyResult<Gauss<$base>> {
            if let Ok(x) = obj.cast::<$py>() {
                return Ok(x.borrow().inner.clone());
            }
            if let Ok(s) = $base_parse(obj) {
                return Ok(Gauss::<$base>::from_base(s));
            }
            if let Ok((num, den)) = obj.extract::<(Vec<Bound<'_, PyAny>>, Vec<Bound<'_, PyAny>>)>()
            {
                let mut parsed_num = Vec::with_capacity(num.len());
                for coeff in &num {
                    parsed_num.push($base_parse(coeff)?);
                }
                let mut parsed_den = Vec::with_capacity(den.len());
                for coeff in &den {
                    parsed_den.push($base_parse(coeff)?);
                }
                return checked_gauss($name, parsed_num, parsed_den);
            }
            if let Ok(items) = obj.extract::<Vec<Bound<'_, PyAny>>>() {
                let mut coeffs = Vec::with_capacity(items.len());
                for item in &items {
                    coeffs.push($base_parse(item)?);
                }
                return Ok(Gauss::<$base>::new(coeffs, vec![<$base>::one()]));
            }
            Err(PyTypeError::new_err(concat!(
                "expected ",
                $name,
                ", base scalar, coefficient list, or (num, den) coefficient lists"
            )))
        }

        pub(crate) fn $wrap(inner: Gauss<$base>) -> $py {
            $py { inner }
        }

        #[pymethods]
        impl $py {
            #[new]
            #[pyo3(signature = (num, den=None))]
            fn new(num: Vec<Bound<'_, PyAny>>, den: Option<Vec<Bound<'_, PyAny>>>) -> PyResult<Self> {
                Self::from_coeffs(num, den)
            }
            #[staticmethod]
            fn zero() -> Self {
                $wrap(Gauss::<$base>::zero())
            }
            #[staticmethod]
            fn one() -> Self {
                $wrap(Gauss::<$base>::one())
            }
            #[staticmethod]
            fn from_base(s: &Bound<'_, PyAny>) -> PyResult<Self> {
                Ok($wrap(Gauss::<$base>::from_base($base_parse(s)?)))
            }
            #[staticmethod]
            fn teichmuller(residue: &Bound<'_, PyAny>) -> PyResult<Self> {
                Ok($wrap(<Gauss<$base> as ResidueField>::teichmuller(
                    $res_parse(residue)?,
                )))
            }
            #[staticmethod]
            fn t() -> Self {
                $wrap(Gauss::<$base>::t())
            }
            #[staticmethod]
            fn uniformizer() -> Self {
                $wrap(<Gauss<$base> as Valued>::uniformizer())
            }
            #[staticmethod]
            fn characteristic() -> u128 {
                Gauss::<$base>::characteristic()
            }
            #[staticmethod]
            #[pyo3(signature = (num, den=None))]
            fn from_coeffs(
                num: Vec<Bound<'_, PyAny>>,
                den: Option<Vec<Bound<'_, PyAny>>>,
            ) -> PyResult<Self> {
                let mut parsed_num = Vec::with_capacity(num.len());
                for coeff in &num {
                    parsed_num.push($base_parse(coeff)?);
                }
                let parsed_den = if let Some(den) = den {
                    let mut out = Vec::with_capacity(den.len());
                    for coeff in &den {
                        out.push($base_parse(coeff)?);
                    }
                    out
                } else {
                    vec![<$base>::one()]
                };
                checked_gauss($name, parsed_num, parsed_den).map($wrap)
            }
            fn parts(&self) -> (Vec<$base_py>, Vec<$base_py>) {
                let (num, den) = self.inner.parts();
                (
                    num.iter().cloned().map($base_wrap).collect(),
                    den.iter().cloned().map($base_wrap).collect(),
                )
            }
            fn valuation(&self) -> Option<i128> {
                self.inner.valuation()
            }
            fn residue(&self) -> Option<$res_py> {
                <Gauss<$base> as ResidueField>::residue(&self.inner).map($res_wrap)
            }
            fn residue_unit(&self) -> Option<$res_py> {
                <Gauss<$base> as ResidueField>::residue_unit(&self.inner).map($res_wrap)
            }
            fn angular_component(&self) -> Option<$res_py> {
                self.residue_unit()
            }
            fn is_integral(&self) -> bool {
                self.inner.is_integral()
            }
            fn is_zero(&self) -> bool {
                self.inner.is_zero()
            }
            fn inv(&self) -> PyResult<Self> {
                self.inner
                    .inv()
                    .map($wrap)
                    .ok_or_else(|| PyValueError::new_err("0 has no inverse in this Gauss field"))
            }
            fn __add__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                Ok($wrap(self.inner.add(&$parse(other)?)))
            }
            fn __radd__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                self.__add__(other)
            }
            fn __sub__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                Ok($wrap(self.inner.sub(&$parse(other)?)))
            }
            fn __rsub__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                Ok($wrap($parse(other)?.sub(&self.inner)))
            }
            fn __neg__(&self) -> Self {
                $wrap(self.inner.neg())
            }
            fn __mul__(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
                match $parse(other) {
                    Ok(o) => $wrap(self.inner.mul(&o)).into_py_any(py),
                    Err(_) => Ok(py.NotImplemented()),
                }
            }
            fn __rmul__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                Ok($wrap(self.inner.mul(&$parse(other)?)))
            }
            fn __truediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                let o = $parse(other)?;
                let oi = o
                    .inv()
                    .ok_or_else(|| PyValueError::new_err("division by 0 in this Gauss field"))?;
                Ok($wrap(self.inner.mul(&oi)))
            }
            fn __rtruediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
                let si = self
                    .inner
                    .inv()
                    .ok_or_else(|| PyValueError::new_err("division by 0 in this Gauss field"))?;
                Ok($wrap($parse(other)?.mul(&si)))
            }
            fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
                matches!($parse(other), Ok(x) if x == self.inner)
            }
            fn __repr__(&self) -> String {
                format!("{:?}", self.inner)
            }
        }
    };
}

gauss_pyclass!(
    PyGaussQp2_4,
    "GaussQp2_4",
    parse_gauss_qp2_4,
    wrap_gauss_qp2_4,
    Qp<2, 4>,
    PyQp2_4,
    parse_qp2_4,
    wrap_qp2_4,
    PyFp2RationalFunction,
    parse_fp2_rational_function,
    wrap_fp2_rational_function
);
gauss_pyclass!(
    PyGaussQp3_4,
    "GaussQp3_4",
    parse_gauss_qp3_4,
    wrap_gauss_qp3_4,
    Qp<3, 4>,
    PyQp3_4,
    parse_qp3_4,
    wrap_qp3_4,
    PyFp3RationalFunction,
    parse_fp3_rational_function,
    wrap_fp3_rational_function
);
gauss_pyclass!(
    PyGaussQp5_4,
    "GaussQp5_4",
    parse_gauss_qp5_4,
    wrap_gauss_qp5_4,
    Qp<5, 4>,
    PyQp5_4,
    parse_qp5_4,
    wrap_qp5_4,
    PyFp5RationalFunction,
    parse_fp5_rational_function,
    wrap_fp5_rational_function
);
gauss_pyclass!(
    PyGaussQp7_4,
    "GaussQp7_4",
    parse_gauss_qp7_4,
    wrap_gauss_qp7_4,
    Qp<7, 4>,
    PyQp7_4,
    parse_qp7_4,
    wrap_qp7_4,
    PyFp7RationalFunction,
    parse_fp7_rational_function,
    wrap_fp7_rational_function
);
gauss_pyclass!(
    PyGaussQp11_4,
    "GaussQp11_4",
    parse_gauss_qp11_4,
    wrap_gauss_qp11_4,
    Qp<11, 4>,
    PyQp11_4,
    parse_qp11_4,
    wrap_qp11_4,
    PyFp11RationalFunction,
    parse_fp11_rational_function,
    wrap_fp11_rational_function
);
gauss_pyclass!(
    PyGaussQp13_4,
    "GaussQp13_4",
    parse_gauss_qp13_4,
    wrap_gauss_qp13_4,
    Qp<13, 4>,
    PyQp13_4,
    parse_qp13_4,
    wrap_qp13_4,
    PyFp13RationalFunction,
    parse_fp13_rational_function,
    wrap_fp13_rational_function
);

fn poly_min_coeff_valuation<S: Valued>(coeffs: Vec<S>) -> Option<i128> {
    Poly::new(coeffs).min_coeff_valuation()
}

/// Rust-name Gauss valuation of a valued-scalar polynomial's coefficient list.
///
/// Python does not carry a generic `Poly<S: Valued>` type parameter at runtime,
/// so the dispatcher requires a homogeneous typed coefficient list. The empty
/// list is the zero polynomial and returns `None`.
#[pyfunction]
fn min_coeff_valuation(coeffs: Vec<Bound<'_, PyAny>>) -> PyResult<Option<i128>> {
    if coeffs.is_empty() {
        return Ok(None);
    }

    macro_rules! try_valued_coeffs {
        ($py:ty) => {{
            let mut parsed = Vec::with_capacity(coeffs.len());
            let mut all_match = true;
            for coeff in &coeffs {
                if let Ok(value) = coeff.cast::<$py>() {
                    parsed.push(value.borrow().inner.clone());
                } else {
                    all_match = false;
                    break;
                }
            }
            if all_match {
                return Ok(poly_min_coeff_valuation(parsed));
            }
        }};
    }

    try_valued_coeffs!(PyQp2_4);
    try_valued_coeffs!(PyQp3_4);
    try_valued_coeffs!(PyQp5_4);
    try_valued_coeffs!(PyQp7_4);
    try_valued_coeffs!(PyQp11_4);
    try_valued_coeffs!(PyQp13_4);

    try_valued_coeffs!(PyQq2_4_2);
    try_valued_coeffs!(PyQq2_4_3);
    try_valued_coeffs!(PyQq2_4_4);
    try_valued_coeffs!(PyQq3_4_2);
    try_valued_coeffs!(PyQq5_4_2);
    try_valued_coeffs!(PyQq3_4_3);

    try_valued_coeffs!(PyLaurentRational6);
    try_valued_coeffs!(PyLaurentFp3_6);
    try_valued_coeffs!(PyLaurentFp5_6);
    try_valued_coeffs!(PyLaurentFp7_6);
    try_valued_coeffs!(PyLaurentFp11_6);
    try_valued_coeffs!(PyLaurentFp13_6);
    try_valued_coeffs!(PyLaurentF9_6);
    try_valued_coeffs!(PyLaurentF25_6);
    try_valued_coeffs!(PyLaurentF27_6);

    try_valued_coeffs!(PyRamifiedQp2_4E2);
    try_valued_coeffs!(PyRamifiedQp3_4E2);
    try_valued_coeffs!(PyRamifiedQp5_4E2);
    try_valued_coeffs!(PyRamifiedQp7_4E2);
    try_valued_coeffs!(PyRamifiedQp11_4E2);
    try_valued_coeffs!(PyRamifiedQp13_4E2);
    try_valued_coeffs!(PyRamifiedQp2_4E3);
    try_valued_coeffs!(PyRamifiedQp3_4E3);
    try_valued_coeffs!(PyRamifiedQp5_4E3);
    try_valued_coeffs!(PyRamifiedQp7_4E3);
    try_valued_coeffs!(PyRamifiedQp11_4E3);
    try_valued_coeffs!(PyRamifiedQp13_4E3);

    try_valued_coeffs!(PyGaussQp2_4);
    try_valued_coeffs!(PyGaussQp3_4);
    try_valued_coeffs!(PyGaussQp5_4);
    try_valued_coeffs!(PyGaussQp7_4);
    try_valued_coeffs!(PyGaussQp11_4);
    try_valued_coeffs!(PyGaussQp13_4);

    Err(PyTypeError::new_err(
        "min_coeff_valuation expects a homogeneous list of typed Qp/Qq/Laurent/Ramified/Gauss coefficients",
    ))
}

#[pyclass(name = "Rational", module = "pleroma", from_py_object)]
#[derive(Clone)]
pub(crate) struct PyRational {
    inner: Rational,
}

#[pymethods]
impl PyRational {
    #[new]
    #[pyo3(signature = (num, den=1))]
    fn new(num: i128, den: i128) -> PyResult<Self> {
        Rational::try_new(num, den)
            .map(|inner| PyRational { inner })
            .ok_or_else(|| {
                PyValueError::new_err(
                    "Rational needs a nonzero denominator and i128-safe normalization",
                )
            })
    }
    #[staticmethod]
    #[pyo3(signature = (num, den=1))]
    fn try_new(num: i128, den: i128) -> Option<Self> {
        Rational::try_new(num, den).map(|inner| PyRational { inner })
    }
    #[staticmethod]
    fn zero() -> Self {
        wrap_rational(Rational::zero())
    }
    #[staticmethod]
    fn one() -> Self {
        wrap_rational(Rational::one())
    }
    #[staticmethod]
    fn integer(n: i128) -> Self {
        wrap_rational(Rational::int(n))
    }
    #[staticmethod]
    fn characteristic() -> u128 {
        Rational::characteristic()
    }
    #[getter]
    fn numerator(&self) -> i128 {
        self.inner.numer()
    }
    #[getter]
    fn denominator(&self) -> i128 {
        self.inner.denom()
    }
    fn numer(&self) -> i128 {
        self.inner.numer()
    }
    fn denom(&self) -> i128 {
        self.inner.denom()
    }
    fn is_integer(&self) -> bool {
        self.inner.is_integer()
    }
    fn is_zero(&self) -> bool {
        self.inner.is_zero()
    }
    /// Sign of this rational: `-1`, `0`, or `1`.
    fn sign(&self) -> i8 {
        ordering_to_i8(self.inner.sign())
    }
    fn floor(&self) -> i128 {
        self.inner.floor()
    }
    fn __add__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyRational> {
        Ok(PyRational {
            inner: self.inner.add(&parse_rational(other)?),
        })
    }
    fn __radd__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyRational> {
        self.__add__(other)
    }
    fn __sub__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyRational> {
        Ok(PyRational {
            inner: self.inner.sub(&parse_rational(other)?),
        })
    }
    fn __rsub__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyRational> {
        Ok(PyRational {
            inner: parse_rational(other)?.sub(&self.inner),
        })
    }
    fn __mul__(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        match parse_rational(other) {
            Ok(o) => PyRational {
                inner: self.inner.mul(&o),
            }
            .into_py_any(py),
            Err(_) => Ok(py.NotImplemented()),
        }
    }
    fn __rmul__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyRational> {
        Ok(PyRational {
            inner: self.inner.mul(&parse_rational(other)?),
        })
    }
    fn __neg__(&self) -> PyRational {
        PyRational {
            inner: self.inner.neg(),
        }
    }
    fn inv(&self) -> PyResult<PyRational> {
        self.inner
            .inv()
            .map(|inner| PyRational { inner })
            .ok_or_else(|| PyValueError::new_err("0 has no inverse"))
    }
    fn __truediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyRational> {
        let o = parse_rational(other)?;
        let oi = o
            .inv()
            .ok_or_else(|| PyValueError::new_err("division by zero"))?;
        Ok(PyRational {
            inner: self.inner.mul(&oi),
        })
    }
    fn __rtruediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyRational> {
        let si = self
            .inner
            .inv()
            .ok_or_else(|| PyValueError::new_err("division by zero"))?;
        Ok(PyRational {
            inner: parse_rational(other)?.mul(&si),
        })
    }
    fn sqrt(&self) -> Option<PyRational> {
        self.inner.sqrt().map(|inner| PyRational { inner })
    }
    fn nth_root(&self, k: u128) -> Option<PyRational> {
        self.inner.nth_root(k).map(|inner| PyRational { inner })
    }
    fn __richcmp__(&self, other: &Bound<'_, PyAny>, op: CompareOp) -> PyResult<bool> {
        Ok(op.matches(self.inner.cmp(&parse_rational(other)?)))
    }
    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}

pub(crate) fn parse_rational(obj: &Bound<'_, PyAny>) -> PyResult<Rational> {
    if let Ok(q) = obj.cast::<PyRational>() {
        return Ok(q.borrow().inner.clone());
    }
    if let Ok(i) = obj.cast::<PyInteger>() {
        return Ok(Rational::int(i.borrow().inner.0));
    }
    if let Ok((num, den)) = obj.extract::<(i128, i128)>() {
        return Rational::try_new(num, den).ok_or_else(|| {
            PyValueError::new_err("rational tuple has zero denominator or overflowed i128")
        });
    }
    if let Ok(v) = obj.extract::<i128>() {
        return Ok(Rational::int(v));
    }
    Err(PyTypeError::new_err(
        "expected Rational, Integer, int, or (num, den) tuple",
    ))
}

#[pyclass(name = "Surreal", module = "pleroma", from_py_object)]
#[derive(Clone)]
pub(crate) struct PySurreal {
    inner: Surreal,
}

#[pymethods]
impl PySurreal {
    #[staticmethod]
    fn zero() -> PySurreal {
        PySurreal {
            inner: Surreal::zero(),
        }
    }
    #[staticmethod]
    fn one() -> PySurreal {
        PySurreal {
            inner: Surreal::one(),
        }
    }
    #[staticmethod]
    fn from_int(n: i128) -> PySurreal {
        PySurreal {
            inner: Surreal::from_int(n),
        }
    }
    #[staticmethod]
    fn from_rational(num: i128, den: i128) -> PyResult<PySurreal> {
        Ok(PySurreal {
            inner: Surreal::from_rational(rational_from_pair(num, den)?),
        })
    }
    #[staticmethod]
    fn characteristic() -> u128 {
        Surreal::characteristic()
    }
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
    fn is_zero(&self) -> bool {
        self.inner.is_zero()
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
    fn __pow__(&self, n: u128, modulo: Option<&Bound<'_, PyAny>>) -> PyResult<PySurreal> {
        if modulo.is_some() {
            return Err(PyValueError::new_err(
                "surreal exponentiation does not take a modulus",
            ));
        }
        let mut acc = Surreal::one();
        for _ in 0..n {
            acc = acc.mul(&self.inner);
        }
        Ok(PySurreal { inner: acc })
    }
    fn __richcmp__(&self, other: &Bound<'_, PyAny>, op: CompareOp) -> PyResult<bool> {
        Ok(op.matches(self.inner.cmp(&parse_surreal(other)?)))
    }
    /// Sign of the leading term: `-1`, `0`, or `1`.
    fn sign(&self) -> i8 {
        ordering_to_i8(self.inner.sign())
    }
    /// The finite-support Hahn/CNF terms `(exponent, coefficient)`.
    #[getter]
    fn terms(&self) -> Vec<(PySurreal, PyRational)> {
        self.inner
            .terms()
            .iter()
            .map(|(e, c)| (PySurreal::from_inner(e.clone()), wrap_rational(c.clone())))
            .collect()
    }
    /// This surreal as a finite rational, if it has no `ω`-term.
    fn as_rational(&self) -> Option<PyRational> {
        self.inner.as_rational().map(wrap_rational)
    }
    /// This surreal as a dyadic rational `num / 2^k`, if it is dyadic.
    fn as_dyadic(&self) -> Option<(i128, u128)> {
        self.inner.as_dyadic()
    }
    /// A single Hahn monomial `coeff * ω^exp`.
    #[staticmethod]
    fn monomial(exp: &Bound<'_, PyAny>, coeff: &Bound<'_, PyAny>) -> PyResult<PySurreal> {
        Ok(PySurreal {
            inner: Surreal::monomial(parse_surreal(exp)?, parse_rational(coeff)?),
        })
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
    /// The surreal reconstructed from transfinite sign-expansion runs
    /// `(sign, length)`, or `None` outside the represented subclass.
    #[staticmethod]
    fn from_transfinite_sign_expansion(runs: Vec<(bool, PyOrdinal)>) -> Option<PySurreal> {
        let se = SignExpansion::from_runs(
            runs.into_iter()
                .map(|(sign, len)| (sign, len.as_ordinal().clone()))
                .collect(),
        );
        Surreal::from_transfinite_sign_expansion(&se).map(|inner| PySurreal { inner })
    }
    /// The surreal reconstructed from a `SignExpansion` value object.
    #[staticmethod]
    fn from_sign_expansion_record(se: &PySignExpansion) -> Option<PySurreal> {
        Surreal::from_transfinite_sign_expansion(&se.inner).map(|inner| PySurreal { inner })
    }
    /// Interpret this surreal as an ordinal, if it is ordinal-valued.
    fn as_ordinal(&self) -> Option<PyOrdinal> {
        self.inner.as_ordinal().map(PyOrdinal::from_inner)
    }
    /// Embed an ordinal as the corresponding surreal ordinal.
    #[staticmethod]
    fn from_ordinal(o: &PyOrdinal) -> PySurreal {
        PySurreal {
            inner: Surreal::from_ordinal(o.as_ordinal()),
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
    fn sqrt_to_terms(&self, n: usize) -> Option<PySurreal> {
        self.inner.sqrt_to_terms(n).map(|inner| PySurreal { inner })
    }
    /// The **truncated real `k`-th root** to `n` leading terms (same ℚ-power scope).
    fn nth_root_to_terms(&self, k: u128, n: usize) -> Option<PySurreal> {
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
    /// `√2` is `None`. The exact companion to [`sqrt_to_terms`](Self::sqrt_to_terms).
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
    /// The transfinite sign expansion as a `SignExpansion` value object.
    fn transfinite_sign_expansion_record(&self) -> Option<PySignExpansion> {
        self.inner
            .transfinite_sign_expansion()
            .map(PySignExpansion::from_inner)
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
    #[staticmethod]
    fn zero() -> PySurcomplex {
        PySurcomplex {
            inner: Surcomplex::zero(),
        }
    }
    #[staticmethod]
    fn one() -> PySurcomplex {
        PySurcomplex {
            inner: Surcomplex::one(),
        }
    }
    #[staticmethod]
    fn characteristic() -> u128 {
        Surcomplex::<Surreal>::characteristic()
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
    fn __rsub__(&self, other: &Bound<'_, PyAny>) -> PyResult<PySurcomplex> {
        Ok(PySurcomplex {
            inner: parse_surcomplex(other)?.sub(&self.inner),
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
    fn is_zero(&self) -> bool {
        self.inner.is_zero()
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
    fn __rtruediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<PySurcomplex> {
        let si = self
            .inner
            .inv()
            .ok_or_else(|| PyValueError::new_err("divisor has no representable inverse"))?;
        Ok(PySurcomplex {
            inner: parse_surcomplex(other)?.mul(&si),
        })
    }
    fn __pow__(&self, n: u128, modulo: Option<&Bound<'_, PyAny>>) -> PyResult<PySurcomplex> {
        if modulo.is_some() {
            return Err(PyValueError::new_err(
                "surcomplex exponentiation does not take a modulus",
            ));
        }
        let mut acc = Surcomplex::<Surreal>::one();
        for _ in 0..n {
            acc = acc.mul(&self.inner);
        }
        Ok(PySurcomplex { inner: acc })
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
    #[staticmethod]
    fn zero() -> Self {
        wrap_integer(Integer::zero())
    }
    #[staticmethod]
    fn one() -> Self {
        wrap_integer(Integer::one())
    }
    #[staticmethod]
    fn characteristic() -> u128 {
        Integer::characteristic()
    }
    #[getter]
    fn value(&self) -> i128 {
        self.inner.0
    }
    fn is_zero(&self) -> bool {
        self.inner.is_zero()
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
    fn __rsub__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyInteger> {
        Ok(PyInteger {
            inner: parse_integer(other)?.sub(&self.inner),
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
    fn inv(&self) -> PyResult<PyInteger> {
        self.inner
            .inv()
            .map(wrap_integer)
            .ok_or_else(|| PyValueError::new_err("Z is a ring: only ±1 are invertible"))
    }
    fn __truediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyInteger> {
        let rhs = parse_integer(other)?;
        let rinv = rhs
            .inv()
            .ok_or_else(|| PyValueError::new_err("integer divisor is not a unit"))?;
        Ok(wrap_integer(self.inner.mul(&rinv)))
    }
    fn __rtruediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyInteger> {
        let si = self
            .inner
            .inv()
            .ok_or_else(|| PyValueError::new_err("integer divisor is not a unit"))?;
        Ok(wrap_integer(parse_integer(other)?.mul(&si)))
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

pub(crate) fn wrap_rational(q: Rational) -> PyRational {
    PyRational { inner: q }
}

pub(crate) fn wrap_nimber(n: Nimber) -> PyNimber {
    PyNimber { inner: n }
}

pub(crate) fn wrap_ordinal(o: Ordinal) -> PyOrdinal {
    PyOrdinal { inner: o }
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
    #[staticmethod]
    fn zero() -> Self {
        wrap_omnific(Omnific::zero())
    }
    #[staticmethod]
    fn one() -> Self {
        wrap_omnific(Omnific::one())
    }
    #[staticmethod]
    fn omega() -> Self {
        wrap_omnific(Omnific::omega())
    }
    #[staticmethod]
    fn from_surreal(s: &PySurreal) -> PyResult<Self> {
        Omnific::from_surreal(s.inner.clone())
            .map(wrap_omnific)
            .ok_or_else(|| PyValueError::new_err("surreal is not an omnific integer"))
    }
    #[staticmethod]
    fn floor(s: &PySurreal) -> Self {
        wrap_omnific(Omnific::floor(&s.inner))
    }
    #[staticmethod]
    fn characteristic() -> u128 {
        Omnific::characteristic()
    }
    /// The underlying surreal value.
    fn surreal(&self) -> PySurreal {
        PySurreal {
            inner: self.inner.inner().clone(),
        }
    }
    fn is_zero(&self) -> bool {
        self.inner.is_zero()
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
    fn __rsub__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyOmnific> {
        Ok(PyOmnific {
            inner: parse_omnific(other)?.sub(&self.inner),
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
    fn __truediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyOmnific> {
        let rhs = parse_omnific(other)?;
        let rinv = rhs
            .inv()
            .ok_or_else(|| PyValueError::new_err("omnific divisor is not a unit"))?;
        Ok(wrap_omnific(self.inner.mul(&rinv)))
    }
    fn __rtruediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyOmnific> {
        let si = self
            .inner
            .inv()
            .ok_or_else(|| PyValueError::new_err("omnific divisor is not a unit"))?;
        Ok(wrap_omnific(parse_omnific(other)?.mul(&si)))
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

/// Rust-name finite-place precision policy used by [`Adele`].
#[pyfunction]
fn adele_prec(p: u128) -> PyResult<u128> {
    adele_precision_for_prime(p)
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
    #[staticmethod]
    fn from_i128(p: u128, k: u128, value: i128) -> PyResult<Self> {
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
    fn __rtruediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyLocalQp> {
        let si = self
            .inner
            .inv()
            .ok_or_else(|| PyValueError::new_err("division by 0 in LocalQp"))?;
        let lhs = parse_local_qp_in_world(other, self.inner.prime(), self.inner.precision())?;
        Ok(PyLocalQp {
            inner: lhs.mul(&si),
        })
    }
    fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
        matches!(other.cast::<PyLocalQp>(), Ok(x) if x.borrow().inner == self.inner)
    }
    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}

pub(crate) fn parse_adele(obj: &Bound<'_, PyAny>) -> PyResult<Adele> {
    if let Ok(a) = obj.cast::<PyAdele>() {
        return Ok(a.borrow().inner.clone());
    }
    if let Ok(v) = obj.extract::<i128>() {
        return Ok(Adele::from_rational(&Rational::int(v)));
    }
    Err(PyTypeError::new_err("expected Adele or int"))
}

pub(crate) fn wrap_adele(inner: Adele) -> PyAdele {
    PyAdele { inner }
}

#[pyclass(name = "AdelePlace", module = "pleroma", from_py_object)]
#[derive(Clone, Copy)]
struct PyAdelePlace {
    inner: AdelePlace,
}

fn adele_place_name(place: AdelePlace) -> String {
    match place {
        AdelePlace::Real => "R".to_string(),
        AdelePlace::Prime(p) => format!("Q_{p}"),
    }
}

fn wrap_adele_place(inner: AdelePlace) -> PyAdelePlace {
    PyAdelePlace { inner }
}

fn parse_adele_place(obj: &Bound<'_, PyAny>) -> PyResult<AdelePlace> {
    if let Ok(place) = obj.cast::<PyAdelePlace>() {
        return Ok(place.borrow().inner);
    }
    Err(PyTypeError::new_err("expected AdelePlace"))
}

#[pymethods]
impl PyAdelePlace {
    #[staticmethod]
    fn real() -> Self {
        wrap_adele_place(AdelePlace::Real)
    }

    #[staticmethod]
    fn prime(p: u128) -> Self {
        wrap_adele_place(AdelePlace::Prime(p))
    }

    #[getter]
    fn name(&self) -> String {
        adele_place_name(self.inner)
    }

    #[getter]
    fn is_real(&self) -> bool {
        self.inner == AdelePlace::Real
    }

    #[getter]
    fn is_prime(&self) -> bool {
        matches!(self.inner, AdelePlace::Prime(_))
    }

    #[getter]
    fn prime_value(&self) -> Option<u128> {
        match self.inner {
            AdelePlace::Real => None,
            AdelePlace::Prime(p) => Some(p),
        }
    }

    fn __str__(&self) -> String {
        adele_place_name(self.inner)
    }

    fn __repr__(&self) -> String {
        match self.inner {
            AdelePlace::Real => "AdelePlace.Real".to_string(),
            AdelePlace::Prime(p) => format!("AdelePlace.Prime({p})"),
        }
    }

    fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
        matches!(other.cast::<PyAdelePlace>(), Ok(place) if place.borrow().inner == self.inner)
    }
}

#[pyclass(name = "Adele", module = "pleroma", from_py_object)]
#[derive(Clone)]
pub(crate) struct PyAdele {
    inner: Adele,
}

#[pymethods]
impl PyAdele {
    #[staticmethod]
    fn from_rational(num: i128, den: i128) -> PyResult<PyAdele> {
        let q = rational_from_pair(num, den)?;
        Ok(PyAdele {
            inner: Adele::from_rational(&q),
        })
    }
    #[staticmethod]
    fn zero() -> PyAdele {
        PyAdele {
            inner: Adele::zero(),
        }
    }
    #[staticmethod]
    fn one() -> PyAdele {
        PyAdele {
            inner: Adele::one(),
        }
    }
    #[staticmethod]
    fn characteristic() -> u128 {
        Adele::characteristic()
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
    fn is_zero(&self) -> bool {
        self.inner.is_zero()
    }
    fn is_idele(&self) -> bool {
        self.inner.is_idele()
    }
    fn is_integral(&self) -> bool {
        self.inner.is_integral()
    }
    fn absolute_value_at(&self, place: &Bound<'_, PyAny>) -> PyResult<(i128, i128)> {
        let place = parse_adele_place(place)?;
        if let AdelePlace::Prime(p) = place {
            adele_precision_for_prime(p)?;
        }
        Ok(rational_pair(&self.inner.absolute_value_at(place)))
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
    fn __rtruediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyAdele> {
        let si = self
            .inner
            .inv()
            .ok_or_else(|| PyValueError::new_err("Adele divisor is not an idele"))?;
        Ok(PyAdele {
            inner: parse_adele(other)?.mul(&si),
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
            fn int(n: i128) -> $py {
                $py {
                    inner: Tropical::<$conv>::int(n),
                }
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
            fn is_zero(&self) -> bool {
                self.inner.is_zero()
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

/// Whether a represented surreal belongs to Oz, the omnific-integer subring.
#[pyfunction]
fn is_omnific_integer(x: &Bound<'_, PyAny>) -> PyResult<bool> {
    Ok(crate::scalar::is_omnific_integer(&parse_surreal(x)?))
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
fn surreal(n: i128) -> PySurreal {
    PySurreal {
        inner: Surreal::from_int(n),
    }
}
// ---------------------------------------------------------------------------
// Nim field operations (the Artin–Schreier ↔ Arf bridge)
// ---------------------------------------------------------------------------

/// Nim addition, i.e. xor on the represented finite nimbers.
#[pyfunction]
fn nim_add(a: u128, b: u128) -> u128 {
    crate::scalar::nim_add(a, b)
}

/// Conway nim multiplication on the represented `F_{2^128}` backend.
#[pyfunction]
fn nim_mul(a: u128, b: u128) -> u128 {
    crate::scalar::nim_mul(a, b)
}

/// Exponentiation by repeated nim multiplication.
#[pyfunction]
fn nim_pow(base: u128, exp: u128) -> u128 {
    crate::scalar::nim_pow(base, exp)
}

/// Frobenius square `x -> x^2` in the nim field.
#[pyfunction]
fn nim_square(x: u128) -> u128 {
    crate::scalar::nim_square(x)
}

/// Iterate Frobenius `k` times: `x -> x^(2^k)`.
#[pyfunction]
fn nim_frobenius_iter(x: u128, k: usize) -> u128 {
    crate::scalar::nim_frobenius_iter(x, k)
}

/// Multiplicative inverse in the nim field (`None` for `*0`).
#[pyfunction]
fn nim_inv(x: u128) -> Option<u128> {
    crate::scalar::nim_inv(x)
}

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
fn nim_degree(x: u128) -> u128 {
    crate::scalar::nim_degree(x)
}

/// The Galois conjugates `x, x², x⁴, …` over F₂.
#[pyfunction]
fn nim_conjugates(x: u128) -> Vec<u128> {
    crate::scalar::nim_conjugates(x)
}

/// Minimal polynomial of `x` over F₂: coefficients `{0,1}` from the constant up.
#[pyfunction]
fn nim_min_poly(x: u128) -> Vec<u128> {
    crate::scalar::nim_min_poly(x)
}

/// Relative trace `Tr_{F_{2^m}/F_{2^e}}(x)` (the `e=1` case is `nim_trace`).
#[pyfunction]
fn nim_relative_trace(x: u128, m: u128, e: u128) -> u128 {
    crate::scalar::nim_relative_trace(x, m, e)
}

/// Relative norm `N_{F_{2^m}/F_{2^e}}(x)` (norm to the prime field is trivial).
#[pyfunction]
fn nim_relative_norm(x: u128, m: u128, e: u128) -> u128 {
    crate::scalar::nim_relative_norm(x, m, e)
}

/// Multiplicative order of `x` in F_{2^128}* (`None` for `*0`).
#[pyfunction]
fn nim_multiplicative_order(x: u128) -> Option<u128> {
    crate::scalar::nim_multiplicative_order(x)
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

#[pyclass(name = "SignExpansion", module = "pleroma", from_py_object)]
#[derive(Clone)]
struct PySignExpansion {
    inner: SignExpansion,
}

impl PySignExpansion {
    fn from_inner(inner: SignExpansion) -> PySignExpansion {
        PySignExpansion { inner }
    }
}

#[pymethods]
impl PySignExpansion {
    #[new]
    fn new(runs: Vec<(bool, PyOrdinal)>) -> Self {
        Self::from_runs(runs)
    }

    /// Build and normalize from `(sign, ordinal_length)` runs.
    #[staticmethod]
    fn from_runs(runs: Vec<(bool, PyOrdinal)>) -> Self {
        PySignExpansion {
            inner: SignExpansion::from_runs(
                runs.into_iter()
                    .map(|(sign, len)| (sign, len.inner))
                    .collect(),
            ),
        }
    }

    /// Run-length-encode a finite sign sequence (`True = +`).
    #[staticmethod]
    fn from_finite(signs: Vec<bool>) -> Self {
        PySignExpansion {
            inner: SignExpansion::from_finite(&signs),
        }
    }

    /// The normalized runs `(sign, length)`, left to right.
    fn runs(&self) -> Vec<(bool, PyOrdinal)> {
        self.inner
            .runs()
            .iter()
            .map(|(sign, len)| (*sign, PyOrdinal { inner: len.clone() }))
            .collect()
    }

    /// The total ordinal length, i.e. the birthday.
    fn length(&self) -> PyOrdinal {
        PyOrdinal {
            inner: self.inner.length(),
        }
    }

    /// The flat finite sign sequence, if no run has transfinite length.
    fn as_finite(&self) -> Option<Vec<bool>> {
        self.inner.as_finite()
    }

    fn __repr__(&self) -> String {
        format!("SignExpansion({:?})", self.inner.runs())
    }
}

pub(crate) fn parse_ordinal(obj: &Bound<'_, PyAny>) -> PyResult<Ordinal> {
    if let Ok(o) = obj.cast::<PyOrdinal>() {
        return Ok(o.borrow().inner.clone());
    }
    if let Ok(v) = obj.extract::<u128>() {
        return Ok(Ordinal::from_u128(v));
    }
    Err(PyTypeError::new_err("expected Ordinal or non-negative int"))
}

#[pymethods]
impl PyOrdinal {
    #[new]
    fn new(n: u128) -> Self {
        PyOrdinal {
            inner: Ordinal::from_u128(n),
        }
    }
    /// The ordinal/nimber zero.
    #[staticmethod]
    fn zero() -> PyOrdinal {
        PyOrdinal {
            inner: Ordinal::zero(),
        }
    }
    /// The ordinal/nimber one.
    #[staticmethod]
    fn one() -> PyOrdinal {
        PyOrdinal {
            inner: Ordinal::one(),
        }
    }
    #[staticmethod]
    fn characteristic() -> u128 {
        Ordinal::characteristic()
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
    /// Build `ω²·c₂ + ω·c₁ + c₀` from `[c₀, c₁, c₂]`.
    #[staticmethod]
    fn from_omega3_coeffs(coeffs: Vec<u128>) -> PyResult<PyOrdinal> {
        let coeffs: [u128; 3] = coeffs.try_into().map_err(|_| {
            PyValueError::new_err(
                "Ordinal.from_omega3_coeffs expects exactly three coefficients [c0, c1, c2]",
            )
        })?;
        Ok(PyOrdinal {
            inner: Ordinal::from_omega3_coeffs(coeffs),
        })
    }
    fn __add__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyOrdinal> {
        Ok(PyOrdinal {
            inner: self.inner.add(&parse_ordinal(other)?),
        })
    }
    fn __radd__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyOrdinal> {
        self.__add__(other)
    }
    fn __sub__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyOrdinal> {
        Ok(PyOrdinal {
            inner: self.inner.sub(&parse_ordinal(other)?),
        })
    }
    fn __rsub__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyOrdinal> {
        Ok(PyOrdinal {
            inner: parse_ordinal(other)?.sub(&self.inner),
        })
    }
    fn __neg__(&self) -> PyOrdinal {
        PyOrdinal {
            inner: self.inner.neg(),
        }
    }
    /// Checked nim-field scalar multiplication.
    fn __mul__(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        match parse_ordinal(other) {
            Ok(rhs) => self
                .inner
                .nim_mul(&rhs)
                .map(|inner| PyOrdinal { inner })
                .ok_or_else(|| {
                    PyValueError::new_err(
                        "ordinal nim product escaped the source-verified Kummer boundary",
                    )
                })?
                .into_py_any(py),
            Err(_) => Ok(py.NotImplemented()),
        }
    }
    fn __rmul__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyOrdinal> {
        parse_ordinal(other)?
            .nim_mul(&self.inner)
            .map(|inner| PyOrdinal { inner })
            .ok_or_else(|| {
                PyValueError::new_err(
                    "ordinal nim product escaped the source-verified Kummer boundary",
                )
            })
    }
    fn inv(&self) -> PyResult<PyOrdinal> {
        self.inner
            .checked_inv()
            .map(|inner| PyOrdinal { inner })
            .ok_or_else(|| PyValueError::new_err("ordinal has no represented inverse"))
    }
    fn __truediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyOrdinal> {
        let rhs = parse_ordinal(other)?;
        let rinv = rhs
            .checked_inv()
            .ok_or_else(|| PyValueError::new_err("ordinal divisor has no represented inverse"))?;
        self.inner
            .nim_mul(&rinv)
            .map(|inner| PyOrdinal { inner })
            .ok_or_else(|| {
                PyValueError::new_err(
                    "ordinal quotient escaped the source-verified Kummer boundary",
                )
            })
    }
    fn __rtruediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyOrdinal> {
        let si = self
            .inner
            .checked_inv()
            .ok_or_else(|| PyValueError::new_err("ordinal divisor has no represented inverse"))?;
        parse_ordinal(other)?
            .nim_mul(&si)
            .map(|inner| PyOrdinal { inner })
            .ok_or_else(|| {
                PyValueError::new_err(
                    "ordinal quotient escaped the source-verified Kummer boundary",
                )
            })
    }
    /// Nim-addition (complete and exact): XOR of like-`ω`-power coefficients.
    fn nim_add(&self, other: &PyOrdinal) -> PyOrdinal {
        PyOrdinal {
            inner: self.inner.nim_add(&other.inner),
        }
    }
    /// Nim-multiplication (partial): exact on the verified Kummer window,
    /// including finite operands and staged transfinite products such as
    /// `ω ⊗ ω = ω²`; `None` beyond that represented boundary.
    fn nim_mul(&self, other: &PyOrdinal) -> Option<PyOrdinal> {
        self.inner
            .nim_mul(&other.inner)
            .map(|o| PyOrdinal { inner: o })
    }
    /// Alias for the represented inverse boundary.
    fn checked_inv(&self) -> Option<PyOrdinal> {
        self.inner.checked_inv().map(|inner| PyOrdinal { inner })
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
    /// Coefficients `[c₀, c₁, c₂]` if this ordinal is below `ω³`.
    fn as_below_omega3(&self) -> Option<Vec<u128>> {
        self.inner.as_below_omega3().map(Vec::from)
    }
    /// The CNF terms `(exponent, coefficient)`.
    #[getter]
    fn terms(&self) -> Vec<(PyOrdinal, u128)> {
        self.inner
            .terms()
            .iter()
            .map(|(e, c)| (PyOrdinal::from_inner(e.clone()), *c))
            .collect()
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
    m.add_class::<PyNimberPoly>()?;
    m.add_class::<PyNimberRationalFunction>()?;
    m.add_class::<PyFp2>()?;
    m.add_class::<PyFp3>()?;
    m.add_class::<PyFp5>()?;
    m.add_class::<PyFp7>()?;
    m.add_class::<PyFp11>()?;
    m.add_class::<PyFp13>()?;
    m.add_class::<PyReductionPolynomialKind>()?;
    m.add_class::<PyF4>()?;
    m.add_class::<PyF8>()?;
    m.add_class::<PyF16>()?;
    m.add_class::<PyF9>()?;
    m.add_class::<PyF25>()?;
    m.add_class::<PyF27>()?;
    m.add_class::<PyFp2Poly>()?;
    m.add_class::<PyFp2RationalFunction>()?;
    m.add_class::<PyFp3Poly>()?;
    m.add_class::<PyFp3RationalFunction>()?;
    m.add_class::<PyFp5Poly>()?;
    m.add_class::<PyFp5RationalFunction>()?;
    m.add_class::<PyFp7Poly>()?;
    m.add_class::<PyFp7RationalFunction>()?;
    m.add_class::<PyFp11Poly>()?;
    m.add_class::<PyFp11RationalFunction>()?;
    m.add_class::<PyFp13Poly>()?;
    m.add_class::<PyFp13RationalFunction>()?;
    m.add_class::<PyZp2_4>()?;
    m.add_class::<PyZp3_4>()?;
    m.add_class::<PyZp5_4>()?;
    m.add_class::<PyZp7_4>()?;
    m.add_class::<PyZp11_4>()?;
    m.add_class::<PyZp13_4>()?;
    m.add_class::<PyQp2_4>()?;
    m.add_class::<PyQp3_4>()?;
    m.add_class::<PyQp5_4>()?;
    m.add_class::<PyQp7_4>()?;
    m.add_class::<PyQp11_4>()?;
    m.add_class::<PyQp13_4>()?;
    m.add_class::<PyWittVec2_4_2>()?;
    m.add_class::<PyWittVec2_4_3>()?;
    m.add_class::<PyWittVec2_4_4>()?;
    m.add_class::<PyWittVec3_4_2>()?;
    m.add_class::<PyWittVec5_4_2>()?;
    m.add_class::<PyWittVec3_4_3>()?;
    m.add_class::<PyQq2_4_2>()?;
    m.add_class::<PyQq2_4_3>()?;
    m.add_class::<PyQq2_4_4>()?;
    m.add_class::<PyQq3_4_2>()?;
    m.add_class::<PyQq5_4_2>()?;
    m.add_class::<PyQq3_4_3>()?;
    m.add_class::<PyLaurentRational6>()?;
    m.add_class::<PyLaurentFp3_6>()?;
    m.add_class::<PyLaurentFp5_6>()?;
    m.add_class::<PyLaurentFp7_6>()?;
    m.add_class::<PyLaurentFp11_6>()?;
    m.add_class::<PyLaurentFp13_6>()?;
    m.add_class::<PyLaurentF9_6>()?;
    m.add_class::<PyLaurentF25_6>()?;
    m.add_class::<PyLaurentF27_6>()?;
    m.add_class::<PyRamifiedQp2_4E2>()?;
    m.add_class::<PyRamifiedQp3_4E2>()?;
    m.add_class::<PyRamifiedQp5_4E2>()?;
    m.add_class::<PyRamifiedQp7_4E2>()?;
    m.add_class::<PyRamifiedQp11_4E2>()?;
    m.add_class::<PyRamifiedQp13_4E2>()?;
    m.add_class::<PyRamifiedQp2_4E3>()?;
    m.add_class::<PyRamifiedQp3_4E3>()?;
    m.add_class::<PyRamifiedQp5_4E3>()?;
    m.add_class::<PyRamifiedQp7_4E3>()?;
    m.add_class::<PyRamifiedQp11_4E3>()?;
    m.add_class::<PyRamifiedQp13_4E3>()?;
    m.add_class::<PyGaussQp2_4>()?;
    m.add_class::<PyGaussQp3_4>()?;
    m.add_class::<PyGaussQp5_4>()?;
    m.add_class::<PyGaussQp7_4>()?;
    m.add_class::<PyGaussQp11_4>()?;
    m.add_class::<PyGaussQp13_4>()?;
    m.add_class::<PyRational>()?;
    m.add_class::<PySurreal>()?;
    m.add_class::<PySurcomplex>()?;
    m.add_class::<PyInteger>()?;
    m.add_class::<PyOmnific>()?;
    m.add_class::<PyLocalQp>()?;
    m.add_class::<PyAdelePlace>()?;
    m.add_class::<PyAdele>()?;
    m.add_class::<PyMaxPlusTropical>()?;
    m.add_class::<PyMinPlusTropical>()?;
    m.add_class::<PyOrdinal>()?;
    m.add_class::<PySignExpansion>()?;
    m.add_function(wrap_pyfunction!(omnific, m)?)?;
    m.add_function(wrap_pyfunction!(omnific_omega, m)?)?;
    m.add_function(wrap_pyfunction!(is_omnific_integer, m)?)?;
    m.add_function(wrap_pyfunction!(omega, m)?)?;
    m.add_function(wrap_pyfunction!(epsilon, m)?)?;
    m.add_function(wrap_pyfunction!(omega_pow, m)?)?;
    m.add_function(wrap_pyfunction!(surreal, m)?)?;
    m.add_function(wrap_pyfunction!(nim_add, m)?)?;
    m.add_function(wrap_pyfunction!(nim_mul, m)?)?;
    m.add_function(wrap_pyfunction!(nim_pow, m)?)?;
    m.add_function(wrap_pyfunction!(nim_square, m)?)?;
    m.add_function(wrap_pyfunction!(nim_frobenius_iter, m)?)?;
    m.add_function(wrap_pyfunction!(nim_inv, m)?)?;
    m.add_function(wrap_pyfunction!(nim_sqrt, m)?)?;
    m.add_function(wrap_pyfunction!(nim_trace, m)?)?;
    m.add_function(wrap_pyfunction!(nim_solve_artin_schreier, m)?)?;
    m.add_function(wrap_pyfunction!(nim_is_artin_schreier_solvable, m)?)?;
    m.add_function(wrap_pyfunction!(nim_degree, m)?)?;
    m.add_function(wrap_pyfunction!(nim_conjugates, m)?)?;
    m.add_function(wrap_pyfunction!(nim_min_poly, m)?)?;
    m.add_function(wrap_pyfunction!(nim_relative_trace, m)?)?;
    m.add_function(wrap_pyfunction!(nim_relative_norm, m)?)?;
    m.add_function(wrap_pyfunction!(nim_multiplicative_order, m)?)?;
    m.add_function(wrap_pyfunction!(nim_is_primitive, m)?)?;
    m.add_function(wrap_pyfunction!(nim_primitive_element, m)?)?;
    m.add_function(wrap_pyfunction!(nim_discrete_log, m)?)?;
    m.add_function(wrap_pyfunction!(adele_prec, m)?)?;
    m.add_function(wrap_pyfunction!(min_coeff_valuation, m)?)?;
    Ok(())
}
