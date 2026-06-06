//! PyO3 bindings — per-backend classes.
//!
//! Each scalar world (nimber / surreal / surcomplex) gets its own scalar type,
//! `<World>Algebra`, and `<World>MV` multivector. The Algebra/MV pair is
//! stamped out by the `backend!` macro, monomorphising the same verified
//! generic engine to the concrete scalar type — so there is no runtime
//! dispatch and no way to mix scalar worlds in one algebra.

use crate::clifford::{CliffordAlgebra, Metric, Multivector};
use crate::nimber::Nimber;
use crate::scalar::{Rational, Scalar};
use crate::surcomplex::Surcomplex;
use crate::surreal::Surreal;
use pyo3::basic::CompareOp;
use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyDict;
use pyo3::IntoPyObjectExt;
use std::collections::BTreeMap;
use std::sync::Arc;

// ---------------------------------------------------------------------------
// Scalar pyclasses + parsers
// ---------------------------------------------------------------------------

#[pyclass(name = "Nimber", module = "pleroma", from_py_object)]
#[derive(Clone)]
struct PyNimber {
    inner: Nimber,
}

#[pymethods]
impl PyNimber {
    #[new]
    fn new(value: u64) -> Self {
        PyNimber { inner: Nimber(value) }
    }
    #[getter]
    fn value(&self) -> u64 {
        self.inner.0
    }
    fn __add__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyNimber> {
        Ok(PyNimber { inner: self.inner.add(&parse_nimber(other)?) })
    }
    fn __radd__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyNimber> {
        self.__add__(other)
    }
    fn __mul__(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        // defer to the other operand (e.g. a multivector's __rmul__) if it isn't a scalar
        match parse_nimber(other) {
            Ok(o) => PyNimber { inner: self.inner.mul(&o) }.into_py_any(py),
            Err(_) => Ok(py.NotImplemented()),
        }
    }
    fn __rmul__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyNimber> {
        Ok(PyNimber { inner: self.inner.mul(&parse_nimber(other)?) })
    }
    fn inv(&self) -> PyResult<PyNimber> {
        self.inner
            .inv()
            .map(|n| PyNimber { inner: n })
            .ok_or_else(|| PyValueError::new_err("*0 has no inverse"))
    }
    fn __truediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyNimber> {
        let o = parse_nimber(other)?;
        let oi = o.inv().ok_or_else(|| PyValueError::new_err("division by *0"))?;
        Ok(PyNimber { inner: self.inner.mul(&oi) })
    }
    fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
        matches!(parse_nimber(other), Ok(n) if n == self.inner)
    }
    fn __hash__(&self) -> u64 {
        self.inner.0
    }
    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}

fn parse_nimber(obj: &Bound<'_, PyAny>) -> PyResult<Nimber> {
    if let Ok(n) = obj.cast::<PyNimber>() {
        return Ok(n.borrow().inner);
    }
    if let Ok(v) = obj.extract::<u64>() {
        return Ok(Nimber(v));
    }
    Err(PyTypeError::new_err("expected Nimber or non-negative int"))
}

#[pyclass(name = "Surreal", module = "pleroma", from_py_object)]
#[derive(Clone)]
struct PySurreal {
    inner: Surreal,
}

#[pymethods]
impl PySurreal {
    fn __add__(&self, other: &Bound<'_, PyAny>) -> PyResult<PySurreal> {
        Ok(PySurreal { inner: self.inner.add(&parse_surreal(other)?) })
    }
    fn __radd__(&self, other: &Bound<'_, PyAny>) -> PyResult<PySurreal> {
        self.__add__(other)
    }
    fn __sub__(&self, other: &Bound<'_, PyAny>) -> PyResult<PySurreal> {
        Ok(PySurreal { inner: self.inner.sub(&parse_surreal(other)?) })
    }
    fn __rsub__(&self, other: &Bound<'_, PyAny>) -> PyResult<PySurreal> {
        Ok(PySurreal { inner: parse_surreal(other)?.sub(&self.inner) })
    }
    fn __mul__(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        match parse_surreal(other) {
            Ok(o) => PySurreal { inner: self.inner.mul(&o) }.into_py_any(py),
            Err(_) => Ok(py.NotImplemented()),
        }
    }
    fn __rmul__(&self, other: &Bound<'_, PyAny>) -> PyResult<PySurreal> {
        Ok(PySurreal { inner: self.inner.mul(&parse_surreal(other)?) })
    }
    fn __neg__(&self) -> PySurreal {
        PySurreal { inner: self.inner.neg() }
    }
    fn inv(&self) -> PyResult<PySurreal> {
        self.inner.inv().map(|s| PySurreal { inner: s }).ok_or_else(|| {
            PyValueError::new_err("only monomials (coeff·ω^e) have a finite-support surreal inverse")
        })
    }
    fn __truediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<PySurreal> {
        let o = parse_surreal(other)?;
        let oi = o
            .inv()
            .ok_or_else(|| PyValueError::new_err("divisor has no finite-support inverse"))?;
        Ok(PySurreal { inner: self.inner.mul(&oi) })
    }
    fn __pow__(&self, n: u32, _modulo: Option<&Bound<'_, PyAny>>) -> PySurreal {
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

fn parse_surreal(obj: &Bound<'_, PyAny>) -> PyResult<Surreal> {
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
struct PySurcomplex {
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
        Ok(PySurcomplex { inner: Surcomplex::new(r, i) })
    }
    #[staticmethod]
    fn i() -> PySurcomplex {
        PySurcomplex { inner: Surcomplex::i() }
    }
    #[getter]
    fn re(&self) -> PySurreal {
        PySurreal { inner: self.inner.re.clone() }
    }
    #[getter]
    fn im(&self) -> PySurreal {
        PySurreal { inner: self.inner.im.clone() }
    }
    fn conj(&self) -> PySurcomplex {
        PySurcomplex { inner: self.inner.conj() }
    }
    fn __add__(&self, other: &Bound<'_, PyAny>) -> PyResult<PySurcomplex> {
        Ok(PySurcomplex { inner: self.inner.add(&parse_surcomplex(other)?) })
    }
    fn __radd__(&self, other: &Bound<'_, PyAny>) -> PyResult<PySurcomplex> {
        self.__add__(other)
    }
    fn __sub__(&self, other: &Bound<'_, PyAny>) -> PyResult<PySurcomplex> {
        Ok(PySurcomplex { inner: self.inner.sub(&parse_surcomplex(other)?) })
    }
    fn __mul__(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        match parse_surcomplex(other) {
            Ok(o) => PySurcomplex { inner: self.inner.mul(&o) }.into_py_any(py),
            Err(_) => Ok(py.NotImplemented()),
        }
    }
    fn __rmul__(&self, other: &Bound<'_, PyAny>) -> PyResult<PySurcomplex> {
        Ok(PySurcomplex { inner: self.inner.mul(&parse_surcomplex(other)?) })
    }
    fn __neg__(&self) -> PySurcomplex {
        PySurcomplex { inner: self.inner.neg() }
    }
    fn inv(&self) -> PyResult<PySurcomplex> {
        self.inner.inv().map(|s| PySurcomplex { inner: s }).ok_or_else(|| {
            PyValueError::new_err("inverse needs an invertible norm a²+b² (a monomial surreal)")
        })
    }
    fn __truediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<PySurcomplex> {
        let o = parse_surcomplex(other)?;
        let oi = o
            .inv()
            .ok_or_else(|| PyValueError::new_err("divisor has no representable inverse"))?;
        Ok(PySurcomplex { inner: self.inner.mul(&oi) })
    }
    fn __pow__(&self, n: u32, _modulo: Option<&Bound<'_, PyAny>>) -> PySurcomplex {
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

fn parse_surcomplex(obj: &Bound<'_, PyAny>) -> PyResult<Surcomplex<Surreal>> {
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

fn wrap_nimber(n: Nimber) -> PyNimber {
    PyNimber { inner: n }
}
fn wrap_surreal(s: Surreal) -> PySurreal {
    PySurreal { inner: s }
}
fn wrap_surcomplex(s: Surcomplex<Surreal>) -> PySurcomplex {
    PySurcomplex { inner: s }
}

// ---------------------------------------------------------------------------
// Algebra + multivector, one pair per backend
// ---------------------------------------------------------------------------

macro_rules! backend {
    ($alg:ident, $alg_name:literal, $mv:ident, $mv_name:literal, $scalar:ty, $parse:path, $scalar_py:ty, $wrap:path) => {
        #[pyclass(name = $alg_name, module = "pleroma", from_py_object)]
        #[derive(Clone)]
        struct $alg {
            inner: Arc<CliffordAlgebra<$scalar>>,
        }

        #[pymethods]
        impl $alg {
            #[new]
            #[pyo3(signature = (q, b=None))]
            fn new(q: Vec<Bound<'_, PyAny>>, b: Option<Bound<'_, PyDict>>) -> PyResult<Self> {
                let mut qv: Vec<$scalar> = Vec::with_capacity(q.len());
                for item in &q {
                    qv.push($parse(item)?);
                }
                let mut bm: BTreeMap<(usize, usize), $scalar> = BTreeMap::new();
                if let Some(d) = b {
                    for (k, v) in d.iter() {
                        let (i, j): (usize, usize) = k.extract()?;
                        let key = if i < j { (i, j) } else { (j, i) };
                        bm.insert(key, $parse(&v)?);
                    }
                }
                let dim = qv.len();
                let metric = Metric { q: qv, b: bm };
                Ok($alg { inner: Arc::new(CliffordAlgebra::new(dim, metric)) })
            }

            #[getter]
            fn dim(&self) -> usize {
                self.inner.dim
            }
            fn gen(&self, i: usize) -> $mv {
                $mv { alg: self.inner.clone(), mv: self.inner.gen(i) }
            }
            fn blade(&self, gens: Vec<usize>) -> $mv {
                $mv { alg: self.inner.clone(), mv: self.inner.blade(&gens) }
            }
            fn scalar(&self, s: &Bound<'_, PyAny>) -> PyResult<$mv> {
                Ok($mv { alg: self.inner.clone(), mv: self.inner.scalar($parse(s)?) })
            }
            fn zero(&self) -> $mv {
                $mv { alg: self.inner.clone(), mv: self.inner.zero() }
            }
            fn pseudoscalar(&self) -> $mv {
                $mv { alg: self.inner.clone(), mv: self.inner.pseudoscalar() }
            }
            fn __repr__(&self) -> String {
                format!("{}(dim={})", $alg_name, self.inner.dim)
            }
        }

        #[pyclass(name = $mv_name, module = "pleroma", from_py_object)]
        #[derive(Clone)]
        struct $mv {
            alg: Arc<CliffordAlgebra<$scalar>>,
            mv: Multivector<$scalar>,
        }

        #[pymethods]
        impl $mv {
            fn __add__(&self, other: &$mv) -> $mv {
                $mv { alg: self.alg.clone(), mv: self.alg.add(&self.mv, &other.mv) }
            }
            fn __sub__(&self, other: &$mv) -> $mv {
                let neg_one = <$scalar as Scalar>::one().neg();
                let neg = self.alg.scalar_mul(&neg_one, &other.mv);
                $mv { alg: self.alg.clone(), mv: self.alg.add(&self.mv, &neg) }
            }
            fn __neg__(&self) -> $mv {
                let neg_one = <$scalar as Scalar>::one().neg();
                $mv { alg: self.alg.clone(), mv: self.alg.scalar_mul(&neg_one, &self.mv) }
            }
            fn __mul__(&self, other: &Bound<'_, PyAny>) -> PyResult<$mv> {
                if let Ok(o) = other.cast::<$mv>() {
                    return Ok($mv {
                        alg: self.alg.clone(),
                        mv: self.alg.mul(&self.mv, &o.borrow().mv),
                    });
                }
                let s = $parse(other)?;
                Ok($mv { alg: self.alg.clone(), mv: self.alg.scalar_mul(&s, &self.mv) })
            }
            fn __rmul__(&self, other: &Bound<'_, PyAny>) -> PyResult<$mv> {
                let s = $parse(other)?;
                Ok($mv { alg: self.alg.clone(), mv: self.alg.scalar_mul(&s, &self.mv) })
            }
            fn __pow__(&self, n: u32, _modulo: Option<&Bound<'_, PyAny>>) -> $mv {
                let mut acc = self.alg.scalar(<$scalar as Scalar>::one());
                for _ in 0..n {
                    acc = self.alg.mul(&acc, &self.mv);
                }
                $mv { alg: self.alg.clone(), mv: acc }
            }
            /// Exterior (wedge) product; also bound to the `^` operator.
            fn wedge(&self, other: &$mv) -> $mv {
                $mv { alg: self.alg.clone(), mv: self.alg.wedge(&self.mv, &other.mv) }
            }
            fn __xor__(&self, other: &$mv) -> $mv {
                self.wedge(other)
            }
            fn reverse(&self) -> $mv {
                $mv { alg: self.alg.clone(), mv: self.alg.reverse(&self.mv) }
            }
            /// `~v` is reversion.
            fn __invert__(&self) -> $mv {
                self.reverse()
            }
            fn grade(&self, k: u32) -> $mv {
                $mv { alg: self.alg.clone(), mv: self.alg.grade_part(&self.mv, k) }
            }
            fn grade_involution(&self) -> $mv {
                $mv { alg: self.alg.clone(), mv: self.alg.grade_involution(&self.mv) }
            }
            /// Versor inverse v⁻¹ = ṽ/(v ṽ); errors if v isn't an invertible versor.
            fn inverse(&self) -> PyResult<$mv> {
                self.alg
                    .versor_inverse(&self.mv)
                    .map(|mv| $mv { alg: self.alg.clone(), mv })
                    .ok_or_else(|| PyValueError::new_err("not an invertible versor"))
            }
            /// Sandwich self · x · self⁻¹ (rotor/versor action).
            fn sandwich(&self, x: &$mv) -> PyResult<$mv> {
                self.alg
                    .sandwich(&self.mv, &x.mv)
                    .map(|mv| $mv { alg: self.alg.clone(), mv })
                    .ok_or_else(|| PyValueError::new_err("not an invertible versor"))
            }
            /// Reflect x in the hyperplane ⊥ self (self must be an invertible vector).
            fn reflect(&self, x: &$mv) -> PyResult<$mv> {
                self.alg
                    .reflect(&self.mv, &x.mv)
                    .map(|mv| $mv { alg: self.alg.clone(), mv })
                    .ok_or_else(|| PyValueError::new_err("not an invertible vector"))
            }
            fn left_contract(&self, other: &$mv) -> $mv {
                $mv { alg: self.alg.clone(), mv: self.alg.left_contract(&self.mv, &other.mv) }
            }
            fn right_contract(&self, other: &$mv) -> $mv {
                $mv { alg: self.alg.clone(), mv: self.alg.right_contract(&self.mv, &other.mv) }
            }
            /// `<<` is left contraction, `>>` is right contraction.
            fn __lshift__(&self, other: &$mv) -> $mv {
                self.left_contract(other)
            }
            fn __rshift__(&self, other: &$mv) -> $mv {
                self.right_contract(other)
            }
            fn dual(&self) -> PyResult<$mv> {
                self.alg
                    .dual(&self.mv)
                    .map(|mv| $mv { alg: self.alg.clone(), mv })
                    .ok_or_else(|| PyValueError::new_err("pseudoscalar not invertible (degenerate metric)"))
            }
            fn norm2(&self) -> $scalar_py {
                $wrap(self.alg.norm2(&self.mv))
            }
            fn scalar_part(&self) -> $scalar_py {
                $wrap(self.alg.scalar_part(&self.mv))
            }
            /// Division: by a scalar, or by a versor (multiply by its inverse).
            fn __truediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<$mv> {
                if let Ok(o) = other.cast::<$mv>() {
                    let oinv = self
                        .alg
                        .versor_inverse(&o.borrow().mv)
                        .ok_or_else(|| PyValueError::new_err("divisor not an invertible versor"))?;
                    return Ok($mv { alg: self.alg.clone(), mv: self.alg.mul(&self.mv, &oinv) });
                }
                let s = $parse(other)?;
                let sinv = <$scalar as Scalar>::inv(&s)
                    .ok_or_else(|| PyValueError::new_err("scalar has no representable inverse"))?;
                Ok($mv { alg: self.alg.clone(), mv: self.alg.scalar_mul(&sinv, &self.mv) })
            }
            fn is_zero(&self) -> bool {
                self.mv.is_zero()
            }
            fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
                if let Ok(o) = other.cast::<$mv>() {
                    self.mv == o.borrow().mv
                } else {
                    false
                }
            }
            fn __repr__(&self) -> String {
                self.mv.display()
            }
        }
    };
}

backend!(
    NimberAlgebra, "NimberAlgebra", NimberMV, "NimberMV",
    Nimber, parse_nimber, PyNimber, wrap_nimber
);
backend!(
    SurrealAlgebra, "SurrealAlgebra", SurrealMV, "SurrealMV",
    Surreal, parse_surreal, PySurreal, wrap_surreal
);
backend!(
    SurcomplexAlgebra, "SurcomplexAlgebra", SurcomplexMV, "SurcomplexMV",
    Surcomplex<Surreal>, parse_surcomplex, PySurcomplex, wrap_surcomplex
);

// ---------------------------------------------------------------------------
// Surreal builders
// ---------------------------------------------------------------------------

#[pyfunction]
fn omega() -> PySurreal {
    PySurreal { inner: Surreal::omega() }
}

#[pyfunction]
fn epsilon() -> PySurreal {
    PySurreal { inner: Surreal::epsilon() }
}

#[pyfunction]
fn omega_pow(exp: &Bound<'_, PyAny>) -> PyResult<PySurreal> {
    Ok(PySurreal { inner: Surreal::omega_pow(parse_surreal(exp)?) })
}

#[pyfunction]
fn rational(num: i128, den: i128) -> PyResult<PySurreal> {
    if den == 0 {
        return Err(PyValueError::new_err("zero denominator"));
    }
    Ok(PySurreal { inner: Surreal::from_rational(Rational::new(num, den)) })
}

#[pyfunction]
fn surreal(n: i128) -> PySurreal {
    PySurreal { inner: Surreal::from_int(n) }
}

#[pyclass(name = "ArfResult", module = "pleroma")]
struct PyArfResult {
    inner: crate::arf::ArfResult,
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
fn arf_invariant(alg: &NimberAlgebra) -> PyArfResult {
    PyArfResult { inner: crate::arf::arf_invariant(&alg.inner.metric) }
}

#[pyfunction]
fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

#[pymodule]
fn pleroma(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyNimber>()?;
    m.add_class::<NimberAlgebra>()?;
    m.add_class::<NimberMV>()?;
    m.add_class::<PySurreal>()?;
    m.add_class::<SurrealAlgebra>()?;
    m.add_class::<SurrealMV>()?;
    m.add_class::<PySurcomplex>()?;
    m.add_class::<SurcomplexAlgebra>()?;
    m.add_class::<SurcomplexMV>()?;
    m.add_function(wrap_pyfunction!(omega, m)?)?;
    m.add_function(wrap_pyfunction!(epsilon, m)?)?;
    m.add_function(wrap_pyfunction!(omega_pow, m)?)?;
    m.add_function(wrap_pyfunction!(rational, m)?)?;
    m.add_function(wrap_pyfunction!(surreal, m)?)?;
    m.add_class::<PyArfResult>()?;
    m.add_function(wrap_pyfunction!(arf_invariant, m)?)?;
    m.add_function(wrap_pyfunction!(version, m)?)?;
    Ok(())
}
