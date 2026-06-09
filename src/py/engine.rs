//! Python bindings for the GA engine: the `backend!` macro that stamps out one
//! `<World>Algebra` / `<World>MV` pyclass pair per scalar backend, the runtime
//! invocations, and conformal GA (`Cga`). The generated structs and their
//! fields are `pub(crate)` so the classifier bindings in [`super::forms`] and the
//! game-exterior binding in [`super::games`] can read `.inner` / `.mv`.

use super::scalars::{
    parse_adele, parse_f16, parse_f25, parse_f27, parse_f4, parse_f8, parse_f9, parse_fp11,
    parse_fp11_poly, parse_fp11_rational_function, parse_fp13, parse_fp13_poly,
    parse_fp13_rational_function, parse_fp2, parse_fp2_poly, parse_fp2_rational_function,
    parse_fp3, parse_fp3_poly, parse_fp3_rational_function, parse_fp5, parse_fp5_poly,
    parse_fp5_rational_function, parse_fp7, parse_fp7_poly, parse_fp7_rational_function,
    parse_gauss_qp11_4, parse_gauss_qp13_4, parse_gauss_qp2_4, parse_gauss_qp3_4,
    parse_gauss_qp5_4, parse_gauss_qp7_4, parse_integer, parse_laurent_f25_6, parse_laurent_f27_6,
    parse_laurent_f9_6, parse_laurent_fp11_6, parse_laurent_fp13_6, parse_laurent_fp3_6,
    parse_laurent_fp5_6, parse_laurent_fp7_6, parse_laurent_rational_6, parse_nimber,
    parse_nimber_poly, parse_nimber_rational_function, parse_omnific, parse_ordinal, parse_qp11_4,
    parse_qp13_4, parse_qp2_4, parse_qp3_4, parse_qp5_4, parse_qp7_4, parse_qq2_4_2, parse_qq2_4_3,
    parse_qq2_4_4, parse_qq3_4_2, parse_qq3_4_3, parse_qq5_4_2, parse_ramified_qp11_4_e2,
    parse_ramified_qp11_4_e3, parse_ramified_qp13_4_e2, parse_ramified_qp13_4_e3,
    parse_ramified_qp2_4_e2, parse_ramified_qp2_4_e3, parse_ramified_qp3_4_e2,
    parse_ramified_qp3_4_e3, parse_ramified_qp5_4_e2, parse_ramified_qp5_4_e3,
    parse_ramified_qp7_4_e2, parse_ramified_qp7_4_e3, parse_rational, parse_surcomplex,
    parse_surreal, parse_witt_vec2_4_2, parse_witt_vec2_4_3, parse_witt_vec2_4_4,
    parse_witt_vec3_4_2, parse_witt_vec3_4_3, parse_witt_vec5_4_2, parse_zp11_4, parse_zp13_4,
    parse_zp2_4, parse_zp3_4, parse_zp5_4, parse_zp7_4, wrap_adele, wrap_f16, wrap_f25, wrap_f27,
    wrap_f4, wrap_f8, wrap_f9, wrap_fp11, wrap_fp11_poly, wrap_fp11_rational_function, wrap_fp13,
    wrap_fp13_poly, wrap_fp13_rational_function, wrap_fp2, wrap_fp2_poly,
    wrap_fp2_rational_function, wrap_fp3, wrap_fp3_poly, wrap_fp3_rational_function, wrap_fp5,
    wrap_fp5_poly, wrap_fp5_rational_function, wrap_fp7, wrap_fp7_poly, wrap_fp7_rational_function,
    wrap_gauss_qp11_4, wrap_gauss_qp13_4, wrap_gauss_qp2_4, wrap_gauss_qp3_4, wrap_gauss_qp5_4,
    wrap_gauss_qp7_4, wrap_integer, wrap_laurent_f25_6, wrap_laurent_f27_6, wrap_laurent_f9_6,
    wrap_laurent_fp11_6, wrap_laurent_fp13_6, wrap_laurent_fp3_6, wrap_laurent_fp5_6,
    wrap_laurent_fp7_6, wrap_laurent_rational_6, wrap_nimber, wrap_nimber_poly,
    wrap_nimber_rational_function, wrap_omnific, wrap_ordinal, wrap_qp11_4, wrap_qp13_4,
    wrap_qp2_4, wrap_qp3_4, wrap_qp5_4, wrap_qp7_4, wrap_qq2_4_2, wrap_qq2_4_3, wrap_qq2_4_4,
    wrap_qq3_4_2, wrap_qq3_4_3, wrap_qq5_4_2, wrap_ramified_qp11_4_e2, wrap_ramified_qp11_4_e3,
    wrap_ramified_qp13_4_e2, wrap_ramified_qp13_4_e3, wrap_ramified_qp2_4_e2,
    wrap_ramified_qp2_4_e3, wrap_ramified_qp3_4_e2, wrap_ramified_qp3_4_e3, wrap_ramified_qp5_4_e2,
    wrap_ramified_qp5_4_e3, wrap_ramified_qp7_4_e2, wrap_ramified_qp7_4_e3, wrap_rational,
    wrap_surcomplex, wrap_surreal, wrap_witt_vec2_4_2, wrap_witt_vec2_4_3, wrap_witt_vec2_4_4,
    wrap_witt_vec3_4_2, wrap_witt_vec3_4_3, wrap_witt_vec5_4_2, wrap_zp11_4, wrap_zp13_4,
    wrap_zp2_4, wrap_zp3_4, wrap_zp5_4, wrap_zp7_4, PyAdele, PyF16, PyF25, PyF27, PyF4, PyF8, PyF9,
    PyFp11, PyFp11Poly, PyFp11RationalFunction, PyFp13, PyFp13Poly, PyFp13RationalFunction, PyFp2,
    PyFp2Poly, PyFp2RationalFunction, PyFp3, PyFp3Poly, PyFp3RationalFunction, PyFp5, PyFp5Poly,
    PyFp5RationalFunction, PyFp7, PyFp7Poly, PyFp7RationalFunction, PyGaussQp11_4, PyGaussQp13_4,
    PyGaussQp2_4, PyGaussQp3_4, PyGaussQp5_4, PyGaussQp7_4, PyInteger, PyLaurentF25_6,
    PyLaurentF27_6, PyLaurentF9_6, PyLaurentFp11_6, PyLaurentFp13_6, PyLaurentFp3_6,
    PyLaurentFp5_6, PyLaurentFp7_6, PyLaurentRational6, PyNimber, PyNimberPoly,
    PyNimberRationalFunction, PyOmnific, PyOrdinal, PyQp11_4, PyQp13_4, PyQp2_4, PyQp3_4, PyQp5_4,
    PyQp7_4, PyQq2_4_2, PyQq2_4_3, PyQq2_4_4, PyQq3_4_2, PyQq3_4_3, PyQq5_4_2, PyRamifiedQp11_4E2,
    PyRamifiedQp11_4E3, PyRamifiedQp13_4E2, PyRamifiedQp13_4E3, PyRamifiedQp2_4E2,
    PyRamifiedQp2_4E3, PyRamifiedQp3_4E2, PyRamifiedQp3_4E3, PyRamifiedQp5_4E2, PyRamifiedQp5_4E3,
    PyRamifiedQp7_4E2, PyRamifiedQp7_4E3, PyRational, PySurcomplex, PySurreal, PyWittVec2_4_2,
    PyWittVec2_4_3, PyWittVec2_4_4, PyWittVec3_4_2, PyWittVec3_4_3, PyWittVec5_4_2, PyZp11_4,
    PyZp13_4, PyZp2_4, PyZp3_4, PyZp5_4, PyZp7_4,
};
use crate::clifford::{
    Cga, CliffordAlgebra, DividedPowerAlgebra, DpVector, LinearMap, Metric, Multivector,
    MAX_BASIS_DIM,
};
use crate::scalar::{
    Adele, Fp, Fpn, Gauss, Integer, Laurent, Nimber, Omnific, Ordinal, Poly, Qp, Qq, Ramified,
    Rational, RationalFunction, Scalar, Surcomplex, Surreal, WittVec, Zp,
};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use pyo3::IntoPyObjectExt;
use std::collections::BTreeMap;
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::sync::Arc;
use std::sync::Mutex;

static PANIC_HOOK_LOCK: Mutex<()> = Mutex::new(());

fn panic_payload_message(payload: Box<dyn std::any::Any + Send>) -> String {
    if let Some(s) = payload.downcast_ref::<&str>() {
        (*s).to_string()
    } else if let Some(s) = payload.downcast_ref::<String>() {
        s.clone()
    } else {
        "Rust operation panicked".to_string()
    }
}

fn scalar_boundary<T>(f: impl FnOnce() -> T) -> PyResult<T> {
    let _guard = PANIC_HOOK_LOCK
        .lock()
        .map_err(|_| PyValueError::new_err("panic hook lock poisoned"))?;
    let old_hook = std::panic::take_hook();
    std::panic::set_hook(Box::new(|_| {}));
    let result = catch_unwind(AssertUnwindSafe(f));
    std::panic::set_hook(old_hook);
    result.map_err(|payload| {
        PyValueError::new_err(format!(
            "operation escaped the represented scalar boundary: {}",
            panic_payload_message(payload)
        ))
    })
}

#[pyclass(name = "SpinorRep", module = "pleroma")]
struct PySpinorRep {
    idempotent: Py<PyAny>,
    basis: Py<PyAny>,
    gen_matrices: Py<PyAny>,
    is_left_regular: bool,
    diagonalized_metric: Py<PyAny>,
    orthogonal_basis_in_original: Py<PyAny>,
    basis_dim: usize,
    generator_count: usize,
}

#[pymethods]
impl PySpinorRep {
    #[getter]
    fn idempotent(&self, py: Python<'_>) -> Py<PyAny> {
        self.idempotent.clone_ref(py)
    }
    #[getter]
    fn basis(&self, py: Python<'_>) -> Py<PyAny> {
        self.basis.clone_ref(py)
    }
    #[getter]
    fn gen_matrices(&self, py: Python<'_>) -> Py<PyAny> {
        self.gen_matrices.clone_ref(py)
    }
    #[getter]
    fn is_left_regular(&self) -> bool {
        self.is_left_regular
    }
    #[getter]
    fn diagonalized_metric(&self, py: Python<'_>) -> Py<PyAny> {
        self.diagonalized_metric.clone_ref(py)
    }
    #[getter]
    fn orthogonal_basis_in_original(&self, py: Python<'_>) -> Py<PyAny> {
        self.orthogonal_basis_in_original.clone_ref(py)
    }
    fn __repr__(&self) -> String {
        format!(
            "SpinorRep(basis_dim={}, generators={}, is_left_regular={})",
            self.basis_dim, self.generator_count, self.is_left_regular
        )
    }
}

#[pyclass(name = "LazySpinorRep", module = "pleroma")]
struct PyLazySpinorRep {
    algebra: Py<PyAny>,
}

#[pymethods]
impl PyLazySpinorRep {
    #[getter]
    fn algebra(&self, py: Python<'_>) -> Py<PyAny> {
        self.algebra.clone_ref(py)
    }

    /// Apply left multiplication by generator `e_i` to a sparse module vector.
    fn apply_generator(
        &self,
        py: Python<'_>,
        i: usize,
        v: Bound<'_, PyAny>,
    ) -> PyResult<Py<PyAny>> {
        Ok(self
            .algebra
            .bind(py)
            .call_method1("apply_generator", (i, v))?
            .unbind())
    }

    /// Apply left multiplication by the vector `Σ coeffs[i] e_i`.
    fn apply_vector(
        &self,
        py: Python<'_>,
        coeffs: Vec<Bound<'_, PyAny>>,
        v: Bound<'_, PyAny>,
    ) -> PyResult<Py<PyAny>> {
        Ok(self
            .algebra
            .bind(py)
            .call_method1("apply_vector", (coeffs, v))?
            .unbind())
    }

    fn __repr__(&self) -> String {
        "LazySpinorRep()".to_string()
    }
}

#[pyclass(name = "VersorClass", module = "pleroma")]
struct PyVersorClass {
    spinor_norm: Py<PyAny>,
    dickson: u128,
}

#[pymethods]
impl PyVersorClass {
    #[getter]
    fn spinor_norm(&self, py: Python<'_>) -> Py<PyAny> {
        self.spinor_norm.clone_ref(py)
    }
    #[getter]
    fn dickson(&self) -> u128 {
        self.dickson
    }
    fn __repr__(&self) -> String {
        format!("VersorClass(dickson={})", self.dickson)
    }
}

fn prime_field_identity_linear_map(py: Python<'_>, p: u128) -> PyResult<Py<PyAny>> {
    match p {
        2 => Fp2LinearMap {
            inner: LinearMap::<Fp<2>>::identity(1),
        }
        .into_py_any(py),
        3 => Fp3LinearMap {
            inner: LinearMap::<Fp<3>>::identity(1),
        }
        .into_py_any(py),
        5 => Fp5LinearMap {
            inner: LinearMap::<Fp<5>>::identity(1),
        }
        .into_py_any(py),
        7 => Fp7LinearMap {
            inner: LinearMap::<Fp<7>>::identity(1),
        }
        .into_py_any(py),
        11 => Fp11LinearMap {
            inner: LinearMap::<Fp<11>>::identity(1),
        }
        .into_py_any(py),
        13 => Fp13LinearMap {
            inner: LinearMap::<Fp<13>>::identity(1),
        }
        .into_py_any(py),
        _ => Err(PyValueError::new_err(
            "unsupported prime field; expected p in {2,3,5,7,11,13}",
        )),
    }
}

/// Rust-name fixed-dispatch constructor for the base-field Galois `LinearMap`.
#[pyfunction]
#[pyo3(signature = (p, degree, power=1))]
fn galois_linear_map(py: Python<'_>, p: u128, degree: usize, power: usize) -> PyResult<Py<PyAny>> {
    match (p, degree) {
        (_, 1) => prime_field_identity_linear_map(py, p),
        (2, 2) => Fp2LinearMap {
            inner: crate::clifford::galois_linear_map::<Fpn<2, 2>>(power),
        }
        .into_py_any(py),
        (2, 3) => Fp2LinearMap {
            inner: crate::clifford::galois_linear_map::<Fpn<2, 3>>(power),
        }
        .into_py_any(py),
        (2, 4) => Fp2LinearMap {
            inner: crate::clifford::galois_linear_map::<Fpn<2, 4>>(power),
        }
        .into_py_any(py),
        (3, 2) => Fp3LinearMap {
            inner: crate::clifford::galois_linear_map::<Fpn<3, 2>>(power),
        }
        .into_py_any(py),
        (3, 3) => Fp3LinearMap {
            inner: crate::clifford::galois_linear_map::<Fpn<3, 3>>(power),
        }
        .into_py_any(py),
        (5, 2) => Fp5LinearMap {
            inner: crate::clifford::galois_linear_map::<Fpn<5, 2>>(power),
        }
        .into_py_any(py),
        _ => Err(PyValueError::new_err(
            "unsupported finite field; expected one of F_p, F4, F8, F16, F9, F25, F27",
        )),
    }
}

/// Rust-name fixed-dispatch constructor for the base-field Frobenius `LinearMap`.
#[pyfunction]
fn frobenius_linear_map(py: Python<'_>, p: u128, degree: usize) -> PyResult<Py<PyAny>> {
    galois_linear_map(py, p, degree, 1)
}

/// Rust-name constructor for the represented nimber-subfield Frobenius `LinearMap`.
#[pyfunction]
#[pyo3(signature = (m, power=1))]
fn nimber_subfield_frobenius_linear_map(
    py: Python<'_>,
    m: usize,
    power: usize,
) -> PyResult<Py<PyAny>> {
    if !m.is_power_of_two() || m > 128 {
        return Err(PyValueError::new_err(
            "nimber subfield degree m must be a power of two <= 128",
        ));
    }
    Fp2LinearMap {
        inner: crate::clifford::nimber_subfield_frobenius_linear_map(m, power),
    }
    .into_py_any(py)
}

/// Ascending generator indices in a blade mask.
#[pyfunction]
fn bits(mask: u128) -> Vec<usize> {
    crate::clifford::bits(mask)
}

/// The grade of a blade mask.
#[pyfunction]
fn grade(mask: u128) -> usize {
    crate::clifford::grade(mask)
}

// ---------------------------------------------------------------------------
// Algebra + multivector, one pair per backend
// ---------------------------------------------------------------------------

// Linear-map wrapper for one concrete scalar backend.
macro_rules! backend_linear_map {
    (
        $alg:ident,
        $alg_name:literal,
        $mv:ident,
        $mv_name:literal,
        $lm:ident,
        $lm_name:literal,
        $scalar:ty,
        $parse:path,
        $scalar_py:ty,
        $wrap:path
    ) => {
        #[pyclass(name = $lm_name, module = "pleroma", from_py_object)]
        #[derive(Clone)]
        pub(crate) struct $lm {
            pub(crate) inner: LinearMap<$scalar>,
        }

        #[pymethods]
        impl $lm {
            /// Rust-name constructor for a column-major `LinearMap`.
            #[staticmethod]
            fn from_columns(cols: Vec<Vec<Bound<'_, PyAny>>>) -> PyResult<Self> {
                Ok($lm {
                    inner: $lm::parse_columns(cols)?,
                })
            }

            /// Rust-name constructor for the identity map on `n` generators.
            #[staticmethod]
            fn identity(n: usize) -> Self {
                $lm {
                    inner: LinearMap::<$scalar>::identity(n),
                }
            }

            #[getter]
            fn n(&self) -> usize {
                self.inner.n
            }

            #[getter]
            fn cols(&self) -> Vec<Vec<$scalar_py>> {
                self.columns_py()
            }

            /// Rust-name `LinearMap::image`: return `f(e_i)` as a grade-1
            /// multivector in the given algebra.
            fn image(&self, alg: &$alg, i: usize) -> PyResult<$mv> {
                alg.ensure_linear_map(&self.inner)?;
                if i >= alg.inner.dim {
                    return Err(PyValueError::new_err("linear-map image index out of range"));
                }
                Ok($mv {
                    alg: alg.inner.clone(),
                    mv: scalar_boundary(|| self.inner.image(&alg.inner, i))?,
                })
            }

            /// Rust-name `LinearMap::compose`: `self ∘ inner`.
            fn compose(&self, inner: &$lm) -> PyResult<$lm> {
                if self.inner.n != inner.inner.n {
                    return Err(PyValueError::new_err("dimension mismatch in compose"));
                }
                Ok($lm {
                    inner: scalar_boundary(|| self.inner.compose(&inner.inner))?,
                })
            }

            fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
                if let Ok(o) = other.cast::<$lm>() {
                    self.inner == o.borrow().inner
                } else {
                    false
                }
            }

            fn __repr__(&self) -> String {
                format!("{}(n={})", $lm_name, self.inner.n)
            }
        }

        impl $lm {
            fn parse_columns(cols: Vec<Vec<Bound<'_, PyAny>>>) -> PyResult<LinearMap<$scalar>> {
                let n = cols.len();
                let mut parsed: Vec<Vec<$scalar>> = Vec::with_capacity(n);
                for col in &cols {
                    if col.len() != n {
                        return Err(PyValueError::new_err(
                            "LinearMap must be square: n columns of length n",
                        ));
                    }
                    let mut out = Vec::with_capacity(n);
                    for x in col {
                        out.push($parse(x)?);
                    }
                    parsed.push(out);
                }
                Ok(LinearMap::from_columns(parsed))
            }

            fn columns_py(&self) -> Vec<Vec<$scalar_py>> {
                self.inner
                    .cols
                    .iter()
                    .map(|col| col.iter().cloned().map($wrap).collect())
                    .collect()
            }
        }
    };
}

// Algebra wrapper and algebra-level operations for one backend.
macro_rules! backend_algebra {
    (
        $alg:ident,
        $alg_name:literal,
        $mv:ident,
        $mv_name:literal,
        $lm:ident,
        $lm_name:literal,
        $scalar:ty,
        $parse:path,
        $scalar_py:ty,
        $wrap:path
    ) => {
        #[pyclass(name = $alg_name, module = "pleroma", from_py_object)]
        #[derive(Clone)]
        pub(crate) struct $alg {
            pub(crate) inner: Arc<CliffordAlgebra<$scalar>>,
        }


        #[pymethods]
        impl $alg {
            #[new]
            #[pyo3(signature = (q, b=None, a=None))]
            fn new(
                q: Vec<Bound<'_, PyAny>>,
                b: Option<Bound<'_, PyDict>>,
                a: Option<Bound<'_, PyDict>>,
            ) -> PyResult<Self> {
                let mut qv: Vec<$scalar> = Vec::with_capacity(q.len());
                for item in &q {
                    qv.push($parse(item)?);
                }
                let dim = qv.len();
                if dim > MAX_BASIS_DIM {
                    return Err(PyValueError::new_err(format!(
                        "algebra dimension must be <= {MAX_BASIS_DIM}"
                    )));
                }
                let mut bm: BTreeMap<(usize, usize), $scalar> = BTreeMap::new();
                if let Some(d) = b {
                    for (k, v) in d.iter() {
                        let (i, j): (usize, usize) = k.extract()?;
                        if i == j {
                            return Err(PyValueError::new_err("b-keys must be off-diagonal"));
                        }
                        if i >= dim || j >= dim {
                            return Err(PyValueError::new_err("b-key index out of range"));
                        }
                        let key = if i < j { (i, j) } else { (j, i) };
                        bm.insert(key, $parse(&v)?);
                    }
                }
                // `a` (the in-order / asymmetric contraction) is keyed (i,j) with
                // i<j; it promotes the algebra to a general bilinear form.
                let mut am: BTreeMap<(usize, usize), $scalar> = BTreeMap::new();
                if let Some(d) = a {
                    for (k, v) in d.iter() {
                        let (i, j): (usize, usize) = k.extract()?;
                        if i >= j {
                            return Err(PyValueError::new_err("a-keys must satisfy i < j"));
                        }
                        if j >= dim {
                            return Err(PyValueError::new_err("a-key index out of range"));
                        }
                        am.insert((i, j), $parse(&v)?);
                    }
                }
                let metric = Metric::general(qv, bm, am);
                Ok($alg {
                    inner: Arc::new(CliffordAlgebra::new(dim, metric)),
                })
            }

            #[getter]
            fn dim(&self) -> usize {
                self.inner.dim
            }

            /// Rust-name constructor for a general-bilinear metric algebra.
            #[staticmethod]
            #[pyo3(signature = (q, b=None, a=None))]
            fn general(
                q: Vec<Bound<'_, PyAny>>,
                b: Option<Bound<'_, PyDict>>,
                a: Option<Bound<'_, PyDict>>,
            ) -> PyResult<Self> {
                Self::new(q, b, a)
            }

            /// Rust-name constructor for the fully-null Grassmann/exterior metric.
            #[staticmethod]
            fn grassmann(n: usize) -> PyResult<Self> {
                if n > MAX_BASIS_DIM {
                    return Err(PyValueError::new_err(format!(
                        "algebra dimension must be <= {MAX_BASIS_DIM}"
                    )));
                }
                let metric = Metric::grassmann(n);
                Ok($alg {
                    inner: Arc::new(CliffordAlgebra::new(n, metric)),
                })
            }

            /// Diagonal quadratic entries `q[i] = e_i^2`.
            fn q(&self) -> Vec<$scalar_py> {
                self.inner.metric.q.iter().cloned().map($wrap).collect()
            }

            /// Nonzero polar entries `(i, j, value)` with `i < j`.
            fn b_terms(&self) -> Vec<(usize, usize, $scalar_py)> {
                self.inner
                    .metric
                    .b
                    .iter()
                    .filter(|(_, v)| !v.is_zero())
                    .map(|(&(i, j), v)| (i, j, $wrap(v.clone())))
                    .collect()
            }

            /// Nonzero upper/in-order contraction entries `(i, j, value)` with `i < j`.
            fn a_terms(&self) -> Vec<(usize, usize, $scalar_py)> {
                self.inner
                    .metric
                    .a
                    .iter()
                    .filter(|(_, v)| !v.is_zero())
                    .map(|(&(i, j), v)| (i, j, $wrap(v.clone())))
                    .collect()
            }

            /// Rust-name metric map, restricted to this same Python backend.
            fn map(&self, py: Python<'_>, f: Bound<'_, PyAny>) -> PyResult<$alg> {
                let apply = |coeff: &$scalar| -> PyResult<$scalar> {
                    let py_coeff = $wrap(coeff.clone()).into_py_any(py)?;
                    let mapped = f.call1((py_coeff,))?;
                    $parse(&mapped)
                };
                let q = self
                    .inner
                    .metric
                    .q
                    .iter()
                    .map(&apply)
                    .collect::<PyResult<Vec<_>>>()?;
                let mut b = BTreeMap::new();
                for (&key, coeff) in &self.inner.metric.b {
                    b.insert(key, apply(coeff)?);
                }
                let mut a = BTreeMap::new();
                for (&key, coeff) in &self.inner.metric.a {
                    a.insert(key, apply(coeff)?);
                }
                let metric = Metric::general(q, b, a);
                Ok($alg {
                    inner: Arc::new(CliffordAlgebra::new(self.inner.dim, metric)),
                })
            }

            /// Rust-name helper: `q[i]`, or zero outside the represented diagonal.
            fn q_val(&self, i: usize) -> $scalar_py {
                $wrap(self.inner.metric.q_val(i))
            }

            /// Rust-name helper: whether the metric has any upper/in-order
            /// contraction terms and therefore needs the general product path.
            fn has_upper(&self) -> bool {
                self.inner.metric.has_upper()
            }

            /// Rust-name helper: whether this basis is orthogonal.
            fn is_orthogonal(&self) -> bool {
                self.inner.metric.is_orthogonal()
            }

            /// The graded (super) tensor product self ⊗̂ other ≅ Cl(self ⟂ other).
            fn graded_tensor(&self, other: &$alg) -> PyResult<$alg> {
                if self.inner.dim + other.inner.dim > MAX_BASIS_DIM {
                    return Err(PyValueError::new_err(format!(
                        "graded tensor dimension exceeds {MAX_BASIS_DIM}"
                    )));
                }
                Ok($alg {
                    inner: Arc::new(self.inner.graded_tensor(&other.inner)),
                })
            }

            /// The tensor square `Cl ⊗̂ Cl`, used by the exterior Hopf coproduct.
            fn tensor_square(&self) -> PyResult<$alg> {
                if self.inner.dim * 2 > MAX_BASIS_DIM {
                    return Err(PyValueError::new_err(format!(
                        "tensor square dimension exceeds {MAX_BASIS_DIM}"
                    )));
                }
                Ok($alg {
                    inner: Arc::new(crate::clifford::tensor_square(&self.inner)),
                })
            }

            /// Embed a multivector of the first graded-tensor factor into this
            /// target algebra.
            fn embed_first(&self, mv: &$mv) -> PyResult<$mv> {
                if mv.alg.dim > self.inner.dim {
                    return Err(PyValueError::new_err(
                        "source multivector dimension exceeds target algebra dimension",
                    ));
                }
                Ok($mv {
                    alg: self.inner.clone(),
                    mv: self.inner.embed_first(&mv.mv),
                })
            }

            /// Embed a multivector of the second graded-tensor factor into this
            /// target algebra by shifting its blade masks by `shift`.
            fn embed_second(&self, mv: &$mv, shift: usize) -> PyResult<$mv> {
                if shift + mv.alg.dim > self.inner.dim {
                    return Err(PyValueError::new_err(
                        "shifted source multivector dimension exceeds target algebra dimension",
                    ));
                }
                Ok($mv {
                    alg: self.inner.clone(),
                    mv: scalar_boundary(|| self.inner.embed_second(&mv.mv, shift))?,
                })
            }

            /// Tensor product of diagonal quadratic-form representatives:
            /// `<a_i> tensor <b_j> = <a_i b_j>`. This is the Witt-ring
            /// multiplication on representatives, distinct from the Clifford
            /// graded tensor product.
            fn tensor_form(&self, other: &$alg) -> PyResult<$alg> {
                let metric = scalar_boundary(|| {
                    crate::forms::tensor_form(&self.inner.metric, &other.inner.metric)
                })?
                .ok_or_else(|| {
                    PyValueError::new_err(
                        "tensor_form needs diagonal form representatives (empty b and a)",
                    )
                })?;
                Ok($alg {
                    inner: Arc::new(CliffordAlgebra::new(metric.q.len(), metric)),
                })
            }

            /// Membership in the fundamental ideal I: for a diagonal representative,
            /// the nondegenerate rank is even.
            fn in_fundamental_ideal(&self) -> PyResult<bool> {
                scalar_boundary(|| crate::forms::in_fundamental_ideal(&self.inner.metric))?
                    .ok_or_else(|| {
                        PyValueError::new_err(
                            "in_fundamental_ideal needs a diagonal form representative",
                        )
                    })
            }

            /// The 1-fold Pfister form `<<a>> = <1, -a>` over this scalar backend.
            #[staticmethod]
            fn pfister1(scale: &Bound<'_, PyAny>) -> PyResult<$alg> {
                let scale = $parse(scale)?;
                let metric = scalar_boundary(|| crate::forms::pfister1(&scale))?;
                Ok($alg {
                    inner: Arc::new(CliffordAlgebra::new(metric.q.len(), metric)),
                })
            }

            /// The n-fold Pfister form `<<a_1,...,a_n>>` over this scalar backend.
            /// The empty product is `<1>`.
            #[staticmethod]
            fn pfister(scales: Vec<Bound<'_, PyAny>>) -> PyResult<$alg> {
                let mut parsed = Vec::with_capacity(scales.len());
                for scale in &scales {
                    parsed.push($parse(scale)?);
                }
                let metric = scalar_boundary(|| crate::forms::pfister(&parsed))?;
                Ok($alg {
                    inner: Arc::new(CliffordAlgebra::new(metric.q.len(), metric)),
                })
            }

            /// Projective geometric algebra `Cl(n,0,1)`: one null ideal/projective
            /// direction followed by `n` unit directions.
            #[staticmethod]
            fn pga(n: usize) -> PyResult<$alg> {
                if n >= MAX_BASIS_DIM {
                    return Err(PyValueError::new_err(format!(
                        "PGA total dimension must be <= {MAX_BASIS_DIM}"
                    )));
                }
                Ok($alg {
                    inner: Arc::new(crate::clifford::pga::<$scalar>(n)),
                })
            }

            /// The even subalgebra as a Clifford algebra one dimension smaller
            /// (orthogonal metrics with a non-null generator only).
            fn even_subalgebra(&self) -> PyResult<$alg> {
                self.inner
                    .even_subalgebra()
                    .map(|a| $alg { inner: Arc::new(a) })
                    .ok_or_else(|| {
                        PyValueError::new_err(
                            "even subalgebra needs an orthogonal metric with a non-null generator",
                        )
                    })
            }
            fn gen(&self, i: usize) -> PyResult<$mv> {
                if i >= self.inner.dim {
                    return Err(PyValueError::new_err("generator index out of range"));
                }
                Ok($mv {
                    alg: self.inner.clone(),
                    mv: self.inner.gen(i),
                })
            }
            fn blade(&self, gens: Vec<usize>) -> PyResult<$mv> {
                let mut seen = std::collections::BTreeSet::new();
                for &g in &gens {
                    if g >= self.inner.dim {
                        return Err(PyValueError::new_err("blade generator index out of range"));
                    }
                    if !seen.insert(g) {
                        return Err(PyValueError::new_err("blade expects distinct generators"));
                    }
                }
                Ok($mv {
                    alg: self.inner.clone(),
                    mv: self.inner.blade(&gens),
                })
            }
            fn scalar(&self, s: &Bound<'_, PyAny>) -> PyResult<$mv> {
                Ok($mv {
                    alg: self.inner.clone(),
                    mv: self.inner.scalar($parse(s)?),
                })
            }
            fn zero(&self) -> $mv {
                $mv {
                    alg: self.inner.clone(),
                    mv: self.inner.zero(),
                }
            }
            fn pseudoscalar(&self) -> $mv {
                $mv {
                    alg: self.inner.clone(),
                    mv: self.inner.pseudoscalar(),
                }
            }

            /// The symmetric Gram matrix of the quadratic form, using
            /// `b/2` off-diagonal. Undefined in characteristic 2.
            fn gram(&self) -> PyResult<Vec<Vec<$scalar_py>>> {
                crate::forms::gram(&self.inner.metric)
                    .map(|rows| {
                        rows.into_iter()
                            .map(|row| row.into_iter().map($wrap).collect())
                            .collect()
                    })
                    .ok_or_else(|| {
                        PyValueError::new_err(
                            "Gram matrix needs 2 invertible; in characteristic 2 use the polar form directly",
                        )
                    })
            }

            /// Congruence-diagonalize the symmetric form, if possible. Returns
            /// `None` in characteristic 2 or when a needed pivot is a nonunit.
            fn diagonalize(&self) -> PyResult<$alg> {
                let metric = scalar_boundary(|| crate::forms::diagonalize(&self.inner.metric))?
                    .ok_or_else(|| {
                        PyValueError::new_err(
                            "metric is not diagonalizable in this scalar world",
                        )
                    })?;
                Ok($alg {
                    inner: Arc::new(CliffordAlgebra::new(metric.q.len(), metric)),
                })
            }

            /// Return this metric unchanged if already diagonal, otherwise
            /// congruence-diagonalize it.
            fn as_diagonal(&self) -> PyResult<$alg> {
                let metric = scalar_boundary(|| crate::forms::as_diagonal(&self.inner.metric))?
                    .ok_or_else(|| {
                        PyValueError::new_err(
                            "metric is not diagonalizable in this scalar world",
                        )
                    })?;
                Ok($alg {
                    inner: Arc::new(CliffordAlgebra::new(metric.q.len(), metric)),
                })
            }

            /// The determinant of a `LinearMap`: the scalar by which its
            /// outermorphism scales the pseudoscalar. Char-faithful (the char-2
            /// determinant over nimbers).
            fn determinant(&self, lm: &$lm) -> PyResult<$scalar_py> {
                self.ensure_linear_map(&lm.inner)?;
                Ok($wrap(scalar_boundary(|| {
                    crate::clifford::determinant(&self.inner, &lm.inner)
                })?))
            }

            /// The trace of a `LinearMap` (`= tr Λ¹f`).
            fn trace(&self, lm: &$lm) -> PyResult<$scalar_py> {
                self.ensure_linear_map(&lm.inner)?;
                Ok($wrap(scalar_boundary(|| {
                    crate::clifford::trace(&self.inner, &lm.inner)
                })?))
            }

            /// The trace of the exterior power `Λ^k f`.
            fn exterior_power_trace(&self, lm: &$lm, k: usize) -> PyResult<$scalar_py> {
                self.ensure_linear_map(&lm.inner)?;
                Ok($wrap(scalar_boundary(|| {
                    crate::clifford::exterior_power_trace(&self.inner, &lm.inner, k)
                })?))
            }

            /// The characteristic polynomial `det(t·I − f)` via exterior-power
            /// traces, as coefficients in descending degree `[1, −c₁, …, (−1)ⁿcₙ]`
            /// (`cₖ = tr Λᵏf`). Char-faithful.
            fn char_poly(&self, lm: &$lm) -> PyResult<Vec<$scalar_py>> {
                self.ensure_linear_map(&lm.inner)?;
                Ok(scalar_boundary(|| crate::clifford::char_poly(&self.inner, &lm.inner))?
                    .into_iter()
                    .map($wrap)
                    .collect())
            }

            /// The inverse `LinearMap`, if it is invertible over this scalar world.
            fn inverse_outermorphism(&self, lm: &$lm) -> PyResult<Option<$lm>> {
                self.ensure_linear_map(&lm.inner)?;
                Ok(scalar_boundary(|| crate::clifford::inverse_outermorphism(&lm.inner))?
                    .map(|inner| $lm { inner }))
            }

            /// Apply the outermorphism of a `LinearMap` to a multivector:
            /// `f(a∧b) = f(a)∧f(b)`.
            fn apply_outermorphism(&self, lm: &$lm, mv: &$mv) -> PyResult<$mv> {
                self.ensure_mv(mv)?;
                self.ensure_linear_map(&lm.inner)?;
                Ok($mv {
                    alg: self.inner.clone(),
                    mv: scalar_boundary(|| {
                        crate::clifford::apply_outermorphism(&self.inner, &lm.inner, &mv.mv)
                    })?,
                })
            }

            /// Full concrete spinor data as a named `SpinorRep` record.
            /// Supports nondegenerate characteristic-0 metrics and nonsingular
            /// characteristic-2 nimber metrics; rejects general-bilinear metrics.
            /// `diagonalized_metric` is returned as `(q, b_terms)` when present,
            /// where `b_terms` contains `(i, j, value)` entries.
            fn spinor_rep(&self, py: Python<'_>) -> PyResult<PySpinorRep> {
                let rep = scalar_boundary(|| crate::clifford::spinor_rep(&self.inner))?.ok_or_else(|| {
	                    PyValueError::new_err(
	                        "spinor_rep needs a supported nondegenerate metric with no general-bilinear a-part",
	                    )
                })?;
                let is_left_regular = rep.is_left_regular;
                let diagonalized_metric: Option<(
                    Vec<$scalar_py>,
                    Vec<(usize, usize, $scalar_py)>,
                )> = rep.diagonalized_metric.map(|metric| {
                    (
                        metric.q.into_iter().map($wrap).collect(),
                        metric
                            .b
                            .into_iter()
                            .map(|((i, j), coeff)| (i, j, $wrap(coeff)))
                            .collect(),
                    )
                });
                let orthogonal_basis_in_original: Option<Vec<Vec<$scalar_py>>> =
                    rep.orthogonal_basis_in_original.map(|matrix| {
                        matrix
                            .into_iter()
                            .map(|row| row.into_iter().map($wrap).collect())
                            .collect()
                    });
                let idempotent = $mv {
                    alg: self.inner.clone(),
                    mv: rep.idempotent,
                };
                let basis: Vec<$mv> = rep
                    .basis
                    .into_iter()
                    .map(|mv| $mv {
                        alg: self.inner.clone(),
                        mv,
                    })
                    .collect();
                let gen_matrices: Vec<Vec<Vec<$scalar_py>>> = rep
                    .gen_matrices
                    .into_iter()
                    .map(|m| {
                        m.into_iter()
                            .map(|row| row.into_iter().map($wrap).collect())
                            .collect()
                    })
                    .collect();
                let basis_dim = basis.len();
                let generator_count = gen_matrices.len();
                Ok(PySpinorRep {
                    idempotent: idempotent.into_py_any(py)?,
                    basis: basis.into_py_any(py)?,
                    gen_matrices: gen_matrices.into_py_any(py)?,
                    is_left_regular,
                    diagonalized_metric: diagonalized_metric.into_py_any(py)?,
                    orthogonal_basis_in_original: orthogonal_basis_in_original.into_py_any(py)?,
                    basis_dim,
                    generator_count,
                })
            }

            /// Apply the lazy left-regular spinor action of generator `e_i` to a
            /// sparse module vector. This reaches dimensions where explicit
            /// `spinor_rep()` matrices are intentionally capped.
            fn apply_generator(&self, i: usize, v: &$mv) -> PyResult<$mv> {
                self.ensure_mv(v)?;
                let rep = scalar_boundary(|| crate::clifford::lazy_spinor_rep(&self.inner))?
                    .ok_or_else(|| {
                        PyValueError::new_err(
                            "lazy_spinor_rep needs a supported nondegenerate metric with no general-bilinear a-part",
                        )
                    })?;
                let mv = scalar_boundary(|| rep.apply_generator(i, &v.mv))?
                    .ok_or_else(|| PyValueError::new_err("generator index out of range"))?;
                Ok($mv {
                    alg: self.inner.clone(),
                    mv,
                })
            }

            /// Build the Rust `LazySpinorRep` façade for this backend.
            fn lazy_spinor_rep(&self, py: Python<'_>) -> PyResult<PyLazySpinorRep> {
                scalar_boundary(|| crate::clifford::lazy_spinor_rep(&self.inner))?
                    .ok_or_else(|| {
                        PyValueError::new_err(
                            "lazy_spinor_rep needs a supported nondegenerate metric with no general-bilinear a-part",
                        )
                    })?;
                Ok(PyLazySpinorRep {
                    algebra: self.clone().into_py_any(py)?,
                })
            }

            /// Apply the lazy left-regular spinor action of a vector
            /// `Σ coeffs[i] e_i` to a sparse module vector.
            fn apply_vector(
                &self,
                coeffs: Vec<Bound<'_, PyAny>>,
                v: &$mv,
            ) -> PyResult<$mv> {
                self.ensure_mv(v)?;
                let mut parsed = Vec::with_capacity(coeffs.len());
                for coeff in &coeffs {
                    parsed.push($parse(coeff)?);
                }
                let rep = scalar_boundary(|| crate::clifford::lazy_spinor_rep(&self.inner))?
                    .ok_or_else(|| {
                        PyValueError::new_err(
                            "lazy_spinor_rep needs a supported nondegenerate metric with no general-bilinear a-part",
                        )
                    })?;
                let mv = scalar_boundary(|| rep.apply_vector(&parsed, &v.mv))?
                    .ok_or_else(|| PyValueError::new_err("coefficient length must equal algebra dimension"))?;
                Ok($mv {
                    alg: self.inner.clone(),
                    mv,
                })
            }

            fn __repr__(&self) -> String {
                format!("{}(dim={})", $alg_name, self.inner.dim)
            }
        }

        impl $alg {
            fn ensure_mv(&self, mv: &$mv) -> PyResult<()> {
                if self.inner.as_ref() == mv.alg.as_ref() {
                    Ok(())
                } else {
                    Err(PyValueError::new_err(
                        "multivector belongs to a different Clifford algebra",
                    ))
                }
            }

            fn ensure_linear_map(&self, lm: &LinearMap<$scalar>) -> PyResult<()> {
                if lm.n != self.inner.dim {
                    return Err(PyValueError::new_err(format!(
                        "linear-map dimension {} does not match algebra dimension {}",
                        lm.n,
                        self.inner.dim
                    )));
                }
                Ok(())
            }
        }

    };
}

// Multivector wrapper and multivector-level operations for one backend.
macro_rules! backend_multivector {
    (
        $alg:ident,
        $alg_name:literal,
        $mv:ident,
        $mv_name:literal,
        $lm:ident,
        $lm_name:literal,
        $scalar:ty,
        $parse:path,
        $scalar_py:ty,
        $wrap:path
    ) => {
        #[pyclass(name = $mv_name, module = "pleroma", from_py_object)]
        #[derive(Clone)]
        pub(crate) struct $mv {
            pub(crate) alg: Arc<CliffordAlgebra<$scalar>>,
            pub(crate) mv: Multivector<$scalar>,
        }

        impl $mv {
            fn ensure_same_algebra(&self, other: &$mv) -> PyResult<()> {
                if self.alg.as_ref() == other.alg.as_ref() {
                    Ok(())
                } else {
                    Err(PyValueError::new_err(
                        "multivectors belong to different Clifford algebras",
                    ))
                }
            }
        }

        #[pymethods]
        impl $mv {
            fn __add__(&self, other: &$mv) -> PyResult<$mv> {
                self.ensure_same_algebra(other)?;
                Ok($mv {
                    alg: self.alg.clone(),
                    mv: self.alg.add(&self.mv, &other.mv),
                })
            }
            fn __sub__(&self, other: &$mv) -> PyResult<$mv> {
                self.ensure_same_algebra(other)?;
                let neg_one = <$scalar as Scalar>::one().neg();
                let neg = scalar_boundary(|| self.alg.scalar_mul(&neg_one, &other.mv))?;
                Ok($mv {
                    alg: self.alg.clone(),
                    mv: scalar_boundary(|| self.alg.add(&self.mv, &neg))?,
                })
            }
            fn __neg__(&self) -> PyResult<$mv> {
                let neg_one = <$scalar as Scalar>::one().neg();
                Ok($mv {
                    alg: self.alg.clone(),
                    mv: scalar_boundary(|| self.alg.scalar_mul(&neg_one, &self.mv))?,
                })
            }
            fn __mul__(&self, other: &Bound<'_, PyAny>) -> PyResult<$mv> {
                if let Ok(o) = other.cast::<$mv>() {
                    let other = o.borrow();
                    self.ensure_same_algebra(&other)?;
                    return Ok($mv {
                        alg: self.alg.clone(),
                        mv: scalar_boundary(|| self.alg.mul(&self.mv, &other.mv))?,
                    });
                }
                let s = $parse(other)?;
                Ok($mv {
                    alg: self.alg.clone(),
                    mv: scalar_boundary(|| self.alg.scalar_mul(&s, &self.mv))?,
                })
            }
            fn __rmul__(&self, other: &Bound<'_, PyAny>) -> PyResult<$mv> {
                let s = $parse(other)?;
                Ok($mv {
                    alg: self.alg.clone(),
                    mv: scalar_boundary(|| self.alg.scalar_mul(&s, &self.mv))?,
                })
            }
            fn __pow__(&self, n: u128, _modulo: Option<&Bound<'_, PyAny>>) -> PyResult<$mv> {
                let acc = scalar_boundary(|| {
                    let mut acc = self.alg.scalar(<$scalar as Scalar>::one());
                    let mut base = self.mv.clone();
                    let mut e = n;
                    while e > 0 {
                        if e & 1 == 1 {
                            acc = self.alg.mul(&acc, &base);
                        }
                        e >>= 1;
                        if e > 0 {
                            base = self.alg.mul(&base, &base);
                        }
                    }
                    acc
                })?;
                Ok($mv {
                    alg: self.alg.clone(),
                    mv: acc,
                })
            }
            /// Exterior (wedge) product; also bound to the `^` operator.
            fn wedge(&self, other: &$mv) -> PyResult<$mv> {
                self.ensure_same_algebra(other)?;
                Ok($mv {
                    alg: self.alg.clone(),
                    mv: scalar_boundary(|| self.alg.wedge(&self.mv, &other.mv))?,
                })
            }
            fn __xor__(&self, other: &$mv) -> PyResult<$mv> {
                self.wedge(other)
            }
            fn reverse(&self) -> PyResult<$mv> {
                Ok($mv {
                    alg: self.alg.clone(),
                    mv: scalar_boundary(|| self.alg.reverse(&self.mv))?,
                })
            }
            /// `~v` is reversion.
            fn __invert__(&self) -> PyResult<$mv> {
                self.reverse()
            }
            fn grade_part(&self, k: usize) -> $mv {
                $mv {
                    alg: self.alg.clone(),
                    mv: self.alg.grade_part(&self.mv, k),
                }
            }
            fn grade_involution(&self) -> PyResult<$mv> {
                Ok($mv {
                    alg: self.alg.clone(),
                    mv: scalar_boundary(|| self.alg.grade_involution(&self.mv))?,
                })
            }
            /// Versor inverse v⁻¹ = ṽ/(v ṽ); errors if v isn't an invertible versor.
            fn versor_inverse(&self) -> PyResult<$mv> {
                scalar_boundary(|| self.alg.versor_inverse(&self.mv))?
                    .map(|mv| $mv {
                        alg: self.alg.clone(),
                        mv,
                    })
                    .ok_or_else(|| PyValueError::new_err("not an invertible versor"))
            }
            /// The **general multivector inverse** (any invertible element, not
            /// just a versor) via the left-multiplication matrix. Errors if the
            /// element is a zero divisor / non-invertible.
            fn multivector_inverse(&self) -> PyResult<$mv> {
                scalar_boundary(|| self.alg.multivector_inverse(&self.mv))?
                    .map(|mv| $mv {
                        alg: self.alg.clone(),
                        mv,
                    })
                    .ok_or_else(|| PyValueError::new_err("not invertible (zero divisor)"))
            }
            /// The **Cayley transform** `(1−B)(1+B)⁻¹` of this bivector — the exact
            /// rational map from the Lie algebra (bivectors) to the Spin group
            /// (rotors). Errors if `1+B` is not invertible.
            fn cayley(&self) -> PyResult<$mv> {
                scalar_boundary(|| self.alg.cayley(&self.mv))?
                    .map(|mv| $mv {
                        alg: self.alg.clone(),
                        mv,
                    })
                    .ok_or_else(|| PyValueError::new_err("1+B not invertible"))
            }
            /// The inverse Cayley transform — a rotor back to its bivector
            /// generator (same involutive formula). Errors if `1+R` is singular.
            fn cayley_inverse(&self) -> PyResult<$mv> {
                scalar_boundary(|| self.alg.cayley_inverse(&self.mv))?
                    .map(|mv| $mv {
                        alg: self.alg.clone(),
                        mv,
                    })
                    .ok_or_else(|| PyValueError::new_err("1+R not invertible"))
            }
            /// Sandwich self · x · self⁻¹ (rotor/versor action; untwisted).
            fn sandwich(&self, x: &$mv) -> PyResult<$mv> {
                self.ensure_same_algebra(x)?;
                scalar_boundary(|| self.alg.sandwich(&self.mv, &x.mv))?
                    .map(|mv| $mv {
                        alg: self.alg.clone(),
                        mv,
                    })
                    .ok_or_else(|| PyValueError::new_err("not an invertible versor"))
            }
            /// Twisted adjoint (Pin/Spin action) α(self) · x · self⁻¹ — the correct
            /// versor action; for an odd versor it gives a genuine reflection.
            fn twisted_sandwich(&self, x: &$mv) -> PyResult<$mv> {
                self.ensure_same_algebra(x)?;
                scalar_boundary(|| self.alg.twisted_sandwich(&self.mv, &x.mv))?
                    .map(|mv| $mv {
                        alg: self.alg.clone(),
                        mv,
                    })
                    .ok_or_else(|| PyValueError::new_err("not an invertible versor"))
            }
            /// Projection onto the even subalgebra (sum of even-grade blades).
            fn even_part(&self) -> PyResult<$mv> {
                Ok($mv {
                    alg: self.alg.clone(),
                    mv: scalar_boundary(|| self.alg.even_part(&self.mv))?,
                })
            }
            /// The exterior-Hopf coproduct Δ, returned as a multivector over the
            /// graded tensor square `Cl ⊗̂ Cl` (a tensor `e_T ⊗ e_U` is the blade
            /// `T | (U << dim)`).
            fn coproduct(&self) -> PyResult<$mv> {
                if self.alg.dim * 2 > MAX_BASIS_DIM {
                    return Err(PyValueError::new_err(format!(
                        "coproduct tensor encoding needs 2*dim <= {MAX_BASIS_DIM}"
                    )));
                }
                let tensor = self.alg.graded_tensor(&self.alg);
                let co = scalar_boundary(|| crate::clifford::coproduct(&self.alg, &self.mv))?;
                Ok($mv {
                    alg: Arc::new(tensor),
                    mv: co,
                })
            }
            /// The exterior-Hopf antipode (the grade involution `(−1)^k`).
            fn antipode(&self) -> PyResult<$mv> {
                Ok($mv {
                    alg: self.alg.clone(),
                    mv: scalar_boundary(|| crate::clifford::antipode(&self.alg, &self.mv))?,
                })
            }
            /// The exterior-Hopf counit (the scalar part).
            fn counit(&self) -> PyResult<$scalar_py> {
                Ok($wrap(scalar_boundary(|| {
                    crate::clifford::counit(&self.alg, &self.mv)
                })?))
            }
            /// `exp(self)` for a nilpotent multivector — the terminating series
            /// `Σ selfᵏ/k!`. Errors if `self` is not nilpotent (a rotational motor,
            /// needing transcendental cos/sin).
            fn exp_nilpotent(&self) -> PyResult<$mv> {
                scalar_boundary(|| crate::clifford::exp_nilpotent(&self.alg, &self.mv))?
                    .map(|mv| $mv {
                        alg: self.alg.clone(),
                        mv,
                    })
                    .ok_or_else(|| {
                        PyValueError::new_err("not nilpotent — would need a transcendental exp")
                    })
            }
            /// Reflect x in the hyperplane ⊥ self (self must be an invertible vector).
            fn reflect(&self, x: &$mv) -> PyResult<$mv> {
                self.ensure_same_algebra(x)?;
                scalar_boundary(|| self.alg.reflect(&self.mv, &x.mv))?
                    .map(|mv| $mv {
                        alg: self.alg.clone(),
                        mv,
                    })
                    .ok_or_else(|| PyValueError::new_err("not an invertible vector"))
            }
            fn left_contract(&self, other: &$mv) -> PyResult<$mv> {
                self.ensure_same_algebra(other)?;
                Ok($mv {
                    alg: self.alg.clone(),
                    mv: scalar_boundary(|| self.alg.left_contract(&self.mv, &other.mv))?,
                })
            }
            fn right_contract(&self, other: &$mv) -> PyResult<$mv> {
                self.ensure_same_algebra(other)?;
                Ok($mv {
                    alg: self.alg.clone(),
                    mv: scalar_boundary(|| self.alg.right_contract(&self.mv, &other.mv))?,
                })
            }
            /// `<<` is left contraction, `>>` is right contraction.
            fn __lshift__(&self, other: &$mv) -> PyResult<$mv> {
                self.left_contract(other)
            }
            fn __rshift__(&self, other: &$mv) -> PyResult<$mv> {
                self.right_contract(other)
            }
            fn dual(&self) -> PyResult<$mv> {
                scalar_boundary(|| self.alg.dual(&self.mv))?
                    .map(|mv| $mv {
                        alg: self.alg.clone(),
                        mv,
                    })
                    .ok_or_else(|| {
                        PyValueError::new_err("pseudoscalar not invertible (degenerate metric)")
                    })
            }
            /// The undual v ↦ v·I (inverse of `dual`).
            fn undual(&self) -> PyResult<$mv> {
                Ok($mv {
                    alg: self.alg.clone(),
                    mv: scalar_boundary(|| self.alg.undual(&self.mv))?,
                })
            }
            /// The Clifford (main) conjugate: reversion ∘ grade involution.
            fn clifford_conjugate(&self) -> PyResult<$mv> {
                Ok($mv {
                    alg: self.alg.clone(),
                    mv: scalar_boundary(|| self.alg.clifford_conjugate(&self.mv))?,
                })
            }
            /// The scalar product ⟨a b⟩₀ (grade-0 part of the geometric product).
            fn scalar_product(&self, other: &$mv) -> PyResult<$scalar_py> {
                self.ensure_same_algebra(other)?;
                Ok($wrap(scalar_boundary(|| {
                    self.alg.scalar_product(&self.mv, &other.mv)
                })?))
            }
            /// The commutator product [a,b] = ab − ba (no ½; char-faithful).
            fn commutator(&self, other: &$mv) -> PyResult<$mv> {
                self.ensure_same_algebra(other)?;
                Ok($mv {
                    alg: self.alg.clone(),
                    mv: scalar_boundary(|| self.alg.commutator(&self.mv, &other.mv))?,
                })
            }
            /// The anticommutator product {a,b} = ab + ba (no ½; char-faithful).
            fn anticommutator(&self, other: &$mv) -> PyResult<$mv> {
                self.ensure_same_algebra(other)?;
                Ok($mv {
                    alg: self.alg.clone(),
                    mv: scalar_boundary(|| self.alg.anticommutator(&self.mv, &other.mv))?,
                })
            }
            /// The regressive (meet) product a ∨ b — intersection dual to the
            /// wedge. Errors if the pseudoscalar is not invertible.
            fn meet(&self, other: &$mv) -> PyResult<$mv> {
                self.ensure_same_algebra(other)?;
                scalar_boundary(|| self.alg.meet(&self.mv, &other.mv))?
                    .map(|mv| $mv {
                        alg: self.alg.clone(),
                        mv,
                    })
                    .ok_or_else(|| {
                        PyValueError::new_err("pseudoscalar not invertible (degenerate metric)")
                    })
            }
            /// Whether this multivector is a blade (a decomposable homogeneous
            /// element — a wedge of vectors).
            fn is_blade(&self) -> bool {
                crate::clifford::is_blade(&self.alg, &self.mv)
            }
            /// A basis of the blade subspace `{x : x ∧ A = 0}`, as coefficient
            /// rows over the algebra generators `e0, e1, ...`. Scalars return an
            /// empty basis; errors for zero or mixed-grade multivectors.
            fn blade_subspace(&self) -> PyResult<Vec<Vec<$scalar_py>>> {
                scalar_boundary(|| crate::clifford::blade_subspace(&self.alg, &self.mv))?
                    .map(|basis| {
                        basis
                            .into_iter()
                            .map(|row| row.into_iter().map($wrap).collect())
                            .collect()
                    })
                    .ok_or_else(|| {
                        PyValueError::new_err(
                            "blade_subspace needs a nonzero homogeneous multivector",
                        )
                    })
            }
            /// Factor a blade into the grade-1 vectors whose wedge is it; errors
            /// if it is not a blade.
            fn factor_blade(&self) -> PyResult<Vec<$mv>> {
                scalar_boundary(|| crate::clifford::factor_blade(&self.alg, &self.mv))?
                    .map(|vs| {
                        vs.into_iter()
                            .map(|mv| $mv {
                                alg: self.alg.clone(),
                                mv,
                            })
                            .collect()
                    })
                    .ok_or_else(|| PyValueError::new_err("not a blade (not decomposable)"))
            }
            fn norm2(&self) -> PyResult<$scalar_py> {
                Ok($wrap(scalar_boundary(|| self.alg.norm2(&self.mv))?))
            }
            /// The Dickson / grade parity of a homogeneous-parity versor candidate:
            /// `0` for even, `1` for odd, `None` for zero or mixed parity.
            fn versor_grade_parity(&self) -> Option<u128> {
                crate::clifford::versor_grade_parity(&self.mv)
            }
            /// Raw spinor norm `<v reverse(v)>_0`; errors when `v` is not an
            /// invertible simple versor. Reduce this scalar modulo squares (char != 2)
            /// or Artin-Schreier (char 2) in the caller's field when needed.
            fn spinor_norm(&self) -> PyResult<$scalar_py> {
                scalar_boundary(|| self.alg.spinor_norm(&self.mv))?
                    .map($wrap)
                    .ok_or_else(|| PyValueError::new_err("not an invertible simple versor"))
            }
            /// Classify a versor as a named `VersorClass` record.
            fn classify_versor(&self, py: Python<'_>) -> PyResult<PyVersorClass> {
                let class = scalar_boundary(|| self.alg.classify_versor(&self.mv))?
                    .ok_or_else(|| PyValueError::new_err("not an invertible simple versor"))?;
                Ok(PyVersorClass {
                    spinor_norm: $wrap(class.spinor_norm).into_py_any(py)?,
                    dickson: class.dickson,
                })
            }
            fn scalar_part(&self) -> PyResult<$scalar_py> {
                Ok($wrap(scalar_boundary(|| self.alg.scalar_part(&self.mv))?))
            }
            /// Division: by a scalar, or by a versor (multiply by its inverse).
            fn __truediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<$mv> {
                if let Ok(o) = other.cast::<$mv>() {
                    let other = o.borrow();
                    self.ensure_same_algebra(&other)?;
                    let oinv = scalar_boundary(|| self.alg.versor_inverse(&other.mv))?
                        .ok_or_else(|| PyValueError::new_err("divisor not an invertible versor"))?;
                    return Ok($mv {
                        alg: self.alg.clone(),
                        mv: scalar_boundary(|| self.alg.mul(&self.mv, &oinv))?,
                    });
                }
                let s = $parse(other)?;
                let sinv = <$scalar as Scalar>::inv(&s)
                    .ok_or_else(|| PyValueError::new_err("scalar has no representable inverse"))?;
                Ok($mv {
                    alg: self.alg.clone(),
                    mv: scalar_boundary(|| self.alg.scalar_mul(&sinv, &self.mv))?,
                })
            }
            #[getter]
            fn terms(&self) -> Vec<(u128, $scalar_py)> {
                self.mv
                    .terms
                    .iter()
                    .map(|(&mask, coeff)| (mask, $wrap(coeff.clone())))
                    .collect()
            }
            fn is_zero(&self) -> bool {
                self.mv.is_zero()
            }
            fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
                if let Ok(o) = other.cast::<$mv>() {
                    let other = o.borrow();
                    self.alg.as_ref() == other.alg.as_ref() && self.mv == other.mv
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

macro_rules! backend {
    (
        $alg:ident,
        $alg_name:literal,
        $mv:ident,
        $mv_name:literal,
        $lm:ident,
        $lm_name:literal,
        $scalar:ty,
        $parse:path,
        $scalar_py:ty,
        $wrap:path
    ) => {
        backend_linear_map!(
            $alg, $alg_name, $mv, $mv_name, $lm, $lm_name, $scalar, $parse, $scalar_py, $wrap
        );
        backend_algebra!(
            $alg, $alg_name, $mv, $mv_name, $lm, $lm_name, $scalar, $parse, $scalar_py, $wrap
        );
        backend_multivector!(
            $alg, $alg_name, $mv, $mv_name, $lm, $lm_name, $scalar, $parse, $scalar_py, $wrap
        );
    };
}
py_engine_backends!(backend);

macro_rules! divided_power_backend {
    ($alg:ident, $alg_name:literal, $vec:ident, $vec_name:literal, $scalar:ty, $parse:path, $scalar_py:ty, $wrap:path) => {
        #[pyclass(name = $alg_name, module = "pleroma", from_py_object)]
        #[derive(Clone)]
        struct $alg {
            inner: Arc<DividedPowerAlgebra>,
        }

        #[pyclass(name = $vec_name, module = "pleroma", from_py_object)]
        #[derive(Clone)]
        struct $vec {
            alg: Arc<DividedPowerAlgebra>,
            vec: DpVector<$scalar>,
        }

        impl $alg {
            fn wrap(&self, vec: DpVector<$scalar>) -> $vec {
                $vec {
                    alg: self.inner.clone(),
                    vec,
                }
            }

            fn ensure_vec(&self, x: &$vec) -> PyResult<()> {
                if self.inner.as_ref() == x.alg.as_ref() {
                    Ok(())
                } else {
                    Err(PyValueError::new_err(
                        "divided-power vector belongs to a different algebra",
                    ))
                }
            }
        }

        #[pymethods]
        impl $alg {
            #[new]
            fn new(dim: usize) -> Self {
                $alg {
                    inner: Arc::new(DividedPowerAlgebra::new(dim)),
                }
            }
            #[getter]
            fn dim(&self) -> usize {
                self.inner.dim
            }
            fn zero(&self) -> $vec {
                self.wrap(self.inner.zero::<$scalar>())
            }
            fn one(&self) -> $vec {
                self.wrap(self.inner.one::<$scalar>())
            }
            fn scalar(&self, s: &Bound<'_, PyAny>) -> PyResult<$vec> {
                Ok(self.wrap(self.inner.scalar::<$scalar>($parse(s)?)))
            }
            fn divided_power(&self, i: usize, k: u128) -> PyResult<$vec> {
                if i >= self.inner.dim {
                    return Err(PyValueError::new_err("generator index out of range"));
                }
                Ok(self.wrap(self.inner.divided_power::<$scalar>(i, k)))
            }
            fn gen(&self, i: usize) -> PyResult<$vec> {
                if i >= self.inner.dim {
                    return Err(PyValueError::new_err("generator index out of range"));
                }
                Ok(self.wrap(self.inner.gen::<$scalar>(i)))
            }
            fn monomial(&self, alpha: Vec<u128>, coeff: &Bound<'_, PyAny>) -> PyResult<$vec> {
                if alpha.len() > self.inner.dim {
                    return Err(PyValueError::new_err("multidegree longer than dim"));
                }
                Ok(self.wrap(self.inner.monomial::<$scalar>(&alpha, $parse(coeff)?)))
            }
            fn add(&self, x: &$vec, y: &$vec) -> PyResult<$vec> {
                self.ensure_vec(x)?;
                self.ensure_vec(y)?;
                Ok(self.wrap(self.inner.add(&x.vec, &y.vec)))
            }
            fn scalar_mul(&self, s: &Bound<'_, PyAny>, x: &$vec) -> PyResult<$vec> {
                self.ensure_vec(x)?;
                let s = $parse(s)?;
                Ok(self.wrap(scalar_boundary(|| self.inner.scalar_mul(&s, &x.vec))?))
            }
            fn mul(&self, x: &$vec, y: &$vec) -> PyResult<$vec> {
                self.ensure_vec(x)?;
                self.ensure_vec(y)?;
                Ok(self.wrap(scalar_boundary(|| self.inner.mul(&x.vec, &y.vec))?))
            }
            fn coproduct(&self, x: &$vec) -> PyResult<Vec<(Vec<u128>, Vec<u128>, $scalar_py)>> {
                self.ensure_vec(x)?;
                Ok(self
                    .inner
                    .coproduct(&x.vec)
                    .into_iter()
                    .map(|((left, right), coeff)| (left, right, $wrap(coeff)))
                    .collect())
            }
            fn counit(&self, x: &$vec) -> PyResult<$scalar_py> {
                self.ensure_vec(x)?;
                Ok($wrap(self.inner.counit(&x.vec)))
            }
            fn antipode(&self, x: &$vec) -> PyResult<$vec> {
                self.ensure_vec(x)?;
                Ok(self.wrap(scalar_boundary(|| self.inner.antipode(&x.vec))?))
            }
            fn __repr__(&self) -> String {
                format!("{}(dim={})", $alg_name, self.inner.dim)
            }
        }

        impl $vec {
            fn ensure_same_algebra(&self, other: &$vec) -> PyResult<()> {
                if self.alg.as_ref() == other.alg.as_ref() {
                    Ok(())
                } else {
                    Err(PyValueError::new_err(
                        "divided-power vectors belong to different algebras",
                    ))
                }
            }

            fn wrap(&self, vec: DpVector<$scalar>) -> $vec {
                $vec {
                    alg: self.alg.clone(),
                    vec,
                }
            }
        }

        #[pymethods]
        impl $vec {
            #[getter]
            fn terms(&self) -> Vec<(Vec<u128>, $scalar_py)> {
                self.vec
                    .terms
                    .iter()
                    .map(|(degree, coeff)| (degree.clone(), $wrap(coeff.clone())))
                    .collect()
            }
            fn is_zero(&self) -> bool {
                self.vec.terms.is_empty()
            }
            fn __add__(&self, other: &$vec) -> PyResult<$vec> {
                self.ensure_same_algebra(other)?;
                Ok(self.wrap(self.alg.add(&self.vec, &other.vec)))
            }
            fn __sub__(&self, other: &$vec) -> PyResult<$vec> {
                self.ensure_same_algebra(other)?;
                let neg_one = <$scalar as Scalar>::one().neg();
                let neg = scalar_boundary(|| self.alg.scalar_mul(&neg_one, &other.vec))?;
                Ok(self.wrap(scalar_boundary(|| self.alg.add(&self.vec, &neg))?))
            }
            fn __neg__(&self) -> PyResult<$vec> {
                let neg_one = <$scalar as Scalar>::one().neg();
                Ok(self.wrap(scalar_boundary(|| {
                    self.alg.scalar_mul(&neg_one, &self.vec)
                })?))
            }
            fn __mul__(&self, other: &$vec) -> PyResult<$vec> {
                self.ensure_same_algebra(other)?;
                Ok(self.wrap(scalar_boundary(|| self.alg.mul(&self.vec, &other.vec))?))
            }
            fn scale(&self, s: &Bound<'_, PyAny>) -> PyResult<$vec> {
                let s = $parse(s)?;
                Ok(self.wrap(scalar_boundary(|| self.alg.scalar_mul(&s, &self.vec))?))
            }
            fn __rmul__(&self, s: &Bound<'_, PyAny>) -> PyResult<$vec> {
                self.scale(s)
            }
            fn coproduct(&self) -> Vec<(Vec<u128>, Vec<u128>, $scalar_py)> {
                self.alg
                    .coproduct(&self.vec)
                    .into_iter()
                    .map(|((left, right), coeff)| (left, right, $wrap(coeff)))
                    .collect()
            }
            fn counit(&self) -> $scalar_py {
                $wrap(self.alg.counit(&self.vec))
            }
            fn antipode(&self) -> PyResult<$vec> {
                Ok(self.wrap(scalar_boundary(|| self.alg.antipode(&self.vec))?))
            }
            fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
                if let Ok(o) = other.cast::<$vec>() {
                    let other = o.borrow();
                    self.alg.as_ref() == other.alg.as_ref() && self.vec == other.vec
                } else {
                    false
                }
            }
            fn __repr__(&self) -> String {
                format!("{:?}", self.vec)
            }
        }
    };
}

py_divided_power_backends!(divided_power_backend);

// ---------------------------------------------------------------------------
// Conformal geometric algebra over bound characteristic-zero worlds
// ---------------------------------------------------------------------------

macro_rules! cga_backend {
    ($py:ident, $name:literal, $scalar:ty, $mv:ident, $scalar_py:ty, $parse:path, $wrap:path) => {
        #[pyclass(name = $name, module = "pleroma")]
        struct $py {
            inner: Cga<$scalar>,
        }

        impl $py {
            fn wrap(&self, mv: Multivector<$scalar>) -> $mv {
                $mv {
                    alg: Arc::new(self.inner.alg.clone()),
                    mv,
                }
            }
        }

        #[pymethods]
        impl $py {
            #[new]
            fn new(n: usize) -> Self {
                $py { inner: Cga::new(n) }
            }
            #[getter]
            fn n(&self) -> usize {
                self.inner.n
            }
            #[getter]
            fn dim(&self) -> usize {
                self.inner.alg.dim
            }
            fn n_o(&self) -> $mv {
                self.wrap(self.inner.n_o())
            }
            fn n_inf(&self) -> $mv {
                self.wrap(self.inner.n_inf())
            }
            /// Lift a Euclidean point to the null cone: `up(p) = n_o + p + ½|p|² n_∞`.
            fn up(&self, p: Vec<Bound<'_, PyAny>>) -> PyResult<$mv> {
                let mut pv = Vec::with_capacity(p.len());
                for x in &p {
                    pv.push($parse(x)?);
                }
                Ok(self.wrap(self.inner.up(&pv)))
            }
            /// Recover a Euclidean point from a null vector (`None` if not normalizable).
            fn down(&self, x: &$mv) -> Option<Vec<$scalar_py>> {
                self.inner
                    .down(&x.mv)
                    .map(|v| v.into_iter().map($wrap).collect())
            }
            /// The conformal inner product `x · y` (= `−½|p−q|²` on lifted points).
            fn inner(&self, x: &$mv, y: &$mv) -> $scalar_py {
                $wrap(self.inner.inner(&x.mv, &y.mv))
            }
            /// The sphere of squared radius `r2` about center `c`.
            fn sphere(&self, c: Vec<Bound<'_, PyAny>>, r2: &Bound<'_, PyAny>) -> PyResult<$mv> {
                let mut cv = Vec::with_capacity(c.len());
                for x in &c {
                    cv.push($parse(x)?);
                }
                Ok(self.wrap(self.inner.sphere(&cv, &$parse(r2)?)))
            }
            /// The plane `{x : x·normal = d}`.
            fn plane(&self, normal: Vec<Bound<'_, PyAny>>, d: &Bound<'_, PyAny>) -> PyResult<$mv> {
                let mut nv = Vec::with_capacity(normal.len());
                for x in &normal {
                    nv.push($parse(x)?);
                }
                Ok(self.wrap(self.inner.plane(&nv, &$parse(d)?)))
            }
            /// The point pair / oriented join `a ∧ b`.
            fn point_pair(&self, a: &$mv, b: &$mv) -> $mv {
                self.wrap(self.inner.point_pair(&a.mv, &b.mv))
            }
            /// The IPNS meet (intersection) `x ∧ y`.
            fn meet(&self, x: &$mv, y: &$mv) -> $mv {
                self.wrap(self.inner.meet(&x.mv, &y.mv))
            }
        }
    };
}

py_cga_backends!(cga_backend);

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PySpinorRep>()?;
    m.add_class::<PyLazySpinorRep>()?;
    m.add_class::<PyVersorClass>()?;
    macro_rules! register_backend {
        ($alg:ident, $alg_name:literal, $mv:ident, $mv_name:literal, $lm:ident, $lm_name:literal, $scalar:ty, $parse:path, $scalar_py:ty, $wrap:path) => {
            m.add_class::<$alg>()?;
            m.add_class::<$mv>()?;
            m.add_class::<$lm>()?;
        };
    }
    macro_rules! register_divided_power_backend {
        ($alg:ident, $alg_name:literal, $vec:ident, $vec_name:literal, $scalar:ty, $parse:path, $scalar_py:ty, $wrap:path) => {
            m.add_class::<$alg>()?;
            m.add_class::<$vec>()?;
        };
    }
    macro_rules! register_cga_backend {
        ($py:ident, $name:literal, $scalar:ty, $mv:ident, $scalar_py:ty, $parse:path, $wrap:path) => {
            m.add_class::<$py>()?;
        };
    }
    py_engine_backends!(register_backend);
    py_divided_power_backends!(register_divided_power_backend);
    py_cga_backends!(register_cga_backend);
    m.add_function(wrap_pyfunction!(galois_linear_map, m)?)?;
    m.add_function(wrap_pyfunction!(frobenius_linear_map, m)?)?;
    m.add_function(wrap_pyfunction!(nimber_subfield_frobenius_linear_map, m)?)?;
    m.add_function(wrap_pyfunction!(bits, m)?)?;
    m.add_function(wrap_pyfunction!(grade, m)?)?;
    Ok(())
}
