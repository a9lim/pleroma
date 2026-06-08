//! Python bindings for the GA engine: the `backend!` macro that stamps out one
//! `<World>Algebra` / `<World>MV` pyclass pair per scalar backend, the five
//! invocations, and conformal GA (`Cga`). The generated structs and their
//! fields are `pub(crate)` so the classifier bindings in [`super::forms`] and the
//! game-exterior binding in [`super::games`] can read `.inner` / `.mv`.

use super::scalars::{
    parse_integer, parse_nimber, parse_omnific, parse_surcomplex, parse_surreal, wrap_integer,
    wrap_nimber, wrap_omnific, wrap_surcomplex, wrap_surreal, PyInteger, PyNimber, PyOmnific,
    PySurcomplex, PySurreal,
};
use crate::clifford::{Cga, CliffordAlgebra, Metric, Multivector, MAX_BASIS_DIM};
use crate::scalar::{Integer, Nimber, Omnific, Scalar, Surcomplex, Surreal};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::collections::BTreeMap;
use std::sync::Arc;

// ---------------------------------------------------------------------------
// Algebra + multivector, one pair per backend
// ---------------------------------------------------------------------------

macro_rules! backend {
    ($alg:ident, $alg_name:literal, $mv:ident, $mv_name:literal, $scalar:ty, $parse:path, $scalar_py:ty, $wrap:path) => {
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

            /// The determinant of a linear map given column-major (`matrix[i]` =
            /// the image of `e_i`): the scalar by which its outermorphism scales
            /// the pseudoscalar. Char-faithful (the char-2 determinant over nimbers).
            fn determinant(&self, matrix: Vec<Vec<Bound<'_, PyAny>>>) -> PyResult<$scalar_py> {
                let lm = self.parse_linear_map(matrix)?;
                Ok($wrap(crate::clifford::determinant(&self.inner, &lm)))
            }

            /// The trace of a (column-major) linear map (`= tr Λ¹f`).
            fn trace(&self, matrix: Vec<Vec<Bound<'_, PyAny>>>) -> PyResult<$scalar_py> {
                let lm = self.parse_linear_map(matrix)?;
                Ok($wrap(crate::clifford::trace(&self.inner, &lm)))
            }

            /// The characteristic polynomial `det(t·I − f)` via exterior-power
            /// traces, as coefficients in descending degree `[1, −c₁, …, (−1)ⁿcₙ]`
            /// (`cₖ = tr Λᵏf`). Char-faithful.
            fn char_poly(&self, matrix: Vec<Vec<Bound<'_, PyAny>>>) -> PyResult<Vec<$scalar_py>> {
                let lm = self.parse_linear_map(matrix)?;
                Ok(crate::clifford::char_poly(&self.inner, &lm)
                    .into_iter()
                    .map($wrap)
                    .collect())
            }

            /// Apply the outermorphism of a (column-major) linear map to a
            /// multivector: `f(a∧b) = f(a)∧f(b)`.
            fn outermorphism(&self, matrix: Vec<Vec<Bound<'_, PyAny>>>, mv: &$mv) -> PyResult<$mv> {
                self.ensure_mv(mv)?;
                let lm = self.parse_linear_map(matrix)?;
                Ok($mv {
                    alg: self.inner.clone(),
                    mv: crate::clifford::apply_outermorphism(&self.inner, &lm, &mv.mv),
                })
            }

            /// A concrete spinor representation: `(idempotent, basis, gen_matrices)`
            /// realizing the classification on column spinors. Nondegenerate
            /// orthogonal char-0 metrics only.
            #[allow(clippy::type_complexity)]
            fn spinor_rep(&self) -> PyResult<($mv, Vec<$mv>, Vec<Vec<Vec<$scalar_py>>>)> {
                let rep = crate::clifford::spinor_rep(&self.inner).ok_or_else(|| {
                    PyValueError::new_err(
                        "spinor_rep needs a nondegenerate orthogonal characteristic-0 metric",
                    )
                })?;
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
                Ok((idempotent, basis, gen_matrices))
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

            fn parse_linear_map(
                &self,
                matrix: Vec<Vec<Bound<'_, PyAny>>>,
            ) -> PyResult<crate::clifford::LinearMap<$scalar>> {
                let n = matrix.len();
                if n != self.inner.dim {
                    return Err(PyValueError::new_err(format!(
                        "matrix dimension {n} does not match algebra dimension {}",
                        self.inner.dim
                    )));
                }
                let mut cols: Vec<Vec<$scalar>> = Vec::with_capacity(n);
                for col in &matrix {
                    if col.len() != n {
                        return Err(PyValueError::new_err(
                            "matrix must be square (n columns of length n)",
                        ));
                    }
                    let mut c = Vec::with_capacity(n);
                    for x in col {
                        c.push($parse(x)?);
                    }
                    cols.push(c);
                }
                Ok(crate::clifford::LinearMap::from_columns(cols))
            }
        }

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
                let neg = self.alg.scalar_mul(&neg_one, &other.mv);
                Ok($mv {
                    alg: self.alg.clone(),
                    mv: self.alg.add(&self.mv, &neg),
                })
            }
            fn __neg__(&self) -> $mv {
                let neg_one = <$scalar as Scalar>::one().neg();
                $mv {
                    alg: self.alg.clone(),
                    mv: self.alg.scalar_mul(&neg_one, &self.mv),
                }
            }
            fn __mul__(&self, other: &Bound<'_, PyAny>) -> PyResult<$mv> {
                if let Ok(o) = other.cast::<$mv>() {
                    let other = o.borrow();
                    self.ensure_same_algebra(&other)?;
                    return Ok($mv {
                        alg: self.alg.clone(),
                        mv: self.alg.mul(&self.mv, &other.mv),
                    });
                }
                let s = $parse(other)?;
                Ok($mv {
                    alg: self.alg.clone(),
                    mv: self.alg.scalar_mul(&s, &self.mv),
                })
            }
            fn __rmul__(&self, other: &Bound<'_, PyAny>) -> PyResult<$mv> {
                let s = $parse(other)?;
                Ok($mv {
                    alg: self.alg.clone(),
                    mv: self.alg.scalar_mul(&s, &self.mv),
                })
            }
            fn __pow__(&self, n: u128, _modulo: Option<&Bound<'_, PyAny>>) -> $mv {
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
                $mv {
                    alg: self.alg.clone(),
                    mv: acc,
                }
            }
            /// Exterior (wedge) product; also bound to the `^` operator.
            fn wedge(&self, other: &$mv) -> PyResult<$mv> {
                self.ensure_same_algebra(other)?;
                Ok($mv {
                    alg: self.alg.clone(),
                    mv: self.alg.wedge(&self.mv, &other.mv),
                })
            }
            fn __xor__(&self, other: &$mv) -> PyResult<$mv> {
                self.wedge(other)
            }
            fn reverse(&self) -> $mv {
                $mv {
                    alg: self.alg.clone(),
                    mv: self.alg.reverse(&self.mv),
                }
            }
            /// `~v` is reversion.
            fn __invert__(&self) -> $mv {
                self.reverse()
            }
            fn grade(&self, k: usize) -> $mv {
                $mv {
                    alg: self.alg.clone(),
                    mv: self.alg.grade_part(&self.mv, k),
                }
            }
            fn grade_involution(&self) -> $mv {
                $mv {
                    alg: self.alg.clone(),
                    mv: self.alg.grade_involution(&self.mv),
                }
            }
            /// Versor inverse v⁻¹ = ṽ/(v ṽ); errors if v isn't an invertible versor.
            fn inverse(&self) -> PyResult<$mv> {
                self.alg
                    .versor_inverse(&self.mv)
                    .map(|mv| $mv {
                        alg: self.alg.clone(),
                        mv,
                    })
                    .ok_or_else(|| PyValueError::new_err("not an invertible versor"))
            }
            /// The **general multivector inverse** (any invertible element, not
            /// just a versor) via the left-multiplication matrix. Errors if the
            /// element is a zero divisor / non-invertible.
            fn inverse_general(&self) -> PyResult<$mv> {
                self.alg
                    .multivector_inverse(&self.mv)
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
                self.alg
                    .cayley(&self.mv)
                    .map(|mv| $mv {
                        alg: self.alg.clone(),
                        mv,
                    })
                    .ok_or_else(|| PyValueError::new_err("1+B not invertible"))
            }
            /// The inverse Cayley transform — a rotor back to its bivector
            /// generator (same involutive formula). Errors if `1+R` is singular.
            fn cayley_inverse(&self) -> PyResult<$mv> {
                self.alg
                    .cayley_inverse(&self.mv)
                    .map(|mv| $mv {
                        alg: self.alg.clone(),
                        mv,
                    })
                    .ok_or_else(|| PyValueError::new_err("1+R not invertible"))
            }
            /// Sandwich self · x · self⁻¹ (rotor/versor action; untwisted).
            fn sandwich(&self, x: &$mv) -> PyResult<$mv> {
                self.ensure_same_algebra(x)?;
                self.alg
                    .sandwich(&self.mv, &x.mv)
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
                self.alg
                    .twisted_sandwich(&self.mv, &x.mv)
                    .map(|mv| $mv {
                        alg: self.alg.clone(),
                        mv,
                    })
                    .ok_or_else(|| PyValueError::new_err("not an invertible versor"))
            }
            /// Projection onto the even subalgebra (sum of even-grade blades).
            fn even_part(&self) -> $mv {
                $mv {
                    alg: self.alg.clone(),
                    mv: self.alg.even_part(&self.mv),
                }
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
                let co = crate::clifford::coproduct(&self.alg, &self.mv);
                Ok($mv {
                    alg: Arc::new(tensor),
                    mv: co,
                })
            }
            /// The exterior-Hopf antipode (the grade involution `(−1)^k`).
            fn antipode(&self) -> $mv {
                $mv {
                    alg: self.alg.clone(),
                    mv: crate::clifford::antipode(&self.alg, &self.mv),
                }
            }
            /// The exterior-Hopf counit (the scalar part).
            fn counit(&self) -> $scalar_py {
                $wrap(crate::clifford::counit(&self.alg, &self.mv))
            }
            /// `exp(self)` for a nilpotent multivector — the terminating series
            /// `Σ selfᵏ/k!`. Errors if `self` is not nilpotent (a rotational motor,
            /// needing transcendental cos/sin).
            fn exp_nilpotent(&self) -> PyResult<$mv> {
                crate::clifford::exp_nilpotent(&self.alg, &self.mv)
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
                self.alg
                    .reflect(&self.mv, &x.mv)
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
                    mv: self.alg.left_contract(&self.mv, &other.mv),
                })
            }
            fn right_contract(&self, other: &$mv) -> PyResult<$mv> {
                self.ensure_same_algebra(other)?;
                Ok($mv {
                    alg: self.alg.clone(),
                    mv: self.alg.right_contract(&self.mv, &other.mv),
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
                self.alg
                    .dual(&self.mv)
                    .map(|mv| $mv {
                        alg: self.alg.clone(),
                        mv,
                    })
                    .ok_or_else(|| {
                        PyValueError::new_err("pseudoscalar not invertible (degenerate metric)")
                    })
            }
            /// The undual v ↦ v·I (inverse of `dual`).
            fn undual(&self) -> $mv {
                $mv {
                    alg: self.alg.clone(),
                    mv: self.alg.undual(&self.mv),
                }
            }
            /// The Clifford (main) conjugate: reversion ∘ grade involution.
            fn clifford_conjugate(&self) -> $mv {
                $mv {
                    alg: self.alg.clone(),
                    mv: self.alg.clifford_conjugate(&self.mv),
                }
            }
            /// The scalar product ⟨a b⟩₀ (grade-0 part of the geometric product).
            fn scalar_product(&self, other: &$mv) -> PyResult<$scalar_py> {
                self.ensure_same_algebra(other)?;
                Ok($wrap(self.alg.scalar_product(&self.mv, &other.mv)))
            }
            /// The commutator product [a,b] = ab − ba (no ½; char-faithful).
            fn commutator(&self, other: &$mv) -> PyResult<$mv> {
                self.ensure_same_algebra(other)?;
                Ok($mv {
                    alg: self.alg.clone(),
                    mv: self.alg.commutator(&self.mv, &other.mv),
                })
            }
            /// The anticommutator product {a,b} = ab + ba (no ½; char-faithful).
            fn anticommutator(&self, other: &$mv) -> PyResult<$mv> {
                self.ensure_same_algebra(other)?;
                Ok($mv {
                    alg: self.alg.clone(),
                    mv: self.alg.anticommutator(&self.mv, &other.mv),
                })
            }
            /// The regressive (meet) product a ∨ b — intersection dual to the
            /// wedge. Errors if the pseudoscalar is not invertible.
            fn meet(&self, other: &$mv) -> PyResult<$mv> {
                self.ensure_same_algebra(other)?;
                self.alg
                    .meet(&self.mv, &other.mv)
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
            /// Factor a blade into the grade-1 vectors whose wedge is it; errors
            /// if it is not a blade.
            fn factor_blade(&self) -> PyResult<Vec<$mv>> {
                crate::clifford::factor_blade(&self.alg, &self.mv)
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
            fn norm2(&self) -> $scalar_py {
                $wrap(self.alg.norm2(&self.mv))
            }
            fn scalar_part(&self) -> $scalar_py {
                $wrap(self.alg.scalar_part(&self.mv))
            }
            /// Division: by a scalar, or by a versor (multiply by its inverse).
            fn __truediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<$mv> {
                if let Ok(o) = other.cast::<$mv>() {
                    let other = o.borrow();
                    self.ensure_same_algebra(&other)?;
                    let oinv = self
                        .alg
                        .versor_inverse(&other.mv)
                        .ok_or_else(|| PyValueError::new_err("divisor not an invertible versor"))?;
                    return Ok($mv {
                        alg: self.alg.clone(),
                        mv: self.alg.mul(&self.mv, &oinv),
                    });
                }
                let s = $parse(other)?;
                let sinv = <$scalar as Scalar>::inv(&s)
                    .ok_or_else(|| PyValueError::new_err("scalar has no representable inverse"))?;
                Ok($mv {
                    alg: self.alg.clone(),
                    mv: self.alg.scalar_mul(&sinv, &self.mv),
                })
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

backend!(
    NimberAlgebra,
    "NimberAlgebra",
    NimberMV,
    "NimberMV",
    Nimber,
    parse_nimber,
    PyNimber,
    wrap_nimber
);
backend!(
    SurrealAlgebra,
    "SurrealAlgebra",
    SurrealMV,
    "SurrealMV",
    Surreal,
    parse_surreal,
    PySurreal,
    wrap_surreal
);
backend!(
    SurcomplexAlgebra,
    "SurcomplexAlgebra",
    SurcomplexMV,
    "SurcomplexMV",
    Surcomplex<Surreal>,
    parse_surcomplex,
    PySurcomplex,
    wrap_surcomplex
);
// ℤ-coefficient backend: the home of the exterior algebra of the game group.
backend!(
    IntegerAlgebra,
    "IntegerAlgebra",
    IntegerMV,
    "IntegerMV",
    Integer,
    parse_integer,
    PyInteger,
    wrap_integer
);
// Omnific-integer backend: the surreal mirror of ℤ — exterior algebra over a
// transfinite ring (ω-scale coefficients).
backend!(
    OmnificAlgebra,
    "OmnificAlgebra",
    OmnificMV,
    "OmnificMV",
    Omnific,
    parse_omnific,
    PyOmnific,
    wrap_omnific
);

// ---------------------------------------------------------------------------
// Conformal geometric algebra over the surreals
// ---------------------------------------------------------------------------

#[pyclass(name = "Cga", module = "pleroma")]
struct PyCga {
    inner: Cga<Surreal>,
}

impl PyCga {
    fn wrap(&self, mv: Multivector<Surreal>) -> SurrealMV {
        SurrealMV {
            alg: Arc::new(self.inner.alg.clone()),
            mv,
        }
    }
}

#[pymethods]
impl PyCga {
    #[new]
    fn new(n: usize) -> Self {
        PyCga { inner: Cga::new(n) }
    }
    #[getter]
    fn n(&self) -> usize {
        self.inner.n
    }
    #[getter]
    fn dim(&self) -> usize {
        self.inner.alg.dim
    }
    fn n_o(&self) -> SurrealMV {
        self.wrap(self.inner.n_o())
    }
    fn n_inf(&self) -> SurrealMV {
        self.wrap(self.inner.n_inf())
    }
    /// Lift a Euclidean point to the null cone: `up(p) = n_o + p + ½|p|² n_∞`.
    fn up(&self, p: Vec<Bound<'_, PyAny>>) -> PyResult<SurrealMV> {
        let mut pv = Vec::with_capacity(p.len());
        for x in &p {
            pv.push(parse_surreal(x)?);
        }
        Ok(self.wrap(self.inner.up(&pv)))
    }
    /// Recover a Euclidean point from a null vector (`None` if not normalizable).
    fn down(&self, x: &SurrealMV) -> Option<Vec<PySurreal>> {
        self.inner
            .down(&x.mv)
            .map(|v| v.into_iter().map(wrap_surreal).collect())
    }
    /// The conformal inner product `x · y` (= `−½|p−q|²` on lifted points).
    fn inner(&self, x: &SurrealMV, y: &SurrealMV) -> PySurreal {
        wrap_surreal(self.inner.inner(&x.mv, &y.mv))
    }
    /// The sphere of squared radius `r2` about center `c`.
    fn sphere(&self, c: Vec<Bound<'_, PyAny>>, r2: &Bound<'_, PyAny>) -> PyResult<SurrealMV> {
        let mut cv = Vec::with_capacity(c.len());
        for x in &c {
            cv.push(parse_surreal(x)?);
        }
        Ok(self.wrap(self.inner.sphere(&cv, &parse_surreal(r2)?)))
    }
    /// The plane `{x : x·normal = d}`.
    fn plane(&self, normal: Vec<Bound<'_, PyAny>>, d: &Bound<'_, PyAny>) -> PyResult<SurrealMV> {
        let mut nv = Vec::with_capacity(normal.len());
        for x in &normal {
            nv.push(parse_surreal(x)?);
        }
        Ok(self.wrap(self.inner.plane(&nv, &parse_surreal(d)?)))
    }
    /// The point pair / oriented join `a ∧ b`.
    fn point_pair(&self, a: &SurrealMV, b: &SurrealMV) -> SurrealMV {
        self.wrap(self.inner.point_pair(&a.mv, &b.mv))
    }
    /// The IPNS meet (intersection) `x ∧ y`.
    fn meet(&self, x: &SurrealMV, y: &SurrealMV) -> SurrealMV {
        self.wrap(self.inner.meet(&x.mv, &y.mv))
    }
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<NimberAlgebra>()?;
    m.add_class::<NimberMV>()?;
    m.add_class::<SurrealAlgebra>()?;
    m.add_class::<SurrealMV>()?;
    m.add_class::<SurcomplexAlgebra>()?;
    m.add_class::<SurcomplexMV>()?;
    m.add_class::<IntegerAlgebra>()?;
    m.add_class::<IntegerMV>()?;
    m.add_class::<OmnificAlgebra>()?;
    m.add_class::<OmnificMV>()?;
    m.add_class::<PyCga>()?;
    Ok(())
}
