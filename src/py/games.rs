//! Python bindings for combinatorial game theory: partizan games, the exterior
//! algebra of the game group (over the `Integer` backend), and nim-mult via the
//! Turning-Corners game recurrence.

use super::engine::IntegerMV;
use crate::clifford::CliffordAlgebra;
use crate::games::{Game, GameExterior};
use crate::scalar::Integer;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::sync::Arc;

/// Nim-multiplication via Conway's Turning-Corners game recurrence (the
/// game-theoretic definition; equals the algebraic nim-product).
#[pyfunction]
fn nim_mul_mex(x: u128, y: u128) -> u128 {
    crate::games::nim_mul_mex(x, y)
}
// ---------------------------------------------------------------------------
// Partizan games + the exterior algebra of the game group
// ---------------------------------------------------------------------------

#[pyclass(name = "Game", module = "pleroma", from_py_object)]
#[derive(Clone)]
struct PyGame {
    inner: Game,
}

#[pymethods]
impl PyGame {
    #[staticmethod]
    fn zero() -> PyGame {
        PyGame {
            inner: Game::zero(),
        }
    }
    #[staticmethod]
    fn star() -> PyGame {
        PyGame {
            inner: Game::star(),
        }
    }
    #[staticmethod]
    fn up() -> PyGame {
        PyGame { inner: Game::up() }
    }
    #[staticmethod]
    fn integer(n: i128) -> PyGame {
        PyGame {
            inner: Game::integer(n),
        }
    }
    #[staticmethod]
    fn switch(a: i128, b: i128) -> PyGame {
        PyGame {
            inner: Game::switch(a, b),
        }
    }
    fn __add__(&self, other: &PyGame) -> PyGame {
        PyGame {
            inner: self.inner.add(&other.inner),
        }
    }
    fn __neg__(&self) -> PyGame {
        PyGame {
            inner: self.inner.neg(),
        }
    }
    fn __sub__(&self, other: &PyGame) -> PyGame {
        PyGame {
            inner: self.inner.add(&other.inner.neg()),
        }
    }
    fn le(&self, other: &PyGame) -> bool {
        self.inner.le(&other.inner)
    }
    fn __eq__(&self, other: &PyGame) -> bool {
        self.inner.eq(&other.inner)
    }
    fn fuzzy(&self, other: &PyGame) -> bool {
        self.inner.fuzzy(&other.inner)
    }
    fn birthday(&self) -> u128 {
        self.inner.birthday()
    }
    fn is_number(&self) -> bool {
        self.inner.is_number()
    }
    fn times_int(&self, n: i128) -> PyGame {
        PyGame {
            inner: self.inner.times_int(n),
        }
    }
    fn __repr__(&self) -> String {
        self.inner.display()
    }
}

#[pyclass(name = "GameExterior", module = "pleroma")]
struct PyGameExterior {
    inner: GameExterior,
    alg: Arc<CliffordAlgebra<Integer>>,
}

#[pymethods]
impl PyGameExterior {
    #[new]
    fn new(gens: Vec<PyGame>) -> Self {
        let games: Vec<Game> = gens.iter().map(|g| g.inner.clone()).collect();
        PyGameExterior::from_inner(GameExterior::new(games))
    }
    #[staticmethod]
    fn free(gens: Vec<PyGame>) -> Self {
        let games: Vec<Game> = gens.iter().map(|g| g.inner.clone()).collect();
        PyGameExterior::from_inner(GameExterior::free(games))
    }
    #[staticmethod]
    fn with_relation_bound(gens: Vec<PyGame>, bound: i128) -> Self {
        let games: Vec<Game> = gens.iter().map(|g| g.inner.clone()).collect();
        PyGameExterior::from_inner(GameExterior::with_relation_search(games, bound))
    }
    #[getter]
    fn dim(&self) -> usize {
        self.inner.algebra().dim
    }
    fn relations(&self) -> Vec<Vec<i128>> {
        self.inner
            .relations()
            .iter()
            .map(|r| r.coeffs.clone())
            .collect()
    }
    /// The grade-1 generator e_i (an `IntegerMV`) standing for game g_i.
    fn generator(&self, i: usize) -> IntegerMV {
        IntegerMV {
            alg: self.alg.clone(),
            mv: self.inner.generator(i),
        }
    }
    /// The game g_i a generator stands for.
    fn game(&self, i: usize) -> PyGame {
        PyGame {
            inner: self.inner.game(i).clone(),
        }
    }
    fn reduce(&self, mv: &IntegerMV) -> PyResult<IntegerMV> {
        self.ensure_mv(mv)?;
        Ok(IntegerMV {
            alg: self.alg.clone(),
            mv: self.inner.reduce(&mv.mv),
        })
    }
    fn add(&self, a: &IntegerMV, b: &IntegerMV) -> PyResult<IntegerMV> {
        self.ensure_mv(a)?;
        self.ensure_mv(b)?;
        Ok(IntegerMV {
            alg: self.alg.clone(),
            mv: self.inner.add(&a.mv, &b.mv),
        })
    }
    fn wedge(&self, a: &IntegerMV, b: &IntegerMV) -> PyResult<IntegerMV> {
        self.ensure_mv(a)?;
        self.ensure_mv(b)?;
        Ok(IntegerMV {
            alg: self.alg.clone(),
            mv: self.inner.wedge(&a.mv, &b.mv),
        })
    }
    fn is_zero(&self, mv: &IntegerMV) -> PyResult<bool> {
        self.ensure_mv(mv)?;
        Ok(self.inner.is_zero(&mv.mv))
    }
    /// Map a grade-1 element Σ c_i e_i back to the game Σ c_i·g_i (the module map
    /// Λ¹ → game group). Errors if the multivector is not purely grade 1.
    fn value_of_grade1(&self, mv: &IntegerMV) -> PyResult<PyGame> {
        self.ensure_mv(mv)?;
        let reduced = self.inner.reduce(&mv.mv);
        if reduced.terms.keys().any(|blade| blade.count_ones() != 1) {
            return Err(PyValueError::new_err("expected a grade-1 element"));
        }
        Ok(PyGame {
            inner: self.inner.value_of_grade1(&reduced),
        })
    }
}

impl PyGameExterior {
    fn from_inner(inner: GameExterior) -> Self {
        let alg = Arc::new(inner.algebra().clone());
        PyGameExterior { inner, alg }
    }

    fn ensure_mv(&self, mv: &IntegerMV) -> PyResult<()> {
        if self.alg.as_ref() == mv.alg.as_ref() {
            Ok(())
        } else {
            Err(PyValueError::new_err(
                "multivector belongs to a different GameExterior algebra",
            ))
        }
    }
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyGame>()?;
    m.add_class::<PyGameExterior>()?;
    m.add_function(wrap_pyfunction!(nim_mul_mex, m)?)?;
    Ok(())
}
