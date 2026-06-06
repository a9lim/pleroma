//! Python bindings for combinatorial game theory: partizan games, the exterior
//! algebra of the game group (over the `Integer` backend), and nim-mult via the
//! Turning-Corners game recurrence.

use super::engine::IntegerMV;
use crate::clifford::{CliffordAlgebra, Metric};
use crate::games::Game;
use crate::scalar::Integer;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::sync::Arc;

/// Nim-multiplication via Conway's Turning-Corners game recurrence (the
/// game-theoretic definition; equals the algebraic nim-product).
#[pyfunction]
fn nim_mul_mex(x: u64, y: u64) -> u64 {
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
    fn integer(n: i64) -> PyGame {
        PyGame {
            inner: Game::integer(n),
        }
    }
    #[staticmethod]
    fn switch(a: i64, b: i64) -> PyGame {
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
    fn birthday(&self) -> u32 {
        self.inner.birthday()
    }
    fn is_number(&self) -> bool {
        self.inner.is_number()
    }
    fn times_int(&self, n: i64) -> PyGame {
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
    alg: Arc<CliffordAlgebra<Integer>>,
    gens: Vec<Game>,
}

#[pymethods]
impl PyGameExterior {
    #[new]
    fn new(gens: Vec<PyGame>) -> Self {
        let games: Vec<Game> = gens.iter().map(|g| g.inner.clone()).collect();
        let n = games.len();
        PyGameExterior {
            alg: Arc::new(CliffordAlgebra::new(n, Metric::grassmann(n))),
            gens: games,
        }
    }
    #[getter]
    fn dim(&self) -> usize {
        self.gens.len()
    }
    /// The grade-1 generator e_i (an `IntegerMV`) standing for game g_i.
    fn generator(&self, i: usize) -> IntegerMV {
        IntegerMV {
            alg: self.alg.clone(),
            mv: self.alg.gen(i),
        }
    }
    /// The game g_i a generator stands for.
    fn game(&self, i: usize) -> PyGame {
        PyGame {
            inner: self.gens[i].clone(),
        }
    }
    /// Map a grade-1 element Σ c_i e_i back to the game Σ c_i·g_i (the module map
    /// Λ¹ → game group). Errors if the multivector is not purely grade 1.
    fn value_of_grade1(&self, mv: &IntegerMV) -> PyResult<PyGame> {
        let mut acc = Game::zero();
        for (&blade, coeff) in &mv.mv.terms {
            if blade.count_ones() != 1 {
                return Err(PyValueError::new_err("expected a grade-1 element"));
            }
            let idx = blade.trailing_zeros() as usize;
            acc = acc.add(&self.gens[idx].times_int(coeff.0));
        }
        Ok(PyGame { inner: acc })
    }
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyGame>()?;
    m.add_class::<PyGameExterior>()?;
    m.add_function(wrap_pyfunction!(nim_mul_mex, m)?)?;
    Ok(())
}
