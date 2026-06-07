//! Python bindings for combinatorial game theory: partizan games, the exterior
//! algebra of the game group (over the `Integer` backend), and nim-mult via the
//! Turning-Corners game recurrence.

use super::engine::IntegerMV;
use super::scalars::{parse_surreal, PySurreal};
use crate::clifford::CliffordAlgebra;
use crate::games::{thermography, Color, Game, GameExterior, Hackenbush, NumberGame};
use crate::scalar::{Integer, Rational, Surreal};
use pyo3::basic::CompareOp;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::sync::Arc;

/// Wrap a dyadic `Rational` (a thermograph coordinate) as a `Surreal` for Python.
fn rat_to_py(r: Rational) -> PySurreal {
    PySurreal::from_inner(Surreal::from_rational(r))
}

/// Nim-multiplication via Conway's Turning-Corners game recurrence (the
/// game-theoretic definition; equals the algebraic nim-product).
#[pyfunction]
fn nim_mul_mex(x: u128, y: u128) -> u128 {
    crate::games::nim_mul_mex(x, y)
}

/// Sprague–Grundy values of a finite **acyclic** impartial game graph given as
/// adjacency lists (`succ[v]` = positions reachable from `v`). Errors on a cycle
/// (Grundy values are undefined on loopy games). A position is a P-position iff
/// its value is 0.
#[pyfunction]
fn grundy_graph(succ: Vec<Vec<usize>>) -> PyResult<Vec<u128>> {
    crate::games::grundy_graph(&succ)
        .ok_or_else(|| PyValueError::new_err("graph has a cycle — Grundy value is undefined"))
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
    /// A general game `{ left | right }` from explicit option lists.
    #[staticmethod]
    fn of(left: Vec<PyGame>, right: Vec<PyGame>) -> PyGame {
        PyGame {
            inner: Game::new(
                left.into_iter().map(|g| g.inner).collect(),
                right.into_iter().map(|g| g.inner).collect(),
            ),
        }
    }
    /// The Left options.
    fn left(&self) -> Vec<PyGame> {
        self.inner
            .left()
            .iter()
            .map(|g| PyGame { inner: g.clone() })
            .collect()
    }
    /// The Right options.
    fn right(&self) -> Vec<PyGame> {
        self.inner
            .right()
            .iter()
            .map(|g| PyGame { inner: g.clone() })
            .collect()
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
    /// The Nim-heap `⋆n` (the remote/far star of the atomic-weight calculus).
    #[staticmethod]
    fn star_n(n: u128) -> PyGame {
        PyGame {
            inner: Game::nim_heap(n),
        }
    }
    /// Whether this game is **all-small** (a Left move iff a Right move at every
    /// position) — the domain of the atomic weight.
    fn is_all_small(&self) -> bool {
        self.inner.is_all_small()
    }
    /// The **atomic weight** as a `Game` (`None` if not all-small).
    fn atomic_weight(&self) -> Option<PyGame> {
        crate::games::atomic_weight(&self.inner).map(|inner| PyGame { inner })
    }
    /// The **atomic weight as an integer** (`None` if not all-small or its atomic
    /// weight is a genuine non-integer game). `aw(↑)=1`, `aw(⋆)=0`, `aw(⇑)=2`.
    fn atomic_weight_int(&self) -> Option<i128> {
        crate::games::atomic_weight_int(&self.inner)
    }
    fn times_int(&self, n: i128) -> PyGame {
        PyGame {
            inner: self.inner.times_int(n),
        }
    }
    /// The canonical form: the unique simplest game equal in value (dominated
    /// options removed, reversible options bypassed).
    fn canonical(&self) -> PyGame {
        PyGame {
            inner: self.inner.canonical(),
        }
    }
    /// Whether this game is already in canonical form.
    fn is_canonical(&self) -> bool {
        self.inner.is_canonical()
    }
    /// An order-independent canonical string `{L|R}` — equal iff the games are
    /// equal in value.
    fn canonical_string(&self) -> String {
        self.inner.canonical_string()
    }
    /// The surreal value of a number-valued game (`None` for non-numbers like
    /// `⋆`, `↑`, switches).
    fn number_value(&self) -> Option<PySurreal> {
        self.inner.number_value().map(PySurreal::from_inner)
    }
    /// The canonical game of a dyadic surreal (or int); errors for non-dyadics.
    #[staticmethod]
    fn from_surreal(s: &Bound<'_, PyAny>) -> PyResult<PyGame> {
        let s = parse_surreal(s)?;
        Game::from_surreal(&s)
            .map(|inner| PyGame { inner })
            .ok_or_else(|| PyValueError::new_err("surreal is not a dyadic rational"))
    }
    /// The ordinal sum `G : H` (play in `H`; a move in the base `G` discards `H`).
    fn ordinal_sum(&self, h: &PyGame) -> PyGame {
        PyGame {
            inner: self.inner.ordinal_sum(&h.inner),
        }
    }
    /// Temperature `t(G)` as a surreal (`−1` for a number); `None` for the rare
    /// degenerate positions outside temperature theory.
    fn temperature(&self) -> Option<PySurreal> {
        thermography::temperature(&self.inner).map(rat_to_py)
    }
    /// Mean (mast) value as a surreal.
    fn mean_value(&self) -> Option<PySurreal> {
        thermography::mean_value(&self.inner).map(rat_to_py)
    }
    /// Left stop `LS(G)` (left wall at temperature 0).
    fn left_stop(&self) -> Option<PySurreal> {
        thermography::left_stop(&self.inner).map(rat_to_py)
    }
    /// Right stop `RS(G)` (right wall at temperature 0).
    fn right_stop(&self) -> Option<PySurreal> {
        thermography::right_stop(&self.inner).map(rat_to_py)
    }
    /// The thermograph as `(mean, temperature, left_wall, right_wall)`, where each
    /// wall is a list of `(t, value)` breakpoints. `None` if undefined.
    #[allow(clippy::type_complexity)]
    fn thermograph(
        &self,
    ) -> Option<(
        PySurreal,
        PySurreal,
        Vec<(PySurreal, PySurreal)>,
        Vec<(PySurreal, PySurreal)>,
    )> {
        let th = thermography::thermograph(&self.inner)?;
        let wall = |w: &thermography::Pl| {
            w.points()
                .iter()
                .map(|(t, v)| (rat_to_py(t.clone()), rat_to_py(v.clone())))
                .collect::<Vec<_>>()
        };
        Some((
            rat_to_py(th.mast.clone()),
            rat_to_py(th.temperature.clone()),
            wall(&th.left_wall),
            wall(&th.right_wall),
        ))
    }
    fn __repr__(&self) -> String {
        self.inner.display()
    }
}

/// Parse a colour name (`"blue"`/`"red"`/`"green"`, case-insensitive) or its
/// initial (`"b"`/`"r"`/`"g"`).
fn parse_color(s: &str) -> PyResult<Color> {
    match s.trim().to_lowercase().as_str() {
        "blue" | "b" | "l" | "left" => Ok(Color::Blue),
        "red" | "r" => Ok(Color::Red),
        "green" | "g" | "e" => Ok(Color::Green),
        other => Err(PyValueError::new_err(format!(
            "unknown colour {other:?} (expected blue/red/green)"
        ))),
    }
}

#[pyclass(name = "Hackenbush", module = "pleroma")]
struct PyHackenbush {
    inner: Hackenbush,
}

#[pymethods]
impl PyHackenbush {
    /// A position from `(u, v, colour)` edges; vertex `0` is the ground.
    #[new]
    fn new(edges: Vec<(usize, usize, String)>) -> PyResult<Self> {
        let edges = edges
            .into_iter()
            .map(|(u, v, c)| Ok((u, v, parse_color(&c)?)))
            .collect::<PyResult<Vec<_>>>()?;
        Ok(PyHackenbush {
            inner: Hackenbush::new(edges),
        })
    }
    /// A stalk `0—1—2—…` from the ground, edge `i` coloured `colors[i]`.
    #[staticmethod]
    fn string(colors: Vec<String>) -> PyResult<Self> {
        let cs = colors
            .iter()
            .map(|c| parse_color(c))
            .collect::<PyResult<Vec<_>>>()?;
        Ok(PyHackenbush {
            inner: Hackenbush::string(&cs),
        })
    }
    /// The partizan game value (the universal evaluator).
    fn to_game(&self) -> PyGame {
        PyGame {
            inner: self.inner.to_game(),
        }
    }
    /// The surreal number value (`None` if the value is not a number).
    fn value(&self) -> Option<PySurreal> {
        self.inner.value().map(PySurreal::from_inner)
    }
    /// The Sprague–Grundy / nim value (`Some` only for all-green positions).
    fn grundy(&self) -> Option<u128> {
        self.inner.grundy()
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

/// A transfinite **number-valued** game, carried by its surreal value (e.g. the
/// game `ω = {0,1,2,…|}`). The numbers-only honoring of transfinite birthdays —
/// value/birthday/sum/order all delegate to the surreal.
#[pyclass(name = "NumberGame", module = "pleroma", from_py_object)]
#[derive(Clone)]
struct PyNumberGame {
    inner: NumberGame,
}

#[pymethods]
impl PyNumberGame {
    #[staticmethod]
    fn from_surreal(s: &Bound<'_, PyAny>) -> PyResult<PyNumberGame> {
        Ok(PyNumberGame {
            inner: NumberGame::from_surreal(&parse_surreal(s)?),
        })
    }
    /// The exact surreal value.
    fn value(&self) -> PySurreal {
        PySurreal::from_inner(self.inner.value().clone())
    }
    /// The birthday as a finite ordinal value, if finite.
    fn birthday_finite(&self) -> Option<u128> {
        self.inner.birthday().and_then(|o| o.as_finite())
    }
    /// The birthday as an ordinal string (`None` outside the representable
    /// subclass, e.g. `√ω`).
    fn birthday_repr(&self) -> Option<String> {
        self.inner.birthday().map(|o| format!("{o:?}"))
    }
    /// The short `Game`, if the value is dyadic; `None` for transfinite numbers.
    fn to_finite_game(&self) -> Option<PyGame> {
        self.inner.to_finite_game().map(|inner| PyGame { inner })
    }
    fn __add__(&self, other: &PyNumberGame) -> PyNumberGame {
        PyNumberGame {
            inner: self.inner.add(&other.inner),
        }
    }
    fn __neg__(&self) -> PyNumberGame {
        PyNumberGame {
            inner: self.inner.neg(),
        }
    }
    fn __richcmp__(&self, other: &PyNumberGame, op: CompareOp) -> bool {
        op.matches(self.inner.cmp(&other.inner))
    }
    fn __repr__(&self) -> String {
        format!("NumberGame({:?})", self.inner.value())
    }
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyGame>()?;
    m.add_class::<PyNumberGame>()?;
    m.add_class::<PyGameExterior>()?;
    m.add_class::<PyHackenbush>()?;
    m.add_function(wrap_pyfunction!(nim_mul_mex, m)?)?;
    m.add_function(wrap_pyfunction!(grundy_graph, m)?)?;
    Ok(())
}
