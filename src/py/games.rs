//! Python bindings for combinatorial game theory: partizan games, the exterior
//! algebra of the game group (over the `Integer` backend), and nim-mult via the
//! Turning-Corners game recurrence.

use super::engine::{IntegerAlgebra, IntegerMV};
use super::forms::{wrap_quadric_fit, PyQuadricFit};
use super::scalars::{
    parse_rational, parse_surreal, wrap_rational, PyOrdinal, PyRational, PySurreal,
};
use crate::clifford::CliffordAlgebra;
use crate::games::{
    thermography, AbstractGame, Color, Game, GameExterior, GameRelation, Hackenbush, LoopyGraph,
    LoopyNimCertificate, LoopyNimber, LoopyValue, NimberGame, NumberGame, Outcome, PartizanOutcome,
    Quotient,
};
use crate::scalar::{Integer, Rational, SignExpansion, Surreal};
use pyo3::basic::CompareOp;
use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;
use std::collections::{HashMap, HashSet};
use std::sync::Arc;

/// Validate an adjacency list: every successor index must be `< n = succ.len()`.
/// Returns `PyValueError` with the offending position and index on failure.
/// This is the shared bounds-check for all list-input Python entry points that
/// forward directly to Rust kernels (which would otherwise panic on out-of-range
/// indices in the predecessor/successor arrays).
fn check_succ_bounds(succ: &[Vec<usize>]) -> PyResult<()> {
    let n = succ.len();
    for (v, neighbors) in succ.iter().enumerate() {
        for &w in neighbors {
            if w >= n {
                return Err(PyValueError::new_err(format!(
                    "adjacency list out of range: succ[{v}] contains index {w}, \
                     but the graph has only {n} positions (0..{n})"
                )));
            }
        }
    }
    Ok(())
}

/// Wrap a dyadic `Rational` (a thermograph coordinate) as a `Surreal` for Python.
fn rat_to_py(r: Rational) -> PySurreal {
    PySurreal::from_inner(Surreal::from_rational(r))
}

fn wrap_pl(inner: thermography::Pl) -> PyPl {
    PyPl { inner }
}

fn wrap_thermograph(inner: thermography::Thermograph) -> PyThermograph {
    PyThermograph { inner }
}

#[pyclass(name = "GameRelation", module = "ogdoad", from_py_object)]
#[derive(Clone)]
struct PyGameRelation {
    inner: GameRelation,
}

#[pymethods]
impl PyGameRelation {
    #[new]
    fn new(coeffs: Vec<i128>) -> Self {
        PyGameRelation {
            inner: GameRelation::new(coeffs),
        }
    }
    #[getter]
    fn coeffs(&self) -> Vec<i128> {
        self.inner.coeffs.clone()
    }
    fn __richcmp__(&self, other: &Bound<'_, PyAny>, op: CompareOp) -> PyResult<bool> {
        let Ok(rel) = other.extract::<PyRef<'_, PyGameRelation>>() else {
            return Ok(matches!(op, CompareOp::Ne));
        };
        match op {
            CompareOp::Eq => Ok(self.inner.coeffs == rel.inner.coeffs),
            CompareOp::Ne => Ok(self.inner.coeffs != rel.inner.coeffs),
            CompareOp::Lt | CompareOp::Le | CompareOp::Gt | CompareOp::Ge => Err(
                PyValueError::new_err("GameRelation only supports equality comparisons"),
            ),
        }
    }
    fn __repr__(&self) -> String {
        format!("GameRelation(coeffs={:?})", self.inner.coeffs)
    }
}

fn wrap_game_relation(inner: GameRelation) -> PyGameRelation {
    PyGameRelation { inner }
}

#[pyclass(name = "GameRelationCertificate", module = "ogdoad")]
struct PyGameRelationCertificate {
    inner: crate::games::GameRelationCertificate,
}

#[pymethods]
impl PyGameRelationCertificate {
    #[getter]
    fn coeffs(&self) -> Vec<i128> {
        self.inner.coeffs.clone()
    }
    #[getter]
    fn value_key(&self) -> String {
        self.inner.value_key.clone()
    }
    #[getter]
    fn independent(&self) -> bool {
        self.inner.independent
    }
    fn __repr__(&self) -> String {
        format!(
            "GameRelationCertificate(coeffs={:?}, value_key={:?}, independent={})",
            self.inner.coeffs, self.inner.value_key, self.inner.independent
        )
    }
}

fn wrap_game_relation_certificate(
    inner: crate::games::GameRelationCertificate,
) -> PyGameRelationCertificate {
    PyGameRelationCertificate { inner }
}

#[pyclass(name = "RelationSearchCertificate", module = "ogdoad")]
struct PyRelationSearchCertificate {
    inner: crate::games::RelationSearchCertificate,
}

#[pymethods]
impl PyRelationSearchCertificate {
    #[getter]
    fn bound(&self) -> i128 {
        self.inner.bound
    }
    #[getter]
    fn exhaustive(&self) -> bool {
        self.inner.exhaustive
    }
    #[getter]
    fn candidate_count(&self) -> Option<usize> {
        self.inner.candidate_count
    }
    #[getter]
    fn relations(&self) -> Vec<PyGameRelationCertificate> {
        self.inner
            .relations
            .iter()
            .cloned()
            .map(wrap_game_relation_certificate)
            .collect()
    }
    fn __repr__(&self) -> String {
        let relations: Vec<_> = self
            .inner
            .relations
            .iter()
            .map(|r| {
                format!(
                    "GameRelationCertificate(coeffs={:?}, value_key={:?}, independent={})",
                    r.coeffs, r.value_key, r.independent
                )
            })
            .collect();
        format!(
            "RelationSearchCertificate(bound={}, exhaustive={}, candidate_count={:?}, relations={:?})",
            self.inner.bound,
            self.inner.exhaustive,
            self.inner.candidate_count,
            relations,
        )
    }
}

fn wrap_relation_search_certificate(
    inner: crate::games::RelationSearchCertificate,
) -> PyRelationSearchCertificate {
    PyRelationSearchCertificate { inner }
}

#[pyclass(name = "Pl", module = "ogdoad", from_py_object)]
#[derive(Clone)]
struct PyPl {
    inner: thermography::Pl,
}

#[pymethods]
impl PyPl {
    #[staticmethod]
    fn constant(value: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(wrap_pl(thermography::Pl::constant(parse_rational(value)?)))
    }
    fn points(&self) -> Vec<(PyRational, PyRational)> {
        self.inner
            .points()
            .iter()
            .map(|(t, v)| (wrap_rational(t.clone()), wrap_rational(v.clone())))
            .collect()
    }
    fn value_at(&self, t: &Bound<'_, PyAny>) -> PyResult<PyRational> {
        Ok(wrap_rational(self.inner.value_at(&parse_rational(t)?)))
    }
    fn oplus_max(&self, other: &PyPl) -> PyPl {
        wrap_pl(self.inner.oplus_max(&other.inner))
    }
    fn oplus_min(&self, other: &PyPl) -> PyPl {
        wrap_pl(self.inner.oplus_min(&other.inner))
    }
    fn otimes(&self, other: &PyPl) -> PyPl {
        wrap_pl(self.inner.otimes(&other.inner))
    }
    fn __repr__(&self) -> String {
        format!("Pl({:?})", self.inner.points())
    }
}

#[pyclass(name = "Thermograph", module = "ogdoad", from_py_object)]
#[derive(Clone)]
struct PyThermograph {
    inner: thermography::Thermograph,
}

#[pymethods]
impl PyThermograph {
    #[getter]
    fn mean(&self) -> PyRational {
        wrap_rational(self.inner.mean())
    }
    #[getter]
    fn temperature(&self) -> PyRational {
        wrap_rational(self.inner.temperature.clone())
    }
    #[getter]
    fn left_wall(&self) -> PyPl {
        wrap_pl(self.inner.left_wall.clone())
    }
    #[getter]
    fn right_wall(&self) -> PyPl {
        wrap_pl(self.inner.right_wall.clone())
    }
    fn left_stop(&self) -> PyRational {
        wrap_rational(self.inner.left_stop())
    }
    fn right_stop(&self) -> PyRational {
        wrap_rational(self.inner.right_stop())
    }
    fn cooled_stops(&self, t: &Bound<'_, PyAny>) -> PyResult<(PyRational, PyRational)> {
        let (l, r) = self.inner.cooled_stops(&parse_rational(t)?);
        Ok((wrap_rational(l), wrap_rational(r)))
    }
    fn __repr__(&self) -> String {
        format!(
            "Thermograph(mean={:?}, temperature={:?})",
            self.inner.mast, self.inner.temperature
        )
    }
}

/// Nim-multiplication via Conway's Turning-Corners game recurrence (the
/// game-theoretic definition; equals the algebraic nim-product).
#[pyfunction]
fn nim_mul_mex(x: u128, y: u128) -> u128 {
    crate::games::nim_mul_mex(x, y)
}

#[derive(Clone, Copy)]
enum CoinFamily {
    Singleton,
    Turtles,
}

type CompanionFn = fn(u128) -> Vec<u128>;

impl CoinFamily {
    fn companions(self) -> CompanionFn {
        match self {
            CoinFamily::Singleton => crate::games::singleton_companions,
            CoinFamily::Turtles => crate::games::turtles_companions,
        }
    }
}

fn parse_coin_family(name: &str) -> PyResult<CoinFamily> {
    match name.trim().to_ascii_lowercase().as_str() {
        "singleton" | "singletons" | "turning-corners" | "turning_corners" | "corners" => {
            Ok(CoinFamily::Singleton)
        }
        "turtles" | "turning-turtles" | "turning_turtles" => Ok(CoinFamily::Turtles),
        other => Err(PyValueError::new_err(format!(
            "unknown coin-turning family {other:?}; expected 'singleton' or 'turtles'"
        ))),
    }
}

fn check_coin_index(n: u128, label: &str) -> PyResult<()> {
    if n >= 128 {
        Err(PyValueError::new_err(format!(
            "{label} must be < 128 because companion sets are u128 bitmasks"
        )))
    } else {
        Ok(())
    }
}

fn lower_coin_mask(n: u128, label: &str) -> PyResult<u128> {
    check_coin_index(n, label)?;
    Ok(if n == 0 { 0 } else { (1u128 << n) - 1 })
}

fn call_u128_moves(callback: &Bound<'_, PyAny>, pos: u128) -> PyResult<Vec<u128>> {
    callback.call1((pos,))?.extract()
}

fn call_usize_moves(callback: &Bound<'_, PyAny>, pos: usize) -> PyResult<Vec<usize>> {
    callback.call1((pos,))?.extract()
}

fn grundy_u128_inner(
    pos: u128,
    moves: &Bound<'_, PyAny>,
    memo: &mut HashMap<u128, u128>,
) -> PyResult<u128> {
    if let Some(&v) = memo.get(&pos) {
        return Ok(v);
    }
    let nexts = call_u128_moves(moves, pos)?;
    let mut values = Vec::with_capacity(nexts.len());
    for next in nexts {
        values.push(grundy_u128_inner(next, moves, memo)?);
    }
    let g = crate::games::mex(values);
    memo.insert(pos, g);
    Ok(g)
}

fn misere_is_n_u128_inner(
    pos: u128,
    moves: &Bound<'_, PyAny>,
    memo: &mut HashMap<u128, bool>,
    visiting: &mut HashSet<u128>,
) -> PyResult<Option<bool>> {
    if let Some(&v) = memo.get(&pos) {
        return Ok(Some(v));
    }
    if !visiting.insert(pos) {
        return Ok(None);
    }
    let nexts = call_u128_moves(moves, pos)?;
    let mut result = nexts.is_empty();
    if !result {
        for next in nexts {
            match misere_is_n_u128_inner(next, moves, memo, visiting)? {
                Some(false) => {
                    result = true;
                    break;
                }
                Some(true) => {}
                None => {
                    visiting.remove(&pos);
                    return Ok(None);
                }
            }
        }
    }
    visiting.remove(&pos);
    memo.insert(pos, result);
    Ok(Some(result))
}

fn loopy_succ_from_callback(n: usize, moves: &Bound<'_, PyAny>) -> PyResult<Vec<Vec<usize>>> {
    let mut succ = Vec::with_capacity(n);
    for v in 0..n {
        let nexts = call_usize_moves(moves, v)?;
        if nexts.iter().any(|&w| w >= n) {
            return Err(PyValueError::new_err(format!(
                "move callback for position {v} returned an out-of-range target"
            )));
        }
        succ.push(nexts);
    }
    Ok(succ)
}

/// Legal companion masks for the named 1-D coin-turning family at coin `n`.
/// Families: `"singleton"` (turn exactly one lower coin) and `"turtles"`
/// (turn optionally one lower coin).
#[pyfunction]
fn coin_companions(kind: &str, n: u128) -> PyResult<Vec<u128>> {
    check_coin_index(n, "n")?;
    Ok((parse_coin_family(kind)?.companions())(n))
}

#[pyfunction]
fn singleton_companions(n: u128) -> PyResult<Vec<u128>> {
    check_coin_index(n, "n")?;
    Ok(crate::games::singleton_companions(n))
}

#[pyfunction]
fn turtles_companions(n: u128) -> PyResult<Vec<u128>> {
    check_coin_index(n, "n")?;
    Ok(crate::games::turtles_companions(n))
}

/// Single-coin Grundy value for a Python companion-mask callback.
/// `companions(n)` must return legal bitmasks over lower coins `{0, …, n-1}`.
#[pyfunction]
fn grundy_1d(companions: Bound<'_, PyAny>, n: u128) -> PyResult<u128> {
    fn inner(
        companions: &Bound<'_, PyAny>,
        n: u128,
        memo: &mut HashMap<u128, u128>,
    ) -> PyResult<u128> {
        let allowed = lower_coin_mask(n, "n")?;
        if let Some(&v) = memo.get(&n) {
            return Ok(v);
        }
        let mut seen = HashSet::new();
        for s in call_u128_moves(companions, n)? {
            if s & !allowed != 0 {
                return Err(PyValueError::new_err(
                    "companion mask contains a coin that is not strictly lower than n",
                ));
            }
            let mut acc = 0u128;
            let mut ss = s;
            while ss != 0 {
                let i = ss.trailing_zeros() as u128;
                ss &= ss - 1;
                acc ^= inner(companions, i, memo)?;
            }
            seen.insert(acc);
        }
        let g = crate::games::mex(seen);
        memo.insert(n, g);
        Ok(g)
    }

    check_coin_index(n, "n")?;
    inner(&companions, n, &mut HashMap::new())
}

/// Single-coin Grundy value of a named 1-D coin-turning family.
#[pyfunction]
fn coin_turning_grundy(kind: &str, n: u128) -> PyResult<u128> {
    check_coin_index(n, "n")?;
    let mut memo = HashMap::new();
    let companions = parse_coin_family(kind)?.companions();
    Ok(crate::games::grundy_1d(&companions, n, &mut memo))
}

/// Single-cell Grundy value of the Tartan product of two named 1-D coin-turning
/// families, computed directly from the 2-D excludant.
#[pyfunction]
fn coin_turning_tartan_grundy(kind_a: &str, kind_b: &str, x: u128, y: u128) -> PyResult<u128> {
    check_coin_index(x, "x")?;
    check_coin_index(y, "y")?;
    let comp_a = parse_coin_family(kind_a)?.companions();
    let comp_b = parse_coin_family(kind_b)?.companions();
    let mut memo = HashMap::new();
    Ok(crate::games::tartan_grundy(
        &comp_a, &comp_b, x, y, &mut memo,
    ))
}

/// Single-cell Grundy value of the Tartan product for two Python companion-mask
/// callbacks.
#[pyfunction]
fn tartan_grundy(
    comp_a: Bound<'_, PyAny>,
    comp_b: Bound<'_, PyAny>,
    x: u128,
    y: u128,
) -> PyResult<u128> {
    fn inner(
        comp_a: &Bound<'_, PyAny>,
        comp_b: &Bound<'_, PyAny>,
        x: u128,
        y: u128,
        memo: &mut HashMap<(u128, u128), u128>,
    ) -> PyResult<u128> {
        let allowed_x = lower_coin_mask(x, "x")?;
        let allowed_y = lower_coin_mask(y, "y")?;
        if let Some(&v) = memo.get(&(x, y)) {
            return Ok(v);
        }
        let mut seen = HashSet::new();
        for ta in call_u128_moves(comp_a, x)? {
            if ta & !allowed_x != 0 {
                return Err(PyValueError::new_err(
                    "row companion mask contains a coin that is not strictly lower than x",
                ));
            }
            let acells = ta | (1u128 << x);
            for tb in call_u128_moves(comp_b, y)? {
                if tb & !allowed_y != 0 {
                    return Err(PyValueError::new_err(
                        "column companion mask contains a coin that is not strictly lower than y",
                    ));
                }
                let bcells = tb | (1u128 << y);
                let mut acc = 0u128;
                let mut aa = acells;
                while aa != 0 {
                    let a = aa.trailing_zeros() as u128;
                    aa &= aa - 1;
                    let mut bb = bcells;
                    while bb != 0 {
                        let b = bb.trailing_zeros() as u128;
                        bb &= bb - 1;
                        if a != x || b != y {
                            acc ^= inner(comp_a, comp_b, a, b, memo)?;
                        }
                    }
                }
                seen.insert(acc);
            }
        }
        let g = crate::games::mex(seen);
        memo.insert((x, y), g);
        Ok(g)
    }

    check_coin_index(x, "x")?;
    check_coin_index(y, "y")?;
    inner(&comp_a, &comp_b, x, y, &mut HashMap::new())
}

/// Grundy value of an acyclic impartial game with `u128` positions and a Python
/// move callback `moves(pos) -> list[u128]`.
#[pyfunction]
fn grundy(pos: u128, moves: Bound<'_, PyAny>) -> PyResult<u128> {
    grundy_u128_inner(pos, &moves, &mut HashMap::new())
}

/// Sprague–Grundy values of a finite **acyclic** impartial game graph given as
/// adjacency lists (`succ[v]` = positions reachable from `v`). Errors on a cycle
/// (Grundy values are undefined on loopy games). A position is a P-position iff
/// its value is 0.
#[pyfunction]
fn grundy_graph(succ: Vec<Vec<usize>>) -> PyResult<Vec<u128>> {
    check_succ_bounds(&succ)?;
    crate::games::grundy_graph(&succ)
        .ok_or_else(|| PyValueError::new_err("graph has a cycle — Grundy value is undefined"))
}

/// The minimal excludant `mex(S)` — the least non-negative integer not in the
/// multiset `values`. The core of the Sprague–Grundy recurrence.
#[pyfunction]
fn mex(values: Vec<u128>) -> u128 {
    crate::games::mex(values)
}

// ---------------------------------------------------------------------------
// Outcomes of a finite game graph (kernel): Win / Loss / Draw + scoring
// ---------------------------------------------------------------------------

fn outcome_name(o: Outcome) -> String {
    match o {
        Outcome::Loss => "Loss",
        Outcome::Win => "Win",
        Outcome::Draw => "Draw",
    }
    .to_string()
}

#[pyclass(name = "Outcome", module = "ogdoad", from_py_object)]
#[derive(Clone)]
struct PyOutcome {
    inner: Outcome,
}

fn wrap_outcome(inner: Outcome) -> PyOutcome {
    PyOutcome { inner }
}

fn parse_outcome(obj: &Bound<'_, PyAny>) -> PyResult<Outcome> {
    if let Ok(outcome) = obj.cast::<PyOutcome>() {
        return Ok(outcome.borrow().inner);
    }
    Err(PyTypeError::new_err("expected Outcome"))
}

#[pymethods]
impl PyOutcome {
    #[staticmethod]
    fn loss() -> Self {
        wrap_outcome(Outcome::Loss)
    }
    #[staticmethod]
    fn win() -> Self {
        wrap_outcome(Outcome::Win)
    }
    #[staticmethod]
    fn draw() -> Self {
        wrap_outcome(Outcome::Draw)
    }
    fn name(&self) -> String {
        outcome_name(self.inner)
    }
    fn is_loss(&self) -> bool {
        self.inner == Outcome::Loss
    }
    fn is_win(&self) -> bool {
        self.inner == Outcome::Win
    }
    fn is_draw(&self) -> bool {
        self.inner == Outcome::Draw
    }
    fn __richcmp__(&self, other: &Bound<'_, PyAny>, op: CompareOp) -> PyResult<bool> {
        match op {
            CompareOp::Eq => Ok(parse_outcome(other).is_ok_and(|o| o == self.inner)),
            CompareOp::Ne => Ok(parse_outcome(other).is_ok_and(|o| o != self.inner)),
            CompareOp::Lt | CompareOp::Le | CompareOp::Gt | CompareOp::Ge => Err(
                PyValueError::new_err("Outcome only supports equality comparisons"),
            ),
        }
    }
    fn __str__(&self) -> String {
        self.name()
    }
    fn __repr__(&self) -> String {
        format!("Outcome.{}", outcome_name(self.inner))
    }
}

fn partizan_outcome_name(o: PartizanOutcome) -> String {
    match o {
        PartizanOutcome::P => "P",
        PartizanOutcome::N => "N",
        PartizanOutcome::L => "L",
        PartizanOutcome::R => "R",
        PartizanOutcome::Draw => "Draw",
    }
    .to_string()
}

#[pyclass(name = "PartizanOutcome", module = "ogdoad", from_py_object)]
#[derive(Clone)]
struct PyPartizanOutcome {
    inner: PartizanOutcome,
}

fn wrap_partizan_outcome(inner: PartizanOutcome) -> PyPartizanOutcome {
    PyPartizanOutcome { inner }
}

fn parse_partizan_outcome(obj: &Bound<'_, PyAny>) -> PyResult<PartizanOutcome> {
    if let Ok(outcome) = obj.cast::<PyPartizanOutcome>() {
        return Ok(outcome.borrow().inner);
    }
    Err(PyTypeError::new_err("expected PartizanOutcome"))
}

#[pymethods]
impl PyPartizanOutcome {
    #[staticmethod]
    fn p() -> Self {
        wrap_partizan_outcome(PartizanOutcome::P)
    }
    #[staticmethod]
    fn n() -> Self {
        wrap_partizan_outcome(PartizanOutcome::N)
    }
    #[staticmethod]
    fn l() -> Self {
        wrap_partizan_outcome(PartizanOutcome::L)
    }
    #[staticmethod]
    fn r() -> Self {
        wrap_partizan_outcome(PartizanOutcome::R)
    }
    #[staticmethod]
    fn draw() -> Self {
        wrap_partizan_outcome(PartizanOutcome::Draw)
    }
    fn name(&self) -> String {
        partizan_outcome_name(self.inner)
    }
    fn __richcmp__(&self, other: &Bound<'_, PyAny>, op: CompareOp) -> PyResult<bool> {
        match op {
            CompareOp::Eq => Ok(parse_partizan_outcome(other).is_ok_and(|o| o == self.inner)),
            CompareOp::Ne => Ok(parse_partizan_outcome(other).is_ok_and(|o| o != self.inner)),
            CompareOp::Lt | CompareOp::Le | CompareOp::Gt | CompareOp::Ge => Err(
                PyValueError::new_err("PartizanOutcome only supports equality comparisons"),
            ),
        }
    }
    fn __str__(&self) -> String {
        self.name()
    }
    fn __repr__(&self) -> String {
        format!("PartizanOutcome.{}", partizan_outcome_name(self.inner))
    }
}

#[pyclass(name = "LoopyNimber", module = "ogdoad", from_py_object)]
#[derive(Clone)]
struct PyLoopyNimber {
    inner: LoopyNimber,
}

fn wrap_loopy_nimber(inner: LoopyNimber) -> PyLoopyNimber {
    PyLoopyNimber { inner }
}

fn parse_loopy_nimber(obj: &Bound<'_, PyAny>) -> PyResult<LoopyNimber> {
    if let Ok(value) = obj.cast::<PyLoopyNimber>() {
        return Ok(value.borrow().inner);
    }
    Err(PyTypeError::new_err("expected LoopyNimber"))
}

#[pymethods]
impl PyLoopyNimber {
    #[staticmethod]
    fn value(n: u128) -> Self {
        wrap_loopy_nimber(LoopyNimber::Value(n))
    }
    #[staticmethod]
    fn side() -> Self {
        wrap_loopy_nimber(LoopyNimber::Side)
    }
    fn to_u128(&self) -> PyResult<u128> {
        match self.inner {
            LoopyNimber::Value(n) => Ok(n),
            LoopyNimber::Side => Err(PyValueError::new_err("LoopyNimber.Side has no u128 value")),
        }
    }
    fn is_side(&self) -> bool {
        self.inner == LoopyNimber::Side
    }
    fn is_value(&self) -> bool {
        matches!(self.inner, LoopyNimber::Value(_))
    }
    fn __richcmp__(&self, other: &Bound<'_, PyAny>, op: CompareOp) -> PyResult<bool> {
        match op {
            CompareOp::Eq => Ok(parse_loopy_nimber(other).is_ok_and(|x| x == self.inner)),
            CompareOp::Ne => Ok(parse_loopy_nimber(other).is_ok_and(|x| x != self.inner)),
            CompareOp::Lt | CompareOp::Le | CompareOp::Gt | CompareOp::Ge => Err(
                PyValueError::new_err("LoopyNimber only supports equality comparisons"),
            ),
        }
    }
    fn __repr__(&self) -> String {
        match self.inner {
            LoopyNimber::Value(n) => format!("LoopyNimber.Value({n})"),
            LoopyNimber::Side => "LoopyNimber.Side".to_string(),
        }
    }
}

/// Normal-play typed [`Outcome`] of every position of a finite
/// game graph given as adjacency lists (`succ[v]` = positions reachable from `v`).
/// Retrograde analysis; `Loss` = P-position. Cyclic graphs are fine (→ `Draw`).
/// Raises `ValueError` if any successor index is out of range.
#[pyfunction]
fn outcomes(succ: Vec<Vec<usize>>) -> PyResult<Vec<PyOutcome>> {
    check_succ_bounds(&succ)?;
    Ok(crate::games::outcomes(&succ)
        .into_iter()
        .map(wrap_outcome)
        .collect())
}

/// The P-positions (Loss positions) of a finite game graph, as node indices.
/// Raises `ValueError` if any successor index is out of range.
#[pyfunction]
fn p_positions(succ: Vec<Vec<usize>>) -> PyResult<Vec<usize>> {
    check_succ_bounds(&succ)?;
    Ok(crate::games::p_positions(&succ))
}

#[pyclass(name = "ScoreInterval", module = "ogdoad")]
struct PyScoreInterval {
    inner: crate::games::ScoreInterval,
}

#[pymethods]
impl PyScoreInterval {
    #[getter]
    fn left(&self) -> i128 {
        self.inner.left
    }
    #[getter]
    fn right(&self) -> i128 {
        self.inner.right
    }
    fn __repr__(&self) -> String {
        format!(
            "ScoreInterval(left={}, right={})",
            self.inner.left, self.inner.right
        )
    }
}

fn wrap_score_interval(inner: crate::games::ScoreInterval) -> PyScoreInterval {
    PyScoreInterval { inner }
}

/// Milnor scoring-game minimax on a finite **acyclic** graph: the `(left, right)`
/// value interval of every position (`left` = optimal score with Left/maximizer
/// to move, `right` with Right/minimizer), where `terminal_score[v]` scores each
/// move-less position. Errors on a cycle (loopy scoring is out of scope) or if any
/// successor index is out of range.
#[pyfunction]
fn scoring_values(
    succ: Vec<Vec<usize>>,
    terminal_score: Vec<i128>,
) -> PyResult<Vec<PyScoreInterval>> {
    check_succ_bounds(&succ)?;
    crate::games::scoring_values(&succ, &terminal_score)
        .map(|v| v.into_iter().map(wrap_score_interval).collect())
        .ok_or_else(|| PyValueError::new_err("graph has a cycle — scoring value is undefined"))
}

// ---------------------------------------------------------------------------
// Misère play: Nim witnesses, octal moves, and the indistinguishability quotient
// ---------------------------------------------------------------------------

/// Normalize a Nim heap-multiset: drop empty heaps and sort ascending.
#[pyfunction]
fn nim_canonical(heaps: Vec<u128>) -> Vec<u128> {
    crate::games::nim_canonical(heaps)
}

/// The closed-form misère-Nim P-position rule: with every nonzero heap equal to
/// `1`, a P iff an odd number of nonzero heaps; otherwise (some heap `≥ 2`) a P
/// iff the XOR is `0`.
#[pyfunction]
fn misere_nim_p_predicted(heaps: Vec<u128>) -> bool {
    crate::games::misere_nim_p_predicted(&heaps)
}

/// Misère outcome of a finite acyclic impartial game with `u128` positions and a
/// Python move callback. Returns `None` if the callback-defined graph has a cycle.
#[pyfunction]
fn try_misere_is_n(pos: u128, moves: Bound<'_, PyAny>) -> PyResult<Option<bool>> {
    misere_is_n_u128_inner(pos, &moves, &mut HashMap::new(), &mut HashSet::new())
}

/// Misère N/P outcome for an acyclic Python callback game; errors on cycles.
#[pyfunction]
fn misere_is_n(pos: u128, moves: Bound<'_, PyAny>) -> PyResult<bool> {
    try_misere_is_n(pos, moves)?
        .ok_or_else(|| PyValueError::new_err("misere_is_n requires an acyclic move graph"))
}

/// `True` iff a callback-defined acyclic game is a misère P-position.
#[pyfunction]
fn misere_is_p(pos: u128, moves: Bound<'_, PyAny>) -> PyResult<bool> {
    Ok(!misere_is_n(pos, moves)?)
}

/// The moves of Nim: reduce any one heap to any strictly smaller size.
#[pyfunction]
fn nim_moves(pos: Vec<u128>) -> Vec<Vec<u128>> {
    crate::games::nim_moves(&pos)
}

/// The moves of an octal game `0.d₁d₂…` (`code[k-1] = dₖ`) on a heap-multiset:
/// remove `k` tokens leaving the heap empty (`dₖ & 1`), one nonempty heap
/// (`dₖ & 2`), or two nonempty heaps (`dₖ & 4`). Nim is `0.333…`.
#[pyfunction]
fn octal_moves(code: Vec<u128>, pos: Vec<u128>) -> Vec<Vec<u128>> {
    crate::games::octal_moves(&code, &pos)
}

/// The bounded misère indistinguishability quotient of an octal game, over single
/// heaps `1..=max_heap` as atoms (elements are sums up to `elem_bound`, separated
/// by tests up to `test_bound`).
#[pyfunction]
fn octal_misere_quotient(
    code: Vec<u128>,
    max_heap: usize,
    elem_bound: usize,
    test_bound: usize,
) -> PyQuotient {
    PyQuotient {
        inner: crate::games::octal_misere_quotient(&code, max_heap, elem_bound, test_bound),
    }
}

/// Loopy impartial nim-values of a (possibly cyclic) game graph: each position is
/// a typed `LoopyNimber.Value(n)`, or `LoopyNimber.Side` for a Draw position.
/// Errors when a cyclic non-Draw subgraph has no unique bounded sidling solution,
/// or when any successor index is out of range.
#[pyfunction]
fn loopy_nim_values(succ: Vec<Vec<usize>>) -> PyResult<Vec<PyLoopyNimber>> {
    check_succ_bounds(&succ)?;
    crate::games::loopy_nim_values(&succ)
        .map(|vs| vs.into_iter().map(wrap_loopy_nimber).collect())
        .ok_or_else(|| {
            PyValueError::new_err("cyclic non-Draw subgraph has no unique bounded sidling solution")
        })
}

/// `(loss_set, draw_set)` for a cyclic Python move rule on positions `0..n`.
#[pyfunction]
fn loopy_decision_sets(n: usize, moves: Bound<'_, PyAny>) -> PyResult<(Vec<usize>, Vec<usize>)> {
    let succ = loopy_succ_from_callback(n, &moves)?;
    let g = LoopyGraph::new(succ);
    Ok((g.loss_set(), g.draw_set()))
}

/// Fit F₂ quadrics to the loss and draw sets of a callback-defined loopy rule on
/// `F_2^k` (`positions = 0..2^k`).
#[pyfunction]
fn loopy_quadric_probe(
    k: usize,
    moves: Bound<'_, PyAny>,
) -> PyResult<(Option<PyQuadricFit>, Option<PyQuadricFit>)> {
    const MAX_ANF_DIM: usize = 20;
    if k > MAX_ANF_DIM {
        return Err(PyValueError::new_err(format!(
            "loopy_quadric_probe is exponential in k; max supported k is {MAX_ANF_DIM}"
        )));
    }
    let n = 1usize << k;
    let (loss, draw) = loopy_decision_sets(n, moves)?;
    let loss_u: Vec<u128> = loss.into_iter().map(|v| v as u128).collect();
    let draw_u: Vec<u128> = draw.into_iter().map(|v| v as u128).collect();
    Ok((
        crate::forms::fit_f2_quadratic(&loss_u, k).map(wrap_quadric_fit),
        crate::forms::fit_f2_quadratic(&draw_u, k).map(wrap_quadric_fit),
    ))
}

#[pyclass(name = "LoopyNimCertificate", module = "ogdoad")]
struct PyLoopyNimCertificate {
    inner: LoopyNimCertificate,
}

#[pymethods]
impl PyLoopyNimCertificate {
    #[getter]
    fn outcomes(&self) -> Vec<PyOutcome> {
        self.inner
            .outcomes
            .iter()
            .copied()
            .map(wrap_outcome)
            .collect()
    }
    #[getter]
    fn side_positions(&self) -> Vec<usize> {
        self.inner.side_positions.clone()
    }
    #[getter]
    fn used_sidling_solver(&self) -> bool {
        self.inner.used_sidling_solver
    }
    #[getter]
    fn sidling_assignments_examined(&self) -> usize {
        self.inner.sidling_assignments_examined
    }
    fn __repr__(&self) -> String {
        format!(
            "LoopyNimCertificate(side_positions={:?}, used_sidling_solver={}, sidling_assignments_examined={})",
            self.inner.side_positions,
            self.inner.used_sidling_solver,
            self.inner.sidling_assignments_examined
        )
    }
}

/// Loopy nim-values plus a certificate explaining Draw/Side promotion and
/// whether the bounded sidling solver was needed.
/// Raises `ValueError` if any successor index is out of range.
#[pyfunction]
fn loopy_nim_values_certified(
    succ: Vec<Vec<usize>>,
) -> PyResult<(Vec<PyLoopyNimber>, PyLoopyNimCertificate)> {
    check_succ_bounds(&succ)?;
    crate::games::loopy_nim_values_certified(&succ)
        .map(|(vs, inner)| {
            let values = vs.into_iter().map(wrap_loopy_nimber).collect();
            (values, PyLoopyNimCertificate { inner })
        })
        .ok_or_else(|| {
            PyValueError::new_err("cyclic non-Draw subgraph has no unique bounded sidling solution")
        })
}

/// Exact thermograph object for a short game, with `Rational` wall coordinates.
#[pyfunction]
fn thermograph(game: &PyGame) -> Option<PyThermograph> {
    thermography::thermograph(&game.inner).map(wrap_thermograph)
}

/// The same exact thermograph, computed through the named tropical wall folds.
#[pyfunction]
fn thermograph_via_tropical(game: &PyGame) -> Option<PyThermograph> {
    crate::games::tropical_thermography::thermograph_via_tropical(&game.inner).map(wrap_thermograph)
}

#[pyfunction]
fn temperature(game: &PyGame) -> Option<PySurreal> {
    crate::games::temperature(&game.inner).map(rat_to_py)
}

#[pyfunction]
fn mean_value(game: &PyGame) -> Option<PySurreal> {
    crate::games::mean_value(&game.inner).map(rat_to_py)
}

#[pyfunction]
fn left_stop(game: &PyGame) -> Option<PySurreal> {
    crate::games::left_stop(&game.inner).map(rat_to_py)
}

#[pyfunction]
fn right_stop(game: &PyGame) -> Option<PySurreal> {
    crate::games::right_stop(&game.inner).map(rat_to_py)
}

/// The atomic weight as a `Game` (`None` if the input is not all-small).
#[pyfunction]
fn atomic_weight(game: &PyGame) -> Option<PyGame> {
    crate::games::atomic_weight(&game.inner).map(|inner| PyGame { inner })
}

/// The atomic weight as an integer (`None` if undefined or genuinely non-integer).
#[pyfunction]
fn atomic_weight_int(game: &PyGame) -> Option<i128> {
    crate::games::atomic_weight_int(&game.inner)
}

// ---------------------------------------------------------------------------
// Partizan games + the exterior algebra of the game group
// ---------------------------------------------------------------------------

#[pyclass(name = "Game", module = "ogdoad", from_py_object)]
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
    /// The Nim-heap `⋆n`, matching Rust's `Game::nim_heap` name.
    #[staticmethod]
    fn nim_heap(n: u128) -> PyGame {
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
    /// An order-independent structural string of the game tree as given; no
    /// canonicalization.
    fn structural_string(&self) -> String {
        self.inner.structural_string()
    }
    /// Structural equality of the game trees as given, without canonicalization.
    fn structural_eq(&self, other: &PyGame) -> bool {
        self.inner.structural_eq(&other.inner)
    }
    /// A readable structural form: `0` for `{|}`, else `{L|R}` recursively.
    fn display(&self) -> String {
        self.inner.display()
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
    /// The thermograph as a first-class exact object with `Pl` walls.
    fn thermograph(&self) -> Option<PyThermograph> {
        thermography::thermograph(&self.inner).map(wrap_thermograph)
    }
    /// The same thermograph, routed through the named tropical max-plus/min-plus
    /// wall folds and pinned equal to `thermograph` in the Rust tests.
    fn thermograph_via_tropical(&self) -> Option<PyThermograph> {
        crate::games::tropical_thermography::thermograph_via_tropical(&self.inner)
            .map(wrap_thermograph)
    }
    /// Cooled stops `(LS(G_t), RS(G_t))` at the rational temperature `num/den`.
    #[pyo3(signature = (num, den=1))]
    fn cooled_stops(&self, num: i128, den: i128) -> PyResult<Option<(PySurreal, PySurreal)>> {
        let t = Rational::try_new(num, den)
            .ok_or_else(|| PyValueError::new_err("zero denominator or bounded i128 overflow"))?;
        Ok(thermography::thermograph(&self.inner).map(|th| {
            let (l, r) = th.cooled_stops(&t);
            (rat_to_py(l), rat_to_py(r))
        }))
    }
    fn __repr__(&self) -> String {
        self.inner.display()
    }
}

fn color_name(c: Color) -> String {
    match c {
        Color::Blue => "blue",
        Color::Red => "red",
        Color::Green => "green",
    }
    .to_string()
}

#[pyclass(name = "Color", module = "ogdoad", from_py_object)]
#[derive(Clone)]
struct PyColor {
    inner: Color,
}

fn wrap_color(inner: Color) -> PyColor {
    PyColor { inner }
}

fn parse_color_obj(obj: &Bound<'_, PyAny>) -> PyResult<Color> {
    if let Ok(color) = obj.cast::<PyColor>() {
        return Ok(color.borrow().inner);
    }
    Err(PyTypeError::new_err("expected Color"))
}

#[pymethods]
impl PyColor {
    #[staticmethod]
    fn blue() -> Self {
        wrap_color(Color::Blue)
    }
    #[staticmethod]
    fn red() -> Self {
        wrap_color(Color::Red)
    }
    #[staticmethod]
    fn green() -> Self {
        wrap_color(Color::Green)
    }
    fn name(&self) -> String {
        color_name(self.inner)
    }
    fn __richcmp__(&self, other: &Bound<'_, PyAny>, op: CompareOp) -> PyResult<bool> {
        match op {
            CompareOp::Eq => Ok(parse_color_obj(other).is_ok_and(|c| c == self.inner)),
            CompareOp::Ne => Ok(parse_color_obj(other).is_ok_and(|c| c != self.inner)),
            CompareOp::Lt | CompareOp::Le | CompareOp::Gt | CompareOp::Ge => Err(
                PyValueError::new_err("Color only supports equality comparisons"),
            ),
        }
    }
    fn __str__(&self) -> String {
        self.name()
    }
    fn __repr__(&self) -> String {
        match self.inner {
            Color::Blue => "Color.Blue".to_string(),
            Color::Red => "Color.Red".to_string(),
            Color::Green => "Color.Green".to_string(),
        }
    }
}

#[pyclass(name = "Hackenbush", module = "ogdoad")]
struct PyHackenbush {
    inner: Hackenbush,
}

#[pymethods]
impl PyHackenbush {
    /// A position from `(u, v, Color)` edges; vertex `0` is the ground.
    #[new]
    fn new(edges: Vec<(usize, usize, PyColor)>) -> Self {
        let edges = edges.into_iter().map(|(u, v, c)| (u, v, c.inner)).collect();
        PyHackenbush {
            inner: Hackenbush::new(edges),
        }
    }
    /// A stalk `0—1—2—…` from the ground, edge `i` coloured `colors[i]`.
    #[staticmethod]
    fn string(colors: Vec<PyColor>) -> Self {
        let cs = colors.into_iter().map(|c| c.inner).collect::<Vec<_>>();
        PyHackenbush {
            inner: Hackenbush::string(&cs),
        }
    }
    /// The edges `(u, v, Color)` as typed Python `Color` values.
    fn edges(&self) -> Vec<(usize, usize, PyColor)> {
        self.inner
            .edges()
            .iter()
            .map(|&(u, v, c)| (u, v, wrap_color(c)))
            .collect()
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

#[pyclass(name = "GameExterior", module = "ogdoad")]
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
    #[staticmethod]
    fn with_relation_search(gens: Vec<PyGame>, bound: i128) -> Self {
        Self::with_relation_bound(gens, bound)
    }
    #[staticmethod]
    fn with_relations(gens: Vec<PyGame>, relations: Vec<PyGameRelation>) -> Self {
        let games: Vec<Game> = gens.iter().map(|g| g.inner.clone()).collect();
        let relations = relations.into_iter().map(|rel| rel.inner).collect();
        PyGameExterior::from_inner(GameExterior::with_relations(games, relations))
    }
    #[getter]
    fn dim(&self) -> usize {
        self.inner.algebra().dim
    }
    /// The underlying free Grassmann algebra before quotienting by game-group
    /// relations. Use `reduce`/`wedge`/`add` on `GameExterior` for quotient-aware
    /// operations.
    fn algebra(&self) -> IntegerAlgebra {
        IntegerAlgebra {
            inner: self.alg.clone(),
        }
    }
    fn relations(&self) -> Vec<PyGameRelation> {
        self.inner
            .relations()
            .iter()
            .cloned()
            .map(wrap_game_relation)
            .collect()
    }
    /// Whether the automatic bounded relation search exhausted its coefficient
    /// box. Explicit relations always report true.
    fn relation_search_complete(&self) -> bool {
        self.inner.relation_search_complete()
    }
    /// Full relation-search certificate as a named record.
    fn relation_search_certificate(&self) -> PyRelationSearchCertificate {
        wrap_relation_search_certificate(self.inner.relation_search_certificate().clone())
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
    fn scalar_mul(&self, s: i128, mv: &IntegerMV) -> PyResult<IntegerMV> {
        self.ensure_mv(mv)?;
        Ok(IntegerMV {
            alg: self.alg.clone(),
            mv: self.inner.scalar_mul(s, &mv.mv),
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
#[pyclass(name = "NumberGame", module = "ogdoad", from_py_object)]
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
    /// The birthday as an `Ordinal`, matching Rust's `NumberGame::birthday` name.
    fn birthday(&self) -> Option<PyOrdinal> {
        self.inner.birthday().map(PyOrdinal::from_inner)
    }
    /// The birthday as an `Ordinal`, when the value is in the representable
    /// sign-expansion subclass.
    fn birthday_ordinal(&self) -> Option<PyOrdinal> {
        self.inner.birthday().map(PyOrdinal::from_inner)
    }
    /// The birthday as an ordinal string (`None` outside the representable
    /// subclass, e.g. `√ω`).
    fn birthday_repr(&self) -> Option<String> {
        self.inner.birthday().map(|o| format!("{o:?}"))
    }
    /// The transfinite sign expansion as runs `(sign, length)` (`True = +`,
    /// length an `Ordinal`), the finite encoding of the number-game tree.
    fn sign_expansion(&self) -> Option<Vec<(bool, PyOrdinal)>> {
        self.inner.sign_expansion().map(|se| {
            se.runs()
                .iter()
                .map(|(s, l)| (*s, PyOrdinal::from_inner(l.clone())))
                .collect()
        })
    }
    /// Reconstruct a number game from transfinite sign-expansion runs
    /// `(sign, length)`.
    #[staticmethod]
    fn from_sign_expansion(runs: Vec<(bool, PyOrdinal)>) -> Option<PyNumberGame> {
        let se = SignExpansion::from_runs(
            runs.into_iter()
                .map(|(sign, len)| (sign, len.as_ordinal().clone()))
                .collect(),
        );
        NumberGame::from_sign_expansion(&se).map(|inner| PyNumberGame { inner })
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

/// A transfinite **nimber-valued** (impartial) game — the Nim heap `⋆α` (e.g.
/// `⋆ω`) — carried by its ordinal Grundy value. The char-2 mirror of `NumberGame`:
/// Grundy value / disjunctive sum (XOR) / Turning-Corners product all delegate to
/// the `Ordinal` (`On₂`) backend, completing the `No ↔ On₂` symmetry at the games
/// layer.
#[pyclass(name = "NimberGame", module = "ogdoad", from_py_object)]
#[derive(Clone)]
struct PyNimberGame {
    inner: NimberGame,
}

#[pymethods]
impl PyNimberGame {
    /// The transfinite Nim heap `⋆α` of a given ordinal Grundy value.
    #[staticmethod]
    fn from_ordinal(o: &PyOrdinal) -> PyNimberGame {
        PyNimberGame {
            inner: NimberGame::from_ordinal(o.as_ordinal()),
        }
    }
    /// The finite Nim heap `⋆n`.
    #[staticmethod]
    fn nim_heap(n: u128) -> PyNimberGame {
        PyNimberGame {
            inner: NimberGame::nim_heap(n),
        }
    }
    /// The exact Grundy value (a transfinite nimber). P-position ⟺ this is `0`.
    fn grundy(&self) -> PyOrdinal {
        PyOrdinal::from_inner(self.inner.grundy().clone())
    }
    /// The Grundy value as a finite nimber, if the heap is finite (`< ω`).
    fn grundy_finite(&self) -> Option<u128> {
        self.inner.grundy().as_finite()
    }
    /// The short `Game`, if the heap is finite; `None` for transfinite heaps.
    fn to_finite_game(&self) -> Option<PyGame> {
        self.inner.to_finite_game().map(|inner| PyGame { inner })
    }
    /// The **Turning-Corners product** (nim-multiplication); `None` only past the
    /// verified Kummer table or at `≥ ω^(ω^ω)`.
    fn turning_corners(&self, other: &PyNimberGame) -> Option<PyNimberGame> {
        self.inner
            .turning_corners(&other.inner)
            .map(|inner| PyNimberGame { inner })
    }
    fn __add__(&self, other: &PyNimberGame) -> PyNimberGame {
        PyNimberGame {
            inner: self.inner.add(&other.inner),
        }
    }
    fn __neg__(&self) -> PyNimberGame {
        PyNimberGame {
            inner: self.inner.neg(),
        }
    }
    fn __richcmp__(&self, other: &PyNimberGame, op: CompareOp) -> bool {
        op.matches(self.inner.cmp(&other.inner))
    }
    fn __repr__(&self) -> String {
        format!("NimberGame(⋆{:?})", self.inner.grundy())
    }
}

/// A bounded misère indistinguishability quotient: the elements (atom-multisets),
/// their class ids, the per-class representatives, and which classes are P.
#[pyclass(name = "Quotient", module = "ogdoad")]
struct PyQuotient {
    inner: Quotient,
}

#[pymethods]
impl PyQuotient {
    /// The enumerated elements (sorted atom-multisets, up to `elem_bound`).
    #[getter]
    fn elements(&self) -> Vec<Vec<usize>> {
        self.inner.elements.clone()
    }
    /// The test positions used to distinguish bounded quotient classes.
    #[getter]
    fn test_positions(&self) -> Vec<Vec<usize>> {
        self.inner.test_positions.clone()
    }
    /// Outcome signatures parallel to `elements`: True means N-position.
    #[getter]
    fn signatures(&self) -> Vec<Vec<bool>> {
        self.inner.signatures.clone()
    }
    /// The class id of each element (parallel to `elements`).
    #[getter]
    fn class_of(&self) -> Vec<usize> {
        self.inner.class_of.clone()
    }
    /// The number of distinct classes (the order of the bounded quotient monoid).
    #[getter]
    fn num_classes(&self) -> usize {
        self.inner.num_classes
    }
    /// A representative multiset for each class.
    #[getter]
    fn class_rep(&self) -> Vec<Vec<usize>> {
        self.inner.class_rep.clone()
    }
    /// P-status of each class (`True` = a misère P-position / second-player win).
    #[getter]
    fn class_is_p(&self) -> Vec<bool> {
        self.inner.class_is_p.clone()
    }
    /// Class multiplication table, if represented at the current bounds.
    #[getter]
    fn multiplication(&self) -> Option<Vec<Vec<usize>>> {
        self.inner.multiplication.clone()
    }
    /// Whether every represented product agrees with the multiplication table.
    #[getter]
    fn multiplication_consistent(&self) -> bool {
        self.inner.multiplication_consistent
    }
    /// Whether the bounded element set is closed under disjunctive sum.
    #[getter]
    fn elements_closed_under_sum(&self) -> bool {
        self.inner.elements_closed_under_sum
    }
    /// Whether represented classes carry a complete sampled monoid table.
    #[getter]
    fn has_complete_bounded_monoid(&self) -> bool {
        self.inner.has_complete_bounded_monoid()
    }
    /// Product of quotient classes, if represented at the current bounds.
    fn class_product(&self, a: usize, b: usize) -> Option<usize> {
        self.inner.class_product(a, b)
    }
    /// Exact outcome signature for an enumerated element.
    fn signature_of_element(&self, element_index: usize) -> Option<Vec<bool>> {
        self.inner
            .signature_of_element(element_index)
            .map(<[bool]>::to_vec)
    }
    fn __repr__(&self) -> String {
        format!(
            "Quotient(num_classes={}, elements={})",
            self.inner.num_classes,
            self.inner.elements.len()
        )
    }
}

/// An abstract finite impartial game for misère analysis: position `0` is the
/// empty game (no moves); position `p` carries the option list `moves[p]` (each
/// option is a position index; `0` = move to empty).
#[pyclass(name = "AbstractGame", module = "ogdoad")]
struct PyAbstractGame {
    inner: AbstractGame,
}

#[pymethods]
impl PyAbstractGame {
    #[new]
    fn new(moves: Vec<Vec<usize>>) -> Self {
        PyAbstractGame {
            inner: AbstractGame { moves },
        }
    }
    /// Misère outcome of a disjunctive sum (a multiset of component positions):
    /// `True` = N (next player / first-player win), `False` = P.
    /// Raises `ValueError` if the move graph has a cycle (e.g. a position whose
    /// option list references itself or forms an indirect cycle).
    fn misere_outcome(&self, pos: Vec<usize>) -> PyResult<bool> {
        let mut memo = std::collections::HashMap::new();
        self.inner.misere_outcome(&pos, &mut memo).ok_or_else(|| {
            PyValueError::new_err("misere_outcome: move graph has a cycle — outcome is undefined")
        })
    }
    /// The bounded misère indistinguishability quotient over the generating
    /// `atoms` (elements are sums up to `elem_bound`, tests up to `test_bound`).
    fn misere_quotient(
        &self,
        atoms: Vec<usize>,
        elem_bound: usize,
        test_bound: usize,
    ) -> PyQuotient {
        PyQuotient {
            inner: crate::games::misere_quotient(&self.inner, &atoms, elem_bound, test_bound),
        }
    }
}

/// Rust-name module-level wrapper for `games::misere_quotient`; Python passes
/// the `AbstractGame` value explicitly.
#[pyfunction]
fn misere_quotient(
    game: &PyAbstractGame,
    atoms: Vec<usize>,
    elem_bound: usize,
    test_bound: usize,
) -> PyQuotient {
    PyQuotient {
        inner: crate::games::misere_quotient(&game.inner, &atoms, elem_bound, test_bound),
    }
}

/// A loopy game as a finite move graph (`succ[v]` = positions reachable from `v`);
/// the graph may be cyclic. Outcomes come from the retrograde kernel analysis
/// (Win / Loss / Draw, where Loss = P-position and Draw is the loopy escape).
#[pyclass(name = "LoopyValue", module = "ogdoad", from_py_object)]
#[derive(Clone)]
struct PyLoopyValue {
    inner: LoopyValue,
}

#[pymethods]
impl PyLoopyValue {
    #[staticmethod]
    fn zero() -> Self {
        PyLoopyValue {
            inner: LoopyValue::Zero,
        }
    }
    #[staticmethod]
    fn star() -> Self {
        PyLoopyValue {
            inner: LoopyValue::Star,
        }
    }
    #[staticmethod]
    fn on() -> Self {
        PyLoopyValue {
            inner: LoopyValue::On,
        }
    }
    #[staticmethod]
    fn off() -> Self {
        PyLoopyValue {
            inner: LoopyValue::Off,
        }
    }
    #[staticmethod]
    fn over() -> Self {
        PyLoopyValue {
            inner: LoopyValue::Over,
        }
    }
    #[staticmethod]
    fn under() -> Self {
        PyLoopyValue {
            inner: LoopyValue::Under,
        }
    }
    #[staticmethod]
    fn dud() -> Self {
        PyLoopyValue {
            inner: LoopyValue::Dud,
        }
    }
    fn name(&self) -> &'static str {
        self.inner.name()
    }
    fn form(&self) -> &'static str {
        self.inner.form()
    }
    fn outcome(&self) -> PyPartizanOutcome {
        wrap_partizan_outcome(self.inner.outcome())
    }
    fn __neg__(&self) -> PyLoopyValue {
        PyLoopyValue {
            inner: self.inner.neg(),
        }
    }
    fn is_stopper(&self) -> bool {
        self.inner.is_stopper()
    }
    fn __add__(&self, other: &PyLoopyValue) -> Option<PyLoopyValue> {
        self.inner
            .add(&other.inner)
            .map(|inner| PyLoopyValue { inner })
    }
    fn __richcmp__(&self, other: &PyLoopyValue, op: CompareOp) -> bool {
        match op {
            CompareOp::Eq => self.inner == other.inner,
            CompareOp::Ne => self.inner != other.inner,
            CompareOp::Lt | CompareOp::Le | CompareOp::Gt | CompareOp::Ge => self
                .inner
                .partial_cmp(&other.inner)
                .is_some_and(|ordering| op.matches(ordering)),
        }
    }
    fn __repr__(&self) -> String {
        format!("LoopyValue({:?})", self.inner)
    }
}

#[pyclass(name = "LoopyGraph", module = "ogdoad")]
struct PyLoopyGraph {
    inner: LoopyGraph,
}

#[pymethods]
impl PyLoopyGraph {
    #[new]
    fn new(succ: Vec<Vec<usize>>) -> PyResult<Self> {
        check_succ_bounds(&succ)?;
        Ok(PyLoopyGraph {
            inner: LoopyGraph::new(succ),
        })
    }
    #[staticmethod]
    fn from_rule(n: usize, moves: Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(PyLoopyGraph {
            inner: LoopyGraph::new(loopy_succ_from_callback(n, &moves)?),
        })
    }
    /// The adjacency lists.
    fn succ(&self) -> Vec<Vec<usize>> {
        self.inner.succ().to_vec()
    }
    /// Typed `Outcome` of every position.
    fn outcomes(&self) -> Vec<PyOutcome> {
        self.inner
            .outcomes()
            .into_iter()
            .map(wrap_outcome)
            .collect()
    }
    /// The Loss positions = P-positions (the player to move loses).
    fn loss_set(&self) -> Vec<usize> {
        self.inner.loss_set()
    }
    /// The Win positions = N-positions (the player to move wins).
    fn win_set(&self) -> Vec<usize> {
        self.inner.win_set()
    }
    /// The Draw positions — the loopy degree of freedom.
    fn draw_set(&self) -> Vec<usize> {
        self.inner.draw_set()
    }
    /// A coarse catalogue reading of position `v` via its impartial outcome:
    /// `LoopyValue.Zero` for a Loss, `LoopyValue.Dud` for a Draw, `None` for a Win (a nonzero loopy
    /// nimber — use `loopy_nim_values`).
    fn classify(&self, v: usize) -> Option<PyLoopyValue> {
        self.inner.classify(v).map(|inner| PyLoopyValue { inner })
    }
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyGame>()?;
    m.add_class::<PyOutcome>()?;
    m.add_class::<PyPartizanOutcome>()?;
    m.add_class::<PyLoopyNimber>()?;
    m.add_class::<PyColor>()?;
    m.add_class::<PyPl>()?;
    m.add_class::<PyThermograph>()?;
    m.add_class::<PyNumberGame>()?;
    m.add_class::<PyNimberGame>()?;
    m.add_class::<PyGameExterior>()?;
    m.add_class::<PyGameRelation>()?;
    m.add_class::<PyGameRelationCertificate>()?;
    m.add_class::<PyRelationSearchCertificate>()?;
    m.add_class::<PyScoreInterval>()?;
    m.add_class::<PyHackenbush>()?;
    m.add_class::<PyQuotient>()?;
    m.add_class::<PyAbstractGame>()?;
    m.add_class::<PyLoopyValue>()?;
    m.add_class::<PyLoopyGraph>()?;
    m.add_class::<PyLoopyNimCertificate>()?;
    m.add_function(wrap_pyfunction!(nim_mul_mex, m)?)?;
    m.add_function(wrap_pyfunction!(coin_companions, m)?)?;
    m.add_function(wrap_pyfunction!(singleton_companions, m)?)?;
    m.add_function(wrap_pyfunction!(turtles_companions, m)?)?;
    m.add_function(wrap_pyfunction!(grundy_1d, m)?)?;
    m.add_function(wrap_pyfunction!(coin_turning_grundy, m)?)?;
    m.add_function(wrap_pyfunction!(coin_turning_tartan_grundy, m)?)?;
    m.add_function(wrap_pyfunction!(tartan_grundy, m)?)?;
    m.add_function(wrap_pyfunction!(grundy, m)?)?;
    m.add_function(wrap_pyfunction!(grundy_graph, m)?)?;
    m.add_function(wrap_pyfunction!(mex, m)?)?;
    m.add_function(wrap_pyfunction!(outcomes, m)?)?;
    m.add_function(wrap_pyfunction!(p_positions, m)?)?;
    m.add_function(wrap_pyfunction!(scoring_values, m)?)?;
    m.add_function(wrap_pyfunction!(nim_canonical, m)?)?;
    m.add_function(wrap_pyfunction!(misere_nim_p_predicted, m)?)?;
    m.add_function(wrap_pyfunction!(try_misere_is_n, m)?)?;
    m.add_function(wrap_pyfunction!(misere_is_n, m)?)?;
    m.add_function(wrap_pyfunction!(misere_is_p, m)?)?;
    m.add_function(wrap_pyfunction!(nim_moves, m)?)?;
    m.add_function(wrap_pyfunction!(octal_moves, m)?)?;
    m.add_function(wrap_pyfunction!(octal_misere_quotient, m)?)?;
    m.add_function(wrap_pyfunction!(misere_quotient, m)?)?;
    m.add_function(wrap_pyfunction!(loopy_nim_values, m)?)?;
    m.add_function(wrap_pyfunction!(loopy_decision_sets, m)?)?;
    m.add_function(wrap_pyfunction!(loopy_quadric_probe, m)?)?;
    m.add_function(wrap_pyfunction!(loopy_nim_values_certified, m)?)?;
    m.add_function(wrap_pyfunction!(thermograph, m)?)?;
    m.add_function(wrap_pyfunction!(thermograph_via_tropical, m)?)?;
    m.add_function(wrap_pyfunction!(temperature, m)?)?;
    m.add_function(wrap_pyfunction!(mean_value, m)?)?;
    m.add_function(wrap_pyfunction!(left_stop, m)?)?;
    m.add_function(wrap_pyfunction!(right_stop, m)?)?;
    m.add_function(wrap_pyfunction!(atomic_weight, m)?)?;
    m.add_function(wrap_pyfunction!(atomic_weight_int, m)?)?;
    Ok(())
}
