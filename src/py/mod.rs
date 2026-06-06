//! PyO3 bindings, split along the same pillars as the math core.
//!
//! Each scalar world (nimber / surreal / surcomplex / integer / omnific) gets
//! its own scalar type plus an `<World>Algebra` / `<World>MV` multivector pair,
//! stamped out by the `backend!` macro in [`engine`] — monomorphising the one
//! verified generic engine to a concrete scalar type, so there is no runtime
//! dispatch and no way to mix scalar worlds in one algebra.
//!
//!   - [`scalars`] — the scalar types, their constructors, nim-field ops.
//!   - [`engine`]  — the `backend!` macro, the algebra/MV pairs, conformal GA.
//!   - [`forms`]   — the classifier / invariant bindings (trichotomy + Witt).
//!   - [`games`]   — partizan games and the game-group exterior algebra.
//!
//! Each submodule registers its own classes and functions through a
//! `pub(crate) fn register`, which the `#[pymodule]` entry point chains
//! together.

use pyo3::prelude::*;

mod engine;
mod forms;
mod games;
mod scalars;

#[pyfunction]
fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

#[pymodule]
fn pleroma(m: &Bound<'_, PyModule>) -> PyResult<()> {
    scalars::register(m)?;
    engine::register(m)?;
    forms::register(m)?;
    games::register(m)?;
    m.add_function(wrap_pyfunction!(version, m)?)?;
    Ok(())
}
