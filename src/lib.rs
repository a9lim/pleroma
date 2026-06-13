//! ogdoad — Clifford algebras (with nilpotents) over the field-like
//! subclasses of combinatorial games.
//!
//! Pure-Rust math core (generic over the `Scalar` trait), with optional PyO3
//! bindings behind the `python` feature (abi3). The source is organised into
//! four pillars plus the bindings:
//!
//! - [`scalar`] — the coefficient worlds: the `Scalar` trait, exact
//!   `Rational`/`Integer`, game-adjacent nimber/surreal backends, finite
//!   fields, p-adic/local functors, and the adelic precision model.
//! - [`clifford`] — the multivector engine (Metric + general bilinear form +
//!   geometric product), generic over `Scalar`, plus the GA layer:
//!   outermorphisms, the exterior Hopf algebra, conformal/projective GA, and
//!   spinor modules.
//! - [`forms`] — quadratic forms and their invariants across the characteristic
//!   trichotomy: char-0 / odd-char / char-2 classifiers, Witt/Brauer-Wall
//!   utilities, Springer decompositions, and rational local-global helpers.
//! - [`games`] — combinatorial game theory: coin-turning & Tartan products,
//!   normal-, misère-, and loopy finite-game probes, plus short partizan games
//!   and the exterior algebra of the game group.
//! - `py` — PyO3 per-backend bindings (feature = "python").
//!
//! See `AGENTS.md` for the mathematical layout and `docs/OPEN.md` for the open problems.

// This crate is matrix/algebra-heavy throughout: linalg solves, Gram matrices,
// Witt/carry formulas, Dickson/symplectic reductions, and spinor reps all walk
// index-parallel arrays where explicit `for i in 0..n` reads clearer than the
// iterator-adapter rewrite (the body indexes several arrays by the same `i`,
// or reads `out[i-1]` while writing `out[i]`). `needless_range_loop` is a false
// positive at every one of those sites, so it is allowed crate-wide here rather
// than suppressed piecemeal at a dozen matrix modules.
#![allow(clippy::needless_range_loop)]

pub mod clifford;
pub mod forms;
pub mod games;
pub(crate) mod linalg;
pub mod ogham;
pub mod scalar;

#[cfg(feature = "python")]
mod py;
