//! Quadratic forms and their invariants, organised by the characteristic
//! trichotomy of the underlying scalar field.
//!
//! The classification of a quadratic form (equivalently, of the Clifford
//! algebra it builds) is *one* theory split three ways by `char F`:
//!
//!   * [`char0`]   — real-closed / algebraically-closed char 0: the 8-fold
//!                   (real) and 2-fold (complex) periodicity tables.
//!   * [`oddchar`] — odd characteristic: discriminant + Hasse invariant.
//!   * [`char2`]   — characteristic 2: the Arf invariant (and Dickson).
//!
//! [`witt`] packages the Witt group across all three legs ([`WittClassG`]),
//! and [`springer`] is the non-Archimedean valuation decomposition over the
//! surreals.
//!
//! [`classify`] is the façade over the trichotomy: which leg classifies a form
//! is a fact about the field, so [`ClassifyForm`] resolves it from the scalar
//! type — call `metric.classify()` / `algebra.classify()` (and `witt_class()`)
//! and the right leg is selected at compile time, no manual char-dispatch.

pub mod brauer_wall;
pub mod char0;
pub mod char2;
pub mod classify;
pub mod diagonalize;
pub mod equivalence;
pub mod hermitian;
pub mod invariants;
pub mod oddchar;
pub mod padic;
pub mod quadric_fit;
pub mod springer;
pub mod springer_padic;
pub mod witt;
pub mod witt_ring;

pub use brauer_wall::*;
pub use char0::*;
pub use char2::*;
pub use classify::*;
pub use diagonalize::*;
pub use equivalence::*;
pub use hermitian::*;
pub use invariants::*;
pub use oddchar::*;
pub use padic::*;
pub use quadric_fit::*;
pub use springer::*;
pub use springer_padic::*;
pub use witt::*;
pub use witt_ring::*;
