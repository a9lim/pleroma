//! **Big** — the transfinite worlds, where the number is allowed to be infinite.
//! Conway normal form / Hahn series `Σ ω^{exp}·coeff` with recursive exponents.
//!
//!   * [`surreal`] — `No`, the real-closed char-0 field. The transfinite mirror of
//!     ℚ/ℝ; coefficients are ℚ (the honest finite truncation), exponents are fully
//!     recursive surreals.
//!   * [`omnific`] — `Oz ⊂ No`, the omnific *integers*: the ring of integers of the
//!     surreals, the transfinite mirror of ℤ (and the surreal mirror of `Z_p`).
//!   * [`onag`] — transfinite (ordinal) **nimbers**: the char-2 sibling of
//!     [`surreal`]. Same CNF representation, but coefficients combine by XOR
//!     (nim-addition) — `surreal : nimber :: No : On₂` extended to the ordinals.
//!
//! `surreal` and `onag` share the descending-CNF *shape* (a `Vec<(exponent,
//! coeff)>` recursing on exponents) and exactly one piece of *code*: the
//! [`cnf::merge_descending`] canonicalizer, into which both feed the three
//! primitives where they differ — the exponent order (`No`'s value order vs the
//! ordinal lexicographic order), the like-coefficient merge (ordinary `+` vs nim
//! `XOR`), and the zero test. Everything else (comparison, multiplication,
//! negation, the field structure) is backend-specific and stays so: `No` is a
//! real-closed field, `On₂` a characteristic-2 world with no negation. See
//! [`cnf`] for why this is a shared *function*, not a shared `Cnf<C>` *type*.

pub(crate) mod cnf;
pub mod omnific;
pub mod onag;
pub mod surreal;

pub use omnific::*;
pub use onag::*;
pub use surreal::*;
