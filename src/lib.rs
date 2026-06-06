//! pleroma — Clifford algebras (with nilpotents) over the field-like
//! subclasses of combinatorial games.
//!
//! Stage 1 (this file's current scope): the pure-Rust math core.
//!   - `scalar`  : the Scalar trait + an exact Rational for engine validation
//!   - `nimber`  : On_2 (characteristic 2) — exact, the novel backend
//!   - `clifford`: the multivector engine, generic over Scalar  [next]
//!   - `surreal` : Conway normal form scalars (characteristic 0) [next]
//!   - `surcomplex`: adjoin i over any backend                    [next]
//!
//! Stage 2: PyO3 bindings + maturin packaging (abi3, robust on Python 3.14).

pub mod arf;
pub mod clifford;
pub mod nimber;
pub mod scalar;
pub mod surcomplex;
pub mod surreal;

#[cfg(feature = "python")]
mod py;
