//! **Exact base** — the Archimedean char-0 anchor of the number table: the field
//! ℚ and its ring of integers ℤ.
//!
//! Neither is a game backend. [`Rational`](rational::Rational) (ℚ over `i128`) is
//! the char-0 yardstick that validates the geometric product against the known
//! `Cl(p,q)` classification before the exotic backends are trusted;
//! [`Integer`](integer::Integer) (ℤ) is the coefficient ring of the exterior
//! algebra of the game group (games are a ℤ-module). They are the (field, ring of
//! integers) pair every other place in the table mirrors — `(No, Oz)`,
//! `(Q_p, Z_p)`, and the unramified mixed-characteristic pair
//! `(Q_q, W_n(F_q))`.

pub mod integer;
pub mod rational;

pub use integer::*;
pub use rational::*;
