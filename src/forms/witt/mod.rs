//! The invariant *groups* of quadratic forms: the Witt group, the Witt ring, and
//! the Brauer–Wall group — the three abelian groups the classifiers land in.
//!
//! * `class` — the Witt **group** `W_q(F)`: [`WittClass`] (the order-2 group of
//!   a finite nim-field, Arf-classified) and [`WittClassG`], the Char0/OddChar/Char2
//!   trichotomy enum that the classifier façade returns. (The Char2 leg is a
//!   *module*, not a ring — its `mul` panics; see `ring` for why.)
//! * `ring` — the Witt **ring**: [`tensor_form`], Pfister forms, the fundamental
//!   ideal `Iⁿ`, and the `eₙ` staircase (`e0 = dim`, `e1 = disc`, `e2 = Hasse`),
//!   with per-field stabilisation (`I² = 0` over `F_q`; the infinite ℝ tower).
//! * `brauer_wall` — the Brauer–Wall group `BW(F)`: [`bw_class_real`] (the Bott
//!   index `(q−p) mod 8`, so `BW(ℝ) ≅ ℤ/8`), [`bw_class_complex`] (`ℤ/2`),
//!   [`bw_class_finite_odd`] (order-4, `≅ W(F_q)`), and [`bw_class_nimber`] (the
//!   char-2 Arf/Witt class `ℤ/2`, nonsingular metrics only). The law is the
//!   graded tensor product.
//! * `brauer_rational` — the **ungraded** rational 2-torsion Brauer class
//!   ([`Brauer2Class`]) as a set of ramified places: the Hasse–Witt invariant
//!   ([`hasse_brauer_class`]) and the Clifford invariant ([`clifford_brauer_class`])
//!   of a `ℚ`-form, which differ by the explicit `n mod 8` / discriminant correction
//!   (Lam). The char-0/odd mirror of the char-2 Bridge B; kept strictly distinct
//!   from the graded `brauer_wall` class.
//! * `cyclic` — Bridge K: the **full `ℚ/ℤ`** ungraded Brauer class ([`BrauerClass`])
//!   and the cyclic-algebra local invariant [`cyclic_algebra_invariant`]
//!   (`inv = v(a)/n mod ℤ`, the unramified class). Lifts `brauer_rational`'s 2-torsion
//!   surface to the full local Brauer group, with [`Brauer2Class`] embedding as the
//!   `½`-slice ([`BrauerClass::from_two_torsion`]).
//!
//! The mod-8 spine lives here: `BW(ℝ) ≅ ℤ/8` is the same periodicity as the char-0
//! 8-fold Clifford table, Bott periodicity, and `E₈` as the rank-8 even unimodular
//! lattice (see [`integral`](crate::forms::integral)).
//!
//! Children are private modules re-exported flat, so the public API stays shallow
//! (`forms::WittClassG`, `forms::tensor_form`, `forms::bw_class_real`, …). The
//! numeric field invariants (level, u-invariant) the ring *implies* live separately
//! in [`field_invariants`](crate::forms::field_invariants).

mod brauer_rational;
mod brauer_wall;
mod class;
mod cyclic;
mod milnor;
mod ring;

pub use brauer_rational::*;
pub use brauer_wall::*;
pub use class::*;
pub use cyclic::*;
pub use milnor::*;
pub use ring::*;
