//! The scalar interface every Clifford backend implements.
//!
//! A Clifford algebra needs a *commutative ring* of scalars. The whole point of
//! this project is that combinatorial games only supply such a ring on their
//! field-like subclasses — nimbers, surreals, surcomplex — so each of those is a
//! `Scalar` impl, and the multivector engine in `clifford/` is written once,
//! generic over this trait.
//!
//! This module is just the trait; every coefficient world is a sibling module,
//! re-exported flat (`scalar::Nimber`, `scalar::Surreal`, …). Two of them are
//! not game-backends but exact rings the engine needs: [`rational`] (ℚ, to
//! validate the geometric product in char 0) and [`integer`] (ℤ, the coefficient
//! ring for the exterior algebra of the game group).

// The coefficient worlds, each a commutative-ring `Scalar` backend, re-exported
// flat so call sites read `scalar::Nimber`, `scalar::Rational`, etc.
pub mod fp;
pub mod fpn;
pub mod integer;
pub mod nimber;
pub mod omnific;
pub mod onag;
pub mod rational;
pub mod surcomplex;
pub mod surreal;
pub mod wittvec;
pub mod zp;

pub use fp::*;
pub use fpn::*;
pub use integer::*;
pub use nimber::*;
pub use omnific::*;
pub use onag::*;
pub use rational::*;
pub use surcomplex::*;
pub use surreal::*;
pub use wittvec::*;
pub use zp::*;

use std::fmt::Debug;

pub trait Scalar: Clone + PartialEq + Debug {
    fn zero() -> Self;
    fn one() -> Self;
    fn add(&self, rhs: &Self) -> Self;
    fn neg(&self) -> Self;
    fn mul(&self, rhs: &Self) -> Self;

    /// Field characteristic: 0 for surreal / surcomplex / rational, 2 for
    /// nimbers. The engine reads this because char 2 is genuinely different:
    /// `-1 == 1` so blade-reordering signs vanish, and the quadratic form Q
    /// (the squares) must be carried independently of the alternating
    /// off-diagonal bilinear form B — see the notes in `clifford/engine.rs`.
    fn characteristic() -> u32;

    /// Multiplicative inverse, or `None` if not invertible (zero) or not
    /// finitely representable in this backend (e.g. a non-monomial surreal,
    /// whose inverse is an infinite Hahn series).
    fn inv(&self) -> Option<Self>;

    fn is_zero(&self) -> bool {
        *self == Self::zero()
    }

    fn sub(&self, rhs: &Self) -> Self {
        self.add(&rhs.neg())
    }
}
