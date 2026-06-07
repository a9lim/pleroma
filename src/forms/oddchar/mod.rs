//! The odd-characteristic Clifford / quadratic-form classifier — the third leg
//! of the trichotomy, companion to [`char0`](crate::forms::char0) and
//! [`char2`](crate::forms::char2).
//!
//! Over a finite field `F_q` of odd characteristic a nondegenerate quadratic
//! form is classified completely by **dimension + discriminant** (det mod
//! squares): for each dimension there are exactly two classes, distinguished by
//! whether the discriminant is a square. So the classifier is essentially
//! `(dim, disc-class)`.
//!
//! We also compute the **Hasse–Witt / Clifford invariant** (a product of Hilbert
//! symbols). Over a finite field this is *always* `+1` — finite fields have
//! trivial Brauer group, so there are no nontrivial quaternion algebras and the
//! Hilbert symbol of any two nonzero elements is `+1`. We compute it the honest
//! way (search for a representing vector, which always exists by
//! Chevalley–Warning) precisely to *exhibit* that triviality, and to make the
//! structural parallel with the Arf invariant explicit — not because it adds
//! classifying power over a finite field. The group-theoretic home of all this
//! is `witt::WittClassG`.

mod field;
mod invariants;

pub use field::*;
pub use invariants::*;

#[cfg(test)]
mod tests;
