//! Quadratic forms and their invariants, organised by the characteristic
//! trichotomy of the underlying scalar field.
//!
//! The classification of a quadratic form (equivalently, of the Clifford
//! algebra it builds) is *one* theory split three ways by `char F`:
//!
//! * [`char0`] — char 0: the 8-fold real and 2-fold complex tables on the
//!   exact-square subdomains represented by the scalar backends.
//! * [`oddchar`] — odd characteristic: discriminant + Hasse invariant.
//! * [`char2`] — characteristic 2: the Arf invariant (and Dickson).
//!
//! Three clusters of cross-cutting machinery sit beside the trichotomy, each in
//! its own shelf so the flat public API stays shallow:
//!
//! * [`witt`] — the **invariant groups**: the Witt group across all three legs
//!   ([`WittClassG`]), the Witt *ring* (`Iⁿ`, Pfister forms, the `eₙ`
//!   staircase), and the Brauer–Wall group ([`bw_class_real`], …).
//! * [`springer`] — the **valuation-graded (local↔global) decomposition** across
//!   the complete valued fields. The discretely-valued legs share **one** engine,
//!   [`springer_decompose_local`], keyed off the
//!   [`ResidueField`](crate::scalar::ResidueField) trait: `Q_p`/`Q_q`
//!   ([`springer_decompose_qp`]/[`springer_decompose_qq`]) and `F_q((t))`
//!   ([`springer_decompose_laurent`]); the surreal entry point
//!   ([`springer_decompose`]) is the one that does *not* fit — its value group is
//!   2-divisible, so the second residue map collapses — and keeps its own engine;
//!   that mismatch *is* the local–global symmetry, not a gap. The char-2 mirror
//!   ([`springer_decompose_local_char2`]) is the Aravire–Jacob three-layer story.
//! * [`integral`] — the **arithmetic view**: integral lattices, genus, mass,
//!   and the Leech lattice.
//!
//! [`classify`] is the façade over the trichotomy: which leg classifies a form
//! is a fact about the field, so [`ClassifyForm`] resolves it from the scalar
//! type — call `metric.classify()` / `algebra.classify()` (and `witt_class()`)
//! and the right leg is selected at compile time, no manual char-dispatch.
//!
//! Alongside the symmetric bilinear forms sit the other two members of the
//! "form + involution" family: [`symplectic`] alternating forms (rank is the
//! complete invariant, char-uniform) and [`hermitian`] forms over the surcomplex
//! field (Sylvester signature; [`HermitianForm::from_skew`] handles the
//! skew-Hermitian case via multiplication by `i`).
//!
//! The local–global layer is unified by [`global_field`] ([`GlobalField`]): the
//! local–global principle (places, Hilbert symbol, reciprocity `∏_v (a,b)_v = +1`,
//! Hasse–Minkowski) written **once** over the two kinds of global field, `ℚ`
//! ([`Rational`](crate::scalar::Rational)) and `F_q(t)`
//! ([`RationalFunction`](crate::scalar::RationalFunction)) — the `forms` mirror of
//! what [`ResidueField`](crate::scalar::ResidueField) did for the discrete Springer
//! engine. And [`trace_form`] bridges *into* this pillar from the field-growing
//! side: the Frobenius-twisted trace form `Tr_{E/F}(x·σ^k(x))` of a
//! [`CyclicGaloisExtension`](crate::scalar::CyclicGaloisExtension) is classified by
//! the same façade — the norm form over `Surcomplex`, the Gold form over the
//! nim-fields.
//!
//! [`field_invariants`] holds the numeric field invariants the Witt ring implies
//! (level/Stufe, Pythagoras number, u-invariant); [`quadric_fit`] is the
//! "is this P-set a quadric?" research bench fed by the game probes.

pub mod char0;
pub mod char2;
pub mod classify;
pub mod diagonalize;
pub mod equivalence;
pub mod field_invariants;
pub mod hermitian;
pub mod integral;
pub mod local_global;
pub mod oddchar;
pub(crate) mod poly_factor;
pub mod quadric_fit;
pub mod springer;
pub mod symplectic;
pub mod trace_form;
pub mod witt;

pub use char0::*;
pub use char2::*;
pub use classify::*;
pub use diagonalize::*;
pub use equivalence::*;
pub use field_invariants::*;
pub use hermitian::*;
pub use integral::*;
pub use local_global::*;
pub use oddchar::*;
pub use quadric_fit::*;
pub use springer::*;
pub use symplectic::*;
pub use trace_form::*;
pub use witt::*;
