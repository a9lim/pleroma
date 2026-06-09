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
//! [`witt`] packages the Witt group across all three legs ([`WittClassG`]),
//! and the Springer decomposition is the valuation-graded decomposition across the
//! complete valued fields. The discretely-valued legs share **one** engine,
//! [`springer_local`] ([`springer_decompose_local`]), keyed off the
//! [`ResidueField`](crate::scalar::ResidueField) trait: [`springer_padic`] over
//! `Q_p`/`Q_q` (char 0, residue `F_p`/`F_q`) and [`springer_laurent`] over
//! `F_q((t))` (char `p`, residue `F_q`). [`springer`] over the surreals (char 0,
//! residue ℝ) is the one that does *not* fit — its value group is 2-divisible, so
//! the second residue map collapses — and keeps its own engine; that mismatch *is*
//! the local–global symmetry, not a gap.
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
//! The arithmetic side sits in two named shelves while preserving the flat public
//! API: [`local_global`] for Hilbert-symbol/Hasse–Minkowski machinery and
//! [`integral`] for lattices, genus, mass, and Leech.

pub mod brauer_wall;
pub mod char0;
pub mod char2;
pub mod classify;
pub mod diagonalize;
pub mod equivalence;
pub mod hermitian;
pub mod integral;
pub mod invariants;
pub mod local_global;
pub mod oddchar;
pub(crate) mod poly_factor;
pub mod quadric_fit;
pub mod springer;
pub mod springer_char2;
pub mod springer_laurent;
pub mod springer_local;
pub mod springer_padic;
pub mod symplectic;
pub mod trace_form;
pub mod witt;
pub mod witt_ring;

pub use brauer_wall::*;
pub use char0::*;
pub use char2::*;
pub use classify::*;
pub use diagonalize::*;
pub use equivalence::*;
pub use hermitian::*;
pub use integral::*;
pub use invariants::*;
pub use local_global::*;
pub use oddchar::*;
pub use quadric_fit::*;
pub use springer::*;
pub use springer_char2::*;
pub use springer_laurent::*;
pub use springer_local::*;
pub use springer_padic::*;
pub use symplectic::*;
pub use trace_form::*;
pub use witt::*;
pub use witt_ring::*;
