//! The **exterior algebra of the game group**: `Λ` over `ℤ` on a chosen tuple of
//! games. This is the Clifford-adjacent structure that lives on *all* of
//! game-world — not just the field-like numbers — because the partizan games form
//! an abelian group (a `ℤ`-module), and the Grassmann algebra is the exterior
//! algebra of that module.
//!
//! Three layers, re-exported flat so every public path is unchanged:
//!
//!   * [`relations`] — [`GameRelation`], [`GameRelationCertificate`],
//!     [`RelationSearchCertificate`]: the relation and certificate record types.
//!   * [`lambda`] — [`GameExterior`]: the free Grassmann engine quotiented by
//!     integer game relations such as `2⋆=0`.
//!   * [`clifford`] — [`GameCliffordError`] and [`GameClifford`]: the checked
//!     integer-valued Clifford deformation surface; constructors verify that every
//!     game relation is null and polar-radical before accepting the metric.
//!
//! Generators may be non-numbers (`⋆`, `↑`, switches) — exactly where the
//! Clifford/scalar story cannot go — which is the point: the
//! [`Game`](crate::games::Game) group is not a ring, but it *is* a `ℤ`-module,
//! and that is enough for `Λ`. The stronger question of a natural game-native
//! source for the quadratic data remains open in `OPEN.md`.

pub mod clifford;
pub mod lambda;
pub mod relations;

pub use clifford::*;
pub use lambda::*;
pub use relations::*;
