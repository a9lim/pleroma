//! The Clifford / geometric-algebra pillar.
//!
//! Two layers, deliberately separated:
//!
//!   * [`engine`] — the associative-algebra core: the `Metric` (carrying the
//!     quadratic form `q`, the alternating polar form `b`, and the optional
//!     asymmetric contraction `a` independently), the `Multivector` /
//!     `CliffordAlgebra` types, and the geometric product itself. This is the
//!     "associative algebra from a general bilinear form" primitive.
//!   * [`versor`] — the geometry built on top: versors and the Pin sandwich
//!     action, reflections, contractions, the pseudoscalar dual, grade
//!     involution, the spinor norm, and the even subalgebra.
//!
//! On top of those sit the structured-algebra modules: [`outermorphism`]
//! (lift a linear map to all grades; determinant via the pseudoscalar),
//! [`hopf`] (the exterior Hopf algebra), [`cga`] (conformal & projective GA),
//! and [`spinor`] (concrete minimal left ideals / spinor modules).
//!
//! Everything is re-exported flat, so downstream code reads `clifford::Metric`,
//! `clifford::sandwich`, `clifford::coproduct`, `clifford::up`, … regardless of
//! which sub-module an item lives in.

pub mod blade;
pub mod cga;
pub mod engine;
pub mod hopf;
pub mod outermorphism;
pub mod spinor;
pub mod spinor_norm;
pub mod versor;

pub use blade::*;
pub use cga::*;
pub use engine::*;
pub use hopf::*;
pub use outermorphism::*;
pub use spinor::*;
pub use spinor_norm::*;
// `versor` adds only inherent methods to `CliffordAlgebra` (reachable through
// the type itself), so there is nothing to glob-re-export from it.
