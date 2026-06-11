//! Loopy combinatorial games — games whose move graph may contain cycles, so
//! play need not terminate. This is the third escape (beside the interactive
//! [`kernel`](crate::games::kernel) route and the [`misere`](crate::games::misere)
//! route) from the XOR-linear P-sets of normal-play disjunctive sums: a cyclic
//! rule admits a **Draw** outcome — a position from which neither player can force
//! a win — and the Draw-set is a genuinely new degree of freedom to test against
//! the Gold quadric `{Q=0}` (see `OPEN.md`, the Tier-2 open question).
//!
//! Four layers, re-exported flat so every public path is unchanged:
//!
//!   * [`catalogue`] — [`LoopyWinner`], [`LoopyPartizanOutcome`],
//!     [`PartizanOutcome`], and the [`LoopyValue`] stopper catalogue
//!     (on/off/over/under/dud/±/tis/tisn/∗/0/`s&t` with outcome/neg/partial
//!     order/partial sum).
//!   * [`graph`] — [`LoopyGraph`], the computable wrapper over
//!     [`kernel::outcomes`](crate::games::outcomes) (Win / Loss / Draw retrograde
//!     analysis).
//!   * [`partizan`] — [`LoopyPartizanGraph`]: the two-sided Left/Right retrograde
//!     solver returning exact [`LoopyPartizanOutcome`] pairs, projecting to the
//!     classical five-class [`PartizanOutcome`] only when honest.
//!   * [`nim_values`] — [`LoopyNimber`], [`LoopyNimCertificate`],
//!     [`loopy_nim_values`], and [`loopy_nim_values_certified`]: impartial loopy
//!     nim-values with certificates (including the checked recovery condition).
//!   * [`research`] — [`loopy_decision_sets`] and [`loopy_quadric_probe`]: the
//!     Loss-set / Draw-set research instrument.
//!
//! Deliberately **out of scope** here: [`Game`](crate::games::Game) stays an acyclic
//! `Arc` tree (it cannot represent cycles, by construction), and
//! [`thermography`](crate::games::thermography) stays finite-game-only — loopy games
//! never freeze to a number, so classical temperature does not apply. The sidling
//! support is finite and certified: over-budget or non-canonical fixed-point
//! systems return `None` rather than pretending to be full loopy-game equality.

pub mod catalogue;
pub mod graph;
pub mod nim_values;
pub mod partizan;
pub mod research;

pub use catalogue::*;
pub use graph::*;
pub use nim_values::*;
pub use partizan::*;
pub use research::*;
