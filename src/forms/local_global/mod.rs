//! Local-global quadratic-form machinery.
//!
//! The scalar pillar names the local and global coefficient worlds. This
//! submodule keeps the corresponding form-theoretic local-global layer together:
//! Hilbert symbols and Hasse-Minkowski over `Q`, the adelic rational facade, and
//! the odd- and characteristic-2 function-field mirrors. The parent `forms`
//! module re-exports these modules and their public items flat for the existing
//! API.

pub mod adelic;
pub mod function_field;
pub mod function_field_char2;
pub mod global_field;
pub mod padic;

pub use adelic::*;
pub use function_field::*;
pub use function_field_char2::*;
pub use global_field::*;
pub use padic::*;
