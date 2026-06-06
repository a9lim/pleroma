//! Combinatorial game theory: the second column of the project, mostly
//! independent of the scalar/Clifford stack.
//!
//!   * [`coin_turning`] — nim-multiplication as Conway's Turning-Corners mex
//!                        recurrence, general 1-D coin-turning, and the 2-D
//!                        Tartan product.
//!   * [`grundy`]       — Sprague–Grundy values of any finite impartial game
//!                        (the normal-play impartial center; P-position ⟺ g = 0).
//!   * [`kernel`]       — normal-play Win/Loss/Draw outcomes of a finite game
//!                        graph (retrograde analysis); P-positions = Loss.
//!   * [`misere`]       — misère-play outcomes, indistinguishability quotients,
//!                        and octal games.
//!   * [`partizan`]     — short partizan games (sum, order, canonical form, the
//!                        surreal-value bridge) plus the exterior algebra of the
//!                        game group (the Clifford-adjacent structure that lives
//!                        on all of game-world, not just the numbers).

pub mod coin_turning;
pub mod grundy;
pub mod kernel;
pub mod misere;
pub mod partizan;

pub use coin_turning::*;
pub use grundy::*;
pub use kernel::*;
pub use misere::*;
pub use partizan::*;
