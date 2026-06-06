//! Combinatorial game theory: the second column of the project, mostly
//! independent of the scalar/Clifford stack.
//!
//!   * [`coin_turning`] — nim-multiplication as Conway's Turning-Corners mex
//!                        recurrence, general 1-D coin-turning, and the 2-D
//!                        Tartan product.
//!   * [`kernel`]       — normal-play Win/Loss/Draw outcomes of a finite game
//!                        graph (retrograde analysis); P-positions = Loss.
//!   * [`misere`]       — misère-play outcomes, indistinguishability quotients,
//!                        and octal games.
//!   * [`partizan`]     — short partizan games plus the exterior algebra of the
//!                        game group (the Clifford-adjacent structure that lives
//!                        on all of game-world, not just the numbers).

pub mod coin_turning;
pub mod kernel;
pub mod misere;
pub mod partizan;

pub use coin_turning::*;
pub use kernel::*;
pub use misere::*;
pub use partizan::*;
