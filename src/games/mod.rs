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
//!   * [`thermography`] — temperature theory: stops, cooling, and the
//!                        thermograph (mean value + temperature) of a short game.
//!   * [`hackenbush`]   — red/blue/green Hackenbush: the one structure whose value
//!                        reads out as surreal (blue–red), nimber (green), or a
//!                        general partizan game (mixed) — the unifier.

pub mod atomic_weight;
pub mod coin_turning;
pub mod game_exterior;
pub mod grundy;
pub mod hackenbush;
pub mod kernel;
pub mod misere;
pub mod partizan;
pub mod thermography;

pub use atomic_weight::*;
pub use coin_turning::*;
pub use game_exterior::*;
pub use grundy::*;
pub use hackenbush::*;
pub use kernel::*;
pub use misere::*;
pub use partizan::*;
pub use thermography::*;
