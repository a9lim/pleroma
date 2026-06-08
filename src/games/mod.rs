//! Combinatorial game theory: the second column of the project, mostly
//! independent of the scalar/Clifford stack.
//!
//! * [`coin_turning`] — nim-multiplication as Conway's Turning-Corners mex
//!   recurrence, general 1-D coin-turning, and the 2-D Tartan product.
//! * [`grundy`] — Sprague–Grundy values of any finite impartial game (the
//!   normal-play impartial center; P-position ⟺ g = 0).
//! * [`kernel`] — normal-play Win/Loss/Draw outcomes of a finite game graph
//!   (retrograde analysis); P-positions = Loss.
//! * [`loopy`] — loopy (cyclic) games: the canonical stoppers
//!   (on/off/over/under/dud), impartial loopy nim-values, and the
//!   Loss-set/Draw-set quadric research instrument.
//! * [`misere`] — misère-play outcomes, indistinguishability quotients, and
//!   octal games.
//! * [`partizan`] — short partizan games (sum, order, canonical form, the
//!   surreal-value bridge).
//! * [`number_game`] — transfinite number-valued games carried by their surreal
//!   value, without materializing infinite options.
//! * [`nimber_game`] — its char-2 mirror: transfinite nimber-valued (impartial)
//!   games — Nim heaps `⋆α` — carried by their ordinal Grundy value (`No ↔ On₂` at
//!   the games layer).
//! * [`thermography`] — temperature theory: stops, cooling, and the thermograph
//!   (mean value + temperature) of a short game.
//! * [`piecewise`] — the piecewise-linear rational scaffold machinery used by
//!   thermography.
//! * [`hackenbush`] — red/blue/green Hackenbush: the one structure whose value
//!   reads out as surreal (blue–red), nimber (green), or a general partizan game
//!   (mixed) — the unifier.

pub mod atomic_weight;
pub mod coin_turning;
pub mod game_exterior;
pub mod grundy;
pub mod hackenbush;
pub mod kernel;
pub mod loopy;
pub mod misere;
pub mod nimber_game;
pub mod number_game;
pub mod partizan;
pub mod piecewise;
pub mod thermography;
pub mod tropical_thermography;

pub use atomic_weight::*;
pub use coin_turning::*;
pub use game_exterior::*;
pub use grundy::*;
pub use hackenbush::*;
pub use kernel::*;
pub use loopy::*;
pub use misere::*;
pub use nimber_game::*;
pub use number_game::*;
pub use partizan::*;
pub use piecewise::*;
pub use thermography::*;
pub use tropical_thermography::*;
