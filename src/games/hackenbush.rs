//! Hackenbush — the unifier.
//!
//! A Hackenbush position is a graph of coloured edges standing on the *ground*
//! (vertex `0`). Players alternately delete an edge of their colour (Left:
//! **blue**, Right: **red**, either player: **green**); any edge no longer
//! connected to the ground falls off. Last player to move wins (normal play).
//!
//! The single evaluator [`Hackenbush::to_game`] — build the partizan game by the
//! move-and-prune recursion — reads out as **all three** of pleroma's game-value
//! worlds at once, which is the whole point:
//!
//! | position             | value world        | bridge                       |
//! |----------------------|--------------------|------------------------------|
//! | blue / red only      | surreal **number** | [`Hackenbush::value`]        |
//! | blue–red string      | dyadic surreal     | = its **sign expansion**     |
//! | green only           | **nimber** (Nim)   | [`Hackenbush::grundy`]       |
//! | mixed                | general partizan   | the `Game` itself            |
//!
//! A blue–red *string* is exactly an [ordinal sum](crate::games::Game::ordinal_sum)
//! of single edges, and Berlekamp's rule says its value's
//! [sign expansion](crate::scalar::Surreal::sign_expansion) is the colour
//! sequence read from the ground up (blue `+`, red `−`). A green *string* of `n`
//! edges is the Nim heap `*n`. Both fall out of the one recursion below.

use crate::games::Game;
use crate::scalar::Surreal;
use std::collections::{BTreeSet, HashSet};

/// An edge colour: who may remove it.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Color {
    /// Left's edge.
    Blue,
    /// Right's edge.
    Red,
    /// Either player's edge (impartial).
    Green,
}

/// A Hackenbush position: undirected coloured edges over the ground vertex `0`.
#[derive(Clone, Debug)]
pub struct Hackenbush {
    edges: Vec<(usize, usize, Color)>,
}

impl Hackenbush {
    /// A position from an explicit edge list `(u, v, colour)`. Vertex `0` is the
    /// ground.
    pub fn new(edges: Vec<(usize, usize, Color)>) -> Hackenbush {
        Hackenbush { edges }
    }

    /// A **stalk** rooted at the ground: `0 — 1 — 2 — …`, edge `i` joining
    /// vertices `i` and `i+1` with colour `colors[i]`.
    pub fn string(colors: &[Color]) -> Hackenbush {
        let edges = colors
            .iter()
            .enumerate()
            .map(|(i, &c)| (i, i + 1, c))
            .collect();
        Hackenbush { edges }
    }

    /// The edges `(u, v, colour)`.
    pub fn edges(&self) -> &[(usize, usize, Color)] {
        &self.edges
    }

    /// The vertices connected to the ground (`0`) through the current edges.
    fn grounded(&self) -> HashSet<usize> {
        let mut reach = HashSet::new();
        reach.insert(0usize);
        let mut changed = true;
        while changed {
            changed = false;
            for &(u, v, _) in &self.edges {
                let (ur, vr) = (reach.contains(&u), reach.contains(&v));
                if ur ^ vr {
                    reach.insert(if ur { v } else { u });
                    changed = true;
                }
            }
        }
        reach
    }

    /// Remove edge `i`, then drop every edge that has fallen off the ground.
    fn remove_edge(&self, i: usize) -> Hackenbush {
        let mut edges = self.edges.clone();
        edges.remove(i);
        let pruned = Hackenbush { edges };
        let grounded = pruned.grounded();
        Hackenbush {
            edges: pruned
                .edges
                .into_iter()
                .filter(|&(u, v, _)| grounded.contains(&u) && grounded.contains(&v))
                .collect(),
        }
    }

    /// The partizan game value, as a [`Game`] — the universal evaluator. Left
    /// options are the blue/green deletions, Right options the red/green ones,
    /// each followed by pruning.
    pub fn to_game(&self) -> Game {
        let mut left = Vec::new();
        let mut right = Vec::new();
        for (i, &(_, _, c)) in self.edges.iter().enumerate() {
            let sub = self.remove_edge(i).to_game();
            match c {
                Color::Blue => left.push(sub),
                Color::Red => right.push(sub),
                Color::Green => {
                    left.push(sub.clone());
                    right.push(sub);
                }
            }
        }
        Game::new(left, right)
    }

    /// The **surreal number** value — `Some` exactly when the position's value is
    /// a number (every blue/red position is). `None` for values carrying an
    /// infinitesimal or switch (green edges, `↑`, `⋆`, …).
    pub fn value(&self) -> Option<Surreal> {
        self.to_game().number_value()
    }

    /// The **Sprague–Grundy (nim) value** — `Some` only for an all-green
    /// (impartial) position, where Hackenbush *is* Nim. `None` if any edge is
    /// blue or red.
    pub fn grundy(&self) -> Option<u128> {
        if self.edges.iter().any(|&(_, _, c)| c != Color::Green) {
            return None;
        }
        Some(self.grundy_green())
    }

    fn grundy_green(&self) -> u128 {
        let reachable: BTreeSet<u128> = (0..self.edges.len())
            .map(|i| self.remove_edge(i).grundy_green())
            .collect();
        let mut m = 0u128;
        while reachable.contains(&m) {
            m += 1;
        }
        m
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scalar::{Rational, Surreal};

    fn blue_red(colors: &[Color]) -> Hackenbush {
        Hackenbush::string(colors)
    }

    /// The canonical Nim-heap game `*n`.
    fn nim_heap(n: u128) -> Game {
        let opts: Vec<Game> = (0..n).map(nim_heap).collect();
        Game::new(opts.clone(), opts)
    }

    #[test]
    fn blue_and_red_strings_are_integers() {
        use Color::*;
        for n in 0u128..5 {
            let blue = Hackenbush::string(&vec![Blue; n as usize]);
            assert_eq!(blue.value(), Some(Surreal::from_int(n as i128)));
            assert!(blue.to_game().eq(&Game::integer(n as i128)));

            let red = Hackenbush::string(&vec![Red; n as usize]);
            assert_eq!(red.value(), Some(Surreal::from_int(-(n as i128))));
        }
    }

    #[test]
    fn green_strings_are_nim_heaps() {
        use Color::*;
        for n in 0u128..6 {
            let g = Hackenbush::string(&vec![Green; n as usize]);
            assert_eq!(g.grundy(), Some(n)); // mex recursion = Nim heap n
            assert!(g.to_game().eq(&nim_heap(n))); // and the game value agrees
            if n >= 1 {
                assert_eq!(g.value(), None); // *n (n≥1) is not a number
            }
        }
    }

    #[test]
    fn blue_red_strings_are_their_sign_expansion() {
        use Color::*;
        // Berlekamp's rule: value's sign expansion = colours ground→top (B=+, R=−).
        let cases: [&[Color]; 6] = [
            &[Blue, Red],            // +−  = 1/2
            &[Blue, Red, Blue],      // +−+ = 3/4
            &[Red, Blue],            // −+  = −1/2
            &[Blue, Blue, Red],      // ++− = 3/2
            &[Blue, Red, Red],       // +−− = 1/4
            &[Red, Blue, Red, Blue], // −+−+ = −5/8
        ];
        for colors in cases {
            let signs: Vec<bool> = colors.iter().map(|&c| c == Blue).collect();
            let expected = Surreal::from_sign_expansion(&signs);
            assert_eq!(
                blue_red(colors).value(),
                Some(expected),
                "colors {:?}",
                colors
            );
        }
    }

    #[test]
    fn the_unifier_one_structure_three_worlds() {
        use Color::*;
        // surreal integer
        assert_eq!(
            Hackenbush::string(&[Blue, Blue, Blue]).value(),
            Some(Surreal::from_int(3))
        );
        // nimber
        assert_eq!(Hackenbush::string(&[Green, Green]).grundy(), Some(2));
        // dyadic surreal via sign expansion
        assert_eq!(
            Hackenbush::string(&[Blue, Red]).value(),
            Some(Surreal::from_rational(Rational::new(1, 2)))
        );
    }

    #[test]
    fn green_cycle_and_mixed() {
        use Color::*;
        // a green triangle hung from the ground by vertex 0. Removing a rim edge
        // leaves a 2-edge path (*2); removing the far edge leaves two separate
        // stalks (*1 ⊕ *1 = *0). So options are {*0, *2} and the triangle is
        // mex{0,2} = *1 — exactly the fusion principle (a green cycle = one edge).
        let triangle = Hackenbush::new(vec![(0, 1, Green), (1, 2, Green), (2, 0, Green)]);
        assert_eq!(triangle.grundy(), Some(1));

        // a blue edge atop a green edge is a partizan infinitesimal: not a number.
        let mixed = Hackenbush::new(vec![(0, 1, Green), (1, 2, Blue)]);
        assert_eq!(mixed.value(), None);
        assert!(mixed.grundy().is_none()); // has a coloured edge
    }
}
