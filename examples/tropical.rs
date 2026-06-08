//! The tropical semiring, and thermography as tropical arithmetic. Run with:
//!   cargo run --example tropical
//!
//! Two takeaways:
//!   1. `(max, +)` and `(min, +)` are dual semirings — same `⊗` (ordinary `+`),
//!      mirror-image `⊕` (max vs min) and mirror-image identities (`−∞` vs `+∞`).
//!   2. A short game's thermograph IS a pair of folds in those dual semirings:
//!      the left wall is a `(max, +)` ⊕-fold, the right wall a `(min, +)` one.
//!      `thermograph_via_tropical` names the folds and comes out identical to the
//!      golden `thermograph`.

use pleroma::games::{thermograph, thermograph_via_tropical, Game};
use pleroma::scalar::{MaxPlus, MinPlus, Rational, Tropical};

fn rule(title: &str) {
    println!("\n── {title} ──");
}

fn main() {
    rule("tropical semiring — (max,+) vs (min,+), the dual conventions");
    let (a, b) = (Tropical::<MaxPlus>::int(2), Tropical::<MaxPlus>::int(5));
    println!("  max-plus:  2 ⊕ 5 = {}   2 ⊗ 5 = {}", a.add(&b), a.mul(&b));
    let (c, d) = (Tropical::<MinPlus>::int(2), Tropical::<MinPlus>::int(5));
    println!("  min-plus:  2 ⊕ 5 = {}   2 ⊗ 5 = {}", c.add(&d), c.mul(&d));
    println!(
        "  identities: max-plus 0 = {} (⊕), 1 = {} (⊗);  min-plus 0 = {}, 1 = {}",
        Tropical::<MaxPlus>::zero(),
        Tropical::<MaxPlus>::one(),
        Tropical::<MinPlus>::zero(),
        Tropical::<MinPlus>::one(),
    );
    // ∞ is the ⊕-identity and absorbs under ⊗.
    let inf = Tropical::<MaxPlus>::infinity();
    println!(
        "  max-plus ∞: 7 ⊕ {} = {}   7 ⊗ {} = {}",
        inf,
        Tropical::<MaxPlus>::int(7).add(&inf),
        inf,
        Tropical::<MaxPlus>::int(7).mul(&inf),
    );

    rule("thermography IS tropical — left wall = (max,+) fold, right = (min,+)");
    // {3 | {1|−1}}: a hot game with a hot right option.
    let g = Game::new(vec![Game::integer(3)], vec![Game::switch(1, -1)]);
    let golden = thermograph(&g).expect("thermograph");
    let named = thermograph_via_tropical(&g).expect("thermograph_via_tropical");
    println!("  game        {}", g.display());
    println!(
        "  mast {:?}   temperature {:?}   stops (LS,RS) = ({:?}, {:?})",
        named.mast,
        named.temperature,
        named.left_stop(),
        named.right_stop(),
    );
    // The whole point: the named-tropical thermograph equals the golden one.
    let same_mast = golden.mast == named.mast;
    let same_temp = golden.temperature == named.temperature;
    let same_walls = golden.left_wall.points() == named.left_wall.points()
        && golden.right_wall.points() == named.right_wall.points();
    assert!(
        same_mast && same_temp && same_walls,
        "naming must be faithful"
    );
    println!("  thermograph_via_tropical == thermograph  ✓");

    rule("a switch {a|b}: the two stops fall out of the two dual semirings");
    for (a, b) in [(5i128, 1i128), (2, -2)] {
        let th = thermograph_via_tropical(&Game::switch(a, b)).unwrap();
        // LS = a comes from the (max,+) left fold; RS = b from the (min,+) right fold.
        assert!(th.left_stop() == Rational::int(a) && th.right_stop() == Rational::int(b));
        println!(
            "  {{{a}|{b}}}: LS = {:?} (max-plus),  RS = {:?} (min-plus),  mean {:?}",
            th.left_stop(),
            th.right_stop(),
            th.mast,
        );
    }
    println!();
}
