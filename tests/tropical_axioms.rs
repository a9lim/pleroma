//! Property-based **semiring** axioms for the tropical worlds, run in both the
//! `(max, +)` and `(min, +)` conventions.
//!
//! The mirror of [`scalar_axioms`](../scalar_axioms.rs) for the sibling
//! [`Semiring`] structure: where that suite fuzzes the commutative-*ring* laws
//! the Clifford engine assumes, this one fuzzes the strictly weaker semiring laws
//! (no additive inverse — `⊕` is idempotent instead), which is exactly why
//! `Tropical` is not a `Scalar`.

use pleroma::scalar::{MaxPlus, MinPlus, Rational, Semiring, Tropical, TropicalConvention};
use proptest::prelude::*;

mod common;
use common::proptest_config;

// Default CI/local runs are smoke-sized. Set `PLEROMA_PROPTEST_CASES=N` for
// deeper semiring fuzzing.
const CASES: u32 = 2;

/// Every commutative-semiring law, checked on one triple `(a, b, c)`.
fn semiring_axioms<T: Semiring>(a: &T, b: &T, c: &T) {
    // (R, ⊕) is a commutative, *idempotent* monoid (no inverse — that is the
    // whole point of a semiring vs a ring).
    assert!(a.add(b).add(c) == a.add(&b.add(c)), "⊕ associative");
    assert!(a.add(b) == b.add(a), "⊕ commutative");
    assert!(a.add(&T::zero()) == *a, "0 is the ⊕-identity");
    assert!(a.add(a) == *a, "⊕ idempotent");

    // (R, ⊗) is a commutative monoid
    assert!(a.mul(b).mul(c) == a.mul(&b.mul(c)), "⊗ associative");
    assert!(a.mul(b) == b.mul(a), "⊗ commutative");
    assert!(a.mul(&T::one()) == *a, "1 is the ⊗-identity");

    // 0 absorbs under ⊗ (a semiring axiom; here ∞ ⊗ a = ∞)
    assert!(a.mul(&T::zero()).is_zero(), "0 absorbs under ⊗ (right)");
    assert!(T::zero().mul(a).is_zero(), "0 absorbs under ⊗ (left)");

    // ⊗ distributes over ⊕, both sides
    assert!(
        a.mul(&b.add(c)) == a.mul(b).add(&a.mul(c)),
        "⊗ over ⊕ (left)"
    );
    assert!(
        a.add(b).mul(c) == a.mul(c).add(&b.mul(c)),
        "⊗ over ⊕ (right)"
    );
}

/// Small tropical values: mostly finite small rationals, with the ⊕-identity
/// (`∞`) weighted in so the absorbing/identity laws actually get exercised.
fn tropicals<C: TropicalConvention>() -> impl Strategy<Value = Tropical<C>> {
    prop_oneof![
        1 => Just(Tropical::<C>::infinity()),
        9 => (-40i128..40, 1i128..40).prop_map(|(n, d)| Tropical::<C>::finite(Rational::new(n, d))),
    ]
}

macro_rules! semiring_suite {
    ($name:ident, $ty:ty, $strat:expr) => {
        proptest! {
            #![proptest_config(proptest_config(CASES))]
            #[test]
            fn $name(a in $strat, b in $strat, c in $strat) {
                semiring_axioms::<$ty>(&a, &b, &c);
            }
        }
    };
}

semiring_suite!(
    max_plus_semiring_axioms,
    Tropical<MaxPlus>,
    tropicals::<MaxPlus>()
);
semiring_suite!(
    min_plus_semiring_axioms,
    Tropical<MinPlus>,
    tropicals::<MinPlus>()
);
