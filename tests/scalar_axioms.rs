//! Property-based commutative-ring axioms, run across every `Scalar` backend.
//!
//! The Clifford engine is generic over `Scalar` and *assumes* a commutative
//! ring; these proptests are the safety net under that assumption. One generic
//! [`ring_axioms`] checker is fed randomized triples from each backend's own
//! strategy, so a regression in any backend's arithmetic surfaces here rather
//! than as a mysterious geometric-product failure.

use ogdoad::scalar::{
    Fp, Integer, Nimber, Ordinal, Poly, Rational, RationalFunction, Scalar, Surcomplex, Surreal,
};
use proptest::prelude::*;

mod common;
use common::proptest_config;

// Default CI/local runs are smoke-sized; deterministic sentinels below own the
// expensive boundary cases. Set `OGDOAD_PROPTEST_CASES=N` for deeper fuzzing.
const FAST_CASES: u32 = 2;
const HEAVY_CASES: u32 = 1;

/// Every commutative-ring law, checked on one triple `(a, b, c)`.
fn ring_axioms<S: Scalar>(a: &S, b: &S, c: &S) {
    // (R, +) is an abelian group
    assert!(a.add(b).add(c) == a.add(&b.add(c)), "+ associative");
    assert!(a.add(b) == b.add(a), "+ commutative");
    assert!(a.add(&S::zero()) == *a, "0 is the additive identity");
    assert!(a.add(&a.neg()).is_zero(), "−a is the additive inverse");

    // (R, ·) is a commutative monoid
    assert!(a.mul(b).mul(c) == a.mul(&b.mul(c)), "· associative");
    assert!(a.mul(b) == b.mul(a), "· commutative");
    assert!(a.mul(&S::one()) == *a, "1 is the multiplicative identity");

    // distributivity, both sides (· need not be symmetric in the engine sense,
    // but the scalar ring is genuinely commutative)
    assert!(
        a.mul(&b.add(c)) == a.mul(b).add(&a.mul(c)),
        "left distributive"
    );
    assert!(
        a.add(b).mul(c) == a.mul(c).add(&b.mul(c)),
        "right distributive"
    );

    // derived subtraction is consistent with negate-then-add
    assert!(a.sub(b) == a.add(&b.neg()), "a − b = a + (−b)");

    // inverse round-trips wherever it exists
    if let Some(ai) = a.inv() {
        assert!(a.mul(&ai) == S::one(), "a · a⁻¹ = 1");
        assert!(ai.mul(a) == S::one(), "a⁻¹ · a = 1");
    }
}

// --- per-backend element strategies (small, to keep arithmetic exact) ---

fn nimbers() -> impl Strategy<Value = Nimber> {
    prop_oneof![
        // Most ring-law fuzz only needs cheap subfields; the dedicated nimber
        // unit tests pin the Conway multiplication table and field operations,
        // while the sentinel test below keeps high-bit representation paths alive.
        8 => 0u128..256,
        // Wider, still cheap finite fuzz catches byte-boundary mistakes.
        2 => any::<u16>().prop_map(u128::from),
    ]
    .prop_map(Nimber)
}

fn integers() -> impl Strategy<Value = Integer> {
    (-1000i128..1000).prop_map(Integer)
}

fn rationals() -> impl Strategy<Value = Rational> {
    (-40i128..40, 1i128..40).prop_map(|(n, d)| Rational::new(n, d))
}

fn fp7() -> impl Strategy<Value = Fp<7>> {
    any::<i128>().prop_map(Fp::<7>::new)
}

/// Small surreals: a handful of monomials `ω^e · (p/q)` with `e ∈ [−2,2]`.
fn surreals() -> impl Strategy<Value = Surreal> {
    prop::collection::vec((-2i128..=2, -4i128..=4, 1i128..=4), 0..3).prop_map(|terms| {
        terms.into_iter().fold(Surreal::zero(), |acc, (e, p, q)| {
            acc.add(&Surreal::monomial(
                Surreal::from_int(e),
                Rational::new(p, q),
            ))
        })
    })
}

fn surcomplexes() -> impl Strategy<Value = Surcomplex<Surreal>> {
    (surreals(), surreals()).prop_map(|(re, im)| Surcomplex::new(re, im))
}

/// Small rational functions over `F_7`: `num/den` with `num, den` of degree < 3,
/// the denominator forced nonzero. `F_q(t)` is exact, so — unlike the local
/// precision models — it belongs in this exact-ring fuzz.
fn rational_functions() -> impl Strategy<Value = RationalFunction<Fp<7>>> {
    let coeffs = || prop::collection::vec((0i128..7).prop_map(Fp::<7>::new), 0..3);
    (coeffs(), coeffs()).prop_map(|(num, den)| {
        let den = if Poly::new(den.clone()).is_zero() {
            vec![Fp::<7>::new(1)]
        } else {
            den
        };
        RationalFunction::new(num, den)
    })
}

macro_rules! axiom_suite {
    ($name:ident, $ty:ty, $strat:expr, $cases:expr) => {
        proptest! {
            #![proptest_config(proptest_config($cases))]
            #[test]
            fn $name(a in $strat, b in $strat, c in $strat) {
                ring_axioms::<$ty>(&a, &b, &c);
            }
        }
    };
}

axiom_suite!(nimber_ring_axioms, Nimber, nimbers(), FAST_CASES);
axiom_suite!(integer_ring_axioms, Integer, integers(), FAST_CASES);
axiom_suite!(rational_ring_axioms, Rational, rationals(), FAST_CASES);
axiom_suite!(fp7_field_axioms, Fp<7>, fp7(), FAST_CASES);
axiom_suite!(surreal_ring_axioms, Surreal, surreals(), HEAVY_CASES);
axiom_suite!(
    surcomplex_ring_axioms,
    Surcomplex<Surreal>,
    surcomplexes(),
    HEAVY_CASES
);
axiom_suite!(
    rational_function_field_axioms,
    RationalFunction<Fp<7>>,
    rational_functions(),
    HEAVY_CASES
);

#[test]
fn nimber_ring_axioms_on_representation_sentinels() {
    let triples = [
        (0, 1, 1u128 << 127),
        (1u128 << 32, 1u128 << 64, 1u128 << 96),
        ((1u128 << 127) ^ 1, (1u128 << 96) ^ 255, (1u128 << 64) ^ 17),
    ];
    for (a, b, c) in triples {
        ring_axioms::<Nimber>(&Nimber(a), &Nimber(b), &Nimber(c));
    }
}

// --- transfinite ordinal nimbers On₂: Scalar on the verified tower, checked partial field ---
//
// `Ordinal` implements `Scalar` with panic-on-escape multiplication for the
// Clifford engine, but the non-panicking mathematical surface is still
// `nim_mul -> Option`. This bespoke checker fuzzes that checked surface:
// nim-addition is total and always checked; nim-multiplication is checked where
// defined, with full commutative-ring laws on the `< ω^ω` segment and
// opportunistic associativity past it.

/// True iff every CNF exponent is finite — i.e. the ordinal is `< ω^ω`, the region
/// where nim-multiplication is implemented (the degree-3 cube-root tower).
fn below_omega_omega(o: &Ordinal) -> bool {
    o.terms().iter().all(|(e, _)| e.as_finite().is_some())
}

/// Small finite ordinal exponents, so random products stay in the closed
/// `< ω^ω` tower. The explicit sentinel test below owns the `≥ ω^ω` boundary.
fn ordinal_exp() -> impl Strategy<Value = Ordinal> {
    (0u128..6).prop_map(Ordinal::from_u128)
}

/// A small ordinal: the nim-sum (XOR) of up to three monomials `ω^exp · coeff`,
/// coefficients in `F_4` (`1..=3`).
fn ordinals() -> impl Strategy<Value = Ordinal> {
    prop::collection::vec((ordinal_exp(), 1u128..4), 0..3).prop_map(|terms| {
        terms.into_iter().fold(Ordinal::zero(), |acc, (e, c)| {
            acc.nim_add(&Ordinal::monomial(e, c))
        })
    })
}

/// The transfinite char-2 field laws on one triple, partiality handled explicitly.
fn ordinal_field_axioms(a: &Ordinal, b: &Ordinal, c: &Ordinal) {
    // (On₂, ⊕) is an abelian group of characteristic 2 — total, always checked.
    assert!(
        a.nim_add(b).nim_add(c) == a.nim_add(&b.nim_add(c)),
        "⊕ associative"
    );
    assert!(a.nim_add(b) == b.nim_add(a), "⊕ commutative");
    assert!(a.nim_add(&Ordinal::zero()) == *a, "0 is the ⊕-identity");
    assert!(a.nim_add(a).is_zero(), "α ⊕ α = 0 (char 2)");

    // ⊗ is commutative and its definedness is symmetric — checkable on every triple
    // regardless of the boundary (both sides `None`, or both equal `Some`).
    for (x, y) in [(a, b), (b, c), (a, c)] {
        assert_eq!(
            x.nim_mul(y),
            y.nim_mul(x),
            "⊗ commutative / symmetric domain"
        );
    }

    // The `< ω^ω` segment is **totally closed**: there every product is defined and
    // the full commutative-ring laws hold. (Pinned, so the suite never degenerates
    // into vacuously skipping multiplication.)
    if below_omega_omega(a) && below_omega_omega(b) && below_omega_omega(c) {
        let one = Ordinal::from_u128(1);
        let ab = a.nim_mul(b).expect("< ω^ω is closed under ⊗");
        assert!(
            ab.nim_mul(c).unwrap() == a.nim_mul(&b.nim_mul(c).unwrap()).unwrap(),
            "⊗ associative"
        );
        assert!(a.nim_mul(&one).unwrap() == *a, "1 is the ⊗-identity");
        assert!(
            a.nim_mul(&Ordinal::zero()).unwrap().is_zero(),
            "0 absorbs under ⊗"
        );
        assert!(
            a.nim_mul(&b.nim_add(c)).unwrap() == ab.nim_add(&a.nim_mul(c).unwrap()),
            "⊗ distributes over ⊕"
        );
    }

    // Past `ω^ω` the engine is partial; where a whole associativity triangle is
    // defined, the law must still hold.
    if let (Some(ab), Some(bc)) = (a.nim_mul(b), b.nim_mul(c)) {
        if let (Some(l), Some(r)) = (ab.nim_mul(c), a.nim_mul(&bc)) {
            assert_eq!(l, r, "⊗ associative where defined past ω^ω");
        }
    }
}

proptest! {
    #![proptest_config(proptest_config(HEAVY_CASES))]
    #[test]
    fn ordinal_partial_field_axioms(a in ordinals(), b in ordinals(), c in ordinals()) {
        ordinal_field_axioms(&a, &b, &c);
    }
}

#[test]
fn ordinal_partial_field_axioms_on_boundary_sentinels() {
    let omega = Ordinal::omega();
    let omega_squared = Ordinal::monomial(Ordinal::from_u128(2), 1);
    let omega_to_omega = Ordinal::monomial(omega.clone(), 1);
    let triples = [
        (
            Ordinal::from_u128(0),
            Ordinal::from_u128(1),
            Ordinal::from_u128(3),
        ),
        (
            omega.clone(),
            omega.nim_add(&Ordinal::from_u128(1)),
            omega_squared,
        ),
        (
            omega_to_omega.clone(),
            omega_to_omega.nim_add(&Ordinal::from_u128(1)),
            Ordinal::monomial(omega.nim_add(&Ordinal::from_u128(1)), 1),
        ),
    ];
    for (a, b, c) in triples {
        ordinal_field_axioms(&a, &b, &c);
    }
}

// --- surreal-eq-cost pin ---
//
// `PartialEq for Surreal` was switched from value-based comparison (which
// allocates a subtraction) to a structural term-vector walk, justified by CNF
// uniqueness (see the `PartialEq` impl's doc comment).  This proptest is the
// permanent pin: for all random canonical `Surreal`s, structural equality
// agrees with value-based comparison.

proptest! {
    #![proptest_config(proptest_config(HEAVY_CASES))]
    #[test]
    fn surreal_structural_eq_matches_value_eq(a in surreals(), b in surreals()) {
        let structural = a == b;
        let value_based = a.cmp(&b) == std::cmp::Ordering::Equal;
        prop_assert_eq!(
            structural,
            value_based,
            "structural PartialEq and value cmp disagree for {:?} vs {:?}",
            a,
            b
        );
    }
}
