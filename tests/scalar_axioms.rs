//! Property-based commutative-ring axioms, run across every `Scalar` backend.
//!
//! The Clifford engine is generic over `Scalar` and *assumes* a commutative
//! ring; these proptests are the safety net under that assumption. One generic
//! [`ring_axioms`] checker is fed randomized triples from each backend's own
//! strategy, so a regression in any backend's arithmetic surfaces here rather
//! than as a mysterious geometric-product failure.

use pleroma::scalar::{
    Fp, Integer, Nimber, Ordinal, Poly, Rational, RationalFunction, Scalar, Surcomplex, Surreal,
};
use proptest::prelude::*;

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
    // any element of F_{2^64} ⊂ F_{2^128}; spans many nim-subfields
    any::<u64>().prop_map(|x| Nimber(x as u128))
}

fn integers() -> impl Strategy<Value = Integer> {
    (-1000i128..1000).prop_map(Integer)
}

fn rationals() -> impl Strategy<Value = Rational> {
    (-40i128..40, 1i128..40).prop_map(|(n, d)| Rational::new(n, d))
}

fn fp7() -> impl Strategy<Value = Fp<7>> {
    any::<i64>().prop_map(|x| Fp::<7>::new(x as i128))
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
    ($name:ident, $ty:ty, $strat:expr) => {
        proptest! {
            #![proptest_config(ProptestConfig::with_cases(256))]
            #[test]
            fn $name(a in $strat, b in $strat, c in $strat) {
                ring_axioms::<$ty>(&a, &b, &c);
            }
        }
    };
}

axiom_suite!(nimber_ring_axioms, Nimber, nimbers());
axiom_suite!(integer_ring_axioms, Integer, integers());
axiom_suite!(rational_ring_axioms, Rational, rationals());
axiom_suite!(fp7_field_axioms, Fp<7>, fp7());
axiom_suite!(surreal_ring_axioms, Surreal, surreals());
axiom_suite!(surcomplex_ring_axioms, Surcomplex<Surreal>, surcomplexes());
axiom_suite!(
    rational_function_field_axioms,
    RationalFunction<Fp<7>>,
    rational_functions()
);

// --- transfinite ordinal nimbers On₂: a PARTIAL field (nim_mul is `Option`) ---
//
// `Ordinal` is deliberately NOT a `Scalar`: nim-multiplication is partial — `None`
// at ω^ω and above, the staged tower boundary — so it can't ride the `axiom_suite!`
// macro. This bespoke checker fuzzes the transfinite char-2 field laws instead.
// nim-addition is total and always checked. nim-multiplication is partial (the
// prime-power tower reaches every ordinal `< ω^(ω²)` and free combinations beyond,
// returning `None` only at the staged non-scalar-Kummer boundary), so its laws are
// checked where defined: commutativity (and symmetric definedness) always; the full
// commutative-ring laws on the `< ω^ω` segment, which is pinned **totally closed** so
// the check never silently degenerates; and associativity opportunistically where the
// whole triangle is defined past `ω^ω`.

/// True iff every CNF exponent is finite — i.e. the ordinal is `< ω^ω`, the region
/// where nim-multiplication is implemented (the degree-3 cube-root tower).
fn below_omega_omega(o: &Ordinal) -> bool {
    o.terms().iter().all(|(e, _)| e.as_finite().is_some())
}

/// Small ordinal exponents: mostly finite (so products land inside the tower),
/// occasionally infinite (`ω`, `ω+1`) to exercise the `≥ ω^ω` boundary.
fn ordinal_exp() -> impl Strategy<Value = Ordinal> {
    prop_oneof![
        8 => (0u128..6).prop_map(Ordinal::from_u128),
        1 => Just(Ordinal::omega()),
        1 => Just(Ordinal::omega().nim_add(&Ordinal::from_u128(1))),
    ]
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
    #![proptest_config(ProptestConfig::with_cases(256))]
    #[test]
    fn ordinal_partial_field_axioms(a in ordinals(), b in ordinals(), c in ordinals()) {
        ordinal_field_axioms(&a, &b, &c);
    }
}
