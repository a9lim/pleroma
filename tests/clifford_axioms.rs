//! Property-based associativity/distributivity of the geometric product, over
//! random metrics and random multivectors, in both characteristic 0 (`Rational`)
//! and characteristic 2 (`Nimber`).
//!
//! The unit suite pins associativity on fixed cases; this fuzzes it. A bug in
//! `geom_product_blades` (sign handling and the polar/quadratic split) shows up
//! as a random associativity failure here, with a shrunk counterexample, instead
//! of as a downstream classifier mystery. The general bilinear `a` term is pinned
//! by fixed unit tests in the engine module.
//!
//! Keep this suite light: the scalar nimber tests pin field arithmetic and
//! high-bit representation paths separately. Here we need enough coefficient
//! variety to exercise the Clifford product law, not hundreds of maximal-width
//! nim products per case.

use ogdoad::clifford::{bits, CliffordAlgebra, Metric, Multivector};
use ogdoad::scalar::{Nimber, Rational, Scalar};
use proptest::prelude::*;
use std::collections::BTreeMap;

mod common;
use common::proptest_config;

const DIM: usize = 3;
const BLADES: usize = 1 << DIM; // 8

// Default CI/local runs are smoke-sized; the explicit sentinel below owns the
// high-bit nimber path. Set `OGDOAD_PROPTEST_CASES=N` for deeper fuzzing.
const CASES: u32 = 2;

/// Build a multivector from coefficients indexed by blade mask `0..2^DIM`.
fn build_mv<S: Scalar>(alg: &CliffordAlgebra<S>, coeffs: &[S]) -> Multivector<S> {
    let mut mv = alg.zero();
    for (mask, c) in coeffs.iter().enumerate() {
        let blade = alg.blade(&bits(mask as u128));
        mv = alg.add(&mv, &alg.scalar_mul(c, &blade));
    }
    mv
}

/// Off-diagonal polar form `b` keyed `(i,j)` with `i<j`, from the three pair
/// values (dim 3).
fn b_map<S: Scalar>(v: [S; 3]) -> BTreeMap<(usize, usize), S> {
    let [v01, v02, v12] = v;
    BTreeMap::from([((0, 1), v01), ((0, 2), v02), ((1, 2), v12)])
}

fn nimber_coeff() -> impl Strategy<Value = u128> {
    prop_oneof![
        // F_16-sized values are enough to exercise nontrivial char-2 products,
        // trace-like cancellations, and q/b independence cheaply. The explicit
        // sentinel test below owns the high-bit representation boundary.
        8 => 0u128..16,
        // A little wider finite fuzz catches table-boundary mistakes.
        2 => any::<u8>().prop_map(u128::from),
    ]
}

fn nimber_mv_coeff() -> impl Strategy<Value = u128> {
    prop_oneof![
        // Sparse multivectors keep each associativity check small while still
        // sampling scalars, vectors, bivectors, and pseudoscalars independently.
        3 => Just(0u128),
        7 => nimber_coeff(),
    ]
}

fn rational_mv_coeff() -> impl Strategy<Value = i128> {
    prop_oneof![3 => Just(0i128), 7 => -3i128..=3]
}

fn check_associative_distributive<S: Scalar>(
    alg: &CliffordAlgebra<S>,
    a: &Multivector<S>,
    b: &Multivector<S>,
    c: &Multivector<S>,
) {
    // (ab)c = a(bc)
    let lhs = alg.mul(&alg.mul(a, b), c);
    let rhs = alg.mul(a, &alg.mul(b, c));
    assert_eq!(lhs, rhs, "geometric product not associative");
    // a(b+c) = ab + ac  and  (a+b)c = ac + bc
    let left = alg.mul(a, &alg.add(b, c));
    let left_expanded = alg.add(&alg.mul(a, b), &alg.mul(a, c));
    assert_eq!(left, left_expanded, "left distributivity");
    let right = alg.mul(&alg.add(a, b), c);
    let right_expanded = alg.add(&alg.mul(a, c), &alg.mul(b, c));
    assert_eq!(right, right_expanded, "right distributivity");
}

proptest! {
    #![proptest_config(proptest_config(CASES))]

    /// Characteristic 2, with an independent quadratic form `q` and polar form
    /// `b` (the nimber-backend point: `q ≠ b` must stay faithful).
    #[test]
    fn nimber_geometric_product_is_a_ring(
        q in prop::array::uniform3(nimber_coeff()),
        bvals in prop::array::uniform3(nimber_coeff()),
        ca in prop::array::uniform::<_, BLADES>(nimber_mv_coeff()),
        cb in prop::array::uniform::<_, BLADES>(nimber_mv_coeff()),
        cc in prop::array::uniform::<_, BLADES>(nimber_mv_coeff()),
    ) {
        let metric = Metric::new(
            q.iter().map(|&x| Nimber(x)).collect(),
            b_map(bvals.map(Nimber)),
        );
        let alg = CliffordAlgebra::new(DIM, metric);
        let mk = |c: [u128; BLADES]| build_mv(&alg, &c.map(Nimber));
        check_associative_distributive(&alg, &mk(ca), &mk(cb), &mk(cc));
    }

    /// Characteristic 0, diagonal metric, small rational coefficients.
    #[test]
    fn rational_geometric_product_is_a_ring(
        q in prop::array::uniform3(-3i128..=3),
        ca in prop::array::uniform::<_, BLADES>(rational_mv_coeff()),
        cb in prop::array::uniform::<_, BLADES>(rational_mv_coeff()),
        cc in prop::array::uniform::<_, BLADES>(rational_mv_coeff()),
    ) {
        let metric = Metric::diagonal(q.iter().map(|&x| Rational::int(x)).collect());
        let alg = CliffordAlgebra::new(DIM, metric);
        let mk = |c: [i128; BLADES]| build_mv(&alg, &c.map(Rational::int));
        check_associative_distributive(&alg, &mk(ca), &mk(cb), &mk(cc));
    }
}

#[test]
fn nimber_geometric_product_sentinel_case() {
    let metric = Metric::new(
        vec![
            Nimber(1u128 << 32),
            Nimber(1u128 << 64),
            Nimber((1u128 << 127) ^ 1),
        ],
        b_map([
            Nimber(1u128 << 96),
            Nimber((1u128 << 127) ^ 3),
            Nimber((1u128 << 96) ^ 255),
        ]),
    );
    let alg = CliffordAlgebra::new(DIM, metric);
    let a = build_mv(
        &alg,
        &[
            Nimber(0),
            Nimber(1),
            Nimber(1u128 << 32),
            Nimber(3),
            Nimber(1u128 << 64),
            Nimber(5),
            Nimber(1u128 << 96),
            Nimber(7),
        ],
    );
    let b = build_mv(
        &alg,
        &[
            Nimber(1),
            Nimber(2),
            Nimber(3),
            Nimber(1u128 << 127),
            Nimber(5),
            Nimber(8),
            Nimber(13),
            Nimber((1u128 << 127) ^ 5),
        ],
    );
    let c = build_mv(
        &alg,
        &[
            Nimber((1u128 << 127) ^ 7),
            Nimber(0),
            Nimber(1u128 << 96),
            Nimber(11),
            Nimber(1u128 << 64),
            Nimber(17),
            Nimber(1u128 << 32),
            Nimber(23),
        ],
    );
    check_associative_distributive(&alg, &a, &b, &c);
}
