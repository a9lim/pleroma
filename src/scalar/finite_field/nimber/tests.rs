use super::galois::ORDER_FACTORS;
use super::*;

#[test]
fn add_is_xor_and_self_inverse() {
    for a in 0u128..64 {
        for b in 0u128..64 {
            assert_eq!(nim_add(a, b), a ^ b);
        }
        assert_eq!(nim_add(a, a), 0); // own inverse
    }
}

#[test]
fn known_small_products() {
    // F_4 = {0,1,2,3}: 2 is a generator with 2^2 = 3.
    assert_eq!(nim_mul(2, 2), 3);
    assert_eq!(nim_mul(2, 3), 1);
    assert_eq!(nim_mul(3, 3), 2);
    // Fermat powers: 4 (x) 4 = 6, distinct powers 4 (x) 2 = 8 (ordinary).
    assert_eq!(nim_mul(4, 4), 6);
    assert_eq!(nim_mul(2, 4), 8);
    assert_eq!(nim_mul(16, 16), 24); // F_2 (x) F_2 = (3/2)*16
                                     // identity / zero
    assert_eq!(nim_mul(1, 37), 37);
    assert_eq!(nim_mul(0, 37), 0);
}

#[test]
fn field_axioms_over_f16() {
    // {0..15} = F_16 should be a field under nim ops.
    for a in 0u128..16 {
        for b in 0u128..16 {
            // commutativity
            assert_eq!(nim_mul(a, b), nim_mul(b, a));
            // closure within F_16
            assert!(nim_mul(a, b) < 16, "{a} (x) {b} left F_16");
            for c in 0u128..16 {
                // associativity
                assert_eq!(
                    nim_mul(nim_mul(a, b), c),
                    nim_mul(a, nim_mul(b, c)),
                    "assoc {a} {b} {c}"
                );
                // distributivity over XOR
                assert_eq!(
                    nim_mul(a, b ^ c),
                    nim_mul(a, b) ^ nim_mul(a, c),
                    "distrib {a} {b} {c}"
                );
            }
        }
    }
}

#[test]
fn every_nonzero_has_inverse_in_f16() {
    for a in 1u128..16 {
        let inv = (1u128..16).find(|&x| nim_mul(a, x) == 1);
        assert!(inv.is_some(), "no inverse for {a} in F_16");
    }
}

#[test]
fn inverse_round_trips() {
    // x ⊗ x^{-1} = 1 for a spread of nonzero nimbers across several fields.
    for x in [
        1u128,
        2,
        3,
        4,
        5,
        15,
        16,
        255,
        256,
        65535,
        65536,
        1_000_000,
        u128::MAX,
    ] {
        let inv = nim_inv(x).unwrap();
        assert_eq!(nim_mul(x, inv), 1, "inverse of {x}");
    }
    assert_eq!(nim_inv(0), None);
    // matches the brute-forced inverses inside F_16
    for x in 1u128..16 {
        let brute = (1u128..16).find(|&y| nim_mul(x, y) == 1).unwrap();
        assert_eq!(nim_inv(x).unwrap(), brute, "F_16 inverse of {x}");
    }
}

#[test]
fn sqrt_is_inverse_frobenius() {
    // √ is the unique inverse of squaring in char 2: (√x)² = x and √(x²) = x.
    for x in [
        0u128,
        1,
        2,
        3,
        5,
        15,
        16,
        255,
        256,
        65535,
        65536,
        1 << 40,
        u128::MAX,
    ] {
        assert_eq!(nim_square(nim_sqrt(x)), x, "(√{x})² ≠ {x}");
        assert_eq!(nim_sqrt(nim_square(x)), x, "√({x}²) ≠ {x}");
    }
    // a square root stays inside the subfield its argument lives in (F_16).
    for x in 0u128..16 {
        assert!(nim_sqrt(x) < 16, "√{x} left F_16");
    }
}

#[test]
fn trace_is_in_f2_and_is_additive() {
    // Tr lands in {0,1} and is F₂-linear (additive) over F_16.
    for x in 0u128..16 {
        assert!(nim_trace(x, 4) <= 1);
        for y in 0u128..16 {
            assert_eq!(nim_trace(x ^ y, 4), nim_trace(x, 4) ^ nim_trace(y, 4));
        }
    }
}

#[test]
fn artin_schreier_solvable_iff_trace_zero() {
    // The unification: y²+y=c is solvable exactly when Tr(c)=0, and the solver
    // returns a genuine root when it is. Checked exhaustively on F_16.
    let m = 4;
    for c in 0u128..16 {
        let solvable = nim_trace(c, m) == 0;
        assert_eq!(nim_is_artin_schreier_solvable(c, m), solvable);
        match nim_solve_artin_schreier(c, m) {
            Some(y) => {
                assert!(solvable, "solver returned a root for trace-1 c={c}");
                assert_eq!(nim_square(y) ^ y, c, "y²+y ≠ c for c={c}");
                assert!(y < 16, "root left F_16");
            }
            None => assert!(!solvable, "solver gave up on trace-0 c={c}"),
        }
    }
    // Exactly half of F_16 is trace-zero (the image is a hyperplane).
    let solvable_count = (0u128..16).filter(|&c| nim_trace(c, m) == 0).count();
    assert_eq!(solvable_count, 8);
}

#[test]
fn artin_schreier_over_f256() {
    // larger field: solver agrees with the trace obstruction on a sample.
    let m = 8;
    for c in (0u128..256).step_by(7) {
        let y = nim_solve_artin_schreier(c, m);
        assert_eq!(y.is_some(), nim_trace(c, m) == 0, "c={c}");
        if let Some(y) = y {
            assert_eq!(nim_square(y) ^ y, c);
        }
    }
}

#[test]
fn associativity_spot_check_large() {
    // a few larger triples to exercise multi-Fermat recursion
    for &(a, b, c) in &[(255u128, 256, 257), (1000, 999, 7), (65535, 65536, 3)] {
        assert_eq!(nim_mul(nim_mul(a, b), c), nim_mul(a, nim_mul(b, c)));
    }
}

// ----- finite-field analysis toolkit -----

fn brute_order(x: u128) -> u128 {
    let mut k = 1u128;
    let mut cur = x;
    while cur != 1 {
        cur = nim_mul(cur, x);
        k += 1;
    }
    k
}

/// Evaluate `Σ poly[i]·x^{⊗i}` in nim arithmetic (poly over F₂).
fn eval_poly_f2(poly: &[u8], x: u128) -> u128 {
    let mut acc = 0u128;
    let mut xpow = 1u128;
    for &c in poly {
        if c == 1 {
            acc ^= xpow;
        }
        xpow = nim_mul(xpow, x);
    }
    acc
}

#[test]
fn order_factors_are_2_128_minus_1() {
    let mut prod = 1u128;
    for &p in &ORDER_FACTORS {
        prod = prod.checked_mul(p).expect("ORDER_FACTORS overflow");
    }
    assert_eq!(prod, u128::MAX); // 2^128 − 1, squarefree
}

#[test]
fn degree_is_smallest_containing_subfield() {
    assert_eq!(nim_degree(0), 1);
    assert_eq!(nim_degree(1), 1);
    assert_eq!(nim_degree(2), 2); // F_4 \ F_2
    assert_eq!(nim_degree(3), 2);
    for x in 4u128..16 {
        assert_eq!(nim_degree(x), 4, "deg {x}"); // F_16 \ F_4
    }
    assert_eq!(nim_degree(16), 8); // F_256 \ F_16
}

#[test]
fn conjugates_and_min_poly() {
    for x in 0u128..16 {
        let conj = nim_conjugates(x);
        assert_eq!(conj.len() as u32, nim_degree(x));
        let mut s = conj.clone();
        s.sort_unstable();
        s.dedup();
        assert_eq!(s.len(), conj.len(), "conjugates of {x} not distinct");

        let mp = nim_min_poly(x);
        assert_eq!(mp.len() as u32, nim_degree(x) + 1);
        assert_eq!(*mp.last().unwrap(), 1, "min poly of {x} not monic");
        assert!(mp.iter().all(|&c| c <= 1));
        for &c in &conj {
            assert_eq!(eval_poly_f2(&mp, c), 0, "min poly of {x}: root {c}");
        }
    }
}

#[test]
fn relative_trace_and_norm() {
    // the e=1 relative trace is the existing F₂ trace
    for x in 0u128..16 {
        assert_eq!(nim_relative_trace(x, 4, 1), nim_trace(x, 4));
    }
    // relative trace/norm land in the target subfield F_16
    for x in 0u128..256 {
        assert!(nim_relative_trace(x, 8, 4) < 16);
        assert!(nim_relative_norm(x, 8, 4) < 16);
    }
    // norm to the prime field is 1 for every nonzero element
    for x in 1u128..16 {
        assert_eq!(nim_relative_norm(x, 4, 1), 1);
    }
    // the relative norm is multiplicative
    for a in 1u128..16 {
        for b in 1u128..16 {
            assert_eq!(
                nim_relative_norm(nim_mul(a, b), 4, 2),
                nim_mul(nim_relative_norm(a, 4, 2), nim_relative_norm(b, 4, 2)),
                "norm({a}⊗{b})"
            );
        }
    }
}

#[test]
fn order_matches_brute_force_in_subfields() {
    for x in 1u128..16 {
        assert_eq!(nim_order(x), Some(brute_order(x)), "order of {x}");
    }
    assert_eq!(nim_order(0), None);
    assert_eq!(nim_order(2), Some(3)); // 2 generates F_4*
    for x in 1u128..16 {
        assert!(!nim_is_primitive(x)); // all sit in a proper subfield
    }
}

#[test]
fn discrete_log_round_trips() {
    // ⟨2⟩ = {1,2,3} ⊂ F_4 (order 3)
    assert_eq!(nim_discrete_log(2, 1), Some(0));
    assert_eq!(nim_discrete_log(2, 2), Some(1));
    assert_eq!(nim_discrete_log(2, 3), Some(2));
    assert_eq!(nim_discrete_log(2, 4), None); // 4 ∉ ⟨2⟩

    // a generator of F_256* (order 255 = 3·5·17): exercises Pohlig–Hellman + CRT
    let g = (16u128..256).find(|&g| nim_order(g) == Some(255)).unwrap();
    for e in 0u128..255 {
        assert_eq!(nim_discrete_log(g, nim_pow(g, e)), Some(e), "log_{g}");
    }
    // a non-generator base (order 51)
    let h = nim_pow(g, 5);
    assert_eq!(nim_order(h), Some(51));
    let target = nim_pow(h, 7);
    let e = nim_discrete_log(h, target).unwrap();
    assert!(e < 51 && nim_pow(h, e) == target);
}

#[test]
fn primitive_element_generates_full_group() {
    let g = nim_primitive_element();
    assert!(nim_is_primitive(g));
    assert_eq!(nim_order(g), Some(u128::MAX)); // order 2^128 − 1
}
