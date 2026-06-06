//! Nimbers: the ordinals under nim-addition and nim-multiplication, Conway's
//! Field On_2 of characteristic 2. Restricted here to `u64`, i.e. nimbers
//! below 2^64 — which *is* exactly the finite nim-field F_{2^64}, and contains
//! every smaller F_{2^{2^k}} (k <= 5: F_4, F_16, F_256, ... F_{2^32}).
//!
//! nim-addition is XOR (trivially correct).
//!
//! nim-multiplication is built from three facts, each individually checkable:
//!   1. distributes over XOR;
//!   2. 2^i = nim-product of the Fermat 2-powers F_n = 2^(2^n) for n in the
//!      binary expansion of i  (distinct Fermat powers nim-multiply to their
//!      ordinary product);
//!   3. F_n (x) F_n = (3/2) F_n = F_n XOR (F_n >> 1).
//! A Fermat power appearing in both factors gets squared via (3), which emits
//! a bit at position 2^n - 1 < 2^n, so every squaring strictly lowers the
//! largest Fermat index in play — that is what makes the recursion terminate.
//!
//! CORRECTNESS STATUS: verified by the `tests` module below (known F_4 / F_16
//! products + field axioms over a range) once `cargo test` can run.

use crate::scalar::Scalar;
use std::cell::RefCell;
use std::collections::HashMap;

#[inline]
pub fn nim_add(a: u64, b: u64) -> u64 {
    a ^ b
}

thread_local! {
    // memo for 2^i (x) 2^j, keyed (min(i,j), max(i,j)); bounded 64x64.
    static POW2_MEMO: RefCell<HashMap<(u32, u32), u64>> = RefCell::new(HashMap::new());
}

/// 2^i (x) 2^j.
fn nim_mul_pow2(i: u32, j: u32) -> u64 {
    let key = if i <= j { (i, j) } else { (j, i) };
    if let Some(v) = POW2_MEMO.with(|m| m.borrow().get(&key).copied()) {
        return v;
    }

    // Fermat indices that appear once (clean product) vs twice (must square).
    let single = i ^ j;
    let common = i & j;

    // distinct Fermat powers nim-multiply to their ordinary product = 2^single
    let clean: u64 = 1u64 << single;

    let result = if common == 0 {
        clean
    } else {
        // fold the squared Fermat powers together, then multiply by `clean`
        let mut squared: u64 = 1; // nim multiplicative identity
        let mut c = common;
        while c != 0 {
            let n = c.trailing_zeros();
            c &= c - 1;
            let f = 1u64 << (1u64 << n); // F_n = 2^(2^n)
            let factor = f ^ (f >> 1); // F_n (x) F_n = (3/2) F_n
            squared = nim_mul(squared, factor);
        }
        nim_mul(clean, squared)
    };

    POW2_MEMO.with(|m| m.borrow_mut().insert(key, result));
    result
}

/// nim-multiplication, by distributing over the bits of both arguments.
pub fn nim_mul(a: u64, b: u64) -> u64 {
    let mut acc = 0u64;
    let mut aa = a;
    while aa != 0 {
        let i = aa.trailing_zeros();
        aa &= aa - 1;
        let mut bb = b;
        while bb != 0 {
            let j = bb.trailing_zeros();
            bb &= bb - 1;
            acc ^= nim_mul_pow2(i, j);
        }
    }
    acc
}

/// Nim-exponentiation by an ordinary integer exponent (square-and-multiply
/// in the multiplicative group, using nim-multiplication).
pub fn nim_pow(mut base: u64, mut exp: u64) -> u64 {
    let mut acc = 1u64; // nim multiplicative identity
    while exp > 0 {
        if exp & 1 == 1 {
            acc = nim_mul(acc, base);
        }
        base = nim_mul(base, base);
        exp >>= 1;
    }
    acc
}

/// Nim-square (the Frobenius endomorphism x ↦ x² of On₂). F₂-linear, and a
/// *bijection* on every finite nim-field F_{2^m} — char-2 squaring has no kernel.
#[inline]
pub fn nim_square(x: u64) -> u64 {
    nim_mul(x, x)
}

/// Nim-square-root: the inverse Frobenius. In F_{2^64} every element is a unique
/// square, and `√x = x^{2^{63}}` because `x^{2^{64}} = x` there, so
/// `(x^{2^{63}})² = x`. The root lands in whatever subfield `x` lives in (the
/// global Frobenius restricts to each subfield), so this is also the square root
/// in any F_{2^{2^k}} ⊆ F_{2^64}. Always defined — no `Option`.
pub fn nim_sqrt(x: u64) -> u64 {
    nim_pow(x, 1u64 << 63)
}

/// Field trace F_{2^m} → F₂:  `Tr(x) = x + x² + x⁴ + … + x^{2^{m-1}} ∈ {0,1}`.
/// This is the canonical map realising k/℘(k) ≅ F₂ that the Arf invariant is read
/// through (see `arf.rs`); it is *also* the obstruction to solving the
/// Artin–Schreier equation `y² + y = c` (solvable iff `Tr(c) = 0`). One trace,
/// both roles — that is the unification. `m` must be the degree of a nim-subfield
/// (a power of two: 1, 2, 4, …, 64).
pub fn nim_trace(x: u64, m: u32) -> u64 {
    let mut acc = x;
    let mut t = x;
    for _ in 1..m {
        t = nim_square(t);
        acc ^= t;
    }
    acc
}

/// Insert `val` (with its associated y-combination `yc`) into an XOR pivot table
/// keyed by highest set bit. Used by the Artin–Schreier solver.
fn xor_basis_insert(table: &mut [Option<(u64, u64)>; 64], mut val: u64, mut yc: u64) {
    while val != 0 {
        let h = (63 - val.leading_zeros()) as usize;
        match table[h] {
            Some((pv, pc)) => {
                val ^= pv;
                yc ^= pc;
            }
            None => {
                table[h] = Some((val, yc));
                return;
            }
        }
    }
}

/// Solve the Artin–Schreier equation `y² + y = c` in F_{2^m}. The map
/// `L(y) = y² + y` is F₂-linear with kernel {0,1}, and its image is exactly the
/// trace-zero hyperplane — so a solution exists **iff `nim_trace(c, m) = 0`**, and
/// when it does there are exactly two (`y` and `y+1`). Returns one solution, or
/// `None` when `c` is not in the image. Solved by Gaussian elimination over F₂ on
/// the bit-basis of F_{2^m} (exact; no fragile closed-form).
pub fn nim_solve_artin_schreier(c: u64, m: u32) -> Option<u64> {
    let mut table: [Option<(u64, u64)>; 64] = [None; 64];
    for k in 0..m {
        let e = 1u64 << k;
        xor_basis_insert(&mut table, nim_square(e) ^ e, e);
    }
    let (mut val, mut yc) = (c, 0u64);
    while val != 0 {
        let h = (63 - val.leading_zeros()) as usize;
        match table[h] {
            Some((pv, pc)) => {
                val ^= pv;
                yc ^= pc;
            }
            None => return None, // c ∉ image(L)  ⇔  Tr(c) ≠ 0
        }
    }
    Some(yc)
}

/// Whether `y² + y = c` is solvable in F_{2^m} — i.e. `Tr(c) = 0`. The same
/// trace, hence the same answer, as the Arf-reduction path.
pub fn nim_is_artin_schreier_solvable(c: u64, m: u32) -> bool {
    nim_trace(c, m) == 0
}

/// Nim-multiplicative inverse in F_{2^64}. The group F_{2^64}^* is cyclic of
/// order 2^64 − 1, so x^(2^64 − 2) = x^{-1}; and the inverse in the big field
/// agrees with the inverse in whatever subfield x actually lives in.
pub fn nim_inv(x: u64) -> Option<u64> {
    if x == 0 {
        None
    } else {
        Some(nim_pow(x, u64::MAX - 1)) // u64::MAX - 1 = 2^64 - 2
    }
}

/// A nimber, i.e. an element of On_2 truncated to F_{2^64}.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct Nimber(pub u64);

impl std::fmt::Debug for Nimber {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "*{}", self.0)
    }
}

impl Scalar for Nimber {
    fn zero() -> Self {
        Nimber(0)
    }
    fn one() -> Self {
        Nimber(1)
    }
    fn add(&self, rhs: &Self) -> Self {
        Nimber(nim_add(self.0, rhs.0))
    }
    fn neg(&self) -> Self {
        // characteristic 2: every element is its own additive inverse
        *self
    }
    fn mul(&self, rhs: &Self) -> Self {
        Nimber(nim_mul(self.0, rhs.0))
    }
    fn characteristic() -> u32 {
        2
    }
    fn inv(&self) -> Option<Self> {
        nim_inv(self.0).map(Nimber)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn add_is_xor_and_self_inverse() {
        for a in 0u64..64 {
            for b in 0u64..64 {
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
        for a in 0u64..16 {
            for b in 0u64..16 {
                // commutativity
                assert_eq!(nim_mul(a, b), nim_mul(b, a));
                // closure within F_16
                assert!(nim_mul(a, b) < 16, "{a} (x) {b} left F_16");
                for c in 0u64..16 {
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
        for a in 1u64..16 {
            let inv = (1u64..16).find(|&x| nim_mul(a, x) == 1);
            assert!(inv.is_some(), "no inverse for {a} in F_16");
        }
    }

    #[test]
    fn inverse_round_trips() {
        // x ⊗ x^{-1} = 1 for a spread of nonzero nimbers across several fields.
        for x in [
            1u64,
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
            u64::MAX,
        ] {
            let inv = nim_inv(x).unwrap();
            assert_eq!(nim_mul(x, inv), 1, "inverse of {x}");
        }
        assert_eq!(nim_inv(0), None);
        // matches the brute-forced inverses inside F_16
        for x in 1u64..16 {
            let brute = (1u64..16).find(|&y| nim_mul(x, y) == 1).unwrap();
            assert_eq!(nim_inv(x).unwrap(), brute, "F_16 inverse of {x}");
        }
    }

    #[test]
    fn sqrt_is_inverse_frobenius() {
        // √ is the unique inverse of squaring in char 2: (√x)² = x and √(x²) = x.
        for x in [
            0u64,
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
            u64::MAX,
        ] {
            assert_eq!(nim_square(nim_sqrt(x)), x, "(√{x})² ≠ {x}");
            assert_eq!(nim_sqrt(nim_square(x)), x, "√({x}²) ≠ {x}");
        }
        // a square root stays inside the subfield its argument lives in (F_16).
        for x in 0u64..16 {
            assert!(nim_sqrt(x) < 16, "√{x} left F_16");
        }
    }

    #[test]
    fn trace_is_in_f2_and_is_additive() {
        // Tr lands in {0,1} and is F₂-linear (additive) over F_16.
        for x in 0u64..16 {
            assert!(nim_trace(x, 4) <= 1);
            for y in 0u64..16 {
                assert_eq!(nim_trace(x ^ y, 4), nim_trace(x, 4) ^ nim_trace(y, 4));
            }
        }
    }

    #[test]
    fn artin_schreier_solvable_iff_trace_zero() {
        // The unification: y²+y=c is solvable exactly when Tr(c)=0, and the solver
        // returns a genuine root when it is. Checked exhaustively on F_16.
        let m = 4;
        for c in 0u64..16 {
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
        let solvable_count = (0u64..16).filter(|&c| nim_trace(c, m) == 0).count();
        assert_eq!(solvable_count, 8);
    }

    #[test]
    fn artin_schreier_over_f256() {
        // larger field: solver agrees with the trace obstruction on a sample.
        let m = 8;
        for c in (0u64..256).step_by(7) {
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
        for &(a, b, c) in &[(255u64, 256, 257), (1000, 999, 7), (65535, 65536, 3)] {
            assert_eq!(nim_mul(nim_mul(a, b), c), nim_mul(a, nim_mul(b, c)));
        }
    }
}
