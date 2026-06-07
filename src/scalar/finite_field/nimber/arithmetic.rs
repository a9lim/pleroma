//! Core nim-field arithmetic.
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

use std::cell::RefCell;
use std::collections::HashMap;

#[inline]
pub fn nim_add(a: u128, b: u128) -> u128 {
    a ^ b
}

thread_local! {
    // memo for 2^i (x) 2^j, keyed (min(i,j), max(i,j)); bounded 128x128.
    static POW2_MEMO: RefCell<HashMap<(usize, usize), u128>> = RefCell::new(HashMap::new());
}

/// 2^i (x) 2^j.
fn nim_mul_pow2(i: usize, j: usize) -> u128 {
    let key = if i <= j { (i, j) } else { (j, i) };
    if let Some(v) = POW2_MEMO.with(|m| m.borrow().get(&key).copied()) {
        return v;
    }

    // Fermat indices that appear once (clean product) vs twice (must square).
    let single = i ^ j;
    let common = i & j;

    // distinct Fermat powers nim-multiply to their ordinary product = 2^single
    let clean: u128 = 1u128 << single;

    let result = if common == 0 {
        clean
    } else {
        // fold the squared Fermat powers together, then multiply by `clean`
        let mut squared: u128 = 1; // nim multiplicative identity
        let mut c = common;
        while c != 0 {
            let n = c.trailing_zeros() as usize;
            c &= c - 1;
            let f = 1u128 << (1u128 << n); // F_n = 2^(2^n)
            let factor = f ^ (f >> 1); // F_n (x) F_n = (3/2) F_n
            squared = nim_mul(squared, factor);
        }
        nim_mul(clean, squared)
    };

    POW2_MEMO.with(|m| m.borrow_mut().insert(key, result));
    result
}

/// nim-multiplication, by distributing over the bits of both arguments.
pub fn nim_mul(a: u128, b: u128) -> u128 {
    let mut acc = 0u128;
    let mut aa = a;
    while aa != 0 {
        let i = aa.trailing_zeros() as usize;
        aa &= aa - 1;
        let mut bb = b;
        while bb != 0 {
            let j = bb.trailing_zeros() as usize;
            bb &= bb - 1;
            acc ^= nim_mul_pow2(i, j);
        }
    }
    acc
}

/// Nim-exponentiation by an ordinary integer exponent (square-and-multiply
/// in the multiplicative group, using nim-multiplication).
pub fn nim_pow(mut base: u128, mut exp: u128) -> u128 {
    let mut acc = 1u128; // nim multiplicative identity
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
pub fn nim_square(x: u128) -> u128 {
    nim_mul(x, x)
}

/// Nim-square-root: the inverse Frobenius. In F_{2^128} every element is a unique
/// square, and `√x = x^{2^127}` because `x^{2^128} = x` there, so
/// `(x^{2^127})² = x`. The root lands in whatever subfield `x` lives in (the
/// global Frobenius restricts to each subfield), so this is also the square root
/// in any F_{2^{2^k}} ⊆ F_{2^128}. Always defined — no `Option`.
pub fn nim_sqrt(x: u128) -> u128 {
    nim_pow(x, 1u128 << 127)
}

/// Nim-multiplicative inverse in F_{2^128}. The group F_{2^128}^* is cyclic of
/// order 2^128 − 1, so x^(2^128 − 2) = x^{-1}; and the inverse in the big field
/// agrees with the inverse in whatever subfield x actually lives in.
pub fn nim_inv(x: u128) -> Option<u128> {
    if x == 0 {
        None
    } else {
        Some(nim_pow(x, u128::MAX - 1)) // u128::MAX - 1 = 2^128 - 2
    }
}
