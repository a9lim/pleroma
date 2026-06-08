//! The transfinite nim-multiplication tower: ordinals `< ω^(ω^ω)` (Conway's
//! algebraically closed field `\bar{F_2}`) as monomials in the prime-power
//! generators `χ_{u^n}`, multiplied by generator-power-vector addition with the
//! per-prime carry relations.
//!
//! This generalizes the degree-3 cube-root tower — the prime-3, exponent-place-`ω^0`
//! special case — to all primes and places (Conway *ONAG* ch. 6; Lenstra *Nim
//! multiplication* 1978; DiMuro *On Onp* arXiv:1108.0962, Thm 3.1 + Table 1).
//!
//! ## The uniform representation
//!
//! Every ordinal `< ω^(ω^ω)` is `Σ ω^E · c` in CNF, with each exponent `E < ω^ω`
//! and `c` a finite nimber. An `E < ω^ω` in turn has CNF `Σ ω^m · V_m` with **finite**
//! exponents `m` and finite coefficients `V_m`. The generator picture:
//!
//! - the exponent-place `ω^m` of `E` is governed by the prime `p(m)` = the
//!   `(m+2)`-th prime (`p(0)=3`, `p(1)=5`, `p(2)=7`, …);
//! - reading `V_m` in base `p(m)` gives the digit `d_{m,k}` = the power of the
//!   generator `χ_{p(m)^{k+1}}` (each `d_{m,k} < p(m)`);
//! - so a pure monomial `ω^E` is `⊗_{m,k} χ_{p(m)^{k+1}}^{⊗ d_{m,k}}`, and a general
//!   element is a finite nim-sum of such monomials with finite-nimber coefficients.
//!
//! We therefore key a monomial by its exponent's decomposition,
//! [`GenKey`] = `place m ↦ base-p(m) digit vector`.
//!
//! ## Multiplication
//!
//! `ω^{E1} ⊗ ω^{E2}`: the generator powers **add** (`χ^i ⊗ χ^j = χ^{i+j}`), i.e. the
//! two keys' digit vectors add componentwise per `(m,k)`; then each place reduces by
//! the tower carries (high→low, exactly as the cube tower's `reduce_key` did):
//!
//! - `χ_{u^{k+1}}^u = χ_{u^k}` (`k ≥ 1`): a digit `≥ u` at level `k` carries one down
//!   to level `k-1`;
//! - `χ_u^u = α_u` (`k = 0`, the **Kummer** relation): a digit `≥ u` at level 0
//!   consumes a factor of the *excess* `α_u`.
//!
//! ## Staging (honest scope)
//!
//! `α_u` is a finite nimber for some primes (`α_3=2`, `α_5=4`, `α_17=16`, …) and a
//! genuine transfinite ordinal for others (`α_7=ω+1`, `α_11=ω^ω+1`, `α_13=ω+4`, …).
//! **Stage 1 (this module)** handles the *scalar* `α_u`: a level-0 carry there just
//! multiplies the coefficient, so the product stays a single reduced monomial. When a
//! level-0 carry needs a *non-scalar* `α_u` it returns `None` — the self-limiting
//! boundary. This already closes every ordinal `< ω^(ω²)` (only primes 3,5 appear
//! there, both with scalar `α`) plus all higher products that never trigger a
//! non-scalar Kummer carry (e.g. `(ω^(ω²))^{⊗k}` for `k < 7`). The non-scalar
//! branching expansion (`α_7 = ω+1` ⇒ a carry splits a monomial into a sum) is Stage 2.

use super::Ordinal;
use crate::scalar::nim_mul;
use std::collections::{BTreeMap, BTreeSet};

/// A monomial's exponent `E < ω^ω`, decomposed per `ω`-place and per prime:
/// `place m ↦ base-p(m) digit vector`, where digit `k` is the power of the generator
/// `χ_{p(m)^{k+1}}` (`0 ≤ digit < p(m)`). An absent place / empty vector is all-zero;
/// the empty map is the exponent `0` (the monomial `1`).
type GenKey = BTreeMap<u128, Vec<u8>>;

/// Whether `n` is prime (trial division; `n` is a small prime index in practice).
fn is_prime(n: u128) -> bool {
    if n < 2 {
        return false;
    }
    let mut d = 2u128;
    while d * d <= n {
        if n.is_multiple_of(d) {
            return false;
        }
        d += 1;
    }
    true
}

/// The prime governing exponent-place `ω^m`: `p(m)` = the `(m+2)`-th prime
/// (`p(0)=3`, `p(1)=5`, `p(2)=7`, …). Prime 2 is excluded — the prime-2 (Fermat)
/// tower is the finite nimber field, handled by [`crate::scalar::nim_mul`].
fn place_prime(m: u128) -> u128 {
    let mut count = 0u128;
    let mut n = 2u128; // skip the prime 2
    loop {
        n += 1;
        if is_prime(n) {
            count += 1;
            if count == m + 1 {
                return n;
            }
        }
    }
}

/// The excess `α_u` as a finite nimber, or `None` if `α_u` is a genuine transfinite
/// ordinal (the Stage-1 boundary — the non-scalar Kummer reduction needs the branching
/// expansion, not yet implemented). Verified On₂ values (DiMuro Table 1; see `NOTES.md`).
fn alpha_scalar(u: u128) -> Option<u128> {
    match u {
        3 => Some(2),
        5 => Some(4),
        17 => Some(16),
        // u = 7, 11, 13, 19, 23, … have non-scalar α_u (Stage 2).
        _ => None,
    }
}

/// Base-`base` digit vector of `v` (least-significant first, no trailing zeros).
fn base_digits(mut v: u128, base: u128) -> Vec<u8> {
    let mut d = Vec::new();
    while v > 0 {
        d.push((v % base) as u8);
        v /= base;
    }
    d
}

/// Decompose an exponent `E` into its [`GenKey`], or `None` if `E ≥ ω^ω` (some CNF
/// exponent of `E` is itself infinite — the whole ordinal is then `≥ ω^(ω^ω)`, beyond
/// the algebraically-closed segment this tower represents).
fn decompose_exp(e: &Ordinal) -> Option<GenKey> {
    let mut key = GenKey::new();
    for (exp, c) in e.terms() {
        let m = exp.as_finite()?; // infinite place ⇒ E ≥ ω^ω ⇒ out of range
        key.insert(m, base_digits(*c, place_prime(m)));
    }
    Some(key)
}

/// Rebuild the exponent ordinal `E = Σ ω^m · V_m` from a [`GenKey`] (`V_m` read back
/// from its base-`p(m)` digits).
fn recompose_exp(key: &GenKey) -> Ordinal {
    key.iter().fold(Ordinal::zero(), |acc, (&m, digits)| {
        let u = place_prime(m);
        let mut v: u128 = 0;
        let mut pw: u128 = 1;
        for &d in digits {
            v += d as u128 * pw;
            pw *= u;
        }
        if v == 0 {
            acc
        } else {
            acc.nim_add(&Ordinal::monomial(Ordinal::from_u128(m), v))
        }
    })
}

/// Reduce one place's raw (post-addition) generator-power digits to canonical digits
/// `< u`, returning `(canonical digits, accumulated scalar)` or `None` if a level-0
/// (Kummer) carry needs a non-scalar `α_u`. Processes high→low: a digit `≥ u` at level
/// `k ≥ 1` carries one to level `k-1` (`χ_{u^{k+1}}^u = χ_{u^k}`); at level 0 it
/// consumes a factor `α_u` (`χ_u^u = α_u`).
fn reduce_place(raw: &[u32], u: u128) -> Option<(Vec<u8>, u128)> {
    let mut d = raw.to_vec();
    let mut scalar = 1u128;
    let uu = u as u32;
    for k in (0..d.len()).rev() {
        while d[k] >= uu {
            d[k] -= uu;
            if k == 0 {
                scalar = nim_mul(scalar, alpha_scalar(u)?);
            } else {
                d[k - 1] += 1;
            }
        }
    }
    let mut digits: Vec<u8> = d.iter().map(|&x| x as u8).collect();
    while digits.last() == Some(&0) {
        digits.pop();
    }
    Some((digits, scalar))
}

/// Multiply two generator monomials by exponent: add their digit vectors per `(m,k)`
/// and reduce each place. Returns the reduced [`GenKey`] and the scalar factor
/// produced by the Kummer carries, or `None` at the Stage-1 boundary.
fn mul_keys(a: &GenKey, b: &GenKey) -> Option<(GenKey, u128)> {
    let mut out = GenKey::new();
    let mut scalar = 1u128;
    let places: BTreeSet<u128> = a.keys().chain(b.keys()).copied().collect();
    for m in places {
        let da = a.get(&m).map(Vec::as_slice).unwrap_or(&[]);
        let db = b.get(&m).map(Vec::as_slice).unwrap_or(&[]);
        let len = da.len().max(db.len());
        let raw: Vec<u32> = (0..len)
            .map(|i| *da.get(i).unwrap_or(&0) as u32 + *db.get(i).unwrap_or(&0) as u32)
            .collect();
        let (red, s) = reduce_place(&raw, place_prime(m))?;
        scalar = nim_mul(scalar, s);
        if !red.is_empty() {
            out.insert(m, red);
        }
    }
    Some((out, scalar))
}

/// Nim-multiply two ordinals `< ω^(ω^ω)`, or `None` outside that range / at the
/// Stage-1 (non-scalar `α_u`) boundary. Distributes over CNF: each monomial pair's
/// coefficients nim-multiply, the exponents multiply via [`mul_keys`], and like
/// monomials XOR-accumulate (char 2).
pub(super) fn mul(a: &Ordinal, b: &Ordinal) -> Option<Ordinal> {
    let mut acc: BTreeMap<GenKey, u128> = BTreeMap::new();
    for (ea, ca) in a.terms() {
        let ka = decompose_exp(ea)?;
        for (eb, cb) in b.terms() {
            let kb = decompose_exp(eb)?;
            let (rk, s) = mul_keys(&ka, &kb)?;
            let coeff = nim_mul(nim_mul(*ca, *cb), s);
            *acc.entry(rk).or_insert(0) ^= coeff;
        }
    }
    Some(acc.into_iter().fold(Ordinal::zero(), |out, (k, c)| {
        if c == 0 {
            out
        } else {
            out.nim_add(&Ordinal::monomial(recompose_exp(&k), c))
        }
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fin(n: u128) -> Ordinal {
        Ordinal::from_u128(n)
    }
    fn w() -> Ordinal {
        Ordinal::omega()
    }
    fn ww() -> Ordinal {
        Ordinal::omega_pow(Ordinal::omega()) // ω^ω
    }

    #[test]
    fn place_primes_are_the_odd_primes() {
        assert_eq!(place_prime(0), 3);
        assert_eq!(place_prime(1), 5);
        assert_eq!(place_prime(2), 7);
        assert_eq!(place_prime(3), 11);
        assert_eq!(place_prime(4), 13);
        assert_eq!(place_prime(5), 17);
    }

    #[test]
    fn reproduces_cube_tower_below_omega_omega() {
        // ω ⊗ ω = ω², ω⊗³ = 2, ω²⊗ω² = ω⁴ = 2⊗ω — the prime-3, place-0 behavior.
        let wsq = mul(&w(), &w()).unwrap();
        assert_eq!(wsq, Ordinal::omega_pow(fin(2)));
        assert_eq!(mul(&wsq, &w()).unwrap(), fin(2)); // ω³ = 2
        assert_eq!(mul(&wsq, &wsq).unwrap(), Ordinal::monomial(fin(1), 2)); // ω⁴ = ω·2
    }

    #[test]
    fn quintic_landmarks_from_dimuro() {
        // ω^ω = χ_5, the degree-5 generator: free powers are ordinary ordinal powers,
        // then the Kummer reduction (ω^ω)⊗⁵ = α_5 = 4.
        let w2 = mul(&ww(), &ww()).unwrap();
        assert_eq!(w2, Ordinal::omega_pow(Ordinal::monomial(fin(1), 2))); // ω^{ω·2}
        let w3 = mul(&w2, &ww()).unwrap();
        assert_eq!(w3, Ordinal::omega_pow(Ordinal::monomial(fin(1), 3))); // ω^{ω·3}
        let w4 = mul(&w3, &ww()).unwrap();
        assert_eq!(w4, Ordinal::omega_pow(Ordinal::monomial(fin(1), 4))); // ω^{ω·4}
        let w5 = mul(&w4, &ww()).unwrap();
        assert_eq!(w5, fin(4)); // (ω^ω)⊗⁵ = α_5 = 4  ← the headline
    }

    #[test]
    fn cross_place_products() {
        // ω^ω ⊗ ω = ω^{ω+1} (exponents ω and 1 add, no carry).
        assert_eq!(
            mul(&ww(), &w()).unwrap(),
            Ordinal::omega_pow(Ordinal::omega().nim_add(&fin(1)))
        );
        // ω^ω ⊗ 2 = ω^ω·2 (a finite-nimber coefficient).
        assert_eq!(
            mul(&ww(), &fin(2)).unwrap(),
            Ordinal::monomial(Ordinal::omega(), 2)
        );
    }

    #[test]
    fn quintic_stage_field_axioms() {
        // The decisive Stage-1 check: the commutative-ring axioms on a sample of
        // ordinals < ω^(ω²) spanning BOTH the prime-3 (place ω^0) and prime-5 (place
        // ω^1) towers, with coefficients in F_4 — every product is defined here, and
        // associativity is what a digit-carry bug would break. (The F_64/g_0 level
        // stays exhaustively pinned by `nim::tests::f4_adjoin_omega_is_a_field`.)
        let w = Ordinal::omega();
        // ω·n = ω^1·n (exponent the FINITE ordinal 1), distinct from ω^ω·n.
        let wn = |n| Ordinal::monomial(fin(1), n);
        let elems: Vec<Ordinal> = vec![
            fin(1),
            fin(2),
            fin(3),                                                         // F_4 scalars
            w.clone(),                                                      // ω = χ_3
            Ordinal::omega_pow(fin(2)),                                     // ω²
            Ordinal::omega_pow(fin(3)),                                     // ω³ = χ_9 (g_1)
            ww(),                                                           // ω^ω = χ_5
            Ordinal::omega_pow(wn(2)),                                      // ω^(ω·2)
            Ordinal::omega_pow(wn(5)),                                      // ω^(ω·5) = χ_25
            Ordinal::omega_pow(w.nim_add(&fin(1))),                         // ω^(ω+1)
            ww().nim_add(&w).nim_add(&fin(1)),                              // ω^ω + ω + 1
            wn(3).nim_add(&fin(2)),                                         // ω·3 + 2
            Ordinal::omega_pow(wn(2)).nim_add(&Ordinal::omega_pow(fin(3))), // ω^(ω·2)+ω³
        ];
        let one = fin(1);
        for a in &elems {
            for b in &elems {
                let ab = a.nim_mul(b).expect("< ω^(ω²) is closed under ⊗");
                assert_eq!(ab, b.nim_mul(a).unwrap(), "non-commutative");
                assert_eq!(a.nim_mul(&one).unwrap(), *a, "identity");
                for c in &elems {
                    let l = ab.nim_mul(c).unwrap();
                    let r = a.nim_mul(&b.nim_mul(c).unwrap()).unwrap();
                    assert_eq!(l, r, "× not associative");
                    let l = a.nim_mul(&b.nim_add(c)).unwrap();
                    let r = ab.nim_add(&a.nim_mul(c).unwrap());
                    assert_eq!(l, r, "× not distributive over ⊕");
                }
            }
        }
    }

    #[test]
    fn boundary_returns_none_on_nonscalar_alpha() {
        // χ_7 = ω^(ω²): free powers are fine (no carry) …
        let w_w2 = Ordinal::omega_pow(Ordinal::omega_pow(fin(2))); // ω^(ω²)
        assert!(mul(&w_w2, &w_w2).is_some()); // (ω^(ω²))⊗² = ω^(ω²·2), free
                                              // … but the 7th power needs the Kummer carry α_7 = ω+1 (non-scalar) ⇒ Stage 2.
        let mut p = w_w2.clone();
        let mut hit_none = false;
        for _ in 0..6 {
            match mul(&p, &w_w2) {
                Some(q) => p = q,
                None => {
                    hit_none = true;
                    break;
                }
            }
        }
        assert!(hit_none, "(ω^(ω²))⊗⁷ must hit the non-scalar-α boundary");
        // and anything ≥ ω^(ω^ω) (an infinite exponent place) is out of range.
        let w_ww = Ordinal::omega_pow(ww()); // ω^(ω^ω)
        assert_eq!(mul(&w_ww, &w()), None);
    }
}
