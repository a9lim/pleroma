//! Finite-field analysis for the finite nimber backend.
//!
//! The shared Galois engine is the [`FiniteField`] trait
//! (`scalar::finite_field`) — `impl FiniteField for Nimber` supplies the
//! Frobenius and field shape, and the free `nim_*` functions here are the
//! `u128`-keyed surface the Python layer binds. Everything lives in F_{2^128}
//! (= `Nimber(u128)`), whose subfields are exactly F_{2^d} for
//! d | 128 = {1,2,4,…,128}.

use super::{nim_inv, nim_mul, nim_pow, nim_square, Nimber};
use crate::scalar::finite_field::FiniteField;
use std::collections::HashMap;

/// The degree of `x` over F₂: the least `d ∈ {1,2,4,…,128}` (a divisor of 128)
/// with `x^{2^d} = x`, i.e. the dimension of the smallest nim-subfield F_{2^d}
/// containing `x`. `nim_degree(0) = nim_degree(1) = 1`.
pub fn nim_degree(x: u128) -> u128 {
    Nimber(x).degree() as u128
}

/// The distinct Galois conjugates of `x` over F₂: `x, x², x⁴, …, x^{2^{d-1}}`
/// where `d = nim_degree(x)`. These are exactly the roots of `x`'s minimal
/// polynomial, each appearing once.
pub fn nim_conjugates(x: u128) -> Vec<u128> {
    Nimber(x).conjugates().into_iter().map(|n| n.0).collect()
}

/// The minimal polynomial of `x` over F₂, as coefficients in `{0,1}` from the
/// constant term up. Monic of degree [`nim_degree`], so the last entry is the
/// leading `1`. The shared `∏(X + xᵢ)` construction is
/// [`FiniteField::min_poly_monic`]; this projects each coefficient (Galois
/// closure guarantees it lands in F₂) to its bit value.
pub fn nim_min_poly(x: u128) -> Vec<u128> {
    Nimber(x)
        .min_poly_monic()
        .into_iter()
        .map(|c| {
            debug_assert!(c.0 <= 1, "minimal-polynomial coefficient left F₂");
            c.0
        })
        .collect()
}

/// Relative trace `Tr_{F_{2^m}/F_{2^e}}(x) = Σ_{i=0}^{m/e−1} x^{2^{ei}}` — the
/// F_{2^e}-linear surjection onto the subfield. [`super::nim_trace`] is the
/// `e = 1` case (target F₂). Requires `e | m`.
pub fn nim_relative_trace(x: u128, m: u128, e: u128) -> u128 {
    Nimber(x).relative_trace_over(m as usize, e as usize).0
}

/// Relative norm `N_{F_{2^m}/F_{2^e}}(x) = ∏_{i=0}^{m/e−1} x^{2^{ei}}
/// = x^{(2^m−1)/(2^e−1)}` — the multiplicative companion of the relative trace.
/// The norm to the *prime* field (`e = 1`) is always `1` for nonzero `x` (F₂* is
/// trivial), so the relative norm to a larger subfield is the informative one.
/// Requires `e | m`.
pub fn nim_relative_norm(x: u128, m: u128, e: u128) -> u128 {
    Nimber(x).relative_norm_over(m as usize, e as usize).0
}

/// The distinct prime factors of `2^128 − 1` (which is squarefree):
/// `3 · 5 · 17 · 257 · 641 · 65537 · 274177 · 6700417 · 67280421310721`.
/// (The Fermat-number factorizations: 2^32+1 = 641·6700417, 2^64+1 =
/// 274177·67280421310721.) Every multiplicative order in F_{2^128}* divides this.
pub(super) const ORDER_FACTORS: [u128; 9] =
    [3, 5, 17, 257, 641, 65537, 274177, 6700417, 67280421310721];

/// The multiplicative order of `x` in F_{2^128}* — the least `k > 0` with the
/// `k`-fold nim-power `x^{⊗k} = 1`. `None` for `x = 0`. Always divides `2^128−1`.
pub fn nim_order(x: u128) -> Option<u128> {
    Nimber(x).multiplicative_order()
}

/// Whether `x` generates the *full* group F_{2^128}* (order `2^128 − 1`). An
/// element can generate a proper subfield's group without being primitive here;
/// a primitive element necessarily lies outside every proper subfield (so
/// `x ≥ 2^64`).
pub fn nim_is_primitive(x: u128) -> bool {
    Nimber(x).is_primitive()
}

/// A primitive element of F_{2^128}* (a generator of the whole multiplicative
/// group). Searches upward from `2^64` — the floor below which everything sits
/// in the proper subfield F_{2^64}. Deterministic; primitive elements have
/// density `φ(2^128−1)/(2^128−1) ≈ 0.50`, so this returns quickly.
pub fn nim_primitive_element() -> u128 {
    let mut x = 1u128 << 64;
    loop {
        if nim_is_primitive(x) {
            return x;
        }
        x += 1;
    }
}

/// `(a·b) mod m` for `a, b < m ≤ 2^128−1`'s largest prime factor (`≈ 6.7e13`),
/// so `a·b < 4.5e27 < u128::MAX` and the direct product is exact.
#[inline]
fn mulmod(a: u128, b: u128, m: u128) -> u128 {
    (a * b) % m
}

/// Modular inverse `a⁻¹ mod m` by extended Euclid (`a, m` fit comfortably in
/// `i128`). `None` iff `gcd(a,m) ≠ 1`.
fn mod_inv(a: u128, m: u128) -> Option<u128> {
    let (mut old_r, mut r) = (a as i128, m as i128);
    let (mut old_s, mut s) = (1i128, 0i128);
    while r != 0 {
        let q = old_r / r;
        old_r -= q * r;
        std::mem::swap(&mut old_r, &mut r);
        old_s -= q * s;
        std::mem::swap(&mut old_s, &mut s);
    }
    if old_r != 1 {
        return None;
    }
    Some(old_s.rem_euclid(m as i128) as u128)
}

/// Garner CRT: the unique `e ∈ [0, ∏mᵢ)` with `e ≡ rᵢ (mod mᵢ)` for pairwise
/// coprime `mᵢ`. The partial products stay below `2^128−1`, so plain `u128`
/// arithmetic suffices.
fn crt(residues: &[u128], moduli: &[u128]) -> Option<u128> {
    let mut e: u128 = 0;
    let mut radix: u128 = 1;
    for (&r, &m) in residues.iter().zip(moduli) {
        let diff = (r % m + m - e % m) % m;
        let inv = mod_inv(radix % m, m)?;
        let coeff = mulmod(diff, inv, m);
        e += coeff * radix; // < ∏ so far ≤ 2^128−1
        radix = radix.checked_mul(m)?;
    }
    Some(e)
}

/// Baby-step/giant-step in a cyclic subgroup of *prime* order `p`: the `k ∈
/// [0,p)` with the `k`-fold nim-power `g^{⊗k} = h`, or `None` if `h ∉ ⟨g⟩`. `g`
/// must have order exactly `p`. Cost ≈ `√p` time and memory.
fn bsgs_prime_order(g: u128, h: u128, p: u128) -> Option<u128> {
    if h == 1 {
        return Some(0);
    }
    let mut m = (p as f64).sqrt() as u128;
    while m * m < p {
        m += 1;
    }
    m = m.max(1);
    let mut table: HashMap<u128, u128> = HashMap::with_capacity(m as usize);
    let mut cur = 1u128;
    for j in 0..m {
        table.entry(cur).or_insert(j);
        cur = nim_mul(cur, g);
    }
    let factor = nim_inv(nim_pow(g, m))?; // g^{−m}
    let mut gamma = h;
    for i in 0..m {
        if let Some(&j) = table.get(&gamma) {
            return Some(i * m + j);
        }
        gamma = nim_mul(gamma, factor);
    }
    None
}

/// Discrete logarithm in F_{2^128}*: the least `e ≥ 0` with the `e`-fold
/// nim-power `base^{⊗e} = x`, or `None` if `x ∉ ⟨base⟩` (or `base = 0`). Solved
/// by Pohlig–Hellman over the (squarefree) factorization of `ord(base)` with a
/// baby-step/giant-step per prime, recombined by CRT.
///
/// Cost is dominated by the largest prime `p | ord(base)`, at `≈ √p`. For a
/// primitive `base` that prime is `67280421310721`, so the table is `≈ 8.2·10⁶`
/// entries — feasible but heavy; logs inside a proper subfield (small `ord`) are
/// effectively instant.
pub fn nim_discrete_log(base: u128, x: u128) -> Option<u128> {
    if base == 0 {
        return None;
    }
    if x == 1 {
        return Some(0);
    }
    if x == base {
        return Some(1); // cheap shortcut: avoids full Pohlig–Hellman for the trivial log
    }
    let n = nim_order(base)?;
    let mut moduli = Vec::new();
    let mut residues = Vec::new();
    for &p in &ORDER_FACTORS {
        if n % p != 0 {
            continue;
        }
        let g_p = nim_pow(base, n / p);
        let h_p = nim_pow(x, n / p);
        residues.push(bsgs_prime_order(g_p, h_p, p)?); // None ⟹ x ∉ ⟨base⟩
        moduli.push(p);
    }
    let e = crt(&residues, &moduli)?;
    (nim_pow(base, e) == x).then_some(e)
}

/// `Nimber` plugs into the shared [`FiniteField`] engine: the Frobenius is
/// nim-squaring, the field is `F_{2^128}` (`ext_degree = 128`), and the
/// multiplicative group has order `2^128 − 1` with the known squarefree
/// factorization `ORDER_FACTORS`. Two methods are overridden with the sharper
/// char-2 / large-field algorithms: [`is_primitive`](FiniteField::is_primitive)
/// (a direct subgroup check, avoiding a full order computation) and
/// [`discrete_log`](FiniteField::discrete_log) (Pohlig–Hellman, vs the trait's
/// brute force — essential for the `≈ 6.7·10¹³` largest prime factor).
impl FiniteField for Nimber {
    fn frobenius(&self) -> Self {
        Nimber(nim_square(self.0))
    }

    fn pow(&self, e: u128) -> Self {
        Nimber(nim_pow(self.0, e))
    }

    fn ext_degree() -> usize {
        128
    }

    fn group_order() -> u128 {
        u128::MAX // 2^128 − 1
    }

    fn group_order_factors() -> Vec<u128> {
        ORDER_FACTORS.to_vec()
    }

    fn is_primitive(&self) -> bool {
        self.0 != 0
            && ORDER_FACTORS
                .iter()
                .all(|&p| nim_pow(self.0, u128::MAX / p) != 1)
    }

    fn discrete_log(&self, x: Nimber) -> Option<u128> {
        nim_discrete_log(self.0, x.0)
    }
}
