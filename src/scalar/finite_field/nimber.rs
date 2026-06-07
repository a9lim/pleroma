//! Nimbers: the ordinals under nim-addition and nim-multiplication, Conway's
//! Field On_2 of characteristic 2. Restricted here to `u128`, i.e. nimbers
//! below 2^128 — which *is* exactly the finite nim-field F_{2^128}, and contains
//! every smaller F_{2^{2^k}} (k <= 7: F_4, F_16, F_256, ... F_{2^128}).
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

use super::FiniteField;
use crate::scalar::Scalar;
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

/// Field trace F_{2^m} → F₂:  `Tr(x) = x + x² + x⁴ + … + x^{2^{m-1}} ∈ {0,1}`.
/// This is the canonical map realising k/℘(k) ≅ F₂ that the Arf invariant is read
/// through (see `arf.rs`); it is *also* the obstruction to solving the
/// Artin–Schreier equation `y² + y = c` (solvable iff `Tr(c) = 0`). One trace,
/// both roles — that is the unification. `m` must be the degree of a nim-subfield
/// (a power of two: 1, 2, 4, …, 128).
pub fn nim_trace(x: u128, m: u128) -> u128 {
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
fn xor_basis_insert(table: &mut [Option<(u128, u128)>; 128], mut val: u128, mut yc: u128) {
    while val != 0 {
        let h = (127 - val.leading_zeros()) as usize;
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
pub fn nim_solve_artin_schreier(c: u128, m: u128) -> Option<u128> {
    let mut table: [Option<(u128, u128)>; 128] = [None; 128];
    for k in 0..m {
        let e = 1u128 << k;
        xor_basis_insert(&mut table, nim_square(e) ^ e, e);
    }
    let (mut val, mut yc) = (c, 0u128);
    while val != 0 {
        let h = (127 - val.leading_zeros()) as usize;
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
pub fn nim_is_artin_schreier_solvable(c: u128, m: u128) -> bool {
    nim_trace(c, m) == 0
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

// ===========================================================================
// Finite-field analysis: conjugates, degree, minimal polynomial, relative
// trace/norm, multiplicative order, primitive elements, discrete log.
//
// The shared Galois engine is the `FiniteField` trait (scalar/finite_field) —
// `impl FiniteField for Nimber` below supplies the Frobenius and the field
// shape, and these free `nim_*` functions (the `u128`-keyed surface the Python
// layer binds) delegate to it. Everything lives in F_{2^128} (= `Nimber(u128)`),
// whose subfields are exactly F_{2^d} for d | 128 = {1,2,4,…,128}.
// ===========================================================================

/// The degree of `x` over F₂: the least `d ∈ {1,2,4,…,128}` (a divisor of 128)
/// with `x^{2^d} = x`, i.e. the dimension of the smallest nim-subfield F_{2^d}
/// containing `x`. `nim_degree(0) = nim_degree(1) = 1`.
pub fn nim_degree(x: u128) -> u32 {
    Nimber(x).degree() as u32
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
pub fn nim_min_poly(x: u128) -> Vec<u8> {
    Nimber(x)
        .min_poly_monic()
        .into_iter()
        .map(|c| {
            debug_assert!(c.0 <= 1, "minimal-polynomial coefficient left F₂");
            c.0 as u8
        })
        .collect()
}

/// Relative trace `Tr_{F_{2^m}/F_{2^e}}(x) = Σ_{i=0}^{m/e−1} x^{2^{ei}}` — the
/// F_{2^e}-linear surjection onto the subfield. [`nim_trace`] is the `e = 1`
/// case (target F₂). Requires `e | m`.
pub fn nim_relative_trace(x: u128, m: u32, e: u32) -> u128 {
    Nimber(x).relative_trace_over(m as usize, e as usize).0
}

/// Relative norm `N_{F_{2^m}/F_{2^e}}(x) = ∏_{i=0}^{m/e−1} x^{2^{ei}}
/// = x^{(2^m−1)/(2^e−1)}` — the multiplicative companion of the relative trace.
/// The norm to the *prime* field (`e = 1`) is always `1` for nonzero `x` (F₂* is
/// trivial), so the relative norm to a larger subfield is the informative one.
/// Requires `e | m`.
pub fn nim_relative_norm(x: u128, m: u32, e: u32) -> u128 {
    Nimber(x).relative_norm_over(m as usize, e as usize).0
}

/// The distinct prime factors of `2^128 − 1` (which is squarefree):
/// `3 · 5 · 17 · 257 · 641 · 65537 · 274177 · 6700417 · 67280421310721`.
/// (The Fermat-number factorizations: 2^32+1 = 641·6700417, 2^64+1 =
/// 274177·67280421310721.) Every multiplicative order in F_{2^128}* divides this.
const ORDER_FACTORS: [u128; 9] = [3, 5, 17, 257, 641, 65537, 274177, 6700417, 67280421310721];

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
/// density `φ(2^128−1)/(2^128−1) ≈ 0.30`, so this returns quickly.
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

/// A nimber, i.e. an element of On_2 truncated to F_{2^128}.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct Nimber(pub u128);

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
    fn characteristic() -> u128 {
        2
    }
    fn inv(&self) -> Option<Self> {
        nim_inv(self.0).map(Nimber)
    }
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
        self.0 != 0 && ORDER_FACTORS.iter().all(|&p| nim_pow(self.0, u128::MAX / p) != 1)
    }

    fn discrete_log(&self, x: Nimber) -> Option<u128> {
        nim_discrete_log(self.0, x.0)
    }
}

#[cfg(test)]
mod tests {
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
}
