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
//! the tower carries (high→low):
//!
//! - `χ_{u^{k+1}}^u = χ_{u^k}` (`k ≥ 1`): a digit `≥ u` at level `k` carries one down
//!   to level `k-1` — exact, keeps a single monomial;
//! - `χ_u^u = α_u` (`k = 0`, the **Kummer** relation): a digit `≥ u` at level 0
//!   consumes a factor of the *excess* `α_u`.
//!
//! ## The branching reduction (non-scalar `α_u`)
//!
//! The excess `α_u` is a finite nimber for some primes (`α_3=2`, `α_5=4`, `α_17=16`)
//! and a genuine transfinite ordinal for others (`α_7=ω+1`, `α_11=ω^ω+1`, `α_13=ω+4`,
//! …). A scalar `α_u` keeps a level-0 carry inside the coefficient — the product stays
//! one monomial. A non-scalar `α_u` is a *sum*, so the carry **branches** the monomial:
//! `χ_u^u = α_u` replaces a generator power by `α_u`, and the (reduced) monomial must be
//! nim-multiplied by that sum, mixing across exponent places.
//!
//! This recursion **descends by place**: every `α_{p(m)}` is built from generators at
//! places strictly `< m` (`α_7 = ω+1` uses `ω = χ_3`, place 0 < 2; `α_11 = ω^ω+1` uses
//! `χ_5`, place 1 < 3; verified from DiMuro Table 1). So `base ⊗ excess` can never
//! re-trigger a carry at the place that produced it, and the recursion bottoms out at
//! `α_3 = 2` in the finite field — the crate's "recurse only on strictly-simpler
//! exponents" discipline. Termination depth is bounded by the largest place index.
//!
//! ## Staging (honest scope)
//!
//! We carry the DiMuro Table 1 excesses through `α_43` plus the locally verified
//! `α_47 = ω^(ω^7)+1` from `experiments/ordinal_excess_probe.py` (see `OPEN.md`).
//! The product of any two ordinals `< ω^(ω^ω)` is therefore exact whenever every
//! Kummer carry it triggers is at a prime `≤ 47`; a carry needing `α_53` or beyond
//! returns `None` — the honest operational boundary, moved up from the earlier
//! "any non-scalar `α_u`" cut. (Anything `≥ ω^(ω^ω)`, an infinite exponent place,
//! is out of range outright.)

use super::Ordinal;
use crate::scalar::{is_prime_u128, nim_mul};
use std::collections::{BTreeMap, BTreeSet};

/// A monomial's exponent `E < ω^ω`, decomposed per `ω`-place and per prime:
/// `place m ↦ base-p(m) digit vector`, where digit `k` is the power of the generator
/// `χ_{p(m)^{k+1}}` (`0 ≤ digit < p(m)`). An absent place / empty vector is all-zero;
/// the empty map is the exponent `0` (the monomial `1`).
type GenKey = BTreeMap<u128, Vec<u128>>;

/// The prime governing exponent-place `ω^m`: `p(m)` = the `(m+2)`-th prime
/// (`p(0)=3`, `p(1)=5`, `p(2)=7`, …). Prime 2 is excluded — the prime-2 (Fermat)
/// tower is the finite nimber field, handled by [`crate::scalar::nim_mul`].
fn place_prime(m: u128) -> u128 {
    let mut count = 0u128;
    let mut n = 2u128; // skip the prime 2
    loop {
        n += 1;
        if is_prime_u128(n) {
            count += 1;
            if count == m + 1 {
                return n;
            }
        }
    }
}

/// The excess `α_u` (`χ_u^u = α_u`, the Kummer relation) as an ordinal, or `None` for
/// primes beyond the verified table (`u > 47` — the staged boundary). Every
/// `α_u` is built from generators at strictly-lower places than `χ_u`'s own, which is
/// what makes the branching reduction descend and terminate. Values through `43`:
/// DiMuro Table 1; value `47`: local fixed-base finite-field oracle (see `OPEN.md`);
/// square brackets in the source table are ordinary ordinal exponentiation, already
/// resolved (`[2^ω]=ω`, `[2^{ω²}]=ω^ω`, …).
fn alpha_ordinal(u: u128) -> Option<Ordinal> {
    let fin = Ordinal::from_u128;
    let wpow = Ordinal::omega_pow;
    let w = Ordinal::omega;
    let val = match u {
        3 => fin(2),
        5 => fin(4),
        7 => w().nim_add(&fin(1)),        // ω + 1
        11 => wpow(w()).nim_add(&fin(1)), // ω^ω + 1
        13 => w().nim_add(&fin(4)),       // ω + 4
        17 => fin(16),
        19 => wpow(fin(3)).nim_add(&fin(4)),       // ω³ + 4
        23 => wpow(wpow(fin(3))).nim_add(&fin(1)), // ω^(ω³) + 1
        29 => wpow(wpow(fin(2))).nim_add(&fin(4)), // ω^(ω²) + 4
        31 => wpow(w()).nim_add(&fin(1)),          // ω^ω + 1
        37 => wpow(fin(3)).nim_add(&fin(4)),       // ω³ + 4
        41 => wpow(w()).nim_add(&fin(1)),          // ω^ω + 1
        43 => wpow(wpow(fin(2))).nim_add(&fin(1)), // ω^(ω²) + 1
        47 => wpow(wpow(fin(7))).nim_add(&fin(1)), // ω^(ω⁷) + 1
        _ => return None,
    };
    Some(val)
}

/// Base-`base` digit vector of `v` (least-significant first, no trailing zeros).
fn base_digits(mut v: u128, base: u128) -> Vec<u128> {
    let mut d = Vec::new();
    while v > 0 {
        d.push(v % base);
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
fn recompose_exp(key: &GenKey) -> Option<Ordinal> {
    key.iter().try_fold(Ordinal::zero(), |acc, (&m, digits)| {
        let u = place_prime(m);
        let mut v: u128 = 0;
        let mut pw: u128 = 1;
        for &d in digits {
            let term = d.checked_mul(pw)?;
            v = v.checked_add(term)?;
            pw = pw.checked_mul(u)?;
        }
        Some(if v == 0 {
            acc
        } else {
            acc.nim_add(&Ordinal::monomial(Ordinal::from_u128(m), v))
        })
    })
}

/// Reduce one place's raw (post-addition) generator-power digits to canonical digits
/// `< u`, returning the canonical digits and the number of **level-0 (Kummer) carries**
/// `q` (each owes one factor of the excess `α_u`, resolved by the caller). Processes
/// high→low: a digit `≥ u` at level `k ≥ 1` carries one to level `k-1`
/// (`χ_{u^{k+1}}^u = χ_{u^k}`); at level 0 it is removed and counted (`χ_u^u = α_u`).
fn reduce_place(raw: &[u128], u: u128) -> Option<(Vec<u128>, u128)> {
    let mut d = raw.to_vec();
    for k in (0..d.len()).rev() {
        let carry = d[k] / u;
        d[k] %= u;
        if carry == 0 {
            continue;
        }
        if k == 0 {
            while d.last() == Some(&0) {
                d.pop();
            }
            return Some((d, carry));
        }
        d[k - 1] = d[k - 1].checked_add(carry)?;
    }
    let mut digits: Vec<u128> = d;
    while digits.last() == Some(&0) {
        digits.pop();
    }
    Some((digits, 0))
}

/// Add two generator monomials' digit vectors per `(m,k)` and reduce each place,
/// returning the canonical base [`GenKey`] and the per-place count of level-0 Kummer
/// carries (the excess `α_{p(m)}` owed). Pure digit bookkeeping — no `α` resolution.
fn reduce_keys(a: &GenKey, b: &GenKey) -> Option<(GenKey, BTreeMap<u128, u128>)> {
    let mut base = GenKey::new();
    let mut overflow: BTreeMap<u128, u128> = BTreeMap::new();
    let places: BTreeSet<u128> = a.keys().chain(b.keys()).copied().collect();
    for m in places {
        let da = a.get(&m).map(Vec::as_slice).unwrap_or(&[]);
        let db = b.get(&m).map(Vec::as_slice).unwrap_or(&[]);
        let len = da.len().max(db.len());
        let raw: Vec<u128> = (0..len)
            .map(|i| {
                da.get(i)
                    .copied()
                    .unwrap_or(0)
                    .checked_add(db.get(i).copied().unwrap_or(0))
            })
            .collect::<Option<_>>()?;
        let (red, q) = reduce_place(&raw, place_prime(m))?;
        if q > 0 {
            overflow.insert(m, q);
        }
        if !red.is_empty() {
            base.insert(m, red);
        }
    }
    Some((base, overflow))
}

/// The product of two generator monomials `ω^{E_a}·c_a` and `ω^{E_b}·c_b`, as a full
/// ordinal (a *sum*, once a non-scalar Kummer carry branches it). Adds the generator
/// powers, reduces, then nim-multiplies in the excess `α` factors the level-0 carries
/// owe — recursively, since `α_u` is itself a (strictly-lower-place) ordinal. `None` if
/// some owed `α_u` is past the verified table (`u > 47`).
fn mul_mono(ka: &GenKey, ca: u128, kb: &GenKey, cb: u128) -> Option<Ordinal> {
    let (base_key, overflow) = reduce_keys(ka, kb)?;
    let coeff = nim_mul(ca, cb);
    let base = if base_key.is_empty() {
        Ordinal::from_u128(coeff)
    } else {
        Ordinal::monomial(recompose_exp(&base_key)?, coeff)
    };
    if overflow.is_empty() {
        return Some(base);
    }
    // Excess factor `∏_m α_{p(m)}^{⊗ q_m}`. Each `α` lives at places `< m` (DiMuro), so
    // both this fold and `base ⊗ excess` descend in place and terminate.
    let mut excess = Ordinal::from_u128(1);
    for (&m, &q) in &overflow {
        let alpha = alpha_ordinal(place_prime(m))?;
        for _ in 0..q {
            excess = mul(&excess, &alpha)?;
        }
    }
    mul(&base, &excess)
}

/// Nim-multiply two ordinals `< ω^(ω^ω)`, or `None` outside that range / when a Kummer
/// carry needs an excess `α_u` past the verified table (`u > 47`). Distributes over CNF
/// (char-2 field addition = nim-add); each monomial pair is handled by [`mul_mono`].
pub(super) fn mul(a: &Ordinal, b: &Ordinal) -> Option<Ordinal> {
    let mut acc = Ordinal::zero();
    for (ea, ca) in a.terms() {
        let ka = decompose_exp(ea)?;
        for (eb, cb) in b.terms() {
            let kb = decompose_exp(eb)?;
            let term = mul_mono(&ka, *ca, &kb, *cb)?;
            acc = acc.nim_add(&term);
        }
    }
    Some(acc)
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
    fn chi7() -> Ordinal {
        Ordinal::omega_pow(Ordinal::omega_pow(fin(2))) // ω^(ω²) = χ_7
    }
    /// `χ_7^⊗n` by repeated nim-multiplication.
    fn chi7_pow(n: u128) -> Ordinal {
        let mut p = fin(1);
        for _ in 0..n {
            p = mul(&p, &chi7()).unwrap();
        }
        p
    }

    #[test]
    fn place_primes_are_the_odd_primes() {
        for (m, p) in [
            (0, 3),
            (1, 5),
            (2, 7),
            (3, 11),
            (4, 13),
            (5, 17),
            (6, 19),
            (7, 23),
            (8, 29),
            (9, 31),
            (10, 37),
            (11, 41),
            (12, 43),
            (13, 47),
            (14, 53),
        ] {
            assert_eq!(place_prime(m), p);
        }
    }

    #[test]
    fn alpha_excesses_descend_in_place() {
        // The termination invariant: every verified α_u is built from generators at
        // places strictly below χ_u's own place (the (u)-index in the odd primes). This
        // is what makes `base ⊗ excess` descend; a typo that violated it would loop.
        for m in 0..=13u128 {
            let u = place_prime(m);
            let alpha = alpha_ordinal(u).unwrap();
            // the highest place appearing anywhere in α_u must be < m.
            let mut hi: Option<u128> = None;
            for (exp, _) in alpha.terms() {
                if let Some(sub) = decompose_exp(exp) {
                    if let Some(&mx) = sub.keys().last() {
                        hi = Some(hi.map_or(mx, |h| h.max(mx)));
                    }
                }
            }
            if let Some(h) = hi {
                assert!(h < m, "α_{u} reaches place {h} ≥ its own place {m}");
            }
        }
    }

    #[test]
    fn reproduces_cube_tower_below_omega_omega() {
        // ω ⊗ ω = ω², ω⊗³ = 2, ω⊗⁴ = ω·2 — the prime-3, place-0 behavior.
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
        assert_eq!(w5, fin(4)); // (ω^ω)⊗⁵ = α_5 = 4
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
    fn septic_kummer_landmark() {
        // THE Stage-2 headline, from DiMuro Table 1 (NOT from the engine — non-circular):
        // χ_7 = ω^(ω²), and χ_7^⊗7 = α_7 = ω + 1. The 7th power is the first non-scalar
        // Kummer carry; it branches the monomial into the sum ω + 1.
        assert_eq!(chi7_pow(7), w().nim_add(&fin(1))); // ω + 1

        // χ_7^⊗8 = α_7 ⊗ χ_7 = (ω+1)⊗ω^(ω²) = ω^(ω²+1) + ω^(ω²).
        let e_w2_1 = Ordinal::omega_pow(fin(2)).nim_add(&fin(1)); // ω² + 1
        let w2 = Ordinal::omega_pow(fin(2)); // ω²
        assert_eq!(
            chi7_pow(8),
            Ordinal::omega_pow(e_w2_1).nim_add(&Ordinal::omega_pow(w2))
        );

        // χ_7^⊗9 = α_7 ⊗ χ_7^⊗2 = (ω+1)⊗ω^(ω²·2) = ω^(ω²·2+1) + ω^(ω²·2). Hand-verified
        // both ways (= χ_7^⊗8 ⊗ χ_7), so it also pins associativity through the carry.
        let w2_2 = Ordinal::monomial(fin(2), 2); // ω²·2  (exponent)
        let w2_2_1 = w2_2.nim_add(&fin(1)); // ω²·2 + 1
        assert_eq!(
            chi7_pow(9),
            Ordinal::omega_pow(w2_2_1).nim_add(&Ordinal::omega_pow(w2_2))
        );
        assert_eq!(chi7_pow(9), mul(&chi7_pow(8), &chi7()).unwrap());
    }

    #[test]
    fn locally_verified_alpha_47_landmark() {
        // `experiments/ordinal_excess_probe.py` independently verifies Lenstra excess
        // m_47 = 1 by a fixed-base finite-field power test. Since f(47)=23, this gives
        // α_47 = κ_23 + 1 = ω^(ω^7) + 1.
        let chi47 = Ordinal::omega_pow(Ordinal::omega_pow(fin(13)));
        let mut pow = fin(1);
        for _ in 0..47 {
            pow = mul(&pow, &chi47).unwrap();
        }
        assert_eq!(
            pow,
            Ordinal::omega_pow(Ordinal::omega_pow(fin(7))).nim_add(&fin(1))
        );
    }

    #[test]
    fn quintic_stage_field_axioms() {
        // The prime-3/prime-5 (scalar-α) commutative-ring sweep on a sample of ordinals
        // < ω^(ω²) spanning both the place-ω⁰ and place-ω¹ towers, coeffs in F_4. Every
        // product is defined here; associativity is what a digit-carry bug would break.
        let w = Ordinal::omega();
        let wn = |n| Ordinal::monomial(fin(1), n); // ω·n = ω^1·n (finite exponent 1)
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
        check_field_axioms(&elems);
    }

    #[test]
    fn septic_stage_field_axioms() {
        // The decisive Stage-2 check: the commutative-ring axioms on a sample built from
        // χ_7 = ω^(ω²) (the first non-scalar-α generator), its powers 1..6, ω (= χ_3,
        // which α_7 = ω+1 drags in), F_4 scalars, and mixed sums. Every pairwise product
        // stays within primes {3, 7} (well inside the verified range ⇒ all `Some`), and associativity /
        // distributivity through the α_7 branching is exactly what a mis-mixed carry
        // would break. The α_7 = ω+1 *value* is source-pinned in `septic_kummer_landmark`.
        let mut elems: Vec<Ordinal> = vec![fin(1), fin(2), fin(3), w(), w().nim_add(&fin(1))];
        for n in 1..=6u128 {
            elems.push(chi7_pow(n));
        }
        elems.push(chi7().nim_add(&w())); // χ_7 + ω
        elems.push(Ordinal::monomial(Ordinal::omega_pow(fin(2)), 2).nim_add(&fin(1))); // χ_7·2 + 1
        elems.push(chi7_pow(3).nim_add(&w()).nim_add(&fin(1))); // χ_7³ + ω + 1
        check_field_axioms(&elems);
    }

    /// Commutativity, identity, associativity, and distributivity over `⊕`, on every
    /// triple of a sample whose pairwise products are all defined.
    fn check_field_axioms(elems: &[Ordinal]) {
        let one = fin(1);
        for a in elems {
            for b in elems {
                let ab = a.nim_mul(b).expect("sample is closed under ⊗");
                assert_eq!(ab, b.nim_mul(a).unwrap(), "non-commutative");
                assert_eq!(a.nim_mul(&one).unwrap(), *a, "identity");
                for c in elems {
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
    fn boundary_returns_none_past_prime_47() {
        // Everything through prime 47 is defined: e.g. χ_47 = ω^(ω^13), free powers fine.
        let chi47 = Ordinal::omega_pow(Ordinal::omega_pow(fin(13)));
        assert!(mul(&chi47, &chi47).is_some()); // (ω^(ω^13))⊗² — free, no carry

        // But a Kummer carry at place 14 (prime 53) is past the verified table ⇒ None.
        // ω^(ω^14·50) = χ_53^⊗50; squaring drives the place-14 digit to 100 ≥ 53, owing
        // the unverified α_53.
        let big = Ordinal::omega_pow(Ordinal::monomial(fin(14), 50)); // ω^(ω^14·50)
        assert_eq!(mul(&big, &big), None);

        // And anything ≥ ω^(ω^ω) (an infinite exponent place) is out of range outright.
        let w_ww = Ordinal::omega_pow(ww()); // ω^(ω^ω)
        assert_eq!(mul(&w_ww, &w()), None);
    }

    #[test]
    fn identity_preserves_large_valid_prime_digits() {
        // p(53)=257, so digit 256 is legal. This used to truncate through `u128`
        // storage and collapse the monomial to 1 even when multiplying by 1.
        assert_eq!(place_prime(53), 257);
        let exp = Ordinal::monomial(fin(53), 256);
        let x = Ordinal::omega_pow(exp);
        assert_eq!(mul(&x, &fin(1)), Some(x.clone()));
        assert_eq!(mul(&fin(1), &x), Some(x));
    }
}
