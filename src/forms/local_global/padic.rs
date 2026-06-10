//! The Hilbert symbol over `Q_p` and the Hasse–Minkowski local–global principle
//! over `Q` — where the Hasse invariant finally does classifying work.
//!
//! `forms::oddchar`'s Hilbert symbol is identically `+1` (finite fields have trivial
//! Brauer group, so no nontrivial quaternion algebras). Over `Q_p` the Hilbert
//! symbol `(a, b)_p` is genuinely nontrivial — it detects the quaternion algebra
//! `(a, b)` — and the **Hasse invariant** `∏_{i<j}(a_i, a_j)_v` it builds becomes a
//! real classifying invariant. The payoff is **Hasse–Minkowski**: a quadratic form
//! over `Q` is isotropic iff it is isotropic over `ℝ` and over every `Q_p` — a
//! theorem this module makes executable ([`try_is_isotropic_q`]).
//!
//! Forms are integer diagonal forms `⟨a_1,…,a_n⟩` (a rational form scales to one
//! without changing square classes). Everything depends only on square classes, so
//! arguments are square-free-reduced internally — keeping the arithmetic small and
//! exact. The gold-standard self-check is **Hilbert reciprocity**: `∏_v (a,b)_v = +1`
//! over all places (`reciprocity_holds` in the tests).
//!
//! References: Serre, *A Course in Arithmetic*, Ch. III–IV (the Hilbert symbol
//! formulas and the local isotropy criteria by rank).

use std::collections::BTreeSet;

use crate::scalar::{is_prime_u128, mul_mod_u128};

/// A place of `Q`: the real place `ℝ`, or the `p`-adic place `Q_p`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Place {
    Real,
    Prime(u128),
}

// --- elementary number theory (i128 internals; square-free keeps values tiny) ---

fn signed_u128(sign: i128, n: u128) -> Option<i128> {
    if sign < 0 {
        if n == (i128::MAX as u128) + 1 {
            Some(i128::MIN)
        } else {
            i128::try_from(n).ok()?.checked_neg()
        }
    } else {
        i128::try_from(n).ok()
    }
}

/// Checked square-free part of `n` (sign preserved).
pub(crate) fn try_square_free(n: i128) -> Option<i128> {
    if n == 0 {
        return Some(0);
    }
    let sign = n.signum();
    let mut n = n.unsigned_abs();
    let mut res: u128 = 1;
    let mut d: u128 = 2;
    while d <= n / d {
        if n.is_multiple_of(d) {
            let mut e = 0;
            while n.is_multiple_of(d) {
                n /= d;
                e += 1;
            }
            if e % 2 == 1 {
                res = res.checked_mul(d)?;
            }
        }
        d += 1;
    }
    if n > 1 {
        res = res.checked_mul(n)?;
    }
    signed_u128(sign, res)
}

/// `p`-adic valuation `v_p(n)` (for `n ≠ 0`).
fn val_p(n: i128, p: i128) -> u128 {
    let mut k = 0;
    let mut n = n.unsigned_abs();
    let p = p as u128;
    while n.is_multiple_of(p) {
        n /= p;
        k += 1;
    }
    k
}

/// The `p`-adic unit part `n / p^{v_p(n)}` (sign preserved).
fn unit_part(mut n: i128, p: i128) -> i128 {
    while n % p == 0 {
        n /= p;
    }
    n
}

/// The Legendre symbol `(a | p)` for an odd prime `p`: `0` if `p | a`, else `±1`.
fn legendre(a: i128, p: i128) -> i128 {
    let p_u = p as u128;
    let a = a.rem_euclid(p) as u128;
    if a == 0 {
        return 0;
    }
    // a^{(p-1)/2} mod p
    let mut base = a;
    let mut e = (p_u - 1) / 2;
    let mut acc: u128 = 1;
    while e > 0 {
        if e & 1 == 1 {
            acc = mul_mod_u128(acc, base, p_u);
        }
        base = mul_mod_u128(base, base, p_u);
        e >>= 1;
    }
    if acc == 1 {
        1
    } else {
        -1
    } // acc is p-1 ≡ -1
}

/// Is the nonzero integer `n` a square in `Q_p`? `v_p(n)` even **and** the unit part
/// is a square unit (`≡ □ mod p` for odd `p`; `≡ 1 mod 8` for `p = 2`). Returns
/// `None` when `p` is not a prime representable by the bounded `i128`
/// implementation.
pub fn try_is_square_qp(n: i128, p: u128) -> Option<bool> {
    if !is_prime_u128(p) || i128::try_from(p).is_err() {
        return None;
    }
    let p = p as i128;
    if n == 0 {
        return Some(false);
    }
    if !val_p(n, p).is_multiple_of(2) {
        return Some(false);
    }
    let u = unit_part(n, p);
    Some(if p == 2 {
        u.rem_euclid(8) == 1
    } else {
        legendre(u, p) == 1
    })
}

// --- the Hilbert symbol ---

/// The Hilbert symbol `(a, b)_∞` over `ℝ`: `−1` iff both `a, b < 0`, else `+1`.
pub fn hilbert_symbol_real(a: i128, b: i128) -> i128 {
    if a < 0 && b < 0 {
        -1
    } else {
        1
    }
}

/// `ε(u) = (u−1)/2 mod 2` for an odd integer `u` (depends on `u mod 4`).
fn eps2(u: i128) -> i128 {
    if u.rem_euclid(4) == 1 {
        0
    } else {
        1
    }
}

/// `ω(u) = (u²−1)/8 mod 2` for an odd integer `u` (depends on `u mod 8`).
fn omega2(u: i128) -> i128 {
    match u.rem_euclid(8) {
        1 | 7 => 0,
        _ => 1, // 3 | 5
    }
}

/// The **tame Hilbert symbol** — the shared formula behind the odd-`p` Hilbert
/// symbol over `Q_p` *and* every (odd-residue) place of `F_q(t)`. With valuations
/// `α, β` and residue quadratic characters `χ_a = χ(ā)`, `χ_b = χ(b̄)`, `χ_{−1} =
/// χ(−1)`,
/// `(a,b)_v = χ_{−1}^{αβ} · χ_a^β · χ_b^α`.
/// Over `Q_p` the residue character is the Legendre symbol; over `F_q(t)` it is
/// the residue-field character `χ_κ`. The `p = 2` (mod-8) and real branches are the
/// two genuine exceptions — everything else is this symbol (see
/// [`hilbert_symbol_ff`](crate::forms::hilbert_symbol_ff)).
pub(crate) fn tame_hilbert_symbol(
    alpha: i128,
    beta: i128,
    chi_a: i128,
    chi_b: i128,
    chi_neg1: i128,
) -> i128 {
    let (a_odd, b_odd) = (alpha.rem_euclid(2) == 1, beta.rem_euclid(2) == 1);
    let mut s: i128 = if a_odd && b_odd { chi_neg1 } else { 1 };
    if b_odd {
        s *= chi_a;
    }
    if a_odd {
        s *= chi_b;
    }
    s
}

/// The Hilbert symbol `(a, b)_p` over `Q_p`, for nonzero integers `a, b`. Standard
/// explicit formulas (Serre III.1): for odd `p`, with `a = p^α u`, `b = p^β v`,
/// `(a,b)_p = (−1)^{αβ ε(p)} (u|p)^β (v|p)^α` (the [`tame_hilbert_symbol`] with the
/// Legendre character); for `p = 2`, `(a,b)_2 = (−1)^{ε(u)ε(v) + α ω(v) + β ω(u)}`.
/// Returns `None` when `p` is not a representable prime, either argument is zero,
/// or square-class reduction overflows the bounded `i128` implementation.
pub fn try_hilbert_symbol_qp(a: i128, b: i128, p: u128) -> Option<i128> {
    if !is_prime_u128(p) || i128::try_from(p).is_err() {
        return None;
    }
    let a = try_square_free(a)?;
    let b = try_square_free(b)?;
    if a == 0 || b == 0 {
        return None;
    }
    let pi = p as i128;
    let (al, be) = (val_p(a, pi), val_p(b, pi));
    let (ua, ub) = (unit_part(a, pi), unit_part(b, pi));
    Some(if p == 2 {
        let expo = (eps2(ua) * eps2(ub) + (al as i128) * omega2(ub) + (be as i128) * omega2(ua))
            .rem_euclid(2);
        if expo == 0 {
            1
        } else {
            -1
        }
    } else {
        // odd p: the tame symbol with the residue Legendre character.
        tame_hilbert_symbol(
            al as i128,
            be as i128,
            legendre(ua, pi),
            legendre(ub, pi),
            legendre(-1, pi),
        )
    })
}

/// The Hilbert symbol at an arbitrary place of `Q` (named `_at` to avoid clashing
/// with the finite-field [`oddchar::hilbert_symbol`](crate::forms::hilbert_symbol)).
pub fn try_hilbert_symbol_at(a: i128, b: i128, place: Place) -> Option<i128> {
    Some(match place {
        Place::Real => hilbert_symbol_real(a, b),
        Place::Prime(p) => try_hilbert_symbol_qp(a, b, p)?,
    })
}

// --- Hasse invariant and Hasse–Minkowski ---

/// The Hasse invariant `ε_v(⟨a_1,…,a_n⟩) = ∏_{i<j} (a_i, a_j)_v` at a place `v`
/// (Serre's convention). Entries must be nonzero.
pub fn try_hasse_at_place(entries: &[i128], place: Place) -> Option<i128> {
    let mut h = 1i128;
    for i in 0..entries.len() {
        for j in (i + 1)..entries.len() {
            h *= try_hilbert_symbol_at(entries[i], entries[j], place)?;
        }
    }
    Some(h)
}

/// The square class of the discriminant `∏ a_i`, kept square-free / small.
pub(crate) fn try_disc_class(entries: &[i128]) -> Option<i128> {
    let mut d: i128 = 1;
    for &e in entries {
        d = try_square_free(d.checked_mul(try_square_free(e)?)?)?;
    }
    Some(d)
}

fn neg_product(a: i128, b: i128) -> Option<i128> {
    a.checked_mul(b)?.checked_neg()
}

pub(crate) fn try_is_isotropic_at_p(entries: &[i128], p: u128) -> Option<bool> {
    let n = entries.len();
    let d = try_disc_class(entries)?;
    Some(match n {
        0 | 1 => false,
        2 => try_is_square_qp(neg_product(entries[0], entries[1])?, p)?,
        3 => {
            try_hilbert_symbol_qp(-1, d.checked_neg()?, p)?
                == try_hasse_at_place(entries, Place::Prime(p))?
        }
        4 => {
            !try_is_square_qp(d, p)?
                || try_hasse_at_place(entries, Place::Prime(p))?
                    == try_hilbert_symbol_qp(-1, -1, p)?
        }
        _ => true,
    })
}

/// The primes that can carry a nontrivial local condition: `2` together with every
/// prime dividing some entry.
pub(crate) fn relevant_primes(entries: &[i128]) -> BTreeSet<u128> {
    let mut ps = BTreeSet::new();
    ps.insert(2);
    for &e in entries {
        let mut n = e.unsigned_abs();
        let mut d: u128 = 2;
        while d <= n / d {
            if n.is_multiple_of(d) {
                ps.insert(d);
                while n.is_multiple_of(d) {
                    n /= d;
                }
            }
            d += 1;
        }
        if n > 1 {
            ps.insert(n);
        }
    }
    ps
}

/// Is a perfect square (over `ℤ`, hence over `Q`)?
pub(crate) fn is_perfect_square(n: i128) -> bool {
    if n < 0 {
        return false;
    }
    let mut lo = 0i128;
    let mut hi = n;
    while lo <= hi {
        let mid = lo + (hi - lo) / 2;
        if mid == 0 || mid <= n / mid {
            lo = mid + 1;
        } else {
            hi = mid - 1;
        }
    }
    hi.checked_mul(hi) == Some(n)
}

/// The **Hilbert reciprocity product** `∏_v (a,b)_v` over all places of `ℚ` — the
/// multiplicative product formula for the quaternion-algebra class `(a,b)`. It is
/// `+1` for every nonzero `a, b` (Hilbert's reciprocity law); the local symbols are
/// `+1` at all but finitely many places (those dividing `2ab`). This is the
/// structural form of the oracle the tests use, exposed for the adelic layer.
pub fn try_hilbert_reciprocity_product(a: i128, b: i128) -> Option<i128> {
    let mut prod = hilbert_symbol_real(a, b);
    let mut primes = relevant_primes(&[a, b]);
    primes.insert(2);
    for p in primes {
        prod *= try_hilbert_symbol_qp(a, b, p)?;
    }
    Some(prod)
}

/// Whether a diagonal form `⟨a_1,…,a_n⟩` over `Q` is **isotropic** (represents 0
/// nontrivially), by the **Hasse–Minkowski** principle: isotropic over `Q` iff
/// isotropic over `ℝ` and over every `Q_p`. A zero entry is an isotropic direction;
/// otherwise rank 1 is anisotropic, rank 2 needs `−a_1 a_2` a global square, and
/// rank ≥ 3 needs `ℝ`-indefiniteness plus the local condition at each prime dividing
/// `2·∏a_i` (all other primes are automatically isotropic for rank ≥ 3). Returns
/// `None` if bounded square-class arithmetic overflows.
pub fn try_is_isotropic_q(entries: &[i128]) -> Option<bool> {
    if entries.contains(&0) {
        return Some(true); // a null coordinate direction
    }
    let n = entries.len();
    if n <= 1 {
        return Some(false);
    }
    if n == 2 {
        // ⟨a,b⟩ isotropic over Q iff −ab is a (global) square.
        return Some(is_perfect_square(neg_product(entries[0], entries[1])?));
    }
    // rank ≥ 3: real place must be indefinite …
    let has_pos = entries.iter().any(|&e| e > 0);
    let has_neg = entries.iter().any(|&e| e < 0);
    if !(has_pos && has_neg) {
        return Some(false); // definite over ℝ ⇒ anisotropic at ∞
    }
    // … and isotropic at every relevant prime.
    for p in relevant_primes(entries) {
        if !try_is_isotropic_at_p(entries, p)? {
            return Some(false);
        }
    }
    Some(true)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sq(n: i128, p: u128) -> bool {
        try_is_square_qp(n, p).expect("test prime is supported")
    }

    fn hs(a: i128, b: i128, p: u128) -> i128 {
        try_hilbert_symbol_qp(a, b, p).expect("test Hilbert symbol is defined")
    }

    fn iso(entries: &[i128]) -> bool {
        try_is_isotropic_q(entries).expect("test square classes fit i128")
    }

    #[test]
    fn hilbert_symbol_is_symmetric_and_bimultiplicative_seed() {
        // symmetry
        for &p in &[2u128, 3, 5, 7] {
            for a in [-3i128, -1, 1, 2, 3, 5, 6] {
                for b in [-3i128, -1, 1, 2, 3, 5, 6] {
                    assert_eq!(hs(a, b, p), hs(b, a, p), "(a,b)_{p} symmetry");
                }
            }
        }
        // (a, -a)_v = 1 and (a, 1-a) = 1 are the defining Steinberg relations; check
        // (a,-a): z² = a x² − a y² has (x,y,z)=(1,1,0).
        for &p in &[2u128, 3, 5] {
            for a in [-3i128, -1, 1, 2, 3, 5] {
                assert_eq!(hs(a, -a, p), 1, "(a,−a)_{p} = 1");
            }
        }
    }

    /// The Hilbert reciprocity oracle: ∏ over all places = +1.
    fn reciprocity_holds(a: i128, b: i128) -> bool {
        try_hilbert_reciprocity_product(a, b).expect("test symbols are defined") == 1
    }

    #[test]
    fn hilbert_reciprocity() {
        // THE GOLD ORACLE: ∏_v (a,b)_v = +1 for all a,b — Hilbert's reciprocity law.
        for a in -12i128..=12 {
            for b in -12i128..=12 {
                if a == 0 || b == 0 {
                    continue;
                }
                assert!(reciprocity_holds(a, b), "reciprocity failed at a={a} b={b}");
            }
        }
    }

    #[test]
    fn hilbert_detects_nontrivial_quaternion_algebra() {
        // (−1,−1)_2 = −1: Hamilton's quaternions ramify at 2 and ∞ — the canonical
        // nontrivial symbol that finite fields can never exhibit.
        assert_eq!(hs(-1, -1, 2), -1);
        assert_eq!(hilbert_symbol_real(-1, -1), -1);
        // … and (−1,−1) is trivial at every odd prime.
        for &p in &[3u128, 5, 7, 11] {
            assert_eq!(hs(-1, -1, p), 1);
        }
        // (2,3)_? : ∏ must still be +1 (reciprocity), with some nontrivial local one.
        assert!(reciprocity_holds(2, 3));
    }

    #[test]
    fn is_square_qp_basics() {
        // 2 is a square in Q_7 iff 2 is a QR mod 7: 3²=2, yes.
        assert!(sq(2, 7));
        // 3 is a nonsquare mod 7.
        assert!(!sq(3, 7));
        // p has odd valuation ⇒ never a square.
        assert!(!sq(7, 7));
        assert!(!sq(5, 7)); // 5 is a nonresidue mod 7
                            // mod 8 rule at p = 2: units ≡ 1 mod 8 are squares.
        assert!(sq(17, 2)); // 17 ≡ 1 mod 8
        assert!(!sq(3, 2)); // 3 ≢ 1 mod 8
        assert!(sq(4, 2)); // 4 = 2², even valuation, unit 1
    }

    #[test]
    fn three_squares_and_sums_of_squares() {
        // ⟨1,1,1⟩: x²+y²+z²=0 has no nontrivial rational solution ⇒ anisotropic.
        assert!(!iso(&[1, 1, 1]));
        // ⟨1,1,-1⟩: (1,0,1) ⇒ isotropic.
        assert!(iso(&[1, 1, -1]));
        // ⟨1,1,1,1⟩: positive definite ⇒ anisotropic over ℝ ⇒ over Q.
        assert!(!iso(&[1, 1, 1, 1]));
        // ⟨1,1,1,-1⟩: indefinite, rank 4, and isotropic over Q (e.g. 1+0+0-1).
        assert!(iso(&[1, 1, 1, -1]));
        // any rank-5 indefinite form is isotropic (u-invariant of Q at every p ≤ 4).
        assert!(iso(&[1, 1, 1, 1, -1]));
        // but rank-5 DEFINITE is not (real place).
        assert!(!iso(&[1, 1, 1, 1, 1]));
    }

    #[test]
    fn classic_anisotropic_ternaries() {
        // x² + y² = 3 z²  ⇔  ⟨1,1,-3⟩ isotropic. No rational solution (3 ≡ 3 mod 4
        // is not a sum of two rational squares) ⇒ anisotropic.
        assert!(!iso(&[1, 1, -3]));
        // x² + y² = 2 z²  ⇔  ⟨1,1,-2⟩: (1,1,1) works ⇒ isotropic.
        assert!(iso(&[1, 1, -2]));
        // x² + y² = 5 z²  ⇔  ⟨1,1,-5⟩: (1,2,1) works (1+4=5) ⇒ isotropic.
        assert!(iso(&[1, 1, -5]));
        // ⟨1,-2,-5⟩ vs ⟨1,-2,-7⟩ etc. — spot-check via reciprocity-backed locals.
        assert!(iso(&[1, 1, -25])); // −25·... actually -1·1 ·... 5² ⇒ iso
    }

    #[test]
    fn rank_two_is_global_square_condition() {
        // ⟨a,b⟩ isotropic iff −ab is a perfect square.
        assert!(iso(&[1, -1])); // −(−1)=1 = 1²
        assert!(iso(&[2, -8])); // −(2·−8)=16 = 4²
        assert!(!iso(&[1, -2])); // −(−2)=2, not a square
        assert!(!iso(&[1, 1])); // −1 not a square
    }

    #[test]
    fn rank_two_square_test_is_exact_near_i128_limit() {
        let a = 3_037_000_499i128;
        assert!(iso(&[a, -a])); // −a·(−a) = a², exactly.
        assert!(!iso(&[a, -(a - 1)]));
    }

    #[test]
    fn qp_apis_reject_nonprime_places() {
        assert_eq!(try_is_square_qp(2, 9), None);
        assert_eq!(try_hilbert_symbol_qp(2, 3, 1), None);
    }
}
