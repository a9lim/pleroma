//! The Hilbert symbol over `Q_p` and the Hasse–Minkowski local–global principle
//! over `Q` — where the Hasse invariant finally does classifying work.
//!
//! `oddchar.rs`'s Hilbert symbol is identically `+1` (finite fields have trivial
//! Brauer group, so no nontrivial quaternion algebras). Over `Q_p` the Hilbert
//! symbol `(a, b)_p` is genuinely nontrivial — it detects the quaternion algebra
//! `(a, b)` — and the **Hasse invariant** `∏_{i<j}(a_i, a_j)_v` it builds becomes a
//! real classifying invariant. The payoff is **Hasse–Minkowski**: a quadratic form
//! over `Q` is isotropic iff it is isotropic over `ℝ` and over every `Q_p` — a
//! theorem this module makes executable ([`is_isotropic_q`]).
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

/// A place of `Q`: the real place `ℝ`, or the `p`-adic place `Q_p`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Place {
    Real,
    Prime(u64),
}

// --- elementary number theory (i128 internals; square-free keeps values tiny) ---

/// The square-free part of `n` (sign preserved): `n` with every squared prime
/// factor removed. The canonical representative of `n`'s class in `Q*/Q*²`.
fn square_free(mut n: i128) -> i128 {
    if n == 0 {
        return 0;
    }
    let sign = n.signum();
    n = n.abs();
    let mut res: i128 = 1;
    let mut d: i128 = 2;
    while d * d <= n {
        if n % d == 0 {
            let mut e = 0;
            while n % d == 0 {
                n /= d;
                e += 1;
            }
            if e % 2 == 1 {
                res *= d;
            }
        }
        d += 1;
    }
    if n > 1 {
        res *= n;
    }
    sign * res
}

/// `p`-adic valuation `v_p(n)` (for `n ≠ 0`).
fn val_p(mut n: i128, p: i128) -> u32 {
    let mut k = 0;
    n = n.abs();
    while n % p == 0 {
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
fn legendre(a: i128, p: i128) -> i8 {
    let a = a.rem_euclid(p);
    if a == 0 {
        return 0;
    }
    // a^{(p-1)/2} mod p
    let mut base = a;
    let mut e = (p - 1) / 2;
    let mut acc: i128 = 1;
    while e > 0 {
        if e & 1 == 1 {
            acc = (acc * base) % p;
        }
        base = (base * base) % p;
        e >>= 1;
    }
    if acc == 1 {
        1
    } else {
        -1
    } // acc is p-1 ≡ -1
}

/// Is the nonzero integer `n` a square in `Q_p`? `v_p(n)` even **and** the unit part
/// is a square unit (`≡ □ mod p` for odd `p`; `≡ 1 mod 8` for `p = 2`).
pub fn is_square_qp(n: i64, p: u64) -> bool {
    let n = n as i128;
    let p = p as i128;
    if n == 0 {
        return false;
    }
    if val_p(n, p) % 2 != 0 {
        return false;
    }
    let u = unit_part(n, p);
    if p == 2 {
        u.rem_euclid(8) == 1
    } else {
        legendre(u, p) == 1
    }
}

// --- the Hilbert symbol ---

/// The Hilbert symbol `(a, b)_∞` over `ℝ`: `−1` iff both `a, b < 0`, else `+1`.
pub fn hilbert_symbol_real(a: i64, b: i64) -> i8 {
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

/// The Hilbert symbol `(a, b)_p` over `Q_p`, for nonzero integers `a, b`. Standard
/// explicit formulas (Serre III.1): for odd `p`, with `a = p^α u`, `b = p^β v`,
/// `(a,b)_p = (−1)^{αβ ε(p)} (u|p)^β (v|p)^α`; for `p = 2`,
/// `(a,b)_2 = (−1)^{ε(u)ε(v) + α ω(v) + β ω(u)}`.
pub fn hilbert_symbol_qp(a: i64, b: i64, p: u64) -> i8 {
    let a = square_free(a as i128);
    let b = square_free(b as i128);
    assert!(a != 0 && b != 0, "Hilbert symbol needs nonzero arguments");
    let pi = p as i128;
    let (al, be) = (val_p(a, pi), val_p(b, pi));
    let (ua, ub) = (unit_part(a, pi), unit_part(b, pi));
    if p == 2 {
        let expo = (eps2(ua) * eps2(ub) + (al as i128) * omega2(ub) + (be as i128) * omega2(ua))
            .rem_euclid(2);
        if expo == 0 {
            1
        } else {
            -1
        }
    } else {
        let eps = ((pi - 1) / 2).rem_euclid(2);
        let mut s: i8 = if (al as i128 * be as i128) % 2 == 1 && eps == 1 {
            -1
        } else {
            1
        };
        if be % 2 == 1 {
            s *= legendre(ua, pi);
        }
        if al % 2 == 1 {
            s *= legendre(ub, pi);
        }
        s
    }
}

/// The Hilbert symbol at an arbitrary place of `Q` (named `_at` to avoid clashing
/// with the finite-field [`oddchar::hilbert_symbol`](crate::forms::hilbert_symbol)).
pub fn hilbert_symbol_at(a: i64, b: i64, place: Place) -> i8 {
    match place {
        Place::Real => hilbert_symbol_real(a, b),
        Place::Prime(p) => hilbert_symbol_qp(a, b, p),
    }
}

// --- Hasse invariant and Hasse–Minkowski ---

/// The Hasse invariant `ε_v(⟨a_1,…,a_n⟩) = ∏_{i<j} (a_i, a_j)_v` at a place `v`
/// (Serre's convention). Entries must be nonzero.
pub fn hasse_at_place(entries: &[i64], place: Place) -> i8 {
    let mut h = 1i8;
    for i in 0..entries.len() {
        for j in (i + 1)..entries.len() {
            h *= hilbert_symbol_at(entries[i], entries[j], place);
        }
    }
    h
}

/// The square class of the discriminant `∏ a_i`, kept square-free / small.
fn disc_class(entries: &[i64]) -> i64 {
    let mut d: i128 = 1;
    for &e in entries {
        d = square_free(d * square_free(e as i128));
    }
    d as i64
}

/// The primes that can carry a nontrivial local condition: `2` together with every
/// prime dividing some entry.
fn relevant_primes(entries: &[i64]) -> BTreeSet<u64> {
    let mut ps = BTreeSet::new();
    ps.insert(2);
    for &e in entries {
        let mut n = (e as i128).abs();
        let mut d: i128 = 2;
        while d * d <= n {
            if n % d == 0 {
                ps.insert(d as u64);
                while n % d == 0 {
                    n /= d;
                }
            }
            d += 1;
        }
        if n > 1 {
            ps.insert(n as u64);
        }
    }
    ps
}

/// Is a perfect square (over `ℤ`, hence over `Q`)?
fn is_perfect_square(n: i128) -> bool {
    if n < 0 {
        return false;
    }
    let r = (n as f64).sqrt() as i128;
    (r - 1..=r + 1).any(|k| k >= 0 && k * k == n)
}

/// Local isotropy of a nondegenerate integer diagonal form over `Q_p`, by rank
/// (Serre IV.2.2): n=1 never; n=2 iff `−d` is a square; n=3 iff `(−1,−d)_p = ε_p`;
/// n=4 iff `d` is a nonsquare or `ε_p = (−1,−1)_p`; n≥5 always.
fn is_isotropic_at_p(entries: &[i64], p: u64) -> bool {
    let n = entries.len();
    let d = disc_class(entries);
    match n {
        0 | 1 => false,
        2 => is_square_qp((-(entries[0] as i128 * entries[1] as i128)) as i64, p),
        3 => hilbert_symbol_qp(-1, -d, p) == hasse_at_place(entries, Place::Prime(p)),
        4 => {
            !is_square_qp(d, p)
                || hasse_at_place(entries, Place::Prime(p)) == hilbert_symbol_qp(-1, -1, p)
        }
        _ => true,
    }
}

/// Whether a diagonal form `⟨a_1,…,a_n⟩` over `Q` is **isotropic** (represents 0
/// nontrivially), by the **Hasse–Minkowski** principle: isotropic over `Q` iff
/// isotropic over `ℝ` and over every `Q_p`. A zero entry is an isotropic direction;
/// otherwise rank 1 is anisotropic, rank 2 needs `−a_1 a_2` a global square, and
/// rank ≥ 3 needs `ℝ`-indefiniteness plus the local condition at each prime dividing
/// `2·∏a_i` (all other primes are automatically isotropic for rank ≥ 3).
pub fn is_isotropic_q(entries: &[i64]) -> bool {
    if entries.iter().any(|&e| e == 0) {
        return true; // a null coordinate direction
    }
    let n = entries.len();
    if n <= 1 {
        return false;
    }
    if n == 2 {
        // ⟨a,b⟩ isotropic over Q iff −ab is a (global) square.
        return is_perfect_square(-(entries[0] as i128 * entries[1] as i128));
    }
    // rank ≥ 3: real place must be indefinite …
    let has_pos = entries.iter().any(|&e| e > 0);
    let has_neg = entries.iter().any(|&e| e < 0);
    if !(has_pos && has_neg) {
        return false; // definite over ℝ ⇒ anisotropic at ∞
    }
    // … and isotropic at every relevant prime.
    relevant_primes(entries)
        .into_iter()
        .all(|p| is_isotropic_at_p(entries, p))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hilbert_symbol_is_symmetric_and_bimultiplicative_seed() {
        // symmetry
        for &p in &[2u64, 3, 5, 7] {
            for a in [-3i64, -1, 1, 2, 3, 5, 6] {
                for b in [-3i64, -1, 1, 2, 3, 5, 6] {
                    assert_eq!(
                        hilbert_symbol_qp(a, b, p),
                        hilbert_symbol_qp(b, a, p),
                        "(a,b)_{p} symmetry"
                    );
                }
            }
        }
        // (a, -a)_v = 1 and (a, 1-a) = 1 are the defining Steinberg relations; check
        // (a,-a): z² = a x² − a y² has (x,y,z)=(1,1,0).
        for &p in &[2u64, 3, 5] {
            for a in [-3i64, -1, 1, 2, 3, 5] {
                assert_eq!(hilbert_symbol_qp(a, -a, p), 1, "(a,−a)_{p} = 1");
            }
        }
    }

    /// The Hilbert reciprocity oracle: ∏ over all places = +1.
    fn reciprocity_holds(a: i64, b: i64) -> bool {
        let mut prod = hilbert_symbol_real(a, b);
        // nontrivial only at primes dividing 2ab
        let mut primes = relevant_primes(&[a, b]);
        primes.insert(2);
        for p in primes {
            prod *= hilbert_symbol_qp(a, b, p);
        }
        prod == 1
    }

    #[test]
    fn hilbert_reciprocity() {
        // THE GOLD ORACLE: ∏_v (a,b)_v = +1 for all a,b — Hilbert's reciprocity law.
        for a in -12i64..=12 {
            for b in -12i64..=12 {
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
        assert_eq!(hilbert_symbol_qp(-1, -1, 2), -1);
        assert_eq!(hilbert_symbol_real(-1, -1), -1);
        // … and (−1,−1) is trivial at every odd prime.
        for &p in &[3u64, 5, 7, 11] {
            assert_eq!(hilbert_symbol_qp(-1, -1, p), 1);
        }
        // (2,3)_? : ∏ must still be +1 (reciprocity), with some nontrivial local one.
        assert!(reciprocity_holds(2, 3));
    }

    #[test]
    fn is_square_qp_basics() {
        // 2 is a square in Q_7 iff 2 is a QR mod 7: 3²=2, yes.
        assert!(is_square_qp(2, 7));
        // 3 is a nonsquare mod 7.
        assert!(!is_square_qp(3, 7));
        // p has odd valuation ⇒ never a square.
        assert!(!is_square_qp(7, 7));
        assert!(!is_square_qp(5, 7)); // 5 is a nonresidue mod 7
                                      // mod 8 rule at p = 2: units ≡ 1 mod 8 are squares.
        assert!(is_square_qp(17, 2)); // 17 ≡ 1 mod 8
        assert!(!is_square_qp(3, 2)); // 3 ≢ 1 mod 8
        assert!(is_square_qp(4, 2)); // 4 = 2², even valuation, unit 1
    }

    #[test]
    fn three_squares_and_sums_of_squares() {
        // ⟨1,1,1⟩: x²+y²+z²=0 has no nontrivial rational solution ⇒ anisotropic.
        assert!(!is_isotropic_q(&[1, 1, 1]));
        // ⟨1,1,-1⟩: (1,0,1) ⇒ isotropic.
        assert!(is_isotropic_q(&[1, 1, -1]));
        // ⟨1,1,1,1⟩: positive definite ⇒ anisotropic over ℝ ⇒ over Q.
        assert!(!is_isotropic_q(&[1, 1, 1, 1]));
        // ⟨1,1,1,-1⟩: indefinite, rank 4, and isotropic over Q (e.g. 1+0+0-1).
        assert!(is_isotropic_q(&[1, 1, 1, -1]));
        // any rank-5 indefinite form is isotropic (u-invariant of Q at every p ≤ 4).
        assert!(is_isotropic_q(&[1, 1, 1, 1, -1]));
        // but rank-5 DEFINITE is not (real place).
        assert!(!is_isotropic_q(&[1, 1, 1, 1, 1]));
    }

    #[test]
    fn classic_anisotropic_ternaries() {
        // x² + y² = 3 z²  ⇔  ⟨1,1,-3⟩ isotropic. No rational solution (3 ≡ 3 mod 4
        // is not a sum of two rational squares) ⇒ anisotropic.
        assert!(!is_isotropic_q(&[1, 1, -3]));
        // x² + y² = 2 z²  ⇔  ⟨1,1,-2⟩: (1,1,1) works ⇒ isotropic.
        assert!(is_isotropic_q(&[1, 1, -2]));
        // x² + y² = 5 z²  ⇔  ⟨1,1,-5⟩: (1,2,1) works (1+4=5) ⇒ isotropic.
        assert!(is_isotropic_q(&[1, 1, -5]));
        // ⟨1,-2,-5⟩ vs ⟨1,-2,-7⟩ etc. — spot-check via reciprocity-backed locals.
        assert!(is_isotropic_q(&[1, 1, -25])); // −25·... actually -1·1 ·... 5² ⇒ iso
    }

    #[test]
    fn rank_two_is_global_square_condition() {
        // ⟨a,b⟩ isotropic iff −ab is a perfect square.
        assert!(is_isotropic_q(&[1, -1])); // −(−1)=1 = 1²
        assert!(is_isotropic_q(&[2, -8])); // −(2·−8)=16 = 4²
        assert!(!is_isotropic_q(&[1, -2])); // −(−2)=2, not a square
        assert!(!is_isotropic_q(&[1, 1])); // −1 not a square
    }
}
