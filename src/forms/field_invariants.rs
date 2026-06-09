//! Numeric **field invariants** the Witt machinery implies but did not yet
//! expose: the **level** (Stufe), the **Pythagoras number**, and the
//! **u-invariant**. Computed honestly over the finite prime fields `F_p`
//! (the textbook constants for the infinite fields are noted in the docs).
//!
//! These round out the Witt-ring picture (`witt_ring.rs` has `Iⁿ` and the `eₙ`
//! staircase): where those are the *cohomological* invariants, these are the
//! classical *numeric* ones, all reducible to sum-of-squares / isotropy
//! questions, which over a finite field are finite computations.
//!
//! Reference values to check against: every finite field has level `1` or `2`
//! and Pythagoras number `≤ 2`; `u(F_q) = 2` (odd `q`). For comparison,
//! formally-real ℝ has level `∞` (no finite `n`), `u(ℝ) = ∞`, Pythagoras number
//! `1`; and `u(Q_p) = 4`.

use crate::scalar::Fp;
use std::collections::BTreeSet;

/// The squares of `F_P` (as residues in `[0, P)`).
fn squares_mod<const P: u128>() -> Vec<u128> {
    (0..P).map(|x| (x * x) % P).collect()
}

/// The set of elements of `F_P` that are sums of exactly `n` squares
/// (`n = 0` is `{0}`).
fn sums_of_n_squares<const P: u128>(n: usize) -> BTreeSet<u128> {
    if n == 0 {
        return BTreeSet::from([0]);
    }
    let squares = squares_mod::<P>();
    let mut cur: BTreeSet<u128> = squares.iter().copied().collect();
    for _ in 1..n {
        let mut next = BTreeSet::new();
        for &a in &cur {
            for &s in &squares {
                next.insert((a + s) % P);
            }
        }
        cur = next;
    }
    cur
}

/// Is `x` a sum of exactly `n` squares in `F_P`?
pub fn is_sum_of_n_squares<const P: u128>(x: Fp<P>, n: usize) -> bool {
    sums_of_n_squares::<P>(n).contains(&(x.value() % P))
}

/// The **level (Stufe)** `s(F_P)`: the least `n` with `−1` a sum of `n` squares,
/// or `None` if `P` is not prime. A finite field has level `1` (iff `−1` is a
/// square: char 2, or `p ≡ 1 mod 4`) or `2`. (ℝ has level `∞` — no finite `n`.)
pub fn level<const P: u128>() -> Option<usize> {
    if !Fp::<P>::modulus_is_prime() {
        return None;
    }
    // −1 in F_P. Level ≤ 2 for any finite field; the search to 4 is a margin.
    let minus_one = (P - 1) % P;
    (1..=4).find(|&n| sums_of_n_squares::<P>(n).contains(&minus_one))
}

/// The **Pythagoras number** `p(F_P)`: the least `n` such that every sum of
/// squares is already a sum of `n` squares (the sum-of-squares set stabilizes).
/// `None` if `P` is not prime. `≤ 2` for finite fields.
pub fn pythagoras_number<const P: u128>() -> Option<usize> {
    if !Fp::<P>::modulus_is_prime() {
        return None;
    }
    let mut prev = sums_of_n_squares::<P>(1);
    for n in 1..=(P as usize + 1) {
        let next = sums_of_n_squares::<P>(n + 1);
        if next == prev {
            return Some(n);
        }
        prev = next;
    }
    Some(P as usize)
}

/// Whether some `code`-indexed nonzero vector isotropes the diagonal form `qs`
/// over `F_P` (brute force over `F_P^dim`).
fn is_anisotropic<const P: u128>(qs: &[u128]) -> bool {
    let dim = qs.len();
    let mut total = 1u128;
    for _ in 0..dim {
        total *= P;
    }
    for code in 1..total {
        // skip the all-zero vector (code 0)
        let mut c = code;
        let mut s = 0u128;
        for &q in qs {
            let xi = c % P;
            c /= P;
            s = (s + q * ((xi * xi) % P)) % P;
        }
        if s == 0 {
            return false; // a nontrivial zero ⇒ isotropic
        }
    }
    true
}

/// Does some diagonal form of dimension `dim` with entries in `F_P*` stay
/// anisotropic?
fn exists_anisotropic_form<const P: u128>(dim: usize) -> bool {
    let mut total = 1u128;
    for _ in 0..dim {
        total *= P - 1;
    }
    for code in 0..total {
        let mut c = code;
        let mut qs = Vec::with_capacity(dim);
        for _ in 0..dim {
            qs.push(1 + c % (P - 1)); // an entry in [1, P-1] = F_P*
            c /= P - 1;
        }
        if is_anisotropic::<P>(&qs) {
            return true;
        }
    }
    false
}

/// The **u-invariant** `u(F_P)`: the largest dimension of an anisotropic form
/// (computed by exhausting diagonal forms up to dim 4). `Some(2)` for every odd
/// prime field (`None` for `P = 2` — char-2 forms are not diagonal — or non-prime
/// `P`). For comparison `u(ℝ) = ∞`, `u(Q_p) = 4`.
pub fn u_invariant<const P: u128>() -> Option<usize> {
    if P == 2 || !Fp::<P>::modulus_is_prime() {
        return None;
    }
    let mut u = 0;
    for dim in 1..=4 {
        if exists_anisotropic_form::<P>(dim) {
            u = dim;
        } else {
            break; // finite-field anisotropic dimensions are an initial segment
        }
    }
    Some(u)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn level_of_finite_fields() {
        // level 1 ⇔ −1 is a square (char 2, or p ≡ 1 mod 4); else level 2.
        assert_eq!(level::<2>(), Some(1)); // char 2: −1 = 1 is a square
        assert_eq!(level::<3>(), Some(2)); // 3 ≡ 3 mod 4
        assert_eq!(level::<5>(), Some(1)); // 5 ≡ 1 mod 4, −1 = 4 = 2²
        assert_eq!(level::<7>(), Some(2)); // 7 ≡ 3 mod 4
        assert_eq!(level::<13>(), Some(1)); // 13 ≡ 1 mod 4
        assert_eq!(level::<9>(), None); // not prime
    }

    #[test]
    fn pythagoras_number_of_finite_fields() {
        assert_eq!(pythagoras_number::<2>(), Some(1)); // every element a square
        assert_eq!(pythagoras_number::<3>(), Some(2));
        assert_eq!(pythagoras_number::<5>(), Some(2));
        assert_eq!(pythagoras_number::<7>(), Some(2));
    }

    #[test]
    fn u_invariant_of_finite_fields_is_two() {
        assert_eq!(u_invariant::<3>(), Some(2));
        assert_eq!(u_invariant::<5>(), Some(2));
        assert_eq!(u_invariant::<7>(), Some(2));
        assert_eq!(u_invariant::<2>(), None); // char 2 out of scope
    }

    #[test]
    fn sum_of_squares_spot_checks() {
        assert!(is_sum_of_n_squares::<3>(Fp::<3>::from_u128(2), 2)); // 2 = 1 + 1
        assert!(!is_sum_of_n_squares::<3>(Fp::<3>::from_u128(2), 1)); // 2 is not a square mod 3
        assert!(is_sum_of_n_squares::<5>(Fp::<5>::from_u128(4), 1)); // 4 = 2²
    }
}
