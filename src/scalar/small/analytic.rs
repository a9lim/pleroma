//! The **analytic layer** for the non-Archimedean local worlds — the p-adic
//! mirror of [`surreal/analytic.rs`](crate::scalar::Surreal#method.sqrt).
//!
//! Where [`Surreal::sqrt`](crate::scalar::Surreal::sqrt) extracts real roots from
//! a Hahn series and [`Surreal::inv_to_terms`](crate::scalar::Surreal::inv_to_terms)
//! sums a Neumann series, the local fields/rings get the operations that make
//! "these worlds can take roots and name canonical lifts" true on the
//! non-Archimedean side:
//!
//!   * [`is_square`](Zp::is_square) / [`sqrt`](Zp::sqrt) — checked squarehood and
//!     root construction. For odd `p` these are **Hensel-lifted** square roots: a
//!     unit is a square iff its residue is a square in the residue field (the seed),
//!     and Newton's iteration `y ← (y + u/y)/2` doubles the correct precision each
//!     step until a fixed point. For the *fields* a `p^{2k}·u` splits off `p^k`;
//!     an odd valuation is never a square. At `p = 2`, squarehood reports only
//!     facts determined by the retained representation, and the root itself remains
//!     unknown unless it is zero or definitely absent.
//!   * [`teichmuller`](Zp::teichmuller) — the **Teichmüller representative** `τ(a)`,
//!     the unique `(q−1)`-th root of unity lifting a residue `a`, via the power
//!     iteration `t ← t^p`. [`WittVec`](crate::scalar::WittVec::teichmuller)
//!     already carries this (its Witt coordinates need it); this adds the same lift
//!     to `Zp`/`Qp`/`Qq`, closing the asymmetry that `Zp = W_K(F_p)` lacked it.
//!
//! ## Scope (honest boundary)
//!
//! The checked root APIs return an outer `None` when the represented precision is
//! insufficient or the construction is not implemented, and `Some(None)` when the
//! value is definitely not a square. `teichmuller` works for every `p` (no division).
//! The natural next operations — `nth_root` (Hensel for `gcd(k, p) = 1`) and the
//! p-adic `log`/`exp` (convergent on `v ≥ 1` / `1 + p𝒪`) — are deliberately left
//! for a follow-up, the same way the surreal layer grew incrementally.

use crate::scalar::analytic::{fp_is_square, fp_sqrt, fq_sqrt};
use crate::scalar::{Fp, Fpn, Qp, Qq, Scalar, WittVec, Zp};

// The residue-field square roots (`fp_is_square` / `fp_sqrt` / `fq_sqrt`) that
// seed the Hensel lifts here now live at the analytic root in
// [`crate::scalar::analytic`], shared with the `ExactRoots` impls for the finite
// fields. This module keeps the p-adic-specific Newton lift and Teichmüller rep.

// ───────────────────────── generic lift helpers ─────────────────────────

/// `base^e` in any [`Scalar`] ring, by square-and-multiply.
fn spow<R: Scalar>(base: &R, mut e: u128) -> R {
    let mut acc = R::one();
    let mut b = base.clone();
    while e > 0 {
        if e & 1 == 1 {
            acc = acc.mul(&b);
        }
        b = b.mul(&b);
        e >>= 1;
    }
    acc
}

/// Hensel/Newton square-root lift in a local ring or field `R`: from a `seed`
/// congruent to `√u` modulo the maximal ideal, iterate `y ← (y + u·y⁻¹)·two_inv`
/// (which doubles the correct precision each step) to a fixed point. `u` must be a
/// unit and `two_inv = 1/2` (a unit iff the residue characteristic is odd).
fn newton_sqrt<R: Scalar>(u: &R, seed: R, two_inv: &R) -> R {
    let mut y = seed;
    for _ in 0..64 {
        let yi = y.inv().expect("Newton sqrt: the unit seed must invert");
        let next = u.mul(&yi).add(&y).mul(two_inv);
        if next == y {
            return next;
        }
        y = next;
    }
    y
}

/// `p^e` as a plain integer (the moduli here fit in `u128`).
fn ipow(p: u128, e: u128) -> u128 {
    let mut acc = 1u128;
    for _ in 0..e {
        acc = acc.checked_mul(p).expect("ipow exceeds u128");
    }
    acc
}

/// Square predicate for the represented residue class in `Z/2^k`. For
/// `a = 2^v u`, `v < k`, a square has even `v`, and odd squares modulo `2^m` are
/// exactly `1 (mod 2)`, `1 (mod 4)`, or `1 (mod 8)` for `m = 1`, `2`, or `≥ 3`.
fn is_square_mod_two_power(a: u128, k: u128) -> bool {
    if a == 0 {
        return true;
    }
    let mut u = a;
    let mut v = 0u128;
    while u.is_multiple_of(2) {
        u /= 2;
        v += 1;
    }
    if !v.is_multiple_of(2) {
        return false;
    }
    let m = k - v;
    match m {
        0 | 1 => true,
        2 => u % 4 == 1,
        _ => u % 8 == 1,
    }
}

/// Square predicate for a capped-relative `Q_2` unit. Returns `None` when the
/// retained unit precision is too short to decide the positive case.
fn q2_unit_is_square(unit: u128, k: u128) -> Option<bool> {
    match k {
        0 => None,
        1 => None,
        2 => {
            if unit % 4 == 3 {
                Some(false)
            } else {
                None
            }
        }
        _ => Some(unit % 8 == 1),
    }
}

// ───────────────────────── Zp = Z/p^k ─────────────────────────

impl<const P: u128, const K: u128> Zp<P, K> {
    /// `1/2` in `Z/p^k`, or `None` at `p = 2` (where `2` is a non-unit).
    fn two_inv() -> Option<Self> {
        Self::one().add(&Self::one()).inv()
    }

    fn is_square_odd(&self) -> bool {
        debug_assert!(P != 2);
        if self.0 == 0 {
            return true; // 0 = 0²
        }
        let v = self.valuation();
        if !v.is_multiple_of(2) {
            return false;
        }
        let unit = self.0 / ipow(P, v);
        fp_is_square(unit % P, P)
    }

    fn sqrt_odd(&self) -> Option<Self> {
        debug_assert!(P != 2);
        if self.0 == 0 {
            return Some(Zp(0));
        }
        let v = self.valuation();
        if !v.is_multiple_of(2) {
            return None;
        }
        let unit_val = self.0 / ipow(P, v);
        let seed_res = fp_sqrt(unit_val % P, P)?;
        let two_inv = Self::two_inv().expect("odd p ⇒ 2 is a unit");
        let root_unit = newton_sqrt(
            &Zp::new(unit_val as i128),
            Zp::new(seed_res as i128),
            &two_inv,
        );
        // reattach p^{v/2}
        Some(Zp::new(ipow(P, v / 2) as i128).mul(&root_unit))
    }

    /// Checked square predicate. For odd `p`, this is the exact Hensel predicate.
    /// For `p = 2`, it decides squarehood in the represented quotient `Z/2^K`.
    pub fn is_square(&self) -> Option<bool> {
        if P != 2 {
            return Some(self.is_square_odd());
        }
        Self::assert_supported_ring();
        Some(is_square_mod_two_power(self.0, K))
    }

    /// Checked square-root entry point. The outer `None` means that the
    /// root-construction algorithm is not implemented for this represented case;
    /// `Some(None)` means the value is definitely not a square.
    pub fn sqrt(&self) -> Option<Option<Self>> {
        if P != 2 {
            return Some(self.sqrt_odd());
        }
        Self::assert_supported_ring();
        if self.0 == 0 {
            return Some(Some(Zp(0)));
        }
        if !is_square_mod_two_power(self.0, K) {
            return Some(None);
        }
        None
    }

    /// The **Teichmüller representative** `τ(a) ∈ Z/p^k` of `a ∈ F_p`: the unique
    /// `(p−1)`-th root of unity with `τ(a) ≡ a (mod p)`, via `t ← t^p`. (`Z/p^k`
    /// is `W_k(F_p)`, so this is the prime-field instance of
    /// [`WittVec::teichmuller`](crate::scalar::WittVec::teichmuller).)
    pub fn teichmuller(a: Fp<P>) -> Self {
        let mut t = Zp::new(a.value() as i128);
        for _ in 0..K {
            t = spow(&t, P);
        }
        t
    }
}

// ───────────────────────── Qp ─────────────────────────

impl<const P: u128, const K: u128> Qp<P, K> {
    fn is_square_odd(&self) -> bool {
        debug_assert!(P != 2);
        match self.valuation() {
            None => true, // 0
            Some(v) if v % 2 != 0 => false,
            Some(_) => fp_is_square(self.unit() % P, P),
        }
    }

    fn sqrt_odd(&self) -> Option<Self> {
        debug_assert!(P != 2);
        let Some(v) = self.valuation() else {
            return Some(Qp::zero());
        };
        if v % 2 != 0 {
            return None;
        }
        let seed_res = fp_sqrt(self.unit() % P, P)?;
        let two_inv = Qp::from_i128(2).inv().expect("odd p ⇒ 2 invertible");
        let unit = Qp::from_i128(self.unit() as i128); // the val-0 unit part
        let root_unit = newton_sqrt(&unit, Qp::from_i128(seed_res as i128), &two_inv);
        Some(Qp::from_p_power(v / 2).mul(&root_unit))
    }

    /// Checked square predicate. At `p = 2`, a nonzero element is a square
    /// iff its valuation is even and its unit is `1 mod 8`; the capped-relative
    /// model can only decide the positive case when at least three unit digits are
    /// retained.
    pub fn is_square(&self) -> Option<bool> {
        if P != 2 {
            return Some(self.is_square_odd());
        }
        Self::assert_supported_field();
        match self.valuation() {
            None => Some(true),
            Some(v) if v % 2 != 0 => Some(false),
            Some(_) => q2_unit_is_square(self.unit(), K),
        }
    }

    /// Checked square-root entry point. Odd `p` delegates to the Hensel
    /// lift. At `p = 2`, the predicate can still reject nonsquares, but the root
    /// construction itself is deliberately not faked.
    pub fn sqrt(&self) -> Option<Option<Self>> {
        if P != 2 {
            return Some(self.sqrt_odd());
        }
        Self::assert_supported_field();
        if self.is_zero() {
            return Some(Some(Qp::zero()));
        }
        match self.is_square()? {
            false => Some(None),
            true => None,
        }
    }

    /// The **Teichmüller representative** `τ(a) ∈ Q_p` of a residue `a ∈ F_p`
    /// (a unit of valuation 0), via `t ← t^p`.
    pub fn teichmuller(a: Fp<P>) -> Self {
        let mut t = Qp::from_i128(a.value() as i128);
        for _ in 0..K {
            t = spow(&t, P);
        }
        t
    }
}

// ───────────────────────── WittVec (sqrt; teichmuller already exists) ─────

impl<const P: u128, const N: usize, const F: usize> WittVec<P, N, F> {
    /// `1/2` in `W_N(F_q)`, or `None` at `p = 2`.
    fn two_inv() -> Option<Self> {
        Self::one().add(&Self::one()).inv()
    }

    /// The unit part `u` with `self = p^{v}·u` (or `None` if `self ≡ 0` to the
    /// retained precision), alongside the valuation `v`.
    fn split_unit(&self) -> Option<(usize, Self)> {
        let v = self.p_valuation();
        if v >= N {
            return None;
        }
        let mut u = *self;
        for _ in 0..v {
            u = u.try_divide_by_p().expect("valuation says p | self");
        }
        Some((v, u))
    }

    fn is_square_odd(&self) -> bool {
        debug_assert!(P != 2);
        match self.split_unit() {
            None => true, // ≡ 0
            Some((v, _)) if v % 2 != 0 => false,
            Some((_, u)) => u.residue().is_square(),
        }
    }

    fn sqrt_odd(&self) -> Option<Self> {
        debug_assert!(P != 2);
        let Some((v, u)) = self.split_unit() else {
            return Some(WittVec::zero());
        };
        if v % 2 != 0 {
            return None;
        }
        let seed_res = fq_sqrt(u.residue())?;
        let two_inv = Self::two_inv().expect("odd p ⇒ 2 is a unit");
        let root_unit = newton_sqrt(&u, WittVec(seed_res.into_coeffs()), &two_inv);
        let mut acc = root_unit;
        let p = WittVec::<P, N, F>::from_int(P as i128);
        for _ in 0..(v / 2) {
            acc = acc.mul(&p);
        }
        Some(acc)
    }

    /// Checked square predicate. Odd `p` delegates to the Hensel predicate.
    /// At `p = 2`, the unramified dyadic unit criterion is not implemented here;
    /// the method reports only the valuation obstruction and otherwise returns
    /// `None` rather than guessing.
    pub fn is_square(&self) -> Option<bool> {
        if P != 2 {
            return Some(self.is_square_odd());
        }
        let v = self.p_valuation();
        if v >= N {
            return Some(true);
        }
        if !v.is_multiple_of(2) {
            return Some(false);
        }
        None
    }

    /// Checked square-root entry point. Odd `p` delegates to the Hensel
    /// lift; dyadic unramified roots return only the zero root or a definite
    /// nonsquare.
    pub fn sqrt(&self) -> Option<Option<Self>> {
        if P != 2 {
            return Some(self.sqrt_odd());
        }
        if self.is_zero() {
            return Some(Some(WittVec::zero()));
        }
        match self.is_square()? {
            false => Some(None),
            true => None,
        }
    }
}

// ───────────────────────── Qq ─────────────────────────

impl<const P: u128, const N: usize, const F: usize> Qq<P, N, F> {
    fn is_square_odd(&self) -> bool {
        debug_assert!(P != 2);
        match self.valuation() {
            None => true, // 0
            Some(v) if v % 2 != 0 => false,
            Some(_) => self.unit_residue().expect("nonzero ⇒ residue").is_square(),
        }
    }

    fn sqrt_odd(&self) -> Option<Self> {
        debug_assert!(P != 2);
        let Some(v) = self.valuation() else {
            return Some(Qq::zero());
        };
        if v % 2 != 0 {
            return None;
        }
        let seed_res = fq_sqrt(self.unit_residue().expect("nonzero ⇒ residue"))?;
        let two_inv = Qq::from_int(2).inv().expect("odd p ⇒ 2 invertible");
        let unit = Qq::from_witt(self.unit()); // val-0 unit part
        let root_unit = newton_sqrt(
            &unit,
            Qq::from_witt(WittVec(seed_res.into_coeffs())),
            &two_inv,
        );
        Some(Qq::from_p_power(v / 2).mul(&root_unit))
    }

    /// Checked square predicate. Odd `p` delegates to the Hensel predicate.
    /// At `p = 2`, the valuation obstruction is exact; the unramified dyadic unit
    /// criterion is left as unknown.
    pub fn is_square(&self) -> Option<bool> {
        if P != 2 {
            return Some(self.is_square_odd());
        }
        match self.valuation() {
            None => Some(true),
            Some(v) if v % 2 != 0 => Some(false),
            Some(_) => None,
        }
    }

    /// Checked square-root entry point. Odd `p` delegates to the Hensel
    /// lift. Dyadic unramified roots return only the zero root or a definite
    /// nonsquare.
    pub fn sqrt(&self) -> Option<Option<Self>> {
        if P != 2 {
            return Some(self.sqrt_odd());
        }
        if self.is_zero() {
            return Some(Some(Qq::zero()));
        }
        match self.is_square()? {
            false => Some(None),
            true => None,
        }
    }

    /// The **Teichmüller representative** `τ(a) ∈ Q_q` of a residue `a ∈ F_q`,
    /// delegating to [`WittVec::teichmuller`](crate::scalar::WittVec::teichmuller).
    pub fn teichmuller(a: Fpn<P, F>) -> Self {
        Qq::from_witt(WittVec::<P, N, F>::teichmuller(a))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ---------- Zp / Qp over F_p ----------

    #[test]
    fn zp_sqrt_round_trips_and_matches_is_square() {
        // Z/3^4, Z/5^4, Z/7^3: every element's square has a root squaring back to it,
        // and is_square exactly classifies the squares.
        fn check<const P: u128, const K: u128>() {
            let m = Zp::<P, K>::modulus();
            for a in 0..m {
                let x = Zp::<P, K>(a);
                let is_sq = x.is_square().expect("odd p squarehood is decidable");
                match x.sqrt().expect("odd p root construction is implemented") {
                    Some(r) => {
                        assert!(is_sq, "sqrt(Some) but !is_square for {x:?}");
                        assert_eq!(r.mul(&r), x, "({r:?})² ≠ {x:?}");
                    }
                    None => assert!(!is_sq, "sqrt(None) but is_square for {x:?}"),
                }
            }
        }
        check::<3, 4>();
        check::<5, 4>(); // p ≡ 1 (mod 4): exercises full Tonelli–Shanks
        check::<7, 3>(); // p ≡ 3 (mod 4): the fast branch
    }

    #[test]
    fn qp_sqrt_handles_valuations() {
        type Q = Qp<5, 5>;
        // a unit square
        let four = Q::from_i128(4);
        let r = four
            .sqrt()
            .expect("odd p root construction is implemented")
            .expect("4 is a square in Q_5");
        assert_eq!(r.mul(&r), four);
        // p^2 · square ⇒ square, root has half the valuation
        let x = Q::from_p_power(2).mul(&four);
        let rx = x
            .sqrt()
            .expect("odd p root construction is implemented")
            .expect("5²·4 is a square");
        assert_eq!(rx.mul(&rx), x);
        assert_eq!(rx.valuation(), Some(1));
        // odd valuation ⇒ never a square
        assert_eq!(Q::from_i128(5).is_square(), Some(false));
        assert_eq!(Q::from_i128(5).sqrt(), Some(None));
        // 2 is a non-residue mod 5 ⇒ not a square in Q_5
        assert_eq!(Q::from_i128(2).is_square(), Some(false));
        assert_eq!(Q::zero().sqrt(), Some(Some(Q::zero())));
    }

    #[test]
    fn teichmuller_is_a_root_of_unity_lifting_the_residue() {
        // τ(a) ≡ a (mod p), and τ(a)^{p-1} = 1 (it is a (p−1)-th root of unity),
        // equivalently τ(a) is Frobenius-fixed: τ(a)^p = τ(a).
        type Z = Zp<7, 4>;
        for a in 1..7u128 {
            let t = Z::teichmuller(Fp::<7>::from_u128(a));
            assert_eq!(t.0 % 7, a, "τ lifts the residue");
            assert_eq!(spow(&t, 7), t, "τ is Frobenius-fixed (τ^p = τ)");
            // a (p−1)-th root of unity
            assert_eq!(spow(&t, 6), Z::one(), "τ^{{p-1}} = 1");
        }
        // Qp agrees with Zp on the lift.
        for a in 1..7u128 {
            let tq = Qp::<7, 4>::teichmuller(Fp::<7>::from_u128(a));
            assert_eq!(tq.unit(), Zp::<7, 4>::teichmuller(Fp::<7>::from_u128(a)).0);
        }
    }

    // ---------- WittVec / Qq over F_q ----------

    #[test]
    fn wittvec_sqrt_over_f9() {
        // W_3(F_9): residue F_9 (odd), exercises the field Tonelli seed + Hensel lift.
        type W = WittVec<3, 3, 2>;
        let q = W::residue_order(); // 9
        let modu = W::modulus();
        // a spread of units (residue ≠ 0) and their squares
        for c0 in 0..modu {
            for c1 in 0..q {
                let x = WittVec::<3, 3, 2>([c0, c1]);
                if x.residue().is_zero() {
                    continue; // focus on units here
                }
                let sq = x.mul(&x);
                assert_eq!(sq.is_square(), Some(true), "a square must read as a square");
                let r = sq
                    .sqrt()
                    .expect("odd p root construction is implemented")
                    .expect("a square has a root");
                assert_eq!(r.mul(&r), sq, "({r:?})² ≠ {sq:?}");
            }
        }
    }

    #[test]
    fn qq_sqrt_and_teichmuller_over_f9() {
        type Q = Qq<3, 3, 2>;
        // unit square in the genuine extension
        let g = WittVec::<3, 3, 2>([0, 1]); // residue = F_9 generator
        let x = Q::from_witt(g);
        let sq = x.mul(&x);
        let r = sq
            .sqrt()
            .expect("odd p root construction is implemented")
            .expect("a square inverts");
        assert_eq!(r.mul(&r), sq);
        // p² · square ⇒ square at half valuation
        let y = Q::from_p_power(2).mul(&sq);
        let ry = y
            .sqrt()
            .expect("odd p root construction is implemented")
            .expect("3²·square is a square");
        assert_eq!(ry.mul(&ry), y);
        assert_eq!(ry.valuation(), Some(1));
        // odd valuation ⇒ not a square
        assert_eq!(Q::from_p_power(1).mul(&sq).is_square(), Some(false));
        // Teichmüller: τ(a) lifts a and is a (q−1)-th root of unity.
        let a = g.residue();
        let t = Q::teichmuller(a);
        assert_eq!(t.unit_residue(), Some(a));
        assert_eq!(t.valuation(), Some(0));
    }

    #[test]
    fn qq_f1_matches_qp() {
        // Q_q with F = 1 IS Q_p: sqrt must agree on the unit residue.
        let x = Qq::<5, 5, 1>::from_int(4);
        let r = x.sqrt().unwrap().unwrap();
        assert_eq!(r.mul(&r), x);
        assert_eq!(Qq::<5, 5, 1>::from_int(2).sqrt(), Some(None)); // 2 a non-residue mod 5
    }

    #[test]
    fn dyadic_square_apis_are_checked_and_honest() {
        type Z = Zp<2, 5>;
        assert_eq!(Zp::<2, 5>(0).is_square(), Some(true));
        assert_eq!(Zp::<2, 5>(1).is_square(), Some(true));
        assert_eq!(Zp::<2, 5>(4).is_square(), Some(true));
        assert_eq!(Zp::<2, 5>(2).is_square(), Some(false)); // odd 2-adic valuation
        assert_eq!(Zp::<2, 5>(3).is_square(), Some(false)); // odd unit not 1 mod 8
        assert_eq!(Zp::<2, 5>(0).sqrt(), Some(Some(Z::zero())));
        assert_eq!(Zp::<2, 5>(3).sqrt(), Some(None));
        assert_eq!(Zp::<2, 5>(1).sqrt(), None); // square known, dyadic root not constructed

        type Q = Qp<2, 5>;
        assert_eq!(Q::zero().is_square(), Some(true));
        assert_eq!(Q::from_i128(1).is_square(), Some(true));
        assert_eq!(Q::from_i128(2).is_square(), Some(false)); // odd valuation
        assert_eq!(Q::from_i128(3).is_square(), Some(false)); // unit 3 mod 8
        assert_eq!(Qp::<2, 2>::from_i128(1).is_square(), None); // not enough unit digits
        assert_eq!(Q::zero().sqrt(), Some(Some(Q::zero())));
        assert_eq!(Q::from_i128(3).sqrt(), Some(None));
        assert_eq!(Q::from_i128(1).sqrt(), None);

        type W = WittVec<2, 4, 2>;
        assert_eq!(W::zero().is_square(), Some(true));
        assert_eq!(W::from_int(2).is_square(), Some(false)); // valuation obstruction
        assert_eq!(W::one().is_square(), None); // dyadic unramified unit path unknown

        type R = Qq<2, 4, 2>;
        assert_eq!(R::zero().is_square(), Some(true));
        assert_eq!(R::from_p_power(1).is_square(), Some(false));
        assert_eq!(R::one().is_square(), None);
        assert_eq!(R::from_p_power(1).sqrt(), Some(None));
    }
}
