//! The **analytic layer** for the non-Archimedean local worlds — the p-adic
//! mirror of [`surreal/analytic.rs`](crate::scalar::Surreal#method.sqrt).
//!
//! Where [`Surreal::sqrt`](crate::scalar::Surreal::sqrt) extracts real roots from
//! a Hahn series and [`Surreal::inv_to_terms`](crate::scalar::Surreal::inv_to_terms)
//! sums a Neumann series, the local fields/rings get the operations that make
//! "these worlds can take roots and name canonical lifts" true on the
//! non-Archimedean side:
//!
//!   * [`is_square`](Zp::is_square) / [`sqrt`](Zp::sqrt) — **Hensel-lifted** square
//!     roots. A unit is a square iff its residue is a square in the residue field
//!     (the seed); Newton's iteration `y ← (y + u/y)/2` then doubles the correct
//!     precision each step until a fixed point. For the *fields* a `p^{2k}·u`
//!     splits off `p^k`; an odd valuation is never a square.
//!   * [`teichmuller`](Zp::teichmuller) — the **Teichmüller representative** `τ(a)`,
//!     the unique `(q−1)`-th root of unity lifting a residue `a`, via the power
//!     iteration `t ← t^p`. [`WittVec`](crate::scalar::WittVec::teichmuller)
//!     already carries this (its Witt coordinates need it); this adds the same lift
//!     to `Zp`/`Qp`/`Qq`, closing the asymmetry that `Zp = W_K(F_p)` lacked it.
//!
//! ## Scope (honest boundary)
//!
//! `is_square`/`sqrt` require **odd residue characteristic** and assert it: the
//! dyadic `p = 2` square root is the mod-8 story that lives in the forms layer, not
//! a Newton lift (`2` is not a unit, so the iteration's `1/2` does not exist).
//! `teichmuller` works for every `p` (no division). The natural next operations —
//! `nth_root` (Hensel for `gcd(k, p) = 1`) and the p-adic `log`/`exp` (convergent
//! on `v ≥ 1` / `1 + p𝒪`) — are deliberately left for a follow-up, the same way
//! the surreal layer grew incrementally.

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

// ───────────────────────── Zp = Z/p^k ─────────────────────────

impl<const P: u128, const K: u128> Zp<P, K> {
    /// `1/2` in `Z/p^k`, or `None` at `p = 2` (where `2` is a non-unit).
    fn two_inv() -> Option<Self> {
        Self::one().add(&Self::one()).inv()
    }

    /// Whether this is a square in `Z_p` (odd `p`): a unit is a square iff its
    /// residue is a square in `F_p`; a `p^{2k}·u` is a square iff `u` is; an odd
    /// valuation is never a square. Asserts odd `p`.
    pub fn is_square(&self) -> bool {
        assert!(
            P != 2,
            "Zp::is_square requires odd p (dyadic squares are the forms mod-8 story)"
        );
        if self.0 == 0 {
            return true; // 0 = 0²
        }
        let v = self.valuation();
        if v % 2 != 0 {
            return false;
        }
        let unit = self.0 / ipow(P, v);
        fp_is_square(unit % P, P)
    }

    /// A square root in `Z_p` (the one congruent to the residue root), or `None` if
    /// this is not a square. Asserts odd `p`.
    pub fn sqrt(&self) -> Option<Self> {
        assert!(
            P != 2,
            "Zp::sqrt requires odd p (dyadic squares are the forms mod-8 story)"
        );
        if self.0 == 0 {
            return Some(Zp(0));
        }
        let v = self.valuation();
        if v % 2 != 0 {
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

    /// The **Teichmüller representative** `τ(a) ∈ Z/p^k` of `a ∈ F_p`: the unique
    /// `(p−1)`-th root of unity with `τ(a) ≡ a (mod p)`, via `t ← t^p`. (`Z/p^k`
    /// is `W_k(F_p)`, so this is the prime-field instance of
    /// [`WittVec::teichmuller`](crate::scalar::WittVec::teichmuller).)
    pub fn teichmuller(a: Fp<P>) -> Self {
        let mut t = Zp::new(a.0 as i128);
        for _ in 0..K {
            t = spow(&t, P);
        }
        t
    }
}

// ───────────────────────── Qp ─────────────────────────

impl<const P: u128, const K: u128> Qp<P, K> {
    /// Whether this is a square in `Q_p` (odd `p`): even valuation and a residue
    /// square. Asserts odd `p`.
    pub fn is_square(&self) -> bool {
        assert!(
            P != 2,
            "Qp::is_square requires odd p (dyadic squares are the forms mod-8 story)"
        );
        match self.valuation() {
            None => true, // 0
            Some(v) if v % 2 != 0 => false,
            Some(_) => fp_is_square(self.unit() % P, P),
        }
    }

    /// A square root in `Q_p`, or `None` if not a square. Asserts odd `p`.
    pub fn sqrt(&self) -> Option<Self> {
        assert!(
            P != 2,
            "Qp::sqrt requires odd p (dyadic squares are the forms mod-8 story)"
        );
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

    /// The **Teichmüller representative** `τ(a) ∈ Q_p` of a residue `a ∈ F_p`
    /// (a unit of valuation 0), via `t ← t^p`.
    pub fn teichmuller(a: Fp<P>) -> Self {
        let mut t = Qp::from_i128(a.0 as i128);
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

    /// Whether this is a square in `W_N(F_q)` (odd `p`): even valuation and a
    /// residue square in `F_q`. Asserts odd `p`.
    pub fn is_square(&self) -> bool {
        assert!(
            P != 2,
            "WittVec::is_square requires odd p (dyadic squares are the forms mod-8 story)"
        );
        match self.split_unit() {
            None => true, // ≡ 0
            Some((v, _)) if v % 2 != 0 => false,
            Some((_, u)) => u.residue().is_square(),
        }
    }

    /// A square root in `W_N(F_q)`, or `None` if not a square. Asserts odd `p`.
    pub fn sqrt(&self) -> Option<Self> {
        assert!(
            P != 2,
            "WittVec::sqrt requires odd p (dyadic squares are the forms mod-8 story)"
        );
        let Some((v, u)) = self.split_unit() else {
            return Some(WittVec::zero());
        };
        if v % 2 != 0 {
            return None;
        }
        let seed_res = fq_sqrt(u.residue())?;
        let two_inv = Self::two_inv().expect("odd p ⇒ 2 is a unit");
        let root_unit = newton_sqrt(&u, WittVec(seed_res.0), &two_inv);
        let mut acc = root_unit;
        let p = WittVec::<P, N, F>::from_int(P as i128);
        for _ in 0..(v / 2) {
            acc = acc.mul(&p);
        }
        Some(acc)
    }
}

// ───────────────────────── Qq ─────────────────────────

impl<const P: u128, const N: usize, const F: usize> Qq<P, N, F> {
    /// Whether this is a square in `Q_q` (odd `p`): even valuation and a residue
    /// square in `F_q`. Asserts odd `p`.
    pub fn is_square(&self) -> bool {
        assert!(
            P != 2,
            "Qq::is_square requires odd p (dyadic squares are the forms mod-8 story)"
        );
        match self.valuation() {
            None => true, // 0
            Some(v) if v % 2 != 0 => false,
            Some(_) => self.unit_residue().expect("nonzero ⇒ residue").is_square(),
        }
    }

    /// A square root in `Q_q`, or `None` if not a square. Asserts odd `p`.
    pub fn sqrt(&self) -> Option<Self> {
        assert!(
            P != 2,
            "Qq::sqrt requires odd p (dyadic squares are the forms mod-8 story)"
        );
        let Some(v) = self.valuation() else {
            return Some(Qq::zero());
        };
        if v % 2 != 0 {
            return None;
        }
        let seed_res = fq_sqrt(self.unit_residue().expect("nonzero ⇒ residue"))?;
        let two_inv = Qq::from_int(2).inv().expect("odd p ⇒ 2 invertible");
        let unit = Qq::from_witt(self.unit()); // val-0 unit part
        let root_unit = newton_sqrt(&unit, Qq::from_witt(WittVec(seed_res.0)), &two_inv);
        Some(Qq::from_p_power(v / 2).mul(&root_unit))
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
                let is_sq = x.is_square();
                match x.sqrt() {
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
        let r = four.sqrt().expect("4 is a square in Q_5");
        assert_eq!(r.mul(&r), four);
        // p^2 · square ⇒ square, root has half the valuation
        let x = Q::from_p_power(2).mul(&four);
        let rx = x.sqrt().expect("5²·4 is a square");
        assert_eq!(rx.mul(&rx), x);
        assert_eq!(rx.valuation(), Some(1));
        // odd valuation ⇒ never a square
        assert!(!Q::from_i128(5).is_square());
        assert_eq!(Q::from_i128(5).sqrt(), None);
        // 2 is a non-residue mod 5 ⇒ not a square in Q_5
        assert!(!Q::from_i128(2).is_square());
        assert_eq!(Q::zero().sqrt(), Some(Q::zero()));
    }

    #[test]
    fn teichmuller_is_a_root_of_unity_lifting_the_residue() {
        // τ(a) ≡ a (mod p), and τ(a)^{p-1} = 1 (it is a (p−1)-th root of unity),
        // equivalently τ(a) is Frobenius-fixed: τ(a)^p = τ(a).
        type Z = Zp<7, 4>;
        for a in 1..7u128 {
            let t = Z::teichmuller(Fp::<7>(a));
            assert_eq!(t.0 % 7, a, "τ lifts the residue");
            assert_eq!(spow(&t, 7), t, "τ is Frobenius-fixed (τ^p = τ)");
            // a (p−1)-th root of unity
            assert_eq!(spow(&t, 6), Z::one(), "τ^{{p-1}} = 1");
        }
        // Qp agrees with Zp on the lift.
        for a in 1..7u128 {
            let tq = Qp::<7, 4>::teichmuller(Fp::<7>(a));
            assert_eq!(tq.unit(), Zp::<7, 4>::teichmuller(Fp::<7>(a)).0);
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
                assert!(sq.is_square(), "a square must read as a square");
                let r = sq.sqrt().expect("a square has a root");
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
        let r = sq.sqrt().expect("a square inverts");
        assert_eq!(r.mul(&r), sq);
        // p² · square ⇒ square at half valuation
        let y = Q::from_p_power(2).mul(&sq);
        let ry = y.sqrt().expect("3²·square is a square");
        assert_eq!(ry.mul(&ry), y);
        assert_eq!(ry.valuation(), Some(1));
        // odd valuation ⇒ not a square
        assert!(!Q::from_p_power(1).mul(&sq).is_square());
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
        let r = x.sqrt().unwrap();
        assert_eq!(r.mul(&r), x);
        assert!(Qq::<5, 5, 1>::from_int(2).sqrt().is_none()); // 2 a non-residue mod 5
    }

    #[test]
    fn dyadic_sqrt_is_rejected() {
        // p = 2 sqrt asserts — the dyadic case is the forms mod-8 story, not a lift.
        assert!(std::panic::catch_unwind(|| Qp::<2, 5>::from_i128(1).sqrt()).is_err());
        assert!(std::panic::catch_unwind(|| Zp::<2, 5>(1).is_square()).is_err());
    }
}
