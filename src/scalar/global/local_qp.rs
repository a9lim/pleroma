//! A **runtime-prime** `p`-adic cell — the gap-filler that makes the
//! [`Adele`](crate::scalar::Adele) ring possible.
//!
//! Every other p-adic backend ([`Qp`](crate::scalar::Qp), [`Zp`](crate::scalar::Zp),
//! [`Qq`](crate::scalar::Qq)) carries its prime as a **const generic** — so `Qp<3,_>`
//! and `Qp<5,_>` are distinct *types* and cannot share a `BTreeMap`. An adele is
//! indexed by a runtime prime, so it needs a cell whose prime is a *field*, not a
//! type parameter. `LocalQp` is exactly that: a near-verbatim port of `Qp`'s
//! capped-relative arithmetic with `(p, k)` moved into the struct.
//!
//! Because its `(p, k)` arrive at construction, `LocalQp` deliberately does **not**
//! implement [`Scalar`](crate::scalar::Scalar) — that trait's zero-argument
//! `zero()`/`one()` can't supply them. It carries the same operations as inherent
//! methods (`add`/`neg`/`mul`/`inv`/`is_zero`), and mixing two different primes is a
//! bug (`assert`-guarded), exactly as the rest of the table forbids mixing scalar
//! worlds.
//!
//! Precision is the same **capped-relative** model as `Qp` (mul/inv exact, addition
//! non-associative across precision boundaries) — see `scalar/small/qp.rs`.

use crate::scalar::{is_prime_u128, mod_inverse_u128, Rational};
use std::fmt;

/// `p^e`, checked against `u128` overflow.
fn p_pow(p: u128, e: u128) -> u128 {
    let mut acc = 1u128;
    for _ in 0..e {
        acc = acc.checked_mul(p).expect("LocalQp: p-power exceeds u128");
    }
    acc
}

/// An element of `Q_p` to precision `k`, with `p` and `k` carried at runtime:
/// `p^{val} · unit` with `p ∤ unit` mod `p^k`, or the sentinel for the field zero.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct LocalQp {
    p: u128,
    k: u128,
    /// The unit mantissa in `[0, p^k)` (`0` only for the field zero).
    unit: u128,
    /// The signed `p`-adic valuation (`0` for the field zero).
    val: i128,
}

impl LocalQp {
    fn check(p: u128, k: u128) {
        assert!(
            is_prime_u128(p) && k > 0,
            "LocalQp needs prime p and positive precision k, got p={p}, k={k}"
        );
        let mut acc = 1u128;
        for _ in 0..k {
            acc = acc.checked_mul(p).expect("LocalQp modulus exceeds u128");
            assert!(
                acc <= i128::MAX as u128,
                "LocalQp modulus must fit i128-backed embeddings, got p={p}, k={k}"
            );
        }
    }

    fn same_world(&self, other: &LocalQp) {
        assert!(
            self.p == other.p && self.k == other.k,
            "LocalQp: cannot mix primes/precisions ({},{}) vs ({},{})",
            self.p,
            self.k,
            other.p,
            other.k
        );
    }

    /// The mantissa modulus `p^k`.
    pub fn modulus(&self) -> u128 {
        p_pow(self.p, self.k)
    }

    /// The prime.
    pub fn prime(&self) -> u128 {
        self.p
    }

    /// The relative precision `k`.
    pub fn precision(&self) -> u128 {
        self.k
    }

    /// Build `p^{val} · unit`, normalizing: factor `p` out of `unit` into the
    /// valuation, reduce mod `p^k`.
    fn normalized(p: u128, k: u128, unit_raw: u128, val: i128) -> Self {
        let m = p_pow(p, k);
        let mut u = unit_raw % m;
        if u == 0 {
            return LocalQp {
                p,
                k,
                unit: 0,
                val: 0,
            };
        }
        let mut v = val;
        while u.is_multiple_of(p) {
            u /= p;
            v += 1;
        }
        LocalQp {
            p,
            k,
            unit: u,
            val: v,
        }
    }

    /// The field zero of `Q_p` at precision `k`.
    pub fn zero(p: u128, k: u128) -> Self {
        Self::check(p, k);
        LocalQp {
            p,
            k,
            unit: 0,
            val: 0,
        }
    }

    /// The field one.
    pub fn one(p: u128, k: u128) -> Self {
        Self::check(p, k);
        LocalQp {
            p,
            k,
            unit: 1 % p_pow(p, k),
            val: 0,
        }
    }

    /// Embed a signed integer, extracting its `p`-adic valuation.
    pub fn from_i128(p: u128, k: u128, n: i128) -> Self {
        Self::check(p, k);
        if n == 0 {
            return LocalQp {
                p,
                k,
                unit: 0,
                val: 0,
            };
        }
        let pp = p as i128;
        let mut w = 0i128;
        let mut nn = n;
        while nn % pp == 0 {
            nn /= pp;
            w += 1;
        }
        let m = p_pow(p, k) as i128;
        let unit = (((nn % m) + m) % m) as u128;
        LocalQp { p, k, unit, val: w }
    }

    /// `p^v`, mantissa `1`. `from_p_power(p, k, -1)` is `1/p`.
    pub fn from_p_power(p: u128, k: u128, v: i128) -> Self {
        Self::check(p, k);
        LocalQp {
            p,
            k,
            unit: 1 % p_pow(p, k),
            val: v,
        }
    }

    /// Embed a rational into `Q_p`: `from_i128(num) · from_i128(den)^{-1}`. The
    /// valuation is `v_p(num) − v_p(den)`.
    pub fn from_rational(p: u128, k: u128, q: &Rational) -> Self {
        let num = LocalQp::from_i128(p, k, q.numer());
        let den = LocalQp::from_i128(p, k, q.denom());
        // den > 0 ⇒ nonzero ⇒ invertible in the field.
        num.mul(
            &den.inv()
                .expect("LocalQp::from_rational: nonzero denominator"),
        )
    }

    /// The valuation, or `None` for zero.
    pub fn valuation(&self) -> Option<i128> {
        if self.unit == 0 {
            None
        } else {
            Some(self.val)
        }
    }

    /// The unit mantissa in `[0, p^k)`.
    pub fn unit(&self) -> u128 {
        self.unit
    }

    /// Whether this is the field zero.
    pub fn is_zero(&self) -> bool {
        self.unit == 0
    }

    /// Addition (capped-relative, like `Qp`).
    pub fn add(&self, rhs: &Self) -> Self {
        self.same_world(rhs);
        if self.unit == 0 {
            return *rhs;
        }
        if rhs.unit == 0 {
            return *self;
        }
        let m = self.modulus();
        let (lo, hi) = if self.val <= rhs.val {
            (self, rhs)
        } else {
            (rhs, self)
        };
        let d = (hi.val - lo.val) as u128;
        let shifted = if d >= self.k {
            0
        } else {
            crate::scalar::mul_mod_u128(p_pow(self.p, d), hi.unit, m)
        };
        let b = lo
            .unit
            .checked_add(shifted)
            .expect("LocalQp addition mantissa sum exceeds u128")
            % m;
        if b == 0 {
            return LocalQp {
                p: self.p,
                k: self.k,
                unit: 0,
                val: 0,
            };
        }
        Self::normalized(self.p, self.k, b, lo.val)
    }

    /// Negation.
    pub fn neg(&self) -> Self {
        if self.unit == 0 {
            return *self;
        }
        LocalQp {
            p: self.p,
            k: self.k,
            unit: self.modulus() - self.unit,
            val: self.val,
        }
    }

    /// Multiplication (exact: valuations add, mantissa is a genuine unit).
    pub fn mul(&self, rhs: &Self) -> Self {
        self.same_world(rhs);
        if self.unit == 0 || rhs.unit == 0 {
            return LocalQp {
                p: self.p,
                k: self.k,
                unit: 0,
                val: 0,
            };
        }
        let m = self.modulus();
        LocalQp {
            p: self.p,
            k: self.k,
            // mul_mod_u128, not checked_mul: p^k can approach i128::MAX, so a
            // schoolbook unit×unit product overflows u128 on in-range inputs.
            unit: crate::scalar::mul_mod_u128(self.unit, rhs.unit, m),
            val: self
                .val
                .checked_add(rhs.val)
                .expect("LocalQp multiplication valuation exceeds i128"),
        }
    }

    /// Multiplicative inverse — total on nonzero (the field property).
    pub fn inv(&self) -> Option<Self> {
        if self.unit == 0 {
            return None;
        }
        let uinv = mod_inverse_u128(self.unit, self.modulus())?;
        Some(LocalQp {
            p: self.p,
            k: self.k,
            unit: uinv,
            val: -self.val,
        })
    }
}

impl fmt::Display for LocalQp {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.unit == 0 {
            return write!(f, "0 (Q_{})", self.p);
        }
        if self.val == 0 {
            write!(f, "{} (mod {}^{})", self.unit, self.p, self.k)
        } else {
            write!(
                f,
                "{}·{}^{} (mod {}^{})",
                self.unit, self.p, self.val, self.p, self.k
            )
        }
    }
}

impl fmt::Debug for LocalQp {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(self, f)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scalar::{Qp, Scalar};

    // The oracle: LocalQp must agree element-for-element with the trusted
    // const-generic Qp<P,K>.
    macro_rules! oracle {
        ($P:literal, $K:literal) => {{
            let p: u128 = $P;
            let k: u128 = $K;
            for n in -40i128..=40 {
                let q = Qp::<$P, $K>::from_i128(n);
                let l = LocalQp::from_i128(p, k, n);
                assert_eq!(q.valuation(), l.valuation(), "val from_i128 {n}");
                assert_eq!(q.unit(), l.unit(), "unit from_i128 {n}");
            }
            for a in -20i128..=20 {
                for b in -20i128..=20 {
                    let (qa, qb) = (Qp::<$P, $K>::from_i128(a), Qp::<$P, $K>::from_i128(b));
                    let (la, lb) = (LocalQp::from_i128(p, k, a), LocalQp::from_i128(p, k, b));
                    let qs = qa.add(&qb);
                    let ls = la.add(&lb);
                    assert_eq!(qs.valuation(), ls.valuation(), "val {a}+{b}");
                    assert_eq!(qs.unit(), ls.unit(), "unit {a}+{b}");
                    let qm = qa.mul(&qb);
                    let lm = la.mul(&lb);
                    assert_eq!(qm.valuation(), lm.valuation(), "val {a}*{b}");
                    assert_eq!(qm.unit(), lm.unit(), "unit {a}*{b}");
                    if a != 0 {
                        let qi = qa.inv().unwrap();
                        let li = la.inv().unwrap();
                        assert_eq!(qi.valuation(), li.valuation(), "val 1/{a}");
                        assert_eq!(qi.unit(), li.unit(), "unit 1/{a}");
                    }
                }
            }
        }};
    }

    #[test]
    fn matches_qp_oracle_p3() {
        oracle!(3, 3);
    }
    #[test]
    fn matches_qp_oracle_p5() {
        oracle!(5, 4);
    }
    #[test]
    fn matches_qp_oracle_p2() {
        oracle!(2, 6);
    }

    #[test]
    fn one_over_p_and_field_property() {
        let p = LocalQp::from_i128(7, 4, 7);
        let pinv = p.inv().unwrap();
        assert_eq!(pinv.valuation(), Some(-1));
        assert_eq!(p.mul(&pinv), LocalQp::one(7, 4));
        assert_eq!(LocalQp::zero(7, 4).inv(), None);
    }

    #[test]
    fn from_rational_valuation() {
        // 50/3 in Q_5: v_5(50) = 2, v_5(3) = 0 ⇒ valuation 2.
        let x = LocalQp::from_rational(5, 4, &Rational::new(50, 3));
        let xq = Qp::<5, 4>::from_rational(&Rational::new(50, 3));
        assert_eq!(x.unit(), xq.unit());
        assert_eq!(x.valuation(), xq.valuation());
        assert_eq!(x.valuation(), Some(2));
        // 3/50 ⇒ valuation -2, and it is the inverse of the above.
        let y = LocalQp::from_rational(5, 4, &Rational::new(3, 50));
        let yq = Qp::<5, 4>::from_rational(&Rational::new(3, 50));
        assert_eq!(y.unit(), yq.unit());
        assert_eq!(y.valuation(), yq.valuation());
        assert_eq!(y.valuation(), Some(-2));
        assert_eq!(x.mul(&y), LocalQp::one(5, 4));
    }

    #[test]
    #[should_panic(expected = "needs prime p")]
    fn invalid_runtime_world_is_rejected_in_release_too() {
        let _ = LocalQp::one(4, 3);
    }

    #[test]
    #[should_panic(expected = "modulus must fit")]
    fn oversized_runtime_modulus_is_rejected() {
        let _ = LocalQp::one(2, 127);
    }

    #[test]
    #[should_panic(expected = "cannot mix primes")]
    fn mixed_runtime_worlds_are_rejected_in_release_too() {
        let x = LocalQp::one(3, 4);
        let y = LocalQp::one(5, 4);
        let _ = x.add(&y);
    }
}
