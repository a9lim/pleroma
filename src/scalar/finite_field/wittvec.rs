//! Witt vectors `W_N(F_q)` — the canonical characteristic-`p` → characteristic-`0`
//! lift, as a `Scalar` backend.
//!
//! **Disambiguation:** this is the ring of *Witt vectors*, unrelated to the
//! quadratic-form *Witt group* in [`forms::witt`](crate::forms::witt) /
//! [`WittClass`](crate::forms::WittClass).
//!
//! The length-`N` `p`-typical Witt vectors over `F_q` form a finite commutative
//! ring of order `q^N`, the degree-`N` truncation of `W(F_q)`, the ring of integers
//! of the **unramified** extension of `Q_p` with residue field `F_q`. That gives an
//! exact, manifestly-correct model:
//!
//! > `W_N(F_q) ≅ (Z/p^N)[t] / (f̃(t))`,
//!
//! where `f̃` is any monic lift to `Z` of the `F_p`-irreducible polynomial defining
//! `F_q` (Hensel: `f̃` stays irreducible over `Z_p`, and the extension is unramified
//! because `f̃ mod p` is separable). We reuse [`Fpn`]'s reduction polynomial as that
//! lift. So a Witt vector is a degree-`F` polynomial over `Z/p^N` — arithmetic is
//! exactly `Fpn`'s, with the coefficient field `F_p` swapped for the coefficient
//! *ring* `Z/p^N`. This avoids hand-deriving the ghost-inversion (Witt addition)
//! polynomials, whose division-by-`p` is the classic correctness trap.
//!
//! The layer accessors ([`WittVec::witt_components`],
//! [`WittVec::from_witt_components`]) use the Teichmüller expansion
//! `a = Σ τ(x_i)·p^i`. Over `F_p` these are the classical `p`-typical Witt
//! coordinates. Over `F_{p^f}` with `f > 1` they are the Frobenius-twisted
//! Teichmüller digits for this unramified-ring representation: the classical
//! coordinates differ by the usual `p^{-i}` Frobenius twists. The round-trip
//! tests pin this representation without claiming the untwisted coordinate
//! polynomials.
//!
//! ## On-brand hook: Artin–Schreier–Witt
//!
//! `W(F₂)` length-`N` is `Z/2^N`, and the additive Frobenius/`℘` on Witt vectors
//! generalises the `y² + y = c` Artin–Schreier solver in `scalar::nim_solve_artin_schreier` to
//! `Z/p^n`-extensions (Artin–Schreier–Witt theory) — extending the Arf↔Artin–
//! Schreier thread. (Documented as motivation; the solver itself is future work.)

use crate::scalar::finite_field::fpn::reduction;
use crate::scalar::{add_mod_u128, mul_mod_u128, reduce_i128_mod_u128, Fpn, Scalar};
use std::fmt;

/// A length-`N` `p`-typical Witt vector over `F_q = F_{p^F}`, realised in
/// unramified-ring coordinates: the `F` coefficients of a polynomial over `Z/p^N`,
/// each a residue in `[0, p^N)`.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct WittVec<const P: u128, const N: usize, const F: usize>(pub [u128; F]);

impl<const P: u128, const N: usize, const F: usize> WittVec<P, N, F> {
    /// The coefficient-ring modulus `p^N` (the precision).
    pub fn modulus() -> u128 {
        let mut acc = 1u128;
        for _ in 0..N {
            acc = acc.checked_mul(P).expect("WittVec modulus exceeds u128");
        }
        acc
    }

    /// The residue field order `q = p^F`.
    pub fn residue_order() -> u128 {
        let mut acc = 1u128;
        for _ in 0..F {
            acc = acc
                .checked_mul(P)
                .expect("WittVec residue order exceeds u128");
        }
        acc
    }

    /// Embed a `Z/p^N` integer as the degree-0 (constant) Witt vector.
    pub fn from_int(n: i128) -> Self {
        let m = Self::modulus();
        let mut out = [0u128; F];
        if F > 0 {
            out[0] = reduce_i128_mod_u128(n, m);
        }
        WittVec(out)
    }

    /// Reduce mod `p`: the **residue** in `F_q` (the ghost/Teichmüller bottom layer).
    pub fn residue(&self) -> Fpn<P, F> {
        let mut c = [0u128; F];
        for i in 0..F {
            c[i] = self.0[i] % P;
        }
        Fpn::from_coeffs(&c)
    }

    /// Multiply two ring elements (polynomials over `Z/p^N` mod the lifted `f̃`).
    fn ring_mul(a: &[u128; F], b: &[u128; F]) -> [u128; F] {
        let m = Self::modulus();
        let mut scratch = vec![0u128; 2 * F - 1];
        for i in 0..F {
            if a[i] == 0 {
                continue;
            }
            let ai = a[i];
            for j in 0..F {
                let term = mul_mod_u128(ai, b[j], m);
                scratch[i + j] = add_mod_u128(scratch[i + j], term, m);
            }
        }
        // reduce mod f̃ (the same reduction poly as Fpn<P,F>, lifted to Z/p^N).
        if F > 1 {
            let red = reduction::<P, F>();
            for k in (F..2 * F - 1).rev() {
                let c = scratch[k];
                if c == 0 {
                    continue;
                }
                scratch[k] = 0;
                for i in 0..F {
                    let term = mul_mod_u128(c, red[i], m);
                    scratch[k - F + i] = add_mod_u128(scratch[k - F + i], term, m);
                }
            }
        }
        let mut out = [0u128; F];
        for i in 0..F {
            out[i] = scratch[i] as u128;
        }
        out
    }

    /// `self^e` in the ring, by square-and-multiply.
    fn pow(&self, mut e: u128) -> Self {
        let mut base = *self;
        let mut acc = Self::one();
        while e > 0 {
            if e & 1 == 1 {
                acc = acc.mul(&base);
            }
            base = base.mul(&base);
            e >>= 1;
        }
        acc
    }

    /// The **Teichmüller representative** `τ(x) ∈ W_N(F_q)` of `x ∈ F_q`: the unique
    /// multiplicative lift with `τ(x) ≡ x mod p`. Computed as `x̃^{q^{N-1}}` for any
    /// lift `x̃` (the power iteration stabilises modulo `p^N`).
    pub fn teichmuller(x: Fpn<P, F>) -> Self {
        let mut y = WittVec::<P, N, F>(x.into_coeffs()); // naive lift (coeffs already in [0,p))
        for _ in 0..N.saturating_sub(1) {
            y = y.pow(Self::residue_order());
        }
        y
    }

    /// Build a Witt vector from its Teichmüller digits `(x₀,…,x_{N-1}) ∈ F_q^N`:
    /// `Σ_i τ(x_i)·p^i`. The inverse of [`witt_components`](Self::witt_components).
    pub fn from_witt_components(xs: &[Fpn<P, F>]) -> Self {
        assert_eq!(xs.len(), N, "need exactly N Teichmuller digits");
        let mut acc = Self::zero();
        let mut pk = Self::one();
        let p_elt = Self::from_int(P as i128);
        for (i, &x) in xs.iter().enumerate() {
            acc = acc.add(&Self::teichmuller(x).mul(&pk));
            if i + 1 < N {
                pk = pk.mul(&p_elt);
            }
        }
        acc
    }

    /// The `p`-adic valuation, capped at the precision `N` (`v_p(0)` reads as `N`,
    /// the precision floor — exactly the [`Zp`](crate::scalar::Zp) discipline). A
    /// unit (residue `≠ 0`) has valuation `0`.
    pub fn p_valuation(&self) -> usize {
        let mut w = *self;
        let mut v = 0;
        while v < N && w.residue().is_zero() {
            w = w.divide_by_p();
            v += 1;
        }
        v
    }

    /// If this vector is divisible by `p` (residue `0`), divide by `p` (precision
    /// drops one digit); else `None`. The hook the fraction field [`Qq`](crate::scalar::Qq)
    /// uses to peel a Witt unit out of an arbitrary vector.
    pub fn try_divide_by_p(&self) -> Option<Self> {
        if self.residue().is_zero() {
            Some(self.divide_by_p())
        } else {
            None
        }
    }

    /// Divide by `p` an element all of whose coefficients are `≡ 0 mod p`.
    fn divide_by_p(&self) -> Self {
        let mut out = [0u128; F];
        for i in 0..F {
            debug_assert_eq!(self.0[i] % P, 0, "divide_by_p on a non-divisible element");
            out[i] = self.0[i] / P;
        }
        WittVec(out)
    }

    /// The Teichmüller digits `(x₀,…,x_{N-1}) ∈ F_q^N`: peel the layers
    /// `x_i = ((a − Σ_{j<i} τ(x_j)p^j)/p^i) mod p`. Over extension residue
    /// fields these are Frobenius-twisted relative to the classical `p`-typical
    /// Witt coordinates.
    pub fn witt_components(&self) -> Vec<Fpn<P, F>> {
        let mut a = *self;
        let mut out = Vec::with_capacity(N);
        for _ in 0..N {
            let r = a.residue();
            out.push(r);
            let t = Self::teichmuller(r);
            a = a.sub(&t).divide_by_p();
        }
        out
    }
}

impl<const P: u128, const N: usize, const F: usize> fmt::Display for WittVec<P, N, F> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // unramified-ring coordinates: coefficients of 1, t, …, t^{F-1} over Z/p^N
        write!(f, "W_{}(F_{}^{}){:?}", N, P, F, self.0)
    }
}

impl<const P: u128, const N: usize, const F: usize> fmt::Debug for WittVec<P, N, F> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(self, f)
    }
}

impl<const P: u128, const N: usize, const F: usize> Scalar for WittVec<P, N, F> {
    fn zero() -> Self {
        WittVec([0u128; F])
    }

    fn one() -> Self {
        let mut out = [0u128; F];
        if F > 0 {
            out[0] = 1 % Self::modulus();
        }
        WittVec(out)
    }

    fn add(&self, rhs: &Self) -> Self {
        let m = Self::modulus();
        let mut out = [0u128; F];
        for i in 0..F {
            out[i] = add_mod_u128(self.0[i], rhs.0[i], m);
        }
        WittVec(out)
    }

    fn neg(&self) -> Self {
        let m = Self::modulus();
        let mut out = [0u128; F];
        for i in 0..F {
            out[i] = if self.0[i] == 0 { 0 } else { m - self.0[i] };
        }
        WittVec(out)
    }

    fn mul(&self, rhs: &Self) -> Self {
        WittVec(Self::ring_mul(&self.0, &rhs.0))
    }

    fn characteristic() -> u128 {
        // The length-N truncation W_N(F_q) is a finite quotient of the
        // characteristic-0 Witt ring, with p^N · 1 = 0 and no smaller positive
        // multiple of 1 vanishing.
        Self::modulus()
    }

    fn inv(&self) -> Option<Self> {
        // Local ring: a unit iff the residue is a unit in F_q (residue ≠ 0). Invert
        // by Hensel/Newton lifting from the residue inverse: b ← b(2 − a·b) doubles
        // the precision each step. `None` for non-units (the Omnific discipline).
        let r = self.residue();
        let rinv = r.inv()?; // None iff residue is 0 ⇒ non-unit
        let mut b = WittVec::<P, N, F>(rinv.into_coeffs()); // lift of the F_q inverse
        let two = Self::one().add(&Self::one());
        let mut prec = 1usize;
        while prec < N {
            // b ← b·(2 − a·b)
            let correction = two.sub(&self.mul(&b));
            b = b.mul(&correction);
            prec *= 2;
        }
        Some(b)
    }
    /// Faster direct construction; semantically identical to the default double-and-add.
    fn from_int(n: i128) -> Self {
        WittVec::<P, N, F>::from_int(n)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scalar::Zp;

    #[test]
    fn q_equals_p_is_z_mod_pn() {
        // ORACLE: W_N(F_p) ≅ Z/p^N — ring ops must match the Zp backend exactly.
        // W_3(F_2) vs Z/8.
        for a in 0..8u128 {
            for b in 0..8u128 {
                let (wa, wb) = (WittVec::<2, 3, 1>([a]), WittVec::<2, 3, 1>([b]));
                let (za, zb) = (Zp::<2, 3>(a), Zp::<2, 3>(b));
                assert_eq!(wa.add(&wb).0[0], za.add(&zb).0);
                assert_eq!(wa.mul(&wb).0[0], za.mul(&zb).0);
                assert_eq!(wa.neg().0[0], za.neg().0);
                assert_eq!(wa.inv().map(|w| w.0[0]), za.inv().map(|z| z.0));
            }
        }
    }

    #[test]
    fn ring_axioms_w2_f4() {
        // W_2(F_4): order 4² = 16, the truncated unramified quadratic extension of Z_2.
        let elems: Vec<WittVec<2, 2, 2>> = (0..16u128)
            .map(|code| WittVec([code & 3, (code >> 2) & 3]))
            .collect();
        let one = WittVec::<2, 2, 2>::one();
        for &a in &elems {
            assert_eq!(a.mul(&one), a);
            assert_eq!(a.add(&a.neg()), WittVec::zero());
            // unit iff residue ≠ 0 in F_4
            if a.residue() != Fpn::<2, 2>::zero() {
                let ai = a.inv().expect("residue-unit must invert");
                assert_eq!(a.mul(&ai), one);
            } else {
                assert!(a.inv().is_none());
            }
            for &b in &elems {
                assert_eq!(a.add(&b), b.add(&a));
                assert_eq!(a.mul(&b), b.mul(&a));
                for &c in &elems {
                    assert_eq!(a.mul(&b).mul(&c), a.mul(&b.mul(&c)));
                    assert_eq!(a.mul(&b.add(&c)), a.mul(&b).add(&a.mul(&c)));
                }
            }
        }
    }

    #[test]
    fn residue_is_the_field_fq() {
        // Reducing W_N(F_q) mod p recovers F_q = Fpn<P,F>, while the finite
        // truncated ring itself has characteristic p^N.
        assert_eq!(
            WittVec::<3, 2, 2>([4, 7]).residue(),
            Fpn::<3, 2>::from_coeffs(&[1, 1])
        );
        assert_eq!(WittVec::<2, 3, 1>::characteristic(), 8);
    }

    #[test]
    fn witt_coordinate_roundtrip() {
        // from_witt_components ∘ witt_components = id, over W_3(F_2) = Z/8.
        for code in 0..8u128 {
            let w = WittVec::<2, 3, 1>([code]);
            let comps = w.witt_components();
            assert_eq!(WittVec::<2, 3, 1>::from_witt_components(&comps), w);
        }
        // and over W_2(F_4).
        for code in 0..16u128 {
            let w = WittVec::<2, 2, 2>([code & 3, (code >> 2) & 3]);
            let comps = w.witt_components();
            assert_eq!(WittVec::<2, 2, 2>::from_witt_components(&comps), w);
        }
    }

    #[test]
    fn large_modulus_arithmetic_reduces_without_wrapping() {
        type W = WittVec<5, 55, 1>;
        let m = W::modulus();
        assert!(m > i128::MAX as u128);
        let minus_one = WittVec::<5, 55, 1>([m - 1]);
        assert_eq!(W::from_int(-1), minus_one);
        assert_eq!(minus_one.add(&minus_one), WittVec::<5, 55, 1>([m - 2]));
        assert_eq!(minus_one.mul(&minus_one), WittVec::<5, 55, 1>([1]));
    }

    #[test]
    fn p2_witt_addition_carry_formula() {
        // The classical p=2 Witt addition: z₀ = x₀+y₀, z₁ = x₁+y₁−x₀y₀ (in F₂,
        // −1=1 so z₁ = x₁+y₁+x₀y₀). Verified against the ring sum, pinning the
        // ghost-coordinate semantics without hand-deriving the polynomial.
        for x0 in 0..2u128 {
            for x1 in 0..2u128 {
                for y0 in 0..2u128 {
                    for y1 in 0..2u128 {
                        let a = WittVec::<2, 2, 1>::from_witt_components(&[
                            Fpn::<2, 1>::from_coeffs(&[x0]),
                            Fpn::<2, 1>::from_coeffs(&[x1]),
                        ]);
                        let b = WittVec::<2, 2, 1>::from_witt_components(&[
                            Fpn::<2, 1>::from_coeffs(&[y0]),
                            Fpn::<2, 1>::from_coeffs(&[y1]),
                        ]);
                        let z = a.add(&b).witt_components();
                        assert_eq!(z[0].coeff(0), (x0 + y0) % 2, "z₀ = x₀+y₀");
                        assert_eq!(z[1].coeff(0), (x1 + y1 + x0 * y0) % 2, "z₁ = x₁+y₁+x₀y₀");
                    }
                }
            }
        }
    }

    #[test]
    fn extension_components_are_teichmuller_digits() {
        type F4 = Fpn<2, 2>;
        type W = WittVec<2, 2, 2>;
        let zero = F4::zero();
        let one = F4::one();
        let t = F4::from_coeffs(&[0, 1]);
        let sqrt_t = F4::from_coeffs(&[1, 1]);
        let a = W::from_witt_components(&[one, zero]);
        let b = W::from_witt_components(&[t, zero]);
        let z = a.add(&b).witt_components();
        assert_eq!(z[0], one.add(&t));
        assert_eq!(
            z[1], sqrt_t,
            "the second Teichmuller digit carries sqrt(x0*y0), not x0*y0"
        );
        assert_ne!(z[1], t);
    }
}
