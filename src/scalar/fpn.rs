//! Finite extension fields `F_{p^n}` — completing the field tower in every
//! characteristic.
//!
//! The odd-characteristic leg of the crate only had the *prime* fields `Fp<P>`;
//! characteristic 2 had the whole nimber tower (`F_{2^{2^k}}`). `Fpn<P, N>` closes
//! that asymmetry: it is `F_{p^n}` for any supported `(p, n)`, the odd-characteristic
//! analogue of the nimber tower. It also supplies the **char-2 odd-degree** fields
//! the nimbers cannot reach — the finite nimbers realise only `F_{2^{2^k}}` (degrees
//! that are powers of two), so `F_8` (degree 3), `F_32` (degree 5), … are *not*
//! nimber subfields; `Fpn<2, 3>` is the only way to get `F_8` here.
//!
//! ## The const-generic modulus, two parameters
//!
//! Like `Fp<P>`, the modulus lives in the **type** (`Scalar::zero()/one()` take no
//! `self`). A field is `Fpn<const P: u128, const N: usize>` = `F_{p^N}`, carried as the
//! `N` coefficients of `c_0 + c_1 x + … + c_{N-1} x^{N-1}` with each `c_i ∈ [0, P)`.
//! A different `(P, N)` is a different type — the same no-mixing discipline the rest
//! of the crate uses. `Fpn<2, 2>` is "the polynomial-basis `F_4`", a *different type*
//! from (but isomorphic to) the nimber `F_4`; the value-add over the nimbers is the
//! odd-degree char-2 layers and the odd-`p` extensions.
//!
//! ## The reduction polynomial
//!
//! Arithmetic is in `F_p[x] / (m(x))` for a monic irreducible `m` of degree `N`.
//! [`reduction`] returns the low coefficients `r` of the reduction rule
//! `x^N = Σ_i r_i x^i` (i.e. `m(x) = x^N − Σ_i r_i x^i`). The polynomials shipped here
//! are verified irreducible by the exhaustive field-axiom tests below; they can be
//! swapped for the canonical **Conway polynomials** later (which additionally give
//! compatible embeddings `F_{p^n} ↪ F_{p^{nm}}`) without touching anything else.
//! `mul` is schoolbook multiply-then-reduce — the degree-`N`, odd-`p` generalisation
//! of `onag.rs`'s "reduce mod `ω³ = 2`".

use crate::scalar::Scalar;
use std::fmt;

/// An element of `F_{p^N}`: the coefficients of `c_0 + c_1 x + … + c_{N-1} x^{N-1}`,
/// each reduced into `[0, P)`.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct Fpn<const P: u128, const N: usize>(pub [u128; N]);

/// Low coefficients `r` of the reduction rule `x^N = Σ_i r_i x^i` for the supported
/// `(P, N)` fields. Each returned slice has length `N`. Unsupported pairs are a
/// compile-time error (the `panic!` fires in a `const`-evaluable position when the
/// field is monomorphised through the engine, and at first use otherwise).
///
/// The chosen reduction polynomials (all verified irreducible by the tests):
///   * `F_4  = F_2[x]/(x²+x+1)`   → `x² = x + 1`
///   * `F_8  = F_2[x]/(x³+x+1)`   → `x³ = x + 1`
///   * `F_9  = F_3[x]/(x²+1)`     → `x² = 2`
///   * `F_25 = F_5[x]/(x²−2)`     → `x² = 2`
///   * `F_27 = F_3[x]/(x³−x+1)`   → `x³ = x + 2`
pub(crate) const fn reduction<const P: u128, const N: usize>() -> &'static [u128] {
    match (P, N) {
        (_, 1) => &[0],       // degree 1: F_p itself, no reduction needed
        (2, 2) => &[1, 1],    // x² = 1 + x
        (2, 3) => &[1, 1, 0], // x³ = 1 + x
        (3, 2) => &[2, 0],    // x² = 2
        (5, 2) => &[2, 0],    // x² = 2
        (3, 3) => &[2, 1, 0], // x³ = 2 + x
        _ => panic!("Fpn: unsupported (P, N) finite field — add its reduction polynomial"),
    }
}

impl<const P: u128, const N: usize> Fpn<P, N> {
    /// The field order `p^N`.
    pub fn order() -> u128 {
        let mut acc = 1u128;
        for _ in 0..N {
            acc = acc.checked_mul(P).expect("Fpn order exceeds u128");
        }
        acc
    }

    /// Embed a base-field constant `c ∈ F_p` as the degree-0 element.
    pub fn constant(c: u128) -> Self {
        let mut out = [0u128; N];
        if N > 0 {
            out[0] = c % P;
        }
        Fpn(out)
    }

    /// Build from a coefficient slice (low-to-high), reducing each entry mod `P`.
    /// Extra trailing coefficients beyond `N` must be zero (else it is not an
    /// element of this field); they are ignored here, so prefer length `≤ N`.
    pub fn from_coeffs(cs: &[u128]) -> Self {
        let mut out = [0u128; N];
        for (i, slot) in out.iter_mut().enumerate() {
            if i < cs.len() {
                *slot = cs[i] % P;
            }
        }
        Fpn(out)
    }

    /// Is this element a square in `F_{p^N}`? In characteristic 2 the Frobenius
    /// `x ↦ x²` is a bijection, so *every* element is a square; in odd
    /// characteristic this is Euler's criterion `x^{(q−1)/2} = 1` (with `0` a
    /// square). The square-class is the `H¹` / discriminant datum the odd-char
    /// classifier reads — so this is what lets the invariant theory run over a
    /// genuine extension field, not just a prime field.
    pub fn is_square(&self) -> bool {
        if self.is_zero() {
            return true;
        }
        if P == 2 {
            return true; // Frobenius is onto in char 2
        }
        // a^{(q−1)/2} == 1
        let mut e = (Self::order() - 1) / 2;
        let mut base = *self;
        let mut acc = Self::one();
        while e > 0 {
            if e & 1 == 1 {
                acc = acc.mul(&base);
            }
            base = base.mul(&base);
            e >>= 1;
        }
        acc == Self::one()
    }

    /// The generator `x` (the class of the indeterminate), i.e. `[0, 1, 0, …]`.
    pub fn generator() -> Self {
        let mut out = [0u128; N];
        if N > 1 {
            out[1] = 1 % P;
        } else if N == 1 {
            // degree-1: the "field" is F_p and x = 0 in it; this is a degenerate case.
            out[0] = 0;
        }
        Fpn(out)
    }
}

impl<const P: u128, const N: usize> fmt::Debug for Fpn<P, N> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut parts: Vec<String> = Vec::new();
        for i in (0..N).rev() {
            let c = self.0[i];
            if c == 0 {
                continue;
            }
            let term = match i {
                0 => format!("{c}"),
                1 if c == 1 => "x".to_string(),
                1 => format!("{c}x"),
                _ if c == 1 => format!("x^{i}"),
                _ => format!("{c}x^{i}"),
            };
            parts.push(term);
        }
        if parts.is_empty() {
            write!(f, "0")
        } else {
            write!(f, "{}", parts.join(" + "))
        }
    }
}

impl<const P: u128, const N: usize> Scalar for Fpn<P, N> {
    fn zero() -> Self {
        Fpn([0u128; N])
    }

    fn one() -> Self {
        let mut out = [0u128; N];
        if N > 0 {
            out[0] = 1 % P;
        }
        Fpn(out)
    }

    fn add(&self, rhs: &Self) -> Self {
        let mut out = [0u128; N];
        for i in 0..N {
            out[i] = ((self.0[i] as u128 + rhs.0[i] as u128) % P as u128) as u128;
        }
        Fpn(out)
    }

    fn neg(&self) -> Self {
        let mut out = [0u128; N];
        for i in 0..N {
            out[i] = if self.0[i] == 0 { 0 } else { P - self.0[i] };
        }
        Fpn(out)
    }

    fn mul(&self, rhs: &Self) -> Self {
        let p = P as u128;
        // Schoolbook product into a degree-(2N-2) scratch, then reduce mod m(x).
        let mut scratch = vec![0u128; 2 * N - 1];
        for i in 0..N {
            if self.0[i] == 0 {
                continue;
            }
            let ai = self.0[i] as u128;
            for j in 0..N {
                scratch[i + j] = (scratch[i + j] + ai * rhs.0[j] as u128) % p;
            }
        }
        // x^k = x^{k-N} · x^N = x^{k-N} · Σ_i red_i x^i, folding top down. (Degree 1 =
        // F_p needs no reduction — the scratch is already a single coefficient.)
        if N > 1 {
            let red = reduction::<P, N>();
            for k in (N..2 * N - 1).rev() {
                let c = scratch[k];
                if c == 0 {
                    continue;
                }
                scratch[k] = 0;
                for i in 0..N {
                    scratch[k - N + i] = (scratch[k - N + i] + c * red[i] as u128) % p;
                }
            }
        }
        let mut out = [0u128; N];
        for i in 0..N {
            out[i] = scratch[i] as u128;
        }
        Fpn(out)
    }

    fn characteristic() -> u128 {
        // The *characteristic* is the prime p, not the order p^N.
        P as u128
    }

    fn inv(&self) -> Option<Self> {
        if self.is_zero() {
            return None;
        }
        // Fermat: a^{p^N − 2} = a^{−1} in F_{p^N}. Square-and-multiply with `mul`.
        let mut e = Self::order() - 2;
        let mut base = *self;
        let mut result = Self::one();
        while e > 0 {
            if e & 1 == 1 {
                result = result.mul(&base);
            }
            base = base.mul(&base);
            e >>= 1;
        }
        Some(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::clifford::{CliffordAlgebra, Metric};

    /// Every element of `F_{p^N}`, enumerated by base-`P` digits.
    fn elems<const P: u128, const N: usize>() -> Vec<Fpn<P, N>> {
        let order = Fpn::<P, N>::order();
        (0..order)
            .map(|mut code| {
                let mut coeffs = [0u128; N];
                for slot in coeffs.iter_mut() {
                    *slot = (code % P as u128) as u128;
                    code /= P as u128;
                }
                Fpn(coeffs)
            })
            .collect()
    }

    fn check_field_axioms<const P: u128, const N: usize>() {
        let es = elems::<P, N>();
        let zero = Fpn::<P, N>::zero();
        let one = Fpn::<P, N>::one();
        assert_eq!(es.len(), Fpn::<P, N>::order() as usize);
        for &a in &es {
            // additive identity / inverse
            assert_eq!(a.add(&zero), a);
            assert_eq!(a.add(&a.neg()), zero);
            // multiplicative identity
            assert_eq!(a.mul(&one), a);
            // inverse: every nonzero element is a unit (THIS is what catches a
            // reducible reduction polynomial — a zero divisor would have no inverse).
            if a.is_zero() {
                assert!(a.inv().is_none());
            } else {
                let ai = a.inv().expect("nonzero element of a field is invertible");
                assert_eq!(a.mul(&ai), one, "a·a⁻¹ = 1");
            }
            for &b in &es {
                assert_eq!(a.add(&b), b.add(&a), "add commutes");
                assert_eq!(a.mul(&b), b.mul(&a), "mul commutes");
                for &c in &es {
                    assert_eq!(a.add(&b).add(&c), a.add(&b.add(&c)), "add assoc");
                    assert_eq!(a.mul(&b).mul(&c), a.mul(&b.mul(&c)), "mul assoc");
                    assert_eq!(a.mul(&b.add(&c)), a.mul(&b).add(&a.mul(&c)), "distrib");
                }
            }
        }
    }

    #[test]
    fn field_axioms_f4_f8_f9_f25_f27() {
        check_field_axioms::<2, 2>(); // F_4
        check_field_axioms::<2, 3>(); // F_8
        check_field_axioms::<3, 2>(); // F_9
        check_field_axioms::<5, 2>(); // F_25
        check_field_axioms::<3, 3>(); // F_27
    }

    #[test]
    fn characteristic_is_p_not_order() {
        assert_eq!(Fpn::<2, 3>::characteristic(), 2); // F_8 has characteristic 2
        assert_eq!(Fpn::<2, 3>::order(), 8);
        assert_eq!(Fpn::<3, 3>::characteristic(), 3); // F_27 has characteristic 3
        assert_eq!(Fpn::<3, 3>::order(), 27);
    }

    #[test]
    fn generator_satisfies_its_minimal_polynomial() {
        // F_8: x³ = x + 1, so x³ + x + 1 = 0 (and −1 = 1 in char 2 ⇒ x³ = x + 1).
        let x = Fpn::<2, 3>::generator();
        let x3 = x.mul(&x).mul(&x);
        assert_eq!(x3, Fpn::<2, 3>::from_coeffs(&[1, 1, 0])); // x + 1
                                                              // F_27: x³ = x + 2.
        let y = Fpn::<3, 3>::generator();
        let y3 = y.mul(&y).mul(&y);
        assert_eq!(y3, Fpn::<3, 3>::from_coeffs(&[2, 1, 0])); // x + 2
    }

    #[test]
    fn frobenius_is_an_automorphism() {
        // x ↦ x^p is additive (the Frobenius) in characteristic p.
        let pow_p = |a: Fpn<3, 3>| {
            let mut r = Fpn::<3, 3>::one();
            for _ in 0..3 {
                r = r.mul(&a);
            }
            r
        };
        for a in elems::<3, 3>() {
            for b in elems::<3, 3>() {
                assert_eq!(pow_p(a.add(&b)), pow_p(a).add(&pow_p(b)));
            }
        }
    }

    #[test]
    fn clifford_over_f9_monomorphises() {
        // Cl over F_9 with q = [x, 1]: the engine runs on the extension field exactly
        // as on a prime field; antisymmetry signs are genuine (−1 = 2 in F_3 ⊂ F_9).
        let x = Fpn::<3, 2>::generator();
        let one = Fpn::<3, 2>::one();
        let alg = CliffordAlgebra::new(2, Metric::diagonal(vec![x, one]));
        let (e0, e1) = (alg.gen(0), alg.gen(1));
        assert_eq!(alg.mul(&e0, &e0), alg.scalar(x));
        assert_eq!(alg.mul(&e1, &e1), alg.scalar(one));
        // e0 e1 = −(e1 e0)
        let neg_one = Fpn::<3, 2>::one().neg();
        assert_eq!(
            alg.mul(&e0, &e1),
            alg.scalar_mul(&neg_one, &alg.mul(&e1, &e0))
        );
    }
}
