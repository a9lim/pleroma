//! The Witt group of quadratic forms over a nim-field — the abstraction that
//! sits behind the `A ⊕ A ≅ H ⊕ H` fact `arf.rs` checks pointwise.
//!
//! Two nonsingular quadratic forms are **Witt-equivalent** if they become
//! isomorphic after adding hyperbolic planes; the equivalence classes form an
//! abelian group `W_q(F)` under orthogonal sum `⊥`, with the hyperbolic plane as
//! identity. Over a *finite* field of characteristic 2 the anisotropic forms are
//! just two — the zero form (Arf 0) and the unique anisotropic plane (Arf 1) —
//! so `W_q(F_{2^m}) ≅ ℤ/2`, **classified completely by the Arf invariant**, and
//! the group law is XOR of Arf invariants. (Over the full algebraically-closed
//! On₂, or other fields, `W_q` can be richer; for the finite nim-subfields this
//! engine targets, Arf is the whole story.)
//!
//! So `WittClass` makes the additivity executable as a group: `w(A) + w(A) = 0`
//! is the same statement as `A ⊕ A ≅ H ⊕ H`, now a one-liner.

use crate::clifford::Metric;
use crate::forms::arf_invariant;
use crate::scalar::Nimber;

/// A class in the Witt group `W_q(F) ≅ ℤ/2` of a finite nim-field: the Arf
/// invariant of a form's anisotropic core (hyperbolic planes are the identity).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct WittClass {
    /// The class, 0 or 1 — equivalently the Arf invariant of the nonsingular core.
    pub arf: u8,
}

impl WittClass {
    /// The identity: the class of the hyperbolic plane (and of the zero form).
    pub fn zero() -> Self {
        WittClass { arf: 0 }
    }

    /// The Witt class of a nimber Clifford metric — the Arf invariant of its
    /// nonsingular core. (Arf is a Witt invariant: it ignores hyperbolic summands
    /// and the polar-form radical, so this is well defined on the class.)
    pub fn from_metric(metric: &Metric<Nimber>) -> Self {
        WittClass {
            arf: arf_invariant(metric).arf,
        }
    }

    /// The group operation: the class of the orthogonal sum `⊥` of two forms.
    /// Arf is additive, hyperbolics vanish, so this is XOR of the Arf invariants.
    pub fn add(&self, other: &WittClass) -> WittClass {
        WittClass {
            arf: self.arf ^ other.arf,
        }
    }

    /// In `ℤ/2` every element is its own inverse (`w + w = 0`).
    pub fn neg(&self) -> WittClass {
        *self
    }

    /// Whether this is the identity class — i.e. the form is hyperbolic (its
    /// anisotropic core is zero).
    pub fn is_hyperbolic(&self) -> bool {
        self.arf == 0
    }

    /// Dimension of the anisotropic core: 0 (hyperbolic) or 2 (the plane).
    pub fn anisotropic_dim(&self) -> usize {
        if self.arf == 0 {
            0
        } else {
            2
        }
    }

    pub fn display(&self) -> String {
        if self.arf == 0 {
            "0 (hyperbolic class)".to_string()
        } else {
            "[anisotropic plane] (Arf 1)".to_string()
        }
    }
}

/// The Witt class across **all three characteristics** — the group-theoretic
/// home of the classifier trichotomy (char-0 signature / odd-char
/// discriminant / char-2 Arf), mirroring the Artin–Schreier↔Arf unification.
///
/// * `Char0`: over a real-closed field `W(ℝ) ≅ ℤ`, classified by the signature
///   `p − q` (the only invariant; the group law adds signatures).
/// * `OddChar`: over a finite field `F_q` of odd characteristic `W(F_q)` has
///   order 4. Its invariants are `e0 = dim mod 2` and `sclass` = the
///   **signed discriminant** `(−1)^{m(m−1)/2}·det` mod squares (a genuine Witt
///   invariant, unlike the ordinary det when `−1` is a nonsquare). The group is
///   `ℤ/4` when `−1` is a nonsquare (`q ≡ 3 mod 4`, `kappa = 1`) and `ℤ/2 × ℤ/2`
///   when `−1` is a square (`q ≡ 1 mod 4`, `kappa = 0`). The `(−1)^{mn}`
///   correction in the signed-discriminant sum is exactly the `kappa` term in
///   `add`, and is what produces the `ℤ/4` when `kappa = 1`.
/// * `Char2`: over a finite nim-field `W ≅ ℤ/2`, classified by the Arf invariant
///   (this is the existing `WittClass`, re-homed as a variant).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WittClassG {
    Char0 {
        signature: isize,
    },
    OddChar {
        /// nonsquareness of `−1`: 0 if `−1` is a square (`q≡1 mod 4`), else 1.
        kappa: u8,
        /// dimension mod 2.
        e0: u8,
        /// signed-discriminant square-class: 0 if a square, 1 if a nonsquare.
        sclass: u8,
    },
    Char2 {
        arf: u8,
    },
}

impl WittClassG {
    /// Char-0 Witt class from a signature `(p, q)`.
    pub fn char0(p: usize, q: usize) -> Self {
        WittClassG::Char0 {
            signature: p as isize - q as isize,
        }
    }

    /// Char-2 Witt class from a nimber metric (the Arf invariant).
    pub fn char2_from_metric(metric: &Metric<Nimber>) -> Self {
        WittClassG::Char2 {
            arf: arf_invariant(metric).arf,
        }
    }

    /// The identity of the odd-char group with the given `kappa`.
    pub fn oddchar_zero(kappa: u8) -> Self {
        WittClassG::OddChar {
            kappa,
            e0: 0,
            sclass: 0,
        }
    }

    /// The group operation `⊥`. Panics if the two classes are in different
    /// characteristic regimes (you cannot add across characteristics).
    pub fn add(&self, other: &WittClassG) -> WittClassG {
        match (*self, *other) {
            (WittClassG::Char0 { signature: a }, WittClassG::Char0 { signature: b }) => {
                WittClassG::Char0 { signature: a + b }
            }
            (WittClassG::Char2 { arf: a }, WittClassG::Char2 { arf: b }) => {
                WittClassG::Char2 { arf: a ^ b }
            }
            (
                WittClassG::OddChar {
                    kappa: ka,
                    e0: e0a,
                    sclass: sa,
                },
                WittClassG::OddChar {
                    kappa: kb,
                    e0: e0b,
                    sclass: sb,
                },
            ) => {
                assert_eq!(ka, kb, "odd-char Witt classes from different fields");
                // signed-disc multiplies with a (−1)^{mn} = (−1)^{e0a·e0b} twist:
                let twist = if e0a & e0b == 1 { ka } else { 0 };
                WittClassG::OddChar {
                    kappa: ka,
                    e0: e0a ^ e0b,
                    sclass: sa ^ sb ^ twist,
                }
            }
            _ => panic!("cannot add Witt classes across characteristics"),
        }
    }

    /// The Witt **ring** multiplication (tensor product of forms, descended to
    /// classes), making `W(F)` a ring with `⟨1⟩` as the unit. Defined for the two
    /// legs where `W` is genuinely a ring:
    ///
    /// * `Char0`: `W(ℝ) ≅ ℤ`, signatures multiply.
    /// * `OddChar`: `W(F_q)` is the order-4 ring — `ℤ/4` when `−1` is a nonsquare
    ///   (`kappa = 1`), via `z = e0 + 2·sclass`; and `F₂[ℤ/2]` when `−1` is a square
    ///   (`kappa = 0`), via `(a,b) = (e0 ⊕ sclass, sclass)` with `t² = 1`. (Both ring
    ///   laws are pinned by `witt_ring`'s test against the concrete `tensor_form`.)
    ///
    /// **Panics on `Char2`:** in characteristic 2 the *quadratic* Witt group `W_q` is
    /// a **module over** the bilinear Witt ring, not a ring — so there is no
    /// quadratic-form ring product to return. (See `forms/witt_ring.rs` for why.)
    pub fn mul(&self, other: &WittClassG) -> WittClassG {
        match (*self, *other) {
            (WittClassG::Char0 { signature: a }, WittClassG::Char0 { signature: b }) => {
                WittClassG::Char0 { signature: a * b }
            }
            (
                WittClassG::OddChar {
                    kappa: ka,
                    e0: e0a,
                    sclass: sa,
                },
                WittClassG::OddChar {
                    kappa: kb,
                    e0: e0b,
                    sclass: sb,
                },
            ) => {
                assert_eq!(ka, kb, "odd-char Witt classes from different fields");
                if ka == 1 {
                    // ℤ/4 via z = e0 + 2·sclass; multiply mod 4.
                    let za = (e0a + 2 * sa) as i32;
                    let zb = (e0b + 2 * sb) as i32;
                    let z = (za * zb).rem_euclid(4);
                    WittClassG::OddChar {
                        kappa: 1,
                        e0: (z & 1) as u8,
                        sclass: ((z >> 1) & 1) as u8,
                    }
                } else {
                    // F₂[ℤ/2] = F₂[t]/(t²−1): (a,b) = (e0⊕sclass, sclass), t² = 1.
                    let (a1, b1) = (e0a ^ sa, sa);
                    let (a2, b2) = (e0b ^ sb, sb);
                    let ar = (a1 & a2) ^ (b1 & b2);
                    let br = (a1 & b2) ^ (a2 & b1);
                    WittClassG::OddChar {
                        kappa: 0,
                        e0: ar ^ br,
                        sclass: br,
                    }
                }
            }
            (WittClassG::Char2 { .. }, WittClassG::Char2 { .. }) => {
                panic!(
                    "char-2 quadratic Witt classes form a MODULE over the bilinear \
                     Witt ring, not a ring — there is no quadratic ring product"
                )
            }
            _ => panic!("cannot multiply Witt classes across characteristics"),
        }
    }

    /// The ring unit `⟨1⟩` of the odd-char Witt ring with the given `kappa`
    /// (`e0 = 1`, `sclass = 0`). The identity for [`mul`](Self::mul).
    pub fn oddchar_one(kappa: u8) -> Self {
        WittClassG::OddChar {
            kappa,
            e0: 1,
            sclass: 0,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeMap;

    fn metric(qs: &[u64], bs: &[((usize, usize), u64)]) -> Metric<Nimber> {
        let q = qs.iter().map(|&x| Nimber(x)).collect();
        let mut b = BTreeMap::new();
        for &((i, j), v) in bs {
            b.insert((i, j), Nimber(v));
        }
        Metric::new(q, b)
    }

    #[test]
    fn hyperbolic_is_identity_anisotropic_is_order_two() {
        let h = WittClass::from_metric(&metric(&[0, 0], &[((0, 1), 1)])); // Arf 0
        let a = WittClass::from_metric(&metric(&[1, 1], &[((0, 1), 1)])); // Arf 1
        assert!(h.is_hyperbolic());
        assert!(!a.is_hyperbolic());
        assert_eq!(h, WittClass::zero());
        assert_eq!(a.anisotropic_dim(), 2);
        // self-inverse: a + a = 0  ⟺  A ⊕ A ≅ H ⊕ H
        assert_eq!(a.add(&a), WittClass::zero());
        assert_eq!(a.add(&h), a); // identity
    }

    #[test]
    fn group_law_is_xor_of_arf() {
        let h = WittClass { arf: 0 };
        let a = WittClass { arf: 1 };
        assert_eq!(a.add(&a), h);
        assert_eq!(a.add(&h), a);
        assert_eq!(h.add(&h), h);
        // direct_sum of the underlying forms agrees with the abstract group law.
        let am = metric(&[1, 1], &[((0, 1), 1)]);
        let combined = WittClass::from_metric(&am.direct_sum(&am));
        assert_eq!(combined, a.add(&a)); // both are 0
    }

    #[test]
    fn witt_class_over_f4() {
        // From arf.rs's F₄ facts: q=[2,2],b=1 is anisotropic (Arf 1); q=[2,3],b=1
        // is hyperbolic-class (Arf 0). Their Witt classes add to the nonzero class.
        let aniso = WittClass::from_metric(&metric(&[2, 2], &[((0, 1), 1)]));
        let split = WittClass::from_metric(&metric(&[2, 3], &[((0, 1), 1)]));
        assert_eq!(aniso.arf, 1);
        assert_eq!(split.arf, 0);
        assert_eq!(aniso.add(&split), aniso);
    }
}
