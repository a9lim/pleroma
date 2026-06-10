//! The Witt group of quadratic forms over a nim-field — the abstraction that
//! sits behind the `A ⊕ A ≅ H ⊕ H` fact the char-2 Arf tests check pointwise.
//!
//! Two nonsingular quadratic forms are **Witt-equivalent** if they become
//! isomorphic after adding hyperbolic planes; the equivalence classes form an
//! abelian group `W_q(F)` under orthogonal sum `⊥`, with the hyperbolic plane as
//! identity. Over a *finite* field of characteristic 2 the anisotropic forms are
//! just two — the zero form (Arf 0) and the unique anisotropic plane (Arf 1) —
//! so `W_q(F_{2^m}) ≅ ℤ/2`, **classified completely by the Arf invariant at that
//! fixed field**, and the group law is XOR of Arf invariants. (Over the full
//! algebraically-closed On₂, `W_q` is trivial; over other fields such as `F₂(t)`,
//! `W_q` can be richer. For each finite field this engine targets, Arf is the
//! whole story.)
//!
//! So `WittClass` makes the additivity executable as a group: `w(A) + w(A) = 0`
//! is the same statement as `A ⊕ A ≅ H ⊕ H`, now a one-liner.

use crate::clifford::Metric;
use crate::forms::arf_invariant;
use crate::scalar::{nim_degree, Nimber};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WittClassError {
    GeneralBilinearMetric,
    Singular {
        radical_dim: usize,
        radical_anisotropic: bool,
    },
}

/// A class in the Witt group `W_q(F) ≅ ℤ/2` of a finite nim-field: the Arf
/// invariant of a form's anisotropic core (hyperbolic planes are the identity).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct WittClass {
    /// Degree `m` of the finite char-2 field `F_{2^m}` this class lives over.
    pub field_degree: u128,
    /// The class, 0 or 1 — equivalently the Arf invariant of the nonsingular core.
    pub arf: u128,
}

impl WittClass {
    /// The identity over `F₂`: the class of the hyperbolic plane (and of the
    /// zero form). Use [`zero_over`](Self::zero_over) when the ground field is
    /// a larger finite char-2 field.
    pub fn zero() -> Self {
        WittClass {
            field_degree: 1,
            arf: 0,
        }
    }

    /// The identity over `F_{2^field_degree}`.
    pub fn zero_over(field_degree: u128) -> Self {
        assert!(field_degree > 0, "char-2 field degree must be positive");
        WittClass {
            field_degree,
            arf: 0,
        }
    }

    /// Checked Witt class of a nimber Clifford metric. The Witt group here is the
    /// group of nonsingular quadratic forms, so a nonzero polar-form radical is
    /// rejected instead of being silently erased.
    pub fn try_from_metric(metric: &Metric<Nimber>) -> Result<Self, WittClassError> {
        let arf = arf_invariant(metric).ok_or(WittClassError::GeneralBilinearMetric)?;
        if arf.radical_dim != 0 {
            return Err(WittClassError::Singular {
                radical_dim: arf.radical_dim,
                radical_anisotropic: arf.radical_anisotropic,
            });
        }
        Ok(WittClass {
            field_degree: nimber_metric_field_degree(metric),
            arf: arf.arf,
        })
    }

    /// The group operation: the class of the orthogonal sum `⊥` of two forms,
    /// checked to stay over the same finite field. Arf is additive only after
    /// the base field is fixed; cross-field sums must first be re-evaluated over
    /// a common extension.
    pub fn try_add(&self, other: &WittClass) -> Result<WittClass, &'static str> {
        if self.field_degree != other.field_degree {
            return Err("char-2 Witt classes are from different finite fields");
        }
        Ok(WittClass {
            field_degree: self.field_degree,
            arf: self.arf ^ other.arf,
        })
    }

    /// Infallible convenience wrapper for callers that already know both
    /// operands live over the same finite field.
    pub fn add(&self, other: &WittClass) -> WittClass {
        self.try_add(other)
            .expect("char-2 Witt classes are from different finite fields")
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
        let field = format!("F_2^{}", self.field_degree);
        if self.arf == 0 {
            format!("0 (hyperbolic class over {field})")
        } else {
            format!("[anisotropic plane] (Arf 1 over {field})")
        }
    }
}

fn nimber_metric_field_degree(metric: &Metric<Nimber>) -> u128 {
    metric
        .q
        .iter()
        .map(|x| nim_degree(x.0))
        .chain(metric.b.values().map(|x| nim_degree(x.0)))
        .max()
        .unwrap_or(1)
}

impl std::ops::Add for WittClass {
    type Output = WittClass;

    fn add(self, rhs: WittClass) -> WittClass {
        WittClass::add(&self, &rhs)
    }
}

impl std::ops::Neg for WittClass {
    type Output = WittClass;

    fn neg(self) -> WittClass {
        WittClass::neg(&self)
    }
}

/// The Witt class across **all three characteristics** — the group-theoretic
/// home of the classifier trichotomy (char-0 signature / odd-char
/// discriminant / char-2 Arf), mirroring the Artin–Schreier↔Arf unification.
///
/// * `Char0`: over the exact-square surreal subdomain, the real-table Witt class
///   is classified by the signature `p − q`; forms outside that subdomain are
///   rejected by the classifier instead of being collapsed to a false real class.
/// * `OddChar`: over a finite field `F_q` of odd characteristic `W(F_q)` has
///   order 4. Its invariants are `e0 = dim mod 2` and `sclass` = the
///   **signed discriminant** `(−1)^{m(m−1)/2}·det` mod squares (a genuine Witt
///   invariant, unlike the ordinary det when `−1` is a nonsquare). The group is
///   `ℤ/4` when `−1` is a nonsquare (`q ≡ 3 mod 4`, `kappa = 1`) and `ℤ/2 × ℤ/2`
///   when `−1` is a square (`q ≡ 1 mod 4`, `kappa = 0`). The `(−1)^{mn}`
///   correction in the signed-discriminant sum is exactly the `kappa` term in
///   `add`, and is what produces the `ℤ/4` when `kappa = 1`.
/// * `Char2`: over a fixed finite char-2 field `W ≅ ℤ/2`, classified by the Arf
///   invariant together with the field degree.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WittClassG {
    Char0 {
        signature: i128,
    },
    OddChar {
        /// Field order `q`; finite fields of the same order are canonically unique.
        field_order: u128,
        /// nonsquareness of `−1`: 0 if `−1` is a square (`q≡1 mod 4`), else 1.
        kappa: u128,
        /// dimension mod 2.
        e0: u128,
        /// signed-discriminant square-class: 0 if a square, 1 if a nonsquare.
        sclass: u128,
    },
    Char2 {
        /// Field degree `m` for `F_{2^m}`.
        field_degree: u128,
        arf: u128,
    },
}

impl WittClassG {
    /// Char-0 Witt class from a signature `(p, q)`.
    pub fn char0(p: usize, q: usize) -> Self {
        WittClassG::Char0 {
            signature: p as i128 - q as i128,
        }
    }

    /// Checked char-2 Witt class from a nonsingular nimber metric.
    pub fn try_char2_from_metric(metric: &Metric<Nimber>) -> Result<Self, WittClassError> {
        let class = WittClass::try_from_metric(metric)?;
        Ok(WittClassG::Char2 {
            field_degree: class.field_degree,
            arf: class.arf,
        })
    }

    /// The identity of the odd-char group with the given `kappa`.
    pub fn oddchar_zero(field_order: u128, kappa: u128) -> Self {
        WittClassG::OddChar {
            field_order,
            kappa,
            e0: 0,
            sclass: 0,
        }
    }

    /// The group operation `⊥`, checked because classes from different
    /// characteristic regimes cannot be added.
    pub fn try_add(&self, other: &WittClassG) -> Result<WittClassG, &'static str> {
        match (*self, *other) {
            (WittClassG::Char0 { signature: a }, WittClassG::Char0 { signature: b }) => {
                Ok(WittClassG::Char0 { signature: a + b })
            }
            (
                WittClassG::Char2 {
                    field_degree: ma,
                    arf: a,
                },
                WittClassG::Char2 {
                    field_degree: mb,
                    arf: b,
                },
            ) => {
                if ma != mb {
                    return Err("char-2 Witt classes are from different finite fields");
                }
                Ok(WittClassG::Char2 {
                    field_degree: ma,
                    arf: a ^ b,
                })
            }
            (
                WittClassG::OddChar {
                    field_order: qa,
                    kappa: ka,
                    e0: e0a,
                    sclass: sa,
                },
                WittClassG::OddChar {
                    field_order: qb,
                    kappa: kb,
                    e0: e0b,
                    sclass: sb,
                },
            ) => {
                if qa != qb || ka != kb {
                    return Err("odd-char Witt classes are from different finite fields");
                }
                // signed-disc multiplies with a (−1)^{mn} = (−1)^{e0a·e0b} twist:
                let twist = if e0a & e0b == 1 { ka } else { 0 };
                Ok(WittClassG::OddChar {
                    field_order: qa,
                    kappa: ka,
                    e0: e0a ^ e0b,
                    sclass: sa ^ sb ^ twist,
                })
            }
            _ => Err("cannot add Witt classes across characteristics"),
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
    ///   laws are pinned by `witt::ring`'s test against the concrete `tensor_form`.)
    ///
    /// In characteristic 2 the *quadratic* Witt group `W_q` is a **module over**
    /// the bilinear Witt ring, not a ring, so char-2 operands are rejected instead
    /// of forcing an infallible product.
    pub fn try_mul(&self, other: &WittClassG) -> Result<WittClassG, &'static str> {
        match (*self, *other) {
            (WittClassG::Char0 { signature: a }, WittClassG::Char0 { signature: b }) => {
                Ok(WittClassG::Char0 { signature: a * b })
            }
            (
                WittClassG::OddChar {
                    field_order: qa,
                    kappa: ka,
                    e0: e0a,
                    sclass: sa,
                },
                WittClassG::OddChar {
                    field_order: qb,
                    kappa: kb,
                    e0: e0b,
                    sclass: sb,
                },
            ) => {
                if qa != qb || ka != kb {
                    return Err("odd-char Witt classes are from different finite fields");
                }
                if ka == 1 {
                    // ℤ/4 via z = e0 + 2·sclass; multiply mod 4.
                    let za = (e0a + 2 * sa) as i128;
                    let zb = (e0b + 2 * sb) as i128;
                    let z = (za * zb).rem_euclid(4);
                    Ok(WittClassG::OddChar {
                        field_order: qa,
                        kappa: 1,
                        e0: (z & 1) as u128,
                        sclass: ((z >> 1) & 1) as u128,
                    })
                } else {
                    // F₂[ℤ/2] = F₂[t]/(t²−1): (a,b) = (e0⊕sclass, sclass), t² = 1.
                    let (a1, b1) = (e0a ^ sa, sa);
                    let (a2, b2) = (e0b ^ sb, sb);
                    let ar = (a1 & a2) ^ (b1 & b2);
                    let br = (a1 & b2) ^ (a2 & b1);
                    Ok(WittClassG::OddChar {
                        field_order: qa,
                        kappa: 0,
                        e0: ar ^ br,
                        sclass: br,
                    })
                }
            }
            (WittClassG::Char2 { .. }, WittClassG::Char2 { .. }) => Err(
                "char-2 quadratic Witt classes form a module over the bilinear Witt ring, not a ring",
            ),
            _ => Err("cannot multiply Witt classes across characteristics"),
        }
    }

    /// The ring unit `⟨1⟩` of the odd-char Witt ring with the given `kappa`
    /// (`e0 = 1`, `sclass = 0`). The identity for [`try_mul`](Self::try_mul).
    pub fn oddchar_one(field_order: u128, kappa: u128) -> Self {
        WittClassG::OddChar {
            field_order,
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

    fn metric(qs: &[u128], bs: &[((usize, usize), u128)]) -> Metric<Nimber> {
        let q = qs.iter().map(|&x| Nimber(x)).collect();
        let mut b = BTreeMap::new();
        for &((i, j), v) in bs {
            b.insert((i, j), Nimber(v));
        }
        Metric::new(q, b)
    }

    #[test]
    fn hyperbolic_is_identity_anisotropic_is_order_two() {
        let h = WittClass::try_from_metric(&metric(&[0, 0], &[((0, 1), 1)]))
            .expect("hyperbolic plane is nonsingular"); // Arf 0
        let a = WittClass::try_from_metric(&metric(&[1, 1], &[((0, 1), 1)]))
            .expect("anisotropic plane is nonsingular"); // Arf 1
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
        let h = WittClass::zero();
        let a = WittClass {
            field_degree: 1,
            arf: 1,
        };
        assert_eq!(a.add(&a), h);
        assert_eq!(a.add(&h), a);
        assert_eq!(h.add(&h), h);
        assert_eq!(a + a, h);
        assert_eq!(-a, a);
        // direct_sum of the underlying forms agrees with the abstract group law.
        let am = metric(&[1, 1], &[((0, 1), 1)]);
        let combined = WittClass::try_from_metric(&am.direct_sum(&am))
            .expect("orthogonal sum of nonsingular planes is nonsingular");
        assert_eq!(combined, a.add(&a)); // both are 0
    }

    #[test]
    fn witt_class_over_f4() {
        // From the char-2 F₄ facts: q=[2,2],b=1 is anisotropic (Arf 1); q=[2,3],b=1
        // is hyperbolic-class (Arf 0). Their Witt classes add to the nonzero class.
        let aniso = WittClass::try_from_metric(&metric(&[2, 2], &[((0, 1), 1)]))
            .expect("F4 anisotropic plane is nonsingular");
        let split = WittClass::try_from_metric(&metric(&[2, 3], &[((0, 1), 1)]))
            .expect("F4 split plane is nonsingular");
        assert_eq!(aniso.field_degree, 2);
        assert_eq!(split.field_degree, 2);
        assert_eq!(aniso.arf, 1);
        assert_eq!(split.arf, 0);
        assert!(split.is_hyperbolic());
        assert_eq!(aniso.add(&split), aniso);
    }

    #[test]
    fn cross_field_char2_witt_addition_is_rejected() {
        let f2_aniso_metric = metric(&[1, 1], &[((0, 1), 1)]);
        let f4_aniso_metric = metric(&[2, 2], &[((0, 1), 1)]);
        let f2_aniso = WittClass::try_from_metric(&f2_aniso_metric)
            .expect("F2 anisotropic plane is nonsingular");
        let f4_aniso = WittClass::try_from_metric(&f4_aniso_metric)
            .expect("F4 anisotropic plane is nonsingular");
        assert_eq!(
            (f2_aniso.field_degree, f2_aniso.arf),
            (1, 1),
            "the F2 plane has Arf 1 over F2"
        );
        assert_eq!(
            (f4_aniso.field_degree, f4_aniso.arf),
            (2, 1),
            "the F4 plane has Arf 1 over F4"
        );

        assert!(f2_aniso.try_add(&f4_aniso).is_err());
        assert!(WittClassG::try_char2_from_metric(&f2_aniso_metric)
            .unwrap()
            .try_add(&WittClassG::try_char2_from_metric(&f4_aniso_metric).unwrap())
            .is_err());

        let summed = f2_aniso_metric.direct_sum(&f4_aniso_metric);
        let re_evaluated = arf_invariant(&summed).expect("direct sum is a nimber metric");
        assert_eq!(
            re_evaluated.arf, 1,
            "bare XOR would predict 0, but the sum re-evaluated over F4 has Arf 1"
        );
    }

    #[test]
    fn singular_forms_are_not_silently_projected_to_witt_classes() {
        let defective = metric(&[0, 0, 1], &[((0, 1), 1)]);
        assert_eq!(
            WittClass::try_from_metric(&defective),
            Err(WittClassError::Singular {
                radical_dim: 1,
                radical_anisotropic: true,
            })
        );
        assert!(WittClassG::try_char2_from_metric(&defective).is_err());
    }
}
