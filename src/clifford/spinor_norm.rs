//! The spinor norm `N: O(Q) → F*/F*²` and the Dickson grade parity — the two
//! invariants that, together, pin down a versor's place in the Pin/Spin tower.
//!
//! `forms::char2` already classifies isometries in characteristic 2 by the
//! **Dickson invariant** (`SO(Q) = ker D`). The companion across the *other*
//! characteristics is the **spinor norm**: the Pin map sends a versor
//! `v = v₁⋯v_k` (a product of vectors) to the composite of the reflections in the
//! `v_i`; the cokernel over a general field is measured by the spinor norm map
//! `O(Q) → F*/F*²`. Concretely `N(v) = ∏ q(v_i) = ⟨v ṽ⟩₀`, read modulo squares.
//! The pair
//! `(Dickson parity, spinor norm)` is what separates the four cosets
//! `Pin/Spin × ±` of `O(Q)`.
//!
//! ## The characteristic-2 caveat (pinned)
//!
//! In characteristic 2 the codomain is **not** `F*/F*²`. There `x ↦ x²` is the
//! Frobenius (a bijection on a perfect field), so every element is a square and
//! `F*/F*²` is trivial — the multiplicative spinor norm collapses. The correct
//! char-2 spinor norm is **additive**, valued in `F/℘(F)` (the Artin–Schreier
//! group, `℘(x) = x² + x`) — the very group the Arf invariant is pushed into by the
//! field trace (`scalar::nim_trace`). So in char 2 the right companion
//! to Dickson is that additive invariant, not this multiplicative one; we expose
//! the **raw** norm `⟨v ṽ⟩₀` generically (correct as an element of `F`) and leave
//! the "mod squares" / "mod ℘" reduction to the caller's field, where the square /
//! Artin–Schreier test lives.

use crate::clifford::{CliffordAlgebra, Multivector};
use crate::scalar::Scalar;

/// The **Dickson invariant** of a versor: its ℤ₂-grade parity. `0` for an even
/// versor (a rotor, in `SO`), `1` for an odd versor (an odd number of reflections).
/// `None` if the multivector is not of homogeneous grade parity — hence not a
/// versor — or is zero. Generic over the scalar (the char-2 `Nimber` specialisation
/// `forms::dickson_of_versor` delegates here).
pub fn versor_grade_parity<S: Scalar>(v: &Multivector<S>) -> Option<u128> {
    let mut parity: Option<u128> = None;
    for &blade in v.terms.keys() {
        let p = (blade.count_ones() % 2) as u128;
        match parity {
            None => parity = Some(p),
            Some(q) if q != p => return None,
            _ => {}
        }
    }
    parity
}

/// The classification of a versor: its raw spinor norm `⟨v ṽ⟩₀ ∈ F` together with
/// its Dickson grade parity.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VersorInvariants<S: Scalar> {
    /// The raw spinor norm `N(v) = ⟨v ṽ⟩₀ = ∏ q(vᵢ)`. Its class in `F*/F*²`
    /// (char ≠ 2) or `F/℘(F)` (char 2) is the spinor-norm invariant.
    pub spinor_norm: S,
    /// The Dickson invariant (grade parity): `0` in `SO`, `1` an odd reflection.
    pub dickson: u128,
}

impl<S: Scalar> CliffordAlgebra<S> {
    /// The **raw spinor norm** `N(v) = ⟨v ṽ⟩₀` of a versor `v`, returned as a field
    /// element. `Some(N)` iff `v ṽ` is a pure invertible scalar (the same
    /// invertibility gate as [`versor_inverse`](CliffordAlgebra::versor_inverse));
    /// `None` if `v` is not a simple invertible versor. For `v = v₁⋯v_k` this
    /// equals `∏ q(vᵢ)`; reduce it modulo squares (char ≠ 2) or modulo `℘` (char 2)
    /// to get the invariant in the appropriate quotient.
    pub fn spinor_norm(&self, v: &Multivector<S>) -> Option<S> {
        let rev = self.reverse(v);
        let vrev = self.mul(v, &rev);
        let n = self.scalar_part(&vrev);
        if self.scalar(n.clone()) != vrev {
            return None; // v ṽ is not a pure scalar ⇒ not a simple versor
        }
        n.inv()?;
        Some(n)
    }

    /// Classify a versor by `(spinor norm, Dickson parity)`. `None` if `v` is not a
    /// versor (mixed grade parity, or `v ṽ` not scalar).
    pub fn classify_versor(&self, v: &Multivector<S>) -> Option<VersorInvariants<S>> {
        let dickson = versor_grade_parity(v)?;
        let spinor_norm = self.spinor_norm(v)?;
        Some(VersorInvariants {
            spinor_norm,
            dickson,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::clifford::Metric;
    use crate::scalar::{Nimber, Rational};

    fn cl3() -> CliffordAlgebra<Rational> {
        // Cl(3,0) over ℚ: q = [1,1,1], orthonormal.
        CliffordAlgebra::new(
            3,
            Metric::diagonal(vec![Rational::one(), Rational::one(), Rational::one()]),
        )
    }

    #[test]
    fn spinor_norm_of_a_reflection_is_q_of_the_vector() {
        let alg = cl3();
        // a unit reflection vector e0: N = q0 = 1.
        assert_eq!(alg.spinor_norm(&alg.gen(0)), Some(Rational::one()));
        // a non-unit vector v = e0 + e1: N = q0 + q1 = 2 (a nonsquare class in ℚ).
        let v = alg.add(&alg.gen(0), &alg.gen(1));
        assert_eq!(alg.spinor_norm(&v), Some(Rational::int(2)));
    }

    #[test]
    fn spinor_norm_is_multiplicative_on_versors() {
        let alg = cl3();
        let v = alg.add(&alg.gen(0), &alg.gen(1)); // N = 2
        let w = alg.gen(2); // N = 1
        let vw = alg.mul(&v, &w);
        let nv = alg.spinor_norm(&v).unwrap();
        let nw = alg.spinor_norm(&w).unwrap();
        let nvw = alg.spinor_norm(&vw).unwrap();
        assert_eq!(nvw, nv.mul(&nw)); // N(vw) = N(v)·N(w)
    }

    #[test]
    fn dickson_parity_counts_reflections_mod_two() {
        let alg = cl3();
        let scalar_one = alg.scalar(Rational::one());
        let e0 = alg.gen(0);
        let e0e1 = alg.mul(&alg.gen(0), &alg.gen(1));
        let e0e1e2 = alg.mul(&e0e1, &alg.gen(2));
        assert_eq!(versor_grade_parity(&scalar_one), Some(0)); // identity rotor
        assert_eq!(versor_grade_parity(&e0), Some(1)); // 1 reflection
        assert_eq!(versor_grade_parity(&e0e1), Some(0)); // 2 reflections (rotor)
        assert_eq!(versor_grade_parity(&e0e1e2), Some(1)); // 3 reflections
                                                           // mixed grade parity ⇒ not a versor
        let mixed = alg.add(&e0, &e0e1);
        assert_eq!(versor_grade_parity(&mixed), None);
        assert_eq!(alg.classify_versor(&mixed), None);
    }

    #[test]
    fn classify_versor_bundles_both() {
        let alg = cl3();
        let e0e1 = alg.mul(&alg.gen(0), &alg.gen(1));
        let c = alg.classify_versor(&e0e1).unwrap();
        assert_eq!(c.dickson, 0); // a rotor
        assert_eq!(c.spinor_norm, Rational::one()); // q0·q1 = 1
    }

    #[test]
    fn null_homogeneous_elements_are_not_versors() {
        let alg = CliffordAlgebra::<Rational>::new(1, Metric::grassmann(1));
        let e0 = alg.gen(0);
        assert_eq!(versor_grade_parity(&e0), Some(1));
        assert_eq!(alg.spinor_norm(&e0), None);
        assert_eq!(alg.classify_versor(&e0), None);
    }

    #[test]
    fn generic_parity_agrees_with_char2_dickson() {
        // The generic versor_grade_parity reproduces forms::dickson_of_versor on the
        // Nimber backend — the char-2 Dickson is this same grade parity.
        let alg = CliffordAlgebra::new(2, Metric::diagonal(vec![Nimber(1), Nimber(1)]));
        let e0 = alg.gen(0);
        let e0e1 = alg.mul(&alg.gen(0), &alg.gen(1));
        assert_eq!(
            versor_grade_parity(&e0),
            crate::forms::dickson_of_versor(&alg, &e0)
        );
        assert_eq!(
            versor_grade_parity(&e0e1),
            crate::forms::dickson_of_versor(&alg, &e0e1)
        );
    }
}
