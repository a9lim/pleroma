//! The geometric-algebra layer built on the [`super::engine`] core: versors
//! and the sandwich / twisted-sandwich (Pin) action, reflections, the
//! left/right contractions, the pseudoscalar Hodge dual, grade involution,
//! the spinor norm, and the even-grade subalgebra projection.

use super::engine::*;
use crate::scalar::Scalar;
use std::collections::BTreeMap;

impl<S: Scalar> CliffordAlgebra<S> {
    /// Projection onto the even subalgebra (the sum of even-grade blades). The
    /// even part is closed under the geometric product — it is a subalgebra.
    pub fn even_part(&self, v: &Multivector<S>) -> Multivector<S> {
        let terms = v
            .terms
            .iter()
            .filter(|(&blade, _)| grade(blade) & 1 == 0)
            .map(|(&blade, c)| (blade, c.clone()))
            .collect();
        Multivector { terms }
    }

    /// The even subalgebra Cl⁰ presented as a Clifford algebra one dimension
    /// smaller. For a diagonal (orthogonal) metric with a non-null generator
    /// e_p, the map `f_i = e_i e_p` (i ≠ p) is an algebra isomorphism
    /// Cl(Q)⁰ ≅ Cl(Q′) with `f_i² = −q_i q_p` — the classical
    /// `Cl(p,q)⁰ ≅ Cl(p, q−1) ≅ Cl(q, p−1)`. Returns the smaller algebra, or
    /// `None` if the metric is non-orthogonal (`b`/`a` nonempty) or has no
    /// non-null generator to pivot on.
    pub fn even_subalgebra(&self) -> Option<CliffordAlgebra<S>> {
        if !self.metric.b.is_empty() || self.metric.has_upper() {
            return None; // only the orthogonal case has this clean presentation
        }
        let p = (0..self.dim)
            .rev()
            .find(|&i| !self.metric.q_val(i).is_zero())?;
        let qp = self.metric.q_val(p);
        let qprime: Vec<S> = (0..self.dim)
            .filter(|&i| i != p)
            .map(|i| self.metric.q_val(i).mul(&qp).neg())
            .collect();
        Some(CliffordAlgebra::new(self.dim - 1, Metric::diagonal(qprime)))
    }

    /// The spinor norm ⟨v ṽ⟩₀ (scalar part of `v * reverse(v)`).
    pub fn norm2(&self, v: &Multivector<S>) -> S {
        let rev = self.reverse(v);
        self.scalar_part(&self.mul(v, &rev))
    }

    /// Grade involution: negate every odd-grade blade.
    pub fn grade_involution(&self, v: &Multivector<S>) -> Multivector<S> {
        let mut terms = BTreeMap::new();
        for (&blade, coeff) in &v.terms {
            let c = if grade(blade) & 1 == 1 {
                coeff.neg()
            } else {
                coeff.clone()
            };
            if !c.is_zero() {
                terms.insert(blade, c);
            }
        }
        Multivector { terms }
    }

    /// Inverse of a versor (a product of invertible vectors): v⁻¹ = ṽ / (v ṽ),
    /// valid exactly when `v * reverse(v)` is a nonzero invertible scalar.
    /// Returns `None` otherwise (null vector, non-versor, or scalar norm not
    /// invertible in the backend).
    pub fn versor_inverse(&self, v: &Multivector<S>) -> Option<Multivector<S>> {
        let rev = self.reverse(v);
        let vrev = self.mul(v, &rev);
        let n = self.scalar_part(&vrev);
        if self.scalar(n.clone()) != vrev {
            return None; // v ṽ is not a pure scalar ⇒ not a simple versor
        }
        let ninv = n.inv()?;
        Some(self.scalar_mul(&ninv, &rev))
    }

    /// The (untwisted) sandwich product v x v⁻¹ — the rotor action. Correct for
    /// *even* versors (rotors); for odd versors use `twisted_sandwich`. `None`
    /// if v isn't invertible as a versor.
    pub fn sandwich(&self, v: &Multivector<S>, x: &Multivector<S>) -> Option<Multivector<S>> {
        let vinv = self.versor_inverse(v)?;
        Some(self.mul(&self.mul(v, x), &vinv))
    }

    /// The **twisted adjoint** (Pin/Spin) action: α(v) x v⁻¹, where α is the
    /// grade involution. This is the representation-theoretically correct versor
    /// action on Pin(Q): for an *odd* versor (e.g. a single vector) the α-sign is
    /// exactly what turns it into a genuine reflection in every signature; for an
    /// *even* versor (rotor) α(v)=v, so it coincides with `sandwich`. `None` if v
    /// isn't an invertible versor.
    pub fn twisted_sandwich(
        &self,
        v: &Multivector<S>,
        x: &Multivector<S>,
    ) -> Option<Multivector<S>> {
        let vinv = self.versor_inverse(v)?;
        let av = self.grade_involution(v);
        Some(self.mul(&self.mul(&av, x), &vinv))
    }

    /// Reflection of x in the hyperplane orthogonal to vector n. This is just the
    /// twisted adjoint by the vector n: α(n) x n⁻¹ = −(n x n⁻¹), since n is odd.
    /// Routing through `twisted_sandwich` keeps the single sign convention.
    pub fn reflect(&self, n: &Multivector<S>, x: &Multivector<S>) -> Option<Multivector<S>> {
        self.twisted_sandwich(n, x)
    }

    /// Left contraction a ⌟ b = Σ_{r≤s} ⟨⟨a⟩_r ⟨b⟩_s⟩_{s−r}.
    pub fn left_contract(&self, a: &Multivector<S>, b: &Multivector<S>) -> Multivector<S> {
        let mut out = self.zero();
        let d = self.dim as u32;
        for r in 0..=d {
            let ar = self.grade_part(a, r);
            if ar.is_zero() {
                continue;
            }
            for s in r..=d {
                let bs = self.grade_part(b, s);
                if bs.is_zero() {
                    continue;
                }
                let prod = self.mul(&ar, &bs);
                out = self.add(&out, &self.grade_part(&prod, s - r));
            }
        }
        out
    }

    /// Right contraction a ⌞ b = Σ_{r≥s} ⟨⟨a⟩_r ⟨b⟩_s⟩_{r−s}.
    pub fn right_contract(&self, a: &Multivector<S>, b: &Multivector<S>) -> Multivector<S> {
        let mut out = self.zero();
        let d = self.dim as u32;
        for s in 0..=d {
            let bs = self.grade_part(b, s);
            if bs.is_zero() {
                continue;
            }
            for r in s..=d {
                let ar = self.grade_part(a, r);
                if ar.is_zero() {
                    continue;
                }
                let prod = self.mul(&ar, &bs);
                out = self.add(&out, &self.grade_part(&prod, r - s));
            }
        }
        out
    }

    /// The unit pseudoscalar I = e₀e₁…e_{dim−1}.
    pub fn pseudoscalar(&self) -> Multivector<S> {
        let mask = if self.dim >= 32 {
            u32::MAX
        } else {
            (1u32 << self.dim) - 1
        };
        let mut terms = BTreeMap::new();
        terms.insert(mask, S::one());
        Multivector { terms }
    }

    /// Hodge-style dual v ↦ v I⁻¹. `None` if the pseudoscalar isn't invertible
    /// (a degenerate metric).
    pub fn dual(&self, v: &Multivector<S>) -> Option<Multivector<S>> {
        let i_inv = self.versor_inverse(&self.pseudoscalar())?;
        Some(self.mul(v, &i_inv))
    }
}
