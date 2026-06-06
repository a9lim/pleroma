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
        let d = self.dim;
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
        let d = self.dim;
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
        let mask = if self.dim >= MAX_BASIS_DIM {
            u128::MAX
        } else {
            (1u128 << self.dim) - 1
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

    /// The undual v ↦ v I — the inverse of [`dual`](Self::dual) (`dual` then
    /// `undual` is the identity). Always defined (no inversion needed).
    pub fn undual(&self, v: &Multivector<S>) -> Multivector<S> {
        self.mul(v, &self.pseudoscalar())
    }

    /// The **Clifford (main) conjugate** `x̄ = α(x̃)` — reversion composed with
    /// grade involution. The third standard involution alongside
    /// [`reverse`](Self::reverse) and [`grade_involution`](Self::grade_involution);
    /// on a grade-`k` blade it is the sign `(−1)^{k(k+1)/2}`.
    pub fn clifford_conjugate(&self, v: &Multivector<S>) -> Multivector<S> {
        self.grade_involution(&self.reverse(v))
    }

    /// The **scalar product** ⟨a b⟩₀ — the grade-0 part of the geometric product.
    pub fn scalar_product(&self, a: &Multivector<S>, b: &Multivector<S>) -> S {
        self.scalar_part(&self.mul(a, b))
    }

    /// The **commutator product** `[a,b] = ab − ba`. No `½` factor, so it is
    /// char-faithful (in characteristic 2 it coincides with the anticommutator,
    /// as it must — there is no sign to separate them).
    pub fn commutator(&self, a: &Multivector<S>, b: &Multivector<S>) -> Multivector<S> {
        let ab = self.mul(a, b);
        let ba = self.mul(b, a);
        self.add(&ab, &self.scalar_mul(&S::one().neg(), &ba))
    }

    /// The **anticommutator product** `{a,b} = ab + ba` (no `½` factor). On two
    /// grade-1 vectors this is the polar form `2B(a,b)` carried by the metric.
    pub fn anticommutator(&self, a: &Multivector<S>, b: &Multivector<S>) -> Multivector<S> {
        self.add(&self.mul(a, b), &self.mul(b, a))
    }

    /// The **regressive (meet) product** `a ∨ b = J⁻¹(J(a) ∧ J(b))`, with `J` the
    /// pseudoscalar dual — the geometric *intersection* of the subspaces `a` and
    /// `b` represent (dual to the wedge, which is their join/span). `None` if the
    /// pseudoscalar is not invertible (a degenerate metric).
    pub fn meet(&self, a: &Multivector<S>, b: &Multivector<S>) -> Option<Multivector<S>> {
        let da = self.dual(a)?;
        let db = self.dual(b)?;
        Some(self.undual(&self.wedge(&da, &db)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scalar::Rational;

    fn r(n: i128) -> Rational {
        Rational::int(n)
    }
    fn euclid(n: usize) -> CliffordAlgebra<Rational> {
        CliffordAlgebra::new(n, Metric::diagonal(vec![r(1); n]))
    }

    #[test]
    fn clifford_conjugate_signs_by_grade() {
        // (−1)^{k(k+1)/2}: scalar +, vector −, bivector −, trivector +.
        let alg = euclid(3);
        let s = alg.scalar(r(1));
        let e0 = alg.gen(0);
        let e01 = alg.wedge(&alg.gen(0), &alg.gen(1));
        let e012 = alg.pseudoscalar();
        assert_eq!(alg.clifford_conjugate(&s), s);
        assert_eq!(alg.clifford_conjugate(&e0), alg.scalar_mul(&r(-1), &e0));
        assert_eq!(alg.clifford_conjugate(&e01), alg.scalar_mul(&r(-1), &e01));
        assert_eq!(alg.clifford_conjugate(&e012), e012);
    }

    #[test]
    fn scalar_and_commutator_products() {
        let alg = euclid(3);
        let (e0, e1) = (alg.gen(0), alg.gen(1));
        // ⟨e0 e0⟩₀ = q0 = 1; ⟨e0 e1⟩₀ = 0 (orthogonal).
        assert_eq!(alg.scalar_product(&e0, &e0), r(1));
        assert_eq!(alg.scalar_product(&e0, &e1), r(0));
        // orthogonal ⇒ {e0,e1} = 0 and [e0,e1] = 2 e0e1.
        assert!(alg.anticommutator(&e0, &e1).is_zero());
        let two_e01 = alg.scalar_mul(&r(2), &alg.wedge(&e0, &e1));
        assert_eq!(alg.commutator(&e0, &e1), two_e01);
    }

    #[test]
    fn meet_of_two_planes_is_their_common_line() {
        // In Cl(3,0): the planes e0∧e1 and e1∧e2 meet in the line e1. The result
        // must be a nonzero grade-1 vector contained in both planes.
        let alg = euclid(3);
        let p1 = alg.wedge(&alg.gen(0), &alg.gen(1));
        let p2 = alg.wedge(&alg.gen(1), &alg.gen(2));
        let line = alg.meet(&p1, &p2).unwrap();
        assert!(!line.is_zero());
        assert_eq!(alg.grade_part(&line, 1), line); // pure grade 1
                                                    // contained in both planes ⇔ wedge with each vanishes
        assert!(alg.wedge(&line, &p1).is_zero());
        assert!(alg.wedge(&line, &p2).is_zero());
    }

    #[test]
    fn dual_undual_round_trip() {
        let alg = euclid(3);
        let v = alg.add(&alg.gen(0), &alg.wedge(&alg.gen(1), &alg.gen(2)));
        let back = alg.undual(&alg.dual(&v).unwrap());
        assert_eq!(back, v);
    }
}
