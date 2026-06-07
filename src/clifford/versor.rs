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

    /// The **Cayley transform** of a bivector `B`: the rotor `(1−B)(1+B)⁻¹`. This
    /// is the missing *exact* arrow between the Lie algebra (bivectors, with the
    /// [`commutator`](Self::commutator)) and the Spin group (rotors): it is
    /// **rational** — no `cos`/`sin`, no `½` — so it works over any char-0 field
    /// backend (Rational/Surreal/Surcomplex), where [`cga::exp_nilpotent`] needs
    /// a terminating series and a general `exp` needs transcendentals. The result
    /// is an even, unit-spinor-norm versor. `None` if `1+B` is not invertible.
    ///
    /// The transform is an **involution** (`cayley∘cayley = id`) for `2`
    /// invertible, so [`cayley_inverse`](Self::cayley_inverse) maps a rotor back
    /// to its bivector generator by the same formula. Degenerate in char 2
    /// (`1−B = 1+B`), where it is identically `1` — intended for char ≠ 2.
    pub fn cayley(&self, b: &Multivector<S>) -> Option<Multivector<S>> {
        let one = self.scalar(S::one());
        let neg_b = self.scalar_mul(&S::one().neg(), b);
        let one_minus_b = self.add(&one, &neg_b);
        let one_plus_b = self.add(&one, b);
        let inv = self.multivector_inverse(&one_plus_b)?;
        Some(self.mul(&one_minus_b, &inv))
    }

    /// The inverse Cayley transform — a rotor `R` back to its bivector generator
    /// `(1−R)(1+R)⁻¹`. Identical formula to [`cayley`](Self::cayley) (the map is
    /// an involution); named for intent. `None` if `1+R` is not invertible.
    pub fn cayley_inverse(&self, r: &Multivector<S>) -> Option<Multivector<S>> {
        self.cayley(r)
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

    #[test]
    fn general_multivector_inverse() {
        let alg = euclid(3);
        // A vector: the general inverse matches versor_inverse and v·v⁻¹ = 1.
        let v = alg.add(&alg.gen(0), &alg.scalar_mul(&r(2), &alg.gen(1)));
        let inv = alg.multivector_inverse(&v).unwrap();
        assert_eq!(inv, alg.versor_inverse(&v).unwrap());
        assert_eq!(alg.mul(&v, &inv), alg.scalar(r(1)));
        // 1 + e0 + e1 : NOT a simple versor (v·ṽ = 3 + 2e0 + 2e1 is not scalar),
        // so versor_inverse declines — but the general inverse succeeds two-sided.
        let x = alg.add(&alg.add(&alg.scalar(r(1)), &alg.gen(0)), &alg.gen(1));
        assert!(alg.versor_inverse(&x).is_none());
        let xi = alg.multivector_inverse(&x).unwrap();
        assert_eq!(alg.mul(&x, &xi), alg.scalar(r(1)));
        assert_eq!(alg.mul(&xi, &x), alg.scalar(r(1)));
        // Zero never inverts.
        assert!(alg.multivector_inverse(&alg.zero()).is_none());
    }

    #[test]
    fn multivector_inverse_in_char_two() {
        use crate::scalar::Nimber;
        // Over a nimber field (char 2, neg = id) the Gauss–Jordan pivots exercise
        // Scalar::sub = add. NB `1+e0` is nilpotent here ((1+e0)² = 1+1 = 0), so it
        // is correctly NON-invertible — use a genuine unit `1 + e0 + e1` (it has
        // odd augmentation, so it inverts in this commutative char-2 algebra).
        let alg = CliffordAlgebra::new(2, Metric::diagonal(vec![Nimber(1), Nimber(1)]));
        assert!(alg
            .multivector_inverse(&alg.add(&alg.scalar(Nimber(1)), &alg.gen(0)))
            .is_none()); // 1 + e0 is nilpotent ⇒ no inverse
        let x = alg.add(&alg.add(&alg.scalar(Nimber(1)), &alg.gen(0)), &alg.gen(1));
        let xi = alg.multivector_inverse(&x).unwrap();
        assert_eq!(alg.mul(&x, &xi), alg.scalar(Nimber(1)));
        assert_eq!(alg.mul(&xi, &x), alg.scalar(Nimber(1)));
    }

    #[test]
    fn cayley_bivector_to_rotor() {
        let alg = euclid(3);
        let b = alg.wedge(&alg.gen(0), &alg.gen(1)); // a bivector generator
        let rotor = alg.cayley(&b).unwrap();
        // The rotor is even and unit spinor norm (R ~R = 1).
        assert_eq!(alg.even_part(&rotor), rotor);
        assert_eq!(alg.norm2(&rotor), r(1));
        // Involution: cayley back to the bivector.
        assert_eq!(alg.cayley_inverse(&rotor).unwrap(), b);
        // The rotor's sandwich preserves length.
        let x = alg.gen(0);
        let rx = alg.sandwich(&rotor, &x).unwrap();
        assert_eq!(alg.norm2(&rx), alg.norm2(&x));
    }
}
