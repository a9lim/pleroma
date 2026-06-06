//! The multivector engine, generic over any `Scalar` backend.
//!
//! ## Metric data — characteristic-faithful by design
//!
//! A blade is a `u32` bitmask over basis generators e_0..e_31. The algebra is
//! defined by two independent pieces of data, *not* a single bilinear form:
//!
//!   * `q[i]`      = e_i²                      (the quadratic form / squares)
//!   * `b[(i,j)]`  = e_i e_j + e_j e_i  (i<j)  (the polar / anticommutator form)
//!
//! In characteristic ≠ 2 these are linked (`b = 2·offdiag`, `q = diag`), so an
//! orthogonal basis just sets `b = 0`. In characteristic 2 they are genuinely
//! independent: the polar form is *alternating* (`b(i,i)=0`) yet `q[i]` can be
//! nonzero, and a nonzero off-diagonal `b[(i,j)]` is exactly what makes the
//! nim-Clifford algebra *non-commutative*. Carrying both is the faithful thing.
//!
//! "With nilpotents": set `q[i] = 0` and you get a null generator, e_i² = 0.
//! All `q = 0`, all `b = 0` ⇒ the exterior/Grassmann algebra.
//!
//! ## Product
//!
//! Two canonical blades multiply by concatenating their (ascending) generator
//! lists into a word and reducing to canonical form with the rules
//!   e_i e_i  → q[i]                            (equal adjacent: contract)
//!   e_i e_j  → b[(j,i)] − e_j e_i   (i>j)      (out of order: swap, emit polar)
//! The `−` goes through the scalar's own `neg()`, so in characteristic 2 it is
//! `+` automatically and signs vanish — no special-casing. Termination: each
//! step lowers (word length, inversion count) lexicographically.

use crate::scalar::Scalar;
use std::collections::BTreeMap;

/// Ascending list of set-bit indices of a blade mask.
fn bits(mask: u32) -> Vec<usize> {
    let mut v = Vec::new();
    let mut m = mask;
    while m != 0 {
        let i = m.trailing_zeros() as usize;
        v.push(i);
        m &= m - 1;
    }
    v
}

fn grade(mask: u32) -> u32 {
    mask.count_ones()
}

/// Sign (+1/-1 as a Scalar) of reordering two disjoint ascending blades when
/// concatenated — i.e. the number of (i in a, j in b) with i > j, mod 2.
fn wedge_sign<S: Scalar>(a: u32, b: u32) -> S {
    let mut swaps = 0u32;
    let mut aa = a;
    while aa != 0 {
        let i = aa.trailing_zeros();
        aa &= aa - 1;
        // count bits of b strictly below i
        let below = b & ((1u32 << i) - 1);
        swaps += below.count_ones();
    }
    if swaps & 1 == 0 {
        S::one()
    } else {
        S::one().neg()
    }
}

/// The metric: squares `q` and anticommutators `b` (keyed (i,j) with i<j).
#[derive(Clone, Debug)]
pub struct Metric<S: Scalar> {
    pub q: Vec<S>,
    pub b: BTreeMap<(usize, usize), S>,
}

impl<S: Scalar> Metric<S> {
    /// Orthogonal metric from a list of squares (b = 0). `Cl(p,q,r)` style.
    pub fn diagonal(q: Vec<S>) -> Self {
        Metric { q, b: BTreeMap::new() }
    }

    /// The fully-null metric: exterior/Grassmann algebra on `n` generators.
    pub fn grassmann(n: usize) -> Self {
        Metric { q: vec![S::zero(); n], b: BTreeMap::new() }
    }

    fn q_val(&self, i: usize) -> S {
        self.q.get(i).cloned().unwrap_or_else(S::zero)
    }

    fn b_val(&self, i: usize, j: usize) -> S {
        let key = if i < j { (i, j) } else { (j, i) };
        self.b.get(&key).cloned().unwrap_or_else(S::zero)
    }

    /// Reduce a generator word to canonical multivector terms.
    fn reduce_word(&self, word: &[usize]) -> BTreeMap<u32, S> {
        for p in 0..word.len().saturating_sub(1) {
            let (a, c) = (word[p], word[p + 1]);
            if a == c {
                // e_a e_a = q[a]
                let q = self.q_val(a);
                let mut rest = Vec::with_capacity(word.len() - 2);
                rest.extend_from_slice(&word[..p]);
                rest.extend_from_slice(&word[p + 2..]);
                return scale(self.reduce_word(&rest), &q);
            } else if a > c {
                // e_a e_c = b[(c,a)] - e_c e_a
                let bv = self.b_val(a, c);
                let mut removed = Vec::with_capacity(word.len() - 2);
                removed.extend_from_slice(&word[..p]);
                removed.extend_from_slice(&word[p + 2..]);
                let mut out = scale(self.reduce_word(&removed), &bv);

                let mut swapped = word.to_vec();
                swapped.swap(p, p + 1);
                let neg = S::one().neg();
                merge(&mut out, scale(self.reduce_word(&swapped), &neg));
                return out;
            }
        }
        // strictly increasing & distinct → a single canonical blade
        let mut mask = 0u32;
        for &g in word {
            mask |= 1 << g;
        }
        let mut m = BTreeMap::new();
        m.insert(mask, S::one());
        m
    }
}

fn scale<S: Scalar>(mut terms: BTreeMap<u32, S>, s: &S) -> BTreeMap<u32, S> {
    if s.is_zero() {
        return BTreeMap::new();
    }
    for v in terms.values_mut() {
        *v = v.mul(s);
    }
    terms.retain(|_, v| !v.is_zero());
    terms
}

fn merge<S: Scalar>(into: &mut BTreeMap<u32, S>, other: BTreeMap<u32, S>) {
    for (blade, coeff) in other {
        let e = into.entry(blade).or_insert_with(S::zero);
        *e = e.add(&coeff);
        if e.is_zero() {
            into.remove(&blade);
        }
    }
}

/// A multivector: blade-mask → coefficient (zeros never stored).
#[derive(Clone, Debug, PartialEq)]
pub struct Multivector<S: Scalar> {
    pub terms: BTreeMap<u32, S>,
}

/// A Clifford algebra: dimension + metric. Produces and combines multivectors.
#[derive(Clone, Debug)]
pub struct CliffordAlgebra<S: Scalar> {
    pub dim: usize,
    pub metric: Metric<S>,
}

impl<S: Scalar> CliffordAlgebra<S> {
    pub fn new(dim: usize, metric: Metric<S>) -> Self {
        CliffordAlgebra { dim, metric }
    }

    pub fn zero(&self) -> Multivector<S> {
        Multivector { terms: BTreeMap::new() }
    }

    pub fn scalar(&self, s: S) -> Multivector<S> {
        let mut terms = BTreeMap::new();
        if !s.is_zero() {
            terms.insert(0u32, s);
        }
        Multivector { terms }
    }

    /// The basis vector e_i.
    pub fn gen(&self, i: usize) -> Multivector<S> {
        let mut terms = BTreeMap::new();
        terms.insert(1u32 << i, S::one());
        Multivector { terms }
    }

    /// A single basis blade from a set of generators, coefficient 1.
    pub fn blade(&self, gens: &[usize]) -> Multivector<S> {
        let mut mask = 0u32;
        for &g in gens {
            mask |= 1 << g;
        }
        let mut terms = BTreeMap::new();
        terms.insert(mask, S::one());
        Multivector { terms }
    }

    pub fn add(&self, a: &Multivector<S>, b: &Multivector<S>) -> Multivector<S> {
        let mut terms = a.terms.clone();
        merge(&mut terms, b.terms.clone());
        Multivector { terms }
    }

    pub fn scalar_mul(&self, s: &S, a: &Multivector<S>) -> Multivector<S> {
        Multivector { terms: scale(a.terms.clone(), s) }
    }

    /// Geometric (Clifford) product.
    pub fn mul(&self, a: &Multivector<S>, b: &Multivector<S>) -> Multivector<S> {
        let mut out: BTreeMap<u32, S> = BTreeMap::new();
        for (&ba, ca) in &a.terms {
            for (&bb, cb) in &b.terms {
                let mut word = bits(ba);
                word.extend(bits(bb));
                let reduced = self.metric.reduce_word(&word);
                let coeff = ca.mul(cb);
                merge(&mut out, scale(reduced, &coeff));
            }
        }
        Multivector { terms: out }
    }

    /// Exterior (wedge) product — metric-independent.
    pub fn wedge(&self, a: &Multivector<S>, b: &Multivector<S>) -> Multivector<S> {
        let mut out: BTreeMap<u32, S> = BTreeMap::new();
        for (&ba, ca) in &a.terms {
            for (&bb, cb) in &b.terms {
                if ba & bb != 0 {
                    continue; // shared generator ⇒ wedge is 0
                }
                let sign = wedge_sign::<S>(ba, bb);
                let coeff = ca.mul(cb).mul(&sign);
                if coeff.is_zero() {
                    continue;
                }
                let e = out.entry(ba | bb).or_insert_with(S::zero);
                *e = e.add(&coeff);
                if e.is_zero() {
                    out.remove(&(ba | bb));
                }
            }
        }
        Multivector { terms: out }
    }

    /// Reversion: reverse the order of generators in every blade.
    /// On a grade-k blade this is (-1)^{k(k-1)/2}.
    pub fn reverse(&self, a: &Multivector<S>) -> Multivector<S> {
        let mut terms = BTreeMap::new();
        for (&blade, coeff) in &a.terms {
            let k = grade(blade);
            let c = if (k * (k.wrapping_sub(1)) / 2) & 1 == 1 {
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

    /// Grade-k projection.
    pub fn grade_part(&self, a: &Multivector<S>, k: u32) -> Multivector<S> {
        let terms = a
            .terms
            .iter()
            .filter(|(&blade, _)| grade(blade) == k)
            .map(|(&blade, c)| (blade, c.clone()))
            .collect();
        Multivector { terms }
    }

    /// The grade-0 (scalar) coefficient.
    pub fn scalar_part(&self, v: &Multivector<S>) -> S {
        v.terms.get(&0).cloned().unwrap_or_else(S::zero)
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
            let c = if grade(blade) & 1 == 1 { coeff.neg() } else { coeff.clone() };
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

    /// The sandwich product v x v⁻¹ (rotor / versor action). `None` if v isn't
    /// invertible as a versor.
    pub fn sandwich(&self, v: &Multivector<S>, x: &Multivector<S>) -> Option<Multivector<S>> {
        let vinv = self.versor_inverse(v)?;
        Some(self.mul(&self.mul(v, x), &vinv))
    }

    /// Reflection of x in the hyperplane orthogonal to vector n: −(n x n⁻¹).
    pub fn reflect(&self, n: &Multivector<S>, x: &Multivector<S>) -> Option<Multivector<S>> {
        let ninv = self.versor_inverse(n)?;
        let nxni = self.mul(&self.mul(n, x), &ninv);
        Some(self.scalar_mul(&S::one().neg(), &nxni))
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

impl<S: Scalar> Multivector<S> {
    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    /// Human-readable form, e.g. `3 + 2*e0 + 1*e0e1`.
    pub fn display(&self) -> String {
        if self.terms.is_empty() {
            return "0".to_string();
        }
        let one = S::one();
        let neg_one = S::one().neg();
        let mut parts = Vec::new();
        for (&blade, coeff) in &self.terms {
            if blade == 0 {
                parts.push(format!("{:?}", coeff));
                continue;
            }
            let label: String = bits(blade).iter().map(|i| format!("e{}", i)).collect();
            if *coeff == one {
                parts.push(label);
            } else if *coeff == neg_one {
                parts.push(format!("-{}", label));
            } else {
                parts.push(format!("{:?}*{}", coeff, label));
            }
        }
        parts.join(" + ")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::nimber::Nimber;
    use crate::scalar::Rational;
    use crate::surreal::Surreal;

    fn r(n: i128) -> Rational {
        Rational::int(n)
    }

    #[test]
    fn complex_numbers_cl01() {
        // Cl(0,1): one generator with e0^2 = -1, the complex numbers.
        let alg = CliffordAlgebra::new(1, Metric::diagonal(vec![r(-1)]));
        let e0 = alg.gen(0);
        let sq = alg.mul(&e0, &e0);
        assert_eq!(sq, alg.scalar(r(-1)));
    }

    #[test]
    fn cl20_bivector_squares_to_minus_one() {
        // Cl(2,0): e0^2 = e1^2 = 1; e0e1 anticommutes and (e0e1)^2 = -1.
        let alg = CliffordAlgebra::new(2, Metric::diagonal(vec![r(1), r(1)]));
        let e0 = alg.gen(0);
        let e1 = alg.gen(1);
        // anticommute: e0e1 = -(e1e0)
        let e0e1 = alg.mul(&e0, &e1);
        let e1e0 = alg.mul(&e1, &e0);
        assert_eq!(e0e1, alg.scalar_mul(&r(-1), &e1e0));
        // (e0e1)^2 = -1
        let sq = alg.mul(&e0e1, &e0e1);
        assert_eq!(sq, alg.scalar(r(-1)));
    }

    #[test]
    fn grassmann_generators_are_nilpotent() {
        // q = 0 ⇒ e_i^2 = 0, and the wedge matches the product off-diagonal.
        let alg = CliffordAlgebra::new(3, Metric::grassmann(3));
        for i in 0..3 {
            let ei = alg.gen(i);
            assert!(alg.mul(&ei, &ei).is_zero(), "e{i}^2 should be 0");
        }
        let (e0, e1) = (alg.gen(0), alg.gen(1));
        assert_eq!(alg.mul(&e0, &e1), alg.wedge(&e0, &e1));
        // antisymmetry
        assert_eq!(alg.mul(&e0, &e1), alg.scalar_mul(&r(-1), &alg.mul(&e1, &e0)));
    }

    #[test]
    fn nimber_orthogonal_is_commutative() {
        // char 2, b = 0 ⇒ e_i e_j = e_j e_i (the genuine char-2-orthogonal fact).
        let alg = CliffordAlgebra::new(2, Metric::diagonal(vec![Nimber(2), Nimber(3)]));
        let e0 = alg.gen(0);
        let e1 = alg.gen(1);
        assert_eq!(alg.mul(&e0, &e1), alg.mul(&e1, &e0)); // commute
        // e0^2 = q0 = 2 (a nimber!), not ±1
        assert_eq!(alg.mul(&e0, &e0), alg.scalar(Nimber(2)));
    }

    #[test]
    fn nimber_offdiagonal_is_noncommutative() {
        // char 2 with b[(0,1)] = t ⇒ e0 e1 + e1 e0 = t ≠ 0 ⇒ non-commutative.
        let mut b = BTreeMap::new();
        b.insert((0usize, 1usize), Nimber(1));
        let alg = CliffordAlgebra::new(2, Metric { q: vec![Nimber(0), Nimber(0)], b });
        let e0 = alg.gen(0);
        let e1 = alg.gen(1);
        let anti = alg.add(&alg.mul(&e0, &e1), &alg.mul(&e1, &e0));
        assert_eq!(anti, alg.scalar(Nimber(1))); // {e0,e1} = 1
        assert_ne!(alg.mul(&e0, &e1), alg.mul(&e1, &e0)); // not commutative
    }

    // The real stress test of reduce_word: associativity on a nontrivial,
    // non-orthogonal metric, in both characteristics.
    fn assert_associative<S: Scalar>(alg: &CliffordAlgebra<S>, gens: &[Multivector<S>]) {
        for a in gens {
            for b in gens {
                for c in gens {
                    let l = alg.mul(&alg.mul(a, b), c);
                    let r = alg.mul(a, &alg.mul(b, c));
                    assert_eq!(l, r, "associativity failed");
                }
            }
        }
    }

    #[test]
    fn associativity_rational_nonorthogonal() {
        let mut b = BTreeMap::new();
        b.insert((0usize, 1usize), r(1)); // non-orthogonal
        b.insert((1usize, 2usize), r(-1));
        let alg = CliffordAlgebra::new(3, Metric { q: vec![r(1), r(-1), r(2)], b });
        let gens = [
            alg.gen(0),
            alg.gen(1),
            alg.gen(2),
            alg.mul(&alg.gen(0), &alg.gen(1)),
            alg.add(&alg.gen(0), &alg.scalar(r(3))),
        ];
        assert_associative(&alg, &gens);
    }

    #[test]
    fn vector_inverse() {
        // Cl(3,0): unit vector squares to 1, so v⁻¹ = v.
        let alg = CliffordAlgebra::new(3, Metric::diagonal(vec![r(1), r(1), r(1)]));
        let v = alg.gen(0);
        let vi = alg.versor_inverse(&v).unwrap();
        assert_eq!(alg.mul(&v, &vi), alg.scalar(r(1)));
        assert_eq!(vi, v);
        // q=2: e0⁻¹ = e0/2
        let alg2 = CliffordAlgebra::new(1, Metric::diagonal(vec![r(2)]));
        let e0 = alg2.gen(0);
        assert_eq!(alg2.mul(&e0, &alg2.versor_inverse(&e0).unwrap()), alg2.scalar(r(1)));
        // a null vector has no inverse
        let alg0 = CliffordAlgebra::new(1, Metric::<Rational>::grassmann(1));
        assert!(alg0.versor_inverse(&alg0.gen(0)).is_none());
    }

    #[test]
    fn reflection_fixes_and_negates() {
        // reflect in the hyperplane ⊥ e1: fixes e0, negates e1.
        let alg = CliffordAlgebra::new(2, Metric::diagonal(vec![r(1), r(1)]));
        let (e0, e1) = (alg.gen(0), alg.gen(1));
        assert_eq!(alg.reflect(&e1, &e0).unwrap(), e0);
        assert_eq!(alg.reflect(&e1, &e1).unwrap(), alg.scalar_mul(&r(-1), &e1));
    }

    #[test]
    fn rotor_preserves_norm() {
        // A rotor (product of two unit vectors) is norm-preserving.
        let alg = CliffordAlgebra::new(2, Metric::diagonal(vec![r(1), r(1)]));
        let rotor = alg.mul(&alg.gen(0), &alg.gen(1));
        let x = alg.add(
            &alg.scalar_mul(&r(3), &alg.gen(0)),
            &alg.scalar_mul(&r(4), &alg.gen(1)),
        );
        let rx = alg.sandwich(&rotor, &x).unwrap();
        assert_eq!(alg.norm2(&rx), alg.norm2(&x)); // 25 either way
    }

    #[test]
    fn left_contraction_lowers_grade() {
        let alg = CliffordAlgebra::new(3, Metric::diagonal(vec![r(1), r(1), r(1)]));
        let e0 = alg.gen(0);
        let e0e1 = alg.mul(&alg.gen(0), &alg.gen(1));
        assert_eq!(alg.left_contract(&e0, &e0e1), alg.gen(1)); // e0 ⌟ (e0∧e1) = e1
        let three = alg.scalar(r(3));
        assert_eq!(alg.left_contract(&three, &e0e1), alg.scalar_mul(&r(3), &e0e1));
    }

    #[test]
    fn dual_of_vector_is_bivector_in_3d() {
        let alg = CliffordAlgebra::new(3, Metric::diagonal(vec![r(1), r(1), r(1)]));
        let d = alg.dual(&alg.gen(0)).unwrap();
        assert!(!d.is_zero());
        assert_eq!(alg.grade_part(&d, 2), d); // purely grade 2
    }

    #[test]
    fn grade_involution_signs() {
        let alg = CliffordAlgebra::new(2, Metric::diagonal(vec![r(1), r(1)]));
        // 5 + e0 + e0e1  ↦  5 − e0 + e0e1
        let v = alg.add(
            &alg.scalar(r(5)),
            &alg.add(&alg.gen(0), &alg.mul(&alg.gen(0), &alg.gen(1))),
        );
        let expect = alg.add(
            &alg.scalar(r(5)),
            &alg.add(&alg.scalar_mul(&r(-1), &alg.gen(0)), &alg.mul(&alg.gen(0), &alg.gen(1))),
        );
        assert_eq!(alg.grade_involution(&v), expect);
    }

    #[test]
    fn versor_over_surreal_metric() {
        // e0² = ω (a monomial ⇒ invertible). e0⁻¹ = ε·e0, and sandwich works.
        let alg = CliffordAlgebra::new(2, Metric::diagonal(vec![Surreal::omega(), Surreal::epsilon()]));
        let e0 = alg.gen(0);
        let inv = alg.versor_inverse(&e0).unwrap();
        assert_eq!(alg.mul(&e0, &inv), alg.scalar(Surreal::one()));
    }

    #[test]
    fn associativity_nimber_nonorthogonal() {
        let mut b = BTreeMap::new();
        b.insert((0usize, 1usize), Nimber(1));
        b.insert((0usize, 2usize), Nimber(3));
        let alg = CliffordAlgebra::new(3, Metric { q: vec![Nimber(2), Nimber(1), Nimber(0)], b });
        let gens = [
            alg.gen(0),
            alg.gen(1),
            alg.gen(2),
            alg.mul(&alg.gen(0), &alg.gen(1)),
            alg.add(&alg.gen(2), &alg.scalar(Nimber(5))),
        ];
        assert_associative(&alg, &gens);
    }
}
