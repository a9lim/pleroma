use super::basis::{bits, grade, wedge_sign, MAX_BASIS_DIM};
use super::metric::Metric;
use super::multivector::Multivector;
use super::terms::{merge, scale};
use crate::scalar::Scalar;
use std::collections::BTreeMap;

/// A Clifford algebra: dimension + metric. Produces and combines multivectors.
#[derive(Clone, Debug, PartialEq)]
pub struct CliffordAlgebra<S: Scalar> {
    pub dim: usize,
    pub metric: Metric<S>,
}

impl<S: Scalar> CliffordAlgebra<S> {
    pub fn new(dim: usize, metric: Metric<S>) -> Self {
        metric.validate_for_dim(dim);
        CliffordAlgebra { dim, metric }
    }

    /// The graded (super) tensor product Cl(self) ⊗̂ Cl(other) ≅ Cl(self ⟂ other).
    pub fn graded_tensor(&self, other: &CliffordAlgebra<S>) -> CliffordAlgebra<S> {
        CliffordAlgebra::new(self.dim + other.dim, self.metric.direct_sum(&other.metric))
    }

    /// Embed a multivector of the first factor into `self ⊗̂ other`.
    pub fn embed_first(&self, v: &Multivector<S>) -> Multivector<S> {
        Multivector {
            terms: v.terms.clone(),
        }
    }

    /// Embed a multivector of the second factor into `first ⊗̂ self`.
    pub fn embed_second(&self, v: &Multivector<S>, shift: usize) -> Multivector<S> {
        assert!(shift <= MAX_BASIS_DIM, "basis shift out of range");
        let terms = v
            .terms
            .iter()
            .map(|(&blade, c)| {
                if blade != 0 {
                    let highest = (u128::BITS - 1 - blade.leading_zeros()) as usize;
                    assert!(
                        highest + shift < MAX_BASIS_DIM,
                        "embedded blade exceeds {MAX_BASIS_DIM} generators"
                    );
                }
                let shifted = if blade == 0 { 0 } else { blade << shift };
                (shifted, c.clone())
            })
            .collect();
        Multivector { terms }
    }

    pub fn zero(&self) -> Multivector<S> {
        Multivector {
            terms: BTreeMap::new(),
        }
    }

    pub fn scalar(&self, s: S) -> Multivector<S> {
        let mut terms = BTreeMap::new();
        if !s.is_zero() {
            terms.insert(0u128, s);
        }
        Multivector { terms }
    }

    /// The basis vector e_i.
    pub fn gen(&self, i: usize) -> Multivector<S> {
        assert!(i < self.dim, "generator index {i} out of range");
        assert!(i < MAX_BASIS_DIM, "generator index {i} exceeds blade mask");
        let mut terms = BTreeMap::new();
        terms.insert(1u128 << i, S::one());
        Multivector { terms }
    }

    /// A single basis blade from a set of generators, coefficient 1.
    pub fn blade(&self, gens: &[usize]) -> Multivector<S> {
        let mut mask = 0u128;
        for &g in gens {
            assert!(g < self.dim, "blade generator index {g} out of range");
            assert!(g < MAX_BASIS_DIM, "blade generator index {g} exceeds mask");
            assert!(
                mask & (1u128 << g) == 0,
                "blade expects a set of distinct generators"
            );
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
        Multivector {
            terms: scale(a.terms.clone(), s),
        }
    }

    /// Geometric (Clifford) product.
    pub fn mul(&self, a: &Multivector<S>, b: &Multivector<S>) -> Multivector<S> {
        let mut out: BTreeMap<u128, S> = BTreeMap::new();
        for (&ba, ca) in &a.terms {
            for (&bb, cb) in &b.terms {
                let reduced = self.metric.geom_product_blades(ba, bb);
                let coeff = ca.mul(cb);
                merge(&mut out, scale(reduced, &coeff));
            }
        }
        Multivector { terms: out }
    }

    /// Exterior (wedge) product — metric-independent.
    pub fn wedge(&self, a: &Multivector<S>, b: &Multivector<S>) -> Multivector<S> {
        let mut out: BTreeMap<u128, S> = BTreeMap::new();
        for (&ba, ca) in &a.terms {
            for (&bb, cb) in &b.terms {
                if ba & bb != 0 {
                    continue;
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

    /// Reversion: reverse the order of generators in every blade and reduce that
    /// reversed word through the Clifford product.
    pub fn reverse(&self, a: &Multivector<S>) -> Multivector<S> {
        let mut out = self.zero();
        for (&blade, coeff) in &a.terms {
            let mut rev_blade = self.scalar(S::one());
            let mut gens = bits(blade);
            gens.reverse();
            for g in gens {
                rev_blade = self.mul(&rev_blade, &self.gen(g));
            }
            out = self.add(&out, &self.scalar_mul(coeff, &rev_blade));
        }
        out
    }

    /// Grade-k projection.
    pub fn grade_part(&self, a: &Multivector<S>, k: usize) -> Multivector<S> {
        let terms = a
            .terms
            .iter()
            .filter(|&(&blade, _)| grade(blade) == k)
            .map(|(&blade, c)| (blade, c.clone()))
            .collect();
        Multivector { terms }
    }

    /// The grade-0 (scalar) coefficient.
    pub fn scalar_part(&self, v: &Multivector<S>) -> S {
        v.terms.get(&0).cloned().unwrap_or_else(S::zero)
    }
}
