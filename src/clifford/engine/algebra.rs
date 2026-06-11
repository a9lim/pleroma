use super::basis::{bits, grade, MAX_BASIS_DIM};
use super::metric::Metric;
use super::multivector::Multivector;
use super::terms::{merge, scale, wedge_terms};
use crate::scalar::Scalar;
use std::collections::BTreeMap;

/// A Clifford algebra: metric + derived dimension. Produces and combines multivectors.
///
/// ## Operator vs context-method policy
///
/// Metric-free additive operations (`+`, `-`, unary `-`, `&` for wedge) are
/// implemented directly on [`Multivector`] as operators. The geometric product
/// and all metric-dependent operations are methods on this type, which provides
/// the metric as context. Use `a + b` / `a & b` for the metric-free ops;
/// `alg.mul(&a, &b)` / `alg.wedge(&a, &b)` for the metric-dependent ones.
/// This mirrors the scalar policy: operators on the concrete type require no
/// extra context; everything that needs context goes through the algebra.
///
/// **Note:** `^` is reserved for scalar power (`x ^ k: u128`); `&` is wedge
/// (`∧` in ogham). See [`Multivector`]'s `BitAnd` impl for the precedence
/// caveat (Rust `&` binds looser than `+` and `*`).
#[derive(Clone, Debug, PartialEq)]
pub struct CliffordAlgebra<S: Scalar> {
    pub(crate) metric: Metric<S>,
}

impl<S: Scalar> CliffordAlgebra<S> {
    pub fn new(dim: usize, metric: Metric<S>) -> Self {
        metric.validate_for_dim(dim);
        CliffordAlgebra { metric }
    }

    /// The number of generators, i.e. the dimension of the underlying vector space.
    /// Derived from the metric (not stored separately); always equal to `metric.dim()`.
    pub fn dim(&self) -> usize {
        self.metric.dim()
    }

    /// Read-only access to the metric of this algebra.
    pub fn metric(&self) -> &Metric<S> {
        &self.metric
    }

    /// The graded (super) tensor product Cl(self) ⊗̂ Cl(other) ≅ Cl(self ⟂ other).
    pub fn graded_tensor(&self, other: &CliffordAlgebra<S>) -> CliffordAlgebra<S> {
        CliffordAlgebra::new(
            self.dim() + other.dim(),
            self.metric.direct_sum(&other.metric),
        )
    }

    /// Embed a multivector of the first factor into `self ⊗̂ other`.
    ///
    /// `self` is ignored — this is a clone of the term map since first-factor
    /// blade masks need no shift. It exists as a method on the algebra for API
    /// symmetry with [`embed_second`](Self::embed_second).
    pub fn embed_first(&self, v: &Multivector<S>) -> Multivector<S> {
        Multivector {
            terms: v.terms.clone(),
        }
    }

    /// Embed a multivector of the second (right) graded-tensor factor into
    /// `left ⊗̂ self` by shifting blade masks left by `left.dim()`.
    ///
    /// The caller passes the left algebra so the shift is read from it directly:
    /// `product_alg.embed_second(&right_mv, &left_alg)`.
    pub fn embed_second(&self, v: &Multivector<S>, left: &CliffordAlgebra<S>) -> Multivector<S> {
        let shift = left.dim();
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

    /// The basis vector `e_i` — named for the `e0∧e1` display language (and to
    /// stay clear of the `gen` keyword reserved in Rust 2024). Python keeps
    /// exposing this as `gen(i)`.
    pub fn e(&self, i: usize) -> Multivector<S> {
        assert!(i < self.dim(), "generator index {i} out of range");
        assert!(i < MAX_BASIS_DIM, "generator index {i} exceeds blade mask");
        let mut terms = BTreeMap::new();
        terms.insert(1u128 << i, S::one());
        Multivector { terms }
    }

    /// A single basis blade from a set of generators, coefficient 1.
    pub fn blade(&self, gens: &[usize]) -> Multivector<S> {
        let mut mask = 0u128;
        for &g in gens {
            assert!(g < self.dim(), "blade generator index {g} out of range");
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
        Multivector {
            terms: wedge_terms(&a.terms, &b.terms),
        }
    }

    /// Reversion: the anti-automorphism `ẽᵢ₁⋯ẽᵢₖ = eᵢₖ⋯eᵢ₁`. Implemented by
    /// reversing the generator list of each wedge-basis blade and reducing through
    /// the Clifford product.
    ///
    /// # Panics
    ///
    /// Panics if the metric has a non-zero `a` (in-order / general-bilinear)
    /// component. In a general-bilinear metric the algebra relations are
    /// asymmetric (`e_i e_j ≠ e_j e_i + symmetric-part`), so blade-by-blade
    /// word reversal is **not** an anti-automorphism of the algebra —
    /// `reverse(xy) ≠ reverse(y)*reverse(x)` in general. Use a symmetric
    /// metric (`Metric::new`/`::diagonal`/`::grassmann`) for operations
    /// that depend on reversion.
    pub fn reverse(&self, a: &Multivector<S>) -> Multivector<S> {
        assert!(
            !self.metric.has_upper(),
            "reverse() is not an anti-automorphism on general-bilinear (a≠0) metrics; \
             use a symmetric metric (Metric::new/diagonal/grassmann)"
        );
        let mut out = self.zero();
        for (&blade, coeff) in &a.terms {
            let mut rev_blade = self.scalar(S::one());
            let mut gens = bits(blade);
            gens.reverse();
            for g in gens {
                rev_blade = self.mul(&rev_blade, &self.e(g));
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

    /// Raise a multivector to a non-negative integer power using square-and-multiply.
    ///
    /// `pow(v, 0)` returns the scalar `1` (the algebra's multiplicative identity),
    /// `pow(v, 1)` returns `v.clone()`, and higher powers are computed via repeated
    /// geometric product — i.e. `self.mul`.
    ///
    /// **Why no `^` operator on `Multivector`?** The geometric product needs the
    /// metric (stored here on the algebra), so iterated geometric multiplication is
    /// not metric-free and cannot live as a bare operator on the `Multivector` type.
    /// Scalar power (`x ^ k: u128` via `impl BitXor<u128>`) is total-product only,
    /// so it CAN live on the scalar type without a metric context. Ogham's `a ↑ k`
    /// desugars to this method for multivectors.
    ///
    /// **Precedence caveat (§5 `spec/ogham.md`):** Rust's `^` binds looser than `*`.
    /// When using scalar `x ^ k`, parenthesize if the intended precedence differs
    /// from ogham's power-tighter-than-product table.
    pub fn pow(&self, v: &Multivector<S>, k: u128) -> Multivector<S> {
        if k == 0 {
            return self.scalar(S::one());
        }
        let mut acc = self.scalar(S::one());
        let mut base = v.clone();
        let mut exp = k;
        // square-and-multiply (binary exponentiation)
        loop {
            if exp & 1 == 1 {
                acc = self.mul(&acc, &base);
            }
            exp >>= 1;
            if exp == 0 {
                break;
            }
            base = self.mul(&base, &base);
        }
        acc
    }
}
