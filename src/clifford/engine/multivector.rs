use super::basis::bits;
use super::terms::{merge, wedge_terms};
use crate::scalar::Scalar;
use std::collections::BTreeMap;
use std::fmt;
use std::ops::{Add, BitXor, Neg, Sub};

/// A multivector: blade-mask → coefficient (zeros never stored).
///
/// ## Operator vs context-method policy
///
/// `Multivector` implements `+`, `-`, unary `-`, and `^` (wedge) as
/// *context-free* operators — no metric needed, so they live on the type.
/// The geometric product (`*`) and all metric-dependent operations
/// (`mul`, `wedge`, `reverse`, `grade_part`, …) live as methods on
/// [`super::algebra::CliffordAlgebra`] and require an algebra context:
///
/// ```text
/// // correct: metric-free additive ops use operators
/// let sum = a + b;
/// let w   = a ^ b;   // exterior/wedge product (metric-independent)
///
/// // correct: metric-dependent ops use the algebra context
/// let prod = alg.mul(&a, &b);
/// let rev  = alg.reverse(&a);
/// ```
///
/// This mirrors the scalar policy from `impl_scalar_ops!`: operators on the
/// concrete type are for the operations that need no extra context; everything
/// else goes through the context object (the algebra) to make the dependency
/// on the metric explicit.
#[derive(Clone, Debug, PartialEq)]
pub struct Multivector<S: Scalar> {
    pub(crate) terms: BTreeMap<u128, S>,
}

impl<S: Scalar> Multivector<S> {
    /// Read-only access to the term map (blade mask → coefficient).
    /// Zeros are never stored; the map is empty iff the multivector is zero.
    pub fn terms(&self) -> &BTreeMap<u128, S> {
        &self.terms
    }

    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    /// Human-readable form, e.g. `3 + 2*e0 + 1*e0e1` (uses `Debug` rendering for
    /// coefficients; works for any `S: Scalar`). The Python binding calls this.
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

/// `fmt::Display` for `Multivector<S>` when `S: fmt::Display` — uses `{}`
/// (Display) for coefficients rather than `{:?}`. Scalars that implement
/// `Display` (e.g. `Fp`, `Fpn`, `Rational` if it did) get clean output.
/// The Python `__repr__` and `display()` method both call the Display-independent
/// path above; this impl is for Rust code that explicitly formats with `{}`.
impl<S: Scalar + fmt::Display> fmt::Display for Multivector<S> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.terms.is_empty() {
            return write!(f, "0");
        }
        let one = S::one();
        let neg_one = S::one().neg();
        let mut parts = Vec::new();
        for (&blade, coeff) in &self.terms {
            if blade == 0 {
                parts.push(format!("{}", coeff));
                continue;
            }
            let label: String = bits(blade).iter().map(|i| format!("e{}", i)).collect();
            if *coeff == one {
                parts.push(label);
            } else if *coeff == neg_one {
                parts.push(format!("-{}", label));
            } else {
                parts.push(format!("{}*{}", coeff, label));
            }
        }
        write!(f, "{}", parts.join(" + "))
    }
}

impl<S: Scalar> Add for Multivector<S> {
    type Output = Multivector<S>;

    fn add(self, rhs: Multivector<S>) -> Multivector<S> {
        let mut terms = self.terms;
        merge(&mut terms, rhs.terms);
        Multivector { terms }
    }
}

impl<S: Scalar> Neg for Multivector<S> {
    type Output = Multivector<S>;

    fn neg(self) -> Multivector<S> {
        let terms = self
            .terms
            .into_iter()
            .map(|(blade, coeff)| (blade, coeff.neg()))
            .filter(|(_, coeff)| !coeff.is_zero())
            .collect();
        Multivector { terms }
    }
}

impl<S: Scalar> Sub for Multivector<S> {
    type Output = Multivector<S>;

    fn sub(self, mut rhs: Multivector<S>) -> Multivector<S> {
        for coeff in rhs.terms.values_mut() {
            *coeff = coeff.neg();
        }
        let mut terms = self.terms;
        merge(&mut terms, rhs.terms);
        Multivector { terms }
    }
}

impl<S: Scalar> BitXor for Multivector<S> {
    type Output = Multivector<S>;

    fn bitxor(self, rhs: Multivector<S>) -> Multivector<S> {
        Multivector {
            terms: wedge_terms(&self.terms, &rhs.terms),
        }
    }
}
