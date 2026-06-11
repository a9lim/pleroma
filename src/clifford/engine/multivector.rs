use super::basis::bits;
use super::terms::{merge, wedge_terms};
use crate::scalar::Scalar;
use std::collections::BTreeMap;
use std::fmt;
use std::ops::{Add, BitAnd, Neg, Sub};

/// A multivector: blade-mask → coefficient (zeros never stored).
///
/// ## Operator vs context-method policy
///
/// `Multivector` implements `+`, `-`, unary `-`, and `&` (wedge / exterior
/// product) as *context-free* operators — no metric needed, so they live on
/// the type. The geometric product (`*`) and all metric-dependent operations
/// (`mul`, `wedge`, `reverse`, `grade_part`, …) live as methods on
/// [`super::algebra::CliffordAlgebra`] and require an algebra context:
///
/// ```text
/// // correct: metric-free additive ops use operators
/// let sum = a + b;
/// let w   = a & b;   // exterior/wedge product (metric-independent)
///                    // ogham ∧; `^` is reserved for power
///
/// // correct: metric-dependent ops use the algebra context
/// let prod = alg.mul(&a, &b);
/// let rev  = alg.reverse(&a);
/// ```
///
/// **Why `&` and not `^` for wedge?** In ogham, `∧`/`&` is the wedge and
/// `↑`/`^` is power. On a type like `Nimber`, element-element `^` would read
/// as XOR = nim-*addition* — not wedge. Using `&` for wedge and reserving `^`
/// for power on scalars (via `impl BitXor<u128>` with a `u128` RHS) makes the
/// type system enforce the distinction: `x ^ y` never compiles when both sides
/// are scalars of the same type (no `BitXor<Self>` impl), preventing the
/// Nimber XOR confusion.
///
/// **Precedence caveat (§5 `spec/ogham.md`):** Rust's `&` binds looser than
/// `+` (and looser than `*`), unlike ogham's wedge-tighter-than-product table.
/// Host code that mixes `+`/`*` and `&` must parenthesize explicitly.
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

    /// Human-readable form, e.g. `3 + 2⋅e0 + e0∧e1` (canonical ogham, Display
    /// v2 §9). A thin alias for the [`fmt::Display`] impl (kept because the
    /// Python binding calls it).
    pub fn display(&self) -> String {
        self.to_string()
    }
}

/// `fmt::Display` for any `Multivector<S>` — `Display` is part of the `Scalar`
/// contract, so coefficients render in their canonical human form (`*n`
/// nimbers, CNF surreals, …). `Debug` on every scalar delegates here-compatible
/// output, so `{}` and `{:?}` agree crate-wide.
impl<S: Scalar> fmt::Display for Multivector<S> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // Display v2 (§9) zero rule: the empty multivector renders as the
        // scalar zero's own display (`*0` in nim-worlds, `0` elsewhere) — bare
        // `0` would not round-trip where bare integers are `E_BareInt`.
        if self.terms.is_empty() {
            return write!(f, "{}", S::zero());
        }
        let one = S::one();
        let neg_one = S::one().neg();
        let mut parts = Vec::new();
        for (&blade, coeff) in &self.terms {
            if blade == 0 {
                parts.push(format!("{coeff}"));
                continue;
            }
            // Display v2 (§9): blades are wedge expressions `e0∧e1` (a single
            // basis vector stays `e0`).
            let label: String = bits(blade)
                .iter()
                .map(|i| format!("e{i}"))
                .collect::<Vec<_>>()
                .join("∧");
            if *coeff == one {
                parts.push(label); // coefficient 1 elided
            } else if *coeff == neg_one {
                parts.push(format!("-{label}")); // -1 → `-label` (via S::one().neg())
            } else {
                // `coeff⋅label`, coefficient parenthesized only when non-atomic.
                parts.push(crate::scalar::poly::attach_coeff(coeff, &label));
            }
        }
        // Display v2 (§9) join rule: a term whose rendering starts with `-`
        // joins with ` - ` (the `-` stripped), string-level and char-agnostic
        // (no sign predicate on `Scalar` exists or is wanted).
        let mut out = String::new();
        for (idx, part) in parts.iter().enumerate() {
            if let Some(stripped) = part.strip_prefix('-') {
                if idx == 0 {
                    out.push('-');
                    out.push_str(stripped);
                } else {
                    out.push_str(" - ");
                    out.push_str(stripped);
                }
            } else {
                if idx != 0 {
                    out.push_str(" + ");
                }
                out.push_str(part);
            }
        }
        write!(f, "{out}")
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

impl<S: Scalar> BitAnd for Multivector<S> {
    type Output = Multivector<S>;

    /// Exterior (wedge) product `a & b` — ogham `a ∧ b`.
    ///
    /// This is metric-independent: it computes the exterior product of the
    /// two term maps directly. In ogham and here, wedge is `∧`/`&`; `^` is
    /// reserved for power. On `Nimber`, an element-element `^` would read as
    /// XOR = nim-*addition*, which is why `BitXor<Self>` does not exist on any
    /// backend: the type system enforces the disambiguation.
    ///
    /// **Precedence caveat (§5 `spec/ogham.md`):** Rust's `&` binds looser
    /// than `+` and `*`. Parenthesize when mixing: `(a + b) & c`, not
    /// `a + b & c`.
    fn bitand(self, rhs: Multivector<S>) -> Multivector<S> {
        Multivector {
            terms: wedge_terms(&self.terms, &rhs.terms),
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::clifford::{CliffordAlgebra, Metric};
    use crate::scalar::{Integer, Nimber};

    #[test]
    fn char2_wedge_blade_and_coefficients() {
        // Nimber world, orthogonal char-2 plane: e0∧e1 (wedge between factors).
        let alg = CliffordAlgebra::new(2, Metric::diagonal(vec![Nimber(1), Nimber(1)]));
        let e0e1 = alg.wedge(&alg.e(0), &alg.e(1));
        assert_eq!(e0e1.to_string(), "e0∧e1");
        // A non-unit nimber coefficient: *3⋅e0∧e1 (atomic, attaches bare).
        let three_e0e1 = alg.scalar_mul(&Nimber(3), &e0e1);
        assert_eq!(three_e0e1.to_string(), "*3⋅e0∧e1");
        // Zero rule: the empty multivector renders as the scalar zero's display
        // — `*0` in nim-worlds, not bare `0`.
        let zero = e0e1.clone() - e0e1;
        assert!(zero.is_zero());
        assert_eq!(zero.to_string(), "*0");
    }

    #[test]
    fn integer_grassmann_negative_and_join_rule() {
        // Integer grassmann world: -2⋅e0∧e1 (negative coefficient attaches bare;
        // the join rule lifts the leading `-`).
        let alg = CliffordAlgebra::new(2, Metric::<Integer>::grassmann(2));
        let e0e1 = alg.wedge(&alg.e(0), &alg.e(1));
        let neg2 = alg.scalar_mul(&Integer(-2), &e0e1);
        assert_eq!(neg2.to_string(), "-2⋅e0∧e1");
        // Join rule on a multi-term element: 3⋅e0 - 2⋅e1.
        let mixed = alg.scalar_mul(&Integer(3), &alg.e(0)) - alg.scalar_mul(&Integer(2), &alg.e(1));
        assert_eq!(mixed.to_string(), "3⋅e0 - 2⋅e1");
        // Coefficient-1 elision (`e0`) and -1 → `-e0` elision.
        assert_eq!(alg.e(0).to_string(), "e0");
        let neg_e0 = alg.scalar_mul(&Integer(-1), &alg.e(0));
        assert_eq!(neg_e0.to_string(), "-e0");
        // grassmann zero is char-0: renders `0`, not `*0`.
        let z = alg.e(0) - alg.e(0);
        assert_eq!(z.to_string(), "0");
    }
}
