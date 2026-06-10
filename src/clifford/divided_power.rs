//! The divided power algebra `Γ(V)` — the **char-faithful symmetric mirror** of
//! the exterior Hopf algebra in [`hopf`](crate::clifford::hopf).
//!
//! Where [`hopf`](crate::clifford::hopf) is the *antisymmetric* side (the wedge
//! product, the unshuffle coproduct with reordering signs), this is the
//! *symmetric* side. The two algebras dual to the symmetric world are the
//! **symmetric algebra** `Sym(V) = k[x_0,…,x_{n−1}]` and its graded dual the
//! **divided power algebra** `Γ(V)`. In characteristic 0 they coincide
//! (`γ^{[α]} ↔ x^α/α!`); in characteristic `p` they do **not**, and `Γ` is the
//! honest dual. We build `Γ` precisely because it is the construction that stays
//! correct in char `p`.
//!
//! A divided-power monomial is `γ^{[α]} = Π_i γ_i^{[α_i]}` indexed by a
//! multidegree `α ∈ ℕ^dim`. The Hopf structure mirrors the exterior one term for
//! term, with the signs replaced by char-sensitive integer coefficients:
//!
//! | | exterior `Λ` ([`hopf`](crate::clifford::hopf)) | divided power `Γ` (here) |
//! |---|---|---|
//! | product   | wedge, reorder **sign** | **binomial** `γ^{[α]}γ^{[β]} = Π_i \binom{α_i+β_i}{α_i} γ^{[α+β]}` |
//! | coproduct | unshuffle (signed) | **deconcatenation** `Δγ^{[α]} = Σ_{β+γ=α} γ^{[β]} ⊗ γ^{[γ]}` (sign-free) |
//! | antipode  | grade involution `(−1)^k` | grade involution `(−1)^{|α|}` |
//! | generators| primitive `e_i` | primitive `γ_i^{[1]}` |
//!
//! ## Where characteristic `p` bites (`Γ ≠ Sym`)
//!
//! The product's binomial coefficients are embedded into the scalar through the
//! ring's own arithmetic (repeated `+`, never a literal), so they reduce mod the
//! characteristic. Hence `(γ_i^{[1]})² = \binom{2}{1} γ_i^{[2]} = 2·γ_i^{[2]}`
//! **vanishes in characteristic 2** even though `γ_i^{[2]} ≠ 0`: the divided power
//! is a genuine new element, not a square of the generator. This is the exact
//! char-faithful analogue of the exterior `e_i² = 0`, and the reason `Γ` (not
//! `Sym`, where `x_i²` is a nonzero basis element) is the right object in char `p`.
//! (See the `divided_square_vanishes_in_char_two` test.)

use crate::scalar::Scalar;
use std::collections::BTreeMap;

/// A multidegree `α = (α_0,…,α_{dim−1})` — the exponent vector of a divided-power
/// monomial `γ^{[α]}`. Always stored at length `dim`.
type Multidegree = Vec<u128>;

/// One tensor term `γ^{[β]} ⊗ γ^{[γ]}` of `Γ ⊗ Γ`, keyed by the pair of
/// multidegrees.
pub type DpTensorKey = (Multidegree, Multidegree);

/// The divided power algebra `Γ(V)` on `dim` generators, the context object
/// (mirroring [`CliffordAlgebra`](crate::clifford::CliffordAlgebra)).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DividedPowerAlgebra {
    pub dim: usize,
}

/// An element of `Γ(V)`: a finite linear combination of divided-power monomials.
#[derive(Clone, Debug, PartialEq)]
pub struct DpVector<S: Scalar> {
    /// multidegree → coefficient (no zero coefficients, no off-length keys).
    pub terms: BTreeMap<Multidegree, S>,
}

/// Integer binomial coefficient `\binom{n}{k}` (exact, small arguments).
fn binom(n: u128, k: u128) -> u128 {
    if k > n {
        return 0;
    }
    let k = k.min(n - k);
    let mut acc = 1u128;
    for i in 0..k {
        acc = acc * (n - i) / (i + 1);
    }
    acc
}

/// Embed a non-negative integer into the scalar ring through its own arithmetic
/// (repeated addition of `one`), so the value reduces modulo the characteristic —
/// the char-faithful integer map `ℤ_{≥0} → S`.
fn embed_int<S: Scalar>(n: u128) -> S {
    let one = S::one();
    let mut acc = S::zero();
    for _ in 0..n {
        acc = acc.add(&one);
    }
    acc
}

impl DividedPowerAlgebra {
    pub fn new(dim: usize) -> Self {
        DividedPowerAlgebra { dim }
    }

    fn empty_degree(&self) -> Multidegree {
        vec![0u128; self.dim]
    }

    /// The zero element.
    pub fn zero<S: Scalar>(&self) -> DpVector<S> {
        DpVector {
            terms: BTreeMap::new(),
        }
    }

    /// The scalar `s` as an element of `Γ` (the empty-degree term).
    pub fn scalar<S: Scalar>(&self, s: S) -> DpVector<S> {
        let mut terms = BTreeMap::new();
        if !s.is_zero() {
            terms.insert(self.empty_degree(), s);
        }
        DpVector { terms }
    }

    /// The identity `1 = γ^{[0]}`.
    pub fn one<S: Scalar>(&self) -> DpVector<S> {
        self.scalar(S::one())
    }

    /// The divided power `γ_i^{[k]}`.
    pub fn divided_power<S: Scalar>(&self, i: usize, k: u128) -> DpVector<S> {
        assert!(i < self.dim, "generator index out of range");
        let mut deg = self.empty_degree();
        deg[i] = k;
        let mut terms = BTreeMap::new();
        if k == 0 {
            return self.one();
        }
        terms.insert(deg, S::one());
        DpVector { terms }
    }

    /// The generator `γ_i = γ_i^{[1]}` (primitive).
    pub fn gen<S: Scalar>(&self, i: usize) -> DpVector<S> {
        self.divided_power(i, 1)
    }

    /// The monomial `coeff · γ^{[α]}` from a multidegree (padded / checked).
    pub fn monomial<S: Scalar>(&self, alpha: &[u128], coeff: S) -> DpVector<S> {
        assert!(alpha.len() <= self.dim, "multidegree longer than dim");
        let mut deg = self.empty_degree();
        deg[..alpha.len()].copy_from_slice(alpha);
        let mut terms = BTreeMap::new();
        if !coeff.is_zero() {
            terms.insert(deg, coeff);
        }
        DpVector { terms }
    }

    pub fn add<S: Scalar>(&self, x: &DpVector<S>, y: &DpVector<S>) -> DpVector<S> {
        let mut terms = x.terms.clone();
        for (deg, c) in &y.terms {
            let e = terms.entry(deg.clone()).or_insert_with(S::zero);
            *e = e.add(c);
            if e.is_zero() {
                terms.remove(deg);
            }
        }
        DpVector { terms }
    }

    pub fn scalar_mul<S: Scalar>(&self, s: &S, x: &DpVector<S>) -> DpVector<S> {
        let mut terms = BTreeMap::new();
        for (deg, c) in &x.terms {
            let v = s.mul(c);
            if !v.is_zero() {
                terms.insert(deg.clone(), v);
            }
        }
        DpVector { terms }
    }

    /// The **binomial product** `γ^{[α]} · γ^{[β]} = Π_i \binom{α_i+β_i}{α_i} γ^{[α+β]}`.
    pub fn mul<S: Scalar>(&self, x: &DpVector<S>, y: &DpVector<S>) -> DpVector<S> {
        let mut terms: BTreeMap<Multidegree, S> = BTreeMap::new();
        for (a, ca) in &x.terms {
            for (b, cb) in &y.terms {
                let mut sum = self.empty_degree();
                let mut mult = 1u128;
                for i in 0..self.dim {
                    sum[i] = a[i] + b[i];
                    mult *= binom(a[i] + b[i], a[i]);
                }
                let coeff = ca.mul(cb).mul(&embed_int::<S>(mult));
                if coeff.is_zero() {
                    continue;
                }
                let e = terms.entry(sum).or_insert_with(S::zero);
                *e = e.add(&coeff);
            }
        }
        terms.retain(|_, c| !c.is_zero());
        DpVector { terms }
    }

    /// The **deconcatenation coproduct** `Δγ^{[α]} = Σ_{β+γ=α} γ^{[β]} ⊗ γ^{[γ]}`
    /// (sign-free — the symmetric mirror of the exterior unshuffle). Returned as a
    /// map `(β, γ) → coeff` over `Γ ⊗ Γ`.
    pub fn coproduct<S: Scalar>(&self, x: &DpVector<S>) -> BTreeMap<DpTensorKey, S> {
        let mut out: BTreeMap<DpTensorKey, S> = BTreeMap::new();
        for (a, c) in &x.terms {
            for beta in sub_multidegrees(a) {
                let gamma: Multidegree = a.iter().zip(&beta).map(|(ai, bi)| ai - bi).collect();
                let key = (beta, gamma);
                let e = out.entry(key.clone()).or_insert_with(S::zero);
                *e = e.add(c);
                if e.is_zero() {
                    out.remove(&key);
                }
            }
        }
        out
    }

    /// The counit `ε: Γ → S` — projection to the empty-degree (scalar) part.
    pub fn counit<S: Scalar>(&self, x: &DpVector<S>) -> S {
        x.terms
            .get(&self.empty_degree())
            .cloned()
            .unwrap_or_else(S::zero)
    }

    /// The antipode `S(γ^{[α]}) = (−1)^{|α|} γ^{[α]}` — the grade involution by
    /// total degree (sign through the scalar's own `neg`, so char 2 ⇒ identity).
    pub fn antipode<S: Scalar>(&self, x: &DpVector<S>) -> DpVector<S> {
        let mut terms = BTreeMap::new();
        for (deg, c) in &x.terms {
            let total: u128 = deg.iter().sum();
            let v = if total % 2 == 1 { c.neg() } else { c.clone() };
            terms.insert(deg.clone(), v);
        }
        DpVector { terms }
    }
}

/// All multidegrees `β` with `0 ≤ β_i ≤ α_i` componentwise (the sub-monomials).
fn sub_multidegrees(alpha: &[u128]) -> Vec<Multidegree> {
    let mut acc = vec![Vec::new()];
    for &ai in alpha {
        let mut next = Vec::new();
        for prefix in &acc {
            for bi in 0..=ai {
                let mut p = prefix.clone();
                p.push(bi);
                next.push(p);
            }
        }
        acc = next;
    }
    acc
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scalar::{Nimber, Rational};

    fn r(n: i128) -> Rational {
        Rational::int(n)
    }

    /// `(ε⊗id)∘Δ = id = (id⊗ε)∘Δ`.
    fn check_counit_law<S: Scalar>(g: &DividedPowerAlgebra, x: &DpVector<S>) {
        let empty = vec![0u128; g.dim];
        let cop = g.coproduct(x);
        let mut left = g.zero::<S>();
        let mut right = g.zero::<S>();
        for ((b, c), coeff) in &cop {
            if *b == empty {
                left = g.add(&left, &g.monomial(c, coeff.clone()));
            }
            if *c == empty {
                right = g.add(&right, &g.monomial(b, coeff.clone()));
            }
        }
        assert_eq!(&left, x, "(ε⊗id)∘Δ ≠ id");
        assert_eq!(&right, x, "(id⊗ε)∘Δ ≠ id");
    }

    /// `(Δ⊗id)∘Δ = (id⊗Δ)∘Δ`.
    fn check_coassociativity<S: Scalar>(g: &DividedPowerAlgebra, x: &DpVector<S>) {
        let cop = g.coproduct(x);
        let mut lhs: BTreeMap<(Multidegree, Multidegree, Multidegree), S> = BTreeMap::new();
        let mut rhs: BTreeMap<(Multidegree, Multidegree, Multidegree), S> = BTreeMap::new();
        for ((b, c), coeff) in &cop {
            for ((b1, b2), d) in &g.coproduct(&g.monomial(b, S::one())) {
                let key = (b1.clone(), b2.clone(), c.clone());
                let e = lhs.entry(key.clone()).or_insert_with(S::zero);
                *e = e.add(&coeff.mul(d));
                if e.is_zero() {
                    lhs.remove(&key);
                }
            }
            for ((c1, c2), d) in &g.coproduct(&g.monomial(c, S::one())) {
                let key = (b.clone(), c1.clone(), c2.clone());
                let e = rhs.entry(key.clone()).or_insert_with(S::zero);
                *e = e.add(&coeff.mul(d));
                if e.is_zero() {
                    rhs.remove(&key);
                }
            }
        }
        assert_eq!(lhs, rhs, "coproduct is not coassociative");
    }

    /// `m∘(S⊗id)∘Δ = η∘ε`.
    fn check_antipode_axiom<S: Scalar>(g: &DividedPowerAlgebra, x: &DpVector<S>) {
        let cop = g.coproduct(x);
        let mut acc = g.zero::<S>();
        for ((b, c), coeff) in &cop {
            let sb = g.antipode(&g.monomial(b, S::one()));
            let term = g.mul(&sb, &g.monomial(c, S::one()));
            acc = g.add(&acc, &g.scalar_mul(coeff, &term));
        }
        let expect = g.scalar(g.counit(x));
        assert_eq!(acc, expect, "antipode axiom failed");
    }

    fn run_axioms<S: Scalar>(g: &DividedPowerAlgebra, elts: &[DpVector<S>]) {
        for x in elts {
            check_counit_law(g, x);
            check_coassociativity(g, x);
            check_antipode_axiom(g, x);
        }
    }

    fn sample<S: Scalar>(g: &DividedPowerAlgebra) -> Vec<DpVector<S>> {
        vec![
            g.one(),
            g.gen(0),
            g.divided_power(0, 2),
            g.divided_power(1, 3),
            g.mul(&g.gen(0), &g.divided_power(1, 2)),
            g.add(&g.gen(0), &g.divided_power(1, 2)),
        ]
    }

    #[test]
    fn hopf_axioms_rational() {
        let g = DividedPowerAlgebra::new(2);
        run_axioms(&g, &sample::<Rational>(&g));
    }

    #[test]
    fn hopf_axioms_nimber() {
        // char 2: every sign is +, antipode = identity — the axioms still hold.
        let g = DividedPowerAlgebra::new(2);
        run_axioms(&g, &sample::<Nimber>(&g));
    }

    #[test]
    fn generators_are_primitive() {
        // Δ(γ_i^{[1]}) = γ_i^{[1]} ⊗ 1 + 1 ⊗ γ_i^{[1]}.
        let g = DividedPowerAlgebra::new(2);
        let cop = g.coproduct(&g.gen::<Rational>(0));
        assert_eq!(cop.len(), 2);
        assert_eq!(cop.get(&(vec![1, 0], vec![0, 0])), Some(&r(1)));
        assert_eq!(cop.get(&(vec![0, 0], vec![1, 0])), Some(&r(1)));
    }

    #[test]
    fn binomial_product_over_rationals() {
        // γ_0^{[1]} · γ_0^{[2]} = C(3,1) γ_0^{[3]} = 3 γ_0^{[3]}.
        let g = DividedPowerAlgebra::new(1);
        let prod = g.mul(&g.gen::<Rational>(0), &g.divided_power(0, 2));
        assert_eq!(prod.terms.get(&vec![3]), Some(&r(3)));
        // (γ_0^{[1]})² = C(2,1) γ_0^{[2]} = 2 γ_0^{[2]} ≠ 0 in char 0.
        let sq = g.mul(&g.gen::<Rational>(0), &g.gen(0));
        assert_eq!(sq.terms.get(&vec![2]), Some(&r(2)));
    }

    #[test]
    fn divided_square_vanishes_in_char_two() {
        // THE char-faithful signature, mirroring exterior e_i² = 0:
        // (γ_0^{[1]})² = 2·γ_0^{[2]} = 0 over a char-2 scalar...
        let g = DividedPowerAlgebra::new(1);
        let sq = g.mul(&g.gen::<Nimber>(0), &g.gen(0));
        assert!(sq.terms.is_empty(), "(γ^[1])² must vanish in char 2");
        // ...yet γ_0^{[2]} itself is a nonzero element — so Γ ≠ Sym in char 2.
        let dp2 = g.divided_power::<Nimber>(0, 2);
        assert_eq!(dp2.terms.get(&vec![2]), Some(&Nimber(1)));
        // and γ_0^{[2]} · γ_0^{[2]} = C(4,2) γ_0^{[4]} = 6 γ_0^{[4]} = 0 in char 2 too.
        let sq2 = g.mul(&dp2, &dp2);
        assert!(sq2.terms.is_empty(), "C(4,2)=6 ≡ 0 mod 2");
    }
}
