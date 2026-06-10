//! **Ordinary (Cantor) ordinal arithmetic** — emphatically *not* nim. These are
//! the classical ordinal sum and product: non-commutative and left-absorbing
//! (`1 + ω = ω`, `2·ω = ω`), with coefficients combining as **natural numbers**
//! (`ω + ω = ω·2`), where the [nim](super::nim) operations would XOR them to `0`.
//!
//! This is the operation the surreal birthday needs (`Omnific::floor` and the
//! sign-expansion length): the length of a concatenated sign expansion is the
//! ordinary ordinal sum of the run lengths. Kept apart from the nim arithmetic
//! because the two share only the CNF representation, not the algebra.

use super::Ordinal;
use std::cmp::Ordering;

impl Ordinal {
    /// **Ordinary** (Cantor) ordinal addition `self + other` — distinct from
    /// [`nim_add`](Self::nim_add). Non-commutative and absorbing: `1 + ω = ω`
    /// (the `1` is swallowed) but `ω + 1 ≠ ω`, and crucially `ω + ω = ω·2`
    /// (coefficients add as **natural numbers**, not XOR — so this is *not*
    /// nim-addition). This is the operation the surreal birthday needs: the
    /// length of a concatenated sign expansion is the ordinary ordinal sum of
    /// the run lengths. Recurses only on `cmp` of exponents ⇒ terminates.
    pub fn ord_add(&self, other: &Ordinal) -> Ordinal {
        if other.is_zero() {
            return self.clone();
        }
        if self.is_zero() {
            return other.clone();
        }
        let beta0 = &other.terms[0].0; // other's leading exponent
        let b0 = other.terms[0].1; // other's leading coefficient
        let mut terms: Vec<(Ordinal, u128)> = Vec::new();
        let mut self_coeff_at_beta0 = 0u128;
        // Keep self's terms strictly above β₀; record a coincident β₀ term;
        // everything of self below β₀ is absorbed by other's leading power.
        for (e, c) in &self.terms {
            match e.cmp(beta0) {
                Ordering::Greater => terms.push((e.clone(), *c)),
                Ordering::Equal => self_coeff_at_beta0 = *c,
                Ordering::Less => break, // descending ⇒ all remaining are smaller
            }
        }
        // Merge at β₀: ordinary natural-number coefficient addition.
        let coeff = self_coeff_at_beta0
            .checked_add(b0)
            .expect("ordinary ordinal addition coefficient exceeds u128");
        terms.push((beta0.clone(), coeff));
        // Then the rest of other (all exponents < β₀).
        terms.extend(other.terms[1..].iter().cloned());
        Ordinal { terms }
    }

    /// **Ordinary** (Cantor) ordinal multiplication `self · other` — distinct
    /// from [`nim_mul`](Self::nim_mul). Left-distributive (`a·(b+c)=a·b+a·c`)
    /// and absorbing on the left (`2·ω = ω`, but `ω·2 = ω+ω`). Built from the
    /// standard CNF rules: `a·(ω^β·m) = ω^{α₀+β}·m` for `β>0`, and `a·n` (finite
    /// `n`) scales `a`'s leading coefficient. Uses [`ord_add`](Self::ord_add) on
    /// exponents (strictly simpler) ⇒ terminates.
    pub fn ord_mul(&self, other: &Ordinal) -> Ordinal {
        if self.is_zero() || other.is_zero() {
            return Ordinal::zero();
        }
        let alpha0 = self.terms[0].0.clone(); // leading exponent of self
        let a0 = self.terms[0].1; // leading coefficient of self
        let mut result = Ordinal::zero();
        for (beta, b) in &other.terms {
            let contribution = if beta.is_zero() {
                // Finite factor n = b:  a·n = ω^{α₀}·(a₀·n) ⊕ (rest of a).
                let mut terms = Vec::with_capacity(self.terms.len());
                terms.push((
                    alpha0.clone(),
                    a0.checked_mul(*b)
                        .expect("ordinary ordinal multiplication coefficient exceeds u128"),
                ));
                terms.extend(self.terms[1..].iter().cloned());
                Ordinal { terms }
            } else {
                // a·(ω^β·b) = ω^{α₀+β}·b  (the lower part of a vanishes).
                Ordinal::monomial(alpha0.ord_add(beta), *b)
            };
            result = result.ord_add(&contribution);
        }
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fin(n: u128) -> Ordinal {
        Ordinal::from_u128(n)
    }

    #[test]
    fn ordinary_ordinal_addition_is_not_nim() {
        let omega = Ordinal::omega();
        let one = fin(1);
        // Absorption on the left: 1 + ω = ω, but ω + 1 ≠ ω.
        assert_eq!(one.ord_add(&omega), omega);
        assert_ne!(omega.ord_add(&one), omega);
        // ω + 1 has both terms (and equals the nim-sum here, disjoint powers).
        assert_eq!(omega.ord_add(&one), omega.nim_add(&one));
        // THE distinction from nim: ω + ω = ω·2 (coeffs add as naturals),
        // whereas nim_add(ω,ω) = 0.
        assert_eq!(omega.ord_add(&omega), Ordinal::monomial(fin(1), 2));
        assert!(omega.nim_add(&omega).is_zero());
        // (ω·2 + 1) + (ω + 5): the lower `+1` of the left is absorbed by ω,
        // coefficients 2+1 = 3 at ω, then +5 ⇒ ω·3 + 5.
        let left = Ordinal::monomial(fin(1), 2).ord_add(&one); // ω·2 + 1
        let right = omega.ord_add(&fin(5)); // ω + 5
        assert_eq!(
            left.ord_add(&right),
            Ordinal::monomial(fin(1), 3).ord_add(&fin(5))
        );
        // associativity on a few triples.
        let a = omega.ord_add(&fin(2));
        let b = Ordinal::omega_pow(fin(2)).ord_add(&one);
        let c = Ordinal::monomial(fin(1), 5);
        assert_eq!(a.ord_add(&b).ord_add(&c), a.ord_add(&b.ord_add(&c)));
    }

    #[test]
    fn ordinary_ordinal_multiplication() {
        let omega = Ordinal::omega();
        // ω·2 = ω + ω = monomial(1, 2); 2·ω = ω (left absorption).
        assert_eq!(omega.ord_mul(&fin(2)), Ordinal::monomial(fin(1), 2));
        assert_eq!(fin(2).ord_mul(&omega), omega);
        // ω·ω = ω².
        assert_eq!(omega.ord_mul(&omega), Ordinal::omega_pow(fin(2)));
        // (ω+1)·2 = ω·2 + 1  (= ω+1+ω+1 = ω+ω+1).
        let w_plus_1 = omega.ord_add(&fin(1));
        assert_eq!(
            w_plus_1.ord_mul(&fin(2)),
            Ordinal::monomial(fin(1), 2).ord_add(&fin(1))
        );
        // (ω+1)·ω = ω²  (the finite tail vanishes against the ω factor).
        assert_eq!(w_plus_1.ord_mul(&omega), Ordinal::omega_pow(fin(2)));
        // associativity: (ω·ω)·ω = ω·(ω·ω) = ω³.
        let lhs = omega.ord_mul(&omega).ord_mul(&omega);
        let rhs = omega.ord_mul(&omega.ord_mul(&omega));
        assert_eq!(lhs, rhs);
        assert_eq!(lhs, Ordinal::omega_pow(fin(3)));
    }

    #[test]
    fn ordinary_ordinal_coefficients_do_not_wrap() {
        let half = 1u128 << 127;
        let a = Ordinal::monomial(fin(1), half);
        let b = Ordinal::monomial(fin(1), half);
        assert!(std::panic::catch_unwind(|| a.ord_add(&b)).is_err());
        assert!(
            std::panic::catch_unwind(|| { Ordinal::monomial(fin(1), half).ord_mul(&fin(4)) })
                .is_err()
        );
    }
}
