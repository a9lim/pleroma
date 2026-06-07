//! Char-2 **nim arithmetic** on ordinals: the transfinite additive group
//! (nim-addition = XOR of like-`ω`-power coefficients) and the field product
//! across `φ_{ω+1}` (ordinals `< ω³`) via the DiMuro/Conway tower. The CNF
//! canonicalizer lives here because its like-term merge *is* the nim addition
//! (XOR); the ordinary-ordinal merge in [`ordinal`](super::ordinal) builds its
//! terms directly instead. See the [module overview](super) for the field tower.

use super::Ordinal;
use crate::scalar::nim_mul;

/// Sort a raw term list into descending CNF and merge like `ω`-powers by **XOR**
/// (nim-addition of coefficients), dropping zeros. Exponents order by the ordinal
/// *lexicographic* order (coefficients are positive naturals, so structure and
/// value agree — unlike the surreals). The descending-merge recipe is shared with
/// the surreal backend via [`cnf::merge_descending`](super::super::cnf::merge_descending);
/// the XOR merge is exactly what makes the coefficient ring char 2.
fn canonicalize(raw: Vec<(Ordinal, u128)>) -> Vec<(Ordinal, u128)> {
    super::super::cnf::merge_descending(raw, |a, b| a.cmp(b), |x, y| x ^ y, |c| *c == 0)
}

impl Ordinal {
    /// Nim-addition: XOR the coefficients of like `ω`-powers. Complete and exact.
    pub fn nim_add(&self, other: &Ordinal) -> Ordinal {
        let mut raw = self.terms.clone();
        raw.extend(other.terms.iter().cloned());
        Ordinal {
            terms: canonicalize(raw),
        }
    }

    /// View this ordinal as an element of the field `φ_{ω+1}` (ordinals `< ω³`
    /// Cantor) — i.e. `ω²·c₂ + ω·c₁ + c₀` with each `cᵢ` finite. Returns the
    /// coefficient vector `[c₀, c₁, c₂]`, or `None` if any CNF exponent is `≥ 3`
    /// (the ordinal lives in a higher, still-staged field).
    pub fn as_below_omega3(&self) -> Option<[u128; 3]> {
        let mut coeffs = [0u128; 3];
        for (exp, c) in &self.terms {
            let e = exp.as_finite()?;
            if e >= 3 {
                return None;
            }
            coeffs[e as usize] = *c;
        }
        Some(coeffs)
    }

    /// Build the ordinal `ω²·c₂ + ω·c₁ + c₀` from its `φ_{ω+1}` coefficients.
    pub fn from_omega3_coeffs(c: [u128; 3]) -> Self {
        let mut raw = Vec::new();
        for (i, &v) in c.iter().enumerate() {
            if v != 0 {
                raw.push((Ordinal::from_u128(i as u128), v));
            }
        }
        Ordinal {
            terms: canonicalize(raw),
        }
    }

    /// Nim-multiplication. Exact across the whole field `φ_{ω+1}` (ordinals
    /// `< ω³`), via polynomial multiplication in `(finite nimbers)[ω]` modulo
    /// `ω³ − 2` — the field-tower construction of Conway / DiMuro (see the
    /// module docs). Returns `None` only when either operand has a CNF exponent
    /// `≥ 3`, i.e. lives in a higher field whose general algorithm is staged.
    pub fn nim_mul(&self, other: &Ordinal) -> Option<Ordinal> {
        // Zero is absorbing in any field.
        if self.is_zero() || other.is_zero() {
            return Some(Ordinal::zero());
        }
        // Fast path: finite × finite is the proven `nimber::nim_mul`.
        if let (Some(a), Some(b)) = (self.as_finite(), other.as_finite()) {
            return Some(Ordinal::from_u128(nim_mul(a, b)));
        }
        // Field path: both ordinals live in `φ_{ω+1}` (below ω³ Cantor).
        if let (Some(a), Some(b)) = (self.as_below_omega3(), other.as_below_omega3()) {
            // Polynomial product in ω over the finite nimbers, degree ≤ 4.
            let mut p = [0u128; 5];
            for (i, &ai) in a.iter().enumerate() {
                if ai == 0 {
                    continue;
                }
                for (j, &bj) in b.iter().enumerate() {
                    if bj != 0 {
                        p[i + j] ^= nim_mul(ai, bj);
                    }
                }
            }
            // Reduce mod ω³ − 2:  ω³ → 2,  ω⁴ → 2⊗ω.
            let c0 = p[0] ^ nim_mul(2, p[3]);
            let c1 = p[1] ^ nim_mul(2, p[4]);
            let c2 = p[2];
            return Some(Ordinal::from_omega3_coeffs([c0, c1, c2]));
        }
        // Higher fields (CNF exponent ≥ 3) — staged.
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fin(n: u128) -> Ordinal {
        Ordinal::from_u128(n)
    }

    #[test]
    fn nim_add_is_xor_below_omega() {
        for a in 0..16u128 {
            for b in 0..16u128 {
                assert_eq!(fin(a).nim_add(&fin(b)), fin(a ^ b));
            }
        }
    }

    #[test]
    fn self_inverse_and_cancellation() {
        let omega = Ordinal::omega();
        // ω ⊕ ω = 0
        assert!(omega.nim_add(&omega).is_zero());
        // (ω·3) ⊕ (ω·3) = 0
        let w3 = Ordinal::monomial(fin(1), 3);
        assert!(w3.nim_add(&w3).is_zero());
        // (ω + 1) ⊕ 1 = ω
        let w_plus_1 = omega.nim_add(&fin(1));
        assert_eq!(w_plus_1.nim_add(&fin(1)), omega);
        // ω·2 ⊕ ω = ω·3  (coefficients XOR: 2 ⊕ 1 = 3)
        let w2 = Ordinal::monomial(fin(1), 2);
        assert_eq!(w2.nim_add(&omega), Ordinal::monomial(fin(1), 3));
    }

    #[test]
    fn additive_group_axioms_with_infinite_terms() {
        let a = Ordinal::omega().nim_add(&fin(2)); // ω + 2
        let b = Ordinal::omega_pow(fin(2)).nim_add(&fin(1)); // ω² + 1
        let c = Ordinal::monomial(fin(1), 5); // ω·5
                                              // associativity + commutativity
        assert_eq!(a.nim_add(&b).nim_add(&c), a.nim_add(&b.nim_add(&c)));
        assert_eq!(a.nim_add(&b), b.nim_add(&a));
        // identity + self-inverse
        assert_eq!(a.nim_add(&Ordinal::zero()), a);
        assert!(a.nim_add(&a).is_zero());
    }

    #[test]
    fn finite_nim_mul_agrees_with_nimber() {
        for a in 0..16u128 {
            for b in 0..16u128 {
                assert_eq!(fin(a).nim_mul(&fin(b)), Some(fin(nim_mul(a, b))));
            }
        }
    }

    #[test]
    fn omega_squared_is_omega_squared() {
        // The minimal computation: ω ⊗ ω = ω² (just polynomial multiplication
        // before any reduction kicks in).
        let omega = Ordinal::omega();
        assert_eq!(omega.nim_mul(&omega).unwrap(), Ordinal::omega_pow(fin(2)));
    }

    #[test]
    fn omega_cubed_is_two() {
        // The headline (Conway/DiMuro): ω is the nim cube root of 2. This is the
        // identity that makes F_2(ω) ≅ F_8 — the cube root x³ = 2 has no
        // solution in any finite F_{2^{2^k}}, so ω supplies it.
        let omega = Ordinal::omega();
        let omega_sq = omega.nim_mul(&omega).unwrap();
        let omega_cubed = omega_sq.nim_mul(&omega).unwrap();
        assert_eq!(omega_cubed, fin(2));
        // And ω² ⊗ ω² = ω⁴ = 2⊗ω.
        assert_eq!(
            omega_sq.nim_mul(&omega_sq).unwrap(),
            Ordinal::monomial(fin(1), 2)
        );
    }

    #[test]
    fn omega_plus_one_squared_and_cubed_by_hand() {
        // (ω+1)² in characteristic 2 = ω² + 1 (cross terms vanish since 1+1=0).
        let w_plus_1 = Ordinal::omega().nim_add(&fin(1));
        let sq = w_plus_1.nim_mul(&w_plus_1).unwrap();
        assert_eq!(sq, Ordinal::omega_pow(fin(2)).nim_add(&fin(1)));
        // (ω+1)³ = (ω+1)·(ω²+1) = ω³ + ω² + ω + 1 = 2 + ω² + ω + 1 = ω² + ω + 3,
        // since nim_add(2, 1) = 2 ⊕ 1 = 3.
        let cubed = sq.nim_mul(&w_plus_1).unwrap();
        let expected = Ordinal::from_omega3_coeffs([3, 1, 1]); // ω² + ω + 3
        assert_eq!(cubed, expected);
    }

    #[test]
    fn f4_adjoin_omega_is_a_field() {
        // The decisive check: F_4(ω) = F_64 (a genuine degree-3 extension of F_4
        // by ω, with ω³ = 2) is closed under the new nim-multiplication and
        // satisfies every field axiom. 64 elements ⇒ 64² × associativity, etc.
        let elems: Vec<Ordinal> = (0..64u128)
            .map(|i| Ordinal::from_omega3_coeffs([i & 3, (i >> 2) & 3, (i >> 4) & 3]))
            .collect();
        let zero = Ordinal::zero();
        let one = fin(1);

        // closure + commutativity (and incidentally that all 64 are distinct).
        for a in &elems {
            for b in &elems {
                let ab = a.nim_mul(b).expect("F_4(ω) is closed");
                assert!(elems.iter().any(|e| e == &ab), "product escaped F_4(ω)");
                assert_eq!(ab, b.nim_mul(a).unwrap(), "non-commutative");
            }
        }

        // associativity of × + distributivity over ⊕.
        for a in &elems {
            for b in &elems {
                for c in &elems {
                    let lhs = a.nim_mul(b).unwrap().nim_mul(c).unwrap();
                    let rhs = a.nim_mul(&b.nim_mul(c).unwrap()).unwrap();
                    assert_eq!(lhs, rhs, "× not associative");
                    let lhs = a.nim_mul(&b.nim_add(c)).unwrap();
                    let rhs = a.nim_mul(b).unwrap().nim_add(&a.nim_mul(c).unwrap());
                    assert_eq!(lhs, rhs, "× not distributive over ⊕");
                }
            }
        }

        // every nonzero element has a multiplicative inverse (search the field).
        for a in elems.iter().filter(|e| !e.is_zero()) {
            let inv = elems
                .iter()
                .find(|b| a.nim_mul(b).unwrap() == one)
                .unwrap_or_else(|| panic!("no inverse for {a:?}"));
            assert_eq!(a.nim_mul(inv).unwrap(), one);
        }

        // and zero is absorbing (sanity).
        for a in &elems {
            assert_eq!(zero.nim_mul(a).unwrap(), zero);
        }
    }

    #[test]
    fn higher_field_nim_mul_is_staged() {
        // The honest boundary: any ordinal whose CNF has an exponent ≥ 3 is in a
        // higher field (degree-5 extension next, then the Lenstra/DiMuro tower
        // climbing through α_p elements) and is **not** implemented. ω³ itself,
        // ω^ω, and any product involving them returns `None`.
        let omega = Ordinal::omega();
        let omega_3 = Ordinal::omega_pow(fin(3)); // the ordinal [ω³]
        let omega_omega = Ordinal::omega_pow(omega.clone()); // [ω^ω]
        assert_eq!(omega.nim_mul(&omega_3), None);
        assert_eq!(omega_3.nim_mul(&omega), None);
        assert_eq!(omega_omega.nim_mul(&omega), None);
        assert_eq!(omega_omega.nim_mul(&omega_omega), None);
    }
}
