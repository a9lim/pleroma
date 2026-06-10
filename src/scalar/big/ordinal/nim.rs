//! Char-2 **nim arithmetic** on ordinals: the transfinite additive group
//! (nim-addition = XOR of like-`ω`-power coefficients) and the field product. The
//! multiplication tower (the prime-power generators `χ_{u^n}` and their carries)
//! lives in [`tower`](super::tower); this module keeps nim-addition, the finite and
//! `φ_{ω+1}` (`< ω³`) helpers used as regression oracles, and the CNF canonicalizer
//! (its like-term merge *is* nim addition / XOR; the ordinary-ordinal merge in
//! [`cantor`](super::cantor) builds its terms directly instead). See the
//! [module overview](super) for the field tower.

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

    /// Nim-multiplication across the prime-power generator tower (Conway / Lenstra /
    /// DiMuro; see the `tower` module). The non-scalar excesses `α_u` (`α_7 = ω+1`,
    /// `α_11 = ω^ω+1`, …) branch a Kummer carry into a *sum*, which is nim-multiplied
    /// back in recursively — descending by place, since every `α_{p(m)}` lives at places
    /// `< m`. Exact for every pair of ordinals `< ω^(ω^ω)` whose product triggers Kummer
    /// carries only at primes `≤ 47` (DiMuro through `43`, plus locally verified `47`).
    ///
    /// Returns `None` when an operand reaches `≥ ω^(ω^ω)` (an infinite exponent place,
    /// outside the algebraically-closed segment) **or** when a level-0 carry would need
    /// an `α_u` past the verified table (`α_53` and beyond) — the staged boundary.
    pub fn nim_mul(&self, other: &Ordinal) -> Option<Ordinal> {
        // Zero is absorbing in any field.
        if self.is_zero() || other.is_zero() {
            return Some(Ordinal::zero());
        }
        // Fast path: finite × finite is the proven `nimber::nim_mul`.
        if let (Some(a), Some(b)) = (self.as_finite(), other.as_finite()) {
            return Some(Ordinal::from_u128(nim_mul(a, b)));
        }
        // The generator tower handles the transfinite case (and its own boundary).
        super::tower::mul(self, other)
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
        // by ω, with ω³ = 2) is closed under the new nim-multiplication. The
        // independent overlap test pins the whole 64×64 product table against the
        // old reduction oracle, so this keeps exhaustive pair facts and checks the
        // triple laws on structured witnesses instead of doing a 64³ sweep.
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

        let witnesses = [
            Ordinal::zero(),
            fin(1),
            fin(2),
            fin(3),
            Ordinal::omega(),
            Ordinal::omega_pow(fin(2)),
            Ordinal::omega().nim_add(&fin(1)),
            Ordinal::omega_pow(fin(2)).nim_add(&Ordinal::omega()),
            Ordinal::from_omega3_coeffs([3, 2, 1]),
            Ordinal::from_omega3_coeffs([1, 3, 2]),
        ];
        for a in &witnesses {
            for b in &witnesses {
                let ab = a.nim_mul(b).unwrap();
                for c in &witnesses {
                    let lhs = ab.nim_mul(c).unwrap();
                    let rhs = a.nim_mul(&b.nim_mul(c).unwrap()).unwrap();
                    assert_eq!(lhs, rhs, "× not associative");
                    let lhs = a.nim_mul(&b.nim_add(c)).unwrap();
                    let rhs = ab.nim_add(&a.nim_mul(c).unwrap());
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
    fn cube_root_tower_relations() {
        // The generators gₙ = ω^(3ⁿ) and their cube-root relations gₙ³ = g_{n-1}.
        let omega = Ordinal::omega(); // g_0
        let w3 = Ordinal::omega_pow(fin(3)); // g_1 = ω^3
        let w9 = Ordinal::omega_pow(fin(9)); // g_2 = ω^9
                                             // (ω^3)² = ω^6, (ω^3) ⊗ ω = ω^4 (the worked examples)
        assert_eq!(w3.nim_mul(&w3).unwrap(), Ordinal::omega_pow(fin(6)));
        assert_eq!(w3.nim_mul(&omega).unwrap(), Ordinal::omega_pow(fin(4)));
        // g_1³ = g_0:  (ω^3)⊗³ = ω
        let w3_cubed = w3.nim_mul(&w3).unwrap().nim_mul(&w3).unwrap();
        assert_eq!(w3_cubed, omega);
        // g_2³ = g_1:  (ω^9)⊗³ = ω^3
        let w9_cubed = w9.nim_mul(&w9).unwrap().nim_mul(&w9).unwrap();
        assert_eq!(w9_cubed, w3);
    }

    #[test]
    fn consistency_with_below_omega3_path() {
        // The new tower path must agree, element-for-element, with the old
        // [c₀,c₁,c₂]-mod-(ω³−2) reduction on every pair of φ_{ω+1} elements — the
        // proof the generalization is faithful on the overlap.
        let elems: Vec<Ordinal> = (0..64u128)
            .map(|i| Ordinal::from_omega3_coeffs([i & 3, (i >> 2) & 3, (i >> 4) & 3]))
            .collect();
        for a in &elems {
            for b in &elems {
                let (ca, cb) = (a.as_below_omega3().unwrap(), b.as_below_omega3().unwrap());
                let mut p = [0u128; 5];
                for (i, &ai) in ca.iter().enumerate() {
                    for (j, &bj) in cb.iter().enumerate() {
                        p[i + j] ^= nim_mul(ai, bj);
                    }
                }
                let old = Ordinal::from_omega3_coeffs([
                    p[0] ^ nim_mul(2, p[3]),
                    p[1] ^ nim_mul(2, p[4]),
                    p[2],
                ]);
                assert_eq!(a.nim_mul(b).unwrap(), old, "tower path disagrees with old");
            }
        }
    }

    #[test]
    fn tower_multiplication_ring_axioms() {
        // The field generated by ω^3 (= g_1) is F_2(ω,ω^3) = F_{2^18} — far too big
        // to enumerate (g_0³=2 already drags in F_4, so it is *not* the naive
        // 0/1-combination of ω^0..ω^8). So the decisive Stage-A check is the
        // commutative-ring axioms on a varied sample of ordinals < ω^ω spanning
        // several generators (exponents up to 27 = 3³, i.e. g_3) and coeffs in
        // F_4 — exercising the digit-carry reduction across the whole tower.
        // Inverses/closure at the g_0 level remain exhaustively pinned by the
        // F_64 test above.
        let mut elems: Vec<Ordinal> = Vec::new();
        for &e in &[0u128, 1, 2, 3, 4, 5, 6, 8, 9, 10, 18, 27] {
            for c in 1..=3u128 {
                elems.push(Ordinal::monomial(fin(e), c));
            }
        }
        // a few genuinely multi-term ordinals.
        elems.push(Ordinal::omega().nim_add(&fin(1))); // ω + 1
        elems.push(
            Ordinal::omega_pow(fin(3))
                .nim_add(&Ordinal::omega())
                .nim_add(&fin(2)),
        ); // ω^3 + ω + 2
        elems.push(Ordinal::omega_pow(fin(9)).nim_add(&Ordinal::omega_pow(fin(3)))); // ω^9 + ω^3

        for a in &elems {
            for b in &elems {
                // every product is defined (all exponents finite ⇒ < ω^ω) …
                let ab = a.nim_mul(b).expect("< ω^ω is closed under ⊗");
                // … and commutative.
                assert_eq!(ab, b.nim_mul(a).unwrap(), "non-commutative");
                for c in &elems {
                    let l = ab.nim_mul(c).unwrap();
                    let r = a.nim_mul(&b.nim_mul(c).unwrap()).unwrap();
                    assert_eq!(l, r, "× not associative");
                    let l = a.nim_mul(&b.nim_add(c)).unwrap();
                    let r = ab.nim_add(&a.nim_mul(c).unwrap());
                    assert_eq!(l, r, "× not distributive over ⊕");
                }
            }
        }
    }

    #[test]
    fn multiplication_reaches_past_omega_omega() {
        // The prime-power tower (`tower.rs`) reaches every ordinal whose Kummer carries
        // stay at primes ≤ 47 — well past ω^ω. Spot-checks of verified landmarks
        // (full coverage lives in `tower::tests`):
        let omega = Ordinal::omega();
        let ww = Ordinal::omega_pow(omega.clone()); // ω^ω = χ_5
                                                    // ω^ω ⊗ ω = ω^(ω+1).
        assert_eq!(
            ww.nim_mul(&omega).unwrap(),
            Ordinal::omega_pow(omega.nim_add(&fin(1)))
        );
        // (ω^ω)⊗⁵ = α_5 = 4 (the quintic, scalar Kummer reduction).
        let mut p = ww.clone();
        for _ in 0..4 {
            p = p.nim_mul(&ww).unwrap();
        }
        assert_eq!(p, fin(4));
        // The Stage-2 branching: χ_7 = ω^(ω²), and (χ_7)⊗⁷ = α_7 = ω + 1 (non-scalar —
        // a single monomial branches into a sum).
        let chi7 = Ordinal::omega_pow(Ordinal::omega_pow(fin(2))); // ω^(ω²)
        let mut q = fin(1);
        for _ in 0..7 {
            q = q.nim_mul(&chi7).unwrap();
        }
        assert_eq!(q, omega.nim_add(&fin(1))); // ω + 1
                                               // The boundary now sits past prime 47, and ≥ ω^(ω^ω) is out of range.
        assert_eq!(Ordinal::omega_pow(ww.clone()).nim_mul(&omega), None); // ω^(ω^ω)
    }
}
