//! Transfinite (ordinal) nimbers — the char-2 mirror of the surreal backend,
//! and the closure the shipped `Nimber(u64)` backend cannot reach.
//!
//! The finite nimbers form `⋃ₙ F_{2^{2^n}}` — the quadratic closure of `F₂` — but
//! this is **not** algebraically closed: it contains `F_{2^d}` only for `d` a
//! power of two, so it misses `F₈` (degree 3), `F₃₂` (degree 5), …. Conway's
//! theorem (ONAG ch. 6) is that the proper class of *all ordinals* under
//! nim-addition and nim-multiplication is an algebraically closed field of
//! characteristic 2, and the algebraic closure of `F₂` already appears among the
//! ordinals below `ω^{ω^ω}`. The first infinite ordinal `ω` supplies the missing
//! cube roots: **`ω³ = 2`** (ω is the nim-cube-root of the nimber 2), which has
//! no solution in any finite layer, so `F₂(ω)` jumps past the 2-power tower and
//! brings in `F₈`.
//!
//! An `Ordinal` is stored in Cantor normal form `Σ ω^{βᵢ}·cᵢ` (`βᵢ` descending
//! ordinals, `cᵢ` finite), mirroring `surreal.rs` — and like there, every
//! operation recurses only on the strictly-simpler *exponents*, which is the
//! termination argument.
//!
//! ## Status (honest scope)
//!
//! * **nim-addition is complete and exact**: like-`ω`-power coefficients combine
//!   by XOR (so `α ⊕ α = 0`, `ω ⊕ 1 = ω+1`), giving the genuine transfinite
//!   characteristic-2 additive group.
//! * **nim-multiplication is implemented across the whole field `φ_{ω+1}`** —
//!   every ordinal strictly below `ω³` (Cantor). Following DiMuro
//!   (*arXiv:1108.0962*, extending Conway *ONAG* ch. 6 and Lenstra 1977 "On
//!   the algebraic closure of two"): the field tower is `φ_n = F_{2^{2^n}}`
//!   (finite), `φ_ω = ω = ⋃ F_{2^{2^n}}` (still not algebraically closed —
//!   missing degree 3), and the next field `φ_{ω+1}` is obtained by adjoining
//!   `ω` as the root of the lex-earliest irreducible `x³ − 2`. DiMuro
//!   Lemma 1.1: the Cantor ordinal `[ω²·a + ω·b + c]` equals the field element
//!   `ω²⊗a ⊕ ω⊗b ⊕ c`. So nim-multiplication of any pair of ordinals
//!   `< ω³` reduces to polynomial multiplication in `(finite nimbers)[ω]` with
//!   the relations `ω³ = 2`, `ω⁴ = 2⊗ω`. The headline `ω ⊗ ω ⊗ ω = 2` and the
//!   full F_4(ω) ≅ F_64 field axioms (exhaustively checked) fall out of this.
//! * **Above `ω³` it is still staged.** The next field `φ_{ω+2}` would adjoin
//!   a degree-5 root over `φ_{ω+1}`, and the general construction climbs the
//!   Lenstra/DiMuro tower via α_p elements that require nontrivial computation
//!   in successively larger finite fields. An ordinal with any CNF exponent
//!   `≥ 3` returns `None`.

use crate::scalar::nim_mul;
use std::cmp::Ordering;
use std::fmt;

/// An ordinal `< ε₀`-ish in Cantor normal form: `Σ ω^{exp}·coeff`, exponents
/// strictly descending, coefficients nonzero finite naturals.
#[derive(Clone, PartialEq, Eq)]
pub struct Ordinal {
    terms: Vec<(Ordinal, u64)>,
}

fn canonicalize(mut raw: Vec<(Ordinal, u64)>) -> Vec<(Ordinal, u64)> {
    raw.sort_by(|a, b| b.0.cmp(&a.0)); // descending by exponent
    let mut out: Vec<(Ordinal, u64)> = Vec::new();
    for (exp, coeff) in raw {
        if let Some(last) = out.last_mut() {
            if last.0 == exp {
                last.1 ^= coeff; // nim-addition of coefficients = XOR
                continue;
            }
        }
        out.push((exp, coeff));
    }
    out.retain(|(_, c)| *c != 0);
    out
}

impl Ordinal {
    /// The ordinal `0`.
    pub fn zero() -> Self {
        Ordinal { terms: Vec::new() }
    }

    /// A finite ordinal / nimber `n`.
    pub fn from_u64(n: u64) -> Self {
        if n == 0 {
            Ordinal::zero()
        } else {
            Ordinal {
                terms: vec![(Ordinal::zero(), n)],
            }
        }
    }

    /// A single monomial `ω^exp · coeff`.
    pub fn monomial(exp: Ordinal, coeff: u64) -> Self {
        if coeff == 0 {
            Ordinal::zero()
        } else {
            Ordinal {
                terms: vec![(exp, coeff)],
            }
        }
    }

    /// `ω^exp` (coefficient 1).
    pub fn omega_pow(exp: Ordinal) -> Self {
        Ordinal::monomial(exp, 1)
    }

    /// `ω`, the first infinite ordinal.
    pub fn omega() -> Self {
        Ordinal::omega_pow(Ordinal::from_u64(1))
    }

    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    pub fn terms(&self) -> &[(Ordinal, u64)] {
        &self.terms
    }

    /// The ordinal order (lexicographic on descending CNF terms).
    pub fn cmp(&self, other: &Ordinal) -> Ordering {
        for ((e1, c1), (e2, c2)) in self.terms.iter().zip(other.terms.iter()) {
            match e1.cmp(e2) {
                Ordering::Equal => {}
                ord => return ord,
            }
            match c1.cmp(c2) {
                Ordering::Equal => {}
                ord => return ord,
            }
        }
        // shared prefix equal: the longer CNF is the larger ordinal
        self.terms.len().cmp(&other.terms.len())
    }

    /// Nim-addition: XOR the coefficients of like `ω`-powers. Complete and exact.
    pub fn nim_add(&self, other: &Ordinal) -> Ordinal {
        let mut raw = self.terms.clone();
        raw.extend(other.terms.iter().cloned());
        Ordinal {
            terms: canonicalize(raw),
        }
    }

    /// True iff this ordinal is finite (a single `ω^0` term, or zero), returning
    /// the finite nimber value.
    pub fn as_finite(&self) -> Option<u64> {
        match self.terms.as_slice() {
            [] => Some(0),
            [(exp, c)] if exp.is_zero() => Some(*c),
            _ => None,
        }
    }

    /// View this ordinal as an element of the field `φ_{ω+1}` (ordinals `< ω³`
    /// Cantor) — i.e. `ω²·c₂ + ω·c₁ + c₀` with each `cᵢ` finite. Returns the
    /// coefficient vector `[c₀, c₁, c₂]`, or `None` if any CNF exponent is `≥ 3`
    /// (the ordinal lives in a higher, still-staged field).
    pub fn as_below_omega3(&self) -> Option<[u64; 3]> {
        let mut coeffs = [0u64; 3];
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
    pub fn from_omega3_coeffs(c: [u64; 3]) -> Self {
        let mut raw = Vec::new();
        for (i, &v) in c.iter().enumerate() {
            if v != 0 {
                raw.push((Ordinal::from_u64(i as u64), v));
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
            return Some(Ordinal::from_u64(nim_mul(a, b)));
        }
        // Field path: both ordinals live in `φ_{ω+1}` (below ω³ Cantor).
        if let (Some(a), Some(b)) = (self.as_below_omega3(), other.as_below_omega3()) {
            // Polynomial product in ω over the finite nimbers, degree ≤ 4.
            let mut p = [0u64; 5];
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

fn fmt_exp(e: &Ordinal) -> String {
    if e.is_zero() {
        String::new()
    } else if *e == Ordinal::from_u64(1) {
        "ω".to_string()
    } else if e.terms.len() == 1 && e.terms[0].0.is_zero() {
        format!("ω^{}", e.terms[0].1) // ω^k for a finite exponent k
    } else {
        format!("ω^({:?})", e)
    }
}

impl fmt::Debug for Ordinal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.terms.is_empty() {
            return write!(f, "0");
        }
        let parts: Vec<String> = self
            .terms
            .iter()
            .map(|(e, c)| {
                let base = fmt_exp(e);
                if base.is_empty() {
                    format!("{}", c) // finite term
                } else if *c == 1 {
                    base
                } else {
                    format!("{}·{}", base, c)
                }
            })
            .collect();
        write!(f, "{}", parts.join(" + "))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fin(n: u64) -> Ordinal {
        Ordinal::from_u64(n)
    }

    #[test]
    fn cantor_normal_form_ordering() {
        let one = fin(1);
        let omega = Ordinal::omega(); // ω
        let omega_times_2 = Ordinal::monomial(one.clone(), 2); // ω·2
        let omega_sq = Ordinal::omega_pow(fin(2)); // ω²
        let omega_omega = Ordinal::omega_pow(Ordinal::omega()); // ω^ω
        assert_eq!(one.cmp(&omega), Ordering::Less);
        assert_eq!(omega.cmp(&omega_times_2), Ordering::Less);
        assert_eq!(omega_times_2.cmp(&omega_sq), Ordering::Less);
        assert_eq!(omega_sq.cmp(&omega_omega), Ordering::Less);
        // ω^ω dominates every ω^n
        assert_eq!(
            omega_omega.cmp(&Ordinal::omega_pow(fin(100))),
            Ordering::Greater
        );
    }

    #[test]
    fn nim_add_is_xor_below_omega() {
        for a in 0..16u64 {
            for b in 0..16u64 {
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
        for a in 0..16u64 {
            for b in 0..16u64 {
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
        let elems: Vec<Ordinal> = (0..64u64)
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

    #[test]
    fn display_reads_as_cnf() {
        assert_eq!(format!("{:?}", Ordinal::omega()), "ω");
        assert_eq!(format!("{:?}", Ordinal::monomial(fin(1), 3)), "ω·3");
        assert_eq!(format!("{:?}", Ordinal::omega_pow(fin(2))), "ω^2");
        assert_eq!(format!("{:?}", Ordinal::omega().nim_add(&fin(1))), "ω + 1");
        assert_eq!(format!("{:?}", fin(5)), "5");
    }
}
