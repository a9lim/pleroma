//! Transfinite (ordinal) nimbers — the char-2 mirror of the surreal backend,
//! and the closure the shipped `Nimber(u128)` backend cannot reach.
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
//! ordinals, `cᵢ` finite), mirroring `surreal/` — and like there, every
//! operation recurses only on the strictly-simpler *exponents*, which is the
//! termination argument. This `mod.rs` is that CNF core (representation,
//! constructors, ordering, display); the two arithmetics live beside it:
//!
//!   * [`nim`] — the char-2 nim arithmetic: nim-addition (XOR of like-power
//!     coefficients) and the `φ_{ω+1}` field product (the DiMuro tower).
//!   * [`cantor`] — the *ordinary* (Cantor) ordinal arithmetic `ord_add`/
//!     `ord_mul` (`ω + ω = ω·2`, `1 + ω = ω`), a genuinely different operation
//!     from nim, used by the surreal birthday's run-length sums.
//!
//! ## Status (honest scope)
//!
//! * **nim-addition is complete and exact** ([`nim`]): like-`ω`-power
//!   coefficients combine by XOR (so `α ⊕ α = 0`, `ω ⊕ 1 = ω+1`), giving the
//!   genuine transfinite characteristic-2 additive group.
//! * **nim-multiplication is implemented across the whole degree-3ⁿ tower** —
//!   every ordinal strictly below **`ω^ω`** (all CNF exponents finite). Following
//!   DiMuro (*arXiv:1108.0962*, extending Conway *ONAG* ch. 6 and Lenstra 1977 "On
//!   the algebraic closure of two"): the finite layers are `F_{2^{2^n}}`, then `ω`
//!   supplies the missing cube roots (`ω³ = 2`), and the tower of cube-root
//!   generators
//!   `g₀ = ω, g₁ = ω³, g₂ = ω⁹, …, gₙ = ω^(3ⁿ)`  with  `g₀³ = 2,  gₙ³ = g_{n-1}`
//!   climbs to `ω^ω`. Every ordinal `< ω^ω` is a multivariate monomial in the `gₙ`
//!   read off the **base-3 digits** of its exponents (`ω^e = ⊗ₖ gₖ^{dₖ}`,
//!   `e = Σ dₖ·3ᵏ`), so nim-multiplication is digit-vector addition with the
//!   cube-root carries `gₖ³ = g_{k-1}`, `g₀³ = 2` (`nim::tower_mul`). This strictly
//!   subsumes the old `< ω³`, `(ω³−2)`-reduction path (the one-generator,
//!   single-digit case) — the `f4_adjoin_omega_is_a_field` (F₆₄) and
//!   `omega_cubed_is_two` checks remain green as regression. New worked relations:
//!   `(ω³)⊗³ = ω`, `(ω⁹)⊗³ = ω³`.
//! * **At `ω^ω` and above it is staged.** `ω^ω` is the first ordinal with an
//!   *infinite* CNF exponent; any such ordinal returns `None`. Reaching the full
//!   algebraic closure (the ordinals `< ω^{ω^ω}`) additionally requires the other
//!   primes: for the `(k+1)`-th prime `u` the generator is `χ_u = ω^(ω^(k-1))`
//!   (`χ_3 = ω`, `χ_5 = ω^ω`, `χ_7 = ω^(ω²)`, …), defined by the **clean Kummer
//!   relation** `(χ_u)^u = α_u` and `(χ_{u^{n+1}})^u = χ_{u^n}`, where `α_u` is
//!   Lenstra's *excess* — the smallest ordinal below `χ_u` with no `u`-th root there
//!   (DiMuro Thm 3.1.4 / Table 1: `α_3 = 2`, `α_5 = 4`, `α_7 = ω+1`,
//!   `α_11 = ω^ω+1`, `α_13 = ω+4`). The root-finding is in *computing* `α_u`, not in
//!   the generator relation. (The Artin–Schreier `x²+x+1` relation is instead the
//!   `p = u = 2` Fermat-tower case — DiMuro Thm 3.1.7 / Cor 3.11 — already handled
//!   inside the finite nimber field
//!   [`finite_field::nimber`](crate::scalar::finite_field).) The transfinite levels
//!   `ω^ω → ω^(ω²) → … → ω^(ω^ω)` (where the CNF exponents become infinite ordinals)
//!   climb prime by prime; the CNF recursion already supports them structurally. A
//!   multi-stage climb, shipped incrementally rather than speculatively; the `α_u`
//!   table is primary-source-verified — see `NOTES.md`.

mod cantor;
mod nim;

use std::cmp::Ordering;
use std::fmt;

/// An ordinal `< ε₀`-ish in Cantor normal form: `Σ ω^{exp}·coeff`, exponents
/// strictly descending, coefficients nonzero finite naturals.
#[derive(Clone, PartialEq, Eq)]
pub struct Ordinal {
    terms: Vec<(Ordinal, u128)>,
}

impl Ordinal {
    /// The ordinal `0`.
    pub fn zero() -> Self {
        Ordinal { terms: Vec::new() }
    }

    /// A finite ordinal / nimber `n`.
    pub fn from_u128(n: u128) -> Self {
        if n == 0 {
            Ordinal::zero()
        } else {
            Ordinal {
                terms: vec![(Ordinal::zero(), n)],
            }
        }
    }

    /// A single monomial `ω^exp · coeff`.
    pub fn monomial(exp: Ordinal, coeff: u128) -> Self {
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
        Ordinal::omega_pow(Ordinal::from_u128(1))
    }

    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    pub fn terms(&self) -> &[(Ordinal, u128)] {
        &self.terms
    }

    /// The ordinal order (lexicographic on descending CNF terms).
    // Inherent value-order, deliberately kept off `std::cmp::Ord`: orders and
    // operators are opt-in here, not blanket trait impls (see AGENTS.md). The
    // ordinal (lex) order also differs from the nim-value structure on the same
    // CNF, so a single std `Ord` impl would be ambiguous.
    #[allow(clippy::should_implement_trait)]
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

    /// True iff this ordinal is finite (a single `ω^0` term, or zero), returning
    /// the finite nimber value.
    pub fn as_finite(&self) -> Option<u128> {
        match self.terms.as_slice() {
            [] => Some(0),
            [(exp, c)] if exp.is_zero() => Some(*c),
            _ => None,
        }
    }
}

fn fmt_exp(e: &Ordinal) -> String {
    if e.is_zero() {
        String::new()
    } else if *e == Ordinal::from_u128(1) {
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

    fn fin(n: u128) -> Ordinal {
        Ordinal::from_u128(n)
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
    fn display_reads_as_cnf() {
        assert_eq!(format!("{:?}", Ordinal::omega()), "ω");
        assert_eq!(format!("{:?}", Ordinal::monomial(fin(1), 3)), "ω·3");
        assert_eq!(format!("{:?}", Ordinal::omega_pow(fin(2))), "ω^2");
        assert_eq!(format!("{:?}", Ordinal::omega().nim_add(&fin(1))), "ω + 1");
        assert_eq!(format!("{:?}", fin(5)), "5");
    }
}
