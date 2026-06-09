//! The **local–global** layer over the adele ring — where the reciprocity facts
//! scattered through [`padic`](crate::forms::padic) become structural statements
//! about [`A_Q`](crate::scalar::Adele).
//!
//! Three theorems, one carrier (the adele/idele):
//!
//!   * **Hilbert reciprocity** `∏_v (a,b)_v = +1` ([`hilbert_product`]) — the
//!     *multiplicative* product formula, the exact mirror of the additive
//!     `∏_v |x|_v = 1` that [`Adele::idele_norm`](crate::scalar::Adele::idele_norm)
//!     already carries. It was an oracle buried in a `padic` test; here it is a
//!     function on the global field.
//!   * **Adelic Hasse–Minkowski** ([`isotropy_over_adeles`]) — a quadratic form
//!     over `ℚ` is isotropic iff it is isotropic over *every* completion `Q_v`. The
//!     flat [`try_is_isotropic_q`](crate::forms::try_is_isotropic_q) verdict, factored into
//!     its per-place breakdown so the "⇔ over `A_Q`" structure is visible.
//!   * **The Brauer fundamental exact sequence** `0 → Br(ℚ) → ⊕_v Br(Q_v) → ℚ/ℤ → 0`
//!     ([`brauer_local_invariants`] / [`brauer_invariant_sum`]) — the local
//!     invariants of a quaternion algebra `(a,b)` are `0` or `½ ∈ ℚ/ℤ`, and their
//!     **sum is `0`**. That sum-zero law *is* Hilbert reciprocity rephrased
//!     additively (`∏ → +1` becomes `∑ → 0 mod ℤ`): an even number of ramified
//!     places.
//!
//! Reuses the `padic` Hilbert-symbol / local-isotropy machinery verbatim; nothing
//! here re-implements a symbol.

use std::collections::BTreeMap;

use crate::forms::padic::{
    relevant_primes, try_hilbert_reciprocity_product, try_hilbert_symbol_at, try_is_isotropic_at_p,
    Place,
};
use crate::scalar::{Rational, Scalar};

/// The integer representative of a rational's class in `ℚ^*/(ℚ^*)²`: `num·den`
/// (since `1/den ~ den` modulo squares). Hilbert symbols depend only on this class.
fn square_class(q: &Rational) -> Option<i128> {
    q.numer().checked_mul(q.denom())
}

// ---------------------------------------------------------------------------
// Hilbert reciprocity — the multiplicative product formula.
// ---------------------------------------------------------------------------

/// The Hilbert symbol product `∏_v (a,b)_v` over all places of `ℚ`, for the
/// quaternion algebra `(a,b)` with `a, b ∈ ℚ^*`. Equal to `+1` for all `a, b`
/// (Hilbert reciprocity) — the multiplicative analogue of the adelic product
/// formula `∏_v |x|_v = 1`.
pub fn hilbert_product(a: &Rational, b: &Rational) -> Option<i128> {
    try_hilbert_reciprocity_product(square_class(a)?, square_class(b)?)
}

// ---------------------------------------------------------------------------
// Adelic Hasse–Minkowski.
// ---------------------------------------------------------------------------

/// The per-place isotropy breakdown of a `ℚ`-form: its verdict at `ℝ` and at every
/// prime that can carry a nontrivial local condition. A form is isotropic over `ℚ`
/// iff [`is_global`](Self::is_global) — isotropic at *every* place (Hasse–Minkowski).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AdelicIsotropy {
    /// Isotropy over the Archimedean completion `ℝ`.
    pub real: bool,
    /// Isotropy over `Q_p` at each relevant prime `p`.
    pub local: BTreeMap<u128, bool>,
}

impl AdelicIsotropy {
    /// Isotropic over `ℚ` iff isotropic at every place — the local–global principle.
    pub fn is_global(&self) -> bool {
        self.real && self.local.values().all(|&b| b)
    }
}

/// The adelic Hasse–Minkowski decomposition of a diagonal integer form of **rank
/// ≥ 3**: isotropy at `ℝ` and at each relevant prime. `is_global()` then equals the
/// flat [`try_is_isotropic_q`](crate::forms::try_is_isotropic_q).
///
/// Rank ≤ 2 is excluded: there isotropy is a *global square* condition (`−ab ∈
/// (ℚ^*)²`), which is not detected by checking only the finitely many relevant
/// primes — use [`try_is_isotropic_q`](crate::forms::try_is_isotropic_q) directly.
pub fn isotropy_over_adeles(entries: &[i128]) -> Option<AdelicIsotropy> {
    assert!(
        entries.len() >= 3,
        "adelic isotropy decomposition needs rank ≥ 3 (rank ≤ 2 is a global-square \
         condition; use try_is_isotropic_q)"
    );
    // Real place: a diagonal form is isotropic over ℝ iff it has a zero entry or is
    // indefinite (entries of both signs).
    let has_zero = entries.contains(&0);
    let has_pos = entries.iter().any(|&e| e > 0);
    let has_neg = entries.iter().any(|&e| e < 0);
    let real = has_zero || (has_pos && has_neg);
    // Finite places: the local isotropy at each relevant prime (rank ≥ 3 ⇒ every
    // other prime is automatically isotropic).
    let nz: Vec<i128> = entries.iter().copied().filter(|&e| e != 0).collect();
    let local = if has_zero {
        BTreeMap::new() // a null direction makes it isotropic everywhere
    } else {
        relevant_primes(&nz)
            .into_iter()
            .map(|p| try_is_isotropic_at_p(&nz, p).map(|b| (p, b)))
            .collect::<Option<BTreeMap<_, _>>>()?
    };
    Some(AdelicIsotropy { real, local })
}

// ---------------------------------------------------------------------------
// The Brauer fundamental exact sequence.
// ---------------------------------------------------------------------------

/// The local invariants `inv_v(a,b) ∈ {0, ½} ⊂ ℚ/ℤ` of the quaternion algebra
/// `(a,b)` at every place `v`: `inv_v = ½` iff `(a,b)_v = −1` (the algebra ramifies
/// at `v`), else `0`. The image of the class in `⊕_v Br(Q_v)` of the Brauer exact
/// sequence.
pub fn brauer_local_invariants(a: &Rational, b: &Rational) -> Option<Vec<(Place, Rational)>> {
    let (ai, bi) = (square_class(a)?, square_class(b)?);
    let mut out = Vec::new();
    let zero = Rational::zero();
    let half = Rational::try_new(1, 2).expect("1/2 is a valid rational");
    let inv = |sym: i128| if sym == 1 { zero.clone() } else { half.clone() };
    out.push((
        Place::Real,
        inv(try_hilbert_symbol_at(ai, bi, Place::Real)?),
    ));
    for p in relevant_primes(&[ai, bi]) {
        out.push((
            Place::Prime(p),
            inv(try_hilbert_symbol_at(ai, bi, Place::Prime(p))?),
        ));
    }
    Some(out)
}

/// The sum `∑_v inv_v(a,b)` of the local invariants. It is always an **integer**
/// (`≡ 0` in `ℚ/ℤ`) — the exactness of the Brauer sequence, which is exactly
/// Hilbert reciprocity additively: an even number of ramified places.
pub fn brauer_invariant_sum(a: &Rational, b: &Rational) -> Option<Rational> {
    Some(
        brauer_local_invariants(a, b)?
            .into_iter()
            .fold(Rational::zero(), |acc, (_, inv)| acc.add(&inv)),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forms::try_is_isotropic_q;

    fn q(n: i128, d: i128) -> Rational {
        Rational::try_new(n, d).expect("test rational is valid")
    }

    #[test]
    fn hilbert_product_is_plus_one_reciprocity() {
        // ∏_v (a,b)_v = +1 for all rationals — the multiplicative product formula.
        for an in -6i128..=6 {
            for ad in 1i128..=4 {
                for bn in -6i128..=6 {
                    for bd in 1i128..=4 {
                        if an == 0 || bn == 0 {
                            continue;
                        }
                        assert_eq!(
                            hilbert_product(&q(an, ad), &q(bn, bd))
                                .expect("test square classes fit i128"),
                            1,
                            "reciprocity failed at ({an}/{ad}, {bn}/{bd})"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn adelic_hasse_minkowski_matches_is_isotropic_q() {
        // is_global() (isotropic at every place) ⇔ is_isotropic_q, on rank-≥3 forms.
        let forms: &[&[i128]] = &[
            &[1, 1, 1],
            &[1, 1, -1],
            &[1, 1, -3],
            &[1, 1, -2],
            &[1, 1, 1, 1],
            &[1, 1, 1, -1],
            &[1, 1, 1, 1, -1],
            &[1, 1, 1, 1, 1],
            &[2, 3, -1],
            &[1, -2, -5],
            &[3, 5, 7, -1],
        ];
        for f in forms {
            assert_eq!(
                isotropy_over_adeles(f)
                    .expect("test square classes fit i128")
                    .is_global(),
                try_is_isotropic_q(f).expect("test square classes fit i128"),
                "adelic vs global isotropy disagree on {f:?}"
            );
        }
    }

    #[test]
    fn brauer_invariant_sum_is_zero_in_q_mod_z() {
        // The sum of local invariants is an integer (≡ 0 mod ℤ) — reciprocity.
        for an in -8i128..=8 {
            for bn in -8i128..=8 {
                if an == 0 || bn == 0 {
                    continue;
                }
                let s = brauer_invariant_sum(&q(an, 1), &q(bn, 1))
                    .expect("test square classes fit i128");
                assert!(
                    s.is_integer(),
                    "Σ inv_v = {s:?} is not ≡ 0 mod ℤ for ({an},{bn})"
                );
            }
        }
    }

    #[test]
    fn hamilton_quaternions_ramify_at_2_and_infinity() {
        // (−1,−1): the canonical nontrivial class, ramified exactly at {2, ∞}, each
        // with local invariant ½ — sum = 1 ≡ 0. The discover-don't-assert version of
        // "Hamilton's quaternions ramify at 2 and ∞".
        let invs =
            brauer_local_invariants(&q(-1, 1), &q(-1, 1)).expect("test square classes fit i128");
        let ramified: Vec<Place> = invs
            .iter()
            .filter(|(_, inv)| *inv == Rational::try_new(1, 2).expect("1/2 is valid"))
            .map(|(pl, _)| *pl)
            .collect();
        assert_eq!(ramified, vec![Place::Real, Place::Prime(2)]);
        // and the real local invariant tracks the real Hilbert symbol: ½ iff a,b<0.
        let real_inv = &invs[0];
        assert_eq!(real_inv.0, Place::Real);
        assert_eq!(real_inv.1, Rational::try_new(1, 2).expect("1/2 is valid"));
        // a split example: (1, anything) is unramified everywhere ⇒ all invariants 0.
        assert_eq!(
            brauer_invariant_sum(&q(1, 1), &q(-7, 1)).expect("test square classes fit i128"),
            Rational::zero()
        );
    }
}
