//! The one piece of machinery the two transfinite backends genuinely share.
//!
//! `surreal` (`No`) and `onag` (`On₂`) both store a number as a descending
//! Conway-normal-form / Hahn series — `Vec<(exponent, coeff)>` with *recursive*
//! exponents, kept strictly descending with like powers merged and zero
//! coefficients dropped. That merge is [`merge_descending`], parameterized by the
//! **three** primitives where the two worlds actually differ:
//!
//!   1. **how exponents are ordered** — `No`'s *value* order (`a < b ⇔ b−a > 0`,
//!      a field operation: `ω−1 < ω` even though it is structurally *longer*) vs
//!      the ordinal *lexicographic* order (coefficients are positive naturals, so
//!      structure and value agree);
//!   2. **how like coefficients combine** — ordinary `ℚ` addition vs nim `XOR`;
//!   3. **which coefficients are zero**.
//!
//! That is the whole of the `surreal : nimber :: No : On₂` analogy at the code
//! level. It is deliberately *not* a shared `Cnf<C>` type: the exponent ordering
//! is field-dependent for `No` and lexicographic for `On₂`, and everything built
//! on top of it (comparison, equality, multiplication, negation) diverges
//! accordingly — `No` is a real-closed *field*, `On₂` a characteristic-2 world
//! with no negation. Sharing the *shape* without falsely sharing the *algebra*
//! is the honest unification; this function is its locus.

use std::cmp::Ordering;

/// Sort a raw `(exponent, coeff)` list into canonical descending CNF: order by
/// exponent (descending) via `exp_cmp`, merge adjacent like exponents with
/// `coeff_merge`, and drop terms whose coefficient `coeff_is_zero`.
///
/// `exp_cmp` must be a total order consistent with equality of exponents
/// (`exp_cmp(a, b) == Equal` ⟺ `a` and `b` are the same exponent), so the
/// post-sort adjacency check correctly groups like powers.
pub(crate) fn merge_descending<E, C>(
    mut raw: Vec<(E, C)>,
    exp_cmp: impl Fn(&E, &E) -> Ordering,
    coeff_merge: impl Fn(&C, &C) -> C,
    coeff_is_zero: impl Fn(&C) -> bool,
) -> Vec<(E, C)> {
    raw.sort_by(|a, b| exp_cmp(&b.0, &a.0)); // descending by exponent
    let mut out: Vec<(E, C)> = Vec::new();
    for (exp, coeff) in raw {
        if let Some(last) = out.last_mut() {
            if exp_cmp(&last.0, &exp) == Ordering::Equal {
                last.1 = coeff_merge(&last.1, &coeff);
                continue;
            }
        }
        out.push((exp, coeff));
    }
    out.retain(|(_, c)| !coeff_is_zero(c));
    out
}
