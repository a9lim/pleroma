//! The [`GlobalField`] trait — the local–global principle written **once**, over
//! the two kinds of global field this crate carries: the number field `ℚ`
//! ([`Rational`]) and the function field `F_q(t)`
//! ([`RationalFunction`]`<S>`).
//!
//! [`forms::padic`](crate::forms::padic)+[`adelic`](crate::forms::adelic) (over
//! `ℚ`) and [`forms::function_field`](crate::forms::function_field) (over
//! `F_q(t)`) were near-line-for-line parallel — the `_ff` suffix on the latter
//! existed only to dodge name collisions with the former. That parallelism is not
//! a coincidence: `ℚ` and `F_q(t)` are *the two kinds of global field*, and the
//! local–global package (places, the Hilbert symbol, the Hasse invariant,
//! reciprocity, Hasse–Minkowski) is **one** theory. This trait states it once, the
//! same move [`ResidueField`](crate::scalar::ResidueField) made for the discrete
//! Springer decomposition.
//!
//! # What is shared vs. what is per-field
//!
//! The **theorem package** — [`hasse_at_place`](GlobalField::hasse_at_place),
//! [`reciprocity_product`](GlobalField::reciprocity_product),
//! [`ramified_places`](GlobalField::ramified_places), and
//! [`is_isotropic_global`](GlobalField::is_isotropic_global) (Hasse–Minkowski) —
//! are **default methods**, written once and identical across both fields. That is
//! the symmetry made executable.
//!
//! The five **arithmetic primitives** stay per-field and delegate to the existing,
//! tested code, because the underlying arithmetic genuinely differs: `ℚ` is
//! `i128` number theory with square-free reduction and an **archimedean place**;
//! `F_q(t)` is `F_q[t]` polynomial factorization with **no archimedean place**
//! (the degree place `∞` is just another tame place). That asymmetry — the missing
//! real place in equal characteristic — is the content, not a gap. This trait is
//! deliberately **not** a [`Valued`](crate::scalar::Valued) abstraction: a global
//! field carries *all* its places at once (the same reason `RationalFunction` and
//! `Adele` are not `Valued`), so per-place residue data stays here in `forms/`.

use crate::forms::FiniteOddField;
use crate::scalar::{Rational, RationalFunction, Scalar};

/// A global field: a field with a family of places, a Hilbert symbol at each, and
/// the local–global principle (reciprocity + Hasse–Minkowski) relating them.
///
/// The implementors are the two kinds of global field:
/// [`Rational`] (`ℚ`, a number field) and [`RationalFunction`]`<S>` (`F_q(t)`, a
/// function field).
pub trait GlobalField: Scalar {
    /// A place of the field: `ℝ`/`Q_p` for `ℚ`, or `∞`/finite-`π` for `F_q(t)`.
    type Place: Clone + std::fmt::Debug + PartialEq;

    // ───────────────────── the five per-field primitives ─────────────────────

    /// The places that can carry a nontrivial local condition for `entries`
    /// (every other place sees only units): the archimedean place(s) **plus** the
    /// finite places dividing some entry.
    fn relevant_places(entries: &[Self]) -> Vec<Self::Place>;

    /// The Hilbert symbol `(a, b)_v ∈ {+1, −1}` over the completion at `place`.
    fn hilbert_symbol_at(a: &Self, b: &Self, place: &Self::Place) -> i128;

    /// Whether a **nonzero** `x` is a square in the local field at `place`.
    fn is_local_square(x: &Self, place: &Self::Place) -> bool;

    /// Whether a **nonzero** `x` is a square in the global field.
    fn is_global_square(x: &Self) -> bool;

    /// Local isotropy of the nondegenerate diagonal form `⟨a_1,…,a_n⟩` over the
    /// completion at `place`, by rank. The archimedean branch (definiteness) and
    /// the finite branch (the Serre rank conditions) live here because the
    /// archimedean place exists only over `ℚ`.
    fn is_isotropic_at_place(entries: &[Self], place: &Self::Place) -> bool;

    // ───────────────────── the local↔global theorem (defaults) ─────────────────────

    /// The Hasse invariant `ε_v(⟨a_1,…,a_n⟩) = ∏_{i<j} (a_i, a_j)_v` at `place`.
    fn hasse_at_place(entries: &[Self], place: &Self::Place) -> i128 {
        let mut h = 1i128;
        for i in 0..entries.len() {
            for j in (i + 1)..entries.len() {
                h *= Self::hilbert_symbol_at(&entries[i], &entries[j], place);
            }
        }
        h
    }

    /// The **Hilbert reciprocity product** `∏_v (a,b)_v` over all places — the
    /// product formula for the quaternion-algebra class `(a,b)`. It is `+1` for
    /// every nonzero `a, b` (Hilbert/Weil reciprocity), the gold oracle.
    fn reciprocity_product(a: &Self, b: &Self) -> i128 {
        assert!(
            !a.is_zero() && !b.is_zero(),
            "reciprocity_product needs nonzero arguments"
        );
        let pair = [a.clone(), b.clone()];
        Self::relevant_places(&pair)
            .iter()
            .fold(1i128, |acc, pl| acc * Self::hilbert_symbol_at(a, b, pl))
    }

    /// The places where the quaternion algebra `(a, b)` **ramifies** (symbol
    /// `−1`). The count is always **even** — reciprocity, additively.
    fn ramified_places(a: &Self, b: &Self) -> Vec<Self::Place> {
        assert!(
            !a.is_zero() && !b.is_zero(),
            "ramified_places needs nonzero arguments"
        );
        let pair = [a.clone(), b.clone()];
        Self::relevant_places(&pair)
            .into_iter()
            .filter(|pl| Self::hilbert_symbol_at(a, b, pl) == -1)
            .collect()
    }

    /// Whether `⟨a_1,…,a_n⟩` is **isotropic** over the global field, by
    /// **Hasse–Minkowski**: isotropic globally iff isotropic at every place. A zero
    /// entry is a null direction; rank ≤ 1 is anisotropic; rank 2 needs `−a_1a_2` a
    /// global square; rank ≥ 3 needs local isotropy at every relevant place.
    fn is_isotropic_global(entries: &[Self]) -> bool {
        if entries.iter().any(|e| e.is_zero()) {
            return true;
        }
        match entries.len() {
            0 | 1 => false,
            2 => Self::is_global_square(&entries[0].mul(&entries[1]).neg()),
            _ => Self::relevant_places(entries)
                .iter()
                .all(|pl| Self::is_isotropic_at_place(entries, pl)),
        }
    }
}

// ───────────────────────────── ℚ (number field) ─────────────────────────────

/// The integer representative of a rational's class in `ℚ*/(ℚ*)²`: `num·den`
/// (since `1/den ~ den` modulo squares). All the local/global symbols depend only
/// on this class.
fn try_rat_square_class(q: &Rational) -> Option<i128> {
    q.numer().checked_mul(q.denom())
}

fn rat_square_class(q: &Rational) -> i128 {
    try_rat_square_class(q).expect("rational square-class representative overflowed i128")
}

impl GlobalField for Rational {
    type Place = crate::forms::padic::Place;

    fn relevant_places(entries: &[Self]) -> Vec<Self::Place> {
        use crate::forms::padic::{relevant_primes, Place};
        assert!(
            entries.iter().all(|x| !x.is_zero()),
            "relevant_places over Q needs nonzero entries"
        );
        let classes: Vec<i128> = entries.iter().map(rat_square_class).collect();
        let mut places = vec![Place::Real];
        places.extend(relevant_primes(&classes).into_iter().map(Place::Prime));
        places
    }

    fn hilbert_symbol_at(a: &Self, b: &Self, place: &Self::Place) -> i128 {
        crate::forms::padic::try_hilbert_symbol_at(rat_square_class(a), rat_square_class(b), *place)
            .expect("GlobalField over Q needs bounded nonzero square classes")
    }

    fn is_local_square(x: &Self, place: &Self::Place) -> bool {
        use crate::forms::padic::{try_is_square_qp, Place};
        if x.is_zero() {
            return false;
        }
        let c = rat_square_class(x);
        match place {
            // a real number is a square iff it is ≥ 0 (its square-class sign).
            Place::Real => c >= 0,
            Place::Prime(p) => {
                try_is_square_qp(c, *p).expect("GlobalField over Q needs a bounded prime place")
            }
        }
    }

    fn is_global_square(x: &Self) -> bool {
        if x.is_zero() {
            return false;
        }
        crate::forms::padic::is_perfect_square(rat_square_class(x))
    }

    fn is_isotropic_at_place(entries: &[Self], place: &Self::Place) -> bool {
        use crate::forms::padic::{try_is_isotropic_at_p, Place};
        match place {
            // archimedean place: isotropic iff a null direction or indefinite.
            Place::Real => {
                let classes: Vec<i128> = entries.iter().map(rat_square_class).collect();
                entries.iter().any(|e| e.is_zero())
                    || (classes.iter().any(|&c| c > 0) && classes.iter().any(|&c| c < 0))
            }
            Place::Prime(p) => {
                let nz: Vec<i128> = entries
                    .iter()
                    .map(rat_square_class)
                    .filter(|&c| c != 0)
                    .collect();
                try_is_isotropic_at_p(&nz, *p)
                    .expect("GlobalField over Q needs bounded local square classes")
            }
        }
    }
}

// ───────────────────────── F_q(t) (function field) ─────────────────────────

impl<S: FiniteOddField> GlobalField for RationalFunction<S> {
    type Place = crate::forms::function_field::FFPlace<S>;

    fn relevant_places(entries: &[Self]) -> Vec<Self::Place> {
        crate::forms::function_field::relevant_places(entries)
    }

    fn hilbert_symbol_at(a: &Self, b: &Self, place: &Self::Place) -> i128 {
        crate::forms::function_field::hilbert_symbol_ff(a, b, place)
    }

    fn is_local_square(x: &Self, place: &Self::Place) -> bool {
        crate::forms::function_field::is_local_square(x, place)
    }

    fn is_global_square(x: &Self) -> bool {
        crate::forms::function_field::is_global_square_ff(x)
    }

    fn is_isotropic_at_place(entries: &[Self], place: &Self::Place) -> bool {
        crate::forms::function_field::is_isotropic_at_place(entries, place)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scalar::Fp;

    // The symmetry made executable: ONE generic test body, instantiated at both
    // kinds of global field, asserting the local↔global theorem package.

    /// Reciprocity `∏_v (a,b)_v = +1` and even ramification, over any global field.
    fn reciprocity_and_even_ramification<G: GlobalField>(samples: &[G]) {
        for a in samples {
            for b in samples {
                assert_eq!(
                    G::reciprocity_product(a, b),
                    1,
                    "reciprocity ∏_v (a,b)_v = +1 failed at a={a:?} b={b:?}"
                );
                assert_eq!(
                    G::ramified_places(a, b).len() % 2,
                    0,
                    "ramified place count must be even at a={a:?} b={b:?}"
                );
            }
        }
    }

    #[test]
    fn reciprocity_over_q() {
        let samples: Vec<Rational> = [-3, -1, 1, 2, 3, 5, 6]
            .iter()
            .map(|&n| Rational::int(n))
            .collect();
        reciprocity_and_even_ramification(&samples);
    }

    #[test]
    fn reciprocity_over_function_field() {
        type F = RationalFunction<Fp<5>>;
        let rf = |num: &[i128], den: &[i128]| -> F {
            RationalFunction::new(
                num.iter().map(|&n| Fp::<5>::new(n)).collect(),
                den.iter().map(|&n| Fp::<5>::new(n)).collect(),
            )
        };
        let samples = [
            rf(&[0, 1], &[1]),    // t
            rf(&[1, 1], &[1]),    // t+1
            rf(&[2], &[1]),       // nonsquare constant 2
            rf(&[0, 1], &[1, 1]), // t/(t+1)
            rf(&[2, 0, 1], &[1]), // t²+2 (irreducible)
        ];
        reciprocity_and_even_ramification(&samples);
    }

    #[test]
    fn global_isotropy_matches_q_field_facade() {
        // The generic Hasse-Minkowski route and the Q-specific facade are the
        // same theorem package, exposed at different abstraction levels.
        use crate::forms::try_is_isotropic_q;
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
            &[1, -1],
            &[2, -8],
            &[1, -2],
        ];
        for f in forms {
            let rats: Vec<Rational> = f.iter().map(|&n| Rational::int(n)).collect();
            assert_eq!(
                Rational::is_isotropic_global(&rats),
                try_is_isotropic_q(f).expect("test square classes fit i128"),
                "generic vs Q-specific isotropy disagree on {f:?}"
            );
        }
    }

    #[test]
    fn global_isotropy_matches_function_field_facade() {
        use crate::forms::is_isotropic_ff;
        type F = RationalFunction<Fp<5>>;
        let rf = |num: &[i128], den: &[i128]| -> F {
            RationalFunction::new(
                num.iter().map(|&n| Fp::<5>::new(n)).collect(),
                den.iter().map(|&n| Fp::<5>::new(n)).collect(),
            )
        };
        let forms: Vec<Vec<F>> = vec![
            vec![rf(&[1], &[1]), rf(&[1], &[1]), rf(&[4], &[1])],
            vec![rf(&[1], &[1]), rf(&[0, 1], &[1]), rf(&[0, 4], &[1])],
            vec![
                rf(&[1], &[1]),
                rf(&[0, 4], &[1]),
                rf(&[3], &[1]),
                rf(&[0, 2], &[1]),
            ],
            vec![rf(&[1], &[1]), rf(&[0, 4], &[1])], // rank 2 anisotropic
            vec![rf(&[1], &[1]), rf(&[0, 0, 4], &[1])], // rank 2 isotropic
        ];
        for f in &forms {
            assert_eq!(
                RationalFunction::is_isotropic_global(f),
                is_isotropic_ff(f),
                "generic vs function-field-specific isotropy disagree on {f:?}"
            );
        }
    }
}
