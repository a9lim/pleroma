//! The [`ResidueField`] trait: a discretely-valued field that knows its **residue
//! field** `k = 𝒪/𝔪` and the reduction maps onto it.
//!
//! A non-Archimedean local field carries a whole package of data:
//!
//! ```text
//! K  — the field itself              (the Scalar)
//! 𝒪  — its ring of integers          (HasRingOfIntegers, valuation ≥ 0)
//! 𝔪  — the maximal ideal             (valuation ≥ 1)
//! Γ  — the value group               (Valued::valuation)
//! ϖ  — a uniformizer                 (Valued::uniformizer)
//! k = 𝒪/𝔪  — the residue field       ← this module
//! ```
//!
//! Until now `k` lived only in the "any number" table's *residue* column as a doc
//! comment, while every other entry of the package had been promoted to the type
//! system: the (field, ring-of-integers) pairing by
//! [`integrality`](crate::scalar::integrality), the valuation + uniformizer by
//! [`valued`](crate::scalar::valued), root-taking by
//! [`analytic`](crate::scalar::analytic). This trait closes the column the same
//! way, and it is what lets the discrete-valuation Springer decomposition
//! ([`forms::springer_decompose_local`](crate::forms::springer_decompose_local)) be written **once**,
//! generic over the residue field, instead of once per local field.
//!
//! # The two maps
//!
//! There are two natural reductions onto `k`, and the Springer second-residue map
//! needs both:
//!
//!   * [`residue`](ResidueField::residue) — the canonical reduction `𝒪 → k`,
//!     `x ↦ x mod 𝔪`. A unit reduces to its residue, a uniformizer to `0`, and a
//!     non-integral element (valuation `< 0`) has **no** residue (`None`).
//!   * [`residue_unit`](ResidueField::residue_unit) — the **angular component**
//!     `ac(x) = residue(x · ϖ^{-v(x)}) ∈ k*`: the residue of the *unit part*,
//!     defined for every nonzero `x` regardless of valuation (`None` only for `0`).
//!     This is the per-layer square-class carrier the valuation filtration reads.
//!
//! # Honest boundaries
//!
//!   * Like [`Valued`](crate::scalar::Valued), this is **not** a [`Scalar`]
//!     supertrait and excludes the globals: [`Adele`](crate::scalar::Adele) and
//!     [`RationalFunction`](crate::scalar::RationalFunction) carry *all* their
//!     places at once, so they have a residue field *per place*, not one — their
//!     residues live at the forms layer
//!     ([`forms::function_field`](crate::forms::function_field)), per place.
//!   * The impls delegate to inherent methods of the same name (inherent shadows
//!     trait in method-call position, so the delegation does not recurse), the same
//!     discipline as [`valued`](crate::scalar::valued).

use crate::scalar::{ExactFieldScalar, Fp, Fpn, Laurent, Qp, Qq, Scalar, Valued};

/// A discretely-valued field with a residue field `k = 𝒪/𝔪` and the canonical and
/// angular-component reductions onto it. The piece of the local-field package
/// `(K, 𝒪, 𝔪, k, Γ, ϖ)` that the rest of the trait layer left in doc comments.
pub trait ResidueField: Valued {
    /// The residue field `k = 𝒪/𝔪`.
    type Residue: ExactFieldScalar;

    /// The canonical reduction `𝒪 → k`, `x ↦ x mod 𝔪`. `None` for a non-integral
    /// element (valuation `< 0`); a uniformizer reduces to `0`, a unit to its
    /// residue, and `0` to `0`.
    fn residue(&self) -> Option<Self::Residue>;

    /// The **angular component** `ac(x) = residue(x · ϖ^{-v(x)}) ∈ k*`: the residue
    /// of the unit part, defined for every nonzero element (`None` only for `0`).
    /// The datum the Springer second-residue map reads at each valuation layer.
    fn residue_unit(&self) -> Option<Self::Residue>;

    /// A canonical lift `τ : k → 𝒪` of the residue map. For the local-field
    /// residue legs (`Qp`, `Qq`, `Laurent`, `Ramified`) this is the usual
    /// multiplicative Teichmüller section: `residue(τ(a)) = a` and
    /// `τ(ab) = τ(a)τ(b)` with `τ(0)=0`, not additive in general. For the
    /// residue-transcendental Gauss functor, the represented operation is the
    /// coefficientwise lift of `k(tbar)`; it is a residue section, but over a
    /// mixed-characteristic base it is not multiplicative on nonconstant
    /// rational functions.
    fn teichmuller(residue: Self::Residue) -> Self;
}

// ───────────────────────── Q_p → F_p ─────────────────────────

impl<const P: u128, const K: u128> ResidueField for Qp<P, K> {
    type Residue = Fp<P>;
    fn residue(&self) -> Option<Fp<P>> {
        match self.valuation() {
            None => Some(Fp::<P>::zero()), // 0 ↦ 0
            Some(v) if v < 0 => None,      // not integral
            Some(0) => Some(Fp::<P>::from_u128(self.unit() % P)),
            Some(_) => Some(Fp::<P>::zero()), // in 𝔪
        }
    }
    fn residue_unit(&self) -> Option<Fp<P>> {
        // `unit()` is the unit mantissa u of x = p^v·u; None only for 0.
        self.valuation()
            .map(|_| Fp::<P>::from_u128(self.unit() % P))
    }
    fn teichmuller(residue: Fp<P>) -> Self {
        Qp::teichmuller(residue)
    }
}

// ───────────────────────── Q_q → F_q ─────────────────────────
//
// The unramified leg: residue field F_q = F_{p^F}, a genuine extension residue
// (F_p only when F = 1, where Q_q is Q_p). Q_q already carries the angular
// component as `unit_residue()` — this just names it through the trait.

impl<const P: u128, const N: usize, const F: usize> ResidueField for Qq<P, N, F> {
    type Residue = Fpn<P, F>;
    fn residue(&self) -> Option<Fpn<P, F>> {
        match self.valuation() {
            None => Some(Fpn::<P, F>::zero()),
            Some(v) if v < 0 => None,
            Some(0) => self.unit_residue(),
            Some(_) => Some(Fpn::<P, F>::zero()),
        }
    }
    fn residue_unit(&self) -> Option<Fpn<P, F>> {
        self.unit_residue()
    }
    fn teichmuller(residue: Fpn<P, F>) -> Self {
        Qq::teichmuller(residue)
    }
}

// ───────────────────────── F_q((t)) → F_q ─────────────────────────
//
// The equal-characteristic leg: the residue field is the coefficient field S
// itself (reduction = "evaluate at t = 0"), and the angular component is the
// leading coefficient u₀.

impl<S: ExactFieldScalar, const K: usize> ResidueField for Laurent<S, K> {
    type Residue = S;
    fn residue(&self) -> Option<S> {
        match self.valuation() {
            None => Some(S::zero()),
            Some(v) if v < 0 => None,
            Some(0) => self.leading_coeff(), // coeff at t⁰ = the constant term
            Some(_) => Some(S::zero()),
        }
    }
    fn residue_unit(&self) -> Option<S> {
        self.leading_coeff()
    }
    fn teichmuller(residue: S) -> Self {
        Laurent::from_scalar(residue)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scalar::{Fpn, Rational};

    #[test]
    fn qp_residue_and_angular_component() {
        type Q5 = Qp<5, 4>;
        // a unit (val 0): residue = angular component = the residue digit.
        let u = Q5::from_i128(7); // 7 ≡ 2 mod 5
        assert_eq!(u.residue(), Some(Fp::<5>::new(2)));
        assert_eq!(u.residue_unit(), Some(Fp::<5>::new(2)));
        // a uniformizer multiple (val 1): residue 0, but angular component 1.
        let p = Q5::from_i128(5); // 5 = 5¹·1
        assert_eq!(p.valuation(), Some(1));
        assert_eq!(p.residue(), Some(Fp::<5>::zero()));
        assert_eq!(p.residue_unit(), Some(Fp::<5>::new(1)));
        // 1/p (val −1): not integral ⇒ no residue, but still an angular component.
        let inv_p = Q5::from_p_power(-1);
        assert_eq!(inv_p.residue(), None);
        assert_eq!(inv_p.residue_unit(), Some(Fp::<5>::new(1)));
        // zero: residue 0, no angular component.
        assert_eq!(Q5::zero().residue(), Some(Fp::<5>::zero()));
        assert_eq!(Q5::zero().residue_unit(), None);
    }

    #[test]
    fn qq_residue_lands_in_the_extension_field() {
        // Q_9 = Frac(W(F_9)): the residue field is F_9, not F_3 — a genuine
        // extension residue carried by the angular component.
        type Q9 = Qq<3, 3, 2>;
        let g = Q9::from_witt(crate::scalar::WittVec::<3, 3, 2>([0, 1])); // residue = F_9 generator
        assert_eq!(g.residue(), Some(Fpn::<3, 2>::from_coeffs(&[0, 1])));
        assert_eq!(g.residue_unit(), Some(Fpn::<3, 2>::from_coeffs(&[0, 1])));
        // scaled by p: residue collapses to 0, angular component survives.
        let pg = g.mul(&Q9::from_int(3));
        assert_eq!(pg.valuation(), Some(1));
        assert_eq!(pg.residue(), Some(Fpn::<3, 2>::zero()));
        assert_eq!(pg.residue_unit(), Some(Fpn::<3, 2>::from_coeffs(&[0, 1])));
    }

    #[test]
    fn laurent_residue_is_evaluation_at_zero() {
        type L = Laurent<Rational, 6>;
        let r = |n: i128| Rational::int(n);
        // 3 + 2t (val 0): residue = constant term 3 = angular component.
        let a = Laurent::<Rational, 6>::from_coeffs(vec![r(3), r(2)], 0);
        assert_eq!(a.residue(), Some(r(3)));
        assert_eq!(a.residue_unit(), Some(r(3)));
        // t·(…) (val 1): residue 0, angular component = leading coeff.
        let b = Laurent::<Rational, 6>::from_coeffs(vec![r(5), r(2)], 1);
        assert_eq!(b.valuation(), Some(1));
        assert_eq!(b.residue(), Some(Rational::zero()));
        assert_eq!(b.residue_unit(), Some(r(5)));
        // t⁻¹ (val −1): not integral.
        assert_eq!(L::t().inv().unwrap().residue(), None);
        assert_eq!(L::t().inv().unwrap().residue_unit(), Some(r(1)));
    }

    /// The trait is usable generically — the whole point of promoting `k` to the
    /// type system (mirrors `analytic`'s `exact_roots_is_generic`).
    #[test]
    fn residue_field_is_generic() {
        fn angular_is_some_for_nonzero<K: ResidueField>(x: &K) {
            if !x.is_zero() {
                assert!(x.residue_unit().is_some());
            }
        }
        angular_is_some_for_nonzero(&Qp::<7, 3>::from_i128(14));
        angular_is_some_for_nonzero(&Laurent::<Rational, 6>::t());
    }

    #[test]
    fn teichmuller_section_lifts_residues() {
        for a in 0..5u128 {
            let r = Fp::<5>::from_u128(a);
            let t = <Qp<5, 4> as ResidueField>::teichmuller(r);
            assert_eq!(t.residue(), Some(r));
        }

        type Q9 = Qq<3, 3, 2>;
        let g = Fpn::<3, 2>::from_coeffs(&[0, 1]);
        let tg = <Q9 as ResidueField>::teichmuller(g);
        assert_eq!(tg.residue(), Some(g));

        let a = Rational::new(3, 2);
        let ta = <Laurent<Rational, 6> as ResidueField>::teichmuller(a.clone());
        assert_eq!(ta.residue(), Some(a));
    }

    #[test]
    fn teichmuller_section_is_multiplicative_not_additive() {
        type Q5 = Qp<5, 4>;
        let a = Fp::<5>::from_u128(1);
        let b = Fp::<5>::from_u128(1);
        let ta = <Q5 as ResidueField>::teichmuller(a);
        let tb = <Q5 as ResidueField>::teichmuller(b);
        assert_eq!(
            ta.mul(&tb).residue(),
            Some(a.mul(&b)),
            "τ(ab) and τ(a)τ(b) reduce to the same residue"
        );

        let lhs = <Q5 as ResidueField>::teichmuller(a.add(&b));
        let rhs = ta.add(&tb);
        assert_eq!(lhs.residue(), rhs.residue());
        assert_ne!(lhs, rhs, "Teichmuller lifts are not generally additive");
    }
}
