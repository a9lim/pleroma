//! Gauss extensions: the rational function field `S(t)` with the **Gauss
//! valuation**, over a [`Valued`] base `S`.
//!
//! This fills the **fourth and last corner** of the functor square (see
//! [`functor`](crate::scalar::functor)) — the **transcendental,
//! residue-extending** one, the twin the table was missing:
//!
//!   * [`Surcomplex`](crate::scalar::Surcomplex): algebraic, residue-extending.
//!   * [`Ramified`](crate::scalar::Ramified): algebraic, value-group-extending.
//!   * [`Laurent`](crate::scalar::Laurent): transcendental, *value*-group-extending
//!     (a fresh valuation, `v(t) = 1`).
//!   * `Gauss<S>`: transcendental, *residue*-extending — adjoin a transcendental
//!     `t` of valuation **0** whose residue `t̄` is transcendental over the residue
//!     field. The Gauss valuation `v(Σ aᵢ tⁱ) = minᵢ v_S(aᵢ)` keeps the value group
//!     of `S` unchanged and grows the residue field `k → k(t̄)`.
//!
//! So `Laurent` and `Gauss` are the two transcendental adjunctions, distinguished
//! by where the new generator lands: `Laurent`'s `t` is a uniformizer (extends the
//! value group), `Gauss`'s `t` is a unit with transcendental residue (extends the
//! residue field).
//!
//! ## Representation (no gcd — `inv` is `den/num`)
//!
//! An element is a quotient `num(t) / den(t)` of [`Poly`]s over `S`, the
//! denominator normalized monic. **No common-factor reduction is performed**: this
//! is a field by `inv(num/den) = den/num` (total on nonzero), so reduction is never
//! needed for invertibility, and gcd over a precision-model base would be unstable.
//! Equality is therefore by **cross-multiplication** (`a/b = c/d ⇔ a·d = c·b`), not
//! structural — degrees can grow under repeated operations (the polynomial analogue
//! of the relative-precision window, acceptable for a precision-model functor). The
//! polynomial arithmetic lives in the shared [`Poly`] type, not here.
//!
//! ## Precision contract
//!
//! Every in-crate [`Valued`] base (`Qp`/`Qq`/`Laurent`) is a capped-relative
//! precision model, so `Gauss` over it inherits that contract and is **excluded
//! from the exact-ring fuzz suite** (by omission, like `Laurent`/`Ramified`). The
//! Gauss *valuation* is exact whenever the base valuation is. Its ring of integers
//! — the `v ≥ 0` subring `S⟨t⟩` ([`Gauss::is_integral`]) — is the same-type
//! valuation subring, so `Gauss` stays out of the [`HasRingOfIntegers`] pairing,
//! the same honesty as `Laurent`/`Ramified`.
//!
//! [`Poly`]: crate::scalar::Poly
//! [`HasRingOfIntegers`]: crate::scalar::HasRingOfIntegers

use crate::scalar::{Poly, RationalFunction, ResidueField, Scalar, Valued};
use std::fmt;

/// An element of the rational function field `S(t)`: `num(t) / den(t)` with `den`
/// normalized monic.
#[derive(Clone)]
pub struct Gauss<S: Valued> {
    num: Poly<S>,
    den: Poly<S>,
}

impl<S: Valued> Gauss<S> {
    /// Assemble `num / den` (already-`Poly`), normalizing the denominator to monic.
    fn from_polys(num: Poly<S>, den: Poly<S>) -> Self {
        assert!(!den.is_zero(), "Gauss: zero denominator");
        if num.is_zero() {
            return Gauss {
                num: Poly::zero(),
                den: Poly::one(),
            };
        }
        // Make the denominator monic: divide both by its leading coefficient (a
        // nonzero element of the field S, hence invertible).
        let lead_inv = den
            .leading()
            .unwrap()
            .inv()
            .expect("a field's nonzero leading coefficient inverts");
        Gauss {
            num: num.scale(&lead_inv),
            den: den.scale(&lead_inv),
        }
    }

    /// Build `num / den` from low-degree-first coefficient vectors, normalizing the
    /// denominator to monic. The denominator must be nonzero. A zero numerator
    /// collapses to the canonical zero `0 / 1`.
    pub fn new(num: Vec<S>, den: Vec<S>) -> Self {
        Gauss::from_polys(Poly::new(num), Poly::new(den))
    }

    /// Embed a base scalar as the constant rational function `s / 1`.
    pub fn from_base(s: S) -> Self {
        Gauss::from_polys(Poly::constant(s), Poly::one())
    }

    /// The indeterminate `t` (a unit of valuation 0 with transcendental residue).
    pub fn t() -> Self {
        Gauss::from_polys(Poly::x(), Poly::one())
    }

    /// The numerator / denominator coefficient slices (low-degree first).
    pub fn parts(&self) -> (&[S], &[S]) {
        (self.num.coeffs(), self.den.coeffs())
    }

    /// The **Gauss valuation** `v(num) − v(den)`, or `None` for zero. Exact
    /// whenever the base valuation is (it never depends on additive cancellation).
    pub fn valuation(&self) -> Option<i128> {
        let vn = self.num.min_coeff_valuation()?; // None ⇒ zero
        let vd = self
            .den
            .min_coeff_valuation()
            .expect("denominator is nonzero");
        Some(vn - vd)
    }

    /// Whether this lies in the ring of integers `S⟨t⟩` (the `v ≥ 0` subring) — the
    /// same-type valuation subring, exactly like [`Laurent::is_integral`]. So
    /// `Gauss` stays out of the [`HasRingOfIntegers`](crate::scalar::HasRingOfIntegers)
    /// pairing.
    ///
    /// [`Laurent::is_integral`]: crate::scalar::Laurent::is_integral
    pub fn is_integral(&self) -> bool {
        self.valuation().is_none_or(|v| v >= 0)
    }
}

impl<S: Valued> PartialEq for Gauss<S> {
    /// Cross-multiplication: `a/b = c/d ⇔ a·d = c·b` (no reduced canonical form).
    fn eq(&self, other: &Self) -> bool {
        self.num.mul(&other.den) == other.num.mul(&self.den)
    }
}

impl<S: Valued> fmt::Debug for Gauss<S> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fn fmt_poly<S: Scalar>(p: &[S]) -> String {
            if p.is_empty() {
                return "0".to_string();
            }
            let mut parts = Vec::new();
            for (i, c) in p.iter().enumerate() {
                if c.is_zero() {
                    continue;
                }
                parts.push(match i {
                    0 => format!("{c:?}"),
                    1 => format!("({c:?})·t"),
                    _ => format!("({c:?})·t^{i}"),
                });
            }
            parts.join(" + ")
        }
        if self.den == Poly::one() {
            write!(f, "{}", fmt_poly(self.num.coeffs()))
        } else {
            write!(
                f,
                "[{}] / [{}]",
                fmt_poly(self.num.coeffs()),
                fmt_poly(self.den.coeffs())
            )
        }
    }
}

impl<S: Valued> Scalar for Gauss<S> {
    fn zero() -> Self {
        Gauss {
            num: Poly::zero(),
            den: Poly::one(),
        }
    }

    fn one() -> Self {
        Gauss {
            num: Poly::one(),
            den: Poly::one(),
        }
    }

    fn add(&self, rhs: &Self) -> Self {
        // a/b + c/d = (a·d + c·b) / (b·d)
        let num = self.num.mul(&rhs.den).add(&rhs.num.mul(&self.den));
        let den = self.den.mul(&rhs.den);
        Gauss::from_polys(num, den)
    }

    fn neg(&self) -> Self {
        Gauss {
            num: self.num.neg(),
            den: self.den.clone(),
        }
    }

    fn mul(&self, rhs: &Self) -> Self {
        Gauss::from_polys(self.num.mul(&rhs.num), self.den.mul(&rhs.den))
    }

    fn characteristic() -> u128 {
        S::characteristic()
    }

    fn inv(&self) -> Option<Self> {
        if self.num.is_zero() {
            return None; // zero
        }
        // (num/den)⁻¹ = den/num — total on nonzero, no gcd needed.
        Some(Gauss::from_polys(self.den.clone(), self.num.clone()))
    }

    fn is_zero(&self) -> bool {
        self.num.is_zero()
    }
}

impl<S: Valued> Valued for Gauss<S> {
    fn valuation(&self) -> Option<i128> {
        Gauss::valuation(self)
    }

    /// The base uniformizer embedded as a constant — `v(ϖ) = 1`, unchanged value
    /// group (the residue-extending signature).
    fn uniformizer() -> Self {
        Gauss::from_base(S::uniformizer())
    }
}

fn reduce_poly_at_min<S: ResidueField>(p: &Poly<S>, min_val: i128) -> Poly<S::Residue> {
    Poly::new(
        p.coeffs()
            .iter()
            .map(|c| match c.valuation() {
                Some(v) if v == min_val => c
                    .residue_unit()
                    .expect("a nonzero coefficient has an angular component"),
                _ => S::Residue::zero(),
            })
            .collect(),
    )
}

fn gauss_angular_component<S: ResidueField>(x: &Gauss<S>) -> Option<RationalFunction<S::Residue>> {
    let vn = x.num.min_coeff_valuation()?;
    let vd = x.den.min_coeff_valuation().expect("denominator is nonzero");
    let num = reduce_poly_at_min(&x.num, vn);
    let den = reduce_poly_at_min(&x.den, vd);
    Some(RationalFunction::new(
        num.coeffs().to_vec(),
        den.coeffs().to_vec(),
    ))
}

impl<S: ResidueField> ResidueField for Gauss<S> {
    type Residue = RationalFunction<S::Residue>;

    fn residue(&self) -> Option<Self::Residue> {
        match self.valuation() {
            None => Some(RationalFunction::zero()),
            Some(v) if v < 0 => None,
            Some(0) => gauss_angular_component(self),
            Some(_) => Some(RationalFunction::zero()),
        }
    }

    fn residue_unit(&self) -> Option<Self::Residue> {
        gauss_angular_component(self)
    }

    fn teichmuller(residue: Self::Residue) -> Self {
        let num = residue
            .num()
            .coeffs()
            .iter()
            .cloned()
            .map(S::teichmuller)
            .collect();
        let den = residue
            .den()
            .coeffs()
            .iter()
            .cloned()
            .map(S::teichmuller)
            .collect();
        Gauss::new(num, den)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scalar::{Fp, Laurent, Qp};

    // Q_3(t): Gauss valuation residue-extends F_3 → F_3(t̄), value group ℤ unchanged.
    type G = Gauss<Qp<3, 6>>;

    fn c(n: i128) -> Qp<3, 6> {
        Qp::from_i128(n)
    }

    #[test]
    fn t_is_a_unit_p_is_the_uniformizer() {
        // The defining contrast with Laurent: t has valuation 0 (a unit), while the
        // base prime p has valuation 1 (the value group is unchanged).
        assert_eq!(G::t().valuation(), Some(0));
        assert_eq!(<G as Valued>::uniformizer().valuation(), Some(1));
        assert_eq!(G::from_base(c(3)).valuation(), Some(1)); // v(3) = 1
        assert_eq!(G::zero().valuation(), None);
    }

    #[test]
    fn gauss_valuation_is_min_of_coefficients() {
        // v(3 + t) = min(v(3), v(1)) = min(1, 0) = 0.
        let three_plus_t = G::new(vec![c(3), c(1)], vec![c(1)]);
        assert_eq!(three_plus_t.valuation(), Some(0));
        // v((9 + 3t) / t) = min(2,1) − 0 = 1.
        let x = G::new(vec![c(9), c(3)], vec![c(0), c(1)]);
        assert_eq!(x.valuation(), Some(1));
        // valuation is additive under multiplication.
        let a = G::from_base(c(3)); // v = 1
        let b = G::new(vec![c(1), c(1)], vec![c(1)]); // 1 + t, v = 0
        assert_eq!(a.mul(&b).valuation(), Some(1));
    }

    #[test]
    fn is_a_field_inv_total_on_nonzero() {
        // A spread of rational functions; every nonzero one inverts to 1.
        let samples = [
            G::t(),
            G::from_base(c(2)),
            G::new(vec![c(1), c(1)], vec![c(1)]), // 1 + t
            G::new(vec![c(1)], vec![c(0), c(1)]), // 1/t
            G::new(vec![c(2), c(0), c(1)], vec![c(1), c(1)]), // (2 + t²)/(1 + t)
        ];
        for x in &samples {
            let xi = x.inv().expect("nonzero inverts in a field");
            assert_eq!(x.mul(&xi), G::one(), "x·x⁻¹ ≠ 1 for {x:?}");
        }
        assert_eq!(G::zero().inv(), None);
    }

    #[test]
    fn cross_multiplication_equality() {
        // t/t = 1 even though it is not structurally reduced.
        let t_over_t = G::new(vec![c(0), c(1)], vec![c(0), c(1)]);
        assert_eq!(t_over_t, G::one());
        // (2t)/(2) = t.
        let two_t_over_two = G::new(vec![c(0), c(2)], vec![c(2)]);
        assert_eq!(two_t_over_two, G::t());
        // distinct elements compare unequal.
        assert_ne!(G::t(), G::one());
    }

    #[test]
    fn ring_axioms_on_a_sample() {
        let es = [
            G::zero(),
            G::one(),
            G::t(),
            G::from_base(c(2)),
            G::new(vec![c(1), c(1)], vec![c(1)]), // 1 + t
        ];
        for a in &es {
            assert_eq!(a.add(&G::zero()), *a);
            assert_eq!(a.add(&a.neg()), G::zero());
            assert_eq!(a.mul(&G::one()), *a);
            for b in &es {
                assert_eq!(a.add(b), b.add(a));
                assert_eq!(a.mul(b), b.mul(a));
                for d in &es {
                    assert_eq!(a.add(b).add(d), a.add(&b.add(d)));
                    assert_eq!(a.mul(b).mul(d), a.mul(&b.mul(d)));
                    assert_eq!(a.mul(&b.add(d)), a.mul(b).add(&a.mul(d)));
                }
            }
        }
    }

    #[test]
    fn integrality_is_the_valuation_subring() {
        assert!(G::t().is_integral()); // v = 0
        assert!(G::from_base(c(3)).is_integral()); // v = 1
                                                   // t is a Gauss unit, so 1/t is too (v = 0) — the residue-extending signature.
        assert!(G::new(vec![c(1)], vec![c(0), c(1)]).is_integral());
        // 1/p has valuation −1 ⇒ not integral.
        assert!(!G::from_base(c(3)).inv().unwrap().is_integral());
    }

    #[test]
    fn composes_over_a_laurent_base() {
        // Gauss over an equal-characteristic local field: F_5((s))(t), residue
        // F_5(s̄)(t̄). Smoke test that the functor stacks on any Valued base.
        type GL = Gauss<Laurent<Fp<5>, 6>>;
        let t = GL::t();
        let s = GL::from_base(Laurent::<Fp<5>, 6>::t()); // the base uniformizer s
        assert_eq!(t.valuation(), Some(0)); // t is a Gauss unit
        assert_eq!(s.valuation(), Some(1)); // s carries the base valuation
        let x = t.add(&s);
        let xi = x.inv().expect("nonzero inverts");
        assert_eq!(x.mul(&xi), GL::one());
    }

    #[test]
    fn residue_field_is_rational_function_over_base_residue() {
        type R = RationalFunction<Fp<3>>;

        assert_eq!(G::t().residue(), Some(R::t()));
        assert_eq!(G::from_base(c(3)).residue(), Some(R::zero()));
        assert_eq!(
            G::from_base(c(3)).residue_unit(),
            Some(R::from_base(Fp::<3>::one()))
        );

        // 3 + 2t reduces to 2 t because the constant coefficient has higher
        // valuation; the Gauss residue field remembers tbar.
        let x = G::new(vec![c(3), c(2)], vec![c(1)]);
        assert_eq!(
            x.residue(),
            Some(RationalFunction::from_poly(Poly::monomial(
                1,
                Fp::<3>::new(2)
            )))
        );
    }

    #[test]
    fn teichmuller_lifts_rational_function_residue() {
        type R = RationalFunction<Fp<3>>;
        let r = R::new(
            vec![Fp::<3>::new(1), Fp::<3>::new(2)],
            vec![Fp::<3>::new(1)],
        );
        let tau = <G as ResidueField>::teichmuller(r.clone());
        assert_eq!(tau.residue(), Some(r));
    }
}
