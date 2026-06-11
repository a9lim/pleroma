//! Ramified (totally ramified) local extensions: adjoin a uniformizer `π` with
//! `πᴱ = ϖ` (the base uniformizer) to a [`Valued`] field — `π` is a root of an
//! Eisenstein polynomial, the mnemonic the old type name carried.
//!
//! This is the **third root-level functor**, completing the square of ways to
//! grow a field beside [`Surcomplex`](crate::scalar::Surcomplex) and
//! [`Laurent`](crate::scalar::Laurent):
//!
//!   * `Surcomplex<S>` adjoins an **algebraic, residue-extending** root (of
//!     `x² + 1`) — the *unramified* flavour.
//!   * `Laurent<S, K>` adjoins a **transcendental, value-group-extending** element
//!     `t` with a fresh valuation.
//!   * `Ramified<S, E>` adjoins an **algebraic, value-group-extending** root of
//!     the Eisenstein polynomial `xᴱ − ϖ` — the *ramified* flavour, the cell the
//!     table was missing. It refines the value group: `v(π) = 1`, `v(ϖ) = E`.
//!
//! Applied to [`Qp`](crate::scalar::Qp) it is the totally ramified extension
//! `Q_p(p^{1/E})` (the ramified twin of the unramified `Q_q`); applied to
//! `Laurent<Fpn, K>` it is `F_q((t^{1/E}))`.
//!
//! ## Always a field (Eisenstein's criterion)
//!
//! `xᴱ − ϖ` is Eisenstein at the prime `ϖ` (`v(ϖ) = 1`, so `ϖ² ∤ ϖ`), hence
//! **irreducible** over the base field for *every* `E`. So `Ramified<S, E>` is a
//! field whenever `S` is — and `inv` is total on nonzero. Over a char-0 base
//! (`Qp`, `Qq`) the extension is also separable for all `E`. Over an
//! equal-characteristic base (`Laurent` of char `p`) it is **inseparable** when
//! `p | E` (wild ramification): still a field, still computable here, but with a
//! trivial automorphism group and a degenerate trace form — a Galois-theory
//! subtlety that does not touch the ring arithmetic. (Note we cannot guard `p | E`
//! generically: the residue characteristic is hidden in the base's const params,
//! and `Qp::characteristic()` reports `0`.)
//!
//! ## Precision contract
//!
//! Every [`Valued`] base in this crate (`Qp`/`Qq`/`Laurent`) is a *capped-relative
//! precision model*, so `Ramified` over it inherits that contract: additive
//! cancellation below the retained window reads as `0`. The `E = 2` inverse is
//! exact (norm/conjugate closed form: a single scalar division, no Gaussian
//! elimination). For `E ≥ 3` the inverse is computed via Gaussian elimination over
//! the base field, so it is correct only to the retained relative precision —
//! `x · x⁻¹ = 1` up to a residual of valuation `≫ K`, not bit-exactly. Like its
//! bases it is therefore **excluded from the exact-ring fuzz suite**. The
//! *valuation* is nonetheless always exact (see [`Ramified::valuation`]).

use crate::scalar::{ResidueField, Scalar, Valued};
use std::fmt;

/// An element of `L = S(π)` with `πᴱ = ϖ`: the coordinate vector
/// `a₀ + a₁π + … + a_{E−1}π^{E−1}` over the base field `S`. The `coeffs` vector is
/// always exactly length `E` (a fixed basis — there is nothing to canonicalize but
/// the length), all-zero being the field zero.
#[derive(Clone, PartialEq, Eq, Hash)]
pub struct Ramified<S: Valued, const E: usize> {
    coeffs: Vec<S>,
}

impl<S: Valued, const E: usize> Ramified<S, E> {
    /// Build from components, padding with zeros (or truncating) to length `E`.
    pub fn new(mut coeffs: Vec<S>) -> Self {
        coeffs.resize(E, S::zero());
        Ramified { coeffs }
    }

    /// Embed a base scalar `s` as the constant `a₀ = s`.
    pub fn from_base(s: S) -> Self {
        let mut coeffs = vec![S::zero(); E];
        if E > 0 {
            coeffs[0] = s;
        }
        Ramified { coeffs }
    }

    /// The uniformizer `π` (the basis element `a₁ = 1`).
    pub fn pi() -> Self {
        debug_assert!(E >= 2, "Ramified needs E >= 2 to be a proper extension");
        let mut coeffs = vec![S::zero(); E];
        coeffs[1] = S::one();
        Ramified { coeffs }
    }

    /// The basis power `π^k` for `k < E` — i.e. the unit basis vector `e_k`.
    fn pi_basis(k: usize) -> Self {
        let mut coeffs = vec![S::zero(); E];
        coeffs[k] = S::one();
        Ramified { coeffs }
    }

    /// The (extension-normalized, `v(π) = 1`) valuation: `min_i (E·v_S(a_i) + i)`
    /// over the nonzero components, or `None` for zero.
    ///
    /// This is **exact** even over a precision-model base: the `E` quantities
    /// `E·v_S(a_i) + i` lie in distinct residue classes mod `E` (the `i` are a
    /// complete residue system), so the minimum is attained uniquely and the
    /// leading term can never cancel.
    pub fn valuation(&self) -> Option<i128> {
        let mut best: Option<i128> = None;
        for (i, a) in self.coeffs.iter().enumerate() {
            if let Some(v) = a.valuation() {
                let val = E as i128 * v + i as i128;
                best = Some(best.map_or(val, |b| b.min(val)));
            }
        }
        best
    }

    fn leading_component(&self) -> Option<&S> {
        let target = self.valuation()?;
        self.coeffs.iter().enumerate().find_map(|(i, a)| {
            a.valuation()
                .filter(|v| E as i128 * *v + i as i128 == target)
                .map(|_| a)
        })
    }

    /// Whether this lies in the ring of integers `O_S[π]` — the same-type
    /// valuation subring (valuation `≥ 0`), exactly as [`Laurent::is_integral`].
    /// So `Ramified` stays out of the [`HasRingOfIntegers`] pairing (it is a
    /// valuation functor, not an algebraic one).
    ///
    /// [`Laurent::is_integral`]: crate::scalar::Laurent::is_integral
    /// [`HasRingOfIntegers`]: crate::scalar::HasRingOfIntegers
    pub fn is_integral(&self) -> bool {
        self.valuation().is_none_or(|v| v >= 0)
    }

    /// The coordinate components `a₀ … a_{E−1}`.
    pub fn components(&self) -> &[S] {
        &self.coeffs
    }
}

impl<S: Valued, const E: usize> fmt::Display for Ramified<S, E> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_zero() {
            return write!(f, "0 (π^{E}=ϖ)");
        }
        let mut first = true;
        for (i, c) in self.coeffs.iter().enumerate() {
            if c.is_zero() {
                continue;
            }
            if !first {
                write!(f, " + ")?;
            }
            first = false;
            match i {
                0 => write!(f, "{c}")?,
                1 => write!(f, "({c})·π")?,
                _ => write!(f, "({c})·π^{i}")?,
            }
        }
        Ok(())
    }
}

impl<S: Valued, const E: usize> fmt::Debug for Ramified<S, E> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(self, f)
    }
}

impl<S: Valued, const E: usize> Scalar for Ramified<S, E> {
    fn zero() -> Self {
        Ramified {
            coeffs: vec![S::zero(); E],
        }
    }

    fn one() -> Self {
        let mut coeffs = vec![S::zero(); E];
        if E > 0 {
            coeffs[0] = S::one();
        }
        Ramified { coeffs }
    }

    fn add(&self, rhs: &Self) -> Self {
        Ramified {
            coeffs: self
                .coeffs
                .iter()
                .zip(&rhs.coeffs)
                .map(|(a, b)| a.add(b))
                .collect(),
        }
    }

    fn neg(&self) -> Self {
        Ramified {
            coeffs: self.coeffs.iter().map(|a| a.neg()).collect(),
        }
    }

    fn mul(&self, rhs: &Self) -> Self {
        // Polynomial multiplication mod (xᴱ − ϖ): an exponent k ≥ E folds once via
        // x^E = ϖ to ϖ·x^{k−E} (single pass — the largest k is 2E−2, so k−E ≤ E−2).
        let w = S::uniformizer();
        let mut out = vec![S::zero(); E];
        for (i, a) in self.coeffs.iter().enumerate() {
            if a.is_zero() {
                continue;
            }
            for (j, b) in rhs.coeffs.iter().enumerate() {
                if b.is_zero() {
                    continue;
                }
                let prod = a.mul(b);
                let k = i + j;
                if k < E {
                    out[k] = out[k].add(&prod);
                } else {
                    out[k - E] = out[k - E].add(&prod.mul(&w));
                }
            }
        }
        Ramified { coeffs: out }
    }

    fn characteristic() -> u128 {
        // Ramification preserves the characteristic.
        S::characteristic()
    }

    fn inv(&self) -> Option<Self> {
        if self.is_zero() {
            return None;
        }
        if E == 2 {
            // Norm/conjugate closed form (Surcomplex with −1 → ϖ): for α = a + bπ,
            // N = α·(a − bπ) = a² − ϖ·b² ∈ S, and α⁻¹ = (a − bπ)·N⁻¹. (Valid in
            // every characteristic: N is a scalar regardless of whether π ↦ −π is
            // an automorphism.)
            let a = &self.coeffs[0];
            let b = &self.coeffs[1];
            let w = S::uniformizer();
            let norm = a.mul(a).sub(&w.mul(&b.mul(b)));
            let ninv = norm.inv()?;
            return Some(Ramified {
                coeffs: vec![a.mul(&ninv), b.neg().mul(&ninv)],
            });
        }
        // General E: the regular representation `L_α : x ↦ α·x` is S-linear; its
        // matrix column `col` is `α·π^col` in the basis. Solve `M·c = e₀` for the
        // coordinates of α⁻¹.
        let mut m = vec![vec![S::zero(); E]; E];
        for col in 0..E {
            let prod = self.mul(&Self::pi_basis(col));
            for row in 0..E {
                m[row][col] = prod.coeffs[row].clone();
            }
        }
        let mut e0 = vec![S::zero(); E];
        e0[0] = S::one();
        let c = crate::linalg::field::solve(m, e0)?;
        Some(Ramified { coeffs: c })
    }

    fn is_zero(&self) -> bool {
        self.coeffs.iter().all(|a| a.is_zero())
    }
}

impl<S: Valued, const E: usize> Valued for Ramified<S, E> {
    fn valuation(&self) -> Option<i128> {
        Ramified::valuation(self)
    }

    /// The extension uniformizer `π`, normalized so `v_L(π)=1`.
    fn uniformizer() -> Self {
        Ramified::pi()
    }
}

impl<S: ResidueField, const E: usize> ResidueField for Ramified<S, E> {
    type Residue = S::Residue;

    fn residue(&self) -> Option<Self::Residue> {
        match self.valuation() {
            None => Some(S::Residue::zero()),
            Some(v) if v < 0 => None,
            Some(0) => self.residue_unit(),
            Some(_) => Some(S::Residue::zero()),
        }
    }

    fn residue_unit(&self) -> Option<Self::Residue> {
        self.leading_component()
            .and_then(|component| component.residue_unit())
    }

    fn teichmuller(residue: Self::Residue) -> Self {
        Ramified::from_base(S::teichmuller(residue))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scalar::{Fp, Laurent, Qp};

    // Q_3(√3): ramified quadratic, the ramified twin of Q_9.
    type E2 = Ramified<Qp<3, 6>, 2>;
    // Q_2(2^{1/3}): exercises the E ≥ 3 matrix-solve inverse.
    type E3 = Ramified<Qp<2, 8>, 3>;

    fn q3(n: i128) -> Qp<3, 6> {
        Qp::from_i128(n)
    }

    #[test]
    fn pi_power_e_is_the_uniformizer() {
        // π² = ϖ = 3.
        let pi = E2::pi();
        assert_eq!(pi.mul(&pi), E2::from_base(q3(3)));
        // π³ = ϖ·π = 3π  (fold once).
        assert_eq!(pi.mul(&pi).mul(&pi), E2::new(vec![q3(0), q3(3)]));
    }

    #[test]
    fn ramified_valuation_is_exact_and_unique() {
        // v(π) = 1, v(ϖ) = v(3) = E·1 + 0 = 2.
        assert_eq!(E2::pi().valuation(), Some(1));
        assert_eq!(E2::from_base(q3(3)).valuation(), Some(2));
        // a + bπ with v(a)=0, v(bπ)=1 ⇒ min(0,1) = 0, no tie (distinct mod 2).
        assert_eq!(E2::new(vec![q3(1), q3(1)]).valuation(), Some(0));
        // 3 + π : v(3)=2, v(π)=1 ⇒ 1.
        assert_eq!(E2::new(vec![q3(3), q3(1)]).valuation(), Some(1));
        assert_eq!(E2::zero().valuation(), None);
    }

    #[test]
    fn valuation_is_additive_under_multiplication() {
        let a = E2::new(vec![q3(1), q3(1)]); // v = 0
        let b = E2::pi(); // v = 1
        assert_eq!(a.mul(&b).valuation(), Some(1));
        let three = E2::from_base(q3(3)); // v = 2
        assert_eq!(b.mul(&three).valuation(), Some(3));
    }

    #[test]
    fn e2_inverse_round_trips_via_norm() {
        // Spread of nonzero elements over the ramified quadratic.
        for a in 1..5i128 {
            for b in 0..5i128 {
                let x = E2::new(vec![q3(a), q3(b)]);
                let xi = x.inv().expect("nonzero inverts in a field");
                assert_eq!(x.mul(&xi), E2::one(), "x·x⁻¹ ≠ 1 for {x:?}");
            }
        }
        assert_eq!(E2::zero().inv(), None);
        // 1/π = π/ϖ = π·(1/3).
        let pinv = E2::pi().inv().unwrap();
        assert_eq!(E2::pi().mul(&pinv), E2::one());
        assert_eq!(pinv.valuation(), Some(-1));
    }

    #[test]
    fn e3_inverse_round_trips_via_matrix_solve() {
        fn q2(n: i128) -> Qp<2, 8> {
            Qp::from_i128(n)
        }
        // π³ = 2 (the mul reduction, exact).
        let pi = E3::pi();
        assert_eq!(pi.mul(&pi).mul(&pi), E3::from_base(q2(2)));
        // The E ≥ 3 inverse is the regular-representation matrix solve. Unlike the
        // E = 2 norm form (a single exact division), Gaussian elimination over the
        // capped-relative Qp model only recovers the inverse to *relative
        // precision* — for a UNIT the determinant is a unit (no precision loss), so
        // x·x⁻¹ = 1 up to a residual of valuation ≫ K. We check that, not bit
        // equality (which a precision model cannot promise).
        const K: i128 = 8;
        for a in [1i128, 3] {
            // odd ⇒ component-0 is a unit ⇒ the whole element is a unit (v = 0)
            for b in 0..3i128 {
                for c in 0..3i128 {
                    let x = E3::new(vec![q2(a), q2(b), q2(c)]);
                    assert_eq!(x.valuation(), Some(0), "expected a unit: {x:?}");
                    let xi = x.inv().expect("a unit inverts");
                    let residual = x.mul(&xi).sub(&E3::one());
                    assert!(
                        residual.is_zero() || residual.valuation().unwrap() >= K,
                        "x·x⁻¹ not ≈ 1 to precision for {x:?}: residual {residual:?}"
                    );
                }
            }
        }
    }

    #[test]
    fn ring_of_integers_is_the_valuation_subring() {
        assert!(E2::pi().is_integral()); // v = 1 ≥ 0
        assert!(E2::one().is_integral());
        assert!(!E2::pi().inv().unwrap().is_integral()); // v = −1
    }

    #[test]
    fn valued_trait_uses_pi_as_uniformizer() {
        assert_eq!(<E2 as Valued>::uniformizer(), E2::pi());
        assert_eq!(<E2 as Valued>::uniformizer().valuation(), Some(1));
    }

    #[test]
    fn residue_field_is_the_base_residue() {
        assert_eq!(E2::pi().residue(), Some(Fp::<3>::zero()));
        assert_eq!(E2::pi().residue_unit(), Some(Fp::<3>::one()));
        assert_eq!(E2::from_base(q3(2)).residue(), Some(Fp::<3>::from_u128(2)));
        assert_eq!(E2::from_base(q3(3)).residue(), Some(Fp::<3>::zero()));
        assert_eq!(E2::from_base(q3(3)).residue_unit(), Some(Fp::<3>::one()));
    }

    #[test]
    fn teichmuller_lifts_base_residue_through_ramification() {
        let r = Fp::<3>::from_u128(2);
        let tau = <E2 as ResidueField>::teichmuller(r);
        assert_eq!(tau.residue(), Some(r));
        assert_eq!(tau.valuation(), Some(0));
    }

    #[test]
    fn characteristic_is_inherited_from_the_base() {
        assert_eq!(E2::characteristic(), 0); // over Q_3
        type EW = Ramified<Laurent<Fp<2>, 6>, 2>;
        assert_eq!(EW::characteristic(), 2); // over F_2((t))
    }

    #[test]
    fn wild_ramification_is_still_a_field() {
        // E = 2 over a char-2 base is WILD (p | E): x² − t is inseparable but
        // still irreducible (Eisenstein), so the extension is a field and inv is
        // total on nonzero — the separability subtlety is invisible to arithmetic.
        type EW = Ramified<Laurent<Fp<2>, 8>, 2>;
        let t = Laurent::<Fp<2>, 8>::t();
        // π² = t.
        let pi = EW::pi();
        assert_eq!(pi.mul(&pi), EW::from_base(t.clone()));
        // a spread of nonzero elements all invert (it is a genuine field).
        for a in 0..2u128 {
            for b in 0..2u128 {
                if a == 0 && b == 0 {
                    continue;
                }
                // a + t  in the first slot (a unit when a = 1, else valuation 1),
                // b in the second.
                let c0 = Laurent::from_scalar(Fp::<2>::new(a as i128)).add(&t);
                let c1 = Laurent::from_scalar(Fp::<2>::new(b as i128));
                let x = EW::new(vec![c0, c1]);
                let xi = x.inv().expect("nonzero inverts: wild extension is a field");
                assert_eq!(x.mul(&xi), EW::one(), "x·x⁻¹ ≠ 1 for {x:?}");
            }
        }
        // 1/π exists (π is a uniformizer, never a zero divisor).
        let pinv = pi.inv().expect("π inverts");
        assert_eq!(pi.mul(&pinv), EW::one());
    }
}
