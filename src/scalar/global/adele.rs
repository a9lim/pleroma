//! The **adele ring** `A_Q` — modeled here as a scalar with components at every
//! rational place at once.
//!
//! Where the rest of the "any number" table picks one place — `Rational` is the
//! Archimedean place `ℝ`, each [`Qp`](crate::scalar::Qp) is one prime place — the
//! adele ring is the **restricted product** `∏'_v Q_v` over *all* places of `ℚ`
//! simultaneously (every prime `p`, plus `ℝ`), the elements `(x_∞; (x_p)_p)` with
//! `x_p` a `p`-adic integer for all but finitely many `p`.
//!
//! Representation (the *global-diagonal + corrections* model): a global rational
//! `principal` embedded **diagonally** (it is the local component at almost every
//! place), a finite map of `p`-adic **deviations** from that diagonal, and an
//! independent Archimedean component `real`. The diagonal copy of `ℚ` is the global
//! field, and the almost-all-integral condition holds for free — a rational has
//! nonzero valuation at only finitely many primes. So the image of `q ∈ ℚ` is just
//! `{ principal: q, real: q, finite: ∅ }`, and **principal adeles cost nothing**.
//!
//! `Adele` implements [`Scalar`] with characteristic 0; its `inv` is intended for
//! the **idele** group `A_Q^*` (units at every place). Like
//! [`Qp`](crate::scalar::Qp), it inherits capped-relative precision at each finite
//! place, so addition is not associative across precision boundaries. It is
//! therefore a *precision model*, **excluded from the exact-ring fuzz suite**.
//! The tested facts are the diagonal embedding, the finite-place bookkeeping, the
//! multiplicative idele behavior in represented cases, and the local–global
//! routines in [`forms::adelic`](crate::forms::adelic).
//!
//! Deliberately **not** [`Valued`](crate::scalar::Valued) (an adele has a whole
//! family of valuations, no single canonical one — use [`Adele::local_at`] and
//! [`Adele::absolute_value_at`]) and **not** in the integrality pairing (the
//! integral adeles `∏ Z_p × ℝ` are a ring but not a separate backend type — use
//! [`Adele::is_integral`]), honest gaps in the same spirit as `Laurent`.

use std::collections::{BTreeMap, BTreeSet};

use crate::scalar::{LocalQp, Rational, Scalar};

/// The nominal relative precision for an adele's finite-place cells.
pub(crate) const ADELE_PREC_NOMINAL: u128 = 16;

/// The effective relative precision at prime `p`: the nominal precision, capped so
/// the schoolbook mantissa product `(p^k)²` never overflows `u128`. A deterministic
/// function of `p`, so all cells at a given prime share one precision (and
/// [`LocalQp`] arithmetic never mixes precisions). Large primes get less precision —
/// the same `i128`-scale limitation as [`Rational`].
pub(crate) fn adele_prec(p: u128) -> u128 {
    let mut k = ADELE_PREC_NOMINAL;
    while k > 1
        && p.checked_pow(
            (2 * k)
                .try_into()
                .expect("adele precision exponent fits the platform exponent type"),
        )
        .is_none()
    {
        k -= 1;
    }
    k
}

/// A place of `ℚ`, for [`Adele::absolute_value_at`]. (Mirrors `forms::padic::Place`;
/// kept here so the scalar layer has no dependency on `forms`.)
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum AdelePlace {
    /// The Archimedean place `ℝ`.
    Real,
    /// The `p`-adic place `Q_p`.
    Prime(u128),
}

/// The primes dividing `n` (trivial trial division; `n` small in practice).
fn primes_dividing(n: i128) -> BTreeSet<u128> {
    let mut ps = BTreeSet::new();
    let mut m = n.abs();
    let mut d = 2i128;
    while d * d <= m {
        if m % d == 0 {
            ps.insert(d as u128);
            while m % d == 0 {
                m /= d;
            }
        }
        d += 1;
    }
    if m > 1 {
        ps.insert(m as u128);
    }
    ps
}

/// `p^e` as a `Rational` (`e` may be negative ⇒ `1/p^{|e|}`).
fn p_pow_rational(p: u128, e: i128) -> Rational {
    let mut acc = 1i128;
    for _ in 0..e.unsigned_abs() {
        acc = acc.checked_mul(p as i128).expect("p-power exceeds i128");
    }
    if e >= 0 {
        Rational::int(acc)
    } else {
        Rational::new(1, acc)
    }
}

/// An element of the adele ring `A_Q`.
#[derive(Clone, Debug, PartialEq)]
pub struct Adele {
    /// The global/diagonal rational — the local component at almost every place.
    principal: Rational,
    /// The Archimedean component `x_∞` (independently overridable; defaults to
    /// `principal`).
    real: Rational,
    /// Finite-place **deviations** from the diagonal: at `p ∈ finite`, the local
    /// component is `(principal in Q_p) + finite[p]`. Invariant: only genuine
    /// (nonzero) deviations are stored.
    finite: BTreeMap<u128, LocalQp>,
}

impl Adele {
    /// The diagonal embedding `ℚ ↪ A_Q`.
    pub fn from_rational(q: &Rational) -> Adele {
        Adele {
            principal: q.clone(),
            real: q.clone(),
            finite: BTreeMap::new(),
        }
    }

    /// Attach a finite-place deviation (for building non-principal adeles/ideles).
    /// A zero deviation is dropped (keeps the representation canonical).
    pub fn with_correction(mut self, p: u128, dev: LocalQp) -> Adele {
        assert_eq!(
            dev.prime(),
            p,
            "Adele correction key p={p} must match LocalQp prime {}",
            dev.prime()
        );
        assert_eq!(
            dev.precision(),
            adele_prec(p),
            "Adele correction at p={p} must use precision {}",
            adele_prec(p)
        );
        if dev.is_zero() {
            self.finite.remove(&p);
        } else {
            self.finite.insert(p, dev);
        }
        self
    }

    /// Override the Archimedean component independently of the diagonal.
    pub fn with_archimedean(mut self, real: Rational) -> Adele {
        self.real = real;
        self
    }

    /// The global/diagonal rational.
    pub fn principal(&self) -> &Rational {
        &self.principal
    }

    /// The Archimedean component `x_∞`.
    pub fn archimedean(&self) -> &Rational {
        &self.real
    }

    /// The diagonal image of `principal` at the prime `p`.
    fn diag_at(&self, p: u128) -> LocalQp {
        LocalQp::from_rational(p, adele_prec(p), &self.principal)
    }

    /// The local component `x_p ∈ Q_p` (diagonal plus any deviation).
    pub fn local_at(&self, p: u128) -> LocalQp {
        let d = self.diag_at(p);
        match self.finite.get(&p) {
            Some(dev) => d.add(dev),
            None => d,
        }
    }

    /// Whether this is a **principal** adele — the image of a global rational
    /// (no deviations, Archimedean component equal to the diagonal).
    pub fn is_principal(&self) -> bool {
        self.finite.is_empty() && self.real == self.principal
    }

    /// The finite set of places that can carry a nontrivial local condition: the
    /// primes dividing `principal`, together with the deviation primes.
    fn active_primes(&self) -> BTreeSet<u128> {
        let mut ps = primes_dividing(self.principal.numer());
        ps.extend(primes_dividing(self.principal.denom()));
        ps.extend(self.finite.keys().copied());
        ps
    }

    /// Whether this adele is an **idele** (a unit of `A_Q`): nonzero at `ℝ` and a
    /// unit (nonzero) at every finite place. Outside `active_primes` the local
    /// component is the diagonal `principal`, automatically a unit there.
    pub fn is_idele(&self) -> bool {
        if self.real.numer() == 0 || self.principal.numer() == 0 {
            return false;
        }
        self.active_primes()
            .into_iter()
            .all(|p| !self.local_at(p).is_zero())
    }

    /// Whether this adele is **integral**: a `p`-adic integer (valuation ≥ 0) at
    /// every finite place (the Archimedean component is unconstrained). The
    /// integral-adele predicate, mirroring `Laurent::is_integral`.
    pub fn is_integral(&self) -> bool {
        self.active_primes().into_iter().all(|p| {
            let x = self.local_at(p);
            x.is_zero() || x.valuation().map(|v| v >= 0).unwrap_or(true)
        })
    }

    /// The normalized absolute value `|x|_v` at a place.
    pub fn absolute_value_at(&self, place: AdelePlace) -> Rational {
        match place {
            AdelePlace::Real => {
                if self.real.sign() == std::cmp::Ordering::Less {
                    self.real.neg()
                } else {
                    self.real.clone()
                }
            }
            AdelePlace::Prime(p) => match self.local_at(p).valuation() {
                None => Rational::zero(),         // x_p = 0 ⇒ |x|_p = 0
                Some(v) => p_pow_rational(p, -v), // |x|_p = p^{-v_p(x)}
            },
        }
    }

    /// The **idele norm** `‖x‖ = ∏_v |x|_v` over all places. Identically `1` on
    /// principal ideles (the image of `ℚ^*`) — that triviality **is** the product
    /// formula (see [`Adele::satisfies_product_formula`]).
    pub fn idele_norm(&self) -> Rational {
        let mut prod = self.absolute_value_at(AdelePlace::Real);
        for p in self.active_primes() {
            prod = prod.mul(&self.absolute_value_at(AdelePlace::Prime(p)));
        }
        prod
    }

    /// Whether the **product formula** `∏_v |x|_v = 1` holds — true exactly for
    /// ideles of norm 1, in particular every principal idele.
    pub fn satisfies_product_formula(&self) -> bool {
        self.is_idele() && self.idele_norm() == Rational::one()
    }
}

impl Scalar for Adele {
    fn zero() -> Self {
        Adele::from_rational(&Rational::zero())
    }

    fn one() -> Self {
        Adele::from_rational(&Rational::one())
    }

    fn add(&self, rhs: &Self) -> Self {
        // Deviations are additive: dev_{a+b}(p) = dev_a(p) + dev_b(p).
        let principal = self.principal.add(&rhs.principal);
        let real = self.real.add(&rhs.real);
        let mut finite = BTreeMap::new();
        let keys: BTreeSet<u128> = self
            .finite
            .keys()
            .chain(rhs.finite.keys())
            .copied()
            .collect();
        for p in keys {
            let da = self
                .finite
                .get(&p)
                .copied()
                .unwrap_or_else(|| LocalQp::zero(p, adele_prec(p)));
            let db = rhs
                .finite
                .get(&p)
                .copied()
                .unwrap_or_else(|| LocalQp::zero(p, adele_prec(p)));
            let dev = da.add(&db);
            if !dev.is_zero() {
                finite.insert(p, dev);
            }
        }
        Adele {
            principal,
            real,
            finite,
        }
    }

    fn neg(&self) -> Self {
        Adele {
            principal: self.principal.neg(),
            real: self.real.neg(),
            finite: self.finite.iter().map(|(&p, d)| (p, d.neg())).collect(),
        }
    }

    fn mul(&self, rhs: &Self) -> Self {
        let principal = self.principal.mul(&rhs.principal);
        let real = self.real.mul(&rhs.real);
        let mut finite = BTreeMap::new();
        let keys: BTreeSet<u128> = self
            .finite
            .keys()
            .chain(rhs.finite.keys())
            .copied()
            .collect();
        for p in keys {
            // new local = x_p · y_p; deviation = that − (principal·principal in Q_p).
            let prod = self.local_at(p).mul(&rhs.local_at(p));
            let diag = LocalQp::from_rational(p, adele_prec(p), &principal);
            let dev = prod.add(&diag.neg());
            if !dev.is_zero() {
                finite.insert(p, dev);
            }
        }
        Adele {
            principal,
            real,
            finite,
        }
    }

    fn characteristic() -> u128 {
        0
    }

    fn inv(&self) -> Option<Self> {
        if !self.is_idele() {
            return None;
        }
        // idele ⇒ principal ≠ 0 and real ≠ 0.
        let principal = self.principal.inv()?;
        let real = self.real.inv()?;
        let mut finite = BTreeMap::new();
        for &p in self.finite.keys() {
            let lx = self.local_at(p);
            let lx_inv = lx.inv().expect("idele ⇒ nonzero at every finite place");
            let diag_inv = LocalQp::from_rational(p, adele_prec(p), &principal);
            let dev = lx_inv.add(&diag_inv.neg());
            if !dev.is_zero() {
                finite.insert(p, dev);
            }
        }
        Some(Adele {
            principal,
            real,
            finite,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn q(n: i128, d: i128) -> Rational {
        Rational::new(n, d)
    }

    #[test]
    fn diagonal_embedding_is_a_ring_homomorphism() {
        for a in -6i128..=6 {
            for b in -6i128..=6 {
                for d in 1i128..=4 {
                    let (ra, rb) = (q(a, d), q(b, d));
                    let (ea, eb) = (Adele::from_rational(&ra), Adele::from_rational(&rb));
                    assert_eq!(ea.add(&eb), Adele::from_rational(&ra.add(&rb)));
                    assert_eq!(ea.mul(&eb), Adele::from_rational(&ra.mul(&rb)));
                }
            }
        }
        assert_eq!(Adele::zero(), Adele::from_rational(&Rational::zero()));
        assert_eq!(Adele::one(), Adele::from_rational(&Rational::one()));
    }

    #[test]
    fn inv_is_total_on_principal_ideles() {
        for n in -10i128..=10 {
            for d in 1i128..=6 {
                if n == 0 {
                    continue;
                }
                let x = Adele::from_rational(&q(n, d));
                assert!(x.is_idele());
                let xi = x.inv().expect("nonzero rational is an idele");
                assert_eq!(xi, Adele::from_rational(&q(d, n)));
                assert_eq!(x.mul(&xi), Adele::one());
            }
        }
        // zero and a non-idele are not invertible.
        assert_eq!(Adele::zero().inv(), None);
    }

    #[test]
    fn product_formula_on_principal_ideles() {
        // ∏_v |q|_v = 1 for every rational q ∈ ℚ*.
        for n in -12i128..=12 {
            for d in 1i128..=8 {
                if n == 0 {
                    continue;
                }
                let x = Adele::from_rational(&q(n, d));
                assert_eq!(x.idele_norm(), Rational::one(), "‖{n}/{d}‖ ≠ 1");
                assert!(x.satisfies_product_formula());
            }
        }
    }

    #[test]
    fn absolute_values_factor_the_rational() {
        // |12/5|_∞ = 12/5; |12/5|_2 = 1/4 (v=2); |12/5|_3 = 1/3 (v=1);
        // |12/5|_5 = 5 (v=-1); product = 1.
        let x = Adele::from_rational(&q(12, 5));
        assert_eq!(x.absolute_value_at(AdelePlace::Real), q(12, 5));
        assert_eq!(x.absolute_value_at(AdelePlace::Prime(2)), q(1, 4));
        assert_eq!(x.absolute_value_at(AdelePlace::Prime(3)), q(1, 3));
        assert_eq!(x.absolute_value_at(AdelePlace::Prime(5)), Rational::int(5));
        assert_eq!(x.absolute_value_at(AdelePlace::Prime(7)), Rational::one());
        assert_eq!(x.idele_norm(), Rational::one());
    }

    #[test]
    fn a_nonprincipal_idele_and_its_inverse() {
        // Perturb at p = 7: a genuine deviation that keeps the element a unit there.
        let dev = LocalQp::from_i128(7, adele_prec(7), 1); // small unit deviation
        let x = Adele::from_rational(&q(2, 3)).with_correction(7, dev);
        assert!(!x.is_principal());
        assert!(x.is_idele());
        let xi = x.inv().expect("idele inverts");
        assert_eq!(x.mul(&xi), Adele::one());
        // and the local components are genuine inverses in Q_7.
        assert_eq!(
            x.local_at(7).mul(&xi.local_at(7)),
            LocalQp::one(7, adele_prec(7))
        );
    }

    #[test]
    fn a_correction_to_zero_is_not_an_idele() {
        // Build a deviation that cancels the diagonal at p = 5 ⇒ x_5 = 0.
        let diag5 = LocalQp::from_rational(5, adele_prec(5), &q(2, 1));
        let x = Adele::from_rational(&q(2, 1)).with_correction(5, diag5.neg());
        assert!(x.local_at(5).is_zero());
        assert!(!x.is_idele());
        assert_eq!(x.inv(), None);
    }

    #[test]
    #[should_panic(expected = "must match LocalQp prime")]
    fn correction_prime_must_match_key() {
        let dev = LocalQp::zero(5, adele_prec(5));
        let _ = Adele::one().with_correction(7, dev);
    }

    #[test]
    #[should_panic(expected = "must use precision")]
    fn correction_precision_must_match_adele_policy() {
        let dev = LocalQp::zero(5, adele_prec(5) + 1);
        let _ = Adele::one().with_correction(5, dev);
    }

    #[test]
    fn additive_group_facts_hold_exactly() {
        let xs: Vec<Adele> = (-4i128..=4)
            .flat_map(|n| (1i128..=3).map(move |d| Adele::from_rational(&q(n, d))))
            .collect();
        let zero = Adele::zero();
        for a in &xs {
            assert_eq!(a.add(&zero), *a);
            assert_eq!(a.add(&a.neg()), zero);
            for b in &xs {
                assert_eq!(a.add(b), b.add(a)); // commutative
            }
        }
    }
}
