//! Conformal (CGA) and projective (PGA) geometric algebra layers, generic over
//! the scalar — so the conformal model can run over the **surreals**, where a
//! point can sit at `ω`-scale and a sphere can have an *infinitesimal* radius
//! `ε`, exactly (impossible with floating point).
//!
//! ## CGA
//!
//! The conformal model of Euclidean `ℝⁿ` lives in `Cl(n+1, 1)`. We build it with
//! a **null basis**: the `n` Euclidean generators `e_0..e_{n−1}` (`q = +1`) plus
//! two null vectors `n_o` (origin) and `n_∞` (infinity), `n_o² = n_∞² = 0` with
//! `{n_o, n_∞} = −2` so `n_o · n_∞ = −1`. A Euclidean point `p` lifts to the null
//! vector `up(p) = n_o + p + ½|p|² n_∞`, and `up(p) · up(q) = −½|p − q|²` turns
//! the conformal inner product into Euclidean distance.
//!
//! Because `½` is needed throughout, **CGA is a characteristic-0 feature** (in
//! char 2 the null pair degenerates, `−2 = 0`); the constructor asserts this.
//!
//! ## PGA and the exact (nilpotent) motor exponential
//!
//! Projective GA uses the degenerate metric `Cl(n,0,1)` (one null `e_0`). The
//! exponential of a **nilpotent** bivector — the translational part of a motor —
//! is a *terminating* polynomial (`B² = 0 ⇒ exp B = 1 + B`), so it is exact over
//! any backend with no transcendentals. The general Euclidean motor carries a
//! rotation (`B² < 0`) and needs `cos`/`sin`, which is out of scope:
//! `exp_nilpotent` returns `None` if the series does not terminate.

use crate::clifford::{CliffordAlgebra, Metric, Multivector};
use crate::scalar::Scalar;
use std::collections::BTreeMap;

/// `k` as a scalar (`1 + 1 + … + 1`).
fn s_int<S: Scalar>(k: usize) -> S {
    let one = S::one();
    let mut acc = S::zero();
    for _ in 0..k {
        acc = acc.add(&one);
    }
    acc
}

/// The conformal geometric algebra of Euclidean `ℝⁿ`: `Cl(n+1, 1)` in a null
/// basis, with helpers for the conformal embedding and round/flat primitives.
pub struct Cga<S: Scalar> {
    pub alg: CliffordAlgebra<S>,
    pub n: usize,
    /// generator index of `n_o` (origin).
    pub no: usize,
    /// generator index of `n_∞` (infinity).
    pub ninf: usize,
}

impl<S: Scalar> Cga<S> {
    /// Build the CGA of `ℝⁿ`. Panics unless the backend has characteristic 0
    /// and `2` is invertible (CGA needs `½`).
    pub fn new(n: usize) -> Self {
        assert_eq!(
            S::characteristic(),
            0,
            "CGA is a characteristic-0 Euclidean construction"
        );
        let two = S::one().add(&S::one());
        assert!(
            two.inv().is_some(),
            "CGA needs 1/2, so 2 must be invertible in the scalar backend"
        );
        let mut q = vec![S::one(); n];
        q.push(S::zero()); // n_o
        q.push(S::zero()); // n_∞
        let mut b = BTreeMap::new();
        // {n_o, n_∞} = −2  ⇒  n_o · n_∞ = −1
        b.insert((n, n + 1), S::one().add(&S::one()).neg());
        let alg = CliffordAlgebra::new(n + 2, Metric::new(q, b));
        Cga {
            alg,
            n,
            no: n,
            ninf: n + 1,
        }
    }

    fn half(&self) -> S {
        S::one()
            .add(&S::one())
            .inv()
            .expect("½ exists in characteristic 0")
    }

    pub fn n_o(&self) -> Multivector<S> {
        self.alg.gen(self.no)
    }
    pub fn n_inf(&self) -> Multivector<S> {
        self.alg.gen(self.ninf)
    }

    /// The conformal (symmetric) inner product `x · y = ½⟨xy + yx⟩₀`. Note the
    /// symmetrization is essential: the engine carries the polar form `b` in the
    /// *anticommutator*, so `⟨xy⟩₀` alone is the asymmetric in-order contraction,
    /// not the bilinear form.
    pub fn inner(&self, x: &Multivector<S>, y: &Multivector<S>) -> S {
        let xy = self.alg.mul(x, y);
        let yx = self.alg.mul(y, x);
        self.half()
            .mul(&self.alg.scalar_part(&self.alg.add(&xy, &yx)))
    }

    /// Lift a Euclidean point `p ∈ ℝⁿ` to the null vector `n_o + p + ½|p|² n_∞`.
    pub fn up(&self, p: &[S]) -> Multivector<S> {
        assert_eq!(p.len(), self.n, "point dimension mismatch");
        let mut acc = self.n_o();
        let mut s = S::zero();
        for (i, pi) in p.iter().enumerate() {
            acc = self
                .alg
                .add(&acc, &self.alg.scalar_mul(pi, &self.alg.gen(i)));
            s = s.add(&pi.mul(pi));
        }
        let coeff = self.half().mul(&s);
        self.alg
            .add(&acc, &self.alg.scalar_mul(&coeff, &self.n_inf()))
    }

    /// Recover a Euclidean point from a (possibly unnormalized) null vector.
    /// `None` if it cannot be normalized (e.g. a direction / point at infinity).
    pub fn down(&self, x: &Multivector<S>) -> Option<Vec<S>> {
        let f = self.inner(x, &self.n_inf()); // = −1 for a normalized point
        let factor = f.neg();
        let inv = factor.inv()?;
        let norm = self.alg.scalar_mul(&inv, x);
        Some(
            (0..self.n)
                .map(|i| {
                    norm.terms
                        .get(&(1u128 << i))
                        .cloned()
                        .unwrap_or_else(S::zero)
                })
                .collect(),
        )
    }

    /// The sphere of squared radius `r2` about center `c`: `up(c) − ½ r² n_∞`.
    pub fn sphere(&self, c: &[S], r2: &S) -> Multivector<S> {
        let coeff = self.half().mul(r2).neg();
        self.alg
            .add(&self.up(c), &self.alg.scalar_mul(&coeff, &self.n_inf()))
    }

    /// The plane `{x : x·normal = d}`: `Σ normalᵢ eᵢ + d n_∞`.
    pub fn plane(&self, normal: &[S], d: &S) -> Multivector<S> {
        assert_eq!(normal.len(), self.n, "normal dimension mismatch");
        let mut acc = self.alg.scalar_mul(d, &self.n_inf());
        for (i, ni) in normal.iter().enumerate() {
            acc = self
                .alg
                .add(&acc, &self.alg.scalar_mul(ni, &self.alg.gen(i)));
        }
        acc
    }

    /// The point pair / oriented join `a ∧ b`.
    pub fn point_pair(&self, a: &Multivector<S>, b: &Multivector<S>) -> Multivector<S> {
        self.alg.wedge(a, b)
    }

    /// The meet (intersection) of two IPNS objects — the outer product `x ∧ y`.
    /// In the inner-product-null-space convention used here (a point is *on* `X`
    /// iff `up(p) · X = 0`), intersection is the wedge; this needs no pseudoscalar
    /// inverse, so it works over every char-0 backend including the surreals.
    pub fn meet(&self, x: &Multivector<S>, y: &Multivector<S>) -> Multivector<S> {
        self.alg.wedge(x, y)
    }
}

/// The projective geometric algebra `Cl(n,0,1)` of `ℝⁿ`: one degenerate
/// generator `e_0` (`q = 0`) and `n` Euclidean ones.
pub fn pga<S: Scalar>(n: usize) -> CliffordAlgebra<S> {
    let mut q = vec![S::zero()]; // e_0 null (the ideal/projective direction)
    q.extend(std::iter::repeat_n(S::one(), n));
    CliffordAlgebra::new(n + 1, Metric::diagonal(q))
}

/// `exp(B)` for a **nilpotent** multivector, as the terminating power series
/// `Σ Bᵏ/k!`. Returns `None` if the series does not terminate within this
/// implementation's conservative bound, or if some `k!` is not invertible in the
/// backend. For the intended PGA translator/motor nilpotents the bound is ample;
/// in large Clifford algebras, `None` is only a refusal, not proof that `B` is
/// non-nilpotent.
pub fn exp_nilpotent<S: Scalar>(
    alg: &CliffordAlgebra<S>,
    b: &Multivector<S>,
) -> Option<Multivector<S>> {
    let cap = 2 * alg.dim() + 2;
    let mut acc = alg.scalar(S::one());
    let mut power = alg.scalar(S::one()); // B^0
    let mut fact = S::one(); // 0!
    for k in 1..=cap {
        power = alg.mul(&power, b);
        if power.is_zero() {
            return Some(acc); // nilpotent: series terminated
        }
        fact = fact.mul(&s_int::<S>(k));
        let finv = fact.inv()?;
        acc = alg.add(&acc, &alg.scalar_mul(&finv, &power));
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scalar::Fp;
    use crate::scalar::Integer;
    use crate::scalar::Rational;
    use crate::scalar::Surreal;

    fn r(n: i128) -> Rational {
        Rational::int(n)
    }
    fn rs(num: i128, den: i128) -> Rational {
        Rational::new(num, den)
    }

    #[test]
    fn cga_rejects_rings_without_one_half() {
        assert!(std::panic::catch_unwind(|| Cga::<Integer>::new(1)).is_err());
    }

    #[test]
    fn cga_rejects_positive_characteristic() {
        assert!(std::panic::catch_unwind(|| Cga::<Fp<3>>::new(1)).is_err());
    }

    #[test]
    fn up_is_null() {
        let cga = Cga::<Rational>::new(2);
        for p in [[r(3), r(4)], [r(0), r(0)], [r(-2), r(5)]] {
            assert_eq!(cga.inner(&cga.up(&p), &cga.up(&p)), r(0));
        }
    }

    #[test]
    fn inner_product_is_euclidean_distance() {
        let cga = Cga::<Rational>::new(2);
        let p = [r(1), r(0)];
        let q = [r(4), r(4)]; // |p−q|² = 9 + 16 = 25
        assert_eq!(cga.inner(&cga.up(&p), &cga.up(&q)), rs(-25, 2));
    }

    #[test]
    fn down_inverts_up() {
        let cga = Cga::<Rational>::new(3);
        let p = vec![r(2), r(-3), r(5)];
        assert_eq!(cga.down(&cga.up(&p)).unwrap(), p);
        // also for a rescaled representative
        let scaled = cga.alg.scalar_mul(&r(7), &cga.up(&p));
        assert_eq!(cga.down(&scaled).unwrap(), p);
    }

    #[test]
    fn point_lies_on_sphere() {
        let cga = Cga::<Rational>::new(2);
        let c = [r(0), r(0)];
        let s = cga.sphere(&c, &r(25)); // radius² = 25
                                        // p at distance 5 from the origin is on the sphere
        assert_eq!(cga.inner(&cga.up(&[r(3), r(4)]), &s), r(0));
        // a point inside is not
        assert_ne!(cga.inner(&cga.up(&[r(1), r(1)]), &s), r(0));
    }

    #[test]
    fn point_lies_on_plane() {
        let cga = Cga::<Rational>::new(2);
        // plane x·(1,0) = 3, i.e. the line x=3
        let pl = cga.plane(&[r(1), r(0)], &r(3));
        assert_eq!(cga.inner(&cga.up(&[r(3), r(9)]), &pl), r(0));
        assert_ne!(cga.inner(&cga.up(&[r(2), r(9)]), &pl), r(0));
    }

    #[test]
    fn meet_of_planes_is_nonzero() {
        let cga = Cga::<Rational>::new(3);
        let a = cga.plane(&[r(1), r(0), r(0)], &r(0));
        let b = cga.plane(&[r(0), r(1), r(0)], &r(0));
        let m = cga.meet(&a, &b);
        assert!(!m.is_zero());
        assert_eq!(cga.alg.grade_part(&m, 2), m); // a line is a grade-2 IPNS blade
    }

    #[test]
    fn surreal_point_at_infinite_scale_is_still_null() {
        // A point with an ω-scale coordinate: the conformal embedding still
        // lands on the null cone, exactly.
        let cga = Cga::<Surreal>::new(2);
        let p = [Surreal::omega(), Surreal::zero()];
        assert_eq!(cga.inner(&cga.up(&p), &cga.up(&p)), Surreal::zero());
    }

    #[test]
    fn surreal_sphere_of_infinitesimal_radius() {
        // A sphere of radius ε: a point at distance ε is exactly on it, a point
        // at distance 2ε is exactly off it — distinguishing infinitesimal radii.
        let cga = Cga::<Surreal>::new(2);
        let eps = Surreal::epsilon();
        let eps2 = eps.mul(&eps);
        let s = cga.sphere(&[Surreal::zero(), Surreal::zero()], &eps2);
        let on = cga.up(&[eps.clone(), Surreal::zero()]);
        assert_eq!(cga.inner(&on, &s), Surreal::zero());
        let off = cga.up(&[eps.mul(&Surreal::from_int(2)), Surreal::zero()]);
        assert_ne!(cga.inner(&off, &s), Surreal::zero());
    }

    #[test]
    fn pga_nilpotent_exp_is_exact() {
        // In Cl(2,0,1), B = e0∧e1 is nilpotent (e0²=0), so exp(B) = 1 + B exactly.
        let alg = pga::<Rational>(2);
        let (e0, e1) = (alg.gen(0), alg.gen(1));
        let b = alg.wedge(&e0, &e1);
        assert!(alg.mul(&b, &b).is_zero());
        assert_eq!(
            exp_nilpotent(&alg, &b).unwrap(),
            alg.add(&alg.scalar(r(1)), &b)
        );
        // scaling: exp(3B) = 1 + 3B
        let b3 = alg.scalar_mul(&r(3), &b);
        assert_eq!(
            exp_nilpotent(&alg, &b3).unwrap(),
            alg.add(&alg.scalar(r(1)), &b3)
        );
    }

    #[test]
    fn pga_motor_translates_exactly() {
        // The motor M = 1 + B (B = e0e1) is a versor; its sandwich translates
        // e1 ↦ e1 + 2 e0 exactly (a translation along the ideal direction).
        let alg = pga::<Rational>(2);
        let (e0, e1) = (alg.gen(0), alg.gen(1));
        let b = alg.wedge(&e0, &e1);
        let motor = exp_nilpotent(&alg, &b).unwrap(); // 1 + B
        let moved = alg.sandwich(&motor, &e1).unwrap();
        let expect = alg.add(&e1, &alg.scalar_mul(&r(2), &e0));
        assert_eq!(moved, expect);
    }

    #[test]
    fn non_nilpotent_exp_returns_none() {
        // A Euclidean rotation bivector squares to −1 (not nilpotent) ⇒ the
        // series never terminates ⇒ None (would need transcendental cos/sin).
        let alg = CliffordAlgebra::new(2, Metric::diagonal(vec![r(1), r(1)]));
        let b = alg.wedge(&alg.gen(0), &alg.gen(1));
        assert!(alg.mul(&b, &b) == alg.scalar(r(-1)));
        assert!(exp_nilpotent(&alg, &b).is_none());
    }
}
