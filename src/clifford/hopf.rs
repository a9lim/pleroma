//! The exterior (Grassmann) Hopf algebra structure on the multivector engine:
//! coproduct, counit, antipode, with the Hopf axioms verified over both
//! characteristics.
//!
//! On the exterior algebra `Λ(V)` (the `grassmann` metric) the generators are
//! **primitive** — `Δ(e_i) = e_i ⊗ 1 + 1 ⊗ e_i` — and `Δ` extends as an algebra
//! map for the wedge product. In closed form this is the **unshuffle
//! coproduct** on blades:
//!
//!   `Δ(e_S) = Σ_{T ⊆ S} sign(T, S∖T) · (e_T ⊗ e_{S∖T})`,
//!
//! where `sign(T, S∖T)` is exactly the reordering sign of `e_T ∧ e_{S∖T} = e_S`.
//! We read that sign straight off `alg.wedge(e_T, e_{S∖T})` — so it is never
//! hand-written, and in characteristic 2 (`−1 = 1`) all signs collapse to `+`
//! automatically.
//!
//! ## The antipode is the grade involution
//!
//! For this (primitively generated) coproduct the antipode is `S(x) = (−1)^k x`
//! on `Λ^k` — the **grade involution**, *not* `reverse ∘ grade_involution`.
//! Check on `v ∧ w`: `m(S⊗id)Δ(vw) = S(vw) − v∧w` must equal `ε(vw)·1 = 0`, so
//! `S(vw) = +v∧w = (−1)^2 vw`; the reversion-twisted sign `(−1)^{k(k+1)/2}` would
//! give `−vw` and fail the axiom. The tests pin this down.
//!
//! A tensor element `e_T ⊗ e_U` of `Cl ⊗̂ Cl` is encoded as the blade
//! `T | (U << dim)` of `tensor_square(alg)` (the low block is the left factor,
//! the high block the right) — matching `embed_first` / `embed_second(·, dim)`.

use crate::clifford::{bits, CliffordAlgebra, Multivector};
use crate::scalar::Scalar;
use std::collections::BTreeMap;

/// The graded tensor square `Cl ⊗̂ Cl`, the codomain of the coproduct.
pub fn tensor_square<S: Scalar>(alg: &CliffordAlgebra<S>) -> CliffordAlgebra<S> {
    alg.graded_tensor(alg)
}

fn blade_of<S: Scalar>(alg: &CliffordAlgebra<S>, mask: u32) -> Multivector<S> {
    alg.blade(&bits(mask))
}

/// The unshuffle coproduct `Δ: Cl → Cl ⊗̂ Cl`, returned as a multivector over
/// `tensor_square(alg)` (a tensor `e_T ⊗ e_U` is the blade `T | (U << dim)`).
pub fn coproduct<S: Scalar>(alg: &CliffordAlgebra<S>, mv: &Multivector<S>) -> Multivector<S> {
    let dim = alg.dim;
    let mut out: BTreeMap<u32, S> = BTreeMap::new();
    for (&mask_s, coeff) in &mv.terms {
        // iterate every submask T of mask_s (including 0 and mask_s)
        let mut t = mask_s;
        loop {
            let u = mask_s ^ t;
            // sign such that e_T ∧ e_U = sign · e_S — read it off `wedge`.
            let w = alg.wedge(&blade_of(alg, t), &blade_of(alg, u));
            let sign = w.terms.get(&mask_s).cloned().unwrap_or_else(S::zero);
            if !sign.is_zero() {
                let tens = t | (u << dim);
                let e = out.entry(tens).or_insert_with(S::zero);
                *e = e.add(&coeff.mul(&sign));
                if e.is_zero() {
                    out.remove(&tens);
                }
            }
            if t == 0 {
                break;
            }
            t = (t - 1) & mask_s;
        }
    }
    Multivector { terms: out }
}

/// The counit `ε: Cl → S` — projection to the scalar part.
pub fn counit<S: Scalar>(alg: &CliffordAlgebra<S>, mv: &Multivector<S>) -> S {
    alg.scalar_part(mv)
}

/// The antipode `S: Cl → Cl` — the grade involution `(−1)^k` on `Λ^k`.
pub fn antipode<S: Scalar>(alg: &CliffordAlgebra<S>, mv: &Multivector<S>) -> Multivector<S> {
    alg.grade_involution(mv)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::clifford::{grade, Metric};
    use crate::scalar::Nimber;
    use crate::scalar::Rational;

    fn r(n: i128) -> Rational {
        Rational::int(n)
    }

    /// Split a tensor-square multivector into a (left-mask, right-mask) → coeff map.
    fn pairs<S: Scalar>(alg: &CliffordAlgebra<S>, x: &Multivector<S>) -> BTreeMap<(u32, u32), S> {
        let dim = alg.dim;
        let low = if dim >= 32 {
            u32::MAX
        } else {
            (1u32 << dim) - 1
        };
        let dtens = coproduct(alg, x);
        dtens
            .terms
            .into_iter()
            .map(|(mask, c)| ((mask & low, mask >> dim), c))
            .collect()
    }

    /// Counit law: `(ε⊗id)∘Δ = id = (id⊗ε)∘Δ`.
    fn check_counit_law<S: Scalar>(alg: &CliffordAlgebra<S>, x: &Multivector<S>) {
        let p = pairs(alg, x);
        // (ε⊗id): keep terms with empty left leg, value is the right blade.
        let mut left = alg.zero();
        let mut right = alg.zero();
        for (&(t, u), c) in &p {
            if t == 0 {
                left = alg.add(&left, &alg.scalar_mul(c, &alg.blade(&bits(u))));
            }
            if u == 0 {
                right = alg.add(&right, &alg.scalar_mul(c, &alg.blade(&bits(t))));
            }
        }
        assert_eq!(&left, x, "(ε⊗id)∘Δ ≠ id");
        assert_eq!(&right, x, "(id⊗ε)∘Δ ≠ id");
    }

    /// Coassociativity: `(Δ⊗id)∘Δ = (id⊗Δ)∘Δ` as triple-tensor maps.
    fn check_coassociativity<S: Scalar>(alg: &CliffordAlgebra<S>, x: &Multivector<S>) {
        let p = pairs(alg, x);
        let mut lhs: BTreeMap<(u32, u32, u32), S> = BTreeMap::new();
        let mut rhs: BTreeMap<(u32, u32, u32), S> = BTreeMap::new();
        for (&(t, u), c) in &p {
            // (Δ⊗id): split the left leg
            for (&(t1, t2), d) in &pairs(alg, &alg.blade(&bits(t))) {
                let key = (t1, t2, u);
                let e = lhs.entry(key).or_insert_with(S::zero);
                *e = e.add(&c.mul(d));
                if e.is_zero() {
                    lhs.remove(&key);
                }
            }
            // (id⊗Δ): split the right leg
            for (&(u1, u2), d) in &pairs(alg, &alg.blade(&bits(u))) {
                let key = (t, u1, u2);
                let e = rhs.entry(key).or_insert_with(S::zero);
                *e = e.add(&c.mul(d));
                if e.is_zero() {
                    rhs.remove(&key);
                }
            }
        }
        assert_eq!(lhs, rhs, "coproduct is not coassociative");
    }

    /// Antipode axiom: `m∘(S⊗id)∘Δ = η∘ε`.
    fn check_antipode_axiom<S: Scalar>(alg: &CliffordAlgebra<S>, x: &Multivector<S>) {
        let p = pairs(alg, x);
        let mut acc = alg.zero();
        for (&(t, u), c) in &p {
            let st = antipode(alg, &alg.blade(&bits(t)));
            let term = alg.mul(&st, &alg.blade(&bits(u)));
            acc = alg.add(&acc, &alg.scalar_mul(c, &term));
        }
        let expect = alg.scalar(counit(alg, x));
        assert_eq!(acc, expect, "antipode axiom failed");
    }

    fn run_axioms<S: Scalar>(alg: &CliffordAlgebra<S>, elts: &[Multivector<S>]) {
        for x in elts {
            check_counit_law(alg, x);
            check_coassociativity(alg, x);
            check_antipode_axiom(alg, x);
        }
    }

    #[test]
    fn hopf_axioms_grassmann_rational() {
        let alg = CliffordAlgebra::new(3, Metric::<Rational>::grassmann(3));
        let elts = [
            alg.scalar(r(1)),
            alg.gen(0),
            alg.gen(1),
            alg.wedge(&alg.gen(0), &alg.gen(1)),
            alg.wedge(&alg.wedge(&alg.gen(0), &alg.gen(1)), &alg.gen(2)),
            alg.add(&alg.gen(0), &alg.wedge(&alg.gen(1), &alg.gen(2))),
        ];
        run_axioms(&alg, &elts);
    }

    #[test]
    fn hopf_axioms_grassmann_nimber() {
        // char 2: every sign is +, antipode = identity — the axioms still hold.
        let alg = CliffordAlgebra::new(3, Metric::<Nimber>::grassmann(3));
        let elts = [
            alg.scalar(Nimber(1)),
            alg.gen(0),
            alg.wedge(&alg.gen(0), &alg.gen(1)),
            alg.wedge(&alg.wedge(&alg.gen(0), &alg.gen(1)), &alg.gen(2)),
        ];
        run_axioms(&alg, &elts);
    }

    #[test]
    fn antipode_is_grade_involution_not_reversion_twist() {
        let alg = CliffordAlgebra::new(3, Metric::<Rational>::grassmann(3));
        for mask in 0u32..8 {
            let blade = alg.blade(&bits(mask));
            let k = grade(mask);
            let expect = if k & 1 == 1 {
                alg.scalar_mul(&r(-1), &blade)
            } else {
                blade.clone()
            };
            assert_eq!(antipode(&alg, &blade), expect, "mask {mask:#b}");
            // the reversion-twist (−1)^{k(k+1)/2} differs at k=2 (gives −1): a
            // grade-2 blade's antipode is +blade, proving we use the right sign.
            if k == 2 {
                assert_eq!(antipode(&alg, &blade), blade);
            }
        }
    }

    #[test]
    fn antipode_is_identity_over_nimber() {
        let alg = CliffordAlgebra::new(3, Metric::<Nimber>::grassmann(3));
        for mask in 0u32..8 {
            let blade = alg.blade(&bits(mask));
            assert_eq!(antipode(&alg, &blade), blade);
        }
    }
}
