//! The Witt **ring**: the multiplicative structure on top of the additive Witt
//! group, the fundamental ideal `Iⁿ`, Pfister forms, and the cohomological
//! invariant staircase `eₙ`.
//!
//! `witt.rs` carries the additive group `W` (forms mod hyperbolics, under `⊥`).
//! Tensor product of forms makes `W` a **ring**, and the powers of the
//! **fundamental ideal** `I = ker(e₀)` (the even-dimensional classes) filter it:
//! `W ⊇ I ⊇ I² ⊇ …`. The Milnor conjecture (Voevodsky) identifies
//! `Iⁿ/Iⁿ⁺¹ ≅ Hⁿ(F, ℤ/2)`, and the **n-fold Pfister forms**
//! `⟨⟨a₁,…,aₙ⟩⟩ = ⊗ᵢ⟨1, −aᵢ⟩` generate `Iⁿ`. The first three quotient maps are
//! invariants the rest of this crate already computes:
//!
//! | n | `eₙ` | reused from |
//! |---|------|-------------|
//! | 0 | `dim mod 2`            | here (trivial) |
//! | 1 | signed **discriminant**| [`oddchar::finite_odd_witt`] (the `sclass`) |
//! | 2 | **Hasse**/Clifford     | [`oddchar::hasse_invariant_finite_odd`] |
//!
//! So this module is the *retro-unification*: discriminant and Hasse stop being
//! separate functions and become `e₁`, `e₂` — successive steps of one staircase.
//!
//! ## Stabilization — where the staircase stops, per field
//!
//! The tower's length is the field's `u`-invariant story:
//!
//!   * **Finite `F_q`** (`u = 2`): `I² = 0`. Every 2-fold Pfister form is
//!     hyperbolic (a 4-dim form over `F_q` is isotropic, and an isotropic Pfister
//!     form is hyperbolic). So the staircase is just `(e₀, e₁)`, and
//!     [`two_fold_pfister_is_hyperbolic`](self) demonstrates `I² = 0` directly.
//!     `e₂` is identically trivial — the same fact `oddchar`'s Hasse `≡ +1` records.
//!   * **`Q_p`** (`u = 4`): `I³ = 0`. The staircase reaches `(e₀, e₁, e₂)` with
//!     **`e₂` = Hasse genuinely nontrivial** — the payoff that `forms/padic.rs`
//!     supplies (over a finite field the Brauer group is trivial so `e₂` carries
//!     nothing; over `Q_p` it does).
//!   * **`ℝ`** (surreal backend, `u = ∞`): the tower is **infinite**. `Iⁿ = 2ⁿℤ`
//!     via the signature, and `eₙ` reads the 2-adic expansion of the signature —
//!     see [`e_real`].
//!
//! ## The characteristic-2 caveat (pinned, not asserted)
//!
//! In characteristic 2 the staircase does **not** index-match the above. The Witt
//! group of *quadratic* forms `W_q(F)` is a **module over** the Witt ring of
//! *bilinear* forms, not a ring; its filtration is Kato's (built from differential
//! forms `Ωⁿ`), not the Milnor `Iⁿ`. The **Arf invariant is the leading
//! invariant** of `W_q` — and for a nonsingular char-2 form `dim` is always even so
//! `e₀ ≡ 0` — but its cohomological *degree* is a convention we deliberately do not
//! force onto the char-0/odd indexing. So: Arf is exposed as *the* char-2 invariant
//! ([`char2::arf_invariant`]); we do **not** claim "Arf = e₂".

use crate::clifford::Metric;
use crate::forms::{finite_odd_witt, hasse_invariant_finite_odd, FiniteOddField, WittClassG};
use crate::scalar::Scalar;

// ---------------------------------------------------------------------------
// The ring multiplication on representatives: the tensor product of forms.
// ---------------------------------------------------------------------------

/// The tensor product `φ ⊗ ψ` of two **diagonal** forms: `⟨a_i⟩ ⊗ ⟨b_j⟩ = ⟨a_i b_j⟩`.
/// This is the Witt-ring multiplication realised on representatives (it descends to
/// classes because tensoring with a hyperbolic plane stays hyperbolic). `None` if
/// either metric is non-diagonal. Generic over the scalar field, so it runs over
/// `Fp`, `Fpn`, `Surreal`, … alike.
pub fn tensor_form<S: Scalar>(a: &Metric<S>, b: &Metric<S>) -> Option<Metric<S>> {
    if !a.b.is_empty() || !a.a.is_empty() || !b.b.is_empty() || !b.a.is_empty() {
        return None;
    }
    let mut q = Vec::with_capacity(a.q.len() * b.q.len());
    for ai in &a.q {
        for bj in &b.q {
            q.push(ai.mul(bj));
        }
    }
    Some(Metric::diagonal(q))
}

/// The 1-fold Pfister form `⟨⟨a⟩⟩ = ⟨1, −a⟩`.
pub fn pfister1<S: Scalar>(a: &S) -> Metric<S> {
    Metric::diagonal(vec![S::one(), a.neg()])
}

/// The n-fold Pfister form `⟨⟨a₁,…,aₙ⟩⟩ = ⊗ᵢ ⟨1, −aᵢ⟩`, the generator of `Iⁿ`.
/// The empty product is the rank-1 unit form `⟨1⟩` (the ring identity).
pub fn pfister<S: Scalar>(scales: &[S]) -> Metric<S> {
    let mut acc = Metric::diagonal(vec![S::one()]);
    for a in scales {
        // tensor_form on two diagonal metrics never returns None.
        acc = tensor_form(&acc, &pfister1(a)).expect("Pfister factors are diagonal");
    }
    acc
}

/// Membership in the **fundamental ideal** `I = ker(e₀)`: a (diagonal) form is in
/// `I` iff its nondegenerate dimension is even.
pub fn in_fundamental_ideal<S: Scalar>(metric: &Metric<S>) -> bool {
    metric.q.iter().filter(|x| !x.is_zero()).count() % 2 == 0
}

// ---------------------------------------------------------------------------
// The cohomological invariant staircase.
// ---------------------------------------------------------------------------

/// The low cohomological invariants `(e₀, e₁, e₂)` of an odd-characteristic form,
/// with the field's stabilization recorded. `e₀ = dim mod 2`, `e₁ =` signed-disc
/// square-class (the genuine `H¹` invariant, reused from [`finite_odd_witt`]), and
/// `e₂ =` the Hasse invariant — `+1` over a finite field, where `I² = 0`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EnStaircase {
    /// `e₀ ∈ H⁰ = ℤ/2`: the dimension mod 2.
    pub e0: u8,
    /// `e₁ ∈ H¹ = F*/F*²`: the signed-discriminant square-class (0 square / 1 not).
    pub e1: u8,
    /// `e₂ ∈ H² = Br(F)[2]`: the Hasse/Clifford invariant. `+1` over a finite field.
    pub e2: i8,
    /// `Iⁿ = 0` for `n ≥ this`. For a finite field, `2`.
    pub stabilizes_at: usize,
}

/// The `(e₀, e₁, e₂)` staircase of a diagonal form over any finite field of odd
/// characteristic. `None` if non-diagonal. Over a finite field `I² = 0`, so
/// `stabilizes_at = 2` and `e₂` is always `+1`; the genuine content is `(e₀, e₁)`.
pub fn e_staircase_finite_odd<F: FiniteOddField>(metric: &Metric<F>) -> Option<EnStaircase> {
    let (e0, e1) = match finite_odd_witt(metric)? {
        WittClassG::OddChar { e0, sclass, .. } => (e0, sclass),
        _ => unreachable!("finite_odd_witt returns the OddChar variant"),
    };
    Some(EnStaircase {
        e0,
        e1,
        e2: hasse_invariant_finite_odd(metric)?,
        stabilizes_at: 2,
    })
}

/// The real cohomological invariant `eₙ` of a form of signature `σ` over `ℝ`
/// (the surreal backend's leg). Over `ℝ`, `W(ℝ) ≅ ℤ` via the signature and
/// `Iⁿ = 2ⁿℤ`, so the form is in `Iⁿ` iff `2ⁿ | σ`, and then
/// `eₙ = (σ / 2ⁿ) mod 2` — the staircase reads the 2-adic expansion of the
/// signature. Returns `None` when the form is **not** in `Iⁿ` (so `eₙ` is undefined).
/// This is the infinite tower the finite-field and `Q_p` legs truncate.
pub fn e_real(signature: i128, n: usize) -> Option<u8> {
    if n >= 128 {
        return None;
    }
    let modulus = 1i128 << n;
    if signature % modulus != 0 {
        return None;
    }
    Some((signature / modulus).rem_euclid(2) as u8)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scalar::{Fp, Fpn};

    fn diag<const P: u128>(qs: &[u128]) -> Metric<Fp<P>> {
        Metric::diagonal(qs.iter().map(|&x| Fp::<P>(x)).collect())
    }

    #[test]
    fn tensor_form_multiplies_entries() {
        // ⟨1,2⟩ ⊗ ⟨1,2⟩ = ⟨1,2,2,4⟩ over F_5.
        let t = tensor_form(&diag::<5>(&[1, 2]), &diag::<5>(&[1, 2])).unwrap();
        assert_eq!(t.q, vec![Fp::<5>(1), Fp::<5>(2), Fp::<5>(2), Fp::<5>(4)]);
    }

    #[test]
    fn pfister_shapes() {
        // ⟨⟨a⟩⟩ = ⟨1,−a⟩; ⟨⟨a,b⟩⟩ = ⟨1,−a,−b,ab⟩.
        let p1 = pfister(&[Fp::<7>(3)]);
        assert_eq!(p1.q, vec![Fp::<7>::one(), Fp::<7>::new(-3)]);
        let p2 = pfister(&[Fp::<7>(3), Fp::<7>(5)]);
        // ⟨1, −3⟩ ⊗ ⟨1, −5⟩ = ⟨1, −5, −3, 15⟩
        assert_eq!(
            p2.q,
            vec![
                Fp::<7>::one(),
                Fp::<7>::new(-5),
                Fp::<7>::new(-3),
                Fp::<7>(15 % 7),
            ]
        );
    }

    #[test]
    fn two_fold_pfister_is_hyperbolic() {
        // THE I²=0 ORACLE: every 2-fold Pfister form over a finite field is
        // Witt-trivial (a 4-dim form over F_q is isotropic ⇒ the Pfister form is
        // hyperbolic). So I² = 0 and the odd-char staircase stops at (e₀, e₁).
        for a in 1..5u128 {
            for b in 1..5u128 {
                let p = pfister(&[Fp::<5>(a), Fp::<5>(b)]);
                assert_eq!(
                    finite_odd_witt(&p).unwrap(),
                    WittClassG::oddchar_zero(0), // F_5: −1 square, identity has kappa 0
                    "2-fold Pfister ⟨⟨{a},{b}⟩⟩ over F_5 must be hyperbolic"
                );
            }
        }
        for a in 1..3u128 {
            for b in 1..3u128 {
                let p = pfister(&[Fp::<3>(a), Fp::<3>(b)]);
                assert_eq!(finite_odd_witt(&p).unwrap(), WittClassG::oddchar_zero(1));
            }
        }
    }

    #[test]
    fn i_squared_is_zero_over_extension_field_f9() {
        // The same I²=0 fact over a genuine EXTENSION field F_9 = F_3², exercising
        // the Fpn backend through the invariant theory: a 2-fold Pfister form has
        // square signed-discriminant (so it is Witt-trivial). Generic tensor_form +
        // pfister + Fpn::is_square — no Fp-specific code.
        let elems: Vec<Fpn<3, 2>> = (1..9u128)
            .map(|code| {
                let mut c = code;
                let mut coeffs = [0u128; 2];
                for slot in coeffs.iter_mut() {
                    *slot = (c % 3) as u128;
                    c /= 3;
                }
                Fpn(coeffs)
            })
            .filter(|x| !x.is_zero())
            .collect();
        for &a in &elems {
            for &b in &elems {
                let p = pfister(&[a, b]);
                // det = 1·(−a)·(−b)·(ab) = (ab)², always a square ⇒ signed disc square.
                let det = p.q.iter().fold(Fpn::<3, 2>::one(), |acc, x| acc.mul(x));
                assert!(
                    det.is_square(),
                    "2-fold Pfister over F_9 disc must be square"
                );
            }
        }
    }

    #[test]
    fn e_staircase_reuses_disc_and_hasse() {
        // e₁ is exactly oddchar's signed-disc class; e₂ is exactly the Hasse
        // invariant (+1 over a finite field).
        let m = diag::<5>(&[1, 2, 3]);
        let s = e_staircase_finite_odd(&m).unwrap();
        assert_eq!(s.e0, 1); // dim 3
        assert_eq!(s.e2, hasse_invariant_finite_odd(&m).unwrap());
        assert_eq!(s.e2, 1); // trivial over a finite field
        assert_eq!(s.stabilizes_at, 2);
        // e₁ matches the Witt class's sclass.
        if let WittClassG::OddChar { sclass, .. } = finite_odd_witt(&m).unwrap() {
            assert_eq!(s.e1, sclass);
        }
    }

    #[test]
    fn class_ring_mul_matches_concrete_tensor_form() {
        // THE RING-LAW ORACLE: the derived class-level product WittClassG::mul must
        // agree with classifying the actual tensor product of forms, for every pair
        // of small nondegenerate forms. Verifies both the ℤ/4 (F_3) and F₂[ℤ/2] (F_5)
        // ring laws against ground truth, and that they are well-defined on classes.
        fn check<const P: u128>() {
            // all 1- and 2-dimensional nondegenerate diagonal forms
            let mut forms: Vec<Metric<Fp<P>>> = Vec::new();
            for e in 1..P {
                forms.push(diag::<P>(&[e]));
            }
            for e in 1..P {
                for f in 1..P {
                    forms.push(diag::<P>(&[e, f]));
                }
            }
            for a in &forms {
                for b in &forms {
                    let lhs = finite_odd_witt(&tensor_form(a, b).unwrap()).unwrap();
                    let rhs = finite_odd_witt(a)
                        .unwrap()
                        .mul(&finite_odd_witt(b).unwrap());
                    assert_eq!(lhs, rhs, "P={P}: ring law disagrees with tensor_form");
                }
            }
        }
        check::<3>(); // kappa = 1, W ≅ ℤ/4
        check::<5>(); // kappa = 0, W ≅ F₂[ℤ/2]
    }

    #[test]
    fn ring_unit_is_neutral() {
        // ⟨1⟩ is the multiplicative identity in both odd-char ring structures.
        let one3 = WittClassG::oddchar_one(1);
        let one5 = WittClassG::oddchar_one(0);
        for m in [diag::<3>(&[1]), diag::<3>(&[2]), diag::<3>(&[1, 2])] {
            let c = finite_odd_witt(&m).unwrap();
            assert_eq!(c.mul(&one3), c);
        }
        for m in [diag::<5>(&[1]), diag::<5>(&[2]), diag::<5>(&[1, 2])] {
            let c = finite_odd_witt(&m).unwrap();
            assert_eq!(c.mul(&one5), c);
        }
        // and Char0: signatures multiply, ⟨1⟩ is unit 1.
        let sig3 = WittClassG::char0(3, 0);
        assert_eq!(sig3.mul(&WittClassG::char0(1, 0)), sig3);
        assert_eq!(
            WittClassG::char0(2, 0).mul(&WittClassG::char0(5, 0)),
            WittClassG::Char0 { signature: 10 }
        );
    }

    #[test]
    fn real_staircase_reads_the_signature_in_binary() {
        // ⟨1,1,1,1⟩ over ℝ has signature 4 = 0b100: in I⁰, I¹, I² but not I³.
        assert_eq!(e_real(4, 0), Some(0)); // 4 even ⇒ in I, e₀ = 0
        assert_eq!(e_real(4, 1), Some(0)); // 4/2 = 2 even ⇒ in I², e₁ = 0
        assert_eq!(e_real(4, 2), Some(1)); // 4/4 = 1 odd  ⇒ e₂ = 1, the leading bit
        assert_eq!(e_real(4, 3), None); //   8 ∤ 4 ⇒ not in I³
                                        // signature 6 = 0b110: e₁ and e₂ both fire (in I² but 6/4 non-integer).
        assert_eq!(e_real(6, 1), Some(1)); // 6/2 = 3 odd
        assert_eq!(e_real(6, 2), None); // 4 ∤ 6
                                        // negative signatures (indefinite-ish) read the same way.
        assert_eq!(e_real(-8, 3), Some(1)); // −8/8 = −1 ≡ 1 mod 2
    }
}
