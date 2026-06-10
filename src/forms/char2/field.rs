//! Char-2 finite-field capability trait — the **additive** mirror of
//! [`FiniteOddField`](crate::forms::FiniteOddField).
//!
//! In odd characteristic the form-theoretic local datum is the **multiplicative**
//! square class `F*/(F*)²`, read by Euler's criterion. In characteristic 2 that
//! datum is *trivial*: the Frobenius `x ↦ x²` is a bijection of a finite field, so
//! **every** element is a square (`F*/(F*)² = 1`). The invariant that does the work
//! is instead **additive** — the Artin–Schreier class `F/℘(F) ≅ F₂`, where
//! `℘(x) = x² + x`, read by the absolute trace `Tr_{F_q/F₂}` (the kernel of which is
//! exactly `℘(F_q)`). This is the same trace that the [`char2`](crate::forms::char2)
//! Arf reduction already pushes through, and the same obstruction that decides
//! Artin–Schreier solvability — one map, both roles.
//!
//! So `FiniteChar2Field` carries [`artin_schreier_class`](FiniteChar2Field::artin_schreier_class)
//! where [`FiniteOddField`](crate::forms::FiniteOddField) carries `is_square_value`.
//! It is the capability the (future) char-2 function-field layer over `F_{2^m}(t)`
//! is generic over, exactly as the odd layer is generic over `FiniteOddField`.
//!
//! **Scope (honest):** implemented for the prime field [`Fp<2>`](crate::scalar::Fp)
//! and its extensions [`Fpn<2, N>`](crate::scalar::Fpn) — the practical coefficient
//! fields of a char-2 function field `F_{2^m}(t)`. [`Nimber`](crate::scalar::Nimber)
//! is deliberately **excluded**, the same boundary `FiniteOddField` draws: it is the
//! *direct* Arf backend (`Metric<Nimber>`), never the coefficient field of a function
//! field, and its order `2^128` does not fit the `u128` field-order metadata the
//! place layer needs. Its absolute trace is still available through
//! [`FieldExtension::trace`](crate::scalar::FieldExtension).

use crate::scalar::{ExactFieldScalar, Fp, Fpn};

/// Finite fields of characteristic 2, with the operations char-2 form theory needs:
/// field-order metadata, an enumeration, and the **Artin–Schreier class** (the
/// additive analogue of the odd-characteristic square class). Intentionally narrower
/// than [`Scalar`] — a form-theory façade, not a new scalar-world requirement, the
/// mirror of [`FiniteOddField`](crate::forms::FiniteOddField).
pub trait FiniteChar2Field: ExactFieldScalar + Copy {
    /// Characteristic prime — always `2` (provided; the trait is char-2 only).
    fn characteristic_prime() -> u128 {
        2
    }

    /// Field order `q = 2^m`.
    fn field_order() -> u128;

    /// Whether this type is a supported finite field of characteristic 2.
    fn is_supported_char2_field() -> bool;

    /// Embed an ordinary integer through the prime subfield `F₂` (so `n ↦ n mod 2`).
    fn from_i128(n: i128) -> Self;

    /// Enumerate the field: index `i ∈ [0, field_order())` ↦ a distinct element,
    /// covering all of `F_q` exactly once (base-2 digits of `i` are the
    /// polynomial-basis coordinates). The char-2 mirror of
    /// [`FiniteOddField::from_index`](crate::forms::FiniteOddField::from_index).
    fn from_index(i: u128) -> Self;

    /// The **Artin–Schreier class** `Tr_{F_q/F₂}(x) ∈ {0, 1}` — the additive
    /// analogue of the odd-characteristic square class. `x ∈ ℘(F_q)` (the image of
    /// `y ↦ y² + y`) **iff** this is `0`; equivalently `y² + y = x` is solvable iff
    /// the class is `0`. `F₂`-linear in `x`.
    fn artin_schreier_class(x: Self) -> u128;

    /// `Some(())` iff this is a supported char-2 finite field (the char-2 mirror of
    /// [`FiniteOddField::ensure_supported`](crate::forms::FiniteOddField::ensure_supported)).
    fn ensure_supported() -> Option<()> {
        Self::is_supported_char2_field().then_some(())
    }
}

impl FiniteChar2Field for Fp<2> {
    fn field_order() -> u128 {
        2
    }

    fn is_supported_char2_field() -> bool {
        Fp::<2>::modulus_is_prime()
    }

    fn from_i128(n: i128) -> Self {
        Fp::<2>::new(n)
    }

    fn from_index(i: u128) -> Self {
        Fp::<2>::from_u128(i)
    }

    fn artin_schreier_class(x: Self) -> u128 {
        // Tr_{F₂/F₂} is the identity; the class is the bit itself.
        x.value() & 1
    }
}

impl<const N: usize> FiniteChar2Field for Fpn<2, N> {
    fn field_order() -> u128 {
        Fpn::<2, N>::field_order()
    }

    fn is_supported_char2_field() -> bool {
        Fpn::<2, N>::is_supported_field()
    }

    fn from_i128(n: i128) -> Self {
        Fpn::<2, N>::constant(n.rem_euclid(2) as u128)
    }

    fn from_index(i: u128) -> Self {
        // base-2 digits of `i` are the polynomial-basis coordinates of the element.
        let mut digits = [0u128; N];
        let mut x = i;
        for d in digits.iter_mut() {
            *d = x & 1;
            x >>= 1;
        }
        Fpn::<2, N>::from_coeffs(&digits)
    }

    fn artin_schreier_class(x: Self) -> u128 {
        use crate::scalar::FieldExtension;
        // The absolute trace `F_{2^N} → F₂` (= the relative trace to the prime
        // subfield), which realises `F/℘(F) ≅ F₂`.
        x.trace().value()
    }
}

/// The Artin–Schreier class over any supported char-2 finite field — the additive
/// analogue of [`is_square_finite`](crate::forms::is_square_finite). Returns
/// `Tr_{F_q/F₂}(x) ∈ {0, 1}`; `0` iff `x ∈ ℘(F_q)`.
pub fn artin_schreier_class_finite<F: FiniteChar2Field>(x: F) -> u128 {
    assert!(
        F::is_supported_char2_field(),
        "characteristic-2 finite-field form theory needs a supported char-2 field"
    );
    F::artin_schreier_class(x)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scalar::Scalar;

    #[test]
    fn f2_class_is_the_identity() {
        // Tr_{F₂/F₂} = id, and ℘(F₂) = {0} (℘(0)=0, ℘(1)=1²+1=0 ⇒ only 0 has class 0).
        assert_eq!(Fp::<2>::artin_schreier_class(Fp::<2>::new(0)), 0);
        assert_eq!(Fp::<2>::artin_schreier_class(Fp::<2>::new(1)), 1);
        assert_eq!(Fp::<2>::field_order(), 2);
        assert!(Fp::<2>::is_supported_char2_field());
    }

    #[test]
    fn f4_artin_schreier_class_matches_solvability() {
        // F₄ = {0, 1, α, α+1}, α² = α+1.  ℘(z) = z²+z:  ℘(0)=0, ℘(1)=0, ℘(α)=1,
        // ℘(α+1)=1 ⇒ ℘(F₄) = {0,1}.  So class 0 on {0,1} (the trace-zero hyperplane),
        // class 1 on {α, α+1}.
        type F4 = Fpn<2, 2>;
        let expect = [(0u128, 0u128), (1, 0), (2, 1), (3, 1)]; // index → class
        for (i, c) in expect {
            assert_eq!(F4::artin_schreier_class(F4::from_index(i)), c, "index {i}");
        }
        assert_eq!(F4::field_order(), 4);

        // class(x) = 0  ⟺  y²+y = x has a solution in F₄ (Artin–Schreier).
        for i in 0..4u128 {
            let x = F4::from_index(i);
            let solvable = (0..4u128)
                .map(F4::from_index)
                .any(|y| y.mul(&y).add(&y) == x);
            assert_eq!(
                F4::artin_schreier_class(x) == 0,
                solvable,
                "AS solvability {i}"
            );
        }
    }

    #[test]
    fn class_is_f2_linear() {
        // The trace is F₂-linear: class(x+y) = class(x) ⊕ class(y). Exhaustive on F₈.
        type F8 = Fpn<2, 3>;
        for i in 0..8u128 {
            for j in 0..8u128 {
                let (x, y) = (F8::from_index(i), F8::from_index(j));
                assert_eq!(
                    F8::artin_schreier_class(x.add(&y)),
                    F8::artin_schreier_class(x) ^ F8::artin_schreier_class(y),
                    "additivity at ({i},{j})"
                );
            }
        }
        // The trace is surjective onto F₂: exactly half of F₈ has class 1.
        let ones = (0..8u128)
            .filter(|&i| F8::artin_schreier_class(F8::from_index(i)) == 1)
            .count();
        assert_eq!(ones, 4);
    }

    #[test]
    fn generic_helper_is_usable() {
        // The capability is usable generically, the point of the trait. The trace
        // is a surjective F₂-linear map, so its kernel ℘(F_q) is an index-2
        // hyperplane: exactly q/2 elements have class 1 (= 1 for F₂, where ℘(F₂)={0}).
        fn class_one_count<F: FiniteChar2Field>() -> u128 {
            (0..F::field_order())
                .filter(|&i| artin_schreier_class_finite(F::from_index(i)) == 1)
                .count() as u128
        }
        assert_eq!(class_one_count::<Fp<2>>(), Fp::<2>::field_order() / 2);
        assert_eq!(
            class_one_count::<Fpn<2, 2>>(),
            Fpn::<2, 2>::field_order() / 2
        );
        assert_eq!(
            class_one_count::<Fpn<2, 3>>(),
            Fpn::<2, 3>::field_order() / 2
        );
    }
}
