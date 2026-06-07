//! Odd finite-field capability trait and scalar-level helpers.

use crate::scalar::{Fp, Fpn, Scalar};

pub(super) fn assert_odd_prime<const P: u128>() {
    Fp::<P>::assert_prime_modulus();
    assert!(P != 2, "odd-characteristic form theory needs P odd");
}

/// Finite fields of odd characteristic, with the operations the form classifiers
/// actually need: field-order metadata, base-field constants, and square classes.
/// This is intentionally narrower than [`Scalar`]: it is a form-theory façade, not
/// a new scalar-world requirement.
pub trait FiniteOddField: Scalar + Copy {
    /// Characteristic prime `p`.
    fn characteristic_prime() -> u128;

    /// Field order `q = p^n`.
    fn field_order() -> u128;

    /// Whether this type is a supported finite field of odd characteristic.
    fn is_supported_odd_field() -> bool;

    /// Embed an ordinary integer through the prime subfield.
    fn from_i128(n: i128) -> Self;

    /// Square-class test in the field. `0` counts as a square.
    fn is_square_value(x: Self) -> bool;

    fn ensure_supported() -> Option<()> {
        Self::is_supported_odd_field().then_some(())
    }
}

impl<const P: u128> FiniteOddField for Fp<P> {
    fn characteristic_prime() -> u128 {
        P
    }

    fn field_order() -> u128 {
        P
    }

    fn is_supported_odd_field() -> bool {
        Fp::<P>::modulus_is_prime() && P != 2
    }

    fn from_i128(n: i128) -> Self {
        Fp::<P>::new(n)
    }

    fn is_square_value(x: Self) -> bool {
        is_square(x)
    }
}

impl<const P: u128, const N: usize> FiniteOddField for Fpn<P, N> {
    fn characteristic_prime() -> u128 {
        P
    }

    fn field_order() -> u128 {
        Fpn::<P, N>::order()
    }

    fn is_supported_odd_field() -> bool {
        Fp::<P>::modulus_is_prime() && P != 2 && N > 0
    }

    fn from_i128(n: i128) -> Self {
        let m = P as i128;
        let v = ((n % m) + m) % m;
        Fpn::<P, N>::constant(v as u128)
    }

    fn is_square_value(x: Self) -> bool {
        x.is_square()
    }
}

/// `base^e` in `F_P` by square-and-multiply.
fn fp_pow<const P: u128>(mut base: Fp<P>, mut e: u128) -> Fp<P> {
    let mut acc = Fp::<P>::one();
    while e > 0 {
        if e & 1 == 1 {
            acc = acc.mul(&base);
        }
        base = base.mul(&base);
        e >>= 1;
    }
    acc
}

/// Euler's criterion: is `x` a square in `F_P`? (`0` counts as a square.)
pub fn is_square<const P: u128>(x: Fp<P>) -> bool {
    assert_odd_prime::<P>();
    if x.is_zero() {
        return true;
    }
    fp_pow(x, (P - 1) / 2) == Fp::<P>::one()
}

/// Square-class predicate over any supported finite field of odd characteristic.
pub fn is_square_finite<F: FiniteOddField>(x: F) -> bool {
    assert!(
        F::is_supported_odd_field(),
        "odd-characteristic finite-field form theory needs odd finite fields"
    );
    F::is_square_value(x)
}

/// The Hilbert symbol `(a, b)` over `F_P`: `+1` iff `z² = a x² + b y²` has a
/// nontrivial solution. Over a finite field this is identically `+1` for nonzero
/// `a, b` (computed by an honest search, which always succeeds).
pub fn hilbert_symbol<const P: u128>(a: Fp<P>, b: Fp<P>) -> i8 {
    assert_odd_prime::<P>();
    for x in 0..P {
        for y in 0..P {
            for z in 0..P {
                if x == 0 && y == 0 && z == 0 {
                    continue;
                }
                let (fx, fy, fz) = (Fp::<P>(x), Fp::<P>(y), Fp::<P>(z));
                let rhs = a.mul(&fx.mul(&fx)).add(&b.mul(&fy.mul(&fy)));
                if fz.mul(&fz) == rhs {
                    return 1;
                }
            }
        }
    }
    -1
}
