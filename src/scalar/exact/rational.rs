//! Exact rational ℚ — *not* a game backend, just the char-0 scalar used to
//! validate the geometric-product engine against the known Cl(p,q)
//! classification before trusting the exotic backends. (The surreal backend is
//! the real char-0 home.)

use crate::scalar::Scalar;
use std::cmp::Ordering;
use std::fmt;

/// Exact rational over i128, used only for engine validation. Overflow is a
/// known limitation — fine for the small forms in the test suite, not meant
/// for serious arithmetic. (The surreal backend is the real char-0 home.)
#[derive(Clone)]
pub struct Rational {
    num: i128,
    den: i128, // always > 0, gcd(num, den) == 1
}

fn gcd_u128(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        let t = b;
        b = a % b;
        a = t;
    }
    a
}

/// Exact integer square root of `n ≥ 0`, or `None` if `n` is not a perfect
/// square.
fn isqrt_exact(n: i128) -> Option<i128> {
    if n < 0 {
        return None;
    }
    if n == 0 {
        return Some(0);
    }
    let mut x = (n as f64).sqrt() as i128;
    while x > 0 && x.checked_mul(x).is_none_or(|v| v > n) {
        x -= 1;
    }
    while (x + 1).checked_mul(x + 1).is_some_and(|v| v <= n) {
        x += 1;
    }
    if x.checked_mul(x) == Some(n) {
        Some(x)
    } else {
        None
    }
}

/// Exact integer `k`-th root of `n` (allowing negative `n` for odd `k`), or
/// `None` if `n` is not a perfect `k`-th power.
fn inth_root_exact(n: i128, k: u128) -> Option<i128> {
    if k == 0 {
        return None;
    }
    if k == 1 {
        return Some(n);
    }
    if n == 0 {
        return Some(0);
    }
    let neg = n < 0;
    if neg && k.is_multiple_of(2) {
        return None; // no real even root of a negative
    }
    let a = n.abs();
    let k_pow = k.try_into().ok()?;
    let pw = |b: i128| -> Option<i128> { b.checked_pow(k_pow) };
    let mut x = (a as f64).powf(1.0 / k as f64) as i128;
    while x > 0 && pw(x).is_none_or(|v| v > a) {
        x -= 1;
    }
    while pw(x + 1).is_some_and(|v| v <= a) {
        x += 1;
    }
    if pw(x) == Some(a) {
        Some(if neg { -x } else { x })
    } else {
        None
    }
}

impl Rational {
    pub fn try_new(num: i128, den: i128) -> Option<Self> {
        if den == 0 {
            return None;
        }
        let (num, den) = if den < 0 {
            (num.checked_neg()?, den.checked_neg()?)
        } else {
            (num, den)
        };
        let g = gcd_u128(num.unsigned_abs(), den as u128).max(1);
        let g = i128::try_from(g).ok()?;
        Some(Rational {
            num: num / g,
            den: den / g,
        })
    }

    pub fn new(num: i128, den: i128) -> Self {
        Self::try_new(num, den).expect("Rational::new received zero denominator or overflowed i128")
    }

    /// The integer `n` as an exact rational. This is the ℤ-embedding for `Rational`.
    ///
    /// Kept as a doc'd alias for `Rational::from_int(n)` — a future sweep retires
    /// this spelling once all call sites migrate.
    pub fn int(n: i128) -> Self {
        Rational::from_int(n)
    }

    /// Sign as an Ordering relative to zero (den is always > 0).
    pub fn sign(&self) -> Ordering {
        self.num.cmp(&0)
    }

    /// True iff this rational is a (rational) integer, i.e. its denominator is 1.
    /// Used by the omnific-integer backend to test the constant CNF term.
    pub fn is_integer(&self) -> bool {
        self.den == 1
    }

    /// The numerator (in lowest terms; carries the sign).
    pub fn numer(&self) -> i128 {
        self.num
    }

    /// The denominator (in lowest terms; always > 0).
    pub fn denom(&self) -> i128 {
        self.den
    }

    /// Total order on values (denominator is always positive).
    // Inherent value-order, deliberately kept off `std::cmp::Ord`: orders and
    // operators are opt-in here, not blanket trait impls (see AGENTS.md).
    #[allow(clippy::should_implement_trait)]
    pub fn cmp(&self, other: &Self) -> Ordering {
        self.sub(other).sign()
    }

    /// The greatest integer ≤ this rational.
    pub fn floor(&self) -> i128 {
        self.num.div_euclid(self.den)
    }

    /// The exact rational square root, or `None` if it is not a perfect square
    /// in ℚ (numerator and denominator both perfect squares, and `self ≥ 0`).
    /// `√2` is `None` here on purpose: it is not rational. This is what bounds
    /// the surreal `sqrt` to the ℚ-coefficient subclass.
    pub fn sqrt(&self) -> Option<Rational> {
        let sn = isqrt_exact(self.num)?;
        let sd = isqrt_exact(self.den)?;
        Some(Rational::new(sn, sd))
    }

    /// The exact rational `k`-th root, or `None` if it is not a perfect `k`-th
    /// power in ℚ (even `k` requires `self ≥ 0`).
    pub fn nth_root(&self, k: u128) -> Option<Rational> {
        if k == 0 {
            return None;
        }
        let rn = inth_root_exact(self.num, k)?;
        let rd = inth_root_exact(self.den, k)?; // den > 0
        Some(Rational::new(rn, rd))
    }
}

impl fmt::Display for Rational {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.den == 1 {
            write!(f, "{}", self.num)
        } else {
            write!(f, "{}/{}", self.num, self.den)
        }
    }
}

impl fmt::Debug for Rational {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // delegate to Display: assert_eq! failure output matches Display everywhere
        fmt::Display::fmt(self, f)
    }
}

impl PartialEq for Rational {
    fn eq(&self, other: &Self) -> bool {
        // both are in lowest terms with positive denominator
        self.num == other.num && self.den == other.den
    }
}

impl From<i128> for Rational {
    /// The ℤ-embedding: the unique unital ring homomorphism ℤ → ℚ.
    fn from(n: i128) -> Self {
        Rational::from_int(n)
    }
}

impl Scalar for Rational {
    fn zero() -> Self {
        Rational { num: 0, den: 1 }
    }
    fn one() -> Self {
        Rational { num: 1, den: 1 }
    }
    /// Faster direct construction; semantically identical to the default double-and-add.
    fn from_int(n: i128) -> Self {
        Rational { num: n, den: 1 }
    }
    fn add(&self, rhs: &Self) -> Self {
        let g = gcd_u128(self.den as u128, rhs.den as u128).max(1);
        let g = i128::try_from(g).expect("Rational denominator gcd overflowed i128");
        let lhs_scale = rhs.den / g;
        let rhs_scale = self.den / g;
        let lhs = self
            .num
            .checked_mul(lhs_scale)
            .and_then(|x| x.checked_add(rhs.num.checked_mul(rhs_scale)?));
        let den = self.den.checked_mul(lhs_scale);
        Rational::try_new(
            lhs.expect("Rational addition overflowed i128"),
            den.expect("Rational addition denominator overflowed i128"),
        )
        .expect("Rational addition normalization overflowed i128")
    }
    fn neg(&self) -> Self {
        Rational {
            num: self
                .num
                .checked_neg()
                .expect("Rational negation overflowed i128"),
            den: self.den,
        }
    }
    fn mul(&self, rhs: &Self) -> Self {
        let mut lhs_num = self.num;
        let mut lhs_den = self.den;
        let mut rhs_num = rhs.num;
        let mut rhs_den = rhs.den;

        let g1 = gcd_u128(lhs_num.unsigned_abs(), rhs_den as u128);
        if g1 > 1 {
            let g1 = i128::try_from(g1).expect("Rational multiplication gcd overflowed i128");
            lhs_num /= g1;
            rhs_den /= g1;
        }
        let g2 = gcd_u128(rhs_num.unsigned_abs(), lhs_den as u128);
        if g2 > 1 {
            let g2 = i128::try_from(g2).expect("Rational multiplication gcd overflowed i128");
            rhs_num /= g2;
            lhs_den /= g2;
        }

        Rational::try_new(
            lhs_num
                .checked_mul(rhs_num)
                .expect("Rational multiplication numerator overflowed i128"),
            lhs_den
                .checked_mul(rhs_den)
                .expect("Rational multiplication denominator overflowed i128"),
        )
        .expect("Rational multiplication normalization overflowed i128")
    }
    fn characteristic() -> u128 {
        0
    }
    fn inv(&self) -> Option<Self> {
        if self.num == 0 {
            None
        } else {
            Some(Rational::new(self.den, self.num))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scalar::Scalar;

    #[test]
    fn rational_arithmetic() {
        let half = Rational::new(1, 2);
        let third = Rational::new(1, 3);
        assert_eq!(half.add(&third), Rational::new(5, 6));
        assert_eq!(half.mul(&third), Rational::new(1, 6));
        assert_eq!(half.sub(&half), Rational::zero());
        assert_eq!(half.add(&half), Rational::one());
        assert_eq!(Rational::new(2, 4), Rational::new(1, 2)); // reduction
    }

    #[test]
    fn rational_adds_before_denominator_product_overflows() {
        let huge_den = 1i128 << 100;
        let x = Rational::new(1, huge_den);
        assert_eq!(x.add(&x), Rational::new(1, 1i128 << 99));
    }
}
