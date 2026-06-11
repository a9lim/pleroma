//! Dense univariate polynomials `S[x]` over a [`Scalar`] base, low-degree-first.
//!
//! The crate's shared polynomial primitive. It backs two clients:
//!
//!   * [`Gauss`](crate::scalar::Gauss) — the rational function field `S(t)` stores
//!     `num/den` as a pair of `Poly`s (this type absorbed the private helpers that
//!     used to live in `functor/gauss.rs`).
//!   * the global **function field** `F_q(t)` and its place/Hilbert-symbol layer
//!     in [`forms::function_field`](crate::forms) — which additionally needs
//!     division, gcd, and modular powers (the residue quadratic character is
//!     Euler's criterion `u^{(|κ|−1)/2}` computed in `F_q[t]/(π)`).
//!
//! Representation is **trimmed** (no trailing zero coefficients; the zero
//! polynomial is the empty vector), so `PartialEq` is structural and exact. The
//! division-flavoured methods (`divrem`, `rem`, `make_monic`, `gcd`, `*_mod`)
//! assume the base is a **field** — they invert the divisor's leading coefficient
//! and panic if it is not invertible. Both clients are fields, so this is the same
//! honesty as `Gauss`'s `inv = den/num`.

use crate::scalar::{Scalar, Valued};

/// A dense univariate polynomial over `S`, coefficients low-degree-first and
/// trimmed (leading coefficient nonzero; the zero polynomial is empty).
#[derive(Clone, PartialEq)]
pub struct Poly<S: Scalar> {
    coeffs: Vec<S>,
}

/// Display v2 (§9) operational atomicity: a coefficient rendering attaches bare
/// iff it contains no spaces and no operator character (`⋅ ∧ ↑ / + -`) outside
/// balanced parentheses; otherwise it is wrapped so `coeff⋅t↑i` stays
/// unambiguous (`(x + 1)⋅t↑2`, but `x⋅t↑2`).
pub(crate) fn atomic(s: &str) -> bool {
    let mut depth: i32 = 0;
    for ch in s.chars() {
        match ch {
            '(' => depth += 1,
            ')' => depth -= 1,
            ' ' if depth == 0 => return false,
            '⋅' | '∧' | '↑' | '/' | '+' | '-' if depth == 0 => return false,
            _ => {}
        }
    }
    true
}

/// Attach a scalar coefficient to a label as `coeff⋅label`, parenthesizing the
/// coefficient only when its rendering is non-atomic (§9). A single leading `-`
/// is a unary sign, not an internal operator, so it is checked separately and
/// carried through bare (`-2⋅e0∧e1`); the Multivector join rule then lifts it to
/// a ` - ` separator. Only a `-`/operator/space *inside* the magnitude forces
/// parens (`(x + 1)⋅e0∧e1`).
pub(crate) fn attach_coeff<S: Scalar>(c: &S, label: &str) -> String {
    let cs = c.to_string();
    let (sign, mag) = match cs.strip_prefix('-') {
        Some(rest) => ("-", rest),
        None => ("", cs.as_str()),
    };
    if atomic(mag) {
        format!("{sign}{mag}⋅{label}")
    } else {
        format!("({cs})⋅{label}")
    }
}

impl<S: Scalar> std::fmt::Display for Poly<S> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.coeffs.is_empty() {
            return write!(f, "0");
        }
        let mut parts = Vec::new();
        for (i, c) in self.coeffs.iter().enumerate() {
            if c.is_zero() {
                continue;
            }
            // Display v2 (§9): variable `t` (matches `F_q[t]`), explicit `⋅`,
            // coefficient parens only when the coefficient renders non-atomically.
            parts.push(match i {
                0 => format!("{c}"),
                1 => attach_coeff(c, "t"),
                _ => attach_coeff(c, &format!("t↑{i}")),
            });
        }
        write!(f, "{}", parts.join(" + "))
    }
}

impl<S: Scalar> std::fmt::Debug for Poly<S> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        std::fmt::Display::fmt(self, f)
    }
}

/// Drop trailing zero coefficients so the leading term is nonzero.
fn trim<S: Scalar>(mut p: Vec<S>) -> Vec<S> {
    while p.last().map(|c| c.is_zero()).unwrap_or(false) {
        p.pop();
    }
    p
}

impl<S: Scalar> Poly<S> {
    /// Build a polynomial from low-degree-first coefficients (trimmed).
    pub fn new(coeffs: Vec<S>) -> Self {
        Poly {
            coeffs: trim(coeffs),
        }
    }

    /// The zero polynomial.
    pub fn zero() -> Self {
        Poly { coeffs: Vec::new() }
    }

    /// The constant polynomial `1`.
    pub fn one() -> Self {
        Poly::constant(S::one())
    }

    /// The constant polynomial `s`.
    pub fn constant(s: S) -> Self {
        Poly::new(vec![s])
    }

    /// The indeterminate `x`.
    pub fn x() -> Self {
        Poly::new(vec![S::zero(), S::one()])
    }

    /// `coeff · x^deg`.
    pub fn monomial(deg: usize, coeff: S) -> Self {
        let mut c = vec![S::zero(); deg];
        c.push(coeff);
        Poly::new(c)
    }

    /// The coefficient slice (low-degree-first; empty iff zero).
    pub fn coeffs(&self) -> &[S] {
        &self.coeffs
    }

    pub fn is_zero(&self) -> bool {
        self.coeffs.is_empty()
    }

    /// The degree, or `None` for the zero polynomial.
    pub fn degree(&self) -> Option<usize> {
        self.coeffs.len().checked_sub(1)
    }

    /// The leading coefficient, or `None` for the zero polynomial.
    pub fn leading(&self) -> Option<&S> {
        self.coeffs.last()
    }

    /// The coefficient of `x^i` (zero past the degree).
    pub fn coeff(&self, i: usize) -> S {
        self.coeffs.get(i).cloned().unwrap_or_else(S::zero)
    }

    pub fn add(&self, rhs: &Self) -> Self {
        let n = self.coeffs.len().max(rhs.coeffs.len());
        let mut out = Vec::with_capacity(n);
        for i in 0..n {
            out.push(self.coeff(i).add(&rhs.coeff(i)));
        }
        Poly::new(out)
    }

    pub fn neg(&self) -> Self {
        Poly {
            coeffs: self.coeffs.iter().map(|c| c.neg()).collect(),
        }
    }

    pub fn sub(&self, rhs: &Self) -> Self {
        self.add(&rhs.neg())
    }

    pub fn mul(&self, rhs: &Self) -> Self {
        if self.is_zero() || rhs.is_zero() {
            return Poly::zero();
        }
        let mut out = vec![S::zero(); self.coeffs.len() + rhs.coeffs.len() - 1];
        for (i, x) in self.coeffs.iter().enumerate() {
            if x.is_zero() {
                continue;
            }
            for (j, y) in rhs.coeffs.iter().enumerate() {
                out[i + j] = out[i + j].add(&x.mul(y));
            }
        }
        Poly::new(out)
    }

    /// Multiply every coefficient by `s`.
    pub fn scale(&self, s: &S) -> Self {
        Poly::new(self.coeffs.iter().map(|c| c.mul(s)).collect())
    }

    /// Evaluate at `x` by Horner's rule.
    pub fn eval(&self, x: &S) -> S {
        let mut acc = S::zero();
        for c in self.coeffs.iter().rev() {
            acc = acc.mul(x).add(c);
        }
        acc
    }

    /// Scale to a monic polynomial (divide through by the leading coefficient).
    /// Panics on the zero polynomial; requires the base to be a field.
    pub fn make_monic(&self) -> Self {
        let lead = self.leading().expect("make_monic of the zero polynomial");
        let inv = lead
            .inv()
            .expect("a field's nonzero leading coefficient inverts");
        self.scale(&inv)
    }

    /// Euclidean division `self = q·divisor + r` with `deg r < deg divisor`,
    /// returning `(q, r)`. Requires `divisor` nonzero over a field.
    pub fn divrem(&self, divisor: &Self) -> (Self, Self) {
        let dd = divisor
            .degree()
            .expect("polynomial division by the zero polynomial");
        let dlead_inv = divisor
            .leading()
            .unwrap()
            .inv()
            .expect("a field's nonzero leading coefficient inverts");
        let mut rem = self.coeffs.clone();
        let mut quot = vec![S::zero(); self.coeffs.len().saturating_sub(dd).max(1)];
        loop {
            rem = trim(rem);
            let rdeg = match rem.len().checked_sub(1) {
                Some(d) if d >= dd => d,
                _ => break,
            };
            let shift = rdeg - dd;
            let factor = rem[rdeg].mul(&dlead_inv);
            quot[shift] = factor.clone();
            for (i, dc) in divisor.coeffs.iter().enumerate() {
                rem[shift + i] = rem[shift + i].sub(&factor.mul(dc));
            }
        }
        (Poly::new(quot), Poly::new(rem))
    }

    /// The remainder `self mod divisor`.
    pub fn rem(&self, divisor: &Self) -> Self {
        self.divrem(divisor).1
    }

    /// Whether `divisor` divides `self` exactly.
    pub fn divides(&self, multiple: &Self) -> bool {
        !self.is_zero() && multiple.rem(self).is_zero()
    }

    /// The monic gcd (the zero polynomial's gcd partner is returned monic).
    pub fn gcd(&self, other: &Self) -> Self {
        let mut a = self.clone();
        let mut b = other.clone();
        while !b.is_zero() {
            let r = a.rem(&b);
            a = b;
            b = r;
        }
        if a.is_zero() {
            a
        } else {
            a.make_monic()
        }
    }

    /// `self · other mod modulus`.
    pub fn mul_mod(&self, other: &Self, modulus: &Self) -> Self {
        self.mul(other).rem(modulus)
    }

    /// `self^e mod modulus` by square-and-multiply.
    pub fn pow_mod(&self, mut e: u128, modulus: &Self) -> Self {
        let mut acc = Poly::one().rem(modulus);
        let mut base = self.rem(modulus);
        while e > 0 {
            if e & 1 == 1 {
                acc = acc.mul_mod(&base, modulus);
            }
            base = base.mul_mod(&base, modulus);
            e >>= 1;
        }
        acc
    }
}

/// `S[t]` is itself a commutative ring — the **ring of integers** of the rational
/// function field [`RationalFunction`](crate::scalar::RationalFunction)`<S> = S(t)`.
/// Its units are the nonzero constants (so `inv` is partial), exactly as `ℤ` sits
/// inside `ℚ`. The trait methods delegate to the inherent ones (inherent shadows
/// trait at the receiver, so this delegates rather than recurses).
impl<S: Scalar> Scalar for Poly<S> {
    fn zero() -> Self {
        Self::constant(S::zero()) // trims to the empty polynomial
    }
    fn one() -> Self {
        Self::constant(S::one())
    }
    fn add(&self, rhs: &Self) -> Self {
        self.add(rhs)
    }
    fn neg(&self) -> Self {
        self.neg()
    }
    fn mul(&self, rhs: &Self) -> Self {
        self.mul(rhs)
    }
    fn characteristic() -> u128 {
        S::characteristic()
    }
    fn inv(&self) -> Option<Self> {
        // units of S[t] are the nonzero constants.
        match self.degree() {
            Some(0) => self.coeff(0).inv().map(Self::constant),
            _ => None,
        }
    }
    fn is_zero(&self) -> bool {
        self.coeffs.is_empty()
    }
}

impl<S: Valued> Poly<S> {
    /// The minimum coefficient valuation (the Gauss valuation of the polynomial),
    /// or `None` for the zero polynomial.
    pub fn min_coeff_valuation(&self) -> Option<i128> {
        self.coeffs.iter().filter_map(|c| c.valuation()).min()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scalar::Fp;

    type P5 = Poly<Fp<5>>;

    fn p(coeffs: &[i128]) -> P5 {
        Poly::new(coeffs.iter().map(|&n| Fp::<5>::from_int(n)).collect())
    }

    #[test]
    fn arithmetic_basics() {
        // (1 + x)(1 + x) = 1 + 2x + x²
        let one_plus_x = p(&[1, 1]);
        assert_eq!(one_plus_x.mul(&one_plus_x), p(&[1, 2, 1]));
        // (1 + x) + (4 + 4x) = 5 + 5x ≡ 0 in F_5
        assert_eq!(p(&[1, 1]).add(&p(&[4, 4])), P5::zero());
        assert_eq!(p(&[1, 1]).neg(), p(&[4, 4]));
        assert_eq!(P5::x().eval(&Fp::<5>::from_int(3)), Fp::<5>::from_int(3));
        assert_eq!(
            p(&[1, 1, 1]).eval(&Fp::<5>::from_int(2)),
            Fp::<5>::from_int(7)
        ); // 1+2+4=7
    }

    #[test]
    fn euclidean_division() {
        // x² − 1 = (x − 1)(x + 1) over F_5  (−1 ≡ 4)
        let x2m1 = p(&[4, 0, 1]);
        let xm1 = p(&[4, 1]); // x − 1
        let (q, r) = x2m1.divrem(&xm1);
        assert_eq!(q, p(&[1, 1])); // x + 1
        assert!(r.is_zero());
        assert!(xm1.divides(&x2m1));
        // a remainder that is genuinely nonzero
        let (_, r2) = p(&[1, 0, 1]).divrem(&xm1); // x² + 1 at x=1 → 2
        assert_eq!(r2, p(&[2]));
    }

    #[test]
    fn gcd_and_monic() {
        // gcd(x² − 1, x² + 2x + 1) = x + 1 (monic)
        let g = p(&[4, 0, 1]).gcd(&p(&[1, 2, 1]));
        assert_eq!(g, p(&[1, 1]));
        // make_monic divides through by the leading coeff: 2x + 2 → x + 1
        assert_eq!(p(&[2, 2]).make_monic(), p(&[1, 1]));
    }

    #[test]
    fn display_v2_canonical_ogham() {
        use crate::scalar::Fpn;
        // Atomic coefficients (a single F_5 residue) attach bare: `1 + 2⋅t`.
        assert_eq!(p(&[1, 2]).to_string(), "1 + 2⋅t");
        assert_eq!(p(&[0, 0, 3]).to_string(), "3⋅t↑2");
        assert_eq!(P5::zero().to_string(), "0");
        // Non-atomic coefficients (an F_8 element `x + 1`) parenthesize.
        // (Poly renders low-degree-first and does not suppress a coefficient-1
        // term — the §9 delta only changes `x→t`, `·→⋅`, and the paren rule.)
        type Q = Poly<Fpn<2, 3>>;
        let xp1 = Fpn::<2, 3>::from_coeffs(&[1, 1]); // x + 1 (non-atomic)
        let x = Fpn::<2, 3>::from_coeffs(&[0, 1]); // x (atomic)
        let one = Fpn::<2, 3>::one();
        // 1 + x⋅t + (x + 1)⋅t↑2
        let poly = Q::new(vec![one, x, xp1]);
        assert_eq!(poly.to_string(), "1 + x⋅t + (x + 1)⋅t↑2");
    }

    #[test]
    fn atomicity_rule() {
        assert!(atomic("42"));
        assert!(atomic("*5"));
        assert!(atomic("*ω"));
        assert!(atomic("x"));
        assert!(atomic("*(ω⋅7)")); // operators only inside balanced parens
        assert!(!atomic("x + 1"));
        assert!(!atomic("ω↑-1"));
        assert!(!atomic("3⋅x")); // bare `⋅`
    }

    #[test]
    fn modular_powers_for_eulers_criterion() {
        // In F_5[x]/(x² + 2) (x² ≡ −2 ≡ 3), the residue field is F_25.
        let modulus = p(&[2, 0, 1]); // x² + 2, irreducible over F_5 (−2=3 is a nonsquare)
                                     // x^(25−1) ≡ 1 (Fermat in F_25*), and x is a nonsquare ⇒ x^((25−1)/2) ≡ −1.
        assert_eq!(P5::x().pow_mod(24, &modulus), P5::one());
        assert_eq!(P5::x().pow_mod(12, &modulus), p(&[4])); // −1 ≡ 4
    }
}
