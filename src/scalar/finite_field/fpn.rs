//! Finite extension fields `F_{p^n}` — completing the field tower in every
//! characteristic.
//!
//! The odd-characteristic leg of the crate only had the *prime* fields `Fp<P>`;
//! characteristic 2 had the whole nimber tower (`F_{2^{2^k}}`). `Fpn<P, N>` closes
//! that asymmetry: it is `F_{p^n}` for any prime `P` and positive `N` whose order
//! fits in the crate's `u128` payload model. It also supplies the **char-2
//! odd-degree** fields the nimbers cannot reach — the finite nimbers realise only
//! `F_{2^{2^k}}` (degrees that are powers of two), so `F_8` (degree 3) is not a
//! nimber subfield; `Fpn<2, 3>` is the way to get it here.
//!
//! ## The const-generic modulus, two parameters
//!
//! Like `Fp<P>`, the modulus lives in the **type** (`Scalar::zero()/one()` take no
//! `self`). A field is `Fpn<const P: u128, const N: usize>` = `F_{p^N}`, carried as the
//! `N` coefficients of `c_0 + c_1 x + … + c_{N-1} x^{N-1}` with each `c_i ∈ [0, P)`.
//! A different `(P, N)` is a different type — the same no-mixing discipline the rest
//! of the crate uses. `Fpn<2, 2>` is "the polynomial-basis `F_4`", a *different type*
//! from (but isomorphic to) the nimber `F_4`; the value-add over the nimbers is the
//! odd-degree char-2 layers and the odd-`p` extensions.
//!
//! ## The reduction polynomial
//!
//! Arithmetic is in `F_p[x] / (m(x))` for a monic irreducible `m` of degree `N`.
//! `reduction` returns the low coefficients `r` of the reduction rule
//! `x^N = Σ_i r_i x^i` (i.e. `m(x) = x^N − Σ_i r_i x^i`). Extension fields are opened
//! by a deterministic search for the first monic irreducible polynomial, certified by
//! Rabin's irreducibility test and cached per `(P,N)`. The old small Conway rows are
//! retained only as test oracles; the runtime model is an honest generated
//! "irreducible polynomial" model, not a compatible Conway embedding. `mul` is
//! schoolbook multiply-then-reduce — the degree-`N`, odd-`p` generalisation of
//! `big::ordinal`'s "reduce mod `ω³ = 2`".

use super::fp::{add_mod, mul_mod};
use super::FiniteField;
use crate::scalar::{add_mod_u128, is_prime_u128, mod_inverse_u128, sub_mod_u128, Fp, Scalar};
use std::collections::BTreeMap;
use std::fmt;
use std::sync::{Mutex, OnceLock};

/// An element of `F_{p^N}`: the coefficients of `c_0 + c_1 x + … + c_{N-1} x^{N-1}`,
/// each reduced into `[0, P)`.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct Fpn<const P: u128, const N: usize>([u128; N]);

/// Provenance of the shipped reduction polynomial for an `Fpn<P,N>` backend.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ReductionPolynomialKind {
    /// Degree-1 prime field, so no extension polynomial is needed.
    PrimeField,
    /// A curated table entry is the Conway polynomial in this polynomial basis.
    /// Production `Fpn` generation no longer returns this tag; old rows use it only
    /// as test-oracle vocabulary.
    Conway,
    /// A curated table entry is verified irreducible, but not claimed as Conway data.
    /// Production `Fpn` generation no longer returns this tag; old rows use it only
    /// as test-oracle vocabulary.
    Irreducible,
    /// The entry was generated deterministically and verified irreducible by Rabin's test.
    GeneratedIrreducible,
}

/// Low coefficients `r` of the reduction rule `x^N = Σ_i r_i x^i`. Each returned
/// slice has length `N`. Degree `1` has the vacuous rule `[0]`; every extension
/// degree is generated deterministically and cached.
pub(crate) fn reduction<const P: u128, const N: usize>() -> &'static [u128] {
    if N == 1 {
        return &[0];
    }
    generated_reduction(P, N)
}

/// Metadata companion to [`reduction`].
pub(crate) fn reduction_kind<const P: u128, const N: usize>() -> ReductionPolynomialKind {
    if N == 1 {
        return ReductionPolynomialKind::PrimeField;
    }
    assert_generated_domain(P, N);
    ReductionPolynomialKind::GeneratedIrreducible
}

type ReductionCache = BTreeMap<(u128, usize), &'static [u128]>;

fn generated_reductions() -> &'static Mutex<ReductionCache> {
    static CACHE: OnceLock<Mutex<ReductionCache>> = OnceLock::new();
    CACHE.get_or_init(|| Mutex::new(BTreeMap::new()))
}

fn generated_reduction(p: u128, n: usize) -> &'static [u128] {
    assert_generated_domain(p, n);
    let cache = generated_reductions();
    let mut guard = cache.lock().expect("Fpn reduction cache poisoned");
    if let Some(&rule) = guard.get(&(p, n)) {
        return rule;
    }
    let rule = deterministic_irreducible_reduction(p, n);
    let leaked: &'static [u128] = Box::leak(rule.into_boxed_slice());
    guard.insert((p, n), leaked);
    leaked
}

fn assert_generated_domain(p: u128, n: usize) {
    assert!(is_prime_u128(p), "Fpn<{p},{n}> needs prime P");
    assert!(n > 0, "Fpn<{p},{n}> needs N > 0");
    assert!(
        field_order_for(p, n).is_some(),
        "Fpn<{p},{n}> field order exceeds u128"
    );
}

fn field_order_for(p: u128, n: usize) -> Option<u128> {
    if n == 0 {
        return None;
    }
    let mut acc = 1u128;
    for _ in 0..n {
        acc = acc.checked_mul(p)?;
    }
    Some(acc)
}

fn deterministic_irreducible_reduction(p: u128, n: usize) -> Vec<u128> {
    let candidates = field_order_for(p, n).expect("generated Fpn domain checked");
    for code in 0..candidates {
        let rule = decode_reduction_code(code, p, n);
        if rule[0] == 0 {
            continue; // monic irreducible degree > 1 cannot have zero constant term
        }
        let modulus = reduction_rule_to_polynomial(&rule, p);
        if is_irreducible_monic(&modulus, p) {
            return rule;
        }
    }
    panic!("Fpn<{p},{n}>: no irreducible polynomial found");
}

fn decode_reduction_code(mut code: u128, p: u128, n: usize) -> Vec<u128> {
    let mut rule = vec![0u128; n];
    for slot in &mut rule {
        *slot = code % p;
        code /= p;
    }
    rule
}

fn reduction_rule_to_polynomial(rule: &[u128], p: u128) -> Vec<u128> {
    let mut poly: Vec<u128> = rule.iter().map(|&c| sub_mod_u128(0, c, p)).collect();
    poly.push(1);
    trim_poly(poly)
}

fn is_irreducible_monic(poly: &[u128], p: u128) -> bool {
    let n = match poly_degree(poly) {
        Some(d) if d > 0 && poly[d] == 1 => d,
        _ => return false,
    };
    if n == 1 {
        return true;
    }

    let x = vec![0, 1];
    for prime in distinct_prime_factors_usize(n) {
        let exp = checked_pow_u128(p, n / prime).expect("Fpn Rabin exponent checked");
        let witness = poly_sub(&poly_pow_x_mod(exp, poly, p), &x, p);
        if !poly_is_one(&poly_gcd(poly.to_vec(), witness, p)) {
            return false;
        }
    }
    let exp = checked_pow_u128(p, n).expect("Fpn Rabin exponent checked");
    poly_sub(&poly_pow_x_mod(exp, poly, p), &x, p).is_empty()
}

fn trim_poly(mut poly: Vec<u128>) -> Vec<u128> {
    while poly.last() == Some(&0) {
        poly.pop();
    }
    poly
}

fn poly_degree(poly: &[u128]) -> Option<usize> {
    poly.iter().rposition(|&c| c != 0)
}

fn poly_is_one(poly: &[u128]) -> bool {
    poly == [1]
}

fn poly_coeff(poly: &[u128], i: usize) -> u128 {
    poly.get(i).copied().unwrap_or(0)
}

fn poly_sub(a: &[u128], b: &[u128], p: u128) -> Vec<u128> {
    let len = a.len().max(b.len());
    let mut out = Vec::with_capacity(len);
    for i in 0..len {
        out.push(sub_mod_u128(poly_coeff(a, i), poly_coeff(b, i), p));
    }
    trim_poly(out)
}

fn poly_mul_mod(a: &[u128], b: &[u128], modulus: &[u128], p: u128) -> Vec<u128> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut out = vec![0u128; a.len() + b.len() - 1];
    for (i, &ai) in a.iter().enumerate() {
        if ai == 0 {
            continue;
        }
        for (j, &bj) in b.iter().enumerate() {
            out[i + j] = add_mod_u128(out[i + j], crate::scalar::mul_mod_u128(ai, bj, p), p);
        }
    }
    poly_rem(out, modulus, p)
}

fn poly_pow_x_mod(mut exp: u128, modulus: &[u128], p: u128) -> Vec<u128> {
    let mut acc = vec![1];
    let mut base = poly_rem(vec![0, 1], modulus, p);
    while exp > 0 {
        if exp & 1 == 1 {
            acc = poly_mul_mod(&acc, &base, modulus, p);
        }
        exp >>= 1;
        if exp > 0 {
            base = poly_mul_mod(&base, &base, modulus, p);
        }
    }
    acc
}

fn poly_rem(mut rem: Vec<u128>, modulus: &[u128], p: u128) -> Vec<u128> {
    let md = poly_degree(modulus).expect("polynomial modulus must be nonzero");
    let lead_inv = mod_inverse_u128(modulus[md], p).expect("nonzero finite-field coefficient");
    loop {
        rem = trim_poly(rem);
        let rd = match rem.len().checked_sub(1) {
            Some(d) if d >= md => d,
            _ => break,
        };
        let factor = crate::scalar::mul_mod_u128(rem[rd], lead_inv, p);
        let shift = rd - md;
        if factor != 0 {
            for (i, &mc) in modulus.iter().take(md + 1).enumerate() {
                let term = crate::scalar::mul_mod_u128(factor, mc, p);
                rem[shift + i] = sub_mod_u128(rem[shift + i], term, p);
            }
        }
    }
    trim_poly(rem)
}

fn poly_gcd(mut a: Vec<u128>, mut b: Vec<u128>, p: u128) -> Vec<u128> {
    a = trim_poly(a);
    b = trim_poly(b);
    while !b.is_empty() {
        let r = poly_rem(a, &b, p);
        a = b;
        b = r;
    }
    poly_make_monic(a, p)
}

fn poly_make_monic(poly: Vec<u128>, p: u128) -> Vec<u128> {
    let d = match poly_degree(&poly) {
        Some(d) => d,
        None => return Vec::new(),
    };
    let inv = mod_inverse_u128(poly[d], p).expect("nonzero finite-field coefficient");
    trim_poly(
        poly.into_iter()
            .map(|c| crate::scalar::mul_mod_u128(c, inv, p))
            .collect(),
    )
}

fn distinct_prime_factors_usize(mut n: usize) -> Vec<usize> {
    let mut out = Vec::new();
    let mut d = 2usize;
    while d <= n / d {
        if n.is_multiple_of(d) {
            out.push(d);
            while n.is_multiple_of(d) {
                n /= d;
            }
        }
        d += 1;
    }
    if n > 1 {
        out.push(n);
    }
    out
}

fn checked_pow_u128(base: u128, exp: usize) -> Option<u128> {
    let mut acc = 1u128;
    for _ in 0..exp {
        acc = acc.checked_mul(base)?;
    }
    Some(acc)
}

impl<const P: u128, const N: usize> Fpn<P, N> {
    /// Whether this const-generic pair has a prime base field, positive degree, and
    /// field order fitting the crate's `u128` payload model. Extension polynomials
    /// beyond the curated rows are generated deterministically on first use.
    pub fn is_supported_field() -> bool {
        Fp::<P>::modulus_is_prime() && field_order_for(P, N).is_some()
    }

    pub fn assert_supported_field() {
        assert!(
            Self::is_supported_field(),
            "Fpn<{P},{N}> needs prime P, N>0, and field order fitting u128"
        );
    }

    /// The field order `p^N`, or `None` when it exceeds `u128` (the public payload
    /// model used throughout the crate).
    pub fn field_order_checked() -> Option<u128> {
        if !Fp::<P>::modulus_is_prime() {
            return None;
        }
        field_order_for(P, N)
    }

    /// The field order `p^N`.
    pub fn field_order() -> u128 {
        Self::assert_supported_field();
        field_order_for(P, N).expect("Fpn order checked")
    }

    /// The low coefficients of the reduction rule `x^N = Σ r_i x^i`.
    pub fn reduction_rule() -> &'static [u128] {
        Self::assert_supported_field();
        reduction::<P, N>()
    }

    /// Whether this backend uses a generated irreducible polynomial, or no extension
    /// polynomial because `N = 1`.
    pub fn reduction_polynomial_kind() -> ReductionPolynomialKind {
        Self::assert_supported_field();
        reduction_kind::<P, N>()
    }

    /// `true` exactly when this backend is tagged with Conway polynomial provenance.
    /// The production generator does not currently return Conway-tagged rows; the
    /// method remains a provenance query rather than an irreducibility claim.
    pub fn is_conway_polynomial() -> bool {
        Self::reduction_polynomial_kind() == ReductionPolynomialKind::Conway
    }

    /// Embed a base-field constant `c ∈ F_p` as the degree-0 element.
    pub fn constant(c: u128) -> Self {
        Self::assert_supported_field();
        let mut out = [0u128; N];
        out[0] = c % P;
        Fpn(out)
    }

    /// Build from a coefficient slice (low-to-high), reducing each entry mod `P`.
    /// Extra trailing coefficients beyond `N` must be zero (else it is not an
    /// element of this field).
    pub fn from_coeffs(cs: &[u128]) -> Self {
        Self::assert_supported_field();
        assert!(
            cs.iter().skip(N).all(|&c| c % P == 0),
            "Fpn::from_coeffs received nonzero coefficients beyond degree {N}"
        );
        let mut out = [0u128; N];
        for (i, slot) in out.iter_mut().enumerate() {
            if i < cs.len() {
                *slot = cs[i] % P;
            }
        }
        Fpn(out)
    }

    /// The canonical coefficient array, low degree first.
    pub fn coeffs(&self) -> &[u128; N] {
        &self.0
    }

    /// Consume the field element and return its canonical coefficient array.
    pub fn into_coeffs(self) -> [u128; N] {
        self.0
    }

    /// The coefficient of `x^i`, or zero past the degree.
    pub fn coeff(&self, i: usize) -> u128 {
        self.0.get(i).copied().unwrap_or(0)
    }

    /// Is this element a square in `F_{p^N}`? In characteristic 2 the Frobenius
    /// `x ↦ x²` is a bijection, so *every* element is a square; in odd
    /// characteristic this is Euler's criterion `x^{(q−1)/2} = 1` (with `0` a
    /// square). The square-class is the `H¹` / discriminant datum the odd-char
    /// classifier reads — so this is what lets the invariant theory run over a
    /// genuine extension field, not just a prime field.
    pub fn is_square(&self) -> bool {
        Self::assert_supported_field();
        if self.is_zero() {
            return true;
        }
        if P == 2 {
            return true; // Frobenius is onto in char 2
        }
        // a^{(q−1)/2} == 1
        let mut e = (Self::field_order() - 1) / 2;
        let mut base = *self;
        let mut acc = Self::one();
        while e > 0 {
            if e & 1 == 1 {
                acc = acc.mul(&base);
            }
            base = base.mul(&base);
            e >>= 1;
        }
        acc == Self::one()
    }

    /// The generator `x` (the class of the indeterminate), i.e. `[0, 1, 0, …]`.
    pub fn generator() -> Self {
        Self::assert_supported_field();
        let mut out = [0u128; N];
        if N > 1 {
            out[1] = 1 % P;
        } else if N == 1 {
            // degree-1: the "field" is F_p and x = 0 in it; this is a degenerate case.
            out[0] = 0;
        }
        Fpn(out)
    }

    /// The element with index `code` in `[0, p^N)` (base-`P` digits = coefficients).
    fn from_code(mut code: u128) -> Self {
        Self::assert_supported_field();
        let mut coeffs = [0u128; N];
        for slot in coeffs.iter_mut() {
            *slot = code % P;
            code /= P;
        }
        Fpn(coeffs)
    }

    // ===== The finite-field analysis toolkit =====
    //
    // The shared Galois engine (degree, conjugates, minimal-polynomial product,
    // relative trace/norm, multiplicative order, discrete log) is the
    // `FiniteField` trait below — one algorithm over `Nimber` and `Fpn` both.
    // `Fpn` keeps only the two pieces that are genuinely per-backend: the `F_p`
    // projection of the minimal polynomial, and primitive-element enumeration.

    /// The **minimal polynomial** over `F_p`, as coefficients in `[0, P)` from the
    /// constant term up — monic of degree [`degree`](FiniteField::degree). The
    /// shared `∏ (X − xᵢ)` construction is [`FiniteField::min_poly_monic`]; this
    /// projects each coefficient (Galois-closure guarantees it lies in `F_p`) to
    /// its base-field value.
    pub fn min_poly(&self) -> Vec<u128> {
        Self::assert_supported_field();
        self.min_poly_monic()
            .into_iter()
            .map(|coeff| {
                debug_assert!(
                    coeff.0[1..].iter().all(|&c| c == 0),
                    "minimal-polynomial coefficient left F_p"
                );
                coeff.coeff(0)
            })
            .collect()
    }

    /// A **primitive element** (a generator of `F_{p^N}*`), found by scanning the
    /// field — cheap for the modest orders in this tower.
    pub fn primitive_element() -> Self {
        Self::assert_supported_field();
        let target = Self::field_order() - 1;
        for code in 1..Self::field_order() {
            let el = Self::from_code(code);
            if el.multiplicative_order() == Some(target) {
                return el;
            }
        }
        panic!("Fpn: no primitive element found (unreachable for a field)");
    }
}

/// `Fpn` plugs into the shared [`FiniteField`] engine by supplying only the
/// field shape: the Frobenius `x ↦ x^p`, integer exponentiation, the extension
/// degree `N`, and the multiplicative-group order `p^N − 1` with its factors.
/// Every Galois notion is then a default method. The brute-force discrete log
/// (the trait default) suffices for the small orders here — no Pohlig–Hellman
/// needed, unlike the nimber `F_{2^128}`.
impl<const P: u128, const N: usize> FiniteField for Fpn<P, N> {
    fn frobenius(&self) -> Self {
        Self::assert_supported_field();
        self.pow(P)
    }

    fn pow(&self, mut e: u128) -> Self {
        Self::assert_supported_field();
        let mut base = *self;
        let mut acc = Self::one();
        while e > 0 {
            if e & 1 == 1 {
                acc = acc.mul(&base);
            }
            base = base.mul(&base);
            e >>= 1;
        }
        acc
    }

    fn ext_degree() -> usize {
        Self::assert_supported_field();
        N
    }

    fn group_order() -> u128 {
        Self::assert_supported_field();
        Self::field_order() - 1
    }

    fn group_order_factors() -> Vec<u128> {
        Self::assert_supported_field();
        distinct_primes(Self::field_order() - 1)
    }
}

/// The distinct prime factors of `n` by trial division (small `n = p^N − 1`).
fn distinct_primes(mut n: u128) -> Vec<u128> {
    let mut out = Vec::new();
    let mut d = 2u128;
    while d <= n / d {
        if n.is_multiple_of(d) {
            out.push(d);
            while n.is_multiple_of(d) {
                n /= d;
            }
        }
        d += 1;
    }
    if n > 1 {
        out.push(n);
    }
    out
}

impl<const P: u128, const N: usize> fmt::Display for Fpn<P, N> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut parts: Vec<String> = Vec::new();
        for i in (0..N).rev() {
            let c = self.0[i];
            if c == 0 {
                continue;
            }
            // Display v2 (§9): explicit `⋅` and `↑`, coefficient-1 suppressed.
            let term = match i {
                0 => format!("{c}"),
                1 if c == 1 => "x".to_string(),
                1 => format!("{c}⋅x"),
                _ if c == 1 => format!("x↑{i}"),
                _ => format!("{c}⋅x↑{i}"),
            };
            parts.push(term);
        }
        if parts.is_empty() {
            write!(f, "0")
        } else {
            write!(f, "{}", parts.join(" + "))
        }
    }
}

impl<const P: u128, const N: usize> fmt::Debug for Fpn<P, N> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(self, f)
    }
}

impl<const P: u128, const N: usize> Scalar for Fpn<P, N> {
    fn zero() -> Self {
        Self::assert_supported_field();
        Fpn([0u128; N])
    }

    fn one() -> Self {
        Self::assert_supported_field();
        let mut out = [0u128; N];
        out[0] = 1 % P;
        Fpn(out)
    }

    fn add(&self, rhs: &Self) -> Self {
        Self::assert_supported_field();
        let mut out = [0u128; N];
        for i in 0..N {
            out[i] = add_mod::<P>(self.0[i], rhs.0[i]);
        }
        Fpn(out)
    }

    fn neg(&self) -> Self {
        Self::assert_supported_field();
        let mut out = [0u128; N];
        for i in 0..N {
            out[i] = if self.0[i] == 0 { 0 } else { P - self.0[i] };
        }
        Fpn(out)
    }

    fn mul(&self, rhs: &Self) -> Self {
        Self::assert_supported_field();
        // Schoolbook product into a degree-(2N-2) scratch, then reduce mod m(x).
        let mut scratch = vec![0u128; 2 * N - 1];
        for i in 0..N {
            if self.0[i] == 0 {
                continue;
            }
            let ai = self.0[i];
            for j in 0..N {
                scratch[i + j] = add_mod::<P>(scratch[i + j], mul_mod::<P>(ai, rhs.0[j]));
            }
        }
        // x^k = x^{k-N} · x^N = x^{k-N} · Σ_i red_i x^i, folding top down. (Degree 1 =
        // F_p needs no reduction — the scratch is already a single coefficient.)
        if N > 1 {
            let red = reduction::<P, N>();
            for k in (N..2 * N - 1).rev() {
                let c = scratch[k];
                if c == 0 {
                    continue;
                }
                scratch[k] = 0;
                for i in 0..N {
                    scratch[k - N + i] = add_mod::<P>(scratch[k - N + i], mul_mod::<P>(c, red[i]));
                }
            }
        }
        let mut out = [0u128; N];
        out[..N].copy_from_slice(&scratch[..N]);
        Fpn(out)
    }

    fn characteristic() -> u128 {
        Self::assert_supported_field();
        // The *characteristic* is the prime p, not the order p^N.
        P
    }

    fn inv(&self) -> Option<Self> {
        Self::assert_supported_field();
        if self.is_zero() {
            return None;
        }
        // Fermat: a^{p^N − 2} = a^{−1} in F_{p^N}. Square-and-multiply with `mul`.
        let mut e = Self::field_order() - 2;
        let mut base = *self;
        let mut result = Self::one();
        while e > 0 {
            if e & 1 == 1 {
                result = result.mul(&base);
            }
            base = base.mul(&base);
            e >>= 1;
        }
        Some(result)
    }
    /// Faster direct construction via the constant coefficient; semantically
    /// identical to the default double-and-add (reduction mod p in degree-0).
    fn from_int(n: i128) -> Self {
        Self::assert_supported_field();
        let mut out = [0u128; N];
        if N > 0 {
            out[0] = Fp::<P>::from_int(n).value();
        }
        Fpn(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::clifford::{CliffordAlgebra, Metric};
    use crate::scalar::FiniteField;

    /// Every element of `F_{p^N}`, enumerated by base-`P` digits.
    fn elems<const P: u128, const N: usize>() -> Vec<Fpn<P, N>> {
        let order = Fpn::<P, N>::field_order();
        (0..order)
            .map(|mut code| {
                let mut coeffs = [0u128; N];
                for slot in coeffs.iter_mut() {
                    *slot = code % P;
                    code /= P;
                }
                Fpn::from_coeffs(&coeffs)
            })
            .collect()
    }

    fn check_field_axioms<const P: u128, const N: usize>() {
        let es = elems::<P, N>();
        let zero = Fpn::<P, N>::zero();
        let one = Fpn::<P, N>::one();
        assert_eq!(es.len(), Fpn::<P, N>::field_order() as usize);
        for &a in &es {
            // additive identity / inverse
            assert_eq!(a.add(&zero), a);
            assert_eq!(a.add(&a.neg()), zero);
            // multiplicative identity
            assert_eq!(a.mul(&one), a);
            // inverse: every nonzero element is a unit (THIS is what catches a
            // reducible reduction polynomial — a zero divisor would have no inverse).
            if a.is_zero() {
                assert!(a.inv().is_none());
            } else {
                let ai = a.inv().expect("nonzero element of a field is invertible");
                assert_eq!(a.mul(&ai), one, "a·a⁻¹ = 1");
            }
            for &b in &es {
                assert_eq!(a.add(&b), b.add(&a), "add commutes");
                assert_eq!(a.mul(&b), b.mul(&a), "mul commutes");
                for &c in &es {
                    assert_eq!(a.add(&b).add(&c), a.add(&b.add(&c)), "add assoc");
                    assert_eq!(a.mul(&b).mul(&c), a.mul(&b.mul(&c)), "mul assoc");
                    assert_eq!(a.mul(&b.add(&c)), a.mul(&b).add(&a.mul(&c)), "distrib");
                }
            }
        }
    }

    #[test]
    fn field_axioms_generated_small_fields() {
        check_field_axioms::<2, 2>(); // F_4
        check_field_axioms::<2, 3>(); // F_8
        check_field_axioms::<2, 4>(); // F_16
        check_field_axioms::<2, 5>(); // F_32, generated
        check_field_axioms::<3, 2>(); // F_9
        check_field_axioms::<5, 2>(); // F_25
        check_field_axioms::<3, 3>(); // F_27
    }

    #[test]
    fn generated_rows_match_small_curated_oracles_without_using_them() {
        // These constants are test-only: the production path above always calls the
        // deterministic generator for extension fields. The comparison protects the
        // generator's scan order and keeps the old Conway rows as oracles, not runtime
        // data.
        assert_eq!(Fpn::<2, 2>::reduction_rule(), &[1, 1]);
        assert_eq!(Fpn::<2, 3>::reduction_rule(), &[1, 1, 0]);
        assert_eq!(Fpn::<2, 4>::reduction_rule(), &[1, 1, 0, 0]);
        assert_eq!(Fpn::<3, 2>::reduction_rule(), &[2, 0]);
        assert_eq!(Fpn::<5, 2>::reduction_rule(), &[2, 0]);
        assert_eq!(
            Fpn::<2, 4>::reduction_polynomial_kind(),
            ReductionPolynomialKind::GeneratedIrreducible
        );
        assert!(!Fpn::<2, 2>::is_conway_polynomial());
        assert_eq!(
            Fpn::<7, 1>::reduction_polynomial_kind(),
            ReductionPolynomialKind::PrimeField
        );
    }

    #[test]
    fn generated_metadata_opens_char2_extension_rows() {
        assert!(Fpn::<2, 5>::is_supported_field()); // F_32
        assert!(Fpn::<2, 6>::is_supported_field()); // F_64
        assert!(Fpn::<2, 7>::is_supported_field()); // F_128
        assert_eq!(Fpn::<2, 7>::field_order(), 128);
        assert_eq!(
            Fpn::<2, 5>::reduction_polynomial_kind(),
            ReductionPolynomialKind::GeneratedIrreducible
        );
        assert_eq!(Fpn::<2, 5>::reduction_rule().len(), 5);
        assert!(is_irreducible_monic(
            &reduction_rule_to_polynomial(Fpn::<2, 5>::reduction_rule(), 2),
            2
        ));

        let g = Fpn::<2, 7>::primitive_element();
        assert_eq!(g.multiplicative_order(), Some(127));
        assert!(g.is_primitive());
    }

    #[test]
    fn characteristic_is_p_not_order() {
        assert_eq!(Fpn::<2, 3>::characteristic(), 2); // F_8 has characteristic 2
        assert_eq!(Fpn::<2, 3>::field_order(), 8);
        assert_eq!(Fpn::<3, 3>::characteristic(), 3); // F_27 has characteristic 3
        assert_eq!(Fpn::<3, 3>::field_order(), 27);
    }

    #[test]
    fn unsupported_parameters_are_rejected() {
        assert!(std::panic::catch_unwind(Fpn::<4, 2>::one).is_err());
        assert!(std::panic::catch_unwind(Fpn::<3, 0>::zero).is_err());
        assert!(std::panic::catch_unwind(Fpn::<2, 128>::one).is_err());
    }

    #[test]
    fn from_coeffs_rejects_nonzero_high_terms() {
        assert_eq!(
            Fpn::<2, 3>::from_coeffs(&[1, 0, 1, 0]),
            Fpn::<2, 3>::from_coeffs(&[1, 0, 1])
        );
        assert!(std::panic::catch_unwind(|| Fpn::<2, 3>::from_coeffs(&[1, 0, 0, 1])).is_err());
    }

    #[test]
    fn display_v2_canonical_ogham() {
        // Display v2 (§9): explicit `⋅` and `↑`, coefficient-1 suppressed.
        // The §9 example `3⋅x↑2 + 2⋅x + 1` needs coefficient 3, so it is only
        // realizable in a field whose characteristic exceeds 3 (in F_27 the
        // coefficient 3 reduces to 0). Pin it in F_125.
        let f125 = Fpn::<5, 3>::from_coeffs(&[1, 2, 3]);
        assert_eq!(format!("{f125:?}"), "3⋅x↑2 + 2⋅x + 1");
        // Over F_27 (the menu's `Fpn<3,3>`), pin a realizable element.
        let f27 = Fpn::<3, 3>::from_coeffs(&[1, 1, 2]);
        assert_eq!(format!("{f27:?}"), "2⋅x↑2 + x + 1");
        // Coefficient-1 and bare-`x` suppression: `x↑2`, `x`.
        assert_eq!(format!("{:?}", Fpn::<5, 3>::from_coeffs(&[0, 1, 1])), "x↑2 + x");
        assert_eq!(format!("{:?}", Fpn::<3, 3>::zero()), "0");
    }

    #[test]
    fn generator_satisfies_its_minimal_polynomial() {
        // F_8: x³ = x + 1, so x³ + x + 1 = 0 (and −1 = 1 in char 2 ⇒ x³ = x + 1).
        let x = Fpn::<2, 3>::generator();
        let x3 = x.mul(&x).mul(&x);
        assert_eq!(x3, Fpn::<2, 3>::from_coeffs(&[1, 1, 0])); // x + 1
                                                              // F_16: x⁴ = x + 1.
        let w = Fpn::<2, 4>::generator();
        let w4 = w.mul(&w).mul(&w).mul(&w);
        assert_eq!(w4, Fpn::<2, 4>::from_coeffs(&[1, 1, 0, 0])); // x + 1
                                                                 // F_27: the reduction is generated, not fixed to the old curated row.
        let y = Fpn::<3, 3>::generator();
        let y3 = y.mul(&y).mul(&y);
        assert_eq!(y3, Fpn::<3, 3>::from_coeffs(Fpn::<3, 3>::reduction_rule()));
    }

    #[test]
    fn frobenius_is_an_automorphism() {
        // x ↦ x^p is additive (the Frobenius) in characteristic p.
        let pow_p = |a: Fpn<3, 3>| {
            let mut r = Fpn::<3, 3>::one();
            for _ in 0..3 {
                r = r.mul(&a);
            }
            r
        };
        for a in elems::<3, 3>() {
            for b in elems::<3, 3>() {
                assert_eq!(pow_p(a.add(&b)), pow_p(a).add(&pow_p(b)));
            }
        }
    }

    #[test]
    fn galois_toolkit_f8_f9_f27() {
        // F_8 = F_2[x]/(x³+x+1): the generator has degree 3 and minimal
        // polynomial x³ + x + 1 = [1,1,0,1]; F_8* is cyclic of prime order 7.
        let x = Fpn::<2, 3>::generator();
        assert_eq!(x.degree(), 3);
        assert_eq!(Fpn::<2, 3>::one().degree(), 1);
        assert_eq!(x.conjugates().len(), 3);
        assert_eq!(x.min_poly(), vec![1u128, 1, 0, 1]); // x³ + x + 1
        assert_eq!(x.multiplicative_order(), Some(7));
        assert!(x.is_primitive());
        // primitive element generates the group; discrete log round-trips.
        let g = Fpn::<2, 3>::primitive_element();
        assert_eq!(g.multiplicative_order(), Some(7));
        for e in 0..7u128 {
            assert_eq!(g.discrete_log(g.pow(e)), Some(e % 7));
        }
        // F_16's Conway generator has order 15 for x^4+x+1.
        let c = Fpn::<2, 4>::generator();
        assert_eq!(c.multiplicative_order(), Some(15));
        assert!(c.is_primitive());
        // Absolute trace/norm to F_2 land in the prime field (constant element).
        let tr = x.relative_trace(1);
        let nm = x.relative_norm(1);
        assert!(tr.coeffs()[1..].iter().all(|&c| c == 0), "trace not in F_2");
        assert!(nm.coeffs()[1..].iter().all(|&c| c == 0), "norm not in F_2");
        // F_9: orders divide 8; the primitive element has order exactly 8.
        let h = Fpn::<3, 2>::primitive_element();
        assert_eq!(h.multiplicative_order(), Some(8));
        assert!(h.is_primitive());
        // F_27: the generator has degree 3 and its conjugate orbit closes.
        let z = Fpn::<3, 3>::generator();
        assert_eq!(z.degree(), 3);
        assert_eq!(z.conjugates().len(), 3);
        // every conjugate is a root of the same minimal polynomial.
        let mp = z.min_poly();
        assert_eq!(mp.len(), 4); // monic degree 3
                                 // Frobenius is an automorphism fixing exactly F_p (degree-1 elements).
        assert_eq!(
            Fpn::<3, 3>::constant(2).frobenius(),
            Fpn::<3, 3>::constant(2)
        );
    }

    #[test]
    fn clifford_over_f9_monomorphises() {
        // Cl over F_9 with q = [x, 1]: the engine runs on the extension field exactly
        // as on a prime field; antisymmetry signs are genuine (−1 = 2 in F_3 ⊂ F_9).
        let x = Fpn::<3, 2>::generator();
        let one = Fpn::<3, 2>::one();
        let alg = CliffordAlgebra::new(2, Metric::diagonal(vec![x, one]));
        let (e0, e1) = (alg.e(0), alg.e(1));
        assert_eq!(alg.mul(&e0, &e0), alg.scalar(x));
        assert_eq!(alg.mul(&e1, &e1), alg.scalar(one));
        // e0 e1 = −(e1 e0)
        let neg_one = Fpn::<3, 2>::one().neg();
        assert_eq!(
            alg.mul(&e0, &e1),
            alg.scalar_mul(&neg_one, &alg.mul(&e1, &e0))
        );
    }
}
