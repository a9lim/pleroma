//! Finite-field polynomial factorization shared by the function-field layers.
//!
//! The public place APIs only need the square-free support: the distinct monic
//! irreducible factors of a polynomial. This helper replaces the old monic
//! trial division with the standard finite-field pipeline:
//!
//! 1. square-free support via `gcd(f, f')`, including inseparable `p`-th roots;
//! 2. distinct-degree factorization using `gcd(x^{q^d} - x, f)`;
//! 3. equal-degree splitting by deterministic Cantor-Zassenhaus style tests.

use crate::scalar::{Poly, Scalar};

fn is_one<S: Scalar>(f: &Poly<S>) -> bool {
    *f == Poly::one()
}

fn checked_pow(base: u128, exp: usize) -> u128 {
    let mut acc = 1u128;
    for _ in 0..exp {
        acc = acc
            .checked_mul(base)
            .expect("finite-field polynomial factorization exponent exceeds u128");
    }
    acc
}

fn scalar_pow<S: Scalar + Copy>(mut base: S, mut exp: u128) -> S {
    let mut acc = S::one();
    while exp > 0 {
        if exp & 1 == 1 {
            acc = acc.mul(&base);
        }
        base = base.mul(&base);
        exp >>= 1;
    }
    acc
}

fn div_exact<S: Scalar>(a: &Poly<S>, b: &Poly<S>) -> Poly<S> {
    let (q, r) = a.divrem(b);
    debug_assert!(r.is_zero(), "expected exact polynomial division");
    q
}

fn proper_factor<S: Scalar>(f: &Poly<S>, g: &Poly<S>) -> bool {
    !g.is_zero() && !is_one(g) && g.degree().unwrap_or(0) < f.degree().unwrap_or(0)
}

fn dedup_push<S: Scalar>(out: &mut Vec<Poly<S>>, f: Poly<S>) {
    let f = f.make_monic();
    if !out.contains(&f) {
        out.push(f);
    }
}

fn formal_derivative<S: Scalar + Copy>(f: &Poly<S>, p: u128, from_index: fn(u128) -> S) -> Poly<S> {
    let cs = f.coeffs();
    if cs.len() <= 1 {
        return Poly::zero();
    }
    let mut out = vec![S::zero(); cs.len() - 1];
    for (i, c) in cs.iter().enumerate().skip(1) {
        let factor = (i as u128) % p;
        if factor != 0 {
            out[i - 1] = c.mul(&from_index(factor));
        }
    }
    Poly::new(out)
}

fn pth_root_poly<S: Scalar + Copy>(
    f: &Poly<S>,
    p: u128,
    q: u128,
    _from_index: fn(u128) -> S,
) -> Poly<S> {
    let p_usize = usize::try_from(p).expect("field characteristic exceeds usize");
    let root_exp = q / p;
    let mut out = vec![S::zero(); f.degree().unwrap_or(0) / p_usize + 1];
    for (i, c) in f.coeffs().iter().enumerate() {
        if c.is_zero() {
            continue;
        }
        assert!(
            i % p_usize == 0,
            "zero derivative polynomial has non-p-multiple support"
        );
        out[i / p_usize] = scalar_pow(*c, root_exp);
    }
    Poly::new(out)
}

fn squarefree_parts<S: Scalar + Copy>(
    f: &Poly<S>,
    p: u128,
    q: u128,
    from_index: fn(u128) -> S,
) -> Vec<Poly<S>> {
    if f.degree().unwrap_or(0) == 0 {
        return Vec::new();
    }
    let f = f.make_monic();
    let der = formal_derivative(&f, p, from_index);
    if der.is_zero() {
        return squarefree_parts(&pth_root_poly(&f, p, q, from_index), p, q, from_index);
    }

    let mut c = f.gcd(&der);
    let mut w = div_exact(&f, &c);
    let mut parts = Vec::new();
    while !is_one(&w) {
        let y = w.gcd(&c);
        let z = div_exact(&w, &y);
        if !is_one(&z) {
            parts.push(z.make_monic());
        }
        w = y.clone();
        c = div_exact(&c, &y);
    }
    if !is_one(&c) {
        parts.extend(squarefree_parts(
            &pth_root_poly(&c, p, q, from_index),
            p,
            q,
            from_index,
        ));
    }
    parts
}

fn poly_from_index<S: Scalar + Copy>(
    mut idx: u128,
    len: usize,
    q: u128,
    from_index: fn(u128) -> S,
) -> Poly<S> {
    let mut coeffs = Vec::with_capacity(len);
    for _ in 0..len {
        coeffs.push(from_index(idx % q));
        idx /= q;
    }
    Poly::new(coeffs)
}

fn split_equal_degree<S: Scalar + Copy>(
    f: &Poly<S>,
    d: usize,
    p: u128,
    q: u128,
    from_index: fn(u128) -> S,
) -> Vec<Poly<S>> {
    let n = f
        .degree()
        .expect("equal-degree factorization needs nonzero input");
    if n == d {
        return vec![f.make_monic()];
    }

    let seed_count = checked_pow(q, n);
    for seed in 0..seed_count {
        let h = poly_from_index(seed, n, q, from_index).rem(f);
        if h.degree().unwrap_or(0) == 0 {
            continue;
        }
        let early = f.gcd(&h);
        if proper_factor(f, &early) {
            let rest = div_exact(f, &early);
            let mut out = split_equal_degree(&early, d, p, q, from_index);
            out.extend(split_equal_degree(&rest, d, p, q, from_index));
            return out;
        }

        let splitter = if p == 2 {
            let mut trace = Poly::zero();
            let mut hp = h;
            for _ in 0..d {
                trace = trace.add(&hp).rem(f);
                hp = hp.pow_mod(q, f);
            }
            trace
        } else {
            let exp = (checked_pow(q, d) - 1) / 2;
            h.pow_mod(exp, f).sub(&Poly::one()).rem(f)
        };
        let g = f.gcd(&splitter);
        if proper_factor(f, &g) {
            let rest = div_exact(f, &g);
            let mut out = split_equal_degree(&g, d, p, q, from_index);
            out.extend(split_equal_degree(&rest, d, p, q, from_index));
            return out;
        }
    }

    panic!("deterministic equal-degree factorization failed to split a reducible factor");
}

fn distinct_degree_factors<S: Scalar + Copy>(
    f: &Poly<S>,
    p: u128,
    q: u128,
    from_index: fn(u128) -> S,
) -> Vec<Poly<S>> {
    let mut out = Vec::new();
    let mut g = f.make_monic();
    let x = Poly::<S>::x();
    let mut h = x.clone();
    let mut d = 1usize;
    while g.degree().unwrap_or(0) >= 2 * d {
        h = h.pow_mod(q, &g);
        let factor = g.gcd(&h.sub(&x).rem(&g));
        if !is_one(&factor) {
            for pi in split_equal_degree(&factor, d, p, q, from_index) {
                dedup_push(&mut out, pi);
            }
            g = div_exact(&g, &factor).make_monic();
            if is_one(&g) {
                return out;
            }
            h = h.rem(&g);
        }
        d += 1;
    }
    if !is_one(&g) {
        dedup_push(&mut out, g);
    }
    out
}

/// The distinct monic irreducible factors of `f` over a finite field with
/// characteristic `p`, order `q`, and deterministic element enumeration.
pub(crate) fn monic_irreducible_factor_support<S: Scalar + Copy>(
    f: &Poly<S>,
    p: u128,
    q: u128,
    from_index: fn(u128) -> S,
) -> Vec<Poly<S>> {
    assert!(
        p >= 2 && q.is_multiple_of(p),
        "invalid finite-field metadata"
    );
    let mut out = Vec::new();
    for part in squarefree_parts(f, p, q, from_index) {
        for pi in distinct_degree_factors(&part, p, q, from_index) {
            dedup_push(&mut out, pi);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forms::{FiniteChar2Field, FiniteOddField};
    use crate::scalar::{Fp, Fpn};

    fn factor_odd<S: FiniteOddField>(f: &Poly<S>) -> Vec<Poly<S>> {
        monic_irreducible_factor_support(
            f,
            S::characteristic_prime(),
            S::field_order(),
            S::from_index,
        )
    }

    fn factor_char2<S: FiniteChar2Field>(f: &Poly<S>) -> Vec<Poly<S>> {
        monic_irreducible_factor_support(
            f,
            S::characteristic_prime(),
            S::field_order(),
            S::from_index,
        )
    }

    #[test]
    fn odd_factorization_splits_equal_degree_products() {
        type F = Fp<5>;
        let x2_minus_1 = Poly::new(vec![F::new(-1), F::zero(), F::one()]);
        let fs = factor_odd(&x2_minus_1);
        assert_eq!(fs.len(), 2);
        assert!(fs.contains(&Poly::new(vec![F::new(-1), F::one()])));
        assert!(fs.contains(&Poly::new(vec![F::one(), F::one()])));
    }

    #[test]
    fn repeated_and_inseparable_support_is_deduped() {
        type F2 = Fp<2>;
        let t_plus_1 = Poly::new(vec![F2::one(), F2::one()]);
        let fourth_power = t_plus_1.mul(&t_plus_1).mul(&t_plus_1).mul(&t_plus_1);
        assert_eq!(factor_char2(&fourth_power), vec![t_plus_1]);
    }

    #[test]
    fn extension_field_coefficients_factor_too() {
        type F9 = Fpn<3, 2>;
        let a = F9::generator();
        let f = Poly::new(vec![a.neg(), F9::one()]); // x - a
        assert_eq!(factor_odd(&f), vec![f]);
    }
}
