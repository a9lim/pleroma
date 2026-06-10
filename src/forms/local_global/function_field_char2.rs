//! The **characteristic-2** local–global symbol over the global function field
//! `F_{2^m}(t)` — the equal-characteristic-2 mirror of
//! [`forms::function_field`](crate::forms::function_field) (which is the odd-`q`
//! mirror of [`forms::padic`](crate::forms::padic)).
//!
//! In odd characteristic the quaternion/quadratic-form symbol is the **symmetric**
//! tame Hilbert symbol `(a,b)_v` built from the multiplicative square class. In
//! characteristic 2 it is the **asymmetric Artin–Schreier symbol** `[a, b)`: the
//! first slot `a ∈ F_q(t)` is *additive* (taken mod `℘(K)`, `℘(x)=x²+x`), the second
//! `b ∈ F_q(t)*` is *multiplicative*. The algebra is the cyclic/quaternion algebra
//! `[a,b) = ⟨i,j | i²+i=a, j²=b, jij⁻¹=i+1⟩`; its reduced-norm form is the 2-fold
//! quadratic Pfister form `[1,a] ⊥ b·[1,a]` (where `[1,a]` is the binary form
//! `x²+xy+ay²`), and that form is isotropic/hyperbolic **iff** `[a,b)` splits.
//!
//! # The local symbol — a residue of a differential
//!
//! At a place `v` (residue field `κ`), the local invariant is
//!
//! ```text
//! s_v(a, b) = Tr_{κ/F₂}( Res_v( a · db/b ) ) ∈ F₂,
//! ```
//!
//! the **Schmid formula** (Serre, *Local Fields* XIV §5; Gille–Szamuely §9.2). The
//! algebra `[a,b)` splits at `v` iff `s_v = 0`. The differential `a·db/b = a·dlog b`
//! is a rational differential; its residue is the coefficient of `ϖ⁻¹` in the local
//! Laurent expansion. Globally, **Weil reciprocity** (the residue theorem on `P¹`)
//! gives `Σ_v s_v(a,b) = 0` — so the algebra ramifies (symbol `1`) at an **even**
//! number of places, the exact char-2 analogue of `∏_v (a,b)_v = +1`.
//!
//! # Computing the residue (the Hensel-parametrization formula)
//!
//! At a finite place `P` (monic irreducible, `κ = F_q[t]/(P)`), write the local
//! parameter `u = P(t)`; there is a unique Hensel series `T(u) ∈ κ[[u]]` with
//! `T(0) = t mod P` and `P(T(u)) = u`, and then `dt = P'(T(u))⁻¹ du`. For
//! `g = N/D` with `D = P^m·E`, `(E,P)=1`, set `B ≡ N·E⁻¹ (mod P^m)`; then
//!
//! ```text
//! Res_P(g dt) = [u^{m-1}]( B(T(u)) · P'(T(u))⁻¹ ) ∈ κ.
//! ```
//!
//! (Even-order poles contribute through the higher `u`-coefficients — there is no
//! "only the simple pole matters" shortcut in char 2. The odd-order tail that
//! Hermite reduction *cannot* remove is the same wild `R_π` phenomenon that
//! Aravire–Jacob's Witt decomposition carries; see root AGENTS.md.) At the degree place
//! `∞` (`κ = F_q`) the substitution `u = 1/t`, `dt = u⁻²du` (the char-2 sign
//! vanishes) gives `Res_∞(g dt) = [u⁻¹]( g(1/u)·u⁻² )`.
//!
//! Scope: this layer is the symbol + reciprocity + quaternion-ramification package.
//! The full char-2 Witt/Springer decomposition of an arbitrary form (the wild
//! `R_π` term) is a separate, larger build tracked in root AGENTS.md.

use crate::forms::{artin_schreier_class_finite, FiniteChar2Field};
use crate::scalar::{Poly, RationalFunction, Scalar};

/// A place of `F_q(t)` in characteristic 2: the degree place `∞` (residue field
/// `F_q`), or a finite place given by a monic irreducible `π(t)` (residue field
/// `F_q[t]/(π)`). The char-2 mirror of
/// [`FFPlace`](crate::forms::function_field::FFPlace).
#[derive(Debug, Clone, PartialEq)]
pub enum Char2Place<S: FiniteChar2Field> {
    /// The degree place `∞` (uniformizer `1/t`, residue field `F_q`).
    Infinite,
    /// A finite place: a monic irreducible `π(t)` (residue field `F_q[t]/(π)`).
    Finite(Poly<S>),
}

// ───────────────────────── polynomial helpers ─────────────────────────

/// The char-2 formal derivative `d/dt`: `Σ cᵢ tⁱ ↦ Σ_{i odd} cᵢ t^{i-1}` (in char 2
/// `i·cᵢ = cᵢ` for odd `i` and `0` for even `i`).
fn dpoly<S: Scalar>(p: &Poly<S>) -> Poly<S> {
    let cs = p.coeffs();
    if cs.len() <= 1 {
        return Poly::zero();
    }
    let mut out = vec![S::zero(); cs.len() - 1];
    for (i, c) in cs.iter().enumerate().skip(1) {
        if i & 1 == 1 {
            out[i - 1] = c.clone();
        }
    }
    Poly::new(out)
}

/// The multiplicity of `pi` in `p` and the cofactor `p / pi^mult`.
pub(crate) fn strip_factor<S: Scalar>(mut p: Poly<S>, pi: &Poly<S>) -> (usize, Poly<S>) {
    let mut mult = 0usize;
    if p.is_zero() {
        return (0, p);
    }
    loop {
        let (quot, rem) = p.divrem(pi);
        if rem.is_zero() {
            p = quot;
            mult += 1;
        } else {
            break;
        }
    }
    (mult, p)
}

/// Extended Euclid on polynomials: returns `(g, x, y)` with `x·a + y·b = g = gcd`.
fn egcd<S: Scalar>(a: &Poly<S>, b: &Poly<S>) -> (Poly<S>, Poly<S>, Poly<S>) {
    let (mut old_r, mut r) = (a.clone(), b.clone());
    let (mut old_s, mut s) = (Poly::one(), Poly::zero());
    let (mut old_t, mut t) = (Poly::zero(), Poly::one());
    while !r.is_zero() {
        let (q, _) = old_r.divrem(&r);
        let nr = old_r.sub(&q.mul(&r));
        old_r = std::mem::replace(&mut r, nr);
        let ns = old_s.sub(&q.mul(&s));
        old_s = std::mem::replace(&mut s, ns);
        let nt = old_t.sub(&q.mul(&t));
        old_t = std::mem::replace(&mut t, nt);
    }
    (old_r, old_s, old_t)
}

/// `e⁻¹ mod m`, for `gcd(e, m) = 1` (so the gcd is a nonzero constant unit).
pub(crate) fn inverse_mod<S: Scalar>(e: &Poly<S>, m: &Poly<S>) -> Poly<S> {
    let (g, x, _) = egcd(e, m);
    let unit = g
        .coeff(0)
        .inv()
        .expect("inverse_mod needs gcd(e, m) = 1 (a unit)");
    x.scale(&unit).rem(m)
}

// ───────────────────────── factorization over F_q ─────────────────────────

/// The distinct monic irreducible factors of `f` over `F_q` (square-free support).
/// The char-2 twin of [`function_field::monic_irreducible_factors`]; the distinct
/// name keeps the flat `forms::*` glob re-export unambiguous.
pub(crate) fn char2_monic_irreducible_factors<S: FiniteChar2Field>(f: &Poly<S>) -> Vec<Poly<S>> {
    crate::forms::poly_factor::monic_irreducible_factor_support(
        f,
        S::characteristic_prime(),
        S::field_order(),
        S::from_index,
    )
}

// ───────────────── κ-power-series arithmetic (mod uᵖʳᵉᶜ over κ=F_q[t]/(P)) ─────────────────
//
// A series is a `Vec<Poly<S>>` of length `prec`; entry `k` is the κ-coefficient (a
// `Poly<S>` reduced mod `P`) of `uᵏ`.

/// `a · b` truncated at `u^prec`, coefficients multiplied in `κ`.
fn ps_mul<S: Scalar>(a: &[Poly<S>], b: &[Poly<S>], prec: usize, p: &Poly<S>) -> Vec<Poly<S>> {
    let mut out = vec![Poly::<S>::zero(); prec];
    for (i, ai) in a.iter().enumerate().take(prec) {
        if ai.is_zero() {
            continue;
        }
        for (j, bj) in b.iter().enumerate().take(prec - i) {
            out[i + j] = out[i + j].add(&ai.mul_mod(bj, p));
        }
    }
    out
}

/// `a⁻¹` truncated at `u^prec` (requires `a[0]` invertible in `κ`).
fn ps_inv<S: FiniteChar2Field>(a: &[Poly<S>], prec: usize, p: &Poly<S>) -> Vec<Poly<S>> {
    let mut b = vec![Poly::<S>::zero(); prec];
    let a0_inv = kappa_inv(&a[0], p);
    b[0] = a0_inv.clone();
    for k in 1..prec {
        let mut acc = Poly::<S>::zero();
        for i in 1..=k {
            let ai = a.get(i).cloned().unwrap_or_else(Poly::zero);
            if !ai.is_zero() {
                acc = acc.add(&ai.mul_mod(&b[k - i], p));
            }
        }
        // char 2: bₖ = −a₀⁻¹·acc = a₀⁻¹·acc.
        b[k] = acc.mul_mod(&a0_inv, p);
    }
    b
}

/// Evaluate the `S`-polynomial `poly` at the κ-series `t` (Horner), truncated at
/// `u^prec`.
pub(crate) fn ps_eval_poly<S: Scalar>(
    poly: &Poly<S>,
    t: &[Poly<S>],
    prec: usize,
    p: &Poly<S>,
) -> Vec<Poly<S>> {
    let mut acc = vec![Poly::<S>::zero(); prec];
    for c in poly.coeffs().iter().rev() {
        acc = ps_mul(&acc, t, prec, p);
        acc[0] = acc[0].add(&Poly::constant(c.clone())).rem(p);
    }
    acc
}

/// `z⁻¹` in `κ = F_q[t]/(P)` by Fermat: `z^{|κ|−2}` with `|κ| = q^{deg P}`.
fn kappa_inv<S: FiniteChar2Field>(z: &Poly<S>, p: &Poly<S>) -> Poly<S> {
    let d = p.degree().expect("a place modulus has degree ≥ 1") as u128;
    let order = S::field_order().pow(
        d.try_into()
            .expect("place degree fits the platform exponent type"),
    );
    z.pow_mod(order - 2, p)
}

/// The Hensel series `T(u) ∈ κ[[u]]` with `T(0) = t mod P`, `P(T(u)) = u`, truncated
/// at `u^prec`. Solved by linear lifting: `[u^k] P(T) = [k=1]` determines `t_k`.
pub(crate) fn hensel_series<S: FiniteChar2Field>(p: &Poly<S>, prec: usize) -> Vec<Poly<S>> {
    let alpha = Poly::monomial(1, S::one()).rem(p); // t mod P
    let mut t = vec![Poly::<S>::zero(); prec];
    if prec == 0 {
        return t;
    }
    t[0] = alpha;
    let pp = dpoly(p);
    let pa_inv = kappa_inv(&pp.rem(p), p); // P'(α)⁻¹ ∈ κ*
    for k in 1..prec {
        // [u^k] of P(T) with t_k currently 0; true value is ck + P'(α)·t_k.
        let pt = ps_eval_poly(p, &t, k + 1, p);
        let target = if k == 1 { Poly::one() } else { Poly::zero() };
        let ck = pt[k].clone();
        // char 2: t_k = (target − ck)·P'(α)⁻¹ = (target + ck)·P'(α)⁻¹.
        t[k] = target.add(&ck).mul_mod(&pa_inv, p);
    }
    t
}

// ───────────────────────── residues ─────────────────────────

/// `Res_P(g dt) ∈ κ = F_q[t]/(P)` for `g = num/den`, `P` monic irreducible.
fn residue_finite<S: FiniteChar2Field>(num: &Poly<S>, den: &Poly<S>, p: &Poly<S>) -> Poly<S> {
    let (m, e) = strip_factor(den.clone(), p);
    if m == 0 {
        return Poly::zero(); // g is regular at P — no residue
    }
    let mut pm = Poly::one();
    for _ in 0..m {
        pm = pm.mul(p);
    }
    let e_inv = inverse_mod(&e, &pm);
    let b = num.mul(&e_inv).rem(&pm); // N·E⁻¹ mod Pᵐ
    let t = hensel_series(p, m);
    let bt = ps_eval_poly(&b, &t, m, p);
    let ppt = ps_eval_poly(&dpoly(p), &t, m, p);
    let ppt_inv = ps_inv(&ppt, m, p);
    let val = ps_mul(&bt, &ppt_inv, m, p);
    val[m - 1].clone()
}

/// `Res_∞(g dt) ∈ F_q` for `g = num/den`. Substitute `u = 1/t`, `dt = u⁻²du`
/// (char 2): the residue is `[u^{k}]( Ñ/D̃ )` with `Ñ, D̃` the coefficient-reversed
/// polynomials and `k = deg N − deg D + 1` (`0` if `k < 0`).
fn residue_infinity<S: FiniteChar2Field>(num: &Poly<S>, den: &Poly<S>) -> S {
    if num.is_zero() {
        return S::zero();
    }
    let dn = num.degree().expect("nonzero numerator") as i128;
    let dd = den.degree().expect("nonzero denominator") as i128;
    let k = dn - dd + 1;
    if k < 0 {
        return S::zero();
    }
    let k = k as usize;
    let rev = |p: &Poly<S>| {
        let mut c = p.coeffs().to_vec();
        c.reverse();
        Poly::new(c)
    };
    let (nt, dt) = (rev(num), rev(den)); // Ñ, D̃ (D̃[0] = leading coeff of D ≠ 0)
                                         // Power series over F_q: Ñ · D̃⁻¹ mod u^{k+1}, then the u^k coefficient.
    let prec = k + 1;
    let d0_inv = dt.coeff(0).inv().expect("D̃(0) = lead(D) inverts");
    let mut binv = vec![S::zero(); prec];
    binv[0] = d0_inv;
    for i in 1..prec {
        let mut acc = S::zero();
        for j in 1..=i {
            acc = acc.add(&dt.coeff(j).mul(&binv[i - j]));
        }
        binv[i] = acc.mul(&d0_inv); // char 2: −d0⁻¹·acc = d0⁻¹·acc
    }
    let mut res = S::zero();
    for i in 0..=k {
        res = res.add(&nt.coeff(i).mul(&binv[k - i]));
    }
    res
}

/// `Tr_{κ/F₂}(z)` for `z ∈ κ = F_q[t]/(P)`: the relative trace `κ → F_q`
/// (`Σ_{i<deg P} z^{q^i}`, a constant) composed with `Tr_{F_q/F₂}` (the
/// Artin–Schreier class).
pub(crate) fn trace_kappa_to_f2<S: FiniteChar2Field>(z: &Poly<S>, p: &Poly<S>) -> u128 {
    let d = p.degree().expect("a place modulus has degree ≥ 1");
    let q = S::field_order();
    let mut term = z.rem(p);
    let mut tr = term.clone();
    for _ in 1..d {
        term = term.pow_mod(q, p);
        tr = tr.add(&term);
    }
    // Tr_{κ/F_q}(z) ∈ F_q: the constant term of the reduced sum.
    artin_schreier_class_finite(tr.rem(p).coeff(0))
}

// ───────────────────────── the symbol ─────────────────────────

/// The reduced `g = a · dlog b = num/den ∈ F_q(t)` (lowest terms) of the differential
/// `a·db/b`. `None` iff `a = 0` (then every symbol is `0`).
fn dlog_differential<S: FiniteChar2Field>(
    a: &RationalFunction<S>,
    b: &RationalFunction<S>,
) -> Option<(Poly<S>, Poly<S>)> {
    assert!(!b.is_zero(), "the Artin–Schreier symbol needs b ≠ 0");
    if a.is_zero() {
        return None;
    }
    let (an, ad) = (a.num(), a.den());
    let (bn, bd) = (b.num(), b.den());
    // dlog b = (Bn'·Bd + Bn·Bd') / (Bn·Bd)   (char 2: the cross signs vanish).
    let dlog_num = dpoly(bn).mul(bd).add(&bn.mul(&dpoly(bd)));
    // g = a · dlog b = An·dlog_num / (Ad·Bn·Bd).
    let mut gnum = an.mul(&dlog_num);
    let mut gden = ad.mul(bn).mul(bd);
    let gg = gnum.gcd(&gden);
    if gg.degree().unwrap_or(0) > 0 {
        gnum = gnum.divrem(&gg).0;
        gden = gden.divrem(&gg).0;
    }
    Some((gnum, gden))
}

/// The Artin–Schreier symbol `s_v(a, b) ∈ {0, 1}` at `place`, for `b ≠ 0`. `0` iff
/// the cyclic algebra `[a, b)` splits over the completion at `place`. The char-2
/// mirror of
/// [`try_hilbert_symbol_ff`](crate::forms::function_field::try_hilbert_symbol_ff).
pub fn as_symbol_at<S: FiniteChar2Field>(
    a: &RationalFunction<S>,
    b: &RationalFunction<S>,
    place: &Char2Place<S>,
) -> u128 {
    let Some((gnum, gden)) = dlog_differential(a, b) else {
        return 0; // a ∈ ℘(K) (in particular a = 0): the symbol vanishes
    };
    match place {
        Char2Place::Finite(pi) => trace_kappa_to_f2(&residue_finite(&gnum, &gden, pi), pi),
        Char2Place::Infinite => artin_schreier_class_finite(residue_infinity(&gnum, &gden)),
    }
}

/// The places that can carry a nontrivial symbol for `[a, b)` (`b ≠ 0`): the poles
/// of `a·dlog b` (monic irreducible factors of its reduced denominator) plus the
/// degree place `∞`. Every other place sees a regular differential, residue `0`.
pub fn as_symbol_places<S: FiniteChar2Field>(
    a: &RationalFunction<S>,
    b: &RationalFunction<S>,
) -> Vec<Char2Place<S>> {
    let mut places = vec![Char2Place::Infinite];
    if let Some((_, gden)) = dlog_differential(a, b) {
        for pi in char2_monic_irreducible_factors(&gden) {
            places.push(Char2Place::Finite(pi));
        }
    }
    places
}

/// The **Weil reciprocity sum** `Σ_v s_v(a, b) ∈ F₂` over all places — identically
/// `0` for every `a` and `b ≠ 0` (the residue theorem on `P¹`). The char-2 additive
/// analogue of the odd-char product formula `∏_v (a,b)_v = +1`, the gold oracle.
pub fn as_symbol_reciprocity_sum<S: FiniteChar2Field>(
    a: &RationalFunction<S>,
    b: &RationalFunction<S>,
) -> u128 {
    as_symbol_places(a, b)
        .iter()
        .fold(0u128, |acc, pl| acc ^ as_symbol_at(a, b, pl))
}

/// The places where the cyclic algebra `[a, b)` **ramifies** (symbol `1`), `b ≠ 0`.
/// The count is always **even** (additive reciprocity), mirroring
/// [`try_ramified_places_ff`](crate::forms::function_field::try_ramified_places_ff).
pub fn as_symbol_ramified_places<S: FiniteChar2Field>(
    a: &RationalFunction<S>,
    b: &RationalFunction<S>,
) -> Vec<Char2Place<S>> {
    as_symbol_places(a, b)
        .into_iter()
        .filter(|pl| as_symbol_at(a, b, pl) == 1)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scalar::{Fp, Fpn};

    type F2 = Fp<2>;
    type R2 = RationalFunction<F2>;

    fn p2(c: &[i128]) -> Poly<F2> {
        Poly::new(c.iter().map(|&n| F2::new(n)).collect())
    }
    fn r2(num: &[i128], den: &[i128]) -> R2 {
        RationalFunction::new(
            num.iter().map(|&n| F2::new(n)).collect(),
            den.iter().map(|&n| F2::new(n)).collect(),
        )
    }

    // ── the residue engine, against Codex's source-derived oracles ──
    // P = t²+t+1 over F₂ (irreducible, κ = F₄ = F₂(α), α²+α+1=0; α = "t mod P").

    #[test]
    fn residue_oracles_at_a_degree_two_place() {
        let p = p2(&[1, 1, 1]); // t² + t + 1
        let p2sq = p.mul(&p); // P²
        let p3 = p2sq.mul(&p); // P³
        let alpha = p2(&[0, 1]); // t mod P = α ∈ F₄
        let one = Poly::<F2>::one();

        // Res_P(t/P² dt) = 1   (Hermite: t/P² dt = d(t/P) + 1/P dt)
        assert_eq!(residue_finite(&p2(&[0, 1]), &p2sq, &p), one);
        // Res_P(1/P³ dt) = 0   (even u-coefficient vanishes)
        assert_eq!(residue_finite(&one, &p3, &p), Poly::zero());
        // Res_P(t/P³ dt) = 1 ∈ F₄   (odd j≥3 term is NOT residue-invisible)
        assert_eq!(residue_finite(&p2(&[0, 1]), &p3, &p), one);
        // A simple pole reads the residue field directly: Res_P((α·t)/P dt)=α·(P')⁻¹·t? —
        // sanity: Res_P(1/P dt) = (P')⁻¹ mod P = 1 (P'=1).
        assert_eq!(residue_finite(&one, &p, &p), one);
        // The Hensel parameter really satisfies P(T)=u: [u¹]T = 1, [u²]T = 1.
        let t = hensel_series(&p, 3);
        assert_eq!(t[0], alpha);
        assert_eq!(t[1], one);
        assert_eq!(t[2], one);
    }

    // ── the symbol + reciprocity, against Codex's worked oracles over F₂(t) ──

    #[test]
    fn symbol_oracle_a1_b_t() {
        // a = 1, b = t:  ω = dt/t.  s = 1 at t=0 and at ∞, 0 elsewhere; Σ = 0.
        let (a, b) = (r2(&[1], &[1]), r2(&[0, 1], &[1]));
        assert_eq!(as_symbol_at(&a, &b, &Char2Place::Finite(p2(&[0, 1]))), 1); // t=0
        assert_eq!(as_symbol_at(&a, &b, &Char2Place::Infinite), 1); // ∞
        assert_eq!(as_symbol_at(&a, &b, &Char2Place::Finite(p2(&[1, 1]))), 0); // t+1: regular
        assert_eq!(as_symbol_reciprocity_sum(&a, &b), 0);
        assert_eq!(as_symbol_ramified_places(&a, &b).len(), 2);
    }

    #[test]
    fn symbol_oracle_a_recip_tp1_b_t() {
        // a = 1/(t+1), b = t:  ω = dt/(t(t+1)).  s = 1 at t=0 and t+1=0, 0 at ∞; Σ = 0.
        let (a, b) = (r2(&[1], &[1, 1]), r2(&[0, 1], &[1]));
        assert_eq!(as_symbol_at(&a, &b, &Char2Place::Finite(p2(&[0, 1]))), 1); // t
        assert_eq!(as_symbol_at(&a, &b, &Char2Place::Finite(p2(&[1, 1]))), 1); // t+1
        assert_eq!(as_symbol_at(&a, &b, &Char2Place::Infinite), 0); // ∞
        assert_eq!(as_symbol_reciprocity_sum(&a, &b), 0);
    }

    #[test]
    fn symbol_oracle_a_recip_irreducible_b_t() {
        // a = 1/(t²+t+1), b = t:  ω = dt/(t·P), P=t²+t+1.  s = 1 at t=0 and at the
        // degree-2 place P (Tr_{F₄/F₂}(t+1)=1), 0 at ∞; Σ = 0.
        let (a, b) = (r2(&[1], &[1, 1, 1]), r2(&[0, 1], &[1]));
        assert_eq!(as_symbol_at(&a, &b, &Char2Place::Finite(p2(&[0, 1]))), 1); // t=0
        assert_eq!(as_symbol_at(&a, &b, &Char2Place::Finite(p2(&[1, 1, 1]))), 1); // P (deg 2)
        assert_eq!(as_symbol_at(&a, &b, &Char2Place::Infinite), 0); // ∞
        assert_eq!(as_symbol_reciprocity_sum(&a, &b), 0);
    }

    #[test]
    fn symbol_oracle_over_f4() {
        // Over F₄(t): a = α (the F₄ generator), b = t. s = 1 at t=0 and ∞ (residue α,
        // Tr_{F₄/F₂}(α)=1), Σ = 0. (And [1,t) would be split: Tr_{F₄/F₂}(1)=0.)
        type F4 = Fpn<2, 2>;
        type R4 = RationalFunction<F4>;
        let alpha = F4::from_coeffs(&[0, 1]);
        let a = R4::from_base(alpha);
        let b = R4::new(
            vec![F4::constant(0), F4::constant(1)],
            vec![F4::constant(1)],
        ); // t
        assert_eq!(
            as_symbol_at(
                &a,
                &b,
                &Char2Place::Finite(Poly::new(vec![F4::constant(0), F4::constant(1)]))
            ),
            1
        );
        assert_eq!(as_symbol_at(&a, &b, &Char2Place::Infinite), 1);
        assert_eq!(as_symbol_reciprocity_sum(&a, &b), 0);
        // [1, t): a = 1 has trace 0 over F₄, so the symbol is 0 at every place.
        let one = R4::from_base(F4::constant(1));
        assert_eq!(as_symbol_reciprocity_sum(&one, &b), 0);
        assert!(as_symbol_ramified_places(&one, &b).is_empty());
    }

    // ── reciprocity sweep: the gold oracle, Σ_v s_v = 0 for every (a, b≠0) ──

    #[test]
    fn reciprocity_sweep_over_f2() {
        let samples = [
            r2(&[0, 1], &[1]),       // t
            r2(&[1, 1], &[1]),       // t+1
            r2(&[1, 0, 1, 1], &[1]), // t³+t+1 (irreducible over F₂)
            r2(&[0, 1], &[1, 1]),    // t/(t+1)
            r2(&[1], &[0, 1]),       // 1/t
            r2(&[1, 1], &[0, 1, 1]), // (t+1)/(t²+t)
        ];
        for a in &samples {
            for b in &samples {
                assert_eq!(
                    as_symbol_reciprocity_sum(a, b),
                    0,
                    "reciprocity Σ_v s_v(a,b) = 0 failed at a={a:?} b={b:?}"
                );
                assert_eq!(
                    as_symbol_ramified_places(a, b).len() % 2,
                    0,
                    "ramified-place count must be even at a={a:?} b={b:?}"
                );
            }
        }
    }

    #[test]
    fn reciprocity_sweep_over_f4() {
        type F4 = Fpn<2, 2>;
        type R4 = RationalFunction<F4>;
        let c = |n: u128| F4::from_index(n);
        let rf = |num: Vec<u128>, den: Vec<u128>| -> R4 {
            RationalFunction::new(
                num.into_iter().map(c).collect(),
                den.into_iter().map(c).collect(),
            )
        };
        let samples = [
            rf(vec![0, 1], vec![1]),    // t
            rf(vec![1, 1], vec![1]),    // t+1
            rf(vec![2, 1], vec![1]),    // t+α
            rf(vec![0, 1], vec![1, 1]), // t/(t+1)
            rf(vec![2], vec![0, 1]),    // α/t
        ];
        for a in &samples {
            for b in &samples {
                assert_eq!(
                    as_symbol_reciprocity_sum(a, b),
                    0,
                    "reciprocity at a={a:?} b={b:?}"
                );
            }
        }
    }

    // ── the defining Steinberg-type relations of the symbol ──

    #[test]
    fn symbol_relations() {
        let places = [
            Char2Place::Infinite,
            Char2Place::Finite(p2(&[0, 1])),    // t
            Char2Place::Finite(p2(&[1, 1])),    // t+1
            Char2Place::Finite(p2(&[1, 1, 1])), // t²+t+1
        ];
        let samples = [r2(&[0, 1], &[1]), r2(&[1, 1], &[1]), r2(&[1], &[0, 1])];
        for a in &samples {
            for b in &samples {
                for pl in &places {
                    // s(a, b²) = 0: dlog(b²) = 0 in char 2.
                    let b2 = b.mul(b);
                    assert_eq!(as_symbol_at(a, &b2, pl), 0, "s(a, b²) = 0");
                    // s(a, a) = 0 (a ≠ 0): a·dlog a = da is exact.
                    assert_eq!(as_symbol_at(a, a, pl), 0, "s(a, a) = 0");
                    // s(℘(x), b) = 0: ℘(x) = x²+x ∈ ℘(K).
                    let wp = a.mul(a).add(a); // x²+x
                    assert_eq!(as_symbol_at(&wp, b, pl), 0, "s(℘(x), b) = 0");
                }
                // additive in the first slot: s(a₁+a₂, b) = s(a₁,b) ⊕ s(a₂,b).
                for pl in &places {
                    let lhs = as_symbol_at(&samples[0].add(&samples[1]), b, pl);
                    let rhs = as_symbol_at(&samples[0], b, pl) ^ as_symbol_at(&samples[1], b, pl);
                    assert_eq!(lhs, rhs, "additive in a");
                }
            }
        }
    }
}
