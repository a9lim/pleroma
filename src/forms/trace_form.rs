//! The **twisted trace form** — the bridge from the "grow a field"
//! [`CyclicGaloisExtension`] layer to the "classify a form" trichotomy.
//!
//! Every cyclic Galois extension `E/F` carries a canonical quadratic form on `E`
//! seen as an `F`-vector space. The *naive* trace form `Tr_{E/F}(x²)` is a **trap**
//! in characteristic 2: Frobenius is additive, so
//! `Tr((x+y)²) = Tr(x²) + Tr(y²)` and the polar form **vanishes** — it degenerates
//! exactly where this project lives. The Arf-bearing object is the
//! **Frobenius-twisted** form
//!
//! ```text
//! Q_k(x) = Tr_{E/F}( x · σ^k(x) )
//! ```
//!
//! with polar `B(x,y) = Tr(x σ^k(y) + y σ^k(x))` (genuinely alternating, since
//! `B(x,x) = Tr(2·x σ^k x) = 0` in char 2). This is precisely the family the Gold
//! research thread builds by hand (`experiments/trace_form_arf.py`): with
//! `E = F_{2^m}`, `σ = ` Frobenius, `k = a`, `Q_a(x) = Tr(x^{1+2^a})` is the Gold
//! form. [`gold_form`] lands it in the typed core over the nim subfields, where
//! `.classify()` reads off the Arf invariant (rank, radical, win-bias zero-count).
//!
//! The same construction over `Surcomplex` (`σ = ` conjugation, `k = 1`) gives the
//! **norm form** `Tr(x·x̄) = 2(a²+b²)` — the binary Pfister/norm form; over `Qq` it
//! gives an unramified local trace form using a Teichmuller-lifted residue basis; and
//! over an odd-characteristic `Fpn` it gives an ordinary diagonalizable trace form.
//!
//! Boundary: the form has dimension `[E:F]`, so as a [`Metric`] it is capped at
//! `MAX_BASIS_DIM = 128` — exactly the degree of the full nim-field `F_{2^128}`.

use crate::clifford::{Metric, MAX_BASIS_DIM};
use crate::forms::ArfResult;
use crate::scalar::{
    nim_square, nim_trace, CyclicGaloisExtension, FieldExtension, Fp, Nimber, Scalar,
};
use std::collections::BTreeMap;

/// Assemble the twisted form `Q(x) = trace(x · twist(x))` over a basis: the shared
/// core behind [`trace_twisted_form`] (trait-driven) and [`gold_form`] (nim-native).
/// `twist` is `σ^k` and `trace` is the relative trace `E → F` (= `T`).
fn assemble_twisted_form<E: Scalar, T: Scalar>(
    basis: &[E],
    twist: impl Fn(&E) -> E,
    trace: impl Fn(&E) -> T,
) -> Metric<T> {
    let n = basis.len();
    let tw: Vec<E> = basis.iter().map(&twist).collect();

    // diagonal: q_i = Tr(e_i · σ^k(e_i))
    let q: Vec<T> = basis
        .iter()
        .zip(&tw)
        .map(|(e, te)| trace(&e.mul(te)))
        .collect();

    // polar: b_{ij} = Tr(e_i σ^k(e_j) + e_j σ^k(e_i))  (i < j), sparse
    let mut b = BTreeMap::new();
    for i in 0..n {
        for j in (i + 1)..n {
            let t = trace(&basis[i].mul(&tw[j]).add(&basis[j].mul(&tw[i])));
            if !t.is_zero() {
                b.insert((i, j), t);
            }
        }
    }

    Metric::general(q, b, BTreeMap::new())
}

fn insert_metric_block<S: Scalar>(
    q: &mut [S],
    b: &mut BTreeMap<(usize, usize), S>,
    offset: usize,
    block: Metric<S>,
) {
    let (bq, bb, ba) = block.into_parts();
    debug_assert!(ba.is_empty());
    for (i, qi) in bq.into_iter().enumerate() {
        q[offset + i] = qi;
    }
    for ((i, j), v) in bb {
        b.insert((offset + i, offset + j), v);
    }
}

/// The Frobenius-twisted trace form `Q_k(x) = Tr_{E/F}(x · σ^k(x))` of a cyclic
/// Galois extension `E/F`, as a [`Metric`] over the base `F` in the distinguished
/// [`basis`](CyclicGaloisExtension::basis) `(e_0,…,e_{d-1})`:
///
/// ```text
/// q_i    = Tr(e_i · σ^k(e_i))
/// b_{ij} = Tr(e_i σ^k(e_j) + e_j σ^k(e_i))      (i < j)
/// ```
///
/// `k = 1` is the standard choice (`σ` itself); larger `k` gives the higher Gold
/// exponents `1 + 2^k`.
///
/// **The transfer reading (`k = 0`).** With `σ^0 = id` the twist vanishes and the
/// form is `Tr_{E/F}(x·y)` — the **Scharlau transfer** `s_*(⟨1⟩)` of the unit form
/// `⟨1⟩ ∈ W(E)` along `s = Tr_{E/F}` (Lam, *Introduction to Quadratic Forms over
/// Fields*, GSM 67, Ch. VII; Scharlau, *Quadratic and Hermitian Forms*, Ch. 2).
/// Equivalently `trace_twisted_form::<E>(0) == transfer_diagonal(&[E::one()])`. The
/// general transfer of a diagonal form is [`transfer_diagonal`].
pub fn trace_twisted_form<E>(k: usize) -> Metric<E::Base>
where
    E: CyclicGaloisExtension,
{
    assemble_twisted_form(&E::basis(), |e| e.sigma_power(k), |z| z.trace())
}

/// The cyclic-algebra trace form `T_A(z) = Trd_A(z²)` for the crossed product
/// algebra `A = (E/F, σ, a) = ⊕ᵢ E·uⁱ`, with `uⁿ = a` and `u·x = σ(x)·u`.
/// The basis is ordered by `u`-line: `(e_0, …, e_{n-1})`,
/// `(e_0u, …, e_{n-1}u)`, …, where `(e_i)` is
/// [`CyclicGaloisExtension::basis`].
///
/// Reduced trace sees only the `u⁰` coefficient, so the form splits into the
/// self-line `E`, the middle self-line `E·u^{n/2}` when `n` is even, and pure
/// polar pairings between `E·uⁱ` and `E·u^{n-i}`:
///
/// ```text
/// T_A(Σ x_i u^i) =
///   Tr(x_0²)
///   + [n even] Tr(a · x_{n/2} σ^{n/2}(x_{n/2}))
///   + Σ_{0<i<n-i} Tr(a · (x_i σ^i(x_{n-i}) + x_{n-i} σ^{n-i}(x_i))).
/// ```
///
/// This is the literal quadratic trace form, not the degree-2 reduced norm. For
/// a quaternion algebra the relation is the Cayley-Hamilton identity
/// `Trd(z²) = Trd(z)² - 2·Nrd(z)`.
pub fn cyclic_algebra_trace_form<E>(a: &E::Base) -> Metric<E::Base>
where
    E: CyclicGaloisExtension,
{
    let basis = E::basis();
    let n = basis.len();
    let dim = n
        .checked_mul(n)
        .expect("cyclic algebra trace-form dimension overflowed");
    assert!(
        dim <= MAX_BASIS_DIM,
        "cyclic_algebra_trace_form has dimension [E:F]^2={dim}, exceeding {MAX_BASIS_DIM}"
    );

    let mut q = vec![E::Base::zero(); dim];
    let mut b = BTreeMap::new();

    let line0 = assemble_twisted_form(&basis, |x| x.clone(), |z| z.trace());
    insert_metric_block(&mut q, &mut b, 0, line0);

    if n % 2 == 0 {
        let mid = n / 2;
        let middle = assemble_twisted_form(&basis, |x| x.sigma_power(mid), |z| a.mul(&z.trace()));
        insert_metric_block(&mut q, &mut b, mid * n, middle);
    }

    for i in 1..n {
        let j = n - i;
        if i >= j {
            continue;
        }
        for r in 0..n {
            for s in 0..n {
                let term = basis[r]
                    .mul(&basis[s].sigma_power(i))
                    .add(&basis[s].mul(&basis[r].sigma_power(j)));
                let value = a.mul(&term.trace());
                if !value.is_zero() {
                    b.insert((i * n + r, j * n + s), value);
                }
            }
        }
    }

    Metric::general(q, b, BTreeMap::new())
}

/// The **Scharlau transfer** `s_*(⟨λ_1,…,λ_r⟩)` of a diagonal form over `E`, pushed
/// to `W(F)` along the field trace `s = Tr_{E/F}` (Lam, GSM 67, Ch. VII; Scharlau,
/// *Quadratic and Hermitian Forms*, Ch. 2). Each diagonal entry `λ ∈ E` contributes
/// the `[E:F]`-dimensional block `(x, y) ↦ Tr_{E/F}(λ·x·y)` over the distinguished
/// [`basis`](CyclicGaloisExtension::basis) `(e_0,…,e_{d-1})`:
///
/// ```text
/// q_a      = Tr(λ · e_a · e_a)
/// b_{ab}   = Tr(λ · (e_a e_b + e_b e_a))            (a < b)
/// ```
///
/// and the blocks are orthogonally summed — `s_*` is additive on `⟂`. The unit case
/// `transfer_diagonal(&[E::one()])` is `Tr(x·y)`, i.e. `s_*(⟨1⟩)`, the `k = 0` member
/// of [`trace_twisted_form`].
///
/// `s_*` is a group homomorphism `W(E) → W(F)` (the transfer of a hyperbolic form is
/// hyperbolic), and satisfies **Frobenius reciprocity** `s_*(r*(x)·y) = x·s_*(y)`,
/// where `r* : W(F) → W(E)` is restriction. Restriction itself is **injective** for
/// odd `[E:F]` (Springer's odd-degree theorem) — the companion to the Springer
/// *residue* theorem that drives the local layer.
///
/// **Boundary:** char `≠ 2`. In characteristic 2 the trace form `Tr(x·σ^k(x))`
/// degenerates (the `Tr(x²)` trap this module documents); the char-2 transfer story
/// is the Artin–Schreier route in `function_field_char2.rs`. The total dimension is
/// `r·[E:F]`, so as a [`Metric`] it is capped at `MAX_BASIS_DIM = 128`.
pub fn transfer_diagonal<E>(entries: &[E]) -> Metric<E::Base>
where
    E: CyclicGaloisExtension,
{
    let basis = E::basis();
    let mut result = Metric::diagonal(Vec::new());
    for lambda in entries {
        let block = assemble_twisted_form(&basis, |x| lambda.mul(x), |z| z.trace());
        result = result.direct_sum(&block);
    }
    result
}

/// The Arf invariant of the **char-2** twisted trace form of `E/F_2` — the typed
/// bridge for the finite-field tower. Builds `Q_k` over `F_2`, lifts the
/// coefficients `F_2 ↪ Nimber` (so the char-2 [`ArfResult`] classifier can read the
/// form), and returns its Arf data. For `E = Fpn<2,m>` with `k = a` this is the Gold
/// form `Tr(x^{1+2^a})`; see [`gold_form`] for the nim-native construction that
/// reaches the larger power-of-two fields.
pub fn trace_form_arf<E>(k: usize) -> Option<ArfResult>
where
    E: CyclicGaloisExtension + FieldExtension<Base = Fp<2>>,
{
    trace_twisted_form::<E>(k)
        .map(|x| Nimber(x.value()))
        .classify()
}

/// The **Gold form** `Q_a(x) = Tr_{F_{2^m}/F_2}(x^{1+2^a})` over the nim subfield
/// `F_{2^m} ⊂ Nimber`, as a [`Metric`]`<Nimber>` (already `F_2`-valued, ready for
/// `.classify()` → [`ArfResult`]). This is the central object of the game-built
/// quadratic-form thread (mirrors `experiments/gold_form_from_games.py`): the bit
/// basis `{1, 2, …, 2^{m-1}}` is an `F_2`-basis of `F_{2^m}`, the twist `σ^a` is the
/// `a`-fold nim-Frobenius `x ↦ x^{2^a}`, and the trace is `nim_trace(·, m)`.
///
/// `m` must be a **power of two** `≤ 128`: only then do the nimbers `< 2^m` form a
/// subfield (`F_{2^{2^k}}`) closed under nim-multiplication. The Gold-rank theorem
/// gives `rank = m − gcd(2a, m)`.
pub fn gold_form(m: usize, a: usize) -> Metric<Nimber> {
    assert!(
        m.is_power_of_two() && m <= 128,
        "the nimbers < 2^m form a subfield only for m a power of two ≤ 128"
    );
    let basis: Vec<Nimber> = (0..m).map(|i| Nimber(1u128 << i)).collect();
    assemble_twisted_form(
        &basis,
        |x| {
            // σ^a = the a-fold nim-Frobenius x ↦ x^{2^a}
            let mut t = x.0;
            for _ in 0..a {
                t = nim_square(t);
            }
            Nimber(t)
        },
        |x| Nimber(nim_trace(x.0, m as u128)),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scalar::{Fp, Fpn, Qq, Rational, Surcomplex};

    fn r(n: i128) -> Rational {
        Rational::int(n)
    }

    fn gcd(a: usize, b: usize) -> usize {
        if b == 0 {
            a
        } else {
            gcd(b, a % b)
        }
    }

    fn eval_rational_metric(m: &Metric<Rational>, coords: &[i128]) -> Rational {
        assert_eq!(m.dim(), coords.len());
        let mut total = Rational::zero();
        for (i, &ci) in coords.iter().enumerate() {
            let x = r(ci);
            total = total.add(&m.q[i].mul(&x).mul(&x));
        }
        for (&(i, j), bij) in &m.b {
            total = total.add(&bij.mul(&r(coords[i])).mul(&r(coords[j])));
        }
        total
    }

    #[test]
    fn surcomplex_twist_is_the_norm_form() {
        // E = ℚ(i)/ℚ, σ = conjugation, k = 1: Q(x) = Tr(x·x̄) = 2(a²+b²), the binary
        // norm form ⟨2, 2⟩ (diagonal, no polar term).
        let m = trace_twisted_form::<Surcomplex<Rational>>(1);
        assert_eq!(m.q, vec![Rational::int(2), Rational::int(2)]);
        assert!(m.b.is_empty());
    }

    #[test]
    fn cyclic_trace_form_degree_two_is_literal_trd_square() {
        // For A = (Q(i)/Q, conjugation, a), ordered as (1, i, u, iu),
        // Trd(z^2) is <2, -2, 2a, 2a>. This is adjacent to, but not equal to,
        // the reduced norm <1, 1, -a, -a>.
        for a in [-3i128, -1, 2, 5] {
            let m = cyclic_algebra_trace_form::<Surcomplex<Rational>>(&r(a));
            assert_eq!(m.q, vec![r(2), r(-2), r(2 * a), r(2 * a)]);
            assert!(m.b.is_empty());
        }
    }

    #[test]
    fn cyclic_trace_form_degree_two_satisfies_cayley_hamilton_relation() {
        // The honest degree-2 tie to the shipped norm-form oracle is
        // Trd(z^2) = Trd(z)^2 - 2*Nrd(z), not equality with Nrd.
        for a in [-3i128, 2, 5] {
            let m = cyclic_algebra_trace_form::<Surcomplex<Rational>>(&r(a));
            for p in -1..=1 {
                for q in -1..=1 {
                    for u in -1..=1 {
                        for v in -1..=1 {
                            let lhs = eval_rational_metric(&m, &[p, q, u, v]);
                            let trd = 2 * p;
                            let nrd = p * p + q * q - a * u * u - a * v * v;
                            let rhs = r(trd * trd - 2 * nrd);
                            assert_eq!(lhs, rhs, "a={a}, coords={:?}", [p, q, u, v]);
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn qq_twist_uses_the_unramified_galois_basis() {
        // E = Q_9/Q_3: the same trace-form bridge now reaches the unramified local
        // leg via the Teichmüller-lifted residue basis and the Witt-Frobenius.
        type Q9 = Qq<3, 3, 2>;
        let m = trace_twisted_form::<Q9>(1);
        assert_eq!(m.q.len(), 2);
        assert!(m.q.iter().all(|x| !x.is_zero()));
        assert!(m.q.iter().all(|x| x.valuation().is_some()));
    }

    #[test]
    fn cyclic_trace_form_degree_three_has_hyperbolic_cross_pair() {
        // For n = 3 there is one self-line (u^0) and one pure polar pair
        // (u^1, u^2). The pair is hyperbolic, so the Witt anisotropic kernel
        // matches the u^0 trace-square block.
        let t = cyclic_algebra_trace_form::<Fpn<3, 3>>(&Fp::<3>::one());
        let line0 = trace_twisted_form::<Fpn<3, 3>>(0);
        assert_eq!(t.dim(), 9);
        assert_eq!(&t.q[..3], line0.q());
        assert!(t.q[3..].iter().all(|x| x.is_zero()));
        for (&(i, j), v) in line0.b() {
            assert_eq!(t.b.get(&(i, j)), Some(v));
        }

        let t_dec = t.witt_decompose().expect("F_3 trace form decomposition");
        let line_dec = line0
            .witt_decompose()
            .expect("F_3 line trace form decomposition");
        assert_eq!(t_dec.anisotropic_dim, line_dec.anisotropic_dim);
    }

    #[test]
    fn gold_form_over_small_fpn_matches_rank_formula() {
        // The typed finite-field path: Gold Q_a over Fpn<2,m>, m = 2, 3.
        // F_4 (m=2), a=1: gcd(2,2)=2 ⇒ Q ≡ 0, all radical.
        let f4 = trace_form_arf::<Fpn<2, 2>>(1).unwrap();
        assert_eq!((f4.rank, f4.radical_dim), (0, 2));
        // F_8 (m=3), a=1: gcd(2,3)=1 ⇒ rank 2, radical 1.
        let f8 = trace_form_arf::<Fpn<2, 3>>(1).unwrap();
        assert_eq!((f8.rank, f8.radical_dim), (2, 1));
    }

    #[test]
    fn gold_form_over_nim_subfields_matches_rank_formula() {
        // The nim-native path reaches the power-of-two fields the Gold survey uses
        // (F_16, F_256). arf_nimber computes rank by independent symplectic reduction
        // of the polar form — agreement with m − gcd(2a, m) is a real cross-check.
        for m in [2usize, 4, 8] {
            let a = 1usize;
            let arf = gold_form(m, a).classify().unwrap();
            let g = gcd(2 * a, m);
            assert_eq!(
                (arf.rank, arf.radical_dim),
                (m - g, g),
                "Gold form over F_2^{m} (a={a})"
            );
        }
        // a higher Gold exponent: m=8, a=3 ⇒ gcd(6,8)=2 ⇒ rank 6.
        let arf = gold_form(8, 3).classify().unwrap();
        assert_eq!((arf.rank, arf.radical_dim), (6, 2));
    }

    #[test]
    fn transfer_of_unit_form_is_the_k0_twisted_form() {
        // s_*(⟨1⟩) = Tr(x·y) is the k = 0 member of the twisted-form family.
        let s = transfer_diagonal::<Fpn<3, 2>>(&[Fpn::<3, 2>::one()]);
        let t0 = trace_twisted_form::<Fpn<3, 2>>(0);
        assert_eq!(s.q, t0.q);
        assert_eq!(s.b, t0.b);
    }

    #[test]
    fn transfer_of_a_hyperbolic_form_is_split() {
        // s_* : W(E) → W(F) is a group homomorphism, so the transfer of the
        // hyperbolic form ⟨1, −1⟩ over E is Witt-trivial over F.
        let one = Fpn::<3, 2>::one();
        let hyp = transfer_diagonal::<Fpn<3, 2>>(&[one, one.neg()]);
        let dec = hyp.witt_decompose().expect("Fp<3> Witt decomposition");
        assert_eq!(
            dec.anisotropic_dim, 0,
            "transfer of a hyperbolic form splits"
        );
    }

    #[test]
    fn frobenius_reciprocity_projection_formula() {
        // s_*(r*(⟨c⟩) · ⟨λ⟩) = ⟨c⟩ · s_*(⟨λ⟩):  c ∈ F factors out of the F-linear
        // trace, so the transfer of (c·λ) equals the c-scaling of the transfer of λ.
        let c = Fp::<3>::new(2); // a unit of F_3
        let lam = Fpn::<3, 2>::from_coeffs(&[1, 1]); // 1 + x ∈ F_9
        let lhs = transfer_diagonal::<Fpn<3, 2>>(&[Fpn::<3, 2>::embed(&c).mul(&lam)]);
        let base = transfer_diagonal::<Fpn<3, 2>>(&[lam]);
        let scaled_q: Vec<Fp<3>> = base.q.iter().map(|x| c.mul(x)).collect();
        let scaled_b: BTreeMap<(usize, usize), Fp<3>> =
            base.b.iter().map(|(k, v)| (*k, c.mul(v))).collect();
        assert_eq!(lhs.q, scaled_q);
        assert_eq!(lhs.b, scaled_b);
    }

    #[test]
    fn springer_odd_degree_restriction_is_injective() {
        // r* : W(K) → W(E) is injective for odd [E:K] (Springer's odd-degree
        // theorem). Witness: the anisotropic binary form ⟨1,1⟩/F_3 stays anisotropic
        // over F_27 (degree 3, odd) — its nonzero Witt class does not die.
        let aniso = Metric::<Fp<3>>::diagonal(vec![Fp::<3>::one(), Fp::<3>::one()]);
        let base_dec = aniso.witt_decompose().expect("Fp<3> Witt decomposition");
        assert_eq!(base_dec.anisotropic_dim, 2, "⟨1,1⟩ anisotropic over F_3");

        let restricted =
            Metric::<Fpn<3, 3>>::diagonal(vec![Fpn::<3, 3>::one(), Fpn::<3, 3>::one()]);
        match restricted
            .witt_decompose()
            .expect("F_27 Witt decomposition")
        {
            crate::forms::FiniteFieldWittDecomp::Odd(d) => {
                assert_eq!(
                    d.anisotropic_dim, 2,
                    "still anisotropic over F_27 ⇒ injective"
                );
            }
            other => panic!("expected odd-characteristic decomposition, got {other:?}"),
        }
    }

    #[test]
    fn metric_map_lifts_fp2_to_nimber() {
        // base-change F_2 ↪ Nimber preserves the form's structure.
        let over_f2 = trace_twisted_form::<Fpn<2, 3>>(1);
        let lifted = over_f2.map(|x| Nimber(x.value()));
        assert_eq!(lifted.q.len(), over_f2.q.len());
        for (i, qi) in over_f2.q.iter().enumerate() {
            assert_eq!(lifted.q[i].0, qi.value());
        }
        assert_eq!(lifted.b.len(), over_f2.b.len());
    }
}
