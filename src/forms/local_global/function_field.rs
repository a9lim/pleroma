//! Places, the Hilbert symbol, and Hilbert reciprocity over the **global function
//! field** `F_q(t)` — the equal-characteristic (char `p`) mirror of
//! [`forms::padic`](crate::forms::padic) over `ℚ`.
//!
//! `F_q(t)` is a global field exactly like `ℚ`, with one structural simplification
//! and one structural difference:
//!
//!   * **simpler:** `q` is odd, so *every* residue field `κ(π) = F_q[t]/(π) =
//!     F_{q^{deg π}}` has odd characteristic. The Hilbert symbol is therefore the
//!     **tame symbol** at every place — there is **no `p = 2` branch** (the messy
//!     mod-8 case that `local_global/padic.rs` carries). The residue-characteristic-2 boundary
//!     is the `springer_laurent` boundary; char-2 function fields are out of scope.
//!   * **different:** there is **no archimedean place**. The role of `ℝ` is played
//!     by the **degree place `∞`** (uniformizer `1/t`, residue field `F_q`), which
//!     is just another tame place — so Hasse–Minkowski over `F_q(t)` has no
//!     definiteness condition (see
//!     [`try_is_isotropic_ff`](crate::forms::try_is_isotropic_ff)).
//!
//! The places of `F_q(t)`: the **finite** places are the monic irreducible
//! polynomials `π(t) ∈ F_q[t]` (residue field `F_{q^{deg π}}`), and the one
//! **infinite** place `∞` is the degree valuation `v_∞(f) = deg(den) − deg(num)`.
//! Reciprocity `∏_v (a,b)_v = +1` (Weil) is the gold oracle, exact here.
//!
//! Entries are elements of [`RationalFunction`] `= F_q(t)`; everything reduces to
//! [`Poly`] arithmetic over `F_q`, with the residue quadratic character computed by
//! Euler's criterion `u^{(|κ|−1)/2}` in `F_q[t]/(π)`.

use crate::forms::{is_square_finite, FiniteOddField};
use crate::scalar::{Poly, RationalFunction, Scalar};

/// A place of `F_q(t)`: the degree place `∞`, or a finite place given by a monic
/// irreducible `π(t)`. The mirror of [`Place`](crate::forms::Place)`{Real,Prime}`.
#[derive(Debug, Clone, PartialEq)]
pub enum FFPlace<S: FiniteOddField> {
    /// The degree place `∞` (uniformizer `1/t`, residue field `F_q`).
    Infinite,
    /// A finite place: a monic irreducible `π(t)` (residue field `F_q[t]/(π)`).
    Finite(Poly<S>),
}

// ───────────────────────── factorization over F_q ─────────────────────────

/// The distinct monic irreducible factors of `f` over `F_q` (the square-free
/// support — multiplicities dropped, which is all the place layer needs).
pub fn monic_irreducible_factors<S: FiniteOddField>(f: &Poly<S>) -> Vec<Poly<S>> {
    crate::forms::poly_factor::monic_irreducible_factor_support(
        f,
        S::characteristic_prime(),
        S::field_order(),
        S::from_index,
    )
}

/// The multiplicity of `pi` in `p` (and `p` with all those factors stripped).
fn strip_factor<S: Scalar>(mut p: Poly<S>, pi: &Poly<S>) -> (i128, Poly<S>) {
    let mut mult = 0i128;
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

// ───────────────────────── per-place local data ─────────────────────────

/// The residue field order `|κ| = q^{deg π}` (or `q` at the degree place).
fn try_kappa_order<S: FiniteOddField>(place: &FFPlace<S>) -> Option<u128> {
    let q = S::field_order();
    match place {
        FFPlace::Finite(pi) => {
            let deg = pi
                .degree()
                .expect("an irreducible has degree ≥ 1")
                .try_into()
                .ok()?;
            q.checked_pow(deg)
        }
        FFPlace::Infinite => Some(q),
    }
}

/// The valuation `v_place(a)` of a **nonzero** `a ∈ F_q(t)`.
pub fn try_valuation_at_ff<S: FiniteOddField>(
    a: &RationalFunction<S>,
    place: &FFPlace<S>,
) -> Option<i128> {
    if a.is_zero() {
        return None;
    }
    Some(match place {
        FFPlace::Finite(pi) => {
            let (mn, _) = strip_factor(a.num().clone(), pi);
            let (md, _) = strip_factor(a.den().clone(), pi);
            mn - md
        }
        FFPlace::Infinite => {
            let dn = a.num().degree().expect("nonzero numerator") as i128;
            let dd = a.den().degree().expect("monic nonzero denominator") as i128;
            dd - dn // v_∞ = deg(den) − deg(num)
        }
    })
}

/// The residue unit `(a / ϖ^{v(a)}) mod ϖ ∈ κ*` of a **nonzero** `a`, as an element
/// of the residue field: a [`Poly`] of degree `< deg π` at a finite place, or a
/// constant (an `F_q` element) at the degree place.
fn try_residue_unit_at<S: FiniteOddField>(
    a: &RationalFunction<S>,
    place: &FFPlace<S>,
) -> Option<Poly<S>> {
    if a.is_zero() {
        return None;
    }
    match place {
        FFPlace::Finite(pi) => {
            let (_, num_s) = strip_factor(a.num().clone(), pi);
            let (_, den_s) = strip_factor(a.den().clone(), pi);
            let num_mod = num_s.rem(pi);
            let den_mod = den_s.rem(pi);
            // den_mod⁻¹ in κ* by Fermat: x^{|κ|−2} (κ* is cyclic of order |κ|−1).
            let den_inv = den_mod.pow_mod(try_kappa_order(place)?.checked_sub(2)?, pi);
            Some(num_mod.mul_mod(&den_inv, pi))
        }
        FFPlace::Infinite => {
            // a·t^{v_∞} → (lead num)/(lead den) as t → ∞.
            let ln = *a.num().leading().expect("nonzero numerator");
            let ld = *a.den().leading().expect("monic nonzero denominator");
            Some(Poly::constant(ln.mul(&ld.inv()?)))
        }
    }
}

/// The residue quadratic character `χ_κ(u) ∈ {+1, −1}` of a **nonzero** residue
/// unit `u ∈ κ*` — Euler's criterion `u^{(|κ|−1)/2}` in `F_q[t]/(π)` (or in `F_q`
/// at the degree place).
fn try_chi_kappa<S: FiniteOddField>(unit: &Poly<S>, place: &FFPlace<S>) -> Option<i128> {
    match place {
        FFPlace::Finite(pi) => {
            let e = (try_kappa_order(place)?.checked_sub(1)?) / 2;
            Some(if unit.pow_mod(e, pi) == Poly::one() {
                1
            } else {
                -1 // the unique order-2 element of κ* is −1
            })
        }
        FFPlace::Infinite => Some(if is_square_finite::<S>(unit.coeff(0)) {
            1
        } else {
            -1
        }),
    }
}

/// Whether a **nonzero** `a` is a square in the local field at `place`: the
/// valuation is even **and** the residue unit is a square in `κ`. The mirror of
/// [`try_is_square_qp`](crate::forms::try_is_square_qp).
pub fn try_is_local_square_ff<S: FiniteOddField>(
    a: &RationalFunction<S>,
    place: &FFPlace<S>,
) -> Option<bool> {
    if a.is_zero() {
        return Some(false);
    }
    Some(
        try_valuation_at_ff(a, place)?.rem_euclid(2) == 0
            && try_chi_kappa(&try_residue_unit_at(a, place)?, place)? == 1,
    )
}

// ───────────────────────── the Hilbert symbol ─────────────────────────

/// The Hilbert symbol `(a, b)_v` over the completion of `F_q(t)` at `place`, for
/// **nonzero** `a, b` — the **tame symbol**. With `α = v(a)`, `β = v(b)` and
/// residue units `ā, b̄`,
/// `(a,b)_v = χ_κ((−1)^{αβ}) · χ_κ(ā)^β · χ_κ(b̄)^α`,
/// exactly the odd-`p` branch of [`try_hilbert_symbol_qp`](crate::forms::try_hilbert_symbol_qp)
/// with the residue Legendre symbol replaced by the residue character `χ_κ`. No
/// `p = 2` branch exists because every residue field has odd characteristic.
pub fn try_hilbert_symbol_ff<S: FiniteOddField>(
    a: &RationalFunction<S>,
    b: &RationalFunction<S>,
    place: &FFPlace<S>,
) -> Option<i128> {
    if a.is_zero() || b.is_zero() {
        return None;
    }
    let al = try_valuation_at_ff(a, place)?;
    let be = try_valuation_at_ff(b, place)?;
    let ca = try_chi_kappa(&try_residue_unit_at(a, place)?, place)?;
    let cb = try_chi_kappa(&try_residue_unit_at(b, place)?, place)?;
    // χ_κ(−1): −1 is a square in κ iff |κ| ≡ 1 (mod 4).
    let chi_neg1 = if try_kappa_order(place)? % 4 == 1 {
        1
    } else {
        -1
    };
    // Exactly the shared tame symbol — the same machine as the odd-`p` Q_p branch,
    // with the residue character `χ_κ` in place of the Legendre symbol.
    Some(crate::forms::padic::tame_hilbert_symbol(
        al, be, ca, cb, chi_neg1,
    ))
}

// ───────────────────────── Hasse invariant + reciprocity ─────────────────────────

/// The relevant places of a list of **nonzero** entries: every finite place where
/// some entry has nonzero valuation (the monic irreducible factors of all
/// numerators and denominators), plus the degree place `∞`. At every other place
/// all entries are units, so every symbol is `+1`. Mirror of
/// [`relevant_primes`](crate::forms::padic).
pub fn try_relevant_places_ff<S: FiniteOddField>(
    entries: &[RationalFunction<S>],
) -> Option<Vec<FFPlace<S>>> {
    if entries.iter().any(|a| a.is_zero()) {
        return None;
    }
    let mut polys: Vec<Poly<S>> = Vec::new();
    for a in entries {
        let factors = monic_irreducible_factors(a.num())
            .into_iter()
            .chain(monic_irreducible_factors(a.den()));
        for pi in factors {
            if !polys.contains(&pi) {
                polys.push(pi);
            }
        }
    }
    let mut places: Vec<FFPlace<S>> = polys.into_iter().map(FFPlace::Finite).collect();
    places.push(FFPlace::Infinite);
    Some(places)
}

/// The Hasse invariant `ε_v(⟨a_1,…,a_n⟩) = ∏_{i<j} (a_i, a_j)_v` at `place`. The
/// `_ff` suffix distinguishes it from the `ℚ` [`hasse_at_place`](crate::forms::hasse_at_place).
pub fn try_hasse_at_place_ff<S: FiniteOddField>(
    entries: &[RationalFunction<S>],
    place: &FFPlace<S>,
) -> Option<i128> {
    let mut h = 1i128;
    for i in 0..entries.len() {
        for j in (i + 1)..entries.len() {
            h *= try_hilbert_symbol_ff(&entries[i], &entries[j], place)?;
        }
    }
    Some(h)
}

/// The **Hilbert reciprocity product** `∏_v (a,b)_v` over all places of `F_q(t)`,
/// `+1` for every nonzero `a, b` (Weil reciprocity / the product formula). Exact —
/// the symbols are `+1` at all but the finitely many relevant places.
pub fn try_hilbert_reciprocity_product_ff<S: FiniteOddField>(
    a: &RationalFunction<S>,
    b: &RationalFunction<S>,
) -> Option<i128> {
    <RationalFunction<S> as crate::forms::GlobalField>::try_reciprocity_product(a, b)
}

// ───────────────────────── Hasse–Minkowski over F_q(t) ─────────────────────────

/// `−1 ∈ F_q(t)`.
fn neg_one<S: FiniteOddField>() -> RationalFunction<S> {
    RationalFunction::from_base(S::one().neg())
}

/// The discriminant `∏ aᵢ ∈ F_q(t)` (its square class is the form's discriminant).
fn disc<S: FiniteOddField>(entries: &[RationalFunction<S>]) -> RationalFunction<S> {
    let mut d = RationalFunction::one();
    for e in entries {
        d = d.mul(e);
    }
    d
}

/// Whether a **nonzero** `x ∈ F_q(t)` is a global square. The square class of `x`
/// is represented by `num(x)·den(x) ∈ F_q[t]`, which is a square iff its `F_q`
/// leading coefficient is a square **and** every irreducible factor has even
/// multiplicity. The char-`p` analogue of `is_perfect_square` over `ℤ`.
pub(crate) fn is_global_square_ff<S: FiniteOddField>(x: &RationalFunction<S>) -> bool {
    if x.is_zero() {
        return false;
    }
    let f = x.num().mul(x.den());
    let lead = *f.leading().expect("nonzero square-class representative");
    if !is_square_finite::<S>(lead) {
        return false;
    }
    let mut g = f.make_monic();
    for pi in monic_irreducible_factors(&f) {
        let (mult, rest) = strip_factor(g, &pi);
        if mult.rem_euclid(2) != 0 {
            return false;
        }
        g = rest;
    }
    true
}

/// Local isotropy of a nondegenerate diagonal form `⟨a_1,…,a_n⟩` over the
/// completion of `F_q(t)` at `place`, by rank — the exact mirror of
/// [`try_is_isotropic_at_p`](crate::forms::padic) (`F_q(t)` and `Q_p` share the
/// u-invariant `4`, so the thresholds match): `n≤1` never, `n=2` iff `−a_1a_2` is a
/// local square, `n=3`/`4` the Hilbert conditions, `n≥5` always. Entries nonzero.
pub fn try_is_isotropic_at_place_ff<S: FiniteOddField>(
    entries: &[RationalFunction<S>],
    place: &FFPlace<S>,
) -> Option<bool> {
    if entries.iter().any(|a| a.is_zero()) {
        return Some(true);
    }
    Some(match entries.len() {
        0 | 1 => false,
        2 => try_is_local_square_ff(&entries[0].mul(&entries[1]).neg(), place)?,
        3 => {
            let d = disc(entries);
            try_hilbert_symbol_ff(&neg_one(), &d.neg(), place)?
                == try_hasse_at_place_ff(entries, place)?
        }
        4 => {
            let d = disc(entries);
            !try_is_local_square_ff(&d, place)?
                || try_hasse_at_place_ff(entries, place)?
                    == try_hilbert_symbol_ff(&neg_one::<S>(), &neg_one::<S>(), place)?
        }
        _ => true,
    })
}

/// Whether a diagonal form over `F_q(t)` is **isotropic**, by **Hasse–Minkowski**:
/// isotropic over the global field iff isotropic at every place. Unlike `ℚ` there
/// is **no archimedean place**, so rank ≥ 3 needs only the local conditions at the
/// relevant places (all other places are unit forms, automatically isotropic for
/// rank ≥ 3 over a finite residue field). Rank 2 reduces to `−a_1a_2` being a
/// global square; a zero entry is an isotropic direction.
pub fn try_is_isotropic_ff<S: FiniteOddField>(entries: &[RationalFunction<S>]) -> Option<bool> {
    <RationalFunction<S> as crate::forms::GlobalField>::try_is_isotropic_global(entries)
}

/// The per-place isotropy breakdown of a rank-`≥3` form — the function-field
/// analogue of [`AdelicIsotropy`](crate::forms::AdelicIsotropy) (no real place).
#[derive(Debug, Clone)]
pub struct FFAdelicIsotropy<S: FiniteOddField> {
    /// `(place, is_isotropic_there)` at each relevant place.
    pub local: Vec<(FFPlace<S>, bool)>,
}

impl<S: FiniteOddField> FFAdelicIsotropy<S> {
    /// Globally isotropic ⟺ isotropic at every place (Hasse–Minkowski).
    pub fn is_global(&self) -> bool {
        self.local.iter().all(|(_, iso)| *iso)
    }
}

/// The adelic Hasse–Minkowski breakdown of a rank-`≥3` form over `F_q(t)`.
pub fn try_isotropy_over_ff_adeles<S: FiniteOddField>(
    entries: &[RationalFunction<S>],
) -> Option<FFAdelicIsotropy<S>> {
    if entries.iter().any(|a| a.is_zero()) {
        return Some(FFAdelicIsotropy { local: Vec::new() });
    }
    let mut local = Vec::new();
    for pl in try_relevant_places_ff(entries)? {
        let iso = try_is_isotropic_at_place_ff(entries, &pl)?;
        local.push((pl, iso));
    }
    Some(FFAdelicIsotropy { local })
}

/// The places where the quaternion algebra `(a, b)` over `F_q(t)` **ramifies** —
/// where `(a,b)_v = −1`. The count is always **even** (the additive form of
/// reciprocity / the even-ramification theorem), the function-field mirror of
/// [`brauer_local_invariants`](crate::forms::brauer_local_invariants).
pub fn try_ramified_places_ff<S: FiniteOddField>(
    a: &RationalFunction<S>,
    b: &RationalFunction<S>,
) -> Option<Vec<FFPlace<S>>> {
    <RationalFunction<S> as crate::forms::GlobalField>::try_ramified_places(a, b)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scalar::Fp;

    type F = RationalFunction<Fp<5>>;
    type PolyF = Poly<Fp<5>>;

    fn rf(num: &[i128], den: &[i128]) -> F {
        RationalFunction::new(
            num.iter().map(|&n| Fp::<5>::new(n)).collect(),
            den.iter().map(|&n| Fp::<5>::new(n)).collect(),
        )
    }
    fn poly(c: &[i128]) -> PolyF {
        Poly::new(c.iter().map(|&n| Fp::<5>::new(n)).collect())
    }

    #[test]
    fn factors_into_monic_irreducibles() {
        // x² − 1 = (x − 1)(x + 1) over F_5  (constant −1 ≡ 4)
        let fs = monic_irreducible_factors(&poly(&[4, 0, 1]));
        assert_eq!(fs.len(), 2);
        assert!(fs.contains(&poly(&[4, 1]))); // x − 1
        assert!(fs.contains(&poly(&[1, 1]))); // x + 1
                                              // x² + 2 is irreducible over F_5 (−2 ≡ 3 is a nonsquare): one factor, itself.
        let irr = monic_irreducible_factors(&poly(&[2, 0, 1]));
        assert_eq!(irr, vec![poly(&[2, 0, 1])]);
        // a repeated factor is reported once (square-free support).
        let sq = monic_irreducible_factors(&poly(&[1, 2, 1])); // (x+1)²
        assert_eq!(sq, vec![poly(&[1, 1])]);
    }

    #[test]
    fn valuations_at_places() {
        // a = t / (t + 1)
        let a = rf(&[0, 1], &[1, 1]);
        assert_eq!(
            try_valuation_at_ff(&a, &FFPlace::Finite(poly(&[0, 1]))),
            Some(1)
        ); // at π = t
        assert_eq!(
            try_valuation_at_ff(&a, &FFPlace::Finite(poly(&[1, 1]))),
            Some(-1)
        ); // at π = t+1
        assert_eq!(try_valuation_at_ff(&a, &FFPlace::Infinite), Some(0)); // deg den − deg num = 0
                                                                          // 1/t² has a double pole at ∞? no: v_∞(1/t²) = deg(t²) − deg(1) = 2.
        assert_eq!(
            try_valuation_at_ff(&rf(&[1], &[0, 0, 1]), &FFPlace::Infinite),
            Some(2)
        );
    }

    #[test]
    fn residue_field_order_overflow_returns_none() {
        let pi = Poly::<Fp<5>>::monomial(56, Fp::<5>::one()).add(&Poly::one());
        assert_eq!(try_kappa_order(&FFPlace::Finite(pi)), None);
    }

    #[test]
    fn hilbert_symbol_is_symmetric_and_steinberg() {
        let samples = [
            rf(&[0, 1], &[1]),
            rf(&[2], &[1]),
            rf(&[1, 1], &[1]),
            rf(&[0, 1], &[1, 1]),
        ];
        let places = [
            FFPlace::Infinite,
            FFPlace::Finite(poly(&[0, 1])),    // t
            FFPlace::Finite(poly(&[1, 1])),    // t+1
            FFPlace::Finite(poly(&[2, 0, 1])), // t²+2 (degree-2 place, residue F_25)
        ];
        for a in &samples {
            for b in &samples {
                for pl in &places {
                    assert_eq!(
                        try_hilbert_symbol_ff(a, b, pl),
                        try_hilbert_symbol_ff(b, a, pl),
                        "symmetry"
                    );
                }
                // Steinberg: (a, −a)_v = 1.
                let neg_a = a.mul(&F::from_base(Fp::<5>::new(-1)));
                for pl in &places {
                    assert_eq!(
                        try_hilbert_symbol_ff(a, &neg_a, pl),
                        Some(1),
                        "(a,−a)_v = 1"
                    );
                }
            }
        }
    }

    #[test]
    fn reciprocity_holds_small() {
        // ∏_v (a,b)_v = +1 — the gold oracle, exact over F_q(t).
        let samples = [
            rf(&[0, 1], &[1]),    // t
            rf(&[1, 1], &[1]),    // t+1
            rf(&[2], &[1]),       // the nonsquare constant 2
            rf(&[0, 1], &[1, 1]), // t/(t+1)
            rf(&[2, 0, 1], &[1]), // t²+2 (irreducible)
        ];
        for a in &samples {
            for b in &samples {
                assert_eq!(
                    try_hilbert_reciprocity_product_ff(a, b),
                    Some(1),
                    "reciprocity failed at a={a:?} b={b:?}"
                );
            }
        }
    }

    #[test]
    fn quaternion_ramifies_at_an_even_number_of_places() {
        // The function-field mirror of "Hamilton's quaternions ramify at 2 and ∞":
        // (t, 2) over F_5(t) — with 2 a nonsquare constant — ramifies at exactly the
        // linear place π = t and the degree place ∞ (an even set, as reciprocity
        // forces).
        let a = rf(&[0, 1], &[1]); // t
        let b = rf(&[2], &[1]); // the nonsquare 2
        let ram = try_ramified_places_ff(&a, &b).unwrap();
        assert_eq!(ram.len(), 2, "even number of ramified places");
        assert!(ram.contains(&FFPlace::Finite(poly(&[0, 1])))); // π = t
        assert!(ram.contains(&FFPlace::Infinite)); // ∞
                                                   // a split quaternion (a square second slot) ramifies nowhere.
        assert!(try_ramified_places_ff(&a, &rf(&[4], &[1]))
            .unwrap()
            .is_empty()); // 4 = 2² is a square
    }

    #[test]
    fn hasse_minkowski_global_matches_adelic() {
        // For rank ≥ 3, try_is_isotropic_ff agrees with the per-place adelic breakdown.
        let forms: [Vec<F>; 4] = [
            vec![rf(&[1], &[1]), rf(&[1], &[1]), rf(&[4], &[1])], // ⟨1,1,−1⟩ isotropic
            vec![rf(&[1], &[1]), rf(&[0, 1], &[1]), rf(&[0, 4], &[1])], // ⟨1,t,−t⟩ isotropic
            // ⟨1,−t,−2,2t⟩ = norm form of the division quaternion (t,2): anisotropic
            vec![
                rf(&[1], &[1]),
                rf(&[0, 4], &[1]),
                rf(&[3], &[1]),
                rf(&[0, 2], &[1]),
            ],
            // …extend by ⟨…,1⟩ to rank 5: u-invariant 4 ⇒ isotropic
            vec![
                rf(&[1], &[1]),
                rf(&[0, 4], &[1]),
                rf(&[3], &[1]),
                rf(&[0, 2], &[1]),
                rf(&[1], &[1]),
            ],
        ];
        let expected = [true, true, false, true];
        for (form, &exp) in forms.iter().zip(&expected) {
            assert_eq!(
                try_is_isotropic_ff(form),
                Some(exp),
                "global isotropy of {form:?}"
            );
            assert_eq!(
                try_isotropy_over_ff_adeles(form).unwrap().is_global(),
                exp,
                "adelic isotropy of {form:?}"
            );
        }
    }

    #[test]
    fn cross_checks_springer_laurent_at_a_linear_place() {
        // Independent-oracle check: at a linear place π = t − 1 (residue field F_5),
        // the F_q(t) place layer must read the same valuation/discriminant data as
        // the established Laurent Springer decomposition. Entries are expanded at
        // t = 1 by hand for the Laurent side and kept as polynomials for the F_q(t)
        // side; the two code paths are completely independent.
        use crate::clifford::Metric;
        use crate::forms::springer_decompose_laurent;
        use crate::scalar::Laurent;
        type L5 = Laurent<Fp<5>, 4>;
        let lc = |cs: &[i128], v: i128| {
            L5::from_coeffs(cs.iter().map(|&n| Fp::<5>::new(n)).collect(), v)
        };

        let pi = poly(&[4, 1]); // t − 1 over F_5
        let place = FFPlace::Finite(pi.clone());

        // ⟨ 2,  t−1,  t²+1 ⟩  — valuations 0,1,0 at π; residues 2 (nonsq), 1, 2.
        let ff = [
            rf(&[2], &[1]),       // 2
            rf(&[4, 1], &[1]),    // t − 1
            rf(&[1, 0, 1], &[1]), // t² + 1
        ];
        // the same form localized at s = t − 1: 2 = 2·s⁰; t−1 = s¹; t²+1 = s²+2s+2.
        let laurent = Metric::diagonal(vec![lc(&[2], 0), lc(&[1], 1), lc(&[2, 2, 1], 0)]);
        let decomp = springer_decompose_laurent(&laurent).unwrap();

        for layer in &decomp.graded {
            let at: Vec<&F> = ff
                .iter()
                .filter(|e| try_valuation_at_ff(e, &place) == Some(layer.valuation))
                .collect();
            assert_eq!(at.len(), layer.dim, "dim at valuation {}", layer.valuation);
            // discriminant square class = XNOR of the per-entry residue characters.
            let disc_sq = at.iter().fold(true, |acc, e| {
                acc == (try_chi_kappa(&try_residue_unit_at(e, &place).unwrap(), &place) == Some(1))
            });
            assert_eq!(
                disc_sq, layer.disc_is_square,
                "disc class at valuation {}",
                layer.valuation
            );
        }
    }

    #[test]
    fn rank_two_is_a_global_square_condition() {
        // ⟨1, −t²⟩: −(1·−t²) = t² is a global square ⇒ isotropic.
        assert_eq!(
            try_is_isotropic_ff(&[rf(&[1], &[1]), rf(&[0, 0, 4], &[1])]),
            Some(true)
        );
        // ⟨1, −t⟩: −(1·−t) = t has an odd-multiplicity place ⇒ not a square ⇒ anisotropic.
        assert_eq!(
            try_is_isotropic_ff(&[rf(&[1], &[1]), rf(&[0, 4], &[1])]),
            Some(false)
        );
        // ⟨1, −2⟩: −(−2) = 2 is a nonsquare constant ⇒ anisotropic.
        assert_eq!(
            try_is_isotropic_ff(&[rf(&[1], &[1]), rf(&[3], &[1])]),
            Some(false)
        );
        // ⟨2, −8⟩: −(2·−8) = 16 = 4² is a square ⇒ isotropic.
        assert_eq!(
            try_is_isotropic_ff(&[rf(&[2], &[1]), rf(&[2], &[1])]),
            Some(true)
        ); // 2 and 2: −(4)=−4≡1=1²
    }
}
