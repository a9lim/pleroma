//! Global isotropy over `F_q(t)` (the Hasse–Minkowski verdict in characteristic 2).
//!
//! This module contains only the globally-scoped functions; it calls back into the
//! local engine (`local_is_isotropic_char2`) through the public API in the parent
//! module hub so that no visibility is widened.

use super::{local_is_isotropic_char2, Char2QuadForm};
use crate::forms::function_field_char2::char2_monic_irreducible_factors;
use crate::forms::{Char2Place, FiniteChar2Field};
use crate::scalar::{Poly, RationalFunction, Scalar};

/// Whether `f ∈ F_q(t)` is a **square**, i.e. lies in `K² = F_q(t²)`. Since
/// `[F_q(t) : F_q(t)²] = 2` (`F_q` perfect, basis `{1, t}` over `K²`), `f = N/D` is a
/// square iff `N·D ∈ F_q[t]²`, and a char-2 polynomial over a perfect field is a
/// square iff every **odd-degree** coefficient vanishes (the even ones are squares
/// automatically). The additive `℘`-analogue of this is [`global_is_pe`].
pub fn ff_is_square<S: FiniteChar2Field>(f: &RationalFunction<S>) -> bool {
    if f.is_zero() {
        return true;
    }
    let prod = f.num().mul(f.den());
    prod.coeffs()
        .iter()
        .enumerate()
        .all(|(i, c)| i & 1 == 0 || c.is_zero())
}

/// Whether `f ∈ ℘(F_q(t))` — the **global** Artin–Schreier triviality test
/// (`℘(x) = x² + x`). By the local–global principle for `℘` over the rational
/// function field, `f ∈ ℘(F_q(t))` iff `f ∈ ℘(K_v)` at **every** place; and the only
/// places that can carry an obstruction are the poles of `f` (finite places dividing
/// `den f`) and `∞` (which also sees the leftover constant's `Tr_{F_q/F₂}`). So a
/// finite sweep of `{∞} ∪ {P | den f}` decides it. The additive analogue of the
/// odd-char `is_global_square_ff`.
pub fn global_is_pe<S: FiniteChar2Field>(f: &RationalFunction<S>) -> bool {
    use super::asnf::local_is_pe;
    if f.is_zero() {
        return true;
    }
    if !local_is_pe(f, &Char2Place::Infinite) {
        return false;
    }
    char2_monic_irreducible_factors(f.den())
        .into_iter()
        .all(|p| local_is_pe(f, &Char2Place::Finite(p)))
}

/// The finite set of places of `F_q(t)` that can make `form` anisotropic: `∞` plus
/// every monic irreducible dividing a numerator or denominator of some coefficient.
/// At every **other** place all coefficients are units, so a rank-`≥ 3` form reduces
/// to a `> 2`-variable form over the finite residue field `κ` — isotropic by
/// Chevalley–Warning and liftable by Hensel — and need not be checked.
pub fn relevant_places_char2<S: FiniteChar2Field>(form: &Char2QuadForm<S>) -> Vec<Char2Place<S>> {
    let mut primes: Vec<Poly<S>> = Vec::new();
    let mut push = |g: &Poly<S>| {
        for p in char2_monic_irreducible_factors(g) {
            if !primes.contains(&p) {
                primes.push(p);
            }
        }
    };
    for (a, b) in &form.blocks {
        push(a.num());
        push(a.den());
        push(b.num());
        push(b.den());
    }
    for c in &form.singular {
        push(c.num());
        push(c.den());
    }
    let mut places = vec![Char2Place::Infinite];
    places.extend(primes.into_iter().map(Char2Place::Finite));
    places
}

/// Whether `form` is **isotropic over `F_q(t)`** (the global Hasse–Minkowski verdict
/// in characteristic 2). The dispatch, all source-pinned (Aravire–Jacob;
/// Elman–Karpenko–Merkurjev; Csahók–Kutas–Montessinos–Zábrádi; Tsen–Lang `C₂`):
///
/// * a null coefficient (`⟨0⟩`) or a hyperbolic block (`[0,b]`/`[a,0]`) ⇒ isotropic;
/// * `rank ≥ 5` ⇒ isotropic — `u(F_q(t)) = 4` (`F_q(t)` is a `C₂` field);
/// * **totally singular** part (`℘`-free, quasilinear): `[K : K²] = 2`, so `≥ 3`
///   singular entries are isotropic, and a binary `⟨c₁, c₂⟩` is isotropic iff
///   `c₁c₂ ∈ K²` (`ff_is_square`); an anisotropic binary quasilinear part is
///   *universal*, so it isotropises any form carrying a nonzero block;
/// * **rank 2** `[a, b]`: isotropic iff `ab ∈ ℘(F_q(t))` ([`global_is_pe`]) — *not* a
///   finite bad-place sweep, since the constant-trace obstruction lives at infinitely
///   many odd-degree places (caught by the global `℘` test);
/// * **rank 3/4 non-degenerate**: Hasse–Minkowski — isotropic iff isotropic over
///   `K_v` at every [`relevant_places_char2`] (a finite set).
///
/// The return type stays optional for API symmetry with the local routines; the
/// current local engine covers every shape routed here.
pub fn is_isotropic_global_char2<S: FiniteChar2Field>(form: &Char2QuadForm<S>) -> Option<bool> {
    // A null direction or a hyperbolic block isotropises the whole form.
    if form.singular.iter().any(|c| c.is_zero()) {
        return Some(true);
    }
    if form.blocks.iter().any(|(a, b)| a.is_zero() || b.is_zero()) {
        return Some(true);
    }
    let nb = form.blocks.len();
    let ns = form.singular.len();
    let rank = 2 * nb + ns;
    if rank == 0 {
        return Some(false); // the empty form is anisotropic by convention
    }
    if rank >= 5 {
        return Some(true); // u(F_q(t)) = 4
    }
    // Totally-singular handling (the quasilinear part), elementary over F_q(t).
    if ns >= 3 {
        return Some(true); // ≥ 3 entries are K²-dependent
    }
    if ns == 2 {
        // A binary block present ⇒ the (universal-if-anisotropic) singular pair
        // isotropises it; otherwise it is the pure quasilinear ⟨c₁,c₂⟩.
        if nb >= 1 {
            return Some(true);
        }
        let prod = form.singular[0].mul(&form.singular[1]);
        return Some(ff_is_square(&prod)); // ⟨c₁,c₂⟩ iso ⟺ c₁c₂ ∈ K²
    }
    // Non-degenerate from here (#singular ≤ 1).
    match (nb, ns) {
        (0, 1) => Some(false), // ⟨c⟩, c ≠ 0
        (1, 0) => {
            // rank 2: [a,b] isotropic ⟺ ab ∈ ℘(F_q(t)).
            let (a, b) = &form.blocks[0];
            Some(global_is_pe(&a.mul(b)))
        }
        // rank 3 ([a,b]⊥⟨c⟩) and rank 4 ([a,b]⊥[a,b]): Hasse–Minkowski.
        _ => {
            let mut all_iso = true;
            for place in relevant_places_char2(form) {
                match local_is_isotropic_char2(form, &place) {
                    Some(true) => {}
                    Some(false) => return Some(false),
                    None => all_iso = false,
                }
            }
            all_iso.then_some(true)
        }
    }
}
