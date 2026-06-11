//! The **characteristic-2 local Witt/Springer decomposition** over `F_q(t)` — the
//! Aravire–Jacob structure theorem, and the rank-by-rank local isotropy it powers.
//! This is the equal-characteristic-2 mirror of the odd-`q`
//! [`springer_decompose_local`](crate::forms::springer_decompose_local) engine, but it is **not** the
//! odd story at `p = 2`: in char 2 a third, *wild* summand appears that the
//! `W = W(k) ⊕ W(k)` grading misses (see root AGENTS.md "the wild-term finding").
//!
//! # The decomposition
//!
//! At a place `v` of `F_q(t)` (completion `K_v = κ((π))`, residue field `κ`,
//! canonical uniformizer `π = P` at a finite place, `π = 1/t` at `∞`), every
//! nonsingular quadratic form `φ = ⊥ [a_i, b_i]` (`[a,b] = ax²+xy+by²`) decomposes
//! **in the Witt group** `W_q(K_v)` as
//!
//! ```text
//! φ  ≅  φ₀  ⊥  ψ  ⊥  ⟨π⟩·φ₁,        φ₀, φ₁ ∈ W_q(κ),  ψ ∈ R_π,
//! ```
//!
//! the Aravire–Jacob theorem (*Quadratic forms over `k(x)` in characteristic 2*,
//! Thm 1.3 / Cor 1.7). Here `W_q(κ) ≅ F₂` via the Arf class `Tr_{κ/F₂}`, and the
//! **wild part** lives in
//!
//! ```text
//! R_π  =  ⊕_{n odd ≥ 1}  κ·π⁻ⁿ      (binary forms [1, r], r a sum of odd negative π-powers),
//! ```
//!
//! a coordinate the odd-characteristic theory has no analogue of: a pole-order-`m`
//! element has `℘`-image of pole order `2m` (even), so **odd-order poles survive
//! `℘`** and cannot be absorbed into `φ₀`/`φ₁`. The decomposition depends on the
//! choice of uniformizer (AJ Rem 1.8); we always use the canonical `π` above.
//!
//! ## The algorithm (`K = K² ⊕ πK²`, then AS-normal-form)
//!
//! Splitting a coefficient by Laurent-exponent parity, `a = a_ev + a_odd` with
//! `a_ev ∈ K²`, `a_odd ∈ πK²` (`κ` perfect ⇒ even-power series are squares), the
//! fundamental relation `[a,b] ≅ [1, a_ev·b] ⊥ ⟨π⟩[1, a_odd·b]` reduces each block
//! to two `[1, c]`'s. Each `c` is pushed to its **Artin–Schreier normal form**: drop
//! the (always-`℘`) positive-power tail, and clear even negative poles bottom-up via
//! `c_n ↦ 0`, `c_{n/2} += √c_n` (subtracting `℘(√c_n·π^{n/2})`). What remains is a
//! `κ`-constant (its `Tr_{κ/F₂}` is the `W_q(κ)` bit) plus odd negative poles (the
//! `R_π` coordinate). The `a_odd` half's wild part folds into `ψ` too, since
//! `⟨π⟩[1, r] ≅ [1, r]` for `r ∈ R_π ⊂ πK²`.
//!
//! # Isotropy (rank-by-rank, `u(K_v) = 4`)
//!
//! `[a,b]` is isotropic iff `ab ∈ ℘(K_v)`; the Pfister/norm criterion routes ranks 3
//! and 4 through the Part-A Artin–Schreier symbol [`as_symbol_at`]
//! (`s_v(d, λ) = 0` ⟺ `[d, λ)` splits); `u(K_v) = 4` makes every rank `≥ 5` form
//! isotropic. (Source-pinned to Aravire–Jacob and Elman–Karpenko–Merkurjev §§7, 13;
//! oracles cross-checked via Codex — see the tests.)
//!
//! # Global isotropy over `F_q(t)` (Hasse–Minkowski)
//!
//! [`is_isotropic_global_char2`] decides isotropy over `F_q(t)` itself. Three
//! ingredients beyond the per-place symbol, all source-pinned: the global
//! `℘`-membership test [`global_is_pe`] (`f ∈ ℘(F_q(t))` ⟺ `f ∈ ℘(K_v)` everywhere,
//! a finite sweep of the poles of `f` plus `∞`) settles **rank 2** (`[a,b]` iso ⟺
//! `ab ∈ ℘`); the elementary `[K : K²] = 2` square test [`ff_is_square`] settles the
//! **totally-singular** part; and **Hasse–Minkowski** over the finite
//! [`relevant_places_char2`] set settles **rank 3/4** (good places are isotropic by
//! Chevalley–Warning). `u(F_q(t)) = 4` (`C₂`, Tsen–Lang) caps it: every `rank ≥ 5`
//! form is isotropic. (Local isotropy itself is reported for the ranks the sources
//! pin exactly — `≤ 4` in the standard block shapes, pure totally-singular tails
//! via the two square classes of `K_v`, the rank-4 mixed case consisting of one
//! binary block plus a two-class quasilinear tail, all nonsingular ranks via the AJ
//! kernel, one-class singular tails via the odd-dimensional Clifford invariant, and
//! `≥ 5` always isotropic.)
//!
//! # Module layout
//!
//! - `asnf` — κ-local arithmetic helpers and the Artin–Schreier normal form
//!   (the private crate layer that feeds the decomposition).
//! - `global` — global isotropy over `F_q(t)` ([`global_is_pe`], [`ff_is_square`],
//!   [`relevant_places_char2`], [`is_isotropic_global_char2`]).
//! - This hub — `Char2QuadForm`, `Char2LocalDecomp`, the Aravire–Jacob decomposition
//!   ([`springer_decompose_local_char2`]), and local isotropy
//!   ([`local_anisotropic_dim_char2`], [`local_is_isotropic_char2`]).

pub(super) mod asnf;
mod global;

pub use global::{ff_is_square, global_is_pe, is_isotropic_global_char2, relevant_places_char2};

use crate::forms::{as_symbol_at, Char2Place, FiniteChar2Field};
use crate::scalar::{Poly, RationalFunction, Scalar};
use std::collections::BTreeMap;

use asnf::{asnf, kmul, laurent, local_is_pe, local_is_square, merge_psi, valuation};

/// A characteristic-2 quadratic form over `F_q(t)`: a sum of nonsingular binary
/// blocks `[a_i, b_i] = a_i x² + xy + b_i y²` and a totally-singular part
/// `⟨c_j⟩ = Σ c_j x_j²` (the radical of the polar form). Rank `= 2·blocks + singular`.
#[derive(Debug, Clone, PartialEq)]
pub struct Char2QuadForm<S: FiniteChar2Field> {
    /// The nonsingular binary blocks `[a_i, b_i]`.
    pub blocks: Vec<(RationalFunction<S>, RationalFunction<S>)>,
    /// The totally-singular diagonal entries `⟨c_j⟩` (each `c_j ≠ 0`).
    pub singular: Vec<RationalFunction<S>>,
}

impl<S: FiniteChar2Field> Char2QuadForm<S> {
    /// A form from binary blocks alone (the generic nonsingular case).
    pub fn from_blocks(blocks: Vec<(RationalFunction<S>, RationalFunction<S>)>) -> Self {
        Self {
            blocks,
            singular: Vec::new(),
        }
    }

    /// A form from binary blocks plus a totally-singular tail.
    pub fn new(
        blocks: Vec<(RationalFunction<S>, RationalFunction<S>)>,
        singular: Vec<RationalFunction<S>>,
    ) -> Self {
        Self { blocks, singular }
    }

    /// The dimension (rank) of the form.
    pub fn rank(&self) -> usize {
        2 * self.blocks.len() + self.singular.len()
    }
}

/// The Aravire–Jacob local decomposition `(φ₀, ψ, φ₁)` of a nonsingular form at a
/// place, in the Witt group `W_q(K_v)`. `phi0`/`phi1 ∈ {0,1}` are the `W_q(κ) ≅ F₂`
/// Arf bits of the unramified and `⟨π⟩`-scaled parts; `psi` is the wild `R_π`
/// coordinate, a sparse map `odd pole order n ↦ κ-coefficient` (no zero entries).
#[derive(Debug, Clone, PartialEq)]
pub struct Char2LocalDecomp<S: FiniteChar2Field> {
    /// The `W_q(κ) ≅ F₂` Arf bit of the unramified part `φ₀`.
    pub phi0: u128,
    /// The wild part `ψ ∈ R_π`: `odd pole order ↦ κ-coefficient` (`κ = F_q[t]/(P)`,
    /// or `F_q` at `∞`, stored as a reduced `Poly<S>`).
    pub psi: BTreeMap<usize, Poly<S>>,
    /// The `W_q(κ) ≅ F₂` Arf bit of the `⟨π⟩`-scaled part `φ₁`.
    pub phi1: u128,
}

// ───────────────────────── the decomposition ─────────────────────────

/// The per-block contribution `(φ₀-bit, φ₁-bit, ψ-part)` of a nonsingular block
/// `[a, b]` (`a, b ≠ 0`) to the Aravire–Jacob decomposition at `place`.
fn block_contribution<S: FiniteChar2Field>(
    a: &RationalFunction<S>,
    b: &RationalFunction<S>,
    place: &Char2Place<S>,
) -> (u128, u128, BTreeMap<usize, Poly<S>>) {
    let va = valuation(a, place).expect("a ≠ 0");
    let vb = valuation(b, place).expect("b ≠ 0");
    // Expand a over [va, max(0,-vb)] and b over [vb, max(0,-va)]: this captures every
    // coefficient of the products a_•·b at exponents ≤ 0 (the only ones ASNF reads).
    let a_hi = std::cmp::max(0, -vb);
    let b_hi = std::cmp::max(0, -va);
    let acoeffs = laurent(a, place, va, a_hi);
    let bcoeffs = laurent(b, place, vb, b_hi);
    let plo = va + vb;
    let mut c_ev: BTreeMap<i128, Poly<S>> = BTreeMap::new();
    let mut c_odd: BTreeMap<i128, Poly<S>> = BTreeMap::new();
    let mut n = plo;
    while n <= 0 {
        let mut sev = Poly::<S>::zero();
        let mut sod = Poly::<S>::zero();
        let i_lo = std::cmp::max(va, n - b_hi);
        let i_hi = std::cmp::min(a_hi, n - vb);
        let mut i = i_lo;
        while i <= i_hi {
            let ai = &acoeffs[(i - va) as usize];
            if !ai.is_zero() {
                let bj = &bcoeffs[(n - i - vb) as usize];
                if !bj.is_zero() {
                    let prod = kmul(ai, bj, place);
                    if i & 1 == 0 {
                        sev = sev.add(&prod);
                    } else {
                        sod = sod.add(&prod);
                    }
                }
            }
            i += 1;
        }
        if !sev.is_zero() {
            c_ev.insert(n, sev);
        }
        if !sod.is_zero() {
            c_odd.insert(n, sod);
        }
        n += 1;
    }
    let (e0, r0) = asnf(&c_ev, plo, place);
    let (e1, r1) = asnf(&c_odd, plo, place);
    let mut psi = r0;
    for (k, v) in r1 {
        merge_psi(&mut psi, k, v);
    }
    (e0, e1, psi)
}

/// The Aravire–Jacob local decomposition `(φ₀, ψ, φ₁)` of `form`'s nonsingular part
/// at `place`, in `W_q(K_v)`. Hyperbolic blocks (`a = 0`, `b = 0`) drop out; the
/// totally-singular part is *not* part of this Witt-group invariant.
pub fn springer_decompose_local_char2<S: FiniteChar2Field>(
    form: &Char2QuadForm<S>,
    place: &Char2Place<S>,
) -> Char2LocalDecomp<S> {
    let mut phi0 = 0u128;
    let mut phi1 = 0u128;
    let mut psi: BTreeMap<usize, Poly<S>> = BTreeMap::new();
    for (a, b) in &form.blocks {
        if a.is_zero() || b.is_zero() {
            continue; // [0,b] / [a,0] is hyperbolic — trivial in W_q
        }
        let (e0, e1, part) = block_contribution(a, b, place);
        phi0 ^= e0;
        phi1 ^= e1;
        for (k, v) in part {
            merge_psi(&mut psi, k, v);
        }
    }
    Char2LocalDecomp { phi0, psi, phi1 }
}

// ───────────────────────── local isotropy (rank-by-rank) ─────────────────────────

/// Whether the binary block `[a, b]` is hyperbolic over `K_v` (isotropic, i.e.
/// `ab ∈ ℘(K_v)`, including the degenerate `a = 0` or `b = 0`).
fn binary_is_hyperbolic<S: FiniteChar2Field>(
    a: &RationalFunction<S>,
    b: &RationalFunction<S>,
    place: &Char2Place<S>,
) -> bool {
    if a.is_zero() || b.is_zero() {
        return true;
    }
    local_is_pe(&a.mul(b), place)
}

fn nonsingular_anisotropic_dim<S: FiniteChar2Field>(
    blocks: &[(RationalFunction<S>, RationalFunction<S>)],
    place: &Char2Place<S>,
) -> usize {
    let form = Char2QuadForm::from_blocks(blocks.to_vec());
    let d = springer_decompose_local_char2(&form, place);
    if d.phi0 == 0 && d.phi1 == 0 && d.psi.is_empty() {
        0
    } else if d.phi0 == 1 && d.phi1 == 1 && d.psi.is_empty() {
        4
    } else {
        2
    }
}

fn singular_anisotropic_dim<S: FiniteChar2Field>(
    singular: &[RationalFunction<S>],
    place: &Char2Place<S>,
) -> usize {
    singular_square_representatives(singular, place).len()
}

fn singular_square_representatives<S: FiniteChar2Field>(
    singular: &[RationalFunction<S>],
    place: &Char2Place<S>,
) -> Vec<RationalFunction<S>> {
    let mut reps = Vec::new();
    for c in singular.iter().filter(|c| !c.is_zero()) {
        if reps.is_empty() {
            reps.push(c.clone());
        } else if reps.len() == 1 {
            let ratio = c.mul(&reps[0].inv().expect("nonzero representative inverts"));
            if !local_is_square(&ratio, place) {
                reps.push(c.clone());
            }
        } else {
            break;
        }
    }
    reps
}

/// The local Clifford invariant of the odd-dimensional form
/// `⟨c⟩ ⊥ r`, where `r = ⊥[a_i,b_i]` is nonsingular. In characteristic 2,
/// `clif(⟨c⟩ ⊥ r) = clif(c r)`; with `[a_i,b_i] ≅ a_i[1,a_i b_i]`, this is the
/// Brauer sum `Σ [a_i b_i, c/a_i)` evaluated by the local Artin-Schreier symbol.
fn semisingular_clifford_at<S: FiniteChar2Field>(
    blocks: &[(RationalFunction<S>, RationalFunction<S>)],
    c: &RationalFunction<S>,
    place: &Char2Place<S>,
) -> u128 {
    let mut inv = 0u128;
    for (a, b) in blocks {
        if a.is_zero() || b.is_zero() {
            continue;
        }
        let d = a.mul(b);
        if local_is_pe(&d, place) {
            continue;
        }
        let lambda = c.mul(&a.inv().expect("a ≠ 0"));
        inv ^= as_symbol_at(&d, &lambda, place);
    }
    inv
}

fn semisingular_anisotropic_dim<S: FiniteChar2Field>(
    blocks: &[(RationalFunction<S>, RationalFunction<S>)],
    c: &RationalFunction<S>,
    place: &Char2Place<S>,
) -> usize {
    if semisingular_clifford_at(blocks, c, place) == 0 {
        1
    } else {
        3
    }
}

/// The **local anisotropic dimension** of `form` over `K_v` at `place`, for the form
/// shapes the char-2 theory pins exactly: pure totally-singular forms (read as a
/// `K_v²`-span inside the two-dimensional vector space `K_v = K_v² ⊕ πK_v²`),
/// nonsingular forms of any rank via the AJ kernel, one-class singular tails via
/// the odd-dimensional Clifford invariant, and any binary-block part plus a
/// two-dimensional singular tail. `u(K_v) = 4`.
pub fn local_anisotropic_dim_char2<S: FiniteChar2Field>(
    form: &Char2QuadForm<S>,
    place: &Char2Place<S>,
) -> Option<usize> {
    let bl = &form.blocks;
    let nb = bl.len();
    let singular = singular_square_representatives(&form.singular, place);
    let ns = singular.len();
    let rank = 2 * nb + ns;
    if rank == 0 {
        return Some(0);
    }
    if nb == 0 {
        return Some(singular_anisotropic_dim(&form.singular, place));
    }
    let reduced_blocks: Vec<_> = bl
        .iter()
        .filter(|(a, b)| !binary_is_hyperbolic(a, b, place))
        .cloned()
        .collect();
    if reduced_blocks.len() != nb {
        let reduced = Char2QuadForm::new(reduced_blocks, singular);
        return local_anisotropic_dim_char2(&reduced, place);
    }
    if ns == 0 {
        return Some(nonsingular_anisotropic_dim(bl, place));
    }
    if ns == 2 {
        // The two coefficients form a K²-basis of K. The quasilinear tail is
        // universal as a value set, so it splits off hyperbolic planes from any
        // nonsingular part and leaves precisely its two-dimensional radical kernel.
        return Some(2);
    }
    if ns == 1 {
        return Some(semisingular_anisotropic_dim(bl, &singular[0], place));
    }
    unreachable!("K_v has K_v²-dimension two")
}

/// Whether `form` is **isotropic over the completion** `K_v` at `place`. `Some(true)`
/// for every rank `≥ 5` (`u(K_v) = 4`); otherwise `anisotropic_dim < rank`.
pub fn local_is_isotropic_char2<S: FiniteChar2Field>(
    form: &Char2QuadForm<S>,
    place: &Char2Place<S>,
) -> Option<bool> {
    let rank = form.rank();
    if rank == 0 {
        return Some(false);
    }
    if rank >= 5 {
        return Some(true);
    }
    local_anisotropic_dim_char2(form, place).map(|d| d < rank)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scalar::{Fp, Fpn};

    type F2 = Fp<2>;
    type R2 = RationalFunction<F2>;

    fn r2(num: &[i128], den: &[i128]) -> R2 {
        RationalFunction::new(
            num.iter().map(|&n| F2::new(n)).collect(),
            den.iter().map(|&n| F2::new(n)).collect(),
        )
    }
    fn p2(c: &[i128]) -> Poly<F2> {
        Poly::new(c.iter().map(|&n| F2::new(n)).collect())
    }
    // The place P = t over F₂: κ = F₂, π = t. Codex's oracles are stated over the
    // local field F_q((π)); identifying π = t realises each as a form over F_q(t).
    fn place_t() -> Char2Place<F2> {
        Char2Place::Finite(p2(&[0, 1]))
    }
    // A constant-coefficient κ entry (κ = F₂ at place t) for the R_π map.
    fn k2(n: i128) -> Poly<F2> {
        Poly::constant(F2::new(n))
    }
    fn decomp(blocks: &[(R2, R2)]) -> Char2LocalDecomp<F2> {
        springer_decompose_local_char2(&Char2QuadForm::from_blocks(blocks.to_vec()), &place_t())
    }

    // ── the Aravire–Jacob decomposition, against Codex's source-pinned oracles ──
    // (π = t; "π⁻¹" ↦ R_π map {1: 1}, "1" in φ₀/φ₁ ↦ the bit 1.)

    #[test]
    fn aj_oracle_perp_killed_by_p() {
        // [1, π⁻²+π⁻¹] ↦ (0,0,0): π⁻²+π⁻¹ = ℘(π⁻¹) ∈ ℘(K), so the block is hyperbolic.
        let d = decomp(&[(r2(&[1], &[1]), r2(&[1, 1], &[0, 0, 1]))]); // b = (1+t)/t²
        assert_eq!(d.phi0, 0);
        assert_eq!(d.phi1, 0);
        assert!(d.psi.is_empty());
    }

    #[test]
    fn aj_oracle_wild_pole() {
        // [1, π⁻²] ↦ (0, π⁻¹, 0): the even pole reduces to the odd (wild) pole π⁻¹.
        let d = decomp(&[(r2(&[1], &[1]), r2(&[1], &[0, 0, 1]))]); // b = 1/t²
        assert_eq!(d.phi0, 0);
        assert_eq!(d.phi1, 0);
        assert_eq!(d.psi, BTreeMap::from([(1usize, k2(1))]));
    }

    #[test]
    fn aj_oracle_residue_bit() {
        // [1+π, 1] ↦ (1,0,0): a_ev = 1 gives Arf bit 1; a_odd = π·b = π is positive (℘).
        let d = decomp(&[(r2(&[1, 1], &[1]), r2(&[1], &[1]))]);
        assert_eq!(d.phi0, 1);
        assert_eq!(d.phi1, 0);
        assert!(d.psi.is_empty());
    }

    #[test]
    fn aj_oracle_split_across_phi1_and_psi() {
        // [1+π, π⁻¹] ↦ (0, π⁻¹, 1): a_ev·b = π⁻¹ (wild), a_odd·b = 1 (φ₁ bit).
        let d = decomp(&[(r2(&[1, 1], &[1]), r2(&[1], &[0, 1]))]); // b = 1/t
        assert_eq!(d.phi0, 0);
        assert_eq!(d.phi1, 1);
        assert_eq!(d.psi, BTreeMap::from([(1usize, k2(1))]));
    }

    #[test]
    fn aj_oracle_sum_of_two_residue_bits_cancels() {
        // [1,1] ⊥ [1,1] ↦ (0,0,0): each [1,1] over F₂ is (1,0,0); the bits cancel.
        let one = r2(&[1], &[1]);
        let single = decomp(&[(one.clone(), one.clone())]);
        assert_eq!((single.phi0, single.phi1, single.psi.len()), (1, 0, 0));
        let d = decomp(&[(one.clone(), one.clone()), (one.clone(), one.clone())]);
        assert_eq!(d.phi0, 0);
        assert_eq!(d.phi1, 0);
        assert!(d.psi.is_empty());
    }

    #[test]
    fn aj_oracle_anisotropic_u4() {
        // [1,1] ⊥ [π,π⁻¹] ↦ (1,0,1): the u = 4 anisotropic class.
        let d = decomp(&[
            (r2(&[1], &[1]), r2(&[1], &[1])),
            (r2(&[0, 1], &[1]), r2(&[1], &[0, 1])),
        ]);
        assert_eq!(d.phi0, 1);
        assert_eq!(d.phi1, 1);
        assert!(d.psi.is_empty());
    }

    #[test]
    fn aj_oracle_residue_plus_wild() {
        // [1, π⁻¹] ⊥ [1,1] ↦ (1, π⁻¹, 0).
        let d = decomp(&[
            (r2(&[1], &[1]), r2(&[1], &[0, 1])),
            (r2(&[1], &[1]), r2(&[1], &[1])),
        ]);
        assert_eq!(d.phi0, 1);
        assert_eq!(d.phi1, 0);
        assert_eq!(d.psi, BTreeMap::from([(1usize, k2(1))]));
    }

    // ── local isotropy, against Codex's rank-by-rank oracles (π = t) ──

    fn form(blocks: &[(R2, R2)], singular: &[R2]) -> Char2QuadForm<F2> {
        Char2QuadForm::new(blocks.to_vec(), singular.to_vec())
    }
    fn anis(blocks: &[(R2, R2)], singular: &[R2]) -> Option<usize> {
        local_anisotropic_dim_char2(&form(blocks, singular), &place_t())
    }
    fn iso(blocks: &[(R2, R2)], singular: &[R2]) -> Option<bool> {
        local_is_isotropic_char2(&form(blocks, singular), &place_t())
    }

    #[test]
    fn iso_rank2() {
        // [1, π⁻²] anisotropic (wild pole); [1, π⁻²+π⁻¹] = ℘ hyperbolic.
        assert_eq!(
            anis(&[(r2(&[1], &[1]), r2(&[1], &[0, 0, 1]))], &[]),
            Some(2)
        );
        assert_eq!(
            iso(&[(r2(&[1], &[1]), r2(&[1], &[0, 0, 1]))], &[]),
            Some(false)
        );
        assert_eq!(
            anis(&[(r2(&[1], &[1]), r2(&[1, 1], &[0, 0, 1]))], &[]),
            Some(0)
        );
        assert_eq!(
            iso(&[(r2(&[1], &[1]), r2(&[1, 1], &[0, 0, 1]))], &[]),
            Some(true)
        );
    }

    #[test]
    fn iso_rank3() {
        // [1, π⁻¹] ⊥ ⟨1⟩: s(π⁻¹, 1) = 0 (1 is a square) ⇒ isotropic, anis dim 1.
        assert_eq!(
            anis(&[(r2(&[1], &[1]), r2(&[1], &[0, 1]))], &[r2(&[1], &[1])]),
            Some(1)
        );
        // [1,1] ⊥ ⟨π⟩: s(1, π) = 1 ⇒ anisotropic, anis dim 3.
        assert_eq!(
            anis(&[(r2(&[1], &[1]), r2(&[1], &[1]))], &[r2(&[0, 1], &[1])]),
            Some(3)
        );
        assert_eq!(
            iso(&[(r2(&[1], &[1]), r2(&[1], &[1]))], &[r2(&[0, 1], &[1])]),
            Some(false)
        );
    }

    #[test]
    fn pure_singular_local_dimension_reads_square_classes() {
        let one = r2(&[1], &[1]);
        let t = r2(&[0, 1], &[1]);
        let t2 = r2(&[0, 0, 1], &[1]);
        let one_plus_t = r2(&[1, 1], &[1]);
        // ⟨1,t⟩ spans both K²-classes over F₂((t)), so it is anisotropic of dim 2.
        assert_eq!(anis(&[], &[one.clone(), t.clone()]), Some(2));
        assert_eq!(iso(&[], &[one.clone(), t.clone()]), Some(false));
        // Same valuation parity is not enough: 1+t has an odd π-coefficient, so
        // it is not in K² and ⟨1,1+t⟩ is still anisotropic of dimension 2.
        assert_eq!(anis(&[], &[one.clone(), one_plus_t.clone()]), Some(2));
        assert_eq!(iso(&[], &[one.clone(), one_plus_t]), Some(false));
        // ⟨1,t²⟩ has both coefficients in the even-valuation square class.
        assert_eq!(anis(&[], &[one.clone(), t2.clone()]), Some(1));
        assert_eq!(iso(&[], &[one.clone(), t2.clone()]), Some(true));
        // Three pure singular entries are necessarily dependent over K².
        assert_eq!(anis(&[], &[one, t, t2]), Some(2));
        assert_eq!(
            iso(
                &[],
                &[r2(&[1], &[1]), r2(&[0, 1], &[1]), r2(&[0, 0, 1], &[1])]
            ),
            Some(true)
        );
    }

    #[test]
    fn mixed_singular_tail_collapses_to_local_square_classes() {
        let one = r2(&[1], &[1]);
        let t = r2(&[0, 1], &[1]);
        let t2 = r2(&[0, 0, 1], &[1]);
        let one_plus_t = r2(&[1, 1], &[1]);
        let block = [(one.clone(), one.clone())];

        // Two singular entries in the same K²-class reduce to the rank-3 case.
        assert_eq!(
            anis(&block, &[one.clone(), t2.clone()]),
            anis(&block, std::slice::from_ref(&one))
        );

        // Distinct K²-lines span K over K². The singular tail is a
        // universal quasilinear plane, so one binary block plus that tail is
        // isotropic with the two-class tail left as the anisotropic kernel.
        assert_eq!(anis(&block, &[one.clone(), t.clone()]), Some(2));
        assert_eq!(iso(&block, &[one.clone(), t]), Some(true));
        assert_eq!(anis(&block, &[one.clone(), one_plus_t.clone()]), Some(2));
        assert_eq!(iso(&block, &[one, one_plus_t]), Some(true));
    }

    #[test]
    fn iso_rank4() {
        let one = r2(&[1], &[1]);
        let t = r2(&[0, 1], &[1]);
        let inv_t = r2(&[1], &[0, 1]);
        // [1, π⁻¹] ⊥ [1,1]: Δ = π⁻¹+1 ∉ ℘ ⇒ isotropic, anis dim 2.
        assert_eq!(
            anis(
                &[(one.clone(), inv_t.clone()), (one.clone(), one.clone())],
                &[]
            ),
            Some(2)
        );
        // [1,1] ⊥ [1,1]: Δ = 0, s(1,1) = 0 ⇒ hyperbolic, anis dim 0.
        assert_eq!(
            anis(
                &[(one.clone(), one.clone()), (one.clone(), one.clone())],
                &[]
            ),
            Some(0)
        );
        // [1,1] ⊥ [π, π⁻¹]: Δ = 0, s(1, π) = 1 ⇒ anisotropic, anis dim 4 (realises u = 4).
        assert_eq!(
            anis(
                &[(one.clone(), one.clone()), (t.clone(), inv_t.clone())],
                &[]
            ),
            Some(4)
        );
        assert_eq!(
            iso(
                &[(one.clone(), one.clone()), (t.clone(), inv_t.clone())],
                &[]
            ),
            Some(false)
        );
    }

    #[test]
    fn rank_ge_5_is_isotropic() {
        let one = r2(&[1], &[1]);
        // any rank-5 form (two binary blocks + a singular vector) is isotropic.
        assert_eq!(
            iso(
                &[(one.clone(), one.clone()), (one.clone(), one.clone())],
                std::slice::from_ref(&one)
            ),
            Some(true)
        );
    }

    #[test]
    fn high_rank_anisotropic_dim_strips_explicit_hyperbolic_blocks() {
        let zero = R2::zero();
        let one = r2(&[1], &[1]);
        let t = r2(&[0, 1], &[1]);
        let inv_t = r2(&[1], &[0, 1]);
        // [0,1] is hyperbolic, leaving the rank-4 anisotropic u=4 form.
        assert_eq!(
            anis(
                &[
                    (zero.clone(), one.clone()),
                    (one.clone(), one.clone()),
                    (t.clone(), inv_t.clone())
                ],
                &[]
            ),
            Some(4)
        );
        // Two hyperbolic blocks strip away, leaving [1,1] anisotropic over F_2((t)).
        assert_eq!(
            anis(
                &[
                    (zero.clone(), one.clone()),
                    (one.clone(), zero.clone()),
                    (one.clone(), one.clone())
                ],
                &[]
            ),
            Some(2)
        );
    }

    #[test]
    fn high_rank_nonsingular_dimension_uses_aj_kernel() {
        let one = r2(&[1], &[1]);
        let t = r2(&[0, 1], &[1]);
        let inv_t = r2(&[1], &[0, 1]);

        // Three copies of [1,1] reduce in W_q to one anisotropic binary block.
        assert_eq!(
            anis(
                &[
                    (one.clone(), one.clone()),
                    (one.clone(), one.clone()),
                    (one.clone(), one.clone())
                ],
                &[]
            ),
            Some(2)
        );

        // The u=4 class plus a binary [1,1] is no longer the pure (1,1) no-wild
        // class, hence its anisotropic kernel is binary.
        assert_eq!(
            anis(
                &[
                    (one.clone(), one.clone()),
                    (t.clone(), inv_t.clone()),
                    (one.clone(), one.clone())
                ],
                &[]
            ),
            Some(2)
        );

        // A two-class quasilinear tail leaves exactly its two-dimensional radical
        // kernel even when the nonsingular part has several blocks.
        assert_eq!(
            anis(
                &[
                    (one.clone(), one.clone()),
                    (t.clone(), inv_t.clone()),
                    (one.clone(), one.clone())
                ],
                &[one.clone(), t.clone()]
            ),
            Some(2)
        );
    }

    #[test]
    fn high_rank_one_class_singular_tail_uses_clifford_invariant() {
        let one = r2(&[1], &[1]);
        let t = r2(&[0, 1], &[1]);
        let inv_t = r2(&[1], &[0, 1]);

        // Two [1,1] blocks are hyperbolic in W_q, so adding a one-class radical
        // leaves a one-dimensional anisotropic kernel.
        assert_eq!(
            anis(
                &[(one.clone(), one.clone()), (one.clone(), one.clone())],
                std::slice::from_ref(&one)
            ),
            Some(1)
        );

        // The u=4 nonsingular core has non-split odd Clifford invariant after any
        // one-class radical is added, leaving the anisotropic rank-3 pure subform.
        assert_eq!(
            anis(
                &[(one.clone(), one.clone()), (t.clone(), inv_t.clone())],
                std::slice::from_ref(&one)
            ),
            Some(3)
        );

        // Adding one more [1,1] cancels the unramified bit of the u=4 core. The
        // nonsingular kernel is binary, but it is the pure ramified class and does
        // not represent the even singular square class.
        assert_eq!(
            anis(
                &[
                    (one.clone(), one.clone()),
                    (t.clone(), inv_t.clone()),
                    (one.clone(), one.clone())
                ],
                std::slice::from_ref(&one)
            ),
            Some(3)
        );

        // Three [1,1] blocks reduce to the binary [1,1], which represents 1; the
        // semisingular kernel therefore collapses to the radical line.
        assert_eq!(
            anis(
                &[
                    (one.clone(), one.clone()),
                    (one.clone(), one.clone()),
                    (one.clone(), one.clone())
                ],
                std::slice::from_ref(&one)
            ),
            Some(1)
        );
    }

    #[test]
    fn decomposition_agrees_with_isotropy_on_rank4() {
        // A pure-binary form is hyperbolic (decomposition trivial) iff anisotropic dim 0.
        let one = r2(&[1], &[1]);
        let t = r2(&[0, 1], &[1]);
        let inv_t = r2(&[1], &[0, 1]);
        for (blocks, expect_trivial) in [
            (
                vec![(one.clone(), one.clone()), (one.clone(), one.clone())],
                true,
            ),
            (
                vec![(one.clone(), one.clone()), (t.clone(), inv_t.clone())],
                false,
            ),
        ] {
            let d = springer_decompose_local_char2(
                &Char2QuadForm::from_blocks(blocks.clone()),
                &place_t(),
            );
            let trivial = d.phi0 == 0 && d.phi1 == 0 && d.psi.is_empty();
            assert_eq!(trivial, expect_trivial);
            let dim = local_anisotropic_dim_char2(&Char2QuadForm::from_blocks(blocks), &place_t());
            assert_eq!(trivial, dim == Some(0));
        }
    }

    // ── F₄(t): the residue-trace distinguishes the F₂ case ──

    #[test]
    fn aj_oracle_over_f4() {
        type F4 = Fpn<2, 2>;
        type R4 = RationalFunction<F4>;
        let c = |n: u128| F4::from_index(n);
        let rf = |num: Vec<u128>, den: Vec<u128>| -> R4 {
            RationalFunction::new(
                num.into_iter().map(c).collect(),
                den.into_iter().map(c).collect(),
            )
        };
        let place = Char2Place::Finite(Poly::new(vec![F4::from_index(0), F4::from_index(1)])); // P = t
                                                                                               // [1, α] ⊥ [π, α/π] ↦ (1,0,1): Tr_{F₄/F₂}(α) = 1, the u = 4 anisotropic class.
        let alpha = rf(vec![2], vec![1]); // α (index 2 = the F₄ generator)
        let blocks = vec![
            (rf(vec![1], vec![1]), alpha.clone()),              // [1, α]
            (rf(vec![0, 1], vec![1]), rf(vec![2], vec![0, 1])), // [t, α/t]
        ];
        let d = springer_decompose_local_char2(&Char2QuadForm::from_blocks(blocks.clone()), &place);
        assert_eq!(d.phi0, 1);
        assert_eq!(d.phi1, 1);
        assert!(d.psi.is_empty());
        assert_eq!(
            local_anisotropic_dim_char2(&Char2QuadForm::from_blocks(blocks), &place),
            Some(4)
        );
        // [1,1] over F₄ is hyperbolic (Tr_{F₄/F₂}(1) = 0).
        let h = springer_decompose_local_char2(
            &Char2QuadForm::from_blocks(vec![(rf(vec![1], vec![1]), rf(vec![1], vec![1]))]),
            &place,
        );
        assert_eq!((h.phi0, h.phi1, h.psi.len()), (0, 0, 0));
    }

    // ── a degree-2 place (κ = F₄): exercise the P-adic digit machinery ──

    #[test]
    fn decomposition_at_a_degree_two_place() {
        // P = t²+t+1 (irreducible over F₂, κ = F₄, π = P). [1, 1/P] has a simple wild
        // pole: the P-adic digit at P⁻¹ is the κ-element 1, so ψ = {1: 1}, and the
        // block is anisotropic (rank-2, ab = 1/P ∉ ℘(K_P)).
        let p = Char2Place::Finite(p2(&[1, 1, 1])); // t²+t+1
        let blocks = vec![(r2(&[1], &[1]), r2(&[1], &[1, 1, 1]))]; // [1, 1/(t²+t+1)]
        let d = springer_decompose_local_char2(&Char2QuadForm::from_blocks(blocks.clone()), &p);
        assert_eq!(d.phi0, 0);
        assert_eq!(d.phi1, 0);
        assert_eq!(d.psi, BTreeMap::from([(1usize, Poly::<F2>::one())]));
        assert_eq!(
            local_anisotropic_dim_char2(&Char2QuadForm::from_blocks(blocks), &p),
            Some(2)
        );
    }

    #[test]
    fn degree_two_place_keeps_hensel_carries_in_asnf() {
        // P = t²+t+1. In the completion at P,
        // ℘(t/P) = (t/P)² + t/P = (t³+t)/P², so [1, ℘(t/P)] is hyperbolic.
        // Treating polynomial P-adic digits as κ[[P]] coefficients drops the
        // Hensel carries and leaves a false wild obstruction here.
        let p_poly = p2(&[1, 1, 1]);
        let place = Char2Place::Finite(p_poly.clone());
        let p_sq = p_poly.mul(&p_poly);
        let wp = R2::new(p2(&[0, 1, 0, 1]).coeffs().to_vec(), p_sq.coeffs().to_vec());

        assert!(asnf::local_is_pe(&wp, &place));
        let form = Char2QuadForm::from_blocks(vec![(r2(&[1], &[1]), wp)]);
        let d = springer_decompose_local_char2(&form, &place);
        assert_eq!((d.phi0, d.phi1), (0, 0));
        assert!(d.psi.is_empty());
        assert_eq!(local_anisotropic_dim_char2(&form, &place), Some(0));
        assert_eq!(local_is_isotropic_char2(&form, &place), Some(true));
    }

    // ── global isotropy over F_q(t) (Hasse–Minkowski), source-pinned oracles ──
    // (verified independently, then cross-checked via Codex; oracle 7 below is the
    // corrected one — [1,1]⊥[t,t] is ISOTROPIC, vector (1,0,1,1).)

    fn gi(blocks: &[(R2, R2)], singular: &[R2]) -> Option<bool> {
        is_isotropic_global_char2(&form(blocks, singular))
    }

    #[test]
    fn global_pe_direct() {
        // global_is_pe over F₂(t): ℘(t)=t²+t ∈ ℘; ℘(1/t)=(1+t)/t² ∈ ℘; 0 ∈ ℘.
        assert!(global_is_pe(&r2(&[0, 1, 1], &[1]))); // t²+t = ℘(t)
        assert!(global_is_pe(&r2(&[1, 1], &[0, 0, 1]))); // (1+t)/t² = ℘(1/t)
        assert!(global_is_pe(&R2::zero())); // 0
                                            // not in ℘: t (odd pole at ∞), 1/t (odd pole at 0), 1 (constant Tr=1).
        assert!(!global_is_pe(&r2(&[0, 1], &[1]))); // t
        assert!(!global_is_pe(&r2(&[1], &[0, 1]))); // 1/t
        assert!(!global_is_pe(&r2(&[1], &[1]))); // 1
    }

    #[test]
    fn global_rank2_pe_obstruction() {
        let one = r2(&[1], &[1]);
        // [1, t²+t]: ab = ℘(t) ∈ ℘ ⇒ isotropic.
        assert_eq!(gi(&[(one.clone(), r2(&[0, 1, 1], &[1]))], &[]), Some(true));
        // [1, 1]: ab = 1, Tr_{F₂/F₂}(1) = 1 ∉ ℘ ⇒ anisotropic (the F₄ norm form). NOT a
        // bad-place sweep: 1 is a unit everywhere; the obstruction is the constant trace.
        assert_eq!(gi(&[(one.clone(), one.clone())], &[]), Some(false));
        // [1, t]: ab = t, odd pole at ∞ ⇒ anisotropic.
        assert_eq!(gi(&[(one.clone(), r2(&[0, 1], &[1]))], &[]), Some(false));
    }

    #[test]
    fn global_rank3() {
        let one = r2(&[1], &[1]);
        let t = r2(&[0, 1], &[1]);
        // [1, t] ⊥ ⟨1⟩: isotropic — (1,0,1) gives 1+1=0 (s_v(t,1)=0 ∀v).
        assert_eq!(
            gi(&[(one.clone(), t.clone())], std::slice::from_ref(&one)),
            Some(true)
        );
        // [1, 1] ⊥ ⟨t⟩: anisotropic — s_v(1, t) = 1 at ∞ (t not a norm from F₄(t)).
        assert_eq!(
            gi(&[(one.clone(), one.clone())], std::slice::from_ref(&t)),
            Some(false)
        );
    }

    #[test]
    fn global_rank4() {
        let one = r2(&[1], &[1]);
        let t = r2(&[0, 1], &[1]);
        let inv_t = r2(&[1], &[0, 1]);
        // [1, t] ⊥ [1, t]: isotropic — equal summands, q(v, v) = 0.
        assert_eq!(
            gi(&[(one.clone(), t.clone()), (one.clone(), t.clone())], &[]),
            Some(true)
        );
        // [1,1] ⊥ [t, t]: ISOTROPIC — (1,0,1,1): [1,1](1,0)+[t,t](1,1)=1+1=0. (The norm
        // form of the division algebra [1,t) is [1,1] ⊥ t·[1,1], and t·[1,1] ≇ [t,t].)
        assert_eq!(
            gi(&[(one.clone(), one.clone()), (t.clone(), t.clone())], &[]),
            Some(true)
        );
        // [1,1] ⊥ [t, 1/t]: anisotropic — the u = 4 class, anisotropic already at t=0.
        assert_eq!(
            gi(
                &[(one.clone(), one.clone()), (t.clone(), inv_t.clone())],
                &[]
            ),
            Some(false)
        );
    }

    #[test]
    fn global_rank_ge_5_is_isotropic() {
        let one = r2(&[1], &[1]);
        let t = r2(&[0, 1], &[1]);
        let inv_t = r2(&[1], &[0, 1]);
        // [1,1] ⊥ [t,1/t] ⊥ ⟨1⟩: the dim-4 core is anisotropic, but u(F₂(t)) = 4 ⇒ the
        // dim-5 form is isotropic (C₂, Tsen–Lang).
        assert_eq!(
            gi(
                &[(one.clone(), one.clone()), (t.clone(), inv_t.clone())],
                std::slice::from_ref(&one)
            ),
            Some(true)
        );
    }

    #[test]
    fn global_totally_singular() {
        let one = r2(&[1], &[1]);
        let t = r2(&[0, 1], &[1]);
        let t2 = r2(&[0, 0, 1], &[1]);
        // ⟨1, t⟩: t ∉ K² (odd degree) ⇒ anisotropic.
        assert_eq!(gi(&[], &[one.clone(), t.clone()]), Some(false));
        // ⟨1, t²⟩: t² ∈ K² ⇒ K²-dependent ⇒ isotropic.
        assert_eq!(gi(&[], &[one.clone(), t2.clone()]), Some(true));
        // ⟨1, t, t²⟩: rank 3 ≥ 3 ⇒ K²-dependent ⇒ isotropic.
        assert_eq!(gi(&[], &[one.clone(), t.clone(), t2.clone()]), Some(true));
        // [1,1] ⊥ ⟨1, t⟩: the anisotropic ⟨1,t⟩ is universal ⇒ it isotropises the block.
        assert_eq!(
            gi(&[(one.clone(), one.clone())], &[one.clone(), t.clone()]),
            Some(true)
        );
    }

    #[test]
    fn global_over_f4() {
        type F4 = Fpn<2, 2>;
        type R4 = RationalFunction<F4>;
        let c = |n: u128| F4::from_index(n);
        let rf = |num: Vec<u128>, den: Vec<u128>| -> R4 {
            RationalFunction::new(
                num.into_iter().map(c).collect(),
                den.into_iter().map(c).collect(),
            )
        };
        let one = rf(vec![1], vec![1]);
        let alpha = rf(vec![2], vec![1]); // α (index 2 = the F₄ generator)
                                          // [1, 1] over F₄(t): ab = 1 = α²+α = ℘(α) ⇒ isotropic (Tr_{F₄/F₂}(1) = 0).
        assert_eq!(
            is_isotropic_global_char2(&Char2QuadForm::from_blocks(vec![(
                one.clone(),
                one.clone()
            )])),
            Some(true)
        );
        // [1, α] over F₄(t): ab = α, Tr_{F₄/F₂}(α) = 1 ⇒ anisotropic.
        assert_eq!(
            is_isotropic_global_char2(&Char2QuadForm::from_blocks(vec![(
                one.clone(),
                alpha.clone()
            )])),
            Some(false)
        );
    }
}
