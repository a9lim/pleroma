//! The Minkowski–Siegel mass formula for even unimodular lattices, and the Leech
//! lattice `Λ₂₄` — the payload that closes the loop to the repo's namesake.
//!
//! The **mass** of a genus is `∑_{L} 1/|Aut(L)|` over its isometry classes; the
//! Minkowski–Siegel formula computes it from local data. For the genus of
//! positive-definite **even unimodular** lattices of rank `n ≡ 0 (mod 8)` it has a
//! clean closed form in Bernoulli numbers ([`mass_even_unimodular`]):
//!
//! ```text
//! mass(n) = (|B_{n/2}| / n) · ∏_{j=1}^{n/2−1} |B_{2j}| / (4j).
//! ```
//!
//! At `n = 8` this is `1/696729600`. Since `E₈` is the *unique* even unimodular
//! rank-8 lattice (alone in its genus), the mass there equals `1/|Aut(E₈)|`, so the
//! formula **recovers** `|Aut(E₈)| = |W(E₈)| = 696729600` — a number the
//! brute-force automorphism counter ([`super::lattice::IntegralForm::automorphism_group_order`])
//! deliberately refuses (it is past the node budget). The mass climbs fast
//! (`n = 16`: `691/277667181515243520000`; `n = 24`: the Niemeier mass), and the
//! `i128` rationals here reach to `n = 24` but no further — the honest ceiling.
//!
//! The **Leech lattice** `Λ₂₄` ([`leech`]) is built from the extended binary Golay
//! code `[24, 12, 8]`: a spanning set of `√8·Λ₂₄ ⊂ ℤ²⁴` (twice the Golay generator
//! rows, the `4(e₀+eᵢ)` glue, and one odd `(−3, 1²³)` vector), reduced to a basis
//! `B` by Hermite normal form, with `Gram = B·Bᵀ / 8`. The construction is
//! **validated, not trusted**: a rank-24 even unimodular lattice with **no roots**
//! (no norm-2 vectors) is the Leech lattice by Niemeier's classification, and
//! [`leech`]'s tests check exactly that (`det = 1`, even, `short_vectors(2)` empty).
//! Its automorphism group is the Conway group `Co₀ = 2·Co₁` of order
//! [`LEECH_AUT_ORDER`] — far past any brute-force reach, recorded as a constant.
//!
//! The **Monster** is *not* here on purpose: `Λ₂₄ → Co₀ → Co₁` is a quadratic-form
//! computation, but the Monster is monstrous moonshine (vertex operator algebras,
//! the Griess algebra, the `j`-function) — a different field, not a lattice
//! invariant. It stays a thematic remark, not a build target.
//!
//! References: Conway–Sloane *SPLAG* Ch. 16 (the mass formula) and Ch. 4/10 (Leech
//! and the Golay code); R. Wilson, *The finite simple groups*, Thm 1.4 (the Leech
//! coordinate conditions). Mass values and the Golay generator cross-checked with
//! an independent Codex pass.

use crate::forms::lattice::IntegralForm;
use crate::linalg::integer::normalize_relation_rows;

/// `|Aut(Λ₂₄)| = |Co₀| = 2·|Co₁| = 2²² · 3⁹ · 5⁴ · 7² · 11 · 13 · 23`. Recorded, not
/// computed: it is far past the brute-force automorphism budget.
pub const LEECH_AUT_ORDER: u128 = 8_315_553_613_086_720_000;

/// `|B_k|` for an even `k`, as `(|numerator|, denominator)` in lowest terms — the
/// Bernoulli numbers up to `B₂₄` (all the mass formula needs through rank 24).
fn bernoulli_abs(k: u32) -> (i128, i128) {
    match k {
        2 => (1, 6),
        4 => (1, 30),
        6 => (1, 42),
        8 => (1, 30),
        10 => (5, 66),
        12 => (691, 2730),
        14 => (7, 6),
        16 => (3617, 510),
        18 => (43867, 798),
        20 => (174611, 330),
        22 => (854513, 138),
        24 => (236364091, 2730),
        _ => panic!("bernoulli_abs only tabulated for even k in 2..=24"),
    }
}

fn gcd(a: i128, b: i128) -> i128 {
    let (mut a, mut b) = (a.abs(), b.abs());
    while b != 0 {
        let t = b;
        b = a % b;
        a = t;
    }
    a
}

/// Multiply two non-negative reduced fractions, cross-reducing first to keep the
/// intermediates small. Returns `None` if the (already reduced) product still
/// overflows `i128`.
fn checked_rat_mul(
    (mut a, mut b): (i128, i128),
    (mut c, mut d): (i128, i128),
) -> Option<(i128, i128)> {
    let g1 = gcd(a, d);
    a /= g1;
    d /= g1;
    let g2 = gcd(c, b);
    c /= g2;
    b /= g2;
    let num = a.checked_mul(c)?;
    let den = b.checked_mul(d)?;
    let g = gcd(num, den);
    Some((num / g, den / g))
}

/// The mass of the genus of positive-definite **even unimodular** lattices of rank
/// `n`. `n` must be a positive multiple of 8 (otherwise no such lattice exists).
/// Returns the exact reduced fraction `(numerator, denominator)`, or `None` if the
/// value overflows the `i128` rational model (beyond `n = 24`).
pub fn mass_even_unimodular(n: u32) -> Option<(i128, i128)> {
    if n == 0 || !n.is_multiple_of(8) {
        return None;
    }
    let m = n / 2;
    let mut acc = (1i128, 1i128);
    for j in 1..m {
        acc = checked_rat_mul(acc, bernoulli_abs(2 * j))?;
        acc = checked_rat_mul(acc, (1, 4 * j as i128))?;
    }
    // times |B_{n/2}| / n
    acc = checked_rat_mul(acc, bernoulli_abs(m))?;
    acc = checked_rat_mul(acc, (1, n as i128))?;
    Some(acc)
}

/// The extended binary Golay code generator `G = [I₁₂ | A]` (12×24), as a `0/1`
/// integer matrix. `A` is the standard icosahedron-based block.
fn golay_generator() -> Vec<Vec<i128>> {
    let a: [[i128; 12]; 12] = [
        [1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1],
        [0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1],
        [0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1],
        [0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1],
        [0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1],
        [0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1],
        [1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0],
        [1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0],
        [1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0],
        [1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0],
        [1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0],
        [1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1],
    ];
    (0..12)
        .map(|i| {
            let mut row = vec![0i128; 24];
            row[i] = 1; // I_12 block
            row[12..24].copy_from_slice(&a[i]); // A block
            row
        })
        .collect()
}

/// A spanning set of `√8 · Λ₂₄ ⊂ ℤ²⁴` (Wilson's coordinate model): twice the Golay
/// generator rows, the `4(e₀ + eᵢ)` glue vectors, and the odd `(−3, 1²³)` vector.
/// Every spanning vector has norm 32 (`= 4` after the `1/√8` scaling).
fn leech_spanning_set() -> Vec<Vec<i128>> {
    let mut s: Vec<Vec<i128>> = Vec::new();
    for g in golay_generator() {
        s.push(g.iter().map(|&x| 2 * x).collect());
    }
    for i in 1..24 {
        let mut v = vec![0i128; 24];
        v[0] = 4;
        v[i] = 4;
        s.push(v);
    }
    let mut odd = vec![1i128; 24];
    odd[0] = -3;
    s.push(odd);
    s
}

/// The Leech lattice `Λ₂₄`: the unique even unimodular rank-24 lattice with no
/// roots, scaled to minimum 4. Built from the Golay spanning set reduced to a basis
/// `B` by Hermite normal form, with `Gram = B·Bᵀ / 8`.
pub fn leech() -> IntegralForm {
    let basis = normalize_relation_rows(leech_spanning_set());
    assert_eq!(basis.len(), 24, "Leech spanning set must have rank 24");
    let n = 24;
    let mut gram = vec![vec![0i128; n]; n];
    for i in 0..n {
        for j in 0..n {
            let dot: i128 = basis[i].iter().zip(&basis[j]).map(|(&x, &y)| x * y).sum();
            assert!(
                dot % 8 == 0,
                "√8·Λ inner products must be divisible by 8 (got {dot})"
            );
            gram[i][j] = dot / 8;
        }
    }
    IntegralForm::new(gram).expect("Leech Gram is symmetric")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forms::root_lattices::e_8;

    #[test]
    fn mass_recovers_e8_weyl_group_order() {
        // mass(8) = 1/696729600 = 1/|W(E_8)| — the formula hands back the
        // automorphism count the brute-force counter refuses to compute.
        assert_eq!(mass_even_unimodular(8), Some((1, 696_729_600)));
    }

    #[test]
    fn mass_matches_known_values_through_24() {
        // n = 16: the two even unimodular lattices E_8⊕E_8 and D_16+.
        assert_eq!(
            mass_even_unimodular(16),
            Some((691, 277_667_181_515_243_520_000))
        );
        // n = 24: the Niemeier mass (24 classes). Still inside i128.
        assert_eq!(
            mass_even_unimodular(24),
            Some((
                1_027_637_932_586_061_520_960_267,
                129_477_933_340_026_851_560_636_148_613_120_000_000
            ))
        );
    }

    #[test]
    fn mass_rejects_non_multiples_of_eight() {
        assert_eq!(mass_even_unimodular(0), None);
        assert_eq!(mass_even_unimodular(4), None);
        assert_eq!(mass_even_unimodular(12), None);
        assert_eq!(mass_even_unimodular(7), None);
    }

    #[test]
    fn golay_is_a_self_dual_doubly_even_code() {
        // The extended Golay code is self-dual: every pair of generator rows has
        // even inner product over F_2, and every row has weight ≡ 0 mod 4
        // (doubly even). Weight 12 for the I-row + complement-ish rows here.
        let g = golay_generator();
        for row in &g {
            let wt: i128 = row.iter().sum();
            assert_eq!(wt % 4, 0, "doubly even: weight {wt}");
        }
        for i in 0..12 {
            for j in 0..12 {
                let ip: i128 = g[i].iter().zip(&g[j]).map(|(&a, &b)| a * b).sum();
                assert_eq!(ip % 2, 0, "self-orthogonal rows {i},{j}");
            }
        }
    }

    #[test]
    fn leech_is_the_rootless_even_unimodular_rank24() {
        let l = leech();
        assert_eq!(l.dim(), 24);
        assert_eq!(l.determinant(), 1, "unimodular");
        assert!(l.is_even(), "even");
        // No roots (no norm-2 vectors): even + unimodular + rank 24 + rootless ⇒
        // Leech, by Niemeier's classification. (We check norm ≤ 2 only — the full
        // minimum 4 would enumerate all 196560 minimal vectors, far too expensive.)
        let roots = l.short_vectors(2).expect("positive definite");
        assert!(
            roots.is_empty(),
            "Leech has no roots (found {})",
            roots.len()
        );
    }

    #[test]
    fn leech_aut_order_has_the_co0_factorisation() {
        // |Co_0| = 2^22 · 3^9 · 5^4 · 7^2 · 11 · 13 · 23.
        let expected: u128 =
            2u128.pow(22) * 3u128.pow(9) * 5u128.pow(4) * 7u128.pow(2) * 11 * 13 * 23;
        assert_eq!(LEECH_AUT_ORDER, expected);
        // E_8 anchors the rank-8 mass; Leech the (uncomputable-here) rank-24 one.
        assert!(e_8().is_even() && e_8().is_unimodular());
    }
}
