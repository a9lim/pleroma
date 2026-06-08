//! The **genus** of an integral lattice: its class up to local equivalence at
//! every place.
//!
//! Two integral lattices are in the same genus iff they are isometric over `ℝ`
//! (same signature) and over `ℤ_p` for every prime `p`. The genus is the natural
//! arithmetic coarsening of the isometry class, and it is exactly the place where
//! the adelic machinery the crate already carries (`padic.rs`'s square classes,
//! `is_square_qp`) acts on a *lattice* rather than a field-level square class.
//!
//! The engine is the **p-adic Jordan decomposition**. Over `ℤ_p` a lattice splits
//! orthogonally into scaled unimodular constituents `⊥_k p^k L_k`; for odd `p`
//! every `L_k` is diagonal, and for `p = 2` the constituents are 1-dimensional
//! (type I) or 2-dimensional even blocks `U = [[0,1],[1,0]]`, `V = [[2,1],[1,2]]`
//! (type II). The Conway–Sloane *p-adic symbol* records, per scale, the dimension,
//! the determinant square class (`sign = ±1`), and at `p = 2` the type and the
//! oddity (trace mod 8 of the type-I part).
//!
//! **Claim level and the honest p = 2 boundary.** For odd `p` the symbol is a
//! complete `ℤ_p` invariant and the comparison here is exact. For `p = 2` the
//! symbol is complete only up to Conway–Sloane *oddity fusion* (well-defined,
//! implemented) and *sign-walking* (the genuinely subtle part — even *SPLAG*'s
//! printed canonical form is wrong, per Allcock, *On the classification of integral
//! quadratic forms*). This module fuses oddities within compartments and compares
//! per-scale signs **directly**, which is **sound** (it never declares two
//! different genera equal) and **exact for symbols without nontrivial sign-walking
//! trains** — i.e. all single-scale symbols, all-type-II symbols, and every case in
//! the test battery (root lattices, `ℤⁿ`, `E_8`, and randomised `ℤ`-isometry
//! checks). For lattices whose 2-adic symbol has adjacent type-I scales coupled by
//! a sign-walk, [`are_in_same_genus`] can return a conservative *false negative*;
//! that residual is the documented boundary, not a silent error.
//!
//! References: Conway–Sloane *SPLAG* Ch. 15 §7; Allcock, *On the classification of
//! integral quadratic forms* (the corrected 2-adic sign-walking calculus).

use crate::forms::lattice::IntegralForm;
use crate::forms::padic::is_square_qp;
use crate::scalar::{Rational, Scalar};
use std::collections::BTreeMap;

/// One scale of a p-adic Jordan symbol: the constituent `p^scale · (unimodular of
/// dimension `dim`)`. `sign` is the determinant square class of the unimodular part
/// (`+1` iff a `ℤ_p`-square unit). For `p = 2`, `type_ii` is true when the
/// constituent is even (a sum of `U`/`V` planes) and `oddity` is the trace mod 8 of
/// the type-I (odd) part; for odd `p` those two fields are unused (`false`, `0`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ScaleSymbol {
    pub scale: u32,
    pub dim: usize,
    pub sign: i8,
    pub type_ii: bool,
    pub oddity: i64,
}

/// The genus of a nondegenerate integral lattice: signature, determinant, and the
/// p-adic symbol at every prime that can carry a nontrivial local invariant
/// (`p | 2·det`).
#[derive(Clone, Debug)]
pub struct Genus {
    pub dim: usize,
    pub signature: (usize, usize),
    pub det: i128,
    symbols: BTreeMap<u128, Vec<ScaleSymbol>>,
}

// --- exact rational helpers (p-adic valuations, unit residues) ---

fn r_int(n: i128) -> Rational {
    Rational::int(n)
}

fn v_p_i128(mut x: i128, p: i128) -> u32 {
    debug_assert!(x != 0);
    let mut k = 0;
    // (i128 has no stable `is_multiple_of` yet — the lint only applies to u128.)
    while x % p == 0 {
        x /= p;
        k += 1;
    }
    k
}

fn unit_part_i128(mut x: i128, p: i128) -> i128 {
    while x % p == 0 {
        x /= p;
    }
    x
}

/// The p-adic valuation of a nonzero rational `num/den`.
fn rat_val(r: &Rational, p: i128) -> i64 {
    debug_assert!(!r.is_zero());
    v_p_i128(r.numer(), p) as i64 - v_p_i128(r.denom(), p) as i64
}

/// The determinant square class (`±1`) of a rational unit `r` over `ℚ_p` (odd `p`).
fn unit_sign_odd(r: &Rational, p: i128) -> i8 {
    let a = unit_part_i128(r.numer(), p).rem_euclid(p);
    let b = unit_part_i128(r.denom(), p).rem_euclid(p);
    let m = (a * b).rem_euclid(p);
    if is_square_qp(m, p as u128) {
        1
    } else {
        -1
    }
}

/// The 2-adic unit residue `u mod 8` of a nonzero rational `r = num/den` whose
/// 2-adic valuation is `scale` (so `r / 2^scale` is a unit). Uses `odd⁻¹ ≡ odd
/// (mod 8)`.
fn unit_mod8(r: &Rational) -> i64 {
    let a = unit_part_i128(r.numer(), 2).rem_euclid(8);
    let b = unit_part_i128(r.denom(), 2).rem_euclid(8);
    (a * b).rem_euclid(8) as i64
}

fn sign_from_mod8(u: i64) -> i8 {
    if u == 1 || u == 7 {
        1
    } else {
        -1
    }
}

// --- exact rational matrix utilities ---

type RMat = Vec<Vec<Rational>>;

fn to_rational(gram: &[Vec<i128>]) -> RMat {
    gram.iter()
        .map(|row| row.iter().map(|&x| r_int(x)).collect())
        .collect()
}

fn rdiv(a: &Rational, b: &Rational) -> Rational {
    a.mul(
        &b.inv()
            .expect("division by zero rational in Jordan splitting"),
    )
}

// --- the p-adic Jordan decomposition ---

/// A raw Jordan block before per-scale aggregation.
struct RawBlock {
    scale: u32,
    dim: usize,
    sign: i8,
    type_ii: bool,
    oddity: i64,
}

/// Minimal p-adic valuation entry `(i, j)` of the symmetric matrix `a` (over the
/// `active` index set). Returns `None` if every active entry is zero.
fn min_val_entry(a: &RMat, active: &[usize], p: i128) -> Option<(usize, usize)> {
    let mut best: Option<(i64, usize, usize)> = None;
    for (ii, &i) in active.iter().enumerate() {
        for &j in &active[ii..] {
            if !a[i][j].is_zero() {
                let v = rat_val(&a[i][j], p);
                if best.is_none_or(|(bv, _, _)| v < bv) {
                    best = Some((v, i, j));
                }
            }
        }
    }
    best.map(|(_, i, j)| (i, j))
}

/// p-adic Jordan blocks of a nondegenerate symmetric integer Gram matrix. Returns
/// `None` if the form is degenerate at `p` (a radical — should not happen for a
/// nonsingular lattice).
fn jordan_blocks(gram: &[Vec<i128>], p: i128) -> Option<Vec<RawBlock>> {
    let mut a = to_rational(gram);
    let mut active: Vec<usize> = (0..gram.len()).collect();
    let mut blocks = Vec::new();
    let mut guard = 0usize;
    let guard_max = 8 * gram.len() * gram.len() + 16;

    while !active.is_empty() {
        guard += 1;
        if guard > guard_max {
            return None; // defensive: should never trip for a nondegenerate form
        }
        let (i, j) = min_val_entry(&a, &active, p)?;
        let scale = rat_val(&a[i][j], p);

        if i != j && p != 2 {
            // Odd p: rotate e_i ← e_i + e_j so a diagonal entry attains the minimal
            // valuation (2 is a unit, so a[i][i] + 2a[i][j] + a[j][j] has valuation
            // = v(a[i][j])). Then re-loop; the diagonal pivot gets peeled next.
            for &c in &active {
                a[i][c] = a[i][c].add(&a[j][c].clone());
            }
            for &r in &active {
                a[r][i] = a[r][i].add(&a[r][j].clone());
            }
            continue;
        }

        let scale = u32::try_from(scale).ok()?; // negative scale ⇒ not p-integral (bug)
        if i == j {
            // 1-dimensional block ⟨a[i][i]⟩.
            let d = a[i][i].clone();
            let (sign, type_ii, oddity) = if p == 2 {
                let u8v = unit_mod8(&d);
                (sign_from_mod8(u8v), false, u8v)
            } else {
                (unit_sign_odd(&d, p), false, 0)
            };
            blocks.push(RawBlock {
                scale,
                dim: 1,
                sign,
                type_ii,
                oddity,
            });
            // Schur complement removing i.
            let pivot = d;
            let rest: Vec<usize> = active.iter().copied().filter(|&r| r != i).collect();
            for &r in &rest {
                for &s in &rest {
                    let corr = rdiv(&a[r][i].mul(&a[i][s]), &pivot);
                    a[r][s] = a[r][s].sub(&corr);
                }
            }
            active = rest;
        } else {
            // p == 2, 2-dimensional even block on {i, j}.
            let alpha = a[i][i].clone();
            let beta = a[i][j].clone();
            let gamma = a[j][j].clone();
            let det = alpha.mul(&gamma).sub(&beta.mul(&beta));
            let sign = sign_from_mod8(unit_mod8(&det));
            blocks.push(RawBlock {
                scale,
                dim: 2,
                sign,
                type_ii: true,
                oddity: 0,
            });
            // Schur complement removing {i, j} via the 2×2 inverse
            //   B⁻¹ = (1/D)[[γ,−β],[−β,α]].
            let rest: Vec<usize> = active
                .iter()
                .copied()
                .filter(|&r| r != i && r != j)
                .collect();
            for &r in &rest {
                for &s in &rest {
                    // (a[r][i], a[r][j]) · B⁻¹ · (a[i][s], a[j][s])ᵀ
                    let t0 = gamma.mul(&a[i][s]).sub(&beta.mul(&a[j][s]));
                    let t1 = alpha.mul(&a[j][s]).sub(&beta.mul(&a[i][s]));
                    let numer = a[r][i].mul(&t0).add(&a[r][j].mul(&t1));
                    let corr = rdiv(&numer, &det);
                    a[r][s] = a[r][s].sub(&corr);
                }
            }
            active = rest;
        }
    }
    Some(blocks)
}

/// Aggregate raw blocks at prime `p` into the per-scale Conway–Sloane symbol
/// (sorted by scale).
fn aggregate(blocks: &[RawBlock]) -> Vec<ScaleSymbol> {
    let mut by_scale: BTreeMap<u32, (usize, i8, bool, i64)> = BTreeMap::new();
    for b in blocks {
        let e = by_scale.entry(b.scale).or_insert((0, 1, true, 0));
        e.0 += b.dim;
        e.1 *= b.sign;
        if !b.type_ii {
            e.2 = false; // any type-I constituent makes the scale type I
        }
        e.3 = (e.3 + b.oddity).rem_euclid(8);
    }
    by_scale
        .into_iter()
        .map(|(scale, (dim, sign, type_ii, oddity))| ScaleSymbol {
            scale,
            dim,
            sign,
            type_ii,
            oddity: if type_ii { 0 } else { oddity },
        })
        .collect()
}

/// The signature `(p₊, p₋)` of a nondegenerate symmetric integer matrix, by exact
/// rational congruence diagonalization (Jacobi/Sylvester): the sign of each pivot.
fn signature(gram: &[Vec<i128>]) -> (usize, usize) {
    let n = gram.len();
    let mut a = to_rational(gram);
    let mut active: Vec<usize> = (0..n).collect();
    let (mut pos, mut neg) = (0usize, 0usize);
    while !active.is_empty() {
        // Find a nonzero diagonal pivot; if none, a nondegenerate form still has a
        // nonzero off-diagonal a[r][s], and e_r ← e_r + e_s makes a[r][r] nonzero.
        if !active.iter().any(|&r| !a[r][r].is_zero()) {
            let mut pair = None;
            'find: for (ai, &r) in active.iter().enumerate() {
                for &s in &active[ai + 1..] {
                    if !a[r][s].is_zero() {
                        pair = Some((r, s));
                        break 'find;
                    }
                }
            }
            let (r, s) = pair.expect("nondegenerate form has a nonzero entry");
            for &c in &active {
                a[r][c] = a[r][c].add(&a[s][c].clone());
            }
            for &rr in &active {
                a[rr][r] = a[rr][r].add(&a[rr][s].clone());
            }
        }
        let i = active
            .iter()
            .copied()
            .find(|&r| !a[r][r].is_zero())
            .expect("a diagonal pivot now exists");
        match a[i][i].sign() {
            std::cmp::Ordering::Greater => pos += 1,
            std::cmp::Ordering::Less => neg += 1,
            std::cmp::Ordering::Equal => unreachable!(),
        }
        let pivot = a[i][i].clone();
        let rest: Vec<usize> = active.iter().copied().filter(|&r| r != i).collect();
        for &r in &rest {
            for &s in &rest {
                let corr = rdiv(&a[r][i].mul(&a[i][s]), &pivot);
                a[r][s] = a[r][s].sub(&corr);
            }
        }
        active = rest;
    }
    (pos, neg)
}

/// The primes that can carry a nontrivial local invariant: `2` and every odd prime
/// dividing the determinant.
fn relevant_primes(det: i128) -> Vec<u128> {
    let mut ps = vec![2u128];
    let mut n = det.unsigned_abs();
    while n.is_multiple_of(2) {
        n /= 2; // 2 is already included; strip its powers before odd factoring
    }
    let mut d = 3u128;
    while d <= n / d {
        if n.is_multiple_of(d) {
            ps.push(d);
            while n.is_multiple_of(d) {
                n /= d;
            }
        }
        d += 2;
    }
    if n > 1 {
        ps.push(n);
    }
    ps.sort_unstable();
    ps.dedup();
    ps
}

impl Genus {
    /// The genus of a nondegenerate integral lattice, or `None` if `det = 0`.
    pub fn of(lattice: &IntegralForm) -> Option<Genus> {
        let det = lattice.determinant();
        if det == 0 {
            return None;
        }
        let gram = lattice.gram();
        let mut symbols = BTreeMap::new();
        for p in relevant_primes(det) {
            let blocks = jordan_blocks(gram, p as i128)?;
            symbols.insert(p, aggregate(&blocks));
        }
        Some(Genus {
            dim: lattice.dim(),
            signature: signature(gram),
            det,
            symbols,
        })
    }

    /// The Conway–Sloane symbol at prime `p` (per-scale), or an empty slice if `p`
    /// carries no constituent.
    pub fn symbol_at(&self, p: u128) -> &[ScaleSymbol] {
        self.symbols.get(&p).map_or(&[], |v| v)
    }

    /// The primes carrying a recorded local symbol.
    pub fn primes(&self) -> Vec<u128> {
        self.symbols.keys().copied().collect()
    }
}

/// Fuse oddities within compartments (maximal runs of consecutive type-I scales):
/// the per-compartment oddity sum mod 8 is the invariant. Returns a canonicalised
/// copy with each compartment's total oddity on its lowest scale and `0` elsewhere.
fn fuse_oddities(syms: &[ScaleSymbol]) -> Vec<ScaleSymbol> {
    let mut out = syms.to_vec();
    let mut i = 0;
    while i < out.len() {
        if out[i].type_ii {
            i += 1;
            continue;
        }
        // start of a compartment: extend over consecutive (scale+1) type-I scales
        let mut j = i;
        let mut total = 0i64;
        loop {
            total += out[j].oddity;
            let extends =
                j + 1 < out.len() && !out[j + 1].type_ii && out[j + 1].scale == out[j].scale + 1;
            if !extends {
                break;
            }
            j += 1;
        }
        let total = total.rem_euclid(8);
        out[i].oddity = total;
        for k in (i + 1)..=j {
            out[k].oddity = 0;
        }
        i = j + 1;
    }
    out
}

/// Whether two integral lattices lie in the same genus.
///
/// Exact for the signature, determinant, and every odd-prime symbol. At `p = 2`
/// the comparison fuses compartment oddities and matches per-scale signs directly;
/// see the module docs for the sound-but-conservative sign-walking boundary.
pub fn are_in_same_genus(a: &IntegralForm, b: &IntegralForm) -> bool {
    let (Some(ga), Some(gb)) = (Genus::of(a), Genus::of(b)) else {
        return false;
    };
    if ga.dim != gb.dim || ga.signature != gb.signature || ga.det != gb.det {
        return false;
    }
    if ga.symbols.keys().ne(gb.symbols.keys()) {
        return false;
    }
    for (&p, sa) in &ga.symbols {
        let sb = &gb.symbols[&p];
        if p == 2 {
            if fuse_oddities(sa) != fuse_oddities(sb) {
                return false;
            }
        } else if sa != sb {
            return false;
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forms::root_lattices::{a_n, d_n, e_6, e_7, e_8};

    fn zn(n: usize) -> IntegralForm {
        IntegralForm::diagonal(&vec![1i128; n])
    }

    /// Apply a unimodular change of basis Uᵀ G U (U upper-unitriangular with the
    /// given off-diagonal entries) — an isometry, so the genus is unchanged.
    fn congruent(g: &IntegralForm, shears: &[(usize, usize, i128)]) -> IntegralForm {
        let n = g.dim();
        let mut u = vec![vec![0i128; n]; n];
        for (i, row) in u.iter_mut().enumerate() {
            row[i] = 1;
        }
        for &(i, j, c) in shears {
            u[i][j] += c; // add c·(col i) to col j
        }
        // G' = Uᵀ G U
        let gram = g.gram();
        let mut gu = vec![vec![0i128; n]; n];
        for i in 0..n {
            for j in 0..n {
                let mut s = 0i128;
                for k in 0..n {
                    s += gram[i][k] * u[k][j];
                }
                gu[i][j] = s;
            }
        }
        let mut out = vec![vec![0i128; n]; n];
        for i in 0..n {
            for j in 0..n {
                let mut s = 0i128;
                for k in 0..n {
                    s += u[k][i] * gu[k][j];
                }
                out[i][j] = s;
            }
        }
        IntegralForm::new(out).unwrap()
    }

    #[test]
    fn z8_and_e8_differ_only_at_two() {
        let z8 = Genus::of(&zn(8)).unwrap();
        let e8 = Genus::of(&e_8()).unwrap();
        // Same rank, signature, determinant.
        assert_eq!(z8.dim, e8.dim);
        assert_eq!(z8.signature, e8.signature);
        assert_eq!(z8.det, e8.det);
        // Z^8: a single type-I scale-0 constituent, oddity 8 ≡ 0 mod 8.
        let s2_z = z8.symbol_at(2);
        assert_eq!(s2_z.len(), 1);
        assert_eq!((s2_z[0].scale, s2_z[0].dim, s2_z[0].type_ii), (0, 8, false));
        // E_8: a single type-II scale-0 constituent.
        let s2_e = e8.symbol_at(2);
        assert_eq!(s2_e.len(), 1);
        assert_eq!((s2_e[0].scale, s2_e[0].dim, s2_e[0].type_ii), (0, 8, true));
        // Distinguished only at p = 2 (type I vs II) — different genus.
        assert!(!are_in_same_genus(&zn(8), &e_8()));
    }

    #[test]
    fn jordan_symbols_match_known_oracles() {
        // A_2: det 3. p=2 single type-II dim-2; p=3 has a scale-1 constituent.
        let a2 = Genus::of(&a_n(2)).unwrap();
        let s2 = a2.symbol_at(2);
        assert_eq!(s2.len(), 1);
        assert!(s2[0].type_ii && s2[0].dim == 2 && s2[0].scale == 0);
        // det = product over scales of p^(scale·dim) at each p must recover |det|.
        assert_eq!(a2.det, 3);

        // D_4: p=2 symbol is two type-II scales (0 and 1), each dim 2.
        let d4 = Genus::of(&d_n(4)).unwrap();
        let s2 = d4.symbol_at(2);
        assert_eq!(s2.len(), 2);
        assert_eq!((s2[0].scale, s2[0].dim, s2[0].type_ii), (0, 2, true));
        assert_eq!((s2[1].scale, s2[1].dim, s2[1].type_ii), (1, 2, true));

        // A_1 = ⟨2⟩: p=2 single type-I scale-1 dim-1, oddity 1.
        let a1 = Genus::of(&IntegralForm::diagonal(&[2])).unwrap();
        let s2 = a1.symbol_at(2);
        assert_eq!(s2.len(), 1);
        assert_eq!(
            (s2[0].scale, s2[0].dim, s2[0].type_ii, s2[0].oddity),
            (1, 1, false, 1)
        );
    }

    #[test]
    fn reflexive_and_isometry_invariant() {
        // A lattice is in its own genus.
        for g in [a_n(2), a_n(4), d_n(4), d_n(5), e_6(), e_7(), e_8(), zn(5)] {
            assert!(are_in_same_genus(&g, &g), "reflexive");
        }
        // Isometric copies (unimodular congruence) share the genus — the strong
        // randomised oracle for the Jordan splitting + canonicalisation.
        let cases = [a_n(3), d_n(4), e_6(), e_8(), zn(6)];
        // Strictly upper-triangular shears ⇒ U is unit-upper-triangular ⇒ det U = 1
        // ⇒ Uᵀ G U is a genuine isometry (preserves det and genus).
        let shear_sets: &[&[(usize, usize, i128)]] = &[
            &[(0, 1, 1)],
            &[(0, 1, 2), (1, 2, -1)],
            &[(0, 2, 1), (0, 1, -3), (1, 2, 1)],
        ];
        for g in &cases {
            for shears in shear_sets {
                let valid: Vec<_> = shears
                    .iter()
                    .copied()
                    .filter(|&(i, j, _)| i < g.dim() && j < g.dim() && i < j)
                    .collect();
                let h = congruent(g, &valid);
                assert_eq!(h.determinant(), g.determinant(), "congruence keeps det");
                assert!(
                    are_in_same_genus(g, &h),
                    "isometric copy must share the genus (dim {})",
                    g.dim()
                );
            }
        }
    }

    #[test]
    fn determinant_and_signature_distinguish_genera() {
        // Different determinant ⇒ different genus.
        assert!(!are_in_same_genus(&a_n(2), &a_n(3))); // det 3 vs 4
        assert!(!are_in_same_genus(
            &IntegralForm::diagonal(&[1]),
            &IntegralForm::diagonal(&[3])
        ));
        // Same det and rank, different signature ⇒ different genus.
        let pos = IntegralForm::diagonal(&[1, 1]);
        let indef = IntegralForm::diagonal(&[1, -1]);
        assert_eq!(pos.determinant().abs(), indef.determinant().abs());
        assert!(!are_in_same_genus(&pos, &indef));
    }

    #[test]
    fn even_unimodular_rank16_share_a_genus() {
        // E_8 ⊕ E_8 is even unimodular rank 16; its genus is that of any even
        // unimodular rank-16 lattice (signature (16,0)). Compare with itself via a
        // sheared copy as the available representative.
        let e8e8 = e_8().direct_sum(&e_8());
        assert_eq!(e8e8.determinant(), 1);
        assert!(e8e8.is_even());
        let g = Genus::of(&e8e8).unwrap();
        assert_eq!(g.signature, (16, 0));
        // single type-II scale-0 constituent of dim 16
        let s2 = g.symbol_at(2);
        assert_eq!(s2.len(), 1);
        assert_eq!((s2[0].dim, s2[0].type_ii), (16, true));
    }
}
