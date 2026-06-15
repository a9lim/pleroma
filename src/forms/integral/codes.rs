//! Binary linear codes and Construction A lattices.
//!
//! This is the finite-code side of the integral lattice story. A binary code
//! `C <= F_2^n` has three compatible readings here:
//!
//! * as a checked F2 row space, with duals and exact weight enumerators;
//! * as a source of integral lattices through Construction A,
//!   `(1/sqrt(2)){x in Z^n : x mod 2 in C}`;
//! * as an exact theta-series oracle through the Hamming weight enumerator.
//!
//! The `1/sqrt(2)` scale is part of the construction. Since [`IntegralForm`]
//! stores an integer Gram matrix, [`BinaryCode::construction_a`] returns `None`
//! unless the resulting Gram is integral; self-orthogonal codes satisfy that
//! boundary. Type I self-dual codes give odd unimodular lattices, while Type II
//! self-dual codes give even unimodular lattices.

use super::lattice::IntegralForm;
use crate::linalg::integer::normalize_relation_rows;

/// `|Aut(D16+)| = 2^15 * 16!`.
pub const D16_PLUS_AUT_ORDER: u128 = 685_597_979_049_984_000;

/// A binary linear code, stored as a row-reduced F2 generator matrix.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BinaryCode {
    n: usize,
    generators: Vec<Vec<u8>>,
}

fn row_weight(row: &[u8]) -> usize {
    row.iter().map(|&x| x as usize).sum()
}

fn dot_mod2(a: &[u8], b: &[u8]) -> u8 {
    a.iter().zip(b).fold(0u8, |acc, (&x, &y)| acc ^ (x & y))
}

fn normalize_generators(mut rows: Vec<Vec<u8>>, n: usize) -> Option<Vec<Vec<u8>>> {
    if rows
        .iter()
        .any(|row| row.len() != n || row.iter().any(|&x| x > 1))
    {
        return None;
    }
    rows.retain(|row| row.iter().any(|&x| x != 0));
    let mut rank = 0usize;
    for col in 0..n {
        let Some(pivot) = (rank..rows.len()).find(|&r| rows[r][col] != 0) else {
            continue;
        };
        rows.swap(rank, pivot);
        let pivot_row = rows[rank].clone();
        for r in 0..rows.len() {
            if r == rank || rows[r][col] == 0 {
                continue;
            }
            for c in col..n {
                rows[r][c] ^= pivot_row[c];
            }
        }
        rank += 1;
        if rank == rows.len() {
            break;
        }
    }
    rows.truncate(rank);
    Some(rows)
}

fn rows_from_strings(rows: &[&str]) -> Vec<Vec<u8>> {
    rows.iter()
        .map(|row| {
            row.bytes()
                .map(|b| match b {
                    b'0' => 0,
                    b'1' => 1,
                    _ => panic!("binary generator rows must contain only 0/1"),
                })
                .collect()
        })
        .collect()
}

fn binomial(n: usize, k: usize) -> i128 {
    if k > n {
        return 0;
    }
    let k = k.min(n - k);
    let mut out = 1i128;
    for i in 1..=k {
        out = out
            .checked_mul((n - k + i) as i128)
            .expect("binomial coefficient exceeds i128")
            / i as i128;
    }
    out
}

fn convolve_i128(a: &[i128], b: &[i128], terms: usize) -> Vec<i128> {
    let mut out = vec![0i128; terms];
    for (i, &ai) in a.iter().enumerate().take(terms) {
        if ai == 0 {
            continue;
        }
        for (j, &bj) in b.iter().enumerate().take(terms - i) {
            if bj == 0 {
                continue;
            }
            out[i + j] = out[i + j]
                .checked_add(ai.checked_mul(bj).expect("series coefficient exceeds i128"))
                .expect("series coefficient exceeds i128");
        }
    }
    out
}

fn series_pow_i128(base: &[i128], exp: usize, terms: usize) -> Vec<i128> {
    let mut out = vec![0i128; terms];
    if terms == 0 {
        return out;
    }
    out[0] = 1;
    for _ in 0..exp {
        out = convolve_i128(&out, base, terms);
    }
    out
}

fn even_residue_theta(terms: usize) -> Vec<i128> {
    let mut out = vec![0i128; terms];
    if terms == 0 {
        return out;
    }
    out[0] = 1;
    let mut m = 1usize;
    while m * m < terms {
        out[m * m] += 2;
        m += 1;
    }
    out
}

fn odd_residue_theta_without_quarter(terms: usize) -> Vec<i128> {
    let mut out = vec![0i128; terms];
    let mut m = 0usize;
    while m * (m + 1) < terms {
        out[m * (m + 1)] += 2; // m and -m-1.
        m += 1;
    }
    out
}

impl BinaryCode {
    /// Build a binary code from generator rows. The stored basis is row-reduced
    /// over F2, so equivalent generator matrices compare equal.
    pub fn new(n: usize, generators: Vec<Vec<u8>>) -> Option<Self> {
        Some(BinaryCode {
            n,
            generators: normalize_generators(generators, n)?,
        })
    }

    /// The block length `n`.
    pub fn len(&self) -> usize {
        self.n
    }

    /// Whether the code has block length zero.
    pub fn is_empty(&self) -> bool {
        self.n == 0
    }

    /// The dimension `k`.
    pub fn dim(&self) -> usize {
        self.generators.len()
    }

    /// Row-reduced generator rows.
    pub fn generators(&self) -> &[Vec<u8>] {
        &self.generators
    }

    /// The number of codewords, `2^k`, when it fits the crate's `u128` payload.
    pub fn size(&self) -> Option<u128> {
        if self.dim() >= 128 {
            None
        } else {
            Some(1u128 << self.dim())
        }
    }

    fn codewords(&self) -> Vec<Vec<u8>> {
        assert!(
            self.dim() < usize::BITS as usize,
            "codeword enumeration is exponential and exceeds usize masks"
        );
        let mut out = Vec::with_capacity(1usize << self.dim());
        for mask in 0usize..(1usize << self.dim()) {
            let mut word = vec![0u8; self.n];
            for (i, row) in self.generators.iter().enumerate() {
                if (mask >> i) & 1 == 0 {
                    continue;
                }
                for j in 0..self.n {
                    word[j] ^= row[j];
                }
            }
            out.push(word);
        }
        out
    }

    /// The dual code `C^perp = {x : x dot c = 0 for all c in C}`.
    pub fn dual(&self) -> BinaryCode {
        let mut pivot_for_row = Vec::new();
        let mut is_pivot = vec![false; self.n];
        for row in &self.generators {
            if let Some(p) = row.iter().position(|&x| x != 0) {
                pivot_for_row.push(p);
                is_pivot[p] = true;
            }
        }

        let mut dual_rows = Vec::new();
        for free in 0..self.n {
            if is_pivot[free] {
                continue;
            }
            let mut v = vec![0u8; self.n];
            v[free] = 1;
            for (r, &pivot) in pivot_for_row.iter().enumerate() {
                v[pivot] = self.generators[r][free];
            }
            dual_rows.push(v);
        }
        BinaryCode::new(self.n, dual_rows).expect("dual rows have the same length")
    }

    /// The block direct sum `C ⊕ D`.
    pub fn direct_sum(&self, other: &BinaryCode) -> BinaryCode {
        let mut rows = Vec::with_capacity(self.dim() + other.dim());
        for row in &self.generators {
            let mut out = vec![0u8; self.n + other.n];
            out[..self.n].copy_from_slice(row);
            rows.push(out);
        }
        for row in &other.generators {
            let mut out = vec![0u8; self.n + other.n];
            out[self.n..].copy_from_slice(row);
            rows.push(out);
        }
        BinaryCode::new(self.n + other.n, rows).expect("direct-sum rows are binary")
    }

    /// `C = C^perp`.
    pub fn is_self_dual(&self) -> bool {
        self.dim() * 2 == self.n && self.generators == self.dual().generators
    }

    /// `C <= C^perp`.
    pub fn is_self_orthogonal(&self) -> bool {
        (0..self.dim()).all(|i| {
            (i..self.dim()).all(|j| dot_mod2(&self.generators[i], &self.generators[j]) == 0)
        })
    }

    /// Every codeword has Hamming weight divisible by 4.
    pub fn is_doubly_even(&self) -> bool {
        if self
            .generators
            .iter()
            .any(|row| !row_weight(row).is_multiple_of(4))
        {
            return false;
        }
        (0..self.dim()).all(|i| {
            (i + 1..self.dim()).all(|j| dot_mod2(&self.generators[i], &self.generators[j]) == 0)
        })
    }

    /// The minimum nonzero Hamming weight, or `None` for the zero code.
    pub fn minimum_distance(&self) -> Option<usize> {
        self.codewords()
            .into_iter()
            .map(|word| row_weight(&word))
            .filter(|&w| w > 0)
            .min()
    }

    /// The Hamming weight enumerator coefficients:
    /// `out[w] = #{c in C : wt(c) = w}`.
    pub fn weight_enumerator(&self) -> Vec<i128> {
        let mut out = vec![0i128; self.n + 1];
        for word in self.codewords() {
            out[row_weight(&word)] += 1;
        }
        out
    }

    /// The exact MacWilliams transform of the Hamming weight enumerator. The
    /// result is the weight enumerator of `C^perp`.
    pub fn macwilliams_transform(&self) -> Vec<i128> {
        let a = self.weight_enumerator();
        let size = i128::try_from(self.size().expect("code size exceeds u128"))
            .expect("code size exceeds i128");
        let mut out = vec![0i128; self.n + 1];
        for (j, out_j) in out.iter_mut().enumerate() {
            let mut acc = 0i128;
            for (i, &ai) in a.iter().enumerate() {
                if ai == 0 {
                    continue;
                }
                let mut kraw = 0i128;
                for s in 0..=j {
                    let sign = if s % 2 == 0 { 1 } else { -1 };
                    kraw = kraw
                        .checked_add(
                            sign * binomial(i, s)
                                .checked_mul(binomial(self.n - i, j - s))
                                .expect("Krawtchouk coefficient exceeds i128"),
                        )
                        .expect("Krawtchouk coefficient exceeds i128");
                }
                acc = acc
                    .checked_add(ai.checked_mul(kraw).expect("MacWilliams sum exceeds i128"))
                    .expect("MacWilliams sum exceeds i128");
            }
            debug_assert_eq!(acc % size, 0);
            *out_j = acc / size;
        }
        out
    }

    /// Construction A with the standard `1/sqrt(2)` scaling. Returns `None`
    /// exactly when the scaled Gram matrix is not integral.
    pub fn construction_a(&self) -> Option<IntegralForm> {
        let mut rows: Vec<Vec<i128>> = self
            .generators
            .iter()
            .map(|row| row.iter().map(|&x| x as i128).collect())
            .collect();
        for i in 0..self.n {
            let mut row = vec![0i128; self.n];
            row[i] = 2;
            rows.push(row);
        }
        let basis = normalize_relation_rows(rows);
        if basis.len() != self.n {
            return None;
        }
        let mut gram = vec![vec![0i128; self.n]; self.n];
        for i in 0..self.n {
            for j in 0..self.n {
                let dot: i128 = basis[i].iter().zip(&basis[j]).map(|(&x, &y)| x * y).sum();
                if dot % 2 != 0 {
                    return None;
                }
                gram[i][j] = dot / 2;
            }
        }
        IntegralForm::new(gram)
    }

    /// Compute the Construction A theta series from the Hamming weight
    /// enumerator and the one-dimensional even/odd residue theta series.
    ///
    /// This returns `None` outside the doubly-even boundary, where the
    /// quarter-powers from odd residue classes do not assemble into integer
    /// exponents.
    pub fn theta_series_via_weight_enumerator(&self, terms: usize) -> Option<Vec<i128>> {
        if !self.is_doubly_even() {
            return None;
        }
        if terms == 0 {
            return Some(Vec::new());
        }
        let weights = self.weight_enumerator();
        let even = even_residue_theta(terms);
        let odd = odd_residue_theta_without_quarter(terms);
        let mut out = vec![0i128; terms];
        for (w, &count) in weights.iter().enumerate() {
            if count == 0 {
                continue;
            }
            debug_assert_eq!(w % 4, 0);
            let shift = w / 4;
            if shift >= terms {
                continue;
            }
            let even_part = series_pow_i128(&even, self.n - w, terms - shift);
            let odd_part = series_pow_i128(&odd, w, terms - shift);
            let product = convolve_i128(&even_part, &odd_part, terms - shift);
            for (i, &coeff) in product.iter().enumerate() {
                out[i + shift] = out[i + shift]
                    .checked_add(
                        count
                            .checked_mul(coeff)
                            .expect("Construction A theta coefficient exceeds i128"),
                    )
                    .expect("Construction A theta coefficient exceeds i128");
            }
        }
        Some(out)
    }
}

/// The binary Hamming `[7,4,3]` code.
pub fn hamming_code() -> BinaryCode {
    BinaryCode::new(
        7,
        rows_from_strings(&["1000011", "0100101", "0010110", "0001111"]),
    )
    .expect("Hamming generator is binary")
}

/// The binary repetition `[n,1,n]` code, for `n > 0`.
pub fn repetition_code(n: usize) -> Option<BinaryCode> {
    if n == 0 {
        return None;
    }
    BinaryCode::new(n, vec![vec![1u8; n]])
}

/// The Type I self-dual `[2,1,2]` code; Construction A gives an odd unimodular
/// rank-2 lattice isometric to `Z^2`.
pub fn type_i_z2_code() -> BinaryCode {
    repetition_code(2).expect("length-2 repetition code exists")
}

/// A Type I self-dual code whose Construction A lattice is `Z^2 ⊕ E8`.
pub fn type_i_z2_plus_e8_code() -> BinaryCode {
    type_i_z2_code().direct_sum(&extended_hamming_code())
}

/// The extended Hamming `[8,4,4]` Type II code; Construction A gives `E8`.
pub fn extended_hamming_code() -> BinaryCode {
    BinaryCode::new(
        8,
        rows_from_strings(&["11110000", "11001100", "10101010", "11111111"]),
    )
    .expect("extended Hamming generator is binary")
}

/// The direct sum of two extended Hamming codes; Construction A gives
/// `E8 + E8`.
pub fn type_ii_e8_sum_code() -> BinaryCode {
    let mut rows = Vec::new();
    for row in extended_hamming_code().generators() {
        let mut r = vec![0u8; 16];
        r[..8].copy_from_slice(row);
        rows.push(r);
    }
    for row in extended_hamming_code().generators() {
        let mut r = vec![0u8; 16];
        r[8..].copy_from_slice(row);
        rows.push(r);
    }
    BinaryCode::new(16, rows).expect("direct sum generator is binary")
}

/// The indecomposable Type II `[16,8,4]` code whose Construction A lattice is
/// `D16+`.
pub fn type_ii_len16_code() -> BinaryCode {
    // One of the two length-16 Type II generators in the Harada-Munemasa
    // self-dual-code tables; the other is `E8` plus `E8`.
    BinaryCode::new(
        16,
        rows_from_strings(&[
            "1000000000001101",
            "0100001110011110",
            "0010001110011011",
            "0001000010000101",
            "0000100100000101",
            "0000011000000101",
            "0000000001010101",
            "0000000000100111",
        ]),
    )
    .expect("length-16 Type II generator is binary")
}

/// The `D16+` lattice, built by Construction A from the indecomposable Type II
/// length-16 code.
pub fn d16_plus() -> IntegralForm {
    type_ii_len16_code()
        .construction_a()
        .expect("Type II Construction A is integral")
}

/// The extended binary Golay `[24,12,8]` code.
pub fn golay_code() -> BinaryCode {
    BinaryCode::new(24, extended_golay_generator_rows()).expect("Golay generator is binary")
}

pub(crate) fn extended_golay_generator_rows() -> Vec<Vec<u8>> {
    let a: [[u8; 12]; 12] = [
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
            let mut row = vec![0u8; 24];
            row[i] = 1;
            row[12..24].copy_from_slice(&a[i]);
            row
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forms::e_8;

    #[test]
    fn hamming_macwilliams_matches_dual() {
        let h = hamming_code();
        assert_eq!(h.len(), 7);
        assert_eq!(h.dim(), 4);
        assert_eq!(h.minimum_distance(), Some(3));
        assert_eq!(h.macwilliams_transform(), h.dual().weight_enumerator());
        assert!(h.construction_a().is_none());
    }

    #[test]
    fn direct_sum_and_repetition_codes_build_type_i_examples() {
        assert!(repetition_code(0).is_none());
        let r2 = type_i_z2_code();
        assert_eq!(r2.len(), 2);
        assert_eq!(r2.dim(), 1);
        assert!(r2.is_self_dual());
        assert!(r2.is_self_orthogonal());
        assert!(!r2.is_doubly_even());
        assert_eq!(r2.weight_enumerator(), vec![1, 0, 1]);

        let z2 = r2.construction_a().unwrap();
        assert_eq!(z2.determinant(), 1);
        assert!(!z2.is_even());
        assert_eq!(z2.minimum(), Some(1));

        let z2_e8 = type_i_z2_plus_e8_code();
        assert_eq!(z2_e8.len(), 10);
        assert_eq!(z2_e8.dim(), 5);
        assert!(z2_e8.is_self_dual());
        assert!(!z2_e8.is_doubly_even());
        let lattice = z2_e8.construction_a().unwrap();
        assert_eq!(lattice.dim(), 10);
        assert!(lattice.is_unimodular());
        assert!(!lattice.is_even());
    }

    #[test]
    fn golay_macwilliams_and_construction_a_boundary() {
        let g = golay_code();
        assert_eq!(g.len(), 24);
        assert_eq!(g.dim(), 12);
        assert_eq!(g.minimum_distance(), Some(8));
        assert!(g.is_self_dual());
        assert!(g.is_doubly_even());
        assert_eq!(g.macwilliams_transform(), g.weight_enumerator());

        let l = g.construction_a().expect("Golay is Type II");
        assert_eq!(l.dim(), 24);
        assert!(l.is_even());
        assert!(l.is_unimodular());
        assert_eq!(g.theta_series_via_weight_enumerator(2), Some(vec![1, 48]));
    }

    #[test]
    fn type_ii_codes_build_the_rank_8_and_rank_16_lattices() {
        let e8_code = extended_hamming_code();
        assert!(e8_code.is_self_dual());
        assert!(e8_code.is_doubly_even());
        let e8_from_code = e8_code.construction_a().unwrap();
        assert_eq!(e8_from_code.dim(), 8);
        assert!(e8_from_code.is_even());
        assert!(e8_from_code.is_unimodular());
        assert_eq!(e8_from_code.determinant(), e_8().determinant());
        assert_eq!(
            e8_code.theta_series_via_weight_enumerator(2),
            Some(vec![1, 240])
        );

        let split = type_ii_e8_sum_code().construction_a().unwrap();
        assert_eq!(split.dim(), 16);
        assert!(split.is_even());
        assert!(split.is_unimodular());
        assert_eq!(split.determinant(), e_8().direct_sum(&e_8()).determinant());

        let d16 = d16_plus();
        assert_eq!(d16.dim(), 16);
        assert!(d16.is_even());
        assert!(d16.is_unimodular());
        assert_eq!(
            type_ii_len16_code().theta_series_via_weight_enumerator(2),
            Some(vec![1, 480])
        );
        assert_eq!(
            D16_PLUS_AUT_ORDER,
            (1u128 << 15) * (1..=16u128).product::<u128>()
        );
    }

    #[test]
    fn weight_enumerator_theta_matches_construction_a_theta() {
        let e8_code = extended_hamming_code();
        let e8_lattice = e8_code.construction_a().unwrap();
        assert_eq!(
            e8_code.theta_series_via_weight_enumerator(3),
            e8_lattice.theta_series(3)
        );

        assert_eq!(
            type_ii_e8_sum_code().theta_series_via_weight_enumerator(3),
            Some(vec![1, 480, 61920])
        );
        assert_eq!(
            type_ii_len16_code().theta_series_via_weight_enumerator(3),
            Some(vec![1, 480, 61920])
        );
        assert_eq!(
            golay_code().theta_series_via_weight_enumerator(2),
            Some(vec![1, 48])
        );
    }
}
