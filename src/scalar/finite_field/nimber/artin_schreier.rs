//! Trace and Artin-Schreier linear algebra over finite nim-fields.

use super::nim_square;

/// Field trace F_{2^m} → F₂:  `Tr(x) = x + x² + x⁴ + … + x^{2^{m-1}} ∈ {0,1}`.
/// This is the canonical map realising k/℘(k) ≅ F₂ that the Arf invariant is read
/// through (see `forms::char2`); it is *also* the obstruction to solving the
/// Artin–Schreier equation `y² + y = c` (solvable iff `Tr(c) = 0`). One trace,
/// both roles — that is the unification. `m` must be the degree of a nim-subfield
/// (a power of two: 1, 2, 4, …, 128).
pub fn nim_trace(x: u128, m: u128) -> u128 {
    let mut acc = x;
    let mut t = x;
    for _ in 1..m {
        t = nim_square(t);
        acc ^= t;
    }
    acc
}

/// Insert `val` (with its associated y-combination `yc`) into an XOR pivot table
/// keyed by highest set bit. Used by the Artin–Schreier solver.
fn xor_basis_insert(table: &mut [Option<(u128, u128)>; 128], mut val: u128, mut yc: u128) {
    while val != 0 {
        let h = (127 - val.leading_zeros()) as usize;
        match table[h] {
            Some((pv, pc)) => {
                val ^= pv;
                yc ^= pc;
            }
            None => {
                table[h] = Some((val, yc));
                return;
            }
        }
    }
}

/// Solve the Artin–Schreier equation `y² + y = c` in F_{2^m}. The map
/// `L(y) = y² + y` is F₂-linear with kernel {0,1}, and its image is exactly the
/// trace-zero hyperplane — so a solution exists **iff `nim_trace(c, m) = 0`**, and
/// when it does there are exactly two (`y` and `y+1`). Returns one solution, or
/// `None` when `c` is not in the image. Solved by Gaussian elimination over F₂ on
/// the bit-basis of F_{2^m} (exact; no fragile closed-form).
pub fn nim_solve_artin_schreier(c: u128, m: u128) -> Option<u128> {
    let mut table: [Option<(u128, u128)>; 128] = [None; 128];
    for k in 0..m {
        let e = 1u128 << k;
        xor_basis_insert(&mut table, nim_square(e) ^ e, e);
    }
    let (mut val, mut yc) = (c, 0u128);
    while val != 0 {
        let h = (127 - val.leading_zeros()) as usize;
        match table[h] {
            Some((pv, pc)) => {
                val ^= pv;
                yc ^= pc;
            }
            None => return None, // c ∉ image(L)  ⇔  Tr(c) ≠ 0
        }
    }
    Some(yc)
}

/// Whether `y² + y = c` is solvable in F_{2^m} — i.e. `Tr(c) = 0`. The same
/// trace, hence the same answer, as the Arf-reduction path.
pub fn nim_is_artin_schreier_solvable(c: u128, m: u128) -> bool {
    nim_trace(c, m) == 0
}
