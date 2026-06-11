use super::basis::wedge_sign;
use super::metric::Metric;
use super::terms::{add_term, merge, scale};
use crate::scalar::Scalar;
use std::collections::BTreeMap;

impl<S: Scalar> Metric<S> {
    fn a_val(&self, i: usize, j: usize) -> S {
        self.a.get(&(i, j)).cloned().unwrap_or_else(S::zero)
    }

    fn b_val(&self, i: usize, j: usize) -> S {
        let key = if i < j { (i, j) } else { (j, i) };
        self.b.get(&key).cloned().unwrap_or_else(S::zero)
    }

    /// The full bilinear form `B(e_i, e_j)`, reconstructed from (q, b, a).
    fn bil(&self, i: usize, j: usize) -> S {
        use std::cmp::Ordering::*;
        match i.cmp(&j) {
            Equal => self.q_val(i),
            Less => self.a_val(i, j),
            Greater => self.b_val(j, i).sub(&self.a_val(j, i)),
        }
    }

    /// The B-contraction `e_i ⌟_B W_T` only (no wedge term).
    fn contract_vec_blade(&self, i: usize, t: u128) -> BTreeMap<u128, S> {
        let mut out: BTreeMap<u128, S> = BTreeMap::new();
        let mut tt = t;
        let mut k = 0u128;
        while tt != 0 {
            let j = tt.trailing_zeros() as usize;
            tt &= tt - 1;
            let c = self.bil(i, j);
            if !c.is_zero() {
                let coeff = if k & 1 == 0 { c } else { c.neg() };
                add_term(&mut out, t ^ (1 << j), coeff);
            }
            k += 1;
        }
        out
    }

    /// `e_i · W_T` in the wedge basis (Chevalley).
    fn vec_times_blade(&self, i: usize, t: u128) -> BTreeMap<u128, S> {
        let mut out = self.contract_vec_blade(i, t);
        if t & (1 << i) == 0 {
            let sign = wedge_sign::<S>(1 << i, t);
            add_term(&mut out, t | (1 << i), sign);
        }
        out
    }

    fn vec_times_mv(&self, i: usize, mv: &BTreeMap<u128, S>) -> BTreeMap<u128, S> {
        let mut out: BTreeMap<u128, S> = BTreeMap::new();
        for (&t, c) in mv {
            merge(&mut out, scale(self.vec_times_blade(i, t), c));
        }
        out
    }

    /// Fast path for an orthogonal basis: `e_S e_T` is a single blade `e_{S△T}`
    /// times the reordering sign and the product of `q_i` over repeated indices.
    fn geom_product_blades_orthogonal(&self, s: u128, t: u128) -> BTreeMap<u128, S> {
        let mut coeff = wedge_sign::<S>(s, t);
        let mut common = s & t;
        while common != 0 {
            let i = common.trailing_zeros() as usize;
            common &= common - 1;
            coeff = coeff.mul(&self.q_val(i));
            if coeff.is_zero() {
                return BTreeMap::new();
            }
        }
        let mut m = BTreeMap::new();
        m.insert(s ^ t, coeff);
        m
    }

    /// The general-bilinear-form geometric product of two wedge blades.
    pub(super) fn geom_product_blades(&self, s: u128, t: u128) -> BTreeMap<u128, S> {
        if self.is_orthogonal() {
            return self.geom_product_blades_orthogonal(s, t);
        }
        if s == 0 {
            let mut m = BTreeMap::new();
            m.insert(t, S::one());
            return m;
        }
        let i = s.trailing_zeros() as usize;
        let s_rest = s ^ (1 << i);
        let xy = self.geom_product_blades(s_rest, t);
        let part1 = self.vec_times_mv(i, &xy);
        let contraction = self.contract_vec_blade(i, s_rest);
        let mut part2: BTreeMap<u128, S> = BTreeMap::new();
        for (&u, cu) in &contraction {
            merge(&mut part2, scale(self.geom_product_blades(u, t), cu));
        }
        let mut out = part1;
        merge(&mut out, scale(part2, &S::one().neg()));
        out
    }

    /// Independent swap/contract reduction retained as a test oracle.
    #[cfg(test)]
    pub(super) fn reduce_word(&self, word: &[usize]) -> BTreeMap<u128, S> {
        for p in 0..word.len().saturating_sub(1) {
            let (a, c) = (word[p], word[p + 1]);
            if a == c {
                let q = self.q_val(a);
                let mut rest = Vec::with_capacity(word.len() - 2);
                rest.extend_from_slice(&word[..p]);
                rest.extend_from_slice(&word[p + 2..]);
                return scale(self.reduce_word(&rest), &q);
            } else if a > c {
                let bv = self.b_val(a, c);
                let mut removed = Vec::with_capacity(word.len() - 2);
                removed.extend_from_slice(&word[..p]);
                removed.extend_from_slice(&word[p + 2..]);
                let mut out = scale(self.reduce_word(&removed), &bv);

                let mut swapped = word.to_vec();
                swapped.swap(p, p + 1);
                let neg = S::one().neg();
                merge(&mut out, scale(self.reduce_word(&swapped), &neg));
                return out;
            }
        }
        let mut mask = 0u128;
        for &g in word {
            mask |= 1 << g;
        }
        let mut m = BTreeMap::new();
        m.insert(mask, S::one());
        m
    }
}
