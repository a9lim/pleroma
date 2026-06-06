//! The multivector engine, generic over any `Scalar` backend.
//!
//! ## Metric data — characteristic-faithful by design
//!
//! A blade is a `u32` bitmask over basis generators e_0..e_31. The algebra is
//! defined by two independent pieces of data, *not* a single bilinear form:
//!
//!   * `q[i]`      = e_i²                      (the quadratic form / squares)
//!   * `b[(i,j)]`  = e_i e_j + e_j e_i  (i<j)  (the polar / anticommutator form)
//!
//! In characteristic ≠ 2 these are linked (`b = 2·offdiag`, `q = diag`), so an
//! orthogonal basis just sets `b = 0`. In characteristic 2 they are genuinely
//! independent: the polar form is *alternating* (`b(i,i)=0`) yet `q[i]` can be
//! nonzero, and a nonzero off-diagonal `b[(i,j)]` is exactly what makes the
//! nim-Clifford algebra *non-commutative*. Carrying both is the faithful thing.
//!
//! "With nilpotents": set `q[i] = 0` and you get a null generator, e_i² = 0.
//! All `q = 0`, all `b = 0` ⇒ the exterior/Grassmann algebra.
//!
//! ## Product
//!
//! Two canonical blades multiply by concatenating their (ascending) generator
//! lists into a word and reducing to canonical form with the rules
//!   e_i e_i  → q[i]                            (equal adjacent: contract)
//!   e_i e_j  → b[(j,i)] − e_j e_i   (i>j)      (out of order: swap, emit polar)
//! The `−` goes through the scalar's own `neg()`, so in characteristic 2 it is
//! `+` automatically and signs vanish — no special-casing. Termination: each
//! step lowers (word length, inversion count) lexicographically.

use crate::scalar::Scalar;
use std::collections::BTreeMap;


/// Ascending list of set-bit indices of a blade mask.
pub fn bits(mask: u32) -> Vec<usize> {
    let mut v = Vec::new();
    let mut m = mask;
    while m != 0 {
        let i = m.trailing_zeros() as usize;
        v.push(i);
        m &= m - 1;
    }
    v
}

/// The grade (number of generators) of a blade mask.
pub fn grade(mask: u32) -> u32 {
    mask.count_ones()
}

/// Sign (+1/-1 as a Scalar) of reordering two disjoint ascending blades when
/// concatenated — i.e. the number of (i in a, j in b) with i > j, mod 2.
fn wedge_sign<S: Scalar>(a: u32, b: u32) -> S {
    let mut swaps = 0u32;
    let mut aa = a;
    while aa != 0 {
        let i = aa.trailing_zeros();
        aa &= aa - 1;
        // count bits of b strictly below i
        let below = b & ((1u32 << i) - 1);
        swaps += below.count_ones();
    }
    if swaps & 1 == 0 {
        S::one()
    } else {
        S::one().neg()
    }
}

/// The metric of a Clifford algebra of a (possibly general, possibly degenerate)
/// bilinear form. We store the full bilinear form `B(e_i, e_j)` factored into
/// three pieces so the common case is cheap and the general case is reachable:
///
///   * `q[i]   = B(e_i, e_i) = e_i²`                  — the quadratic form (diagonal)
///   * `b[(i,j)] = B(e_i,e_j) + B(e_j,e_i) = {e_i,e_j}`  (i<j) — the *polar* /
///     anticommutator form (always symmetric; carried independently of `q`,
///     which is the whole point in char 2 — see the module docs).
///   * `a[(i,j)] = B(e_i, e_j)` for i<j — the **strictly-upper / in-order
///     contraction** part. This is what promotes the engine from "Clifford of a
///     *symmetric* polar form" to "Clifford of a *general* bilinear form"
///     (Chevalley/Fauser): the geometric product gains an in-order scalar
///     `e_i e_j = e_i∧e_j + a[(i,j)]` for i<j. `a` empty ⇒ the ordinary Clifford
///     algebra (`e_i e_j = e_i∧e_j` for i<j); a nonzero antisymmetric choice of
///     `a` interpolates toward the Weyl algebra. `b` stays the anticommutator
///     regardless of `a` (the lower entry is `B(e_j,e_i) = b - a`).
#[derive(Clone, Debug)]
pub struct Metric<S: Scalar> {
    pub q: Vec<S>,
    pub b: BTreeMap<(usize, usize), S>,
    pub a: BTreeMap<(usize, usize), S>,
}

impl<S: Scalar> Metric<S> {
    /// Orthogonal metric from a list of squares (b = 0). `Cl(p,q,r)` style.
    pub fn diagonal(q: Vec<S>) -> Self {
        Metric {
            q,
            b: BTreeMap::new(),
            a: BTreeMap::new(),
        }
    }

    /// The fully-null metric: exterior/Grassmann algebra on `n` generators.
    pub fn grassmann(n: usize) -> Self {
        Metric {
            q: vec![S::zero(); n],
            b: BTreeMap::new(),
            a: BTreeMap::new(),
        }
    }

    /// A symmetric-polar Clifford metric: squares `q` and anticommutators `b`
    /// (i<j), with no in-order contraction (`a` empty). The ordinary case.
    pub fn new(q: Vec<S>, b: BTreeMap<(usize, usize), S>) -> Self {
        Metric {
            q,
            b,
            a: BTreeMap::new(),
        }
    }

    /// A general-bilinear-form metric: squares `q`, polar form `b` (i<j), and the
    /// in-order contraction `a` (i<j). See the struct docs.
    pub fn general(
        q: Vec<S>,
        b: BTreeMap<(usize, usize), S>,
        a: BTreeMap<(usize, usize), S>,
    ) -> Self {
        Metric { q, b, a }
    }

    /// True iff there is any in-order contraction — i.e. this is a genuinely
    /// general bilinear form and needs the Chevalley product path.
    pub(crate) fn has_upper(&self) -> bool {
        self.a.values().any(|v| !v.is_zero())
    }

    fn a_val(&self, i: usize, j: usize) -> S {
        // a is stored for i<j only.
        self.a.get(&(i, j)).cloned().unwrap_or_else(S::zero)
    }

    /// The full bilinear form `B(e_i, e_j)`, reconstructed from (q, b, a):
    ///   * `i == j`: `q_i`
    ///   * `i  < j`: `a_{ij}`                 (the in-order / upper contraction)
    ///   * `i  > j`: `b_{ji} − a_{ji}`        (the lower entry; b is the symmetric
    ///                                         polar form, so `B(e_j,e_i)=b−B(e_i,e_j)`)
    fn bil(&self, i: usize, j: usize) -> S {
        use std::cmp::Ordering::*;
        match i.cmp(&j) {
            Equal => self.q_val(i),
            Less => self.a_val(i, j),
            Greater => self.b_val(j, i).sub(&self.a_val(j, i)),
        }
    }

    /// The B-contraction `e_i ⌟_B W_T` only (no wedge term): the multivector
    /// `Σ_k (−1)^k B(e_i, e_{j_k}) W_{T∖j_k}` over the set bits `j_k` of `t` in
    /// ascending order. Used by the general (Chevalley) geometric product.
    fn contract_vec_blade(&self, i: usize, t: u32) -> BTreeMap<u32, S> {
        let mut out: BTreeMap<u32, S> = BTreeMap::new();
        let mut tt = t;
        let mut k = 0u32;
        while tt != 0 {
            let j = tt.trailing_zeros() as usize;
            tt &= tt - 1;
            let c = self.bil(i, j);
            if !c.is_zero() {
                let coeff = if k & 1 == 0 { c } else { c.neg() };
                let e = out.entry(t ^ (1 << j)).or_insert_with(S::zero);
                *e = e.add(&coeff);
                if e.is_zero() {
                    out.remove(&(t ^ (1 << j)));
                }
            }
            k += 1;
        }
        out
    }

    /// `e_i · W_T` in the wedge basis (Chevalley): wedge part `e_i ∧ W_T` (when
    /// `i ∉ T`) plus the B-contraction `e_i ⌟_B W_T`.
    fn vec_times_blade(&self, i: usize, t: u32) -> BTreeMap<u32, S> {
        let mut out = self.contract_vec_blade(i, t);
        if t & (1 << i) == 0 {
            let sign = wedge_sign::<S>(1 << i, t);
            let e = out.entry(t | (1 << i)).or_insert_with(S::zero);
            *e = e.add(&sign);
            if e.is_zero() {
                out.remove(&(t | (1 << i)));
            }
        }
        out
    }

    fn vec_times_mv(&self, i: usize, mv: &BTreeMap<u32, S>) -> BTreeMap<u32, S> {
        let mut out: BTreeMap<u32, S> = BTreeMap::new();
        for (&t, c) in mv {
            merge(&mut out, scale(self.vec_times_blade(i, t), c));
        }
        out
    }

    /// The general-bilinear-form geometric product of two wedge blades,
    /// `W_S · W_T`, via the Chevalley recursion
    ///   `(e_i ∧ X) · Y = e_i · (X · Y) − (e_i ⌟_B X) · Y`,
    /// peeling the minimum generator `i` of `S` (so `W_S = e_i ∧ W_{S∖i}` with no
    /// sign). Recurses on a strictly smaller leading blade ⇒ terminates. With `a`
    /// empty (`B` lower-triangular) this reproduces `reduce_word` exactly.
    fn geom_product_blades(&self, s: u32, t: u32) -> BTreeMap<u32, S> {
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
        let mut part2: BTreeMap<u32, S> = BTreeMap::new();
        for (&u, cu) in &contraction {
            merge(&mut part2, scale(self.geom_product_blades(u, t), cu));
        }
        let mut out = part1;
        merge(&mut out, scale(part2, &S::one().neg()));
        out
    }

    /// Orthogonal direct sum `M ⟂ M'`: a block-diagonal metric on the disjoint
    /// union of the two generator sets (the second block's indices are shifted up
    /// by `self`'s dimension). The two blocks are mutually orthogonal — cross
    /// polar form 0 — so cross generators anticommute (char 0) / commute (char 2).
    /// These are exactly the relations of the Clifford algebra of the orthogonal
    /// sum of the two forms, i.e. the **graded tensor product** of the algebras;
    /// and the Arf invariant is additive over it.
    pub fn direct_sum(&self, other: &Metric<S>) -> Metric<S> {
        let n = self.q.len();
        let mut q = self.q.clone();
        q.extend(other.q.iter().cloned());
        let mut b = self.b.clone();
        for (&(i, j), v) in &other.b {
            b.insert((i + n, j + n), v.clone());
        }
        let mut a = self.a.clone();
        for (&(i, j), v) in &other.a {
            a.insert((i + n, j + n), v.clone());
        }
        Metric { q, b, a }
    }

    pub(crate) fn q_val(&self, i: usize) -> S {
        self.q.get(i).cloned().unwrap_or_else(S::zero)
    }

    fn b_val(&self, i: usize, j: usize) -> S {
        let key = if i < j { (i, j) } else { (j, i) };
        self.b.get(&key).cloned().unwrap_or_else(S::zero)
    }

    /// Reduce a generator word to canonical multivector terms — the original
    /// swap/contract reduction for the *symmetric-polar* case (`a` empty). It is
    /// now retained solely as the **independent oracle** the general Chevalley
    /// product (`geom_product_blades`) is cross-validated against in the tests
    /// (`general_product_reproduces_reduce_word_when_a_empty`), so the production
    /// engine has a second, structurally different implementation to agree with.
    #[cfg(test)]
    fn reduce_word(&self, word: &[usize]) -> BTreeMap<u32, S> {
        for p in 0..word.len().saturating_sub(1) {
            let (a, c) = (word[p], word[p + 1]);
            if a == c {
                // e_a e_a = q[a]
                let q = self.q_val(a);
                let mut rest = Vec::with_capacity(word.len() - 2);
                rest.extend_from_slice(&word[..p]);
                rest.extend_from_slice(&word[p + 2..]);
                return scale(self.reduce_word(&rest), &q);
            } else if a > c {
                // e_a e_c = b[(c,a)] - e_c e_a
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
        // strictly increasing & distinct → a single canonical blade
        let mut mask = 0u32;
        for &g in word {
            mask |= 1 << g;
        }
        let mut m = BTreeMap::new();
        m.insert(mask, S::one());
        m
    }
}

fn scale<S: Scalar>(mut terms: BTreeMap<u32, S>, s: &S) -> BTreeMap<u32, S> {
    if s.is_zero() {
        return BTreeMap::new();
    }
    for v in terms.values_mut() {
        *v = v.mul(s);
    }
    terms.retain(|_, v| !v.is_zero());
    terms
}

fn merge<S: Scalar>(into: &mut BTreeMap<u32, S>, other: BTreeMap<u32, S>) {
    for (blade, coeff) in other {
        let e = into.entry(blade).or_insert_with(S::zero);
        *e = e.add(&coeff);
        if e.is_zero() {
            into.remove(&blade);
        }
    }
}

/// A multivector: blade-mask → coefficient (zeros never stored).
#[derive(Clone, Debug, PartialEq)]
pub struct Multivector<S: Scalar> {
    pub terms: BTreeMap<u32, S>,
}

/// A Clifford algebra: dimension + metric. Produces and combines multivectors.
#[derive(Clone, Debug)]
pub struct CliffordAlgebra<S: Scalar> {
    pub dim: usize,
    pub metric: Metric<S>,
}

impl<S: Scalar> CliffordAlgebra<S> {
    pub fn new(dim: usize, metric: Metric<S>) -> Self {
        CliffordAlgebra { dim, metric }
    }

    /// The graded (super) tensor product Cl(self) ⊗̂ Cl(other) ≅ Cl(self ⟂ other):
    /// generators of `self` keep their indices, generators of `other` are shifted
    /// up by `self.dim`, and the two blocks are mutually orthogonal. This is the
    /// composable form of "orthogonal sum of quadratic forms" — the operation the
    /// Arf invariant is additive over (`arf(A ⊗̂ B) = arf(A) + arf(B)`).
    pub fn graded_tensor(&self, other: &CliffordAlgebra<S>) -> CliffordAlgebra<S> {
        CliffordAlgebra::new(self.dim + other.dim, self.metric.direct_sum(&other.metric))
    }

    /// Embed a multivector of the *first* factor into `self ⊗̂ other`: the left
    /// block's generator indices are unchanged, so the blades carry over as-is.
    pub fn embed_first(&self, v: &Multivector<S>) -> Multivector<S> {
        Multivector {
            terms: v.terms.clone(),
        }
    }

    /// Embed a multivector of the *second* factor into `first ⊗̂ self`: shift every
    /// blade's generator bits up by `shift` (= the first factor's dimension).
    pub fn embed_second(&self, v: &Multivector<S>, shift: usize) -> Multivector<S> {
        let terms = v
            .terms
            .iter()
            .map(|(&blade, c)| (blade << shift, c.clone()))
            .collect();
        Multivector { terms }
    }

    pub fn zero(&self) -> Multivector<S> {
        Multivector {
            terms: BTreeMap::new(),
        }
    }

    pub fn scalar(&self, s: S) -> Multivector<S> {
        let mut terms = BTreeMap::new();
        if !s.is_zero() {
            terms.insert(0u32, s);
        }
        Multivector { terms }
    }

    /// The basis vector e_i.
    pub fn gen(&self, i: usize) -> Multivector<S> {
        let mut terms = BTreeMap::new();
        terms.insert(1u32 << i, S::one());
        Multivector { terms }
    }

    /// A single basis blade from a set of generators, coefficient 1.
    pub fn blade(&self, gens: &[usize]) -> Multivector<S> {
        let mut mask = 0u32;
        for &g in gens {
            mask |= 1 << g;
        }
        let mut terms = BTreeMap::new();
        terms.insert(mask, S::one());
        Multivector { terms }
    }

    pub fn add(&self, a: &Multivector<S>, b: &Multivector<S>) -> Multivector<S> {
        let mut terms = a.terms.clone();
        merge(&mut terms, b.terms.clone());
        Multivector { terms }
    }

    pub fn scalar_mul(&self, s: &S, a: &Multivector<S>) -> Multivector<S> {
        Multivector {
            terms: scale(a.terms.clone(), s),
        }
    }

    /// Geometric (Clifford) product. Computed in the wedge basis via the
    /// general-bilinear-form (Chevalley) product, so it is correct for an
    /// arbitrary metric `(q, b, a)` — ordinary Clifford when `a` is empty, and
    /// the deformed quantum-Clifford / Weyl-interpolating product when it isn't.
    pub fn mul(&self, a: &Multivector<S>, b: &Multivector<S>) -> Multivector<S> {
        let mut out: BTreeMap<u32, S> = BTreeMap::new();
        for (&ba, ca) in &a.terms {
            for (&bb, cb) in &b.terms {
                let reduced = self.metric.geom_product_blades(ba, bb);
                let coeff = ca.mul(cb);
                merge(&mut out, scale(reduced, &coeff));
            }
        }
        Multivector { terms: out }
    }

    /// Exterior (wedge) product — metric-independent.
    pub fn wedge(&self, a: &Multivector<S>, b: &Multivector<S>) -> Multivector<S> {
        let mut out: BTreeMap<u32, S> = BTreeMap::new();
        for (&ba, ca) in &a.terms {
            for (&bb, cb) in &b.terms {
                if ba & bb != 0 {
                    continue; // shared generator ⇒ wedge is 0
                }
                let sign = wedge_sign::<S>(ba, bb);
                let coeff = ca.mul(cb).mul(&sign);
                if coeff.is_zero() {
                    continue;
                }
                let e = out.entry(ba | bb).or_insert_with(S::zero);
                *e = e.add(&coeff);
                if e.is_zero() {
                    out.remove(&(ba | bb));
                }
            }
        }
        Multivector { terms: out }
    }

    /// Reversion: reverse the order of generators in every blade.
    /// On a grade-k blade this is (-1)^{k(k-1)/2}.
    pub fn reverse(&self, a: &Multivector<S>) -> Multivector<S> {
        let mut terms = BTreeMap::new();
        for (&blade, coeff) in &a.terms {
            let k = grade(blade);
            let c = if (k * (k.wrapping_sub(1)) / 2) & 1 == 1 {
                coeff.neg()
            } else {
                coeff.clone()
            };
            if !c.is_zero() {
                terms.insert(blade, c);
            }
        }
        Multivector { terms }
    }

    /// Grade-k projection.
    pub fn grade_part(&self, a: &Multivector<S>, k: u32) -> Multivector<S> {
        let terms = a
            .terms
            .iter()
            .filter(|(&blade, _)| grade(blade) == k)
            .map(|(&blade, c)| (blade, c.clone()))
            .collect();
        Multivector { terms }
    }

    /// The grade-0 (scalar) coefficient.
    pub fn scalar_part(&self, v: &Multivector<S>) -> S {
        v.terms.get(&0).cloned().unwrap_or_else(S::zero)
    }
}

impl<S: Scalar> Multivector<S> {
    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    /// Human-readable form, e.g. `3 + 2*e0 + 1*e0e1`.
    pub fn display(&self) -> String {
        if self.terms.is_empty() {
            return "0".to_string();
        }
        let one = S::one();
        let neg_one = S::one().neg();
        let mut parts = Vec::new();
        for (&blade, coeff) in &self.terms {
            if blade == 0 {
                parts.push(format!("{:?}", coeff));
                continue;
            }
            let label: String = bits(blade).iter().map(|i| format!("e{}", i)).collect();
            if *coeff == one {
                parts.push(label);
            } else if *coeff == neg_one {
                parts.push(format!("-{}", label));
            } else {
                parts.push(format!("{:?}*{}", coeff, label));
            }
        }
        parts.join(" + ")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scalar::Nimber;
    use crate::scalar::Rational;
    use crate::scalar::Surreal;

    fn r(n: i128) -> Rational {
        Rational::int(n)
    }

    #[test]
    fn complex_numbers_cl01() {
        // Cl(0,1): one generator with e0^2 = -1, the complex numbers.
        let alg = CliffordAlgebra::new(1, Metric::diagonal(vec![r(-1)]));
        let e0 = alg.gen(0);
        let sq = alg.mul(&e0, &e0);
        assert_eq!(sq, alg.scalar(r(-1)));
    }

    #[test]
    fn cl20_bivector_squares_to_minus_one() {
        // Cl(2,0): e0^2 = e1^2 = 1; e0e1 anticommutes and (e0e1)^2 = -1.
        let alg = CliffordAlgebra::new(2, Metric::diagonal(vec![r(1), r(1)]));
        let e0 = alg.gen(0);
        let e1 = alg.gen(1);
        // anticommute: e0e1 = -(e1e0)
        let e0e1 = alg.mul(&e0, &e1);
        let e1e0 = alg.mul(&e1, &e0);
        assert_eq!(e0e1, alg.scalar_mul(&r(-1), &e1e0));
        // (e0e1)^2 = -1
        let sq = alg.mul(&e0e1, &e0e1);
        assert_eq!(sq, alg.scalar(r(-1)));
    }

    #[test]
    fn grassmann_generators_are_nilpotent() {
        // q = 0 ⇒ e_i^2 = 0, and the wedge matches the product off-diagonal.
        let alg = CliffordAlgebra::new(3, Metric::grassmann(3));
        for i in 0..3 {
            let ei = alg.gen(i);
            assert!(alg.mul(&ei, &ei).is_zero(), "e{i}^2 should be 0");
        }
        let (e0, e1) = (alg.gen(0), alg.gen(1));
        assert_eq!(alg.mul(&e0, &e1), alg.wedge(&e0, &e1));
        // antisymmetry
        assert_eq!(
            alg.mul(&e0, &e1),
            alg.scalar_mul(&r(-1), &alg.mul(&e1, &e0))
        );
    }

    #[test]
    fn nimber_orthogonal_is_commutative() {
        // char 2, b = 0 ⇒ e_i e_j = e_j e_i (the genuine char-2-orthogonal fact).
        let alg = CliffordAlgebra::new(2, Metric::diagonal(vec![Nimber(2), Nimber(3)]));
        let e0 = alg.gen(0);
        let e1 = alg.gen(1);
        assert_eq!(alg.mul(&e0, &e1), alg.mul(&e1, &e0)); // commute
                                                          // e0^2 = q0 = 2 (a nimber!), not ±1
        assert_eq!(alg.mul(&e0, &e0), alg.scalar(Nimber(2)));
    }

    #[test]
    fn nimber_offdiagonal_is_noncommutative() {
        // char 2 with b[(0,1)] = t ⇒ e0 e1 + e1 e0 = t ≠ 0 ⇒ non-commutative.
        let mut b = BTreeMap::new();
        b.insert((0usize, 1usize), Nimber(1));
        let alg = CliffordAlgebra::new(2, Metric::new(vec![Nimber(0), Nimber(0)], b));
        let e0 = alg.gen(0);
        let e1 = alg.gen(1);
        let anti = alg.add(&alg.mul(&e0, &e1), &alg.mul(&e1, &e0));
        assert_eq!(anti, alg.scalar(Nimber(1))); // {e0,e1} = 1
        assert_ne!(alg.mul(&e0, &e1), alg.mul(&e1, &e0)); // not commutative
    }

    // The real stress test of reduce_word: associativity on a nontrivial,
    // non-orthogonal metric, in both characteristics.
    fn assert_associative<S: Scalar>(alg: &CliffordAlgebra<S>, gens: &[Multivector<S>]) {
        for a in gens {
            for b in gens {
                for c in gens {
                    let l = alg.mul(&alg.mul(a, b), c);
                    let r = alg.mul(a, &alg.mul(b, c));
                    assert_eq!(l, r, "associativity failed");
                }
            }
        }
    }

    #[test]
    fn associativity_rational_nonorthogonal() {
        let mut b = BTreeMap::new();
        b.insert((0usize, 1usize), r(1)); // non-orthogonal
        b.insert((1usize, 2usize), r(-1));
        let alg = CliffordAlgebra::new(3, Metric::new(vec![r(1), r(-1), r(2)], b));
        let gens = [
            alg.gen(0),
            alg.gen(1),
            alg.gen(2),
            alg.mul(&alg.gen(0), &alg.gen(1)),
            alg.add(&alg.gen(0), &alg.scalar(r(3))),
        ];
        assert_associative(&alg, &gens);
    }

    #[test]
    fn vector_inverse() {
        // Cl(3,0): unit vector squares to 1, so v⁻¹ = v.
        let alg = CliffordAlgebra::new(3, Metric::diagonal(vec![r(1), r(1), r(1)]));
        let v = alg.gen(0);
        let vi = alg.versor_inverse(&v).unwrap();
        assert_eq!(alg.mul(&v, &vi), alg.scalar(r(1)));
        assert_eq!(vi, v);
        // q=2: e0⁻¹ = e0/2
        let alg2 = CliffordAlgebra::new(1, Metric::diagonal(vec![r(2)]));
        let e0 = alg2.gen(0);
        assert_eq!(
            alg2.mul(&e0, &alg2.versor_inverse(&e0).unwrap()),
            alg2.scalar(r(1))
        );
        // a null vector has no inverse
        let alg0 = CliffordAlgebra::new(1, Metric::<Rational>::grassmann(1));
        assert!(alg0.versor_inverse(&alg0.gen(0)).is_none());
    }

    #[test]
    fn reflection_fixes_and_negates() {
        // reflect in the hyperplane ⊥ e1: fixes e0, negates e1.
        let alg = CliffordAlgebra::new(2, Metric::diagonal(vec![r(1), r(1)]));
        let (e0, e1) = (alg.gen(0), alg.gen(1));
        assert_eq!(alg.reflect(&e1, &e0).unwrap(), e0);
        assert_eq!(alg.reflect(&e1, &e1).unwrap(), alg.scalar_mul(&r(-1), &e1));
    }

    #[test]
    fn rotor_preserves_norm() {
        // A rotor (product of two unit vectors) is norm-preserving.
        let alg = CliffordAlgebra::new(2, Metric::diagonal(vec![r(1), r(1)]));
        let rotor = alg.mul(&alg.gen(0), &alg.gen(1));
        let x = alg.add(
            &alg.scalar_mul(&r(3), &alg.gen(0)),
            &alg.scalar_mul(&r(4), &alg.gen(1)),
        );
        let rx = alg.sandwich(&rotor, &x).unwrap();
        assert_eq!(alg.norm2(&rx), alg.norm2(&x)); // 25 either way
    }

    #[test]
    fn twisted_adjoint_matches_reflect_and_sandwich() {
        let alg = CliffordAlgebra::new(2, Metric::diagonal(vec![r(1), r(1)]));
        let (e0, e1) = (alg.gen(0), alg.gen(1));
        let x = alg.add(&alg.scalar_mul(&r(3), &e0), &alg.scalar_mul(&r(4), &e1));
        // Odd versor (a vector): twisted adjoint == reflection in its ⊥ hyperplane.
        assert_eq!(
            alg.twisted_sandwich(&e1, &x).unwrap(),
            alg.reflect(&e1, &x).unwrap()
        );
        // Even versor (a rotor): α(v)=v, so twisted == untwisted sandwich.
        let rotor = alg.mul(&e0, &e1);
        assert_eq!(
            alg.twisted_sandwich(&rotor, &x).unwrap(),
            alg.sandwich(&rotor, &x).unwrap()
        );
    }

    #[test]
    fn left_contraction_lowers_grade() {
        let alg = CliffordAlgebra::new(3, Metric::diagonal(vec![r(1), r(1), r(1)]));
        let e0 = alg.gen(0);
        let e0e1 = alg.mul(&alg.gen(0), &alg.gen(1));
        assert_eq!(alg.left_contract(&e0, &e0e1), alg.gen(1)); // e0 ⌟ (e0∧e1) = e1
        let three = alg.scalar(r(3));
        assert_eq!(
            alg.left_contract(&three, &e0e1),
            alg.scalar_mul(&r(3), &e0e1)
        );
    }

    #[test]
    fn dual_of_vector_is_bivector_in_3d() {
        let alg = CliffordAlgebra::new(3, Metric::diagonal(vec![r(1), r(1), r(1)]));
        let d = alg.dual(&alg.gen(0)).unwrap();
        assert!(!d.is_zero());
        assert_eq!(alg.grade_part(&d, 2), d); // purely grade 2
    }

    #[test]
    fn grade_involution_signs() {
        let alg = CliffordAlgebra::new(2, Metric::diagonal(vec![r(1), r(1)]));
        // 5 + e0 + e0e1  ↦  5 − e0 + e0e1
        let v = alg.add(
            &alg.scalar(r(5)),
            &alg.add(&alg.gen(0), &alg.mul(&alg.gen(0), &alg.gen(1))),
        );
        let expect = alg.add(
            &alg.scalar(r(5)),
            &alg.add(
                &alg.scalar_mul(&r(-1), &alg.gen(0)),
                &alg.mul(&alg.gen(0), &alg.gen(1)),
            ),
        );
        assert_eq!(alg.grade_involution(&v), expect);
    }

    #[test]
    fn versor_over_surreal_metric() {
        // e0² = ω (a monomial ⇒ invertible). e0⁻¹ = ε·e0, and sandwich works.
        let alg = CliffordAlgebra::new(
            2,
            Metric::diagonal(vec![Surreal::omega(), Surreal::epsilon()]),
        );
        let e0 = alg.gen(0);
        let inv = alg.versor_inverse(&e0).unwrap();
        assert_eq!(alg.mul(&e0, &inv), alg.scalar(Surreal::one()));
    }

    #[test]
    fn even_subalgebra_of_cl30_is_quaternions() {
        // Cl(3,0)⁰ ≅ Cl(0,2) ≅ ℍ: two generators squaring to −1 that anticommute.
        let alg = CliffordAlgebra::new(3, Metric::diagonal(vec![r(1), r(1), r(1)]));
        let even = alg.even_subalgebra().unwrap();
        assert_eq!(even.dim, 2);
        let (f0, f1) = (even.gen(0), even.gen(1));
        assert_eq!(even.mul(&f0, &f0), even.scalar(r(-1))); // i² = −1
        assert_eq!(even.mul(&f1, &f1), even.scalar(r(-1))); // j² = −1
        assert_eq!(
            even.mul(&f0, &f1),
            even.scalar_mul(&r(-1), &even.mul(&f1, &f0))
        ); // ij = −ji
    }

    #[test]
    fn even_part_projection() {
        let alg = CliffordAlgebra::new(2, Metric::diagonal(vec![r(1), r(1)]));
        // 5 + 2 e0 + 3 e0e1  ↦  5 + 3 e0e1
        let v = alg.add(
            &alg.scalar(r(5)),
            &alg.add(
                &alg.scalar_mul(&r(2), &alg.gen(0)),
                &alg.mul(&alg.gen(0), &alg.gen(1)),
            ),
        );
        let expect = alg.add(&alg.scalar(r(5)), &alg.mul(&alg.gen(0), &alg.gen(1)));
        assert_eq!(alg.even_part(&v), expect);
    }

    #[test]
    fn graded_tensor_blocks_are_orthogonal() {
        // Cl(1,0) ⊗̂ Cl(0,1): e0²=+1 (left), e1²=−1 (right), and the two blocks
        // anticommute (cross polar form 0 ⇒ in char 0, e0 e1 = −e1 e0).
        let left = CliffordAlgebra::new(1, Metric::diagonal(vec![r(1)]));
        let right = CliffordAlgebra::new(1, Metric::diagonal(vec![r(-1)]));
        let alg = left.graded_tensor(&right);
        assert_eq!(alg.dim, 2);
        let e0 = alg.gen(0); // from the left factor
        let e1 = alg.gen(1); // from the right factor (shifted)
        assert_eq!(alg.mul(&e0, &e0), alg.scalar(r(1)));
        assert_eq!(alg.mul(&e1, &e1), alg.scalar(r(-1)));
        // cross generators anticommute
        assert_eq!(
            alg.mul(&e0, &e1),
            alg.scalar_mul(&r(-1), &alg.mul(&e1, &e0))
        );
        // the embeddings agree with native generators
        assert_eq!(alg.embed_first(&left.gen(0)), e0);
        assert_eq!(alg.embed_second(&right.gen(0), left.dim), e1);
    }

    #[test]
    fn general_product_reproduces_reduce_word_when_a_empty() {
        // The Chevalley general-B product must agree blade-for-blade with the
        // proven `reduce_word` on every metric whose upper part `a` is empty.
        let mut b = BTreeMap::new();
        b.insert((0usize, 1usize), r(1));
        b.insert((1usize, 2usize), r(-1));
        b.insert((0usize, 2usize), r(2));
        let m = Metric::new(vec![r(1), r(-1), r(2)], b);
        for ba in 0u32..8 {
            for bb in 0u32..8 {
                let word: Vec<usize> = bits(ba).into_iter().chain(bits(bb)).collect();
                assert_eq!(
                    m.geom_product_blades(ba, bb),
                    m.reduce_word(&word),
                    "mismatch on blades {ba:#b}·{bb:#b}"
                );
            }
        }
    }

    #[test]
    fn general_bilinear_in_order_contraction() {
        // q=[1,1], b empty, a[(0,1)] = 5: the in-order product gains the scalar,
        //   e0 e1 = e0∧e1 + 5,   but the anticommutator {e0,e1} = b = 0 is unchanged.
        let mut a = BTreeMap::new();
        a.insert((0usize, 1usize), r(5));
        let alg = CliffordAlgebra::new(2, Metric::general(vec![r(1), r(1)], BTreeMap::new(), a));
        let (e0, e1) = (alg.gen(0), alg.gen(1));
        let blade = alg.wedge(&e0, &e1);
        assert_eq!(alg.mul(&e0, &e1), alg.add(&blade, &alg.scalar(r(5)))); // e0∧e1 + 5
                                                                           // anticommutator depends only on b (here 0), not on a
        assert_eq!(alg.add(&alg.mul(&e0, &e1), &alg.mul(&e1, &e0)), alg.zero());
    }

    #[test]
    fn associativity_general_bilinear_form() {
        // A genuinely general bilinear form: nonzero polar b AND nonzero in-order
        // a, in both characteristics. The Chevalley product must stay associative.
        let mut b = BTreeMap::new();
        b.insert((0usize, 1usize), r(1));
        b.insert((1usize, 2usize), r(2));
        let mut a = BTreeMap::new();
        a.insert((0usize, 1usize), r(3));
        a.insert((0usize, 2usize), r(-1));
        let alg = CliffordAlgebra::new(3, Metric::general(vec![r(2), r(-1), r(1)], b, a));
        let gens = [
            alg.gen(0),
            alg.gen(1),
            alg.gen(2),
            alg.mul(&alg.gen(0), &alg.gen(1)),
            alg.add(&alg.gen(2), &alg.scalar(r(3))),
        ];
        assert_associative(&alg, &gens);

        // char 2 version
        let mut bn = BTreeMap::new();
        bn.insert((0usize, 1usize), Nimber(1));
        let mut an = BTreeMap::new();
        an.insert((0usize, 1usize), Nimber(2));
        an.insert((1usize, 2usize), Nimber(3));
        let algn = CliffordAlgebra::new(
            3,
            Metric::general(vec![Nimber(2), Nimber(1), Nimber(0)], bn, an),
        );
        let gensn = [
            algn.gen(0),
            algn.gen(1),
            algn.gen(2),
            algn.mul(&algn.gen(0), &algn.gen(1)),
        ];
        assert_associative(&algn, &gensn);
    }

    #[test]
    fn associativity_nimber_nonorthogonal() {
        let mut b = BTreeMap::new();
        b.insert((0usize, 1usize), Nimber(1));
        b.insert((0usize, 2usize), Nimber(3));
        let alg = CliffordAlgebra::new(3, Metric::new(vec![Nimber(2), Nimber(1), Nimber(0)], b));
        let gens = [
            alg.gen(0),
            alg.gen(1),
            alg.gen(2),
            alg.mul(&alg.gen(0), &alg.gen(1)),
            alg.add(&alg.gen(2), &alg.scalar(Nimber(5))),
        ];
        assert_associative(&alg, &gens);
    }
}
