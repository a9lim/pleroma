//! Witt classes of finite quadratic modules.
//!
//! This is the Wall/Nikulin finite-module side of the integral pillar: a
//! nonsingular finite quadratic module is reduced, prime by prime, by quotienting
//! isotropic cyclic subgroups until an anisotropic core remains. The terminal core
//! is canonicalised as a finite table, so the output is a full Witt normal form,
//! not just the Milgram/Brown phase.

use crate::forms::integral::discriminant::{phase_mod8_from_q_values, DiscriminantForm, IsoTables};
use crate::forms::padic::try_is_square_qp;
use crate::scalar::{Rational, Scalar};
use std::collections::{BTreeMap, BTreeSet, VecDeque};

const FQM_WITT_GROUP_CAP: usize = 512;
const FQM_WITT_TUPLE_CAP: u128 = 2_000_000;

/// A value-count entry in a finite quadratic module normal form.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FqmValueCount {
    /// Numerator of the canonical rational representative.
    pub numer: i128,
    /// Denominator of the canonical rational representative.
    pub denom: i128,
    /// Number of elements carrying this value.
    pub count: u128,
}

/// One p-primary summand of the finite-quadratic-module Witt class.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FqmPrimaryWittClass {
    /// The prime `p`.
    pub prime: u128,
    /// The order of the original p-primary summand.
    pub order: u128,
    /// The order of the anisotropic Witt core.
    pub core_order: u128,
    /// Invariant factors of the anisotropic core.
    pub core_group: Vec<u128>,
    /// Exponent of the anisotropic core.
    pub core_exponent: u128,
    /// The Milgram/Brown phase of this p-primary Witt class.
    pub phase_mod8: i128,
    /// Value counts on the anisotropic core, useful as readable diagnostics.
    pub q_value_counts: Vec<FqmValueCount>,
    /// Opaque exact normal form of the anisotropic core.
    ///
    /// The label is canonical under finite-module isomorphism; equality of these
    /// labels is the equality test for the p-primary Witt class.
    pub normal_form: Vec<i128>,
}

/// The Witt class of a nonsingular finite quadratic module.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FqmWittClass {
    /// Order of the original finite module.
    pub order: u128,
    /// Total Milgram/Brown phase, i.e. the sum of the p-primary phases in `Z/8`.
    pub phase_mod8: i128,
    /// Prime-local Witt normal forms.
    pub primary: Vec<FqmPrimaryWittClass>,
}

impl FqmWittClass {
    /// Whether the class is Witt-trivial, i.e. every p-primary anisotropic core is
    /// the zero module.
    pub fn is_trivial(&self) -> bool {
        self.primary.iter().all(|p| p.core_order == 1)
    }
}

/// A local condition in Nikulin's even-lattice existence criterion.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NikulinPrimaryExistenceInvariants {
    /// The prime `p`.
    pub prime: u128,
    /// The order of the p-primary summand.
    pub order: u128,
    /// The minimal number of generators `l(A_p)`.
    pub length: usize,
    /// Whether the requested rank is exactly `l(A_p)`, so Nikulin's determinant
    /// side condition is active at this prime.
    pub equality_case: bool,
    /// For `p = 2`, whether the 2-primary quadratic form is even in Nikulin's
    /// sense, i.e. it has no order-2 cyclic summand with q-value odd/2.
    pub even_two_primary: bool,
    /// The p-adic determinant square class `discr K(q_p)` represented by an exact
    /// rational. Present only in equality cases where a determinant check is
    /// required.
    pub p_adic_discriminant: Option<Rational>,
    /// Result of the equality-case determinant check, when one is required.
    pub determinant_condition_holds: Option<bool>,
}

/// The first failed condition in Nikulin's theorem 1.10.1.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum NikulinExistenceObstruction {
    /// `sign(q) != t_+ - t_- (mod 8)`.
    SignatureCongruence {
        required_mod8: i128,
        module_phase_mod8: i128,
    },
    /// `rank < l(A_p)` at one prime.
    RankTooSmall {
        prime: u128,
        rank: usize,
        length: usize,
    },
    /// The odd-prime equality case failed:
    /// `(-1)^{t_-}|A_p| != discr K(q_p)` in `Q_p^*/Q_p^{*2}`.
    OddPrimeDeterminant {
        prime: u128,
        signed_order: i128,
        p_adic_discriminant: Rational,
    },
    /// The 2-adic even equality case failed:
    /// `|A_2| != +/- discr K(q_2)` in `Q_2^*/Q_2^{*2}`.
    TwoAdicDeterminant {
        order: u128,
        p_adic_discriminant: Rational,
    },
}

/// Full bounded report for Nikulin's even-lattice existence criterion.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NikulinExistenceInvariants {
    /// Requested signature `(t_+, t_-)`.
    pub signature: (usize, usize),
    /// Requested rank `t_+ + t_-`.
    pub rank: usize,
    /// The finite quadratic module's Gauss/Milgram phase, `sign(q) mod 8`.
    pub module_phase_mod8: i128,
    /// Prime-local rank and determinant checks.
    pub primary: Vec<NikulinPrimaryExistenceInvariants>,
    /// The first failed condition, or `None` when the lattice exists.
    pub obstruction: Option<NikulinExistenceObstruction>,
}

impl NikulinExistenceInvariants {
    /// Whether Nikulin's theorem decides that an even lattice with the requested
    /// signature and discriminant form exists.
    pub fn exists(&self) -> bool {
        self.obstruction.is_none()
    }
}

/// A native finite quadratic module in a cyclic product presentation.
///
/// The `q_values_mod2` slice is ordered lexicographically over the cyclic factors:
/// for factors `[d0, d1, ...]`, index `((x0*d1 + x1)*d2 + ...)` stores
/// `q(x0, x1, ...)` as a rational in `Q/2Z`. The constructor validates
/// nonsingularity and the quadratic law up to `FQM_WITT_GROUP_CAP`.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteQuadraticModule {
    cyclic_factors: Vec<u128>,
    q_values_mod2: Vec<Rational>,
}

impl FiniteQuadraticModule {
    /// Build a nonsingular finite quadratic module from a cyclic presentation and
    /// all of its quadratic values in lexicographic coordinate order.
    pub fn new(cyclic_factors: Vec<u128>, q_values_mod2: Vec<Rational>) -> Option<Self> {
        if cyclic_factors.iter().any(|&d| d <= 1) {
            return None;
        }
        let order = cyclic_factors
            .iter()
            .try_fold(1usize, |acc, &d| acc.checked_mul(usize::try_from(d).ok()?))?;
        if order == 0 || order > FQM_WITT_GROUP_CAP || q_values_mod2.len() != order {
            return None;
        }
        let q_values_mod2 = q_values_mod2
            .into_iter()
            .map(|q| rational_mod_int(q, 2))
            .collect::<Vec<_>>();
        let module = FiniteQuadraticModule {
            cyclic_factors,
            q_values_mod2,
        };
        let table = FqmTable::from_native(&module)?;
        if table.q[table.zero] != Rational::zero()
            || !table.quadratic_values_are_even()
            || !table.bilinear_form_is_biadditive()
            || !table.is_nondegenerate()
        {
            return None;
        }
        Some(module)
    }

    /// The cyclic module generated by `g` with `q(g) = generator_q`.
    pub fn cyclic(order: u128, generator_q: Rational) -> Option<Self> {
        if order <= 1 || usize::try_from(order).ok()? > FQM_WITT_GROUP_CAP {
            return None;
        }
        let qg = rational_mod_int(generator_q, 2);
        let mut q_values = Vec::with_capacity(usize::try_from(order).ok()?);
        for k in 0..order {
            let kk = i128::try_from(k.checked_mul(k)?).ok()?;
            q_values.push(rational_mod_int(Rational::from_int(kk).mul(&qg), 2));
        }
        Self::new(vec![order], q_values)
    }

    /// Orthogonal direct sum.
    pub fn direct_sum(&self, other: &Self) -> Option<Self> {
        let left = FqmTable::from_native(self)?;
        let right = FqmTable::from_native(other)?;
        let mut factors = self.cyclic_factors.clone();
        factors.extend(other.cyclic_factors.iter().copied());
        let mut q_values = Vec::with_capacity(left.q.len().checked_mul(right.q.len())?);
        for ql in &left.q {
            for qr in &right.q {
                q_values.push(rational_mod_int(ql.add(qr), 2));
            }
        }
        Self::new(factors, q_values)
    }

    /// Order of the finite module.
    pub fn order(&self) -> u128 {
        self.q_values_mod2.len() as u128
    }

    /// Cyclic factors of this presentation.
    pub fn cyclic_factors(&self) -> &[u128] {
        &self.cyclic_factors
    }

    /// Quadratic values in lexicographic coordinate order.
    pub fn q_values_mod2(&self) -> &[Rational] {
        &self.q_values_mod2
    }

    /// The Wall/Nikulin Witt normal form.
    pub fn witt_class(&self) -> Option<FqmWittClass> {
        FqmTable::from_native(self)?.witt_class()
    }

    /// Nikulin's even-lattice existence criterion for this finite quadratic
    /// module and the requested signature `(t_+, t_-)`.
    ///
    /// This implements Nikulin, *Integral symmetric bilinear forms and some of
    /// their applications*, Math. USSR Izv. **14** (1980), Theorem 1.10.1, in the
    /// bounded finite-table model used by [`witt_class`](Self::witt_class).
    /// `None` means the table/determinant computation exceeded that bounded exact
    /// surface, not that the theorem failed.
    pub fn nikulin_existence_report(
        &self,
        signature: (usize, usize),
    ) -> Option<NikulinExistenceInvariants> {
        FqmTable::from_native(self)?.nikulin_existence_report(signature)
    }

    /// Boolean convenience wrapper around [`nikulin_existence_report`](Self::nikulin_existence_report).
    pub fn nikulin_even_lattice_exists(&self, signature: (usize, usize)) -> Option<bool> {
        Some(self.nikulin_existence_report(signature)?.exists())
    }
}

impl DiscriminantForm {
    /// The full Wall/Nikulin finite-quadratic-module Witt class of `(A_L, q_L)`.
    ///
    /// This refines [`fqm_gauss_phase`](Self::fqm_gauss_phase): the phase is kept as
    /// a projection, but equality is decided by the p-primary anisotropic normal
    /// forms. The implementation is exact up to the finite enumeration budget; it
    /// returns `None` instead of truncating when `|A_L| > 512`.
    pub fn fqm_witt_class(&self) -> Option<FqmWittClass> {
        FqmTable::from_iso_tables(self.tables_bounded(FQM_WITT_GROUP_CAP)?).witt_class()
    }

    /// Whether two discriminant forms are Witt-equivalent as finite quadratic
    /// modules.
    pub fn is_fqm_witt_equivalent(&self, other: &Self) -> Option<bool> {
        Some(self.fqm_witt_class()? == other.fqm_witt_class()?)
    }

    /// Nikulin's even-lattice existence criterion for this discriminant form and
    /// the requested signature `(t_+, t_-)`.
    ///
    /// This is the existence companion to [`is_isomorphic`](Self::is_isomorphic):
    /// instead of comparing two already-built lattices, it decides whether the
    /// pair `(signature, q)` is realized by some even lattice. The implementation
    /// follows Nikulin theorem 1.10.1 and returns `None` only past the bounded
    /// finite-table surface (`|A| <= 512` here).
    pub fn nikulin_existence_report(
        &self,
        signature: (usize, usize),
    ) -> Option<NikulinExistenceInvariants> {
        FqmTable::from_iso_tables(self.tables_bounded(FQM_WITT_GROUP_CAP)?)
            .nikulin_existence_report(signature)
    }

    /// Boolean convenience wrapper around [`nikulin_existence_report`](Self::nikulin_existence_report).
    pub fn nikulin_even_lattice_exists(&self, signature: (usize, usize)) -> Option<bool> {
        Some(self.nikulin_existence_report(signature)?.exists())
    }
}

#[derive(Clone, Debug)]
struct FqmTable {
    zero: usize,
    q: Vec<Rational>,
    order: Vec<usize>,
    add: Vec<Vec<usize>>,
}

impl FqmTable {
    fn from_iso_tables(t: IsoTables) -> Self {
        FqmTable {
            zero: t.zero,
            q: t.q,
            order: t.order,
            add: t.add,
        }
    }

    fn from_native(module: &FiniteQuadraticModule) -> Option<Self> {
        let n = module.q_values_mod2.len();
        let mut add = vec![vec![0usize; n]; n];
        for i in 0..n {
            let ci = coords_from_index(i, &module.cyclic_factors)?;
            for j in 0..n {
                let cj = coords_from_index(j, &module.cyclic_factors)?;
                let sum = ci
                    .iter()
                    .zip(&cj)
                    .zip(&module.cyclic_factors)
                    .map(|((&a, &b), &d)| (a + b) % d)
                    .collect::<Vec<_>>();
                add[i][j] = index_from_coords(&sum, &module.cyclic_factors)?;
            }
        }
        let zero = 0;
        let mut out = FqmTable {
            zero,
            q: module.q_values_mod2.clone(),
            order: vec![1; n],
            add,
        };
        out.compute_orders();
        Some(out)
    }

    fn compute_orders(&mut self) {
        let n = self.q.len();
        self.order = vec![1usize; n];
        for i in 0..n {
            let mut cur = i;
            let mut k = 1usize;
            while cur != self.zero {
                cur = self.add[cur][i];
                k += 1;
            }
            self.order[i] = k;
        }
    }

    fn witt_class(&self) -> Option<FqmWittClass> {
        if self.q.len() > FQM_WITT_GROUP_CAP {
            return None;
        }
        let mut primes = BTreeSet::new();
        for &ord in &self.order {
            for p in prime_factors_u128(ord as u128) {
                primes.insert(p);
            }
        }
        let mut primary = Vec::new();
        for p in primes {
            let part = self.primary_subtable(p)?;
            let phase = phase_mod8_from_q_values(part.q.iter(), part.q.len())?;
            let mut memo = BTreeMap::new();
            let core = part.anisotropic_core(&mut memo)?;
            let core_phase = phase_mod8_from_q_values(core.q.iter(), core.q.len())?;
            if core_phase != phase {
                return None;
            }
            primary.push(FqmPrimaryWittClass {
                prime: p,
                order: part.q.len() as u128,
                core_order: core.q.len() as u128,
                core_group: core.primary_invariant_factors(p)?,
                core_exponent: core.order.iter().copied().max().unwrap_or(1) as u128,
                phase_mod8: phase,
                q_value_counts: core.q_value_counts(),
                normal_form: core.canonical_label()?,
            });
        }
        let phase_mod8 = primary
            .iter()
            .map(|p| p.phase_mod8)
            .sum::<i128>()
            .rem_euclid(8);
        Some(FqmWittClass {
            order: self.q.len() as u128,
            phase_mod8,
            primary,
        })
    }

    fn nikulin_existence_report(
        &self,
        signature: (usize, usize),
    ) -> Option<NikulinExistenceInvariants> {
        if self.q.len() > FQM_WITT_GROUP_CAP {
            return None;
        }
        let rank = signature.0.checked_add(signature.1)?;
        let sig_plus = i128::try_from(signature.0).ok()?;
        let sig_minus = i128::try_from(signature.1).ok()?;
        let required_mod8 = (sig_plus - sig_minus).rem_euclid(8);
        let module_phase_mod8 = phase_mod8_from_q_values(self.q.iter(), self.q.len())?;
        let mut obstruction = (required_mod8 != module_phase_mod8).then_some(
            NikulinExistenceObstruction::SignatureCongruence {
                required_mod8,
                module_phase_mod8,
            },
        );

        let mut primes = BTreeSet::new();
        for &ord in &self.order {
            for p in prime_factors_u128(ord as u128) {
                primes.insert(p);
            }
        }

        let mut primary = Vec::new();
        for p in primes {
            let part = self.primary_subtable(p)?;
            let length = part.direct_product_generators()?.len();
            let order = part.q.len() as u128;
            let equality_case = rank == length;
            let even_two_primary = p == 2 && !part.has_odd_two_adic_summand();
            let mut p_adic_discriminant = None;
            let mut determinant_condition_holds = None;

            if rank < length && obstruction.is_none() {
                obstruction = Some(NikulinExistenceObstruction::RankTooSmall {
                    prime: p,
                    rank,
                    length,
                });
            }

            if equality_case && p != 2 {
                let discr = part.p_adic_discriminant()?;
                let signed_order = signed_order_for_odd_prime(order, signature.1)?;
                let signed_order_q = Rational::from_int(signed_order);
                let ok = same_square_class_odd(&signed_order_q, &discr, p)?;
                if !ok && obstruction.is_none() {
                    obstruction = Some(NikulinExistenceObstruction::OddPrimeDeterminant {
                        prime: p,
                        signed_order,
                        p_adic_discriminant: discr.clone(),
                    });
                }
                p_adic_discriminant = Some(discr);
                determinant_condition_holds = Some(ok);
            } else if equality_case && even_two_primary {
                let discr = part.p_adic_discriminant()?;
                let order_q = rational_from_u128(order)?;
                let ok = same_square_class_2_up_to_sign(&order_q, &discr)?;
                if !ok && obstruction.is_none() {
                    obstruction = Some(NikulinExistenceObstruction::TwoAdicDeterminant {
                        order,
                        p_adic_discriminant: discr.clone(),
                    });
                }
                p_adic_discriminant = Some(discr);
                determinant_condition_holds = Some(ok);
            }

            primary.push(NikulinPrimaryExistenceInvariants {
                prime: p,
                order,
                length,
                equality_case,
                even_two_primary,
                p_adic_discriminant,
                determinant_condition_holds,
            });
        }

        Some(NikulinExistenceInvariants {
            signature,
            rank,
            module_phase_mod8,
            primary,
            obstruction,
        })
    }

    fn primary_subtable(&self, p: u128) -> Option<Self> {
        let indices = self
            .order
            .iter()
            .enumerate()
            .filter_map(|(i, &ord)| is_prime_power_order(ord as u128, p).then_some(i))
            .collect::<Vec<_>>();
        self.induced_subtable(&indices)
    }

    fn induced_subtable(&self, indices: &[usize]) -> Option<Self> {
        let mut map = vec![usize::MAX; self.q.len()];
        for (new, &old) in indices.iter().enumerate() {
            map[old] = new;
        }
        let zero = map[self.zero];
        if zero == usize::MAX {
            return None;
        }
        let mut add = vec![vec![0usize; indices.len()]; indices.len()];
        for (i, &old_i) in indices.iter().enumerate() {
            for (j, &old_j) in indices.iter().enumerate() {
                let s = self.add[old_i][old_j];
                let mapped = map[s];
                if mapped == usize::MAX {
                    return None;
                }
                add[i][j] = mapped;
            }
        }
        let mut out = FqmTable {
            zero,
            q: indices.iter().map(|&i| self.q[i].clone()).collect(),
            order: vec![1; indices.len()],
            add,
        };
        out.compute_orders();
        Some(out)
    }

    fn anisotropic_core(&self, memo: &mut BTreeMap<Vec<i128>, FqmTable>) -> Option<Self> {
        let raw = self.raw_label()?;
        if let Some(hit) = memo.get(&raw) {
            return Some(hit.clone());
        }
        let isotropic = (0..self.q.len())
            .filter(|&i| i != self.zero && self.q[i] == Rational::zero())
            .collect::<Vec<_>>();
        if isotropic.is_empty() {
            memo.insert(raw, self.clone());
            return Some(self.clone());
        }

        let mut best: Option<(Vec<i128>, FqmTable)> = None;
        for x in isotropic {
            let h = self.subgroup_generated(&[x]);
            if h.len() <= 1 {
                continue;
            }
            let quotient = self.quotient_by_isotropic_subgroup(&h)?;
            if quotient.q.len() >= self.q.len() {
                return None;
            }
            let core = quotient.anisotropic_core(memo)?;
            let label = core.canonical_label()?;
            if best.as_ref().is_none_or(|(b, _)| label < *b) {
                best = Some((label, core));
            }
        }
        let core = best?.1;
        memo.insert(raw, core.clone());
        Some(core)
    }

    fn quotient_by_isotropic_subgroup(&self, subgroup: &BTreeSet<usize>) -> Option<Self> {
        if !subgroup.contains(&self.zero)
            || !subgroup.iter().all(|&h| self.q[h] == Rational::zero())
        {
            return None;
        }
        let orthogonal = (0..self.q.len())
            .filter(|&x| {
                subgroup
                    .iter()
                    .all(|&h| self.bilinear_value(x, h) == Rational::zero())
            })
            .collect::<BTreeSet<_>>();
        if !subgroup.is_subset(&orthogonal) {
            return None;
        }

        let mut coset_of = vec![usize::MAX; self.q.len()];
        let mut reps = Vec::new();
        for &x in &orthogonal {
            if coset_of[x] != usize::MAX {
                continue;
            }
            let id = reps.len();
            reps.push(x);
            for &h in subgroup {
                let y = self.add[x][h];
                if !orthogonal.contains(&y) {
                    return None;
                }
                coset_of[y] = id;
            }
        }
        let zero = coset_of[self.zero];
        if zero == usize::MAX {
            return None;
        }
        let mut add = vec![vec![0usize; reps.len()]; reps.len()];
        for (i, &x) in reps.iter().enumerate() {
            for (j, &y) in reps.iter().enumerate() {
                let s = self.add[x][y];
                let mapped = coset_of[s];
                if mapped == usize::MAX {
                    return None;
                }
                add[i][j] = mapped;
            }
        }
        let mut out = FqmTable {
            zero,
            q: reps.iter().map(|&i| self.q[i].clone()).collect(),
            order: vec![1; reps.len()],
            add,
        };
        out.compute_orders();
        Some(out)
    }

    fn subgroup_generated(&self, gens: &[usize]) -> BTreeSet<usize> {
        let mut set = BTreeSet::new();
        let mut queue = VecDeque::new();
        set.insert(self.zero);
        queue.push_back(self.zero);
        while let Some(x) = queue.pop_front() {
            for &g in gens {
                let nx = self.add[x][g];
                if set.insert(nx) {
                    queue.push_back(nx);
                }
            }
        }
        set
    }

    fn generator_rank(&self) -> usize {
        let mut gens = Vec::new();
        let mut covered = self.subgroup_generated(&gens);
        while covered.len() < self.q.len() {
            let g = (0..self.q.len())
                .filter(|i| !covered.contains(i))
                .max_by_key(|&i| self.order[i])
                .expect("uncovered finite-module element exists");
            gens.push(g);
            covered = self.subgroup_generated(&gens);
        }
        gens.len()
    }

    fn direct_product_generators(&self) -> Option<Vec<usize>> {
        if self.q.len() == 1 {
            return Some(Vec::new());
        }
        let mut candidates = (0..self.q.len())
            .filter(|&i| i != self.zero)
            .collect::<Vec<_>>();
        candidates.sort_by(|&a, &b| self.order[b].cmp(&self.order[a]).then_with(|| a.cmp(&b)));
        let mut gens = Vec::new();
        let covered = self.subgroup_generated(&gens);
        self.direct_product_generators_rec(&candidates, &mut gens, covered)
    }

    fn direct_product_generators_rec(
        &self,
        candidates: &[usize],
        gens: &mut Vec<usize>,
        covered: BTreeSet<usize>,
    ) -> Option<Vec<usize>> {
        if covered.len() == self.q.len() {
            return Some(gens.clone());
        }
        for &g in candidates {
            if covered.contains(&g) || gens.contains(&g) {
                continue;
            }
            let mut trial = gens.clone();
            trial.push(g);
            let trial_covered = self.subgroup_generated(&trial);
            let expected = covered.len().checked_mul(self.order[g])?;
            if trial_covered.len() != expected {
                continue;
            }
            gens.push(g);
            if let Some(out) = self.direct_product_generators_rec(candidates, gens, trial_covered) {
                return Some(out);
            }
            gens.pop();
        }
        None
    }

    fn p_adic_discriminant(&self) -> Option<Rational> {
        let gens = self.direct_product_generators()?;
        if gens.is_empty() {
            return Some(Rational::one());
        }
        let mut matrix = vec![vec![Rational::zero(); gens.len()]; gens.len()];
        for (i, &x) in gens.iter().enumerate() {
            for (j, &y) in gens.iter().enumerate() {
                matrix[i][j] = self.bilinear_value(x, y);
            }
        }
        let det_pairing = rational_det(matrix)?;
        det_pairing.inv()
    }

    fn has_odd_two_adic_summand(&self) -> bool {
        (0..self.q.len()).any(|i| {
            self.order[i] == 2 && self.q[i].denom() == 2 && self.q[i].numer().rem_euclid(2) == 1
        })
    }

    fn canonical_label(&self) -> Option<Vec<i128>> {
        if self.q.len() == 1 {
            return Some(vec![1, 0]);
        }
        let rank = self.generator_rank();
        let candidates = (0..self.q.len())
            .filter(|&i| i != self.zero)
            .collect::<Vec<_>>();
        let mut tuple = Vec::with_capacity(rank);
        let mut best: Option<Vec<i128>> = None;
        let mut seen = 0u128;
        self.canonical_label_rec(rank, &candidates, &mut tuple, &mut best, &mut seen)?;
        best
    }

    fn canonical_label_rec(
        &self,
        rank: usize,
        candidates: &[usize],
        tuple: &mut Vec<usize>,
        best: &mut Option<Vec<i128>>,
        seen: &mut u128,
    ) -> Option<()> {
        if tuple.len() == rank {
            *seen = seen.checked_add(1)?;
            if *seen > FQM_WITT_TUPLE_CAP {
                return None;
            }
            if let Some(order) = self.ordered_elements_from_generators(tuple) {
                let label = self.label_for_order(&order)?;
                if best.as_ref().is_none_or(|b| label < *b) {
                    *best = Some(label);
                }
            }
            return Some(());
        }
        for &cand in candidates {
            if tuple.contains(&cand) {
                continue;
            }
            tuple.push(cand);
            self.canonical_label_rec(rank, candidates, tuple, best, seen)?;
            tuple.pop();
        }
        Some(())
    }

    fn ordered_elements_from_generators(&self, gens: &[usize]) -> Option<Vec<usize>> {
        let mut order = vec![self.zero];
        let mut seen = vec![false; self.q.len()];
        seen[self.zero] = true;
        let mut cursor = 0usize;
        while cursor < order.len() {
            let x = order[cursor];
            for &g in gens {
                let nx = self.add[x][g];
                if !seen[nx] {
                    seen[nx] = true;
                    order.push(nx);
                }
            }
            cursor += 1;
        }
        (order.len() == self.q.len()).then_some(order)
    }

    fn label_for_order(&self, order: &[usize]) -> Option<Vec<i128>> {
        let mut pos = vec![usize::MAX; self.q.len()];
        for (i, &old) in order.iter().enumerate() {
            pos[old] = i;
        }
        let mut out = Vec::with_capacity(2 + 2 * order.len() + order.len() * order.len());
        out.push(i128::try_from(order.len()).ok()?);
        for &old in order {
            out.push(self.q[old].numer());
            out.push(self.q[old].denom());
        }
        for &x in order {
            for &y in order {
                out.push(i128::try_from(pos[self.add[x][y]]).ok()?);
            }
        }
        Some(out)
    }

    fn raw_label(&self) -> Option<Vec<i128>> {
        let order = (0..self.q.len()).collect::<Vec<_>>();
        self.label_for_order(&order)
    }

    fn q_value_counts(&self) -> Vec<FqmValueCount> {
        let mut counts: BTreeMap<(i128, i128), u128> = BTreeMap::new();
        for q in &self.q {
            *counts.entry((q.numer(), q.denom())).or_default() += 1;
        }
        counts
            .into_iter()
            .map(|((numer, denom), count)| FqmValueCount {
                numer,
                denom,
                count,
            })
            .collect()
    }

    fn primary_invariant_factors(&self, p: u128) -> Option<Vec<u128>> {
        let exponent = self.order.iter().copied().max().unwrap_or(1) as u128;
        let max_power = exact_prime_power_exponent(exponent, p)?;
        let mut killed_log = vec![0u128; usize::try_from(max_power + 1).ok()?];
        killed_log[0] = 0;
        let mut p_to_j = 1u128;
        for j in 1..=max_power {
            p_to_j = p_to_j.checked_mul(p)?;
            let count = self
                .order
                .iter()
                .filter(|&&ord| p_to_j.is_multiple_of(ord as u128))
                .count() as u128;
            killed_log[usize::try_from(j).ok()?] = exact_prime_power_exponent(count, p)?;
        }

        let mut ge = vec![0u128; usize::try_from(max_power + 2).ok()?];
        for j in 1..=max_power {
            let ji = usize::try_from(j).ok()?;
            ge[ji] = killed_log[ji].checked_sub(killed_log[ji - 1])?;
        }
        let mut factors = Vec::new();
        for j in 1..=max_power {
            let ji = usize::try_from(j).ok()?;
            let exact = ge[ji].checked_sub(ge[ji + 1])?;
            let factor = pow_u128(p, j)?;
            for _ in 0..exact {
                factors.push(factor);
            }
        }
        Some(factors)
    }

    fn bilinear_value(&self, x: usize, y: usize) -> Rational {
        let diff = self.q[self.add[x][y]].sub(&self.q[x]).sub(&self.q[y]);
        rational_half_mod1(diff)
    }

    fn quadratic_values_are_even(&self) -> bool {
        (0..self.q.len()).all(|x| {
            let nx = self.neg(x);
            self.q[nx] == self.q[x]
        })
    }

    fn bilinear_form_is_biadditive(&self) -> bool {
        for x in 0..self.q.len() {
            for y in 0..self.q.len() {
                for z in 0..self.q.len() {
                    let yz = self.add[y][z];
                    let lhs = self.bilinear_value(x, yz);
                    let rhs = rational_mod_int(
                        self.bilinear_value(x, y).add(&self.bilinear_value(x, z)),
                        1,
                    );
                    if lhs != rhs {
                        return false;
                    }
                    let xy = self.add[x][y];
                    let lhs = self.bilinear_value(xy, z);
                    let rhs = rational_mod_int(
                        self.bilinear_value(x, z).add(&self.bilinear_value(y, z)),
                        1,
                    );
                    if lhs != rhs {
                        return false;
                    }
                }
            }
        }
        true
    }

    fn is_nondegenerate(&self) -> bool {
        (0..self.q.len()).all(|x| {
            x == self.zero
                || (0..self.q.len()).any(|y| self.bilinear_value(x, y) != Rational::zero())
        })
    }

    fn neg(&self, x: usize) -> usize {
        (0..self.q.len())
            .find(|&y| self.add[x][y] == self.zero)
            .expect("finite abelian group element has an inverse")
    }
}

fn coords_from_index(mut index: usize, factors: &[u128]) -> Option<Vec<u128>> {
    let mut out = vec![0u128; factors.len()];
    for (i, &d) in factors.iter().enumerate().rev() {
        let du = usize::try_from(d).ok()?;
        out[i] = (index % du) as u128;
        index /= du;
    }
    Some(out)
}

fn index_from_coords(coords: &[u128], factors: &[u128]) -> Option<usize> {
    let mut out = 0usize;
    for (&x, &d) in coords.iter().zip(factors) {
        if x >= d {
            return None;
        }
        out = out
            .checked_mul(usize::try_from(d).ok()?)?
            .checked_add(usize::try_from(x).ok()?)?;
    }
    Some(out)
}

fn rational_mod_int(x: Rational, modulus: i128) -> Rational {
    debug_assert!(modulus > 0);
    let den = x.denom();
    let mden = den
        .checked_mul(modulus)
        .expect("rational modulus exceeds i128");
    Rational::new(x.numer().rem_euclid(mden), den)
}

fn rational_half_mod1(x: Rational) -> Rational {
    let den = x
        .denom()
        .checked_mul(2)
        .expect("rational denominator exceeds i128");
    rational_mod_int(Rational::new(x.numer(), den), 1)
}

fn rational_from_u128(n: u128) -> Option<Rational> {
    Some(Rational::from_int(i128::try_from(n).ok()?))
}

fn signed_order_for_odd_prime(order: u128, t_minus: usize) -> Option<i128> {
    let order = i128::try_from(order).ok()?;
    Some(if t_minus.is_multiple_of(2) {
        order
    } else {
        order.checked_neg()?
    })
}

fn v_p_i128(mut x: i128, p: i128) -> i128 {
    debug_assert!(x != 0);
    let mut k = 0i128;
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

fn rat_val(r: &Rational, p: i128) -> i128 {
    v_p_i128(r.numer(), p) - v_p_i128(r.denom(), p)
}

fn odd_unit_residue(r: &Rational, p: i128) -> i128 {
    let a = unit_part_i128(r.numer(), p).rem_euclid(p);
    let b = unit_part_i128(r.denom(), p).rem_euclid(p);
    // For square-class purposes b and b^{-1} have the same Legendre symbol.
    (a * b).rem_euclid(p)
}

fn unit_mod8(r: &Rational) -> i128 {
    let a = unit_part_i128(r.numer(), 2).rem_euclid(8);
    let b = unit_part_i128(r.denom(), 2).rem_euclid(8);
    // Odd units are self-inverse modulo 8.
    (a * b).rem_euclid(8)
}

fn same_square_class_odd(a: &Rational, b: &Rational, p: u128) -> Option<bool> {
    if a.is_zero() || b.is_zero() || p == 2 {
        return None;
    }
    let p_i = i128::try_from(p).ok()?;
    let ratio = a.mul(&b.inv()?);
    if rat_val(&ratio, p_i) % 2 != 0 {
        return Some(false);
    }
    try_is_square_qp(odd_unit_residue(&ratio, p_i), p)
}

fn same_square_class_2_up_to_sign(a: &Rational, b: &Rational) -> Option<bool> {
    if a.is_zero() || b.is_zero() {
        return None;
    }
    let ratio = a.mul(&b.inv()?);
    if rat_val(&ratio, 2) % 2 != 0 {
        return Some(false);
    }
    Some(matches!(unit_mod8(&ratio), 1 | 7))
}

fn rational_det(mut a: Vec<Vec<Rational>>) -> Option<Rational> {
    let n = a.len();
    if a.iter().any(|row| row.len() != n) {
        return None;
    }
    let mut det = Rational::one();
    for i in 0..n {
        let pivot = (i..n).find(|&r| !a[r][i].is_zero())?;
        if pivot != i {
            a.swap(i, pivot);
            det = det.neg();
        }
        let pivot_value = a[i][i].clone();
        det = det.mul(&pivot_value);
        let pivot_inv = pivot_value.inv()?;
        for r in (i + 1)..n {
            if a[r][i].is_zero() {
                continue;
            }
            let factor = a[r][i].mul(&pivot_inv);
            for c in i..n {
                let correction = factor.mul(&a[i][c]);
                a[r][c] = a[r][c].sub(&correction);
            }
        }
    }
    Some(det)
}

fn prime_factors_u128(n: u128) -> Vec<u128> {
    let mut m = n;
    let mut out = Vec::new();
    let mut p = 2u128;
    while p <= m / p {
        if m.is_multiple_of(p) {
            out.push(p);
            while m.is_multiple_of(p) {
                m /= p;
            }
        }
        p += if p == 2 { 1 } else { 2 };
    }
    if m > 1 {
        out.push(m);
    }
    out
}

fn is_prime_power_order(order: u128, p: u128) -> bool {
    if order == 1 {
        return true;
    }
    let mut m = order;
    while m.is_multiple_of(p) {
        m /= p;
    }
    m == 1
}

fn exact_prime_power_exponent(mut n: u128, p: u128) -> Option<u128> {
    if n == 1 {
        return Some(0);
    }
    let mut k = 0u128;
    while n > 1 && n.is_multiple_of(p) {
        n /= p;
        k += 1;
    }
    (n == 1).then_some(k)
}

fn pow_u128(base: u128, exp: u128) -> Option<u128> {
    let mut out = 1u128;
    for _ in 0..exp {
        out = out.checked_mul(base)?;
    }
    Some(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forms::{a_n, e_6, e_7, DiscriminantForm, IntegralForm};

    #[test]
    fn native_cyclic_module_matches_lattice_a1() {
        let native = FiniteQuadraticModule::cyclic(2, Rational::new(1, 2)).unwrap();
        let from_lattice = DiscriminantForm::from_lattice(&a_n(1))
            .unwrap()
            .fqm_witt_class()
            .unwrap();
        assert_eq!(native.witt_class().unwrap(), from_lattice);
    }

    #[test]
    fn fqm_witt_reduces_hyperbolic_two_primary_pair() {
        let a1 = DiscriminantForm::from_lattice(&a_n(1)).unwrap();
        let e7 = DiscriminantForm::from_lattice(&e_7()).unwrap();
        assert_ne!(a1.fqm_witt_class().unwrap(), e7.fqm_witt_class().unwrap());

        let hyperbolic = FiniteQuadraticModule::cyclic(2, Rational::new(1, 2))
            .unwrap()
            .direct_sum(&FiniteQuadraticModule::cyclic(2, Rational::new(3, 2)).unwrap())
            .unwrap();
        let class = hyperbolic.witt_class().unwrap();
        assert!(class.is_trivial());
        assert_eq!(class.phase_mod8, 0);
    }

    #[test]
    fn fqm_witt_reduces_hyperbolic_odd_primary_pair() {
        let a2 = DiscriminantForm::from_lattice(&a_n(2)).unwrap();
        let e6 = DiscriminantForm::from_lattice(&e_6()).unwrap();
        assert_eq!(a2.is_fqm_witt_equivalent(&e6), Some(false));
        assert!(DiscriminantForm::from_lattice(&a_n(2).direct_sum(&e_6()))
            .unwrap()
            .fqm_witt_class()
            .unwrap()
            .is_trivial());

        let sum = FiniteQuadraticModule::cyclic(3, Rational::new(2, 3))
            .unwrap()
            .direct_sum(&FiniteQuadraticModule::cyclic(3, Rational::new(4, 3)).unwrap())
            .unwrap();
        let class = sum.witt_class().unwrap();
        assert!(class.is_trivial());
        assert_eq!(class.phase_mod8, 0);
    }

    #[test]
    fn fqm_witt_refines_phase_projection() {
        let a1 = DiscriminantForm::from_lattice(&a_n(1)).unwrap();
        let class = a1.fqm_witt_class().unwrap();
        let phase = a1.fqm_gauss_phase().unwrap();
        assert_eq!(class.order, phase.order as u128);
        assert_eq!(class.phase_mod8, phase.phase_mod8);
        assert_eq!(class.primary[0].phase_mod8, phase.primary[0].phase_mod8);
        assert_eq!(class.primary[0].core_group, vec![2]);
        assert_eq!(class.primary[0].q_value_counts.len(), 2);
    }

    #[test]
    fn nikulin_existence_accepts_realized_lattice_discriminant_forms() {
        for lattice in [a_n(1), a_n(2), e_6(), e_7()] {
            let signature = lattice.signature();
            let q = DiscriminantForm::from_lattice(&lattice).unwrap();
            let report = q.nikulin_existence_report(signature).unwrap();
            assert!(
                report.exists(),
                "realized lattice should pass Nikulin 1.10.1"
            );
            assert_eq!(q.nikulin_even_lattice_exists(signature), Some(true));
        }
    }

    #[test]
    fn nikulin_existence_keeps_odd_two_primary_boundary() {
        let q = DiscriminantForm::from_lattice(&a_n(1)).unwrap();
        let report = q.nikulin_existence_report((1, 0)).unwrap();
        assert!(report.exists());
        assert_eq!(report.primary.len(), 1);
        assert_eq!(report.primary[0].prime, 2);
        assert_eq!(report.primary[0].length, 1);
        assert!(report.primary[0].equality_case);
        assert!(!report.primary[0].even_two_primary);
        assert_eq!(report.primary[0].determinant_condition_holds, None);

        let blocked = q.nikulin_existence_report((0, 1)).unwrap();
        assert_eq!(
            blocked.obstruction,
            Some(NikulinExistenceObstruction::SignatureCongruence {
                required_mod8: 7,
                module_phase_mod8: 1,
            })
        );
    }

    #[test]
    fn nikulin_existence_checks_odd_primary_borderline() {
        let hyperbolic_three = FiniteQuadraticModule::cyclic(3, Rational::new(2, 3))
            .unwrap()
            .direct_sum(&FiniteQuadraticModule::cyclic(3, Rational::new(4, 3)).unwrap())
            .unwrap();

        let report = hyperbolic_three.nikulin_existence_report((1, 1)).unwrap();
        assert!(report.exists());
        assert_eq!(report.primary.len(), 1);
        assert_eq!(report.primary[0].prime, 3);
        assert_eq!(report.primary[0].length, 2);
        assert!(report.primary[0].equality_case);
        assert_eq!(report.primary[0].determinant_condition_holds, Some(true));

        let too_small = hyperbolic_three.nikulin_existence_report((0, 0)).unwrap();
        assert_eq!(
            too_small.obstruction,
            Some(NikulinExistenceObstruction::RankTooSmall {
                prime: 3,
                rank: 0,
                length: 2,
            })
        );
    }

    #[test]
    fn nikulin_existence_checks_even_two_primary_borderline() {
        let u2 = IntegralForm::new(vec![vec![0, 2], vec![2, 0]]).unwrap();
        let q = DiscriminantForm::from_lattice(&u2).unwrap();
        let report = q.nikulin_existence_report((1, 1)).unwrap();
        assert!(report.exists());
        assert_eq!(report.primary.len(), 1);
        assert_eq!(report.primary[0].prime, 2);
        assert_eq!(report.primary[0].length, 2);
        assert!(report.primary[0].equality_case);
        assert!(report.primary[0].even_two_primary);
        assert_eq!(report.primary[0].determinant_condition_holds, Some(true));
    }
}
