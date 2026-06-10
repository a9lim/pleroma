use super::basis::MAX_BASIS_DIM;
use crate::scalar::Scalar;
use std::collections::BTreeMap;

/// Owned metric storage: `(q, b, a)`.
pub type MetricParts<S> = (
    Vec<S>,
    BTreeMap<(usize, usize), S>,
    BTreeMap<(usize, usize), S>,
);

/// The metric of a Clifford algebra of a (possibly general, possibly degenerate)
/// bilinear form. We store the full bilinear form `B(e_i, e_j)` factored into
/// three pieces so the common case is cheap and the general case is reachable:
///
///   * `q[i]   = B(e_i, e_i) = e_i²`                  — the quadratic form (diagonal)
///   * `b[(i,j)] = B(e_i,e_j) + B(e_j,e_i) = {e_i,e_j}`  (i<j) — the *polar* /
///     anticommutator form.
///   * `a[(i,j)] = B(e_i, e_j)` for i<j — the **strictly-upper / in-order
///     contraction** part.
#[derive(Clone, Debug, PartialEq)]
pub struct Metric<S: Scalar> {
    pub(crate) q: Vec<S>,
    pub(crate) b: BTreeMap<(usize, usize), S>,
    pub(crate) a: BTreeMap<(usize, usize), S>,
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
        let metric = Metric {
            q,
            b,
            a: BTreeMap::new(),
        };
        metric.validate_for_dim(metric.q.len());
        metric
    }

    /// A general-bilinear-form metric: squares `q`, polar form `b` (i<j), and the
    /// in-order contraction `a` (i<j). See the struct docs.
    pub fn general(
        q: Vec<S>,
        b: BTreeMap<(usize, usize), S>,
        a: BTreeMap<(usize, usize), S>,
    ) -> Self {
        let metric = Metric { q, b, a };
        metric.validate_for_dim(metric.q.len());
        metric
    }

    /// The represented dimension, i.e. the length of the quadratic diagonal.
    pub fn dim(&self) -> usize {
        self.q.len()
    }

    /// Diagonal quadratic entries `q[i] = e_i^2`.
    pub fn q(&self) -> &[S] {
        &self.q
    }

    /// Polar/anticommutator entries `b[(i,j)] = {e_i,e_j}` with `i < j`.
    pub fn b(&self) -> &BTreeMap<(usize, usize), S> {
        &self.b
    }

    /// Strictly-upper/in-order contraction entries with `i < j`.
    pub fn a(&self) -> &BTreeMap<(usize, usize), S> {
        &self.a
    }

    /// Consume the metric into its invariant-carrying parts.
    pub fn into_parts(self) -> MetricParts<S> {
        (self.q, self.b, self.a)
    }

    pub(crate) fn validate_for_dim(&self, dim: usize) {
        assert!(
            dim <= MAX_BASIS_DIM,
            "CliffordAlgebra supports at most {MAX_BASIS_DIM} generators"
        );
        assert_eq!(
            self.q.len(),
            dim,
            "metric q length must equal algebra dimension"
        );
        Self::validate_keys("b", &self.b, dim);
        Self::validate_keys("a", &self.a, dim);
    }

    fn validate_keys(name: &str, map: &BTreeMap<(usize, usize), S>, dim: usize) {
        for &(i, j) in map.keys() {
            assert!(i < j, "{name}-keys must satisfy i < j");
            assert!(
                j < dim,
                "{name}-key ({i},{j}) is out of range for dimension {dim}"
            );
        }
    }

    /// True iff there is any in-order contraction — i.e. this is a genuinely
    /// general bilinear form and needs the Chevalley product path.
    pub(crate) fn has_upper(&self) -> bool {
        self.a.values().any(|v| !v.is_zero())
    }

    /// True iff the basis is orthogonal: no off-diagonal polar or upper
    /// contraction terms are present. Then a blade product reduces to one blade.
    pub(crate) fn is_orthogonal(&self) -> bool {
        self.b.values().all(|v| v.is_zero()) && self.a.values().all(|v| v.is_zero())
    }

    /// Orthogonal direct sum `M ⟂ M'`: a block-diagonal metric on the disjoint
    /// union of the two generator sets.
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
        Metric::general(q, b, a)
    }

    pub(crate) fn q_val(&self, i: usize) -> S {
        self.q.get(i).cloned().unwrap_or_else(S::zero)
    }

    /// Base-change the metric by applying `f` to every coefficient (`q`, `b`, `a`).
    /// A form over `S` becomes the same form over `T` — e.g. lifting an `F_2`-valued
    /// trace form (`Metric<Fp<2>>`) into `Metric<Nimber>` so the char-2 Arf
    /// classifier can read it (`m.map(|x| Nimber(x.0))`). The structure (which
    /// `(i,j)` entries are present) is preserved verbatim; the caller is responsible
    /// for `f` being a ring map if the result is to mean anything.
    pub fn map<T: Scalar>(&self, f: impl Fn(&S) -> T) -> Metric<T> {
        Metric {
            q: self.q.iter().map(&f).collect(),
            b: self.b.iter().map(|(&k, v)| (k, f(v))).collect(),
            a: self.a.iter().map(|(&k, v)| (k, f(v))).collect(),
        }
    }
}
