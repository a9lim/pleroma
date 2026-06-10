//! The classifier façade: one entry point that dispatches on the scalar field.
//!
//! The three characteristic legs ([`char0`](crate::forms::char0),
//! [`oddchar`](crate::forms::oddchar), [`char2`](crate::forms::char2)) each ship
//! their own classifier with a leg-specific signature — `classify_surreal`,
//! `classify_finite_odd`, `arf_invariant`, … Choosing the right one is a fact
//! about the field, not the form, so it can be resolved *at compile time* from
//! the scalar type. [`ClassifyForm`] does exactly that: write
//! `metric.classify()` (or `S::classify(metric)`) and the correct leg is
//! selected by the monomorphised `S` — no manual `match` on characteristic.
//!
//! [`WittClassify`] is the same idea for the unified [`WittClassG`], over the
//! three legs where a single Witt class exists (real char 0, odd char, char 2).
//! `Rational`'s Witt invariant is the full Hasse–Minkowski datum and surcomplex's
//! is `W(ℂ) = ℤ/2`; neither is a `WittClassG`, so those two backends implement
//! [`ClassifyForm`] but not [`WittClassify`] — honest, not a gap.

use crate::clifford::{CliffordAlgebra, Metric};
use crate::forms::{
    arf_fpn_char2, arf_invariant, arf_ordinal_finite, bw_class_complex, bw_class_finite_odd,
    bw_class_nimber, bw_class_real, classify_finite_odd, classify_rational, classify_surcomplex,
    classify_surreal, finite_odd_witt, isometric_finite_odd, isometric_fpn_char2, isometric_nimber,
    isometric_ordinal_finite, isometric_rational, isometric_real, isometric_surcomplex,
    witt_decompose_finite_odd, witt_decompose_real, ArfResult, BrauerWallClass, CliffordType,
    OddCharType, OddWittDecomp, RationalCliffordType, RealWittDecomp, WittClassG,
};
use crate::scalar::{nim_degree, Fp, Fpn, Nimber, Ordinal, Rational, Scalar, Surcomplex, Surreal};

/// Classification data for the `Fpn<P,N>` finite-field tower. Odd-characteristic
/// extension fields land in the usual finite-odd invariant; characteristic-2
/// extension fields land in the Arf invariant.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FiniteFieldClass {
    /// Finite field of odd characteristic.
    Odd(OddCharType),
    /// Finite field of characteristic 2.
    Char2(ArfResult),
}

/// Witt-decomposition data for the finite-field tower `Fpn<P,N>`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FiniteFieldWittDecomp {
    /// Odd-characteristic finite field.
    Odd(OddWittDecomp),
    /// Characteristic-2 finite field.
    Char2(Char2WittDecomp),
}

/// Witt-decomposition data for a characteristic-2 finite-field form.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Char2WittDecomp {
    /// Extension degree `m` for `F_{2^m}`.
    pub field_degree: u128,
    /// Number of hyperbolic planes split from the nonsingular core.
    pub witt_index: usize,
    /// Dimension of the anisotropic nonsingular core: `0` for hyperbolic, `2` for
    /// the anisotropic plane.
    pub core_anisotropic_dim: usize,
    /// Dimension of the polar radical.
    pub radical_dim: usize,
    /// Whether the quadratic form is nonzero on the radical.
    pub radical_anisotropic: bool,
    /// Arf bit of the nonsingular core.
    pub arf: u128,
}

impl Char2WittDecomp {
    fn from_arf(field_degree: u128, arf: &ArfResult) -> Self {
        let core_anisotropic_dim = if arf.arf == 0 { 0 } else { 2 };
        let witt_index = arf.rank.saturating_sub(core_anisotropic_dim) / 2;
        Char2WittDecomp {
            field_degree,
            witt_index,
            core_anisotropic_dim,
            radical_dim: arf.radical_dim,
            radical_anisotropic: arf.radical_anisotropic,
            arf: arf.arf,
        }
    }
}

impl FiniteFieldClass {
    pub fn display(&self) -> String {
        match self {
            FiniteFieldClass::Odd(c) => c.display(),
            FiniteFieldClass::Char2(c) => {
                format!(
                    "char 2: Arf {} rank {} radical {}{} ({})",
                    c.arf,
                    c.rank,
                    c.radical_dim,
                    if c.radical_anisotropic {
                        " defective"
                    } else {
                        ""
                    },
                    c.o_type
                )
            }
        }
    }
}

/// Classify the quadratic form data attached to a [`Metric`] over `Self`,
/// dispatched on the scalar field. For real/complex-style legs this is also the
/// implemented Clifford-algebra closure type; for `Rational` it is the complete
/// Hasse-Minkowski quadratic-form invariant, not a full rational Brauer-Wall
/// class. The [`Class`](ClassifyForm::Class) associated type is the leg-specific
/// datum:
///
/// | scalar | `Class` | leg |
/// |---|---|---|
/// | [`Surreal`] | [`CliffordType`] | exact-square char 0 subdomain (8-fold) |
/// | [`Surcomplex<Surreal>`](Surcomplex) | [`CliffordType`] | exact-square char 0 subdomain (2-fold) |
/// | [`Rational`] | [`RationalCliffordType`] | char 0, full Hasse-Minkowski form invariant |
/// | [`Fp<P>`](Fp) | [`OddCharType`] | odd characteristic |
/// | [`Fpn<P,N>`](Fpn) | [`FiniteFieldClass`] | finite extension fields, odd or char 2 |
/// | [`Nimber`] | [`ArfResult`] | characteristic 2 (Arf) |
/// | [`Ordinal`] | [`ArfResult`] | detected finite ordinal-nimber windows only |
///
/// `None` means the metric is outside the classifier's domain (e.g. a non-diagonal
/// char-2 form, or a metric the diagonalizer can't reduce).
pub trait ClassifyForm: Scalar {
    /// The classification datum produced for this field's characteristic leg.
    type Class;

    /// Classify the form carried by `metric`.
    fn classify(metric: &Metric<Self>) -> Option<Self::Class>;
}

/// The unified Witt class [`WittClassG`] of a form, for the three legs where a
/// single Witt class exists. (`Rational` and `Surcomplex` deliberately do not
/// implement this — see the module docs.)
pub trait WittClassify: Scalar {
    /// The Witt class of the form carried by `metric`.
    fn witt_class(metric: &Metric<Self>) -> Option<WittClassG>;
}

/// Isometry comparison for scalar worlds with a complete invariant available.
pub trait IsometryClassify: Scalar {
    /// Whether two forms over the same scalar world are isometric.
    fn isometric(m1: &Metric<Self>, m2: &Metric<Self>) -> Option<bool>;
}

/// Constructive Witt decomposition where the crate has a concrete decomposition
/// datum for that scalar world.
pub trait WittDecompose: Scalar {
    /// The decomposition datum for this scalar world.
    type Decomp;

    /// Split a form into hyperbolic planes plus anisotropic kernel data.
    fn witt_decompose(metric: &Metric<Self>) -> Option<Self::Decomp>;
}

/// Brauer-Wall class of the Clifford algebra attached to a form.
pub trait BrauerWallClassify: Scalar {
    /// The Brauer-Wall class of `Cl(metric)`.
    fn bw_class(metric: &Metric<Self>) -> Option<BrauerWallClass>;
}

impl ClassifyForm for Surreal {
    type Class = CliffordType;
    fn classify(metric: &Metric<Self>) -> Option<CliffordType> {
        classify_surreal(metric)
    }
}

impl ClassifyForm for Surcomplex<Surreal> {
    type Class = CliffordType;
    fn classify(metric: &Metric<Self>) -> Option<CliffordType> {
        classify_surcomplex(metric)
    }
}

impl ClassifyForm for Rational {
    type Class = RationalCliffordType;
    fn classify(metric: &Metric<Self>) -> Option<RationalCliffordType> {
        classify_rational(metric)
    }
}

impl<const P: u128> ClassifyForm for Fp<P> {
    type Class = OddCharType;
    fn classify(metric: &Metric<Self>) -> Option<OddCharType> {
        classify_finite_odd(metric)
    }
}

impl<const P: u128, const N: usize> ClassifyForm for Fpn<P, N> {
    type Class = FiniteFieldClass;
    fn classify(metric: &Metric<Self>) -> Option<FiniteFieldClass> {
        if P == 2 {
            arf_fpn_char2(metric).map(FiniteFieldClass::Char2)
        } else {
            classify_finite_odd(metric).map(FiniteFieldClass::Odd)
        }
    }
}

impl ClassifyForm for Nimber {
    type Class = ArfResult;
    fn classify(metric: &Metric<Self>) -> Option<ArfResult> {
        arf_invariant(metric)
    }
}

impl ClassifyForm for Ordinal {
    type Class = ArfResult;
    fn classify(metric: &Metric<Self>) -> Option<ArfResult> {
        arf_ordinal_finite(metric)
    }
}

impl WittClassify for Surreal {
    fn witt_class(metric: &Metric<Self>) -> Option<WittClassG> {
        let (p, q, _r) = crate::forms::char0::surreal_signature(metric)?;
        Some(WittClassG::char0(p, q))
    }
}

impl<const P: u128> WittClassify for Fp<P> {
    fn witt_class(metric: &Metric<Self>) -> Option<WittClassG> {
        finite_odd_witt(metric)
    }
}

impl<const P: u128, const N: usize> WittClassify for Fpn<P, N> {
    fn witt_class(metric: &Metric<Self>) -> Option<WittClassG> {
        if P == 2 {
            let arf = arf_fpn_char2(metric)?;
            if arf.radical_dim != 0 {
                return None;
            }
            Some(WittClassG::Char2 {
                field_degree: N as u128,
                arf: arf.arf,
            })
        } else {
            finite_odd_witt(metric)
        }
    }
}

impl WittClassify for Nimber {
    fn witt_class(metric: &Metric<Self>) -> Option<WittClassG> {
        WittClassG::try_char2_from_metric(metric).ok()
    }
}

impl WittClassify for Ordinal {
    fn witt_class(metric: &Metric<Self>) -> Option<WittClassG> {
        let arf = arf_ordinal_finite(metric)?;
        if arf.radical_dim != 0 {
            return None;
        }
        Some(WittClassG::Char2 {
            field_degree: ordinal_char2_field_degree(metric)?,
            arf: arf.arf,
        })
    }
}

impl IsometryClassify for Surreal {
    fn isometric(m1: &Metric<Self>, m2: &Metric<Self>) -> Option<bool> {
        isometric_real(m1, m2)
    }
}

impl IsometryClassify for Surcomplex<Surreal> {
    fn isometric(m1: &Metric<Self>, m2: &Metric<Self>) -> Option<bool> {
        isometric_surcomplex(m1, m2)
    }
}

impl IsometryClassify for Rational {
    fn isometric(m1: &Metric<Self>, m2: &Metric<Self>) -> Option<bool> {
        isometric_rational(m1, m2)
    }
}

impl<const P: u128> IsometryClassify for Fp<P> {
    fn isometric(m1: &Metric<Self>, m2: &Metric<Self>) -> Option<bool> {
        isometric_finite_odd(m1, m2)
    }
}

impl<const P: u128, const N: usize> IsometryClassify for Fpn<P, N> {
    fn isometric(m1: &Metric<Self>, m2: &Metric<Self>) -> Option<bool> {
        if P == 2 {
            isometric_fpn_char2(m1, m2)
        } else {
            isometric_finite_odd(m1, m2)
        }
    }
}

impl IsometryClassify for Nimber {
    fn isometric(m1: &Metric<Self>, m2: &Metric<Self>) -> Option<bool> {
        isometric_nimber(m1, m2)
    }
}

impl IsometryClassify for Ordinal {
    fn isometric(m1: &Metric<Self>, m2: &Metric<Self>) -> Option<bool> {
        isometric_ordinal_finite(m1, m2)
    }
}

impl WittDecompose for Surreal {
    type Decomp = RealWittDecomp;
    fn witt_decompose(metric: &Metric<Self>) -> Option<Self::Decomp> {
        witt_decompose_real(metric)
    }
}

impl<const P: u128> WittDecompose for Fp<P> {
    type Decomp = OddWittDecomp;
    fn witt_decompose(metric: &Metric<Self>) -> Option<Self::Decomp> {
        witt_decompose_finite_odd(metric)
    }
}

impl<const P: u128, const N: usize> WittDecompose for Fpn<P, N> {
    type Decomp = FiniteFieldWittDecomp;
    fn witt_decompose(metric: &Metric<Self>) -> Option<Self::Decomp> {
        if P == 2 {
            let arf = arf_fpn_char2(metric)?;
            Some(FiniteFieldWittDecomp::Char2(Char2WittDecomp::from_arf(
                N as u128, &arf,
            )))
        } else {
            witt_decompose_finite_odd(metric).map(FiniteFieldWittDecomp::Odd)
        }
    }
}

impl BrauerWallClassify for Surreal {
    fn bw_class(metric: &Metric<Self>) -> Option<BrauerWallClass> {
        bw_class_real(metric)
    }
}

impl BrauerWallClassify for Surcomplex<Surreal> {
    fn bw_class(metric: &Metric<Self>) -> Option<BrauerWallClass> {
        bw_class_complex(metric)
    }
}

impl<const P: u128> BrauerWallClassify for Fp<P> {
    fn bw_class(metric: &Metric<Self>) -> Option<BrauerWallClass> {
        bw_class_finite_odd(metric)
    }
}

impl<const P: u128, const N: usize> BrauerWallClassify for Fpn<P, N> {
    fn bw_class(metric: &Metric<Self>) -> Option<BrauerWallClass> {
        if P == 2 {
            let arf = arf_fpn_char2(metric)?;
            if arf.radical_dim != 0 {
                return None;
            }
            Some(BrauerWallClass::Char2 {
                field_degree: N as u128,
                arf: arf.arf,
            })
        } else {
            bw_class_finite_odd(metric)
        }
    }
}

impl BrauerWallClassify for Nimber {
    fn bw_class(metric: &Metric<Self>) -> Option<BrauerWallClass> {
        bw_class_nimber(metric)
    }
}

impl BrauerWallClassify for Ordinal {
    fn bw_class(metric: &Metric<Self>) -> Option<BrauerWallClass> {
        let arf = arf_ordinal_finite(metric)?;
        if arf.radical_dim != 0 {
            return None;
        }
        Some(BrauerWallClass::Char2 {
            field_degree: ordinal_char2_field_degree(metric)?,
            arf: arf.arf,
        })
    }
}

fn ordinal_char2_field_degree(metric: &Metric<Ordinal>) -> Option<u128> {
    if metric.q.iter().all(|x| x.as_finite().is_some())
        && metric.b.values().all(|x| x.as_finite().is_some())
    {
        return metric
            .q
            .iter()
            .map(|x| x.as_finite().map(nim_degree))
            .chain(metric.b.values().map(|x| x.as_finite().map(nim_degree)))
            .collect::<Option<Vec<_>>>()
            .map(|degrees| degrees.into_iter().max().unwrap_or(1));
    }

    if metric.q.iter().chain(metric.b.values()).all(|x| {
        x.as_below_omega3()
            .is_some_and(|cs| cs.iter().all(|&c| c < 4))
    }) {
        return Some(6);
    }

    None
}

/// Ergonomic methods so callers can write `metric.classify()` /
/// `algebra.classify()` instead of `S::classify(&metric)`.
impl<S: ClassifyForm> Metric<S> {
    /// Classify the form (see [`ClassifyForm`]).
    pub fn classify(&self) -> Option<S::Class> {
        S::classify(self)
    }
}

impl<S: WittClassify> Metric<S> {
    /// The unified Witt class (see [`WittClassify`]).
    pub fn witt_class(&self) -> Option<WittClassG> {
        S::witt_class(self)
    }
}

impl<S: IsometryClassify> Metric<S> {
    /// Test isometry against another form over the same scalar world.
    pub fn isometric_to(&self, other: &Self) -> Option<bool> {
        S::isometric(self, other)
    }
}

impl<S: WittDecompose> Metric<S> {
    /// Split the form into hyperbolic planes plus anisotropic kernel data.
    pub fn witt_decompose(&self) -> Option<S::Decomp> {
        S::witt_decompose(self)
    }
}

impl<S: BrauerWallClassify> Metric<S> {
    /// The Brauer-Wall class of the attached Clifford algebra.
    pub fn bw_class(&self) -> Option<BrauerWallClass> {
        S::bw_class(self)
    }
}

impl<S: ClassifyForm> CliffordAlgebra<S> {
    /// Classify the algebra's underlying form (see [`ClassifyForm`]).
    pub fn classify(&self) -> Option<S::Class> {
        S::classify(&self.metric)
    }
}

impl<S: WittClassify> CliffordAlgebra<S> {
    /// The unified Witt class of the algebra's form (see [`WittClassify`]).
    pub fn witt_class(&self) -> Option<WittClassG> {
        S::witt_class(&self.metric)
    }
}

impl<S: IsometryClassify> CliffordAlgebra<S> {
    /// Test isometry of the underlying forms.
    pub fn isometric_to(&self, other: &Self) -> Option<bool> {
        S::isometric(&self.metric, &other.metric)
    }
}

impl<S: WittDecompose> CliffordAlgebra<S> {
    /// Witt decomposition of the algebra's underlying form.
    pub fn witt_decompose(&self) -> Option<S::Decomp> {
        S::witt_decompose(&self.metric)
    }
}

impl<S: BrauerWallClassify> CliffordAlgebra<S> {
    /// Brauer-Wall class of the algebra.
    pub fn bw_class(&self) -> Option<BrauerWallClass> {
        S::bw_class(&self.metric)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::clifford::Metric;

    #[test]
    fn classify_dispatches_on_scalar_type() {
        // char 0, real-closed: Cl(2,0) over the surreals matches classify_surreal.
        let m = Metric::diagonal(vec![Surreal::one(), Surreal::one()]);
        assert_eq!(m.classify(), classify_surreal(&m));
        assert!(m.classify().is_some());

        // char 2: Arf via the trait matches arf_invariant, and witt_class agrees.
        let n = Metric::diagonal(vec![Nimber::one(), Nimber::one()]);
        assert_eq!(n.classify(), arf_invariant(&n));
        assert_eq!(n.witt_class(), WittClassG::try_char2_from_metric(&n).ok());
        assert_eq!(n.bw_class(), bw_class_nimber(&n));

        // odd char: F_5 dispatch produces the odd-char datum.
        let f = Metric::diagonal(vec![Fp::<5>::new(1), Fp::<5>::new(2)]);
        assert_eq!(f.classify(), classify_finite_odd(&f));
        assert_eq!(f.witt_class(), finite_odd_witt(&f));

        // finite extension field: the same façade reaches the generic odd-field leg.
        let f9 = Metric::diagonal(vec![Fpn::<3, 2>::constant(1), Fpn::<3, 2>::generator()]);
        assert_eq!(
            f9.classify(),
            classify_finite_odd(&f9).map(FiniteFieldClass::Odd)
        );
        assert_eq!(f9.witt_class(), finite_odd_witt(&f9));

        // finite extension field, characteristic 2: the same façade now reaches
        // the generic Arf leg rather than falling through the odd-char classifier.
        let mut b = std::collections::BTreeMap::new();
        b.insert((0usize, 1usize), Fpn::<2, 3>::one());
        let f8 = Metric::new(vec![Fpn::<2, 3>::generator(), Fpn::<2, 3>::generator()], b);
        assert_eq!(
            f8.classify(),
            arf_fpn_char2(&f8).map(FiniteFieldClass::Char2)
        );
        assert!(matches!(f8.classify(), Some(FiniteFieldClass::Char2(_))));

        // ordinal-nimber coefficients classify only inside detected finite
        // windows; the first transfinite one is F_4(ω) = F_64.
        let mut b = std::collections::BTreeMap::new();
        b.insert((0usize, 1usize), Ordinal::one());
        let omega = Ordinal::omega();
        let ord = Metric::new(vec![omega.clone(), omega], b);
        let arf = arf_ordinal_finite(&ord).unwrap();
        assert_eq!(ord.classify(), Some(arf.clone()));
        assert_eq!(
            ord.witt_class(),
            Some(WittClassG::Char2 {
                field_degree: 6,
                arf: arf.arf
            })
        );
        assert_eq!(
            ord.bw_class(),
            Some(BrauerWallClass::Char2 {
                field_degree: 6,
                arf: arf.arf
            })
        );

        let outside_window = Metric::diagonal(vec![Ordinal::omega_pow(Ordinal::omega())]);
        assert_eq!(outside_window.classify(), None);
        assert_eq!(outside_window.bw_class(), None);
    }

    #[test]
    fn algebra_classify_matches_metric_classify() {
        let alg = CliffordAlgebra::new(
            2,
            Metric::diagonal(vec![Surreal::one(), Surreal::one().neg()]),
        );
        assert_eq!(alg.classify(), alg.metric.classify());
        assert_eq!(alg.witt_class(), alg.metric.witt_class());
        assert_eq!(alg.witt_decompose(), alg.metric.witt_decompose());
        assert_eq!(alg.bw_class(), alg.metric.bw_class());
    }

    #[test]
    fn structural_facades_dispatch() {
        let f = Metric::diagonal(vec![Fp::<5>::new(1), Fp::<5>::new(1)]);
        let g = Metric::diagonal(vec![Fp::<5>::new(2), Fp::<5>::new(3)]);
        assert_eq!(f.isometric_to(&g), isometric_finite_odd(&f, &g));
        assert_eq!(f.witt_decompose(), witt_decompose_finite_odd(&f));
        assert_eq!(f.bw_class(), bw_class_finite_odd(&f));

        let f9 = Metric::diagonal(vec![Fpn::<3, 2>::constant(1), Fpn::<3, 2>::constant(1)]);
        let g9 = Metric::diagonal(vec![Fpn::<3, 2>::constant(2), Fpn::<3, 2>::constant(2)]);
        assert_eq!(f9.isometric_to(&g9), isometric_finite_odd(&f9, &g9));
        assert_eq!(
            f9.witt_decompose(),
            witt_decompose_finite_odd(&f9).map(FiniteFieldWittDecomp::Odd)
        );
        assert_eq!(f9.bw_class(), bw_class_finite_odd(&f9));

        let mut b = std::collections::BTreeMap::new();
        b.insert((0usize, 1usize), Fpn::<2, 3>::one());
        let f8 = Metric::new(vec![Fpn::<2, 3>::zero(), Fpn::<2, 3>::zero()], b);
        assert_eq!(
            f8.witt_decompose(),
            Some(FiniteFieldWittDecomp::Char2(Char2WittDecomp {
                field_degree: 3,
                witt_index: 1,
                core_anisotropic_dim: 0,
                radical_dim: 0,
                radical_anisotropic: false,
                arf: 0,
            }))
        );
        assert_eq!(
            f8.bw_class(),
            Some(BrauerWallClass::Char2 {
                field_degree: 3,
                arf: 0
            })
        );

        let mut b = std::collections::BTreeMap::new();
        b.insert((0usize, 1usize), Nimber::one());
        let n = Metric::new(vec![Nimber::zero(), Nimber::zero()], b);
        assert_eq!(n.bw_class(), bw_class_nimber(&n));

        let mut b = std::collections::BTreeMap::new();
        b.insert((0usize, 1usize), Ordinal::one());
        let ord = Metric::new(vec![Ordinal::omega(), Ordinal::omega()], b);
        assert_eq!(ord.isometric_to(&ord), Some(true));
    }
}
