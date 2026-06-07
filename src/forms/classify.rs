//! The classifier façade: one entry point that dispatches on the scalar field.
//!
//! The three characteristic legs ([`char0`](crate::forms::char0),
//! [`oddchar`](crate::forms::oddchar), [`char2`](crate::forms::char2)) each ship
//! their own classifier with a leg-specific signature — `classify_surreal`,
//! `classify_oddchar::<P>`, `arf_invariant`, … Choosing the right one is a fact
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
    arf_invariant, classify_oddchar, classify_rational, classify_surcomplex, classify_surreal,
    oddchar_witt, ArfResult, CliffordType, OddCharType, RationalCliffordType, WittClassG,
};
use crate::scalar::{Fp, Nimber, Rational, Scalar, Surcomplex, Surreal};

/// Classify the quadratic form / Clifford algebra of a [`Metric`] over `Self`,
/// dispatched on the scalar field. The [`Class`](ClassifyForm::Class) associated
/// type is the leg-specific datum:
///
/// | scalar | `Class` | leg |
/// |---|---|---|
/// | [`Surreal`] | [`CliffordType`] | real-closed char 0 (8-fold) |
/// | [`Surcomplex<Surreal>`](Surcomplex) | [`CliffordType`] | alg-closed char 0 (2-fold) |
/// | [`Rational`] | [`RationalCliffordType`] | char 0, full Hasse–Minkowski |
/// | [`Fp<P>`](Fp) | [`OddCharType`] | odd characteristic |
/// | [`Nimber`] | [`ArfResult`] | characteristic 2 (Arf) |
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
        classify_oddchar(metric)
    }
}

impl ClassifyForm for Nimber {
    type Class = ArfResult;
    fn classify(metric: &Metric<Self>) -> Option<ArfResult> {
        arf_invariant(metric)
    }
}

impl WittClassify for Surreal {
    fn witt_class(metric: &Metric<Self>) -> Option<WittClassG> {
        let (p, q, _r) = crate::forms::char0::signature(metric, |x| x.sign())?;
        Some(WittClassG::char0(p, q))
    }
}

impl<const P: u128> WittClassify for Fp<P> {
    fn witt_class(metric: &Metric<Self>) -> Option<WittClassG> {
        oddchar_witt(metric)
    }
}

impl WittClassify for Nimber {
    fn witt_class(metric: &Metric<Self>) -> Option<WittClassG> {
        WittClassG::try_char2_from_metric(metric).ok()
    }
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

        // odd char: F_5 dispatch produces the odd-char datum.
        let f = Metric::diagonal(vec![Fp::<5>::new(1), Fp::<5>::new(2)]);
        assert_eq!(f.classify(), classify_oddchar(&f));
        assert_eq!(f.witt_class(), oddchar_witt(&f));
    }

    #[test]
    fn algebra_classify_matches_metric_classify() {
        let alg = CliffordAlgebra::new(
            2,
            Metric::diagonal(vec![Surreal::one(), Surreal::one().neg()]),
        );
        assert_eq!(alg.classify(), alg.metric.classify());
        assert_eq!(alg.witt_class(), alg.metric.witt_class());
    }
}
