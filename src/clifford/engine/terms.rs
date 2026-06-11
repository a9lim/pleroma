use super::basis::wedge_sign;
use crate::scalar::Scalar;
use std::collections::BTreeMap;

pub(super) fn scale<S: Scalar>(mut terms: BTreeMap<u128, S>, s: &S) -> BTreeMap<u128, S> {
    if s.is_zero() {
        return BTreeMap::new();
    }
    for v in terms.values_mut() {
        *v = v.mul(s);
    }
    terms.retain(|_, v| !v.is_zero());
    terms
}

pub(super) fn merge<S: Scalar>(into: &mut BTreeMap<u128, S>, other: BTreeMap<u128, S>) {
    for (blade, coeff) in other {
        add_term(into, blade, coeff);
    }
}

/// Insert `coeff` for `blade` into `out`, adding to any existing coefficient.
/// If the result is zero it is removed, preserving the "zeros never stored" invariant.
pub(super) fn add_term<S: Scalar>(out: &mut BTreeMap<u128, S>, blade: u128, coeff: S) {
    let e = out.entry(blade).or_insert_with(S::zero);
    *e = e.add(&coeff);
    if e.is_zero() {
        out.remove(&blade);
    }
}

/// The exterior (wedge) product of two term maps — the shared implementation
/// used by both `Multivector::bitand` (`&` operator) and `CliffordAlgebra::wedge`.
/// Metric-independent.
pub(super) fn wedge_terms<S: Scalar>(
    a: &BTreeMap<u128, S>,
    b: &BTreeMap<u128, S>,
) -> BTreeMap<u128, S> {
    let mut out: BTreeMap<u128, S> = BTreeMap::new();
    for (&ba, ca) in a {
        for (&bb, cb) in b {
            if ba & bb != 0 {
                continue;
            }
            let coeff = ca.mul(cb).mul(&wedge_sign::<S>(ba, bb));
            if !coeff.is_zero() {
                add_term(&mut out, ba | bb, coeff);
            }
        }
    }
    out
}
