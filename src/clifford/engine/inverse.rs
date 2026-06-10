use super::algebra::CliffordAlgebra;
use super::multivector::Multivector;
use crate::linalg::field;
use crate::scalar::Scalar;
use std::collections::BTreeMap;

impl<S: Scalar> CliffordAlgebra<S> {
    /// The **general multivector inverse** `v⁻¹` (two-sided), for any element.
    pub fn multivector_inverse(&self, v: &Multivector<S>) -> Option<Multivector<S>> {
        if v.is_zero() {
            return None;
        }
        if v.terms.len() == 1 {
            if let Some(c) = v.terms.get(&0) {
                return Some(self.scalar(c.inv()?));
            }
        }
        let n = 1usize.checked_shl(self.dim.try_into().ok()?)?;
        let mut mat = vec![vec![S::zero(); n]; n];
        for col in 0..n {
            let mut t = BTreeMap::new();
            t.insert(col as u128, S::one());
            let prod = self.mul(v, &Multivector { terms: t });
            for (&blade, c) in &prod.terms {
                mat[blade as usize][col] = c.clone();
            }
        }
        let mut rhs = vec![S::zero(); n];
        rhs[0] = S::one();
        let x = field::solve(mat, rhs)?;
        let mut terms = BTreeMap::new();
        for (bm, c) in x.into_iter().enumerate() {
            if !c.is_zero() {
                terms.insert(bm as u128, c);
            }
        }
        Some(Multivector { terms })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::clifford::Metric;
    use crate::scalar::Rational;

    #[test]
    fn inverse_refuses_huge_non_scalar_without_shift_overflow() {
        let alg: CliffordAlgebra<Rational> = CliffordAlgebra::new(64, Metric::grassmann(64));
        assert_eq!(alg.multivector_inverse(&alg.gen(0)), None);
    }

    #[test]
    fn scalar_inverse_still_works_at_huge_dimension() {
        let alg = CliffordAlgebra::new(64, Metric::grassmann(64));
        let two = alg.scalar(Rational::int(2));
        assert_eq!(
            alg.multivector_inverse(&two),
            Some(alg.scalar(Rational::new(1, 2)))
        );
    }
}
