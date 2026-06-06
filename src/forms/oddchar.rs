//! The odd-characteristic Clifford / quadratic-form classifier — the third leg
//! of the trichotomy, companion to `classify.rs` (char 0) and `arf.rs` (char 2).
//!
//! Over a finite field `F_q` of odd characteristic a nondegenerate quadratic
//! form is classified completely by **dimension + discriminant** (det mod
//! squares): for each dimension there are exactly two classes, distinguished by
//! whether the discriminant is a square. So the classifier is essentially
//! `(dim, disc-class)`.
//!
//! We also compute the **Hasse–Witt / Clifford invariant** (a product of Hilbert
//! symbols). Over a finite field this is *always* `+1` — finite fields have
//! trivial Brauer group, so there are no nontrivial quaternion algebras and the
//! Hilbert symbol of any two nonzero elements is `+1`. We compute it the honest
//! way (search for a representing vector, which always exists by
//! Chevalley–Warning) precisely to *exhibit* that triviality, and to make the
//! structural parallel with the Arf invariant explicit — not because it adds
//! classifying power over a finite field. The group-theoretic home of all this
//! is `witt::WittClassG`.

use crate::clifford::Metric;
use crate::forms::WittClassG;
use crate::scalar::Fp;
use crate::scalar::Scalar;

/// `base^e` in `F_P` by square-and-multiply.
fn fp_pow<const P: u128>(mut base: Fp<P>, mut e: u128) -> Fp<P> {
    let mut acc = Fp::<P>::one();
    while e > 0 {
        if e & 1 == 1 {
            acc = acc.mul(&base);
        }
        base = base.mul(&base);
        e >>= 1;
    }
    acc
}

/// Euler's criterion: is `x` a square in `F_P`? (`0` counts as a square.)
pub fn is_square<const P: u128>(x: Fp<P>) -> bool {
    if x.is_zero() {
        return true;
    }
    fp_pow(x, (P - 1) / 2) == Fp::<P>::one()
}

/// The Hilbert symbol `(a, b)` over `F_P`: `+1` iff `z² = a x² + b y²` has a
/// nontrivial solution. Over a finite field this is identically `+1` for nonzero
/// `a, b` (computed by an honest search, which always succeeds).
pub fn hilbert_symbol<const P: u128>(a: Fp<P>, b: Fp<P>) -> i8 {
    for x in 0..P {
        for y in 0..P {
            for z in 0..P {
                if x == 0 && y == 0 && z == 0 {
                    continue;
                }
                let (fx, fy, fz) = (Fp::<P>(x), Fp::<P>(y), Fp::<P>(z));
                let rhs = a.mul(&fx.mul(&fx)).add(&b.mul(&fy.mul(&fy)));
                if fz.mul(&fz) == rhs {
                    return 1;
                }
            }
        }
    }
    -1
}

/// The classification of a nondegenerate-plus-radical diagonal form over `F_P`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct OddCharType {
    pub p: u128,
    /// Nondegenerate dimension (number of nonzero diagonal entries).
    pub dim: usize,
    /// Radical (null) dimension.
    pub radical_dim: usize,
    /// Discriminant square-class: `true` if `det` of the nondegenerate part is a
    /// square. With `dim`, a complete isometry invariant over a finite field.
    pub disc_is_square: bool,
    /// The Hasse–Witt invariant — always `+1` over a finite field.
    pub hasse: i8,
}

impl OddCharType {
    pub fn display(&self) -> String {
        let d = if self.disc_is_square { "□" } else { "✶" };
        let rad = if self.radical_dim > 0 {
            format!(" ⊗ Λ(F_{}^{})", self.p, self.radical_dim)
        } else {
            String::new()
        };
        format!(
            "F_{}: dim {} disc {} hasse {:+}{}",
            self.p, self.dim, d, self.hasse, rad
        )
    }
}

/// The Hasse invariant `∏_{i<j} (q_i, q_j)` of a diagonal form (nonzero entries
/// only). `None` if the metric is non-diagonal. Always `+1` over a finite field.
pub fn hasse_invariant<const P: u128>(metric: &Metric<Fp<P>>) -> Option<i8> {
    if !metric.b.is_empty() || !metric.a.is_empty() {
        return None;
    }
    let qs: Vec<Fp<P>> = metric.q.iter().copied().filter(|x| !x.is_zero()).collect();
    let mut h = 1i8;
    for i in 0..qs.len() {
        for j in (i + 1)..qs.len() {
            h *= hilbert_symbol(qs[i], qs[j]);
        }
    }
    Some(h)
}

/// The discriminant (product of the nonzero diagonal entries = det of the
/// nondegenerate part). `None` if non-diagonal.
pub fn discriminant<const P: u128>(metric: &Metric<Fp<P>>) -> Option<Fp<P>> {
    if !metric.b.is_empty() || !metric.a.is_empty() {
        return None;
    }
    let mut d = Fp::<P>::one();
    for x in &metric.q {
        if !x.is_zero() {
            d = d.mul(x);
        }
    }
    Some(d)
}

/// Classify a diagonal odd-characteristic form. `None` if non-diagonal.
pub fn classify_oddchar<const P: u128>(metric: &Metric<Fp<P>>) -> Option<OddCharType> {
    if !metric.b.is_empty() || !metric.a.is_empty() {
        return None;
    }
    let dim = metric.q.iter().filter(|x| !x.is_zero()).count();
    let radical_dim = metric.q.len() - dim;
    let disc = discriminant(metric)?;
    Some(OddCharType {
        p: P,
        dim,
        radical_dim,
        disc_is_square: is_square(disc),
        hasse: hasse_invariant(metric)?,
    })
}

/// The odd-characteristic Witt class of a diagonal form: `(dim mod 2, signed
/// discriminant class)`, with `kappa` = nonsquareness of `−1`. `None` if
/// non-diagonal. The signed discriminant `(−1)^{m(m−1)/2}·det` is the genuine
/// Witt invariant; see `witt::WittClassG`.
pub fn oddchar_witt<const P: u128>(metric: &Metric<Fp<P>>) -> Option<WittClassG> {
    if !metric.b.is_empty() || !metric.a.is_empty() {
        return None;
    }
    let mut det = Fp::<P>::one();
    let mut m = 0usize;
    for x in &metric.q {
        if !x.is_zero() {
            det = det.mul(x);
            m += 1;
        }
    }
    // signed discriminant: twist by (−1)^{m(m−1)/2}
    let signed = if (m * (m.wrapping_sub(1)) / 2) & 1 == 1 {
        det.neg()
    } else {
        det
    };
    let kappa = if is_square(Fp::<P>::new(-1)) { 0 } else { 1 };
    Some(WittClassG::OddChar {
        kappa,
        e0: (m & 1) as u8,
        sclass: if is_square(signed) { 0 } else { 1 },
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn diag<const P: u128>(qs: &[u128]) -> Metric<Fp<P>> {
        Metric::diagonal(qs.iter().map(|&x| Fp::<P>(x)).collect())
    }

    #[test]
    fn euler_criterion_matches_brute_force() {
        for p_elt in 0..7u128 {
            let x = Fp::<7>(p_elt);
            let brute = (0..7).any(|y| Fp::<7>(y).mul(&Fp::<7>(y)) == x);
            assert_eq!(is_square::<7>(x), brute, "x = {p_elt} mod 7");
        }
        for p_elt in 0..5u128 {
            let x = Fp::<5>(p_elt);
            let brute = (0..5).any(|y| Fp::<5>(y).mul(&Fp::<5>(y)) == x);
            assert_eq!(is_square::<5>(x), brute, "x = {p_elt} mod 5");
        }
    }

    #[test]
    fn discriminant_distinguishes_planes_over_f3() {
        // <1,1> has disc 1 (square); <1,2> has disc 2 (nonsquare).
        assert!(
            classify_oddchar(&diag::<3>(&[1, 1]))
                .unwrap()
                .disc_is_square
        );
        assert!(
            !classify_oddchar(&diag::<3>(&[1, 2]))
                .unwrap()
                .disc_is_square
        );
    }

    #[test]
    fn hasse_is_trivial_over_finite_fields() {
        // Every Hilbert symbol of nonzero pairs is +1, so every Hasse invariant is.
        for a in 1..5u128 {
            for b in 1..5u128 {
                assert_eq!(hilbert_symbol::<5>(Fp::<5>(a), Fp::<5>(b)), 1);
            }
        }
        for a in 1..7u128 {
            for b in 1..7u128 {
                assert_eq!(hilbert_symbol::<7>(Fp::<7>(a), Fp::<7>(b)), 1);
            }
        }
        assert_eq!(hasse_invariant(&diag::<5>(&[1, 2, 3, 4])).unwrap(), 1);
        assert_eq!(hasse_invariant(&diag::<7>(&[1, 3, 5])).unwrap(), 1);
    }

    // Independent isometry oracle: brute-force search for a congruence
    // M^T diag(d1) M = diag(d2) with M invertible over F_P.
    fn det_small<const P: u128>(m: &[Vec<Fp<P>>]) -> Fp<P> {
        match m.len() {
            1 => m[0][0],
            2 => m[0][0].mul(&m[1][1]).sub(&m[0][1].mul(&m[1][0])),
            _ => unreachable!("only n ≤ 2 in tests"),
        }
    }

    fn is_isometric<const P: u128>(d1: &[Fp<P>], d2: &[Fp<P>]) -> bool {
        let n = d1.len();
        assert_eq!(n, d2.len());
        let mut total = 1u128;
        for _ in 0..(n * n) {
            total *= P;
        }
        for code in 0..total {
            // decode an n×n matrix in base P
            let mut m = vec![vec![Fp::<P>(0); n]; n];
            let mut c = code;
            for row in m.iter_mut() {
                for entry in row.iter_mut() {
                    *entry = Fp::<P>(c % P);
                    c /= P;
                }
            }
            if det_small(&m).is_zero() {
                continue;
            }
            // C = M^T D1 M ; compare to diag(d2)
            let mut ok = true;
            'pair: for i in 0..n {
                for j in 0..n {
                    let mut c_ij = Fp::<P>(0);
                    for k in 0..n {
                        c_ij = c_ij.add(&m[k][i].mul(&d1[k]).mul(&m[k][j]));
                    }
                    let want = if i == j { d2[i] } else { Fp::<P>(0) };
                    if c_ij != want {
                        ok = false;
                        break 'pair;
                    }
                }
            }
            if ok {
                return true;
            }
        }
        false
    }

    #[test]
    fn dim_plus_disc_is_complete_over_finite_fields() {
        // The odd-char analogue of Arf-completeness: two nondegenerate forms are
        // isometric IFF (dim, disc-class) agree. Verified independently against a
        // brute-force congruence search, for dims 1 and 2 over F_3 and F_5.
        fn check<const P: u128>(dim: usize) {
            // all diagonal forms with nonzero entries
            let mut forms: Vec<Vec<Fp<P>>> = vec![vec![]];
            for _ in 0..dim {
                let mut next = vec![];
                for f in &forms {
                    for e in 1..P {
                        let mut g = f.clone();
                        g.push(Fp::<P>(e));
                        next.push(g);
                    }
                }
                forms = next;
            }
            for a in &forms {
                for b in &forms {
                    let disc_a = a.iter().fold(Fp::<P>::one(), |acc, x| acc.mul(x));
                    let disc_b = b.iter().fold(Fp::<P>::one(), |acc, x| acc.mul(x));
                    let same_class = is_square::<P>(disc_a) == is_square::<P>(disc_b);
                    assert_eq!(is_isometric::<P>(a, b), same_class, "P={P} a={a:?} b={b:?}");
                }
            }
        }
        check::<3>(1);
        check::<3>(2);
        check::<5>(1);
        check::<5>(2);
    }

    #[test]
    fn oddchar_witt_is_a_homomorphism() {
        // oddchar_witt(A ⊥ B) = oddchar_witt(A) + oddchar_witt(B): the abstract
        // group law agrees with actual orthogonal sums of forms.
        let forms = [
            diag::<3>(&[1]),
            diag::<3>(&[2]),
            diag::<3>(&[1, 1]),
            diag::<3>(&[1, 2]),
        ];
        for a in &forms {
            for b in &forms {
                let sum = a.direct_sum(b);
                assert_eq!(
                    oddchar_witt(&sum).unwrap(),
                    oddchar_witt(a).unwrap().add(&oddchar_witt(b).unwrap()),
                    "homomorphism failed"
                );
            }
        }
    }

    #[test]
    fn witt_group_is_z4_when_minus_one_nonsquare() {
        // F_3: −1 = 2 is a nonsquare (q ≡ 3 mod 4) ⇒ W(F_3) ≅ ℤ/4.
        let g = oddchar_witt(&diag::<3>(&[1])).unwrap();
        let id = WittClassG::oddchar_zero(1);
        let g2 = g.add(&g);
        let g3 = g2.add(&g);
        let g4 = g3.add(&g);
        assert_ne!(g, id);
        assert_ne!(g2, id); // order > 2 ⇒ not (ℤ/2)²
        assert_ne!(g3, id);
        assert_eq!(g4, id); // order exactly 4 ⇒ ℤ/4
    }

    #[test]
    fn witt_group_is_z2xz2_when_minus_one_square() {
        // F_5: −1 = 4 is a square (q ≡ 1 mod 4) ⇒ W(F_5) ≅ ℤ/2 × ℤ/2.
        let id = WittClassG::oddchar_zero(0);
        let g = oddchar_witt(&diag::<5>(&[1])).unwrap(); // 1 is a square
        let h = oddchar_witt(&diag::<5>(&[2])).unwrap(); // 2 is a nonsquare
                                                         // every nonidentity element has order 2 (exponent 2)
        assert_eq!(g.add(&g), id);
        assert_eq!(h.add(&h), id);
        let gh = g.add(&h);
        assert_eq!(gh.add(&gh), id);
        // the four elements are distinct ⇒ the full Klein four-group
        let elems = [id, g, h, gh];
        for i in 0..4 {
            for j in (i + 1)..4 {
                assert_ne!(elems[i], elems[j], "elements {i},{j} coincide");
            }
        }
    }
}
