use super::*;
use crate::clifford::Metric;
use crate::forms::WittClassG;
use crate::scalar::{Fp, Fpn, Scalar};

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
        classify_finite_odd(&diag::<3>(&[1, 1]))
            .unwrap()
            .disc_is_square
    );
    assert!(
        !classify_finite_odd(&diag::<3>(&[1, 2]))
            .unwrap()
            .disc_is_square
    );
}

#[test]
fn invalid_moduli_are_rejected() {
    assert!(classify_finite_odd(&diag::<2>(&[1, 1])).is_none());
    assert!(classify_finite_odd(&diag::<9>(&[1, 1])).is_none());
    assert!(std::panic::catch_unwind(|| is_square::<2>(Fp::<2>(1))).is_err());
    assert!(std::panic::catch_unwind(|| is_square::<9>(Fp::<9>(1))).is_err());
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
    assert_eq!(
        hasse_invariant_finite_odd(&diag::<5>(&[1, 2, 3, 4])).unwrap(),
        1
    );
    assert_eq!(
        hasse_invariant_finite_odd(&diag::<7>(&[1, 3, 5])).unwrap(),
        1
    );
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
    // finite_odd_witt(A ⊥ B) = finite_odd_witt(A) + finite_odd_witt(B): the abstract
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
                finite_odd_witt(&sum).unwrap(),
                finite_odd_witt(a)
                    .unwrap()
                    .add(&finite_odd_witt(b).unwrap()),
                "homomorphism failed"
            );
        }
    }
}

#[test]
fn witt_group_is_z4_when_minus_one_nonsquare() {
    // F_3: −1 = 2 is a nonsquare (q ≡ 3 mod 4) ⇒ W(F_3) ≅ ℤ/4.
    let g = finite_odd_witt(&diag::<3>(&[1])).unwrap();
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
    let g = finite_odd_witt(&diag::<5>(&[1])).unwrap(); // 1 is a square
    let h = finite_odd_witt(&diag::<5>(&[2])).unwrap(); // 2 is a nonsquare
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

#[test]
fn extension_fields_use_the_same_trait_path() {
    let f9 = Metric::diagonal(vec![Fpn::<3, 2>::constant(1), Fpn::<3, 2>::generator()]);
    let class = classify_finite_odd(&f9).unwrap();
    assert_eq!(class.field_order, 9);
    assert_eq!(
        finite_odd_witt(&f9)
            .unwrap()
            .add(&WittClassG::oddchar_zero(0)),
        finite_odd_witt(&f9).unwrap()
    );
}
