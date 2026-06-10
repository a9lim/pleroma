//! A quick tour of pleroma's verified core. Run with:
//!   cargo run --example tour

use pleroma::clifford::{CliffordAlgebra, Metric};
use pleroma::forms::classify_surreal;
use pleroma::forms::WittClass;
use pleroma::forms::{dickson_matrix, dickson_of_versor};
use pleroma::games::{Game, GameExterior};
use pleroma::scalar::Surcomplex;
use pleroma::scalar::Surreal;
use pleroma::scalar::{nim_solve_artin_schreier, nim_sqrt, nim_trace, Nimber};
use pleroma::scalar::{Integer, Rational, Scalar};
use std::collections::BTreeMap;

fn rule(title: &str) {
    println!("\n── {title} ──");
}

fn main() {
    rule("nimbers On₂ — char 2, the non-commutative Clifford case");
    // b[(0,1)] = *1  ⇒  e0 e1 + e1 e0 = *1 ≠ 0  ⇒  non-commutative
    let mut b = BTreeMap::new();
    b.insert((0usize, 1usize), Nimber(1));
    let alg = CliffordAlgebra::new(2, Metric::new(vec![Nimber(2), Nimber(3)], b));
    let (e0, e1) = (alg.gen(0), alg.gen(1));
    println!("  e0 e1      = {}", alg.mul(&e0, &e1).display());
    println!("  e1 e0      = {}", alg.mul(&e1, &e0).display());
    println!(
        "  {{e0,e1}}   = {}   (the anticommutator b[(0,1)] = *1)",
        alg.add(&alg.mul(&e0, &e1), &alg.mul(&e1, &e0)).display()
    );
    println!(
        "  e0²        = {}   (a nimber square, not ±1)",
        alg.mul(&e0, &e0).display()
    );

    rule("Grassmann — fully null metric, nilpotent generators");
    let g = CliffordAlgebra::new(3, Metric::<Rational>::grassmann(3));
    println!("  e0²        = {}", g.mul(&g.gen(0), &g.gen(0)).display());
    println!("  e0 e1      = {}", g.mul(&g.gen(0), &g.gen(1)).display());
    println!(
        "  e1 e0      = {}   (antisymmetric)",
        g.mul(&g.gen(1), &g.gen(0)).display()
    );

    rule("surreals — a Clifford metric with NO finite entries");
    // e0² = ω (infinite), e1² = ε = ω⁻¹ (infinitesimal), orthogonal.
    let s = CliffordAlgebra::new(
        2,
        Metric::diagonal(vec![Surreal::omega(), Surreal::epsilon()]),
    );
    let e0e1 = s.mul(&s.gen(0), &s.gen(1));
    println!("  e0²        = {}", s.mul(&s.gen(0), &s.gen(0)).display());
    println!("  e1²        = {}", s.mul(&s.gen(1), &s.gen(1)).display());
    println!(
        "  (e0 e1)²   = {}   (= -(ω·ε) = -1, a unit bivector)",
        s.mul(&e0e1, &e0e1).display()
    );

    rule("surreal arithmetic — recursive exponents");
    let w = Surreal::omega();
    println!("  ω·ε        = {:?}", w.mul(&Surreal::epsilon()));
    println!(
        "  (ω+1)(ω-1) = {:?}",
        w.add(&Surreal::from_int(1))
            .mul(&w.sub(&Surreal::from_int(1)))
    );
    println!("  √ω squared = {:?}", {
        let r = Surreal::omega_pow(Surreal::from_rational(Rational::new(1, 2)));
        r.mul(&r)
    });
    println!("  ω^ω        = {:?}", Surreal::omega_pow(Surreal::omega()));

    rule("surcomplex — why it only works over the surreals");
    type NC = Surcomplex<Nimber>;
    let one_plus_i = NC::new(Nimber(1), Nimber(1));
    println!(
        "  over On₂:  i²        = {:?}   (= -1 = 1 in char 2)",
        NC::i().mul(&NC::i())
    );
    println!(
        "  over On₂:  (1+i)²    = {:?}   (nonzero nilpotent ⇒ not a field)",
        one_plus_i.mul(&one_plus_i)
    );
    type SC = Surcomplex<Surreal>;
    let z = SC::new(Surreal::omega(), Surreal::from_int(1)); // ω + i
    println!(
        "  over No:   (ω+i)(ω-i) = {:?}   (= ω²+1, a genuine norm)",
        z.mul(&z.conj())
    );

    rule("char-0 classifier — Cl(p,q) as a matrix algebra (companion to Arf)");
    let cl = |qs: &[i128]| {
        let q = qs.iter().map(|&x| Surreal::from_int(x)).collect();
        classify_surreal(&Metric::diagonal(q)).unwrap().display()
    };
    println!("  Cl(0,2) = {}   (the quaternions)", cl(&[-1, -1]));
    println!("  Cl(1,3) = {}   (spacetime)", cl(&[1, -1, -1, -1]));
    println!("  Cl(3,1) = {}   (≠ Cl(1,3)!)", cl(&[1, 1, 1, -1]));
    println!("  Cl(4,1) = {}   (conformal GA)", cl(&[1, 1, 1, 1, -1]));

    rule("even subalgebra + graded tensor product");
    let cl30 = CliffordAlgebra::new(3, Metric::diagonal(vec![Surreal::from_int(1); 3]));
    let even = cl30.even_subalgebra().unwrap();
    println!(
        "  Cl(3,0)⁰         = {}   (≅ Cl(0,2))",
        classify_surreal(&even.metric).unwrap().display()
    );
    let l = CliffordAlgebra::new(1, Metric::diagonal(vec![Surreal::from_int(1)]));
    let r = CliffordAlgebra::new(1, Metric::diagonal(vec![Surreal::from_int(-1)]));
    let t = l.graded_tensor(&r);
    println!(
        "  Cl(1,0) ⊗̂ Cl(0,1) = {}   (≅ Cl(1,1))",
        classify_surreal(&t.metric).unwrap().display()
    );

    rule("general bilinear form — the in-order contraction `a` deforms the product");
    let mut a = BTreeMap::new();
    a.insert((0usize, 1usize), Surreal::from_int(5));
    let d = CliffordAlgebra::new(
        2,
        Metric::general(vec![Surreal::from_int(1); 2], BTreeMap::new(), a),
    );
    let e0e1 = d.mul(&d.gen(0), &d.gen(1));
    println!("  e0 e1 = {}   (= e0∧e1 + 5)", e0e1.display());

    rule("Artin–Schreier ↔ Arf — one field trace, two roles");
    println!(
        "  √*2 = *{}   (since (√*2)² = *{})",
        nim_sqrt(2),
        nim_mul_sq(nim_sqrt(2))
    );
    for c in 0u128..4 {
        let tr = nim_trace(c, 2);
        match nim_solve_artin_schreier(c, 2) {
            Some(y) => println!("  y²+y=*{c} in F₄: Tr=*{tr} ⇒ y=*{y}"),
            None => println!("  y²+y=*{c} in F₄: Tr=*{tr} ⇒ no solution"),
        }
    }

    rule("Witt group (ℤ/2) + Dickson invariant (the char-2 determinant)");
    let mut b = BTreeMap::new();
    b.insert((0usize, 1usize), Nimber(1));
    let aplane = Metric::new(vec![Nimber(1), Nimber(1)], b);
    let wa = WittClass::try_from_metric(&aplane).expect("anisotropic plane is nonsingular");
    println!(
        "  w(A) = {} ;  w(A)+w(A) = {}",
        wa.display(),
        wa.try_add(&wa).expect("same finite char-2 field").display()
    );
    println!(
        "  Dickson(swap) = {}   (reflection)",
        dickson_matrix(&[vec![0, 1], vec![1, 0]])
    );
    let nb = CliffordAlgebra::new(2, aplane);
    let rotor = nb.mul(&nb.gen(0), &nb.gen(1));
    println!(
        "  Dickson(versor e0e1) = {:?}   (a rotor ⇒ SO)",
        dickson_of_versor(&nb, &rotor)
    );

    rule("exterior algebra of the GAME group — lives where Clifford can't");
    let ext = GameExterior::new(vec![Game::star(), Game::up()]);
    println!(
        "  generators ⋆,↑ are numbers? {} {}",
        ext.game(0).is_number(),
        ext.game(1).is_number()
    );
    let g0g1 = ext.wedge(&ext.generator(0), &ext.generator(1));
    println!(
        "  g0 ∧ g1 = {}   (nonzero grade-2, with game relations imposed)",
        g0g1.display()
    );
    println!(
        "  2·(g0 ∧ g1) = 0 ? {}   (the relation 2⋆=0 propagates)",
        ext.is_zero(&ext.scalar_mul(2, &g0g1))
    );
    let two_g0 = ext.algebra().scalar_mul(&Integer(2), &ext.generator(0));
    println!(
        "  value(2·g0) = ⋆+⋆ = 0 ? {}",
        ext.value_of_grade1(&two_g0).eq(&Game::zero())
    );
}

/// (helper) nim-square, for the tour's √ check.
fn nim_mul_sq(x: u128) -> u128 {
    pleroma::scalar::nim_square(x)
}
