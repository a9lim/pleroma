"""A tour of ogdoad from Python. Run inside the project venv:

    VIRTUAL_ENV=.venv maturin develop && .venv/bin/python demo.py
"""

import ogdoad as pl


def section(title):
    print(f"\n── {title} ──")


def raises_value_error(fn):
    try:
        fn()
    except ValueError:
        return True
    return False


section("ogham — the expression language over fixed worlds")
print(pl.ogham_eval("nimber 2 q=[*1,*1]", "e0 & e0\n[*1,*2] & [*1,*3]\ne0 . e0"))
print("  bare int rejected in nimber world:", raises_value_error(lambda: pl.ogham_eval("nimber 0", "3")))
print(pl.ogham_eval("polyint", "(5.t + 1)@7\ndeg(t^2 + 1)\ngcd(2.t + 2, 4.t + 4)"))
print(pl.ogham_eval("integer 0", "a := 5; a + 1\nnorm1 := (u, v) ↦ (\n  s := u + v;\n  d := u - v;\n  s.s + d.d\n)\nnorm1@(2, 1)"))
print("  ratfunc pole rejected:", raises_value_error(lambda: pl.ogham_eval("ratfunc5", "(1/(t + 1))@4")))

section("nimbers On₂ — char 2, the non-commutative Clifford case")
# b[(0,1)] = *1  ⇒  e0 e1 + e1 e0 = *1 ≠ 0  ⇒  non-commutative.
A = pl.NimberAlgebra(q=[pl.Nimber(2), pl.Nimber(3)], b={(0, 1): 1})
e0, e1 = A.gen(0), A.gen(1)
print("  e0 e1        =", e0 * e1)
print("  e1 e0        =", e1 * e0)
print("  {e0,e1}      =", e0 * e1 + e1 * e0, "  (the anticommutator)")
print("  e0²          =", e0 ** 2, "  (a nimber square, not ±1)")

section("Grassmann — fully null metric, nilpotent generators")
G = pl.SurrealAlgebra(q=[0, 0, 0])
g0, g1 = G.gen(0), G.gen(1)
print("  g0²          =", g0 ** 2)
print("  g0 ∧ g1      =", g0 & g1)
print("  g0∧g1 == g0 g1:", (g0 & g1) == (g0 * g1))

section("surreals — a Clifford metric with NO finite entries")
# e0² = ω (infinite), e1² = ε = ω⁻¹ (infinitesimal), orthogonal.
S = pl.SurrealAlgebra(q=[pl.omega(), pl.epsilon()])
b = S.gen(0) * S.gen(1)
print("  e0²          =", S.gen(0) ** 2)
print("  e1²          =", S.gen(1) ** 2)
print("  (e0 e1)²     =", b ** 2, "  (= -(ω·ε) = -1, a unit bivector)")

section("surreal arithmetic — recursive exponents")
w = pl.omega()
print("  ω·ε          =", w * pl.epsilon())
print("  (ω+1)(ω-1)   =", (w + 1) * (w - 1))
cnf = pl.surreal(3) * w ** 2 - w + 5
print("  3ω² - ω + 5  =", cnf)
print("  sign/lead/dyadic/monomial:", cnf.sign(), cnf.terms[0],
      pl.Surreal.from_rational(3, 4).as_dyadic(), pl.Surreal.monomial(2, 7))
print("  ω^ω          =", pl.omega_pow(pl.omega()))
print("  ω > 10⁶ ?    ", w > 1_000_000)
print("  0 < ε < 1e-9?",
      pl.surreal(0) < pl.epsilon() < pl.Surreal.from_rational(1, 10**9))

section("surcomplex — only a field over the surreals")
z = pl.Surcomplex(pl.omega(), 1)  # ω + i
print("  (ω+i)(ω-i)   =", z * z.conj(), "  (= ω²+1, a genuine norm)")
# algebraic closure: every represented number has a represented root (ExactRoots)
print("  √(3+4i) = 2+i:", pl.Surcomplex(3, 4).sqrt())
print("  √ω over ℂ    :", pl.Surcomplex(pl.omega(), 0).sqrt(), " √2 declines:",
      pl.Surcomplex(2, 0).sqrt())
# division stays first-class even when the norm a²+b² is a non-monomial surreal
zz = pl.Surcomplex(pl.omega() + 1, 1)
print("  1/(ω+1+i)→6t :", zz.inv_to_terms(6), "  z·z⁻¹ =", zz * zz.inv_to_terms(6))

section("versors — reflections & rotations in Cl(3,0)")
E = pl.SurrealAlgebra(q=[1, 1, 1])
e0, e1, e2 = E.gen(0), E.gen(1), E.gen(2)
x = 3 * e0 + 4 * e1
print("  reflect e0 in ⊥e1     =", e1.reflect(e0), "  (fixed)")
print("  reflect e1 in ⊥e1     =", e1.reflect(e1), "  (negated)")
R = e0 * e1  # rotor = product of two unit vectors
print("  rotor (e0 e1) on x    =", R.sandwich(x))
vc = R.classify_versor()
print("  spinor norm/parity    =", R.spinor_norm(), R.versor_grade_parity(),
      " named:", (vc.spinor_norm, vc.dickson))
print("  norm² preserved       =", x.norm2(), "->", R.sandwich(x).norm2())
print("  ~(e0 e1)  (reversion) =", ~R)
print("  e0 ⌟ (e0∧e1)          =", e0 << (e0 & e1))
print("  dual(e0) in 3D        =", e0.dual(), "  (a bivector)")

section("Arf invariant — the char-2 Clifford classifier (see README.md)")
def arf(qs, bs):
    A = pl.NimberAlgebra(q=[pl.Nimber(x) for x in qs], b={k: 1 for k in bs})
    return pl.arf_nimber(A)
print("  Q = x0·x1        (hyperbolic) :", arf([0, 0], [(0, 1)]))
print("  Q = x0²+x0x1+x1² (anisotropic):", arf([1, 1], [(0, 1)]))
print("  A⊕A ≅ H⊕H (Arf additive)      :", arf([1, 1, 1, 1], [(0, 1), (2, 3)]).o_type)
print("  raw F₂ bitmask Arf           :", pl.arf_f2(2, [True, True], [2, 1]))

section("char-0 classifier — the companion to Arf (Cl(p,q) → matrix algebra)")
def cl(qs):
    return pl.classify_surreal(pl.SurrealAlgebra(q=qs))
print("  Cl(0,2)  =", cl([-1, -1]),       "  (the quaternions ℍ)")
print("  Cl(3,0)  =", cl([1, 1, 1]),      "  (M₂(ℂ))")
print("  Cl(1,3)  =", cl([1, -1, -1, -1]),"  (spacetime, M₂(ℍ))")
print("  Cl(3,1)  =", cl([1, 1, 1, -1]),  "  (≠ Cl(1,3): M₄(ℝ))")
print("  Cl(4,1)  =", cl([1, 1, 1, 1, -1]),"  (conformal GA, M₄(ℂ))")
print("  surcomplex Cl(2,ℂ) =", pl.classify_surcomplex(
    pl.SurcomplexAlgebra([pl.Surcomplex(1), pl.Surcomplex(1)])))

section("even subalgebra + graded tensor product")
cl30 = pl.SurrealAlgebra(q=[1, 1, 1])
print("  Cl(3,0)⁰         =", pl.classify_surreal(cl30.even_subalgebra()), "  (≅ Cl(0,2) = ℍ)")
tens = pl.SurrealAlgebra(q=[1]).graded_tensor(pl.SurrealAlgebra(q=[-1]))
print("  Cl(1,0) ⊗̂ Cl(0,1) =", pl.classify_surreal(tens), "  (≅ Cl(1,1) = M₂(ℝ))")

section("general bilinear form — the in-order contraction `a` deforms the product")
# a[(0,1)] = 5: e0 e1 = e0∧e1 + 5, while the anticommutator {e0,e1}=b stays 0.
D = pl.SurrealAlgebra(q=[1, 1], b=None, a={(0, 1): pl.surreal(5)})
d0, d1 = D.gen(0), D.gen(1)
print("  e0 e1            =", d0 * d1, "  (= e0∧e1 + 5)")
print("  {e0,e1} = b = 0  :", d0 * d1 + d1 * d0)

section("twisted adjoint (Pin) — the correct versor action")
P = pl.SurrealAlgebra(q=[1, 1])
p0, p1 = P.gen(0), P.gen(1)
print("  twisted_sandwich e1 on 3e0+4e1 =", p1.twisted_sandwich(3 * p0 + 4 * p1), " (= reflection)")

section("Artin–Schreier ↔ Arf — the same field trace, two roles")
root = pl.Nimber(pl.nim_sqrt(2))
print("  √*2 in On₂           =", root, " (since (√*2)² =", root * root, ")")
for c in range(4):
    y = pl.nim_solve_artin_schreier(c, 2)
    print(f"  y²+y=*{c} in F₄: Tr=*{pl.nim_trace(c,2)}  ->  "
          + (f"y=*{y}" if y is not None else "no solution"))

section("Witt group (ℤ/2) + Dickson invariant (char-2 determinant)")
A = pl.NimberAlgebra(q=[1, 1], b={(0, 1): 1})  # anisotropic plane
wA = pl.witt_class(A)
print("  w(A) =", wA, "   w(A)+w(A) =", wA + wA, " (A⊕A ≅ H⊕H)")
print("  WittClass constructors:", pl.WittClass.zero(), -pl.WittClass(1), pl.WittClass(1).arf)
print("  WittClass metric constructor:", pl.WittClass.try_from_metric(A) == wA)
print("  Dickson(swap)  =", pl.dickson_matrix([[0, 1], [1, 0]]), " (a reflection)")
print("  Dickson(diag *2,*3 rotation) =", pl.dickson_matrix([[2, 0], [0, 3]]), " (in SO)")

section("exterior algebra of the GAME group — lives where Clifford can't")
# Λ needs only a ℤ-module; the game group is one, even for non-numbers (⋆, ↑).
ext = pl.GameExterior([pl.Game.star(), pl.Game.up()])
g0, g1 = ext.generator(0), ext.generator(1)
print("  generators are non-numbers:", not ext.game(0).is_number(), not ext.game(1).is_number())
cert = ext.relation_search_certificate()
print("  raw free algebra/search cert:", ext.algebra().dim, ext.relation_search_complete(),
      cert.bound, cert.relations)
print("  relation-search constructor:",
      pl.GameExterior.with_relation_search([pl.Game.star(), pl.Game.up()], 1).relation_search_complete())
g0g1 = ext.wedge(g0, g1)
print("  g0 ∧ g1 = -(g1 ∧ g0):", g0g1, "==", -ext.wedge(g1, g0))
print("  2·(g0 ∧ g1) = 0       :", ext.is_zero(ext.scalar_mul(2, g0g1)), " (relation 2⋆=0 propagates)")
print("  value(g0 + g1) = ⋆ + ↑ :", ext.value_of_grade1(g0 + g1))
print("  value(2·g0) = ⋆+⋆ = 0  :", ext.value_of_grade1(2 * g0) == pl.Game.zero())
rel = pl.GameRelation([2, 0])
ext_explicit = pl.GameExterior.with_relations([pl.Game.star(), pl.Game.up()], [rel])
print("  explicit relation 2⋆=0   :", ext_explicit.is_zero(2 * ext_explicit.generator(0)),
      ext_explicit.relations()[0].coeffs, cert.bound)
checked = pl.GameClifford.with_quadratic_data([pl.Game.star(), pl.Game.up()], [rel], [0, 5])
c0, c1 = checked.generator(0), checked.generator(1)
print("  checked Clifford ↑²       :", checked.mul(c1, c1),
      "  2·(⋆↑)=0:", checked.is_zero(checked.scalar_mul(2, checked.mul(c0, c1))))
try:
    pl.GameClifford.with_quadratic_data([pl.Game.star(), pl.Game.up()], [rel], [1, 0])
except ValueError as exc:
    print("  rejects Q(⋆)=1 under 2⋆=0:", "polar pairing" in str(exc))


# ===========================================================================
# The expansion pass: new scalars, GA configurations, deeper invariants
# ===========================================================================

section("Fp — odd characteristic, completing the classification trichotomy")
# char 0: signature → matrix algebra. char 2: Arf. odd char: dim + discriminant.
f3_11 = pl.OddFiniteFieldForm(3, [1, 1])
f3_12 = pl.OddFiniteFieldForm(3, [1, 2])
f5_1234 = pl.OddFiniteFieldForm(5, [1, 2, 3, 4])
print("  F₃ <1,1> :", f3_11.classify())   # disc 1 = square
print("  F₃ <1,2> :", f3_12.classify())   # disc 2 = nonsquare
print("  Hasse always +1 over a finite field:", f5_1234.hasse_invariant())
print("  finite odd helper package     :", f3_12.classify().display(),
      pl.OddFiniteFieldForm(5, [1, 2]).witt_class(), pl.hilbert_symbol(5, 2, 3),
      pl.OddFiniteFieldForm(7, []).is_square(2))
f9 = pl.OddFiniteFieldForm(3, [1, 1], degree=2)
f9_stair = f9.e_staircase()
print("  F₉ BW/e-staircase helpers    :", f9.bw_class(),
      (f9_stair.e0, f9_stair.e1, f9_stair.e2, f9_stair.stabilizes_at))
f9_class = f9.classify_unified()
print("  F₉ <1,1> via OddFiniteFieldForm:", f9.classify(), "class=", f9_class.kind,
      f9_class.display(), "W=", f9.witt_class())
print("  finite form isometric_to     :", f9.isometric_to(f9))
# the odd-char Witt group: ℤ/4 when −1 is a nonsquare (F₃), ℤ/2×ℤ/2 when it is (F₅)
g3 = pl.OddFiniteFieldForm(3, [1]).witt_class(); zero3 = pl.OddFiniteFieldForm(3, []).witt_class()
print("  W(F₃) is ℤ/4 :", g3 + g3 != zero3, "and", g3 + g3 + g3 + g3 == zero3)
g5, h5 = pl.OddFiniteFieldForm(5, [1]).witt_class(), pl.OddFiniteFieldForm(5, [2]).witt_class()
zero5 = pl.OddFiniteFieldForm(5, []).witt_class()
print("  W(F₅) is ℤ/2×ℤ/2 :", g5 + g5 == zero5 and h5 + h5 == zero5)
f2a = pl.Char2FiniteFieldForm([1, 1], {(0, 1): 1})
f8h = pl.Char2FiniteFieldForm([0, 0], {(0, 1): 3}, degree=3)
print("  F₂ char-2 anisotropic plane :", f2a.classify(), f2a.classify_unified().kind, f2a.bw_class())
print("  F₈ char-2 hyperbolic plane  :", f8h.classify(), f8h.is_isometric(
      pl.Char2FiniteFieldForm([0, 0], {(0, 1): 1}, degree=3)))
n2a = pl.NimberAlgebra(q=[1, 1], b={(0, 1): 1})
print("  char-2 form methods          :", pl.arf_nimber(n2a).arf,
      f8h.classify().arf,
      f8h.isometric_to(pl.Char2FiniteFieldForm([0, 0], {(0, 1): 1}, degree=3)),
      f8h.isometric_to(pl.Char2FiniteFieldForm([0, 0], {(0, 1): 1}, degree=3)))
fp = pl.Fp5(2) * pl.Fp5(3)
f8x = pl.F8.generator()
print("  F₅ scalar 2·3 and inverse   :", fp, fp.inv())
print("  F₈ generator minpoly/order  :", f8x.min_poly(), f8x.multiplicative_order())
print("  F₈ from_coeffs/into_coeffs  :", pl.F8.from_coeffs([0, 1, 0]) == f8x,
      f8x.into_coeffs())
print("  F₈ trait helpers             :", pl.F8.ext_degree(), f8x.multiplicative_order(),
      f8x.min_poly_monic(), f8x.relative_trace_over(3, 1), f8x.relative_norm_over(3, 1))
print("  F₅/F₈ roots and reduction   :", pl.Fp5(4).sqrt(), f8x.is_square(), (f8x * f8x).sqrt(),
      pl.F8.reduction_polynomial_kind(), pl.F8.reduction_rule())
print("  finite-field root methods    :", pl.Fp5(4).is_square(), pl.Fp5(4).sqrt(),
      pl.F8.from_coeffs([0, 0, 1]).sqrt())
print("  finite-field validators      :", pl.Fp5.assert_prime_modulus() is None,
      pl.F8.assert_supported_field() is None)
print("  F₈ Frobenius²/group factors :", f8x.frobenius_iter(2), pl.F8.group_order_factors())
print("  Nimber trait helpers         :", pl.Nimber.ext_degree(), pl.Nimber(2).multiplicative_order(),
      pl.Nimber(2).relative_trace_over(4, 1), pl.Nimber(2).relative_norm_over(4, 2))
print("  scalar operators              :", pl.Nimber(2) + 3, pl.Fp5(2) * 3,
      -pl.F8.generator(), pl.Rational(3, 4) - pl.Rational(1, 6))
print("  F₅ disc / F₄ AS class        :", pl.OddFiniteFieldForm(5, [1, 2]).discriminant(),
      pl.artin_schreier_class_finite(2, degree=2))
F3Cl = pl.Fp3Algebra([1, 2])
print("  Cl_F3 <1,2> bivector²       :", (F3Cl.gen(0) * F3Cl.gen(1)) ** 2)
print("  classify Cl_F3 algebra       :", pl.classify_finite_algebra(F3Cl),
      pl.classify_finite_algebra_unified(F3Cl).display(), pl.bw_class_finite_algebra(F3Cl))
F8DP = pl.F8DividedPowerAlgebra(1)
print("  Γ_F8: x·γ₁                  :", F8DP.scalar(f8x) * F8DP.gen(0))

section("Oz — omnific integers: an exterior algebra over a transfinite ring")
Oz = pl.OmnificAlgebra(q=[0, 0, 0])  # Grassmann over Oz
e0, e1, e2 = Oz.gen(0), Oz.gen(1), Oz.gen(2)
w = pl.omnific_omega()
print("  e0² = 0 (nilpotent):", (e0 * e0).is_zero())
print("  (ω·e0) ∧ e1 ∧ e2   :", (w * e0) & e1 & e2, "  (ω-scale coefficient)")
print("  Oz validator ω / ε :", pl.is_omnific_integer(pl.omega()), "/", pl.is_omnific_integer(pl.epsilon()))
print("  ω is not a unit (1/ω=ε ∉ Oz):", end=" ")
try:
    w.inv(); print("?!")
except ValueError:
    print("correctly rejected")

section("ordinal nimbers On₂ — the char-2 mirror of the surreals")
omega = pl.Ordinal.omega()
print("  ω ⊕ ω        =", omega.nim_add(omega), "   (self-inverse)")
print("  ω·2 ⊕ ω      =", pl.Ordinal.monomial(pl.Ordinal(1), 2).nim_add(omega))
print("  ω < ω²       :", omega < pl.Ordinal.omega_pow(pl.Ordinal(2)))
print("  ordinal order ω < ω²:", omega < pl.Ordinal.omega_pow(pl.Ordinal(2)))
print("  2 ⊗ 2 = *3   :", pl.Ordinal(2).nim_mul(pl.Ordinal(2)))
# nim-multiplication: implemented below ω^ω via the current DiMuro/Conway
# degree-3 tower. The old φ_{ω+1} (<ω³) case is the first layer.
print("  ω ⊗ ω        =", omega.nim_mul(omega), "  (just polynomial mult)")
omega_sq = omega.nim_mul(omega)
print("  ω ⊗ ω ⊗ ω    =", omega_sq.nim_mul(omega), "  (the headline: ω³ = 2)")
# (ω + 1)³ in characteristic 2 = ω³ + ω² + ω + 1 = 2 + ω² + ω + 1 = ω² + ω + *3
w1 = omega.nim_add(pl.Ordinal(1))
print("  (ω+1)³       =", w1.nim_mul(w1).nim_mul(w1), "  (= ω² + ω + nim_add(2,1))")
omega3_cell = pl.Ordinal.from_omega3_coeffs([3, 1, 1])
print("  <ω³ coeffs   :", omega3_cell, "↔", omega3_cell.as_below_omega3())
# finite CNF exponents above 3 now stay inside the implemented tower
print("  ω³ ⊗ ω       :", pl.Ordinal.omega_pow(pl.Ordinal(3)).nim_mul(omega))
print("  ω^ω staged   :", pl.Ordinal.omega_pow(omega).nim_mul(omega))
print("  scalar ops/inv in On₂         :", omega + pl.Ordinal(1), pl.Ordinal(2) * pl.Ordinal(3),
      pl.Ordinal(3).inv())
O = pl.OrdinalAlgebra([omega, omega.nim_mul(omega)])
print("  Cl_On₂ e0²/e1²:", O.gen(0) * O.gen(0), O.gen(1) * O.gen(1))
OH = pl.OrdinalAlgebra([0, 0], b={(0, 1): 1})
print("  finite ordinal Arf/Witt/BW:", pl.arf_ordinal_finite(OH), pl.ordinal_witt(OH),
      pl.bw_class_ordinal(OH))
print("  ordinal finite isometry      :", pl.isometric_ordinal_finite(OH, OH))

section("coin-turning games — the game-built nim product")
print("  singleton companions at n=4:", pl.coin_companions("singleton", 4))
print("  exact companion fns match   :", pl.singleton_companions(4) == pl.coin_companions("singleton", 4),
      pl.turtles_companions(4) == pl.coin_companions("turtles", 4))
print("  singleton Grundy g(n):", [pl.coin_turning_grundy("singleton", n) for n in range(6)])
print("  turtles Grundy g(n)  :", [pl.coin_turning_grundy("turtles", n) for n in range(6)])
print("  Tartan singleton×singleton (2,4):",
      pl.coin_turning_tartan_grundy("singleton", "singleton", 2, 4),
      "= nim_mul_mex", pl.nim_mul_mex(2, 4))
print("  Tartan singleton×turtles (3,2):",
      pl.coin_turning_tartan_grundy("singleton", "turtles", 3, 2),
      "=", pl.nim_mul_mex(3, pl.coin_turning_grundy("turtles", 2)))
singleton_cb = lambda n: [1 << i for i in range(n)]
print("  callback g₁D/Tartan          :", pl.grundy_1d(singleton_cb, 5),
      pl.tartan_grundy(singleton_cb, singleton_cb, 2, 4))
print("  callback Grundy path         :", pl.grundy(3, lambda h: list(range(h))))

section("outermorphisms + determinant — Grassmann's def, char-faithful")
R = pl.SurrealAlgebra(q=[1, 1])
lin = [[2, 3], [1, 4]]  # columns f(e_i), i.e. row matrix [[2,1],[3,4]]
lin_map = pl.SurrealLinearMap.from_columns(lin)
lin_inv = R.inverse_outermorphism(lin_map)
print("  det [[2,1],[3,4]] = 2·4−1·3 =", R.determinant(lin_map))
print("  tr Λ²[[2,1],[3,4]]           =", R.exterior_power_trace(lin_map, 2))
print("  inverse outermorphism returns:",
      R.apply_outermorphism(lin_inv, R.apply_outermorphism(lin_map, R.gen(0))) == R.gen(0))
skewed = pl.RationalAlgebra([0, 0], b={(0, 1): pl.Rational(2)})
print("  Gram/diagonalized H over ℚ   :", skewed.gram(), skewed.diagonalize().gram())
general_alg = pl.RationalAlgebra.general([0, 0], a={(0, 1): pl.Rational(5)})
grass_alg = pl.SurrealAlgebra.grassmann(2)
print("  metric Rust helpers          :", R.q(), R.q_val(9), R.is_orthogonal(),
      general_alg.has_upper(), (grass_alg.gen(0) * grass_alg.gen(0)).is_zero())
mapped_R = R.map(lambda x: -x)
print("  metric map helper            :", mapped_R.q()[0] == -pl.surreal(1),
      mapped_R.q()[1] == -pl.surreal(1))
N = pl.NimberAlgebra(q=[1, 1])
nim_lin = pl.NimberLinearMap.from_columns([[2, 1], [3, 1]])
print("  the char-2 determinant (= permanent):", N.determinant(nim_lin))
f8_frob_map = pl.frobenius_linear_map(2, 3)
f16_frob_map = pl.nimber_subfield_frobenius_linear_map(4)
print("  Frobenius spectra F₈/F₁₆     :",
      pl.Fp2Algebra([0, 0, 0]).char_poly(f8_frob_map),
      pl.Fp2Algebra([0, 0, 0, 0]).determinant(f16_frob_map))
print("  Frobenius LinearMap constructors:",
      pl.galois_linear_map(2, 3, 1) == f8_frob_map,
      f16_frob_map.n == 4,
      R.apply_outermorphism(lin_map, R.gen(0)) == lin_map.image(R, 0))
id_map = pl.SurrealLinearMap.identity(R.dim)
print("  LinearMap Rust helpers       :", lin_map.image(R, 0),
      lin_map.compose(id_map) == lin_map,
      R.determinant(id_map) == pl.surreal(1))
print("  MV Rust-name helpers         :", R.gen(0).grade_part(1) == R.gen(0),
      R.gen(0).versor_inverse() == R.gen(0),
      R.scalar(1).multivector_inverse() == R.scalar(1))
TT = R.tensor_square()
print("  tensor embeddings e0/e1      :", TT.embed_first(R.gen(0)), TT.embed_second(R.gen(1), R.dim))

section("exterior Hopf algebra — antipode = grade involution (not reversion-twist)")
H = pl.SurrealAlgebra(q=[0, 0])  # exterior algebra
b = H.gen(0) & H.gen(1)
print("  Δ(e0) primitive (lives in Cl⊗̂Cl):", H.gen(0).coproduct())
print("  S(e0) = −e0          :", H.gen(0).antipode() == -H.gen(0))
print("  S(e0∧e1) = +e0∧e1    :", b.antipode() == b, " (grade 2: (−1)²=+1)")

section("concrete spinor modules — the classification as matrices on spinors")
S = pl.SurrealAlgebra(q=[1, 1, 1])  # Cl(3,0) ≅ M₂(ℂ)
spin = S.spinor_rep()
idem, basis, M = spin.idempotent, spin.basis, spin.gen_matrices
print("  Cl(3,0) minimal ideal real-dim:", len(basis), "(= 2 cols × ℂ)")
print("  spinor metadata says regular?:", spin.is_left_regular,
      " basis:", len(spin.basis))


def _matmul(a, b):
    n = len(a)
    return [[sum((a[i][k] * b[k][j] for k in range(n)), pl.surreal(0))
             for j in range(n)] for i in range(n)]


M0sq = _matmul(M[0], M[0])
holds = all(M0sq[i][j] == pl.surreal(1 if i == j else 0)
            for i in range(len(basis)) for j in range(len(basis)))
print("  M0² = q0·I (the Clifford relation, on the spinor matrices):", holds)
H2 = pl.NimberAlgebra(q=[0, 0], b={(0, 1): 1})
basis2 = H2.spinor_rep().basis
print("  char-2 hyperbolic plane spinor ideal dim:", len(basis2), "(idempotent e0e1)")
Large = pl.RationalAlgebra(q=[1] * 11)
lazy_one = Large.scalar(1)
lazy_e0 = Large.apply_generator(0, lazy_one)
lazy_back = Large.apply_vector([1] + [0] * 10, lazy_e0)
lazy_rep = Large.lazy_spinor_rep()
print("  lazy spinor beyond matrix cap :", lazy_back == lazy_one,
      Large.apply_generator(0, lazy_one) == lazy_e0,
      Large.apply_vector([1] + [0] * 10, lazy_e0) == lazy_back,
      lazy_rep.apply_generator(0, lazy_one) == lazy_e0)

section("conformal GA over the surreals — exact ∞ and infinitesimal radii")
cga = pl.SurrealCga(2)
p_inf = cga.up([pl.omega(), 0])             # a point at ω-scale
print("  up(ω, 0) still null:", cga.inner(p_inf, p_inf) == pl.surreal(0))
eps = pl.epsilon()
sph = cga.sphere([0, 0], eps * eps)         # a sphere of radius ε
on, off = cga.up([eps, 0]), cga.up([2 * eps, 0])
print("  ε-sphere: ε-point on, 2ε-point off:",
      cga.inner(on, sph) == pl.surreal(0), cga.inner(off, sph) != pl.surreal(0))

section("projective GA — exact nilpotent motor (no transcendentals)")
P = pl.SurrealAlgebra.pga(2)                # Cl(2,0,1), e0 the ideal direction
motor = (P.gen(0) & P.gen(1)).exp_nilpotent()  # B² = 0 ⇒ exp(B) = 1 + B
print("  exp(e0∧e1) = 1 + B:", motor == P.scalar(1) + (P.gen(0) & P.gen(1)))
print("  it translates e1 ↦ e1 + 2e0:", motor.sandwich(P.gen(1)) == P.gen(1) + 2 * P.gen(0))

section("non-Archimedean Springer decomposition (surreal)")
form = pl.SurrealAlgebra(q=[pl.omega(), pl.epsilon(), 1, -1])
sp = pl.springer_decompose(form)
print("  valuation filtration of ⟨ω, ε, 1, −1⟩:")
print("   ", sp)
print("  top residue layer:", sp.graded[0].valuation, sp.graded[0].signature)
print("  (W(No)=W(ℝ)=ℤ — the value group is 2-divisible; the filtration is the novelty)")

section("eₙ staircase — discriminant & Hasse as one filtration")
# Over a finite field I²=0, so the staircase is (e₀, e₁) and e₂ is trivial.
stair = pl.OddFiniteFieldForm(5, [1, 2, 3]).e_staircase()
print(f"  ⟨1,2,3⟩/F5: e0={stair.e0} (dim) e1={stair.e1} (disc) e2={stair.e2:+} (Hasse), I^{stair.stabilizes_at}=0")
# Over ℝ the tower is infinite: eₙ reads the 2-adic expansion of the signature.
print("  ⟨1,1,1,1⟩/ℝ (sig 4): eₙ for n=0..3 =", [pl.e_real(4, n) for n in range(4)])
print("  numeric invariants of F5:",
      "level", pl.level(5),
      "pythagoras", pl.pythagoras_number(5),
      "u", pl.u_invariant(5))
print("  WittClassG constructors:",
      pl.WittClassG.char0(3, 1), pl.WittClassG.oddchar_one(5, 0) * pl.WittClassG.oddchar_zero(5, 0),
      pl.WittClassG.char2(1).arf())
print("  WittClassG operators:",
      pl.WittClassG.try_char2_from_metric(A).arf(),
      (pl.WittClassG.char0(2, 0) + pl.WittClassG.char0(1, 0)).signature(),
      (pl.WittClassG.oddchar_one(5, 0) * pl.WittClassG.oddchar_one(5, 0)).kind())
print("  2 is a sum of two squares in F3:", pl.is_sum_of_n_squares(3, 2, 2))

section("symplectic and Hermitian forms — form + involution siblings")
print("  2 hyperbolic symplectic planes:", pl.SymplecticForm.hyperbolic(2).classify())
print("  char-2 alternating [[0,1],[1,0]]:", pl.classify_symplectic_nimber([[0, 1], [1, 0]]))
H = pl.HermitianForm.from_gram([[pl.Surcomplex(2, 0), pl.Surcomplex(0, 1)],
                                [pl.Surcomplex(0, -1), pl.Surcomplex(2, 0)]])
print("  Hermitian [[2,i],[-i,2]]:", H.signature(), "diagonal", H.diagonalize())
print("  diagonal Hermitian ⟨1,-1,0⟩:", pl.HermitianForm.diagonal([1, -1, 0]).signature())
print("  form Rust constructors       :",
      pl.SymplecticForm.from_gram([[0, 1], [-1, 0]]).classify().planes(),
      (lambda sig: (sig.pos, sig.neg, sig.radical))(
          pl.HermitianForm.from_gram([[pl.Surcomplex(1, 0)]]).signature()))

section("p-adic Hilbert symbol + Hasse–Minkowski over Q")
print("  (−1,−1)_2 =", pl.hilbert_symbol_qp(-1, -1, 2),
      " — Hamilton's quaternions ramify at 2 (finite fields can't show this)")
print("  ε₂(⟨1,1,1⟩), (−1,−1)∞:", pl.hasse_at_place([1, 1, 1], pl.RationalPlace.prime(2)),
      pl.hilbert_symbol_at(-1, -1))
print("  checked Qp helpers reject 9 :", raises_value_error(lambda: pl.hilbert_symbol_qp(2, 3, 9)),
      raises_value_error(lambda: pl.is_square_qp(2, 9)), pl.is_isotropic_q([1, 1, 1]))
print("  Qp isotropy checked helpers :",
      pl.is_isotropic_at_p([1, 1, -1], 2),
      raises_value_error(lambda: pl.is_isotropic_at_p([1, 1, 1], 9)))
for f in ([1, 1, 1], [1, 1, -1], [1, 1, -3], [1, 1, 1, 1, -1]):
    print(f"  ⟨{','.join(map(str, f))}⟩ isotropic over Q:", pl.is_isotropic_q(f))

section("loopy impartial games — Side values with a certificate")
values, cert = pl.loopy_nim_values_certified([[1], [0], []])
print("  2-cycle plus terminal:", values, cert, "outcomes", cert.outcomes,
      "recovery:", cert.recovery_condition_holds)

section("Brauer–Wall group — BW(ℝ)=ℤ/8 is the Bott clock")
# walk ⟨−1⟩⊗̂…⊗̂⟨−1⟩: the Bott index cycles mod 8.
g = pl.bw_class_real(pl.SurrealAlgebra(q=[-1]))
walk, cur = [], g
for _ in range(8):
    walk.append(repr(cur)); cur = cur + g
print("  [Cl⟨−1⟩]ⁿ for n=1..8:", " ".join(w.replace("Real(", "").rstrip(")") for w in walk))
print("  BW constructors/zero_like:", pl.BrauerWallClass.real(9), g.zero_like(),
      pl.BrauerWallClass.char2(1) + pl.BrauerWallClass.char2(1))
print("  BW(F_3) of ⟨1⟩:", pl.OddFiniteFieldForm(3, [1]).bw_class(), "(order-4 graded part ≅ W(F_3))")
A2 = pl.NimberAlgebra(q=[1, 1], b={(0, 1): 1})
print("  BW(F_2^m) anisotropic nimber plane:", pl.bw_class_nimber(A2), "(Z/2 Arf class)")


# ===========================================================================
# Arc IV: the CGT/surreal core, forms foundations, and GA depth
# ===========================================================================

section("partizan canonical form — Conway's simplicity theorem")
# Value-preserving reduction: G − G = 0 for any G, so its canonical form is {|}.
up = pl.Game.up()
print("  ↑ − ↑ canonical = 0   :", (up - up).canonical_string() == pl.Game.zero().canonical_string())
print("  ↑ is already canonical:", up.is_canonical())
raw_zero = up - up
print("  raw tree vs value key  :", raw_zero.structural_eq(pl.Game.zero()),
      raw_zero.canonical_string() == pl.Game.zero().canonical_string())
# A messy sum reduces to a single simple game (here ↑ + ⋆ + ⋆ = ↑, since ⋆+⋆=0).
messy = pl.Game.up() + pl.Game.star() + pl.Game.star()
print("  canonical(↑+⋆+⋆) == ↑ :", messy.canonical_string() == pl.Game.up().canonical_string())
print("  operator add/neg/eq checks   :",
      up + pl.Game.star() == up + pl.Game.star(), -up == -up)

section("game ↔ surreal bridge — numbers carry a surreal value")
half = pl.Game.from_surreal(pl.Surreal.simplest_between(0, 1))  # ½ = {0|1}
print("  Game.from_surreal(½)  =", half, " value =", half.number_value(), " birthday =", half.birthday())
print("  ⋆, ↑, ±1 are non-numbers (no value):",
      pl.Game.star().number_value(), pl.Game.up().number_value(), pl.Game.switch(1, -1).number_value())
print("  simplest_between(⅓,⅔) = ½:",
      pl.Surreal.simplest_between(pl.Surreal.from_rational(1, 3), pl.Surreal.from_rational(2, 3)))

section("Sprague–Grundy — the impartial center (P-position ⟺ g = 0)")
# Nim heap of size n as the path n → {n−1,…,0}: g(n) = n.
heap = [[j for j in range(h)] for h in range(6)]
print("  Grundy of nim-heap paths 0..5:", pl.grundy_graph(heap), " (g=0 ⟺ Loss/P)")

section("forms now accept ARBITRARY (non-diagonal) metrics — diagonalization")
# A skew hyperbolic plane whose diagonalization is exactly represented.
H = pl.SurrealAlgebra(q=[0, 0], b={(0, 1): 1})
print("  classify skew-H over the surreals:", pl.classify_surreal(H), " (= M₂(ℝ), as ⟨1,−1⟩)")
print("  is_isometric ⟨1,1⟩≅⟨2,3⟩ over F₅ :",
      pl.OddFiniteFieldForm(5, [1, 1]).isometric_to(pl.OddFiniteFieldForm(5, [2, 3])))
print("  isometric ℝ/ℚ/𝔽₂ mirrors      :",
      pl.isometric_real(H, pl.SurrealAlgebra(q=[1, -1])),
      pl.isometric_rational(pl.RationalAlgebra(q=[1, -1]), pl.RationalAlgebra(q=[-1, 1])),
      pl.isometric_surcomplex(pl.SurcomplexAlgebra(q=[1]), pl.SurcomplexAlgebra(q=[-1])),
      pl.isometric_nimber(pl.NimberAlgebra(q=[0, 0], b={(0, 1): 1}), pl.NimberAlgebra(q=[0, 0], b={(0, 1): 1})))

section("Witt decomposition — k·H ⊥ anisotropic kernel")
rw = pl.witt_decompose_real(pl.SurrealAlgebra(q=[1, 1, 1, -1, -1]))
print("  ⟨1,1,1,−1,−1⟩/ℝ:", rw, "fields",
      (rw.witt_index, rw.anisotropic_pos, rw.anisotropic_neg, rw.radical_dim))
ow5 = pl.OddFiniteFieldForm(5, [1, 1]).witt_decompose()
ow3 = pl.OddFiniteFieldForm(3, [1, 1]).witt_decompose()
print("  ⟨1,1⟩/F₅ = H             :", ow5, "index", ow5.witt_index)
print("  ⟨1,1⟩/F₃ anisotropic plane:", ow3, "aniso dim", ow3.anisotropic_dim)
pf = pl.SurrealAlgebra.pfister([-1, -1])
print("  Pfister <<-1,-1>> in I? dim/type:", pf.in_fundamental_ideal(), pf.dim,
      pl.classify_surreal(pf))

section("characteristic polynomial via exterior powers (cₖ = tr Λᵏf)")
R = pl.SurrealAlgebra(q=[1, 1])
# columns f(e_i): M = [[2,1],[3,4]], trace 6, det 5, char poly t²−6t+5.
lin_map = pl.SurrealLinearMap.from_columns([[2, 3], [1, 4]])
print("  char_poly [[2,1],[3,4]] =", R.char_poly(lin_map), " (t²−6t+5)")
print("  trace =", R.trace(lin_map), " det =", R.determinant(lin_map))

section("GA depth — conjugate, scalar/commutator products, meet, blade factoring")
E = pl.SurrealAlgebra(q=[1, 1, 1])
e0, e1, e2 = E.gen(0), E.gen(1), E.gen(2)
print("  Clifford conjugate of e0∧e1   :", (e0 & e1).clifford_conjugate(), " (sign (−1)^{k(k+1)/2})")
print("  scalar product ⟨e0 e0⟩₀       :", e0.scalar_product(e0))
print("  commutator [e0,e1] = 2 e0e1   :", e0.commutator(e1))
blade = (e0 + e1) & e2
print("  blade subspace dimension       :", len(blade.blade_subspace()), blade.blade_subspace())
print("  factor the blade (e0+e1)∧e2   :", blade.factor_blade())
print("  raw blade term/bits/grade     :",
      blade.terms, pl.bits(blade.terms[0][0]), pl.grade(blade.terms[0][0]))
print("  e0∧e1 + e1∧e2 ... meet(planes):", (e0 & e1).meet(e1 & e2), " (their common line, ±e1)")

section("nimber field toolkit — degree, minimal polynomial, order, discrete log")
print("  *2 over F₂: degree", pl.nim_degree(2), " min_poly", pl.nim_min_poly(2), " (x²+x+1)")
print("  conjugates of *5         :", pl.nim_conjugates(5), " order", pl.nim_multiplicative_order(5))
print("  raw nim arithmetic        :", pl.nim_add(5, 7), pl.nim_mul(2, 3),
      pl.nim_pow(2, 3), pl.nim_square(5), pl.nim_frobenius_iter(5, 2), pl.nim_inv(3))
print("  relative trace F₁₆→F₄(*7):", pl.nim_relative_trace(7, 4, 2),
      " norm:", pl.nim_relative_norm(7, 4, 2))
print("  log_*2(*3) in F₄         :", pl.nim_discrete_log(2, 3), " (2² = 3)")
print("  F₂¹²⁸ metadata          :", pl.Nimber.ext_degree(), pl.Nimber.group_order_factors(),
      pl.Nimber(5).frobenius_iter(3))

section("thermography — temperature, mean value, stops")
for g, name in [(pl.Game.switch(1, -1), "{1|-1}"),
                (pl.Game.switch(3, -1), "{3|-1}"),
                (pl.Game.star(), "*"), (pl.Game.up(), "↑")]:
    print(f"  {name:8} temp={g.temperature()!r:>4}  mean={g.mean_value()!r:>4}"
          f"  stops=({g.left_stop()!r},{g.right_stop()!r})")
hot = pl.Game.switch(3, -1)
def same_thermograph(a, b):
    return (a.mean == b.mean and a.temperature == b.temperature and
            a.left_wall.points() == b.left_wall.points() and
            a.right_wall.points() == b.right_wall.points())

th = hot.thermograph()
th_trop = hot.thermograph_via_tropical()
print("  cooled stops at t=1        :", hot.cooled_stops(1),
      " tropical mirror equal:", same_thermograph(th, th_trop))
print("  exact Pl wall at t=1       :", th.left_wall.value_at(pl.Rational(1)),
      " mean/temp:", th.mean, th.temperature,
      " object mirror equal:", same_thermograph(th, th_trop))
print("  module temp/mean/stops      :", pl.temperature(hot), pl.mean_value(hot),
      pl.left_stop(hot), pl.right_stop(hot))
print("  Pl tropical ⊗ sample       :", th.left_wall.otimes(pl.Pl.constant(1)).value_at(1))

section("surreal sign-expansion & floor (the omnific bridge)")
print("  sign expansion of 3/4    :", pl.Surreal.from_rational(3, 4).sign_expansion(), " (+ − +)")
print("  from [+,−,+]             :", pl.Surreal.from_sign_expansion([True, False, True]))
se = pl.SignExpansion.from_finite([True, True, False])
ser = pl.SignExpansion.from_runs([(True, pl.Ordinal(1)), (True, pl.Ordinal(2)), (False, pl.Ordinal(1))])
print("  SignExpansion object      :", se.as_finite(), ser.runs(), ser.length(),
      pl.Surreal.from_sign_expansion_record(se))
print("  surreal order ω > 1      :", pl.omega() > 1)
print("  floor(ω + ½) = ω, frac = :", (pl.omega() + pl.Surreal.from_rational(1, 2)).floor(),
      "/", (pl.omega() + pl.Surreal.from_rational(1, 2)).frac())

section("HACKENBUSH — one structure, three value worlds")
B, Rcol, G = pl.Color.blue(), pl.Color.red(), pl.Color.green()
print("  blue–blue–blue stalk → surreal number :", pl.Hackenbush.string([B, B, B]).value())
print("  green–green stalk     → nimber (Nim)   : *%d" % pl.Hackenbush.string([G, G]).grundy())
print("  blue–red string  → sign-exp surreal    :", pl.Hackenbush.string([B, Rcol]).value(),
      "(= sign expansion + −)")
tri = pl.Hackenbush([(0, 1, G), (1, 2, G), (2, 0, G)])
print("  triangle edge payload                  :", tri.edges())
print("  green triangle (fusion principle)      : *%d" % tri.grundy())

section("ANY NUMBER, BIG OR SMALL — truncated surreal field arithmetic")
w = pl.omega()
print("  √ω = ω^{1/2}                 :", w.sqrt_to_terms(4),
      "  (squares back to", w.sqrt_to_terms(4) ** 2, ")")
print("  √(ω²+2ω+1) = ω+1             :", (w * w + 2 * w + 1).sqrt_to_terms(2))
print("  √2 is not finite-CNF/ℚ       :", pl.surreal(2).sqrt_to_terms(4))   # None — honest scope
print("  exact √ω (no precision arg)  :", w.exact_sqrt(), " is_square(ω):", w.is_square())
print("  1/(ω+1) Neumann to 3 terms   :", (w + 1).inv_to_terms(3))
print("  surreal root helpers         :",
      w.sqrt_to_terms(4),
      (w ** 3).nth_root_to_terms(3, 2))

section("transfinite birthdays & Gonshor sign expansions")
for s, name in [(w, "ω"), (pl.epsilon(), "ε"), (w + 1, "ω+1"),
                (pl.omega_pow(pl.omega()), "ω^ω"), (w.sqrt_to_terms(4), "√ω")]:
    print(f"  birthday({name:3}) = {s.birthday_ordinal()!r:>5}   sign-exp = {s.transfinite_sign_expansion()!r}")
sur_back = pl.Surreal.from_transfinite_sign_expansion((w + 1).transfinite_sign_expansion())
print("  surreal sign-exp inverse ω+1:", sur_back, " ordinal mirror:", pl.Surreal.from_ordinal(pl.Ordinal.omega()).as_ordinal())

section("ordinary (Cantor) ordinal arithmetic — NOT nim")
o = pl.Ordinal.omega()
print("  1 + ω = ω (absorption)       :", pl.Ordinal(1).ord_add(o))
print("  ω + 1 ≠ ω                    :", o.ord_add(pl.Ordinal(1)))
print("  ω + ω = ω·2 (nim would be 0) :", o.ord_add(o))
print("  CNF terms of ω·2             :", o.ord_add(o).terms)

section("NumberGame — games of transfinite birthday (numbers only)")
ng = pl.NumberGame.from_surreal(w)
print("  value =", ng.value(), " birthday =", ng.birthday(), "/", ng.birthday_repr(),
      " short game? =", ng.to_finite_game())
print("  sign-exp(ω) as runs         :", ng.sign_expansion())
ng_back = pl.NumberGame.from_sign_expansion(pl.NumberGame.from_surreal(w + 1).sign_expansion())
print("  sign-exp round trip ω+1     :", ng_back.value())
print("  ω + 1 (delegated to surreal) :", (ng + pl.NumberGame.from_surreal(pl.surreal(1))).value())
one_game = pl.NumberGame.from_surreal(pl.surreal(1))
print("  operator add/neg/order       :",
      (ng + one_game).value(), (-ng).value(), ng > one_game)

section("NimberGame — the char-2 mirror: transfinite Nim heaps ⋆α (No ↔ On₂)")
wg = pl.NimberGame.from_ordinal(pl.Ordinal.omega())  # the heap ⋆ω
print("  grundy(⋆ω) =", wg.grundy(), " short game? =", wg.to_finite_game())
print("  ⋆ω + ⋆ω = 0 (XOR, a P-position) :", (wg + wg).grundy())
print("  ⋆ω + ⋆1 (disjunctive sum)       :", (wg + pl.NimberGame.nim_heap(1)).grundy())
one_heap = pl.NimberGame.nim_heap(1)
print("  operator add/neg/order          :",
      (wg + one_heap).grundy(), (-wg).grundy(), wg > one_heap)
print("  ⋆ω ⊗ ⋆ω ⊗ ⋆ω = ⋆2 (Conway ω³=2):",
      wg.turning_corners(wg).turning_corners(wg).grundy())

section("Cayley transform — bivector (Lie algebra) ↔ rotor (Spin group)")
G = pl.SurrealAlgebra(q=[1, 1, 1])
B = G.gen(0) & G.gen(1)
R = B.cayley()
print("  cayley(e0∧e1) = rotor        :", R, "  norm² =", R.norm2())
print("  cayley_inverse(rotor) = B    :", R.cayley_inverse())
x = G.scalar(1) + G.gen(0) + G.gen(1)            # NOT a simple versor
print("  general inverse of 1+e0+e1   : x·x⁻¹ =", x * x.multivector_inverse())

section("atomic weight — finishing thermography (all-small games)")
for g, name in [(pl.Game.zero(), "0"), (pl.Game.star(), "⋆"), (pl.Game.nim_heap(2), "⋆2"),
                (pl.Game.up(), "↑"), (pl.Game.up().times_int(2), "⇑"),
                (pl.Game.up().times_int(-1), "↓"), (pl.Game.up() + pl.Game.star(), "↑*")]:
    print(f"  aw({name:3}) = {pl.atomic_weight_int(g)!r:>3}   all-small={g.is_all_small()}  tree={g.display()}")

section("game outcomes — Win/Loss/Draw of a finite move graph (kernel)")
# cycle-with-exit: 0→1, 1→{2,0}, 2 terminal. Retrograde analysis resolves it.
succ = [[1], [2, 0], []]
print("  outcomes (Loss=P)            :", pl.outcomes(succ))
print("  P-positions                  :", pl.p_positions(succ))
# Milnor scoring: 0→{1,2} with terminal scores 3,7 ⇒ first-move advantage (7,3).
score = pl.scoring_values([[1, 2], [], []], [0, 3, 7])[0]
print("  scoring (left,right)         :", score, score.left, score.right)

section("loopy games — the Draw escape a cyclic rule opens up")
g = pl.LoopyGraph([[1], [0], [3], []])           # 0↔1 a drawn 2-cycle; 2→3 terminal
print("  adjacency lists             :", g.succ())
print("  outcomes                     :", g.outcomes())
print("  draw-set (the loopy d.o.f.)  :", g.draw_set())
print("  loopy nim-values             :", pl.loopy_nim_values([[1], [0], [3], []]))
print("  loopy values on+off, over+under:",
      pl.LoopyValue.on() + pl.LoopyValue.off(), pl.LoopyValue.over() + pl.LoopyValue.under())
print("  loopy value dud is stopper?  :", pl.LoopyValue.dud().is_stopper(),
      " outcome:", pl.LoopyValue.dud().outcome())
tis_graph = pl.LoopyPartizanGraph(
    [[2], [0], []],        # tis -> 0, tisn -> tis, 0 terminal
    [[1], [2], []],        # tis -> tisn, tisn -> 0, 0 terminal
)
print("  partizan tis/tisn outcomes   :", tis_graph.outcomes(), " classical:", tis_graph.partizan_outcomes())
print("  tis sides / class            :", pl.LoopyValue.tis().sides(), pl.LoopyValue.tis().partizan_outcome())
loopy_rule = lambda v: [[1], [0], [3], []][v]
print("  LoopyGraph.from_rule         :", pl.LoopyGraph.from_rule(4, loopy_rule).succ())
print("  callback loss/draw sets      :", pl.loopy_decision_sets(4, loopy_rule))
print("  callback quadric loss/draw   :", pl.loopy_quadric_probe(2, loopy_rule))

section("misère play — Nim witness + the octal indistinguishability quotient")
print("  misère-Nim P([1,1,1])        :", pl.misere_nim_p_predicted([1, 1, 1]))
print("  callback misère path         :",
      pl.try_misere_is_n(1, lambda h: list(range(h))),
      pl.misere_is_n(1, lambda h: list(range(h))),
      pl.misere_is_p(1, lambda h: list(range(h))))
q = pl.octal_misere_quotient([3, 3, 3], 3, 3, 3)  # 0.333 = Nim, bounded quotient
ag = pl.AbstractGame([[], [0], [0, 1]])
print("  abstract quotient consistency:",
      pl.misere_quotient(ag, [1, 2], 2, 2).num_classes == ag.misere_quotient([1, 2], 2, 2).num_classes)
print("  Nim misère quotient          :", q, " classes:", q.num_classes)
print("  quotient product/signature   :", q.class_product(0, 0), q.signature_of_element(0))

section("char-0 Clifford from a bare signature (no metric needed)")
print("  Cl(0,3) over ℝ               :", pl.classify_real(0, 3))
print("  Cl(3) over ℂ                 :", pl.classify_complex(3))
print("  raw char-0 helpers           :",
      pl.surreal_signature(pl.SurrealAlgebra(q=[1, -1, 0])),
      pl.surcomplex_rank(pl.SurcomplexAlgebra(q=[pl.Surcomplex(1, 0), pl.Surcomplex(0, 0)])))

section("local–global — Hasse–Minkowski + Hilbert reciprocity over ℚ")
ai = pl.isotropy_over_adeles([1, 1, 1])           # ⟨1,1,1⟩ anisotropic over ℚ
print("  ⟨1,1,1⟩ isotropy by place    :", ai, " global:", ai.is_global())
print("  ∏_v (−1,−1)_v = +1 (recip.)  :", pl.hilbert_product((-1, 1), (-1, 1)))
print("  inv_v(−1,−1), sum            :", pl.brauer_local_invariants((-1, 1), (-1, 1)),
      pl.brauer_invariant_sum((-1, 1), (-1, 1)))

section("runtime p-adic cells + adeles — the scalar side of local–global")
k3 = pl.adele_prec(3)
third_3 = pl.LocalQp.from_rational(3, k3, 1, 3)
two_3 = pl.LocalQp.from_int(3, k3, 2)
print("  1/3 in Q₃ has valuation     :", third_3.valuation(), " unit:", third_3.unit)
print("  Local/fixed Qp from_int     :", pl.LocalQp.from_int(3, k3, 9).valuation(),
      pl.Qp3_4.from_int(9).valuation())
print("  2·(1/3) in Q₃               :", two_3 * third_3)
adelic = pl.Adele.from_rational(2, 3)
print("  diagonal 2/3 local at 3      :", adelic.local_at(3),
      " norm:", adelic.idele_norm(), " product formula:", adelic.satisfies_product_formula())
print("  adelic precision policy       :", pl.adele_prec(3))
print("  adding a 3-adic correction   :", adelic.with_correction(3, pl.LocalQp.from_int(3, k3, 1)).local_at(3))
adele_alg = pl.AdeleAlgebra([adelic])
adele_cga = pl.AdeleCga(1)
adele_pt = adele_cga.up([adelic])
print("  adelic Cl/CGA backends        :",
      (adele_alg.gen(0) * adele_alg.gen(0)).terms[0][1] == adelic,
      adele_cga.inner(adele_pt, adele_pt).is_zero())
z8 = pl.Zp2_4(2)
q5 = pl.Qp5_4.from_rational(1, 5)
Q5 = pl.Qp5_4Algebra(q=[1, pl.Qp5_4.from_p_power(1)])
e0q, e1q = Q5.gen(0), Q5.gen(1)
dp_q5 = pl.Qp5_4DividedPowerAlgebra(1)
def local_springer_rows(sp):
    return [(g.valuation, g.dim, g.disc_is_square) for g in sp.graded]

def same_local_springer(a, b):
    return a.radical_dim == b.radical_dim and local_springer_rows(a) == local_springer_rows(b)

print("  fixed Z/2⁴ and Q₅⁽⁴⁾ cells  :", z8.is_unit(), q5.valuation(), (e0q * e1q) ** 2)
print("  p-adic validators            :", pl.Zp2_4.assert_supported_ring() is None,
      pl.Qp5_4.assert_supported_field() is None)
q5_sp = pl.springer_decompose_qp(5, [(1, 1), (5, 1)])
q5_odd = q5_sp.parity_layer(1)[0]
print("  Q₅ Springer and Γ backend    :",
      q5_sp.graded, (q5_odd.valuation, q5_odd.dim, q5_odd.disc_is_square),
      dp_q5.gen(0) * dp_q5.gen(0))
print("  Q₅ Springer repeat check     :",
      same_local_springer(pl.springer_decompose_qp(5, [(1, 1), (5, 1)]), q5_sp),
      same_local_springer(pl.springer_decompose_local(Q5), q5_sp))
print("  Q₅ residue/integral package  :",
      pl.Qp5_4.uniformizer().valuation(), pl.Qp5_4.teichmuller(pl.Fp5(2)).residue(),
      pl.Qp5_4.from_int(25).to_integer(), q5.residue(), q5.residue_unit(),
      q5.is_integral())
print("  valued polynomial Gauss min  :",
      pl.min_coeff_valuation([pl.Qp5_4.from_p_power(2), pl.Qp5_4.zero(), pl.Qp5_4.from_p_power(1)]))
print("  p-adic checked roots         :",
      pl.Zp2_4(4).is_square(), raises_value_error(lambda: pl.Zp2_4(4).sqrt()),
      pl.Qp5_4.from_int(4).sqrt())
f9_ns = pl.F9.primitive_element()
w9_ns = pl.WittVec3_4_2.teichmuller(f9_ns)
q9_ns = pl.Qq3_4_2.from_witt(w9_ns)
Q9 = pl.Qq3_4_2Algebra(q=[1, q9_ns * pl.Qq3_4_2.from_p_power(1)])
e0u, e1u = Q9.gen(0), Q9.gen(1)
dp_w9 = pl.WittVec3_4_2DividedPowerAlgebra(1)
print("  fixed W₄(F₉) and Q₉ cells    :",
      w9_ns.residue(), q9_ns.valuation(), (e0u * e1u) ** 2)
q9_sp = pl.springer_decompose_qq(3, 2, [([1, 0], 0), (f9_ns.coeffs, 1)])
print("  Q₉ Springer and Γ backend    :",
      q9_sp.graded, q9_sp.radical_dim,
      dp_w9.gen(0) * dp_w9.gen(0))
print("  Q₉ Springer repeat check     :",
      same_local_springer(pl.springer_decompose_qq(3, 2, [([1, 0], 0), (f9_ns.coeffs, 1)]), q9_sp))
print("  Q₉ trace/norm/sigma          :",
      q9_ns.trace(), q9_ns.norm(), q9_ns.sigma(),
      len(pl.Qq3_4_2.basis()), pl.Qq3_4_2.embed(pl.Qp3_4.from_p_power(1)).valuation())
q9_scaled = q9_ns * pl.Qq3_4_2.uniformizer()
print("  Q₉ residue/integral package  :",
      q9_scaled.residue(), q9_scaled.residue_unit(), q9_scaled.to_integer(),
      pl.Qq3_4_2.teichmuller(f9_ns).residue())
lr = pl.LaurentRational_6([1, -1])
l5_t = pl.LaurentFp5_6.t()
L5loc = pl.LaurentFp5_6Algebra(q=[1, l5_t])
e0l, e1l = L5loc.gen(0), L5loc.gen(1)
dp_l5 = pl.LaurentFp5_6DividedPowerAlgebra(1)
print("  ℚ((t)) and F₅((t)) cells     :",
      lr.inv(), l5_t.valuation(), (e0l * e1l) ** 2)
print("  F₅((t)) residue package      :",
      pl.LaurentFp5_6.uniformizer().valuation(), l5_t.residue(), l5_t.residue_unit(),
      pl.LaurentFp5_6.teichmuller(pl.Fp5(3)).residue())
l5_sp = pl.springer_decompose_laurent(5, 1, [([1], 0), ([2], 1)])
print("  F₅((t)) Springer and Γ       :",
      l5_sp.graded, l5_sp.graded[0].disc_is_square,
      dp_l5.gen(0) * dp_l5.gen(0))
print("  F₅((t)) Springer repeat check:",
      same_local_springer(pl.springer_decompose_laurent(5, 1, [([1], 0), ([2], 1)]), l5_sp))
r5_pi = pl.RamifiedQp5_4_E2.pi()
R5ram = pl.RamifiedQp5_4_E2Algebra(q=[1, r5_pi])
e0r, e1r = R5ram.gen(0), R5ram.gen(1)
dp_r5 = pl.RamifiedQp5_4_E2DividedPowerAlgebra(1)
print("  ramified Q₅(π), π²=5 cell    :",
      (r5_pi * r5_pi).valuation(), r5_pi.valuation(), (e0r * e1r) ** 2)
print("  ramified residue package     :",
      pl.RamifiedQp5_4_E2.uniformizer().valuation(), r5_pi.residue(), r5_pi.residue_unit(),
      pl.RamifiedQp5_4_E2.teichmuller(pl.Fp5(2)).residue())
ram_sp = pl.springer_decompose_ramified_qp4_e2(5, [[1], [0, 1]])
print("  ramified Springer and Γ      :",
      ram_sp.graded, ram_sp.graded[1].valuation,
      dp_r5.gen(0) * dp_r5.gen(0))
r5_pi3 = pl.RamifiedQp5_4_E3.pi()
R5ram3 = pl.RamifiedQp5_4_E3Algebra(q=[1, r5_pi3])
e0r3, e1r3 = R5ram3.gen(0), R5ram3.gen(1)
dp_r5_3 = pl.RamifiedQp5_4_E3DividedPowerAlgebra(1)
ram3_sp = pl.springer_decompose_ramified_qp4_e3(5, [[1], [0, 1, 0]])
print("  cubic ramified Q₅ and Γ      :",
      (r5_pi3 * r5_pi3 * r5_pi3).valuation(), (e0r3 * e1r3) ** 2,
      ram3_sp.graded, dp_r5_3.gen(0) * dp_r5_3.gen(0))
g5_t = pl.GaussQp5_4.t()
g5_p = pl.GaussQp5_4.from_base(5)
G5gauss = pl.GaussQp5_4Algebra(q=[1, g5_t])
e0g, e1g = G5gauss.gen(0), G5gauss.gen(1)
dp_g5 = pl.GaussQp5_4DividedPowerAlgebra(1)
print("  Gauss Q₅(t), t residue unit  :",
      g5_t.valuation(), g5_p.valuation(), (e0g * e1g) ** 2)
g5_res = pl.Fp5RationalFunction.t() + pl.Fp5RationalFunction.from_base(pl.Fp5(2))
print("  Gauss Q₅(t) residue package  :",
      g5_t.residue(), (g5_p + 2 * g5_t).residue_unit(),
      pl.GaussQp5_4.teichmuller(g5_res).residue())
print("  Gauss Q₅(t) Γ backend        :",
      dp_g5.gen(0) * dp_g5.gen(0))

section("tropical semirings — the dual walls behind thermography")
mx3, mx5 = pl.MaxPlusTropical.finite(3, 1), pl.MaxPlusTropical.finite(5, 1)
mn3, mn5 = pl.MinPlusTropical.finite(3, 1), pl.MinPlusTropical.finite(5, 1)
print("  max-plus 3 ⊕ 5, 3 ⊗ 5       :", mx3 + mx5, mx3 * mx5)
print("  min-plus 3 ⊕ 5, 3 ⊗ 5       :", mn3 + mn5, mn3 * mn5)
print("  tropical operator checks      :", mx3 + mx5, mn3 * mn5)

section("quadric fitting + Gold forms — the Python research bench")
fit = pl.fit_f2_quadratic([0, 1, 2], 2)            # Q=x0x1 zero set
print("  fit zero-set {00,01,10}     :", fit, " genuine:", fit.is_genuinely_quadratic())
gold = pl.gold_form_arf(8, 1)
print("  Gold Q₁ over F₂⁸             :", gold, " rank/rad:", (gold.rank, gold.radical_dim))
gold_alg = pl.gold_form(4, 1)
print("  same Gold form as Cl metric  :", pl.arf_nimber(gold_alg))
print("  Gold/trace helpers           :", pl.arf_nimber(pl.gold_form(4, 1)).arf,
      pl.trace_form_arf(3).arf,
      pl.classify_finite_algebra(pl.trace_twisted_form(3, 2)))
print("  typed trace forms F₈/F₉      :",
      pl.trace_form_arf(3),
      pl.classify_finite_algebra(pl.trace_twisted_form(3, 2)))

section("integral lattices — ADE, genus, mass, Leech constants")
A2 = pl.a_n(2)
E8 = pl.e_8()
print("  A₂ det/min/kissing/Coxeter   :",
      A2.determinant(), A2.minimum(), A2.kissing_number(), pl.coxeter_number(A2),
      pl.is_root_lattice(A2))
print("  E₈ even unimodular aut order :", E8.is_even(), E8.is_unimodular(), E8.automorphism_group_order())
gen = A2.genus()
print("  genus(A₂) primes/symbols     :",
      gen.primes(), gen.symbol_at(3), gen.canonical_symbol_at(2))
print("  mass rank 8 even unimodular  :", pl.mass_even_unimodular(8),
      " Leech |Aut|:", pl.leech_aut_order())
print("  pinned automorphism constants:", pl.E8_WEYL_GROUP_ORDER, pl.D16_PLUS_AUT_ORDER,
      pl.AUTO_NODE_BUDGET)
d16p = pl.d16_plus()
print("  Leech/D16+ constructors       :", pl.leech().dim, d16p.dim, d16p.determinant(), d16p.is_even())
rq = pl.Rational(3, 4)
print("  ℚ scalar 3/4 + 1/6, √9/16   :", rq + pl.Rational(1, 6), pl.Rational(9, 16).sqrt(),
      " sign:", (-rq).sign(), " > 1/2:", rq > pl.Rational(1, 2),
      " reduced:", pl.Rational.try_new(6, 8).numer(), pl.Rational.try_new(6, 8).denom())
RQ = pl.RationalAlgebra([pl.Rational(1), pl.Rational(-3)])
rq_type = pl.classify_rational(RQ)
print("  Cl_Q <1,-3> invariants      :", rq_type)
print("  rational local Hasse records :", rq_type.local_hasse, rq_type.real_closure.ground)
print("  A₂ signature/theta head      :", A2.signature(), A2.theta_series(4))
print("  A₂ rational Clifford metric  :", pl.classify_rational(A2.clifford_metric()))
A2_f2 = A2.clifford_metric_f2()
print("  A₂ mod-2 Clifford Arf        :", pl.arf_nimber(A2_f2) if A2_f2 else None)
code = pl.BinaryCode.extended_hamming()
print("  [8,4,4] code weight/theta    :", code.weight_enumerator(), code.theta_series_via_weight_enumerator(3))
print("  Construction A kissing       :", code.construction_a().kissing_number())
type_i = pl.BinaryCode.type_i_z2()
z2 = type_i.construction_a()
odd_disc = pl.OddDiscriminantForm.from_lattice(pl.IntegralForm.diagonal([3]))
print("  Type I Construction A        :",
      type_i.weight_enumerator(), z2.is_even(), z2.theta_series_level4(5))
print("  odd discr/Milgram report     :",
      odd_disc.group, odd_disc.quadratic_value_mod1([1]), pl.odd_milgram_report(pl.IntegralForm.diagonal([3])))
print("  Golay raw generator rows     :", len(pl.extended_golay_generator_rows()),
      len(pl.extended_golay_generator_rows()[0]))
disc = pl.DiscriminantForm.from_lattice(A2)
print("  discr(A₂) group/Milgram/Weil :", disc.group, disc.milgram_signature_mod8(),
      pl.genus_signature_mod8(A2), disc.verify_weil_relations())
print("  E4 head via modular helper   :", pl.eisenstein_e4(3))
DP = pl.RationalDividedPowerAlgebra(1)
gamma = DP.gen(0)
print("  Γ_Q: γ·γ and coproduct       :", gamma * gamma, gamma.coproduct())

section("exact function fields — F_{2^128}(t) and fixed F_p(t)")
one_plus_t = pl.NimberPoly([1, 1])
print("  (1+t)^2 in char 2            :", one_plus_t * one_plus_t)
print("  gcd(t²+1, t+1)               :", pl.NimberPoly([1, 0, 1]).gcd(one_plus_t))
print("  F₂¹²⁸[t] operator checks      :",
      one_plus_t + pl.NimberPoly.one(), one_plus_t * one_plus_t)
func = pl.NimberRationalFunction([1, 1], [0, 1])
func_t = pl.NimberRationalFunction.t()
print("  ((1+t)/t)·t                  :", func * func_t, " den:", (func * func_t).den)
print("  F₂¹²⁸(t) operator checks      :",
      func - func_t, -func)
FF = pl.NimberRationalFunctionAlgebra([func_t, 1])
print("  Cl_F(t) e0²/e1²              :", FF.gen(0) * FF.gen(0), FF.gen(1) * FF.gen(1))
FDP = pl.NimberRationalFunctionDividedPowerAlgebra(1)
print("  Γ_F(t): t·γ₁                 :", FDP.scalar(func_t) * FDP.gen(0))
fp5_poly_t = pl.Fp5Poly.x()
FP5 = pl.Fp5PolyAlgebra([fp5_poly_t, 1])
print("  Cl_F5[t] e0²/e1²             :", FP5.gen(0) * FP5.gen(0), FP5.gen(1) * FP5.gen(1))
fp5_func_t = pl.Fp5RationalFunction.t()
FF5 = pl.Fp5RationalFunctionAlgebra([fp5_func_t, 1])
FF5DP = pl.Fp5RationalFunctionDividedPowerAlgebra(1)
print("  Cl_F5(t) and Γ_F5(t)         :",
      FF5.gen(0) * FF5.gen(0), FF5DP.scalar(fp5_func_t) * FF5DP.gen(0))
print("  F₅[t]/F₅(t) operator checks   :",
      fp5_poly_t + 1, fp5_func_t * fp5_func_t)
int_poly_t = pl.IntegerPoly.x()
print("  Z[t] eval/rem/gcd             :",
      (5 * int_poly_t + 1) @ 7,
      (int_poly_t * int_poly_t - 1) % (int_poly_t - 1),
      pl.IntegerPoly([2, 2]).gcd(pl.IntegerPoly([4, 4])))
print("  F₅[t]/F₅(t) @ checks          :",
      (fp5_poly_t * fp5_poly_t + 1) @ (fp5_poly_t + 1),
      (1 / fp5_func_t) @ 2)
t = ([0, 1], [1])          # t in F_5(t)
two = ([2], [1])           # nonsquare constant 2 in F_5
norm_form = [([1], [1]), ([0, 4], [1]), ([3], [1]), ([0, 2], [1])]
print("  F₅(t) factors t²+2          :", pl.monic_irreducible_factors(5, [2, 0, 1]))
print("  F₅(t) (t,2) ramifies        :", pl.try_ramified_places_ff(5, t, two),
      " reciprocity:", pl.try_hilbert_reciprocity_product_ff(5, t, two))
ff_adeles = pl.try_isotropy_over_ff_adeles(5, norm_form)
ff_local0 = ff_adeles.local[0] if ff_adeles is not None else None
print("  F₅(t) norm form isotropic?  :", pl.try_is_isotropic_ff(5, norm_form),
      ff_adeles.is_global() if ff_adeles is not None else None,
      ff_local0.place if ff_local0 is not None else None,
      ff_local0.is_isotropic if ff_local0 is not None else None,
      (ff_local0.place, ff_local0.is_isotropic) if ff_local0 is not None else None)
print("  F₅(t) helper checks         :",
      pl.try_valuation_at_ff(5, t, [0, 1]),
      pl.is_global_square_ff(5, ([1, 2, 1], [1])),
      pl.try_hilbert_symbol_ff(5, t, two, [0, 1]),
      pl.try_hilbert_reciprocity_product_ff(5, t, two),
      pl.try_is_isotropic_ff(5, norm_form))
f2_one = ([1], [1])
f2_t = ([0, 1], [1])
print("  F₂(t) factors t²+t+1       :", pl.char2_monic_irreducible_factors([1, 1, 1]))
print("  F₂(t) [1,t) ramifies       :", pl.as_symbol_ramified_places(f2_one, f2_t),
      " reciprocity:", pl.as_symbol_reciprocity_sum(f2_one, f2_t))
print("  F₂(t) [1,t) at t,∞         :",
      pl.as_symbol_at(f2_one, f2_t, [0, 1]),
      pl.as_symbol_at(f2_one, f2_t))
f2_pe = ([0, 1, 1], [1])    # t²+t = wp(t)
char2_form = pl.Char2FunctionFieldForm([(f2_one, f2_t)], [f2_one])
char2_blocks = pl.Char2FunctionFieldForm.from_blocks([(f2_one, f2_t)])
print("  F₂(t) ℘ tests t²+t, t      :", pl.global_is_pe(f2_pe),
      pl.global_is_pe(f2_t))
print("  F₂(t) [1,t]⊥<1> local/global:",
      char2_form.local_anisotropic_dim([0, 1]),
      char2_form.is_isotropic_at_place([0, 1]),
      char2_form.is_isotropic())
aj = char2_form.decompose_at()
print("  F₂(t) AJ decomp at ∞        :", aj, [(term.pole_order, term.coefficient) for term in aj.psi],
      (aj.phi0, [(term.pole_order, term.coefficient) for term in aj.psi], aj.phi1))
print("  F₂(t) helper checks         :",
      pl.as_symbol_places(f2_one, f2_t),
      pl.as_symbol_at(f2_one, f2_t, [0, 1]),
      pl.as_symbol_reciprocity_sum(f2_one, f2_t),
      pl.global_is_pe(f2_pe),
      len(pl.relevant_places_char2(char2_form)),
      pl.local_is_isotropic_char2(char2_form, [0, 1]),
      pl.is_isotropic_global_char2(char2_form),
      char2_blocks.rank())

section("nimber Galois — Frobenius x↦x² and its inverse, the nim √")
n = pl.Nimber(5)
print("  *5² (Frobenius)              :", n.frobenius(), " == *5**2:", n ** 2)
print("  √*5  (inverse Frobenius)     :", n.sqrt(), " squares back:", n.sqrt() ** 2)

section("py-waves parity — Bridge J/K/M/N/O from Python")
lex = pl.lexicode(7, 3)
nimlex = pl.nim_lexicode_naive(2, 2, 2)
print("  lexicode L(7,3)              :",
      (lex.len(), lex.dim(), lex.minimum_distance(), lex.weight_enumerator()))
print("  nim lexicode base 4          :",
      nimlex.words(), nimlex.is_closed_under_nim_scalars())
print("  Brown β and doubled Arf      :",
      pl.brown_f2(1, [1], [0]).beta,
      pl.double_f2([True, True], [2, 1]).beta)
a1_disc = pl.DiscriminantForm.from_lattice(pl.IntegralForm.a(1))
print("  discriminant form iso/Brown  :",
      a1_disc.is_isomorphic(a1_disc),
      a1_disc.brown_invariant())
np = pl.newton_polygon([pl.Qp5_4.from_int(-5), pl.Qp5_4.zero(), pl.Qp5_4.one()])
print("  Newton roots of x²-5 over Q₅ :", np.root_valuations(),
      " τ(5²)=", pl.tropicalize(pl.Qp5_4.from_p_power(2)))
transfer = pl.transfer_diagonal(3, 2, [1])
print("  Scharlau transfer F₉/F₃ <1> :", transfer.dim, pl.classify_finite_algebra(transfer))
quat = pl.Brauer2Class.quaternion(-1, -1)
full = pl.BrauerClass.from_two_torsion(quat)
print("  Brauer 2-torsion → Q/Z      :", quat.ramified_places(), full.local(), full.invariant_sum())
print("  unramified cyclic invariant  :",
      pl.cyclic_algebra_invariant(5, 2, pl.Qp5_4.from_p_power(1)))
ff_t = ([0, 1], [1])
print("  Milnor residues over Q/F₅(t):",
      pl.global_residues([3, 5]),
      pl.global_residues_ff(5, [ff_t]))
print("  constant-extension reciprocity:",
      pl.constant_extension_invariants(5, 3, ff_t),
      pl.constant_extension_invariant_sum(5, 3, ff_t))
