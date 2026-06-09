"""A tour of pleroma from Python. Run inside the project venv:

    VIRTUAL_ENV=.venv maturin develop && .venv/bin/python demo.py
"""

import pleroma as pl


def section(title):
    print(f"\n── {title} ──")


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
print("  g0 ∧ g1      =", g0 ^ g1, "  (^ is the wedge)")
print("  g0∧g1 == g0 g1:", (g0 ^ g1) == (g0 * g1))

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
print("  3ω² - ω + 5  =", pl.surreal(3) * w ** 2 - w + 5)
print("  ω^ω          =", pl.omega_pow(pl.omega()))
print("  ω > 10⁶ ?    ", w > 1_000_000)
print("  0 < ε < 1e-9?", pl.surreal(0) < pl.epsilon() < pl.rational(1, 10**9))

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
print("  norm² preserved       =", x.norm2(), "->", R.sandwich(x).norm2())
print("  ~(e0 e1)  (reversion) =", ~R)
print("  e0 ⌟ (e0∧e1)          =", e0 << (e0 ^ e1))
print("  dual(e0) in 3D        =", e0.dual(), "  (a bivector)")

section("Arf invariant — the char-2 Clifford classifier (see README.md)")
def arf(qs, bs):
    A = pl.NimberAlgebra(q=[pl.Nimber(x) for x in qs], b={k: 1 for k in bs})
    return pl.arf_invariant(A)
print("  Q = x0·x1        (hyperbolic) :", arf([0, 0], [(0, 1)]))
print("  Q = x0²+x0x1+x1² (anisotropic):", arf([1, 1], [(0, 1)]))
print("  A⊕A ≅ H⊕H (Arf additive)      :", arf([1, 1, 1, 1], [(0, 1), (2, 3)]).o_type)

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
print("  Dickson(swap)  =", pl.dickson_matrix([[0, 1], [1, 0]]), " (a reflection)")
print("  Dickson(diag *2,*3 rotation) =", pl.dickson_matrix([[2, 0], [0, 3]]), " (in SO)")

section("exterior algebra of the GAME group — lives where Clifford can't")
# Λ needs only a ℤ-module; the game group is one, even for non-numbers (⋆, ↑).
ext = pl.GameExterior([pl.Game.star(), pl.Game.up()])
g0, g1 = ext.generator(0), ext.generator(1)
print("  generators are non-numbers:", not ext.game(0).is_number(), not ext.game(1).is_number())
g0g1 = ext.wedge(g0, g1)
print("  g0 ∧ g1 = -(g1 ∧ g0):", g0g1, "==", -ext.wedge(g1, g0))
print("  2·(g0 ∧ g1) = 0       :", ext.is_zero(2 * g0g1), " (relation 2⋆=0 propagates)")
print("  value(g0 + g1) = ⋆ + ↑ :", ext.value_of_grade1(g0 + g1))
print("  value(2·g0) = ⋆+⋆ = 0  :", ext.value_of_grade1(2 * g0) == pl.Game.zero())


# ===========================================================================
# The expansion pass: new scalars, GA configurations, deeper invariants
# ===========================================================================

section("Fp — odd characteristic, completing the classification trichotomy")
# char 0: signature → matrix algebra. char 2: Arf. odd char: dim + discriminant.
print("  F₃ <1,1> :", pl.classify_oddchar(3, [1, 1]))   # disc 1 = square
print("  F₃ <1,2> :", pl.classify_oddchar(3, [1, 2]))   # disc 2 = nonsquare
print("  Hasse always +1 over a finite field:", pl.hasse_invariant(5, [1, 2, 3, 4]))
f9 = pl.FiniteFieldForm(3, [1, 1], degree=2)
print("  F₉ <1,1> via FiniteFieldForm:", f9.classify(), "W=", f9.witt_class())
# the odd-char Witt group: ℤ/4 when −1 is a nonsquare (F₃), ℤ/2×ℤ/2 when it is (F₅)
g3 = pl.oddchar_witt(3, [1]); zero3 = pl.oddchar_witt(3, [])
print("  W(F₃) is ℤ/4 :", g3 + g3 != zero3, "and", g3 + g3 + g3 + g3 == zero3)
g5, h5 = pl.oddchar_witt(5, [1]), pl.oddchar_witt(5, [2])
print("  W(F₅) is ℤ/2×ℤ/2 :", g5 + g5 == pl.oddchar_witt(5, []) and h5 + h5 == pl.oddchar_witt(5, []))

section("Oz — omnific integers: an exterior algebra over a transfinite ring")
Oz = pl.OmnificAlgebra(q=[0, 0, 0])  # Grassmann over Oz
e0, e1, e2 = Oz.gen(0), Oz.gen(1), Oz.gen(2)
w = pl.omnific_omega()
print("  e0² = 0 (nilpotent):", (e0 * e0).is_zero())
print("  (ω·e0) ∧ e1 ∧ e2   :", (w * e0) ^ e1 ^ e2, "  (ω-scale coefficient)")
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
print("  2 ⊗ 2 = *3   :", pl.Ordinal(2).nim_mul(pl.Ordinal(2)))
# nim-multiplication: implemented below ω^ω via the current DiMuro/Conway
# degree-3 tower. The old φ_{ω+1} (<ω³) case is the first layer.
print("  ω ⊗ ω        =", omega.nim_mul(omega), "  (just polynomial mult)")
omega_sq = omega.nim_mul(omega)
print("  ω ⊗ ω ⊗ ω    =", omega_sq.nim_mul(omega), "  (the headline: ω³ = 2)")
# (ω + 1)³ in characteristic 2 = ω³ + ω² + ω + 1 = 2 + ω² + ω + 1 = ω² + ω + *3
w1 = omega.nim_add(pl.Ordinal(1))
print("  (ω+1)³       =", w1.nim_mul(w1).nim_mul(w1), "  (= ω² + ω + nim_add(2,1))")
# finite CNF exponents above 3 now stay inside the implemented tower
print("  ω³ ⊗ ω       :", pl.Ordinal.omega_pow(pl.Ordinal(3)).nim_mul(omega))
print("  ω^ω staged   :", pl.Ordinal.omega_pow(omega).nim_mul(omega))

section("coin-turning games — the game-built nim product")
print("  singleton companions at n=4:", pl.coin_companions("singleton", 4))
print("  singleton Grundy g(n):", [pl.coin_turning_grundy("singleton", n) for n in range(6)])
print("  turtles Grundy g(n)  :", [pl.coin_turning_grundy("turtles", n) for n in range(6)])
print("  Tartan singleton×singleton (2,4):",
      pl.tartan_grundy("singleton", "singleton", 2, 4),
      "= nim_mul_mex", pl.nim_mul_mex(2, 4))
print("  Tartan singleton×turtles (3,2):",
      pl.tartan_grundy("singleton", "turtles", 3, 2),
      "=", pl.nim_mul_mex(3, pl.coin_turning_grundy("turtles", 2)))

section("outermorphisms + determinant — Grassmann's def, char-faithful")
R = pl.SurrealAlgebra(q=[1, 1])
print("  det [[2,1],[3,4]] = 2·4−1·3 =", R.determinant([[2, 3], [1, 4]]))  # columns f(e_i)
N = pl.NimberAlgebra(q=[1, 1])
print("  the char-2 determinant (= permanent):", N.determinant([[2, 1], [3, 1]]))

section("exterior Hopf algebra — antipode = grade involution (not reversion-twist)")
H = pl.SurrealAlgebra(q=[0, 0])  # exterior algebra
b = H.gen(0) ^ H.gen(1)
print("  Δ(e0) primitive (lives in Cl⊗̂Cl):", H.gen(0).coproduct())
print("  S(e0) = −e0          :", H.gen(0).antipode() == -H.gen(0))
print("  S(e0∧e1) = +e0∧e1    :", b.antipode() == b, " (grade 2: (−1)²=+1)")

section("concrete spinor modules — the classification as matrices on spinors")
S = pl.SurrealAlgebra(q=[1, 1, 1])  # Cl(3,0) ≅ M₂(ℂ)
idem, basis, M = S.spinor_rep()
print("  Cl(3,0) minimal ideal real-dim:", len(basis), "(= 2 cols × ℂ)")


def _matmul(a, b):
    n = len(a)
    return [[sum((a[i][k] * b[k][j] for k in range(n)), pl.surreal(0))
             for j in range(n)] for i in range(n)]


M0sq = _matmul(M[0], M[0])
holds = all(M0sq[i][j] == pl.surreal(1 if i == j else 0)
            for i in range(len(basis)) for j in range(len(basis)))
print("  M0² = q0·I (the Clifford relation, on the spinor matrices):", holds)
H2 = pl.NimberAlgebra(q=[0, 0], b={(0, 1): 1})
_idem2, basis2, _M2 = H2.spinor_rep()
print("  char-2 hyperbolic plane spinor ideal dim:", len(basis2), "(idempotent e0e1)")

section("conformal GA over the surreals — exact ∞ and infinitesimal radii")
cga = pl.Cga(2)
p_inf = cga.up([pl.omega(), 0])             # a point at ω-scale
print("  up(ω, 0) still null:", cga.inner(p_inf, p_inf) == pl.surreal(0))
eps = pl.epsilon()
sph = cga.sphere([0, 0], eps * eps)         # a sphere of radius ε
on, off = cga.up([eps, 0]), cga.up([2 * eps, 0])
print("  ε-sphere: ε-point on, 2ε-point off:",
      cga.inner(on, sph) == pl.surreal(0), cga.inner(off, sph) != pl.surreal(0))

section("projective GA — exact nilpotent motor (no transcendentals)")
P = pl.SurrealAlgebra(q=[0, 1, 1])          # Cl(2,0,1), e0 the ideal direction
motor = (P.gen(0) ^ P.gen(1)).exp_nilpotent()  # B² = 0 ⇒ exp(B) = 1 + B
print("  exp(e0∧e1) = 1 + B:", motor == P.scalar(1) + (P.gen(0) ^ P.gen(1)))
print("  it translates e1 ↦ e1 + 2e0:", motor.sandwich(P.gen(1)) == P.gen(1) + 2 * P.gen(0))

section("non-Archimedean Springer decomposition (surreal)")
form = pl.SurrealAlgebra(q=[pl.omega(), pl.epsilon(), 1, -1])
print("  valuation filtration of ⟨ω, ε, 1, −1⟩:")
print("   ", pl.springer_decompose(form))
print("  (W(No)=W(ℝ)=ℤ — the value group is 2-divisible; the filtration is the novelty)")

section("eₙ staircase — discriminant & Hasse as one filtration")
# Over a finite field I²=0, so the staircase is (e₀, e₁) and e₂ is trivial.
e0, e1, e2, stab = pl.e_staircase_oddchar(5, [1, 2, 3])
print(f"  ⟨1,2,3⟩/F5: e0={e0} (dim) e1={e1} (disc) e2={e2:+} (Hasse), I^{stab}=0")
# Over ℝ the tower is infinite: eₙ reads the 2-adic expansion of the signature.
print("  ⟨1,1,1,1⟩/ℝ (sig 4): eₙ for n=0..3 =", [pl.e_real(4, n) for n in range(4)])
print("  numeric invariants of F5:",
      "level", pl.finite_field_level(5),
      "pythagoras", pl.finite_field_pythagoras_number(5),
      "u", pl.finite_field_u_invariant(5))
print("  2 is a sum of two squares in F3:", pl.is_sum_of_n_squares(3, 2, 2))

section("symplectic and Hermitian forms — form + involution siblings")
print("  2 hyperbolic symplectic planes:", pl.SymplecticForm.hyperbolic(2).classify())
print("  char-2 alternating [[0,1],[1,0]]:", pl.classify_symplectic_nimber([[0, 1], [1, 0]]))
H = pl.HermitianForm([[pl.Surcomplex(2, 0), pl.Surcomplex(0, 1)],
                      [pl.Surcomplex(0, -1), pl.Surcomplex(2, 0)]])
print("  Hermitian [[2,i],[-i,2]]:", H.signature(), "diagonal", H.diagonalize())
print("  diagonal Hermitian ⟨1,-1,0⟩:", pl.HermitianForm.diagonal([1, -1, 0]).signature())

section("p-adic Hilbert symbol + Hasse–Minkowski over Q")
print("  (−1,−1)_2 =", pl.hilbert_symbol_qp(-1, -1, 2),
      " — Hamilton's quaternions ramify at 2 (finite fields can't show this)")
for f in ([1, 1, 1], [1, 1, -1], [1, 1, -3], [1, 1, 1, 1, -1]):
    print(f"  ⟨{','.join(map(str, f))}⟩ isotropic over Q:", pl.is_isotropic_q(f))

section("loopy impartial games — Side values with a certificate")
values, cert = pl.loopy_nim_values_certified([[1], [0], []])
print("  2-cycle plus terminal:", values, cert, "outcomes", cert.outcomes)

section("Brauer–Wall group — BW(ℝ)=ℤ/8 is the Bott clock")
# walk ⟨−1⟩⊗̂…⊗̂⟨−1⟩: the Bott index cycles mod 8.
g = pl.bw_class_real(pl.SurrealAlgebra(q=[-1]))
walk, cur = [], g
for _ in range(8):
    walk.append(repr(cur)); cur = cur.add(g)
print("  [Cl⟨−1⟩]ⁿ for n=1..8:", " ".join(w.replace("Real(", "").rstrip(")") for w in walk))
print("  BW(F_3) of ⟨1⟩:", pl.bw_class_oddchar(3, [1]), "(order-4 graded part ≅ W(F_3))")
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
# A messy sum reduces to a single simple game (here ↑ + ⋆ + ⋆ = ↑, since ⋆+⋆=0).
messy = pl.Game.up() + pl.Game.star() + pl.Game.star()
print("  canonical(↑+⋆+⋆) == ↑ :", messy.canonical_string() == pl.Game.up().canonical_string())

section("game ↔ surreal bridge — numbers carry a surreal value")
half = pl.Game.from_surreal(pl.Surreal.simplest_between(0, 1))  # ½ = {0|1}
print("  Game.from_surreal(½)  =", half, " value =", half.number_value(), " birthday =", half.birthday())
print("  ⋆, ↑, ±1 are non-numbers (no value):",
      pl.Game.star().number_value(), pl.Game.up().number_value(), pl.Game.switch(1, -1).number_value())
print("  simplest_between(⅓,⅔) = ½:", pl.Surreal.simplest_between(pl.rational(1, 3), pl.rational(2, 3)))

section("Sprague–Grundy — the impartial center (P-position ⟺ g = 0)")
# Nim heap of size n as the path n → {n−1,…,0}: g(n) = n.
heap = [[j for j in range(h)] for h in range(6)]
print("  Grundy of nim-heap paths 0..5:", pl.grundy_graph(heap), " (g=0 ⟺ Loss/P)")

section("forms now accept ARBITRARY (non-diagonal) metrics — diagonalization")
# A skew hyperbolic plane whose diagonalization is exactly represented.
H = pl.SurrealAlgebra(q=[0, 0], b={(0, 1): 1})
print("  classify skew-H over the surreals:", pl.classify_surreal(H), " (= M₂(ℝ), as ⟨1,−1⟩)")
print("  is_isometric ⟨1,1⟩≅⟨2,3⟩ over F₅ :", pl.is_isometric_oddchar(5, [1, 1], [2, 3]))

section("Witt decomposition — k·H ⊥ anisotropic kernel")
print("  ⟨1,1,1,−1,−1⟩/ℝ (idx,+,−,rad):", pl.witt_decompose_real(pl.SurrealAlgebra(q=[1, 1, 1, -1, -1])))
print("  ⟨1,1⟩/F₅ = H (idx,aniso,□,rad):", pl.witt_decompose_oddchar(5, [1, 1]))
print("  ⟨1,1⟩/F₃ anisotropic plane    :", pl.witt_decompose_oddchar(3, [1, 1]))

section("characteristic polynomial via exterior powers (cₖ = tr Λᵏf)")
R = pl.SurrealAlgebra(q=[1, 1])
# columns f(e_i): M = [[2,1],[3,4]], trace 6, det 5, char poly t²−6t+5.
print("  char_poly [[2,1],[3,4]] =", R.char_poly([[2, 3], [1, 4]]), " (t²−6t+5)")
print("  trace =", R.trace([[2, 3], [1, 4]]), " det =", R.determinant([[2, 3], [1, 4]]))

section("GA depth — conjugate, scalar/commutator products, meet, blade factoring")
E = pl.SurrealAlgebra(q=[1, 1, 1])
e0, e1, e2 = E.gen(0), E.gen(1), E.gen(2)
print("  Clifford conjugate of e0∧e1   :", (e0 ^ e1).clifford_conjugate(), " (sign (−1)^{k(k+1)/2})")
print("  scalar product ⟨e0 e0⟩₀       :", e0.scalar_product(e0))
print("  commutator [e0,e1] = 2 e0e1   :", e0.commutator(e1))
blade = (e0 + e1) ^ e2
print("  factor the blade (e0+e1)∧e2   :", blade.factor_blade())
print("  e0∧e1 + e1∧e2 ... meet(planes):", (e0 ^ e1).meet(e1 ^ e2), " (their common line, ±e1)")

section("nimber field toolkit — degree, minimal polynomial, order, discrete log")
print("  *2 over F₂: degree", pl.nim_degree(2), " min_poly", pl.nim_min_poly(2), " (x²+x+1)")
print("  conjugates of *5         :", pl.nim_conjugates(5), " order", pl.nim_order(5))
print("  relative trace F₁₆→F₄(*7):", pl.nim_relative_trace(7, 4, 2),
      " norm:", pl.nim_relative_norm(7, 4, 2))
print("  log_*2(*3) in F₄         :", pl.nim_discrete_log(2, 3), " (2² = 3)")

section("thermography — temperature, mean value, stops")
for g, name in [(pl.Game.switch(1, -1), "{1|-1}"),
                (pl.Game.switch(3, -1), "{3|-1}"),
                (pl.Game.star(), "*"), (pl.Game.up(), "↑")]:
    print(f"  {name:8} temp={g.temperature()!r:>4}  mean={g.mean_value()!r:>4}"
          f"  stops=({g.left_stop()!r},{g.right_stop()!r})")
hot = pl.Game.switch(3, -1)
print("  cooled stops at t=1        :", hot.cooled_stops(1),
      " tropical mirror equal:", hot.thermograph() == hot.thermograph_via_tropical())

section("surreal sign-expansion & floor (the omnific bridge)")
print("  sign expansion of 3/4    :", pl.rational(3, 4).sign_expansion(), " (+ − +)")
print("  from [+,−,+]             :", pl.Surreal.from_sign_expansion([True, False, True]))
print("  floor(ω + ½) = ω, frac = :", (pl.omega() + pl.rational(1, 2)).floor(),
      "/", (pl.omega() + pl.rational(1, 2)).frac())

section("HACKENBUSH — one structure, three value worlds")
print("  blue–blue–blue stalk → surreal number :", pl.Hackenbush.string(["b", "b", "b"]).value())
print("  green–green stalk     → nimber (Nim)   : *%d" % pl.Hackenbush.string(["g", "g"]).grundy())
print("  blue–red string  → sign-exp surreal    :", pl.Hackenbush.string(["b", "r"]).value(),
      "(= sign expansion + −)")
tri = pl.Hackenbush([(0, 1, "g"), (1, 2, "g"), (2, 0, "g")])
print("  green triangle (fusion principle)      : *%d" % tri.grundy())

section("ANY NUMBER, BIG OR SMALL — truncated surreal field arithmetic")
w = pl.omega()
print("  √ω = ω^{1/2}                 :", w.sqrt(4), "  (squares back to", w.sqrt(4) ** 2, ")")
print("  √(ω²+2ω+1) = ω+1             :", (w * w + 2 * w + 1).sqrt(2))
print("  √2 is not finite-CNF/ℚ       :", pl.surreal(2).sqrt(4))   # None — honest scope
print("  exact √ω (no precision arg)  :", w.exact_sqrt(), " is_square(ω):", w.is_square())
print("  1/(ω+1) Neumann to 3 terms   :", (w + 1).inv_to_terms(3))

section("transfinite birthdays & Gonshor sign expansions")
for s, name in [(w, "ω"), (pl.epsilon(), "ε"), (w + 1, "ω+1"),
                (pl.omega_pow(pl.omega()), "ω^ω"), (w.sqrt(4), "√ω")]:
    print(f"  birthday({name:3}) = {s.birthday_ordinal()!r:>5}   sign-exp = {s.transfinite_sign_expansion()!r}")

section("ordinary (Cantor) ordinal arithmetic — NOT nim")
o = pl.Ordinal.omega()
print("  1 + ω = ω (absorption)       :", pl.Ordinal(1).ord_add(o))
print("  ω + 1 ≠ ω                    :", o.ord_add(pl.Ordinal(1)))
print("  ω + ω = ω·2 (nim would be 0) :", o.ord_add(o))

section("NumberGame — games of transfinite birthday (numbers only)")
ng = pl.NumberGame.from_surreal(w)
print("  value =", ng.value(), " birthday =", ng.birthday_repr(),
      " short game? =", ng.to_finite_game())
print("  sign-exp(ω) as runs         :", ng.sign_expansion())
print("  ω + 1 (delegated to surreal) :", (ng + pl.NumberGame.from_surreal(pl.surreal(1))).value())

section("NimberGame — the char-2 mirror: transfinite Nim heaps ⋆α (No ↔ On₂)")
wg = pl.NimberGame.from_ordinal(pl.Ordinal.omega())  # the heap ⋆ω
print("  grundy(⋆ω) =", wg.grundy(), " short game? =", wg.to_finite_game())
print("  ⋆ω + ⋆ω = 0 (XOR, a P-position) :", (wg + wg).grundy())
print("  ⋆ω + ⋆1 (disjunctive sum)       :", (wg + pl.NimberGame.nim_heap(1)).grundy())
print("  ⋆ω ⊗ ⋆ω ⊗ ⋆ω = ⋆2 (Conway ω³=2):",
      wg.turning_corners(wg).turning_corners(wg).grundy())

section("Cayley transform — bivector (Lie algebra) ↔ rotor (Spin group)")
G = pl.SurrealAlgebra(q=[1, 1, 1])
B = G.gen(0) ^ G.gen(1)
R = B.cayley()
print("  cayley(e0∧e1) = rotor        :", R, "  norm² =", R.norm2())
print("  cayley_inverse(rotor) = B    :", R.cayley_inverse())
x = G.scalar(1) + G.gen(0) + G.gen(1)            # NOT a simple versor
print("  general inverse of 1+e0+e1   : x·x⁻¹ =", x * x.inverse_general())

section("atomic weight — finishing thermography (all-small games)")
for g, name in [(pl.Game.zero(), "0"), (pl.Game.star(), "⋆"), (pl.Game.star_n(2), "⋆2"),
                (pl.Game.up(), "↑"), (pl.Game.up().times_int(2), "⇑"),
                (pl.Game.up().times_int(-1), "↓"), (pl.Game.up() + pl.Game.star(), "↑*")]:
    print(f"  aw({name:3}) = {g.atomic_weight_int()!r:>3}   all-small={g.is_all_small()}")

section("game outcomes — Win/Loss/Draw of a finite move graph (kernel)")
# cycle-with-exit: 0→1, 1→{2,0}, 2 terminal. Retrograde analysis resolves it.
succ = [[1], [2, 0], []]
print("  outcomes (Loss=P)            :", pl.outcomes(succ))
print("  P-positions                  :", pl.p_positions(succ))
# Milnor scoring: 0→{1,2} with terminal scores 3,7 ⇒ first-move advantage (7,3).
print("  scoring (left,right)         :", pl.scoring_values([[1, 2], [], []], [0, 3, 7]))

section("loopy games — the Draw escape a cyclic rule opens up")
g = pl.LoopyGraph([[1], [0], [3], []])           # 0↔1 a drawn 2-cycle; 2→3 terminal
print("  outcomes                     :", g.outcomes())
print("  draw-set (the loopy d.o.f.)  :", g.draw_set())
print("  loopy nim-values (None=Side) :", pl.loopy_nim_values([[1], [0], [3], []]))

section("misère play — Nim witness + the octal indistinguishability quotient")
print("  misère-Nim P([1,1,1])        :", pl.misere_nim_p_predicted([1, 1, 1]))
q = pl.octal_misere_quotient([3, 3, 3], 3, 3, 3)  # 0.333 = Nim, bounded quotient
print("  Nim misère quotient          :", q, " classes:", q.num_classes)

section("char-0 Clifford from a bare signature (no metric needed)")
print("  Cl(0,3) over ℝ               :", pl.classify_real(0, 3))
print("  Cl(3) over ℂ                 :", pl.classify_complex(3))

section("local–global — Hasse–Minkowski + Hilbert reciprocity over ℚ")
ai = pl.isotropy_over_adeles([1, 1, 1])           # ⟨1,1,1⟩ anisotropic over ℚ
print("  ⟨1,1,1⟩ isotropy by place    :", ai, " global:", ai.is_global())
print("  ∏_v (−1,−1)_v = +1 (recip.)  :", pl.hilbert_product((-1, 1), (-1, 1)))

section("runtime p-adic cells + adeles — the scalar side of local–global")
k3 = pl.Adele.finite_precision(3)
third_3 = pl.LocalQp.from_rational(3, k3, 1, 3)
two_3 = pl.LocalQp(3, k3, 2)
print("  1/3 in Q₃ has valuation     :", third_3.valuation(), " unit:", third_3.unit)
print("  2·(1/3) in Q₃               :", two_3 * third_3)
adelic = pl.Adele.from_rational(2, 3)
print("  diagonal 2/3 local at 3      :", adelic.local_at(3),
      " norm:", adelic.idele_norm(), " product formula:", adelic.satisfies_product_formula())
print("  adding a 3-adic correction   :", adelic.with_correction(3, pl.LocalQp(3, k3, 1)).local_at(3))

section("tropical semirings — the dual walls behind thermography")
mx3, mx5 = pl.MaxPlusTropical.finite(3, 1), pl.MaxPlusTropical.finite(5, 1)
mn3, mn5 = pl.MinPlusTropical.finite(3, 1), pl.MinPlusTropical.finite(5, 1)
print("  max-plus 3 ⊕ 5, 3 ⊗ 5       :", mx3 + mx5, mx3 * mx5)
print("  min-plus 3 ⊕ 5, 3 ⊗ 5       :", mn3 + mn5, mn3 * mn5)

section("quadric fitting + Gold forms — the Python research bench")
fit = pl.fit_f2_quadratic([0, 1, 2], 2)            # Q=x0x1 zero set
print("  fit zero-set {00,01,10}     :", fit, " genuine:", fit.is_genuinely_quadratic())
gold = pl.gold_form_arf(8, 1)
print("  Gold Q₁ over F₂⁸             :", gold, " rank/rad:", (gold.rank, gold.radical_dim))
gold_alg = pl.gold_form_algebra(4, 1)
print("  same Gold form as Cl metric  :", pl.arf_invariant(gold_alg))

section("integral lattices — ADE, genus, mass, Leech constants")
A2 = pl.IntegralForm.a(2)
E8 = pl.root_lattice_e8()
print("  A₂ det/min/kissing/Coxeter   :",
      A2.determinant(), A2.minimum(), A2.kissing_number(), A2.coxeter_number())
print("  E₈ even unimodular aut order :", E8.is_even(), E8.is_unimodular(), E8.automorphism_group_order())
gen = A2.genus()
print("  genus(A₂) primes/symbol@3    :", gen.primes(), gen.symbol_at(3))
print("  mass rank 8 even unimodular  :", pl.mass_even_unimodular(8),
      " Leech |Aut|:", pl.leech_aut_order())

section("nimber Galois — Frobenius x↦x² and its inverse, the nim √")
n = pl.Nimber(5)
print("  *5² (Frobenius)              :", n.frobenius(), " == *5**2:", n ** 2)
print("  √*5  (inverse Frobenius)     :", n.sqrt(), " squares back:", n.sqrt() ** 2)
