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

section("Arf invariant — the char-2 Clifford classifier (see NOTES.md)")
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
# nim-multiplication: full across the field φ_{ω+1} (ordinals < ω³), via the
# DiMuro/Conway construction. ω is the nim cube root of 2.
print("  ω ⊗ ω        =", omega.nim_mul(omega), "  (just polynomial mult)")
omega_sq = omega.nim_mul(omega)
print("  ω ⊗ ω ⊗ ω    =", omega_sq.nim_mul(omega), "  (the headline: ω³ = 2)")
# (ω + 1)³ in characteristic 2 = ω³ + ω² + ω + 1 = 2 + ω² + ω + 1 = ω² + ω + *3
w1 = omega.nim_add(pl.Ordinal(1))
print("  (ω+1)³       =", w1.nim_mul(w1).nim_mul(w1), "  (= ω² + ω + nim_add(2,1))")
# higher fields stay staged: any CNF exponent ≥ 3 returns None
print("  ω³·ω staged  :", pl.Ordinal.omega_pow(pl.Ordinal(3)).nim_mul(omega))

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

section("p-adic Hilbert symbol + Hasse–Minkowski over Q")
print("  (−1,−1)_2 =", pl.hilbert_symbol_qp(-1, -1, 2),
      " — Hamilton's quaternions ramify at 2 (finite fields can't show this)")
for f in ([1, 1, 1], [1, 1, -1], [1, 1, -3], [1, 1, 1, 1, -1]):
    print(f"  ⟨{','.join(map(str, f))}⟩ isotropic over Q:", pl.is_isotropic_q(f))

section("Brauer–Wall group — BW(ℝ)=ℤ/8 is the Bott clock")
# walk ⟨−1⟩⊗̂…⊗̂⟨−1⟩: the Bott index cycles mod 8.
g = pl.bw_class_real(pl.SurrealAlgebra(q=[-1]))
walk, cur = [], g
for _ in range(8):
    walk.append(repr(cur)); cur = cur.add(g)
print("  [Cl⟨−1⟩]ⁿ for n=1..8:", " ".join(w.replace("Real(", "").rstrip(")") for w in walk))
print("  BW(F_3) of ⟨1⟩:", pl.bw_class_oddchar(3, [1]), "(order-4 graded part ≅ W(F_3))")


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
# A skew hyperbolic plane: q=[0,0], off-diagonal {e0,e1}=2 ⇒ B(e0,e1)=1.
H = pl.SurrealAlgebra(q=[0, 0], b={(0, 1): 2})
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
