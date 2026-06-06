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
