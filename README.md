# pleroma

A playground for **Clifford algebras (with nilpotents)** over the field-like
subclasses of Conway's combinatorial games — built as a verified Rust core with
Python bindings.

## why this is the interesting shape

Games under disjunctive sum are an abelian group but **not a ring** — Conway's
multiplication is only a congruence on the *numbers*. A Clifford algebra needs a
commutative scalar ring, so "Clifford over games" funnels you into the three
field-like cores of game-world, and the answer forks on which you pick:

| scalar backend | structure | what Clifford does there |
|---|---|---|
| **surreals** `No` | real-closed field, char 0 | exactly the Cl(p,q) classification over ℝ — 8-fold periodicity — but with surreal (infinite/infinitesimal) metric entries |
| **surcomplex** `No[i]` | algebraically closed, char 0 | complex Clifford theory, 2-fold periodicity |
| **nimbers** `On₂` | algebraically closed, char **2** | the genuinely new beast: char-2 Clifford, alternating bilinear form |

"With nilpotents" = the quadratic form may be **degenerate**: a basis vector
with `Q(eᵢ) = 0` squares to zero. `Q ≡ 0` recovers the full Grassmann/exterior
algebra; mixed signatures give `Cl(p,q,r)` with `r` null directions.

## char 2 is not a footnote

Over the nimbers, `eᵢeⱼ + eⱼeᵢ = 2B = 0`, so the bilinear form drops out of the
anticommutator and the algebra over an *orthogonal* basis becomes commutative.
To get the real (non-commutative) char-2 theory the engine carries the quadratic
form `q` (the squares `eᵢ²`) **independently** of the alternating off-diagonal
form `b` (the anticommutators `{eᵢ,eⱼ}`) — they no longer determine each other.
A nonzero off-diagonal `b[(i,j)]` is exactly what makes a nim-Clifford algebra
non-commutative.

## quickstart

Needs Rust (`curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh`)
and a Python ≥ 3.9. From the repo root:

```sh
python3 -m venv .venv
.venv/bin/pip install maturin
VIRTUAL_ENV=.venv .venv/bin/maturin develop   # builds the abi3 extension
.venv/bin/python demo.py
```

```python
import pleroma as pl

# char-2 nimber Clifford, non-orthogonal ⇒ non-commutative
A = pl.NimberAlgebra(q=[pl.Nimber(2), pl.Nimber(3)], b={(0, 1): 1})
e0, e1 = A.gen(0), A.gen(1)
e0 * e1            # e0e1            (* is the geometric product)
e0 * e1 + e1 * e0  # *1              (the anticommutator b[(0,1)])

# surreal metric: infinite & infinitesimal squares, unit bivector anyway
S = pl.SurrealAlgebra(q=[pl.omega(), pl.epsilon()])
(S.gen(0) * S.gen(1)) ** 2          # -1

# surreal scalars: recursive exponents, exact order
w = pl.omega()
(w + 1) * (w - 1)                   # ω^2 - 1
pl.omega_pow(pl.omega())           # ω^(ω)
w > 1_000_000                       # True
```

```python
# versors: reflections and rotations (Cl(3,0))
E = pl.SurrealAlgebra(q=[1, 1, 1])
e0, e1 = E.gen(0), E.gen(1)
e1.reflect(e0)                      # e0   (reflection in the hyperplane ⊥ e1)
(e0 * e1).sandwich(3*e0 + 4*e1)    # rotor action; norm² preserved
e0.inverse()                        # versor inverse
~(e0 * e1)                          # reversion (~ = reverse)
e0 << (e0 ^ e1)                     # left contraction  (e1)
e0.dual()                           # Hodge dual
```

Operators: `*` geometric, `^` wedge, `<<`/`>>` left/right contraction, `**`
power, `/` divide (by scalar or versor), `~` reverse, `+`/`-`, `==`. Multivector
methods: `.inverse()`, `.sandwich(x)`, `.reflect(x)`, `.dual()`, `.norm2()`,
`.grade(k)`, `.grade_involution()`, `.reverse()`, `.scalar_part()`,
`.is_zero()`; algebra has `.pseudoscalar()`. Scalars have `.inv()` and `/`.
Builders: `omega()`, `epsilon()`, `omega_pow(x)`, `rational(p, q)`, `surreal(n)`,
`Nimber(n)`, `Surcomplex(re, im)`.

## architecture

Pure-Rust math core (`cargo test`, no Python in the loop), Python layer on top.

- `scalar.rs` — the `Scalar` trait + an exact `Rational` (engine validation only)
- `nimber.rs` — On₂ in `u64` (= F_{2^64}): nim-add = XOR, nim-mul via Fermat-power recursion
- `clifford.rs` — the multivector engine, generic over `Scalar`, with independent `q`/`b` (characteristic-faithful)
- `surreal.rs` — Conway normal form with recursive surreal exponents (ℚ coefficients)
- `surcomplex.rs` — adjoin `i` over any backend
- `py.rs` — PyO3 per-backend classes (`python` feature; abi3)

Run the Rust tour without Python: `cargo run --example tour`.

## status

**Core + bindings + versor/GA layer complete and verified.** 34 `cargo test`
checks green — nim-field axioms and inverses, Cl(0,1)≅ℂ, Cl(2,0), Grassmann
nilpotents, char-2 commutativity *and* the faithful non-commutative char-2 case,
associativity over non-orthogonal metrics in both characteristics, recursive-
exponent surreal arithmetic (`ω^ω`, `ω·ε=1`, `√ω`), a Clifford metric with
`e0²=ω, e1²=ε`, versor inverse / reflection / rotor / contraction / dual, and
the surcomplex char-2 degeneracy theorem. Python bindings build as an abi3 wheel
and import on CPython 3.14.

## honest limitations / future directions

- **Surreal coefficients are ℚ**, the exact finite stand-in for ℝ (true CNF
  allows any real). Exponents *are* fully recursive surreals.
- **Surreal inverse is exact only for monomials** (`coeff·ω^e`). A genuine sum
  inverts to an infinite Hahn series this finite representation can't hold, so
  `inv()` returns `None`/raises. Versor inverse needs the scalar norm to invert,
  so it works whenever that norm is a monomial (e.g. metrics like `[ω, ε]`).
- **Nimbers cap at `u64`** (the field F_{2^64}); widening to `u128` is mechanical.
- The general char-2 product reduces blade words directly; fine for a playground,
  not optimized for large algebras.
- **Not yet explored:** whether a nim-Clifford form's Arf invariant over the
  finite nim-fields F_{2^{2^n}} carries Sprague–Grundy meaning — the one thread
  that would make this "Clifford that knows it's made of games." Also: Gonshor's
  surreal exponential, and an inverse/division op on the surreal backend.
