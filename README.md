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
`Nimber(n)`, `Surcomplex(re, im)`, `omnific(n)`, `Ordinal(n)`. Module functions:
`arf_invariant(nimber_alg)` (the char-2 Clifford classifier) and `nim_mul_mex(x, y)`
(the Turning-Corners game product).

The expansion pass adds more from Python too: the `OmnificAlgebra` (Oz) and
`Ordinal` (transfinite nimber) backends; `classify_oddchar(p, q)` /
`oddchar_witt(p, q)` (the odd-characteristic third leg of the trichotomy); the GA
methods `.determinant(M)` / `.outermorphism(M, v)` / `.spinor_rep()` /
`.coproduct()` / `.antipode()` / `.exp_nilpotent()`; `Cga(n)` (conformal GA over
the surreals — exact `ω`-scale points and `ε`-radius spheres); and
`springer_decompose(surreal_alg)`. See the second half of `demo.py`.

## architecture

Pure-Rust math core (`cargo test`, no Python in the loop), Python layer on top.
Four pillars under `src/`, each re-exported flat (`scalar::Nimber`,
`clifford::sandwich`, `forms::arf_invariant`, …):

- `scalar/` — the `Scalar` trait + an exact `Rational`/`Integer`, and the
  game-backed coefficient worlds: `nimber` (On₂ in `u64`, nim-add = XOR, nim-mul
  via Fermat-power recursion), `surreal` (Conway normal form with recursive
  exponents), `surcomplex` (adjoin `i`), `omnific` (`Oz`), `onag` (ordinal
  nimbers), `fp` (odd-characteristic prime fields).
- `clifford/` — the multivector engine (`engine`: independent `q`/`b`/`a`,
  characteristic-faithful) with the geometry split out (`versor`), plus
  outermorphisms, the exterior Hopf algebra, conformal/projective GA, and
  spinor modules.
- `forms/` — quadratic-form classifiers across the characteristic trichotomy:
  `char0` (Cl(p,q) → matrix algebra), `oddchar` (discriminant/Hasse), `char2`
  (the Arf invariant), plus the Witt group and the Springer decomposition.
- `games/` — combinatorial game theory: `coin_turning` (nim-mult as Conway's
  Turning-Corners game), normal- and misère-play outcomes, and short partizan
  games with the exterior algebra of the game group.
- `py/` — PyO3 per-backend classes (`python` feature; abi3), split per pillar.

`experiments/` uses the shipped library for the research probe: Arf invariants
of Gold forms over the nim-fields, and the demonstration that those forms are
composites of game operations.

See [`NOTES.md`](NOTES.md) for the math thread: the Arf↔Clifford classification,
the coin-turning↔nim-multiplication bridge to games, and the open question.

Run the Rust tour without Python: `cargo run --example tour`.

## status

**Core + bindings + versor/GA layer + Arf classifier + the expansion pass
complete and verified.** 144 `cargo test` checks green — nim-field axioms and
inverses, Cl(0,1)≅ℂ, Cl(2,0), Grassmann nilpotents, char-2 commutativity *and*
the faithful non-commutative char-2 case, associativity over non-orthogonal
metrics in both characteristics, recursive-exponent surreal arithmetic (`ω^ω`,
`ω·ε=1`, `√ω`), a Clifford metric with `e0²=ω, e1²=ε`, versor inverse /
reflection / rotor / contraction / dual, the Arf invariant (`A⊕A ≅ H⊕H`, and
Gold-function ranks `m−2·gcd(a,m)`), the game definition of nim-multiplication
(Turning Corners), and the surcomplex char-2 degeneracy theorem. The expansion
pass adds: the odd-characteristic trichotomy (`dim+disc` completeness vs a
brute-force congruence oracle; the order-4 Witt group `ℤ/4` vs `ℤ/2×ℤ/2`),
outermorphism determinants (char-faithful), the exterior Hopf axioms (both
characteristics), conformal GA with exact surreal `∞`/`ε`, concrete spinor
modules matching the classifier, omnific-integer exterior algebra, ordinal
nim-arithmetic — including ordinal nim-multiplication across the whole field
`φ_{ω+1}` via the Conway/DiMuro construction, with `ω⊗ω⊗ω = 2` and the
F₄(ω) ≅ F₆₄ field axioms exhaustively verified — and the non-Archimedean
Springer decomposition. Python bindings build as an abi3 wheel and import on
CPython 3.14; `demo.py` tours all of it.

The `experiments/` probes (run on the shipped library) reproduce the
Gold-function ranks, show the Arf-bearing forms are composites of game
operations, and confirm the Arf invariant equals the win-bias (zero-count) of
those forms — see `NOTES.md` for the full thread and the remaining open
question (a *natural* game realizing those forms' P-positions).

## honest limitations / future directions

- **Surreal coefficients are ℚ**, the exact finite stand-in for ℝ (true CNF
  allows any real). Exponents *are* fully recursive surreals.
- **Surreal inverse is exact only for monomials** (`coeff·ω^e`). A genuine sum
  inverts to an infinite Hahn series this finite representation can't hold, so
  `inv()` returns `None`/raises. Versor inverse needs the scalar norm to invert,
  so it works whenever that norm is a monomial (e.g. metrics like `[ω, ε]`).
- **Nimbers cap at `u64`** (the field F_{2^64}); widening to `u128` is mechanical.
- **Ordinal nim-multiplication is complete inside `φ_{ω+1}`** (every ordinal
  `< ω³` Cantor) but **staged above**: the next field would adjoin a degree-5
  root over `φ_{ω+1}`, and the general construction climbs the Lenstra/DiMuro
  tower via `α_p` elements that need real computation in successively larger
  finite fields. Any ordinal whose CNF has an exponent `≥ 3` returns `None`.
- The general char-2 product reduces blade words directly; fine for a playground,
  not optimized for large algebras.
- **The open research question** (explored in `NOTES.md`): the Arf invariant of a
  nim-Clifford form *is* the win-bias of the form in the counting sense, and the
  Arf-bearing forms are built from game operations — but a *natural* game whose
  P-positions are those forms' zero sets is not yet found (normal disjunctive
  play is XOR-linear, so a candidate must be interactive or misère). Also open:
  Gonshor's surreal exponential.
