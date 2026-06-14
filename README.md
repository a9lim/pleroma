# ogdoad

[![CI](https://github.com/a9lim/ogdoad/actions/workflows/ci.yml/badge.svg)](https://github.com/a9lim/ogdoad/actions/workflows/ci.yml)
[![crates.io](https://img.shields.io/crates/v/ogdoad)](https://crates.io/crates/ogdoad)
[![PyPI](https://img.shields.io/pypi/v/ogdoad)](https://pypi.org/project/ogdoad/)
[![docs.rs](https://img.shields.io/docsrs/ogdoad)](https://docs.rs/ogdoad)
[![License: AGPL v3](https://img.shields.io/badge/License-AGPL_v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)

`ogdoad` is a Rust research playground for Clifford algebras, quadratic forms,
and combinatorial-game arithmetic, with optional Python bindings. It is built
around one observation: the exotic number systems it implements — surreals,
nimbers, p-adics, Witt vectors, Laurent series — are not a grab bag. They are
cells of *one table*, and the same structures recur from cell to cell with the
characteristic and the place swapped. The code is organized to make those
symmetries visible.

The central constraint is mathematical, not architectural. Conway games under
disjunctive sum form an abelian group, **not a scalar ring** — Conway
multiplication is defined only on the number/nimber subclasses. A Clifford
algebra needs a commutative scalar ring, so this project does **not** build
Clifford algebras over all games. It builds a generic Clifford engine over the
commutative scalar worlds *adjacent* to game theory, and a forms layer that
classifies the result.

## Two views of one table of numbers

Every backend is a cell in a table with two axes:

- **place** — *where* the number lives (Archimedean, p-adic, finite,
  transfinite), and whether it is a field or its ring of integers. This is how
  `src/scalar/` is organized.
- **characteristic** — *which* classification theory applies (char 0 / odd / 2).
  This is how `src/forms/` is organized.

The axes are independent; the two pillars are complementary readings of the same
objects. The place axis pairs each **field** with its **ring of integers**:

| | field | ring of integers |
| --- | --- | --- |
| Archimedean (char 0) | `Rational` ℚ | `Integer` ℤ |
| transfinite | `Surreal` (No) | `Omnific` (Oz) |
| p-adic (char 0) | `Qp`, `Qq` | `Zp`, `WittVec` |
| function field (char p) | `RationalFunction` F_q(t) | `Poly` F_q[t] |
| finite | `Fp`, `Fpn`, `Nimber` | — |

The pairing is structural, not decorative: the `HasFractionField` /
`HasRingOfIntegers` trait pair makes ℤ⊂ℚ, Oz⊂No, Zp⊂Qp, W_N⊂Qq, and F_q[t]⊂F_q(t)
explicit in the type system (with ℤ[i]⊂ℚ[i] following for free via the surcomplex
transport). The rest of the local-field data is structural too — the valuation
and uniformizer (`Valued`), and the residue field `k = 𝒪/𝔪` with angular
component and Teichmüller section (`ResidueField`) — so the whole package
`(K, 𝒪, 𝔪, k, Γ, ϖ)` lives in the type system rather than the comments.

## The symmetries

**char 0 ↔ char 2.** Classifying a quadratic form is one theory split by
`char F`. Over a real-closed field it is the 8-fold periodic Cl(p,q) table
(`M_n(ℝ/ℂ/ℍ)`); in characteristic 2 the quadratic and polar forms part ways and
the same role is played by the Arf invariant and the Brauer–Wall group. On the
finite char-2 legs (`Nimber`, generated `Fpn<2,N>`, the documented finite ordinal
windows) a nonsingular form carries both the Arf classifier and the
`BW(F_{2^m}) ≅ ℤ/2` class, under the same XOR law. The classifier façade picks
the leg from the scalar type at compile time, so `metric.classify()` /
`.bw_class()` are one call across every implemented leg.

**surreal No ↔ ordinal On₂.** The surreals (a char-0 field) and the ordinal
nimbers (a char-2 non-field) are mirror images: both are Cantor-normal-form towers
over recursive exponents, sharing one canonicalizer. They differ in exactly three
places — the exponent order, the coefficient merge (`+` vs `XOR`), and the zero
test — which is why the shared code is a *function*, not a type. No is where
infinite and infinitesimal Clifford metrics live; On₂ is the proper-class char-2
field. The mirror reads out again at the games layer: `NumberGame` (a transfinite
surreal-valued game) and `NimberGame` (a transfinite Nim heap `⋆α` carried by its
ordinal Grundy value) are the two views, one per characteristic.

**the 2×2 functor table.** Orthogonal to the place table, there are four ways to
grow a field, and all four corners are filled:

| | residue-extending | value-extending |
| --- | --- | --- |
| **algebraic** | `Surcomplex` (adjoin `i`) | `Ramified` (adjoin `π = ϖ^{1/e}`) |
| **transcendental** | `Gauss` (adjoin a unit `t`) | `Laurent` (adjoin a uniformizer `t`) |

`Laurent` over a finite field is the equal-characteristic mirror of `Qp`;
`Ramified` is the ramified twin of the unramified `Qq`. The finite *separable*
extensions among these carry a uniform relative trace/norm (`FieldExtension`):
the algebraic-closure functor `Surcomplex`, the finite tower `Fpn/Fp`, the
unramified `Qq/Qp`, and the nim-field `Nimber/F_2` (= `F_{2^128}`) — one interface
for the norm map that feeds Hilbert symbols, the Brauer–Wall group, and Hermitian
forms. The cyclic-Galois refinement (`CyclicGaloisExtension`, adding a basis and
the generator `σ`) feeds the **twisted trace form** `Tr_{E/F}(x·σ^k(x))`, which
lands back in the classifiers — the binary norm form over `Surcomplex`, trace
forms over `Qq` and `Fpn`, and the **Gold form** `Tr(x^{1+2^a})` over the
nim-fields, Arf-classified. The same Galois data also builds Frobenius linear maps
in `clifford::frobenius`, so the scalar trace maps and the Clifford outermorphism
spectra share one basis-level computation.

**local ↔ global.** The Springer decomposition appears across the complete valued
fields, and the value group controls the answer: over the surreals the value group
is 2-divisible (`W(No) = W(ℝ) = ℤ`), but over `Q_p`, the unramified `Q_q`, and
`F_q((t))` it is `ℤ`, so two residue layers survive (`W(Q_p) = W(F_p)²`). The
discretely-valued legs share **one** generic engine keyed on the `ResidueField`
trait; the surreal leg keeps its own, exactly because its value group is divisible
— that mismatch *is* the symmetry, not a gap. The adelic layer then glues the
local data: Hasse–Minkowski isotropy over ℚ and Hilbert reciprocity
`∏_v (a,b)_v = +1`. Those per-prime residues also assemble into Milnor's exact
sequence `0 → W(ℤ) → W(ℚ) → ⊕_p W(F_p) → 0` — the global Witt group with the
Springer residue as its boundary map and the signature as its kernel. The same
package recurs in **equal characteristic** over the
global function field `F_q(t)`: the tame Hilbert symbol at each monic-irreducible
place plus the degree place `∞`, reciprocity, Hasse–Minkowski, and the split Milnor
map `W(F_q(t)) ≅ W(F_q) ⊕ ⊕_π W(F_q[t]/π)` — and here it is **exact** (no precision
model), the char-`p` mirror of the ℚ stack. Both global
fields answer **one** interface: the `GlobalField` trait states the places, the
local Hilbert symbol, reciprocity, and Hasse–Minkowski once, with `ℚ` and `F_q(t)`
as its two implementors.

The integral leg carries its own local/global echo: even lattices produce
discriminant quadratic modules, p-primary Milgram/Brown phase projections, bounded
exact finite-quadratic-module Witt normal forms, and rational or mod-2 Clifford
metrics, making the lattice signature, the real Brauer–Wall mod-8 cycle, and the
Clifford classifier directly comparable in the core. Conway-Sloane `p`-adic genus
symbols, including the corrected 2-adic train/compartment/oddity calculus, give the
integral genus comparison without discriminant-form search budgets. The same leg crosses
the code/theta boundary — binary codes feed Construction A lattices, exact theta
series are identified inside `ℂ[E4, E6]`, `D16+` and `E8 ⊕ E8` share the `E4²`
theta series, Leech is pinned by rootlessness in weight 12, and the Niemeier catalogue
checks the rank-24 mass and weighted theta average against `E12` with the 691
coefficient. Discriminant forms expose Weil `S`/`T` matrices with the Milgram phase
recovered from the standard conjugate `S` prefactor.

**the games bridge.** Red/blue/green Hackenbush is the one object that reads out
as a surreal (blue − red), a nimber (all-green = Nim), or a general partizan game
— and nim-multiplication itself is realized by Conway's Turning-Corners coin game.
This is the seam where the game pillar meets the scalar pillar. The game pillar even
reaches the lattice world: a greedy binary **lexicode** is built by the **mex** rule,
so the Conway–Sloane codes (the `[7,4,3]` Hamming, the `[24,12,8]` Golay) are
Sprague–Grundy P-sets that feed straight into the Construction A lattices of the
integral leg — `mex → lexicode → Golay → Construction A → theta`, one chain crossing
three pillars. The same file also ships base-`2^k` nim-alphabet lexicodes, verifying
nim-additive closure and witnessing scalar-linearity at Fermat bases (4/16) but not
base 8. And thermography itself **is** tropical arithmetic: the option folds are the
tropical `⊕` and cooling is the tropical `⊗`, with the two scaffold walls living in
the dual `(max,+)`/`(min,+)` semirings — named in `scalar/tropical.rs` (a
`Semiring`, not a `Scalar`: an idempotent `⊕` has no inverse) and machine-checked
equal to the golden thermograph.

## The char-2 point

In characteristic 2 the quadratic form and its polar form carry different data.
The engine stores them separately:

```text
e_i^2             = q_i      # the quadratic form
e_i e_j + e_j e_i = b_ij     # the polar / anticommutator (alternating: b_ii = 0)
```

For nimbers `-1 = 1`, so an orthogonal basis with `b = 0` gives a *commutative*
Clifford product; a nonzero off-diagonal `b[(i,j)]` is what makes a
characteristic-2 example noncommutative. Collapsing `q` and `b` into one symmetric
form would silently throw away the entire point of the nimber backend. (An optional
third field `a` lifts the engine to a general, non-symmetric bilinear form.)

On nonsingular metrics over the finite char-2 legs, the form layer also exposes the
Brauer–Wall class as the same Arf/Witt `ℤ/2` datum: hyperbolic planes are zero, the
anisotropic plane has class one, and orthogonal sum / graded tensor adds by XOR.
The spinor module has a separate characteristic-2 representation path: it never uses
the char-0 `½(1+w)` idempotent, accepts nonsingular polar forms such as the
hyperbolic plane with null-square generators, takes blade idempotents like `e_i e_j`
when they shrink a left ideal, and otherwise falls back honestly to the complete
left-regular action.

## Quickstart

Requires Rust and Python ≥ 3.9.

```sh
python3 -m venv .venv
.venv/bin/pip install maturin
VIRTUAL_ENV=.venv .venv/bin/maturin develop
.venv/bin/python demo.py
```

```python
import ogdoad as pl

# characteristic-2 nimber Clifford: non-orthogonal => noncommutative
A = pl.NimberAlgebra(q=[pl.Nimber(2), pl.Nimber(3)], b={(0, 1): 1})
e0, e1 = A.gen(0), A.gen(1)
e0 * e1 + e1 * e0                   # *1  (the anticommutator b[(0,1)])

# surreal metric: infinite and infinitesimal squares are exact
S = pl.SurrealAlgebra(q=[pl.omega(), pl.epsilon()])
(S.gen(0) * S.gen(1)) ** 2         # -1

# the games bridge: Hackenbush reads out as a surreal OR a nimber
B, G = pl.Color.blue(), pl.Color.green()
pl.Hackenbush.string([B, B]).value()      # a surreal number
pl.Hackenbush.string([G, G]).grundy()     # a nimber (all-green = Nim)

# char 0 <-> char 2: a classification on each leg
pl.classify_real(1, 3)             # Cl(1,3) over R, the 8-fold table
pl.arf_nimber(A)                   # the char-2 mirror invariant
pl.bw_class_nimber(A)              # the char-2 Brauer-Wall class, if nonsingular

# local <-> global: Hasse-Minkowski + Hilbert reciprocity over Q
pl.is_isotropic_q([1, 1, 1])       # False (anisotropic over Q)
pl.hilbert_product((-1, 1), (-1, 1))  # +1  (reciprocity)
```

The Python surface is **runtime-friendly parity**: every backend that is a plain
runtime type is bound, while open-ended const-generic families (arbitrary
`Qp<P,K>`, `Qq<P,N,F>`, …) stay Rust-only unless they get an explicit fixed
dispatch slice. See [`src/py/AGENTS.md`](src/py/AGENTS.md) for the full bound
surface and the binding-scope policy.

Run the Rust tour without Python:

```sh
cargo run --example tour
```

## Layout

A pure Rust math core, generic over a `Scalar` trait, with PyO3 per-backend
bindings on top. Each `src/` pillar has its own `AGENTS.md` with the file-by-file
breakdown:

- `src/scalar/` — the `Scalar` trait and every coefficient world, grouped by place.
- `src/clifford/` — the multivector engine, geometric product, and the GA layer
  (versors, outermorphisms, Hopf/divided-power structures, conformal/projective GA,
  spinors, Frobenius linear maps, including the characteristic-2 nimber spinors).
- `src/forms/` — the quadratic-form classifiers across the characteristic
  trichotomy, plus Witt/Brauer–Wall, the Springer trio, `local_global/` for
  Hasse–Minkowski and Hilbert symbols, and `integral/` for lattices, genus,
  discriminant forms, Weil matrices, codes/Construction A, theta/modular forms,
  `D16+`, Leech, and the Niemeier catalogue.
- `src/games/` — normal-, misère-, and loopy-play impartial games, finite
  loopy-partizan graphs, short partizan games, thermography/atomic weight,
  Hackenbush, the exterior algebra of the game group, and the checked integer
  Clifford deformation surface on game generators.
- `src/ogham/` — the Ogham expression-language core: lexer/parser/AST/unparser,
  fixed-world evaluator (Clifford worlds plus polynomial/ratfunc function worlds),
  error taxonomy, and conformance runner support.
- `src/py/` — the optional PyO3 bindings behind the `python` feature.
- `src/linalg/` — crate-private shared linear algebra (exact integer HNF/Smith,
  F₂/nim-field rank, generic field solves), consumed by the pillars above.

See `AGENTS.md` for the working-notes summary, `docs/OPEN.md` for the genuine research
problems, `docs/` (COMPLETENESS.md and CONTINUATIONS.md for the game-valued ledgers of
buildable work and the deferred stars, DONE.md for the go-forward ledger) for the
cross-pillar work, `docs/ogham/`
for the Ogham language contract and hand-verified corpus, and `writeups/goldarf.tex`
for the draft note on the Gold/Arf game thread.

## The bridges — a traveller's catalog

The construction era left the pillars joined by named bridges (summarized in the
`AGENTS.md` files; the catalog below walks them). Five islands: **S**calar,
**C**lifford, **F**orms (the classifier core), the **I**ntegral wing, **G**ames.
Seventeen crossings — Bridge N is four footbridges — each listed with its banks. A
bridge with both feet on one island is a loop; crossing it counts like any other.

| bridge | banks | what it carries |
|---|---|---|
| A | I–C | even lattice → Clifford metric; bounded FQM Witt class and Milgram phase = signature mod 8 |
| B | C–F | char-2 Arf/Brauer–Wall classification over the `Fpn<2,N>` coefficient fields |
| C | S–C | Frobenius/Galois maps as outermorphisms, with flat exterior spectrum |
| D | S–C | `Ordinal` as a checked Clifford scalar — genuinely transfinite char-2 squares |
| E | I–I | theta series identified in `ℂ[E₄,E₆]`; the Milnor isospectral pair, executable |
| F | C–F | the rational Clifford invariant `c(q) = s(q) + δ(n mod 8, disc)`, corrected |
| H | I–I | Construction A: codes ↔ lattices; MacWilliams ↔ the theta transformation |
| I | I–F | the Weil representation of the discriminant form; a third route to σ mod 8 |
| J | S–F | the valuation as (lax) tropicalization; Newton slopes **are** Springer layers |
| K | S–F | the full `ℚ/ℤ` cyclic-algebra Brauer invariant; reciprocity over `F_q(t)` |
| M | F–I | the Brown `ℤ/8` invariant — the char-2 cell of the mod-8 spine, float-free |
| N.1 | F–I | Milnor's exact sequence: the Springer residues go global over `ℚ` and `F_q(t)` |
| N.2 | S–F | the Scharlau transfer, named and tested |
| N.3 | I–I | Nikulin: genus and existence via signature + discriminant form |
| N.4 | I–I | one Bernoulli source for the Eisenstein constants and the mass formula |
| O | G–I | lexicodes: greedy = mex; the `[24,12,8]` lexicode is Golay |
| `game-clifford-checked` | C–G | checked integer Clifford data on game generators; quotient-compatible, not game-native |

(G and L were never built under those letters — they became the deferred stars
`*1` (spinor genus, `docs/COMPLETENESS.md`) and `*2` (the char-`p` Drinfeld mirror,
`docs/CONTINUATIONS.md`). The alphabet itself still has two pontoons
missing; `game-clifford-checked` is the later unlettered C–G span.)

**The traveller's question** (Euler, 1736): can you cross every bridge exactly
once and end where you began? Count the bridge-ends per island:

| island | S | C | F | I | G |
|---|---|---|---|---|---|
| degree | 5 | 6 | **8** | 13 | 2 |

An Euler circuit needs every island even. **Forms — the island the mod-8 spine
runs through — remains balanced, with degree exactly 8**, and the new checked
C–G span also balances Clifford and Games. Only Scalar and the Integral wing are
odd now, so an open Euler stroll exists, but the closed grand tour still does not.
The integral wing, with its four loops (E, H, N.3, N.4), remains the one place a
traveller may wander in circles.

The remaining pure-building closure is now singular: **`*2` (S–I)**, the
Drinfeld/Carlitz mirror, would make every island even. The other open/deferred
threads still matter on their own terms — `*1` for the spinor genus, `under` for a
constructive thermography ↔ Newton-polygon bridge, and `*4` for the wild local
symbol — but after `game-clifford-checked`, they are no longer parity-paired
solutions to the current round-trip obstruction.

## Research thread

The narrow mathematical thread in `docs/OPEN.md` and `writeups/goldarf.tex` is *not* a
claim of a new Clifford classification theorem. It is an investigation of
game-built quadratic forms in the nimber backend:

1. Turning-Corners games realize nim multiplication.
2. Frobenius squaring and traces are built from nim multiplication and XOR.
3. Gold-style trace forms `Tr(λ · x^{1+2^a})` are therefore expressible from
   game-value operations.
4. The Arf invariant gives the standard zero-count bias for a quadratic zero set.
5. The open question is whether a natural, non-tautological game rule has such a
   zero set as its P-positions. Current probes span normal play, misère quotient,
   interactive (`kernel`), loopy (Draw-set), and bent-form searches; they narrow
   the target but do not solve it.

## Status and limits

This is active research code with tests, examples, and experiments. Treat green
tests as regression evidence, not as a proof of the mathematical program. CI runs
`cargo fmt --check`, `cargo clippy --all-targets` (warning-clean), `cargo test`,
`cargo check --features python`, `cargo check --examples`, and `cargo doc --no-deps`
(intra-doc links kept warning-clean).

Scope boundaries worth stating plainly:

- `Nimber(u128)` is exactly `F_{2^128}`. It contains the nim subfields of degree
  dividing 128; it is not the proper-class field of all nimbers.
- `Ordinal` nim-addition is general on the represented CNF terms, and it implements
  `Scalar` for Clifford experiments inside the checked Kummer boundary.
  Nim-multiplication is implemented below `ω^(ω^ω)` when every carry uses a
  verified finite Lenstra excess row: the OEIS A380496 b-file (126 rows, odd primes
  `3..=709`); a carry needing a prime past that table (the first unknown is `719`)
  returns `None`.
  Finite ordinal-nimber metrics expose their minimal detected `F_{2^m}` and use
  that degree for Arf/Witt/Brauer-Wall classification; genuinely transfinite
  metrics remain outside the classifier.
- `Surreal` uses finite support and rational coefficients — the honest truncation
  of true CNF. Non-monomial inverses are infinite Hahn series and are not
  represented.
- `Qp`, `Qq`, `Laurent`, `Ramified`, `Gauss`, and `Adele` are finite-precision
  (capped-relative) models, not exact infinite-memory local fields. They are useful
  for local/global form experiments and excluded from the exact-ring fuzz.
- `ExactScalar` / `ExactFieldScalar` / `PrecisionScalar` name that exact-vs-capped
  boundary explicitly. They are opt-in markers, not `Scalar` supertraits.
- Fixed-width integer payloads are consistently `u128`/`i128` for arithmetic
  carriers, residues, invariants, counts, and budgets. `usize` is used for indices,
  dimensions, and platform ABI hooks.
- The Gold/Arf game thread is conditional: *if* a game has P-set `{Q = 0}`, Arf
  predicts the win-bias. No non-tautological natural game with that P-set has been
  found.

License: AGPL-3.0-or-later (see `LICENSE`).
