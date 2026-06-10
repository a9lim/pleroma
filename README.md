# pleroma

`pleroma` is a Rust research playground for Clifford algebras, quadratic forms,
and combinatorial-game arithmetic, with optional Python bindings. It is built
around a single observation: the exotic number systems it implements — surreals,
nimbers, p-adics, Witt vectors, Laurent series — are not a grab bag. They are the
cells of *one table*, and the same structures recur across them with the
characteristic and the place swapped. The code is organized to make those
symmetries visible.

The central constraint is mathematical, not architectural: Conway games under
disjunctive sum form an abelian group, **not a scalar ring**. Conway
multiplication is defined only on the number/nimber sides. So this project does
**not** build Clifford algebras over all games. It builds a generic Clifford
engine over the commutative scalar backends that are *adjacent* to game theory,
and a forms layer that classifies the result.

## Two views of one table of numbers

Every backend is a cell in a table with two axes:

- **place** — *where* the number lives (Archimedean, p-adic, finite, transfinite),
  and whether it is a field or its ring of integers. This is how `src/scalar/` is
  organized.
- **characteristic** — *which* classification theory applies (char 0 / odd / 2).
  This is how `src/forms/` is organized.

The two axes are independent, and the two pillars are complementary readings of
the same objects. The place axis pairs each **field** with its **ring of
integers**:

| | field | ring of integers |
| --- | --- | --- |
| Archimedean (char 0) | `Rational` ℚ | `Integer` ℤ |
| transfinite | `Surreal` (No) | `Omnific` (Oz) |
| p-adic (char 0) | `Qp`, `Qq` | `Zp`, `WittVec` |
| function field (char p) | `RationalFunction` F_q(t) | `Poly` F_q[t] |
| finite | `Fp`, `Fpn`, `Nimber` | — |

The pairing is not decorative — it is made structural by the
`HasFractionField` / `HasRingOfIntegers` trait pair (`ℤ⊂ℚ`, `Oz⊂No`, `Zp⊂Qp`,
`W_N⊂Qq`, `F_q[t]⊂F_q(t)`, and `ℤ[i]⊂ℚ[i]` for free via the surcomplex transport).
The rest of the local-field data is structural too — the valuation and uniformizer
(`Valued`) and the residue field `k = 𝒪/𝔪` with angular component plus
Teichmuller section (`ResidueField`: `Qp→F_p`, `Qq→F_q`, `Laurent→k`,
`Ramified→k`, `Gauss→k(tbar)`) — so the whole package `(K, 𝒪, 𝔪, k, Γ, ϖ)`
lives in the type system rather than the comments.

## The symmetries

**char 0 ↔ char 2.** Classifying a quadratic form is one theory split by `char F`.
Over a real-closed field it is the 8-fold periodic Cl(p,q) table (`M_n(ℝ/ℂ/ℍ)`);
in characteristic 2 the quadratic and polar forms part ways and the same role is
played by the Arf invariant and the Brauer–Wall group. On the finite char-2 legs
(`Nimber`, supported `Fpn<2,N>`, and the documented finite ordinal windows),
nonsingular forms now have both the Arf classifier and the `BW(F_{2^m}) ≅ Z/2`
Brauer-Wall class, with the same XOR law. The classifier façade picks the leg from
the scalar type at compile time, so `metric.classify()` / `.bw_class()` are one
call across the implemented legs.

**surreal No ↔ ordinal On₂.** The surreals (a char-0 field) and the ordinal
nimbers (a char-2 non-field) are mirror images: both are Cantor-normal-form towers
over recursive exponents, sharing one canonicalizer. They differ in exactly three
places — the exponent order, the coefficient merge (`+` vs `XOR`), and the zero
test — which is why the shared code is a *function*, not a type. No is where
infinite/infinitesimal Clifford metrics live; On₂ is the proper-class char-2 field.
The same mirror reads out at the *games* layer: `NumberGame` (a transfinite
surreal-valued game) and `NimberGame` (a transfinite Nim heap `⋆α` carried by its
ordinal Grundy value) are the two views, one per characteristic.

**The 2×2 functor table.** Orthogonal to the place table, there are four ways to
grow a field, and all four corners are filled:

| | residue-extending | value-extending |
| --- | --- | --- |
| **algebraic** | `Surcomplex` (adjoin `i`) | `Ramified` (adjoin `π = ϖ^{1/e}`) |
| **transcendental** | `Gauss` (adjoin a unit `t`) | `Laurent` (adjoin a uniformizer `t`) |

`Laurent` over a finite field is the equal-characteristic mirror of `Qp`;
`Ramified` is the ramified twin of the unramified `Qq`. The functor corners now keep
their residue/value signatures too: `Ramified` is valued by the new uniformizer `π`
with the base residue field, and `Gauss` has unchanged value group with residue
`k(tbar)`. The finite *separable* extensions among these carry a uniform relative
trace/norm (`FieldExtension`): the algebraic-closure functor `Surcomplex` (degree 2),
the finite tower `Fpn/Fp`, the unramified `Qq/Qp`, and the nim-field `Nimber/F_2`
(`F_{2^128}`, the main char-2 field) — one interface for the norm map that feeds
Hilbert symbols, the Brauer–Wall group, and Hermitian forms. The Python fixed
`Qq*_4_*` classes expose the same trace/norm surface into their matching `Qp*_4`
base classes. Python also exposes the fixed local-field residue package
(`uniformizer`, residue/angular-component maps, Teichmüller sections, and
`Qp`/`Qq` ring-of-integers projections) wherever the residue field is already a
bound Python scalar; the fixed Gauss slices return typed `Fp*RationalFunction`
residues for the `F_p(tbar)` field. The cyclic-Galois ones
(`CyclicGaloisExtension`, adding a basis and the generator `σ`) feed the **twisted
trace form** `Tr_{E/F}(x·σ^k(x))` (`forms::trace_form`), which lands back in the form
classifiers: the binary norm form over `Surcomplex`, unramified local trace forms
over `Qq`, ordinary finite-field trace forms over `Fpn`, and the **Gold form**
`Tr(x^{1+2^a})` over the nim-fields, Arf-classified — the typed-core home of the
game-built quadratic-form thread. The same cyclic Galois data now also builds
Frobenius linear maps in `clifford::frobenius`, so the scalar trace maps and
Clifford outermorphism spectra share one basis-level computation.

**local ↔ global.** The Springer decomposition appears across the complete valued
fields, and the value group controls the answer: over the surreals the value group
is 2-divisible (`W(No)=W(ℝ)=ℤ`), but over `Q_p`, the unramified `Q_q`, and
`F_q((t))` it is `ℤ`, so two residue layers survive (`W(Q_p)=W(F_p)²`). The
discretely-valued legs share **one** generic engine keyed on the `ResidueField`
trait; the surreal leg keeps its own, exactly because its value group is divisible —
that mismatch is the symmetry, not a gap. The adelic layer then
glues the local data: Hasse–Minkowski isotropy over ℚ and Hilbert reciprocity
`∏_v (a,b)_v = +1` (the multiplicative product formula). The same local–global
package recurs in **equal characteristic** over the global function field `F_q(t)`
(`forms::function_field`): the tame Hilbert symbol at each monic-irreducible place
plus the degree place `∞`, reciprocity `∏_v (a,b)_v = +1`, and Hasse–Minkowski —
and here it is **exact** (no precision model), the char-`p` mirror of the ℚ stack.
The two global fields answer **one** interface: the `GlobalField` trait states the
places, the local Hilbert symbol, reciprocity, and Hasse–Minkowski once, with `ℚ`
and `F_q(t)` as its two implementors — the local↔global mirror of the
`ResidueField`-keyed Springer engine.

The integral leg has its own local/global echo: even lattices now produce
discriminant quadratic modules, Milgram Gauss-sum phases, and rational or mod-2
Clifford metrics. That makes the lattice signature, the real Brauer-Wall mod-8
cycle, and the Clifford classifier directly comparable in the core. The same leg
now carries the code/theta boundary too: binary codes feed Construction A
lattices, exact theta series are identified inside `C[E4,E6]`, `D16+` and
`E8 + E8` have the same `E4^2` theta series, Leech is pinned by rootlessness in
weight 12, and discriminant forms expose Weil `S`/`T` matrices with the Milgram
phase recovered from the standard conjugate `S` prefactor.

**the games bridge.** Red/blue/green Hackenbush is the one object that reads out as
a surreal (blue−red), a nimber (all-green = Nim), or a general partizan game — and
nim-multiplication itself is realized by Conway's Turning-Corners coin game. This
is the seam where the game pillar meets the scalar pillar. And thermography itself
**is** tropical arithmetic — the option folds are the tropical `⊕` and cooling is
the tropical `⊗`, with the two scaffold walls living in the dual `(max,+)`/`(min,+)`
semirings — named in `scalar/tropical.rs` (a `Semiring`, not a `Scalar`: an
idempotent `⊕` has no inverse) and machine-checked equal to the golden thermograph
by `games::thermograph_via_tropical`.

## The char-2 point

In characteristic 2 the quadratic form and its polar form carry different data.
The engine stores them separately:

```text
e_i^2          = q_i        # the quadratic form
e_i e_j + e_j e_i = b_ij    # the polar / anticommutator (alternating: b_ii = 0)
```

For nimbers `-1 = 1`, so an orthogonal basis with `b = 0` gives a *commutative*
Clifford product; a nonzero off-diagonal `b[(i,j)]` is what makes a
characteristic-2 example noncommutative. Collapsing `q` and `b` into one symmetric
form would silently throw away the entire point of the nimber backend. (An optional
third field `a` lifts the engine to a general, non-symmetric bilinear form.)

On nonsingular metrics over `Nimber`, supported `Fpn<2,N>`, and the documented
finite ordinal windows, the form layer also exposes the Brauer-Wall class as the
same Arf/Witt `Z/2` datum: hyperbolic planes are zero, the anisotropic plane has
class one, and orthogonal sum / graded tensor adds by XOR. The spinor module has a
separate characteristic-2 representation path: it never uses the char-0
`1/2(1+w)` idempotent, accepts nonsingular polar forms such as the hyperbolic plane
with null-square generators, uses blade idempotents like `e_i e_j` when they shrink
a left ideal, and otherwise falls back honestly to the complete left-regular
action.

## Quickstart

Requires Rust and Python ≥ 3.9.

```sh
python3 -m venv .venv
.venv/bin/pip install maturin
VIRTUAL_ENV=.venv .venv/bin/maturin develop
.venv/bin/python demo.py
```

```python
import pleroma as pl

# characteristic-2 nimber Clifford: non-orthogonal => noncommutative
A = pl.NimberAlgebra(q=[pl.Nimber(2), pl.Nimber(3)], b={(0, 1): 1})
e0, e1 = A.gen(0), A.gen(1)
e0 * e1 + e1 * e0          # *1  (the anticommutator b[(0,1)])

# surreal metric: infinite and infinitesimal squares are exact
S = pl.SurrealAlgebra(q=[pl.omega(), pl.epsilon()])
(S.gen(0) * S.gen(1)) ** 2          # -1

# the games bridge: Hackenbush reads out as a surreal OR a nimber
pl.Hackenbush.string(["blue", "blue"]).value()      # a surreal number
pl.Hackenbush.string(["green", "green"]).grundy()   # a nimber (all-green = Nim)

# char 0 <-> char 2: a classification on each leg
pl.classify_real(1, 3)              # Cl(1,3) over R, the 8-fold table
pl.arf_nimber(A)                    # the char-2 mirror invariant
pl.bw_class_nimber(A)               # the char-2 Brauer-Wall class, if nonsingular

# local <-> global: Hasse-Minkowski + Hilbert reciprocity over Q
pl.is_isotropic_q([1, 1, 1])        # False (anisotropic over Q)
pl.hilbert_product((-1, 1), (-1, 1))# +1  (reciprocity)
```

Python exposes the Clifford backends (`NimberAlgebra`, the fixed finite-field
backends `Fp2Algebra`/`Fp3Algebra`/…/`F27Algebra`, the fixed p-adic slices
`Zp2_4Algebra`/…/`Zp13_4Algebra` and `Qp2_4Algebra`/…/`Qp13_4Algebra`,
the fixed unramified slices `WittVec2_4_2Algebra`/…/`WittVec3_4_3Algebra` and
`Qq2_4_2Algebra`/…/`Qq3_4_3Algebra`, the fixed functor slices
`Laurent*_6Algebra`, `RamifiedQp*_4_E2Algebra`/`RamifiedQp*_4_E3Algebra`,
and `GaussQp*_4Algebra`,
`RationalAlgebra`,
`NimberPolyAlgebra`, `NimberRationalFunctionAlgebra`,
`Fp*PolyAlgebra`, `Fp*RationalFunctionAlgebra`, `SurrealAlgebra`,
`SurcomplexAlgebra`, `IntegerAlgebra`, `OmnificAlgebra`, `OrdinalAlgebra`), the
same scalar families' divided-power algebras and Witt-ring representative methods
(`tensor_form`, `pfister1`, `pfister`, `in_fundamental_ideal`), the generic GA
versor/spinor-norm and explicit/lazy spinor-representation methods, raw
blade-term/basis utilities, full outermorphism/inverse-outermorphism utilities,
fixed Frobenius/Galois
matrix helpers, the scalar classes (`Nimber`, `Fp2`,
`Fp3`, `Fp5`, `Fp7`, `Fp11`, `Fp13`, `F4`, `F8`, `F16`, `F9`, `F25`, `F27`,
`Zp2_4`/…/`Zp13_4`, `Qp2_4`/…/`Qp13_4`,
`WittVec2_4_2`/…/`WittVec3_4_3`, `Qq2_4_2`/…/`Qq3_4_3`,
`LaurentRational_6`, `LaurentFp*_6`/`LaurentF*_6`,
`RamifiedQp*_4_E2`/`RamifiedQp*_4_E3`,
`GaussQp*_4`, `NimberPoly`, `NimberRationalFunction`, `Fp*Poly`,
`Fp*RationalFunction`, `Rational`, `Surreal`,
`Surcomplex`, `Integer`, `Omnific`, `Ordinal`), the form helpers (`classify_*`, `classify_rational`,
`arf_nimber`, `arf_ordinal_finite`, `OddFiniteFieldForm`, `Char2FiniteFieldForm`,
finite-algebra classifiers, base-field isometry helpers, trace/Gold-form helpers,
`hilbert_symbol_qp`, `hilbert_symbol_at`, `hasse_at_place`, rational Brauer
local-invariant helpers,
`isotropy_over_adeles`, `bw_class_*`, `BinaryCode`, `DiscriminantForm`,
Milgram/genus-signature helpers, modular q-expansion helpers, …),
and the game helpers (`Game`, `GameExterior`,
`Hackenbush`, `Pl`, `Thermograph`, exact thermograph wall operations,
`LoopyValue`, the kernel outcome/scoring functions, the misère and loopy
machinery, Python-callback Grundy/misère/loopy probes, `nim_mul_mex`,
`grundy_graph`). `springer_decompose_qp` and `springer_decompose_qq` expose the
fixed odd-residue p-adic and unramified Springer cells; `springer_decompose_laurent`
and `springer_decompose_ramified_qp4_e2` / `springer_decompose_ramified_qp4_e3`
do the same for the fixed
equal-characteristic and ramified local cells. Open-ended const-generic families
beyond those fixed dispatch sets (`Qp<P,K>`, `Zp<P,K>`, `Qq<P,N,F>`,
`WittVec<P,N,F>`, arbitrary `Laurent<S,K>`, arbitrary `Ramified<S,E>`,
arbitrary `Gauss<S>`, …) remain Rust-only
unless they get another explicit fixed dispatch set or runtime design.
The exact nimber function-field monomorph `F_{2^128}(t)` is bound, the prime-field
`F_p[t]` / `F_p(t)` rows are bound for the fixed Gauss residue fields with
matching Clifford and divided-power backends, and the
odd-characteristic `F_q(t)` local-global form layer is exposed through fixed
finite-field dispatch; the characteristic-2 Artin-Schreier function-field symbol
and the char-2 Springer local/global form-isotropy layer are exposed the same way
over `F_2(t)`, `F_4(t)`, `F_8(t)`, and `F_16(t)`. See `src/py/AGENTS.md` for the
binding-scope policy.

Run the Rust tour without Python:

```sh
cargo run --example tour
```

## Layout

A pure Rust math core, generic over a `Scalar` trait, with PyO3 per-backend
bindings on top. Each `src/` pillar has its own `AGENTS.md` with the detailed
file-by-file breakdown:

- `src/scalar/` — the `Scalar` trait and every coefficient world, grouped by place.
- `src/clifford/` — the multivector engine, geometric product, and the GA layer
  (versors, outermorphisms, Hopf/divided-power structures, conformal/projective GA,
  spinors, Frobenius linear maps, including characteristic-2 regular/left-ideal
  nimber spinors).
- `src/forms/` — the quadratic-form classifiers and invariants across the
  characteristic trichotomy, plus Witt/Brauer–Wall, the Springer trio,
  `local_global/` for Hasse–Minkowski/Hilbert symbols, and `integral/` for
  lattices, genus, discriminant forms, Weil matrices, Milgram checks, mass, binary
  codes/Construction A, theta/modular forms, `D16+`, and Leech.
- `src/games/` — normal-, misère-, and loopy-play impartial games, short partizan
  games, thermography/atomic weight, Hackenbush, and the exterior algebra of the
  game group.
- `src/py/` — the optional PyO3 bindings behind the `python` feature.

See `AGENTS.md` for the working-notes summary, `OPEN.md` for the genuine research
problems, `ROADMAP.md` for the implemented and proposed cross-pillar bridges and
remaining boundaries, and
`writeups/goldarf.tex` for the narrow draft note on the Gold/Arf game thread.

## Research thread

The narrow mathematical thread in `OPEN.md` and `writeups/goldarf.tex` is not a
claim of a new Clifford classification theorem. It is an investigation of
game-built quadratic forms in the nimber backend:

1. Turning-Corners games realize nim multiplication.
2. Frobenius squaring and traces are built from nim multiplication and XOR.
3. Gold-style trace forms `Tr(λ · x^(1+2^a))` are therefore expressible from
   game-value operations.
4. The Arf invariant gives the standard zero-count bias for a quadratic zero set.
5. The open question is whether a natural, non-tautological game rule has such a
   zero set as its P-positions. Current probes span normal-play, misère quotient,
   interactive (`kernel`), loopy (Draw-set), and bent-form searches; they narrow
   the target but do not solve it. The detailed problem statement and probe map are
   in `OPEN.md`.

## Status and limits

This is active research code with tests, examples, and experiments. Treat green
tests as regression evidence, not as a proof of the mathematical program. The CI
checks are `cargo fmt --check`, `cargo clippy --all-targets` (warning-clean),
`cargo test`, `cargo check --features python`, and `cargo check --examples`.

Scope boundaries worth stating plainly:

- `Nimber(u128)` is exactly `F_{2^128}`. It contains the nim subfields of degree
  dividing 128; it is not the proper-class field of all nimbers.
- `Ordinal` nim-addition is general on the represented CNF terms and it implements
  `Scalar` for Clifford experiments inside the checked Kummer boundary.
  Nim-multiplication is implemented below `ω^(ω^ω)` when every carry uses the
  verified excess table: DiMuro through `alpha_u` for `u <= 43`, plus the local
  fixed-base oracle value `alpha_47`; a carry past that table returns `None`.
- `Surreal` uses finite support and rational coefficients — the honest truncation of
  true CNF. Non-monomial inverses are infinite Hahn series and are not represented.
- `Qp`, `Qq`, `Laurent`, `Ramified`, `Gauss`, and `Adele` are finite-precision
  (capped-relative) models, not exact infinite-memory local fields. They are useful
  for local/global form experiments and excluded from the exact-ring fuzz.
- `ExactScalar` / `ExactFieldScalar` / `PrecisionScalar` name that exact-vs-capped
  precision boundary explicitly. They are opt-in markers, not `Scalar` supertraits.
- Fixed-width integer payloads are consistently `u128`/`i128` for arithmetic
  carriers, residues, invariants, counts, and budgets. `usize` is still used for
  indices, dimensions, and platform ABI hooks.
- The Gold/Arf game thread is conditional: *if* a game has P-set `{Q = 0}`, Arf
  predicts the win-bias. No non-tautological natural game with that P-set has been
  found.
