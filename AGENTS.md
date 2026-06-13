# AGENTS.md — ogdoad

Working notes for agents editing this repo. Global rules still apply. This file is
the **summary**; each `src/` pillar has its own `AGENTS.md` with the file-by-file
breakdown and the layer-specific "things that look like bugs".

## What this is

Clifford algebras (with nilpotents) over commutative scalar worlds adjacent to
Conway's combinatorial games. Games under disjunctive sum are an abelian group,
**not a ring** — Conway multiplication is only defined on the number/nimber
subclasses. A Clifford algebra needs a commutative scalar ring, so the direct
game-valued Clifford story only lives on the field-like cores (nimbers, surreals,
surcomplex). The repo also carries comparison scalar worlds (Fp/Fpn, Zp/Qp/Qq,
Laurent, ramified/Gauss functors, an adelic precision model, and the exact global
function field `F_q(t)` = `RationalFunction` over `Poly` = `F_q[t]`) for form-theory
experiments. A pure Rust math core, generic over a `Scalar` trait, with PyO3
per-backend bindings on top. "With nilpotents" = the quadratic form may be
degenerate (`q[i]=0` ⇒ `eᵢ²=0`); all-zero `q` is the exterior/Grassmann algebra.

## Layout — the four pillars + bindings

`src/lib.rs` is the crate root: the four pillars + the (feature-gated) `py` module.
Each pillar's `mod.rs` re-exports its children flat, so public paths stay shallow
(`scalar::Nimber`, `clifford::Metric`, …).

| dir | pillar | detail |
|---|---|---|
| `src/scalar/`   | commutative coefficient worlds (the `Scalar` trait + every backend) | [`src/scalar/AGENTS.md`](src/scalar/AGENTS.md) |
| `src/clifford/` | the multivector engine + the GA layer | [`src/clifford/AGENTS.md`](src/clifford/AGENTS.md) |
| `src/forms/`    | quadratic forms & invariants, by the char trichotomy plus local-global and integral layers | [`src/forms/AGENTS.md`](src/forms/AGENTS.md) (+ [`integral/`](src/forms/integral/AGENTS.md)) |
| `src/games/`    | combinatorial game theory | [`src/games/AGENTS.md`](src/games/AGENTS.md) |
| `src/ogham/`    | the Ogham expression-language core (parser, fixed-world evaluator, polynomial/ratfunc function worlds, errors, conformance support) | root rules |
| `src/py/`       | PyO3 bindings (feature = "python") + the binding-scope policy | [`src/py/AGENTS.md`](src/py/AGENTS.md) |
| `src/linalg/`   | crate-private shared linear algebra | [`src/linalg/AGENTS.md`](src/linalg/AGENTS.md) |

Beyond the library: `examples/` (Rust demos `tour`/`tropical` + the open-question
probes `interactive_kernel`, `octal_hunt`, `loopy_quadric`, `misere_quotient`,
`bent_route`), `experiments/` (Python research probes on top of the shipped
lib), `demo.py` (the Python tour),
`docs/` (OPEN.md — the genuine research problems; COMPLETENESS.md — the game-valued
ledger of buildable items completing symmetries/connections already in the code;
CONTINUATIONS.md — the game-valued ledger of buildable items that are genuinely new
features (the ogham language work, the char-`p` Drinfeld mirror); the deferred stars
`*1`/`*2`/`*4`/`*8` split across those two; DONE.md — the go-forward ledger for new
work; CONSISTENCY.md — the aesthetic/structural ledger; CORRECTNESS.md — the
verification-status ledger (machine-verified / source-pinned / asserted); TABLES.md —
the inventory of curated hardcoded tables),
`docs/ogham/` (ogham.md — the expression-language spec, the shipped Display-v2 +
host-operator contract, backend-helper surface, v1 parser/evaluator contract,
v1.1 polynomial/ratfunc function-world contract, the shipped v2.0 abstraction
layer, the shipped v2.1 program layer, and the pre-contract v3.0 stub, §§17–19;
conformance.txt — the hand-verified corpus the language must pass),
and `writeups/`
(`goldarf.tex` — the consolidated draft note on the Gold/Arf game thread,
including the Tier-2 no-go/construction program; `excess.tex` — the
consolidated note on the transfinite nim excess problem).

## Claim levels and non-claims

Use these labels when changing prose, papers, examples, or comments:

- **Standard math**: external facts such as Sprague-Grundy, the Turning-Corners
  product theorem, Arf classification of nonsingular binary quadratic forms, the
  Gold rank formula, zero-count formulas over `F_2`, local Hilbert/Artin-Schreier
  reciprocity, and Conway-Sloane lattice theory.
- **Implemented and tested**: statements backed by the Rust tests, examples, Python
  experiments, or the `demo.py` tour in this checkout.
- **Interpretation**: bridges such as "Arf is a win-bias" are conditional on a
  game whose P-set is the corresponding quadratic zero set.
- **Open**: the natural Gold-quadric game rule, a genuine game-native quadratic
  deformation of `GameExterior`, and transfinite nim multiplication beyond the
  source-verified excess table. These live in `docs/OPEN.md`.

Scope boundaries to preserve:

- Not a Clifford algebra over arbitrary partizan games. A Clifford algebra needs a
  commutative scalar ring; the full game group is only an abelian group.
- Not a new classification theorem for all characteristic-2 Clifford algebras over
  arbitrary fields. The code computes Arf and Brauer-Wall data for `Nimber`,
  supported `Fpn<2,N>` fields, and the documented finite ordinal windows; it rejects
  singular metrics where a nonsingular Witt/BW class is required, and keeps
  rank/radical data explicit.
- Not a solved game-semantics theorem. Gold forms are built from game operations,
  but no non-tautological natural game is known whose P-set is their zero set.
- Not an algebraically closed finite backend. `Nimber(u128)` is `F_{2^128}` and
  contains only finite nimber subfields whose degrees divide 128.

## The symmetries (the intellectual spine — see README)

The two pillars are complementary views of the *same* table of numbers:
`scalar/` groups **by place** (each field beside its ring of integers); `forms/`
cuts ACROSS by **characteristic** (the one classification theory, three ways). The
recurring symmetries the project is built around:

- **char 0 ↔ char 2** classifiers (the real 8-fold table mirrored by Arf/Brauer–Wall).
- **surreal No ↔ ordinal On₂** (char-0 field and char-2 non-field sharing one CNF core).
- **(field, ring of integers)** pairings, made structural in `scalar/integrality.rs`.
- the **2×2 functor table** (algebraic|transcendental × residue|value-extending),
  with cyclic Galois/Frobenius maps also feeding Clifford outermorphisms.
- **local ↔ global** (the Springer trio + the local-global Hasse–Minkowski layer).

## Implemented mathematical state

The scalar landscape is broad, but not all backends have the same exactness claim:

| backend | role |
|---|---|
| `Nimber(u128)` | finite nim-field `F_{2^128}` with nim add/mul; main char-2 backend |
| `Surreal` | finite-support Hahn/CNF char-0 backend; real-table classification only on represented exact square classes |
| `Surcomplex` | `Surreal[i]`; complex-table classification only on represented exact square classes |
| `Integer`, `Omnific` | coefficient rings for exterior/nilpotent structures |
| `Fp`, `Fpn`, `Zp`, `WittVec` | comparison scalar worlds for the characteristic trichotomy |
| `Qp`, `Qq`, `Laurent`, `Ramified`, `Gauss` | local-field-style backends/functors, mostly capped-relative precision models |
| `Adele`, `LocalQp` | runtime-prime adelic precision model over `Q` |
| `RationalFunction` | exact global function field `F_q(t)` over `Poly = F_q[t]` |
| `Ordinal` | staged transfinite nimbers; a checked/panic-on-escape `Scalar` for Clifford metrics; nim-addition on represented CNF terms, nim-multiplication through Kummer carries `α_u` assembled from `ord_u(2)`, `Q(f(u))`, and the source-verified finite `m_u` rows (DiMuro `u ≤ 43`, plus certified `m_47=1`) below `ω^(ω^ω)` |

The char-2 Clifford point is load-bearing. In characteristic 2, `q` and `b` are
independent:

```text
e_i^2             = q_i
e_i e_j + e_j e_i = b_ij
```

The polar form is alternating, so `b_ii = 0`, while `q_i` can be nonzero. If the
engine collapses `q` and `b`, the char-2 Clifford product becomes the wrong
commutative object. The char-2 form layer reports `ArfInvariants { arf: u128, rank,
radical_dim, radical_anisotropic }` plus the derived `o_type()` →
`OrthogonalType`; for degenerate forms, the Arf of the
nonsingular core is not the whole form. On nonsingular metrics over `Nimber`,
supported `Fpn<2,N>`, and the documented finite ordinal windows, the same Arf bit is
the Brauer-Wall class `BW(F_{2^m}) ≅ Z/2`; hyperbolic planes are `0`, the
anisotropic plane is `1`, and direct sum / graded tensor adds by XOR.
`clifford::spinor` has a separate char-2 route: no `1/2(1+w)`, blade idempotents
such as `e_i e_j` when they shrink a left ideal, and otherwise the full
regular/lazy left action. Singular polar forms and general-bilinear `a` metrics are
rejected.

The cross-pillar bridges live in the Rust core. `IntegralForm` exports rational and even-mod-2 Clifford metrics plus
discriminant Gauss-sum/Milgram checks; finite char-2 `Fpn<2,N>` classification runs
through the façade; cyclic Galois/Frobenius maps have Clifford linear-map constructors;
the **rational 2-torsion Brauer class** `Brauer2Class` (`witt/brauer_rational.rs`:
Hasse–Witt `s(q)` vs the Clifford invariant `c(q) = s(q) + δ(n mod 8, disc)`) and its
**full `ℚ/ℤ` lift** `BrauerClass` (`witt/cyclic.rs`: `cyclic_algebra_invariant = v(a)/n`,
with full-strength reciprocity over `F_q(t)`); the **valuation as (lax) tropicalization**
with `NewtonPolygon` over the valued legs (`scalar/newton.rs`, slope = root valuation =
Springer residue layer); `Ordinal` serves as a Clifford scalar inside the verified Kummer
boundary;
`forms/integral/codes.rs` carries binary codes, MacWilliams, and Construction A
(with the `1/sqrt(2)` scaling and an `Option` boundary when the scaled Gram is not
integral), including the Type II length-16 code whose lattice is `D16+`;
`forms/integral/{theta,modular}.rs` give exact theta coefficients and `E4`/`E6`/`E12`
identification (`theta_E8 = E4`, `theta_{E8+E8} = theta_{D16+} = E4^2`, the rootless
Leech `q^1` oracle), while `forms/integral/niemeier.rs` carries the 24-class
Niemeier root/glue/Aut catalogue and verifies the rank-24 mass plus weighted
Siegel-Weil identity against `E12` and the 691 coefficient; and
`DiscriminantForm` exposes dependency-free `Complex64` Weil
`S`/`T` matrices, with the `S` prefactor the conjugate of the positive Milgram phase
and `verify_weil_relations` checking the honest metaplectic relations (not the
oversimplified `S^4 = I`). The fourth-wave joins are shipped too: Milnor's exact
sequence `W(ℤ)→W(ℚ)→⊕_p W(F_p)` (`witt/milnor.rs::global_residues`, odd `p`), the named
Scharlau transfer (`trace_form::transfer_diagonal`), Nikulin's genus criterion
(`DiscriminantForm::is_isomorphic`) plus the theorem-1.10.1 existence predicate
(`nikulin_existence_report` / `nikulin_even_lattice_exists`), exact Conway-Sloane
`p`-adic genus symbols (`Genus::canonical_symbol_at`, with the corrected 2-adic
train/compartment/oddity reduction), the games↔integral lexicode edge
(`games/lexicode.rs`: greedy = mex, so the `[24,12,8]` lexicode is Golay), and the
Brown `ℤ/8` invariant — the char-2 cell of the mod-8 spine (`char2/brown.rs`:
`brown_f2`/`double_f2`, with `β = 4·Arf`, plus `DiscriminantForm::brown_invariant`
giving the float-free `β ≡ sign(L) mod 8` on 2-elementary discriminant forms). The
fifth-wave Bridge K is shipped too: the full `ℚ/ℤ` ungraded Brauer invariant
(`witt/cyclic.rs`: `BrauerClass` + `cyclic_algebra_invariant` = `v(a)/n` for the
unramified local cyclic class over the `Qq` leg) with full-strength reciprocity over
`F_q(t)` (`constant_extension_invariants`, `Σ_v deg(v)·v(a)/n = 0`); it lifts the
2-torsion `Brauer2Class` (which embeds as its `½`-slice) to the full local Brauer group.

The checked game-Clifford deformation surface is implemented as an engineering
bridge, not as a game-native scalar claim. `GameClifford::with_quadratic_data` accepts
integer `q`/polar data on a chosen game-generator tuple only after verifying every
game relation in the quotient is null and polar-radical for that data; over the
torsion-free target `ℤ`, relations such as `2* = 0` force `Q(*)` and all pairings
with `*` to vanish. The stronger question of a natural game-native source for such
quadratic data remains open in `docs/OPEN.md`.

The game-built Gold-form bridge is implemented, but the play rule is not. The
standard chain is:

```text
x + y      = XOR = disjunctive sum of impartial game values
x * y      = nim product = Turning-Corners product value
x -> x^2   = Frobenius = diagonal product x*x
Tr(x)      = x + x^2 + ... + x^(2^(m-1))
Q_a(x)     = Tr(x * x^(2^a))
```

Implemented probes verify Gold ranks, Arf zero-count bias, literal Turning-Corners
reconstruction on small fields, frame-obstruction experiments, misère-kernel
obstruction examples, loopy Draw/Loss-set experiments, and bent Gold-component route
probes. The conditional statement: if a game has P-positions `{Q = 0}`, Arf gives
the sign and size of the second-player win-bias. The existence of a non-tautological
natural rule with P-set `{Q = 0}` is open (`docs/OPEN.md`), but the σ-valued
echo-fifo+dummy realizer is **verified** (2026-06-10, adversarial review:
`experiments/echo_solver.py`, 391,680/391,680 m=8 checks, zero misses — record in
`writeups/goldarf.tex` §8); the open steps are recasting its forced-charge readout into
normal/misère/loopy outcome semantics and the general-n linking proof. The
realizer's *mechanism* is reduced (2026-06-10 second pass,
`experiments/linking_game.py`, goldarf §8 `sec:linking`): the σ-game is the
odd-close parity game on the support graph, and the linking theorem — an isolated
coin forces flips even, hence exactness for all m — is machine-verified on every
graph class k ≤ 7 with a strictly-verified two-mode defender strategy; only the
general-n induction is open.

Appendix-grade shipped layers that should not be mistaken for new Gold/Arf claims:
tropical thermography (`Semiring` + dual `Tropical<MaxPlus/MinPlus>`), the
source-verified ordinal nim Kummer tower below `ω^(ω^ω)`, the characteristic-2
Artin-Schreier local-global layer over `F_{2^m}(t)` including the Aravire-Jacob wild
summand, and the integral lattice/genus/mass/Leech/Niemeier/theta/code/Weil chain. These are
standard-math implementations and useful infrastructure; cite them as such.

## Commands

```sh
cargo test                                    # the math core (pure Rust, no Python)
cargo clippy --all-targets                    # lint (kept warning-clean)
cargo doc --no-deps                           # rustdoc (intra-doc links warning-clean)
cargo run --example tour                      # Rust demo
cargo run --example tropical                  # tropical-semiring / thermography demo
cargo run --example interactive_kernel        # open-problem probe
cargo run --example loopy_quadric             # open-problem probe
cargo run --example bent_route                # open-problem probe
python3 -m venv .venv && .venv/bin/pip install maturin
VIRTUAL_ENV=.venv .venv/bin/maturin develop   # build + install the abi3 extension
.venv/bin/python demo.py
.venv/bin/python experiments/framing_obstruction.py
.venv/bin/python experiments/gold_family_survey.py
.venv/bin/python experiments/misere_kernel.py
python3 experiments/echo_solver.py selftest   # echo adversarial-review harness (stdlib, no venv)
python3 experiments/linking_game.py all 5     # linking-reduction harness (stdlib, no venv; `all 7` ≈ 75 s)
python3 experiments/exception_column_m4.py    # 2·3^k excess column m=4 certification (stdlib, no venv; ≈ 2 min)
```

`maturin develop` needs `VIRTUAL_ENV` set (or a `.venv` in cwd) and `cargo` on PATH
(`. "$HOME/.cargo/env"`).

## Hard rules

1. **The math core is generic over `Scalar` and pure Rust.** PyO3 lives behind the
   `python` feature; never `use pyo3` outside `src/py/`, never make it non-optional.
   This is what keeps `cargo test` from linking libpython. (Details: `src/py/AGENTS.md`.)
2. **The metric carries `q` and `b` independently — do not collapse them.** In char
   2, `b` is alternating yet `q[i]` can be nonzero; collapsing them makes every
   char-2 algebra commutative. The optional `a[(i,j)]` lifts the engine to a general
   bilinear form. Build with `Metric::new`/`::diagonal`/`::grassmann`/`::general`,
   never the bare struct literal. (Details: `src/clifford/AGENTS.md`.)
3. **Signs go through the scalar's own `neg()`, never a literal `-1` or a
   `characteristic()` branch.** For nimbers `neg` is identity, so char-2 sign-
   vanishing falls out for free; hardcoding signs breaks char 2.
4. **Surreal arithmetic recurses only on exponents** (strictly simpler than the
   number). That is the entire termination argument — never recurse on the number.
5. **Per-backend, no mixing.** Each Python backend monomorphises the engine to one
   concrete scalar type; mixing worlds raises `TypeError` by construction. Intended.
6. **Verify, don't claim.** The `associativity_*` tests catch product bugs and
   `general_product_reproduces_reduce_word_when_a_empty` pins the general engine to
   an independent oracle. Add a test before trusting a new operation.
7. **Fixed-width integer payloads are `u128`/`i128`.** Public arithmetic carriers,
   invariant bits/residues, counts, budgets, and signed integral data use those
   widths consistently. `usize` remains for indexing, dimensions, and platform ABI
   hooks such as Python hashing; otherwise do not introduce narrower fixed-width
   carriers without a hard external signature forcing it.

## Style

- Rust 2021, `cargo fmt` clean, `cargo clippy --all-targets` warning-clean (the one
  crate-level allow — `needless_range_loop` for the matrix code — is justified in
  `lib.rs`; targeted `#[allow]`s carry a one-line reason). License: AGPL-3.0-or-later.
- Numeric payload style is deliberate: non-index fixed-width integers are
  `u128`/`i128` throughout the core, docs, examples, and tests.
- Display is deliberate and canonical (ogham Display v2, `docs/ogham/ogham.md` §9):
  blades render as wedge expressions `e0∧e1` (`∧` = U+2227); coefficients attach
  `coeff⋅label` (`⋅` = U+22C5) with coefficient-`1` elided and `-1` → `-label`
  (compared via `S::one().neg()`, never a literal). A term whose rendering starts
  with `-` joins with ` - ` (string-level, char-agnostic); the empty multivector
  renders as `S::zero()`'s display (`*0` in nim-worlds, `0` elsewhere). Nimbers
  print `*n`; ordinals are star-wrapped (`*5`, `*ω`, else `*(…)`: `*(ω + 1)`,
  `*(ω↑2)`, `*(ω⋅3)`, `*(ω↑(ω))`); surreals print CNF (`3⋅ω↑2 - ω + 5`, `ω↑(ω)`,
  `ω↑-1`, `ω↑(1/2)` — exponent bare iff a signed integer); `Fpn` `3⋅x↑2 + 2⋅x + 1`;
  `Poly` uses variable `t` (`1 + 2⋅t`, parens on a non-atomic coefficient);
  `RationalFunction` `(num)/(den)`. Atomic = no spaces and no `⋅ ∧ ↑ / + -`
  outside balanced parens.
- Rust scalar operators: total-product backends have `+ - *` and unary `-`
  (concrete-only, via `impl_scalar_ops!`), plus `^ u128` for power (`x ^ 3` =
  square-and-multiply via `Scalar::mul`; `x ^ 0 == one()`). The `u128` RHS
  prevents element-element `^` from compiling — on `Nimber`, `x ^ x` would silently
  mean nim-addition (XOR), so no `BitXor<Self>` impl exists on any backend. **Rust
  `^` binds looser than `*`; parenthesize when mixing product and power.**
  `Multivector` has `&` (wedge, ogham `∧`) via `impl BitAnd`; **no `^` operator on
  `Multivector`** — the geometric product needs the metric, so use
  `CliffordAlgebra::pow(&self, v, k)` for repeated geometric multiplication.
  `Ordinal` deliberately omits owned `*` and `^` because transfinite
  nim-multiplication is partial at the verified Kummer boundary; use `nim_mul` and
  `nim_pow` for the checked paths. Generic engine code over `S: Scalar` still calls
  `.add(&x)`/`.mul(&x)` — operators are NOT a supertrait (see `src/scalar/AGENTS.md`).

## Testing

`cargo test` is the source of truth and needs no Python. It does **NOT** compile the
`python` feature — after touching `src/py/` or any core API the bindings call, run
`cargo check --features python` (and `cargo clippy --features python --all-targets`).
After touching `clifford/` or `scalar/big/surreal/`, also rebuild + run `demo.py`
(display changes don't surface in `cargo test`).

After touching any doc comment (`//!` / `///`), run `cargo doc --no-deps` and keep it
warning-clean — intra-doc links to renamed/private items and `[i]`-style brackets in
prose are the common breakages. **Run it cold** (`rm -rf target/doc` first, or
`RUSTDOCFLAGS="-D warnings" cargo doc --no-deps`): an incremental `cargo doc` only
re-checks what it recompiles, so it silently under-reports stale-link warnings in
untouched modules.

Two property suites (dev-dep `proptest`, in `tests/`): `scalar_axioms.rs` fuzzes the
commutative-ring axioms across every backend, `clifford_axioms.rs` fuzzes geometric-
product associativity/distributivity over random metrics in char 0 and char 2. The
capped-relative precision models (Qp/Qq/Laurent/Ramified/Gauss/Adele) are excluded
from the exact-ring fuzz by design; `ExactScalar`/`ExactFieldScalar`/`PrecisionScalar`
mark that boundary without becoming `Scalar` supertraits. (serde is intentionally NOT
shipped — the invariant-carrying types need custom deserialization, not a naive
derive.)

The narrow Gold/Arf game thread and the genuine open problems live in `docs/OPEN.md`; the
draft notes are `writeups/goldarf.tex` (Gold/Arf) and `writeups/excess.tex`
(transfinite excess). Read `docs/OPEN.md` before touching `forms/char2/`,
`forms/quadric_fit.rs`, `forms/char0.rs`, `games/coin_turning.rs`, `games/kernel.rs`,
`games/misere.rs`, `games/loopy/`, `forms/witt/`, `experiments/`, or the
open-question example probes.
