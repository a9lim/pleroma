# AGENTS.md — pleroma

Working notes for agents editing this repo. Global rules still apply. This file is
the **summary**; each `src/` pillar has its own `AGENTS.md` with the detailed
file-by-file breakdown and the layer-specific "things that look like bugs".

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
(`scalar::Nimber`, `clifford::sandwich`, …).

| dir | pillar | detail |
|---|---|---|
| `src/scalar/`   | commutative coefficient worlds (the `Scalar` trait + every backend) | [`src/scalar/AGENTS.md`](src/scalar/AGENTS.md) |
| `src/clifford/` | the multivector engine + the GA layer | [`src/clifford/AGENTS.md`](src/clifford/AGENTS.md) |
| `src/forms/`    | quadratic forms & invariants, by the char trichotomy plus local-global and integral layers | [`src/forms/AGENTS.md`](src/forms/AGENTS.md) |
| `src/games/`    | combinatorial game theory | [`src/games/AGENTS.md`](src/games/AGENTS.md) |
| `src/py/`       | PyO3 bindings (feature = "python") + the binding-scope policy | [`src/py/AGENTS.md`](src/py/AGENTS.md) |
| `src/linalg/`   | crate-private shared linear algebra | [`src/linalg/AGENTS.md`](src/linalg/AGENTS.md) |

Beyond the library: `examples/` (Rust demos + the open-question probes:
`interactive_kernel`, `octal_hunt`, `loopy_quadric`, `misere_quotient`,
`bent_route`, `tour`), `experiments/` (Python research probes on top of the shipped
lib), `demo.py` (the Python tour), `OPEN.md` (the genuine research problems),
`ROADMAP.md` (the implemented cross-pillar bridge map, a proposed second-wave bridge
spec, plus remaining boundaries), and
`writeups/goldarf.tex` (the narrow draft note on the Gold/Arf game thread).

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
- **Open**: the natural Gold-quadric game rule, genuine game-native quadratic
  deformation of `GameExterior`, and transfinite nim multiplication beyond the
  source-verified DiMuro excess table. These live in `OPEN.md`.

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
| `Ordinal` | staged transfinite nimbers; now a checked/panic-on-escape `Scalar` for Clifford metrics; nim-addition on represented CNF terms and nim-multiplication through source-verified Kummer carries `alpha_u`, `u <= 43` |

The char-2 Clifford point is load-bearing. In characteristic 2, `q` and `b` are
independent:

```text
e_i^2              = q_i
e_i e_j + e_j e_i = b_ij
```

The polar form is alternating, so `b_ii = 0`, while `q_i` can be nonzero. If the
engine collapses `q` and `b`, the char-2 Clifford product becomes the wrong
commutative object. The implemented char-2 form layer reports `ArfResult { arf:
u128, rank, radical_dim, radical_anisotropic, o_type }`; for degenerate forms, Arf
of the nonsingular core is not the whole form. On nonsingular metrics over `Nimber`,
supported `Fpn<2,N>`, and the documented finite ordinal windows, the same Arf bit is
also the implemented Brauer-Wall class `BW(F_{2^m}) ~= Z/2`; hyperbolic planes are
`0`, the anisotropic plane is `1`, and direct sum / graded tensor adds by XOR.
`clifford::spinor` has a separate char-2 route: no `1/2(1+w)`, honest blade
idempotents such as `e_i e_j` when they shrink a left ideal, and otherwise the full
regular/lazy left action. Singular polar forms and general-bilinear `a` metrics are
rejected.

The cross-pillar bridge pass in `ROADMAP.md` is implemented in the Rust core:
`IntegralForm` now exports rational and even-mod-2 Clifford metrics plus
discriminant Gauss-sum/Milgram checks; finite char-2 `Fpn<2,N>` classification is
wired through the façade; cyclic Galois/Frobenius maps have Clifford linear-map
constructors; and `Ordinal` can serve as a Clifford scalar inside the verified
Kummer boundary.

The second-wave bridges E/H/I in `ROADMAP.md` are also implemented in the Rust core.
`forms/integral/codes.rs` carries binary codes, MacWilliams, and Construction A
(with the required `1/sqrt(2)` scaling and an `Option` boundary when the scaled Gram
is not integral), including the Type II length-16 code whose lattice is `D16+`.
`forms/integral/theta.rs` and `modular.rs` give exact theta coefficients and
`E4`/`E6` q-expansion identification, pinning `theta_E8 = E4`,
`theta_{E8+E8} = theta_{D16+} = E4^2`, and the rootless Leech `q^1` oracle.
`DiscriminantForm` now exposes dependency-free `Complex64` Weil `S`/`T` matrices;
the `S` prefactor is the conjugate of the positive Milgram phase, and
`verify_weil_relations` checks the honest metaplectic relations rather than the
oversimplified `S^4 = I` statement.

The game-built Gold-form bridge is implemented, but the play rule is not. The
standard chain is:

```text
x + y      = XOR = disjunctive sum of impartial game values
x * y      = nim product = Turning-Corners product value
x -> x^2   = Frobenius = diagonal product x*x
Tr(x)      = x + x^2 + ... + x^(2^(m-1))
Q_a(x)     = Tr(x * x^(2^a))
```

Implemented probes verify Gold ranks, Arf zero-count bias, literal
Turning-Corners reconstruction on small fields, frame-obstruction experiments,
misere-kernel obstruction examples, loopy Draw/Loss-set experiments, and bent
Gold-component route probes. The conditional statement is: if a game has
P-positions `{Q = 0}`, Arf gives the sign and size of the second-player win-bias.
The existence of a non-tautological natural rule is open; details live in `OPEN.md`.

Appendix-grade shipped layers that should not be mistaken for new Gold/Arf claims:
tropical thermography (`Semiring` + dual `Tropical<MaxPlus/MinPlus>`), the
source-verified ordinal nim Kummer tower below `omega^(omega^omega)` with `u <= 43`
excesses, the characteristic-2 Artin-Schreier local-global layer over `F_{2^m}(t)`
including the Aravire-Jacob wild summand, and the integral
lattice/genus/mass/Leech/theta/code/Weil chain. These are standard-math
implementations and useful infrastructure; cite them as such.

## Commands

```sh
cargo test                                    # the math core (pure Rust, no Python)
cargo clippy --all-targets                    # lint (kept warning-clean)
cargo run --example tour                      # Rust demo
cargo run --example interactive_kernel        # open-problem probe
cargo run --example loopy_quadric             # open-problem probe
cargo run --example bent_route                # open-problem probe
python3 -m venv .venv && .venv/bin/pip install maturin
VIRTUAL_ENV=.venv .venv/bin/maturin develop   # build + install the abi3 extension
.venv/bin/python demo.py
.venv/bin/python experiments/framing_obstruction.py
.venv/bin/python experiments/gold_family_survey.py
.venv/bin/python experiments/misere_kernel.py
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
  `lib.rs`; targeted `#[allow]`s carry a one-line reason). License: see `LICENSE`.
- Numeric payload style is deliberate: non-index fixed-width integers are
  `u128`/`i128` throughout the core, docs, examples, and tests.
- Display is deliberate: blades render `e0e1`; coefficients `1`/`-1` elided; nimbers
  print `*n`; surreals print CNF (`3ω^2 - ω + 5`, `ω^(ω)`, `ω^-1`).
- Rust scalar operators: total-product backends have `+ - *` and unary `-`
  (concrete-only, via `impl_scalar_ops!`). `Ordinal` deliberately omits owned `*`
  because transfinite nim-multiplication is partial at the verified Kummer boundary;
  use `nim_mul` for the checked path. Generic engine code over `S: Scalar` still
  calls `.add(&x)`/`.mul(&x)` — operators are NOT a supertrait (see
  `src/scalar/AGENTS.md`).

## Testing

`cargo test` is the source of truth and needs no Python. It does **NOT** compile the
`python` feature — after touching `src/py/` or any core API the bindings call, run
`cargo check --features python` (and `cargo clippy --features python --all-targets`).
After touching `clifford/` or `scalar/big/surreal/`, also rebuild + run `demo.py`
(display changes don't surface in `cargo test`).

Two property suites (dev-dep `proptest`, in `tests/`): `scalar_axioms.rs` fuzzes the
commutative-ring axioms across every backend, `clifford_axioms.rs` fuzzes geometric-
product associativity/distributivity over random metrics in char 0 and char 2. The
capped-relative precision models (Qp/Qq/Laurent/Ramified/Gauss/Adele) are excluded
from the exact-ring fuzz by design; `ExactScalar`/`ExactFieldScalar`/`PrecisionScalar`
mark that boundary without becoming `Scalar` supertraits. (serde is intentionally NOT
shipped — the invariant-carrying types need custom deserialization, not a naive
derive.)

The narrow Gold/Arf game thread and the genuine open problems now live in
`OPEN.md`; the current draft note is `writeups/goldarf.tex`. Read `OPEN.md` before
touching `forms/char2/`, `forms/quadric_fit.rs`, `forms/char0.rs`,
`games/coin_turning.rs`, `games/kernel.rs`, `games/misere.rs`, `games/loopy.rs`,
`forms/witt/`, `experiments/`, or the open-question example probes.
