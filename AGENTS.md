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
lib), `demo.py` (the Python tour), and `NOTES.md` (the mathematical thread).

## The symmetries (the intellectual spine — see README)

The two pillars are complementary views of the *same* table of numbers:
`scalar/` groups **by place** (each field beside its ring of integers); `forms/`
cuts ACROSS by **characteristic** (the one classification theory, three ways). The
recurring symmetries the project is built around:

- **char 0 ↔ char 2** classifiers (the real 8-fold table mirrored by Arf/Brauer–Wall).
- **surreal No ↔ ordinal On₂** (char-0 field and char-2 non-field sharing one CNF core).
- **(field, ring of integers)** pairings, made structural in `scalar/integrality.rs`.
- the **2×2 functor table** (algebraic|transcendental × residue|value-extending).
- **local ↔ global** (the Springer trio + the local-global Hasse–Minkowski layer).

## Commands

```sh
cargo test                                    # the math core (pure Rust, no Python)
cargo clippy --all-targets                    # lint (kept warning-clean)
cargo run --example tour                      # Rust demo
python3 -m venv .venv && .venv/bin/pip install maturin
VIRTUAL_ENV=.venv .venv/bin/maturin develop   # build + install the abi3 extension
.venv/bin/python demo.py
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

## Style

- Rust 2021, `cargo fmt` clean, `cargo clippy --all-targets` warning-clean (the one
  crate-level allow — `needless_range_loop` for the matrix code — is justified in
  `lib.rs`; targeted `#[allow]`s carry a one-line reason). License: see `LICENSE`.
- Display is deliberate: blades render `e0e1`; coefficients `1`/`-1` elided; nimbers
  print `*n`; surreals print CNF (`3ω^2 - ω + 5`, `ω^(ω)`, `ω^-1`).
- Rust scalar operators: every backend has `+ - *` and unary `-` (concrete-only, via
  `impl_scalar_ops!`). Generic engine code over `S: Scalar` still calls `.add(&x)`/
  `.mul(&x)` — operators are NOT a supertrait (see `src/scalar/AGENTS.md`).

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

The narrow math thread (Arf↔Clifford, the games bridge, the char-0/char-2 classifier
symmetry, the open play-semantics question) is written up in `NOTES.md` — read it
before touching `forms/char2/`, `forms/quadric_fit.rs`, `forms/char0.rs`,
`games/coin_turning.rs`, `games/kernel.rs`, `games/misere.rs`, `games/loopy.rs`,
`forms/witt.rs`, `experiments/`, or the open-question example probes.
