# Contributing

Thanks for the interest. ogdoad is a research codebase, so the bar is correctness
first: a new operation lands with a test that pins it to an independent oracle, not
on a "looks right" basis.

## Read the working notes first

`AGENTS.md` is the map: the four pillars (`scalar/`, `clifford/`, `forms/`,
`games/`) plus the PyO3 bindings, and each pillar has its own `AGENTS.md` with the
file-by-file breakdown and the layer-specific "things that look like bugs but
aren't". `docs/OPEN.md` is the genuine open problems — read it before touching
`forms/char2/`, `games/`, the `experiments/`, or the open-question example probes,
so you don't file a research question as a bug or a solved theorem.

## The non-negotiables

These are the invariants the whole thing rests on (full list in AGENTS.md → Hard
rules):

- **The math core is generic over `Scalar` and pure Rust.** PyO3 lives behind the
  `python` feature — never `use pyo3` outside `src/py/`, never make it
  non-optional. This is what keeps `cargo test` from linking libpython.
- **The metric carries `q` and `b` independently — do not collapse them.** In
  characteristic 2 the polar form `b` is alternating yet `q[i]` can be nonzero;
  collapsing them makes every char-2 algebra commutative (the wrong object).
- **Signs go through the scalar's own `neg()`**, never a literal `-1` or a
  `characteristic()` branch — for nimbers `neg` is identity, so char-2
  sign-vanishing falls out for free.
- **Surreal arithmetic recurses only on exponents** (strictly simpler than the
  number). That's the entire termination argument.
- **Verify, don't claim.** Add a test before trusting a new operation.

## Test plan

```sh
cargo test                                  # the math core — source of truth, no Python
cargo clippy --all-targets                  # kept warning-clean
cargo fmt --check
RUSTDOCFLAGS="-D warnings" cargo doc --no-deps   # run COLD (rm -rf target/doc first)
```

`cargo test` does **not** compile the `python` feature. After touching `src/py/` or
any core API the bindings call:

```sh
cargo check --features python
cargo clippy --features python --all-targets
```

After touching `clifford/` or `scalar/big/surreal/`, rebuild and run the tour —
Display changes (`e0e1`, `*n`, CNF) don't surface in `cargo test`:

```sh
VIRTUAL_ENV=.venv .venv/bin/maturin develop
.venv/bin/python demo.py
```

## Claim levels

When you change prose, comments, examples, or the writeup, label the claim:
**standard math** (external fact) · **implemented and tested** (backed by this
checkout) · **interpretation** (a conditional bridge) · **open** (lives in
`docs/OPEN.md`). A new "X is true" statement is backed by a test or a citation, not
asserted.

## Releasing

The version in `Cargo.toml` is the single source of truth (pyproject and the
maturin build inherit it). The release workflow is **dormant** while the version is
`0.0.0`; bumping it arms the pipeline, which on the next push to `main` publishes to
crates.io and PyPI (both via OIDC trusted publishing), tags `vX.Y.Z`, and cuts a
GitHub release. Each target is checked independently, so a partial-failure run
resumes cleanly.

## License

By contributing you agree your contributions are licensed under AGPL-3.0-or-later,
same as the rest of the project.
