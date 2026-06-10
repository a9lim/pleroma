# AUDIT.md — mathematical correctness audit

**Date:** 2026-06-09 · **Tree:** `main` @ `0bbaec6` (clean) · **Baseline:** `cargo test` fully green

A comprehensive mathematical audit of the repository, performed by a 16-agent
review fleet with per-finding adversarial verification. **52 findings confirmed**
(3 critical, 10 major, 25 minor, 14 doc), 4 raw findings refuted in verification,
2 areas left unaudited (see [Coverage](#coverage-map)).

Severity vocabulary (as used throughout):

- **critical** — wrong results from core operations on mainline inputs
- **major** — wrong results on real but narrower inputs, or a wrong invariant/classification
- **minor** — edge-case wrongness, recoverable, or panic-instead-of-answer
- **doc** — false mathematical claim in prose/comments only; code unaffected

---

## Implementation progress

This section tracks fixes landed after the original audit. The original finding
writeups are retained below as the defect record and repro context.

| Batch | Status | Findings | Current evidence |
|---|---|---|---|
| 1 | fixed and validated | L-1, C-2, F-3, F-1, G-1, G-2, G-3, G-5, P-5 | Unit-pivot nullspace now returns `None` on nonunit pivots and callers propagate it; Hermitian off-diagonal pivot uses the conjugate partner; loopy stopper sums/order and misere-Nim zero heaps are corrected. Validation: `cargo test`, `cargo check --all-targets`, `cargo check --features python --all-targets`, `cargo check --examples`, both clippy `-D warnings` gates, `maturin develop`, a targeted Python probe, and `git diff --check`. |
| 2 | fixed and validated | F-7, F-8 | 2-adic Jordan splitting now prefers an odd block when a diagonal entry ties the minimal valuation, and 2-adic trains now continue across one empty scale between type-I constituents. Added the audit's `[[2,1],[1,1]] ~ Z^2` and `diag(1,20) ~ diag(5,4)` counterexamples. Validation: `cargo test`, `cargo check --all-targets`, `cargo check --features python --all-targets`, `cargo check --examples`, both clippy `-D warnings` gates, and `git diff --check`. |
| 3 | fixed and validated | C-1, C-4 | Reversion now reduces the reversed generator word through the Clifford product, so non-orthogonal metrics get true Clifford reversion instead of the exterior grade sign shortcut; versor inverse/norm/classification consumers now use that route. Added non-orthogonal anti-automorphism and `q=[1,1], b01=1` rotor norm-3 regressions. Validation: `cargo test`, `cargo check --all-targets`, `cargo check --features python --all-targets`, `cargo check --examples`, both clippy `-D warnings` gates, and `git diff --check`. |
| 4 | fixed and validated | S-1, S-2, S-3 | Surreal lazy inverse/root series now compute exact finite powers until the requested leading window is stable instead of trusting a fixed `2n+8` truncation; binomial coefficients and series merges are checked so unrepresented deep windows return `None` rather than panic or emit guessed terms. `Rational` addition/multiplication now reduce before multiplying to avoid avoidable pre-normalization overflow. Added dense inverse, deep square-root cancellation, and `ω+1`/`1+ω^-1` exact-square refusal regressions. Validation: `cargo test`, `cargo check --all-targets`, `cargo check --features python --all-targets`, `cargo check --examples`, both clippy `-D warnings` gates, `maturin develop`, `.venv/bin/python demo.py`, and `git diff --check`. |
| 5 | fixed and validated | F-10 | Finite-place char-2 Springer Laurent expansion now evaluates the P-free unit through the shared Hensel parameter `T(u)` from `function_field_char2`, so degree-`>1` places keep κ[[u]] carries instead of treating polynomial digits as coefficients. Local square testing now uses valuation parity plus the char-2 derivative-zero criterion. Added the audit's `P=t^2+t+1`, `℘(t/P)=(t^3+t)/P^2` regression, pinning `[1, ℘(t/P)]` as locally hyperbolic. Validation: `cargo test`, `cargo check --all-targets`, `cargo check --features python --all-targets`, `cargo check --examples`, both clippy `-D warnings` gates, and `git diff --check`. |
| 6 | fixed and validated | F-5, F-6, F-12 | Char-2 `WittClass`, `WittClassG`, and `BrauerWallClass` now carry the finite field degree and reject cross-field addition instead of XORing bare Arf bits. `Fpn<2,N>`, nimber, ordinal, and Python façade paths preserve the tag; Python constructors keep degree-1 defaults and expose `field_degree`. Added the audit's F2/F4 mixed-field regression, including the re-evaluated direct sum with Arf 1. The Arf and Witt prose now states the fixed-rank graded/field-fixed boundary and the trivial algebraically closed On₂ case. Validation: `cargo test`, `cargo check --all-targets`, `cargo check --features python --all-targets`, `cargo check --examples`, both clippy `-D warnings` gates, `maturin develop`, `.venv/bin/python demo.py`, targeted Python cross-field rejection probe, and `git diff --check`. |
| 7 | fixed and validated | S-4, S-5, S-6, S-7, S-8, S-9, S-10, S-11, S-12, S-13, S-14, S-15 | Scalar boundary arithmetic now fails loudly or reduces safely instead of wrapping: ordinary ordinal Cantor coefficients use checked arithmetic; `mod_inverse_u128`, residue-field modular powers, and large-modulus `WittVec` arithmetic use u128-native modular helpers; `Zp`/`Qp`/`LocalQp` reject moduli beyond the i128-backed embedding boundary; `Qp::to_integer` uses modular multiplication; and `Laurent<S,0>` is rejected. The Witt-vector component prose now names the Teichmüller-digit convention, Gauss residue lifting documents and tests the mixed-characteristic non-multiplicative boundary, and the primitive-density / exact-table docs are corrected. Validation: scalar-focused `cargo test scalar:: --lib`, full `cargo test`, `cargo check --all-targets`, `cargo check --features python --all-targets`, `cargo check --examples`, both clippy `-D warnings` gates, `maturin develop`, `.venv/bin/python demo.py`, and `git diff --check`. |
| 8 | fixed and validated | C-3, C-5, C-6, C-7, C-8, F-2, F-4, F-9, F-11, G-4, G-6, G-7, L-2, L-3 | High-dimension Clifford helpers now use checked basis sizing: general inverse refuses non-scalars beyond the materializable matrix range while still inverting scalars, and outermorphism grade masks cover the full 128-bit basis window. Char-0 Clifford matrix dimensions are `u128`, so dim-128 tables report `2^64` instead of wrapping; function-field residue orders use checked exponentiation. Integer row-lattice arithmetic now uses the module's checked helpers at the remaining overflow sites. Explicit game-relation certificates compute row independence, and the remaining prose/doc false claims were narrowed: graded tensor products for radical factors, Cayley chart boundaries, spinor-norm exactness, nilpotent-exp cap refusal, A1 automorphisms, Hackenbush sign expansion, Turning-Corners Kummer boundary, and property-suite `a`-term coverage. Validation: `cargo test --lib`, full `cargo test`, `cargo check --all-targets`, `cargo check --features python --all-targets`, `cargo check --examples`, both clippy `-D warnings` gates, `maturin develop`, `.venv/bin/python demo.py`, and `git diff --check`. |

---

## Methodology

- **Partition.** Sixteen reviewers, one per domain slice: nimber+ordinal,
  surreal/omnific/surcomplex, p-adic/Witt/functor/adele, scalar core + finite
  fields, clifford engine, clifford GA structures, forms char-0/odd-char, forms
  char-2/Witt, forms local–global/Springer, forms integral, games, linalg +
  property suites, py scalars/engine, py forms/games/catalog, docs/prose, and
  examples/experiments. Every reviewer was required to read the root and pillar
  `AGENTS.md` first, so documented-intentional behavior (nimber `neg()` ≡ id,
  independent `q`/`b` in char 2, precision-capped backends, panic/`Option`
  scope boundaries, the Kummer window) was excluded by construction.
- **Independent recomputation.** Reviewers were required to verify anything
  reported with confidence > 0.5 by hand computation, python3 brute force
  (nim tables, zero counts, theta/Eisenstein coefficients, retrograde game
  analysis, Hilbert symbols, Witt addition laws), or scratch crates in `/tmp`
  compiled against this repo. No repo files were modified.
- **Adversarial verification.** Each raw finding then faced two independent
  skeptics instructed to *refute* it (defaulting to refuted when uncertain):
  one recomputing the mathematics from scratch, one checking the docs/tests
  for documented intent or misreading. Only double-confirmed findings appear
  below; the four that died are in [Refuted](#refuted-in-verification).
- **Scale.** 128 agents, ~9.7M tokens, 1,581 tool calls, ~50 min wall clock.
- **A meta-observation up front:** several confirmed bugs are *pinned by
  passing tests* (the loopy-game sum table, the `clifford_axioms.rs` header
  claim, the genus shear test that never disturbs scan order). A green suite
  was treated as evidence of nothing beyond what each assertion literally says.

---

## Headline findings

The three criticals, plus the cross-cutting root causes that account for
about half of everything else.

| # | Severity | Where | What |
|---|---|---|---|
| F-7 | critical | `forms/integral/genus.rs` | 2-adic Jordan splitting misclassifies odd blocks as type II on a valuation tie — `are_in_same_genus` returns `false` for ℤ-isometric lattices |
| F-10 | critical | `forms/springer/char2.rs` | char-2 local engine drops multiplication/squaring carries at every finite place of degree ≥ 2 |
| S-3 | critical | `scalar/big/surreal/analytic.rs` | `is_square`/`sqrt`/Surreal form classification panic with i128 overflow on mainline inputs (`ω+1`) |

**Recurring root causes:**

1. **`reverse()` is not the reversion anti-automorphism on non-orthogonal
   metrics** (C-1, C-4). The engine's per-blade grade sign equals true
   reversion only when `b = a = ∅` — and every nonsingular char-2 metric (the
   flagship nim-Clifford case) has nonzero `b`. Everything downstream of
   `reverse` (`norm2`, `versor_inverse`, `sandwich`, `spinor_norm`,
   `classify_versor`) is wrong or falsely refuses on such metrics.
2. **`linalg::field::unit_pivot_nullspace` has no honest-failure path over
   rings** (L-1, C-2, F-5). Columns with no *unit* pivot are silently treated
   as free, so "kernels" contain non-kernel vectors over `Integer`/`Omnific`.
   The sibling routines (`solve`, `inverse_matrix`) return `None` in exactly
   this situation; this one fabricates an answer.
3. **Fixed-window truncation in the surreal lazy-series layer** (S-1, S-2,
   S-3). The `w = 2n+8` working window is treated as guaranteeing *n exact
   leading terms*, but Hahn-series cancellation makes the required depth
   input-dependent, and i128 `Rational` coefficients cannot hold deep binomial
   coefficients. One root design issue, three findings.
4. **Unchecked u128/i128 arithmetic at the representable boundary** (~9
   minors). `Cargo.toml` carries no `[profile]` override, so release builds
   (including `maturin --release` wheels) *wrap silently* where debug builds
   panic — in Cantor ordinal coefficients, WittVec products, `mod_pow`,
   `Qp::to_integer`, residue-field orders `q^deg`, SNF corner cases. The
   affected modules' own siblings (`tower.rs`, `qp.rs`, `fp.rs`,
   `linalg/integer.rs`) demonstrate the intended checked-arithmetic
   discipline; these sites missed it. A single `overflow-checks = true` in
   the release profile would convert every one of these from silent
   wrongness to a loud panic at near-zero cost for this workload.
5. **The loopy-game catalogue tables were filled from intuition rather than
   the stopper survival criterion** (G-1, G-2, G-3), and the inline tests pin
   the wrong table.

---

## Findings — `src/scalar/`

### S-1 · MAJOR · `scalar/big/surreal/analytic.rs` — `inv_to_terms` returns spurious terms inside the claimed *n* leading terms

`inv_to_terms`, lines 26–56. Doc contract: "the n leading terms of 1/x." The
Neumann loop truncates each power `(-r)^k` and the running series to
`w = 2n+8` terms, silently dropping contributions the true leading terms need
whenever intermediate partial sums are denser than the final answer.
Verified in real Rust (path-dep crate): for `x = 1 + ω⁻¹ + … + ω⁻²⁰`, the
true inverse is `1 − ω⁻¹ + ω⁻²¹ − ω⁻²² + …`, but `inv_to_terms(3)` returns
`1 − ω⁻¹ + ω⁻¹⁴` — the third term has the wrong exponent and a coefficient
whose true value in `1/x` is 0. Not truncation-below-precision: a wrong
coefficient *within* the claimed window. (A faithful Python mirror matches
exact series inversion on 300 random dense polynomials — mainline inputs are
fine; the failure needs cancellation deeper than the window.)

### S-2 · MAJOR · `scalar/big/surreal/analytic.rs` — `sqrt_to_terms` / `nth_root_to_terms`: same fixed-window flaw through `binomial_series`

Lines 106–133. For `y = 1 − ω⁻¹ + ω⁻²⁵` (so `z = y²` is finite-support and
an exact perfect square), `z.sqrt_to_terms(3)` returns
`1 − ω⁻¹ − (52003/8388608)·ω⁻²⁶` instead of `1 − ω⁻¹ + ω⁻²⁵` — wrong
exponent, junk coefficient. Cube-root analogue fails the same way. Gap ≤ ~22
works; gap ≥ 25 fails at n = 3.

### S-3 · CRITICAL · `scalar/big/surreal/analytic.rs` — `is_square` / exact `sqrt` / Surreal form classification panic with i128 overflow on mainline inputs

`binomial_series` accumulates `binom(α, j)` as a reduced i128 `Rational`
(denominators grow like `4^j` for α = 1/2) with the coefficient update *before*
any break check; `ExactRoots::sqrt` retries with windows up to ~104, and the
coefficients overflow i128 around j ≈ 60–65. Verified panics end-to-end:
`is_square(ω+1)` panics (`rational.rs:198`) where the documented answer is
`false`; `Surreal::witt_class(Metric::diagonal(vec![ω+1]))` panics where the
documented behavior (AGENTS: "classification only on represented exact square
classes") is `None`; `is_square(1+ε)` panics; and `is_square(y²)` for
`y = 1 − ω⁻¹ + ω⁻²⁵` panics even though `y²` *is* a represented exact square
(expected `true`). `ω+1` is about the most natural Surreal metric entry after
`ω` itself, hence critical.

*Fix shape for S-1/2/3 (per the reviewer):* cancellation-aware adaptive
window with bigint coefficients, or an honest weakening of the doc contracts
to best-effort plus a `None` on window exhaustion.

### S-4 · MINOR · `scalar/big/ordinal/cantor.rs` — unchecked u128 coefficient arithmetic in Cantor `ord_add`/`ord_mul`

`ord_add` line 43 (`+`), `ord_mul` line 66 (`*`). Release builds wrap mod
2^128: `monomial(1, 2^127).ord_mul(fin(4))` yields the CNF term `ω·0` — a
zero coefficient violating the type's own invariant (mod.rs:67–68), with
`is_zero()` false. Verified in a release-mode scratch crate; debug panics, so
`cargo test` can never see it. Inconsistent with the same pillar's `tower.rs`,
which returns `None` on every analogous overflow. Exposed to Python via
`PyOrdinal.ord_add/ord_mul` (release wheels only; `maturin develop` is debug).

### S-5 · DOC · `scalar/finite_field/nimber/galois.rs` — primitive-element density misstated

`nim_primitive_element` doc, lines 81–84: claims φ(2¹²⁸−1)/(2¹²⁸−1) ≈ 0.30.
Exact computation over the file's own nine `ORDER_FACTORS` (whose product was
verified to equal 2¹²⁸−1, all factors prime): ∏(1−1/p) = 0.4992… ≈ **0.50**.
The function and its "returns quickly" conclusion are correct — more strongly
so, in fact.

### S-6 · MINOR · `scalar/finite_field/wittvec.rs` — `witt_components` are Teichmüller digits, not classical Witt coordinates

Module doc (23–32) claims the coordinates are "the genuine Witt/Teichmüller
coordinates" satisfying the classical p = 2 carry `S₁ = x₁ + y₁ − x₀y₀`. The
code defines digits via `a = Σ τ(x_i)·pⁱ`; the classical Witt isomorphism
(Serre, *Local Fields* II §6) is `a = Σ τ(a_i^{p^{−i}})·pⁱ` — the code's
digits are the Frobenius-twisted coordinates, coinciding only over F_p
(exactly the only case the test covers). Verified over W₂(F₄): code addition
satisfies the twisted law `z₁ = x₁ + y₁ + √(x₀y₀)`, not the classical
`x₀·y₀` carry; counterexample with `x₀ = 1, y₀ = t̄` gives `√t̄ = 1+t̄ ≠ t̄`.

### S-7 · MINOR · `scalar/functor/gauss.rs` — `Gauss::teichmuller` violates the documented multiplicativity contract over mixed-characteristic bases

Lines 266–282. `ResidueField::teichmuller` is documented as a multiplicative
section, but the Gauss impl lifts residues *coefficientwise* through the base
section — multiplicative only when the base section is a ring homomorphism
(equal characteristic). Over `Gauss<Qp<5,6>>` with `r = 1+t̄`, `s = 1+2t̄`:
`τ(rs)` has t-coefficient `τ(3) = 1068` while `τ(r)τ(s)` has
`τ(1)+τ(2) = 14558`; difference has valuation 1, well inside precision 6.
`residue(τ(a)) = a` does hold.

### S-8 · MINOR · `scalar/finite_field/wittvec.rs` — unchecked u128 products in `ring_mul`/`add` for p^N > 2^64

`ring_mul` (92–123) computes `ai·b[j]` with operands < m = p^N, so products
reach ~m² and overflow once p^N > 2^64 (legal const parameters, e.g.
`WittVec<2,65,1>`; `modulus()` accepts anything ≤ u128::MAX). Release wraps
silently; siblings Zp/Qp/LocalQp use `checked_*` precisely to avoid this.
Currently unreachable from the Python bindings (max exposed modulus 5⁴).

### S-9 · MINOR · `scalar/small/zp.rs` (+ `qp.rs`, `local_qp.rs`) — moduli in (i128::MAX, u128::MAX] accepted but break embedding and inversion

`Zp::new`, `Qp::from_i128`, `LocalQp::from_i128` cast the modulus
`as i128`; for p^K > 2¹²⁷ (constructible: `Qp<2,127>`, `Qp<5,55>`) the cast
wraps negative and reduction produces out-of-range mantissas (e.g.
`Qp::<2,127>::from_i128(3)` → unit `2¹²⁷+3`, violating `unit ∈ [0, p^K)` and
breaking `PartialEq`/`Hash`). Independently, `inv()` of genuine units returns
`None` on this band (see S-13). Fix: assert `modulus ≤ i128::MAX` in
`assert_supported_*`, or a u128-native reduction.

### S-10 · MINOR · `scalar/functor/laurent.rs` — `Laurent<S, 0>` unguarded: `one()*one() = 0` and `inv()` panics

The p-adic backends assert K > 0; Laurent never validates. At K = 0,
constructors bypass `normalized()` while arithmetic truncates to length 0:
`one().mul(&one())` = 0 and `one().inv()` indexes `w[0]` of an empty vec
(out-of-bounds panic). A degenerate parameter yields an internally
inconsistent non-ring instead of a clean rejection.

### S-11 · MINOR · `scalar/integrality.rs` — `Qp::to_integer` uses `wrapping_mul`: wrong Zp residue when P^(K+1) > 2^128

Lines 129–143. For valid parameters with P^K < 2^128 < P^(K+1) (e.g.
`Qp<3,80>`), `acc·P` wraps and `(x mod 2^128) mod m ≠ x mod m` for odd P.
Verified: unit `3⁸⁰−1`, v = 1 → correct residue `3⁸⁰−3 =
147808829414345923316083210206383297598`, wrapped path gives
`103144121322099306484875023187381681344`. Notably `qp.rs` itself uses
`checked_mul + expect` everywhere; this is the one site that silences the
check and returns a wrong ring element.

### S-12 · MINOR · `scalar/integrality.rs` — `Zp::to_fraction` casts residue `as i128`: wraps for residues ≥ 2¹²⁷

Lines 116–121. Moduli above 2¹²⁷ are constructible (5⁵⁵, 11³⁷); a residue
r ≥ 2¹²⁷ casts to r − 2^128 and `from_i128` embeds the wrong integer
(verified numerically for `Zp<5,55>`, r = 2¹²⁷+1; 2^128 mod 5⁵⁵ ≠ 0). The
round-trip test only exercises `Zp<3,3>`.

### S-13 · MINOR · `scalar/mod.rs` — `mod_inverse_u128` returns `None` for every modulus > i128::MAX, even for genuine units

Lines 119–134: `i128::try_from(modulus).ok()?`. So `Zp::<5,55>::inv(2)` =
`None` although 2 is a unit — a silent false claim of non-invertibility
against the documented "None if not invertible (zero)" contract, also
affecting runtime-prime `LocalQp`. Verified against an exact copy of the
function (`mod_inverse_u128(2, 5³) = Some(63)` ✓; `(2, 5⁵⁵) = None` ✗).

### S-14 · MINOR · `scalar/analytic.rs` — raw u128 products in `mod_pow`: misclassified quadratic residues for P > 2^64

Lines 111–122 (`mod_pow`), feeding `fp_is_square`/`fp_sqrt` and the Hensel
seeds for Zp/Qp/Qq. For p = 2⁶⁴+13 and a = p−2, the wrapped Euler criterion
returns 1 ("square") while the true Legendre symbol is −1 (verified with
sympy). Debug panics instead. Inconsistent with `fp.rs`, which deliberately
implements overflow-safe double-and-add `mul_mod` so the *field* works for
all P < 2^128 — the type is fine, its `ExactRoots` impl silently is not.

### S-15 · DOC · `scalar/exact/mod.rs` — `(F_{p^n}, W_n)` listed as a (field, ring-of-integers) pair

Module doc, lines 8–10. The fraction field of W_N(F_q) is **Q_q**; F_{p^n} is
its *residue* field. The code itself implements the correct pairing
(`impl HasFractionField for WittVec { type Frac = Qq }`), as do
`scalar/mod.rs`'s table and `scalar/AGENTS.md` — the doc inverts the
residue-vs-fraction direction the whole pillar is organized around.

---

## Findings — `src/clifford/`

### C-1 · MAJOR · `clifford/engine/algebra.rs` — `reverse()` is the per-blade exterior reversal, not Clifford reversion, on any non-orthogonal metric

`CliffordAlgebra::reverse`, lines 142–157. The diagonal sign
`(−1)^{k(k−1)/2}` per wedge-blade equals the reversion anti-automorphism (the
unique anti-automorphism fixing vectors) **only when `b` and `a` vanish**: in
the engine's own basis identification, true reversion is non-diagonal —
`τ(e₀e₁) = e₁e₀ = b₀₁ − e₀e₁`, but `reverse` gives `−e₀e₁`. Verified against
a line-by-line Python transcription of `product.rs`: with q = (1,1), b₀₁ = 1
over ℚ, `reverse(e₀·e₁) = −e₀e₁` while `reverse(e₁)·reverse(e₀) = 1 − e₀e₁` —
not even an anti-homomorphism. Same failure for a-only metrics. **In char 2,
`reverse` is the identity map**, which cannot be the reversion of any
noncommutative algebra — and every nonsingular char-2 metric (the project's
flagship nim-Clifford case) has nonzero off-diagonal `b`. Downstream damage
verified numerically: `norm2(1+e₀e₁)` returns 2 where the true spinor norm of
the genuine rotor `R = e₀(e₀+e₁)` is 3 (see C-4).

### C-2 · MAJOR · `clifford/blade.rs` — `blade_subspace` silently returns vectors *not* in {x : x∧A = 0} over non-field scalars

Lines 125–171; root cause is L-1. Counterexample over `Integer`, found by
exhaustive scan: n = 3, `A = −2·e₀e₁ − 2·e₀e₂ = −2·(e₀ ∧ (e₁+e₂))`, a genuine
2-blade. The wedge-map matrix is the single row `[0, 2, −2]`; no entry is a
unit, so all three columns are declared free and the function returns
`Some([e₀, e₁, e₂])` — three "basis" vectors for a 2-dimensional subspace,
including `e₁` with `e₁ ∧ A = 2·e₀₁₂ ≠ 0`. Both the dimension and the
membership claims are silently wrong (no `None` exists on this path),
violating the unconditional docstring contract.

### C-3 · MINOR · `clifford/engine/inverse.rs` — `multivector_inverse` shift overflow for dim ≥ 64

Line 13, `1usize << self.dim` with dim allowed up to 128. Debug panics at the
shift; release masks the shift amount and then indexes out of bounds (panic)
for any non-scalar input. No wrong result escapes — a robustness gap in a
documented "for any element" API (which is memory-bound far below dim 64
anyway).

### C-4 · MAJOR · `clifford/versor.rs` — grade-sign `reverse` breaks the entire versor layer on b ≠ 0 metrics

`norm2` (45–49), `versor_inverse` (67–80), `sandwich`/`twisted_sandwich`
(85–104), `spinor_norm`, `classify_versor`. Companion to C-1, verified
independently against a Python mirror of `product.rs`: on the nondegenerate
metric q = [1,1], b₀₁ = 1 over ℚ (Gram det 3/4), the genuine rotor
`R = e₀·(e₀+e₁) = 1 + e₀e₁` has true spinor norm `R·τ(R) = q(e₀)q(e₀+e₁) = 3`,
but `norm2(R)` returns **2** (no gate protects it — silently wrong
invariant), and `versor_inverse`/`sandwich`/`spinor_norm`/`classify_versor`
all return `None` on R because `R·reverse(R) = 2 − e₀e₁` is not scalar — a
false refusal of a genuine versor. Together with C-1 this means the
versor/Pin layer is only trustworthy on orthogonal (b = a = ∅) metrics.

*Fix shape for C-1/C-4:* implement reversion honestly — e.g. reduce the
reversed generator word through the engine (`reduce_word` on the reversed
word), or restrict `reverse`-consuming APIs to orthogonal metrics with a
loud boundary.

### C-5 · MINOR · `clifford/outermorphism.rs` — `trace`/`char_poly`/`exterior_power_trace` wrong at dim = 128

`grade_k_masks` (122–142) evaluates `c >> limit_bits` with `limit_bits = 128`
on u128: debug panics ("attempt to shift right with overflow"); in release
the shift is masked and the Gosper loop never runs, so
`exterior_power_trace = 0` for every k ≥ 1 — `trace(identity) = 0` instead of
128, `char_poly = [1, 0, …, 0]` for *any* map, silently. `hopf.rs`'s `pairs`
helper shows the correct guard (`if dim >= MAX_BASIS_DIM { u128::MAX }`);
this site lacks it. (`determinant()` is unaffected via `pseudoscalar()`.)

### C-6 · DOC · `clifford/versor.rs` — Cayley-transform doc claims the result is a versor/Spin element; false from dim 6

Doc at 218–229. Even and unit-spinor-norm are correct, but for the non-simple
bivector `B = e₀₁ + e₂₃ + e₄₅` in Cl(6,0), `R = (1−B)(1+B)⁻¹` satisfies
`R·τ(R) = 1` yet `R e₀ τ(R)` leaks grade 5 (verified with exact Fraction
arithmetic) — R is not in the Lipschitz group. Consistent with the classical
fact that {even, RR̃ = 1} = Spin only for dim ≤ 5. The code (standard Cayley
formula) is fine; the unit test only exercises dim 3.

### C-7 · DOC · `clifford/spinor_norm.rs` — module doc asserts 1 → {±1} → Pin(Q) → O(Q) → 1, which fails over general fields

Lines 6–12. The correct sequence ends `→ O(Q) → F*/F*²` — the last map being
exactly the spinor norm this module computes. The module's own test
(`spinor_norm(e₀+e₁) = 2`, a nonsquare in ℚ) exhibits an isometry not hit by
norm-(±1) versors, contradicting the claimed surjectivity. The invariant is
interesting *because* the stated sequence is not exact at O(Q).

### C-8 · DOC · `clifford/cga.rs` — `exp_nilpotent` doc equates cap-exhaustion with "not nilpotent"

Lines 182–204. The cap `2·dim+2` is below the maximal nilpotency index from
dim 10 up (Cl(5,5) ≅ M₃₂(ℚ) contains nilpotents of index up to 32 > 22), so a
genuinely nilpotent element can get `None` and be misclassified by a caller
trusting the parenthetical. For the intended PGA-motor inputs the cap is
ample. The `None` is a refusal, not a wrong value; only the documented
equivalence is wrong.

---

## Findings — `src/forms/`

### F-1 · MAJOR · `forms/hermitian.rs` — pivot manufacture uses λ instead of conj(λ): panics on valid nondegenerate forms

`ensure_pivot` (79–87), panic surfaces in `diagonalize` at line 189. The
congruence E*HE with E = I + λ·E_{j,k} gives new
`H[k][k] = λ² + conj(λ)² = 2·Re(λ²)` — **not** the `2|λ|²` the in-code
comment claims; the comment describes the correct algorithm (λ = conj(h[k][j]))
and the code implements the wrong one. `2·Re(λ²)` vanishes whenever
`|Re λ| = |Im λ| ≠ 0`. Counterexample, verified by running the crate:
`H = [[0, 1+i], [1−i, 0]]` (conjugate-symmetric, det −2, signature (1,1))
passes `from_gram` and then `.diagonalize()`/`.signature()` **panic** with
"nonzero real pivot inverts in a field". With the conjugate fix the same
input diagonalizes to `[4, −1/2]`, signature (1,1) ✓.

### F-2 · MINOR · `forms/char0.rs` — `classify_real`/`classify_complex` overflow at dimension 128

`p2` (171–173): `matrix_dim = 1usize << matrix_exp` overflows for
`matrix_exp = 64` (reachable: nondegenerate dim-128 form with q−p ≡ 0, 6 mod 8,
or `classify_complex(128)`; the engine accepts dim 128 =
MAX_BASIS_DIM). Verified: debug panics; release silently claims
Cl(128,0) ≅ ℝ where the correct answer is M_{2^64}(ℝ). `usize` cannot hold
2⁶⁴; the type needs an `Option` guard or a log₂ representation.

### F-3 · MINOR · `forms/symplectic.rs` — `classify` returns wrong rank/radical over non-field scalars, silently

Lines 107–114; root cause L-1. Over ℤ the alternating Gram `[[0,2],[−2,0]]`
passes `from_gram`, yet `classify` returns `{rank: 0, radical_dim: 2}`
(verified by running the crate) — the true kernel is 0 (det = 4). The module
doc states the theorem "over any field" but the API neither restricts to
fields nor returns `None` on a non-unit pivot.

### F-4 · DOC · `forms/char0.rs` — "Cl(p,q,r) ≅ Cl(p,q) ⊗ Λ(F^r)" is false for the ungraded tensor product

Module doc 32–33; `CliffordType` display ("C ⊗ Λ(R²)"); mirrored in
`OddCharType::display`. The correct statement needs the **graded** product
⊗̂. Counterexample verified with the crate's own engine: Cl(0,1,1) is
4-dimensional noncommutative (`e₀e₁ = −e₁e₀ ≠ e₁e₀` over `Rational`), while
ℂ ⊗ Λ(ℝ¹) is commutative. The classification *data* (p,q,r) and all
comparisons are unaffected.

### F-5 · MAJOR · `forms/witt/class.rs` — char-2 Witt/Brauer–Wall class addition is not a homomorphism across nim-subfields

`WittClass::add` (60–64), `WittClassG::try_add` char-2 arm (180–182), also
`BrauerWallClass::try_add`. `arf_nimber` evaluates each metric's Arf over its
own minimal field of definition, and **Arf is not stable under field
extension**: q = [1,1], b₀₁ = 1 has Arf 1 over F₂ but Arf 0 over F₄
(Tr_{F₄/F₂}(1) = 0). The char-2 class variants store only the bare arf bit
and XOR it unconditionally — unlike the OddChar variant, which stores
`field_order` and rejects cross-field addition. Verified in the shipped
crate: A = Metric([\*1,\*1], b₀₁=\*1) → arf 1 (over F₂); B = Metric([\*2,\*2],
b₀₁=\*1) → arf 1 (over F₄); XOR predicts 0 for A ⊥ B, but
`arf_invariant(A.direct_sum(B))` = 1 — confirmed by brute-force zero count
over F₄: the rank-4 sum has 52 = 4³ − 3·4 zeros (the Arf-1 count; Arf 0 would
give 76). Fix shape: carry the field of definition in the char-2 classes
exactly as OddChar already does, and evaluate both operands over the
compositum (or reject) in `try_add`.

### F-6 · DOC · `forms/char2/arf.rs` — header misstates the classification theorem (drops "graded" and "same rank")

Lines 1–5: "two such algebras are isomorphic iff their F₂ forms share an Arf
invariant." As *ungraded* algebras, the Clifford algebras of **both** rank-2
nonsingular F₂ forms are M₂(F₂) (Wedderburn: no finite noncommutative
division algebra, Br(F_q) = 0), so different Arf, isomorphic algebras. Also
false across ranks with equal Arf. The true statement — ℤ/2-**graded**
isomorphism at fixed rank — is exactly what the rest of the file and
`witt/brauer_wall.rs` correctly compute.

### F-7 · CRITICAL · `forms/integral/genus.rs` — 2-adic Jordan splitting misclassifies odd blocks as type II on a valuation tie

`jordan_blocks`/`min_val_entry`, ~157–265. `min_val_entry` scans the upper
triangle row-major and keeps the *first* entry of minimal 2-valuation; the
code peels a 2-dim "type II" block whenever that entry is off-diagonal. The
correct dispatch: peel a 1-dim (odd) block if **any diagonal** entry attains
the minimal valuation; the 2-dim even block is only correct when all diagonal
valuations are strictly larger. The scan order lets `a[0][1]` win a tie
against `a[1][1]`. Counterexample, verified against the compiled crate:
`G = [[2,1],[1,1]]` is ℤ-isometric to ℤ² (basis (1,1),(0,1); det 1,
represents 1), but `Genus::of(ℤ²)` gives (scale 0, dim 2, **I**, oddity 2)
while `Genus::of(G)` gives (scale 0, dim 2, **II**, oddity 0), and
`are_in_same_genus(ℤ², G)` returns **false** for isometric lattices. The
randomized isometry-invariance test misses this because its shears are all
strictly upper-triangular, which never disturbs the scan order.

### F-8 · MAJOR · `forms/integral/genus.rs` — 2-adic trains break across a missing scale, losing Conway–Sloane sign-walking identifications

`two_adic_trains`, ~430–448. The code requires strictly consecutive scales,
but in Conway–Sloane (SPLAG ch. 15 §7.5; also Sage's
`canonical_2_adic_trains`) absent scales count as dimension-0 type-II
constituents, so a train *continues* across one empty scale when both
flanking constituents are type I. The code therefore never sign-walks
between scales 0 and 2, splitting one genus into two. Verified two ways:
brute force over ℤ/64 exhibits `U ∈ GL₂` with `Uᵀ·diag(1,20)·U = diag(5,4)`
(the same brute-forcer exactly reproduces the code's equivalence classes on
all 16 adjacent-scale forms, validating both oracle and walking rule), and
the two forms agree in det, signature, and p = 5 symbol — same genus, yet
`are_in_same_genus` returns false.

### F-9 · DOC · `forms/integral/root_lattices.rs` — |Aut(A_n)| = 2(n+1)! claimed unconditionally; false for n = 1

`a_n` doc ~57–59. A₁ = ⟨2⟩ has Aut = {±1}, order 2, not 4 (for n ≥ 2 the
formula is right since −1 ∉ W(A_n); for n = 1, −1 *is* the nontrivial Weyl
element). The implementation is correct — `a_root_automorphism_order`
special-cases n = 1 → 2 — so prose-only.

### F-10 · CRITICAL · `forms/springer/char2.rs` — local engine treats P-adic polynomial-digit expansions as κ((π)) series, dropping carries at every finite place of degree ≥ 2

`laurent_finite` (~198–239), `asnf` (~316–346), `block_contribution`
(~474–527), `local_is_square`. At a place P of degree d ≥ 2, elements are
expanded as Σ g_k·P^k with polynomial digits of degree < d, then the digits
are combined *digit-wise* (`asnf` folds an even pole via `m[n/2] += sqrt(c_n)
mod P`; `block_contribution` multiplies digit sequences with `mul_mod P`;
`local_is_square` splits by digit parity). That is only valid over a genuine
coefficient field: polynomial digit representatives are not multiplicatively
closed for d ≥ 2 — `s² = digit + e·P`, and the dropped carry `e·P^{n+1}`
changes exactly the Laurent coefficients the Artin–Schreier normal form
reads. The sibling engine (`local_global/function_field_char2.rs`) does this
correctly via honest κ[[u]] arithmetic (`hensel_series`/`ps_eval_poly`).
Empirically confirmed against the public API with ground-truth
counterexamples starting at `c = ℘(t/P) = (t³+t)/(t²+t+1)²` over F₂(t).
*Fix shape:* route springer/char2's expansion and folding through the
function_field_char2 power-series helpers.

### F-11 · MINOR · `forms/local_global/function_field.rs` — u128 overflow computing |κ| = q^(deg π) at high-degree places

`kappa_order` (74–85), same pattern in `kappa_inv`.
`S::field_order().pow(deg)` overflows for deg ≥ 128 over F₂ (deg ≥ 56 over
F₅, …): debug panics, release wraps — after which the Euler-criterion
exponent (|κ|−1)/2, the Fermat-inverse exponent |κ|−2, and the sqrt exponent
|κ|/2 are all wrong, producing wrong Hilbert symbols/residues/AS classes with
no error. Extreme but legal inputs; everything below the threshold is exact.

### F-12 · DOC · `forms/witt/class.rs` — W_q "can be richer" over algebraically-closed On₂; it is trivial there

Module header 10–12. Over an algebraically closed char-2 field every binary
form ax² + xy + by² is isotropic (y² + y = ab solvable), so W_q(On₂) = 0 —
strictly *smaller* than the ℤ/2 of the finite subfields. "Richer" points the
wrong way for the named example (genuinely richer fields exist, e.g. F₂(t),
which the file's siblings handle).

---

## Findings — `src/games/`

### G-1 · MAJOR · `games/loopy.rs` — `over + under = 0` is false; the sum is a draw-class value

`LoopyValue::add` line 171, doc line 157, pinned by test `the_closed_sums`
line 598. In over+under each player owns a loop move while every exiting move
loses, so under the module's own normal-play-with-draws convention the
position is a **Draw** with either player to move — and a game equal to 0
must be a P-position. Not `dud` either: brute-force retrograde analysis shows
over+under+1 is a Left-win regardless of who starts, while dud+1 stays drawn
(dud is absorbing). The correct return per the function's own contract
(`None` when the sum leaves the catalogue) is `None`.

### G-2 · MINOR · `games/loopy.rs` — `partial_cmp` declares ⋆ confused with over/under; actually under < ⋆ < over

Lines 179–212. By the standard stopper-comparison criterion (Siegel: G ≥ H
iff Left survives G−H playing second; over and ⋆ are stoppers): in over+⋆
with Right first, Left wins (Right's loop is answered by ⋆→0 leaving over;
Right's ⋆→0 leaves over); conversely ⋆ ≥ over fails. Hence over > ⋆ > under
strictly. Only ⋆ ∥ 0 and the dud-confusions are correct. Pinned by the
passing test.

### G-3 · MINOR · `games/loopy.rs` — over+over, under+under, ⋆+over, ⋆+under do *not* leave the catalogue

`add` comment + `_ => None` arm. By the stopper criterion, over+over = over
and ⋆+over = over (mirrors for under): Left survives both
over+over+under and over+under+under playing second; for ⋆+over both
difference games reduce to the same multiset {over, ⋆, under}. Independently,
outcomes of over+over+X vs over+X agree across 14 test contexts (1, −1, up,
down, star, on, off, …). Returning `None` on sums equal to a catalogue value
contradicts the function's own spec.

*Fix shape for G-1/2/3 (per the reviewer):* `(Over,Under) → None`; add
`(Over,Over) → Over`, `(Under,Under) → Under`, `(Star,Over) → Over`,
`(Star,Under) → Under`; `partial_cmp`: ⋆ vs over = Less, ⋆ vs under =
Greater — and re-pin the tests to the corrected table.

### G-4 · MINOR · `games/game_exterior.rs` — relation certificates from `with_relations` hardcode `independent: true`

Line 131 passes a trailing `true` to `relation_search_certificate`, marking
every explicitly supplied relation independent without checking. Supplying
relations [2,0] and [4,0] over gens [⋆, ↑] yields two rows both flagged
independent — a false certificate. The search path enforces independence
properly; the algebra quotient itself is unaffected
(`reduce_integer_vector` handles dependent rows).

### G-5 · MINOR · `games/misere.rs` — misère-Nim predicate miscounts when zero heaps are present

`misere_nim_p_predicted` (113–121) uses `heaps.len() % 2` in the all-heaps-≤1
branch; Bouton's misère theorem counts **nonzero** heaps. `[1, 0]`: returns N,
truly P (verified by brute-force misère retrograde analysis). All internal
callers canonicalize first (zeros dropped) so tests pass, but the function is
`pub`, takes any `&[u128]`, and documents no nonzero-heap precondition.
Mirrored at the Python surface (P-5).

### G-6 · DOC · `games/hackenbush.rs` — comment value wrong: the stalk Red,Blue,Red,Blue is −5/8, not −3/8

Test comment line 202 (`// −+−+ = −3/8`). The sign expansion −+−+ walks
−1, −1/2, −3/4, **−5/8**; an independent dyadic-simplicity evaluation of the
stalk agrees. The assertion itself pins game value == `from_sign_expansion`
(both correct); only the annotation is wrong — in a file whose stated purpose
is pinning Berlekamp's rule.

### G-7 · DOC · `games/AGENTS.md` — stale Turning-Corners boundary "None ≥ ω^ω" contradicts the implementation

Line 33. The implementation and its passing test compute
⋆ω^ω ⊗ ⋆ω^ω = ⋆ω^(ω·2) (Some); `None` only at ≥ ω^(ω^ω) or Kummer carries
needing a prime > 47 (consistent with α₄₇ in `tower.rs`). The test comment at
`nimber_game.rs:195–205` even says "(was staged under the old ω^ω boundary)".
Same stale sentence copied into the Python docstring (P-7).

---

## Findings — `src/linalg/` and the test suites

### L-1 · MAJOR · `linalg/field.rs` — `unit_pivot_nullspace` silently returns non-null vectors over ring scalars; module doc falsely claims a None boundary

Lines 84–126; module doc 1–6. The doc says these kernels "return None when a
required nonunit pivot appears" — true for `solve`/`inverse_matrix`, **false**
for `unit_pivot_nullspace`, which has no `None` path: a column whose remaining
entries are nonzero non-units is treated as *free*, so the returned "basis"
contains x with Mx ≠ 0. Verified end-to-end through the public API
(`clifford::blade_subspace`, see C-2: row [0, −2, 2] over ℤ → all three
columns "free"). Over fields — the mainline Nimber/Rational/Surreal uses —
the routine is correct. This is the root cause of C-2 and F-3; one honest
`Option`/`None` path fixes all three.

### L-2 · MINOR · `linalg/integer.rs` — three unchecked i128 sites deviate from the module's checked-arithmetic discipline

(1) `reduce_integer_vector`'s `v[i] -= q * row[i]` — the only unchecked row
op in the file; (2) `ext_gcd`'s sign fix `(−r0, −s0, −t0)` overflows for
input `(i128::MIN, 0)`, returning a negative "gcd" in violation of the
documented g ≥ 0 contract; (3) `smith_normal_form`'s final `.abs()` on a
surviving `i128::MIN` diagonal. Release wraps silently; debug panics without
the module's documented overflow messages. Adversarial magnitudes only.

### L-3 · DOC · `tests/clifford_axioms.rs` — header claims the fuzz covers the general bilinear `a` term; no strategy ever sets `a`

Header lines 5–8 vs the strategies: both property tests build metrics
exclusively via `Metric::new(q, b)`/`Metric::diagonal(q)`, which hard-code
`a = ∅`. The only `a`-term coverage in the repo is one fixed unit test.
Either drop the claim or add a `Metric::general` strategy.

---

## Findings — `src/py/` (bindings)

### P-1 · MINOR · `py/scalars.rs` — `relative_trace_over`/`relative_norm_over` return values that are not traces/norms

`validate_relative_degrees` (29–37) checks only e > 0 and e | m; the core
default methods compute Σ/Π of x^(p^(e·i)) with no check that x ∈ F_{p^m} or
m | ext_degree. Verified: `Nimber(4)` (degree 4) `.relative_trace_over(2,1)`
returns ⋆2 — not even an element of the claimed codomain F₂;
`Nimber(2).relative_trace_over(256,1)` accepts m = 256 > 128 and returns ⋆0
(each conjugate counted twice). Partial validation invites the inference that
the rest is validated; instead out-of-subfield inputs produce silent garbage
labelled as a trace. Fix: require m | ext_degree() and x ∈ F_{p^m}.

### P-2 · MINOR · `py/scalars.rs` — `__eq__`/`__hash__` contract violated for Nimber, Fp\*, Zp\*, WittVec\* against Python ints

These families parse ints through the reducing constructor in `__eq__` but
hash the canonical internal value: `Fp13(-1) == -1` is True yet
`{Fp13(-1): 'x'}[-1]` raises KeyError; `Nimber(2**64+5)` hashes the u128
truncated to usize. The reduce-on-parse equality is also non-transitive
(`14 == Fp13(1)`, `Fp13(1) == 27`, `14 != 27`), breaking dict/set semantics a
second way. Classes *without* explicit `__hash__` (Surreal, Rational, Qp, …)
are simply unhashable — the safe behavior — so the bug is confined to these
four families. Fix: restrict `__eq__` to canonical ints, or hash canonically.

### P-3 · MINOR · `py/scalars.rs` — three-argument `pow(x, n, m)` silently ignores the modulus for Surreal, Surcomplex, and every multivector class

`pow(Surreal(2), 3, 5)` returns 8 (verified; same for Surcomplex and
multivectors) instead of raising — inconsistent with the same file's own
convention (`PyNimber.__pow__` and the finite-field classes raise ValueError
on a supplied modulus). Returning the unreduced power is a silently wrong
result for an accepted input.

### P-4 · DOC · `py/scalars.rs` — `Ordinal.nim_mul` docstring claims None for any infinite operand

Lines 5443–5449 vs core `Ordinal::nim_mul` (nim.rs:72–83), which returns Some
throughout the verified Kummer window — `Ordinal.omega().nim_mul(omega)`
returns ω² (verified empirically; Conway: ω³ = 2). The docstring contradicts
the adjacent `__mul__` docstring and the AGENTS/ROADMAP claims.

### P-5 · MINOR · `py/games.rs` — `misere_nim_p_predicted` wrong on heap lists containing zeros

Lines 800–808; the binding surface of G-5, independently verified by
brute-force misère analysis: `[1,0]` truly P, predicted N; `[0]` truly N,
predicted P. Self-contradictory against
`try_misere_is_n([1,0], nim_moves)` in the same module. On canonical
zero-free inputs the predictor is correct.

### P-6 · MINOR · `py/forms.rs` — `arf_f2` accepts n > 128, which the core cannot represent

Lines 163–185: the guard `if n >= u128::BITS { domain_mask = u128::MAX }`
lets any n > 128 through to core `arf_f2`, which builds basis vectors
`1u128 << i` — shift overflow at i = 128 (debug: PanicException; release: the
129th basis vector aliases e₀, producing silently wrong Arf data — verified
construction: n = 129, q anisotropic only on e₁₂₈ reports
`radical_anisotropic = false`). Accept exactly n == 128, reject n > 128
(mirroring the explicit k ≤ 20 cap in `fit_f2_quadratic`).

### P-7 · DOC · `py/games.rs` — `turning_corners` docstring repeats the stale "None at/above ω^ω" boundary

Line 1591; same stale sentence as G-7 (apparently copied from games/AGENTS.md
rather than the current core doc). A Python user is told ⋆ω^ω ⊗ ⋆ω fails when
it succeeds.

---

## Refuted in verification

Raw findings that died under the two-skeptic protocol — kept here because the
refutation reasoning is itself useful.

1. **"q=[1,1], b₀₁=1 over Nimber is called 'the anisotropic plane' but is
   hyperbolic over F_{2^128}"** (`clifford/spinor.rs`, doc). The ambient-field
   arithmetic is right (Tr_{F_{2^128}/F₂}(1) = 0, isotropic vector (ω, 1)),
   but the repo's documented convention — `forms/char2/arf.rs` classifies nim
   metrics **over the field of definition** of their entries, reduced via
   trace — makes the label correct as used. Refuted on documented intent.
2. **u128 overflow in `divided_power.rs`'s binomial helper around total
   degree 126+.** The math skeptic *confirmed* the overflow (central-binomial
   intermediates first exceed u128 at n = 126), but the intent skeptic found
   it documented/bounded; under the both-must-confirm rule it dropped.
   Honest residual: if `divided_power` is ever exercised at degree ≥ 126,
   re-examine. Borderline case worth knowing about.
3. **"ArfResult.arf is not an isometry invariant for defective singular
   forms."** The counterexample is mathematically real (verified), but it is
   the crate's documented, deliberately-handled boundary —
   `equivalence.rs:299–315` has a test using the *exact same metric pair*,
   asserting `a1.arf != a2.arf` and routing isometry through the
   radical-aware path. Refuted on documented intent.
4. **"Float-bounded Fincke–Pohst fallback has no rigorous no-drop
   guarantee"** (`integral/lattice.rs` vs an AGENTS claim). The math verifier
   died (API error) before reaching a verdict, so this defaulted to refuted
   rather than being judged. Honest status: **unresolved**, not refuted —
   worth a five-minute look at the float-bound margin if exact enumeration
   matters to a downstream claim.

---

## Sub-threshold observations

Reviewer notes that didn't rise to findings but are worth recording:

- **Root `AGENTS.md` says the Kummer window is "u ≤ 43"; the implementation
  and `scalar/AGENTS.md` ship u ≤ 47** (the 47 row probe-verified rather than
  source-verified; OPEN.md documents the provenance split precisely). Stale
  summary line.
- **Property-test volume is smoke-sized as configured** (FAST_CASES = 2,
  HEAVY_CASES = 1 per suite per run; documented, with the
  `PLEROMA_PROPTEST_CASES` override) — near-zero random coverage per CI run.
- **The rational Clifford fuzz uses diagonal metrics only**, so char-0
  nonzero-`b` associativity is exercised only by fixed unit tests (compare
  L-3 for the `a` term).
- `versor.rs::even_subalgebra`'s claimed isomorphism Cl(Q)⁰ ≅ Cl(Q′) needs
  `q_p` invertible; over ring backends with q_p = ±2 the doc overstates.
- `mass_formula.rs` has a stale test comment about Golay generator weights
  (all 12 shipped rows have weight 8; the assertion itself is fine).
- `lattice.rs::level()`'s doc is precise only for even lattices.
- `py/engine.rs::backend_algebra::new` lets `b = {(0,1): x, (1,0): y}`
  silently last-wins after key normalization.
- The two char-2 local engines disagree by construction
  (`local_global/function_field_char2.rs` is the honest one); see F-10's fix
  direction.
- `BW(F_{2^m}) ≅ ℤ/2` could not be source-checked offline; the reviewer's
  plausibility argument (char-2 Koszul sign vanishes; Br(F_q) = 0) is
  recorded in the transcripts.

---

## Coverage map

| Area | Files | Status |
|---|---|---|
| nimber-ordinal | `scalar/finite_field/nimber/*`, `scalar/big/ordinal/*`, `big/cnf.rs` | read fully; nim table 0..16 recomputed independently |
| surreal-omnific | `scalar/big/surreal/*`, `big/omnific.rs`, `functor/surcomplex.rs` | read fully + `exact/rational.rs` for the coefficient world |
| padic-local | `scalar/small/*`, `wittvec.rs`, `functor/*`, `global/*`, `valued.rs` | read fully; "in good mathematical shape" overall |
| scalar-core-ff | `scalar/{mod,exactness,integrality,extension,residue,poly,tropical,analytic}.rs`, `exact/*`, `fp/fpn` | read fully; Fpn irreducibles pinned by the exhaustive field-axiom tests |
| clifford-engine | `engine*`, `blade.rs`, `mod.rs` | read fully; sign logic verified against brute-force transposition counts |
| clifford-structures | `spinor*`, `versor`, `cga`, `outermorphism`, `frobenius`, `hopf`, `divided_power` | read fully |
| forms-core-char0 | `char0`, `classify`, `diagonalize`, `equivalence`, `field_invariants`, `hermitian`, `symplectic`, `trace_form`, `poly_factor`, `oddchar/*` | read fully; documented-intentional Nones respected |
| forms-char2-witt | `char2/*`, `witt/*`, `quadric_fit.rs` | read fully; zero counts brute-forced over F₂/F₄ |
| forms-local-global | `local_global/*`, `springer/*` | read fully; Hilbert symbols checked line-by-line against Serre |
| forms-integral | `integral/*` | read fully; Bernoulli/mass constants, E₄/E₆/Δ expansions, Leech/D₁₆⁺ data all recomputed exactly |
| games | `games/*` (15 files) | read fully; partizan recursions checked against standard definitions |
| linalg-axioms | `linalg/*`, `lib.rs`, `tests/*` | read fully; SNF ported to Python and fuzzed against the port |
| py-scalars-engine | `py/{mod,scalars,engine}.rs` | read fully (7.5k lines); installed `.venv` extension verified to match this checkout |
| py-forms-games | `py/{forms,games,catalog}.rs` | read fully (9.2k lines); catalog is a pure name/type manifest, spot-verified clean |
| **docs-claims** | README, ROADMAP, OPEN, OPEN-3, TABLES, `writeups/goldarf.tex` | **NOT AUDITED** — reviewer died (API socket error); left as a gap by decision |
| **examples-experiments** | `demo.py`, `examples/*`, `experiments/*` | **NOT AUDITED** — reviewer died (API socket error); left as a gap by decision |

The two gaps mean: no checkable-claim sweep was performed over the prose
documents or the probe scripts. Findings above that touch them (G-7 on
`games/AGENTS.md`, the stale u ≤ 43 note) came from adjacent reviewers,
not a systematic pass.

---

## Suggested triage order

1. **F-7 / F-8** (genus symbols) — `are_in_same_genus` returning false for
   isometric lattices poisons anything downstream that trusts genus
   equality; both fixes are localized (tie-break dispatch; train continuation
   across one empty scale per SPLAG §7.5).
2. **S-3** (surreal analytic panics) + **S-1/S-2** — one design fix.
3. **C-1 / C-4** (reversion) — decide: honest reversion via the engine, or a
   loud orthogonal-metrics-only boundary on the versor layer. The current
   state silently mis-handles exactly the char-2 metrics the project is
   about.
4. **F-10** (springer/char2 carries) — route through the
   `function_field_char2` power-series helpers.
5. **L-1** (+ C-2, F-3 for free) — add the `None` path.
6. **F-5** (char-2 Witt class field tag) and **G-1/2/3** (loopy table) —
   wrong group elements / wrong CGT identities with tests pinning them.
7. **F-1** (hermitian conj) — one-character-class fix, panic on valid input.
8. The u128/i128 overflow family — consider `overflow-checks = true` in the
   release profile as a blanket mitigation, then pick off the sites that
   should return `None` instead of panicking.
9. The 14 doc items — mostly one-line edits; F-6, C-7, S-15 are the ones a
   reader would actually be misled by.

---

*Generated from the audit run `wf_85c8b24b-348` (2026-06-09). Full per-finding
descriptions, verifier transcripts, and reproduction scripts are preserved in
the session transcript directory; raw structured output in
`/tmp/pleroma-audit-result.json`.*
