# Cross-pillar bridges — DONE (the go-forward ledger)

The running ledger of cross-pillar work **completed from here on**.

The cross-pillar bridge-building era (bridges **A–O** plus **K** — lattice/Clifford/
Brauer–Wall, the char-2 Arf classifier, Frobenius outermorphisms, the transfinite
Clifford engine, theta/modular forms, Construction-A codes, the Weil representation, the
rational and full-`ℚ/ℤ` Brauer invariants, Newton polygons, the Brown invariant, the
unification pass, lexicodes) closed with every non-deferred bridge shipped. Its full
record — mathematical contracts, surfaces, oracles, boundaries, and the formal proofs —
lives in [`roadmap/CODA.md`](CODA.md); the terse working-notes summary is in the
`AGENTS.md` files (root + per-pillar).

What remains unbuilt is tracked in [`roadmap/TODO.md`](TODO.md) — the game-valued
ledger of buildable items plus the deferred stars `*1`/`*2`/`*4`; genuine open
problems stay in [`OPEN.md`](../OPEN.md), loopy-valued: `tis`/`tisn`, `on`/`off`,
`over`/`under` (the old numerals §1–§4 survive as aliases).

## How to use this ledger

When a new piece of cross-pillar work lands, add a short entry here:

```
## completed items

### <date>: `<name>`
**Summary:** <one-line what-it-connects>
**Pillars:** … ↔ …    **Claim level:** standard math / implemented-and-tested / …
- surface: the functions/types that shipped
- oracles: the tests that pin it
- boundaries: the honest non-claims
```

Promote anything proof-heavy or too long for working notes into `roadmap/CODA.md`, and
fold the one-line structural fact into the relevant `AGENTS.md`.

## completed items

### 2026-06-12: `exception-column-m4`
**Summary:** the Lenstra excess on the `f(p) = 2*3^k` column is `m_p = 4` exactly at
every prime current factor tables reach
**Pillars:** scalar (ordinal tower) ↔ excess analysis (`OPEN.md` §3 / `on`)
**Claim level:** proved norm identity + certified computation (analysis-level
A380496-type rows; no new shippable `alpha_u` carries — the Rust boundary at
`alpha_53` is untouched)
- surface: `experiments/exception_column_m4.py` (committed, stdlib-only, ~2 min).
  The key fact is a *corrected* compositum norm: `sigma(4) = 5` (the
  Artin–Schreier conjugate over `F_4`; `writeups/excess.tex` §4.3 had
  `(kappa+4)(kappa+6)`, a Frobenius slip — erratum recorded in the note), so
  `Norm(kappa+4) = (kappa+4)(kappa+5) = kappa^2 + kappa + omega`, collapsing the
  `m = 4` root test into the same trinomial field as the `C_k` certification:
  `m_p = 4 <=> p | ord(M_k)` with `M_k = Nbar/N`, `N = zeta^2 + zeta + zeta^h`.
- results: `m_p = 4` **universally** for `k = 2..6` (fully factored levels; 14
  rows, 11 new — among them `87211`, `135433`, `272010961`, `139483`, two
  ~20-digit primes, four larger, and a P78); **consistent** at `k = 7, 8` (all
  11 known primes). New conjecture `D_k` (prime-to-3 part of `ord(M_k)` full) =
  the column analogue of `C_k`; proved twisted-norm lemma
  `Norm(N_k) = eta^2 + omega^2*eta + 1 != N_{k-1}` (no `gamma`-style
  propagation); observed 3-part parity `v_3(ord M_k) = k+1` (odd `k`) / `k-1`
  (even `k`), unexplained.
- oracles: in-script audits (product reconstruction, primality det-MR/MR-64,
  per-prime order + squarefreeness/Wieferich checks, LTE, sieve re-derivation of
  small `k = 7, 8` factors); anchors `m_19 = m_163 = m_1459 = 4` reproduced
  never assumed; independent term-algebra re-derivation (`p = 19` full `m`
  pattern, `p = 87211`); explicit AS-compositum model cross-checks
  (`sigma(4) = 5`, the norm identity, direct power test per level `k <= 6`).
- boundaries: `k = 7, 8` are per-known-prime only (factordb CF cofactors); an
  `m_p >= 5` example, if any, hides inside those cofactors. Factorization
  provenance: factordb FF entries (fetched 2026-06-12) for `k = 5, 6`,
  re-verified locally by product/primality audit; the P78 and the 43-digit
  `k = 7` prime are PRP-local (factordb marks them proven).

### 2026-06-12: `ogham-2.0`
**Summary:** Ogham first-order abstraction: functions, booleans, ternary, and sorted binders
**Pillars:** scalar ↔ Clifford ↔ Ogham    **Claim level:** engineering — implemented and tested language surface
- surface: `src/ogham/{ast,lex,parse,unparse,error,eval}.rs` now carries the
  v2.0 abstraction layer: `↦`/`~` lambdas as closed AST Function values,
  capture-by-substitution and definition-time beta, tuple application and
  function composition through `@`, per-binder Element/Index/Bool inference
  (including the Gold `a : Index, u : Element` family), Bool literals and
  lazy `and`/`or`/`not`, lazy ternary, Index relations, and ordinary `t`
  bindings outside poly/ratfunc worlds. The live evaluator stores Element,
  Index, Bool, and Function bindings while keeping the pure Rust core free of
  PyO3.
- surface: poly/ratfunc worlds retain the v1.1 element substitution contract
  and add the v2.0 coherence edge: Element ∘ Function yields a Function, and
  Function @ Element evaluates through the same sorted-argument machinery.
  The `;` token now reaches the staged `E_SeqValue` boundary, but sequencing
  remains `ogham-2.1`.
- oracles: the v2.0 slice from `spec/conformance_v2.txt` is merged into
  `spec/conformance.txt`, replacing the four superseded v1.1 reserved-syntax
  vectors. `cargo test ogham_conformance -- --nocapture` passes over the
  expanded corpus, covering capture visibility, self-reference rejection,
  sort errors, Gold traces over `f4`, poly composition, nim-world four-way
  honesty, multivector-valued binders, and definition-time world checks.
- boundaries: no sequencing, let-bodies, continuation lines, recursion,
  game forms, precision literals, or higher-order values. Function equality is
  still `E_FnSort`; Functions may only appear as binding RHS, `@` operands,
  or whole statements.

### 2026-06-12: `nikulin-existence`
**Summary:** Nikulin's even-lattice existence theorem for `(signature, FQM)` pairs
**Pillars:** forms ↔ integral    **Claim level:** standard math made computational + implemented and tested
- surface: `FiniteQuadraticModule::nikulin_existence_report` /
  `nikulin_even_lattice_exists` and the matching `DiscriminantForm` methods implement
  Nikulin theorem 1.10.1 over the bounded finite-table FQM surface. The report type
  (`NikulinExistenceInvariants`) records the signature phase, primary lengths, equality
  cases, determinant square-class checks, and the first obstruction when existence
  fails.
- oracles: tests pin realized ADE discriminant forms, the odd 2-primary boundary
  (`A_1`, where the 2-adic determinant side condition does not fire), an odd
  3-primary equality case, a rank-too-small obstruction, and the even 2-primary
  borderline `U(2)` case.
- boundaries: this decides existence only; it does not enumerate lattices. Like
  `fqm_witt_class`, it is exact up to the finite enumeration budget (`|A| <= 512`) and
  returns `None` rather than guessing past that surface.

### 2026-06-12: `ogham-v1.1`
**Summary:** Ogham function-shaped polynomial and rational-function worlds
**Pillars:** scalar ↔ Ogham ↔ Python    **Claim level:** engineering — implemented and tested calculator surface
- surface: `src/ogham/eval.rs` now dispatches `poly2`/`poly3`/`poly5`/`poly7`,
  `polyint`, and `ratfunc2`/`ratfunc3`/`ratfunc5`/`ratfunc7` as function-shaped
  worlds with reserved `t`, bare constants, bindings, exact equality, live `@`
  substitution/compose, polynomial `%`, exact polynomial `/`, `deg`, and `gcd`.
  The two pending §16 choices are settled: the final names are the literal
  `poly*`/`polyint`/`ratfunc*` rows, and `deg(0)` is `E_Domain`.
- surface: `src/py/scalars.rs` binds `IntegerPoly`, adds `%` to the v1.1
  polynomial classes, and adds `@` eval/compose to the fixed `Fp*Poly` and
  `Fp*RationalFunction` classes; `src/py/catalog.rs` exposes `Poly<Integer>` as
  an engine/divided-power backend.
- oracles: `spec/conformance.txt` carries the new `@world` blocks for `t`
  round-trips, polynomial eval/compose, exact division/remainder, `deg` as an
  Index, primitive `polyint` gcd, monic-divisor errors, ratfunc poles, and
  ratfunc field-boundary errors; `tests/ogham_conformance.rs` passes over the
  expanded corpus.
- boundaries: these are still first-order calculator worlds, not user-defined
  functions. `deg`/`gcd` are polynomial-world functions; ratfunc worlds expose
  equality and substitution but no `%`, order, `deg`, or `gcd`. Precision worlds,
  games mode, and invariant colon-commands remain out.

### 2026-06-12: `ogham-v1`
**Summary:** the v1 Ogham expression language
**Pillars:** scalar ↔ clifford ↔ Python    **Claim level:** engineering — implemented and tested calculator surface
- surface: `src/ogham/` now ships the zero-dependency lexer/parser/AST/unparser,
  typed `OghamError` taxonomy, monomorphised world dispatch over the v1 menu
  (`nimber`, `ordinal`, `surreal`, `omnific`, `integer`, `fp2`/`fp3`/`fp5`/`fp7`,
  `f4`/`f8`/`f16`/`f9`/`f25`/`f27`), metric declarations (`q=`, `b=`, `a=`,
  `grassmann`, and `nimber gold(m,a)`), Element/Index evaluation, bindings,
  relations, `%`, exact integer `/`, factorial, `rev`/`grade`/`even`/`dual`,
  `frob`/`tr`, and checked ordinal Kummer-boundary errors. `examples/ogham_repl.rs`
  provides the colon-command REPL, and the Python module exposes
  `ogham_eval(world, src)`.
- oracles: `tests/ogham_conformance.rs` runs the hand-verified
  `spec/conformance.txt` corpus, including sugar canonicalization, char-2
  Clifford products with independent `q`/`b`, surreal CNF arithmetic, ordinal
  star-literals and Kummer escape, finite-field Frobenius/Wilson checks, exact
  integer division, and relation cells. The Python engine dunder alignment now
  makes multivector `^` raise the Ogham `E_ExpSort` hint instead of delegating
  to wedge.
- boundaries: v1 is a calculator over the original Clifford/scalar world menu;
  the function-shaped `poly*`/`polyint`/`ratfunc*` rows are recorded separately
  in `ogham-v1.1`, while user functions remain `ogham 2.0` work. The Rust
  corpus harness checks the committed hand vectors; corpus expansion/blessing is
  still an operator workflow rather than a rich generated-vector system.

### 2026-06-12: `ogham-backend`
**Summary:** evaluator helper surface for ogham `%`, exact `/`, `@`, `!`, and relations
**Pillars:** scalar (+py bindings)    **Claim level:** engineering — implemented and tested support surface
- surface: `Integer::{divrem,rem,div_exact}` with Euclidean remainders and
  `IntegerDivExactError`; `Surreal::rem` / `Omnific::rem` as the monic-`ω↑e`
  CNF-tail filter; `Poly::compose`; `checked_factorial_i128` plus
  `factorial_in_scalar`; `Ord`/`PartialOrd` for the totally ordered scalar
  worlds (`Integer`, `Rational`, `Surreal`, `Omnific`); `fuzzy()` on `Nimber`
  and `Ordinal`. Python mirrors the new scalar methods where runtime-friendly:
  integer `divrem`/`rem`/`div_exact`/exact `/`/`%`, surreal/omnific `rem`,
  poly `compose`, rich comparisons for `Integer`/`Omnific`, and `fuzzy()` on
  the nim classes.
- oracles: unit tests pin Euclidean signs (`-7 % 3 = 2`), exact-division
  remainder reporting, monic-omega-power rejection/filtering, polynomial
  composition, the `!33`/`!34` i128 boundary, Wilson's `!6 = -1` in `F_7`,
  total-order trait delegation, and nim fuzzy-as-distinctness. Validation:
  `cargo fmt --check`; `git diff --check`; `cargo test`; both clippy passes
  with `-D warnings`; `cargo check --features python`; cold rustdoc with
  `RUSTDOCFLAGS="-D warnings"`; `maturin develop`; `demo.py`; focused Python
  probe over the new methods.
- boundaries: superseded by the `ogham-v1` language entry above. Surreal/omnific
  `%` deliberately rejects non-monic and
  non-omega-power moduli; exact division is not generalized to surreal/omnific
  long division; nim worlds remain unordered (no `PartialOrd`, no `BitOr`
  shorthand).

### 2026-06-11: `ogham-foundations`
**Summary:** the expression-language spec, canonical-ogham Display v2, and host operator alignment
**Pillars:** scalar ↔ clifford (+py and games touchpoints)    **Claim level:** engineering — design contract + behavior-preserving display/operator refactor
- surface: [`spec/ogham.md`](../spec/ogham.md) — the implementation contract
  (canonical-unicode/ASCII-sugar layering `ω ↑ ∧ ⋅` over `w ^ & .`, star as the
  value marker in nim-worlds, structural CNF star-literals, two-sort
  Element/Index typing, the 15-world v1 menu, error taxonomy, work packages) —
  plus [`spec/conformance.txt`](../spec/conformance.txt), the hand-verified
  corpus (incl. `*ω↑3 = *2` and the char-2 polar-pair vector). Display v2
  (spec §9): Ordinal star-wrapped (`*5`, `*ω`, `*(ω + 1)`, `*(ω↑(ω))`),
  Surreal/Fpn explicit `⋅`/`↑`, Poly `x→t`, RationalFunction `(num)/(den)`,
  Multivector `e0∧e1` blades with the ` - ` join rule and the `S::zero()` zero
  rule. Operators (spec §13): `&` is wedge — element-element `BitXor` removed
  and banned on every type (on `Nimber` it would read as XOR = nim-addition);
  `CliffordAlgebra::pow(v, k)`; `x ^ k` (`u128` RHS) scalar power via
  `impl_scalar_ops!`; `Ordinal::nim_pow` checked beside `nim_mul`. Also:
  demo.py `from_i128`→`from_int` repair (missed by the taste-sweep rename).
- oracles: display unit tests pinning the §9 strings per backend;
  `CliffordAlgebra::pow` tests in char 0 and char 2; `nim_pow(ω, 3) = *2`
  (Conway's cube root) + escape-returns-`None`; operator-forwarding tests
  migrated to `&`; full gate green (cargo test 813+16, clippy, cold rustdoc,
  `--features python` check+clippy, demo.py tour end-to-end).
- boundaries: superseded by the `ogham-v1` language entry above. Poly/ratfunc,
  precision worlds, and `{L|R}` game forms are reserved syntax, not shipped.

### 2026-06-11: `taste-sweep`
**Summary:** the taste-audit ledger, played
**Pillars:** scalar ↔ clifford ↔ forms ↔ integral ↔ games (+py touchpoints)
**Claim level:** engineering — behavior-preserving refactor, adversarially reviewed

Thirteen of the fifteen items in [`roadmap/TASTE.md`](TASTE.md) (2026-06-11 audit)
played in one session, on a9's switch decisions: full suffix glossary, outright
Option→Result façade conversion, `e(i)`, IntoIterator-only metric ctors.

- surface: `Scalar::from_int` (the unital ring hom ℤ→R, one spelling, trait-level
  default + per-backend overrides; legacy spellings retired; the duplicate
  `FiniteOddField`/`FiniteChar2Field` embedding methods deleted); `Display` as a
  `Scalar` supertrait with `Debug` delegating (byte-identical output, `{}` now
  works crate-wide); the `…Invariants` record glossary + verb-first façade traits
  (glossary in `src/forms/AGENTS.md`); façade methods return
  `Result<_, ClassifyError>` with leg-faithful variants; `o_type()` →
  `OrthogonalType`; typed Witt/BW arithmetic errors; `WittClass::zero_f2()`;
  `CliffordAlgebra::e(i)` (+`DividedPowerAlgebra::gamma1(i)`; Python keeps `gen`);
  engine encapsulation (`terms`/`metric` pub(crate) + accessors, denormalized
  `dim` field removed); `add_term`/`wedge_terms` dedup; `embed_second` takes the
  left algebra; `Complex64` in its own module; five monoliths split into layered
  directories with public paths frozen (loopy/, game_exterior/, discriminant/
  with phases.rs, lattice/, springer/char2/); Surreal `PartialEq` switched to the
  structural walk (CNF uniqueness argument on the impl).
- oracles: full suite green at every step (797 unit incl. two new `ClassifyError`
  variant pins + the `surreal_structural_eq_matches_value_eq` proptest pin);
  clippy/cold-doc/`--features python` clean on both feature sets; three
  adversarial reviewers (math invariants / motion fidelity / glossary
  completeness) over the full diff, findings resolved.
- boundaries: `py-dunder-pyramid` (2·e_y) and the starred `experiments-as-essays`
  deliberately unplayed — the Python pass is a9's; `Metric::new` keeps its name
  (a9 declined the rename). Python-visible names are stable; six façade bindings
  now surface typed reasons. `demo.py` needs a `maturin develop` rebuild + visual
  pass (display paths changed under it).
- process note (for the next fleet): the worktree-default-branch footgun, bare
  subagent system prompts, and four false fleet self-reports were all caught by
  absolute checks against ground truth (test counts vs dev, clippy vs dev, grep
  vs main), never by re-reading agent summaries. Verify, don't claim — it
  applies to fleets too.

### 2026-06-11: `fqm-witt`
**Summary:** finite-quadratic-module Witt normal forms
**Pillars:** forms ↔ integral    **Claim level:** standard math made computational + implemented and tested
- surface: `FiniteQuadraticModule::{new,cyclic,direct_sum,witt_class}` adds a native
  cyclic-product presentation for nonsingular finite quadratic modules;
  `DiscriminantForm::fqm_witt_class` and `is_fqm_witt_equivalent` expose the
  Wall/Nikulin Witt class of represented discriminant forms. The class is stored
  p-primary as a canonical anisotropic core (`FqmWittClass` /
  `FqmPrimaryWittClass`) with the Milgram/Brown phase retained as a projection.
- oracles: tests pin agreement between native `Z/2, q=1/2` and the lattice `A_1`
  discriminant form, separation of `A_1` and `E_7`, split reduction of
  `⟨1/2⟩ ⊕ ⟨3/2⟩`, split reduction of the actual lattice sum `A_2 ⊕ E_6`, and
  compatibility with the older `fqm_gauss_phase` projection.
- boundaries: this is an exact bounded finite-table normal form (`|A| <= 512`,
  plus an explicit generator-tuple budget); it returns `None` rather than
  truncating. The lattice-existence successor is shipped separately as
  `nikulin-existence`.

### 2026-06-11: `fqm-gauss-phase`
**Summary:** p-primary finite-quadratic-module phase projection
**Pillars:** forms ↔ integral    **Claim level:** standard math made computational + implemented and tested
- surface: `DiscriminantForm::fqm_gauss_phase` returns `FqmGaussPhase` with
  p-primary `FqmPrimaryPhase { prime, order, exponent, phase_mod8 }` factors and the
  total Milgram/Brown `Z/8` phase; `milgram_signature_mod8_fqm` exposes the total.
  The existing `GaussSum`/`milgram_signature_mod8` float route stays as an oracle.
- oracles: tests pin the mixed-primary `A_1 ⊕ A_2` factorization, extend the 2-primary
  phase past the old Brown-only `Z/2` slice on `A_3`, cover odd torsion on `E_6`, and
  cross-check the FQM phase against exact lattice signature, the Conway-Sloane genus
  oddity route, and the legacy float phase across the ADE zoo.
- boundaries: this is the Gauss-sum phase projection of the finite quadratic module,
  not Wall/Nikulin/Kawauchi-Kojima's full generator-and-relation Witt normal form.
  It still enumerates represented discriminant groups (`FQM_GAUSS_GROUP_CAP = 4096`)
  and uses the principal embedding only after an exact cyclotomic shape check chooses
  between the two square-root branches; the FQM normal-form successor is the
  `fqm-witt` entry above.

### 2026-06-11: `game-clifford-checked`
**Summary:** quotient-compatible integer Clifford data on game generators
**Pillars:** clifford ↔ games ↔ Python    **Claim level:** implemented and tested
- surface: `GameClifford::{new,free,with_relation_search,with_quadratic_data}` wraps
  the integer Clifford engine on a chosen game-generator tuple; explicit or discovered
  `GameRelation` rows are accepted only after the game-group relation evaluates to
  zero and the supplied integer `q`/polar data makes that relation null and
  polar-radical. Quotient-aware `reduce`, `add`, `scalar_mul`, `mul`, `wedge`,
  `is_zero`, and `value_of_grade1` mirror `GameExterior`; PyO3 exposes
  `GameClifford` with the same checked constructors.
- oracles: tests pin free Clifford anticommutators, rejection of `Q(*) != 0` and
  nonzero pairings under `2* = 0`, accepted torsion vanishings with `2*(e_* e_up)=0`,
  duplicate-generator compatibility, and bounded relation search. The Rust tour and
  `demo.py` include the accepted/rejected `2* = 0` examples.
- boundaries: this is an integer-valued checked deformation engine, not a Clifford
  algebra over arbitrary games and not a proof that the quadratic data is game-native.
  The natural-source and torsion-target questions remain in `OPEN.md` §2.

### 2026-06-11: `loopy-partizan`
**Summary:** finite Left/Right loopy outcome engine
**Pillars:** games ↔ Python    **Claim level:** implemented and tested
- surface: `LoopyPartizanGraph`, `LoopyPartizanOutcome`, `LoopyWinner`;
  `LoopyValue::{PlusMinus,Tis,Tisn,OnsideOffside}` plus exact starter-pair
  `outcome`, classical `partizan_outcome`, `sides`, negation, conservative
  comparison, and partial sums; `LoopyNimCertificate::{recovery_condition_holds,
  recovery_blockers}` for the checked finite recovery condition behind additive
  finite-nimber claims.
- oracles: loopy unit tests pin the classical `P/N/L/R` embedding, the repo
  `tis={0|tisn}` / `tisn={tis|0}` mixed draw classes, impartial agreement with
  `kernel::outcomes`, and recovery blockers for finite positions with `Side`
  options; Python feature build exposes the new classes.
- boundaries: this is finite retrograde analysis and finite/certified sidling,
  not full loopy-game equality; mixed classes are kept as starter pairs instead
  of being forced into `P/N/L/R/Draw`; sums outside the represented catalogue
  return `None`.

### 2026-06-11: `py-waves`
**Summary:** Python parity for waves J/K/M/N/O
**Pillars:** scalar/forms/games ↔ Python    **Claim level:** implemented and tested
- surface: Python now exposes `lexicode`/`lexicode_naive`/`lexicode_bounded`,
  `NimLexicode`, Brown invariants (`BrownResult`, `brown_f2`, `double_f2`,
  `DiscriminantForm.brown_invariant`), discriminant-form isomorphism checks,
  `NewtonPolygon`/`newton_polygon`/`tropicalize`, Scharlau `transfer_diagonal`,
  rational and full `Q/Z` Brauer classes (`Brauer2Class`, `BrauerClass`,
  `cyclic_algebra_invariant`), Milnor residues (`global_residues`,
  `global_residues_ff`), and function-field constant-extension invariants.
- oracles: the Rust math oracles remain the source of truth; `demo.py` now includes a
  py-waves parity rung covering each newly bound surface, and the Python feature check
  compiles the wrappers.
- boundaries: bindings follow the existing fixed-dispatch policy. Const-generic
  open families remain Rust-only unless represented by an existing Python slice; the
  char-2 transfer story still belongs to the Artin-Schreier layer, not
  `transfer_diagonal`.

### 2026-06-11: `fpn-gen`
**Summary:** generated finite-field reduction polynomials
**Pillars:** scalar    **Claim level:** standard math made computational + implemented and tested
- surface: `Fpn<P,N>` now supports every prime `P` and positive `N` whose field order
  `P^N` fits in `u128`; every extension field uses a cached deterministic search for
  the first monic irreducible reduction polynomial certified by Rabin's irreducibility
  test. The public `ReductionPolynomialKind::GeneratedIrreducible` tag records
  generated provenance, including through the Python enum wrapper.
- oracles: tests keep the curated `F_4/F_8/F_16/F_9/F_25/F_27` field-law sweeps, add
  the old small Conway/irreducible rows as test-only scan-order oracles, add a
  generated `F_32` field-law sweep, verify generated `F_32` irreducibility metadata,
  and check a primitive element in generated `F_128`.
- boundaries: this does not claim Conway-compatible embeddings for generated rows,
  and the `u128` payload model still rejects fields whose order overflows `u128`
  (for example `Fpn<2,128>`).

### 2026-06-11: `ordinal-principled`
**Summary:** Kummer carries from `m_u` only
**Pillars:** scalar    **Claim level:** standard math made computational + implemented and tested
- surface: `scalar/big/ordinal/tower.rs::alpha_ordinal` now reconstructs each shipped
  Kummer carry by computing `f(u)=ord_u(2)`, recursively computing DiMuro's `Q(f(u))`,
  assembling the corresponding `χ`-sum, and nim-adding the finite excess integer
  `m_u`. The only production row data left for the staged tower is
  `finite_excess(u)` through `u=47`.
- oracles: tests separately pin every shipped row's `f(u)`, `Q(f(u))`, finite `m_u`,
  and resulting `α_u`; the existing cube/quintic/septic/`α_47` landmarks and field-law
  sweeps remain green.
- boundaries: the operational boundary is unchanged. A carry needing `m_53` or beyond
  still returns `None`; computing or certifying new finite excess integers remains the
  open/research step.

### 2026-06-11: `nim-lexicodes`
**Summary:** lexicodes over nim alphabets
**Pillars:** games ↔ integral    **Claim level:** standard math made computational + implemented and tested
- surface: `games::nim_lexicode_naive` / `nim_lexicode_naive_bounded` build literal
  greedy lexicodes over the alphabet `{0, ..., 2^k-1}` and return a `NimLexicode`
  record with packed/decoded codewords, nim-add closure verification, F2-dimension,
  Fermat-base detection, and coordinatewise nim-scalar closure checks.
- oracles: tests verify the base-`2^k` repetition lexicodes are closed under
  coordinatewise nim-addition; base `4` and `16` are closed under nim-scalar
  multiplication, while base `8` fails exactly because the alphabet is not a finite
  nim-field.
- boundaries: the binary optimized `lexicode`/`BinaryCode` route is unchanged. The
  q-ary path is intentionally literal/budgeted and does not implement the deeper
  Conway-Sloane turning-game realization, which remains the separate `lexicode-game`
  TODO.

### 2026-06-11: `subfield-detect`
**Summary:** finite ordinal-nimber subfield detection
**Pillars:** scalar ↔ forms    **Claim level:** standard math made computational + implemented and tested
- surface: `Ordinal::finite_subfield_degree`, `scalar::ordinal_finite_subfield_degree`,
  and `scalar::ordinal_common_finite_subfield_degree` detect the minimal represented
  finite field `F_{2^m}` by generator support plus Frobenius minimization. The forms
  side exposes `forms::ordinal_metric_finite_subfield_degree` and routes
  `arf_ordinal_finite`, ordinal Witt classes, ordinal Brauer-Wall classes, and
  ordinal isometry through the detected/common degree.
- oracles: tests pin finite-nimber degrees, `ω` as degree 6, `ω^3` as degree 18,
  `ω^ω` as degree 20, `ω^(ω^2)` as degree 42, common-degree lcm behavior, inversion
  in a detected non-`F_64` field, ordinal Arf classification past the old `F_64`
  window, and rejection at `ω^(ω^ω)`.
- boundaries: detection is limited to the source-verified staged tower and the
  shipped finite Kummer excess table (`m_u` through `47`); genuinely transfinite
  ordinal-nimber metrics still return `None` for finite Arf/Witt/Brauer-Wall
  classification.

### 2026-06-11: `brown-reduce`
**Summary:** Brown invariant by reduction, not enumeration
**Pillars:** forms    **Claim level:** standard math made computational + implemented and tested
- surface: `forms/char2/brown.rs::brown_f2` now splits off the polar radical and
  reduces the nonsingular `ℤ/4`-quadratic core into odd lines and even symplectic
  planes, adding their Brown phases in `ℤ/8`; `BROWN_MAX_ENUM_RANK` is gone, so the
  public route no longer panics above rank 26.
- oracles: the old exact Gauss-sum enumeration route is retained as a test-only
  checker; tests exhaust every four-dimensional input, compare the former rank-26
  budget edge against enumeration, preserve the Arf-doubling/additivity/radical
  checks, and verify a rank-40 form reduces past the old ceiling.
- boundaries: the `u128` bitmask surface still bounds represented dimensions to
  `n ≤ 128`; Brown's symmetric polar category remains distinct from the Clifford
  char-2 alternating polar, with `double_f2` still the explicit bridge.

### 2026-06-11: `cyclic-trace-form`
**Summary:** cyclic algebra `Trd(z²)` trace form
**Pillars:** forms ↔ scalar    **Claim level:** standard math made computational + implemented and tested
- surface: `forms/trace_form.rs::cyclic_algebra_trace_form::<E>(&a)` builds the
  literal cyclic-algebra quadratic trace form `T_A(z)=Trd_A(z²)` for
  `A=(E/F,σ,a)=⊕ E·u^i`, ordered by the `E·u^i` lines. The `u^0` and, when present,
  `u^{n/2}` self-blocks reuse `assemble_twisted_form`; the `i`/`n-i` line pairs are
  pure polar blocks. The existing `MAX_BASIS_DIM=128` cap gives the boundary
  `[E:F]^2 ≤ 128`.
- oracles: over `Q(i)/Q`, the degree-2 trace form is pinned as
  `⟨2,-2,2a,2a⟩`, and a grid check verifies the honest reduced-norm relation
  `Trd(z²)=Trd(z)^2-2·Nrd(z)` against the shipped `⟨1,1,-a,-a⟩` norm form. Over
  `F_27/F_3`, the `u`/`u^2` cross block is checked to be Witt-hyperbolic.
- boundaries: this is not the reduced norm for general `n` (nor equal to it for
  quaternions); it is the quadratic trace companion named in CODA K §6(c).

### 2026-06-11: `milnor-d2`
**Summary:** dyadic cell of Milnor's exact sequence
**Pillars:** forms    **Claim level:** standard math made computational + implemented and tested
- surface: `forms/witt/milnor.rs::global_residues` now includes the `p=2` component
  of Milnor's residue map. The dyadic residue uses Milnor's hand convention, not the
  odd-prime Springer residue: a diagonal line contributes iff its `2`-adic valuation
  is odd, landing in the `W(F_2) ≅ Z/2` carrier
  `WittClassG::Char2 { field_degree: 1, arf }`.
- oracles: tests pin `⟨2⟩`, `⟨1,2⟩`, and `⟨-2⟩` as nonzero dyadic residue classes,
  verify `⟨2,-2⟩` cancels, check mixed support such as `⟨6⟩`, and cross-check
  reconstruction against `try_is_isotropic_q` for `⟨2⟩` vs `⟨8⟩` and `⟨2⟩` vs `⟨1⟩`.
- boundaries: this is only the dyadic `W(ℚ)` cell; the equal-characteristic
  `F_q(t)` twin shipped separately as `milnor-ff`.

### 2026-06-11: `milnor-ff`
**Summary:** split Milnor residues over `F_q(t)`
**Pillars:** forms    **Claim level:** standard math made computational + implemented and tested
- surface: `forms/witt/milnor.rs::global_residues_ff` implements the split
  equal-characteristic Milnor map for odd constant fields:
  `W(F_q(t)) ≅ W(F_q) ⊕ ⊕_π W(F_q[t]/π)`. The first component is the constant-field
  class selected by the even-valuation layer at the degree place `∞`;
  `FunctionFieldMilnorResidues<S>` records that class plus the finite vector of
  nonzero second residues at monic irreducible places, using the exact
  `local_global/function_field.rs` valuation, residue-unit, and residue-character
  helpers.
- oracles: tests pin constant forms, the `t` place, nonsquare constants, a degree-2
  irreducible place over `F_5`, square-multiple invariance, hyperbolic cancellation,
  and radical-entry rejection.
- boundaries: odd constant fields only (`FiniteOddField`); characteristic-2
  function fields stay on the separate Artin-Schreier/Aravire-Jacob layer, and this
  does not claim tame or wild norm-residue symbols beyond the second-residue Witt
  map.

### 2026-06-10: `echo-solver`
**Summary:** the echo-fifo+dummy adversarial review: CONFIRM
**Pillars:** games ↔ forms    **Claim level:** implemented and tested
- surface: `experiments/echo_solver.py` (stdlib-only, maintained) — direct full-state
  σ-explicit solvers for the echo family (`fifo_value`, `ko_value`), the parameterized
  prose-reading sweep that *failed* (all 80+ readings of the old §8.3 prose die at
  m=8 — the prose misdescribed the rule), nim/mex self-test, torsor and
  decision-liveness instrumentation, and the full-sweep stage.
- verdict: the formerly **unverified** echo-fifo+dummy `m = 8` exactness claim is
  **CONFIRMED** — 391,680/391,680 checks (765 scaled Gold forms × 256 positions ×
  both stances), zero misses, re-derived with no decomposition and no isomorphism
  caching. The m=4 family is 30/30 with the dummy and only 15/30 without it; the
  §8.2 echo-ko table was independently reproduced (σ-explicit, miss x=224 included)
  and shown to be the σ=1-stance face of a stance-asymmetric rule.
- oracles: explicit no-memo tree enumeration (m=4 exhaustive + m=8 small supports);
  the original `direct_fifo_value` executed verbatim (1,920 agreeing solves); the
  Turning-Corners mex recurrence for the nim product; a Codex cross-run of every
  stage including nim_mul vs the original probe on all 65,536 pairs.
- boundaries: the verified realizer is **σ-valued** (forced terminal charge, not a
  P-set); normal/misère/loopy recasting and the even-`a` diagonal lemma remain open
  (`OPEN.md` tis (§1)); the bounded-*window* blocker conjecture is untouched (FIFO
  memory is unbounded); goldarf §8.3's old prose rule description was corrected in
  the same pass.

### 2026-06-10: `linking-reduction`
**Summary:** the echo-fifo+dummy mechanism, reduced and screened
**Pillars:** games ↔ forms    **Claim level:** standard math (reductions) + implemented and tested (screens) + open (general n)
- surface: `experiments/linking_game.py` (stdlib-only, maintained) — the abstract
  odd-close parity game, validated against `echo_solver.fifo_value` through the
  `SynthForm` bridge; all-iso-classes rigidity/Bad screens; the strict menu
  verification of the two-mode defender strategy. Writeup: goldarf §8 "the linking
  reduction and the general-m theorem" (`\label{sec:linking}`).
- reductions (proven, machine-validated): FIFO ⇒ no nesting ⇒ linked = overlap;
  `D = σ ⊕ und` flips only on odd-`deg_U(front)` closes (the σ-game IS the
  odd-close parity game); ko localizes to fronts opened onto an empty queue and
  passes occur only after `U = ∅` (irrelevant to the flip fight).
- verified: the **linking theorem** (isolated coin ⇒ flips forced even ⇒
  `σ = |E| mod 2 = Q(x)` on Gold boards, i.e. exactness for all m) holds on every
  graph iso class `k ≤ 7` + dummy, both seats (1,044 classes at k=7) — beyond
  Gold-arising boards. No-dummy Bad census `{3:1, 5:4, 7:34}`: all
  mover-controlled, none with an isolated vertex, 33/34 dominated at n=7; all
  even-n boards rigid. The domination device (queue empty, `v` dominating an
  even nonempty remainder ⇒ forced flip in two plies) is the unique local
  obstruction; the dummy kills it at every root. Two-mode defender strategy
  (prevention/debt menus) strict-verified `k ≤ 7`, both seats —
  menu-existential semantics (a Codex review exhibited a losing in-menu choice).
- boundaries: the general-n proof is **open** — parity-local invariants provably
  insufficient (feature-mining inconsistency); the working architecture
  (Codex spar, thread 019eb4ff-695b-7762-97e8-c0bea66c4e7e) is firewall
  segmentation + mutual no-debt/one-debt induction + certificate-depth
  completeness, with the poison transition as the hard obligation.
