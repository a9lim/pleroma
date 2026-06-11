# Taste audit — TASTE.md (the aesthetic ledger)

A structural/stylistic read of the Rust core (with a broad-stroke coda on the
Python surface), 2026-06-11. This is an *opinion* document — claim level:
**interpretation**, one reviewer's eye — but every observation below is checked
against the actual source, not vibes. Items are valued like
[`TODO.md`](TODO.md): a game value `g` on a pillar blade `e_B`
(`e_s` scalar, `e_c` clifford, `e_f` forms, `e_i` integral, `e_g` games,
`e_y` py). Numbers ≈ focused days; `±n` means the real decision is a9's
(API-churn scope calls, mostly); `↑` is worth less than any number but
strictly positive; `*n` is real, on-thesis, unscheduled.

Renames here are *internal* unless marked otherwise — the Python surface names
its classes independently, so most of this churn never reaches a user.

---

## Status — played 2026-06-11 (same session as the audit)

| item | outcome |
|---|---|
| `int-embedding-babel` | **played** — `Scalar::from_int`, legacy spellings retired, duplicate field-trait embeddings deleted |
| `debug-as-display` | **played full** — `Display` supertrait, `Debug` delegates, byte-identical output |
| `monolith-modules` | **played** — all five splits, paths frozen |
| `record-suffix-zoo` | **played full** (a9's call) — `…Invariants` glossary + verb-first traits, glossary in forms/AGENTS.md |
| `none-reasons` | **played as conversion** (a9's call) — `Result<_, ClassifyError>`, leg-faithful variants, pinned |
| `engine-encapsulation-split` | **played** — pub(crate) fields + accessors, `dim` field removed |
| `terms-idiom-dup` | **played** — `add_term`/`wedge_terms`, twins collapsed |
| `metric-ctor-ergonomics` | **played IntoIterator-only** (a9 declined the `new` rename) |
| `stringly-edges` | **played** — `o_type()`/`OrthogonalType`, typed try_* errors, `zero_f2()` |
| `gen-keyword` | **played** — `e(i)` / `gamma1(i)`; Python keeps `gen` |
| `complex64-placement` | **played** — own module, shadow documented |
| `mv-context-asymmetry` | **played** — `embed_second` takes the left algebra |
| `surreal-eq-cost` | **played** — structural eq + CNF-uniqueness argument + proptest pin |
| `py-dunder-pyramid` | **unplayed** — deferred to a9's Python pass |
| `experiments-as-essays` | **unplayed by design** — the star stays a star |

Process record: `roadmap/DONE.md` → `taste-sweep` (incl. the fleet-failure
archaeology). The audit prose below is left as written — findings text describes
the pre-sweep tree.

---

## What the codebase does well (0 — already played, listed so they don't get refactored away)

These aren't filler; they're the things a cleanup pass could accidentally
destroy, so they go on record first.

- **The narrated-library doc culture.** Module docs are essays with theses
  (`scalar/mod.rs`'s "any number" table, `engine.rs`'s char-2 manifesto), the
  per-pillar AGENTS files carry intent *and* the "looks like a bug, isn't"
  ledgers, and claim levels are enforced in prose. Most research code documents
  what; this documents *why* and *why not*. It is the single best thing about
  the repo.
- **Boundary honesty as a house aesthetic.** Partial operations return `Option`
  and the `None`s mean something (`Surreal::inv` on non-monomials, `Rational::sqrt`
  on non-squares, the p-adic `Option<Option<_>>` precision contract); `Ordinal`
  *omits* owned `*` rather than panic prettily; capability traits
  (`Valued`/`ExactRoots`/`ResidueField`/exactness markers) are opt-in, never
  `Scalar` supertraits. The refusal to fake totality is consistent and lovely.
- **Decisions with recorded failure evidence.** The operator policy
  (`impl_scalar_ops!` concrete-only, "70+ generic sites broke when tried") is
  the model: a tradeoff was *measured*, decided, and documented where the next
  editor will trip over it. Same for `cnf.rs` as a shared function-not-type.
- **Oracle discipline.** `reduce_word` kept as an independent `cfg(test)` oracle
  for the Chevalley product; proptest fuzz with the precision worlds honestly
  excluded; one verified Galois engine shared by `Nimber` and `Fpn` through the
  `FiniteField` trait. "Verify, don't claim" is actually practiced.
- **Literate mathematical naming.** `Metric::grassmann`, `Surreal::omega_pow`,
  `Game::star`/`nim_heap`/`fuzzy`, `springer_decompose_laurent` — objects are
  named what mathematicians call them, and the No↔On₂ mirror keeps constructor
  names parallel across the symmetry (`Surreal::omega()` / `Ordinal::omega()`).
- **Substrate kept in its place.** `linalg/` is `pub(crate)` on purpose;
  `py/catalog.rs` makes the bound-monomorph set *data* in one manifest. Both are
  the right call and both are documented as policy.

The critiques below are sharper *because* the baseline is this high — most of
them are places where the codebase fails to meet its own standard, which is the
most useful kind of finding.

---

## numbers — buildable now

### 1·e_s: `int-embedding-babel`
**The canonical map ℤ → R is spelled six ways across a table whose entire
thesis is uniformity.** `Rational::int`, `Surreal::from_int`,
`Omnific::from_int`, `Qq::from_int`, `WittVec::from_int`, `Qp::from_i128`,
`Fp::new(i128)`, `Zp::new(i128)`, `Ordinal::from_u128`, `Fp::from_u128`,
`Tropical::int`. Every commutative ring has exactly one unital embedding of ℤ;
the scalar pillar's whole point is that these worlds are *one table*, and the
single most uniform mathematical operation is the least uniform name in the
crate. Compounding it: there are **zero** `From`/`TryFrom` impls in the entire
crate — the std conversion vocabulary is simply unused.
Recommendation: add `Scalar::from_int(n: i128) -> Self` as a trait method with
a default double-and-add impl over `one()` (backends override for speed), keep
the ergonomic inherent shortcuts as aliases, and let generic forms/code stop
hand-rolling small-integer lifts. `impl From<i128>` for the concrete worlds
where it's total is then a free one-liner each.

### 1·(e_s∧e_c∧e_f): `debug-as-display`
**The crate's entire human-facing rendering layer lives in hand-written `Debug`
impls; there is exactly one `Display` impl in the whole crate**
(`GameCliffordError`). Seven types carry a bespoke `pub fn display(&self) ->
String` (Multivector, WittClass, CliffordType ×2, OddCharType,
FiniteFieldClass, Game), and `Multivector::display` renders coefficients via
`{:?}`. The root cause is the trait bound: `Scalar: Debug`, not `Display`. I
understand why it happened — `assert_eq!` prints `Debug`, so pretty-Debug makes
test failures readable — but the cost is real: `format!("{x}")` doesn't compile
for any scalar, no *structural* dump exists anywhere (debugging a misbehaving
`Surreal`'s actual term tree requires reading the pretty form back in your
head), and `display()`-returning-String is Java idiom, not Rust (no `{}`
interpolation, no `ToString`, no trait-object printing).
Minimal fix, fully additive: `impl Display` for the scalars by delegating to
the existing `Debug` bodies, add `Display` to the `Scalar` bounds, implement
`Display` for `Multivector` and the record types, keep `display()` as thin
aliases for the Python layer. The *full* inversion (derive `Debug`, move pretty
into `Display`) changes test-failure output crate-wide — that's a `±2` switch
a9 should call deliberately, and "don't" is a defensible answer; the additive
half is just hygiene.

### 1·(e_f∧e_g): `monolith-modules`
**The crate has a demonstrated, excellent splitting pattern — and five files
that ignore it.** `clifford/engine/` splits at ~100–500 lines by concept
(basis/metric/product/algebra/multivector/inverse/terms); `nimber/` likewise
(arithmetic/artin_schreier/galois). Meanwhile `loopy.rs` (1230 lines — AGENTS
itself describes it as "four layers"), `springer/char2.rs` (1424),
`integral/discriminant.rs` (1483: GaussSum + Complex64 + DiscriminantForm +
two phase types), `game_exterior.rs` (1184: Λ-engine + relation certificates +
GameClifford deformation), and `integral/lattice.rs` (1342) each hold 3–5
separable concerns. This isn't line-count pedantry — the small-module
discipline is what makes the engine readable, and the forms/games pillars'
hardest files are exactly the ones that didn't get it. Mechanical, low-risk,
pattern already in-house.

### ±1·e_f: `record-suffix-zoo`
**The forms layer names the same kind of object — "the record a classifier
returns" — under at least five suffix conventions.** `ArfResult` and
`BrownResult`; `CliffordType`, `RationalCliffordType`, `OddCharType`;
`FiniteFieldClass`, `SymplecticClass`, `VersorClass`; `HermitianSignature`;
the `*Decomp` family. The damage: `…Class` sometimes means a literal element of
a classifying group (`WittClass`, `BrauerWallClass`, `BrauerClass` — the *good*
usage, where `try_add` is a group law) and sometimes just a report
(`FiniteFieldClass` is an enum of reports). A reader can't tell from the name
whether a type carries algebra or prose. Proposed glossary: `…Class` =
element of an actual group/pointed set (keeps Witt/Brauer/BW); `…Decomp` =
decomposition (already consistent); everything else — the leg reports — picks
ONE suffix (`…Invariants` reads best to me; `…Type` is the incumbent with
seniority). Fold in the façade-trait word-order wobble while there:
`ClassifyForm` (verb-object) vs `WittClassify`/`IsometryClassify`/
`BrauerWallClassify` (object-verb) vs `WittDecompose` — five traits, one job,
two grammars. A switch because it's rename churn across forms/ + py/forms.rs
docs; the blast radius is internal but wide.

### ±1·e_f: `none-reasons`
**`Option` is doing error-enum work in the classifier façade.** 278
`Option`-returning fns in forms/ vs 17 `Result`s. For genuinely partial math
(`inv`, `sqrt`) bare `None` is the honest house style and should stay. But
`metric.classify() -> None` means *one of*: general-bilinear metric, singular
polar form, outside the represented window, diagonalizer failed — the AGENTS
files then carry the disambiguation in prose that the type system was offered
and declined. The crate already knows the better pattern:
`WittClassError { GeneralBilinearMetric, Singular {…} }` is exemplary and sits
right there in `witt/class.rs`. Recommendation: small `#[non_exhaustive]`
reason enums for the façade entry points where `None` is ≥3-valued
(`classify`/`witt_class`/`bw_class`), nothing else. A switch because
Option-plus-docs is arguably *also* a deliberate aesthetic — but I'd argue the
docs are currently load-bearing in a place types are cheaper.

## halves — an afternoon each

### ½·e_c: `engine-encapsulation-split`
**Three core engine types, three different encapsulation postures, one rule.**
`Metric` is the hard-rule type ("never the bare struct literal") and is
properly guarded: `pub(crate)` fields, validating constructors, accessors,
`into_parts`. But `Multivector.terms` is a bare `pub` field on a type whose
invariant ("zeros never stored" — `is_zero()` *is* `terms.is_empty()`) any
downstream user can silently violate with a struct literal; and
`CliffordAlgebra { pub dim, pub metric }` validates `dim == metric.dim()` in
`new` then lets anyone mutate either field. `dim` is moreover *denormalized* —
it always equals `metric.dim()`; storing both is carrying a proof obligation
for free. Fix: `terms` → `pub(crate)` + accessor (the py layer is in-crate;
nothing breaks), drop the `dim` field for a `dim()` delegating method, and the
crate's stated invariants become structural instead of social. While in there:
state the operator-vs-context-method policy for `Multivector` the way
`impl_scalar_ops!`'s doc states it for scalars — right now `a + b` /
`alg.add(&a,&b)` and `a ^ b` / `alg.wedge(&a,&b)` coexist with no canonical
choice on record.

### ½·e_c: `terms-idiom-dup`
**The "add into entry, remove if zero" dance is hand-inlined five-plus times,
and `BitXor` duplicates `wedge` verbatim.** `terms.rs` exists precisely to hold
this idiom but only ships `merge`/`scale`; the single-term form
(`entry → or_insert(zero) → add → remove-if-zero`) is replicated in
`Multivector::bitxor`, `CliffordAlgebra::wedge`, `contract_vec_blade`,
`vec_times_blade`, and `merge` itself. Worse, `Multivector::bitxor` and
`CliffordAlgebra::wedge` are the *same loop body copy-pasted* — a divergence
bug waiting for whichever one gets fixed first. Fix: `add_term(&mut BTreeMap,
blade, coeff)` in `terms.rs`, a shared `wedge_terms` free fn both call. Pure
deletion; the associativity suite already pins behavior.

### ½·e_c: `metric-ctor-ergonomics`
**Every off-diagonal metric in the repo is built with a three-line `BTreeMap`
ritual.** `let mut b = BTreeMap::new(); b.insert((0,1), x); Metric::new(q, b)` —
in the engine tests alone this appears a dozen times, and it's the documented
public path. Accept `impl IntoIterator<Item = ((usize, usize), S)>` in
`new`/`general` and call sites collapse to
`Metric::new(q, [((0,1), x)])`. Secondary nit, same constructor family: `new`
being the *middle* generality (q,b) while `general` is the actual general case
reads backwards — `new` is conventionally either the primary or the most
general constructor. Renaming `new → with_polar` (keeping `new` as alias a
deprecation-cycle) would make the family self-describing:
`diagonal / grassmann / with_polar / general`. The iterator change is free;
the rename is optional polish.

### ½·e_f: `stringly-edges`
**Three small sharp edges where the strong-typing discipline lapses.**
(1) `ArfResult.o_type: &'static str` — an invariant ("O+"/"O−") carried as a
string in the crate's flagship char-2 record, *and* it's derivable from the
Arf bit; an `enum OType { Plus, Minus }` (or dropping the field for a method)
is strictly better. (2) `try_add(&self, …) -> Result<_, &'static str>` — three
string-error methods in `witt/class.rs` sitting in the same file as the
properly-typed `WittClassError`; two error styles within sixty lines.
(3) `WittClass::zero()` silently means *over F₂* — sum it with a class over
F₄ and you get a runtime error from the identity element; the parameterless
zero of a field-indexed group family is a footgun (`zero_over` exists and is
the honest one; `zero()` could go, or become `zero_f2()`).

## ups — worth less than any number, still strictly positive

### ↑·e_c: `gen-keyword`
`CliffordAlgebra::gen` collides with the `gen` keyword reserved in Rust
edition 2024 — an edition migration rewrites every call site to `r#gen(i)`,
which is ugly enough to count as breakage. The repo is 2021 so nothing is
on fire, but the rename is better done by choice than by `cargo fix`:
`alg.e(i)` is shorter, matches the display language (`e0e1`) and the math
(`e_i`), and frees the keyword. (`blade(&[i])` already exists as the general
form.)

### ↑·e_i: `complex64-placement`
A hand-rolled `pub struct Complex64` lives inside
`forms/integral/discriminant.rs` — general-purpose float-complex machinery
embedded in (and re-exported from) a Weil-representation module, sharing its
name with `num_complex::Complex64` (a deliberate-or-not shadow worth a doc
line either way). Dependency-free is the right call; the placement isn't —
it's substrate, and the crate has a substrate floor (`linalg/`, or a sibling
util) where it would stop looking like part of the discriminant-form theory.

### ↑·e_c: `mv-context-asymmetry`
The graded-tensor embedding API is asymmetric in a way that pushes
bookkeeping onto the caller: `embed_first(&self, v)` ignores `self` entirely
(it's a clone), while `embed_second(&self, v, shift)` takes the *first
factor's dimension* as a raw `usize` the caller must remember
(`alg.embed_second(&right_v, left.dim)`). Either have `graded_tensor` return
a small product type that remembers the split and owns both embeddings, or at
minimum take the left algebra by reference instead of a bare integer.

### ↑·e_s: `surreal-eq-cost`
`PartialEq for Surreal` routes through `self.sub(other).sign()` — a full
subtraction (clone + canonicalize + recursive exponent comparison) per `==`,
although the representation is already canonical by construction (every
constructor canonicalizes; exponents recursively so). If canonical-form
structural equality is provably equivalent — and the constructors' discipline
suggests it is — `==` becomes a cheap structural walk and `canonicalize`'s own
exponent comparisons stop paying the subtraction tax. Worth a proof-comment
either way: right now the value-eq impl silently implies the representation
*isn't* trusted to be canonical, which contradicts the module doc.

---

## the Python-facing side — broader strokes

### 2·e_y: `py-dunder-pyramid`
The binding layer is the largest code in the repo (`scalars.rs` 5823 lines,
`forms.rs` 5458, `games.rs` 2531) and it's *better* organized than its size
suggests — `catalog.rs` as the single manifest is genuinely good architecture.
But the macro layer stops one level short: eleven per-family pyclass macros
(`prime_field_pyclass!`, `qp_pyclass!`, `laurent_pyclass!`, …) **each
re-implement the shared `Scalar` dunder surface** (`__add__`/`__radd__`/
`__sub__`/`__mul__`/`__neg__`/`__eq__`/`__repr__`/`zero`/`one`/`is_zero`/…) —
23 hand-rolled `fn __add__` definitions in one file. The AGENTS file even *names* the
concept ("the shared runtime `Scalar` surface") that the code never factored.
One `scalar_dunders!($ty, $parse, $wrap)` macro invoked *by* the family macros
would delete on the order of a thousand lines and make "add a method to every
scalar" a one-site edit. Same move then splits `scalars.rs` by family into a
`py/scalars/` directory, which the engine/ pattern already licenses.

### *1·e_y: `experiments-as-essays`
The honest read of `experiments/` + `demo.py`: **much better than advertised.**
The module docstrings are publication-grade (the `misere_kernel.py` header is a
small literature review with theorem citations and an honest caveat), the
stdlib-only/no-venv constraint is a stated feature and a good one, and the
probes are self-contained on purpose. The mess is real but shallow: each script
re-grows its own CLI conventions (`selftest` vs `all 5` vs bare argv), its own
mini-engines, and its own output formatting; there's no shared `_lib.py` for
the recurring outcome-solver/table-printing/assert-census patterns; and
`demo.py` at 1004 lines is a tour that can only be run whole. Star-valued, not
numbered, deliberately: consolidating research probes into a framework is the
kind of cleanup that can *destroy* value (self-containedness is why the
adversarial-review harnesses are trustworthy — `echo_solver.py` being one file
IS its audit story). The buildable kernel, if any: a tiny shared CLI/reporting
helper, adopted only by *new* probes, and a `demo.py` section index. Leave the
verified harnesses byte-stable.

---

## the disposition (one paragraph, hat off)

This is a codebase with an unusually strong *macro*-aesthetic — the pillar
symmetries, the honesty boundaries, the narrated docs are all top-percentile —
and a merely good *micro*-aesthetic that hasn't had the same love: names and
encapsulation postures drift between files in ways the math never does. The
fix-shape is correspondingly small: almost everything above is an afternoon-to
two-day item, none of it touches mathematical content, and the two switches
(`record-suffix-zoo`, `none-reasons`) are the only ones where reasonable
people disagree. If I had to play one move first it would be
`int-embedding-babel` — it's cheap, it's on-thesis (the table should *look*
like a table), and it makes several other items (generic forms code, py
parity) slightly easier downstream.
