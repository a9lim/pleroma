# ogham — the ogdoad expression language

Status: **DRAFT v0.1** (2026-06-11, design by a9lim + Claude; inventory pass run
against this checkout). This document is the implementation contract: every
decision below either cashes out as a vector in [`spec/conformance.txt`](conformance.txt)
or it is not really decided. Implementing agents work until the corpus is green;
judgment calls not covered here go back to the spec, not into the code.

ogham is a small calculator language over the ogdoad core: one world (scalar
backend + Clifford metric) per session, expressions over that world's algebra,
bindings, and nothing else. No control flow, no user functions, no floats.
File extension `.og`. The name: og(doad) + the ancient stroke-script — fitting
a language whose operators are strokes and ticks (`*`, `↑`, `∧`, `⋅`, `/`).

---

## 1. Design principles

1. **Weird numbers first.** Scalar literals are the richest part of the
   grammar. `*` belongs to nimbers, not to multiplication.
2. **parse ∘ display = id.** Display (v2, §9) emits canonical ogham; the
   parser's input language is a superset of Display's output language.
3. **Two layers: canonical and sugar.** Canonical uses the unicode math glyphs
   where ASCII is contested (`ω ↑ ∧ ⋅`); ASCII stays canonical where it is
   uncontested (`* e + - / == [ ] ( )`). Sugar is input-only; the REPL echoes
   canonical (the REPL is the tutor).
4. **Unambiguous to the end.** No juxtaposition anywhere — not even as sugar.
   No inference of worlds from literals. No context-sensitive operators.
5. **One world per session** (core rule 5). Mixing is a parse/eval-time error,
   never a coercion.
6. **Errors are mathematical content.** Partiality (Kummer boundary, monomial
   inverses, non-fields) surfaces as typed errors with the math in the message.
7. **Pure Rust, zero deps, no pyo3 outside `src/py/`** (core rule 1). The
   parser/evaluator is a new `src/ogham/` area; the Python `eval` hook lives in
   `src/py/` behind the `python` feature.

## 2. Symbols and codepoints

| meaning | canonical | codepoint | ASCII sugar | notes |
|---|---|---|---|---|
| omega | `ω` | U+03C9 | `w` | atom; also inside star-literals |
| power | `↑` | U+2191 | `^` | right-assoc; Knuth's arrow |
| wedge | `∧` | U+2227 | `&` | exterior product |
| product | `⋅` | U+22C5 | `.` | the algebra's product; U+00B7 `·` also accepted on input |
| nimber prefix | `*` | — | — | value marker in nim-worlds (§6.3) |
| blade prefix | `e` | — | — | `e0`, `e1`, … basis 1-blades |
| neg / sub | `-` | — | — | unary and binary |
| recip / div | `/` | — | — | unary and binary (§7.4) |
| add | `+` | — | — | |
| equality | `==` | — | — | top-level only, prints `true`/`false` |
| binding | `=` | — | — | `name = expr` |
| vector | `[a,b,c]` | — | — | `Σ aᵢ⋅eᵢ`; length must equal world dim |
| comment | `#` | — | — | to end of line |

Reserved for the future, must lex but reject with `E_Reserved`: `↑↑` (towers),
`{` `}` `|` (game forms `{L|R}`, contractions), `O(` (precision tails),
`!` (factorial — genuinely unassigned), `t` (poly/ratfunc variable, §6.8).

**Unary-fill principle**: a unary form of a binary operator fills the left
operand with the operator's identity. `-a = 0 - a`, `/a = 1/a`. Only the two
inverse-taking operators have non-trivial unary forms; no other operator gets
one.

## 3. Lexical structure

- Tokens are self-delimiting; there are **zero juxtaposition / maximal-munch
  rules**. Whitespace separates tokens but is never semantic.
- `INT`: `[0-9]+`, value must fit `u128`. No sign (sign is unary `-`); the one
  exception is a tight signed exponent immediately after `↑` (§5).
- `IDENT`: `[a-z][a-z0-9_]*`, excluding reserved words. Reserved everywhere:
  `w`, `true`, `false`, stdlib function names (§8). Reserved per-world: `x` in
  `f4…f27` worlds (the field generator), `t` in future poly/ratfunc worlds.
- `e` followed immediately by digits lexes as a BLADE token (`e0`, `e12`).
  `e` alone is an error (not an identifier).
- `*` followed by anything lexes as the STAR prefix token; `*` is never an
  infix operator.
- Sugar substitution happens in the lexer: `w→ω`, `^→↑`, `&→∧`, `.→⋅`, `·→⋅`.
  After the lexer, only canonical tokens exist.

## 4. Grammar (EBNF)

Statements (one per line; blank lines and comment-only lines are no-ops):

```ebnf
statement   = binding | expression ;
binding     = IDENT "=" expression ;          (* rebinding allowed *)

expression  = additive [ "==" additive ] ;    (* == not nestable *)
additive    = mulexpr { ("+" | "-") mulexpr } ;
mulexpr     = wedge   { ("⋅" | "/") wedge } ;
wedge       = unary   { "∧" unary } ;
unary       = { "-" | "/" } power ;
power       = atom [ "↑" exponent ] ;         (* right-assoc via recursion *)
exponent    = [ "-" ] INT
            | "(" expression ")" ;            (* Index sort; Scalar iff base is ω in surreal-family worlds *)
atom        = INT | starlit | "ω" | BLADE | vector | call
            | IDENT | "(" expression ")" ;
vector      = "[" expression { "," expression } "]" ;
call        = IDENT "(" [ arglist ] ")" ;
arglist     = arg { "," arg } ;
arg         = expression ;                    (* sort checked per signature *)

starlit     = "*" ( INT | "ω" | "(" cnf ")" ) ;
cnf         = cnfterm { "+" cnfterm } ;       (* strictly descending exponents, else E_CnfOrder *)
cnfterm     = INT
            | "ω" [ "↑" cnfexp ] [ "⋅" INT ] ;
cnfexp      = INT | "ω" | "(" cnf ")" ;
```

Notes:

- **Star-literals are structural, not arithmetic.** Inside `*(…)` the symbols
  `+ ⋅ ↑` build a CNF ordinal *index* (the nimber's address in On₂), they do
  not evaluate. Exponents strictly descend or `E_CnfOrder`. Outside the star,
  all arithmetic is nim arithmetic. `*(ω + 1)` is the nimber at ordinal ω+1;
  `*ω + *1` is a nim-sum that happens to equal it.
- Unparenthesized star applies only to `INT` and bare `ω`: `*5`, `*ω`.
  Everything else takes parens: `*(ω↑2)`, `*(ω⋅3 + 5)`. So `(*ω)↑2` (nim
  square) and `*(ω↑2)` (index ω²) are visibly different, settling the binding
  question: **the star binds tighter than `↑`**, i.e. `*ω↑2 = (*ω)↑2`.
- The surreal-family worlds (`surreal`, `omnific`) allow CNF **at expression
  level, unstarred and live**: `3⋅ω↑2 - ω + 5` is ordinary arithmetic over
  monomials. `ω↑e` with non-integer `e` is the Hahn monomial constructor and
  requires base exactly `ω` (§7.3).

## 5. Precedence (tight → loose)

```text
atoms:  INT, *‹i›, ω, e‹i›, [a,b,c], f(...), (...)
↑           power, right-assoc (2↑3↑2 = 2↑9); tight signed INT exponent ok (ω↑-1)
unary - /   neg, reciprocal
∧           wedge
⋅  /        product, right-division, left-assoc
+  -        add, subtract
==          equality (non-associative, top level only)
```

Wedge tighter than `⋅` follows Hestenes (outer binds tighter than geometric).
Check: `*3⋅e0∧e1` = `*3 ⋅ (e0∧e1)`. Display v2 relies on this: blade terms
print unparenthesized.

**Host-language caveat** (§13): Rust and Python cannot reproduce this table
for the overloaded operators (`&` binds looser than `+` in Python). The
precedence above is ogham's, full stop; host code uses parens.

## 6. Worlds

A session holds exactly one world: a scalar backend monomorphised into a
`CliffordAlgebra<S>` plus environment. Declared by colon-command (REPL) or a
leading directive line (`.og` files use the same syntax without the colon
prompt):

```text
:world ‹name› ‹dim› q=[s0,…,s(n-1)] [b=[(i,j):s, …]] [a=[(i,j):s, …]]
:world ‹name› ‹dim› grassmann
:world nimber gold(m,a)            # dim = m, metric = forms::trace_form::gold_form(m,a)
:world ‹name› 0                    # pure scalar work, no metric
```

`q`/`b`/`a` mirror `Metric::diagonal` / `::new` / `::general`
(src/clifford/engine/metric.rs): `q` dense length-n, `b`/`a` sparse `i<j`
pairs. Values are scalar literals of the world. Declaring `a≠∅` prints a
warning that `rev`, `dual`, and everything built on reverse is unavailable
(the engine panics there; ogham refuses earlier with `E_GeneralMetric`).
`dim ≤ 128` (`MAX_BASIS_DIM`).

**Typing.** Two value sorts only. **Element**: every value in a world is a
`Multivector<S>`; pure scalars are grade-0 elements (in `dim 0` worlds,
everything is grade-0). **Index**: meta-integers (`i128`) used for exponents,
grades, blade indices, stdlib integer args; Index expressions allow
`+ - ⋅ ↑` and parens, nothing else. Position determines sort; there are no
coercions between sorts.

### 6.1 v1 world menu (fixed dispatch table)

Const-generic backends require a compiled-in menu; v1 ships:

| world name(s) | backend | field? | notes |
|---|---|---|---|
| `nimber` | `Nimber` (u128) | yes | F_{2^128} |
| `ordinal` | `Ordinal` | partial | Kummer-checked (§7.5) |
| `surreal` | `Surreal` | partial | monomial inverses only |
| `omnific` | `Omnific` | no (units ±1) | |
| `integer` | `Integer` (i128) | no (units ±1) | |
| `fp2 fp3 fp5 fp7` | `Fp<2|3|5|7>` | yes | |
| `f4 f8 f16` | `Fpn<2,2|3|4>` | yes | char-2 extension fields |
| `f9 f27` | `Fpn<3,2|3>` | yes | |
| `f25` | `Fpn<5,2>` | yes | |

(The six `f*` names match the Python binding classes `F4…F27`,
src/py/scalars.rs. Extending the menu = adding one arm to the dispatch enum.)

Deferred to v1.1: `poly(p)`, `ratfunc(p)` (need the `t`-variable literal
grammar and ride on the Display-v2 `x→t` fix); all precision worlds
(`Qp/Qq/Laurent/Ramified/Gauss/Adele` — `O(p^k)` literal design is its own
iteration); games mode (`{L|R}`).

### 6.2 Integer literals per world (the `from_int` trap)

`Scalar::from_int` is the ℤ-ring map — in char-2 backends `from_int(3) = 1`.
Literal meaning is therefore defined per world and **never** via `from_int`
in nim-worlds:

| world | bare `INT` at Element position |
|---|---|
| `nimber`, `ordinal` | **error `E_BareInt`**, hint: `did you mean *3?` |
| `surreal`, `omnific`, `integer` | exact integer (`from_int`, overridden exactly there) |
| `fp*`, `f*` | residue (`from_u128`-style reduction; `f*` worlds: degree-0 constant) |

Bare `INT` at Index position is always a meta-integer, in every world.

### 6.3 Star-literals per world

- `nimber` world: `*n` with `n` a u128 — `Nimber(n)` (the representation
  constructor, src/scalar/finite_field/nimber/mod.rs). `*` alone is sugar for
  `*1` (CGT star); canonical prints `*1`.
- `ordinal` world: `*n`, `*ω`, `*(cnf)` — assembled from `Ordinal::from_u128`
  / `::monomial` / `::omega_pow` per the structural CNF. The star is the value
  marker; there are no unstarred Element literals in this world.
- All other worlds: `E_WrongWorld`, hint names the world that wanted it.

### 6.4–6.8 Other literal forms

- `ω` (atom): `surreal`/`omnific` worlds — `Surreal::omega()`. In `ordinal`
  world bare `ω` is `E_BareOrdinal` (hint: `*ω`); the glyph appears there only
  inside star-literals.
- Dyadic/rational values are spelled with division: `1/2`, `3/2` (the field
  operation *is* the literal syntax; in non-field worlds it errors honestly).
- `f*` worlds: the generator is the reserved identifier `x`
  (`Fpn::generator()`); elements are reached arithmetically (`x↑2 + x + 1`).
- `e‹digits›` blades: `alg.e(i)`, `E_BladeIndex` if `i ≥ dim`.
- Future `poly`/`ratfunc`: reserved `t`; fractions print as `(num)/(den)` —
  the current `[num] / [den]` display collides with vector syntax and is fixed
  by Display v2 (§9).

## 7. Semantics (desugaring to the engine)

All file:line references are to this checkout.

| ogham | engine call |
|---|---|
| `a + b` | `Multivector::add` (multivector.rs:85) |
| `a - b` | `Multivector::sub` (:109) — scalar `neg()` underneath, never literal −1 (core rule 3) |
| `-a` | `Multivector::neg` (:95) |
| `a ⋅ b` | `alg.mul(&a, &b)` (algebra.rs:141) |
| `a ∧ b` | `alg.wedge(&a, &b)` (algebra.rs:153) |
| `a / b` | `a ⋅ inv(b)` — **right division**; noncommutative worlds beware, documented not hidden |
| `/a` | grade-0: `Scalar::inv` else `alg.multivector_inverse(&a)` (inverse.rs:9); `None → E_NotInvertible` |
| `a ↑ k` (k ≥ 0) | iterated `alg.mul`, left fold; `a↑0 = 1` |
| `a ↑ -k` | `(/a) ↑ k` |
| `ω ↑ s` (surreal world, s an Element) | `Surreal::omega_pow(s)` — Hahn monomial constructor; any other base with Element exponent is `E_ExpSort` |
| `[a0,…,a(n-1)]` | `Σ alg.scalar_mul(&ai, &alg.e(i))`; length ≠ dim → `E_DimMismatch` |
| `a == b` | `PartialEq`, prints `true`/`false` |

Evaluation is strict, left-to-right; bindings live in a per-world environment
(cleared on `:world`). A bare expression statement evaluates and prints the
value's canonical display. If the *input* was not already canonical, the REPL
first echoes the canonical form of the parsed expression (the unparser, §10),
then the value.

### 7.5 Partiality (the honest edges)

| operation | behavior |
|---|---|
| `ordinal` mul/inv escaping the verified Kummer tower | `Ordinal::nim_mul`/`checked_inv` return `None` → `E_KummerEscape` ("beyond the source-verified tower below ω^(ω^ω)"). ogham never calls the panicking `Scalar::mul` path on Ordinal. |
| `surreal` inverse of a non-monomial | `Surreal::inv = None` → `E_NotInvertible` ("only CNF monomials invert exactly; 1/(ω+1) is an infinite Hahn series") |
| `integer`/`omnific` inverse of non-units | `E_NotInvertible` |
| `/0` anywhere | `E_DivisionByZero` |
| grassmann/degenerate inverses | `multivector_inverse → None → E_NotInvertible` |

## 8. Stdlib v1

Six functions, all thin wrappers; signatures are sorted (E = Element,
I = Index):

| call | engine | notes |
|---|---|---|
| `rev(E)` | `alg.reverse` (algebra.rs:173) | `E_GeneralMetric` if `a ≠ ∅` |
| `grade(E, I)` | `alg.grade_part` (:193) | |
| `even(E)` | `alg.even_part` (versor.rs:13) | |
| `dual(E)` | `alg.dual` (versor.rs:183) | `None → E_NotInvertible` (pseudoscalar) |
| `tr(E, I)` | `nim_trace(x, m)` (artin_schreier.rs:11) | nimber world, grade-0 arg; m a power of 2 ≤ 128; `f*` worlds: `FieldExtension::trace` (extension.rs:60), 1-arg form `tr(E)` |
| `frob(E)` | `FiniteField::frobenius` (finite_field/mod.rs:48) | finite-field worlds, grade-0 arg |

Everything else (versors, sandwiches, contractions, meet, spinor norms) is
deliberately out of v1 — reach those from Rust/Python. The Gold chain works
day one: `tr(x ⋅ x↑(2↑a), m)`.

## 9. Display v2 (canonical form) — normative delta

**Shipped 2026-06-11** (ledger: `roadmap/DONE.md` → `ogham-foundations`); the
"v1" column below records the pre-v2 state for the historical delta.

Display v2 makes every `Display` impl in language scope emit canonical ogham.
One rendering path each (the Multivector path is already unified,
multivector.rs:59-83).

| type | v1 (current, cited) | v2 (canonical) |
|---|---|---|
| `Nimber` | `*5` (nimber/mod.rs:36) | unchanged |
| `Ordinal` | `5`, `ω`, `ω·3`, `ω^2`, `ω + 1`, `ω^(ω)` (ordinal/mod.rs:217) | **star-wrapped**: `*5`, `*ω`, `*(ω⋅3)`, `*(ω↑2)`, `*(ω + 1)`, `*(ω↑(ω))`; `·`(U+00B7)→`⋅`(U+22C5) |
| `Surreal` | `3ω^2 - ω + 5`, `ω^-1`, `ω^1/2`, `ω^(ω)` (surreal/mod.rs:247) | `3⋅ω↑2 - ω + 5`, `ω↑-1`, `ω↑(1/2)`, `ω↑(ω)` — explicit `⋅`, `↑`; exponent bare iff (signed) integer, else parenthesized |
| `Omnific` | delegates to Surreal | follows |
| `Integer`, `Fp` | plain int | unchanged |
| `Fpn` | `3x^2 + 2x + 1` (fpn.rs:549) | `3⋅x↑2 + 2⋅x + 1` (in a field where the coefficients exist, e.g. F_125; coefficients reduce mod p) |
| `Poly` | `1 + (2)·x` (poly.rs:29) | `1 + 2⋅t` — variable **x→t** (matches `F_q[t]`, frees `x` for Fpn); parens only when the coefficient renders non-atomically (more than one token) |
| `RationalFunction` | `[num] / [den]` (function_field.rs:115) | `(num)/(den)` — `[…]` belongs to vectors |
| `Multivector` | `e0e1`, `3*e0e1`, joined ` + ` (multivector.rs:59) | blades `e0∧e1`; coefficients `c⋅e0∧e1` with the same atomic-parens rule (`(x + 1)⋅e0∧e1`); elision unchanged (`1` elided; `-1` → `-e0∧e1`); **join rule**: if a term's rendering starts with `-`, strip it and join with ` - ` instead of ` + ` (string-level, char-agnostic — no sign predicate on `Scalar` exists or is wanted); **zero rule**: the empty multivector renders as `S::zero()`'s display (`*0` in nim-worlds, `0` elsewhere) — the current bare `0` would not round-trip where bare integers are `E_BareInt` |

**Atomicity (operational rule)** for the coefficient-parens decisions above: a
rendering is atomic iff it contains no spaces and no operator characters
(`⋅ ∧ ↑ /`, internal `+ -`) outside balanced parentheses; a single *leading*
`-` is a unary sign, not an operator — it attaches bare and is then lifted by
the join rule. So `42`, `-2`, `*5`, `*ω`, `x`, and `*(ω⋅7)` are atomic (the
star-literal is self-delimiting); `x + 1` and `ω↑-1` are not. Atomic coefficients attach bare (`*(ω⋅7)⋅e0∧e2` — never
double-wrapped); non-atomic ones get parens (`(x + 1)⋅e0∧e1`).

Unchanged and out of scope: `CliffordInvariants` names (`M_2(R)` …), tropical
display, game displays, error-message strings.

**Blast radius** (inventoried): assertions at ordinal/mod.rs:286-291 (5),
games/nimber_game.rs:221 (1); doc-comment examples (e.g. multivector.rs:48);
AGENTS.md style bullets ("blades render `e0e1`"; ordinal/surreal display
examples in the root and pillar AGENTS.md files); demo.py prose prints.
`tests/` has zero display assertions. The `*0` strings in py error messages
stay (Nimber display unchanged).

## 10. The unparser

A canonical pretty-printer over the AST, used for (a) the REPL echo of
non-canonical input and (b) conformance `~` vectors. Rules: canonical glyphs;
minimal parens per §5 precedence (re-parsing the output must yield the same
AST); spacing exactly as Display v2 emits: single spaces around `+ - ==` and
after `,`; `⋅ / ∧ ↑` and unary operators tight (`3⋅ω↑2 - ω + 5`,
`*ω⋅e0∧e1`).

## 11. Error taxonomy

Every error is `OghamError { kind, span, message, hint: Option<String> }`.
Kinds and canonical hints (conformance `!` vectors match on kind + message
substring):

| kind | trigger | hint example |
|---|---|---|
| `E_Parse` | token/grammar violation | |
| `E_Reserved` | `{ } \| ! ↑↑ O(` etc. | "reserved for future games/precision syntax" |
| `E_ExpSort` | non-integer exponent, e.g. `e0^e1` | "`↑`/`^` is power; the wedge product is `∧`/`&`" |
| `E_IndexSort` | Element where Index expected, and vice versa | |
| `E_BareInt` | bare integer at Element position in nim-worlds | "did you mean `*3`?" |
| `E_BareOrdinal` | bare `ω` in ordinal world | "values are starred here: `*ω`" |
| `E_WrongWorld` | literal form foreign to the session world | "`*3` is a nimber; this is the `surreal` world" |
| `E_CnfOrder` | star-literal exponents not strictly descending | "CNF indices are structural: write `*(ω + 1)`, not `*(1 + ω)`" |
| `E_KummerEscape` | ordinal mul/inv past the verified tower | "below ω^(ω^ω), primes ≤ 47 — see OPEN.md" |
| `E_NotInvertible` | failed `inv`/`multivector_inverse`/`dual` | per-world math in message (§7.5) |
| `E_DivisionByZero` | `/0` | |
| `E_BladeIndex` | `e‹i›` with i ≥ dim | |
| `E_DimMismatch` | vector length ≠ dim; vector in dim-0 world | |
| `E_GeneralMetric` | `rev`/`dual` with `a ≠ ∅` | "reverse is undefined for the Chevalley construction" |
| `E_Unbound` | unknown identifier | |
| `E_Arity`, `E_UnknownFn` | call errors | |

## 12. REPL

`examples/ogham_repl.rs` (the binary driver; the library lives in
`src/ogham/`). The REPL layer owns a dispatch enum over the §6.1 menu — one
arm per monomorphised `CliffordAlgebra<S>` — which is exactly how rule 5 is
preserved. Colon-commands (REPL only, not in the grammar): `:world …` (§6),
`:env` (bindings + world summary), `:help`, `:quit`. Echo behavior per §7.
Invariant queries (`:arf` etc.) deferred — the colon-command namespace is
where they will land, not the function namespace.

## 13. Host operator alignment (Rust + Python)

**Shipped 2026-06-11** alongside Display v2 (same ledger entry). The
overloads speak the same dialect as the display. Pre-alignment state, for the
record: `Multivector` had `^` = `BitXor` = wedge; scalars had `+ - *` via
`impl_scalar_ops!`; `Ordinal` deliberately has no owned `*` (unchanged).

| op | Rust | Python |
|---|---|---|
| wedge | `impl BitAnd for Multivector` (`a & b`); **remove** `BitXor`-as-wedge | `__and__`; `__xor__` raises `TypeError` with the §11 `E_ExpSort` hint during a deprecation window |
| power | scalars: `impl BitXor<u128>` for total-product backends (`x ^ 3`, square-and-multiply via `mul`); RHS is the meta-integer type, so no clash with any element-element op. **Multivectors get no power operator** — the geometric product needs the metric, so power ships as `CliffordAlgebra::pow(&self, v: &Multivector<S>, k: u128)` (ogham's `↑` desugars to it; negative exponents stay in the evaluator via `multivector_inverse`) | **`**` (`__pow__`)** — Python has a native power operator and the Py multivector classes are algebra-bound, so this works where the Rust operator can't; do not bend `__xor__` into power |
| ordinal power | **no operator.** `Ordinal` omits owned `*` because nim-mul is partial; an `^` that panics through iterated partial mul would contradict that deliberate omission. Add `Ordinal::nim_pow(&self, k: u128) -> Option<Ordinal>` beside `nim_mul` instead | same: `pow()` method returning/raising honestly |
| product | `*` stays `Mul` (Rust has no native power operator to displace it; `⋅` isn't typeable as an operator) | `*` stays |

Two flags, decided here:

- **Nimber `^` danger**: Rust users may expect `Nimber ^ Nimber` = XOR =
  nim-*addition*. The power overload takes `u128` on the right, so
  `Nimber ^ Nimber` simply does not compile — the type system is the
  disambiguation. Never implement element-element `BitXor` on any backend.
- **Precedence mismatch is documented, not fixed**: Python's `&` binds looser
  than `+` (so `a + b & c` ≠ ogham's reading) and Rust's `^` looser than `*`.
  Host code parenthesizes; rustdoc/docstrings on the overloads say so.

## 14. Conformance corpus

`spec/conformance.txt`, UTF-8, line-based:

```text
@world ‹world-decl args, exactly as after ":world"›   # resets bindings
> ‹input line›            # statement, exactly as typed (may use sugar)
~ ‹canonical unparse›     # optional: expected canonical echo of the input
= ‹expected display›      # value line; or:
! ‹E_Kind›: ‹message substring›
```

Blocks separated by blank lines; `@world` persists until the next `@world`.
The harness is `tests/ogham_conformance.rs` (pure Rust, reads the file,
no_python), mirrored by a pytest that drives the Python `eval` hook over the
same file — one corpus, two frontends. The corpus ships with hand-verified
vectors (small nim arithmetic, char-2 wedges, dyadic surreals, Conway's
`(*ω)↑3 = *2`); after WP3 lands, the harness gains a `--bless` mode to
*extend* (never overwrite) the corpus with engine-generated vectors — the
engine is the value-oracle for values, the spec is the oracle for syntax,
sorts, and errors.

## 15. Work packages

Status: **WP1 and WP7 shipped 2026-06-11** (sequentially, live tree — they
share multivector.rs; ledger: `roadmap/DONE.md` → `ogham-foundations`).
Remaining sequencing: WP2 → WP3 → (WP4 ∥ WP5 ∥ WP6). Every agent gets an
explicit `model:` pin. Acceptance for all: `cargo test`, `cargo clippy
--all-targets`, cold `cargo doc --no-deps` warning-clean; WP6/WP7 add
`cargo check --features python` + `clippy --features python --all-targets`;
WP1 adds a `demo.py` rerun (display changes don't surface in `cargo test`).

| WP | scope | model |
|---|---|---|
| **WP1 Display v2** | §9 exactly: ordinal star-wrap + `⋅` codepoint, surreal/Fpn explicit ops, Poly `x→t` + paren rule, RationalFunction `(…)/(…)`, Multivector blade/coefficient/join rules. Fix the §9 blast-radius assertions and doc comments; update AGENTS.md display bullets (root + scalar + clifford). | sonnet |
| **WP2 Lexer / parser / AST / unparser** | `src/ogham/{lex,ast,parse,unparse}.rs`, pure Rust, zero deps, world-independent (literal *forms* parse everywhere; world legality is WP3's). §3–§5, §10. Unit tests: golden token streams, precedence cases from §5, unparse∘parse = id on the corpus's `~` lines. | sonnet |
| **WP3 Worlds + evaluator** | `src/ogham/{world,eval,error}.rs`: the §6.1 dispatch enum, per-world literal mapping (§6.2–6.8), §7 desugaring, §7.5 partiality, §8 stdlib, §11 errors. The judgment-heavy package. | opus |
| **WP4 REPL** | `examples/ogham_repl.rs` + colon commands (§12). | sonnet |
| **WP5 Conformance harness** | `tests/ogham_conformance.rs` + corpus format parser + `--bless` extension mode (§14). | sonnet |
| **WP6 Python eval** | `ogham_eval(world: &str, src: &str)` pyfunction + per-class `__and__`/`__pow__`/`__xor__` alignment (§13); pytest mirror of the corpus. | sonnet |
| **WP7 Host operators (Rust)** | §13: `BitAnd` wedge, remove `BitXor`-as-wedge, `BitXor<u128>` power on Multivector + total scalar backends, `Ordinal::nim_pow`, rustdoc precedence caveats. Migrate the ~3 in-repo uses of `^`-as-wedge. | sonnet |
