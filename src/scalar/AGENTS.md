# AGENTS.md — `src/scalar/`

The PILLAR of commutative coefficient worlds: the `Scalar` trait and every
concrete backend the Clifford engine and forms layer run over. Pure Rust,
generic; the per-backend Python wrappers live in `src/py/scalars.rs`.

## Two orthogonal organizing axes

The directory is grouped **by place** (the "any number" table: each field beside
its ring of integers). A *second* axis — the **characteristic trichotomy** (char
0 / odd / 2) — cuts ACROSS the table and organizes `forms/` instead. Hold both:
the place table says *where a number lives* (Archimedean, p-adic, finite,
transfinite); the trichotomy says *which classification theory applies*.

```
              FIELD                     RING OF INTEGERS
 Archimedean  Rational (ℚ)              Integer (ℤ)               exact/
 transfinite  Surreal (No)              Omnific (Oz)              big/
 p-adic       Qp, Qq                    Zp, WittVec              small/, finite_field/
 fn field     RationalFunction (F_q(t)) Poly (F_q[t])            global/, poly.rs
 finite       Fp, Fpn, Nimber           —                        finite_field/
```

The trait layer makes the table structural:

- `integrality.rs` — the (field, ring-of-integers) pairing (`HasFractionField` /
  `HasRingOfIntegers`).
- `valued.rs` — the local-field view (`Valued`: valuation + uniformizer ϖ).
- `residue.rs` — the residue field column `k = 𝒪/𝔪` (`ResidueField`).
- `analytic.rs` — root-taking.
- `extension.rs` — finite separable extensions / relative trace-norm
  (`FieldExtension`).
- `exactness.rs` — the exact-vs-capped-precision boundary.

Fixed-width mathematical payloads here are `u128`/`i128`: finite-field values,
field orders, p-adic moduli/valuations, nimber data, ordinal coefficients, and
exact rational/integer numerators. `usize` is only for array lengths, dimensions,
and const-generic sizes that are inherently indices.

- **`poly.rs`** — `Poly<S>`, the shared dense-univariate polynomial ring `S[t]`
  (low-degree-first, trimmed). The crate's one polynomial primitive: `Gauss` and
  `RationalFunction` store `num/den` as `Poly` pairs, and the function-field place
  layer (`forms/function_field.rs`) uses its `divrem`/`gcd`/`pow_mod`. As `S[t]` it
  is the **ring of integers** of `S(t)`, so it impls `Scalar` + `HasFractionField`
  (Frac = `RationalFunction<S>`); its units are the nonzero constants, so `inv` is
  partial. Display is canonical ogham (Display v2, `docs/ogham/ogham.md` §9): variable
  `t`, explicit `⋅`, coefficient parens only when non-atomic — and it owns the
  shared `pub(crate)` `atomic`/`attach_coeff` helpers the `Multivector` display
  also uses (atomic = no spaces and no `⋅ ∧ ↑ / + -` outside balanced parens; a
  single leading `-` is a unary sign, carried bare).
- **`newton.rs`** — `NewtonPolygon`: the lower convex hull of `{(i, v(aᵢ))}` for
  `f = Σ aᵢtⁱ` over a `Valued` field (`of`/`vertices`/`degree`/`slopes`/
  `root_valuations`/`zero_root_multiplicity`). The tropicalization of the Springer
  residue filtration — the Newton slope theorem (root valuations = negated hull
  slopes) is the `Valued`-side oracle, tested over `Qp`/`Laurent`/`Ramified`.

## The `Scalar` trait + the trait layer

- **`mod.rs`** — the `Scalar` trait (`add`/`neg`/`mul`/`zero`/`one`/`is_zero`/
  `inv`/`characteristic`/`from_int`) + `Display` as a supertrait + the "any number"
  table doc + the flat re-export hub. Also `impl_scalar_ops!`: total-product backends
  get concrete-type operators (`+ - *` and unary `-`) forwarding to the trait methods.
  `Ordinal` is the deliberate exception: additive operators only, multiplication behind
  the checked `nim_mul` API (the represented Kummer tower has an honest boundary). `/`
  stays a method (inv is partial). **The operators are NOT a `Scalar` supertrait** —
  see "things that look like bugs". Ogham backend helpers live here too:
  `checked_factorial_i128` (host-carrier roof: `33!`, not `34!`) and
  `factorial_in_scalar<S>` (in-world product via `from_int`, zeroing at positive
  characteristic).

  **`Scalar::from_int(n: i128) -> Self`** is the ℤ-embedding (unique unital ring
  homomorphism ℤ → R). The default double-and-add over `one()`/`neg()` is correct for
  every characteristic: char-2 worlds automatically get `n mod 2` because `1+1=0`.
  Backends with a direct construction (`Rational`, `Integer`, `Fp`, `Fpn`, `Zp`,
  `Qp`, `Qq`, `WittVec`, `Surreal`, `Omnific`) override for efficiency. Do NOT
  override in char-2 worlds with a bit-cast constructor — `Nimber(n)` and
  `Ordinal::from_u128(n)` are *representation* constructors (which nimber / which
  ordinal), not ℤ-embeddings.

  **`Display` is a `Scalar` supertrait** (`Scalar: Clone + PartialEq + Debug +
  Display`). Every backend has an `impl fmt::Display`, with `Debug` delegating to
  `Display` (byte-identical output). The legacy ℤ-embedding spellings
  (`Rational::int`, `Fp::new`, `Zp::new`, `Qp::from_i128`) were retired in the
  taste sweep — `from_int` is the one spelling, trait-level with per-backend
  overrides; representation constructors with different semantics keep their
  own names (see above).
- **`integrality.rs`** — `HasFractionField {Frac; to_fraction}` +
  `HasRingOfIntegers {Int; is_integral/to_integer}`, with `Int:
  HasFractionField<Frac=Self>` tying the loop. Impl'd for the **five** distinct-type
  rows (ℤ⊂ℚ, Oz⊂No, Zp⊂Qp, W_N⊂Qq, F_q[t]⊂F_q(t)) PLUS the blanket Surcomplex
  transport (ℤ[i]⊂ℚ[i] falls out). Laurent/Ramified `F_q[[t]]`/`O[π]` are same-type
  valuation subrings, so they stay out (`is_integral` only) — honest.
- **`valued.rs`** — the `Valued` trait: a discrete valuation + canonical uniformizer
  ϖ, impl'd for the local fields/functors (Qp/Qq/Laurent/Gauss, and Ramified where ϖ
  is the adjoined π). The spine of the "local fields" view (cuts across `small/` +
  `functor/`). NOT a `Scalar` supertrait (rings of integers + exact Archimedean
  worlds excluded).
- **`analytic.rs`** — the analytic layer, split on where precision lives.
  `ExactRoots {is_square; sqrt}` (no precision arg) for Rational, Nimber, Fp, Fpn,
  Surreal (exact via the fixed-point bridge over the lazy roots), Laurent, and the
  blanket `Surcomplex<R: ExactRoots+Ordered>` (the algebraic-closure √(a+bi)).
  `SeriesRoots {sqrt_to_terms; nth_root_to_terms; inv_to_terms}` (caller-chosen n)
  is the lazy interface — Surreal-only (the one world with unbounded, not type-fixed,
  precision). `Ordered {sign}` is the branch-picking datum the Surcomplex blanket
  needs. The residue Tonelli roots (`fp_sqrt`/`fq_sqrt`) live here (shared with
  `small/analytic`'s Hensel seed). P-adics use checked inherent APIs (dyadic
  squarehood/root construction can be unknown); Gauss/Ramified are excluded
  honestly. NOT a `Scalar` supertrait, like `Valued`.
- **`residue.rs`** — the `ResidueField: Valued` trait: `k = 𝒪/𝔪` (assoc type
  `Residue`) + two reductions — `residue` (canonical `𝒪 → k`, `None` below the
  integers) and `residue_unit` (the **angular component** `ac(x) ∈ k*`, `None` only
  for 0) + the multiplicative Teichmüller section `τ : k → 𝒪`. Impl'd for the local
  fields/functors (Qp→Fp, Qq→Fq, Laurent→S, Ramified→base residue, Gauss→k(tbar)).
  It is what lets `forms/springer/local.rs` write the discrete Springer
  decomposition once. Globals (Adele/RationalFunction) stay out — per-place residues
  live at the forms layer. NOT a `Scalar` supertrait.
- **`exactness.rs`** — marker traits for the representation contract: `ExactScalar`,
  `ExactFieldScalar`, `PrecisionScalar`. Exact finite/truncated rings (`Zp`,
  `WittVec`) are exact as represented; capped-relative models
  (`Qp`/`Qq`/`Laurent`/`Ramified`/`Gauss`/`Adele`) are marked separately and stay
  out of the exact-ring fuzz. NOT `Scalar` supertraits — generic code opts in only
  when it needs the promise.
- **`extension.rs`** — the `FieldExtension: Scalar` trait: a finite separable
  extension `E/F` with `extension_degree`/`embed`/`trace`/`norm` to a distinguished
  `Base`. The orthogonal view of `FiniteField`'s relative trace/norm (one fixed base
  vs. any subfield). Impl'd for:
  - `Surcomplex<S: Ordered>` (deg 2);
  - `Fpn<P,N>` over `Fp<P>` (deg N, delegating to the tested `FiniteField` relative
    trace/norm);
  - `Qq<P,N,F>` over `Qq<P,N,1>`=`Q_p` (deg F, via the Witt Frobenius; base kept in
    the Qq family to dodge the Qp-`u128`/Qq-`usize` const-kind clash);
  - `Nimber` (`Base = Fp<2>`, degree 128, F_{2^128}/F_2): trace = the absolute
    `nim_trace(·,128)`; norm onto F_2*={1} is trivial (1 for nonzero).

  Ramified (non-Galois, degenerate trace form) and Gauss (transcendental, infinite
  degree) excluded honestly — the same boundary `analytic.rs` draws.
  **`CyclicGaloisExtension: FieldExtension`** adds the Galois data the twisted trace
  form needs: `basis()` (an F-basis) + `sigma()`/`sigma_power(k)` (the cyclic
  generator). Impl'd for `Surcomplex` (σ=conj, basis {1,i}), `Fpn` (σ=Frobenius,
  coordinate basis), `Qq` (σ=Witt-Frobenius, Teichmuller-lifted residue basis), and
  `Nimber` (σ=nim-Frobenius, bit basis {1,2,…,2¹²⁷}). The relative trace stays
  `FieldExtension::trace` — σ/basis are the only new data. Consumed by
  `forms/trace_form.rs` and `clifford/frobenius.rs`. NOT a `Scalar` supertrait.
- **`tropical.rs`** — the `Semiring` trait + `Tropical<C>`, the `(min,+)/(max,+)`
  semiring. A SIBLING structure (like `Valued`), **not** a `Scalar`: a semiring
  drops the additive inverse (tropical `⊕` is idempotent), so it is not a ring and
  never enters `clifford/` — the same boundary the game *group* hits. The convention
  is a sealed compile-time marker (`TropicalConvention`: `MaxPlus`/`MinPlus`), so the
  two dual semirings are distinct, non-interoperating types sharing one impl body.
  The games payoff (thermography IS tropical) lives in
  `games/tropical_thermography.rs`; the laws are fuzzed in both conventions in
  `tests/tropical_axioms.rs`.

## `exact/` — the Archimedean char-0 base

- **`rational.rs`** — exact ℚ over i128, NOT a game backend: the char-0 scalar that
  validates the geometric product against the known Cl(p,q) classification before the
  exotic backends are trusted. (Overflow is a known limit; the surreal backend is the
  real char-0 home.) Implements the standard total-order traits by delegating to its
  inherent value `cmp`.
- **`integer.rs`** — exact ℤ, the coefficient ring for the exterior algebra of the
  game group (`games/game_exterior/`): games are a ℤ-module, not a ring, so Λ over
  ℤ is the structure that lives on all of game-world. Only ±1 invertible. Ogham's
  exact-division support is here: `divrem`/`rem` are Euclidean (`0 <= r < |b|`),
  and `div_exact` returns `IntegerDivExactError::Remainder(r)` on non-exact
  division. It also implements the standard total-order traits.

## `big/` — the transfinite worlds

- **`cnf.rs`** — `merge_descending`, the descending-CNF canonicalizer parameterized
  by the 3 places surreal & ordinal differ (exponent order: No value-order vs ordinal
  lex; coeff merge: + vs XOR; zero test). Deliberately a shared FUNCTION, not a
  `Cnf<C>` TYPE — the orders/algebras diverge (No is a field, On₂ isn't), so a shared
  type would be a false identity.
- **`surreal/`** — finite-support surreal Hahn/CNF backend (char 0):
  - `mod.rs` — CNF core: `Vec<(exponent: Surreal, coeff: Rational)>`, recursive
    exponents, Hahn arithmetic `ω^a·ω^b = ω^{a+b}`, Scalar, Debug, `truncate()`,
    standard total-order traits, and `rem` by a monic omega-power modulus
    (filter terms with exponent strictly below the modulus exponent).
  - `simplicity.rs` — the {L|R}/simplicity bridge (dyadic): `as_rational`/`as_dyadic`/
    `dyadic_birthday` + `simplest_above`/`_below`/`_between`, floor/frac (the Oz
    bridge).
  - `sign_expansion.rs` — exact `sign_expansion`/`from_sign_expansion` (dyadic,
    round-trips, length = birthday) + `as_ordinal`/`from_ordinal` + the transfinite
    (Gonshor) `SignExpansion` + `birthday_ordinal` + the transfinite inverse.
  - `analytic.rs` — the LAZY field layer (the `SeriesRoots` primitives):
    `inv_to_terms` (Neumann series) + `sqrt_to_terms`/`nth_root_to_terms` (real-closed
    roots to n terms; `Some` iff the leading coeff is a perfect ℚ-power).
- **`omnific.rs`** — the omnific integers Oz: `Omnific(Surreal)`, a transfinite
  commutative RING (not field). The surreal mirror of `Integer`; inherits the
  total order and monic-omega-power `rem` from the underlying surreal while
  revalidating the omnific-integral invariant.
- **`ordinal/`** — transfinite (ordinal) NIMBERS On₂, the char-2 mirror of surreal:
  - `mod.rs` — CNF core: `Ordinal = Vec<(exponent: Ordinal, coeff: u128)>`, the lex
    cmp, `as_finite`, `checked_inv`, `fuzzy` (`a != b` as the nimber game-value
    incomparability test), `Scalar`. The `Scalar::mul` route is
    panic-on-escape, matching the Kummer tower boundary; callers needing an explicit
    mathematical boundary use `nim_mul`.
  - `nim.rs` — char-2 NIM arithmetic: `nim_add` (coeff XOR) COMPLETE; `nim_mul`
    dispatches zero / finite×finite / the generator tower.
  - `subfield.rs` — finite-subfield detection for represented ordinal nimbers:
    minimal `F_{2^m}` by generator support plus Frobenius minimization, with
    common-degree helper for the forms façade. `checked_inv` uses this finite-field
    route beyond the old `F_64` window.
  - `tower.rs` — the prime-power generator tower (Conway/Lenstra/DiMuro): a monomial
    `ω^E` keyed by `place m ↦ base-p(m) digit vector`; `⊗` adds digit vectors and
    reduces with the Kummer carries `χ_u^u = α_u`. Non-scalar `α_u` (`α_7=ω+1`, …)
    branch a carry into a *sum*, recursed in by descending place. Carries are assembled
    from `ord_u(2)`, DiMuro's `Q(f(u))`, and the finite `m_u` rows from OEIS A380496
    (the b-file's 126 known rows, odd primes `3..=709`); a carry needing `m_719` (the
    first OEIS-unknown row)+ returns `None`, as does anything `≥ ω^(ω^ω)` (see
    `docs/OPEN.md`). The table extends reach, not feasibility: large primes are in the
    table but their `q_set`/finite-subfield reconstruction over the huge component field
    (`e_p` in the millions) is costly.
  - `cantor.rs` — ORDINARY (Cantor) `ord_add`/`ord_mul` (ω+ω=ω·2, 1+ω=ω) — the
    surreal birthday's run-length arithmetic. A distinct algebra, sharing only CNF.

The surreal↔ordinal **mirror** (No char 0 / On₂ char 2, sharing `cnf.rs`) is one of
the project's central symmetries.

## `small/` — the non-Archimedean (p-adic) local world

- **`qp.rs`** — `Qp<const P, const K>`: the p-adic FIELD Q_p. `p^val·unit`, char 0,
  inv total on nonzero. CAPPED-RELATIVE precision: mul/inv exact, addition NOT
  associative across precision boundaries (a precision model, like float). EXCLUDED
  from the exact-ring fuzz.
- **`zp.rs`** — `Zp<const P, const K>`: the p-adic integers Z_p (= Z/p^k), the ring
  of integers of Q_p. A LOCAL RING (p a non-unit), residue field F_p; Cl over it is
  non-semisimple.
- **`qq.rs`** — `Qq<const P, const N, const F>`: the UNRAMIFIED extension Q_q =
  Frac(W_N(F_q)), residue degree F (residue field F_q). To WittVec what Qp is to Zp;
  Qq with F=1 IS Qp.
- **`analytic.rs`** — the p-adic analytic layer over all four backends (mirror of
  `surreal/analytic`): checked `is_square() -> Option<bool>`, `sqrt() ->
  Option<Option<_>>`, and the Teichmüller rep τ. Odd `p` takes the Hensel/Newton
  path; dyadic inputs return `None` when the represented precision cannot decide
  squarehood or construct the root.

## `finite_field/` — the finite residue worlds

- **`mod.rs`** — the `FiniteField` TRAIT: the shared Galois engine (degree,
  conjugates, min_poly_monic, relative_trace/_norm, multiplicative_order, is_primitive,
  discrete_log) as default methods. An impl supplies only `frobenius`, integer `pow`,
  `ext_degree`, `group_order`, `group_order_factors`. nimber + fpn both impl it — one
  verified algorithm, two backends.
- **`fp.rs`** — `Fp<const P>`: the prime field F_P (any prime P — the odd-char
  comparison backend, and `F_2 = Base` for `Nimber`); the `Qp → Fp` residue field.
- **`fpn.rs`** — `Fpn<const P, const N>`: F_{p^N} via a (P,N)-keyed irreducible
  reduction poly (`reduction`, `reduction_kind` → the public `ReductionPolynomialKind`
  metadata, `is_supported_field`). Completes the odd-char tower AND the char-2
  odd-degree fields nimbers can't reach (F_8); supported `Fpn<2,N>` metrics classify
  through the char-2 Arf façade. (NB the static `field_order()` = field order p^N, ≠
  `multiplicative_order(&self)`.)
- **`nimber/`** — On₂ in u128 (= F_{2^128}), split by layer, re-exporting `nim_*`
  flat: `mod.rs` (wrapper + Scalar), `arithmetic.rs` (`nim_add`=XOR; `nim_mul` via
  Fermat-power recursion; `nim_square`/`nim_sqrt`/`nim_inv`), `artin_schreier.rs`
  (`nim_trace` + y²+y=c solver), `galois.rs` (impl FiniteField, with Pohlig–Hellman +
  BSGS overrides for `is_primitive`/`discrete_log`). `Nimber::fuzzy` is the
  game-value incomparability predicate: exactly `self != other`; do not turn that
  into `PartialOrd`.
- **`wittvec.rs`** — `WittVec<const P, const N, const F>`: Witt vectors W_N(F_q) as
  the truncated unramified ring (Z/p^N)[t]/(f̃). The char-p analogue of Z_p; its field
  of fractions is `small/qq.rs`.

## `functor/` — the root-level functors (ways to GROW a field)

Orthogonal to the place table: a 2×2 of (algebraic|transcendental) ×
(residue|value-extending), **all four corners filled**.

| | residue-extending | value-extending |
|---|---|---|
| **algebraic** | `surcomplex.rs` (root of x²+1) | `ramified.rs` (root of Eisenstein xᴱ−ϖ) |
| **transcendental** | `gauss.rs` (adjoin t as a unit, v(t)=0) | `laurent.rs` (adjoin t as uniformizer, v(t)=1) |

- **`surcomplex.rs`** — `Surcomplex<S>` = adjoin i over ANY backend (carries
  `conj()`). Only meaningful over char-0 worlds (over nimbers i²=1, degenerate).
- **`laurent.rs`** — `Laurent<S, const K>` = S((t)) to relative precision K. Over a
  finite field, the EQUAL-characteristic local cell F_q((t)) (the char-p mirror of
  Qp); ring of integers F_q[[t]] = the val≥0 subring. Capped-relative; EXCLUDED.
- **`ramified.rs`** — `Ramified<S: Valued, const E>` = adjoin a root of xᴱ−ϖ. The
  RAMIFIED local cell Q_p(p^{1/E}), the ramified twin of Qq. Always a field
  (Eisenstein), incl. wild/inseparable p|E. `Valued` with uniformizer π and
  `ResidueField` with the base residue field. EXCLUDED.
- **`gauss.rs`** — `Gauss<S: Valued>` = S(t) with the Gauss valuation (v(t)=0,
  transcendental residue ⇒ residue field k(t̄)). The last corner, Laurent's
  residue-extending twin. `Valued` itself, `ResidueField` as `k(tbar)`; precision
  model, EXCLUDED.

## `global/` — the adelic/global place

`Adele` is a finite-precision restricted-product model over ℚ, with `LocalQp` as the
runtime-prime p-adic cell. Useful for product-formula / Hilbert-reciprocity /
Hasse–Minkowski experiments in `forms/adelic.rs`; not an exact infinite-memory adele.
`LocalQp` (runtime prime, NOT const-generic) deliberately does NOT impl `Scalar` —
its world `(p,k)` is only known at construction — so it is the runtime-only analogue
of `forms`'s runtime `OddFiniteFieldForm`.

`RationalFunction<S>` (in `global/function_field.rs`) is the **equal-characteristic
mirror**: the global function field `F_q(t)`, the char-`p` analogue of `ℚ` as a
global field. Same field-of-fractions arithmetic as `Gauss` (over `Poly`, `inv =
den/num`, cross-mult equality) but a different ROLE — it carries *all* its place
valuations at once, so like `Adele` it is deliberately **not** `Valued`. Unlike the
precision-model functors it is **exact**, so it joins the `scalar_axioms` fuzz and
carries the `ExactScalar`/`ExactFieldScalar` markers. It feeds
`forms/function_field.rs`.

## Things that look like bugs but are not (scalar layer)

- **Scalar `+ - *` operators are concrete-only, NOT a `Scalar` supertrait.**
  Making `Scalar: Add+Sub+Mul+Neg` brings the ops into scope for every generic `S`,
  where `Mul::mul(self, Self)` shadows `Scalar::mul(&self, &Self)` at owned-receiver
  sites and forces clones the borrow-based engine avoids (70+ generic sites broke when
  tried). Don't promote them; don't migrate the engine's `.add()`/`.mul()` to
  operators. Don't add owned `*` back for `Ordinal`; `nim_mul` is the checked product
  boundary.
- **`ExactRoots`/`SeriesRoots`/`Ordered`/`Valued`/`ResidueField`/the exactness markers
  are NOT `Scalar` supertraits.** Not every world takes roots or has a valuation, so
  the bounds stay opt-in. The trait impls *delegate to inherent methods of the same
  name* (inherent-shadows-trait makes that delegate-not-recurse).
- **`Tropical` has no `neg`/`inv` and is deliberately not a `Scalar`.** A semiring's
  `⊕` is idempotent (`a ⊕ a = a`), so there is no additive inverse — that is the
  defining difference from a ring, and why `Semiring` is a sibling trait. `Tropical`
  never reaches `clifford/` (a Clifford algebra needs a commutative *ring*), exactly
  the boundary the game group hits. The two conventions are distinct types on purpose:
  `Tropical<MaxPlus>` and `Tropical<MinPlus>` do not interoperate, because
  thermography's two walls live in dual semirings.
- **`Surreal` has two square roots, by design.** `sqrt_to_terms(n)` is the lazy
  `SeriesRoots` primitive; `ExactRoots::sqrt(&self)` (0 args) is the exact value.
  Different arities, different precision contracts — don't unify them. (Python:
  `Surreal.sqrt(n)` lazy, `Surreal.exact_sqrt()` exact.)
- **P-adic square roots are checked inherent APIs, not `ExactRoots`.** `Zp`/`Qp`/
  `Qq`/`WittVec` use `Option<bool>` for squarehood and `Option<Option<_>>` for roots:
  outer `None` = unknown/unimplemented at the represented precision, inner `None` =
  definitely nonsquare. The finite fields and `Laurent` handle char 2 natively inside
  `ExactRoots`.
- **Surcomplex over nimbers is degenerate.** `i²=1`, `(1+i)²=0`, not a field.
  Surcomplex is only meaningful over char-0 worlds.
- **Surreal coefficients are ℚ, not ℝ.** The honest finite truncation of true CNF.
  Exponents *are* fully recursive surreals. `√2`, `√(2ω)` are honestly `None` (the
  leading coeff must be a perfect ℚ-power); `√ω = ω^{1/2}` IS exact (monomial).
- **`Surreal::inv` returns `None` for any non-monomial.** `1/(ω+1)` is an infinite
  Hahn series; finite support can't hold it.
- **`Surreal::birthday_ordinal`/`transfinite_sign_expansion` are `None` outside the
  representable subclass** (`√ω`, `ω−1`, `½ω`, mixed). Every *ordinal* (incl. ω^ω) is
  handled; `ε` is the one infinitesimal pinned. The honest Gonshor scope boundary.
- **`PartialEq for Surreal` is structural, not value-based.** The previous
  `self.cmp(other) == Ordering::Equal` was correct but allocated a subtraction.
  CNF uniqueness (Hahn series in reduced form — see the inline proof comment on the
  impl) guarantees structural equality and value equality coincide for all canonical
  surreals this module produces. A proptest in `tests/scalar_axioms.rs`
  (`surreal_structural_eq_matches_value_eq`) pins the agreement permanently.
- **`Qp` addition is not associative across precision boundaries.** Capped-relative
  (the standard p-adic model, like float). No finite-memory exact Q_p exists.
- **`nim_mul`'s `1u128 << (1u128 << n)` is not overflow-prone** for valid u128: bit
  positions < 128 ⇒ Fermat indices n ≤ 6, shift ≤ 64.
- **`Fpn::field_order()` is the field order `p^N` (static, no self); the element's
  multiplicative order is `multiplicative_order(&self)`.** Different things.
- **The `nim_*` Galois free fns delegate to the `FiniteField` trait; don't re-add
  inherent `Nimber` Galois methods.** An inherent `Nimber::degree` would shadow and
  recurse forever back through the free fn. To add a Galois op, add a default method
  to the trait (both nimber and fpn get it free). Nimber *overrides* `is_primitive`/
  `discrete_log` for the sharper large-field algorithms — intended, not duplication.
- **`scalar * multivector` works via the scalar's `__mul__` returning
  `NotImplemented`** so Python falls back to the MV's `__rmul__`. Don't make the
  scalar ops raise on a non-scalar operand.
- **`Poly<S>` has BOTH inherent methods and a `Scalar` impl with the same names**
  (`add`/`mul`/`neg`/`zero`/`one`/`is_zero`). Not duplication: the inherent methods
  are what `Gauss`/`RationalFunction`/the place layer call by value, and they SHADOW
  the trait at every concrete `Poly::…` site (so the trait bodies delegate, not
  recurse). The `Scalar` impl exists only so `Poly = F_q[t]` can be the
  `HasFractionField` ring of integers of `RationalFunction = F_q(t)`.

## Math facts worth not re-deriving

- nim-field: `F_{2^{2^k}}` = nimbers `< 2^{2^k}`. `F_n ⊗ F_n = (3/2)F_n` for a
  Fermat 2-power `F_n = 2^{2^n}`; distinct Fermat powers multiply ordinarily.
- Surreal CNF = finite-support Hahn series with ℚ coefficients; the ω-map is the
  monomial map and `ω^a·ω^b = ω^{a+b}` is a group homomorphism on represented
  monomials.
