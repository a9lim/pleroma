# AGENTS.md — `src/scalar/`

The PILLAR of commutative coefficient worlds: the `Scalar` trait and every
concrete backend the Clifford engine and forms layer run over. Pure Rust,
generic; the per-backend Python wrappers live in `src/py/scalars.rs`.

## Two orthogonal organizing axes

The directory is grouped **by place** (the "any number" table: each field beside
its ring of integers). A *second* axis — the **characteristic trichotomy** (char
0 / odd / 2) — cuts ACROSS the table and is what organizes `forms/` instead. Hold
both: the place table says *where a number lives* (Archimedean, p-adic, finite,
transfinite); the trichotomy says *which classification theory applies*.

```
              FIELD                     RING OF INTEGERS
 Archimedean  Rational (ℚ)              Integer (ℤ)               exact/
 transfinite  Surreal (No)              Omnific (Oz)              big/
 p-adic       Qp, Qq                    Zp, WittVec              small/, finite_field/
 fn field     RationalFunction (F_q(t)) Poly (F_q[t])            global/, poly.rs
 finite       Fp, Fpn, Nimber           —                        finite_field/
```

The (field, ring-of-integers) pairing is made **structural** in `integrality.rs`
(`HasFractionField` / `HasRingOfIntegers`); the local-field view is made
structural in `valued.rs` (`Valued`); the **residue field** column `k = 𝒪/𝔪` in
`residue.rs` (`ResidueField`); root-taking in `analytic.rs`; finite separable
**extensions** (relative trace/norm) in `extension.rs` (`FieldExtension`); and the
exact-vs-capped-precision boundary in `exactness.rs`.

Fixed-width mathematical payloads in this pillar are `u128`/`i128`: finite-field
values, field orders, p-adic moduli/valuations, nimber data, ordinal coefficients,
and exact rational/integer numerators all follow that contract. `usize` remains for
array lengths, dimensions, and const-generic sizes that are inherently indices.

- **`poly.rs`** — `Poly<S>`, the shared dense-univariate polynomial ring `S[t]`
  (low-degree-first, trimmed). The crate's one polynomial primitive: `Gauss` and
  `RationalFunction` store `num/den` as `Poly` pairs, and the function-field place
  layer (`forms/function_field.rs`) uses its `divrem`/`gcd`/`pow_mod` (Euler's
  criterion in `F_q[t]/(π)`). As `S[t]` it is the **ring of integers** of `S(t)`,
  so it impls `Scalar` + `HasFractionField` (Frac = `RationalFunction<S>`); its
  units are the nonzero constants, so `inv` is partial.

## The `Scalar` trait + the trait layer (`mod.rs` and friends)

- **`mod.rs`** — the `Scalar` trait (`add`/`neg`/`mul`/`zero`/`one`/`is_zero`/
  `inv`/`characteristic`) + the "any number" table doc + the flat re-export hub.
  Also `impl_scalar_ops!`: every backend gets concrete-type operators (`+ - *`
  and unary `-`) forwarding to the trait methods (so `Surreal + Surreal`,
  `-nimber` work). `/` stays a method (inv is partial). **The operators are NOT a
  `Scalar` supertrait** — see "things that look like bugs".
- **`integrality.rs`** — the (field, ring-of-integers) pairing made structural:
  `HasFractionField {Frac; to_fraction}` + `HasRingOfIntegers {Int; is_integral/
  to_integer}` (with `Int: HasFractionField<Frac=Self>` tying the loop). Impl'd
  for the four distinct-type rows (ℤ⊂ℚ, Oz⊂No, Zp⊂Qp, W_N⊂Qq) PLUS the blanket
  Surcomplex transport (ℤ[i]⊂ℚ[i] falls out). Laurent/Ramified `F_q[[t]]`/`O[π]`
  are same-type valuation subrings, so they stay out (`is_integral` only) — honest.
- **`valued.rs`** — the `Valued` trait: a discrete valuation + canonical
  uniformizer ϖ, impl'd for the local fields/functors (Qp/Qq/Laurent/Gauss and
  Ramified, where the uniformizer is the adjoined π). The spine of the "local
  fields" view (cuts across `small/` + `functor/`). NOT a `Scalar` supertrait
  (rings of integers + exact Archimedean worlds excluded).
- **`analytic.rs`** — the ANALYTIC layer unified as two traits split on where
  precision lives. `ExactRoots {is_square; sqrt}` (no precision arg — exact, or
  exact to the type's K) for Rational, Nimber, Fp, Fpn, Zp, Qp, Qq, WittVec,
  Surreal (exact via the fixed-point bridge over the lazy roots), Laurent, AND the
  blanket `Surcomplex<R: ExactRoots+Ordered>` (the algebraic-closure √(a+bi)).
  `SeriesRoots {sqrt_to_terms; nth_root_to_terms; inv_to_terms}` (caller-chosen n)
  is the lazy interface — Surreal-only (the one world with unbounded, not
  type-fixed, precision). `Ordered {sign}` is the branch-picking datum the
  Surcomplex blanket needs. The residue Tonelli roots (`fp_sqrt`/`fq_sqrt`) live
  here (shared with `small/analytic`'s Hensel seed). Gauss/Ramified excluded
  honestly. NOT a `Scalar` supertrait, like `Valued`.
- **`residue.rs`** — the `ResidueField: Valued` trait: the residue field `k = 𝒪/𝔪`
  (assoc type `Residue`) + two reductions — `residue` (canonical `𝒪 → k`, `None`
  below the integers) and `residue_unit` (the **angular component** `ac(x) ∈ k*`,
  `None` only for 0) + the multiplicative Teichmuller section `τ : k → 𝒪`. Impl'd
  for the local fields/functors (Qp→Fp, Qq→Fq, Laurent→S, Ramified→base residue,
  Gauss→k(tbar)), the last piece of the local-field package `(K,𝒪,𝔪,k,Γ,ϖ)`. It's
  what lets `forms/springer/local.rs` write the discrete Springer decomposition once.
  Globals (Adele/RationalFunction) stay out — per-place residues live at the forms
  layer. NOT a `Scalar` supertrait, like `Valued`.
- **`exactness.rs`** — marker traits for the representation contract: `ExactScalar`,
  `ExactFieldScalar`, and `PrecisionScalar`. Exact finite/truncated rings (`Zp`,
  `WittVec`) are exact as represented; capped-relative models (`Qp`/`Qq`/`Laurent`/
  `Ramified`/`Gauss`/`Adele`) are marked separately and stay out of exact-ring fuzz.
  NOT `Scalar` supertraits — generic code opts in only when it needs the promise.
- **`extension.rs`** — the `FieldExtension: Scalar` trait: a finite separable
  extension `E/F` with `extension_degree`/`embed`/`trace`/`norm` to a distinguished
  `Base`. The orthogonal view of `FiniteField`'s relative trace/norm (one fixed base
  vs. any subfield). Impl'd for `Surcomplex<S: Ordered>` (deg 2), `Fpn<P,N>` over
  `Fp<P>` (deg N, **delegates to** the tested `FiniteField` relative trace/norm), and
  `Qq<P,N,F>` over `Qq<P,N,1>`=`Q_p` (deg F, via the Witt Frobenius `witt_components
  ∘ frobenius ∘ from_witt_components`; base is in the Qq family to dodge the
  Qp-`u128`/Qq-`usize` const-kind clash). Ramified (non-Galois, degenerate trace
  form) and Gauss (transcendental, infinite degree) excluded honestly — the SAME
  boundary `analytic.rs` draws. NOT a `Scalar` supertrait.
  - **`Nimber` impls `FieldExtension`** (`Base = Fp<2>`, degree 128, F_{2^128}/F_2):
    trace = the absolute `nim_trace(·,128)`; norm onto F_2*={1} (1 for nonzero — the
    norm map is trivial). It carried `nim_trace` all along but was absent from the
    interface; this closes "the main char-2 field is missing from `FieldExtension`".
  - **`CyclicGaloisExtension: FieldExtension`** adds the Galois data the twisted
    trace form needs: `basis()` (an F-basis) + `sigma()`/`sigma_power(k)` (the cyclic
    generator). Impl'd for `Surcomplex` (σ=conj, basis {1,i}), `Fpn` (σ=Frobenius,
    coordinate basis), `Qq` (σ=Witt-Frobenius, Teichmuller-lifted residue basis),
    and `Nimber` (σ=nim-Frobenius, bit basis {1,2,…,2¹²⁷}). The relative trace stays
    `FieldExtension::trace` — σ/basis are the only new data, so it's a thin subtrait.
    Consumed by `forms/trace_form.rs` and by `clifford/frobenius.rs`, where σ is
    turned into a basis-level `LinearMap` for outermorphism spectral checks.
- **`tropical.rs`** — the `Semiring` trait + the tropical number type `Tropical<C>`,
  the `(min,+)/(max,+)` semiring. A SIBLING structure (like `Valued`), **not** a
  `Scalar`: a semiring drops the additive inverse (tropical `⊕` is idempotent), so it
  is not a ring and never enters `clifford/` — the same boundary the game *group*
  hits. The convention is a sealed compile-time marker (`MaxPlus`/`MinPlus`), so the
  two dual semirings are distinct, non-interoperating types sharing one impl body (the
  `Surcomplex<S>`/`Laurent<S,K>` move). The games-pillar payoff (thermography IS
  tropical) lives in `games/tropical_thermography.rs`; the semiring laws are fuzzed in
  both conventions in `tests/tropical_axioms.rs`.

## `exact/` — the Archimedean char-0 base (field + ring of integers)

- **`rational.rs`** — exact ℚ over i128, NOT a game backend: the char-0 scalar
  that validates the geometric product against the known Cl(p,q) classification
  before the exotic backends are trusted. (Overflow is a known limit; the surreal
  backend is the real char-0 home.)
- **`integer.rs`** — exact ℤ, the coefficient ring for the exterior algebra of the
  game group (`games/game_exterior.rs`): games are a ℤ-module, not a ring, so Λ
  over ℤ is the structure that lives on all of game-world. Only ±1 invertible.

## `big/` — the transfinite worlds (the number may be infinite)

- **`cnf.rs`** — `merge_descending`, the descending-CNF canonicalizer parameterized
  by the 3 places surreal & ordinal differ (exponent order: No value-order vs
  ordinal lex; coeff merge: + vs XOR; zero test). Deliberately a shared FUNCTION,
  not a `Cnf<C>` TYPE — the orders/algebras diverge (No is a field, On₂ isn't), so
  a shared type would be a false identity.
- **`surreal/`** — finite-support surreal Hahn/CNF backend (char 0), all `impl Surreal`:
  - `mod.rs` — CNF core: `Vec<(exponent: Surreal, coeff: Rational)>`, recursive
    exponents, Hahn arithmetic `ω^a·ω^b = ω^{a+b}`, Scalar, Debug, `truncate()`.
  - `simplicity.rs` — the {L|R}/simplicity bridge (dyadic): `as_rational`/
    `as_dyadic`/`dyadic_birthday` + `simplest_above`/`_below`/`_between`, floor/frac
    (the Oz bridge).
  - `sign_expansion.rs` — exact `sign_expansion`/`from_sign_expansion` (dyadic,
    round-trips, length = birthday) + `as_ordinal`/`from_ordinal` + the transfinite
    (Gonshor) `SignExpansion` + `birthday_ordinal` + the transfinite inverse.
  - `analytic.rs` — the LAZY field layer (the `SeriesRoots` primitives):
    `inv_to_terms` (Neumann series) + `sqrt_to_terms`/`nth_root_to_terms` (real-closed
    roots to n terms; `Some` iff the leading coeff is a perfect ℚ-power).
- **`omnific.rs`** — the omnific integers Oz: `Omnific(Surreal)`, a transfinite
  commutative RING (not field). The surreal mirror of `Integer`.
- **`ordinal/`** — transfinite (ordinal) NIMBERS On₂, the char-2 mirror of surreal:
  - `mod.rs` — CNF core: `Ordinal = Vec<(exponent: Ordinal, coeff: u128)>`, the lex
    cmp, `as_finite`, `checked_mul`, `checked_inv`, and `Scalar`. The `Scalar::mul`
    route is panic-on-escape, matching the checked Kummer tower boundary; callers
    that need an explicit mathematical boundary use the checked methods.
  - `nim.rs` — char-2 NIM arithmetic: `nim_add` (coeff XOR) COMPLETE; `nim_mul`
    dispatches zero / finite×finite / the generator tower.
  - `tower.rs` — the prime-power generator tower (Conway/Lenstra/DiMuro): a monomial
    `ω^E` keyed by `place m ↦ base-p(m) digit vector`; `⊗` adds digit vectors and
    reduces with the Kummer carries `χ_u^u = α_u`. Non-scalar `α_u` (`α_7=ω+1`, …)
    branch a carry into a *sum*, recursed in by descending place. Carries verified
    `α_u` for primes `u ≤ 47`; `None` past that or at `≥ ω^(ω^ω)` (see root `OPEN.md`).
  - `cantor.rs` — ORDINARY (Cantor) `ord_add`/`ord_mul` (ω+ω=ω·2, 1+ω=ω) — the
    surreal birthday's run-length arithmetic. A distinct algebra, sharing only CNF.

The surreal↔ordinal **mirror** (No char 0 / On₂ char 2, sharing `cnf.rs`) is one of
the project's central symmetries.

## `small/` — the non-Archimedean (p-adic) local world

- **`qp.rs`** — `Qp<const P, const K>`: the p-adic FIELD Q_p (the p-adic mirror of
  ℚ / of Omnific⊂Surreal). `p^val·unit`, char 0, inv total on nonzero. CAPPED-
  RELATIVE precision: mul/inv exact, addition NOT associative across precision
  boundaries (a precision model, like float). EXCLUDED from the exact-ring fuzz.
- **`zp.rs`** — `Zp<const P, const K>`: the p-adic integers Z_p (= Z/p^k), the ring
  of integers of Q_p. A LOCAL RING (p a non-unit), residue field F_p; Cl over it is
  non-semisimple.
- **`qq.rs`** — `Qq<const P, const N, const F>`: the UNRAMIFIED extension Q_q =
  Frac(W_N(F_q)), residue degree F (residue field F_q). To WittVec what Qp is to Zp;
  Qq with F=1 IS Qp.
- **`analytic.rs`** — the p-adic ANALYTIC layer over all four backends (mirror of
  `surreal/analytic`): Hensel-lifted `is_square`/`sqrt` (Newton, ODD p only) + the
  Teichmüller rep τ. These inherent methods are what `ExactRoots` delegates to.

## `finite_field/` — the finite residue worlds (the trichotomy's finite leg)

- **`mod.rs`** — the `FiniteField` TRAIT: the shared Galois engine (degree,
  conjugates, min_poly, relative_trace/_norm, multiplicative_order, is_primitive,
  discrete_log) as default methods. An impl supplies only `frobenius`, integer
  `pow`, `ext_degree`, `group_order`, `group_order_factors`. nimber + fpn both
  impl it — one verified algorithm, two backends.
- **`fp.rs`** — `Fp<const P>`: the prime field F_P (odd char), residue field of Zp.
- **`fpn.rs`** — `Fpn<const P, const N>`: F_{p^N} via a (P,N)-keyed irreducible
  reduction poly. Completes the odd-char tower AND the char-2 odd-degree fields
  nimbers can't reach (F_8); supported `Fpn<2,N>` metrics now classify through the
  char-2 Arf façade. (NB the static `order()` = field order p^N, ≠
  `multiplicative_order(&self)`.)
- **`nimber/`** — On₂ in u128 (= F_{2^128}), split by layer, re-exporting `nim_*`
  flat: `mod.rs` (wrapper + Scalar), `arithmetic.rs` (`nim_add`=XOR; `nim_mul` via
  Fermat-power recursion; `nim_square`/`nim_sqrt`/`nim_inv`), `artin_schreier.rs`
  (`nim_trace` + y²+y=c solver), `galois.rs` (impl FiniteField, with Pohlig–Hellman
  + BSGS overrides for `is_primitive`/`discrete_log`).
- **`wittvec.rs`** — `WittVec<const P, const N, const F>`: Witt vectors W_N(F_q) as
  the truncated unramified ring (Z/p^N)[t]/(f̃). The char-p analogue of Z_p; its
  field of fractions is `small/qq.rs`.

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
  Qp); ring of integers F_q[[t]] = the val≥0 subring. Capped-relative; EXCLUDED
  from the fuzz.
- **`ramified.rs`** — `Ramified<S, const E>` = adjoin a root of xᴱ−ϖ over a Valued
  base. The RAMIFIED local cell Q_p(p^{1/E}), the ramified twin of Qq. Always a
  field (Eisenstein), incl. wild/inseparable p|E. `Valued` with uniformizer π and
  `ResidueField` with the base residue field. EXCLUDED from the fuzz.
- **`gauss.rs`** — `Gauss<S>` = S(t) with the Gauss valuation (v(t)=0, transcendental
  residue ⇒ residue field k(t̄)). The last corner, Laurent's residue-extending twin.
  `Valued` itself and `ResidueField` as `k(tbar)`; precision model, EXCLUDED.

## `global/` — the adelic/global place

`Adele` is a finite-precision restricted-product model over ℚ, with `LocalQp` as
the runtime-prime p-adic cell. Useful for product-formula / Hilbert-reciprocity /
Hasse–Minkowski experiments in `forms/adelic.rs`; not an exact infinite-memory
adele. `LocalQp` (runtime prime, NOT const-generic) is the analogue of
`forms`'s runtime `FiniteFieldForm`.

`RationalFunction<S>` (in `global/function_field.rs`) is the **equal-characteristic
mirror**: the global function field `F_q(t)`, the char-`p` analogue of `ℚ` as a
global field. Same field-of-fractions arithmetic as `Gauss` (over `Poly`, `inv =
den/num`, cross-mult equality) but a different ROLE — it carries *all* its place
valuations at once, so like `Adele` it is deliberately **not** `Valued`. Unlike
the precision-model functors it is **exact**, so it *joins* the `scalar_axioms`
fuzz and carries the `ExactScalar`/`ExactFieldScalar` markers. It feeds
`forms/function_field.rs` (the `forms/padic`+`forms/adelic` mirror).

## Things that look like bugs but are not (scalar layer)

- **Scalar `+ - *` operators are concrete-only, NOT a `Scalar` supertrait.**
  Making `Scalar: Add+Sub+Mul+Neg` brings the ops into scope for every generic
  `S`, where `Mul::mul(self, Self)` shadows `Scalar::mul(&self, &Self)` at
  owned-receiver sites and forces clones the borrow-based engine avoids (70+
  generic sites broke when tried). Don't promote them; don't migrate the engine's
  `.add()`/`.mul()` to operators.
- **`ExactRoots`/`SeriesRoots`/`Ordered`/`Valued`/`ResidueField`/the exactness markers
  are NOT `Scalar` supertraits.**
  Not every world takes roots or has a valuation, so the bounds stay opt-in. The
  trait impls *delegate to inherent methods of the same name* (inherent-shadows-
  trait makes that delegate-not-recurse).
- **`Tropical` has no `neg`/`inv` and is deliberately not a `Scalar`.** A semiring's
  `⊕` is idempotent (`a ⊕ a = a`), so there is no additive inverse — that is the
  defining difference from a ring, and the reason `Semiring` is a sibling trait, not a
  `Scalar` supertrait. `Tropical` never reaches `clifford/` (a Clifford algebra needs
  a commutative *ring*), exactly the boundary the game group hits. The two conventions
  are distinct types on purpose: `Tropical<MaxPlus>` and `Tropical<MinPlus>` do not
  interoperate, because thermography's two walls live in dual semirings.
- **`Surreal` has two square roots, by design.** `sqrt_to_terms(n)` is the lazy
  `SeriesRoots` primitive; `ExactRoots::sqrt(&self)` (0 args) is the exact value.
  Different arities, different precision contracts — don't unify them. (Python:
  `Surreal.sqrt(n)` lazy, `Surreal.exact_sqrt()` exact.)
- **`ExactRoots::sqrt`/`is_square` on `Zp`/`Qp`/`Qq`/`WittVec` panic at p=2.** They
  inherit the inherent odd-p assertion (the dyadic case is the forms mod-8 story).
  The finite fields and `Laurent` handle char 2 natively.
- **Surcomplex over nimbers is degenerate.** `i²=1`, `(1+i)²=0`, not a field.
  Surcomplex is only meaningful over char-0 worlds.
- **Surreal coefficients are ℚ, not ℝ.** The honest finite truncation of true CNF.
  Exponents *are* fully recursive surreals. `√2`, `√(2ω)` are honestly `None` (the
  leading coeff must be a perfect ℚ-power); `√ω = ω^{1/2}` IS exact (monomial).
- **`Surreal::inv` returns `None` for any non-monomial.** `1/(ω+1)` is an infinite
  Hahn series; finite support can't hold it.
- **`Surreal::birthday_ordinal`/`transfinite_sign_expansion` are `None` outside the
  representable subclass** (`√ω`, `ω−1`, `½ω`, mixed). Every *ordinal* (incl. ω^ω)
  is handled; `ε` is the one infinitesimal pinned. The honest Gonshor scope boundary.
- **`Qp` addition is not associative across precision boundaries.** Capped-relative
  (the standard p-adic model, like float). No finite-memory exact Q_p exists.
- **`nim_mul`'s `1u128 << (1u128 << n)` is not overflow-prone** for valid u128:
  bit positions < 128 ⇒ Fermat indices n ≤ 6, shift ≤ 64.
- **`Fpn::order()` is the field order `p^N` (static, no self); the element's
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
  (`add`/`mul`/`neg`/`zero`/`one`/`is_zero`). Not duplication to "clean up": the
  inherent methods are what `Gauss`/`RationalFunction`/the place layer call by
  value, and they SHADOW the trait at every concrete `Poly::…` site (so the trait
  bodies delegate, not recurse). The `Scalar` impl exists only so `Poly = F_q[t]`
  can be the `HasFractionField` ring of integers of `RationalFunction = F_q(t)`.

## Math facts worth not re-deriving

- nim-field: `F_{2^{2^k}}` = nimbers `< 2^{2^k}`. `F_n ⊗ F_n = (3/2)F_n` for a
  Fermat 2-power `F_n = 2^{2^n}`; distinct Fermat powers multiply ordinarily.
- Surreal CNF = finite-support Hahn series with ℚ coefficients; the ω-map is the
  monomial map and `ω^a·ω^b = ω^{a+b}` is a group homomorphism on represented
  monomials.
