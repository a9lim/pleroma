# Cross-pillar bridges — TBD (proposed & deferred)

This is the *unbuilt* half of the cross-pillar bridge map: connections whose
mathematics is standard but which are **not yet implemented** — either *proposed* (a
concrete buildable surface) or *deferred* (real and on-thesis, but a larger build not
slated into the current order). It is one of three companion documents:

- **`roadmap/DONE.md`** — the bridges already built and tested (first wave A–D, second
  wave E/F/H/I, third-wave J), each with its formal appendix.
- **`roadmap/TBD.md`** (this file) — the deferred bridge **G** (spinor genus), the
  proposed bridge **K** (the full `ℚ/ℤ` cyclic-algebra Brauer invariant), the
  deferred bridge **L** (the char-`p` Drinfeld/Carlitz mirror of the integral
  pillar), and the proposed **fourth wave**: **M** (the Brown `ℤ/8` invariant —
  the char-2 cell of the mod-8 spine), **N** (a unification pass over
  already-shipped parts), and **O** (lexicodes — the games ↔ integral edge). The
  formalization-pass draft for Bridge K is appended after its entry.
- **`OPEN.md`** — genuine research problems with no known answer.

Claim-level discipline (`AGENTS.md` → "Claim levels and non-claims") still applies:
every piece here is **standard math made computational**, the same status the built
bridges shipped at — not a new theorem. References to "the now-built Bridge F",
Bridge B, Bridge C, etc. point at `roadmap/DONE.md`.

## G — spinor genus (deferred, noted for completeness)

Refining `genus → spinor genus → isometry class` via the spinor norm is classical
(Eichler; Cassels–Hall), and the `clifford/spinor_norm.rs` map is the right
primitive in spirit. But it is **not buildable from the current surface**:
`spinor_norm` computes one versor's norm, whereas the spinor genus needs the local
spinor-norm *images* `θ(O(L ⊗ ℤ_p))` at every prime plus adelic class-group
bookkeeping and the proper/improper class distinction. The one cheap, honest piece
is **Eichler's theorem** as a documented predicate — *indefinite, rank ≥ 3* ⇒ spinor
genus = isometry class — which would let `Genus` upgrade to a class statement in
exactly that regime. The full definite-lattice computation is a larger build; it
stays out of the second wave, adjacent to `OPEN.md` rather than scheduled here.


---

# Third wave — K proposed, L deferred

The third-wave review ("deepen, don't sprawl") produced three bridges; **J** is built
(`roadmap/DONE.md`). The remaining two are here:

```
   CyclicGaloisExt ──cyclic algebra (χ,a)── brauer (full ℚ/ℤ) ──norm form── trace_form     (K)
   F_q[t] ⊂ F_q(t) ──Carlitz / Drinfeld── (char-p mirror of) integral/{theta,modular,codes} (L)
```

Bridge **K** lifts the existing 2-torsion Brauer surface to the full `Br(K_v) = ℚ/ℤ`
image via cyclic algebras built from the Galois data Bridge C already exposes; it
shares a class type with the now-built Bridge F (`Brauer2Class` is its 2-torsion
½-slice). Bridge **L** is the deferred large wing — the char-`p` mirror of the whole
integral pillar — noted for completeness like Bridge G.

## Bridge K — cyclic algebras: the full `ℚ/ℤ` Brauer invariant from the Galois data

**Pillars:** `scalar/…CyclicGaloisExtension` ↔ a new rational/cyclic Brauer class in
`forms/witt/` ↔ `forms/local_global/adelic` (the exact sequence) ↔ `forms/trace_form`
(the norm form).
**Claim level:** PROPOSED — standard math (local class field theory; the cyclic-algebra
invariant map; Serre, *Local Fields*). Lifts the **2-torsion** Brauer surface already in
`adelic.rs` to the full **`Br(K_v) = ℚ/ℤ`** image. The natural completion of the
Brauer thread (and the home Bridge F's rational Clifford invariant sits inside).

### Context: what already exists, and the cap

`local_global/adelic.rs` already builds `brauer_local_invariants` (`inv_v ∈ {0, ½}`),
`brauer_invariant_sum`, and documents the fundamental exact sequence
`0 → Br(ℚ) → ⊕_v Br(ℚ_v) → ℚ/ℤ → 0`. But the local invariant only sees **quaternion**
(degree-2, 2-torsion) classes, so the sequence is realized only in its `½ℤ/ℤ` shadow.

### The mathematics

A cyclic extension `E/K` of degree `n` with a distinguished generator `σ` and an element
`a ∈ K*` defines the **cyclic algebra** `(χ_σ, a) = ⊕_{i<n} E·uⁱ`, with `uⁿ = a` and
`u·x = σ(x)·u`. Its class generates `Br(E/K)`, and when `E/K_v` is **unramified** with `σ`
the arithmetic Frobenius, the local **invariant map** sends `(χ_σ, a) ↦ v(a)/n ∈
(1/n)ℤ/ℤ ⊂ ℚ/ℤ` — the *full* local Brauer group, not just its 2-torsion. So the project
already owns every input — the cyclic Galois data (`σ`, the basis), the local valuations,
the reciprocity sum — and is one constructor away from the full invariant.

Three corrections the formalization pass pinned (full statements in the appendix below):

- **Ramified caveat (load-bearing).** `v(a)/n` is the invariant *only* when `E/K_v` is
  **unramified**; the ramified case needs the general local symbol. Scope the surface to
  unramified-at-`v` data — it suffices for everything below.
- **Where full-strength reciprocity lives.** Over `ℚ`, Minkowski forces every cyclic
  `E/ℚ` of degree `>1` to ramify somewhere, so an `n>2` reciprocity test over `ℚ` needs
  ramified symbols. The clean route is `F_q(t)`: the **constant extension** `F_{qⁿ}(t)`
  is unramified at *every* place, `Frob_v = σ^{deg v}`, and `Σ_v inv_v = (1/n)·deg(div a)
  = 0` — full `ℚ/ℤ` reciprocity reduces to "principal divisors have degree 0", the
  product formula the function-field layer already embodies.
- **The `trace_form` tie is loose as a one-liner.** `Nrd` is degree-`n`, not quadratic;
  the quadratic companion is the algebra trace form `T_A(z) = Trd(z²)`, which
  `assemble_twisted_form` already builds block-by-block. Honest cases: `n=2` char≠2 gives
  `Nrd ≅ ½Q₁ ⟂ (−a/2)Q₁`; `n=2` char 2 *is* the Artin–Schreier symbol Pfister form
  already shipped in `function_field_char2.rs`. So `cyclic_algebra_trace_form` is a
  composition, not new math.

### Proposed surface

- generalize the (proposed Bridge F) `Brauer2Class` to
  `BrauerClass { local: BTreeMap<Place, Rational /* in ℚ/ℤ */> }` with additive
  (mod-`ℤ`) law; the quaternion case is the `½` slice. (`Place` needs an `Ord` derive.)
- `cyclic_algebra_invariant(E, a) -> Rational` `= v(a)/n (mod 1)` for the **unramified**
  local class; `None` on the capped-precision boundary (never a wrong value).
- `constant_extension_invariants(n, a)` over `F_q(t)` — `inv_v = deg(v)·v(a)/n`, the exact
  full-`ℚ/ℤ` reciprocity oracle (everywhere unramified, no ramified symbols needed).
- tie `(χ_σ, a)`'s **trace form** `T_A(z) = Trd(z²)` to `trace_form` as the independent
  oracle (the degree-2 norm-form identity is the cleanest instance).

### Oracles / proposed tests

- Reciprocity at full strength: `Σ_v inv_v ≡ 0 (mod ℤ)` for degree-`n` cyclic classes,
  not only for `½`.
- the degree-2 cyclic class reproduces the existing quaternion `brauer_local_invariants`.
- an unramified cyclic class has `inv_v = 0` at the good places.
- Bridge F's rational Clifford invariant embeds as the 2-torsion part — the two proposed
  bridges share one class type, F supplying the char-0 Clifford correction and K the full
  `ℚ/ℤ` lift.

### Scope / caveats

- **Unramified-at-`v` only** for the `v(a)/n` formula (ramified local symbols are out of
  scope; the `F_q(t)` route delivers full `ℚ/ℤ` strength without them). Reads only `v(a)`,
  `n`, `deg(v)`, so the invariant is **exact** even over the capped-precision local models.
- **Finite legs carry no Brauer content.** Over `Nimber`/`Fpn` every central simple algebra
  splits (Wedderburn), so the Gold forms have no `inv`; their classifier is Arf/Brauer–Wall
  (Bridge B). Bridge K lives only on the local/global legs (`Qq`, `Adele` places, `F_q(t)`, `ℝ`).
- This is the **ungraded** Brauer group; keep it distinct from the graded `BrauerWallClass`
  exactly as the Bridge F section insists. Full lemmas, the convention fix (arithmetic
  Frobenius, `χ_σ(σ)=+1/n`), and the proposed tests are in the appendix below.


---

## Bridge K — formal statements and proposed surface (formalization-pass appendix)

> Moved here from the former `BRIDGES-DRAFT.md` (a parallel formalization front).
> Standard math made computational unless marked; this is the full draft behind the
> Bridge K entry above.

**Status:** PROPOSED. Every theorem below is **standard math** (local/global class field theory); the bridge consists of making it computational on surfaces the crate already ships. The shipped inputs it builds on are labeled **implemented-and-tested** where cited. Nothing here is a new theorem, an Arf/Gold claim, or a graded (Brauer–Wall) statement.

**Pillars:** `scalar/extension.rs` (`CyclicGaloisExtension`: `Surcomplex`, `Fpn<P,N>`, `Qq<P,N,F>`, `Nimber`) ↔ a new ungraded Brauer class in `forms/witt/` ↔ `forms/local_global/adelic.rs` (`brauer_local_invariants`, `brauer_invariant_sum`) ↔ `forms/trace_form.rs` (`trace_twisted_form`) ↔ `forms/local_global/function_field{,_char2}.rs` (places, valuations, the Artin–Schreier symbol).

---

## 1. The cyclic algebra *(standard math)*

Let $E/K$ be a cyclic Galois extension of degree $n$ with a distinguished generator $\sigma$ of $\mathrm{Gal}(E/K)$, and let $\chi_\sigma : \mathrm{Gal}(E/K) \to \frac{1}{n}\mathbb{Z}/\mathbb{Z}$ be the character with $\chi_\sigma(\sigma) = \tfrac1n$. For $a \in K^\times$ the **cyclic algebra** is

$$(\chi_\sigma, a) \;=\; \bigoplus_{i=0}^{n-1} E\,u^i, \qquad u^n = a, \qquad u\,x = \sigma(x)\,u \quad (x \in E),$$

a central simple $K$-algebra of degree $n$ (dimension $n^2$), containing $E$ as a maximal subfield. Standard properties (Gille–Szamuely, *Central Simple Algebras and Galois Cohomology*, Ch. 2):

- $(\chi_\sigma, a) \otimes_K (\chi_\sigma, b) \sim (\chi_\sigma, ab)$ in $\mathrm{Br}(K)$;
- $(\chi_\sigma, a)$ splits $\iff a \in N_{E/K}(E^\times)$; in particular $(\chi_\sigma, N_{E/K}(x))$ splits;
- $a \mapsto [(\chi_\sigma, a)]$ induces an isomorphism $K^\times/N_{E/K}(E^\times) \xrightarrow{\sim} \mathrm{Br}(E/K)$;
- for $n = 2$, $E = K(\sqrt d)$ (char $\neq 2$): $(\chi_\sigma, a)$ **is** the quaternion algebra $(d, a)_K$; in char 2, $E = K(\wp^{-1}(d))$: it is the Artin–Schreier symbol algebra $[d, a)$ already implemented in `function_field_char2.rs`.

The crate's `CyclicGaloisExtension` trait carries exactly the defining data: `basis()` (the $K$-basis of $E$), `sigma()`, `sigma_power(k)`, plus `FieldExtension::{trace, norm, extension_degree}`.

## 2. The local invariant *(standard math)*

Let $K$ be a nonarchimedean local field with normalized valuation $v$, and let $E/K$ be **unramified** of degree $n$ with $\sigma$ the arithmetic Frobenius (inducing $x \mapsto x^{|\kappa|}$ on the residue field). Then the invariant isomorphism $\mathrm{inv}_K : \mathrm{Br}(K) \xrightarrow{\sim} \mathbb{Q}/\mathbb{Z}$ of local class field theory satisfies

$$\boxed{\;\mathrm{inv}_K\big[(\chi_\sigma, a)\big] \;=\; \frac{v(a)}{n} \pmod{\mathbb{Z}}\;}$$

and every class in $\mathrm{Br}(K)$ arises this way (every central simple algebra over a local field has an unramified splitting field). References: Serre, *Local Fields* (GTM 67), Ch. XII; Gille–Szamuely §6.3–6.4; Reiner, *Maximal Orders*, §31. Consequences pinned by the formula: $(\chi_\sigma, a)$ splits at $K$ iff $n \mid v(a)$; the image is the full cyclic group $\frac1n\mathbb{Z}/\mathbb{Z}$, not just its 2-torsion.

**Convention warning.** The sign of $\mathrm{inv}$ depends on choosing the *arithmetic* Frobenius and $\chi_\sigma(\sigma) = +\frac1n$; the geometric-Frobenius convention negates it. The crate's `sigma()` impls (`Fpn::frobenius`, the Witt–Frobenius on `Qq`, nim-squaring on `Nimber`) are all arithmetic, so $+v(a)/n$ is the consistent choice. Reciprocity ($\S3$) is convention-independent; degree-2 compatibility ($\S4$) is not — fix it once, test it.

**Archimedean place.** $\mathrm{Br}(\mathbb{R}) = \frac12\mathbb{Z}/\mathbb{Z}$; for $E = \mathbb{C}$, $\sigma$ = conjugation, $\mathrm{inv}_\mathbb{R}[(\chi_\sigma, a)] = \tfrac12$ iff $a < 0$. There is no valuation to read; this place is special-cased exactly as `brauer_local_invariants` already does via the real Hilbert symbol. $\mathrm{Br}(\mathbb{C}) = 0$.

**Ramified caveat (load-bearing).** If $E/K_v$ is *ramified*, $v(a)/n$ is **not** the invariant; the general local symbol is needed. The proposed surface below is scoped to unramified-at-$v$ data, which suffices for everything in §5–§7.

## 3. Global reciprocity *(standard math)*

For a global field $K$ (number field or function field), the Albert–Brauer–Hasse–Noether exact sequence

$$0 \longrightarrow \mathrm{Br}(K) \longrightarrow \bigoplus_v \mathrm{Br}(K_v) \xrightarrow{\;\sum_v \mathrm{inv}_v\;} \mathbb{Q}/\mathbb{Z} \longrightarrow 0$$

(Reiner §32; Tate, "Global class field theory", in Cassels–Fröhlich, *Algebraic Number Theory*, Ch. VII) gives, for every central simple $K$-algebra $A$:

$$\sum_v \mathrm{inv}_v(A \otimes_K K_v) \;\equiv\; 0 \pmod{\mathbb{Z}},$$

with $\mathrm{inv}_v(A) = 0$ for all but finitely many $v$. For a global cyclic class $(\chi_\sigma, a)$ and a place $v$ unramified in $E$ with $\mathrm{Frob}_v = \sigma^{m_v} \in \mathrm{Gal}(E/K)$, the local term is

$$\mathrm{inv}_v\big[(\chi_\sigma,a)\big] \;=\; \frac{m_v \, v(a)}{n} \pmod{\mathbb{Z}}.$$

**Scope fact, not a gap:** over $\mathbb{Q}$, by Minkowski's theorem every cyclic $E/\mathbb{Q}$ of degree $>1$ ramifies somewhere, so a *full-strength* $n>2$ reciprocity test over $\mathbb{Q}$ would require ramified-place symbols. The crate already owns the clean alternative: over $K = \mathbb{F}_q(t)$ (`RationalFunction` / `FFPlace`), the **constant extension** $E = \mathbb{F}_{q^n}(t)$ is unramified at *every* place (including $\infty$), with $\mathrm{Frob}_v = \sigma^{\deg v}$, so

$$\sum_v \mathrm{inv}_v \;=\; \frac1n \sum_v \deg(v)\, v(a) \;=\; \frac1n \deg\big(\mathrm{div}(a)\big) \;=\; 0,$$

i.e. full $\mathbb{Q}/\mathbb{Z}$-strength reciprocity reduces to "principal divisors have degree 0" — the product formula the function-field layer already embodies. (The Brauer group of $\mathbb{F}_q(t)$ via residues: Faddeev's sequence, Gille–Szamuely §6.4, using $\mathrm{Br}(\mathbb{F}_q) = 0$.)

## 4. How this lifts the shipped 2-torsion surface

**Implemented and tested today** (`forms/local_global/adelic.rs`): `brauer_local_invariants(a, b) -> Option<Vec<(Place, Rational)>>` with values in $\{0, \tfrac12\}$ — the local invariants of the *quaternion* class $(a,b)_\mathbb{Q}$, $\mathrm{inv}_v = \tfrac12 \iff (a,b)_v = -1$ — and `brauer_invariant_sum`, whose vanishing mod $\mathbb{Z}$ is Hilbert reciprocity stated additively. This realizes the exact sequence of §3 only in its $\frac12\mathbb{Z}/\mathbb{Z}$ shadow.

The lift: quaternions are precisely the $n = 2$ cyclic algebras. For $p$ odd and $d$ a nonsquare unit at $p$, $E = \mathbb{Q}_p(\sqrt d)$ is the unramified quadratic extension and

$$\mathrm{inv}_p\big[(\chi_\sigma, a)\big] = \frac{v_p(a)}{2} \equiv \tfrac12\,[\,v_p(a) \text{ odd}\,], \qquad (d,a)_p = \Big(\frac{d}{p}\Big)^{v_p(a)} = (-1)^{v_p(a)},$$

so the degree-2 cyclic invariant reproduces the shipped quaternion invariant place-by-place (at $p = 2$ take $d = 5$; at $\infty$, §2's special case). The new class type replaces "a set of ramified places" by "a $\mathbb{Q}/\mathbb{Z}$-valued divisor of places", and the shipped surface becomes its $\{0,\tfrac12\}$ slice.

## 5. Bridge F as the 2-torsion part

Bridge F's proposed `Brauer2Class { ramified: BTreeSet<Place> }` with symmetric-difference addition embeds via

$$\texttt{ramified} \;\longmapsto\; \Big(v \mapsto \tfrac12\,[\,v \in \texttt{ramified}\,]\Big),$$

a group monomorphism onto the 2-torsion of $\bigoplus_v \mathbb{Q}/\mathbb{Z}$ (XOR of indicator sets $=$ addition of $\tfrac12$'s mod 1). Quadratic-form Brauer classes are 2-torsion, so **all** of Bridge F (Hasse–Witt $s(q)$, the even-Clifford class $c(q)$, and the Lam Prop. V.3.20 $n \bmod 8$/disc correction between them) lands inside the Bridge K class type; K supplies the full-$\mathbb{Q}/\mathbb{Z}$ ambient group and the $n>2$ classes F cannot see. One shared type, two constructors. The reciprocity law specializes correctly: "sum of invariants $\equiv 0$" restricted to the $\tfrac12$-slice is "$|\texttt{ramified}|$ even".

Keep this **ungraded** Brauer class strictly distinct from the graded `BrauerWallClass` in `forms/witt/brauer_wall.rs`, exactly as the Bridge F section insists.

## 6. The tie to `trace_form.rs` *(standard math; the precise statements)*

The Bridge K entry's one-line gloss ("the reduced norm form of $(\chi_\sigma,a)$ *is* the twisted trace form") is loose; the honest statements are:

**(a) $n = 2$, char $\neq 2$.** $\mathrm{Nrd}(x + yu) = N_{E/K}(x) - a\,N_{E/K}(y)$. Since $x\sigma(x) \in K$, the shipped twisted form satisfies $Q_1(x) := \mathrm{Tr}_{E/K}(x\,\sigma(x)) = 2\,N_{E/K}(x)$, hence

$$\mathrm{Nrd} \;\cong\; \tfrac12\,Q_1 \;\perp\; \big(-\tfrac a2\big)\,Q_1 .$$

Pinned instance: `trace_twisted_form::<Surcomplex<Rational>>(1)` $= \langle 2,2\rangle$ (the existing test `surcomplex_twist_is_the_norm_form`), giving $\mathrm{Nrd}\big[(-1,a)_\mathbb{Q}\big] = \langle 1,1,-a,-a\rangle$ — and $(\chi_\sigma,a)$ splits at $v$ iff this form is isotropic over $K_v$ iff $\mathrm{inv}_v = 0$. The norm form is the **independent oracle** for the degree-2 invariant.

**(b) $n = 2$, char 2.** Here $Q_1(x) = \mathrm{Tr}(x\sigma(x)) = 2N(x) = 0$ identically and $\mathrm{Tr}(x^2)$ has vanishing polar — both degenerations `trace_form.rs` already documents as the char-2 trap. The reduced-norm form of $[d, a)$ is instead the 2-fold quadratic Pfister form $[1,d] \perp a\,[1,d]$, **already implemented** in `function_field_char2.rs` with Schmid's residue formula (Serre, *Local Fields*, XIV §5; Gille–Szamuely §9.2) for the local symbol — that layer *is* the char-2, $n=2$ instance of Bridge K, shipped.

**(c) General $n$.** $\mathrm{Nrd}$ is a degree-$n$ form, not quadratic; the quadratic companion is the algebra trace form $T_A(z) = \mathrm{Trd}(z^2)$. Since $\mathrm{Trd}$ kills $E u^i$ for $i \not\equiv 0$ and restricts to $\mathrm{Tr}_{E/K}$ on $E$, $T_A$ decomposes over the lines $Eu^i$ (collecting $i + j \equiv 0 \bmod n$):

$$T_A \;\cong\; Q_0 \;\perp\; \Big(\perp_{0<i<n/2} M_i\Big) \;\perp\; \big[\,n \text{ even}: \; \mathrm{Tr}_{E/K}(a\,x\,\sigma^{n/2}(x))\,\big],$$

where $Q_0(x) = \mathrm{Tr}(x^2)$, the middle block is the $a$-scaled $\sigma^{n/2}$-twist, and $M_i$ is the metabolic pairing $Eu^i \times Eu^{n-i} \to K$, $(x,y) \mapsto \mathrm{Tr}_{E/K}\big(a(x\,\sigma^i(y) + y\,\sigma^{n-i}(x))\big)$. Every block is an instance of the crate's `assemble_twisted_form` core — so `trace_form.rs` already contains the assembler for $T_A$, and a `cyclic_algebra_trace_form` constructor is a composition, not new math.

**(d) Non-tie, for honesty.** Over the finite legs (`Nimber`, `Fpn`) every central simple algebra splits (Wedderburn), so the Gold forms $Q_a$ carry **no** Brauer content; their classifier is Arf/Brauer–Wall (Bridge B), not $\mathrm{inv}$. Bridge K's invariant lives only on the local/global legs (`Qq`, `Adele`-places, $\mathbb{F}_q(t)$, $\mathbb{R}$).

## 7. Proposed surface

```rust
// forms/witt/brauer.rs  (shared with Bridge F)
pub struct BrauerClass {
    /// inv_v ∈ ℚ/ℤ, canonical representative in [0,1); zero entries omitted,
    /// so the split class is the empty map (matching Brauer2Class's ∅).
    pub local: BTreeMap<Place, Rational>,
}
impl BrauerClass {
    pub fn add(&self, other: &Self) -> Self;          // entrywise, mod ℤ, drop zeros
    pub fn invariant_sum(&self) -> Rational;          // ≡ 0 mod ℤ for global classes
    pub fn from_quaternion(ramified: &BTreeSet<Place>) -> Self;   // the ½-slice (Bridge F)
    pub fn two_torsion(&self) -> Option<BTreeSet<Place>>;          // back down, when it is one
}

/// inv = v(a)/n mod ℤ for the unramified local cyclic class (χ_σ, a),
/// E = Qq<P,N,F> over Q_p = Qq<P,N,1>, σ = the Witt–Frobenius, n = F.
/// None on the capped-precision Option boundary (a not invertibly represented).
pub fn cyclic_algebra_invariant<E: CyclicGaloisExtension>(a: &E::Base) -> Option<Rational>
where E::Base: Valued;

/// inv_v = deg(v)·v(a)/n mod ℤ over F_q(t) with E = F_{q^n}(t) (constant extension,
/// everywhere unramified, Frob_v = σ^{deg v}); exact.
pub fn constant_extension_invariants<S: FiniteOddField>(
    n: u128, a: &RationalFunction<S>,
) -> Option<Vec<(FFPlace<S>, Rational)>>;
```

Implementation notes: `Place` (in `padic.rs`) currently derives only `PartialEq, Eq` — keying a `BTreeMap` needs `Ord` (derive it; document that `Real` sorts per declaration order). All invariants are tiny exact `Rational`s ($i128$-backed); the construction reads only $v(a)$, $n$, $\deg v$, so it is **exact even over the capped-precision local models**, with `None` (never a wrong value) when precision loss hides $v(a)$.

## 8. Proposed tests / oracles

1. **Degree-2 compatibility** *(the lift is a lift)*: for $p$ odd, $d$ a nonsquare unit mod $p$ (and $d=5$ at $p=2$), `cyclic_algebra_invariant` over the unramified quadratic equals the entry of the shipped `brauer_local_invariants(d, a)` at $p$, across a sweep of $a$ with $v_p(a) \in \{0,1,2,3\}$.
2. **Splitting law**: $\mathrm{inv} = 0 \iff n \mid v(a)$; in particular $(\chi_\sigma, \text{unit}) $ splits (the "unramified class at good places" oracle) and $(\chi_\sigma, N_{E/K}(x))$ splits for sampled $x$ (norms via the existing `FieldExtension::norm`).
3. **Additivity / $n$-torsion**: $\mathrm{inv}(ab) = \mathrm{inv}(a) + \mathrm{inv}(b) \bmod \mathbb{Z}$; $n \cdot \mathrm{inv}(a) \equiv 0$; the image for fixed $n$ is exactly $\frac1n\mathbb{Z}/\mathbb{Z}$ (full local Brauer group, not 2-torsion).
4. **Full-strength reciprocity** over $\mathbb{F}_q(t)$: for constant extensions of degree $n \in \{2,3,4,5\}$ and random $a \in \mathbb{F}_q(t)^\times$, $\sum_v \deg(v)\,v(a)/n \equiv 0 \bmod \mathbb{Z}$ — discover-don't-assert via the place enumeration of `function_field.rs`, with the independent check $\deg(\mathrm{div}(a)) = 0$.
5. **Reciprocity over $\mathbb{Q}$, degree-2 slice**: the existing `brauer_invariant_sum_is_zero_in_q_mod_z` re-read through `BrauerClass::from_quaternion(…).invariant_sum()` — pins the §5 embedding.
6. **Norm-form oracle** ($n=2$, char $\neq 2$): $\mathrm{inv}_v = 0 \iff \langle 1,-d,-a,da\rangle$ isotropic over $\mathbb{Q}_v$ (`try_is_isotropic_at_p`), tying the invariant to the shipped Hasse–Minkowski layer; plus the $\tfrac12 Q_1 \perp (-\tfrac a2)Q_1$ identity of §6(a) against `trace_twisted_form`.
7. **Char-2 cross-check**: the $\{0,\tfrac12\}$ class of $[d,a)$ from the shipped `as_symbol_places` agrees with `BrauerClass` arithmetic, and `as_symbol_reciprocity_sum` is its reciprocity instance.
8. **Bridge F embedding** (once F lands): `from_quaternion` ∘ XOR $=$ `add` ∘ `from_quaternion`; `two_torsion` round-trips.

## 9. Scope and caveats

- **Unramified-at-$v$ classes only** for the $v(a)/n$ formula; ramified local symbols (needed for full-strength $n>2$ reciprocity over $\mathbb{Q}$, by Minkowski) are out of this bridge's minimal scope — the function-field route (§3, test 4) delivers full $\mathbb{Q}/\mathbb{Z}$ strength without them. Document the boundary; don't fake the ramified case.
- **Ungraded Brauer only.** No contact with `BrauerWallClass` / Arf; the finite legs carry no invariant (Wedderburn, §6(d)).
- **Convention is part of the spec**: arithmetic Frobenius, $\chi_\sigma(\sigma) = +\frac1n$ (§2); a sign flip is invisible to every 2-torsion test and to reciprocity, so pin it with an $n \geq 3$ asymmetric case (e.g. $\mathrm{inv} = \frac13$ vs $\frac23$ distinguished via additivity under $a \mapsto a^2$).
- **Claim levels**: §§1–3, 6 standard math (Serre, *Local Fields*, Ch. XII, XIV §5; Gille–Szamuely Ch. 2, §§6.3–6.4, §9.2; Reiner, *Maximal Orders*, §§31–32; Tate in Cassels–Fröhlich Ch. VII; Lam, *Introduction to Quadratic Forms over Fields*, Ch. III, V); §4's existing surface implemented-and-tested; everything in §§7–8 proposed; no interpretation-level or open-level claims are introduced.

---

## Bridge L — the char-`p` mirror of the integral pillar (deferred, large)

**Pillars:** `scalar/global/function_field` (`F_q(t)`, `F_q[t]`) ↔ a large new
Drinfeld/Carlitz layer ↔ `forms/integral/{theta,modular,codes}`.
**Claim level:** PROPOSED but **large** — standard math (Goss, *Basic Structures of
Function Field Arithmetic*; Gekeler, Drinfeld modular forms; Goppa / AG codes). Noted
like Bridge G: real and on-thesis, **not** scheduled into a build order.

### The mirror

The entire `integral/` wing — even-unimodular `ℤ`-lattices, `θ`-series,
`M_*(SL₂ℤ) = ℂ[E₄, E₆]`, Construction-A codes, Leech — is char-0. The project already
ships **exact** `F_q[t] ⊂ F_q(t)`, the char-`p` global field, and its arithmetic carries
a complete mirror of the integral pillar:

- the **Carlitz module** `C_t(x) = t·x + x^q` is the char-`p` analogue of `exp` / the
  lattice exponential; the mirror of `E₄, E₆` are **Drinfeld modular forms** for
  `GL₂(F_q[t])`, with Goss `ζ`-values mirroring the Eisenstein constants.
- rank-`r` `F_q[t]`-lattices mirror even-unimodular `ℤ`-lattices and their reduction
  theory.
- **Goppa / algebraic-geometry codes** from function fields would tie *straight back into
  the existing `codes.rs`* Construction-A machinery — the same code↔lattice seam, read in
  characteristic `p`.

This is the `No ↔ On₂` / char-0 ↔ char-2 move applied to the richest pillar — the most
*on-thesis* possible "new structure," which is exactly why it earns a mention while
smaller additions do not.

### Why deferred

A genuine new wing (Drinfeld modules, the Carlitz exponential, rank-`r` reduction
theory): weeks of work, specialized, and worth starting only if the goal is a *second
headline pillar* rather than finishing the first. Like G, it sits adjacent to the
roadmap, not inside its build order.


---

# Fourth wave — M/N/O proposed

The fourth-wave review asked where the **symmetry table** itself (README → "The
symmetries") is still uneven, rather than where a new number system could go. Three
answers, three proposals:

1. the mod-8 spine has **four char-0 routes and no char-2 cell** — Bridge **M** (the
   Brown invariant);
2. several shipped computations are pairwise-disconnected instances of **one
   classical statement** — Bridge **N** (a pass of four small unifications);
3. the pillar *graph* closed in the first wave, but one **edge** was never drawn:
   games ↔ integral — Bridge **O** (lexicodes).

```
  char2/arf ──β = 4·Arf── Brown β ∈ ℤ/8 ──β ≡ σ (mod 8)── integral/discriminant      (M)
  springer ∂_p ──Milnor exact sequence── W(ℚ), W(F_q(t));  trace_form = Scharlau s_*;
       genus ⟺ (σ, q_L) (Nikulin);  one Bernoulli source for E₄/E₆ and the mass      (N)
  games mex/Grundy ──greedy lexicode── codes.rs ──Construction A── lattices ── theta  (O)
```

Claim-level discipline still applies: every item below is **standard math made
computational**, the same status the built bridges shipped at — not a new theorem.
Where a statement must be transcribed from a source rather than reconstructed, the
entry says so.

## Bridge M — the Brown invariant: the char-2 cell of the mod-8 spine

**Pillars:** `forms/char2/` (Arf) ↔ `forms/integral/discriminant.rs` (Milgram,
Bridge A) ↔ `forms/witt/brauer_wall.rs` (the mod-8 cycle).
**Claim level:** PROPOSED — standard math (E. H. Brown, *Generalizations of the
Kervaire invariant*, Ann. of Math. 95 (1972); C. T. C. Wall, *Quadratic forms on
finite groups, and related topics*, Topology 2 (1963); Nikulin) made computational.

### The asymmetry it repairs

The mod-8 spine currently lives entirely on the char-0 side: the exact rational
signature, the genus oddity (`genus_signature_mod8`), the Milgram Gauss-sum phase
(`milgram_signature_mod8`, Bridge A), and the Weil `S` prefactor (Bridge I) are four
routes to `σ mod 8`. The char-2 side carries only the `ℤ/2` Arf/BW bit. The
classical object filling the char-2 mod-8 cell is the **Brown invariant** of
`ℤ/4`-valued quadratic refinements.

### The mathematics

A **`ℤ/4`-quadratic form** on a finite-dimensional `F₂`-space `V` is `q : V → ℤ/4`
with

```text
q(x+y) = q(x) + q(y) + 2·b(x,y),
```

where `b : V×V → F₂` is bilinear and `2· : F₂ ↪ ℤ/4`. Setting `y = x` forces
`b(x,x) = q(x) mod 2` — so `b` is symmetric **but not alternating**.

**Category trap (load-bearing).** This `b` is *not* the engine's polar form: the
crate's char-2 `Metric` carries an alternating `b` (`b_ii = 0`) with `q` valued in
the field, while Brown's category has `ℤ/4`-valued `q` with `b_ii = q_i mod 2`.
Hard rule 2 (don't collapse `q` and `b`) has a cousin here: don't identify the two
categories. The doubling map below is the only bridge between them.

For `b` nondegenerate, the Gauss sum is a `ℤ[i]`-integer of absolute value `2^{n/2}`:

```text
Σ_{x ∈ V} i^{q(x)} = 2^{n/2} · ζ₈^β,    ζ₈ = e^{2πi/8},
```

and **`β(q) ∈ ℤ/8` is the Brown invariant**: additive under `⊥`, and a complete
invariant up to adding split planes, making the Witt group of the category cyclic of
order 8, generated by `⟨1⟩` (the 1-dimensional form with `q(x) = 1`). **[Pin the
exact stable-equivalence statement from Brown 1972 / Wall 1963 during the
formalization pass; do not paraphrase it into the prose before then.]**

Three identifications make this the missing cell rather than a fifth pillar:

1. **Arf is the 2-torsion.** For a classical nonsingular char-2 form
   `q′ : V → F₂` (alternating polar), the **doubled** form `2q′ : V → ℤ/4` has
   Gauss sum `Σ (−1)^{q′(x)} = (−1)^{Arf} · 2^{n/2}`, so

   ```text
   β(2q′) = 4 · Arf(q′).
   ```

   The shipped Arf bit embeds as `{0, 4} ⊂ ℤ/8` — the char-2 classifier becomes
   the 2-torsion of a `ℤ/8` invariant, mirroring "the real Witt class is the
   2-torsion shadow of the signature".
2. **Milgram on the 2-elementary slice is Brown.** For an even lattice `L` with
   2-elementary `A_L ≅ (ℤ/2)^k`, the discriminant form `q_L` takes values in
   `½ℤ/2ℤ`, and `t ↦ 2t` identifies `(A_L, 2q_L)` with a `ℤ/4`-quadratic form
   whose Brown sum *is* the Milgram Gauss sum. Milgram/van der Blij then reads

   ```text
   β(2·q_L) ≡ sign(L)   (mod 8)
   ```

   — computed from the **integer value-counts** `(n₀ − n₂) + i(n₁ − n₃)`, i.e.
   exact `ℤ/8` arithmetic. That is a **fifth independent route to `σ mod 8`**, and
   the first with no floating point (the `GaussSum` route is `f64`).
3. **The generators are shipped lattices.** `a_n(1)` (= `A₁ = ⟨2⟩`): `A = ℤ/2`,
   `q = ½`, `β = 1 ≡ σ`. `e_7()`: `q = 3/2`, `β = 7 ≡ σ`. `d_n(4)`: three nonzero
   elements of `q`-value 1, sum `1 + 3i² = −2`, `β = 4 ≡ σ`. The `ℤ/8` generator
   `⟨1⟩` is literally the discriminant form of `A₁`.

### Proposed surface

- `forms/char2/brown.rs`
  - input in the `arf_f2` idiom: `brown_f2(n: usize, q4: &[u128] /* mod 4 */,
    bmat: &[u128]) -> BrownResult`, with the constructor-level check
    `b_ii = q_i mod 2`.
  - `BrownResult { beta: u128 /* of the nonsingular core, mod 8 */, rank: usize,
    radical_dim: usize, radical_anisotropic: bool }` — mirroring `ArfResult`
    field-for-field. On the radical of `b`, `q` takes values in `{0, 2}`; `q ≡ 0`
    there ⇒ `beta` is the core invariant with `radical_dim` reported; any radical
    value `2` ⇒ the full Gauss sum vanishes (`radical_anisotropic`), and `beta`
    still reports the core. Data, not a panic.
  - primary route: reduction to `⟨±1⟩` / split summands (the `arf_char2`-style
    reduction); oracle route: direct `2^n` enumeration of the value distribution
    with exact integer phase recovery — the same enumeration budget
    `DiscriminantForm` already pays.
  - `double_f2(...)` — the embedding from `arf_f2` input data;
    `from_discriminant_form(&DiscriminantForm) -> Option<...>` — `Some` only for
    2-elementary groups (read off the invariant factors).

### Oracles / proposed tests

- `β` additivity under `⊥`: reduction route vs enumeration route, fuzzed.
- `brown_f2(double_f2(q′)).beta == 4 * arf_f2(q′).arf` across nonsingular metrics.
- the split objects: the hyperbolic plane `[q(e)=0, q(f)=0, b(e,f)=1]` and
  `⟨1⟩ ⊥ ⟨−1⟩` both have `β = 0`; `β(⟨1⟩^{⊥8}) = 0` (the order-8 relation).
- the lattice slice: `from_discriminant_form` of `a_n(1)`, `e_7()`, `d_n(4)`,
  `d_n(8)` gives `β ≡ signature mod 8`, cross-checked against
  `milgram_signature_mod8`, `genus_signature_mod8`, and the f64 `GaussSum` phase;
  `e_8()` collapses to the empty form, `β = 0 ≡ 8`.

### Scope / caveats

- The lattice tie is **2-elementary discriminant groups only**. Higher 2-power
  torsion needs `ℤ/2^{k+1}`-valued refinements and odd torsion has its own odd
  Gauss sums — both stay with the shipped f64 `GaussSum`. A full
  finite-quadratic-module Witt group (Nikulin's generators and relations) is a
  further rung, not this bridge.
- No new theorem: Brown 1972 is the source; the bridge is the wiring to Arf
  (shipped) and Milgram (Bridge A).

## Bridge N — the unification pass: four small joins of already-shipped parts

**Pillars:** vary per item — each joins surfaces that already exist.
**Claim level:** PROPOSED — standard math; every item is assembly + verification of
shipped machinery, deliberately smaller than a headline bridge. The four items are
independent; ship in any order.

### N.1 — Milnor's exact sequences: the Springer trio goes global

The shipped Springer engine computes per-place residue buckets, and the
local-global layer decides per-form isotropy; the Witt-*group*-level global
statement is assembled nowhere:

```text
0 → W(ℤ) → W(ℚ)     →∂  ⊕_p W(F_p)   → 0     (exact; Milnor)
0 → W(F_q) → W(F_q(t)) →∂  ⊕_P W(k(P)) → 0     (split exact; q odd)
```

(Milnor–Husemoller, *Symmetric Bilinear Forms*, Ch. IV; Lam, GSM 67, Ch. IX —
Milnor's sequence over `k(t)`, split via the `∞`/leading-coefficient place.) For
odd `p`, `∂_p` **is** the shipped second Springer residue, lifted from
`LocalResidueForm` buckets to Witt classes — and the target type already exists:
`WittClassG::OddChar { field_order, kappa, e0, sclass }`. The kernel `W(ℤ) ≅ ℤ` is
detected by the signature, so the sequence ties **three** pillar surfaces together:
Springer residues, the global field layer, and the integral pillar's headline
invariant. Once again the char-`p` mirror is the better-behaved one (split, exact
arithmetic) — the project's recurring punchline.

- **Surface:** `forms/witt/milnor.rs` — `global_residues(entries: &[i128]) ->
  Option<(i128 /* signature */, BTreeMap<u128, WittClassG>)>` over `ℚ`, computing
  `v_p(aᵢ)` and residue square classes exactly from the `i128` entries (Legendre
  symbols), with `springer_decompose_qp` on the capped model as the cross-check;
  the parity-1 layer of the `p`-buckets gives `e0 = dim mod 2` and `sclass` from
  the bucket discriminants — **matching the shipped signed-discriminant `sclass`
  convention, not the raw determinant**. The `F_q(t)` twin keys on `FFPlace` and
  is exact outright.
- **Oracles:** finite support (`∂_p = 0` for `p ∤ 2·∏aᵢ`); reconstruction — equal
  signature + equal residues at every `p` ⇒ Witt-equivalent over `ℚ`, cross-checked
  against the shipped Hasse–Minkowski layer and `witt_decompose`; split exactness
  over `F_q(t)` discover-don't-assert on random diagonal forms.
- **Boundary (load-bearing):** `∂₂` is **not** Springer's (residue char 2); Milnor
  defines it by hand. Either transcribe Milnor's `∂₂` from Ch. IV during the
  formalization pass, or ship the odd-`p` sequence with the `p = 2` component
  documented as the known boundary and reconstruction tests restricted to
  odd-support forms. **Do not guess the `∂₂` convention.** Char-2 constant fields
  of `F_q(t)` stay excluded (that is the separate Aravire–Jacob layer).

### N.2 — name the Scharlau transfer

`trace_twisted_form::<E>(0)` is not just an example generator — it is `s_*(⟨1⟩)`
for the **Scharlau transfer** `s_* : W(E) → W(K)` along `s = Tr_{E/K}` (Lam, GSM 67,
Ch. VII; Scharlau, *Quadratic and Hermitian Forms*, Ch. 2). Naming it buys two
theorems whose tests are shipped-shaped:

- **Frobenius reciprocity** (projection formula): `s_*(r*(x)·y) = x·s_*(y)` in
  `W(K)` — the product via the shipped `tensor_form`.
- **Springer's odd-degree theorem**: `r* : W(K) → W(E)` is injective for odd
  `[E:K]` — the *other* Springer theorem (the 1952 C. R. note), joining the 1955
  residue theorem that runs the local layer. Two Springer theorems, one crate,
  one per pillar seam.

- **Surface:** `transfer_diagonal<E: CyclicGaloisExtension>(entries: &[E]) ->
  Metric<E::Base>` — the Gram of `(x,y) ↦ Tr(λ·x·y)` per entry `λ`, through the
  existing `assemble_twisted_form` core; plus doc-level naming in `trace_form.rs`
  declaring the `k = 0` twisted form as the transfer of `⟨1⟩`.
- **Oracles:** the transfer is a Witt-group map (transfer of a split form is
  split — checked via the classify façade); Frobenius reciprocity fuzzed over
  `Fpn<p,2>/Fp` and `Qq/Qp` for odd `p`; odd-degree injectivity sampled over
  `Fpn<p,3>/Fp` (anisotropic stays anisotropic — both sides decidable by the
  finite-field classifiers).
- **Boundary:** char ≠ 2 legs only. The char-2 degenerations
  (`Tr(x·σ(x)) = 2N = 0`) are exactly the trap `trace_form.rs` already documents;
  the char-2 transfer story is the Artin–Schreier route already shipped in
  `function_field_char2.rs`.

### N.3 — Nikulin: genus ⟺ (signature, discriminant form)

Bridges A and I compare mod-8 *phases*; Nikulin's criterion upgrades the
relationship to a classification equivalence: two **even** lattices lie in the same
genus iff they have equal signature pairs and isomorphic discriminant quadratic
forms (Nikulin, *Integral symmetric bilinear forms and some of their geometric
applications*, Izv. Akad. Nauk SSSR 43 (1979), Cor. 1.9.4). Both sides are shipped
— `are_in_same_genus` (Conway–Sloane symbols, Allcock-corrected at `p = 2`) and
`DiscriminantForm` — the only missing piece is the isomorphism test between finite
quadratic modules.

- **Surface:** `DiscriminantForm::is_isomorphic(&self, other: &Self) ->
  Option<bool>` — match invariant factors, then search value-preserving group
  isomorphisms; budgeted (`None` past the node budget — the
  `automorphism_group_order` precedent).
- **Oracles:** `are_in_same_genus(a, b) == (equal signatures && q_a ≅ q_b)` across
  the shipped zoo (`a_n`, `d_n`, `e_6`/`e_7`/`e_8`, `d16_plus`, direct sums); the
  Milnor pair pins the unimodular case (trivial `A_L`, equal signature ⇒ same
  genus, non-isometric); easy separations (`A₂` vs `A₁ ⊕ A₁`: `ℤ/3` vs `(ℤ/2)²`).
- **Boundary:** even lattices only (Bridge A's boundary); the brute-force budget is
  honest — this is a *cross-check of two shipped routes*, not a p-adic-symbol
  reimplementation. A small **multi-class genus** example would also de-degenerate
  Bridge E's Siegel–Weil rung; finding one is a good shared errand for this item.

### N.4 — one Bernoulli source for Eisenstein and mass

`mass_formula.rs` already computes Bernoulli numbers exactly (Akiyama–Tanigawa, per
`TABLES.md`); `modular.rs` hardcodes `240`/`504` as literals. They are the same
numbers — `E_{2k} = 1 − (4k/B_{2k}) Σ σ_{2k−1}(n) qⁿ` gives `240 = −8/B₄` and
`−504 = −12/B₆` — and the even-unimodular mass constants are Bernoulli products
whose rank-8 instance lands on `1/696729600 = 1/E8_WEYL_GROUP_ORDER`.

- **Surface:** export the Bernoulli helper from one shared home; derive the
  Eisenstein coefficients from it, keeping the literals as the pinned oracle
  (the `TABLES.md` discipline: derived value asserted equal to curated constant).
- **Oracles:** `B₂ … B₂₄` pinned (von Staudt–Clausen denominators as a free extra
  check); `240`/`504` reproduced; `mass_even_unimodular(8)` through the shared
  helper still equals `1/E8_WEYL_GROUP_ORDER`.
- **Boundary:** nothing new mathematically — this item exists purely so the
  `ζ(1−2k)` thread is one computation instead of two private ones.

## Bridge O — lexicodes: greedy = mex, the games ↔ integral edge

**Pillars:** `games/` (mex/Grundy) ↔ `forms/integral/codes.rs` (Bridge H) → thence
Construction A / theta (Bridges H/E).
**Claim level:** PROPOSED — standard math (Conway–Sloane, *Lexicographic codes:
error-correcting codes from game theory*, IEEE Trans. Inform. Theory 32 (1986)
337–348; independently Brualdi–Pless, *Greedy codes*, 1993 — verify venue in the
formalization pass). Closes the one pillar edge still missing from the bridge
graph: games ↔ integral.

### The mathematics

The **lexicode** `L(n, d)`: scan `F₂ⁿ` in lexicographic order, greedily keeping
every vector at Hamming distance ≥ `d` from all vectors kept so far. Conway–Sloane:

- **Binary lexicodes are linear, and the proof is game theory.** The greedy scan is
  the mex rule — the codewords are the Grundy-value-0 positions of a turning-game
  move structure — and XOR-closure is Sprague–Grundy theory, not coding theory.
- More generally, lexicodes over base `B = 2^k` are closed under nim-addition, and
  are **linear over the fields of size `2^(2^k)`** — i.e. precisely over Conway's
  nim-fields, the same `F_{2^{2^k}}` tower the `Nimber` backend realizes
  (`Nimber(u128) = F_{2^{2^7}}`). The fields where greedy respects linear algebra
  are the fields nim-multiplication builds. **[Transcribe the exact statement and
  the specific game from the 1986 paper during the formalization pass; do not pin
  prose from memory.]**
- The celebrated instances are the shipped codes: Hamming codes are lexicodes, and
  the length-24, `d = 8` lexicode is the extended binary **Golay** code — so
  `hamming_code()`, `extended_hamming_code()`, and `golay_code()` each acquire a
  second, game-theoretic construction.

That makes a full-length chain executable: **mex → lexicode → Golay →
Construction A → even unimodular rank 24 with roots (≠ Leech, Bridge H's pinned
distinction) → theta (Bridge E)**. Every link except the first is already shipped.

### The OPEN #1 tie (interpretation level)

`OPEN.md` §1 records that normal-play P-sets are *linear* in Grundy coordinates —
the polar form is exactly the XOR-closure obstruction. Lexicodes are the classical
demonstration of the solved side of that line: a fixed, natural, non-tautological
greedy rule whose P-set is a rich **linear** code. The Gold-quadric question is
then precisely whether the lexicode phenomenon admits a quadratic refinement.
Bridge O supplies the degree-1 case as executable context; it does **not** touch
the open problem and must not be cited as progress on it.

### Proposed surface

- `games/lexicode.rs` (games already depends on forms — `loopy.rs` imports
  `fit_f2_quadratic` — so importing `BinaryCode` keeps the existing direction):
  - `lexicode_naive(n: usize, d: usize) -> Option<BinaryCode>` — the literal
    greedy scan for small `n`, **discover-don't-assert**: collect greedily, check
    XOR-closure, and return `None` on a closure failure (which would falsify the
    linearity theorem rather than hide it), else hand a basis to
    `BinaryCode::new`.
  - `lexicode(n: usize, d: usize) -> Option<BinaryCode>` — the production route,
    using linearity incrementally: the next generator is the lex-least `v` whose
    coset `v + C` has minimum weight ≥ `d`. For `n = 24` that is 12 generator
    searches against ≤ `2^12`-word coset scans with popcount early-exit, not a
    naive `2^24 × 2^12` sweep; budgeted like `automorphism_group_order`.
  - a doc-level statement (plus a toy-`n` witness test) that the greedy choice at
    each step is a mex over the move structure.

### Oracles / proposed tests

- `lexicode_naive == lexicode` wherever both run (`n ≤ 14` across small `d`).
- `lexicode(7, 3)` vs `hamming_code()`, `lexicode(8, 4)` vs
  `extended_hamming_code()` — expect literal equality under the natural bit order;
  on a convention mismatch, assert equivalence (dim + minimum distance + weight
  enumerator) and document the order.
- `lexicode(24, 8)`: assert `[24, 12, 8]`, doubly-even, self-dual, and
  `weight_enumerator` equal to `golay_code()`'s. Uniqueness of the `[24,12,8]`
  Type II code up to permutation (standard math — MacWilliams–Sloane; Pless)
  upgrades the predicate bundle to "is the Golay code".
- sanity rungs: `d = 1` gives `F₂ⁿ`; `d = 2` gives the even-weight code.
- the chain test: `lexicode(24, 8).construction_a()` is even unimodular rank 24
  **with** roots — re-pinning Bridge H's Golay-is-not-Leech boundary from the
  games side.

### Scope / caveats

- Binary only. The nim-field `2^{2^k}` statement is documented context, not a
  shipped surface — a `Nimber`-coefficient lexicode is a separate rung; note the
  temptation and resist it until the binary layer is solid.
- The greedy order is part of the definition: lexicographic in the standard bit
  order. A permuted coordinate order gives a different (equivalent) code —
  document the convention once.
- The turning-game realization must be transcribed from the 1986 paper, not
  reconstructed. The OPEN #1 paragraph above is interpretation-level framing and
  stays subordinate to `OPEN.md`.


---

## TBD — status snapshot

**K, M, N, and O are proposed; G and L are deferred:**

- **K (proposed):** lifts the shipped 2-torsion Brauer surface (`adelic.rs`) to the
  full `ℚ/ℤ` invariant via cyclic algebras built from the Galois data Bridge C
  exposes; shares a class type with the now-built Bridge F (`roadmap/DONE.md`) —
  `Brauer2Class` is its 2-torsion ½-slice. Full formal draft appended above.
- **M (proposed):** the Brown `ℤ/8` invariant of `ℤ/4`-valued quadratic
  refinements — the char-2 cell of the mod-8 spine. Contains the shipped Arf bit
  as its 2-torsion (`β = 4·Arf`) and computes the Milgram phase exactly (no `f64`)
  on 2-elementary discriminant forms.
- **N (proposed):** the unification pass — four independent small joins: Milnor's
  exact sequences (the Springer residues assembled globally over `ℚ` and
  `F_q(t)`), the Scharlau-transfer naming (+ Frobenius reciprocity and Springer's
  odd-degree theorem), Nikulin's genus ⟺ (signature, discriminant form)
  criterion, and one shared Bernoulli source for the Eisenstein and mass
  constants.
- **O (proposed):** lexicodes — greedy = mex, binary lexicodes are linear by
  Sprague–Grundy theory, and the length-24 `d = 8` lexicode is the Golay code.
  Closes the games ↔ integral edge and documents the solved *linear* shadow of
  `OPEN.md` §1.
- **G (deferred):** the spinor-genus refinement `genus → spinor genus → isometry
  class`; classical but not buildable from the current surface. The cheap honest
  piece is Eichler's theorem as a documented predicate (indefinite, rank ≥ 3 ⇒
  spinor genus = isometry class).
- **L (deferred, large):** the char-`p` Drinfeld/Carlitz mirror of the whole
  `integral/` pillar — a genuine second-headline-pillar build, not a task.

Recommended order: **N**'s items are the cheapest and ship independently (N.2/N.4
are afternoon-sized; N.3 needs one budgeted brute-force routine; N.1 carries the
`∂₂` decision); **M** is the highest mathematical payoff per line of code (one new
file, three shipped identifications); **O** first if feeding the `OPEN.md` §1
framing matters most (its formalization pass needs the 1986 paper in hand). **K**
remains the natural completion of the Brauer thread; **L** is a project-scope
decision. The built bridges are in `roadmap/DONE.md`; the genuine open problems
stay in `OPEN.md`.
