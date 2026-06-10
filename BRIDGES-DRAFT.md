# Bridge & lemma drafts (fable assault, ogdoad)

> Draft formal statements produced by the parallel formalization front. Claim-level: standard math made computational unless marked. Triage before pasting into ROADMAP.md / writeups/goldarf.tex.

## Bridge J — valuation as tropicalization; Newton polygons

All surfaces read (ROADMAP Bridge J, OPEN.md, `scalar/valued.rs`, `scalar/tropical.rs`, `scalar/functor/{ramified,gauss}.rs`, `scalar/poly.rs`, `forms/springer/{mod,local}.rs`, `forms/local_global/function_field.rs`, pillar AGENTS.md files). Drafting now — the deliverable follows.

# Bridge J — formal statements and proofs

**Status.** Everything below is **standard math** (no new theorems), per the third-wave discipline in `ROADMAP.md`. Items marked ⟦implemented⟧ are witnessed by tests in this checkout; items marked ⟦proposed⟧ name the tests that would witness the proposed `NewtonPolygon` surface. Nothing here is at *interpretation* or *open* level.

## 0. Setup and notation

Throughout, $(K, v)$ is a field with a **normalized discrete valuation**: $v : K^\times \twoheadrightarrow \mathbb{Z}$ with $v(xy) = v(x) + v(y)$ and $v(x+y) \ge \min(v(x), v(y))$, extended by $v(0) = +\infty$. Write $\mathcal{O} = \{v \ge 0\}$, $\mathfrak{m} = \{v \ge 1\}$, residue field $k = \mathcal{O}/\mathfrak{m}$, and fix the uniformizer $\varpi$ (so $v(\varpi) = 1$). The **angular component** of $x \ne 0$ is $\mathrm{ac}(x) = \overline{x\varpi^{-v(x)}} \in k^\times$ (it depends on the choice of $\varpi$).

$\mathbb{T}$ denotes the min-plus tropical semiring $(\mathbb{Q} \cup \{+\infty\},\ \oplus = \min,\ \otimes = +)$, with $\oplus$-identity $\infty$ and $\otimes$-identity $0$.

Dictionary to the code (all in `/Users/a9lim/Work/ogdoad`):

| math | code |
|---|---|
| $v$, $\varpi$ | `Valued::valuation` (`None` = $\infty$), `Valued::uniformizer` — `src/scalar/valued.rs` |
| $\mathbb{T}$ | `Tropical<MinPlus>` — `src/scalar/tropical.rs` (`Semiring`; ⟦implemented⟧, fuzzed in `tests/tropical_axioms.rs`) |
| $k$, $\mathrm{ac}$ | `ResidueField::Residue`, `residue_unit` — `src/scalar/residue.rs` |
| discretely-valued legs | `Qp<P,K>` ($v(p){=}1$), `Qq<P,N,F>` (unramified, $v(p){=}1$), `Laurent<S,K>` ($v(t){=}1$), `Ramified<S,E>` (renormalized $v(\pi){=}1$, value group $\mathbb{Z}$), `Gauss<S>` ($v(t){=}0$) |
| $\mathbb{F}_q(t)$ per place | `try_valuation_at_ff`, `FFPlace::{Finite(π), Infinite}` — `src/forms/local_global/function_field.rs` |
| Springer buckets | `springer_decompose_local`, `LocalResidueForm { valuation, dim, disc_is_square }`, `parity_layer` — `src/forms/springer/local.rs` |
| Gauss valuation on $K[y]$ | `Poly::min_coeff_valuation` (`src/scalar/poly.rs`), coefficientwise reduction at the minimum (`reduce_poly_at_min` in `src/scalar/functor/gauss.rs`) |

---

## 1. (a) The valuation is the tropicalization

**Lemma J.1 (valuation–tropical dictionary).** ⟦standard math⟧ Define $\tau : K \to \mathbb{T}$ by $\tau(x) = v(x)$ (so $\tau(0) = \infty$). Then:

$$
\begin{aligned}
\text{(i)}\quad & \tau(xy) \;=\; \tau(x) \otimes \tau(y) \quad\text{for all } x, y \in K \text{ (including } 0\text{, by absorption)};\\
\text{(ii)}\quad & \tau(x+y) \,\oplus\, \bigl(\tau(x) \oplus \tau(y)\bigr) \;=\; \tau(x) \oplus \tau(y) \quad\text{i.e.}\quad v(x+y) \ge \min(v(x), v(y));\\
\text{(iii)}\quad & \tau(x+y) \;=\; \tau(x) \oplus \tau(y) \quad\text{whenever } \tau(x) \neq \tau(y);\\
\text{(iv)}\quad & \tau(1) = 0 = 1_{\mathbb{T}}, \qquad \tau(0) = \infty = 0_{\mathbb{T}}.
\end{aligned}
$$

*Proof.* (i), (ii), (iv) restate the valuation axioms in the $(\min,+)$ dictionary; the $\oplus$-internal phrasing of (ii) uses $a \ge b \iff a \oplus b = b$ in $(\mathbb{Q}\cup\{\infty\}, \min)$. For (iii): note first $v(-1) = 0$ (since $2\,v(-1) = v(1) = 0$ in $\mathbb{Z}$), so $v(-y) = v(y)$. Assume WLOG $v(x) < v(y)$, and suppose $v(x+y) > v(x)$. Then $v(x) = v\bigl((x+y) + (-y)\bigr) \ge \min(v(x+y), v(y)) > v(x)$, a contradiction. $\blacksquare$

**Remark J.2 (how "semiring homomorphism" is meant — a non-claim).** $\tau$ is a homomorphism of multiplicative monoids $(K, \cdot, 1, 0) \to (\mathbb{T}, \otimes, 1_\mathbb{T}, 0_\mathbb{T})$ and is **lax** for addition: (ii) with equality (iii) exactly off the *tropical vanishing locus* (the locus where the minimum is attained at least twice — e.g. $v(1 + (-1)) = \infty \ne 0$). No discretely-valued field admits a *strict* additive homomorphism onto $\mathbb{T}$; strict functoriality is restored by replacing $\mathbb{T}$ with the tropical **hyperfield** [Viro 2010], or by taking Lemma J.1(i)–(iii) as the *definition* of a valuation, as in [Maclagan–Sturmfels, Ch. 2]. ROADMAP's slogan "the valuation **is** the tropicalization" has Lemma J.1 as its precise content; prose should not claim strictness.

**Lemma J.3 (graded ring of the valuation filtration).** ⟦standard math⟧ Let $\mathfrak{m}^\lambda = \{x : v(x) \ge \lambda\}$ for $\lambda \in \mathbb{Z}$ (fractional ideals). The associated graded ring of the filtration,
$$
\mathrm{gr}_v(K) \;=\; \bigoplus_{\lambda \in \mathbb{Z}} \mathfrak{m}^{\lambda}/\mathfrak{m}^{\lambda+1},
$$
is, after the choice of $\varpi$, isomorphic to $k[u, u^{-1}]$ ($u = $ class of $\varpi$), and the leading-form map $\sigma : K^\times \to \mathrm{gr}_v(K)$, $\sigma(x) = x \bmod \mathfrak{m}^{v(x)+1}$, is multiplicative, with
$$
\sigma(x) \;=\; \mathrm{ac}(x)\, u^{v(x)}.
$$

*Proof.* Write $x = \varpi^{v(x)} u_x$ with $u_x \in \mathcal{O}^\times$; then $\mathrm{ac}(x) = \bar{u}_x$, each graded piece is a one-dimensional $k$-vector space spanned by $u^\lambda$, and multiplicativity of $\sigma$ is multiplicativity of $v$ and of the residue map on units ($k$ is a field, so there is no cancellation of leading terms). $\blacksquare$

The two lemmas together say: **the valuation/tropical filtration of $K$ has tropical shadow $\tau$ and graded pieces $k \cdot u^\lambda$** — the "residue layers" of part (c).

**Witness tests (a).**
- ⟦implemented⟧ `src/scalar/valued.rs::tests::{uniformizers_have_valuation_one, zero_valuation_is_none}` (J.1(iv) and the $\infty$ convention); `src/scalar/functor/ramified.rs::tests::valuation_is_additive_under_multiplication` (J.1(i) on the ramified leg); `tests/tropical_axioms.rs` ($\mathbb{T}$ is a semiring, both conventions).
- ⟦proposed⟧ `tests/tropicalization.rs`, with the thin adaptor (the Bridge J surface):
  ```rust
  fn trop<K: Valued>(x: &K) -> Tropical<MinPlus> {
      match x.valuation() { Some(v) => Tropical::int(v), None => Tropical::infinity() }
  }
  ```
  proptest over `Qp<5,8>`, `Qq<3,4,2>`, `Laurent<Fp<7>,8>`, `Ramified<Qp<3,8>,2>`, `Gauss<Qp<5,6>>`:
  - `tropicalize_is_multiplicative`: `trop(x.mul(&y)) == trop(&x).mul(&trop(&y))` — exact, zero included;
  - `tropicalize_is_subadditive`: `let s = trop(&x).add(&trop(&y)); trop(&x.add(&y)).add(&s) == s` — the $\oplus$-internal J.1(ii), **truncation-safe**: if a deep cancellation renders the sum as the represented $0$, the left side is $\infty$ and the identity still holds;
  - `tropicalize_equality_off_vanishing_locus`: `if trop(&x) != trop(&y) { trop(&x.add(&y)) == trop(&x).add(&trop(&y)) }` — exact even in the capped models, since the leading term survives truncation.

---

## 2. (b) The Newton-polygon slope theorem

**Definition J.4 (Newton polygon).** For $f = \sum_{i=0}^{n} a_i x^i \in K[x]$ with $a_0 a_n \ne 0$, the **Newton polygon** $\mathrm{NP}(f)$ is the lower boundary of the convex hull of $\{(i, v(a_i)) : a_i \ne 0\} \subset \mathbb{R}^2$, a convex piecewise-linear chain from $(0, v(a_0))$ to $(n, v(a_n))$ with strictly increasing side slopes in $\mathbb{Q}$. (If $a_0 = 0$, factor out $x^m$ first; those $m$ roots are $0$, "valuation $\infty$".)

*Orientation convention — an implementation trap.* With points $(i, v(a_i))$, a side of slope $-\lambda$ corresponds to roots of valuation $+\lambda$. To keep the public surface matching ROADMAP's "slopes are the valuations of the roots", the proposed type should expose `root_valuations() -> Vec<(Rational, u128)>` (negated slopes with horizontal lengths) rather than asking callers to negate; slopes are `Rational` (ratios of `i128`) since root valuations can be fractional even though $\Gamma = \mathbb{Z}$.

**Theorem J.5 (slope theorem).** ⟦standard math: Koblitz, GTM 58, Ch. IV; Neukirch, Ch. II⟧ Let $K$ be **complete** (henselian suffices) with respect to the discrete valuation $v$, let $f \in K[x]$ with $a_0 a_n \neq 0$, let $L$ be a splitting field of $f$, and let $w$ be the unique extension of $v$ to $L$. If $\mathrm{NP}(f)$ has a side of slope $-\lambda$ with horizontal length $\ell$, then $f$ has **exactly $\ell$ roots $r \in L$ (with multiplicity) with $w(r) = \lambda$**, and every root arises this way. In particular $\sum_{\text{sides}} \ell = n$ and the multiset of root valuations is determined by the coefficient valuations alone.

*Proof.* Existence/uniqueness of $w$ on the finite extension $L/K$ is the standard consequence of completeness, $w = \tfrac{1}{[L:K]}\, v \circ N_{L/K}$ [Neukirch, Ch. II]. Normalize $f$ monic (dividing by $a_n$ translates the polygon vertically; slopes and lengths are unchanged). Write $f = \prod_{j=1}^n (x - r_j)$ with $w(r_1) \le \cdots \le w(r_n)$. The coefficients are signed elementary symmetric functions: $a_{n-m} = \pm e_m(r_1, \dots, r_n)$, so by J.1(ii)–(iii) applied in $(L, w)$:
$$
v(a_{n-m}) \;=\; w(e_m) \;\ge\; \min_{|S| = m} \sum_{j \in S} w(r_j) \;=\; \sum_{j \le m} w(r_j),
$$
with **equality whenever the minimizing $m$-subset is unique**, i.e. whenever $w(r_m) < w(r_{m+1})$, and unconditionally at $m = 0$ and $m = n$ (a unique subset each). Let $h(i) := \sum_{j \le n-i} w(r_j)$ for $i = 0, \dots, n$ (height as a function of the point index $i = n - m$). Its successive slopes are $h(i+1) - h(i) = -w(r_{n-i})$, non-decreasing in $i$ because the $w(r_j)$ are sorted — so the graph of $h$ is convex; it lies on or below every point $(i, v(a_i))$; and it touches them at $i \in \{0, n\}$ and at every index where the sorted valuations jump — exactly the vertices of the graph of $h$. Hence the lower convex hull of the points **is** the graph of $h$, and the side of slope $-\lambda$ spans exactly the indices $i$ with $w(r_{n-i}) = \lambda$, of horizontal length $\#\{j : w(r_j) = \lambda\}$. $\blacksquare$

**Lemma J.6 (additivity; Dumas).** ⟦standard math: Dumas 1906⟧ For $f, g \in K[x]$ with nonzero constant terms, the sides of $\mathrm{NP}(fg)$ are obtained by concatenating the sides of $\mathrm{NP}(f)$ and $\mathrm{NP}(g)$ in increasing slope order; per-slope horizontal lengths add.

*Proof (complete case, which is all the project legs need).* Immediate from Theorem J.5: the root multiset of $fg$ in a common splitting field is the union of the two root multisets. (Dumas's original proof is a direct coefficient estimate and needs no completeness.) $\blacksquare$

**Corollary J.7 (Eisenstein).** ⟦standard math: Serre, *Local Fields*, Ch. I⟧ If $f$ is monic of degree $n$ with $v(a_i) \ge 1$ for $i < n$ and $v(a_0) = 1$, then $\mathrm{NP}(f)$ is the single side from $(0,1)$ to $(n,0)$, so every root has valuation $1/n$; $f$ is irreducible, and a root generates a totally ramified extension of degree $n$.

*Proof.* The polygon claim is immediate (all interior points lie on or above the segment). If $h \mid f$ is monic of degree $d$, then $v(h(0)) = \sum_{d \text{ roots}} w(r) = d/n \in \mathbb{Z}$ forces $d \in \{0, n\}$. The value group of $K(r)$ contains $\tfrac{1}{n}\mathbb{Z}$, so $e = n = [K(r):K]$. $\blacksquare$

This is exactly the project's `Ramified<S, E>` ($x^E - \varpi$): its *renormalized* valuation $\min_i\,(E \cdot v_S(a_i) + i)$ rescales the slope-$\tfrac{1}{E}$ root to $v(\pi) = 1$, restoring $\Gamma = \mathbb{Z}$ — which is why the Newton lattice stays integral on that leg.

**Corollary J.8 (unit roots ⟺ flat polygon).** For monic $f \in \mathcal{O}[x]$: all roots of $f$ are units of (the integral closure of $\mathcal{O}$ in) $L$ $\iff$ $\mathrm{NP}(f)$ is the single horizontal side at height $0$ $\iff$ $v(a_0) = 0$ $\iff$ the residue reduction $\bar{f} \in k[x]$ has $\bar{f}(0) \ne 0$.

*Proof.* $v(a_0) = \sum_j w(r_j)$ with every $w(r_j) \ge 0$ (monic, integral coefficients, J.5), so the sum vanishes iff every term does. $\blacksquare$

**Corollary J.9 (per-place polygons over the global $\mathbb{F}_q(t)$).** ⟦standard math: Stichtenoth, GTM 254, Ch. 1⟧ For $f \in \mathbb{F}_q(t)[x]$ and a place $P$ of $\mathbb{F}_q(t)$ (a monic irreducible $\pi$, or $\infty$ with $v_\infty = \deg \mathrm{den} - \deg \mathrm{num}$), the polygon $\mathrm{NP}_P(f)$ computed from the **exact** valuations $v_P(a_i)$ equals the Newton polygon of $f$ over the completion $\mathbb{F}_q(t)_P \cong \mathbb{F}_{q^{\deg P}}((\pi))$, and Theorem J.5 applies there. (The completion at a degree-1 finite place is literally the `Laurent` backend; coefficient valuations are insensitive to completion, so the global leg's polygon is exact with no precision model at all.)

**Witness tests (b)** — all ⟦proposed⟧, on `NewtonPolygon::of(coeffs: &[K]) -> NewtonPolygon` for `K: Valued`:
- `eisenstein_single_slope`: $\mathrm{NP}(x^E - p)$ over `Qp<5,8>` has one side, `root_valuations() == [(1/E, E)]`; cross-check `Ramified::<Qp<5,8>, E>::pi().valuation() == Some(1)` (J.7 ↔ the renormalization).
- `sqrt_p_slope_half`: $\mathrm{NP}(x^2 - p)$ over `Qp<5,8>` gives root valuation $\tfrac12 \notin \mathbb{Z}$; cross-check `Qp::<5,8>::from_i128(5).is_square() == Some(false)` (odd valuation ⇒ nonsquare; `src/scalar/small/analytic.rs`).
- `dumas_additivity`: for $f, g$ with distinct slopes over `Qp`/`Laurent`, per-slope lengths of $\mathrm{NP}(fg)$ are the sums (J.6).
- `flat_polygon_iff_unit_roots`: monic integral $f$; all-zero slopes $\iff$ `a₀.valuation() == Some(0)` $\iff$ the residue reduction has nonzero constant term (J.8, via `ResidueField::residue`).
- `ff_place_polygon_matches_completion`: $f$ over `RationalFunction<Fp<5>>` at the place $t$: polygon from `try_valuation_at_ff` equals the polygon of the coefficientwise image in `Laurent<Fp<5>, K>` (J.9 — the exact-global vs local-model agreement).

---

## 3. (c) Slopes are the Springer residue layers

**Theorem J.10 (Springer).** ⟦standard math: Springer, Indag. Math. 17 (1955); Lam, GSM 67, Ch. VI⟧ Let $K$ be complete discretely valued with $\operatorname{char} k \ne 2$, and fix $\varpi$. Every nondegenerate diagonal form over $K$ is isometric to $q_0 \perp \varpi\, q_1$ with $q_0, q_1$ having unit diagonal entries, and the two **residue homomorphisms** $\partial_0, \partial_1$ (sending $\langle u \rangle \mapsto \langle \bar{u} \rangle$ and $\langle \varpi u \rangle \mapsto \langle \bar{u} \rangle$ respectively) induce a group isomorphism
$$
(\partial_0, \partial_1) : W(K) \;\xrightarrow{\ \sim\ }\; W(k) \oplus W(k),
$$
where $\partial_1$ (not $\partial_0$) depends on the choice of $\varpi$. The two summands are indexed by $\Gamma/2\Gamma = \mathbb{Z}/2$ — they exist *because* the value group is not 2-divisible: $\langle \varpi^2 a \rangle \cong \langle a \rangle$, while $\langle \varpi a \rangle \not\cong \langle a \rangle$ in general.

This is the theorem behind `springer_decompose_local` + `parity_layer` ⟦implemented: `src/forms/springer/local.rs::tests::*`⟧; the code records, per valuation $\lambda$, the layer $(\lambda, \dim, \mathrm{disc\ square\text{-}class})$, and `parity_layer(ε)` is the data of $\partial_\varepsilon$.

**Definition J.11 ($\lambda$-initial form — the graded/tropical piece).** For $\lambda \in \mathbb{Z}$ and $f = \sum a_i x^i \in K[x]$, let
$$
m_\lambda(f) \;=\; \min_i \bigl(v(a_i) + i\lambda\bigr) \;=\; \bigoplus_i \tau(a_i) \otimes \lambda^{\otimes i} \quad(\text{the tropicalized } f \text{ evaluated at } \lambda),
$$
and define the **initial form** $\mathrm{in}_\lambda(f) \in k[y]$ as the coefficientwise reduction of $\varpi^{-m_\lambda(f)} f(\varpi^\lambda y)$ — i.e. substitute $x = \varpi^\lambda y$, then take the Gauss-valuation angular component (in the code: a $\varpi^\lambda$-shift, `Poly::min_coeff_valuation`, and the reduce-at-the-minimum step that `reduce_poly_at_min` in `src/scalar/functor/gauss.rs` already performs — `Gauss<S>` *is* the Gauss valuation this construction lives in). Two standard facts: $\lambda$ is the negative of a slope of $\mathrm{NP}(f)$ iff $\deg \mathrm{in}_\lambda(f) > \operatorname{ord}_y \mathrm{in}_\lambda(f)$ (the minimum is attained at two distinct $i$ — the **tropical-root** criterion [Maclagan–Sturmfels, Ch. 2–3]); and $\mathrm{in}_\lambda(fg) = \mathrm{in}_\lambda(f)\,\mathrm{in}_\lambda(g)$, since the Gauss valuation is a valuation on $K[y]$ and its angular component into the domain $k[y]$ is multiplicative (Lemma J.3 applied to $\mathrm{Gauss}$).

**Proposition J.12 (slope ⟺ residue layer, for diagonal forms).** ⟦standard math; elementary given J.5/J.6 + J.10⟧ Let $q = \langle a_1, \dots, a_n \rangle$ with all $a_i \in K^\times$ (zero entries are the radical, tracked separately as `radical_dim`), and let $f_q(x) = \prod_{i=1}^n (x - a_i)$. Then:

**(i) (the polygon is the bucket shadow).** $\mathrm{NP}(f_q)$ has a side of slope $-\lambda$ and horizontal length $\ell$ $\iff$ $\#\{i : v(a_i) = \lambda\} = \ell$. Hence the side multiset of $\mathrm{NP}(f_q)$ equals the multiset $\{(\texttt{g.valuation}, \texttt{g.dim})\}$ of the Springer decomposition — every Newton slope **is** a residue layer, and conversely.

**(ii) (the initial form is the residue layer's contents).** For each such $\lambda$,
$$
\mathrm{in}_\lambda(f_q) \;=\; c\, \cdot\, y^{\,\#\{i\,:\,v(a_i) > \lambda\}} \prod_{i\,:\,v(a_i) = \lambda} \bigl(y - \mathrm{ac}(a_i)\bigr), \qquad c = \prod_{i\,:\,v(a_i) < \lambda} \bigl(-\mathrm{ac}(a_i)\bigr) \in k^\times,
$$
so the nonzero roots of $\mathrm{in}_\lambda(f_q)$ in $\bar{k}$ are exactly the angular components of the layer, and the layer discriminant is recovered as $\prod_{v(a_i) = \lambda} \mathrm{ac}(a_i)$, whose $k$-square class is `disc_is_square`.

**(iii) (the Witt-level collapse).** If moreover $\operatorname{char} k \ne 2$, the Witt class of $q$ depends only on the layers grouped by $\lambda \bmod 2$: since $\langle a \rangle \cong \langle \varpi^{\,v(a) \bmod 2}\, u_a \rangle$, one gets $\partial_\varepsilon[q] = \bigl[\bigoplus_{v(a_i) \equiv \varepsilon (2)} \langle \mathrm{ac}(a_i) \rangle\bigr] \in W(k)$, and $(\partial_0, \partial_1)$ is Springer's isomorphism. `parity_layer(ε)` computes exactly the data of $\partial_\varepsilon$.

*Proof.* (i): each factor $(x - a_i)$ has the two-point polygon with the single side of slope $-v(a_i)$ and length 1 (using $v(-a_i) = v(a_i)$); apply Lemma J.6. (ii): $\mathrm{in}_\lambda(x - a) = y - \mathrm{ac}(a)$, $y$, or $-\mathrm{ac}(a)$ according as $v(a) = \lambda$, $> \lambda$, $< \lambda$ (compute $m_\lambda = \min(\lambda, v(a))$ directly); multiply, using multiplicativity of $\mathrm{in}_\lambda$ (Definition J.11). (iii): $a = \bigl(\varpi^{\lfloor v(a)/2 \rfloor}\bigr)^2\, \varpi^{\,v(a) \bmod 2}\, u_a$ and, for units, $\langle u \rangle \cong \langle u' \rangle$ over $K$ iff $\bar{u}/\bar{u}'$ is a square in $k$ (Hensel's lemma lifts residue squares when $\operatorname{char} k \ne 2$); then apply Theorem J.10. $\blacksquare$

**Remark J.13 (the forgetful hierarchy — what each level sees).** The data refine strictly:
$$
\underbrace{\mathrm{NP}(f_q)}_{\text{tropical shadow: } (\lambda, \dim) \text{ per layer}} \;\prec\; \underbrace{\{\mathrm{in}_\lambda(f_q)\}_\lambda}_{\text{graded pieces: } + \text{ angular components, hence } \texttt{disc\_is\_square}} \;\prec\; \underbrace{q \text{ itself}}_{\text{the form}}
$$
The polygon is precisely the image of the Springer decomposition under the tropicalization of Lemma J.1 — it sees valuations and dimensions and forgets the residue square classes. This is the exact sense of ROADMAP's "the Springer layers are the graded pieces of the valuation/tropical filtration"; it is the place-axis twin of the games-side identity (thermography in $\mathbb{T}_{\max}$; the sign mirror `MinPlus`↔`MaxPlus` is a convention flip, not a second semiring — `src/scalar/tropical.rs` already enforces the two-type separation).

**Witness tests (c).**
- ⟦implemented⟧ `src/forms/springer/local.rs::tests::{one_engine_decomposes_every_discrete_leg, unramified_qq_reads_extension_residue, residue_char_two_is_rejected_uniformly}` — the bucket engine, the extension-residue square class, and the char-2 boundary.
- ⟦proposed⟧ `polygon_is_the_springer_shadow`: diagonal $\langle a_i \rangle$ over `Qp<5,8>`, `Qq<3,3,2>`, `Laurent<Fp<7>,8>`; build $f_q = \prod (x - a_i)$ via `Poly`; assert the side multiset `{(root_valuation, length)}` equals `{(g.valuation, g.dim)}` from `springer_decompose_local`, and that grouping sides by slope parity reproduces `parity_layer(0)`/`parity_layer(1)` cardinalities (J.12(i), (iii)).
- ⟦proposed⟧ `initial_form_recovers_layer_discriminant`: compute $\mathrm{in}_\lambda(f_q)$ by the shift + `min_coeff_valuation` + reduce-at-min recipe; assert the product of its nonzero roots (equivalently $\pm$ its lowest nonvanishing coefficient ratio) has `is_square_finite::<K::Residue>` equal to the layer's `disc_is_square` (J.12(ii)).
- ⟦proposed⟧ `polygon_outlives_springer`: over `Qp<2,8>` (residue char 2) and `Gauss<Qp<5,6>>` (infinite residue field), `NewtonPolygon::of` succeeds while `springer_decompose_local` returns `None` — J.12(i)–(ii) need no Witt theory; only (iii) does.

---

## 4. Scope boundaries and non-claims

- **Discretely-valued legs only.** The surreal leg has 2-divisible value group: the second Springer layer collapses ($W(\mathrm{No}) = W(\mathbb{R})$, `springer/surreal.rs`) and there is no integer Newton lattice. Polygons over divisible $\Gamma$ are definable but are *not claimed or scheduled* — the same boundary the Springer engine already documents, and itself an instance of the local↔global symmetry.
- **Char-2 residue fields.** J.5/J.6/J.12(i)–(ii) hold for any residue characteristic; J.10/J.12(iii) require $\operatorname{char} k \ne 2$. The char-2 local Witt theory is the separate Aravire–Jacob layer (`springer/char2.rs`) and is outside Bridge J.
- **Precision.** On the capped-relative models (`Qp`/`Qq`/`Laurent`/`Ramified`/`Gauss`), valuations of *represented nonzero* elements are exact, so polygons of represented coefficients are exact; a coefficient whose true valuation exceeds the precision horizon renders as $0$ (vertex absent). J.1(ii) is truncation-safe; equality claims hold off the vanishing locus. The $\mathbb{F}_q(t)$ leg (Corollary J.9) is exact outright.
- **Choice of $\varpi$.** $\mathrm{ac}$, $\mathrm{in}_\lambda$, and $\partial_1$ depend on it; the code pins it to `Valued::uniformizer` via `residue_unit`. $\partial_0$ and the polygon do not.
- **No strictness claim** for "$v$ is a semiring homomorphism" (Remark J.2). No new theorem anywhere in this bridge: J is standard math made computational, the same status as shipped bridges A–I.

## 5. References

- T. A. Springer, *Quadratic forms over fields with a discrete valuation I*, Indag. Math. **17** (1955).
- T. Y. Lam, *Introduction to Quadratic Forms over Fields*, GSM 67, AMS, 2005 — Ch. VI (residue homomorphisms, Springer's theorem).
- N. Koblitz, *p-adic Numbers, p-adic Analysis, and Zeta-Functions*, GTM 58, Springer, 2nd ed. 1984 — Ch. IV (Newton polygons).
- J. Neukirch, *Algebraic Number Theory*, Grundlehren 322, Springer, 1999 — Ch. II (complete/henselian valued fields, unique extension of valuations).
- G. Dumas, *Sur quelques cas d'irréductibilité des polynômes à coefficients rationnels*, J. Math. Pures Appl., 1906 (polygon additivity; the irreducibility criterion).
- J.-P. Serre, *Local Fields*, GTM 67, Springer, 1979 — Ch. I (Eisenstein polynomials, total ramification).
- D. Maclagan, B. Sturmfels, *Introduction to Tropical Geometry*, GSM 161, AMS, 2015 — Ch. 2–3 (valuations as tropicalization; tropical roots/Kapranov in rank 1).
- O. Viro, *Hyperfields for tropical geometry I. Hyperfields and dequantization*, arXiv:1006.3034, 2010 (strict functoriality via the tropical hyperfield).
- H. Stichtenoth, *Algebraic Function Fields and Codes*, GTM 254, Springer, 2009 — Ch. 1 (places of $\mathbb{F}_q(t)$).

## Extraspecial reframing of the Gold-quadric problem

% ============================================================================
% PASTE 1 — preamble additions for writeups/goldarf.tex
% (place with the existing \newtheorem declarations; the file already has
%  observation/proposition/question in plain style and remark in definition
%  style, so add:)
% ============================================================================

\theoremstyle{plain}
\newtheorem{lemma}{Lemma}
\newtheorem{corollary}{Corollary}
\theoremstyle{definition}
\newtheorem{definition}{Definition}

% ============================================================================
% PASTE 2 — new section. Suggested placement: after the subsection "The
% diagonal framing" (end of "What blocks the game semantics"), before
% "Validation status". Uses existing macros \F, \Arf, \Tr and the existing
% label prop:nogo; new citations Quillen71, Gorenstein80, Winter72 are in
% PASTE 3.
% ============================================================================

\section{The extraspecial reframing}

The diagonal-framing question can be restated in group-theoretic terms, and
the restatement sharpens it. The dictionary is classical -- the
characteristic-$2$ Heisenberg/extraspecial picture of quadratic
forms~\cite{Quillen71,Gorenstein80} -- and this section proves the pieces we
use, to keep the note self-contained. Claim levels: Lemmas~\ref{lem:extdict},
\ref{lem:arftypes}, \ref{lem:abelian} and
Proposition~\ref{prop:between} are standard mathematics; the $R_8$
instantiation in Corollary~\ref{cor:kernel} is implemented and tested
(\path{experiments/misere_kernel.py}); reading $E$-equivariance as the right
naturality criterion is interpretation; the existence of a game-native model
of the extension is open.

Throughout, $V=\F_2^n$ written additively, and $Z=\{1,z\}\cong\mathbb{Z}/2$
written multiplicatively. A \emph{quadratic map} $Q\colon V\to\F_2$ is a
function with $Q(0)=0$ whose polarization
$B(x,y)=Q(x+y)+Q(x)+Q(y)$ is bilinear; $B$ is then alternating (hence
symmetric), and $Q$ may have any diagonal $q_i=Q(e_i)$, exactly the
char-$2$ datum the engine keeps separate from $b$.

\subsection{Quadratic forms are central extensions}

\begin{lemma}[extraspecial dictionary]\label{lem:extdict}
Let $1\to Z\to E\xrightarrow{\ \pi\ }V\to 1$ be a central extension of the
elementary abelian $2$-group $V$ by $Z\cong\mathbb{Z}/2$.
\begin{enumerate}[label=(\roman*),leftmargin=*]
\item For $x\in V$ and any lift $\tilde x\in\pi^{-1}(x)$ one has
$\tilde x^2\in Z$, and $\tilde x^2$ is independent of the lift; the
\emph{squaring form} $Q_E\colon V\to\F_2$ defined by
$\tilde x^2=z^{Q_E(x)}$ is well defined, as is the \emph{commutator pairing}
$B_E$ defined by $[\tilde x,\tilde y]=z^{B_E(x,y)}$, which is bilinear and
alternating. They satisfy
\[
(\tilde x\tilde y)^2=\tilde x^2\,\tilde y^2\,[\tilde x,\tilde y],
\qquad\text{equivalently}\qquad
Q_E(x+y)=Q_E(x)+Q_E(y)+B_E(x,y),
\]
so $Q_E$ is a quadratic map with polar form $B_E$.
\item Conversely, every quadratic map $Q$ arises: fixing a basis
$e_1,\dots,e_n$ of $V$, the bilinear $2$-cocycle
\[
\varphi(x,y)=\sum_i x_iy_i\,Q(e_i)+\sum_{i>j}x_iy_j\,B(e_i,e_j)
\]
defines a group $E_Q=Z\times V$ with
$(s,x)(t,y)=(s+t+\varphi(x,y),\,x+y)$, whose squaring form is $Q$ and whose
commutator pairing is the polar form $B$.
\item The assignment $Q\mapsto[E_Q]$ is a bijection from quadratic maps on
$V$ to equivalence classes of central extensions of $V$ by $\mathbb{Z}/2$.
\item $E_Q$ is abelian iff $B=0$ (iff $Q$ is linear); and for $n\geq 1$,
$E_Q$ is \emph{extraspecial} (i.e.\ $Z(E)=[E,E]=\Phi(E)\cong\mathbb{Z}/2$)
iff $B$ is nondegenerate, which forces $n=2r$ even and $|E_Q|=2^{1+2r}$.
\end{enumerate}
\end{lemma}

\begin{proof}
(i) $\pi(\tilde x^2)=2x=0$, so $\tilde x^2\in Z$; replacing $\tilde x$ by
$z\tilde x$ multiplies the square by $z^2=1$. Commutators of lifts lie in
$Z$ because $V$ is abelian, are central, and are unchanged by central
retagging of the lifts; centrality of commutators gives bimultiplicativity,
and $[\tilde x,\tilde x]=1$ gives $B_E(x,x)=0$. For the identity: with
$c=[\tilde y,\tilde x]$ central,
$\tilde x\tilde y\tilde x\tilde y
=\tilde x(\tilde x\tilde y\,c)\tilde y=\tilde x^2\tilde y^2c$,
and $c=c^{-1}=[\tilde x,\tilde y]$ since $Z$ has exponent $2$. As
$\tilde x\tilde y$ lifts $x+y$, the additive identity follows.

(ii) A bilinear map is a $2$-cocycle, so $E_Q$ is a group with center
containing $Z=\{(0,0),(1,0)\}$. Squares: $(s,x)^2=(\varphi(x,x),0)$ and
$\varphi(x,x)=\sum_i x_iQ(e_i)+\sum_{i>j}x_ix_jB(e_i,e_j)=Q(x)$ by the
polarization expansion of $Q$ in the basis. Commutators:
$(s,x)(t,y)(s,x)^{-1}(t,y)^{-1}=(\varphi(x,y)+\varphi(y,x),0)
=(B(x,y),0)$, using that the symmetrization of $\varphi$ is $B$ (the
diagonal terms cancel and $B$ is symmetric).

(iii) Equivalent extensions have cohomologous cocycles, and a coboundary
$\delta\psi(x,y)=\psi(x)+\psi(y)+\psi(x+y)$ has zero diagonal
($\delta\psi(x,x)=\psi(0)=0$), so the squaring form is an invariant of the
class; (ii) gives surjectivity. Injectivity: if two extensions have equal
squaring form $Q$, choose sections and cocycles $\varphi,\varphi'$; both
have diagonal $Q$, and both have symmetrization equal to the polar form of
$Q$, so $d=\varphi+\varphi'$ is a symmetric cocycle with zero diagonal. The
extension $E_d$ is then abelian of exponent $2$, i.e.\ an $\F_2$-vector
space, so $1\to Z\to E_d\to V\to 1$ splits and $d$ is a coboundary; the two
extensions are equivalent.

(iv) Commutators generate $z^{\,\mathrm{im}\,B}$, so $E_Q$ is abelian iff
$B=0$. In general $Z(E_Q)=\pi^{-1}(\operatorname{rad}B)$, since $\tilde x$
is central iff $B(x,\cdot)=0$; thus $Z(E_Q)=Z$ iff $B$ is nondegenerate. If
so, $B\neq0$ gives $[E,E]=Z$, and $\Phi(E)=E^2[E,E]=Z$ because all squares
lie in $Z$ and $E/Z$ is elementary abelian. A nondegenerate alternating
form exists only in even dimension.
\end{proof}

\begin{remark}
This is the characteristic-$2$ Heisenberg picture: $E_Q$ is the Heisenberg
group of $(V,B)$ and $Q$ selects the central extension; cohomologically,
(iii) is the standard identification of $H^2(V;\F_2)$ with quadratic maps
$V\to\F_2$~\cite{Quillen71}. The same finite-quadratic-module data drives
the Weil $S$/$T$ matrices in the library's integral layer. Note that the
dictionary is exactly the engine's discipline in group clothing: the
diagonal $q_i$ (squares) and the off-diagonal $b_{ij}$ (commutators) are
independent central data, and collapsing them is collapsing $E_Q$ onto an
abelian quotient.
\end{remark}

\subsection{Arf classifies the two extraspecial types}

Write $H$ for the hyperbolic plane ($Q(e_1)=Q(e_2)=0$, $B(e_1,e_2)=1$;
$\Arf=0$) and $A$ for the anisotropic plane ($Q\equiv1$ on $V\setminus\{0\}$;
$\Arf=1$), and $\circ$ for the central product (identify the centers).

\begin{lemma}[Arf = the two types]\label{lem:arftypes}
Let $B$ be nondegenerate, $\dim V=2r\geq2$.
\begin{enumerate}[label=(\roman*),leftmargin=*]
\item $E_{Q\perp Q'}\cong E_Q\circ E_{Q'}$.
\item $E_H\cong D_8$ and $E_A\cong Q_8$.
\item Every extraspecial $2$-group of order $2^{1+2r}$ is isomorphic to
$E_Q$ for a nondegenerate $Q$, and
$E_Q\cong E_{Q'}$ iff $(V,Q)\cong(V',Q')$ iff the ranks and Arf invariants
agree~\cite{Dickson01,Arf41}. Hence there are exactly two extraspecial
groups of order $2^{1+2r}$:
\[
E^+_{2^{1+2r}}\;=\;D_8^{\circ r}\ \ (\Arf=0,\ Q\cong H^{\perp r}),
\qquad
E^-_{2^{1+2r}}\;=\;D_8^{\circ(r-1)}\circ Q_8\ \ (\Arf=1,\ Q\cong
H^{\perp(r-1)}\perp A).
\]
In particular $D_8\circ D_8\cong Q_8\circ Q_8$, the group avatar of
$H\perp H\cong A\perp A$ (additivity of $\Arf$).
\item (Census.) With $\varepsilon=(-1)^{\Arf Q}$,
\[
\#\{g\in E_Q: g^2=1\}=2\,\#\{x:Q(x)=0\}=2^{2r}+\varepsilon\,2^{r},
\]
so the element-order census of $E_Q$ is the zero-count bias of
Section~5 read multiplicatively, and it distinguishes the two types.
\end{enumerate}
\end{lemma}

\begin{proof}
(i) The cocycle $\varphi\oplus\varphi'$ on $V\oplus V'$ is bilinear with
diagonal $Q\perp Q'$, and
$(s,(x,x'))\mapsto[((s,x),(0,x'))]$ is an isomorphism onto
$(E_Q\times E_{Q'})/\langle(z,z')\rangle$.

(ii) $E_H$ has order $8$, is nonabelian ($B\neq0$), and has a noncentral
involution (any lift of $e_1$ squares to $z^{Q(e_1)}=1$), so $E_H\cong D_8$;
in $E_A$ every noncentral element squares to $z$ (as $Q\equiv1$ off $0$), so
$z$ is the unique involution and $E_A\cong Q_8$.

(iii) For $E$ extraspecial, $V=E/Z(E)$ is elementary abelian (as
$\Phi(E)=Z(E)$), so Lemma~\ref{lem:extdict} applies to
$1\to Z(E)\to E\to V\to1$ and gives $E\cong E_{Q_E}$ with $B_{Q_E}$
nondegenerate ($Z(E)$ exactly central). If $\psi\colon E_Q\to E_{Q'}$ is an
isomorphism, it preserves the (characteristic) centers, fixes $z\mapsto z'$
(the unique nontrivial central element), and induces an $\F_2$-linear map
$\bar\psi\colon V\to V'$ with
$Q'(\bar\psi x)$ read from $\psi(\tilde x)^2=\psi(\tilde x^2)$ -- an
isometry. Conversely an isometry $g$ transports the cocycle:
$\varphi'\circ(g\times g)$ has diagonal $Q'\circ g=Q$, hence by
Lemma~\ref{lem:extdict}(iii) yields an extension equivalent to $E_Q$, and
$(s,x)\mapsto(s,gx)$ closes the isomorphism. Dickson's classification of
nonsingular quadratic forms over $\F_2$ (two isometry classes per even
rank, separated by $\Arf$) finishes; the central-product normal forms
follow from (i), (ii) and Witt decomposition.

(iv) The two lifts of $x$ both square to $z^{Q(x)}$, so they are
involutions or the identity iff $Q(x)=0$; now apply the zero-count formula
of Section~5. The counts differ for $\varepsilon=\pm1$, so
$E^+\not\cong E^-$.
\end{proof}

\subsection{The abelian obstruction}

\begin{lemma}[abelian obstruction]\label{lem:abelian}
Let $M$ be a commutative monoid, $Z=\{1,z\}\subseteq M$ a two-element
subgroup, $N\subseteq M$ a submonoid containing $Z$, and
$\pi\colon N\twoheadrightarrow V$ a surjective monoid homomorphism onto an
$\F_2$-vector space with $\pi^{-1}(\pi(m))=mZ$ for all $m\in N$. Then the
squaring form $Q\colon V\to\F_2$, defined by $\tilde m^2=z^{Q(x)}$ for any
$\tilde m\in\pi^{-1}(x)$, is well defined and \emph{$\F_2$-linear}; its
polar form vanishes identically. Consequently no commutative monoid
realizes, through its own squaring map, a quadratic form with nonzero --
a fortiori nondegenerate -- characteristic-$2$ polar form. Equivalently, by
Lemma~\ref{lem:extdict}(iv): if $B\neq0$, the extension $E_Q$ is
nonabelian and admits no model $(N,Z,\pi)$ inside any commutative monoid.
\end{lemma}

\begin{proof}
$\pi(m^2)=2\pi(m)=0$, so $m^2\in\pi^{-1}(0)=Z$, and
$(zm)^2=z^2m^2=m^2$, so $Q$ is well defined. Commutativity gives
$(mn)^2=m^2n^2$; since $\tilde m\tilde n$ is a preimage of $x+y$, this is
$Q(x+y)=Q(x)+Q(y)$, i.e.\ $B\equiv0$.
\end{proof}

\begin{corollary}[misère quotients: the kernel shadow is linear]\label{cor:kernel}
Let $\mathcal{Q}$ be a finite misère quotient -- a finite commutative
monoid -- with kernel $K$, the maximal subgroup of $\mathcal{Q}$ (the
mutual-divisibility class of the product of all idempotents)~\cite{PS}.
\begin{enumerate}[label=(\roman*),leftmargin=*]
\item Every configuration $(N,Z,\pi)$ as in Lemma~\ref{lem:abelian} inside
$\mathcal{Q}$ has linear squaring form: the intrinsic multiplication of a
misère quotient cannot supply a quadratic refinement with $B\neq0$, on $K$
or anywhere else in $\mathcal{Q}$.
\item Under the regularity hypothesis of \cite[Thm.~6.4]{PS} (satisfied by
the regular finite quotients arising in practice), $K\cong(\mathbb{Z}/2)^k$
is the group of normal-play Grundy values and the $P$-portion meets $K$ in
the XOR-linear normal-play set.
\item On the smallest wild quotient
$R_8=\langle a,b,c\mid a^2=1,\ b^3=b,\ bc=ab,\ c^2=b^2\rangle$ with
$P=\{a,b^2\}$: the idempotents are $\{1,b^2\}$, the kernel is
$K=\{b^2,b,ab,ab^2\}\cong(\mathbb{Z}/2)^2$ with identity $b^2$,
$P\cap K=\{b^2\}\mapsto\{0\}$ is linear, and the genuinely misère
$P$-element $a$ lies outside $K$. (Implemented and tested:
\path{experiments/misere_kernel.py}.)
\end{enumerate}
So the linear obstruction observed on $R_8$ is forced by commutativity, not
an accident of the example: a nondegenerate polar form is the commutator
pairing of a nonabelian group, and a misère quotient has none to offer.
\end{corollary}

\begin{proof}
(i) is Lemma~\ref{lem:abelian} applied to submonoids of $\mathcal{Q}$,
which are commutative. (ii) is quoted from \cite{PS}. (iii) is a finite
verification.
\end{proof}

\subsection{A Tier-2 naturality screen: $E$-equivariance}

Proposition~\ref{prop:nogo} kills rules that are blind to everything but
the symplectic structure (Tier~1), while per-position evaluator circuits
for $Q_a$ realize the quadric tautologically (Tier~3). The extraspecial
dictionary suggests a criterion strictly between the two: the rule should
see the extension $E_Q$ -- which carries the diagonal data $q_i$
structurally, as squares -- but only up to its automorphisms, so that no
basis, field structure, or evaluation circuit can be smuggled in.

\begin{definition}[uniform rules; the three tiers]\label{def:tier2}
Let $Q$ be nondegenerate on $V=\F_2^{2r}$ with polar form $B$, and
$E=E_Q$, $\pi\colon E\to V$, $\Sigma=\{x\in V: Q(x)=0\}$.
\begin{enumerate}[label=(\alph*),leftmargin=*]
\item A \emph{uniform rule on a finite set $X$} is a single move relation
$M\subseteq X\times X$, fixed once for all positions; for normal play we
require the move digraph acyclic, and outcomes are computed by the usual
retrograde recursion. (For loopy or misère readings, replace the outcome
recursion; the equivariance constraint below applies verbatim, since
outcome classes are invariants of digraph automorphisms.)
\item A uniform rule on $E$ is \emph{$E$-equivariant} if every
$\alpha\in\operatorname{Aut}(E)$ is an automorphism of the move digraph.
Its \emph{shadow} is $\pi(P)\subseteq V$, where $P$ is its $P$-set; it
\emph{realizes} $T\subseteq V$ if $\pi(P)=T$.
\item \emph{Tier 1} ($\operatorname{Sp}(B)$-blind): a uniform rule on $V$
invariant under $\operatorname{Sp}(B)$, as in
Proposition~\ref{prop:nogo}. \emph{Tier 3} ($Q_a$-evaluation): a family of
games $\{\Gamma_x\}_{x\in V}$ whose structure varies with the designated
input $x$ (e.g.\ the trace-circuit evaluator); this is not a uniform rule.
\emph{Tier 2 screen}: a route to the Gold quadric is Tier-2 natural
\emph{only if} its positions and moves lift to an $E$-equivariant uniform
rule on $E_{Q_a}$ -- with $Q_a$ restricted to its nonsingular core when the
form is degenerate, matching the classifier's radical discipline -- that
realizes $\{Q_a=0\}$.
\end{enumerate}
\end{definition}

\begin{proposition}[the screen sits strictly between the solved tiers]\label{prop:between}
Let $Q$ be nondegenerate on $V=\F_2^{2r}$, $r\geq2$, and $E=E_Q$.
\begin{enumerate}[label=(\roman*),leftmargin=*]
\item The image of $\operatorname{Aut}(E)\to GL(V)$ is exactly the
orthogonal group $O(Q)$, with kernel the central twists
$g\mapsto g\,z^{\ell(\pi g)}$, $\ell\in V^{*}$, which equal
$\operatorname{Inn}(E)\cong V$ because $B$ is nondegenerate; thus
$\operatorname{Out}(E)\cong O(Q)$ (cf.~\cite{Winter72}). The
$\operatorname{Aut}(E)$-orbits on $E$ are
\[
\{1\},\qquad \{z\},\qquad \pi^{-1}(\Sigma\setminus\{0\}),\qquad
\pi^{-1}(V\setminus\Sigma).
\]
\item Hence the $P$-set of any $E$-equivariant uniform rule is a union of
these four orbits, and its shadow is one of the eight unions of $\{0\}$,
$\Sigma\setminus\{0\}$, $V\setminus\Sigma$. The quadric $\Sigma$ is among
them: the Tier-1 exclusion does not apply.
\item (Tier 1 $\subsetneq$ Tier 2.) Every $\operatorname{Sp}(B)$-invariant
uniform rule on $V$ pulls back along $\pi$ to an $E$-equivariant rule with
$P$-set $\pi^{-1}$ of the original, whose shadow is therefore one of
$\varnothing$, $\{0\}$, $V\setminus\{0\}$, $V$ -- never the quadric
(Proposition~\ref{prop:nogo} in this language). The containment is strict:
the \emph{squaring rule}
\[
g\to h\ \text{legal}\iff g^2=z\ \text{and}\ h^2=1
\]
(``move from any order-$4$ position to any position of order at most
$2$'') is $E$-equivariant and acyclic, and its $P$-set is the involution
locus $\{g:g^2=1\}=\pi^{-1}(\Sigma)$, with shadow exactly $\Sigma$.
\item (Tier 2 $\subsetneq$ Tier 3.) Up to the isometry induced by any group
isomorphism, the family of $E$-equivariant rules and their shadows depends
only on $(r,\Arf Q)$ (Lemma~\ref{lem:arftypes}(iii)). In particular the
screen cannot separate Gold forms of equal core rank and Arf invariant, and
an $E$-equivariant rule has no access to $m$, the exponent $a$, the field
multiplication, the coordinate frame $q_i=Q_a(e_i)$, or any evaluation
circuit. Conversely, arbitrary subsets of $V$ -- all realizable by ad hoc
acyclic lookup games and by Tier-3 evaluators (Section~7.1) -- are not
shadows of $E$-equivariant rules once they leave the eight-set list of
(ii).
\end{enumerate}
\end{proposition}

\begin{proof}
(i) Automorphisms preserve $Z=Z(E)$ and fix $z$
($\operatorname{Aut}(\mathbb{Z}/2)=1$), so the induced map preserves the
squaring form: $Q(\bar\alpha x)$ is read from
$\alpha(\tilde x)^2=\alpha(\tilde x^2)=z^{Q(x)}$; the image lies in
$O(Q)$. Surjectivity: an isometry $g$ transports the standard cocycle as in
the proof of Lemma~\ref{lem:arftypes}(iii), producing an automorphism of
$E_Q$ inducing $g$. The kernel consists of maps $g\mapsto g\,z^{\ell(\pi g)}$
with $\ell$ linear (multiplicativity forces additivity of $\ell$); inner
automorphisms are exactly the twists $\ell=B(\bar g,\cdot)$, and
nondegeneracy makes $\bar g\mapsto B(\bar g,\cdot)$ onto $V^{*}$. Orbits:
$O(Q)$ is transitive on $\Sigma\setminus\{0\}$ and on $V\setminus\Sigma$ by
Witt's extension theorem in its characteristic-$2$ quadratic-space form
\cite{Taylor} (both sets are nonempty for $r\geq2$); the two lifts of any
$x\neq0$ are conjugate under conjugation by $\tilde g$ with $B(\bar g,x)=1$;
and $1$, $z$ are fixed by everything.

(ii) Digraph automorphisms preserve the retrograde outcome recursion, so
$P$ is $\operatorname{Aut}(E)$-invariant; apply (i) and project.

(iii) The pullback $M=\{(g,h):(\pi g,\pi h)\in R\}$ of an
$\operatorname{Sp}(B)$-invariant rule $R$ is preserved by
$\operatorname{Aut}(E)$, since the induced action on $V$ lies in
$O(Q)\subseteq\operatorname{Sp}(B)$; it is acyclic when $R$ is, and an easy
induction along the acyclic rank gives
$P(M)=\pi^{-1}(P(R))$. By the transitivity of
$\operatorname{Sp}(B)$ on $V\setminus\{0\}$ for $r\geq2$, $P(R)$ is a union
of $\{0\}$ and $V\setminus\{0\}$. For the squaring rule: automorphisms
preserve element orders, so the rule is $E$-equivariant; positions with
$g^2=1$ are terminal, hence $P$; positions with $g^2=z$ have a move (to
$1$), hence $N$. Its shadow $\Sigma$ satisfies
$1<|\Sigma|=2^{2r-1}+\varepsilon2^{r-1}<2^{2r}-1$ for $r\geq2$, so it is
none of the four Tier-1 shadows.

(iv) If $\psi\colon E_Q\to E_{Q'}$ is an isomorphism (it exists iff ranks
and Arf agree, Lemma~\ref{lem:arftypes}(iii)), then
$M\mapsto(\psi\times\psi)(M)$ is a bijection between $E$-equivariant rules
carrying $P$-sets to $P$-sets and shadows to shadows through the induced
isometry $\bar\psi$; in particular it carries rules realizing $\{Q=0\}$ to
rules realizing $\{Q'=0\}$. The shadow constraint of (ii) excludes all but
eight subsets, while Tier-3 constructions realize any subset; the
containment is proper.
\end{proof}

\begin{remark}[what the screen does and does not settle]
Equivariance is a screen, not a construction. Over the \emph{abstract}
group $E$ the screen is satisfiable -- Proposition~\ref{prop:between}(iii)
exhibits the squaring rule, which is nothing but the extraspecial squaring
map of Lemma~\ref{lem:extdict} read as a one-move relation. That rule
consults only the multiplication of $E$, so it is a $Q$-evaluator exactly
insofar as $E$ itself is fed in by hand. The reframing therefore relocates
the open problem: instead of \emph{find a play rule that computes $Q_a$},
it asks \emph{exhibit a game-native model of the extension $E_{Q_a}$} --
positions built from game values, multiplication from game constructions --
after which the rule layer is canonical. Lemma~\ref{lem:abelian} is the
sharp constraint on that search: no commutative value-world, in particular
no misère-quotient kernel (Corollary~\ref{cor:kernel}), can host $E$ once
$B\neq0$. The noncommutativity must enter from a structurally available
asymmetry, and the one that normal, misère, and partizan play all possess
-- and that the symmetric polar form $B$ discards -- is the
first-/second-player asymmetry of the move relation. Testing the existing
routes (kernel, loopy Draw-set, misère quotient) against the screen --
full $E$-symmetry versus only $\operatorname{Sp}(B)$ or an abelian
quotient -- is recorded as a progress target in \texttt{OPEN.md}. Claim
levels: the proposition is standard mathematics; treating
$E$-equivariance as \emph{the} naturality criterion is interpretation;
the existence of a game-native $E$ is open. (For $r=1$ the screen degenerates:
on the anisotropic plane $\Sigma=\{0\}$ and $O(Q)=\operatorname{Sp}(B)$,
so the statement is vacuous there; the Gold targets of interest have
$r\geq2$ cores.)
\end{remark}

\begin{question}
Does any game-native construction realize the extraspecial extension
$E_{Q_a}$ -- equivalently, produce the diagonal data $q_i=Q_a(e_i)$ as
squares in a noncommutative game-built structure -- without an evaluation
circuit for $Q_a$? By Lemma~\ref{lem:abelian} the source cannot be any
commutative game-value monoid.
\end{question}

% ============================================================================
% PASTE 3 — bibliography additions (alphabetical slots indicated):
%   after \bibitem{Gold68}:
% ============================================================================

\bibitem{Gorenstein80} D.~Gorenstein, \emph{Finite Groups}, 2nd ed.,
Chelsea Publishing, New York, 1980.

%   after \bibitem{PS}:

\bibitem{Quillen71} D.~Quillen, \emph{The mod~$2$ cohomology rings of
extra-special $2$-groups and the spectra of quadratic forms}, Math.\ Ann.\
\textbf{194} (1971), 197--212.

%   after \bibitem{Taylor}:

\bibitem{Winter72} D.~L.~Winter, \emph{The automorphism group of an
extraspecial $p$-group}, Rocky Mountain J.\ Math.\ \textbf{2} (1972),
159--168.

% ============================================================================
% Notes for the integrator (not for the .tex):
% - Sources verified against repo surfaces: OPEN.md (extraspecial paragraph,
%   tier dichotomy, progress targets), experiments/misere_kernel.py (R8
%   presentation, idempotents {1,b^2}, kernel {b^2,b,ab,ab^2}, P∩K={b^2},
%   a ∉ K, PS Thm 6.4 + regularity caveat), writeups/goldarf.tex
%   (prop:nogo label, \F/\Arf/\Tr macros, existing Taylor/PS/Arf41/Dickson01
%   bibitems, zero-count formula in §5, bench list in §7.1).
% - "Section~5"/"Section~7.1" in the prose refer to "Arf as conditional
%   win-bias" and "Test benches" in the current draft; renumber if sections
%   move.
% - Claim-level discipline: all lemmas/propositions are standard math with
%   self-contained proofs; Corollary (iii) is implemented-and-tested;
%   the criterion-as-naturality reading is flagged interpretation; existence
%   of a game-native E is flagged open. No new mathematical claims beyond
%   standard extraspecial-group theory are introduced.
% ============================================================================

## Bridge K — cyclic algebras and the full Q/Z Brauer invariant

# Bridge K — cyclic algebras and the full $\mathbb{Q}/\mathbb{Z}$ local Brauer invariant

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

ROADMAP's one-line gloss ("the reduced norm form of $(\chi_\sigma,a)$ *is* the twisted trace form") is loose; the honest statements are:

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
