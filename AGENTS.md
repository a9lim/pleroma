# AGENTS.md — pleroma

Working notes for agents editing this repo. Global rules still apply.

## What this is

Clifford algebras (with nilpotents) over the field-like subclasses of Conway's
combinatorial games. Games under disjunctive sum are an abelian group, **not a
ring** — Conway multiplication is only a congruence on the numbers. A Clifford
algebra needs a commutative scalar ring, so this only lives on the three
field-like cores of game-world, and each is a backend:

- **nimbers** `On₂` — algebraically closed, characteristic **2**. The only
  backend where Clifford gets a genuinely new flavour (alternating polar form,
  `q ≠ b`).
- **surreals** `No` — real-closed, char 0. Cl(p,q) exactly as over ℝ, but metric
  entries may be infinite/infinitesimal.
- **surcomplex** `No[i]` — algebraically closed, char 0.

A pure Rust math core, generic over a `Scalar` trait, with PyO3 per-backend
bindings on top. "With nilpotents" = the quadratic form may be degenerate
(`q[i]=0` ⇒ `eᵢ²=0`); all-zero `q` is the exterior/Grassmann algebra.

## Layout

```
src/
  lib.rs          # crate root: the four pillars + the (feature-gated) py module.
                # Each pillar's mod.rs re-exports its children flat, so public
                # paths stay shallow (scalar::Nimber, clifford::sandwich, …).

  scalar/         # PILLAR — the commutative coefficient worlds (generic Scalar)
    mod.rs        # just the Scalar trait (add/neg/mul/zero/one/is_zero/inv/
                  # characteristic) + the flat re-export hub for the worlds below.
    rational.rs   # exact ℚ over i128, NOT a game backend — the char-0 scalar that
                  # validates the geometric product against the known Cl(p,q)
                  # classification before we trust the exotic backends.
    integer.rs    # exact ℤ — the coefficient ring for the exterior algebra of the
                  # game group (games/partizan.rs): games are a ℤ-module, not a
                  # ring, so Λ over ℤ is the structure that lives on all of
                  # game-world. Only ±1 invertible (Grassmann never calls inv).
    nimber.rs     # On₂ in u128 (= F_{2^128}): nim_add = XOR; nim_mul via Fermat-
                  # power recursion, memoised on 2^i ⊗ 2^j. Also nim_square /
                  # nim_sqrt (Frobenius & its inverse), nim_trace, and the
                  # Artin–Schreier solver (y²+y=c, solvable ⇔ Tr(c)=0).
    surreal.rs    # Conway normal form: Vec<(exponent: Surreal, coeff: Rational)>
                  # with recursive exponents. Hahn arithmetic: ω^a·ω^b = ω^{a+b}.
                  # Plus the {L|R}/simplicity bridge (dyadic case): as_rational/
                  # as_dyadic/is_dyadic/dyadic_birthday + simplest_above/_below/
                  # simplest_between (the shallowest surreal-tree node in (a,b)).
    surcomplex.rs # Surcomplex<S> = adjoin i over any backend.
    omnific.rs    # the omnific integers Oz: Omnific(Surreal) newtype, a transfinite
                  # commutative RING (not field). Surreal mirror of Integer; the
                  # exterior algebra with ω-scale coefficients.
    onag.rs       # transfinite (ordinal) nimbers: Ordinal in CNF (mirror of
                  # scalar/surreal.rs). nim-add COMPLETE (coeff XOR); nim-mul
                  # COMPLETE across φ_{ω+1} (all ordinals < ω³ Cantor) via DiMuro
                  # Lemma 1.1: poly mult in (finite nimbers)[ω] mod ω³−2.
                  # ω⊗ω⊗ω=2; F₄(ω)≅F₆₄ verified. Above ω³ staged (Lenstra tower).
    fp.rs         # Fp<const P>: the prime field F_P (odd characteristic), a
                  # comparison backend completing the char trichotomy. Genuine neg.
    fpn.rs        # Fpn<const P, const N>: finite extension fields F_{p^N} via a
                  # (P,N)-keyed irreducible reduction poly. Completes the odd-char
                  # tower AND the char-2 odd-DEGREE fields nimbers can't reach (F_8).
                  # Schoolbook mul+reduce, Fermat inv, is_square (Euler/Frobenius).
    zp.rs         # Zp<const P, const K>: the p-adic integers Z_p to precision k
                  # (= Z/p^k). A LOCAL RING (p a non-unit), not Q_p — char()=0,
                  # inv = Omnific pattern (units only). Cl over it is non-semisimple.
    wittvec.rs    # WittVec<const P, const N, const F>: Witt vectors W_N(F_q), as the
                  # truncated unramified ring (Z/p^N)[t]/(f̃) (NOT the forms Witt
                  # group). Witt/Teichmüller coords + carry-formula oracle on top.

  clifford/       # PILLAR — the multivector engine + GA layer (generic Scalar)
    mod.rs        # thin hub: re-exports engine + versor + the structured-algebra
                  # modules flat (clifford::Metric, clifford::sandwich, …).
    engine.rs     # Metric { q, b, a } + CliffordAlgebra<S> + Multivector<S>. The
                  # associative-algebra core: geom_product_blades (general-bilinear
                  # Chevalley product; reduce_word is a #[cfg(test)] oracle it is
                  # cross-validated against), blade arithmetic, grade projection.
                  # (bits/grade are pub; q_val/has_upper are pub(crate).)
    versor.rs     # the GA layer on top: versor_inverse, sandwich, twisted_sandwich
                  # (Pin action), reflect, left/right_contract, dual/undual,
                  # grade_involution, norm2, even_part / even_subalgebra. Plus the
                  # product/involution suite: clifford_conjugate, scalar_product
                  # ⟨ab⟩₀, commutator/anticommutator (½-free, char-faithful), and
                  # the regressive meet a∨b (intersection dual to wedge).
    blade.rs      # blade analysis: blade_subspace {x:x∧A=0}, is_blade (dim⟨A⟩=grade),
                  # factor_blade (decompose a blade into grade-1 vectors). Nullspace
                  # over the field; char-faithful.
    outermorphism.rs # lift a grade-1 LinearMap<S> to all grades (f(a∧b)=f(a)∧f(b));
                  # determinant as the pseudoscalar action f(I)=det·I; compose,
                  # inverse_outermorphism. Plus the char poly via exterior powers:
                  # exterior_power_trace (cₖ=tr Λᵏf), trace, char_poly. Char-faithful
                  # (the char-2 determinant/permanent too).
    hopf.rs       # the exterior Hopf algebra: unshuffle coproduct (sign read off
                  # wedge), counit, antipode = grade involution (NOT reversion-
                  # twist). Hopf axioms tested over Rational AND Nimber.
    cga.rs        # conformal (Cl(n+1,1) null basis: up/down/inner/sphere/plane/
                  # meet) + projective (pga = Cl(n,0,1), exp_nilpotent terminating
                  # motor exp) GA. Char-0 (needs ½); surreal ∞/ε radii are exact.
    spinor.rs     # concrete minimal left ideals Cl·f from a primitive idempotent
                  # ∏½(1+w); basis + gen_matrices realizing M_d(K) on column
                  # spinors. Ideal dim matches classify; Clifford relations hold.
    spinor_norm.rs # the spinor norm N:O(Q)→F*/F*² (= norm2 mod squares) + the
                  # generic versor_grade_parity (Dickson; char2::dickson_of_versor
                  # delegates here) + classify_versor. Char-2 codomain is F/℘(F).

  forms/          # PILLAR — quadratic forms & invariants, by the char trichotomy
    mod.rs        # re-exports the legs + diagonalize/equivalence + witt/witt_ring
                  # + brauer_wall + padic + springer.
    diagonalize.rs # congruence diagonalization (char ≠ 2): gram, diagonalize,
                  # as_diagonal. Returns None in char 2 (nonsingular forms aren't
                  # diagonalizable there — use char2.rs's symplectic Arf reduction).
                  # This is what lets char0/oddchar classify ARBITRARY metrics.
    equivalence.rs # isometry (per backend, via the complete invariant) +
                  # Witt decomposition (k·H ⊥ anisotropic kernel) over ℝ and F_q.
    char0.rs      # (was classify.rs) the char-0 Clifford classifier: Cl(p,q) →
                  # matrix algebra over ℝ/ℂ/ℍ via the 8-fold table (real-closed
                  # surreal/rational) and the 2-fold table (surcomplex). Non-diagonal
                  # metrics are diagonalized first (signature is pub(crate)).
    oddchar.rs    # (was disc.rs) the odd-char classifier: discriminant + is_square
                  # (Euler) + hilbert_symbol/hasse_invariant (≡ +1 over finite
                  # fields) + classify_oddchar + oddchar_witt. dim+disc complete.
                  # Non-diagonal metrics are diagonalized first (as_diagonal).
    char2.rs      # (was arf.rs) the Arf invariant (char-2 classifier): arf_f2 (F₂
                  # bitmask) + arf_nimber (any nim-field, symplectic reduction + the
                  # field trace). Also the Dickson invariant (dickson_matrix =
                  # rank(g−I) mod 2, ker = SO; dickson_of_versor) and
                  # fit_f2_quadratic (is a set a quadric? its Arf?).
    witt.rs       # WittClass: the Witt group W_q(F) ≅ ℤ/2 of a finite nim-field,
                  # Arf-classified. Plus WittClassG: the Char0/OddChar/Char2
                  # trichotomy enum (odd-char is order-4: ℤ/4 or ℤ/2×ℤ/2) with the
                  # ring `mul` (Char2 panics — W_q is a module, not a ring).
    witt_ring.rs  # the Witt RING: tensor_form, Pfister forms, fundamental ideal Iⁿ,
                  # and the eₙ staircase (e0=dim, e1=disc, e2=Hasse — reused from
                  # oddchar). Stabilization per field (I²=0 over F_q; infinite ℝ tower
                  # via e_real). DON'T claim Arf=e2 (char-2 indexing is Kato's, pinned).
    brauer_wall.rs # the Brauer–Wall group BW(F): bw_class_real (= char0's Bott index
                  # (q−p) mod 8 ⇒ BW(ℝ)=ℤ/8), bw_class_complex (ℤ/2), bw_class_oddchar
                  # (order-4 ≅ W(F_q), DISCOVERED not asserted). Law = graded_tensor.
    padic.rs      # the GENUINE Hilbert symbol over Q_p (odd-p + p=2 mod-8) — nontrivial
                  # unlike oddchar's +1 — + Hasse–Minkowski is_isotropic_q over Q.
                  # Oracle: Hilbert reciprocity ∏_v=+1. Where Hasse does real work.
    springer.rs   # non-Archimedean Springer decomposition over the surreals: a
                  # diagonal form's ω-adic valuation filtration into residue ℝ-forms.
                  # Honest: value group 2-divisible ⇒ W(No)=W(ℝ)=ℤ; the filtration
                  # itself is the novelty.

  games/          # PILLAR — combinatorial game theory (mostly Scalar-free)
    mod.rs        # re-exports the modules below flat.
    coin_turning.rs # (was games.rs) nim_mul_mex: nim-mult as Conway's Turning-
                  # Corners mex recurrence (== algebraic nim_mul). Plus general 1-D
                  # coin-turning (grundy_1d) and the 2-D Tartan product
                  # (tartan_grundy), with the Tartan/Product theorem verified.
    grundy.rs     # general Sprague–Grundy (normal-play impartial center): mex,
                  # grundy_graph (DAG; None on a cycle), closure-based grundy.
                  # P-position ⟺ g=0; SG theorem g(G+H)=g(G)⊕g(H) pinned vs Bouton.
    kernel.rs     # normal-play Win/Loss/Draw outcomes of any finite game graph
                  # (retrograde analysis); P-positions = Loss. The interactive route.
                  # Plus scoring_values: the Milnor minimax interval (left,right) on a
                  # DAG — the integer-valued scoring knob for the open question.
    misere.rs     # misère-play outcomes (misere_is_n/_is_p) for any finite
                  # impartial game; misère Nim vs Bouton; the bounded
                  # indistinguishability quotient (misere_quotient); octal games
                  # (octal_moves, octal_misere_quotient). The non-linear route.
    partizan.rs   # short partizan games (sum/neg/order/birthday/is_number) + the
                  # CANONICAL FORM (dominated/reversible reduction; structural_string
                  # vs canonical_string — the latter canonicalizes, a value key) +
                  # the game↔surreal bridge (number_value / from_surreal, numbers
                  # only) + the exterior algebra of the GAME group: Λ over ℤ on game
                  # generators (living on all of game-world, incl. non-numbers ⋆/↑).
                  # NB: distinct from coin_turning.rs — that is coin-turning.

  py/             # PyO3 bindings (feature = "python"), split per pillar
    mod.rs        # the #[pymodule]; chains each submodule's pub(crate) register().
    scalars.rs    # the scalar pyclasses + constructors + nim-field fns; the
                  # parse_*/wrap_* hooks the backend! macro threads in. Surreal also
                  # exposes the simplicity bridge (simplest_between/above/below,
                  # dyadic_birthday, is_dyadic).
    engine.rs     # the backend! macro → <World>Algebra + <World>MV pairs (incl.
                  # the Integer backend) + conformal GA (Cga). MV methods include the
                  # Arc-C suite (clifford_conjugate, scalar_product, commutator,
                  # anticommutator, undual, meet, is_blade, factor_blade); algebra
                  # methods add trace / char_poly alongside determinant.
    forms.rs      # classify / witt / dickson / odd-char / springer bindings, plus
                  # witt_decompose_real / witt_decompose_oddchar / is_isometric_oddchar.
    games.rs      # Game (incl. canonical/is_canonical/canonical_string/
                  # number_value/from_surreal) / GameExterior + nim_mul_mex +
                  # grundy_graph.

examples/tour.rs   # cargo run --example tour   (Rust-only demo)
examples/misere_quotient.rs    # misère quotients + the quadric test on P-sets
examples/interactive_kernel.rs # B-coupled interactive games vs {Q=0}
examples/octal_hunt.rs         # sweep octal games for a (ℤ/2)^k quadric P-set
                               # (cargo run --release --example octal_hunt)
demo.py            # the same tour from Python
experiments/       # research probes ON TOP of the shipped lib: Arf of Gold
                   # forms, the game-built synthesis, the Arf win-bias,
                   # artin_arf (the trace ↔ Arf unification),
                   # open_question_probe (the polar-form obstruction),
                   # tartan_bilinear (B realized by Turning-Corners), and
                   # framing_obstruction (the Sp(B) no-go + the diagonal-framing
                   # ladder for the open question). See NOTES.md.
```

The math thread (Arf↔Clifford, the games bridge, the char-0/char-2 classifier
symmetry, the Artin–Schreier ↔ Arf unification, the open play-semantics
question) is written up in `NOTES.md` — read it before touching `forms/char2.rs`,
`forms/char0.rs`, `games/coin_turning.rs`, `games/kernel.rs`, `games/misere.rs`,
`forms/witt.rs`, `experiments/`, or the `misere_quotient` / `interactive_kernel`
examples.

## Commands

```sh
cargo test                                    # the math core (pure Rust, no Python)
cargo run --example tour                      # Rust demo
python3 -m venv .venv && .venv/bin/pip install maturin
VIRTUAL_ENV=.venv .venv/bin/maturin develop   # build + install the abi3 extension
.venv/bin/python demo.py
```

`maturin develop` needs `VIRTUAL_ENV` set (or a `.venv` in cwd) and `cargo` on
PATH (`. "$HOME/.cargo/env"`).

## Hard rules

1. **The math core is generic over `Scalar` and pure Rust.** PyO3 lives behind
   the `python` feature (`pyo3` is an optional dep; `extension-module` only
   enabled there). This is what keeps `cargo test` from linking libpython.
   Never `use pyo3` outside the `py/` module; never make it non-optional.

2. **The metric carries `q` and `b` independently — do not collapse them.**
   `q[i] = eᵢ²` (quadratic form); `b[(i,j)] = {eᵢ,eⱼ}` (polar/anticommutator,
   i<j). In char ≠ 2 they're linked; in char 2 they are NOT — `b` is alternating
   (`b(i,i)=0`) yet `q[i]` can be nonzero. Collapsing to one symmetric bilinear
   form silently makes every char-2 algebra commutative and throws away the
   entire point of the nimber backend. There is now a THIRD, *optional* field
   `a[(i,j)]` (i<j): the in-order / asymmetric contraction that lifts the engine
   to a general (non-symmetric) bilinear form `B` — `e_i e_j = e_i∧e_j + a_{ij}`
   for i<j; `b` stays the (symmetric) anticommutator regardless. `a` empty ⇒ the
   ordinary Clifford algebra. Build metrics with `Metric::new(q, b)` (a empty),
   `Metric::diagonal`, `Metric::grassmann`, or `Metric::general(q, b, a)` rather
   than the bare struct literal, so the `a` field is handled for you (`a` is keyed
   i<j only).

3. **Signs go through the scalar's own `neg()`, never a literal `-1` or a
   `characteristic()` branch.** The product (`geom_product_blades`, and the
   `#[cfg(test)]` oracle `reduce_word`) emits `S::one().neg()` from the wedge
   antisymmetry. For nimbers `neg` is identity, so `-1 = 1` and char-2
   sign-vanishing falls out for free. Hardcoding signs breaks char 2.

4. **Surreal arithmetic recurses only on exponents.** Every op (add/mul/cmp) on
   a `Surreal` recurses into its *exponents*, which are strictly simpler (lower
   depth) than the number itself. That is the entire termination argument. Never
   write a recursion that calls back on the number.

5. **Per-backend, no mixing.** Each Python backend monomorphises the generic
   engine to one concrete scalar type. Mixing scalar worlds in one algebra is
   impossible by construction (raises `TypeError`) and that's intended — do not
   add a runtime-tagged "any scalar" path.

6. **Verify, don't claim.** Engine + every backend have `cargo test` checks. The
   `associativity_*` tests (incl. `associativity_general_bilinear_form`) are the
   ones that actually catch product bugs, and `general_product_reproduces_reduce_word_when_a_empty`
   pins the general engine to the independent oracle — add a test before trusting
   a new operation. The char-0 classifier is checked against the known low-dim
   table + a dimension-consistency sweep; Dickson against known O(Q) elements;
   the Artin–Schreier solver against the trace obstruction exhaustively on F₁₆.

## Style

- Rust 2021, `cargo fmt` clean, no warnings. License: see `LICENSE`.
- Display is deliberate and should stay readable: blades render `e0e1`;
  coefficients `1`/`-1` are elided; nimbers print `*n`; surreals print CNF
  (`3ω^2 - ω + 5`, `ω^(ω)`, `ω^-1`). Keep `display()` / `Debug` matching this.
- Python operators: `*` geometric, `^` wedge, `<<`/`>>` left/right contraction,
  `~` reverse, `/` divide (scalar or versor), `**` power, `+`/`-`, `==`.

## Testing

`cargo test` is the source of truth and needs no Python. The Python layer is
smoke-tested via `demo.py`. After touching `clifford/` or `scalar/surreal.rs`, run
`cargo test` **and** rebuild + run `demo.py` — display changes don't surface in
`cargo test`.

## Things that look like bugs but are not

- **Char-2 Clifford over an orthogonal basis is commutative.** `e0*e1 == e1*e0`
  when `b` is empty and the scalar is a nimber. Correct: `{e0,e1}=2B=0` and
  `-1=1`. Set an off-diagonal `b[(i,j)]` to get non-commutativity.
- **Surcomplex over nimbers is degenerate.** `i²=1`, `(1+i)²=0`, not a field.
  That's the theorem — On₂ is already algebraically closed, so `i` adjoins
  nothing. Surcomplex is only meaningful over the surreals.
- **Surreal coefficients are ℚ, not ℝ** — the honest finite truncation of true
  CNF. Exponents *are* fully recursive surreals. Don't "fix" this expecting
  irrational coefficients.
- **`Surreal::inv` returns `None` for any non-monomial.** `1/(ω+1)` is an
  infinite Hahn series; finite-support can't hold it. So `versor_inverse`
  succeeds iff the spinor norm `v ṽ` is a scalar *and* a monomial. Intended.
- **`scalar * multivector` works via the scalar's `__mul__` returning
  `NotImplemented`** so Python falls back to the MV's `__rmul__`. Don't make the
  scalar ops raise on a non-scalar operand — that breaks `omega() * e0`.
- **`nim_mul`'s `1u128 << (1u128 << n)` looks overflow-prone.** It isn't for valid
  u128 inputs: bit positions are < 128, so Fermat indices `n ≤ 6` and the shift is
  ≤ 64.
- **`nim_mul_mex` is the slow *game* definition (the mex recurrence), for
  validation and small arguments only.** It's exponential in the argument size —
  fine up to ~48, infeasible over a whole field like F_{2^16}. For real
  computation use the algebraic product (`nim_mul` / `Nimber.__mul__`), which it
  is proven equal to. Experiments use the fast product and only `nim_mul_mex` on
  tiny fields.
- **Pyright flags `import pleroma` as unresolved.** It's installed in `.venv`;
  the editor's interpreter is the system Python. `.venv/bin/python` runs fine.
- **The `neg_one` branch in `Multivector::display` never fires for nimbers.**
  `neg(one)=one` in char 2, so the `coeff==one` branch catches it first.
  Harmless.
- **`Game::canonical_string` canonicalizes; `structural_string` does not.**
  `structural_string` is an order-independent fingerprint of the tree *as given*
  (so `(↑−↑).structural_string() ≠ 0`); `canonical_string` reduces first, so it
  *is* a value key. `structural_eq` compares the as-given structures — use it as
  `a.canonical().structural_eq(&b.canonical())`, or just compare `canonical_string`s.
- **`forms::diagonalize`/`as_diagonal` return `None` in characteristic 2.** Not a
  bug: a nonsingular char-2 form has an alternating polar form and is *not*
  diagonalizable. The char-2 leg classifies via the symplectic Arf reduction
  (`char2.rs`) instead, which takes the full (q, b) metric directly.

## Math facts worth not re-deriving

- nim-field: `F_{2^{2^k}}` = nimbers `< 2^{2^k}`. `F_n ⊗ F_n = (3/2)F_n` for a
  Fermat 2-power `F_n = 2^{2^n}`; distinct Fermat powers multiply ordinarily.
- A real-closed field gives the full Cl(p,q) classification (8-fold periodicity);
  that's why the surreal backend reproduces ℝ-Clifford with exotic scalars.
- Surreal CNF is the Hahn series field ℝ((ω^No)); the ω-map is the monomial map
  and `ω^a·ω^b = ω^{a+b}` is a group homomorphism (No,+) → (No>0,×).
