# AGENTS.md — pleroma

Working notes for agents editing this repo. Global rules still apply.

## What this is

Clifford algebras (with nilpotents) over commutative scalar worlds adjacent to
Conway's combinatorial games. Games under disjunctive sum are an abelian group,
**not a ring** — Conway multiplication is only defined on the number/nimber
subclasses. A Clifford algebra needs a commutative scalar ring, so the direct
game-valued Clifford story only lives on the field-like cores:

- **nimbers** `Nimber(u128)` — the finite field `F_{2^128}` with nim
  arithmetic, characteristic **2**. This is the backend where the char-2
  distinction between quadratic data and polar data (`q ≠ b`) matters most.
  Full `On₂` is the mathematical horizon, not what the fixed-width backend
  stores.
- **surreals** `Surreal` — finite-support Hahn/CNF numbers with rational
  coefficients, char 0. The real-closed Clifford table is exposed only when the
  represented metric entries are exactly square-equivalent to `±1`; entries may
  still be infinite/infinitesimal.
- **surcomplex** `Surcomplex<Surreal>` — adjoin `i` to the implemented surreal
  backend. The algebraically-closed Clifford table is likewise restricted to
  represented exact square classes.

The repo also has comparison scalar worlds (`Fp/Fpn`, `Zp/Qp/Qq`, Laurent,
ramified/Gauss functors, and an adelic precision model) for form-theory
experiments. A pure Rust math core, generic over a `Scalar` trait, with PyO3
per-backend bindings on top. "With nilpotents" = the quadratic form may be degenerate
(`q[i]=0` ⇒ `eᵢ²=0`); all-zero `q` is the exterior/Grassmann algebra.

## Layout

```
src/
  lib.rs          # crate root: the four pillars + the (feature-gated) py module.
                # Each pillar's mod.rs re-exports its children flat, so public
                # paths stay shallow (scalar::Nimber, clifford::sandwich, …).

  scalar/         # PILLAR — the commutative coefficient worlds (generic Scalar).
                  # Grouped BY PLACE into the "any number" table: each field beside
                  # its ring of integers. Subdirs re-export flat, so public paths
                  # stay shallow (scalar::Nimber, scalar::Surreal, …) regardless of
                  # depth. (The char trichotomy that organises forms/ cuts ACROSS
                  # this table — the two pillars are complementary views.)
    mod.rs        # the Scalar trait (add/neg/mul/zero/one/is_zero/inv/
                  # characteristic) + the "any number" table doc + the flat
                  # re-export hub chaining the family subdirs below. Also the
                  # impl_scalar_ops! macro: every backend gets concrete-type
                  # operators (+ - * and unary -) forwarding to the trait methods
                  # (so `Surreal + Surreal`, `-nimber` work). NOT a Scalar
                  # supertrait — see the operators note under "things that look
                  # like bugs". `/` stays a method (inv is partial).
    integrality.rs # the (field, ring-of-integers) pairing made STRUCTURAL: the
                  # HasFractionField {Frac; to_fraction} + HasRingOfIntegers {Int;
                  # is_integral/to_integer} trait pair (Int: HasFractionField<Frac=
                  # Self> ties the loop). Impl'd for the four distinct-type rows
                  # (ℤ⊂ℚ, Oz⊂No, Zp⊂Qp, W_N⊂Qq) PLUS the blanket Surcomplex transport
                  # (ℤ[i]⊂ℚ[i] falls out); the generic frac∘int=id round-trip is the
                  # manifest's test. Laurent/Ramified F_q[[t]]/O[π] are same-type
                  # valuation subrings, so they stay out (is_integral only) — honest.
    valued.rs     # the Valued trait: a discrete valuation + canonical uniformizer ϖ,
                  # impl'd for the local FIELDS (Qp/Qq/Laurent). The spine of the
                  # "local fields" view (cuts across small/ + functor/), and the datum
                  # Ramified folds from its base. NOT a Scalar supertrait (rings of
                  # integers + exact Archimedean worlds excluded), same as the ops.
    analytic.rs   # the ANALYTIC layer UNIFIED: the root-taking surface as two traits
                  # split on where precision lives. ExactRoots {is_square; sqrt} (no
                  # precision arg — exact, or exact to the type's K) for Rational,
                  # Nimber, Fp, Fpn, Zp, Qp, Qq, WittVec, Surreal (exact via the
                  # fixed-point bridge over the lazy roots), Laurent (Newton in odd/0
                  # char, even-exponent inverse-Frobenius in char 2), AND the blanket
                  # Surcomplex<R: ExactRoots+Ordered> — the algebraic-closure √(a+bi)
                  # that used to be a private helper in forms/char0 (the classifier now
                  # calls the trait). SeriesRoots {sqrt_to_terms; nth_root_to_terms;
                  # inv_to_terms} (caller-chosen n) is the lazy interface — Surreal-only
                  # (the one world with unbounded, not type-fixed, precision). Ordered
                  # {sign} is the datum the Surcomplex blanket needs to pick the branch.
                  # The residue Tonelli roots (fp_sqrt/fq_sqrt) live here now (shared
                  # with small/analytic's Hensel seed). Surcomplex<Surreal> also gets a
                  # lazy inv_to_terms (division when the norm a²+b² is non-monomial).
                  # Gauss/Ramified excluded honestly (rational-function/ramified roots
                  # need a different machine). NOT a Scalar supertrait, like valued.
    functor/      # the root-level *functors* — ways to GROW a field, orthogonal to
                  # the place table. 2×2 (algebraic|transcendental × residue|value-
                  # extending), ALL FOUR corners filled (see functor/mod.rs).
      surcomplex.rs # Surcomplex<S> = adjoin i over ANY backend (carries conj()). The
                  #   ALGEBRAIC, residue-extending corner (root of x²+1), not a world.
      laurent.rs  #   Laurent<S, const K> = S((t)), formal Laurent series to relative
                  #   precision K. The TRANSCENDENTAL, VALUE-extending corner (adjoin t
                  #   as a uniformizer, v(t)=1). Over a finite field it fills the EQUAL-
                  #   CHARACTERISTIC local cell: F_q((t)) (char p, the mirror of Qp);
                  #   ring of integers F_q[[t]] = the val≥0 subring (Laurent::is_
                  #   integral). Capped-relative like Qp (mul/inv exact, addition non-
                  #   assoc across precision); field iff S is. EXCLUDED from the fuzz.
      ramified.rs #   Ramified<S, const E> (was Eisenstein) = adjoin a root of the
                  #   Eisenstein polynomial xᴱ−ϖ over a Valued base. The ALGEBRAIC,
                  #   value-extending corner — the RAMIFIED local cell: Q_p(p^{1/E})
                  #   over Qp, the ramified twin of Qq. Exact ramified valuation
                  #   (unique min via distinct residues mod E); E=2 norm-form inverse,
                  #   E≥3 regular-rep solve. Always a field (Eisenstein's criterion),
                  #   incl. wild/inseparable p|E. Capped-relative, EXCLUDED from fuzz.
      gauss.rs    #   Gauss<S> = S(t) with the GAUSS valuation over a Valued base. The
                  #   TRANSCENDENTAL, RESIDUE-extending corner (adjoin t as a UNIT,
                  #   v(t)=0, transcendental residue ⇒ residue field k(t̄), value group
                  #   unchanged) — the last corner, Laurent's residue-extending twin.
                  #   num/den rational functions, NO gcd (inv=den/num; eq by cross-
                  #   mult; monic denom). Valued itself; precision model, EXCLUDED.

    global/       # FAMILY — the adelic/global place. Adele is a finite-precision
                  # restricted-product model over Q, with LocalQp as the runtime-prime
                  # p-adic cell. Useful for product formula / Hilbert reciprocity /
                  # Hasse-Minkowski experiments in forms/adelic.rs; not an exact
                  # infinite-memory adele implementation.

    exact/        # FAMILY — the Archimedean char-0 base (field + ring of integers)
      rational.rs # exact ℚ over i128, NOT a game backend — the char-0 scalar that
                  # validates the geometric product against the known Cl(p,q)
                  # classification before we trust the exotic backends.
      integer.rs  # exact ℤ — the coefficient ring for the exterior algebra of the
                  # game group (games/game_exterior.rs): games are a ℤ-module, not a
                  # ring, so Λ over ℤ is the structure that lives on all of
                  # game-world. Only ±1 invertible (Grassmann never calls inv).

    big/          # FAMILY — the transfinite worlds (the number may be infinite)
      cnf.rs      # the ONE thing surreal & ordinal genuinely share: merge_descending,
                  # the descending-CNF canonicalizer parameterized by the 3 places
                  # they differ (exponent order: No value-order vs ordinal lex;
                  # coeff merge: + vs XOR; zero test). Deliberately a shared
                  # FUNCTION, not a Cnf<C> TYPE — the orders/algebras diverge (No is
                  # a field, On₂ isn't), so a shared type would be a false identity.
      surreal/    # finite-support surreal Hahn/CNF backend (char 0). SPLIT into a subdir, all impl Surreal:
        mod.rs    #   CNF core: Vec<(exponent: Surreal, coeff: Rational)>, recursive
                  #   exponents, Hahn arithmetic ω^a·ω^b = ω^{a+b}, Scalar, Debug,
                  #   the truncate() precision knob, and the (shared) test module.
        simplicity.rs    # the {L|R}/simplicity bridge (dyadic case): as_rational/
                  #   as_dyadic/is_dyadic/dyadic_birthday + simplest_above/_below/
                  #   simplest_between, and floor/frac (the Oz bridge — Omnific::floor
                  #   wraps it). simplest_in_cut is pub(super) for sign_expansion.
        sign_expansion.rs # EXACT sign_expansion/from_sign_expansion (dyadic,
                  #   round-trips, length = birthday) + as_ordinal + from_ordinal (its
                  #   inverse: Ordinal→Surreal CNF) + the TRANSFINITE (Gonshor)
                  #   SignExpansion{runs:(bool,Ordinal)} + birthday_ordinal + the
                  #   transfinite inverse from_transfinite_sign_expansion (closes the
                  #   round trip on the representable subclass). Every ordinal ↦
                  #   all-pluses incl ω^ω; ε=+(−)^ω. None outside the representable
                  #   subclass (√ω, ½ω, ω−1) — honest, not ℝ-trunc.
        analytic.rs      # the LAZY FIELD layer (the SeriesRoots primitives): inv_to_terms
                  #   (Neumann-series inverse to n terms — non-monomials too, the Zp-style
                  #   precision contract) + sqrt_to_terms/nth_root_to_terms (real-closed
                  #   roots to n terms; Some iff the leading coeff is a perfect ℚ-power, so
                  #   √2/√(2ω) honestly None, √ω exact). Named *_to_terms (matching
                  #   inv_to_terms) so the precision-free names sqrt/is_square belong to
                  #   ExactRoots (scalar/analytic.rs), whose Surreal sqrt is the exact
                  #   value recovered by squaring these truncations back.
      omnific.rs  # the omnific integers Oz: Omnific(Surreal) newtype, a transfinite
                  # commutative RING (not field). Surreal mirror of Integer (the ring
                  # of integers of No); the exterior algebra with ω-scale coeffs.
      ordinal/    # transfinite (ordinal) NIMBERS On₂ — the char-2 mirror of surreal
                  # (was onag/, after Conway's ONAG; renamed name-by-object for the
                  # type Ordinal). SPLIT into a subdir like surreal/, all impl Ordinal:
        mod.rs    #   CNF core: Ordinal = Vec<(exponent: Ordinal, coeff: u128)>,
                  #   constructors, the lexicographic cmp, as_finite, Debug.
        nim.rs    #   char-2 NIM arithmetic: nim_add (coeff XOR) COMPLETE; nim_mul
                  #   implemented below ω^ω via the current degree-3 tower. The old
                  #   φ_{ω+1}/<ω³ path is still tested as the one-generator case
                  #   (ω⊗ω⊗ω=2, F₄(ω)≅F₆₄). At ω^ω and above, multiplication returns
                  #   None rather than speculating. The XOR canonicalize (= the
                  #   char-2 coeff merge) lives here.
        cantor.rs #   ORDINARY (Cantor) ord_add/ord_mul (NOT nim: ω+ω=ω·2, 1+ω=ω) —
                  #   the surreal birthday's run-length arithmetic. A distinct
                  #   algebra from nim, sharing only the CNF shape. (Was ordinal.rs.)

    small/        # FAMILY — the non-Archimedean (p-adic) local world
      qp.rs       # Qp<const P, const K>: the p-adic FIELD Q_p (field of fractions of
                  # Zp; the p-adic mirror of ℚ / of Omnific⊂Surreal). p^val·unit,
                  # char()=0, inv TOTAL on nonzero (1/p exists). CAPPED-RELATIVE
                  # precision: mul/inv exact, addition NOT associative across
                  # precision boundaries (a precision model, like float) — used at
                  # the forms layer (valuation/residue), EXCLUDED from the exact-ring
                  # fuzz suite.
      zp.rs       # Zp<const P, const K>: the p-adic integers Z_p to precision k
                  # (= Z/p^k), the ring of integers of Q_p. A LOCAL RING (p a
                  # non-unit) — char()=0, inv = Omnific pattern (units only). Cl over
                  # it is non-semisimple. residue field F_p.
      qq.rs       # Qq<const P, const N, const F>: the UNRAMIFIED extension Q_q =
                  # Frac(W_N(F_q)) of Q_p, residue degree F (residue field F_q). The
                  # field of fractions of WittVec — to WittVec what Qp is to Zp; Qq
                  # with F=1 IS Qp. p^val·(Witt unit), char()=0, inv TOTAL on nonzero.
                  # Capped-relative (excluded from the fuzz suite). Completes the
                  # (field, ring of integers) pairing on the unramified leg.
      analytic.rs # the p-adic ANALYTIC layer over all four backends (mirror of surreal/
                  # analytic.rs): Hensel-lifted is_square/sqrt (Newton lift y←(y+u/y)/2
                  # seeded by the Tonelli residue roots — fp_sqrt/fq_sqrt now live in
                  # scalar/analytic.rs and are imported here; ODD p only — dyadic sqrt
                  # is the forms mod-8 story, asserted) and the Teichmüller rep τ
                  # (the (q−1)th root of unity; added to Zp/Qp/Qq, WittVec already
                  # had it). The inherent is_square/sqrt are what ExactRoots delegates to
                  # (scalar/analytic.rs). nth_root + p-adic log/exp are the next ops.

    finite_field/ # FAMILY — the finite residue worlds (the char trichotomy's finite
                  # leg + the unramified ring of integers)
      mod.rs      # the FiniteField TRAIT: the shared Galois engine (degree,
                  # conjugates, min_poly_monic, relative_trace/_norm[_over],
                  # multiplicative_order, is_primitive, discrete_log) as default
                  # methods. An impl supplies only frobenius, integer pow, ext_degree,
                  # group_order, and group_order_factors. nimber + fpn both impl it —
                  # one verified algorithm, two backends (was duplicated).
      fp.rs       # Fp<const P>: the prime field F_P (odd characteristic), the residue
                  # field of Zp and the base of every extension here. Genuine neg.
      fpn.rs      # Fpn<const P, const N>: finite extension fields F_{p^N} via a
                  # (P,N)-keyed irreducible reduction poly. Completes the odd-char
                  # tower AND the char-2 odd-DEGREE fields nimbers can't reach (F_8).
                  # Schoolbook mul+reduce, Fermat inv, is_square (Euler/Frobenius).
                  # impl FiniteField (frobenius/pow/ext_degree/group_order); keeps
                  # only min_poly (F_p projection) + primitive_element. (NB the static
                  # order() = field order p^N, ≠ multiplicative_order(&self).)
      nimber/     # On₂ in u128 (= F_{2^128}), split by layer while re-exporting the
                  # same nim_* surface flat:
        mod.rs    #   Nimber wrapper + Scalar impl + public re-exports.
        arithmetic.rs # nim_add = XOR; nim_mul via Fermat-power recursion memoised
                  #   on 2^i ⊗ 2^j; nim_square/nim_sqrt/nim_inv.
        artin_schreier.rs # nim_trace + y²+y=c solver (solvable ⇔ Tr(c)=0).
        galois.rs #   impl FiniteField; degree/conjugates/min_poly/relative
                  #   trace/norm/order/primitive/discrete-log. OVERRIDES is_primitive
                  #   and discrete_log with Pohlig–Hellman + BSGS over ORDER_FACTORS.
      wittvec.rs  # WittVec<const P, const N, const F>: Witt vectors W_N(F_q), as the
                  # truncated unramified ring (Z/p^N)[t]/(f̃) (NOT the forms Witt
                  # group). The char-p analogue of Z_p (= W(F_p)) — the ring of
                  # integers of the unramified extension (its FIELD of fractions is
                  # small/qq.rs). Witt/Teichmüller coords + carry-formula oracle, plus
                  # p_valuation/try_divide_by_p (the hooks Qq uses to peel a unit).

  linalg/         # crate-private shared linear algebra, deliberately below the
                  # mathematical pillars rather than a public API.
    field.rs      # Gaussian solve / inverse_matrix / unit-pivot nullspace over a
                  # Scalar field. Used by multivector_inverse, blade analysis, and
                  # inverse_outermorphism.
    f2.rs         # nim-field row rank for F₂/F_{2^k}-style Dickson computations.
    integer.rs    # integer relation row normalization + vector reduction for the
                  # game exterior algebra's lattice quotient.

  clifford/       # PILLAR — the multivector engine + GA layer (generic Scalar)
    mod.rs        # thin hub: re-exports engine + versor + the structured-algebra
                  # modules flat (clifford::Metric, clifford::sandwich, …).
    engine.rs     # thin engine hub + product/regression tests. Public paths stay
                  # clifford::Metric / CliffordAlgebra / Multivector / bits / grade.
    engine/       # the associative-algebra core split by concept:
      basis.rs    # bits / grade / MAX_BASIS_DIM / wedge_sign.
      metric.rs   # Metric {q,b,a}, constructors, direct_sum, q_val/has_upper.
      product.rs  # geom_product_blades (general-bilinear Chevalley product) plus
                  # the cfg(test) reduce_word oracle it is cross-validated against.
      algebra.rs  # CliffordAlgebra<S>: blade arithmetic, grade projection,
                  # wedge/reverse/graded_tensor/embeddings.
      multivector.rs # Multivector<S>: term store, zero/display helpers.
      inverse.rs  # GENERAL multivector_inverse via the shared linalg::field solver.
      terms.rs    # local term-map scale/merge helpers.
    versor.rs     # the GA layer on top: versor_inverse, sandwich, twisted_sandwich
                  # (Pin action), reflect, left/right_contract, dual/undual,
                  # grade_involution, norm2, even_part / even_subalgebra. Plus the
                  # product/involution suite: clifford_conjugate, scalar_product
                  # ⟨ab⟩₀, commutator/anticommutator (½-free, char-faithful), and
                  # the regressive meet a∨b (intersection dual to wedge). Plus the
                  # CAYLEY transform cayley/cayley_inverse = (1−B)(1+B)⁻¹: the exact
                  # RATIONAL bivector↔rotor map (Lie algebra ↔ Spin group, no cos/sin,
                  # char≠2) — uses engine::multivector_inverse since 1+B isn't a versor.
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
    divided_power.rs # the CHAR-FAITHFUL symmetric mirror of hopf.rs: the divided
                  # power algebra Γ(V) (dual of Sym), DividedPowerAlgebra +
                  # DpVector{terms: multidegree→coeff}. BINOMIAL product (the dual of
                  # Sym's free product), DECONCATENATION coproduct (sign-free mirror
                  # of the exterior unshuffle), antipode = (−1)^{|α|} grade
                  # involution. Binomials embed via the scalar's own +, so they
                  # reduce mod char: (γᵢ⁽¹⁾)²=2γᵢ⁽²⁾=0 in char 2 while γᵢ⁽²⁾≠0 — the
                  # honest Γ≠Sym (mirror of exterior eᵢ²=0). Axioms tested over
                  # Rational AND Nimber. Standalone (own monomials, not the blade
                  # engine); no Python binding.
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
    mod.rs        # re-exports the legs + classify + diagonalize/equivalence +
                  # witt/witt_ring + brauer_wall + padic + adelic + springer.
    classify.rs   # the classifier FAÇADE: ClassifyForm + WittClassify +
                  # IsometryClassify + WittDecompose + BrauerWallClassify, keyed on
                  # the scalar so `metric.classify()` / `.witt_class()` /
                  # `.isometric_to()` / `.witt_decompose()` / `.bw_class()` pick the
                  # right leg at compile time (Surreal→CliffordType, Fp/Fpn→
                  # OddCharType, Nimber→ArfResult, …). Rational & Surcomplex impl
                  # ClassifyForm but not WittClassify (their Witt data isn't a single
                  # WittClassG — honest, not a gap).
    diagonalize.rs # congruence diagonalization (char ≠ 2): gram, diagonalize,
                  # as_diagonal. Returns None in char 2 (nonsingular forms aren't
                  # diagonalizable there — use char2/arf.rs's symplectic reduction).
                  # This is what lets char0/oddchar classify ARBITRARY metrics.
    equivalence.rs # isometry (per backend, via the complete invariant) +
                  # Witt decomposition (k·H ⊥ anisotropic kernel) over ℝ and F_q.
    char0.rs      # (was classify.rs) the char-0 Clifford classifier: Cl(p,q) →
                  # matrix algebra over ℝ/ℂ/ℍ via the 8-fold table (real-closed
                  # surreal/rational) and the 2-fold table (surcomplex). Non-diagonal
                  # metrics are diagonalized first (signature is pub(crate)).
    oddchar/      # odd-characteristic forms, re-exported flat:
      mod.rs      #   hub and trichotomy docs.
      field.rs    #   FiniteOddField unifies Fp and Fpn square classes + metadata.
      invariants.rs # classify_finite_odd / finite_odd_witt / discriminant_finite_odd
                  #   / hasse_invariant_finite_odd (≡ +1 over finite fields): ONE
                  #   generic FiniteOddField implementation keyed off the trait, so
                  #   Fp and Fpn share the path (the const-P prime-field wrappers were
                  #   folded away). dim+disc complete. Non-diagonal metrics
                  #   diagonalized first.
    char2/        # characteristic-2 invariants, re-exported flat:
      mod.rs      #   hub: forms::arf_invariant / forms::dickson_matrix stay shallow.
      arf.rs      #   the Arf invariant (char-2 classifier): arf_f2 (F₂ bitmask) +
                  #   arf_nimber (any nim-field, symplectic reduction + trace).
      dickson.rs  #   Dickson invariant: dickson_matrix = rank(g−I) mod 2, ker = SO;
                  #   dickson_of_versor delegates to generic versor grade parity.
    quadric_fit.rs # the "is this P-set a quadric?" research BENCH (split from
                  # the char2 classifier): fit_f2_quadratic (Gaussian elim over the 2^k membership
                  # equations) + QuadricFit{constant,qd,bmat,arf} + is_genuinely_
                  # quadratic. The instrument the game probes / misere_quotient /
                  # octal_hunt feed P-positions into — distinct from the classifier.
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
    adelic.rs     # local-global rational form helpers: Hilbert product over all
                  # places, rank>=3 adelic Hasse-Minkowski breakdown, and Brauer
                  # local invariant sums. Reuses padic.rs; not a new symbol engine.
    springer.rs   # non-Archimedean Springer decomposition over the surreals: a
                  # diagonal form's ω-adic valuation filtration into residue ℝ-forms.
                  # Honest: value group 2-divisible ⇒ W(No)=W(ℝ)=ℤ; the filtration
                  # itself is the novelty.
    springer_padic.rs # the p-adic MIRROR of springer.rs over the Qp field: bucket a
                  # diagonal form by p-adic valuation, residue F_p-forms (dim + disc
                  # square-class). The novelty vs surreal: value group ℤ NOT 2-divisible
                  # ⇒ TWO residue layers (val even/odd) survive (parity_layer) =
                  # W(Q_p)=W(F_p)⊕W(F_p). Odd p, already-diagonal only.
    springer_laurent.rs # the THIRD Springer sibling, over the equal-characteristic
                  # local field F_q((t)) (Laurent backend). t-adic valuation buckets,
                  # residue F_q-forms (dim + disc square-class via is_square_finite).
                  # Novelty vs the Qp twin: EQUAL characteristic (the field is char p)
                  # and a general residue field F_q (extension square-class, not just
                  # F_p). Value group ℤ ⇒ two parity layers = W(F_q((t)))=W(F_q)². Odd
                  # residue char, already-diagonal only; residue char 2 REJECTED (the
                  # Springer second-residue-map boundary — char-2 Witt lives in char2/).
    invariants.rs # numeric FIELD INVARIANTS the Witt ring implies: level/Stufe s(F),
                  # pythagoras_number, u_invariant, is_sum_of_n_squares — computed over
                  # finite F_p (level≤2, u=2); ℝ/Q_p textbook constants documented.
    hermitian.rs  # HERMITIAN forms over Surcomplex (the involution conj() the forms
                  # pillar never used): HermitianForm (conj-symmetric Gram), unitary
                  # (conjugate) congruence diagonalize → real diagonal, signature
                  # (Sylvester, the complete invariant = U(p,q)). Generic sign closure.

  games/          # PILLAR — combinatorial game theory (mostly Scalar-free)
    mod.rs        # re-exports the modules below flat.
    piecewise.rs  # Pl: exact rational piecewise-linear wall arithmetic (max/min,
                  # subtraction, crossings, cleanup) used by thermography.
    thermography.rs # temperature theory: the thermograph of a short game — left/right
                  # scaffolds, stops, cooling (cooled_stops), temperature, and mean
                  # (mast) value. Switches/numbers/↑/⋆ pinned; mean is additive.
                  # (Atomic weight now lives in atomic_weight.rs.)
    atomic_weight.rs # atomic weight of ALL-SMALL games (finishes thermography): the
                  # two-ahead rule (Siegel Constructive Atomic Weight; Larsson–
                  # Nowakowski arXiv:2007.03949 Thm 10). A={aw(G^L)−2|aw(G^R)+2};
                  # non-integer ⇒ aw=A, else far-star ⋆N (N=birthday+1) max/min over
                  # A's RAW option games (handles fractional option weights — a naive
                  # 1+max_R aw is WRONG, Codex-caught). aw IS additive on all-small.
                  # Game::is_all_small / Game::nim_heap live in partizan.rs.
    hackenbush.rs # red/blue/green Hackenbush: Hackenbush{edges, ground=0}, to_game()
                  # (the universal evaluator via move-and-prune), value() → surreal
                  # number (blue–red), grundy() → nimber (all-green = Nim). The
                  # bridge tying surreals+nimbers+sign-expansion through one object.
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
    misere.rs     # misère-play outcomes (misere_is_n/_is_p) for finite acyclic
                  # impartial games; misère Nim vs Bouton; the bounded
                  # indistinguishability quotient (misere_quotient); octal games
                  # (octal_moves, octal_misere_quotient). The non-linear route.
    partizan.rs   # short partizan games (sum/neg/order/birthday/is_number) + the
                  # CANONICAL FORM (dominated/reversible reduction; structural_string
                  # vs canonical_string — the latter canonicalizes, a value key) +
                  # the game↔surreal bridge (number_value / from_surreal, numbers
                  # only). Also Game::ordinal_sum (G:H — Hackenbush strings are
                  # these), Game::nim_heap (⋆n, the far star) + Game::is_all_small
                  # (atomic-weight domain). The Λ-of-the-game-group exterior algebra
                  # split out to game_exterior.rs. NB: distinct from coin_turning.rs.
    number_game.rs # transfinite NUMBER games (ω, ε) carried by their Surreal value —
                  # value/birthday(Ordinal)/sum/cmp delegate to surreal, no infinite
                  # option tree (the finite Game engine is untouched). Plus the FULL
                  # transfinite round trip via sign_expansion/from_sign_expansion: the
                  # run-length sign expansion is the finite encoding of the (infinite)
                  # {L|R} tree (ω={0,1,2,…|} can't be listed but +^ω is finite data),
                  # closing surreal↔game THROUGH the canonical birthday path, not a
                  # stored value — the transfinite analogue of the dyadic Game bridge.
    game_exterior.rs # the exterior algebra of the GAME group: Λ over ℤ on game
                  # generators (living on all of game-world, incl. non-numbers ⋆/↑ —
                  # needs only the ℤ-module structure, not the game product). Split
                  # from partizan.rs: GameExterior (free Grassmann engine quotiented
                  # by integer game relations such as 2⋆=0, propagated through the
                  # exterior ideal) + GameRelation; integer-lattice normalization and
                  # reduction live in linalg/integer.rs.

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
    forms.rs      # classify / witt / dickson / springer bindings, plus
                  # FiniteFieldForm (runtime Fp/Fpn form wrapper) and compatibility
                  # wrappers for classify_oddchar / witt_decompose_oddchar /
                  # is_isometric_oddchar.
    games.rs      # Game (incl. canonical/is_canonical/canonical_string/
                  # number_value/from_surreal) / GameExterior + nim_mul_mex +
                  # grundy_graph.

examples/tour.rs   # cargo run --example tour   (Rust-only demo)
examples/misere_quotient.rs    # misère quotients + the quadric test on P-sets
examples/interactive_kernel.rs # B-coupled interactive games vs {Q=0}
examples/octal_hunt.rs         # sweep octal games for a (ℤ/2)^k quadric P-set
                               # (cargo run --release --example octal_hunt)
examples/loopy_quadric.rs      # cyclic (Draw-set) rules vs {Q=0}; the radical collapse
examples/bent_route.rs         # route probes on one BENT form: B+frame reaches a right-
                               # Arf quadric class there; the local-field (Ising)
                               # completion fails to hit the target Gold zero set
demo.py            # the same tour from Python
experiments/       # research probes ON TOP of the shipped lib: Arf of Gold
                   # forms, the game-built synthesis, the Arf win-bias,
                   # artin_arf (the trace ↔ Arf unification),
                   # open_question_probe (the polar-form obstruction),
                   # tartan_bilinear (B realized by Turning-Corners),
                   # framing_obstruction (the Sp(B) no-go + the diagonal-framing
                   # ladder for the open question), and gold_family_survey (the
                   # trace-family probe Σ Tr(c_i x^{1+2^i}) and where sampled/scaled
                   # cases go BENT — components Tr(λ x^{1+2^a}) bent for 2/3 of λ
                   # in the APN cases tested), and misere_kernel (the misère kernel
                   # obstruction: the canonical group/kernel shadow carries the
                   # XOR-linear normal-play P-set by Plambeck–Siegel Thm 6.4;
                   # checked on R8). See NOTES.md.
```

The math thread (Arf↔Clifford, the games bridge, the char-0/char-2 classifier
symmetry, the Artin–Schreier ↔ Arf unification, the open play-semantics
question) is written up in `NOTES.md` — read it before touching `forms/char2/`,
`forms/quadric_fit.rs`, `forms/char0.rs`, `games/coin_turning.rs`,
`games/kernel.rs`, `games/misere.rs`, `forms/witt.rs`, `experiments/`, or the
`misere_quotient` / `interactive_kernel` examples.

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
- Rust scalar operators: every backend has `+ - *` and unary `-` (via
  `impl_scalar_ops!` in `scalar/mod.rs`), forwarding to the `Scalar` methods. Use
  them on concrete-typed scalar code; generic engine code over `S: Scalar` still
  calls `.add(&x)`/`.mul(&x)` (operators aren't a supertrait — see the note below).

## Testing

`cargo test` is the source of truth and needs no Python. The Python layer is
smoke-tested via `demo.py`. After touching `clifford/` or `scalar/big/surreal/`,
run `cargo test` **and** rebuild + run `demo.py` — display changes don't surface
in `cargo test`.

**`cargo test` does NOT compile the `python` feature** (it's gated, and that's
deliberate — keeps the core libpython-free). So a green `cargo test` can hide a
broken `py/` build: after touching `py/` *or any core API the bindings call*
(e.g. renaming a `Scalar`/`FiniteField` method), run `cargo check --features
python`. The Galois-trait unification once shipped a broken `py/scalars.rs` this
way (the nim Galois methods moved to the trait; the bindings now call the
`nim_*` u128 free fns instead).

Beyond the per-module unit tests there are two **property-based** suites (dev-dep
`proptest`, integration tests in `tests/`): `tests/scalar_axioms.rs` fuzzes the
commutative-ring axioms across every backend, and `tests/clifford_axioms.rs`
fuzzes geometric-product associativity/distributivity over random metrics in
char 0 and char 2. These are the randomized safety net under the
`Scalar`-is-a-commutative-ring assumption and the product engine. (serde
(de)serialization is intentionally NOT shipped yet — the invariant-carrying types
need custom invariant-preserving deserialization, not a naive derive.)

## Things that look like bugs but are not

- **`Surreal` has two square roots, by design.** The inherent `sqrt_to_terms(n)` is
  the lazy `SeriesRoots` primitive (n leading terms); `ExactRoots::sqrt(&self)` (0
  args) is the exact value. They are *different methods with different arities*, so
  on a concrete `Surreal` `x.sqrt_to_terms(4)` is lazy and `x.sqrt()` (with
  `ExactRoots` in scope) is exact. The rename to `*_to_terms` (matching the existing
  `inv_to_terms`) is what frees the bare names `sqrt`/`is_square` for the
  precision-free `ExactRoots` family across every backend. Don't "unify" them back to
  one `sqrt` — the two precision contracts are genuinely different. The Python
  `Surreal.sqrt(n)` stays the lazy one; `Surreal.exact_sqrt()` is the exact one.
- **`ExactRoots`/`SeriesRoots`/`Ordered` are NOT `Scalar` supertraits** (same as the
  operators and `Valued`): not every world takes roots, so the bounds stay opt-in.
  The trait impls *delegate to inherent methods of the same name* (e.g.
  `ExactRoots::sqrt` for `Qp` calls the inherent `Qp::sqrt`); inherent-shadows-trait
  in method position makes that delegate-not-recurse, exactly the `Valued` pattern.
- **`ExactRoots::sqrt`/`is_square` on `Zp`/`Qp`/`Qq`/`WittVec` panic at p=2.** They
  inherit the inherent odd-p assertion (the dyadic case is the forms mod-8 story).
  The finite fields and `Laurent` handle char 2 natively, so the trait is total
  there; the p-adic rings are the documented exception, not a bug.
- **Scalar `+ - *` operators are concrete-only, NOT a `Scalar` supertrait.** This
  is deliberate: making `Scalar: Add+Sub+Mul+Neg` brings the ops into scope for
  every generic `S`, where `Mul::mul(self, Self)` shadows `Scalar::mul(&self,
  &Self)` at owned-receiver sites (`m[i][j].mul(&x)`) and forces clones the
  borrow-based engine avoids (70+ generic sites broke when tried). Concrete-only
  ops give users `Surreal + Surreal` with zero churn to the generic core. Don't
  "promote" them to a supertrait, and don't migrate the engine's `.add()`/`.mul()`
  to operators — the `&self` methods are the right tool there.
- **Char-2 Clifford over an orthogonal basis is commutative.** `e0*e1 == e1*e0`
  when `b` is empty and the scalar is a nimber. Correct: `{e0,e1}=2B=0` and
  `-1=1`. Set an off-diagonal `b[(i,j)]` to get non-commutativity.
- **Surcomplex over nimbers is degenerate.** `i²=1`, `(1+i)²=0`, not a field.
  Full `On₂` is algebraically closed, while the shipped `Nimber(u128)` is
  `F_{2^128}`; either way the char-2 adjunction is not the useful complex
  backend. Surcomplex is only meaningful over char-0 scalar worlds here.
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
  (`forms::char2`) instead, which takes the full (q, b) metric directly.
- **`Qp` addition is not associative across precision boundaries.** Capped-relative
  precision (the standard p-adic model, like float): mul/inv are exact, but
  additive cancellation below the retained window reads as `0`. `Qp` is a *precision
  model*, not an exact ring — that's why it's excluded from `tests/scalar_axioms.rs`
  and used only at the forms layer (valuation/residue). Don't "fix" it to be exact;
  no finite-memory exact `Q_p` exists (`1/(p+1)` has infinite support).
- **`Surreal::sqrt`/`nth_root` return `None` for `√2`, `√(2ω)`, `½ω`'s roots, etc.**
  Intended: the leading coefficient must be a perfect ℚ-power. `√2` is not a
  finite-CNF-with-ℚ-coeffs surreal (it's born at ω), so we honestly decline — same
  ethos as "coefficients are ℚ, not ℝ". `√ω = ω^{1/2}` IS exact (monomial).
- **`Surreal::birthday_ordinal`/`transfinite_sign_expansion` are `None` outside the
  representable subclass** (`√ω`, `ω−1`, `½ω`, mixed ordinal+infinitesimal). Every
  *ordinal* (incl. `ω^ω`) is handled (all-pluses); `ε` is the one infinitesimal
  pinned. Not a gap to close blindly — it's the honest Gonshor scope boundary.
- **`Fpn::order()` is the field order `p^N` (static, no self); the element's
  multiplicative order is `multiplicative_order(&self)`.** Different things; the name
  split is deliberate (the static `order()` predates the Galois toolkit).
- **The `nim_*` Galois free fns delegate to the `FiniteField` trait; don't re-add
  inherent `Nimber` Galois methods.** `nim_degree(x)` etc. call `Nimber(x).degree()`
  (the trait default). An inherent `Nimber::degree` would take precedence and
  recurse forever back through the free fn — that's why the old `impl Nimber`
  wrappers were removed. The trait lives in `scalar/finite_field/mod.rs`; to add a
  Galois op, add a default method there (both nimber and fpn get it free). Nimber
  *overrides* `is_primitive` / `discrete_log` for the sharper large-field algorithms
  — that's intended, not duplication.
- **Atomic weight's integer branch is NOT `1 + max_R aw(G^R)`.** It's a predicate
  over `A`'s raw option *games* (`A^R = aw(G^R)+2`) comparing an integer `n` via
  `le`/`fuzzy`, bounded by the *tightest* right option — so it stays correct when an
  option's atomic weight is a fraction (e.g. `½`). The naive max-of-integers form
  panics/misreads there (Codex-caught; see `integer_branch_handles_fractional_option_weights`).
  And atomic weight IS additive on all-small games — don't reinstate a "not additive" claim.

## Math facts worth not re-deriving

- nim-field: `F_{2^{2^k}}` = nimbers `< 2^{2^k}`. `F_n ⊗ F_n = (3/2)F_n` for a
  Fermat 2-power `F_n = 2^{2^n}`; distinct Fermat powers multiply ordinarily.
- A real-closed field gives the full Cl(p,q) classification (8-fold periodicity);
  the implemented surreal backend reproduces that table only on the exact-square
  subdomain it can represent.
- Surreal CNF is modeled as finite-support Hahn series with rational
  coefficients; the ω-map is the monomial map and `ω^a·ω^b = ω^{a+b}` is a group
  homomorphism on represented monomials.
