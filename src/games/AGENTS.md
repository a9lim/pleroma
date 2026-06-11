# AGENTS.md — `src/games/`

The PILLAR of combinatorial game theory — the second column of the project,
mostly independent of the scalar/Clifford stack (the bridge is the number/nimber
subclasses, where Conway multiplication is defined). Games under disjunctive sum
are an abelian GROUP, not a ring; that constraint is *why* the Clifford story
lives on the scalar backends and not on all games.

> Read root `OPEN.md` before touching `coin_turning.rs`, `kernel.rs`, `misere.rs`, or
> the example probes — they feed the open play-semantics question.

`mod.rs` re-exports every module below flat.

Fixed-width game payloads use `u128`/`i128`: finite nim heaps, octal codes, Grundy
values, scoring integers, and quotient counters. `usize` is for graph nodes, option
indices, and collection lengths.

## Values & arithmetic

- **`partizan.rs`** — short partizan games (sum/neg/order/birthday/is_number) + the
  CANONICAL FORM (dominated/reversible reduction; `structural_string` vs
  `canonical_string` — the latter canonicalizes, a value key) + the game↔surreal
  bridge (`number_value`/`from_surreal`, numbers only). Also `Game::ordinal_sum`
  (G:H — Hackenbush strings are these), `Game::nim_heap` (⋆n), `Game::is_all_small`.
- **`number_game.rs`** — transfinite NUMBER games (ω, ε) carried by their Surreal
  value — value/birthday/sum/cmp delegate to surreal, no infinite option tree. Plus
  the FULL transfinite round trip via sign_expansion/from_sign_expansion (the
  run-length sign expansion is the finite encoding of the infinite {L|R} tree).
- **`nimber_game.rs`** — the CHAR-2 MIRROR of `number_game.rs`: transfinite NIMBER
  games (Nim heaps `⋆α`) carried by their `Ordinal` (On₂) Grundy value — grundy/
  add(=nim-add XOR)/cmp/to_finite_game delegate to `Ordinal`, no infinite option
  tree. `neg` is the identity (char 2: every impartial game is self-inverse);
  `turning_corners` is the nim-product (Conway's coin game, `ω³=2`); `None` only when
  a Kummer carry needs a prime past the verified table or at `≥ ω^(ω^ω)`. This is the
  `No ↔ On₂` symmetry at the games layer (the rest lives at the scalar layer via the
  shared CNF core, reaching Clifford through `Scalar for Ordinal` inside the checked
  Kummer boundary). Bound to Python as `NimberGame`.
- **`game_exterior.rs`** — the exterior algebra of the GAME group: Λ over ℤ on game
  generators (living on all of game-world, incl. non-numbers ⋆/↑ — needs only the
  ℤ-module structure). `GameExterior` — three constructors: `new` (auto-search for
  integer relations), `free` (no quotient), `with_relations` (explicit), and
  `with_relation_search(bound)` — quotients the free Grassmann engine by integer game
  relations such as 2⋆=0. Carries `GameRelation` + the `GameRelationCertificate` /
  `RelationSearchCertificate` evidence records; lattice normalization in
  `linalg/integer.rs`.

## Temperature theory

- **`thermography.rs`** — the thermograph of a short game: left/right scaffolds,
  stops, cooling (`cooled_stops`), temperature, and mean (mast) value.
- **`atomic_weight.rs`** — atomic weight of ALL-SMALL games (finishes thermography):
  the two-ahead rule (Siegel Constructive Atomic Weight; Larsson–Nowakowski
  arXiv:2007.03949 Thm 10). `aw` IS additive on all-small games.
- **`piecewise.rs`** — `Pl`: exact rational piecewise-linear wall arithmetic used by
  thermography. `add_pl` (pointwise sum) is the tropical `⊗`; `sub_pl` is the arithmetic
  difference (`left_raw − right_raw`) in the meeting-temperature recursion, NOT a
  tropical operation.
- **`tropical_thermography.rs`** — names the latent tropical structure in
  thermography and machine-checks it. The option folds are tropical `⊕` in DUAL
  semirings — the left wall a `(max,+)` fold over the Left options' right walls, the
  right wall a `(min,+)` fold over the Right options' left walls — and cooling is
  tropical `⊗`. `Pl::oplus_max`/`oplus_min`/`otimes` name the wall operations;
  `thermograph_via_tropical` is a parallel recursion pinned EQUAL to
  `thermography::thermograph`. It reuses the identical `pub(crate)`
  freeze/meeting-temperature cooling tail — it only renames the folds, it does not
  reimplement cooling. The `Semiring`/`Tropical<C>` algebra it points at lives in
  `scalar/tropical.rs`.

## Impartial / outcome analysis

- **`coin_turning.rs`** — `nim_mul_mex`: nim-mult as Conway's Turning-Corners mex
  recurrence (a different *definition* from the algebraic `nim_mul`, proven equal).
  Plus general 1-D coin-turning (`grundy_1d`) and the 2-D Tartan product
  (`tartan_grundy`), with the Tartan/Product theorem verified.
- **`grundy.rs`** — general Sprague–Grundy (normal-play impartial center): `mex`,
  `grundy_graph` (DAG; None on a cycle), closure-based `grundy`. P-position ⟺ g=0;
  SG theorem `g(G+H)=g(G)⊕g(H)` pinned vs Bouton.
- **`kernel.rs`** — normal-play Win/Loss/Draw outcomes of any finite game graph
  (retrograde analysis); `p_positions` = Loss. The interactive route to the open
  question. Plus `scoring_values`: the Milnor minimax `ScoreInterval { left, right }`
  (`i128`) on a DAG — the integer-valued scoring knob.
- **`loopy.rs`** — loopy (cyclic) games, the third escape from XOR-linear P-sets: a
  cyclic rule admits a **Draw** outcome (a genuinely new degree of freedom). Three
  layers: `LoopyGraph` (a thin computable wrapper over `kernel::outcomes` —
  loss/win/draw sets), `loopy_nim_values`/`loopy_nim_values_certified`
  (+ `LoopyNimCertificate`: Draw ⇒ `Side`/∞, else a nimber; exact on an acyclic
  non-Draw subgraph; **not additive over sums when Draw options are present** —
  values are Grundy values of the Draw-deleted subgraph, not Smith/Conway loopy
  nim-values), and the `LoopyValue` stopper catalogue
  (`Zero`/`Star`/`On`/`Off`/`Over`/`Under`/`Dud`, with `outcome` → `PartizanOutcome
  {P,N,L,R,Draw}`, neg/partial order/partial sum). The payoff is
  `loopy_decision_sets`/`loopy_quadric_probe`: read a cyclic rule's Loss-set AND
  Draw-set, each fit by `fit_f2_quadratic`.
- **`misere.rs`** — checked misère-play outcomes (`try_misere_is_n`/`misere_is_p`)
  for finite acyclic impartial games; cycles return `None`. Covers misère Nim vs
  Bouton; the bounded indistinguishability quotient (`misere_quotient`,
  `AbstractGame`, `Quotient`); octal games (`octal_moves`, `octal_misere_quotient`).
  The non-linear route to the open question.
- **`lexicode.rs`** — **Bridge O**, the games ↔ integral edge: greedy binary
  lexicodes `L(n,d)` (Conway–Sloane 1986). `lexicode`/`lexicode_naive`/
  `lexicode_bounded` (+ `LEXICODE_NODE_BUDGET`, an honest backstop → `None`, not a
  silent cap). The greedy step is exactly `mex(Forbidden)` over radius-`(d−1)` Hamming
  balls (`grundy::mex`); linearity is the Sprague–Grundy theorem, *discovered* not
  assumed. Ships the `[7,4,3]` Hamming, `[8,4,4]` extended Hamming, and `[24,12,8]`
  Golay codes as lexicodes, chaining `mex → lexicode → Golay → Construction A → theta`.
  **Claim level:** the degree-1 (solved, linear) side of `OPEN.md` §1 — explicitly does
  NOT touch the open Gold-quadric question; do not cite as progress on it.

## The bridge object

- **`hackenbush.rs`** — red/blue/green Hackenbush: `Hackenbush { edges }` (vertex 0
  is the ground by convention; edges colored by the `Color {Blue, Red, Green}` enum)
  with the `string` stalk constructor, `to_game()` (the universal evaluator),
  `value()` → surreal (blue–red), `grundy()` → nimber (all-green = Nim). The one
  structure tying surreals + nimbers + sign-expansion through a single object.

## Things that look like bugs but are not (games layer)

- **`Game::canonical_string` canonicalizes; `structural_string` does not.**
  `structural_string` is an order-independent fingerprint of the tree *as given* (so
  `(↑−↑).structural_string() ≠ 0`); `canonical_string` reduces first, so it *is* a
  value key. Compare `a.canonical().structural_eq(&b.canonical())` or just compare
  `canonical_string`s.
- **Atomic weight's integer branch is NOT `1 + max_R aw(G^R)`.** It's a predicate
  over `A`'s raw option *games* (`A^R = aw(G^R)+2`) comparing an integer `n` via
  `le`/`fuzzy`, bounded by the *tightest* right option — so it stays correct when an
  option's atomic weight is a fraction (e.g. ½). The naive max-of-integers form
  misreads there. And atomic weight IS additive on all-small games.
- **`nim_mul_mex` is the slow *game* definition (the mex recurrence), for validation
  and small arguments only** — exponential in the argument size. For real computation
  use the algebraic product (`nim_mul`), which it is proven equal to.
- **`nim_moves` takes `&Vec<u128>` (not `&[u128]`) on purpose** (with a `ptr_arg`
  allow): it is passed as a `fn` matching the generic move-generator bound `Fn(&P)`
  with `P = Vec<u128>` in `misere_is_p`/`grundy`, where a `fn(&[u128])` pointer would
  not unify.
- **`Game` stays an acyclic `Arc` tree by construction** (it cannot represent cycles).
  Loopy games are a separate `LoopyGraph` engine; `thermography` is
  finite-game-only (loopy games never freeze to a number).
- **`Pl` does NOT implement `Semiring`.** A `Pl` wall has no representable ∞-wall
  (the tropical `⊕`-identity), so the semiring law-checking lives on `Tropical<C>`
  (which has `Infinity`), not on `Pl`; `Pl` only gets the named wrappers
  `oplus_max`/`oplus_min`/`otimes` it actually uses. Do NOT fake an ∞-wall with empty
  `pts` — that breaks the `self.pts[0]` invariant `value_at` assumes. The
  `pub(crate)` `walls_with` (which internally calls `freeze`/`meeting_temperature`) is
  shared by `thermograph_via_tropical` so it reuses the identical cooling tail; the
  golden thermography tests pin that this sharing is inert.
