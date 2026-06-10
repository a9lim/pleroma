# AGENTS.md — `src/games/`

The PILLAR of combinatorial game theory — the second column of the project,
mostly independent of the scalar/Clifford stack (the bridge is the
number/nimber subclasses, where Conway multiplication is defined). Games under
disjunctive sum are an abelian GROUP, not a ring; that constraint is *why* the
Clifford story lives on the scalar backends and not on all games.

> Read root `OPEN.md` before touching `coin_turning.rs`, `kernel.rs`, `misere.rs`, or
> the example probes — they feed the open play-semantics question.

`mod.rs` re-exports every module below flat.

Fixed-width game payloads use `u128`/`i128`: finite nim heaps, octal codes, Grundy
values, scoring integers, and quotient counters follow the same repo-wide width
contract. `usize` remains for graph nodes, option indices, and collection lengths.

## Values & arithmetic

- **`partizan.rs`** — short partizan games (sum/neg/order/birthday/is_number) + the
  CANONICAL FORM (dominated/reversible reduction; `structural_string` vs
  `canonical_string` — the latter canonicalizes, a value key) + the game↔surreal
  bridge (`number_value`/`from_surreal`, numbers only). Also `Game::ordinal_sum`
  (G:H — Hackenbush strings are these), `Game::nim_heap` (⋆n), `Game::is_all_small`.
- **`number_game.rs`** — transfinite NUMBER games (ω, ε) carried by their Surreal
  value — value/birthday/sum/cmp delegate to surreal, no infinite option tree. Plus
  the FULL transfinite round trip via sign_expansion/from_sign_expansion (the run-
  length sign expansion is the finite encoding of the infinite {L|R} tree).
- **`nimber_game.rs`** — the CHAR-2 MIRROR of `number_game.rs`: transfinite NIMBER
  games (Nim heaps `⋆α`) carried by their `Ordinal` (On₂) Grundy value — grundy/
  add(=nim-add XOR)/cmp/to_finite_game delegate to `Ordinal`, no infinite option
  tree. `neg` is the identity (char 2: every impartial game is self-inverse);
  `turning_corners` is the nim-product (Conway's coin game, `ω³=2`); `None` only
  when a Kummer carry needs a prime past the verified table or at `≥ ω^(ω^ω)`.
  Closes the `No ↔ On₂` symmetry at the games layer (it had only lived at the scalar
  layer via the shared CNF core, and now reaches Clifford through `Scalar for
  Ordinal` inside the checked Kummer boundary). Bound to Python as `NimberGame`.
- **`game_exterior.rs`** — the exterior algebra of the GAME group: Λ over ℤ on game
  generators (living on all of game-world, incl. non-numbers ⋆/↑ — needs only the
  ℤ-module structure). `GameExterior` (free Grassmann engine quotiented by integer
  game relations such as 2⋆=0) + `GameRelation`; lattice normalization in
  `linalg/integer.rs`.

## Temperature theory

- **`thermography.rs`** — the thermograph of a short game: left/right scaffolds,
  stops, cooling (`cooled_stops`), temperature, and mean (mast) value.
- **`atomic_weight.rs`** — atomic weight of ALL-SMALL games (finishes thermography):
  the two-ahead rule (Siegel Constructive Atomic Weight; Larsson–Nowakowski
  arXiv:2007.03949 Thm 10). `aw` IS additive on all-small games.
- **`piecewise.rs`** — `Pl`: exact rational piecewise-linear wall arithmetic used by
  thermography.
- **`tropical_thermography.rs`** — names the latent tropical structure in
  thermography and machine-checks it. The option folds are tropical `⊕` in DUAL
  semirings — the left wall a `(max,+)` fold over the Left options' right walls, the
  right wall a `(min,+)` fold over the Right options' left walls — and cooling is
  tropical `⊗`. `Pl::oplus_max`/`oplus_min`/`otimes` name the wall operations;
  `thermograph_via_tropical` is a parallel recursion pinned EQUAL to
  `thermography::thermograph` (it reuses the IDENTICAL freeze/meeting-temperature
  tail — it only renames the folds, it does not reimplement cooling). The
  `Semiring`/`Tropical<C>` algebra it points at lives in `scalar/tropical.rs`.

## Impartial / outcome analysis

- **`coin_turning.rs`** — `nim_mul_mex`: nim-mult as Conway's Turning-Corners mex
  recurrence (== algebraic `nim_mul`). Plus general 1-D coin-turning (`grundy_1d`)
  and the 2-D Tartan product (`tartan_grundy`), with the Tartan/Product theorem
  verified. (Distinct from `coin_turning` ≠ the algebraic `nim_mul`.)
- **`grundy.rs`** — general Sprague–Grundy (normal-play impartial center): `mex`,
  `grundy_graph` (DAG; None on a cycle), closure-based `grundy`. P-position ⟺ g=0;
  SG theorem `g(G+H)=g(G)⊕g(H)` pinned vs Bouton.
- **`kernel.rs`** — normal-play Win/Loss/Draw outcomes of any finite game graph
  (retrograde analysis); P-positions = Loss. The interactive route to the open
  question. Plus `scoring_values`: the Milnor minimax interval (left, right) on a
  DAG — the integer-valued scoring knob.
- **`loopy.rs`** — loopy (cyclic) games, the third escape from XOR-linear P-sets: a
  cyclic rule admits a **Draw** outcome (a genuinely new degree of freedom). Three
  layers: `LoopyGraph` (a thin computable wrapper over `kernel::outcomes` —
  loss/win/draw sets), `loopy_nim_values` (Draw ⇒ `Side`/∞, else a nimber; exact on
  an acyclic non-Draw subgraph), and the `LoopyValue` stopper catalogue
  (on/off/over/under/dud with outcome/neg/partial order/partial sum). The payoff is
  `loopy_decision_sets`/`loopy_quadric_probe`: read a cyclic rule's Loss-set AND
  Draw-set, each fit by `fit_f2_quadratic`.
- **`misere.rs`** — checked misère-play outcomes (`try_misere_is_n`/`misere_is_p`)
  for finite acyclic impartial games; cycles return `None`. Covers misère Nim vs
  Bouton; the bounded indistinguishability quotient (`misere_quotient`,
  `AbstractGame`, `Quotient`); octal games (`octal_moves`,
  `octal_misere_quotient`). The non-linear route to the open question.

## The bridge object

- **`hackenbush.rs`** — red/blue/green Hackenbush: `Hackenbush{edges, ground=0}`,
  `to_game()` (the universal evaluator), `value()` → surreal (blue–red), `grundy()`
  → nimber (all-green = Nim). The one structure tying surreals + nimbers + sign-
  expansion through a single object.

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
  misreads there (Codex-caught). And atomic weight IS additive on all-small games.
- **`nim_mul_mex` is the slow *game* definition (the mex recurrence), for validation
  and small arguments only** — exponential in the argument size. For real
  computation use the algebraic product (`nim_mul`), which it is proven equal to.
- **`nim_moves` takes `&Vec<u128>` (not `&[u128]`) on purpose** (with a `ptr_arg`
  allow): it is passed as a `fn` matching the generic move-generator bound `Fn(&P)`
  with `P = Vec<u128>` in `misere_is_p`/`grundy`, where a `fn(&[u128])` pointer would
  not unify.
- **`Game` stays an acyclic `Arc` tree by construction** (it cannot represent
  cycles). Loopy games are a separate `LoopyGraph` engine; `thermography` is finite-
  game-only (loopy games never freeze to a number).
- **`Pl` does NOT implement `Semiring`, and `add_pl` + the `pub(crate)` widening in
  `thermography.rs` are deliberate.** A `Pl` wall has no representable ∞-wall (the
  tropical `⊕`-identity), so the semiring law-checking lives on `Tropical<C>` (which
  has `Infinity`), not on `Pl`; `Pl` only gets the named wrappers
  `oplus_max`/`oplus_min`/`otimes` it actually uses. Do NOT fake an ∞-wall with empty
  `pts` — that breaks the `self.pts[0]` invariant `value_at` assumes. `add_pl` (in
  `piecewise.rs`) is the additive twin of `sub_pl` naming the `⊗`; `freeze` and
  `meeting_temperature` were widened from private to `pub(crate)` (bodies untouched)
  only so `thermograph_via_tropical` reuses the identical cooling tail — the golden
  thermography tests are the proof those visibility edits are inert.
