# Cross-pillar work — CONTINUATIONS (genuinely new features)

The ledger of buildable items that **extend ogdoad past what it currently covers** —
new directions and features, not the completion of a connection already in the code.
The two exemplars are the **ogham** language work (a whole sub-language growing toward
recursion + games) and the **char-`p` Drinfeld/Carlitz mirror** (a candidate second
headline pillar). Items that round out an existing symmetry or bridge — most of the
standing content — live in [`COMPLETENESS.md`](COMPLETENESS.md); newly completed work
goes in [`DONE.md`](DONE.md); genuine research questions in [`OPEN.md`](OPEN.md).

Claim-level discipline (`AGENTS.md` → "Claim levels and non-claims") applies to every
item: each is **standard math** or **engineering** when built — not a new theorem.

Items are valued exactly as in [`COMPLETENESS.md`](COMPLETENESS.md) — a game value `g`
on a pillar blade `e_B` (the "How items are valued" legend is canonical there). Numbers
are cold/buildable, `±n` switches are a9's scope call first, `↑` ups are infinitesimal,
`*n` stars are deferred not-yet-numbers; reference items by **slug**.

---

## numbers — ogham (the language)

### 2·e_o: `ogham-reflect`
**The consolidation pass before release** (a9, 2026-06-12: "consolidate it
into the best version of itself before release") — plays after `ogham 3.0`
(star `*8` below), before any 4.0 design. Scope: (1) rewrite the spec §1 identity:
the principles describe a v1 calculator, and by 3.0 the honest description
is the **lisp-for-games** — the value-rich/computation-thin inversion,
Conway's ontology as the data model, the construct↔math coincidences on
record (four-way relations = outcome classes, `=:` = loopy definition, the
lazy trio = play-one-branch); (2) fold the §17–§19 delta sections into the
main spec body so the language reads as one contract, and merge/reorganize
the conformance corpus; (3) a CONSISTENCY.md-style audit of `src/ogham/` after
three builds of growth — naming, error taxonomy, dispatch-enum shape, REPL
UX; (4) release scoping, **a9's call**: ogham as ogdoad's front door vs an
`ogham` crate re-exporting the core, README/writeup, the public name. Also
worth an hour inside this pass: a CGScript/CGSuite comparison read, for
ideas and for honest differentiation. The refactor is licensed; the
identity questions are the point.

---

## stars (deferred — the not-yet-numbers, confused with zero)

The star numbers are one shared nim-sum scheme across both buildable ledgers; the
sibling stars `*1` (spinor genus) and `*4` (the wild local symbol) live in
[`COMPLETENESS.md`](COMPLETENESS.md).

### *2: `the char-p Drinfeld/Carlitz mirror of the integral pillar` (large)

The entire `integral/` wing — even-unimodular `ℤ`-lattices, `θ`-series,
`M_*(SL₂ℤ) = ℂ[E₄, E₆]`, Construction-A codes, Leech — is char 0. The project already
ships **exact** `F_q[t] ⊂ F_q(t)`, the char-`p` global field, whose arithmetic carries a
complete mirror:

- the **Carlitz module** `C_t(x) = t·x + x^q` is the char-`p` analogue of `exp` / the
  lattice exponential; the mirror of `E₄, E₆` are **Drinfeld modular forms** for
  `GL₂(F_q[t])`, with Goss `ζ`-values mirroring the Eisenstein constants;
- rank-`r` `F_q[t]`-lattices mirror even-unimodular `ℤ`-lattices and their reduction
  theory;
- **Goppa / algebraic-geometry codes** from function fields tie straight back into the
  existing `codes.rs` Construction-A machinery — the same code↔lattice seam in char `p`.

This is the `No ↔ On₂` / char-0 ↔ char-2 move applied to the richest pillar — the most
on-thesis possible "new structure." But it is a genuine new wing (Drinfeld modules, the
Carlitz exponential, rank-`r` reduction theory): weeks of specialized work, worth starting
only as a *second headline pillar* rather than a task. References: Goss, *Basic Structures
of Function Field Arithmetic*; Gekeler, Drinfeld modular forms; Goppa / AG codes.

### *8: `ogham 3.0 — recursion + games`

The `docs/ogham/ogham.md` §19 stub (2026-06-12; the predecessor `*8` — ogham 2.x
functions/programs — was converted to the numbered `ogham-2.0`/`ogham-2.1`
entries, recorded in the DONE ledger now in git history when their sketches landed). The semantic
break and the telos: **totality traded for
attributable partiality** (fuel + `E_Depth`, `:depth`), `=:` fixpoint
bindings (μ — `:=` captures the past, `=:` is an equation the name
satisfies; a9's notation), local `=:` in bodies, and the `game` world —
`{L|R}` as ogham's cons cell: game forms as the recursive data constructor,
CGT's full four-way order live, `⋅` rejected (group-not-ring made an
evaluator fact), Index-based option access, `grundy` via `=:` as the
acceptance example — and **Element-`=:` as loopy games** (§19.4, folded in
at a9's call 2026-06-12: guarded fixpoint equations on game forms *are*
coinductive definitions — `dud =: {dud | dud}` — with outcomes from
`games/loopy/`; the construct and the math object coincide again). Owed to
the real sketch: mutual-recursion groups, fuel default, up/down naming (the
`↑` glyph collision), the loopy stopper/sum envelope, game-form display, and
the sequence-sort/HOF gate. Held as a star until §19 grows into a real sketch.
Nimbered `*8`: every smaller name is a nim-sum of the shipped
stars (`*3 = *1 + *2`, …, `*7 = *1 + *2 + *4`).
