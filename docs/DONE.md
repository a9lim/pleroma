# Cross-pillar bridges — DONE (the go-forward ledger)

The running ledger of cross-pillar work **completed from here on**.

The cross-pillar bridge-building era (bridges **A–O** plus **K** — lattice/Clifford/
Brauer–Wall, the char-2 Arf classifier, Frobenius outermorphisms, the transfinite
Clifford engine, theta/modular forms, Construction-A codes, the Weil representation, the
rational and full-`ℚ/ℤ` Brauer invariants, Newton polygons, the Brown invariant, the
unification pass, lexicodes) closed with every non-deferred bridge shipped, as did the
ogham 1.x–2.x language work and the transfinite-excess thread. The working-notes summary
of all of it is in the `AGENTS.md` files (root + per-pillar); the historical entry-level
ledger is in git history.

What remains unbuilt is tracked in the two buildable ledgers —
[`COMPLETENESS.md`](COMPLETENESS.md) (completing symmetries and connections already in
the code) and [`CONTINUATIONS.md`](CONTINUATIONS.md) (genuinely new features), each
carrying its slice of the deferred stars `*1`/`*2`/`*4`/`*8`; genuine open problems
stay in [`OPEN.md`](OPEN.md), loopy-valued: `tis`/`tisn`, `on`/`off`, `over`/`under`
(the old numerals §1–§4 survive as aliases).

## How to use this ledger

When a new piece of cross-pillar work lands, add a short entry here:

```
## completed items

### <date>: `<name>`
**Summary:** <one-line what-it-connects>
**Pillars:** … ↔ …    **Claim level:** standard math / implemented-and-tested / …
- surface: the functions/types that shipped
- oracles: the tests that pin it
- boundaries: the honest non-claims
```

Fold the one-line structural fact into the relevant `AGENTS.md`; keep any longer
derivation alongside the code or in a `writeups/` note.

## completed items

### 2026-06-12: `padic-symbols`
**Summary:** Conway-Sloane `p`-adic genus symbols now give exact integral-lattice
genus comparison, with the canonical 2-adic train/compartment/oddity reduction
exposed on the Rust and Python `Genus` surface.
**Pillars:** integral    **Claim level:** standard math, implemented and tested
- surface: `Genus::of`, `Genus::symbol_at`, `Genus::canonical_symbol_at`,
  `are_in_same_genus`, and Python `Genus.canonical_symbol_at`.
- oracles: odd-prime determinant-square-class symbols, Sage/Allcock-style 2-adic
  canonical-symbol examples, random unimodular-congruence invariance, `Z^8` vs
  `E8`, `E8⊕E8` vs `D16+`, and Nikulin/discriminant-form agreement across the
  ADE zoo and Milnor pair.
- boundaries: odd-lattice discriminant forms, full spinor-genus computation, and
  level-`N` theta machinery stay on their separate docket items.
