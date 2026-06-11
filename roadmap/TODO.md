# Cross-pillar bridges — TODO (deferred)

Every bridge that was *explicitly on the build order* is done — the full record is in
[`roadmap/CODA.md`](CODA.md), and newly completed work goes in the
[`roadmap/DONE.md`](DONE.md) ledger. What's left are the two **deferred** bridges: real,
standard-math, on-thesis, but each a larger build not slated into the current order.
Fittingly, they're the *star* values — nimbered `*1` and `*2`, the not-yet-numbers of the
roadmap.

Claim-level discipline (`AGENTS.md` → "Claim levels and non-claims") still applies: both
are **standard math made computational** when built — not new theorems. Genuine open
problems (no known answer) live in [`OPEN.md`](../OPEN.md), not here.

## `*1` — spinor genus (was Bridge G)

Refine `genus → spinor genus → isometry class` via the spinor norm (Eichler;
Cassels–Hall). `clifford/spinor_norm.rs` is the right primitive in spirit, but the full
bridge is **not buildable from the current surface**: `spinor_norm` computes one versor's
norm, whereas the spinor genus needs the local spinor-norm *images* `θ(O(L ⊗ ℤ_p))` at
every prime, adelic class-group bookkeeping, and the proper/improper class distinction.

The one cheap, honest piece is **Eichler's theorem** as a documented predicate —
*indefinite, rank ≥ 3* ⇒ spinor genus = isometry class — which would let `Genus` upgrade
to a class statement in exactly that regime. The full definite-lattice computation is the
larger build; it sits adjacent to the roadmap, not inside it.

## `*2` — the char-`p` Drinfeld/Carlitz mirror of the integral pillar (large)

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
