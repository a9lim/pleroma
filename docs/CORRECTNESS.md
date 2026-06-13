# CORRECTNESS.md (the verification-status ledger)

The verification-status ledger: which shipped claims are **machine-verified**, which
are **source-pinned**, and which are **asserted-but-unproven** — valued like
[`COMPLETENESS.md`](COMPLETENESS.md) — a game value `g` on a pillar blade `e_B` (`e_s`
scalar, `e_c` clifford, `e_f` forms, `e_i` integral, `e_g` games, `e_o` ogham, `e_y`
py). Claim level **interpretation/engineering**: each entry is a status call on the
existing verification surface, checked against the actual oracles, not vibes. Numbers
≈ focused days to close a verification gap; `±n` flags an a9 scope call; `↑` is worth
less than any number but strictly positive; `*n` is real, on-thesis, unscheduled.

The standing verification surface is the baseline this ledger reads against: `cargo
test` (the `proptest` suites `tests/scalar_axioms.rs` and `tests/clifford_axioms.rs`,
the `associativity_*` oracles, and `general_product_reproduces_reduce_word_when_a_empty`),
the adversarial stdlib harnesses `experiments/echo_solver.py` and
`experiments/linking_game.py`, and the source-pinned finite tables inventoried in
[`TABLES.md`](TABLES.md). Its aesthetic sibling — structural/stylistic findings rather
than soundness — is [`CONSISTENCY.md`](CONSISTENCY.md).

---

Currently empty. The 2026-06-10 correctness sweep is recorded against the code and the
draft notes (`experiments/audit/`, `writeups/goldarf.tex`, `writeups/excess.tex`); the
next verification audit — a claim resting on bounded observation, an assertion no
oracle yet pins, or a harness not yet wired into CI — lands its findings here.
