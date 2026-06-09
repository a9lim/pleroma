# TODO: Genuine Research Problems

This file is intentionally narrow. It lists directions from the repo audit that
look like genuine new research rather than implementation of known formulas,
standard algorithms, or already-source-pinned theory.

## Natural Gold-quadric game rule

Find, or rule out under a precise naturality condition, a non-tautological game
rule whose P-positions are the zero set `{Q = 0}` of a game-built Gold quadratic
form.

Why this is research:
- The repo already builds the Gold forms and tests several routes, but the
  missing datum is not just code. It is a play rule, or a definition of
  "natural" strong enough to make the question non-ad-hoc.
- Current probes separate frame-blind failures from rules that directly encode
  the quadratic form, leaving an actual definitional gap.

Concrete progress targets:
- Formalize a naturality criterion, such as equivariance, locality,
  encoding-complexity, or basis/framing access.
- Prove a no-go theorem for a class of frame-blind or polar-form-only rules.
- Exhibit a fixed uniform rule, more constrained than an arbitrary lookup game,
  whose P-set or Draw-set is a Gold quadric.

Relevant surfaces:
- `NOTES.md` open-question and naturality-dichotomy sections.
- `experiments/open_question_probe.py`
- `experiments/framing_obstruction.py`
- `examples/interactive_kernel.rs`
- `examples/loopy_quadric.rs`

## Ordinal nim multiplication beyond the verified excess table

Push transfinite nim multiplication beyond the source-verified DiMuro excess
table, especially at the first unverified Kummer carry such as `alpha_47`.

Why this is research:
- Rewriting the current table-driven code to compute the known shape
  `f(u)`, `Q(f(u))`, and hardcode only the finite excess integers is an
  implementation improvement, not research.
- Extending past the verified table is different: the finite excess term has no
  closed form in the cited theorem, and shipping new values would require an
  independent oracle or a new algorithmic proof.

Concrete progress targets:
- Derive or certify finite excess terms beyond the published table.
- Build a verified `u`-th-power/root-search oracle for the transfinite field.
- Prove enough about the search to avoid merely numerological extensions.

Relevant surfaces:
- `src/scalar/big/ordinal/tower.rs`
- `src/scalar/big/ordinal/mod.rs`
- `NOTES.md` section on Lenstra-DiMuro excess elements.

