"""The Gold form, built entirely out of game operations.

src/games.rs verifies that Conway's Turning-Corners game *is* nim-multiplication:
pl.nim_mul_mex(x, y) == Nimber(x) * Nimber(y). Three operations are therefore
game-realizable:

    ⊗   nim-product            =  Turning-Corners Grundy value
    □   Frobenius v ↦ v²        =  v ⊗ v   (the diagonal of Turning Corners)
    ⊕   XOR                     =  disjunctive sum of single-coin positions

and the field trace Tr(x) = x ⊕ x² ⊕ … ⊕ x^{2^{m-1}} is just iterated □ and ⊕.

So the Gold quadratic form Q_a(v) = Tr(v ⊗ v^{2^a}) is a *composite of game
operations* on a position's nimber value v. We confirm two things:

  1. the bridge: the Turning-Corners product equals the algebraic product;
  2. the form computed with literal game products (nim_mul_mex) equals the form
     computed with the fast algebraic product, on the small fields where the
     mex recurrence is tractable.

Then the field-wide form (used by trace_form_arf.py for the Arf invariant) is
computed with the fast product, which (1) shows is the same operation.

Open: Q_a is a derived form, not one position's Grundy value, and the
play-meaning of its Arf invariant is not established. The form is made of games;
the invariant's game-semantics is the remaining gap.
"""

import pleroma as pl


def bridge_holds(limit: int = 24) -> bool:
    """Turning-Corners (game) product == algebraic nim-product."""
    return all(
        pl.nim_mul_mex(x, y) == (pl.Nimber(x) * pl.Nimber(y)).value
        for x in range(limit)
        for y in range(limit)
    )


def gold_fast(v: int, a: int, m: int) -> int:
    """Q_a(v) = Tr(v ⊗ v^{2^a}) via the fast algebraic product (= game product)."""
    x = pl.Nimber(v)
    g = x
    for _ in range(a):
        g = g * g                       # Frobenius^a
    s = x * g                           # v ⊗ v^{2^a}
    acc, t = s, s
    for _ in range(m - 1):
        t = t * t
        acc = acc + t                   # trace = iterated square + XOR
    return acc.value


def gold_literal(v: int, a: int, m: int) -> int:
    """The same form using ONLY literal Turning-Corners products (small fields)."""
    sq = lambda z: pl.nim_mul_mex(z, z)
    g = v
    for _ in range(a):
        g = sq(g)
    s = pl.nim_mul_mex(v, g)
    acc, t = s, s
    for _ in range(m - 1):
        t = sq(t)
        acc ^= t
    return acc


if __name__ == "__main__":
    print("bridge  Turning-Corners product == algebraic nim-product (x,y<24):",
          bridge_holds())
    print()
    print("literal game ops reproduce the Gold form where the mex recurrence is feasible:")
    all_ok = True
    for k in (1, 2):                    # F2^2, F2^4
        m = 1 << k
        for a in range(1, k + 1):
            ok = all(gold_literal(v, a, m) == gold_fast(v, a, m) for v in range(1 << m))
            all_ok = all_ok and ok
            print(f"  F2^{m}  a={a}:  {'✓' if ok else '✗'}")
    print()
    print(f"=> the Gold / Arf-bearing form is a composite of game operations: {all_ok}")
    print("   (Arf invariants of these forms: see experiments/trace_form_arf.py)")
