"""The polar form is game-built: B realized by Turning-Corners (tartan) products.

open_question_probe.py showed the obstruction to {Q_a=0} being a normal-play
P-set is *exactly* the polar form B of the Gold form. This script closes the
"is B game-realizable?" half of that, concretely: it recomputes B(e_i,e_j) using
ONLY the game definition of the nim-product — Conway's Turning Corners,
`nim_mul_mex`, which is the tartan square of the 1-D game g(n)=n (verified in
games.rs) — together with the trace (iterated game-squaring + XOR).

The Gold form is Q_a(v) = Tr(v^{1+2^a}); its polar form bilinearises to

    B(u, v) = Tr( u ⊗ v^{2^a}  ⊕  v ⊗ u^{2^a} ),

every ⊗ here a Turning-Corners product and every squaring the diagonal of one.
So B is a composite of coin-turning games. We check the game-built B agrees with
the algebraic polar form on the bit basis, over the tiny nim-fields F₄ and F₁₆
(the only sizes where the exponential `nim_mul_mex` is feasible — exactly the
regime AGENTS.md reserves it for).

Upshot: of the three layers a P-position game would need — linear (Grundy/XOR),
the bilinear coupling B, and a quadratic play rule — the first two are now both
realized in code from actual games. Only the quadratic play rule is missing.
"""

import pleroma as pl

from common import gold


def nim_sq_game(x: int) -> int:
    return pl.nim_mul_mex(x, x)          # the diagonal of Turning Corners


def frob_game(x: int, a: int) -> int:
    for _ in range(a):                   # Frobenius^a = a game-squarings
        x = nim_sq_game(x)
    return x


def trace_game(x: int, m: int) -> int:
    acc, t = x, x
    for _ in range(m - 1):               # Tr = iterated squaring + XOR
        t = nim_sq_game(t)
        acc ^= t
    return acc


def b_game(i: int, j: int, a: int, m: int) -> int:
    """B(e_i, e_j) via Turning-Corners products and the trace."""
    ei, ej = 1 << i, 1 << j
    term = pl.nim_mul_mex(ei, frob_game(ej, a)) ^ pl.nim_mul_mex(ej, frob_game(ei, a))
    return trace_game(term, m)


def b_form(i: int, j: int, a: int, m: int) -> int:
    """B(e_i, e_j) = Q(e_i⊕e_j) ⊕ Q(e_i) ⊕ Q(e_j), the algebraic polar form."""
    return gold((1 << i) ^ (1 << j), a, m) ^ gold(1 << i, a, m) ^ gold(1 << j, a, m)


if __name__ == "__main__":
    hdr = f"{'field':>7} {'a':>2} {'pairs':>6} {'B game=form?':>13}"
    print(hdr)
    print("-" * len(hdr))
    all_ok = True
    for k in (1, 2):                     # F2^2, F2^4  (nim_mul_mex stays tiny)
        m = 1 << k
        for a in range(1, k + 1):
            pairs = 0
            ok = True
            for i in range(m):
                for j in range(m):
                    if b_game(i, j, a, m) != b_form(i, j, a, m):
                        ok = False
                    pairs += 1
            all_ok = all_ok and ok
            print(f"  F2^{m:<3} {a:>2} {pairs:>6} {str(ok):>13}")
    print("-" * len(hdr))
    print(f"The polar form B is computed by Turning-Corners (tartan) games: {all_ok}")
    print("→ linear (Grundy) and bilinear (B) layers are both game-built;")
    print("  the open piece is a play rule reading out the quadratic Q.")
