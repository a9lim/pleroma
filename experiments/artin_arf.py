"""Artin–Schreier ↔ Arf: the field trace does double duty.

The trace Tr_{F_{2^m}/F2}(x) = Σ_{i<m} x^{2^i} shows up in two apparently
separate places in this repo:

 (1) it pushes the Arf invariant of a char-2 quadratic form down to F2 (the
     canonical k/℘(k) ≅ F2, used in arf.rs);
 (2) it is the exact obstruction to solving the Artin–Schreier equation
     y² + y = c   (solvable in F_{2^m}  ⟺  Tr(c) = 0).

The game-built Gold form is Q_a(v) = Tr(v^{1+2^a}). So its zero set is literally
the set of v whose "twisted square" c = v^{1+2^a} admits an Artin–Schreier root:

    Q_a(v) = 0   ⟺   Tr(v^{1+2^a}) = 0   ⟺   y² + y = v^{1+2^a} is solvable.

Hence the Arf win-bias #{Q_a = 0} (see arf_win_bias.py) is an *Artin–Schreier
solvability count*. This script checks, exhaustively on the nim-subfields
F_{2^{2^k}}:

  • nim_sqrt is the inverse Frobenius (√x squared is x);
  • the Artin–Schreier theorem: solvable ⟺ Tr(c)=0, and the solver returns a
    genuine root y with y²+y = c;
  • the link: Gold-form-zero(v) ⟺ the twisted square v^{1+2^a} is AS-solvable.
"""

import pleroma as pl

from common import gold


def twisted_square(v: int, a: int) -> int:
    """c = v^{1+2^a} = v ⊗ Frobenius^a(v)."""
    x = pl.Nimber(v)
    g = x
    for _ in range(a):
        g = g * g
    return (x * g).value


if __name__ == "__main__":
    # nim_sqrt is the inverse Frobenius across the 64-bit field (sampled).
    assert all(
        pl.Nimber(pl.nim_sqrt(x)) * pl.Nimber(pl.nim_sqrt(x)) == pl.Nimber(x)
        for x in range(512)
    ), "nim_sqrt is not the inverse Frobenius"

    hdr = f"{'field':>7} {'a':>2} {'#Q=0':>6} {'#AS-solv':>9} {'link?':>6} {'AS-thm?':>8}"
    print(hdr)
    print("-" * len(hdr))
    all_ok = True
    for k in range(1, 5):                       # F2^2, F2^4, F2^8, F2^16
        m = 1 << k
        # the Artin–Schreier theorem on the WHOLE field, once per m
        as_theorem = True
        for c in range(1 << m):
            solv = pl.nim_is_artin_schreier_solvable(c, m)
            if solv != (pl.nim_trace(c, m) == 0):
                as_theorem = False
            y = pl.nim_solve_artin_schreier(c, m)
            if solv:
                yy = pl.Nimber(y)
                if (yy * yy + yy).value != c:
                    as_theorem = False
            elif y is not None:
                as_theorem = False

        for a in range(1, k + 1):
            zeros_Q = 0
            solv_AS = 0
            link = True
            for v in range(1 << m):
                q0 = gold(v, a, m) == 0
                as0 = pl.nim_is_artin_schreier_solvable(twisted_square(v, a), m)
                zeros_Q += q0
                solv_AS += as0
                if q0 != as0:                   # the headline equivalence
                    link = False
            ok = link and as_theorem and zeros_Q == solv_AS
            all_ok = all_ok and ok
            print(f"  F2^{m:<3} {a:>2} {zeros_Q:>6} {solv_AS:>9} "
                  f"{('✓' if link else '✗'):>6} {('✓' if as_theorem else '✗'):>8}")
    print("-" * len(hdr))
    print(f"Gold-zero(v) ⟺ v^(1+2^a) is Artin–Schreier-solvable; one trace, both roles: {all_ok}")
