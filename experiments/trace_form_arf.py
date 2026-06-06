"""Empirical probe (extension direction 1): quadratic structure of the nim-field.

The nim-field On₂ is the Grundy-value field of coin-turning games, so it is the
natural substrate to look for quadratic (Arf-bearing) structure. Grundy values
are linear and coin-turning products are bilinear; a quadratic form needs an
independent diagonal. The natural quadratic forms on a char-2 field are

    Q_a(x) = Tr(x · x^{2^a}) = Tr(x^{1+2^a})          (g = Frobenius^a is additive)

— exactly the forms studied in coding theory / Kloosterman-sum land, so their
Arf invariants and ranks have known patterns to check against.

We build each Q_a as an F₂ quadratic form in the bit-basis e_i = 2^i of
F_{2^m} (m = 2^k) and read off its Arf invariant with pleroma's classifier.
Everything here runs on top of the shipped library — this is the tool used for
its own research question.
"""

import pleroma as pl


def nim_trace(x: int, m: int) -> int:
    """Tr_{F_{2^m}/F_2}(x) = x + x^2 + ... + x^{2^{m-1}}, in {0,1}."""
    acc = pl.Nimber(x)
    t = pl.Nimber(x)
    for _ in range(m - 1):
        t = t * t
        acc = acc + t
    assert acc.value in (0, 1), f"trace not in F2: {acc.value}"
    return acc.value


def frob(x: pl.Nimber, a: int) -> pl.Nimber:
    """Frobenius^a: x -> x^{2^a}."""
    for _ in range(a):
        x = x * x
    return x


def trace_form(m: int, a: int):
    """Arf result of Q_a(x) = Tr(x^{1+2^a}) on F_{2^m}, in the bit-basis."""
    n = m

    def Q(xi: int) -> int:
        x = pl.Nimber(xi)
        return nim_trace((x * frob(x, a)).value, m)

    def B(xi: int, yj: int) -> int:
        x, y = pl.Nimber(xi), pl.Nimber(yj)
        return nim_trace((x * frob(y, a) + y * frob(x, a)).value, m)

    q = [pl.Nimber(Q(1 << i)) for i in range(n)]
    b = {}
    for i in range(n):
        for j in range(i + 1, n):
            if B(1 << i, 1 << j):
                b[(i, j)] = pl.Nimber(1)
    return pl.arf_invariant(pl.NimberAlgebra(q=q, b=b))


from math import gcd

if __name__ == "__main__":
    hdr = f"{'field':>7} {'a':>2} {'Q':>10} {'Arf':>3} {'type':>4} {'rank':>4} {'rad':>3} {'2·gcd':>5} {'ok?':>3}"
    print(hdr)
    print("-" * len(hdr))
    all_ok = True
    for k in range(1, 6):           # m = 2,4,8,16,32
        m = 1 << k
        for a in range(1, k + 1):
            r = trace_form(m, a)
            exp = 1 + (1 << a)
            predicted_rad = 2 * gcd(a, m)   # Gold-function rank: rank = m - 2·gcd(a,m)
            ok = (r.radical_dim == predicted_rad) and (r.rank == m - predicted_rad)
            all_ok = all_ok and ok
            print(f"  F2^{m:<3} {a:>2} {('Tr(x^%d)' % exp):>10} "
                  f"{r.arf:>3} {r.o_type:>4} {r.rank:>4} {r.radical_dim:>3} {predicted_rad:>5} {('✓' if ok else '✗'):>3}")
    print("-" * len(hdr))
    print(f"rank formula rank = m - 2·gcd(a,m) matches the Arf classifier: {all_ok}")
