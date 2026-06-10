"""Shared Gold-form helpers for the experiment scripts."""

import ogdoad as pl


def frob(x: pl.Nimber, a: int) -> pl.Nimber:
    """Frobenius^a: x -> x^(2^a)."""
    for _ in range(a):
        x = x * x
    return x


def nim_trace(x: int, m: int) -> int:
    """Trace from F_{2^m} to F_2, returned as 0 or 1."""
    acc = pl.Nimber(x)
    t = pl.Nimber(x)
    for _ in range(m - 1):
        t = t * t
        acc = acc + t
    assert acc.value in (0, 1), f"trace not in F2: {acc.value}"
    return acc.value


def gold(v: int, a: int, m: int) -> int:
    """Gold form Q_a(v) = Tr(v^(1+2^a)) over F_{2^m}."""
    x = pl.Nimber(v)
    return nim_trace((x * frob(x, a)).value, m)


def polar(u: int, v: int, a: int, m: int) -> int:
    """Polar form B(u,v) = Q(u+v) + Q(u) + Q(v)."""
    return gold(u ^ v, a, m) ^ gold(u, a, m) ^ gold(v, a, m)
