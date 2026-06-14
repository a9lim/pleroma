"""Small independent probes for Lenstra excess in the ordinal nim tower.

This is not a replacement for CGSuite's calculator. It is a deliberately small
term-algebra oracle for the first cases where the excess search is subtle.

For a nonzero finite-field element beta in F_{2^E}, beta has no p-th root exactly
when p divides the multiplicative order of beta. Lenstra excess can therefore be
viewed as the least finite translate m for which

    p | ord(kappa_{f(p)} + m).

The cases below independently reproduce:

* the singleton-odd lower bound: m=0 still has a root for p=7;
* the first 2*3^k exception: p=19 needs m=4, not m=1;
* the fact that Q alone is not enough: Q={9} gives p=19 -> 4 but p=73 -> 1;
* a couple of larger singleton witnesses;
* a fixed-base p=47 test, using only lower verified rows, that certifies
  alpha_47 = omega^(omega^7)+1 for the Rust tower;
* a deeper fixed-base dependency rehearsal certifying m_179 = 1 in the
  E=19580 component field (run with --deep).
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from functools import cache


# Enough Lenstra data to build the component fields used by this probe. These are
# the small source-table values, not values inferred from this script.
Q_SET: dict[int, tuple[int, ...]] = {
    2: (),
    3: (2,),
    5: (4,),
    7: (3,),
    11: (5,),
    13: (3, 4),
    17: (8,),
    19: (9,),
    23: (11,),
    # Target rows used only to build kappa_{f(p)} + m.
    47: (23,),
    53: (13,),
    73: (9,),
    # Dependency chain for the first OEIS-unknown row p=719:
    # 719 -> 359 -> 179 -> 89 -> 11 -> 5 -> finite.
    # The p=89 row is fast; p=179 is re-certified by the --deep path below.
    89: (11,),
    179: (89,),
    359: (179,),
}

EXCESS: dict[int, int] = {
    2: 0,
    3: 0,
    5: 0,
    7: 1,
    11: 1,
    13: 0,
    17: 0,
    19: 4,
    23: 1,
    89: 1,
    179: 1,
}

# Factorizations of 2^E - 1 for the small component fields exercised below.
GROUP_ORDER_FACTORS: dict[int, dict[int, int]] = {
    6: {3: 2, 7: 1},
    18: {3: 3, 7: 1, 19: 1, 73: 1},
    36: {3: 3, 5: 1, 7: 1, 13: 1, 19: 1, 37: 1, 73: 1, 109: 1},
    156: {
        3: 2,
        5: 1,
        7: 1,
        13: 2,
        53: 1,
        79: 1,
        157: 1,
        313: 1,
        1249: 1,
        1613: 1,
        2731: 1,
        3121: 1,
        8191: 1,
        21841: 1,
        121369: 1,
        22366891: 1,
    },
    220: {
        3: 1,
        5: 2,
        11: 2,
        23: 1,
        31: 1,
        41: 1,
        89: 1,
        397: 1,
        683: 1,
        881: 1,
        2113: 1,
        2971: 1,
        3191: 1,
        201961: 1,
        48912491: 1,
        415878438361: 1,
        3630105520141: 1,
    },
}


def prime_power(q: int) -> tuple[int, int]:
    for p in range(2, q + 1):
        if q % p != 0:
            continue
        n = 0
        x = q
        while x % p == 0:
            n += 1
            x //= p
        if x == 1:
            return p, n
        raise ValueError(f"{q} is not a prime power")
    raise ValueError(f"{q} is not a prime power")


def factor_product(factors: dict[int, int]) -> int:
    out = 1
    for p, e in factors.items():
        out *= p**e
    return out


def finite_summand(q_set: tuple[int, ...], excess: int) -> int:
    even_q = next((q for q in q_set if q % 2 == 0), None)
    base = 0 if even_q is None else 1 << (even_q // 2)
    return base + excess


def finite_components_for(n: int) -> set[int]:
    if n <= 1:
        return set()
    # To represent all finite terms up to the top bit of n, include
    # kappa_2, kappa_4, ..., kappa_{2^r}.
    return {1 << (k + 1) for k in range((n.bit_length() - 1).bit_length())}


@cache
def primitive_components(q: int) -> frozenset[int]:
    p, _ = prime_power(q)
    components: set[int] = set()
    pn = 1
    while pn < q:
        pn *= p
        components.add(pn)
    for q1 in Q_SET[p]:
        components.update(primitive_components(q1))
    components.update(finite_components_for(finite_summand(Q_SET[p], EXCESS[p])))
    return frozenset(components)


def toggle_terms(dst: set[int], terms: frozenset[int]) -> None:
    for term in terms:
        if term in dst:
            dst.remove(term)
        else:
            dst.add(term)


@dataclass
class TermAlgebra:
    q_components: tuple[int, ...]

    def __post_init__(self) -> None:
        self.q_components = tuple(sorted(set(self.q_components), key=prime_power))
        self.q_degrees = tuple(prime_power(q)[0] for q in self.q_components)
        basis = [1]
        for degree in self.q_degrees:
            basis.append(basis[-1] * degree)
        self.basis = tuple(basis)
        self.term_count = self.basis[-1]
        self.index = {q: i for i, q in enumerate(self.q_components)}
        self.kappa_table = tuple(self._kappa_power(i) for i in range(len(self.q_components)))
        self._q_power_cache: dict[tuple[int, int, int], frozenset[int]] = {}
        self._term_product_cache: dict[tuple[int, int], frozenset[int]] = {}
        self._square_term_table: tuple[frozenset[int], ...] | None = None

    def _kappa_power(self, index: int) -> frozenset[int]:
        q = self.q_components[index]
        p = self.q_degrees[index]
        if p == 2:
            return frozenset((self.basis[index] - 1, self.basis[index]))
        if q == p:
            terms = [self.basis[self.index[q2]] for q2 in Q_SET[p]]
            excess = EXCESS[p]
            if excess:
                terms.append(excess.bit_length() - 1)
            return frozenset(terms)
        return frozenset((self.basis[index - 1],))

    def q_power_times_term(self, q_index: int, q_exponent: int, term: int) -> frozenset[int]:
        key = (q_index, q_exponent, term)
        if key in self._q_power_cache:
            return self._q_power_cache[key]

        p = self.q_degrees[q_index]
        current = (term % self.basis[q_index + 1]) // self.basis[q_index]
        new = current + q_exponent
        if new < p:
            result = frozenset((term + q_exponent * self.basis[q_index],))
        else:
            high = (term // self.basis[q_index + 1]) * self.basis[q_index + 1]
            high += (new % p) * self.basis[q_index]
            low = term % self.basis[q_index]
            terms: set[int] = set()
            for summand in self.kappa_table[q_index]:
                product = frozenset(high + t for t in self.term_times_term(low, summand))
                toggle_terms(terms, product)
            result = frozenset(terms)

        self._q_power_cache[key] = result
        return result

    def term_times_term(self, x: int, y: int) -> frozenset[int]:
        key = (x, y)
        if key in self._term_product_cache:
            return self._term_product_cache[key]

        terms = {y}
        for q_index in range(len(self.q_components) - 1, -1, -1):
            x_exponent = (x % self.basis[q_index + 1]) // self.basis[q_index]
            if x_exponent == 0:
                continue
            new_terms: set[int] = set()
            for term in terms:
                toggle_terms(new_terms, self.q_power_times_term(q_index, x_exponent, term))
            terms = new_terms

        result = frozenset(terms)
        self._term_product_cache[key] = result
        return result

    def multiply(self, a: frozenset[int], b: frozenset[int]) -> frozenset[int]:
        terms: set[int] = set()
        for x in a:
            for y in b:
                toggle_terms(terms, self.term_times_term(x, y))
        return frozenset(terms)

    def square(self, a: frozenset[int]) -> frozenset[int]:
        terms: set[int] = set()
        for x in a:
            toggle_terms(terms, self.term_times_term(x, x))
        return frozenset(terms)

    def square_with_table(self, a: frozenset[int]) -> frozenset[int]:
        if self._square_term_table is None:
            self._square_term_table = tuple(
                self.term_times_term(term, term) for term in range(self.term_count)
            )
        terms: set[int] = set()
        for x in a:
            toggle_terms(terms, self._square_term_table[x])
        return frozenset(terms)

    def power(self, a: frozenset[int], exponent: int) -> frozenset[int]:
        result = frozenset((0,))
        current = a
        while exponent:
            if exponent & 1:
                result = self.multiply(result, current)
            exponent >>= 1
            if exponent:
                current = self.square(current)
        return result

    def fixed_base_power(self, a: frozenset[int], exponent: int) -> frozenset[int]:
        result = frozenset((0,))
        if exponent == 0:
            return result

        term_times_base = tuple(
            self.multiply(frozenset((term,)), a) for term in range(self.term_count)
        )
        for bit in range(exponent.bit_length() - 1, -1, -1):
            result = self.square_with_table(result)
            if (exponent >> bit) & 1:
                terms: set[int] = set()
                for term in result:
                    toggle_terms(terms, term_times_base[term])
                result = frozenset(terms)
        return result

    def order(self, a: frozenset[int]) -> int:
        group_order = (1 << self.term_count) - 1
        factors = GROUP_ORDER_FACTORS[self.term_count]
        assert factor_product(factors) == group_order
        order = group_order
        for prime, exponent in factors.items():
            for _ in range(exponent):
                candidate = order // prime
                if order % prime == 0 and self.power(a, candidate) == frozenset((0,)):
                    order = candidate
                else:
                    break
        return order


def beta_for(target_p: int, excess: int) -> tuple[TermAlgebra, frozenset[int]]:
    q_set = Q_SET[target_p]
    components: set[int] = set()
    for q in q_set:
        components.update(primitive_components(q))
    finite = finite_summand(q_set, excess)
    components.update(finite_components_for(finite))
    algebra = TermAlgebra(tuple(components))

    terms = [bit for bit in range(finite.bit_length()) if (finite >> bit) & 1]
    terms.extend(algebra.basis[algebra.index[q]] for q in q_set if q % 2 == 1)
    return algebra, frozenset(terms)


def has_pth_root(algebra: TermAlgebra, beta: frozenset[int], p: int) -> bool:
    order = algebra.order(beta)
    return order % p != 0


def has_pth_root_by_power(algebra: TermAlgebra, beta: frozenset[int], p: int) -> bool:
    group_order = (1 << algebra.term_count) - 1
    if group_order % p != 0:
        return False
    return algebra.fixed_base_power(beta, group_order // p) == frozenset((0,))


def fixed_base_certificate(p: int, excess: int) -> str:
    algebra, beta = beta_for(p, excess)
    root = has_pth_root_by_power(algebra, beta, p)
    return (
        f"p={p}, m={excess}, Q={Q_SET[p]}, components={algebra.q_components}, "
        f"E={algebra.term_count}, root? {root}"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--deep",
        action="store_true",
        help="also run the E=19580 fixed-base certificate for p=179 (about a minute locally)",
    )
    args = parser.parse_args()

    cases = [
        (7, 0),
        (7, 1),
        (19, 1),
        (19, 4),
        (73, 1),
        (23, 1),
        (53, 1),
        (89, 1),
    ]
    header = (
        f"{'p':>3} {'Q(f(p))':>9} {'m':>2} {'components':>18} {'E':>4} "
        f"{'ord(beta)':>18} {'p|ord?':>7} {'root?':>6}"
    )
    print(header)
    print("-" * len(header))
    for p, m in cases:
        algebra, beta = beta_for(p, m)
        order = algebra.order(beta)
        print(
            f"{p:>3} {str(Q_SET[p]):>9} {m:>2} {str(algebra.q_components):>18} "
            f"{algebra.term_count:>4} {order:>18} {str(order % p == 0):>7} "
            f"{str(has_pth_root(algebra, beta, p)):>6}"
        )

    print()
    print("Reading: excess is the first row for each p where root? becomes False.")
    print("In particular, Q={9} has p=19 needing m=4, but p=73 already works at m=1.")
    print()
    print("Targeted fixed-base test beyond the Rust table:")
    print(fixed_base_certificate(47, 1))
    if args.deep:
        print(fixed_base_certificate(179, 1))
    else:
        print("p=179 deep certificate skipped; rerun with --deep")


if __name__ == "__main__":
    main()
