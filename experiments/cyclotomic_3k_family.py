"""The 3-power excess family: ord(kappa_{3^k} + 1) = 3^(k+1) * (2^(3^k) - 1).

writeups/excess.tex (the 3^k family thread) asks to prove this formula or find its
first failure. This probe carries the June 2026 result:

THE KEY RECOGNITION. The tower relations (Lenstra/DiMuro; independently encoded
in `ordinal_excess_probe.py`) give `kappa_2^2 = kappa_2 + 1`, `kappa_3^3 =
kappa_2`, and `kappa_{3^j}^3 = kappa_{3^(j-1)}`. So `zeta := kappa_{3^k}`
satisfies `zeta^(3^k) = kappa_2`, an element of order 3: **zeta is a primitive
3^(k+1)-th root of unity**. Since 2 is a primitive root mod 3^(k+1), the
cyclotomic polynomial `Phi_{3^(k+1)}(x) = x^(2h) + x^h + 1` (h = 3^k) is
irreducible over F_2, and the component field is `F = F_2(zeta) = F_{2^(2h)}`
with index-2 subfield `L = F_{2^h}`. All arithmetic below happens in the sparse
trinomial model `F_2[x]/(x^(2h) + x^h + 1)`.

PROVED (machine-checked here at every level, plus in the term algebra of
`ordinal_excess_probe.py` for k = 1, 2):

* Half-angle splitting. With `s = (3^(k+1)+1)/2` (the inverse of 2 mod 3^(k+1)),
      kappa_{3^k} + 1  =  zeta^s * (zeta^s + zeta^-s),
  where `zeta^s` lies in the norm-one circle `U` (order exactly 3^(k+1)) and
  `(zeta^s + zeta^-s)^2 = zeta + zeta^-1 =: gamma_k` lies in `L*`. Since
  `F* = L* x U` with coprime orders,
      ord(kappa_{3^k} + 1)  =  3^(k+1) * ord(gamma_k),   ord(gamma_k) | 2^(3^k)-1.
  Corollaries used as machine checks: `(kappa+1)^(2^h-1) = zeta^-1` and
  `Norm_{F/L}(kappa+1) = gamma_k` (the closed-form instance of
  writeups/excess.tex's norm-reduction identity when E/f = 2).
* Translates. `kappa_{3^k} + 2 = zeta + zeta^(3^k)` and `kappa_{3^k} + 3 =
  zeta + zeta^(2*3^k)` split the same way with `gamma_k` replaced by a Galois
  conjugate, so the m = 1, 2, 3 translates ALL have order `3^(k+1) * ord(gamma_k)`,
  and `ord(kappa_{3^k}) = 3^(k+1)`.
* The 2*3^k exception, unconditionally. Any prime p with `f(p) = 2*3^k` divides
  `2^(3^k)+1`, hence divides neither `3^(k+1)` nor `2^(3^k)-1`, hence p never
  divides `ord(kappa_{3^k}+m)` for m in {0,1,2,3}: by the order criterion the
  Lenstra excess has `m_p >= 4`. The *full* formula is NOT needed for the
  exception - the splitting alone forces it. (New reach example: p = 87211,
  f = 54, beyond the current tables, has m_p >= 4.)
* Norm tower. `gamma_{k-1} = gamma_k^3 + gamma_k = Norm_{L_k/L_{k-1}}(gamma_k)`
  (minimal polynomial `X^3 + X + gamma_{k-1}`), so `ord(gamma_{k-1}) | ord(gamma_k)`;
  with LTE (`v_r(2^(3^k)-1) = v_r(2^(f(r))-1)` for old primes r != 3) full order
  parts PROPAGATE upward. The conjecture C_k ("gamma_k primitive in L") reduces
  to C_{k-1} plus full order parts at the NEW primes r | Phi_{3^k}(2).
* Equivalence. For a prime r with `f(r) = 3^k` (exactly the primitive prime
  factors of `Phi_{3^k}(2)`), `Q(f(r)) = {3^k}` and `m_r = 1  <=>  r | ord(gamma_k)`.
  So, when `2^(3^k)-1` is squarefree across its known factors (it is, in range):
      C_k  <=>  m_r = 1 for every prime r with f(r) in {3, 9, ..., 3^k}.
  The family formula IS the candidate 0/1/4 rule restricted to the 3-power
  column. Also structurally forced: on that column `m_r in {0,2,3}` is
  impossible (m = 0 has a root since ord(kappa) = 3^(k+1); m = 2,3 share m = 1's
  order), so any failure of C_k jumps straight to `m_r >= 4`.

VERIFIED STATUS (this script):
  k = 1..6  C_k holds: gamma_k is primitive, ord(kappa_{3^k}+1) = 3^(k+1)*(2^(3^k)-1).
  k = 7, 8  consistent: all KNOWN prime factors of Phi_{3^7}(2), Phi_{3^8}(2)
            divide ord(gamma_k); full certification blocked only by the
            unfactored cofactors (factordb status CF).
Cross-validations: term algebra reproduces k = 1, 2; the recorded calculator
rows m_2593 = 1 (k = 4) and m_487 = 1 (k = 5) are re-derived independently.
Newly certified excess rows (order criterion, independent of the cofactors):
  f = 27:   m = 1 for 262657
  f = 81:   m = 1 for 71119, 97685839
  f = 243:  m = 1 for 16753783618801, 192971705688577, 3712990163251158343
  f = 729:  m = 1 for 80191, 97687, 379081, P42, P90
  f = 2187: m = 1 for 39367, 7606246033, 263196614521, 529063556041
  f = 6561: m = 1 for 209953, 1299079, 70063267397606709277393
These are excess-table science (A380496-type rows), not new shippable alpha_u
carries: the Rust tower's operational boundary at alpha_53 is untouched.

FACTORIZATION PROVENANCE. 2^(3^k)-1 for k <= 4 is classical; the Phi_243(2) and
Phi_729(2) splittings are factordb "FF" entries (2026-06), re-verified here by
exact product reconstruction and primality tests (deterministic Miller-Rabin
below 3.3e24, MR-64 PRP above - only the 42- and 90-digit primes of Phi_729(2)
are PRP-local; factordb marks them proven). The small factors 39367, 209953,
1299079 are also re-derived by direct sieve over r = 2*3^j*t + 1.
"""

from __future__ import annotations

import random

# ---------------------------------------------------------------------------
# GF(2)[x] arithmetic; ints as coefficient bit-vectors, modulus x^(2h)+x^h+1.
# ---------------------------------------------------------------------------


def gf2_mul(a: int, b: int) -> int:
    r = 0
    while b:
        lsb = b & -b
        r ^= a << (lsb.bit_length() - 1)
        b ^= lsb
    return r


def tri_reduce(v: int, h: int) -> int:
    two_h = h << 1
    mask = (1 << two_h) - 1
    while True:
        hi = v >> two_h
        if not hi:
            return v & mask
        v = (v & mask) ^ hi ^ (hi << h)


def fmul(a: int, b: int, h: int) -> int:
    return tri_reduce(gf2_mul(a, b), h)


def fpow(a: int, e: int, h: int) -> int:
    r = 1
    while e:
        if e & 1:
            r = fmul(r, a, h)
        e >>= 1
        if e:
            a = fmul(a, a, h)
    return r


def poly_mod(a: int, f: int) -> int:
    df = f.bit_length() - 1
    while a and a.bit_length() - 1 >= df:
        a ^= f << (a.bit_length() - 1 - df)
    return a


def poly_gcd(a: int, b: int) -> int:
    while b:
        a, b = b, poly_mod(a, b)
    return a


def certify_irreducible(h: int) -> None:
    """Certify x^(2h)+x^h+1 irreducible over F_2 (n = 2h has prime divisors 2, 3)."""
    n = 2 * h
    f = (1 << n) | (1 << h) | 1
    s = 2
    frob = {}
    for d in range(1, n + 1):
        s = fmul(s, s, h)
        frob[d] = s
    assert frob[n] == 2, "x^(2^n) != x"
    for d in (n // 2, n // 3):
        assert poly_gcd(frob[d] ^ 2, f) == 1, f"gcd(x^(2^{d})-x, f) != 1: reducible"


# ---------------------------------------------------------------------------
# Primality: deterministic Miller-Rabin below 3.3e24, MR-64 PRP above.
# ---------------------------------------------------------------------------

SMALL_DET_BASES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
DET_LIMIT = 3317044064679887385961981


def is_prime(n: int, rounds: int = 64) -> tuple[bool, str]:
    if n < 2:
        return False, "det"
    for p in SMALL_DET_BASES:
        if n % p == 0:
            return (n == p), "det"
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1

    def witness(a: int) -> bool:
        x = pow(a, d, n)
        if x in (1, n - 1):
            return False
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                return False
        return True

    if n < DET_LIMIT:
        return (not any(witness(a) for a in SMALL_DET_BASES)), "det"
    rng = random.Random(0xC0FFEE)
    return (not any(witness(rng.randrange(2, n - 1)) for _ in range(rounds))), "prp"


# ---------------------------------------------------------------------------
# Verified factorizations of 2^(3^k)-1 (products re-checked in main()).
# ---------------------------------------------------------------------------

PHI_3J_2_FACTORS: dict[int, tuple[int, ...]] = {
    1: (7,),
    2: (73,),
    3: (262657,),
    4: (2593, 71119, 97685839),
    5: (487, 16753783618801, 192971705688577, 3712990163251158343),
    6: (
        80191,
        97687,
        379081,
        664728004346558283448724389870269691211809,
        101213745778143742250901040788003424950068418098259161142719688891708905138274462262307761,
    ),
    # Known prime factors only; cofactors composite and unfactored (factordb CF).
    7: (39367, 7606246033, 263196614521, 529063556041),
    8: (209953, 1299079, 70063267397606709277393),
}

FULLY_FACTORED_LEVELS = (1, 2, 3, 4, 5, 6)


def factors_through(k: int) -> list[int]:
    out: list[int] = []
    for j in range(1, k + 1):
        out.extend(PHI_3J_2_FACTORS[j])
    return sorted(out)


def sieve_small_factors(j: int, t_max: int) -> list[int]:
    """Primes r = 2*3^j*t + 1 with ord_r(2) = 3^j, i.e. r | Phi_{3^j}(2)."""
    base = 2 * 3**j
    found = []
    for t in range(1, t_max + 1):
        r = base * t + 1
        if pow(2, 3**j, r) == 1 and pow(2, 3 ** (j - 1), r) != 1 and is_prime(r)[0]:
            found.append(r)
    return found


# ---------------------------------------------------------------------------
# Per-level verification.
# ---------------------------------------------------------------------------


def verify_level(k: int, complete: bool) -> list[int]:
    """Check identities at level k; return the new primes r with gamma an r-th power."""
    h = 3**k
    n = (1 << h) - 1
    certify_irreducible(h)
    x = 2  # zeta
    x_inv = (1 << (2 * h - 1)) ^ (1 << (h - 1))
    assert fmul(x, x_inv, h) == 1
    assert fpow(x, 3 ** (k + 1), h) == 1 and fpow(x, 3**k, h) != 1
    beta = x ^ 1
    gamma = x ^ x_inv

    # circle / norm / half-angle identities
    assert fpow(beta, n, h) == x_inv, "circle identity fails"
    assert fmul(beta, fpow(beta, 1 << h, h), h) == gamma, "norm identity fails"
    half = pow(2, -1, 3 ** (k + 1))
    z_half = fpow(x, half, h)
    assert fmul(z_half, z_half ^ fpow(x_inv, half, h), h) == beta, "half-angle fails"
    assert fpow(gamma, 1 << h, h) == gamma, "gamma not in the subfield"
    # translates m = 2, 3: circle parts are primitive 3^(k+1)-th roots of unity
    assert fpow(x ^ (1 << h), n, h) == fpow(x, 2 * h - 1, h), "m=2 circle fails"
    assert fpow(x ^ (1 << h) ^ 1, n, h) == 1 << (h - 1), "m=3 circle fails"

    assert fpow(gamma, n, h) == 1
    new = PHI_3J_2_FACTORS[k]
    old = factors_through(k - 1) if k > 1 else []
    failures = [r for r in new if fpow(gamma, n // r, h) == 1]
    # old primes are guaranteed by the norm tower; spot-check anyway when complete
    if complete:
        assert not [r for r in old if fpow(gamma, n // r, h) == 1], "norm tower broken"

    label = f"k={k} (h=3^{k}={h}, F_2(zeta_{3 ** (k + 1)}) = F_2[x]/(x^{2 * h}+x^{h}+1))"
    if failures:
        print(f"{label}: *** C_{k} FAILS at new primes {failures} ***")
    elif complete:
        print(f"{label}: identities OK; gamma primitive -> ord(kappa+1) = 3^{k + 1}*(2^{h}-1)")
    else:
        print(f"{label}: identities OK; all {len(new)} known new primes divide ord(gamma)")
        print(f"   -> m_r = 1 certified for {new}; C_{k} open pending the unfactored cofactor")
    return failures


def term_algebra_cross_check() -> None:
    """Anchor the cyclotomic model to the independent term algebra (k = 1, 2).

    Verifies, against `ordinal_excess_probe.TermAlgebra` (which knows nothing of
    the zeta-reading), that kappa_{3^k} has order 3^(k+1), that the circle and
    norm identities hold, that the m = 1, 2, 3 translates share one order, that
    the order equals 3^(k+1)*(2^(3^k)-1), and that the norm-tower image
    gamma_2^3 + gamma_2 has order 7 (a conjugate of gamma_1).
    """
    from ordinal_excess_probe import TermAlgebra

    one = frozenset((0,))
    for k, comps, n in ((1, (2, 3), 7), (2, (2, 3, 9), 511)):
        alg = TermAlgebra(comps)
        top = comps[-1]
        kap = frozenset((alg.basis[alg.index[top]],))
        b1 = frozenset(kap | {0})
        b2 = frozenset(kap | {1})
        b3 = frozenset(kap | {0, 1})
        assert alg.power(kap, 3 ** (k + 1)) == one and alg.power(kap, 3**k) != one
        kap_inv = alg.power(kap, 3 ** (k + 1) - 1)
        gamma = frozenset(kap ^ kap_inv)
        assert alg.power(b1, n) == kap_inv, "term-algebra circle identity fails"
        assert alg.power(b1, n + 2) == gamma, "term-algebra norm identity fails"
        want = 3 ** (k + 1) * n
        for b in (b1, b2, b3):
            assert alg.order(b) == want, "term-algebra order mismatch"
        assert alg.order(gamma) == n, "term-algebra gamma not primitive"
        if k == 2:
            g3g = frozenset(alg.power(gamma, 3) ^ gamma)
            assert alg.power(g3g, 7) == one and g3g != one
    print("term-algebra cross-check (k=1,2): orders, circle/norm, translates, norm tower OK")


def main() -> None:
    term_algebra_cross_check()
    # factorization audit: exact products, squarefreeness, primality
    for k in FULLY_FACTORED_LEVELS:
        prod = 1
        for r in factors_through(k):
            prod *= r
        assert prod == (1 << 3**k) - 1, f"k={k}: factor product != 2^(3^{k})-1"
    for k in (7, 8):
        for r in PHI_3J_2_FACTORS[k]:
            assert ((1 << 3**k) - 1) % r == 0
    prp_only = []
    for j, rs in PHI_3J_2_FACTORS.items():
        for r in rs:
            ok, kind = is_prime(r)
            assert ok, f"{r} composite"
            if kind == "prp":
                prp_only.append(r)
            assert ((1 << 3**j) - 1) % (r * r) != 0, f"level-Wieferich {r}"
    # LTE: v_3(2^(3^k)+1) = k+1 exactly
    for k in range(1, 9):
        assert (2 ** 3**k + 1) % 3 ** (k + 1) == 0 and (2 ** 3**k + 1) % 3 ** (k + 2) != 0
    print("factorizations audited: products exact, all factors prime, squarefree;")
    print(f"  PRP-only (no local deterministic proof): {[len(str(r)) for r in prp_only]}-digit primes")
    # independent re-derivation of the small Phi_{3^7}, Phi_{3^8} factors
    assert sieve_small_factors(7, 100000) == [39367]
    assert sieve_small_factors(8, 35000) == [209953, 1299079]
    print("  sieve re-derives 39367 (j=7) and 209953, 1299079 (j=8)\n")

    any_failures = False
    for k in (1, 2, 3, 4, 5, 6):
        any_failures |= bool(verify_level(k, complete=True))
    for k in (7, 8):
        any_failures |= bool(verify_level(k, complete=False))
    print()
    if any_failures:
        print("FAILURE FOUND: see above - this falsifies the 0/1/4 candidate rule too.")
    else:
        print("No failure: formula proved-equal to 3^(k+1)*ord(gamma_k) and verified")
        print("primitive for k <= 6; k = 7, 8 consistent on every known prime factor.")


if __name__ == "__main__":
    main()
