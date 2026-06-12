"""The 2*3^k exception column: the Lenstra excess is m_p = 4 exactly.

writeups/excess.tex section "The m=4 upper bound: the unexplored gap" asks whether
m_p = 4 exactly (not merely >= 4) for every prime p with f(p) = ord_p(2) = 2*3^k.
This probe carries the June 2026 answer: YES at every prime the factor tables can
see — universally for k <= 6 (fully factored levels), and at every known prime for
k = 7, 8.

THE CORRECTED NORM. excess.tex section 4.3 stated
Norm_{F_{2^(4h)}/F_{2^(2h)}}(kappa+4) = (kappa+4)(kappa+6); that is a slip. The
norm conjugate is sigma(4) where sigma = Frob^(2h) restricted to F_16; since
F_16 ∩ F_{2^(2h)} = F_4 (gcd(4, 2*3^k) = 2), sigma|F_16 is the NONTRIVIAL element
of Gal(F_16/F_4), which swaps the two roots {4, 5} of the Artin-Schreier minimal
polynomial y^2 + y + 2 of nimber 4 over F_4. (Equivalently: 2h = 2 mod 4 and the
Frobenius orbit of 4 is 4 -> 6 -> 5 -> 7, so Frob^(2h)(4) = Frob^2(4) = 4^4 = 5;
6 = 4^2 is Frob^1(4), one squaring short.) Hence, writing omega := zeta^h = kappa_2
(the cube root of unity "nimber 2", via the tower equation kappa_{3^k}^(3^k) = 2):

    N := Norm(kappa + 4) = (kappa + 4)(kappa + 5) = kappa^2 + kappa + omega,

by Artin-Schreier itself (4*5 = 4^2 + 4 = 2). The m = 4 translate test therefore
collapses into the SAME sparse trinomial field F_2[x]/(x^(2h) + x^h + 1) that the
C_k certification (cyclotomic_3k_family.py) already uses — no compositum
arithmetic is needed at all.

THE IN-FIELD CRITERION (the chain, each step PROVED/machine-checked here):

  kappa+4 has no p-th root              [excess def.; m_p >= 4 already PROVED]
    <=>  (kappa+4)^((2^(4h)-1)/p) != 1   [power criterion, v_p(2^(4h)-1) = 1]
    <=>  N^((2^(2h)-1)/p) != 1           [norm reduction + the corrected N]
    <=>  M^((2^h+1)/p) != 1              [F* = L* x U; M := Nbar/N = N^(2^h-1)
                                          is the unit-circle part squared, and
                                          p | 2^h+1 since f(p) = 2h]

so: m_p = 4  <=>  p | ord(M_k),  M_k = Nbar/N,  N = zeta^2 + zeta + zeta^h,
all inside F_2[x]/(x^(2h) + x^h + 1). Wieferich guard: the middle equivalences
need v_p(2^(2h)-1) = 1, asserted per prime (cf. excess.tex Remark 3).

RESULTS (this script; claim levels per the excess.tex convention):

  CERTIFIED (k <= 6, fully factored levels — universal on the column):
    m_p = 4 for EVERY prime p with f(p) = 2*3^k, k = 2..6:
      k=2 (f=18):    19                                  [DiMuro anchor row]
      k=3 (f=54):    87211                               [NEW: only >=4 known]
      k=4 (f=162):   163 [calculator anchor], 135433, 272010961        [NEW]
      k=5 (f=486):   1459 [calculator anchor], 139483,
                     10429407431911334611, 918125051602568899753       [NEW]
      k=6 (f=1458):  227862073, 3110690934667, 216892513252489863991753,
                     1102099161075964924744009, P78                    [NEW]
  CONSISTENT (k = 7, 8 — every KNOWN prime of the level; cofactors unfactored):
      k=7 (f=4374):  17497, 5419387, 6049243, 24796936459, 184318815979,
                     1104193455232918700687932947755896164888067,
                     6326666886932800988419273258756815291776881
      k=8 (f=13122): 52489, 28934011, 526084074721, 28905359674006243,
                     72687062849889979
  No counterexample: the candidate 0/1/4 rule's upper bound m_p <= 4 survives
  every prime any current factor table can reach on this column. An m_p >= 5
  example, if one exists, now hides strictly inside the unfactored cofactors of
  Phi_{2*3^7}(2) and beyond.

  PROVED (the twisted norm-tower lemma — why no C_k-style descent): with
  eta = zeta^3 = zeta_{k-1},
      Norm_{F_k/F_{k-1}}(N_k) = eta^2 + omega^2*eta + 1  =: G_{k-1},
  a TWISTED companion of N_{k-1} = eta^2 + eta + omega, not N_{k-1} itself nor a
  conjugate/zeta-shift of it (checked). So old-prime parts do not propagate up
  the way gamma's did, and the per-level conjecture D_k below is genuinely
  per-level. (Derivation: sigma(zeta) = omega^a * zeta for the level-3 Galois
  group, expand the three-factor product; machine-checked at k = 2, 3, 4.)

  CONJECTURE D_k (the column analogue of C_k): the prime-to-3 part of ord(M_k)
  is full, i.e. ((2^h+1)/3^(k+1)) | ord(M_k). Equivalent (given the audits) to
  m_p = 4 for every prime with f(p) = 2*3^k. Status: certified k <= 6,
  consistent k = 7, 8. The 3-part of ord(M_k) is NOT always full and carries no
  column information; measured values are printed per level.

CROSS-VALIDATION. (1) The independent term algebra of ordinal_excess_probe.py
(which knows nothing of the zeta-reading or the norm reduction) re-derives the
root/no-root pattern m = 0..4 at p = 19 and m = 4 at p = 87211. (2) An explicit
Artin-Schreier compositum model F_{2^(4h)} = F[y]/(y^2 + y + omega) verifies
sigma(4) = 5, the norm identity, and (for k <= 6) the full direct power test
(kappa+4)^((2^(4h)-1)/p) != 1 on the smallest prime of each level, against the
in-field M route. (3) The DiMuro/calculator anchor rows m_19 = m_163 = m_1459 = 4
are reproduced, never assumed.

FACTORIZATION PROVENANCE. Phi_{2*3^k}(2) = 2^(2*3^(k-1)) - 2^(3^(k-1)) + 1.
k <= 4 factorizations are classical/small; k = 5, 6 are factordb "FF" entries
(fetched 2026-06-12), re-verified here by exact product reconstruction and
primality tests (deterministic Miller-Rabin below 3.3e24, MR-64 PRP above; the
P78 of Phi_{2*3^6}(2) is PRP-local, factordb marks it proven). k = 7, 8 known
factors (factordb CF entries, same fetch; the second 43-digit k=7 prime added
2026-06-12 after an adversarial review caught the first fetch dropping it) are
verified individually by exact divisibility and order tests, and the small ones
re-derived by direct sieve over p = 2*3^k*t + 1; the unfactored cofactors are
composite (no universal claim is made at those levels). NOTE: these hardcoded
known-factor lists go stale as factordb advances; staleness can only shrink the
k = 7, 8 coverage claim, never break a certified row. Squarefreeness (the Wieferich guard) is asserted per
prime. 3-part: v_3(Phi_{2*3^k}(2)) = 1 exactly (LTE), asserted.

These are excess-table science (A380496-type rows), not new shippable alpha_u
carries: the Rust tower's operational boundary at alpha_53 is untouched.
"""

from __future__ import annotations

from cyclotomic_3k_family import certify_irreducible, fmul, fpow, is_prime, tri_reduce

# ---------------------------------------------------------------------------
# Fast arithmetic helpers on top of the trinomial model: squaring via byte
# spreading (squaring in F_2[x] interleaves zeros), inversion via poly xgcd.
# ---------------------------------------------------------------------------

_SPREAD = tuple(sum(((b >> i) & 1) << (2 * i) for i in range(8)) for b in range(256))


def fsq(a: int, h: int) -> int:
    src = a.to_bytes((a.bit_length() + 7) // 8 or 1, "little")
    dst = bytearray(2 * len(src))
    for i, b in enumerate(src):
        w = _SPREAD[b]
        dst[2 * i] = w & 0xFF
        dst[2 * i + 1] = w >> 8
    return tri_reduce(int.from_bytes(bytes(dst), "little"), h)


def frob(a: int, t: int, h: int) -> int:
    """a^(2^t) by t fast squarings."""
    for _ in range(t):
        a = fsq(a, h)
    return a


def fpow_f(a: int, e: int, h: int) -> int:
    """fpow with fast squaring (matches fpow; pinned by an audit assert)."""
    r = 1
    while e:
        if e & 1:
            r = fmul(r, a, h)
        e >>= 1
        if e:
            a = fsq(a, h)
    return r


def finv(a: int, h: int) -> int:
    """Inverse mod x^(2h)+x^h+1 by extended Euclid in F_2[x]."""
    f = (1 << (2 * h)) | (1 << h) | 1
    r0, r1, s0, s1 = f, a, 0, 1
    while r1:
        d = r0.bit_length() - r1.bit_length()
        if d < 0:
            r0, r1, s0, s1 = r1, r0, s1, s0
            d = -d
        r0 ^= r1 << d
        s0 ^= s1 << d
    assert r0 == 1, "gcd != 1: modulus not irreducible or a = 0"
    return tri_reduce(s0, h)


# ---------------------------------------------------------------------------
# The Artin-Schreier compositum F_{2^(4h)} = F[y]/(y^2 + y + omega), elements
# (a, b) = a + b*y. Independent cross-check route only; never the main path.
# ---------------------------------------------------------------------------


def comp_mul(p: tuple[int, int], q: tuple[int, int], h: int, c: int) -> tuple[int, int]:
    a, b = p
    d, e = q
    be = fmul(b, e, h)
    return (fmul(a, d, h) ^ fmul(be, c, h), fmul(a, e, h) ^ fmul(b, d, h) ^ be)


def comp_sq(p: tuple[int, int], h: int, c: int) -> tuple[int, int]:
    a, b = p
    b2 = fsq(b, h)
    return (fsq(a, h) ^ fmul(b2, c, h), b2)


def comp_pow(p: tuple[int, int], e: int, h: int, c: int) -> tuple[int, int]:
    r = (1, 0)
    while e:
        if e & 1:
            r = comp_mul(r, p, h, c)
        e >>= 1
        if e:
            p = comp_sq(p, h, c)
    return r


# ---------------------------------------------------------------------------
# Column data: the primes with f(p) = ord_p(2) = 2*3^k, i.e. the prime factors
# of Phi_{2*3^k}(2)/3. Provenance in the module docstring; audited in main().
# ---------------------------------------------------------------------------

P78 = 393063301203384521164229656203691748263012766081190297429488962985651210769817

COLUMN_PRIMES: dict[int, tuple[int, ...]] = {
    2: (19,),
    3: (87211,),
    4: (163, 135433, 272010961),
    5: (1459, 139483, 10429407431911334611, 918125051602568899753),
    6: (
        227862073,
        3110690934667,
        216892513252489863991753,
        1102099161075964924744009,
        P78,
    ),
    # Known prime factors only; cofactors composite and unfactored (factordb CF).
    7: (
        17497,
        5419387,
        6049243,
        24796936459,
        184318815979,
        1104193455232918700687932947755896164888067,
        6326666886932800988419273258756815291776881,
    ),
    8: (52489, 28934011, 526084074721, 28905359674006243, 72687062849889979),
}

FULLY_FACTORED_LEVELS = (2, 3, 4, 5, 6)

# Rows already in the DiMuro table / calculator records: reproduced, not assumed.
ANCHOR_ROWS = {19, 163, 1459}


def phi_2x3k(k: int) -> int:
    hp = 3 ** (k - 1)
    return (1 << (2 * hp)) - (1 << hp) + 1


def sieve_column_factors(k: int, t_max: int) -> list[int]:
    """Primes p = 2*3^k*t + 1 with ord_p(2) = 2*3^k exactly."""
    base = 2 * 3**k
    found = []
    for t in range(1, t_max + 1):
        p = base * t + 1
        if (
            pow(2, base, p) == 1
            and pow(2, base // 2, p) != 1
            and pow(2, base // 3, p) != 1
            and is_prime(p)[0]
        ):
            found.append(p)
    return found


# ---------------------------------------------------------------------------
# Audits: factorizations, primality, order conditions, Wieferich, LTE.
# ---------------------------------------------------------------------------


def audit() -> None:
    prp_only = []
    for k, primes in COLUMN_PRIMES.items():
        phi = phi_2x3k(k)
        assert phi % 3 == 0 and phi % 9 != 0, f"k={k}: v_3(Phi) != 1"
        prod = 3
        for p in primes:
            ok, kind = is_prime(p)
            assert ok, f"k={k}: {p} composite"
            if kind == "prp":
                prp_only.append(p)
            assert phi % p == 0, f"k={k}: {p} does not divide Phi_{{2*3^{k}}}(2)"
            base = 2 * 3**k
            assert pow(2, base, p) == 1, f"k={k}: ord_{p}(2) does not divide 2*3^{k}"
            assert pow(2, base // 2, p) != 1 and pow(2, base // 3, p) != 1, (
                f"k={k}: ord_{p}(2) proper divisor of 2*3^{k}"
            )
            # Wieferich guard: v_p(2^(2h)-1) = 1, i.e. p is not a base-2
            # Wieferich-type prime at its own level.
            assert pow(2, base, p * p) != 1, f"k={k}: Wieferich prime {p}"
            prod *= p
        if k in FULLY_FACTORED_LEVELS:
            assert prod == phi, f"k={k}: factor product != Phi_{{2*3^{k}}}(2)"
        else:
            cof = phi // prod
            assert phi % prod == 0 and cof > 1 and not is_prime(cof)[0], (
                f"k={k}: cofactor bookkeeping broken"
            )
    # independent re-derivation of the small k = 7, 8 factors by sieve
    assert sieve_column_factors(7, 1400) == [17497, 5419387, 6049243]
    assert sieve_column_factors(8, 5) == [52489]
    # fast-squaring path pinned to the reference multiply
    assert fsq(0b1011001, 9) == fmul(0b1011001, 0b1011001, 9)
    v = (1 << 17) ^ (1 << 9) ^ 5
    assert fpow_f(v, 12345, 9) == fpow(v, 12345, 9)
    assert fmul(v, finv(v, 9), 9) == 1
    print("audits OK: products exact (k<=6), every factor prime + order-verified,")
    print("  squarefree at its level, v_3(Phi)=1; cofactors composite (k=7,8);")
    print("  sieve re-derives 17497, 5419387, 6049243 (k=7) and 52489 (k=8);")
    print(f"  PRP-only (no local deterministic proof): {[len(str(p)) for p in prp_only]}-digit primes\n")


# ---------------------------------------------------------------------------
# Per-level verification.
# ---------------------------------------------------------------------------


def m_element(h: int) -> int:
    """M = Nbar/N for N = zeta^2 + zeta + zeta^h in F_2[x]/(x^(2h)+x^h+1)."""
    n = 4 ^ 2 ^ (1 << h)  # zeta^2 + zeta + omega, zeta = x
    return fmul(frob(n, h, h), finv(n, h), h)


def compositum_checks(k: int, direct_prime: int | None) -> None:
    """The AS-compositum route: sigma(4) = 5, Norm = N, optional direct test."""
    h = 3**k
    c = 1 << h  # omega
    y = (0, 1)  # nimber 4: y^2 + y + omega = 0
    assert comp_mul(y, y, h, c) == (c, 1), "y^2 != y + omega"
    s = y
    for _ in range(2 * h):
        s = comp_sq(s, h, c)
    assert s == (1, 1), f"sigma(4) != 5 at k={k}"
    kap4 = (2, 1)  # kappa + 4
    n = 4 ^ 2 ^ c
    assert comp_mul(kap4, (2 ^ 1, 1), h, c) == (n, 0), f"Norm != kappa^2+kappa+omega at k={k}"
    if direct_prime is not None:
        e = ((1 << (4 * h)) - 1) // direct_prime
        assert comp_pow(kap4, e, h, c) != (1, 0), (
            f"direct compositum test disagrees with the M route at k={k}"
        )


def twisted_norm_lemma(k: int) -> None:
    """Norm_{F_k/F_{k-1}}(N_k) = eta^2 + omega^2*eta + 1, eta = zeta^3."""
    h = 3**k
    c = 1 << h
    n = 4 ^ 2 ^ c
    t = 2 * 3 ** (k - 1)
    nrm = fmul(fmul(n, frob(n, t, h), h), frob(n, 2 * t, h), h)
    eta = 8  # zeta^3
    g = fsq(eta, h) ^ fmul(fsq(c, h), eta, h) ^ 1
    assert nrm == g, f"twisted norm lemma fails at k={k}"


def verify_level(k: int) -> list[int]:
    h = 3**k
    u = (1 << h) + 1
    certify_irreducible(h)
    m = m_element(h)
    assert fpow_f(m, u, h) == 1, "M not in the unit circle"
    assert frob(m, h, h) == finv(m, h), "Mbar != M^-1: not a circle element"

    label = f"k={k} (f=2*3^{k}={2 * 3 ** k}, F_2[x]/(x^{2 * h}+x^{h}+1))"
    failures = []
    certified = []
    for p in COLUMN_PRIMES[k]:
        if fpow_f(m, u // p, h) == 1:
            failures.append(p)
        else:
            certified.append(p)
    if failures:
        print(f"{label}: *** m_p >= 5 at {failures}: THE 0/1/4 RULE IS BROKEN ***")
        return failures

    anchors = sorted(set(certified) & ANCHOR_ROWS)
    tag = f" (anchors reproduced: {anchors})" if anchors else ""
    if k in FULLY_FACTORED_LEVELS:
        # exact ord(M): strip prime factors of |U| = 3^(k+1) * (column primes <= k)
        factors = [3] * (k + 1)
        for j in range(2, k + 1):
            factors.extend(COLUMN_PRIMES[j])
        prod = 1
        for q in factors:
            prod *= q
        assert prod == u, f"k={k}: 2^h+1 factorization chain broken"
        o = u
        for q in sorted(set(factors)):
            while o % q == 0 and fpow_f(m, o // q, h) == 1:
                o //= q
        v3 = 0
        while o % 3 == 0:
            o //= 3
            v3 += 1
        # the prime-to-3 part of ord(M) must be full for D_k; recompute cleanly
        odd_part_full = all(fpow_f(m, u // p, h) != 1 for j in range(2, k + 1) for p in COLUMN_PRIMES[j])
        assert odd_part_full, f"k={k}: an old-level prime is missing from ord(M)"
        print(f"{label}: m_p = 4 CERTIFIED for ALL primes of the level{tag}")
        print(f"   ord(M_{k}) = (2^{h}+1)/3^{k + 1 - v3}; 3-part 3^{v3} of max 3^{k + 1}; D_{k} holds")
    else:
        print(f"{label}: m_p = 4 at every KNOWN prime: {certified}")
        print(f"   (cofactor of Phi_{{2*3^{k}}}(2) unfactored: level consistent, not complete)")
    return failures


def term_algebra_cross_check() -> None:
    """Independent term-algebra re-derivation (knows nothing of the norm route)."""
    from ordinal_excess_probe import Q_SET, beta_for, has_pth_root_by_power

    # 87211 is beyond the probe's tabulated rows; f(87211) = 54 = 2*27, so its
    # odd component set is {27} (same convention as Q_SET[19] = (9,)).
    Q_SET.setdefault(87211, (27,))
    for m in range(4):
        alg, beta = beta_for(19, m)
        assert has_pth_root_by_power(alg, beta, 19), f"term algebra: kappa_9+{m} lost its 19th root"
    for p in (19, 87211):
        alg, beta = beta_for(p, 4)
        assert not has_pth_root_by_power(alg, beta, p), f"term algebra: kappa+4 has a {p}th root"
    print("term-algebra cross-check: root pattern m=0..3/4 at p=19; m=4 at p=87211 OK")


def main() -> None:
    audit()
    term_algebra_cross_check()
    for k in (2, 3, 4):
        twisted_norm_lemma(k)
    print("twisted norm lemma Norm(N_k) = eta^2 + omega^2*eta + 1 verified (k=2,3,4)\n")

    any_failures = False
    for k in (2, 3, 4, 5, 6, 7, 8):
        direct = min(COLUMN_PRIMES[k]) if k <= 6 else None
        compositum_checks(k, direct)
        any_failures |= bool(verify_level(k))
    print()
    if any_failures:
        print("FAILURE FOUND: an m_p >= 5 row exists — the 0/1/4 rule is refuted; see above.")
    else:
        print("No failure: m_p = 4 exactly on the whole 2*3^k column for k <= 6 (universal),")
        print("and at every known prime for k = 7, 8. The 0/1/4 upper bound m <= 4 survives")
        print("every prime any current factor table can reach on this column.")


if __name__ == "__main__":
    main()
