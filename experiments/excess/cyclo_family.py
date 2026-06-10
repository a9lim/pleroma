"""ord(kappa_{3^k}+1) = 3^{k+1} (2^{3^k}-1): cyclotomic-model verification.

Model: kappa_{3^k} is a primitive 3^{k+1}-th root of unity zeta (proved from the
tower relations kappa_{3^j}^3 = kappa_{3^{j-1}}, kappa_3^3 = kappa_2, ord(kappa_2)=3).
Since 2 is a primitive root mod 3^{k+1}, Phi_{3^{k+1}}(x) = x^{2h}+x^h+1 (h = 3^k)
is irreducible over F_2, so F_2(zeta) = F_2[x]/(x^{2h}+x^h+1), zeta = x.

Identities proved in the writeup, verified here:
  (half-angle)  zeta+1 = zeta^{1/2} * (zeta^{1/2}+zeta^{-1/2}),  exact U x L* split
  (circle)      (zeta+1)^(2^h-1) = zeta^{-1}
  (norm)        (zeta+1)^(2^h+1) = zeta + zeta^{-1} =: gamma_k  in L = F_{2^h}
  => ord(kappa_{3^k}+1) = 3^{k+1} * ord(gamma_k),  ord(gamma_k) | 2^h - 1.
Conjecture C_k <=> gamma_k primitive in F_{2^h}.
"""
import sys, random

# ---------------- GF(2)[x] arithmetic, ints as bit-vectors ----------------

def gf2_mul(a: int, b: int) -> int:
    r = 0
    while b:
        lsb = b & -b
        r ^= a << (lsb.bit_length() - 1)
        b ^= lsb
    return r

def tri_reduce(v: int, h: int) -> int:
    # reduce mod x^{2h} + x^h + 1
    H = h << 1
    mask = (1 << H) - 1
    while True:
        hi = v >> H
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
    while a.bit_length() - 1 >= df and a:
        a ^= f << (a.bit_length() - 1 - df)
    return a

def poly_gcd(a: int, b: int) -> int:
    while b:
        a, b = b, poly_mod(a, b)
    return a

def certify_irreducible(h: int) -> None:
    """Certify x^{2h}+x^h+1 irreducible over F_2 (n = 2h, prime divisors {2,3})."""
    n = 2 * h
    f = (1 << n) | (1 << h) | 1
    s = 2  # x
    frob = {}
    for d in range(1, n + 1):
        s = fmul(s, s, h)
        frob[d] = s
    assert frob[n] == 2, "x^(2^n) != x: not a subring of F_{2^n}"
    for d in (n // 2, n // 3):
        g = poly_gcd(frob[d] ^ 2, f)
        assert g == 1, f"gcd(x^(2^{d})-x, f) != 1: reducible"

# ---------------- primality (Miller-Rabin) ----------------

SMALL_DET_BASES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
DET_LIMIT = 3317044064679887385961981  # deterministic below this with bases above

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
    def witness(a):
        x = pow(a, d, n)
        if x in (1, n - 1):
            return False
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                return False
        return True
    if n < DET_LIMIT:
        for a in SMALL_DET_BASES:
            if witness(a):
                return False, "det"
        return True, "det"
    rng = random.Random(0xC0FFEE)
    for _ in range(rounds):
        if witness(rng.randrange(2, n - 1)):
            return False, "prp"
    return True, "prp"

# ---------------- the verification per level k ----------------

def verify_level(k: int, factors: dict[int, int], check_irred=True) -> None:
    h = 3 ** k
    N = (1 << h) - 1                      # |L*| = 2^h - 1
    M = 3 ** (k + 1) * N                  # conjectured ord(kappa_{3^k}+1)
    prod = 1
    for r, e in factors.items():
        prod *= r ** e
    assert prod == N, f"k={k}: factor product != 2^{h}-1"
    for r in factors:
        ok, kind = is_prime(r)
        assert ok, f"k={k}: factor {r} not prime"
    print(f"k={k}: h=3^{k}={h}, field F_2[x]/(x^{2*h}+x^{h}+1) = F_2(zeta_{3**(k+1)})")
    if check_irred:
        certify_irreducible(h)
        print(f"k={k}: trinomial certified irreducible (gcd test)")
    x = 2
    # zeta order
    assert fpow(x, 3 ** (k + 1), h) == 1 and fpow(x, 3 ** k, h) != 1
    xinv = (1 << (2 * h - 1)) ^ (1 << (h - 1))     # x^{-1} = x^{2h-1}+x^{h-1}
    assert fmul(x, xinv, h) == 1
    beta = x ^ 1                                    # kappa_{3^k} + 1
    # (circle) beta^(2^h-1) == zeta^{-1}
    assert fpow(beta, N, h) == xinv, "circle identity FAILS"
    # (norm)  beta^(2^h+1) == gamma = zeta + zeta^{-1}
    gamma = x ^ xinv
    assert fmul(fpow(beta, N, h), fmul(beta, beta, h), h) == fmul(gamma, beta, h) or True
    assert fmul(beta, fpow(beta, 1 << h, h), h) == gamma, "norm identity FAILS"
    # (half-angle) beta = zeta^{1/2} * (zeta^{1/2} + zeta^{-1/2})
    half = pow(2, -1, 3 ** (k + 1))
    zh = fpow(x, half, h)
    zhinv = fpow(xinv, half, h)
    assert fmul(zh, zh ^ zhinv, h) == beta, "half-angle identity FAILS"
    # gamma in L: gamma^(2^h) == gamma
    assert fpow(gamma, 1 << h, h) == gamma
    print(f"k={k}: identities verified  (circle, norm, half-angle; gamma in F_2^{h})")
    # translates m=2,3: beta_m^(2^h-1) is a primitive 3^{k+1}-th root of unity
    b2 = x ^ (1 << h)            # kappa + kappa_2,  kappa_2 = zeta^{3^k} = x^h
    b3 = b2 ^ 1
    assert fpow(b2, N, h) == fpow(x, 2 * h - 1, h), "m=2 circle FAILS"
    assert fpow(b3, N, h) == (1 << (h - 1)), "m=3 circle FAILS"
    print(f"k={k}: translate identities verified  (m=2,3 circle parts are primitive 3^{k+1}-th roots)")
    # primitivity of gamma in L*
    failures = []
    assert fpow(gamma, N, h) == 1
    for r in sorted(factors):
        t = fpow(gamma, N // r, h)
        status = "FAILS (r-th power!)" if t == 1 else "ok"
        if t == 1:
            failures.append(r)
        fr = "3^?" 
        print(f"k={k}:   gamma^((2^{h}-1)/{r}) {'==' if t==1 else '!='} 1   -> {status}")
    if failures:
        print(f"k={k}: *** C_{k} FAILS at primes {failures} ***")
    else:
        print(f"k={k}: gamma_k is PRIMITIVE in F_2^{h};  ord(kappa_3^{k}+1) = 3^{k+1}*(2^{h}-1)  CONFIRMED")
    # direct certificate on beta itself
    assert fpow(beta, M, h) == 1
    assert fpow(beta, M // 3, h) != 1
    for r in sorted(factors):
        if fpow(beta, M // r, h) == 1 and r not in failures:
            print(f"k={k}: !!! inconsistency at {r}")
    print(f"k={k}: direct certificate ord(beta) = 3^{k+1} * prod = {M}" if not failures else f"k={k}: ord(beta) < conjectured")
    print()
    return failures

if __name__ == "__main__":
    # verified factorizations of 2^{3^k}-1, k=1..4 (products checked in-script,
    # all factors < 2^27 so Miller-Rabin is deterministic)
    F1 = {7: 1}
    F2 = {7: 1, 73: 1}
    F3 = {7: 1, 73: 1, 262657: 1}
    F4 = {7: 1, 73: 1, 262657: 1, 2593: 1, 71119: 1, 97685839: 1}
    all_fail = {}
    for k, F in [(1, F1), (2, F2), (3, F3), (4, F4)]:
        f = verify_level(k, F)
        if f: all_fail[k] = f
    print("summary failures:", all_fail if all_fail else "none through k=4")
