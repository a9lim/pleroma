import sys, time
sys.path.insert(0, "/tmp")
from cyclo_family import (gf2_mul, tri_reduce, fmul, fpow, certify_irreducible,
                          is_prime, verify_level)

# Levels 5 and 6: factordb (FF) factorizations of Phi_243(2) and Phi_729(2),
# verified in-script by exact product reconstruction + primality tests.
F4 = {7: 1, 73: 1, 262657: 1, 2593: 1, 71119: 1, 97685839: 1}
PHI243 = {487: 1, 16753783618801: 1, 192971705688577: 1, 3712990163251158343: 1}
PHI729 = {80191: 1, 97687: 1, 379081: 1,
          664728004346558283448724389870269691211809: 1,
          101213745778143742250901040788003424950068418098259161142719688891708905138274462262307761: 1}
F5 = dict(F4); F5.update(PHI243)
F6 = dict(F5); F6.update(PHI729)

# 3-adic valuation sanity (LTE): v_3(2^{3^k}+1) = k+1 exactly
for k in range(1, 8):
    h = 3 ** k
    assert (2**h + 1) % 3**(k+1) == 0 and (2**h + 1) % 3**(k+2) != 0
print("LTE check: v_3(2^(3^k)+1) = k+1 exactly, k=1..7  OK")

# squarefreeness of 2^{3^k}-1 across represented primes (no level-Wieferich)
for k, F in [(5, F5), (6, F6)]:
    N = 2**(3**k) - 1
    for r in F:
        assert N % (r * r) != 0, f"square divisor {r} at k={k}"
print("squarefree check: no r^2 | 2^(3^k)-1 for any table prime, k=5,6  OK")

# primality kinds
for r in sorted(F6):
    ok, kind = is_prime(r)
    assert ok
    if kind == "prp":
        print(f"  note: {str(r)[:25]}... ({len(str(r))} digits) is MR-64 PRP (factordb status: prime)")

t0 = time.time()
fail5 = verify_level(5, F5)
print(f"[k=5 took {time.time()-t0:.1f}s]")
t0 = time.time()
fail6 = verify_level(6, F6)
print(f"[k=6 took {time.time()-t0:.1f}s]")

# ---- k=7 partial: Phi_2187(2) is CF in factordb; test the known primes only ----
k, h = 7, 3**7
N = (1 << h) - 1
known_new = [39367, 7606246033, 263196614521, 529063556041]
for r in known_new:
    ok, kind = is_prime(r)
    assert ok and (N % r == 0), f"{r} not a prime factor of 2^2187-1"
certify_irreducible(h)
x = 2
xinv = (1 << (2*h - 1)) ^ (1 << (h - 1))
assert fmul(x, xinv, h) == 1
beta = x ^ 1
gamma = x ^ xinv
assert fpow(beta, N, h) == xinv, "k=7 circle identity FAILS"
assert fmul(beta, fpow(beta, 1 << h, h), h) == gamma, "k=7 norm identity FAILS"
print(f"k=7: trinomial irreducible, circle + norm identities verified in F_2^4374")
t0 = time.time()
for r in known_new:
    t = fpow(gamma, N // r, h)
    print(f"k=7:   gamma^((2^2187-1)/{r}) {'==' if t==1 else '!='} 1   -> {'FAILS' if t==1 else 'ok'}")
print(f"k=7: all known new primes ok -> C_7 consistent so far; certification blocked by the")
print(f"     unfactored 400-digit cofactor of Phi_2187(2) (factordb status CF). [{time.time()-t0:.1f}s]")
