import sys
sys.path.insert(0, "/Users/a9lim/Work/ogdoad/experiments")
from ordinal_excess_probe import TermAlgebra

ONE = frozenset((0,))

# ---- k=1: components (2,3), F_{2^6}, kappa_3 ----
alg = TermAlgebra((2, 3))
k3 = frozenset((alg.basis[alg.index[3]],))
beta1 = frozenset((0, alg.basis[alg.index[3]]))      # kappa_3 + 1
beta2 = frozenset((1, alg.basis[alg.index[3]]))      # kappa_3 + 2   (finite 2 = kappa_2 = term 1)
beta3 = frozenset((0, 1, alg.basis[alg.index[3]]))   # kappa_3 + 3
assert alg.power(k3, 9) == ONE and alg.power(k3, 3) != ONE, "ord(kappa_3) != 9"
print("k=1: ord(kappa_3) = 9  [primitive 9th root of unity]")
print("k=1: ord(kappa_3+1) =", alg.order(beta1), " (conjecture: 9*7 = 63)")
print("k=1: ord(kappa_3+2) =", alg.order(beta2))
print("k=1: ord(kappa_3+3) =", alg.order(beta3))
# circle identity: beta1^(2^3-1) == kappa_3^{-1} = kappa_3^8
assert alg.power(beta1, 7) == alg.power(k3, 8), "circle identity fails k=1"
# norm identity: beta1^(2^3+1) == kappa_3 + kappa_3^{-1}
gamma1 = frozenset(k3 ^ alg.power(k3, 8))
assert alg.power(beta1, 9) == gamma1, "norm identity fails k=1"
print("k=1: beta^(2^3-1) == kappa_3^{-1} and Norm(beta) == kappa_3+kappa_3^{-1}  OK")
print("k=1: ord(gamma_1) =", alg.order(gamma1), " (conjecture: 2^3-1 = 7)")
print()

# ---- k=2: components (2,3,9), F_{2^18}, kappa_9 ----
alg = TermAlgebra((2, 3, 9))
k9 = frozenset((alg.basis[alg.index[9]],))
beta1 = frozenset((0, alg.basis[alg.index[9]]))
beta2 = frozenset((1, alg.basis[alg.index[9]]))
beta3 = frozenset((0, 1, alg.basis[alg.index[9]]))
assert alg.power(k9, 27) == ONE and alg.power(k9, 9) != ONE, "ord(kappa_9) != 27"
print("k=2: ord(kappa_9) = 27  [primitive 27th root of unity]")
print("k=2: ord(kappa_9+1) =", alg.order(beta1), " (conjecture: 27*511 = 13797)")
print("k=2: ord(kappa_9+2) =", alg.order(beta2))
print("k=2: ord(kappa_9+3) =", alg.order(beta3))
assert alg.power(beta1, 511) == alg.power(k9, 26), "circle identity fails k=2"
gamma2 = frozenset(k9 ^ alg.power(k9, 26))
assert alg.power(beta1, 513) == gamma2, "norm identity fails k=2"
print("k=2: beta^(2^9-1) == kappa_9^{-1} and Norm(beta) == kappa_9+kappa_9^{-1}  OK")
print("k=2: ord(gamma_2) =", alg.order(gamma2), " (conjecture: 2^9-1 = 511)")
# tower checks: kappa_9^3 == kappa_3, kappa_3^3 == kappa_2, gamma_1 == gamma_2^3 + gamma_2
k3 = frozenset((alg.basis[alg.index[3]],))
k2 = frozenset((alg.basis[alg.index[2]],))
assert alg.power(k9, 3) == k3 and alg.power(k3, 3) == k2
g2cubed = alg.power(gamma2, 3)
gamma1_in_18 = frozenset(g2cubed ^ gamma2)
assert alg.power(gamma1_in_18, 7) == ONE and gamma1_in_18 != ONE, "norm-tower image not of order 7"
print("k=2: kappa_9^3 == kappa_3, kappa_3^3 == kappa_2, and gamma_2^3+gamma_2 has order 7 (== gamma_1 up to conjugacy)  OK")
