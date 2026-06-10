# Teichmuller reps in Z_5 mod 5^6 via t <- t^5 iterated (as Qp::teichmuller does, K=6 iterations)
M = 5**6
def teich(a, p=5, K=6):
    t = a % M
    for _ in range(K):
        t = pow(t, p, M)
    return t

t1, t2, t3 = teich(1), teich(2), teich(3)
print("tau(1)=", t1, "tau(2)=", t2, "tau(3)=", t3)
print("tau(2)^4 mod 5^6 =", pow(t2,4,M), " (should be 1)")
print("tau(3)^4 mod 5^6 =", pow(t3,4,M), " (should be 1)")
# Gauss<Qp<5,6>> teichmuller of r=1+tbar, s=1+2tbar  (coefficientwise lift)
# tau(r)*tau(s) = 1 + (tau(1)+tau(2)) t + tau(1)tau(2) t^2
# tau(r*s) = tau(1 + 3 tbar + 2 tbar^2) = 1 + tau(3) t + tau(2) t^2
lin_prod = (t1 + t2) % M
lin_lift = t3
print("t-coefficient of tau(r)tau(s):", lin_prod)
print("t-coefficient of tau(rs):     ", lin_lift)
print("equal?", lin_prod == lin_lift)
# difference and its 5-adic valuation
d = (lin_prod - lin_lift) % M
v = 0
dd = d
while dd and dd % 5 == 0:
    dd //= 5; v += 1
print("difference:", d, "5-adic valuation:", v, "(nonzero mod 5^6 => detectable at precision 6)")
