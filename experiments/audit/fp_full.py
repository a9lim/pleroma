import math
from fractions import Fraction
exec(open('/tmp/fp_emulate.py').read().split('# scan')[0])  # reuse ldl, fp_search, norm

K, a, b, c = 286203218045748, 7, 7, 10
M = 70
g = [[2*K, -(K-a), -(K-b)],
     [-(K-a), 2*K, -(K-c)],
     [-(K-b), -(K-c), 2*K]]
nrm = 2*(a+b+c)
bound = nrm*M*M
print("K =", K, "(a,b,c) =", (a,b,c), "bound =", bound, "norm(1,1,1) =", norm(g,[1,1,1]))
print("entries exactly f64-representable:", all(abs(x) < 2**53 and float(x) == x for row in g for x in row))

# 1. positive definite (leading minors)
m1 = g[0][0]
m2 = g[0][0]*g[1][1]-g[0][1]**2
m3 = (g[0][0]*(g[1][1]*g[2][2]-g[1][2]*g[2][1])
    - g[0][1]*(g[1][0]*g[2][2]-g[1][2]*g[2][0])
    + g[0][2]*(g[1][0]*g[2][1]-g[1][1]*g[2][0]))
print("minors:", m1>0, m2>0, m3>0, "det =", m3)

# 2. exact box count (mimic short_vectors_exact_bounded): r_i = ceil(sqrt(bound*inv_ii))
import itertools
def cof(i):
    idx = [k for k in range(3) if k != i]
    return g[idx[0]][idx[0]]*g[idx[1]][idx[1]] - g[idx[0]][idx[1]]*g[idx[1]][idx[0]]
count = 1
for i in range(3):
    r2 = Fraction(bound) * Fraction(cof(i), m3)
    r = math.isqrt(r2.numerator // r2.denominator)
    while Fraction(r*r) < r2: r += 1
    count *= 2*r + 1
print("exact box count:", count, "> 2e6:", count > 2_000_000)

# 3. size reduction fixed point: round_div_nearest(g_ij, g_ii) == 0 for i<j, diagonals nondecreasing
def round_div_nearest(p, q):
    # round-half-away-from-zero on p/q, q>0 (check rust impl separately)
    fl = p // q
    rem = p - fl*q
    if 2*rem >= q: fl += 1
    return fl
ks = [round_div_nearest(g[i][j], g[i][i]) for i in range(3) for j in range(i+1,3)]
print("shear ks:", ks, "swaps needed:", any(g[i+1][i+1] < g[i][i] for i in range(2)))

# 4. full fp_search emulation
d, u = ldl(g)
print("d =", d)
det_true = m3
d2_true = Fraction(m3, m2)
print("d2 true =", float(d2_true), " d2 float =", d[2], " abs err =", d[2]-float(d2_true))
eps = 1e-9*max(float(bound),1.0)+1e-9
out = []
x = [0,0,0]
fp_search(g, 3, bound, d, u, eps, 0.0, x, out)
expected = [[m,m,m] for m in range(-M, M+1) if m != 0]
got = sorted(out)
print("vectors found:", len(got), " expected:", len(expected))
missing = [v for v in expected if v not in got]
spurious = [v for v in got if v not in expected]
print("missing:", missing)
print("spurious:", spurious)
print("norm of missing:", [norm(g, v) for v in missing], "<= bound:", bound)
