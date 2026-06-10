import math
from fractions import Fraction

def ldl(gram):
    n = len(gram)
    d = [0.0]*n
    l = [[0.0]*n for _ in range(n)]
    for j in range(n):
        dj = float(gram[j][j])
        for k in range(j):
            dj -= l[j][k]*l[j][k]*d[k]
        d[j] = dj
        l[j][j] = 1.0
        for i in range(j+1, n):
            s = float(gram[i][j])
            for k in range(j):
                s -= l[i][k]*l[j][k]*d[k]
            l[i][j] = s/dj if dj != 0.0 else 0.0
    u = [[0.0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            u[i][j] = l[j][i]
    return d, u

def norm(gram, x):
    n = len(gram)
    return sum(gram[i][j]*x[i]*x[j] for i in range(n) for j in range(n))

def fp_search(gram, i, bound, d, u, eps, tail, x, out):
    if i == 0:
        q = norm(gram, x)
        if 0 < q <= bound:
            out.append(list(x))
        return
    idx = i-1
    center = 0.0
    for j in range(i, len(d)):
        center += u[idx][j]*float(x[j])
    remaining = float(bound) - tail
    if remaining < -eps:
        return
    radius = math.sqrt(max(remaining, 0.0)/d[idx]) + eps
    lo = math.ceil(-center - radius)
    hi = math.floor(-center + radius)
    for xi in range(lo, hi+1):
        x[idx] = xi
        coord = float(xi) + center
        fp_search(gram, idx, bound, d, u, eps, tail + d[idx]*coord*coord, x, out)
    x[idx] = 0

def exact_box_count(gram, bound):
    # mimic short_vectors_exact_bounded box size with exact rationals
    n = len(gram)
    import sympy
    M = sympy.Matrix(gram)
    inv = M.inv()
    count = 1
    for i in range(n):
        radius2 = Fraction(bound) * Fraction(inv[i,i].p, inv[i,i].q)
        # ceil sqrt of rational
        r = math.isqrt(radius2.numerator // radius2.denominator)
        while Fraction(r*r) < radius2:
            r += 1
        count *= (2*r+1)
    return count

def test_K(K, M):
    g = [[2*K, -(K-1), -(K-1)],[-(K-1), 2*K, -(K-1)],[-(K-1), -(K-1), 2*K]]
    bound = 6*M*M
    d, u = ldl(g)
    # exact d2 via Fractions
    fd0 = Fraction(2*K)
    fl10 = Fraction(-(K-1), 2*K)
    fd1 = Fraction(2*K) - fl10*fl10*fd0
    fl20 = fl10
    fl21 = (Fraction(-(K-1)) - fl20*fl10*fd0)/fd1
    fd2 = Fraction(2*K) - fl20*fl20*fd0 - fl21*fl21*fd1
    eps = 1e-9*max(float(bound),1.0) + 1e-9
    radius_top = math.sqrt(float(bound)/d[2]) + eps
    hi = math.floor(radius_top)
    true_radius = math.sqrt(bound/float(fd2))
    return d[2], float(fd2), hi, true_radius

# scan K for d2_float overestimating enough that hi < M
M = 70
found = []
for K in [10**12 + t for t in range(0, 200)]:
    d2f, d2t, hi, tr = test_K(K, M)
    if hi < M:
        found.append((K, d2f, d2t, hi))
print("M =", M, "bound =", 6*M*M)
print("hits:", len(found))
for K, d2f, d2t, hi in found[:5]:
    print(f"K={K}: d2_float={d2f!r} d2_true={d2t!r} hi={hi}")

# diagnostics
for K in [10**12, 10**13, 10**14, 4*10**14, 9*10**14, 2*10**15]:
    d2f, d2t, hi, tr = test_K(K, M)
    eps = 1e-9*6*M*M + 1e-9
    print(f"K={K:.2e} d2_float={d2f:.12f} d2_true={d2t:.12f} err={d2f-d2t:+.3e} hi={hi} true_radius={tr:.9f} eps={eps:.3e}")
