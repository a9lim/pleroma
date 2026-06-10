import math, random
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

def det3(g):
    return (g[0][0]*(g[1][1]*g[2][2]-g[1][2]*g[2][1])
          - g[0][1]*(g[1][0]*g[2][2]-g[1][2]*g[2][0])
          + g[0][2]*(g[1][0]*g[2][1]-g[1][1]*g[2][0]))

def make(K, a, b, c):
    return [[2*K, -(K-a), -(K-b)],
            [-(K-a), 2*K, -(K-c)],
            [-(K-b), -(K-c), 2*K]]

random.seed(42)
hits = []
trials = 0
M = 70
while trials < 300000 and len(hits) < 10:
    trials += 1
    K = random.randint(10**11, 4*10**15)
    a = random.randint(1, 12); b = random.randint(1, 12); c = random.randint(1, 12)
    g = make(K, a, b, c)
    # norm of (1,1,1) = 2(a+b+c); bound for M copies:
    nrm111 = 2*(a+b+c)
    bound = nrm111 * M * M
    # exact d2 = det(G)/det(2x2 leading minor)
    det = det3(g)
    if det <= 0:
        continue
    minor = g[0][0]*g[1][1] - g[0][1]*g[1][0]
    d2_true = Fraction(det, minor)
    d, u = ldl(g)
    eps = 1e-9*max(float(bound), 1.0) + 1e-9
    radius = math.sqrt(float(bound)/d[2]) + eps
    hi = math.floor(radius)
    # true radius^2 = bound/d2_true ; check M is genuinely inside true box
    inside = Fraction(bound) / d2_true >= M*M
    if inside and hi < M:
        hits.append((K, a, b, c, d[2], float(d2_true), hi))
print("trials:", trials, "hits:", len(hits))
for h in hits:
    print(h)
