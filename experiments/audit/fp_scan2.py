import math
from fractions import Fraction

def ldl3(g):
    # faithful transcription for n=3
    d0 = float(g[0][0])
    l10 = float(g[1][0]) / d0
    l20 = float(g[2][0]) / d0
    d1 = float(g[1][1]) - l10*l10*d0
    s21 = float(g[2][1]) - l20*l10*d0
    l21 = s21 / d1
    d2 = float(g[2][2]) - l20*l20*d0 - l21*l21*d1
    return d0, d1, d2, l10, l20, l21

def det3i(g):
    return (g[0][0]*(g[1][1]*g[2][2]-g[1][2]*g[2][1])
          - g[0][1]*(g[1][0]*g[2][2]-g[1][2]*g[2][0])
          + g[0][2]*(g[1][0]*g[2][1]-g[1][1]*g[2][0]))

I128MAX = 2**127 - 1

def bareiss_fits(g):
    # emulate bareiss_det checked muls for 3x3
    a = [row[:] for row in g]
    prev = 1
    for k in range(2):
        for i in range(k+1,3):
            for j in range(k+1,3):
                p1 = a[i][j]*a[k][k]
                p2 = a[i][k]*a[k][j]
                if abs(p1) > I128MAX or abs(p2) > I128MAX: return False
                a[i][j] = (p1-p2)//prev
        prev = a[k][k]
    return True

def check(K, variant, M):
    if variant == 0:   # (a,b,c)=(1,0,0)
        g = [[2*K, -(K-1), -K],[-(K-1), 2*K, -K],[-K, -K, 2*K]]
    elif variant == 1: # (a,b,c)=(0,0,1)
        g = [[2*K, -K, -K],[-K, 2*K, -(K-1)],[-K, -(K-1), 2*K]]
    elif variant == 2: # (a,b,c)=(0,1,0)
        g = [[2*K, -K, -(K-1)],[-K, 2*K, -K],[-(K-1), -K, 2*K]]
    else:              # s3=1, abc=0
        g = [[2*K, -K, -K],[-K, 2*K, -K],[-K, -K, 2*K+2]]
    bound = 2*M*M
    det = det3i(g)
    if det <= 0: return None
    minor = g[0][0]*g[1][1]-g[0][1]*g[1][0]
    d2t = Fraction(det, minor)
    # vector (M,M,M) must have norm exactly bound
    nrm111 = sum(g[i][j] for i in range(3) for j in range(3))
    if nrm111*M*M != bound: return None
    d0,d1,d2,l10,l20,l21 = ldl3(g)
    eps = 1e-9*float(bound) + 1e-9
    hi = math.floor(math.sqrt(float(bound)/d2) + eps)
    if hi >= M: return None
    # confirm exact path would refuse: box > 2e6
    count = 1
    for i in range(3):
        idx = [k for k in range(3) if k != i]
        cof = g[idx[0]][idx[0]]*g[idx[1]][idx[1]] - g[idx[0]][idx[1]]*g[idx[1]][idx[0]]
        r2 = Fraction(bound*cof, det)
        r = math.isqrt(r2.numerator//r2.denominator)
        while Fraction(r*r) < r2: r += 1
        count *= 2*r+1
    if count <= 2_000_000: return None
    if not bareiss_fits(g): return None
    # shears all zero? |g_ij| <= g_ii/2 with round-half-up handling: -K/2K rounds to 0 (half rounds up)
    # round_div_nearest(p,q): q=2K; for p=-K: div_euclid=-1 rem=K, 2K>=2K -> 0 OK; p=-(K-1): 0 OK
    return (g, bound, count, d2, float(d2t), hi)

found = []
tried = 0
for M in (63, 64, 70, 80, 90, 100):
    for variant in range(4):
        K = 1_000_000_000
        while K <= 1_650_000_000:
            tried += 1
            r = check(K, variant, M)
            if r:
                found.append((K, variant, M, r))
                if len(found) >= 3: break
            K += 1
        if found: break
    if found: break
print("tried:", tried, "found:", len(found))
for K, variant, M, (g, bound, count, d2f, d2t, hi) in found:
    print(f"K={K} variant={variant} M={M} bound={bound} box={count} d2f={d2f!r} d2t={d2t!r} hi={hi}")
