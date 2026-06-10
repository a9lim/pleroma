import itertools, random, functools
random.seed(7)

# ---------------- nim arithmetic ----------------
@functools.lru_cache(maxsize=None)
def nim_mul(a, b):
    if a > b: a, b = b, a   # a <= b
    if a == 0: return 0
    if a == 1: return b
    # largest Fermat power F = 2^{2^k} with F <= b < F^2
    F = 2
    while F * F <= b:
        F = F * F
    bh, bl = divmod(b, F)
    if a < F:
        # (bh F + bl) * a = (bh a) F + bl a   (x<F times F is ordinary product)
        return nim_mul(bh, a) * F ^ nim_mul(bl, a)
    ah, al = divmod(a, F)
    ahbh = nim_mul(ah, bh)
    cross = nim_mul(ah, bl) ^ nim_mul(al, bh)
    albl = nim_mul(al, bl)
    # F*F = (3/2) F  -> ahbh * F^2 = ahbh*(F/2)  ^  ahbh * F... careful:
    # F (x) F = F + F/2 (nim), so ahbh (x) F^2 = nim_mul(ahbh, F) ^ nim_mul(ahbh, F//2)
    #   nim_mul(ahbh,F): ahbh < F so ordinary ahbh*F
    return (cross ^ ahbh) * F ^ albl ^ nim_mul(ahbh, F // 2)

# sanity: known small nim table
assert nim_mul(2,2)==3 and nim_mul(2,3)==1 and nim_mul(3,3)==2
assert nim_mul(4,4)==6 and nim_mul(2,4)==8 and nim_mul(4,5)==2

def nim_inv(x):
    # brute force within F_{2^m} containing x
    m = 1
    while x >= (1 << m): m *= 2
    for y in range(1, 1 << m):
        if nim_mul(x, y) == 1: return y
    raise ValueError

def nim_square(x): return nim_mul(x, x)

def nim_trace(x, m):
    acc = x; t = x
    for _ in range(1, m):
        t = nim_square(t); acc ^= t
    return acc

def min_field_degree(maxv):
    m = 1
    while m < 128:
        if maxv < (1 << m): return m
        m <<= 1
    return 128

# ---------------- port of arf_nimber ----------------
def qf(v, q, bmat):
    n = len(v); acc = 0
    for i in range(n):
        acc ^= nim_mul(nim_mul(v[i], v[i]), q[i])
        for j in range(i+1, n):
            acc ^= nim_mul(nim_mul(v[i], v[j]), bmat[i][j])
    return acc

def bf(u, v, bmat):
    n = len(u); acc = 0
    for i in range(n):
        for j in range(i+1, n):
            cross = nim_mul(u[i], v[j]) ^ nim_mul(u[j], v[i])
            acc ^= nim_mul(cross, bmat[i][j])
    return acc

def arf_nimber(q, bpairs):
    n = len(q)
    bmat = [[0]*n for _ in range(n)]
    for (i, j), v in bpairs.items():
        bmat[i][j] = v; bmat[j][i] = v
    maxv = max(q + [0] + [x for row in bmat for x in row])
    m = min_field_degree(maxv)
    vectors = [[1 if k == i else 0 for k in range(n)] for i in range(n)]
    s = 0; pairs = 0; radical_dim = 0; radical_aniso = False
    while vectors:
        a = vectors.pop()
        pos = next((k for k, w in enumerate(vectors) if bf(a, w, bmat) != 0), None)
        if pos is not None:
            braw = vectors[pos]
            vectors[pos] = vectors[-1]; vectors.pop()  # swap_remove
            c = bf(a, braw, bmat)
            ci = nim_inv(c)
            b = [nim_mul(ci, x) for x in braw]
            for idx, w in enumerate(vectors):
                wb = bf(w, b, bmat); wa = bf(w, a, bmat)
                nw = list(w)
                if wb: nw = [x ^ nim_mul(wb, y) for x, y in zip(nw, a)]
                if wa: nw = [x ^ nim_mul(wa, y) for x, y in zip(nw, b)]
                vectors[idx] = nw
            s ^= nim_mul(qf(a, q, bmat), qf(b, q, bmat))
            pairs += 1
        else:
            radical_dim += 1
            if qf(a, q, bmat) != 0: radical_aniso = True
    return dict(arf=nim_trace(s, m), rank=2*pairs, radical_dim=radical_dim,
                radical_anisotropic=radical_aniso, m=m)

# ---------------- ground truth: zero counts ----------------
def zero_count(q, bpairs, m):
    n = len(q)
    bmat = [[0]*n for _ in range(n)]
    for (i, j), v in bpairs.items():
        bmat[i][j] = v; bmat[j][i] = v
    Fq = 1 << m
    cnt = 0
    for v in itertools.product(range(Fq), repeat=n):
        if qf(list(v), q, bmat) == 0: cnt += 1
    return cnt

def expected_zeros(qsize, npairs, arf):
    # nondegenerate rank-2n form over F_q: q^{2n-1} +- (q-1) q^{n-1}
    base = qsize ** (2*npairs - 1)
    corr = (qsize - 1) * qsize ** (npairs - 1)
    return base + corr if arf == 0 else base - corr

# ---- check 1: exhaustive F2 nonsingular forms n=2 and n=4 (subset of polar configs)
print("== F2 exhaustive checks ==")
fails = 0
for n, bconfigs in [(2, [ {(0,1):1} ]),
                    (4, [ {(0,1):1,(2,3):1}, {(0,1):1,(1,2):1,(2,3):1},
                          {(0,1):1,(0,2):1,(0,3):1,(1,2):1,(1,3):1,(2,3):1},
                          {(0,2):1,(1,3):1}, {(0,3):1,(1,2):1,(0,1):1} ])]:
    for bp in bconfigs:
        for qv in itertools.product([0,1], repeat=n):
            r = arf_nimber(list(qv), bp)
            if r['radical_dim'] != 0:   # need nonsingular for the count formula
                continue
            zc = zero_count(list(qv), bp, 1)
            exp = expected_zeros(2, r['rank']//2, r['arf'])
            if zc != exp:
                fails += 1
                print("FAIL F2", qv, bp, r, "zeros", zc, "expected", exp)
print("F2 nonsingular check fails:", fails)

# ---- check 2: random F4 forms, 2 and 4 vars
print("== F4 random checks ==")
fails = 0
for trial in range(200):
    n = random.choice([2, 4])
    q = [random.randrange(4) for _ in range(n)]
    bp = {}
    for i in range(n):
        for j in range(i+1, n):
            v = random.randrange(4)
            if v: bp[(i,j)] = v
    # force at least the field to be F4 (else field-of-def shrinks; that's fine too)
    r = arf_nimber(q, bp)
    if r['radical_dim'] != 0: continue
    m = r['m']
    zc = zero_count(q, bp, m)
    exp = expected_zeros(1 << m, r['rank']//2, r['arf'])
    if zc != exp:
        fails += 1
        print("FAIL F4", q, bp, r, "zeros", zc, "expected", exp)
print("F4 random check fails:", fails)

# ---- check 3: F16 planes (2 vars)
print("== F16 random planes ==")
fails = 0
for trial in range(60):
    q = [random.randrange(16) for _ in range(2)]
    bp = {(0,1): random.randrange(1,16)}
    r = arf_nimber(q, bp)
    if r['radical_dim'] != 0: continue
    m = r['m']
    zc = zero_count(q, bp, m)
    exp = expected_zeros(1 << m, 1, r['arf'])
    if zc != exp:
        fails += 1
        print("FAIL F16", q, bp, r, "zeros", zc, "expected", exp)
print("F16 plane check fails:", fails)

# ---- check 4: the documented F4 test values
r1 = arf_nimber([2,3], {(0,1):1}); r2 = arf_nimber([2,2], {(0,1):1})
print("F4 [2,3] arf:", r1['arf'], " (test says 0)   F4 [2,2] arf:", r2['arf'], " (test says 1)")

# ---- check 5: cross-subfield Witt-class addition
print("== cross-subfield Witt/BW additivity ==")
A = ([1,1], {(0,1):1})          # F2 anisotropic plane: arf over F2
B = ([2,2], {(0,1):1})          # F4 anisotropic plane: arf over F4
ra = arf_nimber(*A); rb = arf_nimber(*B)
# direct sum
qs = A[0] + B[0]
bp = dict(A[1]); bp[(2,3)] = B[1][(0,1)]
rsum = arf_nimber(qs, bp)
print("arf(A) =", ra['arf'], "(field deg", ra['m'], ")  arf(B) =", rb['arf'], "(field deg", rb['m'], ")")
print("arf(A perp B) =", rsum['arf'], "(field deg", rsum['m'], ")   XOR law predicts", ra['arf'] ^ rb['arf'])
# ground truth via zero count of the rank-4 form over F4
zc = zero_count(qs, bp, 2)
print("zero count of A perp B over F4:", zc,
      " Arf0 ->", expected_zeros(4,2,0), " Arf1 ->", expected_zeros(4,2,1))
