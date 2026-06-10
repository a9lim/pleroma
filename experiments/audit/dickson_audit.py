import itertools
from arf_audit import nim_mul, nim_inv

def mat_mul(A, B):
    n = len(A)
    return [[__import__('functools').reduce(lambda a, k: a ^ nim_mul(A[i][k_]
            , B[k_][j]), [], 0) for j in range(n)] for i in range(n)]

def mmul(A, B):
    n = len(A)
    C = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            acc = 0
            for k in range(n):
                acc ^= nim_mul(A[i][k], B[k][j])
            C[i][j] = acc
    return C

def nim_rank(rows):
    rows = [list(r) for r in rows]
    nrows = len(rows)
    if nrows == 0: return 0
    ncols = len(rows[0])
    pr = 0
    for col in range(ncols):
        p = next((r for r in range(pr, nrows) if rows[r][col] != 0), None)
        if p is None: continue
        rows[pr], rows[p] = rows[p], rows[pr]
        inv = nim_inv(rows[pr][col])
        for c in range(col, ncols):
            rows[pr][c] = nim_mul(rows[pr][c], inv)
        for r in range(nrows):
            if r != pr and rows[r][col] != 0:
                f = rows[r][col]
                for c in range(col, ncols):
                    rows[r][c] ^= nim_mul(f, rows[pr][c])
        pr += 1
        if pr == nrows: break
    return pr

def dickson(g):
    n = len(g)
    m = [list(r) for r in g]
    for i in range(n):
        m[i][i] ^= 1
    return nim_rank(m) % 2

def qval(v, q, bmat):
    n = len(v); acc = 0
    for i in range(n):
        acc ^= nim_mul(nim_mul(v[i], v[i]), q[i])
        for j in range(i+1, n):
            acc ^= nim_mul(nim_mul(v[i], v[j]), bmat[i][j])
    return acc

def orthogonal_group(q, bpairs, m):
    """all g in GL_n(F_{2^m}) with Q(gv) = Q(v) for all v  (n small)"""
    n = len(q)
    bmat = [[0]*n for _ in range(n)]
    for (i, j), v in bpairs.items():
        bmat[i][j] = v; bmat[j][i] = v
    Fq = 1 << m
    vecs = list(itertools.product(range(Fq), repeat=n))
    group = []
    # g given by images of basis vectors (columns)
    for cols in itertools.product(vecs, repeat=n):
        # build matrix with columns cols: g[i][j] = cols[j][i]
        g = [[cols[j][i] for j in range(n)] for i in range(n)]
        if nim_rank(g) != n: continue
        ok = True
        for v in vecs:
            gv = [0]*n
            for i in range(n):
                acc = 0
                for j in range(n):
                    acc ^= nim_mul(g[i][j], v[j])
                gv[i] = acc
            if qval(gv, q, bmat) != qval(list(v), q, bmat):
                ok = False; break
        if ok: group.append(g)
    return group

# O+_4(2): Q = x0x1 + x2x3 over F2
q = [0,0,0,0]; bp = {(0,1):1, (2,3):1}
G = orthogonal_group(q, bp, 1)
print("|O+_4(2)| =", len(G), "(theory: 72)")
ker = sum(1 for g in G if dickson(g) == 0)
print("|ker D| =", ker, "(theory: 36)")
# homomorphism check on all pairs
bad = 0
for a in G:
    for b in G:
        if dickson(mmul(a, b)) != (dickson(a) ^ dickson(b)):
            bad += 1
print("homomorphism failures O+_4(2):", bad)

# O-_4(2): Q = x0x1 + x2^2 + x2x3 + x3^2 over F2
q = [0,0,1,1]; bp = {(0,1):1, (2,3):1}
G = orthogonal_group(q, bp, 1)
print("|O-_4(2)| =", len(G), "(theory: 120)")
ker = sum(1 for g in G if dickson(g) == 0)
print("|ker D| =", ker, "(theory: 60)")
bad = 0
for a in G:
    for b in G:
        if dickson(mmul(a, b)) != (dickson(a) ^ dickson(b)):
            bad += 1
print("homomorphism failures O-_4(2):", bad)

# O_2 over F4 (hyperbolic plane, nim-field entries)
q = [0,0]; bp = {(0,1):1}
G = orthogonal_group(q, bp, 2)
print("|O(H) over F4| =", len(G), "(theory: 6 = dihedral)")
ker = sum(1 for g in G if dickson(g) == 0)
print("|ker D| =", ker, "(theory: 3)")
bad = 0
for a in G:
    for b in G:
        if dickson(mmul(a, b)) != (dickson(a) ^ dickson(b)):
            bad += 1
print("homomorphism failures over F4:", bad)
