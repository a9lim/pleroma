# Brute-force Gold form Q_a(x) = Tr(x^{2^a+1}) over F_{2^m}, m = 4, 8, 16.
# Verify: rank of polar form = m - gcd(2a, m), zero counts, Arf via count, radical isotropy.
from math import gcd

IRRED = {2: 0b111, 4: 0b10011, 8: 0b100011011, 16: 0b10000000000101101, 32: 0b100000000000000000000000011000101}

def mk_field(m):
    mod = IRRED[m]
    def mul(a, b):
        r = 0
        while b:
            if b & 1: r ^= a
            b >>= 1
            a <<= 1
            if a >> m & 1: a ^= mod
        return r
    return mul

def field_ops(m):
    mul = mk_field(m)
    def sq(a): return mul(a, a)
    def tr(a):
        t, c = 0, a
        for _ in range(m):
            t ^= c
            c = sq(c)
        # t should be 0 or 1
        return t & 1 if t in (0,1) else None
    return mul, sq, tr

def gold(m, a):
    mul, sq, tr = field_ops(m)
    def Q(x):
        y = x
        for _ in range(a): y = sq(y)   # x^{2^a}
        return tr(mul(x, y))
    N = 1 << m
    qv = [Q(x) for x in range(N)]
    assert all(v in (0,1) for v in qv), "trace not in F2!"
    zeros = qv.count(0)
    # polar form on basis e_i = 1<<i
    B = [[qv[(1<<i)^(1<<j)] ^ qv[1<<i] ^ qv[1<<j] for j in range(m)] for i in range(m)]
    # rank over F2
    M = [int(''.join(map(str, row)), 2) for row in B]
    rank = 0
    rows = M[:]
    for col in range(m):
        piv = None
        for r in range(rank, m):
            if rows[r] >> (m-1-col) & 1: piv = r; break
        if piv is None: continue
        rows[rank], rows[piv] = rows[piv], rows[rank]
        for r in range(m):
            if r != rank and (rows[r] >> (m-1-col) & 1): rows[r] ^= rows[rank]
        rank += 1
    # radical: x with B(x, e_j)=0 for all j; B(x,y) = qv[x^y]^qv[x]^qv[y] (bilinear)
    rad = [x for x in range(N) if all((qv[x^(1<<j)] ^ qv[x] ^ qv[1<<j]) == 0 for j in range(m))]
    rad_dim = len(rad).bit_length() - 1
    q_on_rad = [qv[x] for x in rad]
    rad_isotropic = all(v == 0 for v in q_on_rad)
    # predicted zero count if Arf=arf and radical isotropic:
    # zeros = 2^rad * (2^{rank-1} + (-1)^arf 2^{rank/2-1})
    pred = {}
    for arf in (0,1):
        if rank > 0:
            pred[arf] = (1<<rad_dim) * ((1<<(rank-1)) + (-1)**arf * (1<<(rank//2 - 1)))
    return dict(m=m, a=a, rank=rank, expected_rank=m - gcd(2*a, m), zeros=zeros,
                rad_dim=rad_dim, rad_isotropic=rad_isotropic, pred0=pred.get(0), pred1=pred.get(1),
                bias=zeros - (1<<(m-1)))

for (m,a) in [(4,1),(4,2),(8,1),(8,2),(8,3),(16,1),(16,4)]:
    r = gold(m,a)
    arf = 1 if r['zeros'] == r['pred1'] else (0 if r['zeros'] == r['pred0'] else '??')
    print(f"m={m:2d} a={a} rank={r['rank']:2d} (exp {r['expected_rank']:2d}) rad_dim={r['rad_dim']} "
          f"rad_iso={r['rad_isotropic']} zeros={r['zeros']:6d} bias={r['bias']:6d} arf={arf}")
