# Over F_256: count lambda != 0 with Q_lambda(x) = Tr(lambda x^{2^a+1}) bent (rank 8), for gcd(a,8)=1.
m = 8; mod = 0b100011011
def mul(a,b):
    r=0
    while b:
        if b&1: r^=a
        b>>=1; a<<=1
        if a>>m&1: a^=mod
    return r
def sq(a): return mul(a,a)
def tr(a):
    t,c=0,a
    for _ in range(m): t^=c; c=sq(c)
    return t
def rank_of(lam, a):
    def Q(x):
        y=x
        for _ in range(a): y=sq(y)
        return tr(mul(lam, mul(x,y)))
    qv=[Q(1<<i) for i in range(m)]
    B=[[Q((1<<i)^(1<<j))^qv[i]^qv[j] for j in range(m)] for i in range(m)]
    rows=[int(''.join(map(str,r)),2) for r in B]; rank=0
    for col in range(m):
        piv=None
        for r in range(rank,m):
            if rows[r]>>(m-1-col)&1: piv=r;break
        if piv is None: continue
        rows[rank],rows[piv]=rows[piv],rows[rank]
        for r in range(m):
            if r!=rank and rows[r]>>(m-1-col)&1: rows[r]^=rows[rank]
        rank+=1
    return rank
for a in (1,3):
    bent=sum(1 for lam in range(1,256) if rank_of(lam,a)==m)
    print(f"a={a}: bent components {bent}/255 (classical 2(2^8-1)/3 = {2*255//3})")
