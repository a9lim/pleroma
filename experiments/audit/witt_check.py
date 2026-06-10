# Simulate WittVec<2,2,2> = (Z/4)[t]/(t^2 - (1+t)), i.e. W_2(F_4) per the code.
M = 4  # p^N = 2^2
# F_4 = F_2[t]/(t^2+t+1); elements as (c0,c1) mod 2

def ring_mul(a, b):
    # (a0+a1 t)(b0+b1 t) mod (t^2 = 1+t), coeffs mod 4
    s0 = a[0]*b[0]
    s1 = a[0]*b[1] + a[1]*b[0]
    s2 = a[1]*b[1]
    # t^2 -> 1 + t
    return ((s0 + s2) % M, (s1 + s2) % M)

def ring_add(a, b):
    return ((a[0]+b[0]) % M, (a[1]+b[1]) % M)

def ring_pow(a, e):
    acc = (1, 0)
    base = a
    while e:
        if e & 1:
            acc = ring_mul(acc, base)
        base = ring_mul(base, base)
        e >>= 1
    return acc

def residue(a):
    return (a[0] % 2, a[1] % 2)

def teich(x):
    # naive lift then ^(q^{N-1}) = ^4 once (N-1 = 1 iteration of pow q=4)
    y = (x[0], x[1])
    y = ring_pow(y, 4)
    return y

def divide_by_2(a):
    assert a[0] % 2 == 0 and a[1] % 2 == 0, a
    return (a[0]//2, a[1]//2)

def sub(a, b):
    return ((a[0]-b[0]) % M, (a[1]-b[1]) % M)

def witt_components(a):
    out = []
    for _ in range(2):
        r = residue(a)
        out.append(r)
        t = teich(r)
        a = divide_by_2(sub(a, t))
    return out

def from_witt_components(xs):
    acc = (0,0)
    pk = (1,0)
    for x in xs:
        acc = ring_add(acc, ring_mul(teich(x), pk))
        pk = ring_mul(pk, (2,0))
    return acc

# F_4 arithmetic for the comparison
def f4_mul(x, y):
    s0 = x[0]*y[0]; s1 = x[0]*y[1]+x[1]*y[0]; s2 = x[1]*y[1]
    return ((s0+s2) % 2, (s1+s2) % 2)

def f4_add(x, y):
    return ((x[0]+y[0]) % 2, (x[1]+y[1]) % 2)

def f4_sqrt(x):
    # Frobenius inverse: sqrt = x^2 in F_4
    return f4_mul(x, x)

F4 = [(0,0),(1,0),(0,1),(1,1)]

# sanity: roundtrip
for c0 in range(4):
    for c1 in range(4):
        w = (c0, c1)
        assert from_witt_components(witt_components(w)) == w, w
print("roundtrip ok")

classical_ok = True
twisted_ok = True
bad_classical = []
for x0 in F4:
    for x1 in F4:
        for y0 in F4:
            for y1 in F4:
                a = from_witt_components([x0, x1])
                b = from_witt_components([y0, y1])
                z = witt_components(ring_add(a, b))
                z0, z1 = z
                # classical Witt addition (p=2, char 2): S0 = x0+y0, S1 = x1+y1+x0*y0
                c1 = f4_add(f4_add(x1, y1), f4_mul(x0, y0))
                # twisted (Teichmuller-digit) law: carry = sqrt(x0*y0)
                t1 = f4_add(f4_add(x1, y1), f4_sqrt(f4_mul(x0, y0)))
                if z0 != f4_add(x0, y0):
                    print("S0 mismatch!", x0, y0, z0)
                if z1 != c1:
                    classical_ok = False
                    if len(bad_classical) < 3:
                        bad_classical.append((x0,x1,y0,y1,z1,c1))
                if z1 != t1:
                    twisted_ok = False
print("classical Witt law z1 = x1+y1+x0*y0 holds over F_4:", classical_ok)
print("twisted law z1 = x1+y1+sqrt(x0*y0) holds over F_4:", twisted_ok)
for r in bad_classical:
    print("counterexample x0,x1,y0,y1 -> z1, classical:", r)
