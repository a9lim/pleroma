import random
from math import gcd

def div_euclid(a, b):
    q, r = divmod(a, b)
    return q + 1 if r < 0 else q

# ---- faithful port of normalize_relation_rows (integer.rs) ----
def leading(row):
    for i, x in enumerate(row):
        if x != 0: return i
    return None

def normalize_relation_rows(rows):
    rows = [r[:] for r in rows if any(x != 0 for x in r)]
    width = len(rows[0]) if rows else 0
    rank = 0
    for col in range(width):
        piv = next((r for r in range(rank, len(rows)) if rows[r][col] != 0), None)
        if piv is None: continue
        rows[rank], rows[piv] = rows[piv], rows[rank]
        if rows[rank][col] < 0:
            rows[rank] = [-x for x in rows[rank]]
        while True:
            r = next((r for r in range(rank+1, len(rows)) if rows[r][col] != 0), None)
            if r is None: break
            pv = rows[rank][col]
            q = div_euclid(rows[r][col], pv)
            src = rows[rank][:]
            rows[r] = [t - q*s for t, s in zip(rows[r], src)]
            if rows[r][col] != 0 and abs(rows[r][col]) < abs(rows[rank][col]):
                rows[rank], rows[r] = rows[r], rows[rank]
                if rows[rank][col] < 0:
                    rows[rank] = [-x for x in rows[rank]]
        if rows[rank][col] < 0:
            rows[rank] = [-x for x in rows[rank]]
        pv = rows[rank][col]
        src = rows[rank][:]
        for r in range(len(rows)):
            if r == rank or rows[r][col] == 0: continue
            q = div_euclid(rows[r][col], pv)
            rows[r] = [t - q*s for t, s in zip(rows[r], src)]
        rank += 1
    rows = [r for r in rows if any(x != 0 for x in r)]
    rows.sort(key=lambda r: leading(r) if leading(r) is not None else 10**9)
    return rows

def reduce_integer_vector(v, rows):
    v = v[:]
    for row in normalize_relation_rows(rows):
        l = leading(row)
        if l is None: continue
        p = row[l]
        q = div_euclid(v[l], p)
        if q != 0:
            v = [a - q*b for a, b in zip(v, row)]
    return v

# ---- independent incremental HNF (structurally different algorithm) ----
def ext_gcd(a, b):
    old_r, r = a, b
    old_s, s = 1, 0
    old_t, t = 0, 1
    while r:
        q = old_r // r
        old_r, r = r, old_r - q*r
        old_s, s = s, old_s - q*s
        old_t, t = t, old_t - q*t
    if old_r < 0:
        return -old_r, -old_s, -old_t
    return old_r, old_s, old_t

def incremental_hnf(rows, width):
    basis = {}  # leading col -> row
    def insert(v):
        v = v[:]
        while True:
            l = leading(v)
            if l is None: return
            if l not in basis:
                if v[l] < 0: v = [-x for x in v]
                basis[l] = v
                return
            p = basis[l]
            if v[l] % p[l] == 0:
                q = v[l] // p[l]
                v = [a - q*b for a, b in zip(v, p)]
            else:
                g, x, y = ext_gcd(p[l], v[l])
                new = [x*a + y*b for a, b in zip(p, v)]
                # v' has zero at l: (p[l]/g)*v - (v[l]/g)*p
                v = [(p[l]//g)*a - (v[l]//g)*b for a, b in zip(v, p)]
                basis[l] = new
                # re-insert the replaced pivot row? new has lead l (value g>0), fine.
        # unreachable
    for r in rows:
        insert(r)
    # canonicalize: reduce above-pivot entries, bottom-up
    cols = sorted(basis)
    out = [basis[c] for c in cols]
    for i in range(len(out)-1, -1, -1):
        for j in range(i+1, len(out)):
            cj = cols[j]
            pj = out[j][cj]
            q = div_euclid(out[i][cj], pj)
            if q:
                out[i] = [a - q*b for a, b in zip(out[i], out[j])]
    return out

def shape_check(out):
    leads = [leading(r) for r in out]
    assert all(l is not None for l in leads)
    assert leads == sorted(leads) and len(set(leads)) == len(leads)
    for i, r in enumerate(out):
        p = r[leads[i]]
        assert p > 0
        for k in range(len(out)):
            if k != i:
                # entry of row k at pivot col i: below must be 0, above in [0,p)
                e = out[k][leads[i]]
                if k > i:
                    assert e == 0, (out, k, i)
                else:
                    assert 0 <= e < p, (out, k, i)

random.seed(31337)
for trial in range(3000):
    nrows = random.randint(1, 5)
    width = random.randint(1, 5)
    rows = [[random.randint(-12, 12) for _ in range(width)] for _ in range(nrows)]
    out = normalize_relation_rows(rows)
    if out:
        shape_check(out)
    H = incremental_hnf(rows, width)
    assert out == H, (rows, out, H)
    # reduce_integer_vector: lattice elements reduce to zero
    coeffs = [random.randint(-3, 3) for _ in rows]
    v = [sum(c*r[i] for c, r in zip(coeffs, rows)) for i in range(width)]
    assert reduce_integer_vector(v, rows) == [0]*width, (rows, v)
    # canonicity: congruent vectors reduce identically
    w = [random.randint(-20, 20) for _ in range(width)]
    w2 = [a + b for a, b in zip(w, v)]
    assert reduce_integer_vector(w, rows) == reduce_integer_vector(w2, rows), (rows, w, v)
print("HNF: 3000 random cases OK (shape, canonical-form match vs independent impl, reduce canonicity)")
