"""Normal-frame sweep at m=8 for the four Gold forms, ko=self, both orients.
Also: confirm the (8,2) bit-frame miss x=224 is a bad k=3 pattern."""
import time
from extraspecial_core import *

validate()

def is_normal(beta, m):
    """beta normal <=> {beta^(2^i)} linearly independent over F_2"""
    vecs = []
    b = beta
    for _ in range(m):
        vecs.append(b)
        b = nim_mul(b, b)
    # gaussian elim
    basis = []
    for v in vecs:
        for w in basis:
            v = min(v, v ^ w)
        if v == 0: return False, None
        basis.append(v)
        basis.sort(reverse=True)
    return True, vecs

m = 8
forms = [
    ("(8,1)l1", gold_q(8, 1, 1)),
    ("(8,2)l1", gold_q(8, 2, 1)),
    ("(8,1)l2", gold_q(8, 1, 2)),
    ("(8,1)l3", gold_q(8, 1, 3)),
]

# confirm bit-frame miss pattern for (8,2), support {5,6,7}
Q82 = forms[1][1]; B82 = polar(Q82, 8)
S = [5, 6, 7]
qS = [Q82[1 << c] for c in S]
BS = [(a, b, (B82[S[a]] >> S[b]) & 1) for a in range(3) for b in range(a + 1, 3)]
x = sum(1 << c for c in S)
print(f"(8,2) miss support {S}: q|S={qS} B|S={BS} Q(x)={Q82[x]}")
print(f"  -> bad-pattern test: Q=0 and B-edges>=2: "
      f"{Q82[x] == 0 and sum(e[2] for e in BS) >= 2}")

# enumerate normal elements
normals = []
for beta in range(1, 256):
    ok, vecs = is_normal(beta, m)
    if ok: normals.append((beta, vecs))
print(f"\nnormal elements at m=8: {len(normals)}")

def frame_sweep(Q, vecs, maxfirst):
    """positions = coordinate vectors over frame 'vecs'; value vs Q(field elt)."""
    mm = len(vecs)
    qover = [Q[v] for v in vecs]
    Bover = []
    for i in range(mm):
        row = 0
        for j in range(mm):
            if i == j: continue
            b = Q[vecs[i] ^ vecs[j]] ^ Q[vecs[i]] ^ Q[vecs[j]]
            row |= b << j
        Bover.append(row)
    miss = 0
    for cm in range(1 << mm):
        xf = 0
        for i in range(mm):
            if (cm >> i) & 1: xf ^= vecs[i]
        v, _ = echo_value(cm, None, None, mm, ko='self', maxfirst=maxfirst,
                          qover=qover, Bover=Bover)
        if v != Q[xf]: miss += 1
    return miss

# cheap k<=3 pre-screen per frame: count bad triples
def bad_triples(Q, vecs):
    mm = len(vecs)
    bad = 0
    import itertools
    for tri in itertools.combinations(range(mm), 3):
        x = vecs[tri[0]] ^ vecs[tri[1]] ^ vecs[tri[2]]
        if Q[x] != 0: continue
        edges = 0
        for a in range(3):
            for b in range(a + 1, 3):
                i, j = tri[a], tri[b]
                edges += Q[vecs[i] ^ vecs[j]] ^ Q[vecs[i]] ^ Q[vecs[j]]
        if edges >= 2: bad += 1
    return bad

for name, Q in forms:
    # prescreen all normal frames by k=3 cleanliness
    clean = []
    bt_hist = {}
    for beta, vecs in normals:
        bt = bad_triples(Q, vecs)
        bt_hist[bt] = bt_hist.get(bt, 0) + 1
        if bt == 0: clean.append((beta, vecs))
    print(f"\n{name}: bad-triple histogram over {len(normals)} normal frames: "
          f"{dict(sorted(bt_hist.items()))}")
    print(f"  k=3-clean normal frames: {len(clean)}")
    # full sweep on up to 4 cleanest frames (or least-bad if none clean)
    cands = clean[:4]
    if not cands:
        best = sorted(normals, key=lambda bv: bad_triples(Q, bv[1]))[:2]
        cands = best
    for beta, vecs in cands:
        for mf in (True, False):
            t0 = time.time()
            miss = frame_sweep(Q, vecs, mf)
            print(f"  beta={beta:3d} {'P1max' if mf else 'P1min'}: "
                  f"miss={miss}/256 ({time.time()-t0:.0f}s)")
