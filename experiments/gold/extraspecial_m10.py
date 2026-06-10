"""Scale check: synthetic rank-4 Arf-0/1 forms on F_2^10, adapted frame,
full 1024-position sweep."""
import time
from extraspecial_core import echo_value

def coord_form(m, hyp_pairs, q1pairs=0):
    Q = []
    for x in range(1 << m):
        t = 0
        for p in range(hyp_pairs):
            a, b = (x >> (2 * p)) & 1, (x >> (2 * p + 1)) & 1
            t ^= a & b
            if p < q1pairs: t ^= a ^ b
        Q.append(t)
    return Q

m = 10
for name, Q in (("r4 Arf0 rad6iso m=10", coord_form(10, 2, 0)),
                ("r4 Arf1 rad6iso m=10", coord_form(10, 2, 1))):
    qover = [Q[1 << i] for i in range(m)]
    Bover = []
    for i in range(m):
        row = 0
        for j in range(m):
            if i == j: continue
            row |= (Q[(1 << i) ^ (1 << j)] ^ Q[1 << i] ^ Q[1 << j]) << j
        Bover.append(row)
    t0 = time.time()
    miss = 0; ndg = 0
    for x in range(1 << m):
        v, ch = echo_value(x, None, None, m, ko='self', maxfirst=True,
                           qover=qover, Bover=Bover)
        ndg += ch
        if v != Q[x]: miss += 1
    print(f"{name}: miss={miss}/1024 choice-states={ndg} ({time.time()-t0:.0f}s)")
