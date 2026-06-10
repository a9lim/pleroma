"""Verify the attack's CORRECTED sweep numbers absent from the cited script:
   (a) m=4 exhaustive misere: 32672 affine, 96 deg>2, 0 non-affine quadrics
   (b) random m=8 (seed 0, same generator shape): 183 affine, 217 deg>2,
       0 non-affine quadrics."""
import random

def solve_coin(n, comps, misere):
    N = 1 << n
    out = [False] * N
    for v in range(N):
        has_move = to_p = False
        for r in range(n):
            if not (v >> r) & 1: continue
            for S in comps[r]:
                has_move = True
                if out[v ^ (1 << r) ^ S]:
                    to_p = True; break
            if to_p: break
        out[v] = (not misere) if not has_move else (not to_p)
    return [v for v in range(N) if out[v]]

def is_affine(pts):
    if not pts: return True
    s0 = pts[0]
    sh = {p ^ s0 for p in pts}
    return all((x ^ y) in sh for x in sh for y in sh)

def anf_deg(pset, k):
    N = 1 << k
    c = [1] * N
    for v in pset: c[v] = 0
    for i in range(k):
        bit = 1 << i
        for mask in range(N):
            if mask & bit: c[mask] ^= c[mask ^ bit]
    return max((bin(m).count("1") for m in range(N) if c[m]), default=0)

# (a) exhaustive m=4
per_coin = []
for r in range(4):
    masks = list(range(1 << r))
    per_coin.append([[masks[i] for i in range(len(masks)) if (b >> i) & 1]
                     for b in range(1 << len(masks))])
aff = d2 = quad_nonaff = 0
for f0 in per_coin[0]:
    for f1 in per_coin[1]:
        for f2 in per_coin[2]:
            for f3 in per_coin[3]:
                p = solve_coin(4, [f0, f1, f2, f3], True)
                a = is_affine(p)
                d = anf_deg(p, 4)
                if a: aff += 1
                elif d > 2: d2 += 1
                else: quad_nonaff += 1
print(f"(a) m=4 exhaustive misere: affine={aff} (claim 32672), deg>2 nonaffine={d2} (claim 96), nonaffine quadric={quad_nonaff} (claim 0)")

# (b) random m=8 with the SAME seed/generator shape as the cited script
random.seed(0)
aff = d2 = quad_nonaff = 0
for trial in range(400):
    comps = []
    for r in range(8):
        kfam = random.randint(1, 4)
        comps.append([random.randrange(1 << r) for _ in range(kfam)])
    p = solve_coin(8, comps, True)
    a = is_affine(p)
    d = anf_deg(p, 8)
    if a: aff += 1
    elif d > 2: d2 += 1
    else: quad_nonaff += 1
print(f"(b) random m=8: affine={aff} (claim 183), deg>2 nonaffine={d2} (claim 217), nonaffine quadric={quad_nonaff} (claim 0)")
