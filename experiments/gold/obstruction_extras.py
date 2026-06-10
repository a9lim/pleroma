"""Verify: (a) {h in E : h.I = I} = Z(E) for nondegenerate Q (kill arm);
(b) Frobenius in O(Q_a) for unscaled Gold (under-constrain arm);
(c) Gold radicals are isotropic: Q_a|R(B) = 0 for (8,1),(8,2),(4,1)."""
from extraspecial_core import nim_mul, gold_q, polar, validate

validate()

# (a) on the bent (8,1,lam=2) form: E = V x F2 with triangular cocycle
Q = gold_q(8, 1, 2)   # nondegenerate rank 8
m = 8
q = [Q[1 << i] for i in range(m)]
B = polar(Q, m)
def coc(u, v):
    s = 0
    for i in range(m):
        if (u >> i) & 1 and (v >> i) & 1: s ^= q[i]
    for k in range(m):
        for j in range(k):
            if (u >> k) & 1 and (v >> j) & 1 and ((B[k] >> j) & 1): s ^= 1
    return s
def emul(g, h): return (g[0] ^ h[0] ^ coc(g[1], h[1]), g[1] ^ h[1])
E = [(a, u) for a in (0, 1) for u in range(256)]
I = {g for g in E if emul(g, g) == (0, 0)}
print(f"|I| = {len(I)} (= 2*|Z| = {2 * sum(1 for v in Q if v == 0)})")
stab = [h for h in E if {emul(h, g) for g in I} == I]
print(f"left-translation stabilizer of I: {stab}  (predict [(0,0),(1,0)] = Z(E))")

# (b) Frobenius in O(Q_a) for unscaled Gold (8,1), (8,2)
for a in (1, 2):
    Qa = gold_q(8, a, 1)
    ok = all(Qa[nim_mul(x, x)] == Qa[x] for x in range(256))
    print(f"Frobenius preserves Q_{a} (m=8, lam=1): {ok}")

# (c) Gold radical isotropy
for (mm, a) in ((4, 1), (8, 1), (8, 2)):
    Qa = gold_q(mm, a, 1)
    Ba = polar(Qa, mm)
    rad = [v for v in range(1 << mm)
           if all((Qa[v ^ (1 << i)] ^ Qa[v] ^ Qa[1 << i]) == 0 for i in range(mm))]
    print(f"({mm},{a}): |R(B)|={len(rad)}, Q|rad zero: "
          f"{all(Qa[v] == 0 for v in rad)}")
