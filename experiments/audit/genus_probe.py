# Is diag(1,20) ~ diag(5,4) over Z/64 (hence Z_2 by lifting)?
M = 64
import itertools
A = [[1,0],[0,20]]
B = [[5,0],[0,4]]
found = None
# columns u, v of U: need u^T A u = 5, v^T A v = 4, u^T A v = 0 (mod 64), det U odd
us = [(x,y) for x in range(M) for y in range(M) if (x*x + 20*y*y) % M == 5]
vs = [(x,y) for x in range(M) for y in range(M) if (x*x + 20*y*y) % M == 4]
print("candidates:", len(us), len(vs))
for u in us:
    for v in vs:
        if (u[0]*v[0] + 20*u[1]*v[1]) % M == 0 and (u[0]*v[1] - u[1]*v[0]) % 2 == 1:
            found = (u, v)
            break
    if found: break
print("U columns (mod 64):", found)

# Cross-validate the oracle on adjacent-scale forms diag(a,2b): classes under Z/64-equiv
def equiv(A, B, M=64):
    a11,a22 = A; b11,b22 = B
    us = [(x,y) for x in range(M) for y in range(M) if (a11*x*x + a22*y*y) % M == b11 % M]
    vs = [(x,y) for x in range(M) for y in range(M) if (a11*x*x + a22*y*y) % M == b22 % M]
    for u in us:
        for v in vs:
            if (a11*u[0]*v[0] + a22*u[1]*v[1]) % M == 0 and (u[0]*v[1]-u[1]*v[0]) % 2 == 1:
                return True
    return False

# all diag(a, 2b), a,b odd in {1,3,5,7}: which pairs are Z_2-equivalent?
forms = [(a, 2*b) for a in (1,3,5,7) for b in (1,3,5,7)]
classes = []
for f in forms:
    placed = False
    for c in classes:
        if equiv(f, c[0]):
            c.append(f); placed = True; break
    if not placed:
        classes.append([f])
print("adjacent-scale diag(a,2b) Z/64 classes:")
for c in classes: print("  ", c)

# gap forms diag(a, 4b)
forms4 = [(a, 4*b) for a in (1,3,5,7) for b in (1,3,5,7)]
classes4 = []
for f in forms4:
    placed = False
    for c in classes4:
        if equiv(f, c[0]):
            c.append(f); placed = True; break
    if not placed:
        classes4.append([f])
print("gap-scale diag(a,4b) Z/64 classes:")
for c in classes4: print("  ", c)
