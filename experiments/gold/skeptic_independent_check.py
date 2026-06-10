# Independent adversarial re-verification (final skeptic round).
# Targets the load-bearing symmetry claims S1/S3/S4 + maximality, and probes
# two suspected soft spots:
#   (a) S1 parenthetical "no nontrivial pure translation" in the DEGENERATE case
#   (b) S3 "orbitals = Gram classes" exactness (diagonal/zero degeneracies)
import itertools

def bits(x, m): return [i for i in range(m) if (x >> i) & 1]

def make_Q(q, B, m):
    def Q(x):
        s = 0
        bs = bits(x, m)
        for i in bs: s ^= q[i]
        for a in range(len(bs)):
            for b in range(a+1, len(bs)):
                s ^= B[bs[a]][bs[b]]
        return s
    return Q

def make_Bv(B, m):
    def Bv(x, y):
        s = 0
        for i in bits(x, m):
            for j in bits(y, m):
                s ^= B[i][j]
        return s
    return Bv

m = 4
V = list(range(1 << m))
# B = standard symplectic: e0~e1, e2~e3
B = [[0]*m for _ in range(m)]
B[0][1] = B[1][0] = 1
B[2][3] = B[3][2] = 1
Bv = make_Bv(B, m)

# Arf 0: q = x0x1 + x2x3 (q_i = 0); Arf 1: q0=q1=1 (x0^2+x0x1+x1^2 + x2x3)
forms = {"Arf0": [0,0,0,0], "Arf1": [1,1,0,0]}

# enumerate GL(4,2) as column matrices encoded as tuples of 4 column vectors
def apply(g, x, m):  # g = tuple of column images of e_i
    y = 0
    for i in bits(x, m): y ^= g[i]
    return y

cols = list(itertools.product(V, repeat=m))
GL = []
for g in cols:
    # invertible iff images of all 16 vectors distinct
    img = set(apply(g, x, m) for x in V)
    if len(img) == 1 << m:
        GL.append(g)
print("|GL(4,2)| =", len(GL), "(expect 20160)")

Sp = [g for g in GL if all(Bv(apply(g,x,m), apply(g,y,m)) == Bv(x,y)
                           for x in V for y in V)]
print("|Sp(4,2)| =", len(Sp), "(expect 720)")

for name, q in forms.items():
    Q = make_Q(q, B, m)
    Z = frozenset(x for x in V if Q(x) == 0)
    O = [g for g in Sp if all(Q(apply(g,x,m)) == Q(x) for x in V)]
    # --- S1: Stab_AGL(Z) = AO(Q); pure translations
    stab = []
    for g in GL:
        for c in V:
            img = frozenset(apply(g,x,m) ^ c for x in Z)
            if img == Z:
                stab.append((g,c))
    AO = []
    for g in Sp:
        for c in V:
            if all(Q(apply(g,x,m) ^ c) == Q(x) for x in V):
                AO.append((g,c))
    ident = tuple(1 << i for i in range(m))
    pure_trans = [c for (g,c) in stab if g == ident and c != 0]
    print(f"{name}: |Z|={len(Z)} |O(Q)|={len(O)} |Stab_AGL(Z)|={len(stab)} "
          f"|AO(Q)|={len(AO)} equal={sorted(stab)==sorted(AO)} "
          f"pure_translations={pure_trans}")
    # --- S3: O(Q)-orbitals on VxV vs Gram classes
    # orbitals
    pairs = [(u,w) for u in V for w in V]
    seen = {}
    orbitals = []
    for p in pairs:
        if p in seen: continue
        orb = set()
        stack = [p]
        while stack:
            (u,w) = stack.pop()
            if (u,w) in orb: continue
            orb.add((u,w))
            for g in O:
                t = (apply(g,u,m), apply(g,w,m))
                if t not in orb: stack.append(t)
        for t in orb: seen[t] = len(orbitals)
        orbitals.append(orb)
    # Gram classes (Q(u), Q(w), B(u,w))
    gram = {}
    for (u,w) in pairs:
        gram.setdefault((Q(u), Q(w), Bv(u,w)), set()).add((u,w))
    # extended Gram: + degeneracy pattern (u==0, w==0, u==w)
    egram = {}
    for (u,w) in pairs:
        egram.setdefault((Q(u), Q(w), Bv(u,w), u==0, w==0, u==w), set()).add((u,w))
    bare_match = sorted(map(sorted, orbitals)) == sorted(map(sorted, gram.values()))
    ext_match  = sorted(map(sorted, orbitals)) == sorted(map(sorted, egram.values()))
    print(f"  S3: #orbitals={len(orbitals)} #GramClasses={len(gram)} "
          f"bare_Gram_match={bare_match} #extGram={len(egram)} ext_match={ext_match}")
    # --- maximality: every g in Sp \ O together with O generates Sp
    Oset = set(O)
    outside = [g for g in Sp if g not in Oset]
    def compose(g, h):  # (g.h)(x) = g(h(x))
        return tuple(apply(g, h[i], m) for i in range(m))
    def gen_closure(gens):
        elems = set(gens)
        frontier = list(elems)
        while frontier:
            nf = []
            for a in frontier:
                for b in list(elems):
                    for c in (compose(a,b), compose(b,a)):
                        if c not in elems:
                            elems.add(c); nf.append(c)
            frontier = nf
            if len(elems) == 720: break
        return elems
    allgen = all(len(gen_closure(O + [g])) == 720 for g in outside)
    print(f"  maximality: |Sp\\O|={len(outside)}, all generate Sp: {allgen}")
    # --- S4 counting: N(v) = #{x: Q(x)=0=Q(x+v)} for v != 0
    Ns = {}
    for v in V[1:]:
        Ns.setdefault(sum(1 for x in V if Q(x)==0 and Q(x^v)==0), []).append(v)
    arf = 0 if len(Z) == 10 else 1
    formula = (1 << (2*2-2)) + ((-1)**arf) * (1 << (2-1))
    print(f"  S4: N(v) value -> #v:",
          {n: len(vlist) for n, vlist in sorted(Ns.items())},
          f"formula_claim={formula} min>0: {min(Ns) > 0}")

# --- (a) degenerate-case probe: m=5, core rank 4, isotropic radical e4
m2 = 5
V2 = list(range(1 << m2))
B2 = [[0]*m2 for _ in range(m2)]
B2[0][1] = B2[1][0] = 1
B2[2][3] = B2[3][2] = 1
q2 = [0]*m2  # isotropic radical: q[4] = 0
Q2 = make_Q(q2, B2, m2)
Z2 = frozenset(x for x in V2 if Q2(x) == 0)
c = 1 << 4  # translation by radical vector e4
img = frozenset(x ^ c for x in Z2)
print("\nDegenerate probe (m=5, rank-4 core, isotropic radical):",
      "translation by e4 stabilizes Z:", img == Z2,
      "| Q(x+e4)==Q(x) for all x:", all(Q2(x ^ c) == Q2(x) for x in V2))
