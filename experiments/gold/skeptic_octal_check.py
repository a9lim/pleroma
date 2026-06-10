"""Independent adversarial verification of the 'anisotropic-normal-frame
unwinding game' claims (octal/coin-turning attack). Written from the attack's
STATED rule, not from its code.

  K1. Key lemma: Q_a(beta^{2^i}) = Q_a(beta) on conjugate orbits.
  K2. Normal-generator counts: m=4: 8 normal, all anisotropic (Q_1=1);
      m=8: 128 normal, 64 anisotropic (a=1); m=16: beta=2^15 normal+anisotropic.
  K3. Unwinding game P-set == {(S,eps): eps = Q(x_S)} for anisotropic normal
      frames: m=4 (all), m=8 a=1/a=2, m=8 bent lambda=2,3, m=16 a=1.
      Zero counts: F16:4, F256 a=1:112, a=2:96, F_2^16: 32512.
  K4. Decision-degeneracy: all options of every position share one outcome.
  K5. Isotropic-frame control: P-set is the parity-mixed set, NOT {eps=Q}.
  K6. Pair-move control: pair moves break exactness AND make it nondegenerate.
  K7. Identity Q_lambda(x^2) = Q_sqrt(lambda)(x) on F_256.
  K8. m=4 unscaled Gold zero set = F_4 = {0,1,2,3} (affine, vacuous target).
"""
import sys
from functools import lru_cache

sys.setrecursionlimit(100000)

@lru_cache(maxsize=None)
def nim_mul(a, b):
    if a < b:
        a, b = b, a
    if b == 0:
        return 0
    if b == 1:
        return a
    F = 2
    while F * F <= a:
        F = F * F
    a1, a0 = divmod(a, F)
    b1, b0 = divmod(b, F)
    c2 = nim_mul(a1, b1)
    c1 = nim_mul(a1, b0) ^ nim_mul(a0, b1)
    c0 = nim_mul(a0, b0)
    return ((c1 ^ c2) * F) ^ c0 ^ nim_mul(c2, F >> 1)

def nim_sq(x): return nim_mul(x, x)

def frob(x, a):
    for _ in range(a):
        x = nim_sq(x)
    return x

def trace(x, m):
    acc, t = 0, x
    for _ in range(m):
        acc ^= t
        t = nim_sq(t)
    return acc

def gold(v, lam, a, m):
    return trace(nim_mul(lam, nim_mul(v, frob(v, a))), m)

def polar(u, v, lam, a, m):
    return gold(u ^ v, lam, a, m) ^ gold(u, lam, a, m) ^ gold(v, lam, a, m)

# pinned validation (repo: coin_turning.rs tests + goldarf.tex Table 2)
assert nim_mul(2, 2) == 3 and nim_mul(2, 3) == 1 and nim_mul(4, 4) == 6
assert nim_mul(2, 4) == 8 and nim_mul(16, 16) == 24
print("[ok] nim products match repo-pinned values")

def conjugates(beta, m):
    out, x = [], beta
    for _ in range(m):
        out.append(x)
        x = nim_sq(x)
    return out

def f2_rank_vecs(vecs):
    rows = [v for v in vecs if v]
    rank = 0
    width = max((v.bit_length() for v in rows), default=0)
    for col in range(width):
        piv = next((i for i in range(rank, len(rows)) if (rows[i] >> col) & 1), None)
        if piv is None:
            continue
        rows[rank], rows[piv] = rows[piv], rows[rank]
        for i in range(len(rows)):
            if i != rank and (rows[i] >> col) & 1:
                rows[i] ^= rows[rank]
        rank += 1
    return rank

def is_normal(beta, m):
    return f2_rank_vecs(conjugates(beta, m)) == m

# ---- K1
bad = 0
for beta in range(1, 256):
    if not is_normal(beta, 8):
        continue
    for a in (1, 2):
        if len({gold(c, 1, a, 8) for c in conjugates(beta, 8)}) != 1:
            bad += 1
print(f"[K1] Q_a constant on conjugate orbits (m=8, a=1,2): violations={bad}")

# ---- K2
for m, en, ea in [(4, 8, 8), (8, 128, 64)]:
    normals = [b for b in range(1, 1 << m) if is_normal(b, m)]
    aniso = [b for b in normals if gold(b, 1, 1, m) == 1]
    print(f"[K2] m={m}: normal={len(normals)} (claim {en}), anisotropic a=1: {len(aniso)} (claim {ea})")
b16 = 1 << 15
print(f"[K2] m=16: beta=2^15 normal? {is_normal(b16, 16)}, Q_1(beta)={gold(b16, 1, 1, 16)}")

# ---- K8
z = [v for v in range(16) if gold(v, 1, 1, 4) == 0]
print(f"[K8] m=4 {{Q_1=0}} = {z} (claim [0,1,2,3] = F_4)")

# ---- unwinding game -------------------------------------------------------
def solve_unwind(m, frame, lam, a, pair_moves=False):
    B = [[polar(frame[i], frame[j], lam, a, m) for j in range(m)] for i in range(m)]
    Brow = [sum(B[i][j] << j for j in range(m)) for i in range(m)]
    N = 1 << m
    xs = [0] * N
    for S in range(1, N):
        low = (S & -S).bit_length() - 1
        xs[S] = xs[S & (S - 1)] ^ frame[low]
    out = [False] * (N << 1)
    for S in sorted(range(N), key=lambda s: bin(s).count("1")):
        for eps in (0, 1):
            opts = []
            for i in range(m):
                if not (S >> i) & 1:
                    continue
                flip = bin(Brow[i] & S).count("1") & 1  # B[i][i]=0
                opts.append(((S ^ (1 << i)) << 1) | (eps ^ flip))
            if pair_moves:
                for i in range(m):
                    if not (S >> i) & 1: continue
                    for j in range(i + 1, m):
                        if not (S >> j) & 1: continue
                        fl = (bin(Brow[i] & S).count("1") + bin(Brow[j] & S).count("1")) & 1
                        opts.append(((S ^ (1 << i) ^ (1 << j)) << 1) | (eps ^ fl))
            out[(S << 1) | eps] = (eps == 0) if not opts else not any(out[o] for o in opts)
    return out, xs, Brow

def check_frame(m, beta, lam, a, label, pair_moves=False):
    frame = conjugates(beta, m)
    if f2_rank_vecs(frame) != m:
        return f"  {label}: beta NOT normal, skip"
    out, xs, Brow = solve_unwind(m, frame, lam, a, pair_moves=pair_moves)
    N = 1 << m
    qt = [gold(xs[S], lam, a, m) for S in range(N)]
    exact = all(out[(S << 1) | e] == (e == qt[S]) for S in range(N) for e in (0, 1))
    pm = all(out[(S << 1) | e] == ((qt[S] ^ e ^ (bin(S).count('1') & 1)) == 0)
             for S in range(N) for e in (0, 1))
    degen = True
    for S in range(1, N):
        for e in (0, 1):
            vals = set()
            for i in range(m):
                if not (S >> i) & 1: continue
                flip = bin(Brow[i] & S).count("1") & 1
                vals.add(out[((S ^ (1 << i)) << 1) | (e ^ flip)])
            if pair_moves:
                for i in range(m):
                    if not (S >> i) & 1: continue
                    for j in range(i + 1, m):
                        if not (S >> j) & 1: continue
                        fl = (bin(Brow[i] & S).count("1") + bin(Brow[j] & S).count("1")) & 1
                        vals.add(out[((S ^ (1 << i) ^ (1 << j)) << 1) | (e ^ fl)])
            if len(vals) > 1:
                degen = False
    zc = sum(1 for S in range(N) if qt[S] == 0)
    return (f"  {label}: Q(beta)={gold(beta, lam, a, m)} exact={exact} "
            f"parity-mixed={pm} degenerate={degen} |{{Q=0}}|={zc}")

print("\n[K3/K4] anisotropic normal frames:")
for beta in [b for b in range(1, 16) if is_normal(b, 4)]:
    print(check_frame(4, beta, 1, 1, f"m=4 a=1 beta={beta}"))

done = 0
for beta in range(1, 256):
    if done >= 3: break
    if is_normal(beta, 8) and gold(beta, 1, 1, 8) == 1:
        print(check_frame(8, beta, 1, 1, f"m=8 a=1 beta={beta}")); done += 1
done = 0
for beta in range(1, 256):
    if done >= 3: break
    if is_normal(beta, 8) and gold(beta, 1, 2, 8) == 1:
        print(check_frame(8, beta, 1, 2, f"m=8 a=2 beta={beta}")); done += 1
for lam in (2, 3):
    full = [b for b in range(1, 256) if is_normal(b, 8)
            and all(gold(c, lam, 1, 8) == 1 for c in conjugates(b, 8))]
    print(f"  m=8 bent lam={lam}: fully-anisotropic normal generators: {len(full)}/128 (claim 32)")
    if full:
        print(check_frame(8, full[0], lam, 1, f"m=8 bent lam={lam} beta={full[0]}"))

print("\n[K5] isotropic-frame control (m=8 a=1, Q(beta)=0):")
done = 0
for beta in range(1, 256):
    if done >= 2: break
    if is_normal(beta, 8) and gold(beta, 1, 1, 8) == 0:
        print(check_frame(8, beta, 1, 1, f"m=8 a=1 ISO beta={beta}")); done += 1

print("\n[K6] pair-move control (m=8 a=1, anisotropic frame + pair moves):")
done = 0
for beta in range(1, 256):
    if done >= 1: break
    if is_normal(beta, 8) and gold(beta, 1, 1, 8) == 1:
        print(check_frame(8, beta, 1, 1, f"m=8 a=1 +pairs beta={beta}", pair_moves=True)); done += 1

print("\n[K7] Q_lambda(x^2) == Q_sqrt(lambda)(x) on F_256:")
viol = 0
for lam in range(1, 256):
    s = frob(lam, 7)
    assert nim_sq(s) == lam
    for x in range(256):
        if gold(nim_sq(x), lam, 1, 8) != gold(x, s, 1, 8):
            viol += 1
print(f"  violations: {viol}")

print("\n[K3] m=16 a=1 (the big one):")
print(check_frame(16, 1 << 15, 1, 1, "m=16 a=1 beta=2^15"))
print("done")
