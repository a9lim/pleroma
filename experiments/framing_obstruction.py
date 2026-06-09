"""The framing obstruction probe: why the Gold quadric needs more than B alone.

This probe stress-tests the structure of OPEN.md's open question — "a game whose
moves are built from B (coin-turning) alone, not from Q itself, with P-set the
Gold quadric {Q=0}" — by organizing the obstruction into a symmetry-breaking
ladder, using two classical facts:

  (i)  the quadratic refinements of a fixed symplectic form B form a torsor under the
       linear duals V* (two forms share a polar form iff they differ by a linear
       functional — arXiv:2506.23613); and
  (ii) Sp(2m,2) acts on those refinements with exactly two orbits, distinguished by
       the Arf invariant, the stabilizer of one form being its orthogonal group O(Q)
       (the F2 classification, e.g. arXiv:1001.4751).

THE LADDER (each rung adds exactly one datum, breaking one layer of symmetry):

  sees                       symmetry     reachable P-set
  --------------------------------------------------------------------------------
  abstract B only            Sp(B)        NO quadric  (Sp(B) is transitive on V\\{0},
                                          so the only invariant sets are 0/{0}/V\\0/V)
  B + coordinate frame       O(Q_frame)   the FRAME quadric {Q_frame=0}, via an
                                          explicit oracle rule using Q_frame
  B + frame + diagonal q_i   O(Q_gold)    the GOLD quadric {Q_gold=0}, any Arf

The diagonal q_i = Q(e_i) = Tr(e_i^{1+2^a}) is the single-coin self-Gold value — a
FRAMING in the Arf-Kervaire sense. And empirically (Part 3) it always carries Arf 1:
Q_frame is split (O+) in the standard basis, the framing flips it to the Gold O-.

So the probe supports this structural reading:
  * frame-blind (Sp(B)-equivariant) games provably cannot realize ANY quadric;
  * a deliberately framed B+frame rule realizes the split frame quadric, not by
    natural play semantics but by directly reading Q_frame;
  * the gap to the Gold quadric is exactly the diagonal framing, m bits not derivable
    from B in this model.
"""

import itertools
from math import gcd

import pleroma as pl

# ----------------------------------------------------------------------------- helpers


def gold(v, a, m):
    """Gold form Q_a(v) = Tr(v^{1+2^a}) over F_2^m, valued in {0,1}."""
    x = pl.Nimber(v)
    g = x
    for _ in range(a):
        g = g * g
    s = x * g
    acc, t = s, s
    for _ in range(m - 1):
        t = t * t
        acc = acc + t
    return acc.value


def gold_pairs(a, m):
    """The nonzero pairs (i<j) of the Gold polar form B on the bit basis."""
    return [(i, j) for i in range(m) for j in range(i + 1, m)
            if (gold((1 << i) ^ (1 << j), a, m) ^ gold(1 << i, a, m) ^ gold(1 << j, a, m)) == 1]


def Bform(pairs):
    def B(u, v):
        acc = 0
        for (i, j) in pairs:
            if (u >> i) & 1 and (v >> j) & 1:
                acc ^= 1
            if (u >> j) & 1 and (v >> i) & 1:
                acc ^= 1
        return acc
    return B


def arf(qd, pairs, m):
    q = [pl.Nimber(1 if qd[i] else 0) for i in range(m)]
    b = {(i, j): pl.Nimber(1) for (i, j) in pairs}
    return pl.arf_invariant(pl.NimberAlgebra(q=q, b=b))


def kernel(succ):
    """P-positions (Loss) of a finite DAG game (normal play). Loss = the kernel."""
    label = [""] * len(succ)  # "" = unvisited, else "Win"/"Loss"

    def solve(v):
        if label[v]:
            return label[v]
        label[v] = "Win"  # guard against revisits (graph is a DAG: w<v)
        res = "Loss"
        for w in succ[v]:
            if solve(w) == "Loss":
                res = "Win"
                break
        label[v] = res
        return res

    for v in range(len(succ)):
        solve(v)
    return set(v for v in range(len(succ)) if label[v] == "Loss")


def invertible(cols, m):
    rows = list(cols)
    r = 0
    for bit in range(m):
        piv = next((k for k in range(r, len(rows)) if (rows[k] >> bit) & 1), None)
        if piv is None:
            continue
        rows[r], rows[piv] = rows[piv], rows[r]
        for k in range(len(rows)):
            if k != r and (rows[k] >> bit) & 1:
                rows[k] ^= rows[r]
        r += 1
    return r == m


def apply_map(cols, v, m):
    out = 0
    for i in range(m):
        if (v >> i) & 1:
            out ^= cols[i]
    return out


# ----------------------------------------------------------------------------- part 1


def part1_nogo():
    print("=" * 72)
    print("PART 1 — the no-go: a frame-blind (Sp(B)-equivariant) game has no quadric")
    print("=" * 72)
    m = 4
    N = 1 << m
    pairs = [(0, 1), (2, 3)]  # nondegenerate symplectic form on F_2^4
    B = Bform(pairs)
    basis = [1 << i for i in range(m)]
    sp = [cols for cols in itertools.product(range(N), repeat=m)
          if invertible(cols, m)
          and all(B(cols[i], cols[j]) == B(basis[i], basis[j])
                  for i in range(m) for j in range(i + 1, m))]
    print(f"|Sp(B)| = {len(sp)}   (= |Sp(4,2)| = 720)")

    # Sp(B) is transitive on nonzero vectors ⇒ only invariant sets are ∅,{0},V\0,V.
    orbit_of_e0 = sorted(set(apply_map(g, 1, m) for g in sp))
    print(f"orbit of e_0 under Sp(B): {len(orbit_of_e0)} vectors "
          f"(= all {N-1} nonzero ⇒ transitive)")
    print("⇒ the ONLY Sp(B)-invariant subsets of V are ∅, {0}, V\\{0}, V — no quadric.\n")

    for name, qd in [("split  Q (Arf 0)", [0, 0, 0, 0]),
                     ("Gold-type Q (Arf 1)", [1, 1, 0, 0])]:
        def Q(v, qd=qd):
            acc = sum(1 for i in range(m) if qd[i] and (v >> i) & 1)
            for (i, j) in pairs:
                if (v >> i) & 1 and (v >> j) & 1:
                    acc ^= 1
            return acc & 1
        zeros = set(v for v in range(N) if Q(v) == 0)
        ar = arf(qd, pairs, m)
        stab = [g for g in sp if set(apply_map(g, v, m) for v in zeros) == zeros]
        print(f"  {name}: |Q=0|={len(zeros):>2}, |O(Q)|={len(stab):>3}, "
              f"orbit={len(sp)//len(stab):>2} (= # refinements with Arf {ar.arf}); "
              f"{{Q=0}} moved by Sp(B): {len(stab) < len(sp)}")
    print()


# ----------------------------------------------------------------------------- part 2


def part2_frame_quadric():
    print("=" * 72)
    print("PART 2 — an explicit B+frame rule realizes the FRAME quadric {Q_frame=0}")
    print("=" * 72)
    print("Rule: move v -> any w<v with Q_frame(w) != Q_frame(v). Q_frame(v) =")
    print("Σ_{i<j∈v} B(e_i,e_j). This is a framed oracle rule, not natural semantics.\n")
    import random
    rng = random.Random(0xF00D)
    for m in (4, 6, 8):
        allp = [(i, j) for i in range(m) for j in range(i + 1, m)]
        ok = tot = genuine = 0
        for _ in range(40):
            pairs = [p for p in allp if rng.random() < 0.4]
            if not pairs:
                continue

            def Qf(v, pairs=pairs):
                acc = 0
                for (i, j) in pairs:
                    if (v >> i) & 1 and (v >> j) & 1:
                        acc ^= 1
                return acc
            zeros = set(v for v in range(1 << m) if Qf(v) == 0)
            succ = [[w for w in range(v) if Qf(w) != Qf(v)] for v in range(1 << m)]
            ker = kernel(succ)
            tot += 1
            ok += (ker == zeros)
            if arf([0] * m, pairs, m).rank >= 2:
                genuine += 1
        print(f"  m={m}: kernel == {{Q_frame=0}} in {ok}/{tot} random B/frame samples "
              f"({genuine} genuinely quadratic)")
    print("  ⇒ the explicit framed rule realizes the frame quadric by construction.\n")


# ----------------------------------------------------------------------------- part 3


def part3_framing_flips_arf():
    print("=" * 72)
    print("PART 3 — the framing: Q_gold = Q_frame ⊕ ℓ_diag, and ℓ_diag flips Arf 0→1")
    print("=" * 72)
    print(f"{'m':>3} {'a':>2} {'rank':>4} {'gcd':>3} | {'Arf(gold)':>9} {'Arf(frame)':>10} "
          f"{'decomp ok':>9}")
    print("-" * 50)
    for k in (1, 2, 3, 4):
        m = 1 << k
        for a in range(1, k + 1):
            pairs = gold_pairs(a, m)
            diag = [gold(1 << i, a, m) for i in range(m)]

            def Qframe(v, pairs=pairs):
                acc = 0
                for (i, j) in pairs:
                    if (v >> i) & 1 and (v >> j) & 1:
                        acc ^= 1
                return acc

            def ell(v, diag=diag):
                return sum(diag[i] for i in range(m) if (v >> i) & 1) & 1

            decomp = all(gold(v, a, m) == (Qframe(v) ^ ell(v)) for v in range(1 << m))
            ag = arf([d == 1 for d in diag], pairs, m)
            af = arf([0] * m, pairs, m)
            print(f"{m:>3} {a:>2} {ag.rank:>4} {gcd(a, m):>3} | {ag.arf:>9} {af.arf:>10} "
                  f"{str(decomp):>9}")
    print("\n  Reading for these Gold polar forms: whenever the form is genuinely")
    print("  quadratic (rank≥2), the frame quadric is split (Arf 0) and the diagonal")
    print("  framing flips it to the Gold O- (Arf 1).")


if __name__ == "__main__":
    part1_nogo()
    part2_frame_quadric()
    part3_framing_flips_arf()
