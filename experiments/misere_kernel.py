"""The misère kernel no-go: why the misère route cannot give a genuine quadric.

The misère route looked promising because misère sums are NOT XOR-linear, so the
P-set escapes the normal-play subspace obstruction. But a quadric {Q=0} needs an
ambient F_2-VECTOR SPACE to live in, and the Plambeck–Siegel structure theory says
the only vector-space part of a misère quotient is exactly the place where the
non-linearity is gone. Concretely (Plambeck–Siegel, "Misère quotients for impartial
games", JCTA 2008):

  * Every finite misère quotient Q has a KERNEL K — the mutual-divisibility class of
    the product z of all idempotents. K is the maximal subgroup of Q, the quotient
    map x ↦ zx surjects Q ↠ K, and every homomorphism from Q to a group factors
    through it. So K is the canonical (and essentially only) F_2-vector-space shadow
    of Q. (Tame quotients: T_n = K_n ∪ {1,a}, K_n ≅ (Z/2)^n, |T_n|=2^n+2 — never a
    group; the genuine (Z/2)^n is the kernel, not the whole quotient.)

  * THEOREM 6.4: z·Φ(G) = z·Φ(H)  ⟺  G, H have the same NORMAL-play Grundy value.
    So K ≅ the normal-play nim-value group (Z/2)^k under XOR, and the P-portion
    restricted to K is the normal-play {Grundy = 0} set — XOR-LINEAR, a subspace.

THE NO-GO. A genuine quadric is a nonlinear zero set on a vector space. The only
vector-space part of a misère quotient is its kernel K, and on K the P-structure is
(Theorem 6.4) the linear normal-play one. The genuine misère non-linearity lives
OFF the kernel, among the non-group elements (the "fickle units"), where there is
no ambient vector space and "quadric" is not even defined. So the misère quotient
places its non-linearity exactly where the quadratic-form framing cannot reach it.

This is the misère analog of the frame-blind Sp(B) no-go: a structural reason the
route cannot work, not merely an empty search. (Caveat: Theorem 6.4 carries a
regularity hypothesis on the closed game set A; it holds for the regular finite
quotients that arise in practice.)

This script verifies the no-go on R8, the SMALLEST WILD quotient (the most
promising-looking candidate: order 8, wild, genuinely non-linear): it extracts the
kernel, coordinatises it as F_2^k, projects the P-portion in, and shows P∩K is the
linear normal-play set {0} — while the genuine misère P-element sits outside K.
"""

# --------------------------------------------------------------------- monoid engine

# An element is a^i · m with i∈{0,1} and m ∈ {1, b, b2, c}. The 4×4 product table for
# m·n gives (extra power of a, resulting base); full product multiplies the a-powers.
BASE = ["1", "b", "b2", "c"]
# table[m][n] = (extra_a, base_index)
_MN = {
    ("1", "1"): (0, "1"), ("1", "b"): (0, "b"), ("1", "b2"): (0, "b2"), ("1", "c"): (0, "c"),
    ("b", "b"): (0, "b2"), ("b", "b2"): (0, "b"), ("b", "c"): (1, "b"),
    ("b2", "b2"): (0, "b2"), ("b2", "c"): (1, "b2"),
    ("c", "c"): (0, "b2"),
}


def _mn(m, n):
    return _MN.get((m, n)) or _MN[(n, m)]  # commutative


ELEMENTS = [(i, m) for i in (0, 1) for m in BASE]  # 8 elements


def name(e):
    i, m = e
    a = "a" if i else ""
    body = "" if m == "1" else m
    return (a + body) or "1"


def mul(x, y):
    (i, m), (j, n) = x, y
    extra, base = _mn(m, n)
    return ((i + j + extra) % 2, base)


# --------------------------------------------------------------------- structure tools


def check_axioms(elems):
    assert all(mul(x, y) == mul(y, x) for x in elems for y in elems), "not commutative"
    assert all(mul(mul(x, y), z) == mul(x, mul(y, z))
               for x in elems for y in elems for z in elems), "not associative"
    one = (0, "1")
    assert all(mul(one, x) == x for x in elems), "no identity"


def idempotents(elems):
    return [x for x in elems if mul(x, x) == x]


def kernel(elems):
    """Maximal subgroup K: mutual-divisibility class of z = product of idempotents.
    For a finite commutative monoid this is { x : z·x = x } around the idempotent z,
    which forms a group with identity z."""
    z = (0, "1")
    for e in idempotents(elems):
        z = mul(z, e)
    K = [x for x in elems if mul(z, x) == x]
    return z, K


def coordinatise(z, K):
    """K is an elementary abelian 2-group with identity z. Return a dict element→bitmask
    over a chosen basis, verifying every element has order 2 (x·x = z)."""
    assert all(mul(x, x) == z for x in K), "kernel is not elementary abelian 2-group"
    basis = []
    spanned = {z}
    for x in K:
        if x in spanned:
            continue
        # x is independent of the current span
        new = {mul(x, s) for s in spanned}
        if new & spanned == set():  # genuinely new coset
            basis.append(x)
            spanned |= new
    coord = {}
    for x in K:
        for bits in range(1 << len(basis)):
            acc = z
            for t, bx in enumerate(basis):
                if bits & (1 << t):
                    acc = mul(acc, bx)
            if acc == x:
                coord[x] = bits
                break
    return basis, coord


def is_affine_subspace(points):
    """Is S ⊆ F_2^k an affine subspace (a coset of a linear subspace)? — i.e. linear,
    NOT a genuine quadric. Translate by any point; closed under XOR ⟺ affine."""
    if not points:
        return True
    s0 = points[0]
    shifted = {p ^ s0 for p in points}
    return all((x ^ y) in shifted for x in shifted for y in shifted)


# --------------------------------------------------------------------- the R8 demo

if __name__ == "__main__":
    print("R8 — the smallest wild misère quotient (Plambeck–Siegel)")
    print("  R8 = ⟨a,b,c | a²=1, b³=b, bc=ab, c²=b²⟩,  P = {a, b²}\n")

    check_axioms(ELEMENTS)
    print(f"  elements ({len(ELEMENTS)}): {[name(e) for e in ELEMENTS]}")
    P = [(1, "1"), (0, "b2")]  # {a, b2}
    print(f"  P-portion: {[name(e) for e in P]}")

    idem = idempotents(ELEMENTS)
    print(f"\n  idempotents (x²=x): {[name(e) for e in idem]}")
    z, K = kernel(ELEMENTS)
    print(f"  kernel identity z = {name(z)};  kernel K = {[name(e) for e in K]}  (|K|={len(K)})")

    basis, coord = coordinatise(z, K)
    k = len(basis)
    print(f"  K ≅ (Z/2)^{k}  with basis {[name(b) for b in basis]}  (z ↦ 0)")
    print(f"    coordinates: " + ", ".join(f"{name(e)}↦{coord[e]:0{k}b}" for e in K))

    PK = [coord[e] for e in P if e in K]
    P_outside = [name(e) for e in P if e not in K]
    print(f"\n  P ∩ K  (the only vector-space part) = {PK}  "
          f"= {{0}} = the normal-play {{Grundy=0}} set (Theorem 6.4)")
    print(f"  P elements OUTSIDE the kernel (fickle units, no ambient space): {P_outside}")

    affine = is_affine_subspace(PK)
    print(f"\n  Is P∩K an affine subspace (linear, NOT a genuine quadric)? {affine}")
    print("\nConclusion. On R8 — the smallest, most promising wild quotient — the only")
    print("F_2-vector-space part is the kernel (Z/2)², and there the P-set is the linear")
    print("normal-play {Grundy=0} (Thm 6.4). The genuine misère non-linearity (the element")
    print("a ∈ P) lives among the non-group fickle units, where 'quadric' has no ambient.")
    print("So the misère route cannot realise a genuine Gold quadric — a structural no-go,")
    print("the misère analog of the frame-blind Sp(B) obstruction.")
