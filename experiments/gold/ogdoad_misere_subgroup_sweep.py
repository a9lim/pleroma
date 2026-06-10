"""Sweep bounded misere quotients of octal games; analyze EVERY maximal subgroup.

New surface vs the repo's probes (octal_hunt.rs / misere_kernel.py):
  * octal_hunt only fits P-sets of quotients that are GLOBALLY (Z/2)^k;
  * misere_kernel.py analyzes only the kernel K of R8.
This script, for every bounded quotient (signature-closed commutative monoid):
  1. enumerates all idempotents e and maximal subgroups G_e = unit group of eM;
  2. checks every G_e for exponent 2 (elementary abelian) -- any exponent-4
     element would be news (a place where a squaring map V -> Z could even live);
  3. coordinatizes each elementary-abelian G_e as F_2^k and quadric-fits the
     trace P /\ G_e (ANF/Mobius, degree <= 2, polar rank);
  4. records |P /\ K| (kernel singleton prediction, P-S Thm 6.4);
  5. computes |Aut(M, P)| for small quotients -- the symmetry budget any
     E-equivariant (extraspecial) rule would need to occupy.

Bounded-observation discipline: signature classes w.r.t. a bounded test set,
closed under products up to a class cap; congruence is spot-checked on in-bound
element pairs (violations counted, quotient flagged).  Coarser-than-true
quotients are possible; anything interesting gets re-examined, not trusted.
"""

import sys
from itertools import combinations_with_replacement, permutations

# ---------------------------------------------------------------- octal engine


def octal_moves(code, pos):
    out = set()
    lst = list(pos)
    for idx in range(len(lst)):
        n = lst[idx]
        base = lst[:idx] + lst[idx + 1:]
        for k in range(1, n + 1):
            d = code[k - 1] if k - 1 < len(code) else 0
            rem = n - k
            if rem == 0:
                if d & 1:
                    out.add(tuple(sorted(base)))
            else:
                if d & 2:
                    out.add(tuple(sorted(base + [rem])))
                if d & 4:
                    for a in range(1, rem // 2 + 1):
                        out.add(tuple(sorted(base + [a, rem - a])))
    return out


def make_outcome(code):
    memo = {}

    def is_n(pos):  # misere: terminal => N (cannot-move wins)
        r = memo.get(pos)
        if r is not None:
            return r
        res = True
        nxt = octal_moves(code, pos)
        if nxt:
            res = False
            for q in nxt:
                if not is_n(q):
                    res = True
                    break
        memo[pos] = res
        return res

    return is_n


def multisets(atoms, maxlen):
    out = [()]
    for length in range(1, maxlen + 1):
        out.extend(combinations_with_replacement(atoms, length))
    return out


# ---------------------------------------------------- bounded, closed quotient


def closed_quotient(code, max_heap, elem_bound, test_bound, cap=64):
    """Signature-closed bounded quotient. Returns None if class count > cap.
    Result: (reps, table, is_p, violations) -- a finite commutative monoid
    table on signature classes, P-portion, and congruence spot-check count."""
    atoms = list(range(1, max_heap + 1))
    tests = multisets(atoms, test_bound)
    is_n = make_outcome(code)
    sig_memo = {}

    def sig_of(g):
        g = tuple(sorted(g))
        s = sig_memo.get(g)
        if s is None:
            s = tuple(is_n(tuple(sorted(g + t))) for t in tests)
            sig_memo[g] = s
        return s

    sig_to_class = {}
    reps = []

    def cls(g):
        s = sig_of(g)
        c = sig_to_class.get(s)
        if c is None:
            c = len(reps)
            sig_to_class[s] = c
            reps.append(tuple(sorted(g)))
        return c

    elements = multisets(atoms, elem_bound)
    for g in elements:
        cls(g)
        if len(reps) > cap:
            return None

    # close the class set under products of representatives
    frontier = list(range(len(reps)))
    while frontier:
        new_frontier = []
        n_before = len(reps)
        for i in frontier:
            for j in range(len(reps)):
                cls(reps[i] + reps[j])
                if len(reps) > cap:
                    return None
        for c in range(n_before, len(reps)):
            new_frontier.append(c)
        frontier = new_frontier

    nc = len(reps)
    table = [[0] * nc for _ in range(nc)]
    for i in range(nc):
        for j in range(i, nc):
            c = cls(reps[i] + reps[j])
            table[i][j] = table[j][i] = c
    is_p = [not is_n(reps[i]) for i in range(nc)]

    # congruence spot-check on in-bound element pairs
    violations = 0
    cls_of_elem = {g: cls(g) for g in elements}
    for gi in range(len(elements)):
        for gj in range(gi, len(elements)):
            g, h = elements[gi], elements[gj]
            if len(g) + len(h) <= elem_bound:
                prod = tuple(sorted(g + h))
                if cls_of_elem[prod] != table[cls_of_elem[g]][cls_of_elem[h]]:
                    violations += 1
    return reps, table, is_p, violations


# -------------------------------------------------------- subgroup analysis


def idempotents(table):
    return [e for e in range(len(table)) if table[e][e] == e]


def kernel_idempotent(table):
    idem = idempotents(table)
    z = idem[0]
    for e in idem[1:]:
        z = table[z][e]
    return z


def maximal_subgroup(table, e):
    """G_e: the unit group of eM with identity e."""
    n = len(table)
    eM = sorted({table[e][x] for x in range(n)})
    in_eM = set(eM)
    G = []
    for x in eM:
        if table[e][x] != x:
            continue
        if any(table[x][y] == e for y in in_eM):
            G.append(x)
    return G


def is_exponent_2(table, e, G):
    return all(table[x][x] == e for x in G)


def coordinatize(table, e, G):
    """G elementary abelian with identity e -> dict element -> bitmask."""
    basis = []
    span = {e}
    for x in G:
        if x not in span:
            basis.append(x)
            span |= {table[x][s] for s in span}
    coord = {}
    for x in G:
        for bits in range(1 << len(basis)):
            acc = e
            for t, bx in enumerate(basis):
                if bits & (1 << t):
                    acc = table[acc][bx]
            if acc == x:
                coord[x] = bits
                break
    assert len(coord) == len(G)
    return basis, coord


def is_affine(points):
    if not points:
        return True
    s0 = points[0]
    shifted = {p ^ s0 for p in points}
    return all((x ^ y) in shifted for x in shifted for y in shifted)


def anf_quadric_fit(zero_set, k):
    """Port of forms::quadric_fit::fit_f2_quadratic.  Returns None if the set
    is not the zero set of a degree-<=2 boolean polynomial; else
    (constant, diag, bmat, polar_rank)."""
    n = 1 << k
    coeffs = [True] * n
    for v in zero_set:
        coeffs[v] = False
    for i in range(k):
        bit = 1 << i
        for mask in range(n):
            if mask & bit:
                coeffs[mask] ^= coeffs[mask ^ bit]
    for mask in range(n):
        if coeffs[mask] and bin(mask).count("1") > 2:
            return None
    constant = coeffs[0]
    diag = [coeffs[1 << i] for i in range(k)]
    bmat = [[False] * k for _ in range(k)]
    for i in range(k):
        for j in range(i + 1, k):
            if coeffs[(1 << i) | (1 << j)]:
                bmat[i][j] = bmat[j][i] = True
    # polar rank: F_2 Gaussian elimination on the alternating matrix
    m = [row[:] for row in bmat]
    rank = 0
    rows = list(range(k))
    for col in range(k):
        piv = next((r for r in rows if m[r][col]), None)
        if piv is None:
            continue
        rows.remove(piv)
        rank += 1
        for r in rows:
            if m[r][col]:
                for c in range(k):
                    m[r][c] ^= m[piv][c]
    return constant, diag, bmat, rank


def automorphism_count(table, is_p, limit_n=9):
    """|Aut(M, P)| by brute force (small monoids only)."""
    n = len(table)
    if n > limit_n:
        return None
    idx = list(range(n))
    count = 0
    # identity element: the e with e*x = x for all x
    one = next(e for e in idx if all(table[e][x] == x for x in idx))
    others = [x for x in idx if x != one]
    for perm in permutations(others):
        f = [0] * n
        f[one] = one
        for src, dst in zip(others, perm):
            f[src] = dst
        if any(is_p[x] != is_p[f[x]] for x in idx):
            continue
        if all(f[table[x][y]] == table[f[x]][f[y]] for x in idx for y in idx):
            count += 1
    return count


# ------------------------------------------------------------------ analysis


def analyze_quotient(name, reps, table, is_p, violations, verbose=False,
                     stats=None):
    n = len(table)
    idem = idempotents(table)
    z = kernel_idempotent(table)
    findings = []
    for e in idem:
        G = maximal_subgroup(table, e)
        exp2 = is_exponent_2(table, e, G)
        if not exp2:
            findings.append(("EXPONENT>2", e, len(G)))
            if stats is not None:
                stats["exp_gt2"] += 1
            continue
        basis, coord = coordinatize(table, e, G)
        k = len(basis)
        ptrace = sorted(coord[x] for x in G if is_p[x])
        affine = is_affine(ptrace)
        if stats is not None:
            stats["subgroups"] += 1
            stats["max_rank"] = max(stats["max_rank"], k)
            if e == z:
                stats["kernels"] += 1
                if len(ptrace) == 1:
                    stats["kernel_singleton"] += 1
                else:
                    stats["kernel_nonsingleton"].append(
                        (name, len(G), ptrace))
            if not affine:
                fit = anf_quadric_fit(ptrace, k)
                stats["nonaffine"].append((name, e, k, ptrace, fit))
        if verbose:
            tag = "KERNEL" if e == z else ("UNITS" if e != z and all(
                table[e][x] == x for x in range(n)) else "mid")
            print(f"    e={e} ({tag}) |G_e|={len(G)} rank={k} "
                  f"P-trace={ptrace} affine={affine}")
        findings.append((e, len(G), k, tuple(ptrace), affine))
    return idem, z, findings


# ------------------------------------------------------------------- R8 check

R8_BASE = ["1", "b", "b2", "c"]
_MN = {
    ("1", "1"): (0, "1"), ("1", "b"): (0, "b"), ("1", "b2"): (0, "b2"),
    ("1", "c"): (0, "c"), ("b", "b"): (0, "b2"), ("b", "b2"): (0, "b"),
    ("b", "c"): (1, "b"), ("b2", "b2"): (0, "b2"), ("b2", "c"): (1, "b2"),
    ("c", "c"): (0, "b2"),
}


def r8_tables():
    elems = [(i, m) for i in (0, 1) for m in R8_BASE]

    def mul(x, y):
        (i, m), (j, n) = x, y
        extra, base = _MN.get((m, n)) or _MN[(n, m)]
        return ((i + j + extra) % 2, base)

    index = {e: k for k, e in enumerate(elems)}
    table = [[index[mul(x, y)] for y in elems] for x in elems]
    P = {index[(1, "1")], index[(0, "b2")]}
    is_p = [k in P for k in range(len(elems))]
    names = []
    for (i, m) in elems:
        names.append((("a" if i else "") + ("" if m == "1" else m)) or "1")
    return names, table, is_p


def main():
    print("== R8 (hardcoded Plambeck-Siegel table): all maximal subgroups ==")
    names, table, is_p = r8_tables()
    stats0 = dict(subgroups=0, kernels=0, kernel_singleton=0,
                  kernel_nonsingleton=[], nonaffine=[], exp_gt2=0, max_rank=0)
    analyze_quotient("R8", None, table, is_p, 0, verbose=True, stats=stats0)
    aut = automorphism_count(table, is_p)
    print(f"    |Aut(R8, P)| = {aut}")
    print(f"    R8 stats: {stats0}")

    print("\n== octal sweep ==")
    max_heap = 4
    elem_bound, test_bound = 4, 4
    codes = []
    for d1 in (1, 3, 5, 7):
        codes.append((d1,))
        for d2 in range(8):
            codes.append((d1, d2))
            for d3 in range(8):
                codes.append((d1, d2, d3))
    print(f"codes: {len(codes)}, heap cutoffs 2..{max_heap}, "
          f"bounds elem<={elem_bound}/test<={test_bound}")

    stats = dict(subgroups=0, kernels=0, kernel_singleton=0,
                 kernel_nonsingleton=[], nonaffine=[], exp_gt2=0, max_rank=0)
    n_quot = 0
    n_capped = 0
    n_viol = 0
    order_hist = {}
    aut_max = 0
    for code in codes:
        for k in range(2, max_heap + 1):
            q = closed_quotient(code, k, elem_bound, test_bound, cap=40)
            if q is None:
                n_capped += 1
                continue
            reps, table, is_p, violations = q
            n_quot += 1
            if violations:
                n_viol += 1
                continue  # untrusted table; skip (counted)
            name = "0." + "".join(map(str, code)) + f"@h{k}"
            order_hist[len(table)] = order_hist.get(len(table), 0) + 1
            analyze_quotient(name, reps, table, is_p, violations, stats=stats)
            a = automorphism_count(table, is_p)
            if a is not None:
                aut_max = max(aut_max, a)

    print(f"\nquotients analyzed: {n_quot} (capped: {n_capped}, "
          f"congruence violations skipped: {n_viol})")
    print(f"order histogram: {dict(sorted(order_hist.items()))}")
    print(f"maximal subgroups analyzed: {stats['subgroups']}; "
          f"max rank seen: {stats['max_rank']}; "
          f"exponent>2 subgroups: {stats['exp_gt2']}")
    print(f"kernels: {stats['kernels']}; with |P /\\ K| = 1: "
          f"{stats['kernel_singleton']}")
    if stats["kernel_nonsingleton"]:
        print("KERNEL NON-SINGLETON cases:")
        for item in stats["kernel_nonsingleton"][:20]:
            print("   ", item)
    if stats["nonaffine"]:
        print("NON-AFFINE P-traces (POSSIBLE QUADRIC HOSTS):")
        for item in stats["nonaffine"][:20]:
            print("   ", item)
    else:
        print("non-affine P-traces found: NONE "
              "(every P /\\ G_e affine across the sweep)")
    print(f"max |Aut(M, P)| over small quotients: {aut_max}")


if __name__ == "__main__":
    main()
