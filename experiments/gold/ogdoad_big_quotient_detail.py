"""Detail pass: the order >= 10 bounded quotients, with backtracking Aut."""
import sys
import ogdoad_misere_subgroup_sweep as s


def aut_count_backtrack(table, is_p):
    """|Aut(M, P)| via backtracking on partial maps (preserves table and P)."""
    n = len(table)
    one = next(e for e in range(n) if all(table[e][x] == x for x in range(n)))
    order = sorted(range(n), key=lambda x: x != one)  # map identity first

    count = 0

    def extend(f, depth):
        nonlocal count
        if depth == len(order):
            count += 1
            return
        x = order[depth]
        used = set(v for v in f if v is not None)
        for y in range(n):
            if y in used or is_p[x] != is_p[y]:
                continue
            f[x] = y
            ok = True
            for u in order[:depth + 1]:
                if f[u] is None:
                    continue
                xu = table[x][u]
                if f[xu] is not None and f[xu] != table[y][f[u]]:
                    ok = False
                    break
                ux = table[u][x]
                if f[ux] is not None and f[ux] != table[f[u]][y]:
                    ok = False
                    break
            if ok:
                extend(f, depth + 1)
            f[x] = None

    f = [None] * n
    f[one] = one
    # identity must map to identity; start depth 1
    extend(f, 1)
    return count


codes = []
for d1 in (1, 3, 5, 7):
    codes.append((d1,))
    for d2 in range(8):
        codes.append((d1, d2))
        for d3 in range(8):
            codes.append((d1, d2, d3))

for code in codes:
    for k in range(2, 5):
        q = s.closed_quotient(code, k, 4, 4, cap=40)
        if q is None:
            continue
        reps, table, is_p, viol = q
        if len(table) >= 10:
            name = "0." + "".join(map(str, code)) + f"@h{k}"
            aut = aut_count_backtrack(table, is_p)
            idem = s.idempotents(table)
            z = s.kernel_idempotent(table)
            print(f"{name}: order {len(table)}, idempotents {idem}, "
                  f"z={z}, |Aut(M,P)|={aut}", flush=True)
            stats = dict(subgroups=0, kernels=0, kernel_singleton=0,
                         kernel_nonsingleton=[], nonaffine=[], exp_gt2=0,
                         max_rank=0)
            s.analyze_quotient(name, reps, table, is_p, viol, verbose=True,
                               stats=stats)
