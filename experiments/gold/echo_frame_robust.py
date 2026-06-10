"""Frame-order robustness: does the m=4 exact hit survive permuting the
triangular order of the cocycle? Also: which orientation wins, per order."""

import itertools
import sys

sys.setrecursionlimit(10000)
sys.path.insert(0, "/tmp")

from echo_charge_probe import make_form  # noqa: E402


def solve_perm(m, qd, Bm, order, p1_maximizes):
    """order[i] = priority of coin i in the triangular splitting."""

    def charge(o, i):
        acc = qd[i] if (o >> i) & 1 else 0
        rel = o
        while rel:
            k = (rel & -rel).bit_length() - 1
            rel &= rel - 1
            if order[k] > order[i]:
                acc ^= Bm[k][i]
        return acc

    def forced(x):
        bits = [i for i in range(m) if (x >> i) & 1]
        if not bits:
            return 0
        memo = {}

        def rec(u, o, last, mover_is_p1):
            if u == 0 and o == 0:
                return 0
            key = (u, o, last, mover_is_p1)
            if key in memo:
                return memo[key]
            legal = []
            for i in bits:
                if i == last:
                    continue
                if (u >> i) & 1:
                    legal.append((i, u ^ (1 << i), o ^ (1 << i)))
                elif (o >> i) & 1:
                    legal.append((i, u, o ^ (1 << i)))
            if not legal:
                res = rec(u, o, -1, not mover_is_p1)
                memo[key] = res
                return res
            maximize = mover_is_p1 == p1_maximizes
            best = None
            for (i, u2, o2) in legal:
                v = charge(o, i) ^ rec(u2, o2, i, not mover_is_p1)
                best = v if best is None else (max(best, v) if maximize else min(best, v))
                if (maximize and best == 1) or (not maximize and best == 0):
                    break
            memo[key] = best
            return best

        return rec(x, 0, -1, True)

    return [forced(x) for x in range(1 << m)]


def main():
    for (m, a, lam, label) in [(4, 1, 1, "Gold (4,1) rank 2"),
                               (4, 1, 2, "bent (4,1,lam=2) rank 4 Arf 0"),
                               (4, 1, 3, "component (4,1,lam=3)")]:
        Qf, Bf, qd, Bm = make_form(m, a, lam)
        Qtab = [Qf(v) for v in range(1 << m)]
        nz = sum(1 for t in Qtab if t == 0)
        print(f"\n=== {label}: |Q=0|={nz}, q_diag={qd} ===")
        hits_max = hits_min = tot = 0
        for perm in itertools.permutations(range(4)):
            order = list(perm)
            FA = solve_perm(m, qd, Bm, order, True)
            FB = solve_perm(m, qd, Bm, order, False)
            tot += 1
            hits_max += FA == Qtab
            hits_min += FB == Qtab
        print(f"  P1-max == Q in {hits_max}/{tot} frame orders; "
              f"P1-min == Q in {hits_min}/{tot}")


if __name__ == "__main__":
    main()
