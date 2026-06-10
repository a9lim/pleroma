"""Window-w ko variant of the echo-charge game + failure localization.

Hypothesis: ko with memory window w forces the linked (round-robin) pattern on
popcount <= w+1 positions, so the minimax value equals Q there; adversarial
unlinking re-enters at popcount >= w+2.
"""

import sys

sys.setrecursionlimit(10000)
sys.path.insert(0, "/tmp")

from echo_charge_probe import make_form, charge_move, anf  # noqa: E402


def solve_form_window(m, qd, Bm, w, lower=True):
    above = [(~((1 << (i + 1)) - 1)) & ((1 << m) - 1) for i in range(m)]

    def forced(x, p1_maximizes):
        bits = [i for i in range(m) if (x >> i) & 1]
        if not bits:
            return 0
        memo = {}

        def rec(u, o, window, mover_is_p1):
            if u == 0 and o == 0:
                return 0
            key = (u, o, window, mover_is_p1)
            if key in memo:
                return memo[key]
            legal = []
            for i in bits:
                if i in window:
                    continue
                ui, oi = (u >> i) & 1, (o >> i) & 1
                if ui:
                    legal.append((i, u ^ (1 << i), o ^ (1 << i)))
                elif oi:
                    legal.append((i, u, o ^ (1 << i)))
            if not legal:  # stuck: pass clears the window
                res = rec(u, o, (), not mover_is_p1)
                memo[key] = res
                return res
            maximize = mover_is_p1 == p1_maximizes
            best = None
            for (i, u2, o2) in legal:
                ch = charge_move(o, i, qd, Bm, above, lower)
                nw = (window + (i,))[-w:] if w > 0 else ()
                v = ch ^ rec(u2, o2, nw, not mover_is_p1)
                if best is None:
                    best = v
                elif maximize:
                    best = max(best, v)
                else:
                    best = min(best, v)
                if (maximize and best == 1) or (not maximize and best == 0):
                    break
            memo[key] = best
            return best

        return rec(x, 0, (), True)

    n = 1 << m
    FA = [forced(x, True) for x in range(n)]
    FB = [forced(x, False) for x in range(n)]
    return FA, FB


def report(label, F, Qtab, m):
    n = 1 << m
    agree = sum(1 for v in range(n) if F[v] == Qtab[v])
    mism = [x for x in range(n) if F[x] != Qtab[x]]
    by_pop = {}
    for x in mism:
        by_pop.setdefault(bin(x).count("1"), []).append(x)
    co = anf(F)
    deg = max((bin(k).count("1") for k in range(n) if co[k]), default=0)
    s = f"  {label}: agree {agree}/{n} deg={deg}"
    if not mism:
        s += "  == Q EXACTLY"
    else:
        pops = {k: len(v) for k, v in sorted(by_pop.items())}
        s += f"  mismatch popcounts: {pops}"
        if len(mism) <= 6:
            s += f"  at {mism}"
    print(s)


def main():
    # locate the failing position of the bent (4,1,lam=2) P1-min run, w=1
    print("=== bent (4,1,lam=2), window w=1: localize the P1-min failure ===")
    Qf, Bf, qd, Bm = make_form(4, 1, 2)
    Qtab = [Qf(v) for v in range(16)]
    FA, FB = solve_form_window(4, qd, Bm, 1)
    report("w=1 P1-max", FA, Qtab, 4)
    report("w=1 P1-min", FB, Qtab, 4)

    print("\n=== window w=2, m=4 cases ===")
    for (m, a, lam, label) in [(4, 1, 1, "Gold (4,1)"), (4, 1, 2, "bent (4,1,lam=2)")]:
        Qf, Bf, qd, Bm = make_form(m, a, lam)
        Qtab = [Qf(v) for v in range(1 << m)]
        FA, FB = solve_form_window(m, qd, Bm, 2)
        report(f"{label} w=2 P1-max", FA, Qtab, m)
        report(f"{label} w=2 P1-min", FB, Qtab, m)

    print("\n=== window w=2 and w=3, m=8 (8,2) Gold rank 4 ===")
    Qf, Bf, qd, Bm = make_form(8, 2, 1)
    Qtab = [Qf(v) for v in range(256)]
    for w in (2, 3):
        FA, FB = solve_form_window(8, qd, Bm, w)
        report(f"(8,2) w={w} P1-max", FA, Qtab, 8)
        report(f"(8,2) w={w} P1-min", FB, Qtab, 8)

    print("\n=== window w=2 and w=3, m=8 (8,1) Gold rank 6 ===")
    Qf, Bf, qd, Bm = make_form(8, 1, 1)
    Qtab = [Qf(v) for v in range(256)]
    for w in (2, 3):
        FA, FB = solve_form_window(8, qd, Bm, w)
        report(f"(8,1) w={w} P1-max", FA, Qtab, 8)
        report(f"(8,1) w={w} P1-min", FB, Qtab, 8)


if __name__ == "__main__":
    main()
