"""Echo-charge probe: does alternating play on the extraspecial group E realize Q?

Model (the [asymmetry] angle made precise):
  V = F_2^m, Q = (scaled) Gold form, B its polar. Fix the bit frame and the
  triangular 2-cocycle c(u,v) = sum_i q_i u_i v_i  XOR  sum_{k>j} B_kj u_k v_j,
  so c(v,v) = Q(v) and c(u,v)+c(v,u) = B(u,v). E = V x F_2 with
  (s,u)(t,v) = (s+t+c(u,v), u+v) is the extraspecial-type extension:
  squaring map = Q, commutator = B.

Game ECHO(x): coins = bits(x), each must be touched exactly twice. Players
  alternate; touching coin i when the open set (coins touched an odd number of
  times) is o charges sigma += c(o, e_i)  [= right multiplication of the running
  word by the lift (0,e_i)].  KO variant: may not touch the coin touched in the
  immediately preceding touch; a stuck player passes (pass clears the ko).
  The complete play word lies over 0 in V, hence equals 1 or z in E.
  Readout = the central character: one player wants the word to be 1 (sigma=0),
  the other wants z (sigma=1).

Question: is the forced value of sigma equal to Q(x) for all x?
  popcount(x) <= 2 with ko: provably yes (ko forces the linked pattern i j i j).
  popcount(x) >= 3: this probe computes the minimax.

Theory cross-checks asserted below:
  * c(v,v) = Q(v), c+c^T = B.
  * echo identity: any word w with endpoint x, played twice -> charge Q(x).
  * chord-linking formula: full double-touch play has
      sigma_final = l_diag(x) + sum_{i<j in bits(x)} B_ij * [chords linked].
  * no-ko mirror collapse: P2 forces sigma = l_diag(x).
"""

import random
import sys
from functools import lru_cache

sys.setrecursionlimit(10000)
sys.path.insert(0, "/Users/a9lim/Work/ogdoad/experiments")

import ogdoad as pl  # noqa: E402


def frob(x: "pl.Nimber", a: int) -> "pl.Nimber":
    for _ in range(a):
        x = x * x
    return x


def nim_trace_val(x: int, m: int) -> int:
    acc = pl.Nimber(x)
    t = pl.Nimber(x)
    for _ in range(m - 1):
        t = t * t
        acc = acc + t
    assert acc.value in (0, 1)
    return acc.value


def make_form(m: int, a: int, lam: int):
    """Q(v) = Tr(lam * v^(1+2^a)) over F_{2^m} (nimber subfield)."""
    L = pl.Nimber(lam)

    def Q(v: int) -> int:
        x = pl.Nimber(v)
        return nim_trace_val((L * x * frob(x, a)).value, m)

    qtab = [Q(v) for v in range(1 << m)]

    def Qf(v):
        return qtab[v]

    def Bf(u, v):
        return qtab[u ^ v] ^ qtab[u] ^ qtab[v]

    qd = [Qf(1 << i) for i in range(m)]
    Bm = [[Bf(1 << i, 1 << j) for j in range(m)] for i in range(m)]
    return Qf, Bf, qd, Bm


def cocycle_full(u: int, v: int, qd, Bm, m: int, lower=True) -> int:
    acc = 0
    for i in range(m):
        if (u >> i) & 1 and (v >> i) & 1:
            acc ^= qd[i]
    for k in range(m):
        if not (u >> k) & 1:
            continue
        for j in range(m):
            if not (v >> j) & 1:
                continue
            if (lower and k > j) or (not lower and k < j):
                acc ^= Bm[k][j]
    return acc


def charge_move(o: int, i: int, qd, Brows, lower_mask_above, lower=True) -> int:
    """c(o, e_i): the charge for touching coin i when the open set is o."""
    acc = qd[i] if (o >> i) & 1 else 0
    mask = lower_mask_above[i] if lower else ((1 << i) - 1)
    rel = o & mask
    while rel:
        k = (rel & -rel).bit_length() - 1
        rel &= rel - 1
        acc ^= Brows[k][i]
    return acc


# ----------------------------------------------------------------- game solver


def solve_form(m, Qf, qd, Bm, ko=True, lower=True, verbose=True):
    """Return (F_A, F_B): F_A[x] = forced sigma when P1 maximizes sigma (P2 min),
    F_B[x] = forced sigma when P1 minimizes (P2 max). Both with the same rules."""
    above = [(~((1 << (i + 1)) - 1)) & ((1 << m) - 1) for i in range(m)]

    def forced(x, p1_maximizes):
        bits = [i for i in range(m) if (x >> i) & 1]
        if not bits:
            return 0

        memo = {}

        def rec(u, o, last, mover_is_p1):
            # returns forced FUTURE charge (F_2) under optimal play
            if u == 0 and o == 0:
                return 0
            key = (u, o, last, mover_is_p1)
            if key in memo:
                return memo[key]
            legal = []
            for i in bits:
                if ko and i == last:
                    continue
                ui, oi = (u >> i) & 1, (o >> i) & 1
                if ui:  # untouched -> open
                    legal.append((i, u ^ (1 << i), o ^ (1 << i)))
                elif oi:  # open -> closed
                    legal.append((i, u, o ^ (1 << i)))
            if not legal:  # stuck: pass, ko clears
                res = rec(u, o, -1, not mover_is_p1)
                memo[key] = res
                return res
            maximize = mover_is_p1 == p1_maximizes
            vals = []
            for (i, u2, o2) in legal:
                ch = charge_move(o, i, qd, Bm, above, lower)
                vals.append(ch ^ rec(u2, o2, i, not mover_is_p1))
                # short-circuit
                if maximize and vals[-1] == 1:
                    break
                if not maximize and vals[-1] == 0:
                    break
            res = max(vals) if maximize else min(vals)
            memo[key] = res
            return res

        return rec(x, 0, -1, True)

    n = 1 << m
    FA = [forced(x, True) for x in range(n)]
    FB = [forced(x, False) for x in range(n)]
    return FA, FB


# ----------------------------------------------------------------- ANF fitting


def anf(table):
    """Mobius transform: truth table (list over F_2^k) -> ANF coefficients."""
    n = len(table)
    k = n.bit_length() - 1
    co = list(table)
    for i in range(k):
        bit = 1 << i
        for mask in range(n):
            if mask & bit:
                co[mask] ^= co[mask ^ bit]
    return co


def describe_table(table, m, Bm, Qtab, ltab):
    n = 1 << m
    co = anf(table)
    deg = max((bin(mask).count("1") for mask in range(n) if co[mask]), default=0)
    out = [f"deg={deg}"]
    if table == Qtab:
        out.append("== Q EXACTLY")
        return " ".join(out)
    agree = sum(1 for v in range(n) if table[v] == Qtab[v])
    out.append(f"agree with Q: {agree}/{n}")
    if table == ltab:
        out.append("== l_diag exactly")
    if deg <= 2 and deg > 0:
        # polar of the fitted quadratic part vs B
        same_B = all(
            (co[(1 << i) | (1 << j)] if i != j else 0) == Bm[i][j]
            for i in range(m)
            for j in range(i + 1, m)
        )
        out.append(f"quadratic; polar == B: {same_B}")
        if same_B:
            diag = [co[1 << i] for i in range(m)]
            out.append(f"refinement of B with diagonal {diag}")
    return " ".join(out)


# ----------------------------------------------------------------- cross-checks


def run_checks(m, a, lam):
    Qf, Bf, qd, Bm = make_form(m, a, lam)
    n = 1 << m
    # cocycle identities
    rng = random.Random(2026)
    pairs = (
        [(u, v) for u in range(n) for v in range(n)]
        if m <= 4
        else [(rng.randrange(n), rng.randrange(n)) for _ in range(400)]
    )
    for (u, v) in pairs:
        assert cocycle_full(v, v, qd, Bm, m) == Qf(v)
        assert cocycle_full(u, v, qd, Bm, m) ^ cocycle_full(v, u, qd, Bm, m) == Bf(u, v)
    above = [(~((1 << (i + 1)) - 1)) & ((1 << m) - 1) for i in range(m)]

    # echo identity: random word w (single-bit moves), played twice -> Q(endpoint)
    for _ in range(200):
        k = rng.randrange(1, 9)
        word = [rng.randrange(m) for _ in range(k)]
        sigma, o = 0, 0
        for i in word + word:
            sigma ^= charge_move(o, i, qd, Bm, above)
            o ^= 1 << i
        assert o == 0
        endpoint = 0
        for i in word:
            endpoint ^= 1 << i
        assert sigma == Qf(endpoint), "echo identity failed"

    # chord-linking formula on random full double-touch plays
    for _ in range(200):
        x = rng.randrange(1, n)
        bits = [i for i in range(m) if (x >> i) & 1]
        sched = bits + bits
        rng.shuffle(sched)
        sigma, o, times = 0, 0, {i: [] for i in bits}
        for t, i in enumerate(sched):
            sigma ^= charge_move(o, i, qd, Bm, above)
            o ^= 1 << i
            times[i].append(t)
        link = 0
        for ii in range(len(bits)):
            for jj in range(ii + 1, len(bits)):
                i, j = bits[ii], bits[jj]
                (a1, b1), (a2, b2) = times[i], times[j]
                linked = (a1 < a2 < b1 < b2) or (a2 < a1 < b2 < b1)
                if linked:
                    link ^= Bm[i][j]
        ell = 0
        for i in bits:
            ell ^= qd[i]
        assert sigma == ell ^ link, "chord-linking formula failed"
    print(f"  [checks pass: cocycle identities, echo identity, linking formula]")
    return Qf, Bf, qd, Bm


def main():
    cases = []
    # (m, a, lam, label)
    cases.append((4, 1, 1, "Gold (4,1)  rank 2, Arf 1"))
    # bent scaled component over F_16
    Qf, _, _, _ = make_form(4, 1, 1)
    for lam in range(1, 16):
        Qb, _, _, _ = make_form(4, 1, lam)
        z = sum(1 for v in range(16) if Qb(v) == 0)
        if z in (6, 10):
            cases.append((4, 1, lam, f"bent Gold component (4,1,lam={lam}) rank 4"))
            break
    cases.append((8, 1, 1, "Gold (8,1)  rank 6, Arf 1"))
    cases.append((8, 2, 1, "Gold (8,2)  rank 4, Arf 1"))
    for lam in range(1, 256):
        Qb, _, _, _ = make_form(8, 1, lam)
        z = sum(1 for v in range(256) if Qb(v) == 0)
        if z in (120, 136):
            cases.append((8, 1, lam, f"bent Gold component (8,1,lam={lam}) rank 8"))
            break

    for (m, a, lam, label) in cases:
        n = 1 << m
        print(f"\n=== {label} ===")
        Qf, Bf, qd, Bm = run_checks(m, a, lam)
        Qtab = [Qf(v) for v in range(n)]
        ltab = []
        for v in range(n):
            e = 0
            for i in range(m):
                if (v >> i) & 1:
                    e ^= qd[i]
            ltab.append(e)
        nz = sum(1 for t in Qtab if t == 0)
        print(f"  |{{Q=0}}| = {nz}/{n}   q_diag = {qd}")

        for ko in (False, True):
            for lower in ((True,) if not ko else (True, False)):
                FA, FB = solve_form(m, Qf, qd, Bm, ko=ko, lower=lower)
                tag = f"ko={'on ' if ko else 'off'} cocycle={'lower' if lower else 'upper'}"
                # sanity: popcount<=2 forced to Q when ko on
                if ko:
                    for x in range(n):
                        if bin(x).count("1") <= 2:
                            assert FA[x] == Qtab[x] == FB[x], (
                                f"popcount<=2 not forced to Q at x={x}"
                            )
                print(f"  [{tag}] P1-max: {describe_table(FA, m, Bm, Qtab, ltab)}")
                print(f"  [{tag}] P1-min: {describe_table(FB, m, Bm, Qtab, ltab)}")
                if not ko:
                    # Theorem C: P2 (minimizer in orientation A) forces l_diag:
                    # whenever l(x)=0, FA[x] must be 0.
                    viol = [x for x in range(n) if ltab[x] == 0 and FA[x] == 1]
                    print(f"      mirror-collapse check (l=0 => FA=0): violations {len(viol)}")


if __name__ == "__main__":
    main()
