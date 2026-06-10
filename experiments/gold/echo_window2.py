"""Window-w ko ECHO with the CORRECT (sigma-in-key) solver.

Variants:
  pass_clears=True : stuck => pass, ko window cleared (round-1 rule).
  pass_clears=False: stuck => pass, window kept; double-pass => forced
                     canonical completion (close open coins ascending).
Sweep (8,2) and (8,1) lower side, orientation A (P1 wants 1), w in {1,2}.
"""
import sys
sys.path.insert(0, "/tmp")
from asym2_probe import make_form, make_charge

def solve(x, m, ch, tgt, w, pass_clears):
    bits = [i for i in range(m) if (x >> i) & 1]
    if not bits:
        return 0
    memo = {}

    def forced_finish(u, o, sigma):
        # double-pass: close open coins ascending, then open+close untouched ones
        for i in sorted(i for i in range(len(bits)) if False):
            pass
        oo = o
        s = sigma
        for i in bits:
            if (oo >> i) & 1:
                s ^= ch(oo, i)
                oo ^= 1 << i
        uu = u
        for i in bits:
            if (uu >> i) & 1:
                s ^= ch(oo, i)
                oo ^= 1 << i
                s ^= ch(oo, i)
                oo ^= 1 << i
                uu ^= 1 << i
        return s

    def rec(u, o, win, mover, sigma, passed):
        if u == 0 and o == 0:
            return sigma
        key = (u, o, win, mover, sigma, passed)
        r = memo.get(key)
        if r is not None:
            return r
        blocked = set(win)
        legal = []
        for i in bits:
            if i in blocked:
                continue
            bit = 1 << i
            if u & bit:
                legal.append((i, u ^ bit, o ^ bit))
            elif o & bit:
                legal.append((i, u, o ^ bit))
        if not legal:
            if pass_clears:
                res = rec(u, o, (), 1 - mover, sigma, False)
            else:
                if passed:  # double pass: forced completion
                    res = forced_finish(u, o, sigma)
                else:
                    res = rec(u, o, win, 1 - mover, sigma, True)
            memo[key] = res
            return res
        want = tgt if mover == 0 else 1 - tgt
        res = 1 - want
        for (i, u2, o2) in legal:
            nwin = (tuple(win) + (i,))[-w:]
            r2 = rec(u2, o2, nwin, 1 - mover, sigma ^ ch(o, i), False)
            if r2 == want:
                res = want
                break
        memo[key] = res
        return res

    return rec(x, 0, (), 0, 0, False)

import time
for (m, a, lam) in [(8, 2, 1), (8, 1, 1)]:
    Q, qd, B = make_form(m, a, lam)
    ch = make_charge(qd, B, m, True)
    for w in (1, 2):
        for pc in (True, False):
            t0 = time.time()
            val = [solve(x, m, ch, 1, w, pc) for x in range(1 << m)]
            agree = sum(1 for x in range(1 << m) if val[x] == Q[x])
            misses = [x for x in range(1 << m) if val[x] != Q[x]]
            mtxt = ""
            if 0 < len(misses) <= 6:
                mtxt = "  misses=" + ",".join(
                    f"{x}(pc{bin(x).count('1')},Q={Q[x]},v={val[x]})" for x in misses)
            print(f"(m={m},a={a}) w={w} pass_clears={pc}: {agree}/{1<<m}"
                  f"{' EXACT' if agree == (1<<m) else ''}{mtxt} [{time.time()-t0:.0f}s]")
