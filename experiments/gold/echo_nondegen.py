"""Decision-(non)degeneracy instrument for the corrected ECHO-ko game.

A game is decision-degenerate (a 'clock') if from every reachable state all
legal moves lead to the same outcome for the mover. T2 is provably a clock.
Question: are the ECHO-ko exact hits at m=4 (Gold (4,1,1) and bent (4,1,2))
genuinely adversarial -- do positions exist with both winning and losing moves?
Counts, per starting x: reachable states under the game tree where the mover's
legal options have DIFFERING minimax outcomes (i.e. mistakes are possible).
"""
import sys
sys.path.insert(0, "/tmp")
from asym2_probe import make_form, make_charge, solve_position
from functools import lru_cache

def analyse(m, a, lam, side, tgt):
    Q, qd, B = make_form(m, a, lam)
    ch = make_charge(qd, B, m, side)
    exact = all(solve_position(x, m, ch, tgt) == Q[x] for x in range(1 << m))

    total_choice_states = 0   # reachable states with >= 2 legal moves
    total_mistake_states = 0  # ... where options' values differ for the mover
    for x in range(1, 1 << m):
        bits = [i for i in range(m) if (x >> i) & 1]
        memo = {}

        def val(u, o, last, mover, sigma):
            if u == 0 and o == 0:
                return sigma
            key = (u, o, last, mover, sigma)
            if key in memo:
                return memo[key]
            legal = []
            for i in bits:
                if i == last:
                    continue
                bit = 1 << i
                if u & bit:
                    legal.append((i, u ^ bit, o ^ bit))
                elif o & bit:
                    legal.append((i, u, o ^ bit))
            if not legal:
                res = val(u, o, -1, 1 - mover, sigma)
            else:
                want = tgt if mover == 0 else 1 - tgt
                outs = [val(u2, o2, i, 1 - mover, sigma ^ ch(o, i))
                        for (i, u2, o2) in legal]
                res = want if want in outs else 1 - want
            memo[key] = res
            return res

        # walk reachable states, count choice/mistake states
        seen = set()
        stack = [(x, 0, -1, 0, 0)]
        cs = ms = 0
        while stack:
            (u, o, last, mover, sigma) = stack.pop()
            if (u, o, last, mover, sigma) in seen:
                continue
            seen.add((u, o, last, mover, sigma))
            if u == 0 and o == 0:
                continue
            legal = []
            for i in bits:
                if i == last:
                    continue
                bit = 1 << i
                if u & bit:
                    legal.append((i, u ^ bit, o ^ bit))
                elif o & bit:
                    legal.append((i, u, o ^ bit))
            if not legal:
                stack.append((u, o, -1, 1 - mover, sigma))
                continue
            outs = []
            for (i, u2, o2) in legal:
                s2 = sigma ^ ch(o, i)
                outs.append(val(u2, o2, i, 1 - mover, s2))
                stack.append((u2, o2, i, 1 - mover, s2))
            if len(legal) >= 2:
                cs += 1
                if len(set(outs)) > 1:
                    ms += 1
        total_choice_states += cs
        total_mistake_states += ms
    return exact, total_choice_states, total_mistake_states

for (m, a, lam, side, tgt, name) in [
    (4, 1, 1, True, 1, "Gold (4,1) low/A"),
    (4, 1, 1, True, 0, "Gold (4,1) low/B"),
    (4, 1, 2, True, 1, "bent (4,1,lam=2) low/A"),
    (4, 1, 2, True, 0, "bent (4,1,lam=2) low/B"),
]:
    exact, cs, ms = analyse(m, a, lam, side, tgt)
    print(f"{name:<26} exact={exact}  choice-states={cs}  mistake-states={ms}"
          f"  {'NON-DEGENERATE' if ms > 0 else 'clock'}")
