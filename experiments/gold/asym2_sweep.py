"""Capstone: full lambda sweeps at m=8 (a=1,2,3), all positions incl. full board."""
import sys, time
sys.path.insert(0, "/tmp")
from asym2_probe import make_form
from asym2_dummy import support_value_dummy

sys.setrecursionlimit(1000000)

for a in (1, 2, 3):
    t0 = time.time()
    exact_both = exact_one = fail = 0
    fails = []
    for lam in range(1, 256):
        Q, qd, B = make_form(8, a, lam)
        ok = {}
        for t in (1, 0):
            good = True
            for x in range(256):
                S = tuple(i for i in range(8) if (x >> i) & 1)
                if support_value_dummy(S, qd, B, t, 8) != Q[x]:
                    good = False
                    break
            ok[t] = good
        if ok[0] and ok[1]:
            exact_both += 1
        elif ok[0] or ok[1]:
            exact_one += 1
        else:
            fail += 1
            fails.append(lam)
    print(f"m=8 a={a}: 255 lambda: exact both orientations={exact_both}, "
          f"one={exact_one}, fail={fail} {fails[:10]} [{time.time()-t0:.0f}s]")
