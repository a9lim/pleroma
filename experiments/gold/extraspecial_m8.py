"""m=8 sweep with the correct solver + decision-nondegeneracy counts."""
import time
from extraspecial_core import *

validate()

forms = [
    ("(8,1) lam=1  rank6", gold_q(8, 1, 1)),
    ("(8,2) lam=1  rank4", gold_q(8, 2, 1)),
    ("(8,1) lam=2  bent8", gold_q(8, 1, 2)),
    ("(8,1) lam=3  bent8", gold_q(8, 1, 3)),
]

# first: decision-nondegeneracy on the m=4 exact instances
print("=== m=4 decision-nondegeneracy (ko=self, P1max) ===")
for lam in (1, 2, 12, 13, 14):
    Q = gold_q(4, 1, lam); B = polar(Q, 4)
    tot_choice = 0; per = []
    for x in range(16):
        v, ch = echo_value(x, Q, B, 4, maxfirst=True)
        tot_choice += ch
        per.append(ch)
    print(f"lam={lam:2d}: total choice-states={tot_choice}, per-position={per}")

print("\n=== m=8, ko=self ===")
for name, Q in forms:
    B = polar(Q, 8)
    for mf in (True, False):
        t0 = time.time()
        misses = []
        tot_choice = 0
        for x in range(256):
            v, ch = echo_value(x, Q, B, 8, maxfirst=mf)
            tot_choice += ch
            if v != Q[x]:
                misses.append((x, bin(x).count('1'), Q[x], v))
        dt = time.time() - t0
        ori = 'P1max' if mf else 'P1min'
        print(f"{name} {ori}: agree={256-len(misses)}/256 "
              f"choice-states={tot_choice} ({dt:.0f}s)")
        if misses and len(misses) <= 12:
            print(f"   misses (x, popcount, Q, val): {misses}")
