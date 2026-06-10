"""Bench: m=4 full sweep + m=8 forms, correct solver."""
import sys, time
sys.path.insert(0, "/tmp")
from asym2_probe import make_form, make_charge, solve_position

print("=" * 70)
print("m=4, a=1: all 15 lambda, orientations A (P1 wants 1) / B (P1 wants 0),")
print("cocycle sides lower/upper. 'EXACT' = value(x) == Q(x) for all 16 x.")
print("=" * 70)
hits4 = {}
for lam in range(1, 16):
    Q, qd, B = make_form(4, 1, lam)
    zc = sum(1 for v in Q if v == 0)
    row = []
    for side, sname in [(True, "low"), (False, "up")]:
        ch = make_charge(qd, B, 4, side)
        for tgt, tname in [(1, "A"), (0, "B")]:
            val = [solve_position(x, 4, ch, tgt) for x in range(16)]
            agree = sum(1 for x in range(16) if val[x] == Q[x])
            tag = f"{sname}/{tname}:{agree}"
            if agree == 16:
                tag += "*"
                hits4.setdefault(lam, []).append((side, tgt))
            row.append(tag)
    print(f"lam={lam:2d} |Q=0|={zc:2d}  " + "  ".join(row))
print("EXACT hits at m=4:", sorted(hits4))

print()
print("=" * 70)
print("m=8 forms (lower side), both orientations")
print("=" * 70)
for (m, a, lam) in [(8, 1, 1), (8, 2, 1), (8, 1, 2), (8, 1, 3)]:
    Q, qd, B = make_form(m, a, lam)
    zc = sum(1 for v in Q if v == 0)
    for side, sname in [(True, "low"), (False, "up")]:
        ch = make_charge(qd, B, m, side)
        for tgt, tname in [(1, "A"), (0, "B")]:
            t0 = time.time()
            val = [solve_position(x, m, ch, tgt) for x in range(256)]
            agree = sum(1 for x in range(256) if val[x] == Q[x])
            misses = [x for x in range(256) if val[x] != Q[x]]
            mtxt = ""
            if 0 < len(misses) <= 8:
                mtxt = " misses=" + ",".join(
                    f"x={x}(pc{bin(x).count('1')},Q={Q[x]},v={val[x]})" for x in misses)
            print(f"(m={m},a={a},lam={lam}) |Q=0|={zc:3d} {sname}/{tname}: "
                  f"{agree}/256{' EXACT' if agree == 256 else ''}{mtxt} "
                  f"[{time.time()-t0:.1f}s]")
