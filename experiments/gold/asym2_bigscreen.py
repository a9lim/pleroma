import sys, time, random
sys.path.insert(0, "/tmp")
from asym2_fifo import abstract_value
sys.setrecursionlimit(2000000)
rng = random.Random(31337)
for (kbase, n) in [(7, 200), (8, 40), (9, 10)]:
    pairs = [(i, j) for i in range(kbase) for j in range(i + 1, kbase)]
    fails = 0
    t0 = time.time()
    for _ in range(n):
        edges = frozenset(p for p in pairs if rng.randint(0, 1))
        par = len(edges) & 1
        for want in (0, 1):
            if abstract_value(kbase + 1, edges, want) != par:
                fails += 1
                print(f"FAIL board={kbase+1} edges={sorted(edges)} want={want}", flush=True)
    print(f"board {kbase+1} (k={kbase}+dummy): {n} random patterns x2 wants, "
          f"{fails} failures [{time.time()-t0:.0f}s]", flush=True)
