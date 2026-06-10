#!/usr/bin/env python3
"""Independent cross-check (different algorithm from /tmp/loopy_audit.py):
standard retrograde W/L/D labeling on the (position, mover) graph.
Position = sorted tuple of component names; moving a component to 0 deletes it.
"""

GAMES = {  # name -> (Left options, Right options); "0" means component vanishes
    "0":     ([], []),
    "*":     (["0"], ["0"]),
    "*2":    (["0", "*"], ["0", "*"]),
    "1":     (["0"], []),
    "-1":    ([], ["0"]),
    "up":    (["0"], ["*"]),
    "down":  (["*"], ["0"]),
    "on":    (["on"], []),
    "off":   ([], ["off"]),
    "over":  (["0"], ["over"]),
    "under": (["under"], ["0"]),
    "dud":   (["dud"], ["dud"]),
}
NEG = {"0": "0", "*": "*", "*2": "*2", "1": "-1", "-1": "1", "up": "down",
       "down": "up", "on": "off", "off": "on", "over": "under",
       "under": "over", "dud": "dud"}

def canon(pos):
    return tuple(sorted(pos))

def moves(pos, player):  # player 0 = Left, 1 = Right
    out = set()
    for i, comp in enumerate(pos):
        for o in GAMES[comp][player]:
            nxt = list(pos)
            if o == "0":
                nxt.pop(i)
            else:
                nxt[i] = o
            out.add(canon(nxt))
    return out

def solve(start):
    """node = (pos, mover). Returns label dict: 'W' mover wins, 'L' mover
    loses, 'D' neither forces a win (infinite play)."""
    start = canon(start)
    nodes, stack = set(), [(start, 0), (start, 1)]
    succ = {}
    while stack:
        n = stack.pop()
        if n in nodes:
            continue
        nodes.add(n)
        pos, pl = n
        succ[n] = [(s, 1 - pl) for s in moves(pos, pl)]
        stack.extend(succ[n])
    pred = {n: [] for n in nodes}
    for n, ss in succ.items():
        for s in ss:
            pred[s].append(n)
    label, work = {}, []
    deg = {n: len(succ[n]) for n in nodes}
    for n in nodes:
        if deg[n] == 0:
            label[n] = "L"  # stuck mover loses
            work.append(n)
    while work:
        n = work.pop()
        for p in pred[n]:
            if p in label:
                continue
            if label[n] == "L":
                label[p] = "W"; work.append(p)
            elif label[n] == "W":
                deg[p] -= 1
                if deg[p] == 0:
                    label[p] = "L"; work.append(p)
    for n in nodes:
        label.setdefault(n, "D")
    return label

def outcome(pos):
    lab = solve(pos)
    return lab[(canon(pos), 0)], lab[(canon(pos), 1)]  # (Left first, Right first)

def left_survives_second(pos):
    return solve(pos)[(canon(pos), 1)] != "W"  # Right (first) cannot force a win

def ge(g, h):
    return left_survives_second(list(g) + [NEG[c] for c in h])

def rel(g, h):
    a, b = ge(g, h), ge(h, g)
    return {(True, True): "=", (True, False): ">",
            (False, True): "<", (False, False): "||"}[(a, b)]

print("== survival comparisons (independent solver) ==")
for lhs, rhs in [(("over", "over"), ("over",)), (("under", "under"), ("under",)),
                 (("*", "over"), ("over",)), (("*", "under"), ("under",)),
                 (("over", "under"), ()), (("*", "*"), ()),
                 (("on", "over"), ("on",)), (("over",), ()), (("over",), ("*",))]:
    print(f"  {'+'.join(lhs):<13} vs {'+'.join(rhs) if rhs else '0':<6}: {rel(lhs, rhs)}")

print("\n== standalone outcomes (Left-first, Right-first; W/L/D from mover's view) ==")
for g in [(), ("over",), ("over", "over"), ("*", "over"), ("over", "under"), ("dud",)]:
    print(f"  o({'+'.join(g) if g else '0'}) = {outcome(g)}")

print("\n== full-outcome context battery ==")
ctxs = [(), ("*",), ("*2",), ("1",), ("-1",), ("up",), ("down",), ("on",),
        ("off",), ("over",), ("under",), ("dud",), ("up", "up"), ("1", "under"),
        ("-1", "over"), ("down", "down"), ("*", "under"), ("under", "under")]
pairs = [(("over", "over"), ("over",)), (("*", "over"), ("over",)),
         (("under", "under"), ("under",)), (("*", "under"), ("under",))]
bad = 0
for lhs, rhs in pairs:
    for X in ctxs:
        a, b = outcome(lhs + X), outcome(rhs + X)
        if a != b:
            bad += 1
            print(f"  MISMATCH {'+'.join(lhs)} vs {'+'.join(rhs)} @ X={X}: {a} vs {b}")
print(f"  mismatches: {bad} over {len(pairs) * len(ctxs)} comparisons")

print("\n== over+under vs 0: contexts where outcomes DIFFER (substitution check) ==")
diffs = []
for X in ctxs:
    a, b = outcome(("over", "under") + X), outcome(X)
    if a != b:
        diffs.append((X, a, b))
for X, a, b in diffs:
    print(f"  X={'+'.join(X) if X else '(empty)'}: over+under+X={a}  X={b}")
print(f"  differing contexts: {len(diffs)}")
