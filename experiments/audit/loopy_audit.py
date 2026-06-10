#!/usr/bin/env python3
"""Audit: do over+over, under+under, *+over, *+under stay in the loopy stopper
catalogue (i.e. equal over / under / over / under)?

Components (each a letter); a position is a sorted tuple (multiset).
Moves return either a replacement letter or None (component removed -> 0).

  O = over  = {0 | over}   : L -> 0,    R -> O
  U = under = {under | 0}  : L -> U,    R -> 0
  S = *     = {0 | 0}      : L -> 0,    R -> 0
  N = on    = {on |}       : L -> N,    R: none
  F = off   = {| off}      : L: none,   R -> F
  P = up    = {0 | *}      : L -> 0,    R -> S
  D = down  = {* | 0}      : L -> S,    R -> 0
  1 = {0 |}                : L -> 0,    R: none
  M = -1 = {| 0}           : L: none,   R -> 0
"""
from functools import lru_cache
import itertools

LMOVES = {'O': [None], 'U': ['U'], 'S': [None], 'N': ['N'], 'F': [],
          'P': [None], 'D': ['S'], '1': [None], 'M': []}
RMOVES = {'O': ['O'], 'U': [None], 'S': [None], 'N': [], 'F': ['F'],
          'P': ['S'], 'D': [None], '1': [], 'M': [None]}

def moves(pos, player):
    tbl = LMOVES if player == 'L' else RMOVES
    out = set()
    for i, c in enumerate(pos):
        if i > 0 and pos[i-1] == c:
            continue  # multiset: identical components give identical moves
        rest = pos[:i] + pos[i+1:]
        for r in tbl[c]:
            out.add(tuple(sorted(rest + (r,))) if r is not None else rest)
    return out

def closure(start):
    seen, todo = {start}, [start]
    while todo:
        p = todo.pop()
        for pl in 'LR':
            for q in moves(p, pl):
                if q not in seen:
                    seen.add(q); todo.append(q)
    return seen

def forced_win(start, winner):
    """Positions (p, mover) from which `winner` can force the opponent stuck.
    Least fixed point; infinite play is NOT a win."""
    space = closure(start)
    states = [(p, m) for p in space for m in 'LR']
    win = {s: False for s in states}
    changed = True
    while changed:
        changed = False
        for (p, m) in states:
            if win[(p, m)]:
                continue
            opts = moves(p, m)
            if m == winner:
                v = any(win[(q, 'L' if m == 'R' else 'R')] for q in opts)
            else:
                # opponent to move: stuck => winner wins; else all moves lose
                v = (not opts) or all(win[(q, 'L' if m == 'R' else 'R')] for q in opts)
            if v:
                win[(p, m)] = True
                changed = True
    return win

def left_survives_second(pos):
    """Left survives pos with Right to move <=> Right cannot force Left stuck."""
    w = forced_win(pos, 'R')
    return not w[(pos, 'R')]

def neg(pos):
    m = {'O':'U','U':'O','S':'S','N':'F','F':'N','P':'D','D':'P','1':'M','M':'1'}
    return tuple(sorted(m[c] for c in pos))

def geq(g, h):
    """g >= h via the stopper survival criterion: Left survives g + (-h) second."""
    return left_survives_second(tuple(sorted(g + neg(h))))

def cmp_games(g, h):
    a, b = geq(g, h), geq(h, g)
    return {(True,True):'=', (True,False):'>', (False,True):'<', (False,False):'||'}[(a,b)]

print("== survival-criterion comparisons ==")
print("over+over   vs over :", cmp_games(('O','O'), ('O',)))
print("under+under vs under:", cmp_games(('U','U'), ('U',)))
print("*+over      vs over :", cmp_games(('S','O'), ('O',)))
print("*+under     vs under:", cmp_games(('S','U'), ('U',)))
# sanity: identities the code DOES implement
print("-- sanity (code's closed cases) --")
print("over+under  vs 0    :", cmp_games(('O','U'), ()))
print("*+*         vs 0    :", cmp_games(('S','S'), ()))
print("on+on       vs on   :", cmp_games(('N','N'), ('N',)))
print("on+*        vs on   :", cmp_games(('N','S'), ('N',)))
print("on+over     vs on   :", cmp_games(('N','O'), ('N',)))
print("over        vs 0    :", cmp_games(('O',), ()))   # expect >
print("over vs *           :", cmp_games(('O',), ('S',)))

def outcome(pos):
    wl = forced_win(pos, 'L')
    wr = forced_win(pos, 'R')
    lf = wl[(pos,'L')]; ls = wl[(pos,'R')]   # Left wins moving first / second
    rf = wr[(pos,'R')]; rs = wr[(pos,'L')]   # Right wins moving first / second
    return (lf, ls, rf, rs)

print("\n== outcome cross-check over 16 contexts ==")
ctxs = [(), ('S',), ('P',), ('D',), ('N',), ('F',), ('1',), ('M',),
        ('O',), ('U',), ('S','S'), ('P','P'), ('N','S'), ('1','U'), ('M','O'), ('D','D')]
ok = True
for X in ctxs:
    for (lhs, rhs, nm) in [(('O','O'), ('O',), 'over+over vs over'),
                           (('S','O'), ('O',), '*+over vs over'),
                           (('U','U'), ('U',), 'under+under vs under'),
                           (('S','U'), ('U',), '*+under vs under')]:
        a = outcome(tuple(sorted(lhs + X)))
        b = outcome(tuple(sorted(rhs + X)))
        if a != b:
            ok = False
            print(f"MISMATCH {nm} in context {X}: {a} vs {b}")
print("all context outcomes agree" if ok else "CONTEXT MISMATCH FOUND")

# stopper check: every alternating play from over+over / *+over terminates?
def is_stopper(pos):
    # an infinite alternating play exists iff there's a cycle in the (p,mover)
    # graph reachable with both players always able to move; detect via
    # "can play continue forever": greatest fixed point of "has a move to a
    # state where play can continue forever"
    space = closure(pos)
    states = [(p, m) for p in space for m in 'LR']
    inf = {s: True for s in states}
    changed = True
    while changed:
        changed = False
        for (p, m) in states:
            if not inf[(p, m)]:
                continue
            opts = moves(p, m)
            if not any(inf[(q, 'L' if m == 'R' else 'R')] for q in opts):
                inf[(p, m)] = False
                changed = True
    return not (inf[(pos,'L')] or inf[(pos,'R')])

print("\n== stopper checks (criterion applicability) ==")
for nm, p in [('over', ('O',)), ('over+over', ('O','O')), ('*+over', ('S','O')),
              ('over+under', ('O','U')), ('on+off', ('N','F'))]:
    print(f"{nm}: stopper = {is_stopper(p)}")
