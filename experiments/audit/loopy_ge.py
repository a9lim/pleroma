# Stopper comparison via survival (Siegel CGT, Thm on stoppers):
# G >= H  iff  Left survives G + (-H) playing second (Right moves first).
# "Survives" = play is infinite, or Right eventually cannot move.
# Model: each atomic game = dict pos -> (left_moves, right_moves).

import itertools

def atom(name):
    if name == '0':     return {'0': ([], [])}
    if name == 'star':  return {'*': (['0'], ['0']), '0': ([], [])}
    if name == 'over':  return {'over': (['0'], ['over']), '0': ([], [])}
    if name == 'under': return {'under': (['under'], ['0']), '0': ([], [])}
    if name == 'on':    return {'on': (['on'], []), '0': ([], [])}
    if name == 'off':   return {'off': ([], ['off']), '0': ([], [])}
    raise ValueError(name)

NEG = {'0':'0','star':'star','over':'under','under':'over','on':'off','off':'on'}
START = {'0':'0','star':'*','over':'over','under':'under','on':'on','off':'off'}

def sum_game(g1, g2):
    pos = {}
    for p1, p2 in itertools.product(g1, g2):
        l = [(q,p2) for q in g1[p1][0]] + [(p1,q) for q in g2[p2][0]]
        r = [(q,p2) for q in g1[p1][1]] + [(p1,q) for q in g2[p2][1]]
        pos[(p1,p2)] = (l, r)
    return pos

def left_survives_second(g, start):
    # attractor for Right to states (p, 'L') where Left has no move
    states = [(p,t) for p in g for t in 'LR']
    A = set((p,'L') for p in g if not g[p][0])
    changed = True
    while changed:
        changed = False
        for p in g:
            if (p,'R') not in A and any((q,'L') in A for q in g[p][1]):
                A.add((p,'R')); changed = True
            if (p,'L') not in A and g[p][0] and all((q,'R') in A for q in g[p][0]):
                A.add((p,'L')); changed = True
    return (start,'R') not in A

def ge(a, b):
    g = sum_game(atom(a), atom(NEG[b]))
    return left_survives_second(g, (START[a], START[NEG[b]]))

names = ['off','under','0','star','over','on']
print("G \\ H  : G>=H table")
for a in names:
    print(f"{a:>6}:", {b: ge(a,b) for b in names})
print()
def rel(a,b):
    x, y = ge(a,b), ge(b,a)
    return '=' if x and y else '>' if x else '<' if y else '||'
for a,b in [('over','star'),('star','under'),('star','0'),('over','0'),
            ('under','0'),('on','star'),('off','star'),('over','under')]:
    print(f"{a} {rel(a,b)} {b}")
