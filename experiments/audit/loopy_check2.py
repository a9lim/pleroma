import itertools
def sum_game(g1, g2):
    pos = {}
    for p1, p2 in itertools.product(g1, g2):
        l = [(q,p2) for q in g1[p1][0]] + [(p1,q) for q in g2[p2][0]]
        r = [(q,p2) for q in g1[p1][1]] + [(p1,q) for q in g2[p2][1]]
        pos[(p1,p2)] = (l, r)
    return pos
def left_survives_second(g, start):
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

# up = {0|*}, down = -up = {*|0}, star = {0|0}
up    = {'up': (['0'], ['*']), '*': (['0'], ['0']), '0': ([], [])}
down  = {'down': (['*'], ['0']), '*': (['0'], ['0']), '0': ([], [])}
star  = {'*': (['0'], ['0']), '0': ([], [])}
zero  = {'0': ([], [])}
over  = {'over': (['0'], ['over']), '0': ([], [])}
under = {'under': (['under'], ['0']), '0': ([], [])}

def ge(gA, sA, gB_neg, sB_neg):
    return left_survives_second(sum_game(gA, gB_neg), (sA, sB_neg))

# classical sanity checks on finite games:
print("up >= star :", ge(up,'up', star,'*'))      # expect False (up || *)
print("star >= up :", ge(star,'*', down,'down'))  # expect False
print("up >= 0    :", ge(up,'up', zero,'0'))      # expect True (up > 0)
print("0 >= up    :", ge(zero,'0', down,'down'))  # expect False
# and the loopy facts under audit:
print("over >= star :", ge(over,'over', star,'*'))     # claim: True
print("star >= over :", ge(star,'*', under,'under'))   # claim: False
# bonus: over vs up (standard: over > up)
print("over >= up :", ge(over,'over', down,'down'))
print("up >= over :", ge(up,'up', under,'under'))
