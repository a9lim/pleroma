import itertools, collections

# component games: name -> (Left options, Right options)
G = {
    'zero':  ([], []),
    'star':  (['zero'], ['zero']),
    'one':   (['zero'], []),
    'on':    (['on'], []),
    'off':   ([], ['off']),
    'over':  (['zero'], ['over']),
    'under': (['under'], ['zero']),
    'dud':   (['dud'], ['dud']),
}

def succ(state):
    comps, player = state
    out = []
    for i, c in enumerate(comps):
        opts = G[c][0] if player == 'L' else G[c][1]
        for o in opts:
            nc = tuple(sorted(comps[:i] + (o,) + comps[i+1:]))
            out.append((nc, 'R' if player == 'L' else 'L'))
    return out

def outcomes(start_comps):
    # enumerate reachable states
    seen = set()
    stack = [(tuple(sorted(start_comps)), 'L'), (tuple(sorted(start_comps)), 'R')]
    while stack:
        s = stack.pop()
        if s in seen: continue
        seen.add(s)
        stack.extend(succ(s))
    # retrograde Win/Loss/Draw for the player to move
    succs = {s: succ(s) for s in seen}
    preds = collections.defaultdict(list)
    cnt = {}
    for s, ss in succs.items():
        cnt[s] = len(ss)
        for t in ss:
            preds[t].append(s)
    label = {}
    from collections import deque
    q = deque()
    for s in seen:
        if cnt[s] == 0:
            label[s] = 'Loss'; q.append(s)
    while q:
        s = q.popleft()
        for p in preds[s]:
            if p in label: continue
            if label[s] == 'Loss':
                label[p] = 'Win'; q.append(p)
            else:  # successor is Win for opponent
                cnt[p] -= 1
                if cnt[p] == 0:
                    label[p] = 'Loss'; q.append(p)
    for s in seen:
        label.setdefault(s, 'Draw')
    k = tuple(sorted(start_comps))
    return label[(k,'L')], label[(k,'R')]

def show(comps):
    l, r = outcomes(comps)
    # translate to who-wins
    def who(lab, mover):
        if lab == 'Draw': return 'Draw'
        winner = mover if lab == 'Win' else ('L' if mover=='R' else 'R')
        return winner + ' wins'
    print(f"{'+'.join(comps):20s}  L to move: {who(l,'L'):8s}  R to move: {who(r,'R'):8s}")

print("-- baseline outcomes (sanity, vs LoopyValue::outcome) --")
for c in ['zero','star','on','off','over','under','dud','one']:
    show([c])
print("-- the disputed identity --")
show(['over','under'])         # code claims this == zero, i.e. P (mover loses both ways)
show(['zero'])
print("-- distinguish from dud via context +1 --")
show(['over','under','one'])
show(['dud','one'])
print("-- distinguish from 0 via context +star --")
show(['over','under','star'])
show(['star'])
print("-- check the other closed sums the code claims --")
show(['on','off'])             # claimed dud
show(['on','over']); show(['on','under']); show(['on','star'])   # claimed on
show(['off','over']); show(['off','under']); show(['off','star']) # claimed off
show(['star','star'])          # claimed zero
