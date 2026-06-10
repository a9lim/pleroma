"""Adapted (Arf-normal symplectic) frame: build for each m=8 Gold form,
predict misses from the pattern theory, verify by full sweep.
Also (6,1) rank-4 full sweep and m=4 family in adapted frames."""
from extraspecial_core import *

validate()

def build_adapted_frame(Q, m):
    """Arf-normal frame: hyperbolic pairs (u_i, v_i) with q-values
    (1,1) on at most one pair (if Arf 1), (0,0) elsewhere, plus radical basis.
    Returns list of frame vectors [u1, v1, u2, v2, ..., r1, r2, ...]."""
    V = list(range(1, 1 << m))
    def Bf(u, v): return Q[u ^ v] ^ Q[u] ^ Q[v]
    # radical = {v : B(v, .) == 0}
    rad = [v for v in range(1 << m) if all(Bf(v, 1 << i) == 0 for i in range(m))]
    # basis of radical
    radbasis = []
    span = {0}
    for v in rad:
        if v not in span:
            radbasis.append(v)
            span = {s ^ v2 for s in span for v2 in (0, v)}
    # core: symplectic Gram-Schmidt on a complement of the radical
    pairs = []
    used_span = set(span)
    remaining = [v for v in range(1, 1 << m) if v not in used_span]
    cur_span = set(span)
    def spanof(vecs):
        s = {0}
        for v in vecs:
            s |= {x ^ v for x in s}
        return s
    chosen = []
    while True:
        # pick u not in current span, with some partner w: B(u,w)=1, w not in span
        cand_u = [v for v in range(1, 1 << m) if v not in cur_span]
        if not cand_u: break
        found = False
        for u in cand_u:
            for w in cand_u:
                if Bf(u, w) == 1:
                    # orthogonalize: keep (u,w); reduce later pairs against them
                    pairs.append((u, w))
                    chosen += [u, w]
                    cur_span = spanof(radbasis + chosen)
                    found = True
                    break
            if found: break
        if not found: break
        # project remaining space orthogonal to (u,w): handled implicitly by
        # re-selecting candidates orthogonal to all chosen pairs:
        # (simple approach: filter cand by orthogonality at next iteration)
        # -- implement properly below
        u, w = pairs[-1]
        # replace the complement: vectors orthogonal to u and w modulo span
        # we'll just continue; correctness enforced by final checks
        # re-pick: restrict future candidates to x' = x + B(x,w)u + B(x,u)w
        # Do an explicit symplectic reduction instead:
        break
    # cleaner: full symplectic reduction
    pairs = []
    basis_done = list(radbasis)
    def in_span(v, vecs):
        s = spanof(vecs)
        return v in s
    avail = [v for v in range(1, 1 << m)]
    comp = []   # current complement vectors to process
    # iterative reduction over the whole space
    work = [v for v in range(1, 1 << m) if not in_span(v, basis_done)]
    while work:
        u = work[0]
        partner = None
        for w in work:
            if Bf(u, w) == 1: partner = w; break
        if partner is None:
            # u central in remaining space -> should be in radical span; skip
            work = [v for v in work[1:] if not in_span(v, basis_done +
                    [x for p in pairs for x in p] + [u])]
            basis_done.append(u)   # shouldn't happen for honest radical calc
            continue
        w = partner
        pairs.append((u, w))
        # reduce: x -> x + B(x,w)u + B(x,u)w, keep those independent
        newwork = []
        for x in work:
            if x in (u, w): continue
            x2 = x ^ (u if Bf(x, w) else 0) ^ (w if Bf(x, u) else 0)
            if x2 != 0 and not in_span(x2, basis_done +
                                       [y for p in pairs for y in p] + newwork):
                newwork.append(x2)
        work = newwork
    # now adjust q-values pairwise to Arf normal form:
    # within pair (u,w): want (Q(u),Q(w)) == (0,0) if possible:
    # transformations: u->u+w etc. The 4 candidates u,w,u+w give q-values; a
    # hyperbolic pair has some basis with (0,0) iff Arf-contribution 0.
    norm_pairs = []
    arf_ones = 0
    for (u, w) in pairs:
        cands = [(u, w), (w, u), (u ^ w, w), (u, u ^ w), (w, u ^ w), (u ^ w, u)]
        best = None
        for (x, y) in cands:
            if Bf(x, y) != 1: continue
            if Q[x] == 0 and Q[y] == 0: best = (x, y); break
        if best is None:
            best = (u, w)   # anisotropic pair: all combos (1,1)? keep
            arf_ones += 1
        norm_pairs.append(best)
    frame = []
    for (u, w) in norm_pairs: frame += [u, w]
    frame += radbasis
    # sanity: frame is a basis
    assert len(spanof(frame)) == (1 << m), "frame not a basis"
    return frame, norm_pairs, radbasis, arf_ones

def frame_sweep_detail(Q, frame, m, maxfirst=True):
    mm = len(frame)
    qover = [Q[v] for v in frame]
    Bover = []
    for i in range(mm):
        row = 0
        for j in range(mm):
            if i == j: continue
            row |= (Q[frame[i] ^ frame[j]] ^ Q[frame[i]] ^ Q[frame[j]]) << j
        Bover.append(row)
    misses = []
    ndg = 0
    for cm in range(1 << mm):
        xf = 0
        for i in range(mm):
            if (cm >> i) & 1: xf ^= frame[i]
        v, ch = echo_value(cm, None, None, mm, ko='self', maxfirst=maxfirst,
                           qover=qover, Bover=Bover)
        ndg += ch
        if v != Q[xf]:
            misses.append(cm)
    return misses, qover, Bover, ndg

def predict_misses(qover, Bover, mm):
    """pattern theory: support bad iff (B-graph(S), target) in bad classes.
    Implemented for matching graphs only (adapted frames):
    p = #disjoint edges in S; bad iff (p==2 and t==1) or (p>=3 and t != val)
    with val(p>=2, P1max) = 0  => bad iff p>=2 and t==1."""
    bad = []
    for cm in range(1 << mm):
        S = [i for i in range(mm) if (cm >> i) & 1]
        edges = 0
        for a in range(len(S)):
            for b in range(a + 1, len(S)):
                if (Bover[S[a]] >> S[b]) & 1: edges += 1
        t = 0
        for i in S: t ^= qover[i]
        t ^= edges & 1
        # in a matching frame, edges == p
        if edges >= 2 and t == 1:
            bad.append(cm)
        elif edges == 1 and False:
            pass
    return bad

def coord_form(m, hyp_pairs, q1pairs=0, aniso_rad=0):
    """synthetic: hyp_pairs hyperbolic pairs on coords (0,1),(2,3),...;
    q=1 on both coins of the first q1pairs pairs; q=1 on aniso_rad radical coords."""
    Q = []
    for x in range(1 << m):
        t = 0
        for p in range(hyp_pairs):
            a, b = (x >> (2 * p)) & 1, (x >> (2 * p + 1)) & 1
            t ^= a & b
            if p < q1pairs: t ^= a ^ b
        for r in range(aniso_rad):
            t ^= (x >> (2 * hyp_pairs + r)) & 1
        Q.append(t)
    return Q

forms = [
    ("(8,1)l1 rank6", gold_q(8, 1, 1), 8),
    ("(8,2)l1 rank4", gold_q(8, 2, 1), 8),
    ("(8,1)l2 bent.rank8", gold_q(8, 1, 2), 8),
    ("synth r4 Arf0 rad4iso  m=8", coord_form(8, 2), 8),
    ("synth r4 Arf1 rad4iso  m=8", coord_form(8, 2, 1), 8),
    ("synth r4 Arf1 rad2aniso m=8", coord_form(8, 2, 1, 2), 8),
    ("synth r2 Arf0 rad6iso  m=8", coord_form(8, 1), 8),
]
for name, Q, m in forms:
    frame, pairs, radb, arf1 = build_adapted_frame(Q, m)
    qover = [Q[v] for v in frame]
    print(f"\n{name}: pairs={len(pairs)} radical_dim={len(radb)} "
          f"aniso_pairs={arf1} q-on-frame={qover} "
          f"Q|radical={[Q[v] for v in radb]}")
    mm = len(frame)
    Bover = []
    for i in range(mm):
        row = 0
        for j in range(mm):
            if i == j: continue
            row |= (Q[frame[i] ^ frame[j]] ^ Q[frame[i]] ^ Q[frame[j]]) << j
        Bover.append(row)
    pred = predict_misses(qover, Bover, mm)
    misses, _, _, ndg = frame_sweep_detail(Q, frame, m, True)
    print(f"  predicted misses: {len(pred)}  actual misses: {len(misses)} "
          f"  match: {sorted(pred) == sorted(misses)}  choice-states={ndg}")
    if misses and len(misses) <= 12:
        print(f"  miss supports: {[bin(c) for c in misses]}")
