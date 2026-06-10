"""Play-semantics probe: the Arf invariant as a win-bias.

Dickson (1901): the Arf invariant of a quadratic form over F₂ is the value the
form takes *most often*. Quantitatively, for a nonsingular form on F₂^{2m'},

    #{v : Q(v) = 0} = 2^{2m'-1} + (-1)^{Arf} · 2^{m'-1}.

So IF a game had P-positions (second-player wins) exactly {v : Q(v) = 0}, then
the Arf invariant would be the *sign of the win-bias*: Arf=0 ⇒ the second player
wins from 2^{m'-1} more starting positions than the first; Arf=1 ⇒ the reverse.
The margin is fixed (a Gauss-sum); the Arf invariant chooses its direction.

This script checks that signature directly on the game-built Gold forms
Q_a(v) = Tr(v ⊗ v^{2^a}): brute-force the value distribution and confirm the
zero-count matches the count predicted from the form's Arf classification
(handling the radical). Restricted to nim-subfields F_{2^{2^k}} (the only m for
which the nimbers below 2^m are field-closed).

This does NOT exhibit a game whose P-positions are {Q=0}; disjunctive sums have
XOR-linear (not quadratic) outcomes, so any such game must be interactive, and
whether a natural one exists is the open problem. What it does pin down: the
Arf invariant *is* the win-bias, exactly, in the counting sense.
"""

import ogdoad as pl

from common import gold


def predicted_zeros(arf: int, rank: int, radical_dim: int, radical_aniso: bool, m: int) -> int:
    if radical_aniso:
        return 2 ** (m - 1)               # nonzero linear part balances Q exactly
    bias = 2 ** (m - rank // 2 - 1)
    return 2 ** (m - 1) + (-1) ** arf * bias


if __name__ == "__main__":
    hdr = f"{'field':>7} {'a':>2} {'Arf':>3} {'rank':>4} {'zeros':>7} {'predicted':>9} {'bias':>7} {'ok?':>3}"
    print(hdr)
    print("-" * len(hdr))
    all_ok = True
    for k in range(1, 5):                 # F2^2, F2^4, F2^8, F2^16
        m = 1 << k
        for a in range(1, k + 1):
            # classification, via the shipped Arf classifier
            q = [pl.Nimber(gold_basis) for gold_basis in
                 (gold(1 << i, a, m) for i in range(m))]
            b = {}
            for i in range(m):
                for j in range(i + 1, m):
                    # polar form B(e_i,e_j) = Q(e_i+e_j)+Q(e_i)+Q(e_j)
                    bij = gold((1 << i) ^ (1 << j), a, m) ^ gold(1 << i, a, m) ^ gold(1 << j, a, m)
                    if bij:
                        b[(i, j)] = pl.Nimber(1)
            res = pl.arf_nimber(pl.NimberAlgebra(q=q, b=b))
            # brute-force the value distribution
            zeros = sum(1 for v in range(1 << m) if gold(v, a, m) == 0)
            pred = predicted_zeros(res.arf, res.rank, res.radical_dim,
                                   res.radical_anisotropic, m)
            bias = zeros - 2 ** (m - 1)
            ok = zeros == pred
            all_ok = all_ok and ok
            print(f"  F2^{m:<3} {a:>2} {res.arf:>3} {res.rank:>4} {zeros:>7} {pred:>9} "
                  f"{bias:>+7} {('✓' if ok else '✗'):>3}")
    print("-" * len(hdr))
    print(f"Arf invariant predicts the exact win-bias (zero-count) of the game-built form: {all_ok}")
