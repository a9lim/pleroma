# Open Problem 3: Ordinal Nim Multiplication Beyond the Verified Excess Table

This file records the June 2026 push on `OPEN.md` problem 3: derive, falsify, or
sharpen a closed formula for Lenstra excess in ordinal nim multiplication below
the first transcendental boundary `omega^(omega^omega)`.

The result is progress, not closure. We now have a sharper finite-field
reformulation, an independent local oracle, and a locally certified first new
carry `alpha_47`. The global closed formula is still unproved.

## Current State

Implemented in the Rust tower:

- DiMuro Table 1 rows through `alpha_43`.
- A locally verified row
  `alpha_47 = omega^(omega^7) + 1`.
- The operational boundary is now `alpha_53`: a carry needing `alpha_53` or
  beyond returns `None`.

Current external data, refreshed on 2026-06-09:

- OEIS A380496 has 1417 extended rows: 799 known and 618 unknown.
- The OEIS b-file has 126 initial known rows.
- The first OEIS unknown is row `n=127`, the 127th odd prime `p=719`.
- For `p=719`, `f(719) = ord_719(2) = 359` and `Q(359) = {359}`.
- The transfinite-nim-calculator logs record the direct component exponent as
  `e_719 = 1258230380`, which is the practical wall for direct exponentiation.

## Notation

For an odd prime `p`:

- `f(p) = ord_p(2)`.
- `Q(h)` is Lenstra's set of prime-power components appearing in `kappa_h`.
- The Lenstra excess `m_p` is the least finite `m` such that
  `kappa_{f(p)} + m` has no `p`-th root in the relevant finite component field.
- The Kummer carry is

```text
alpha_p = kappa_{f(p)} + m_p.
```

For the ordinal tower, a row with `Q(f(p)) = {q}` and finite excess `m` gives the
ordinal sum corresponding to `kappa_q + m`.

## Exact Reformulation

Let `beta = kappa_{f(p)} + m` lie in the finite component field `F_{2^E}`.
The multiplicative group is cyclic of order `N = 2^E - 1`.

The equation `x^p = beta` is solvable in `F_{2^E}` iff `beta` lies in the image
of the `p`-power map on this cyclic group. Since `p` is prime, this gives:

```text
beta has no p-th root  <=>  p divides ord(beta).
```

Equivalently, when `p | N`:

```text
beta has no p-th root  <=>  beta^((2^E - 1)/p) != 1.
```

Thus:

```text
m_p = least m such that p | ord(kappa_{f(p)} + m).
```

This is more useful than the original root-search phrasing because it turns the
finite correction into a statement about prime divisors of a specific
multiplicative order.

## Independent Oracle

`experiments/ordinal_excess_probe.py` is a small local term-algebra oracle. It is
not a replacement for CGSuite or the C++ calculator. It exists to verify the
first subtle cases without using the Rust production tower as an oracle.

It implements:

- The impartial term algebra used by the calculators.
- A multiplicative-order test for small component fields.
- A fixed-base exponentiation path, ported from the C++ calculator's strategy,
  for targeted root tests where full order factorization is unnecessary.

Current probe output includes:

```text
p=7,  m=0, Q=(3,),  root? True
p=7,  m=1, Q=(3,),  root? False
p=19, m=1, Q=(9,),  root? True
p=19, m=4, Q=(9,),  root? False
p=73, m=1, Q=(9,),  root? False
p=47, m=1, Q=(23,), root? False
```

The last line certifies `m_47 = 1` using only lower verified rows. Since
`f(47) = 23` and `Q(23) = {23}`,

```text
alpha_47 = kappa_23 + 1 = omega^(omega^7) + 1.
```

That value is now implemented in `src/scalar/big/ordinal/tower.rs`.

## Candidate Formula

The empirical rule that survived the current audit is:

```text
m_p = 0  if Q(f(p)) is not a singleton odd prime-power
m_p = 1  if Q(f(p)) is a singleton odd prime-power
         except:
m_p = 4  when f(p) = 2 * 3^k, k >= 1
```

This matched:

- all 950 calculator records with known `Q`-sets;
- all OEIS-known rows covered by those calculator `Q`-sets.

The only `m=4` records in the calculator logs are:

```text
p=19,   f(p)=18,  Q(f(p))={9}
p=163,  f(p)=162, Q(f(p))={81}
p=1459, f(p)=486, Q(f(p))={243}
```

This is evidence, not a theorem.

## What Is Already Falsified

`Q(f(p))` alone does not determine `m_p`.

Examples:

```text
Q={9}:   m_19=4,    m_73=1
Q={81}:  m_163=4,   m_2593=1
Q={243}: m_1459=4,  m_487=1
```

The order reformulation explains the first split:

```text
ord(kappa_9 + 1) = 3^3 * (2^9 - 1).
```

So `73 | ord(kappa_9 + 1)`, but `19` does not divide it. Adding `4` changes the
order and picks up `19`.

## Why the Candidate Is Still Not Proved

If the `0/1/4` rule were true, it would imply a global bound:

```text
m_p <= 4.
```

Lenstra explicitly left absolute boundedness open after giving lower-bound rules
such as:

- singleton odd `Q(f(p))` forces positive excess;
- `f(p) = 2 * 3^k` forces excess at least `4`.

So proving the candidate is not just table cleanup; it would settle a stronger
boundedness question.

## The p=719 Wall

For the first unknown row:

```text
p = 719
f(p) = 359
Q(359) = {359}
predicted m_719 = 1
```

The component chain is sparse:

```text
359 -> 179 -> 89 -> 11 -> 5 -> finite
```

The component field degree is:

```text
E = 2 * 2 * 5 * 11 * 89 * 179 * 359
  = 1258230380.
```

Directly testing

```text
(kappa_359 + 1)^((2^E - 1)/719)
```

is not locally feasible with the current term-array algorithm.

## Norm Reduction Direction

Since `f(p) = ord_p(2)`, `p | 2^f - 1`. If `beta in F_{2^E}` and `f | E`, then:

```text
beta^((2^E - 1)/p)
  = Norm_{F_{2^E}/F_{2^f}}(beta)^((2^f - 1)/p).
```

So the `p=719` test can be reduced to:

1. Compute
   `Norm_{F_{2^1258230380}/F_{2^359}}(kappa_359 + 1)`.
2. Test whether that norm has order divisible by `719` in `F_{2^359}`.

This is the next algorithmic target. It avoids the final huge target field, but
still requires a structural way to compute a norm over

```text
E / f = 3504820
```

Frobenius conjugates without materializing the giant term algebra.

## Files Updated

Core:

- `src/scalar/big/ordinal/tower.rs`
  - added `alpha_47`;
  - added `locally_verified_alpha_47_landmark`;
  - moved refusal boundary to `alpha_53`.
- `src/scalar/big/ordinal/nim.rs`
  - updated the documented boundary.
- `src/scalar/big/ordinal/mod.rs`
  - updated provenance and boundary docs.
- `src/games/nimber_game.rs`
  - updated the turning-corners boundary.

Docs and local guidance:

- `OPEN.md`
  - records the order criterion, candidate rule, `alpha_47`, and the `p=719`
    pressure point.
- `README.md`
  - updated the ordinal boundary.
- `src/scalar/AGENTS.md`
  - updated the scalar-pillar boundary note.

Experiment:

- `experiments/ordinal_excess_probe.py`
  - independent term-algebra probe;
  - fixed-base `p=47` root test;
  - documents the `Q={9}` split.

## Verification Run

Commands run successfully after the `alpha_47` promotion:

```sh
python3 -m py_compile experiments/ordinal_excess_probe.py
python3 experiments/ordinal_excess_probe.py
cargo fmt --check
cargo test
cargo check --all-targets
cargo check --features python --all-targets
cargo clippy --all-targets -- -D warnings
cargo clippy --features python --all-targets -- -D warnings
git diff --check
```

Focused Rust tests added/passing:

```text
scalar::big::ordinal::tower::tests::locally_verified_alpha_47_landmark
scalar::big::ordinal::tower::tests::boundary_returns_none_past_prime_47
```

## Sources Checked

- OEIS A380496: `https://oeis.org/A380496`
- Lenstra, "On the algebraic closure of two":
  `https://pub.math.leidenuniv.nl/~lenstrahw/PUBLICATIONS/1977e/art.pdf`
- CGSuite `NimFieldCalculator.scala`:
  `https://github.com/aaron-siegel/cgsuite`
- Django Peeters `transfinite-nim-calculator`:
  `https://github.com/DjangoPeeters/transfinite-nim-calculator`

No subagents, gaslamp, or Claude consultation were used.

## Next Concrete Steps

1. Implement a structural norm computation for singleton odd `Q={q}` cases.
2. Apply it to `p=719` and decide whether `m_719=1` is certified or falsified.
3. If the norm route works, test the next OEIS unknown singleton rows before
   promoting any more Rust carries.
4. Try to prove the special family:

```text
ord(kappa_{3^k} + 1) = 3^(k+1) * (2^(3^k) - 1)
```

or find the first failure. This is the clearest route to explaining the observed
`f(p)=2*3^k` exception.
