# Production Hardcoded Tables

This file records production hardcoded tables and finite case tables found in the
runtime library and Python bindings. It excludes tests, examples, and experiment
oracles, and it also excludes trivial enum display strings unless the mapping is a
semantic catalogue.

## Remaining Data Tables

These are still real production tables because the finite data is curated,
sourced, or public API vocabulary rather than a theorem with a simpler closed
form.

| table | source | should stay a table? | note |
|---|---|---|---|
| Prime factors of `2^128 - 1` for nimber multiplicative orders | `src/scalar/finite_field/nimber/galois.rs::ORDER_FACTORS` | Yes. | The coarse identity is `2^128 - 1 = prod_{i=0..6} (2^(2^i)+1)`, but the prime factors of `F_5` and `F_6` are still recorded arithmetic data. |
| Finite Lenstra excess integers `m_u`, odd primes `3..=709` | `src/scalar/big/ordinal/tower.rs::finite_excess` | Yes. | OEIS A380496 ("Lenstra excess of the n-th odd prime"), the b-file's 126 known rows (the first 14 reproduce DiMuro Table 1 + the old `m_47`; first OEIS-unknown row is `p=719`). Indexed by odd-prime place; pinned against the b-file in `excess_table_matches_oeis_a380496`. `alpha_u` is assembled from `ord_u(2)`, `Q(f(u))`, and this finite integer. Provenance: Conway/Lenstra/Le Bruyn/Siegel/Peeters via CGSuite's calculator. |
| Named binary-code generator matrices: Hamming `[7,4,3]`, extended Hamming `[8,4,4]`, the indecomposable Type II `[16,8,4]`, and extended Golay `[24,12,8]` | `src/forms/integral/codes.rs::{hamming_code,extended_hamming_code,type_ii_len16_code,extended_golay_generator_rows}` | Yes. | These are finite named representatives for the Construction A bridge. The split length-16 Type II code is derived from the extended Hamming code; the indecomposable length-16 generator is sourced from the Harada-Munemasa self-dual-code table; Golay is shared with the Leech construction. |
| `E_6`, `E_7`, `E_8` Dynkin edge lists | `src/forms/integral/root_lattices.rs::{e_6,e_7,e_8}` | Yes. | These are exceptional finite diagrams. They could be generated from branch specs, but that would still encode the same exceptional data. |
| Exceptional automorphism-order constants: `E_6/E_7/E_8` orders, the `E_8` Weyl-group order, the Conway-group `Co_0` order (Leech / rootless Niemeier class), and the `D16+` order | `src/forms/integral/lattice/::RootComponentKind::automorphism_order`, `src/forms/integral/root_lattices.rs::E8_WEYL_GROUP_ORDER`, `src/forms/integral/mass_formula.rs::LEECH_AUT_ORDER`, `src/forms/integral/codes.rs::D16_PLUS_AUT_ORDER` | Yes. | The infinite `A_n`/`D_n` families are formulaic; the exceptional orders are curated constants. `E8_WEYL_GROUP_ORDER` and `D16_PLUS_AUT_ORDER = 2^15·16!` anchor the rank-8/rank-16 theta/Siegel-Weil bridge; `LEECH_AUT_ORDER = Co_0 = 2^22·3^9·5^4·7^2·11·13·23` is returned by `Niemeier::automorphism_group_order` for the rootless class. `LEECH_AUT_ORDER` and `D16_PLUS_AUT_ORDER` are also exported as Python module constants. |
| Exceptional Coxeter numbers `h(E_6)=12`, `h(E_7)=18`, `h(E_8)=30` | `src/forms/integral/niemeier.rs::NiemeierComponentKind::coxeter_number` | Yes. | The `A_n` (`n+1`) and `D_n` (`2n-2`) cases are formulaic; the exceptional `E` Coxeter numbers are constants. Used with rank to count roots per Niemeier class for the theta-series weighting. |
| Niemeier root-system, glue-index, and `Aut(N)/W(R)` catalogue | `src/forms/integral/niemeier.rs::NIEMEIER_CLASSES` | Yes. | This is the 24-class rank-24 even-unimodular catalogue from Conway-Sloane's Niemeier table, cross-checked by the glue-square determinant, mass sum, and weight-12 Siegel-Weil identity. The code builds root sublattices and the explicit Leech lattice; it does not encode 23 full glued Gram matrices. |
| Clifford-invariant vs Hasse-Witt correction `delta(n mod 8, d)` | `src/forms/witt/brauer_rational.rs::clifford_correction` | Yes. | The `n mod 8` -> `{(-1,-1), (-1,d)}` correction between `c(q)` and `s(q)` (Bridge F), from Lam GSM 67 pp. 117-119 and cross-checked against SageMath's `clifford_invariant`. The eight residue cases have no simpler closed form. |
| Real Clifford 8-fold (Bott) classification table `s = (q-p) mod 8 -> R/C/H` | `src/forms/char0.rs::real_core` | Yes. | The mod-8 period of the real Clifford classification (`0,6,7 -> R`; `1,5 -> C`; `2..=4 -> H`). Load-bearing for `classify_real`; the char-0 mirror of the `clifford_correction` `delta(n mod 8)` table above, with no simpler closed form than the eight cases. |
| Finite loopy-value catalogue (`0`, `*`, `on`, `off`, `over`, `under`, `±`, `tis`, `tisn`, `dud`, plus integer `s&t` tags) | `src/games/loopy/::LoopyValue` methods | Yes. | The named finite catalogue and onside/offside tag surface are the intended public boundary; full loopy equality remains outside this table. |
| Python finite odd-field dispatch table | `src/py/forms.rs::finite_field_order`, `with_finite_odd_metric`, `with_finite_odd_metrics`, `with_finite_odd_value` | Yes for now. | Rust must monomorphise concrete const-generic types; replacing this needs a dynamic finite-field type or a generated support macro, not a numeric formula. |
| Python prime-field dispatch table | `src/py/forms.rs::with_prime_field`, `is_sum_of_n_squares` | Yes for now. | A formula such as "all primes" would not instantiate Rust types. |
| Python char-2 finite-field dispatch table | `src/py/forms.rs::{with_finite_char2_field, with_finite_char2_metric, with_finite_char2_metrics}` | Yes for now. | Degree in `{1,2,3,4}` dispatch to `Fpn<2,N>` monomorphs; the char-2 companion to the finite odd-field dispatch. Same monomorphisation constraint. |
| Python local-field Springer dispatch tables | `src/py/forms.rs::{springer_decompose_qp, springer_decompose_qq, springer_decompose_laurent, springer_decompose_ramified_qp4_e2, springer_decompose_ramified_qp4_e3, springer_decompose_local}` | Yes for now. | Finite `(p)`, `(p, residue_degree)`, `(p, degree)`, and local-algebra-type dispatch over the supported `Qp/Qq/Laurent/Ramified` monomorphs. Rust must instantiate concrete const-generic field types, so the supported cells are an explicit table, not a formula. |
| Python coin-family string aliases | `src/py/games.rs::parse_coin_family` | Yes. | API vocabulary. |

## Out of scope (deliberately not tables)

Recorded so a future sweep does not re-flag them as gaps:

- **Eisenstein normalization constants** (`240`, `-504`, `65520/691`, ...) are computed
  at runtime from the single Bernoulli source (`forms/integral/mass_formula.rs::bernoulli`);
  the literals live only in `forms/integral/modular.rs` tests as pinned oracles. Keeping
  them derived rather than tabled is a deliberate discipline (see the comment in
  `eisenstein_constants_derive_from_the_shared_bernoulli_source`).
- **Local Hilbert-symbol factors** `eps(u) = (u-1)/2 mod 2` and `omega(u) = (u^2-1)/8 mod 2`
  (`forms/local_global/padic.rs::{eps2, omega2}`) are closed-form number-theoretic
  functions, not curated data.
- **The ogham language surface** — the world catalogue, builtin-function names, and reserved
  keywords (`src/ogham/{eval,parse,lex}.rs`) — is public API vocabulary but is owned by the
  language spec `docs/ogham/ogham.md`, not this inventory.
- **`clifford/` and `linalg/`** carry no curated lookup tables: signs go through `Scalar::neg`
  and blade products / reductions are computed index arithmetic.
