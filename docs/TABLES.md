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
| Finite Lenstra excess integers `m_u`, `u <= 47` | `src/scalar/big/ordinal/tower.rs::finite_excess` | Yes. | Values through `43` are from DiMuro's source table; `47` is locally certified. `alpha_u` itself is now assembled from `ord_u(2)`, `Q(f(u))`, and this finite integer. |
| Named binary-code generator matrices: Hamming `[7,4,3]`, extended Hamming `[8,4,4]`, the indecomposable Type II `[16,8,4]`, and extended Golay `[24,12,8]` | `src/forms/integral/codes.rs::{hamming_code,extended_hamming_code,type_ii_len16_code,extended_golay_generator_rows}` | Yes. | These are finite named representatives for the Construction A bridge. The split length-16 Type II code is derived from the extended Hamming code; the indecomposable length-16 generator is sourced from the Harada-Munemasa self-dual-code table; Golay is shared with the Leech construction. |
| `E_6`, `E_7`, `E_8` Dynkin edge lists | `src/forms/integral/root_lattices.rs::{e_6,e_7,e_8}` | Yes. | These are exceptional finite diagrams. They could be generated from branch specs, but that would still encode the same exceptional data. |
| Exceptional `E_6/E_7/E_8` automorphism orders | `src/forms/integral/lattice/::RootComponentKind::automorphism_order`, `src/forms/integral/root_lattices.rs::E8_WEYL_GROUP_ORDER` | Yes. | The infinite `A_n` and `D_n` families are formulaic; the exceptional `E` orders remain constants. `E8_WEYL_GROUP_ORDER` is exposed separately because the theta/Siegel-Weil bridge records it alongside the `D16+` order. |
| Clifford-invariant vs Hasse-Witt correction `delta(n mod 8, d)` | `src/forms/witt/brauer_rational.rs::clifford_correction` | Yes. | The `n mod 8` -> `{(-1,-1), (-1,d)}` correction between `c(q)` and `s(q)` (Bridge F), from Lam GSM 67 pp. 117-119 and cross-checked against SageMath's `clifford_invariant`. The eight residue cases have no simpler closed form. |
| Finite loopy-value catalogue (`0`, `*`, `on`, `off`, `over`, `under`, `±`, `tis`, `tisn`, `dud`, plus integer `s&t` tags) | `src/games/loopy/::LoopyValue` methods | Yes. | The named finite catalogue and onside/offside tag surface are the intended public boundary; full loopy equality remains outside this table. |
| Python finite odd-field dispatch table | `src/py/forms.rs::finite_field_order`, `with_finite_odd_metric`, `with_finite_odd_metrics`, `with_finite_odd_value` | Yes for now. | Rust must monomorphise concrete const-generic types; replacing this needs a dynamic finite-field type or a generated support macro, not a numeric formula. |
| Python prime-field dispatch table | `src/py/forms.rs::with_prime_field`, `is_sum_of_n_squares` | Yes for now. | A formula such as "all primes" would not instantiate Rust types. |
| Python coin-family string aliases | `src/py/games.rs::parse_coin_family` | Yes. | API vocabulary. |
