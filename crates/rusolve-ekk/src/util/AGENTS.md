# util

The HiGHS utility layer: vectors, sparse matrices, the LU factor `HFactor`, hash tables, sets, sorting, random numbers, extended precision and C++ standard-library stand-ins.

## What lives here
- `hfactor/` (`build.rs`, `kernel.rs`, `refactor.rs`, `solve.rs`, `update.rs`, `extend.rs`, `debug.rs`), `hfactor_const.rs`, `factor_timer.rs`.
- `hvector.rs`, `hset.rs`, `sparse_matrix/` (`HighsSparseMatrix`), `highs_matrix_slice.rs`, `highs_matrix_utils.rs`.
- Hashing and trees: `highs_hash_table.rs`, `highs_hash_tree.rs`, `highs_hash.rs`, `highs_rb_tree.rs`, `highs_splay.rs`, `highs_disjoint_sets.rs`.
- Numerics and order: `highs_cdouble.rs`, `highs_random.rs`, `highs_sort.rs`, `highs_pdqsort.rs`, `highs_sort_set.rs`, `highs_integers.rs`, `highs_linear_sum_bounds.rs`.
- `cpp_std.rs` (`std::sort`, `mt19937`, priority queue stand-ins), `c_printf.rs`, `highs_timer.rs`, `highs_data_stack.rs`, `highs_utils.rs`.

## HiGHS files mirrored
`util/HFactor*.cpp`, `HVectorBase.cpp`, `HSet.cpp`, `HighsSparseMatrix.cpp`, `HighsHash.h`, `HighsHashTree.h`, `HighsRandom.h`, `HighsCDouble.h` and the other `util/*.h`. Table: `docs/agent/codemap.md`.

## Invariants
- Order is behavior: `HVector` index order, hash-table iteration order, `std::sort` tie order and `mt19937` output all affect pivots and branching. Match them exactly (`cpp_std.rs` exists for that; `RUSOLVE_STD_LIBSTDCXX` selects the libstdc++ variant).
- `HighsCDouble` arithmetic must keep its error-free transformations; no FMA, no reassociation.
- Hot kernels (`hfactor/`) must not change results; any speedup is proven bit-identical on the factor and solve tests first.

## How to test
- `cargo test --release -p rusolve-ekk --test t_hfactor_build --test t_hfactor_refactor --test t_hfactor_solve --test t_hfactor_update --test t_hvector --test t_sparse_matrix --test t_sort --test t_random --test t_cdouble --test t_gksort`.
- Fixtures: `tests/fixtures/factor/`, `tests/fixtures/unit/`, `tests/fixtures/random.jsonl`; C++ standard-library cases (`std::sort`, `mt19937`, hash tree) are programs with expected output in `tests/fixtures/mip/unit/`. `HighsSort` has no oracle fixture (see the header of `tests/t_sort.rs`).
- Identity after any change in a hot kernel.

## Change recipes
1. Speed a factor kernel: change `hfactor/kernel.rs`, rerun the `t_hfactor_*` tests (bit-exact), then identity, then a measured benchmark.
2. Port a missing utility: `docs/skills/contrib-port-function/SKILL.md`; add a unit golden under `tests/fixtures/unit/`.
3. Add a C++ stand-in: put it in `cpp_std.rs` with a test comparing against oracle output.

## Pitfalls
- The C++ standard library is not specified bit for bit; `std::sort` and hash behavior differ between libstdc++ and libc++, so `cpp_std.rs` carries the implementations the port uses (`RUSOLVE_STD_LIBSTDCXX` selects the libstdc++ variant).
