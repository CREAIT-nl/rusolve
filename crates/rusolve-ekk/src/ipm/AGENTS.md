# ipm

The interior point solvers: IPX (with BASICLU), HiPO (with its orderings and BLAS-level kernels), and the wrappers that `Highs::run` calls.

## What lives here
- `ipx/`: IPX LP solver, model, basis, KKT solvers, crossover, power method; `ipx_wrapper.rs`: `IpxWrapper.cpp`.
- `basiclu/`: the BASICLU LU kernel (`lu_*.rs`, `basiclu_*.rs`).
- `hipo/`: `factorhighs/`, `ipm/`, `orderings/` (AMD, METIS, RCM), `auxiliary/`, `blas.rs`; `hipo_wrapper/`: the HiPO entry point and report.
- `oracle_trace.rs` files under `ipx/` and `hipo/`: trace emitters for the oracle comparison.

## HiGHS files mirrored
`ipm/IpxWrapper.cpp`, `ipm/ipx/*`, `ipm/basiclu/*`, `ipm/hipo/**`, `extern/amd`, `extern/rcm`, `extern/metis`. Per-file mapping: `docs/agent/codemap.md`.

## Invariants
- Bit-exact with the HiGHS 1.16 oracle: same iteration trace, same ordering, same crossover pivots. No FMA; keep summation order, because the dense kernels are order sensitive.
- The MIP analytic centre and the LP IPM fallback run through this IPX (`src/mip/ipm_replay.rs`); changing numerics here can move identity.
- Default behavior must match HiGHS; a rusolve variant needs a default-off option.

## How to test
- IPX: `cargo test --release -p rusolve-ekk --test ipx_basiclu --test ipx_linalg --test ipx_model --test ipx_run`; fixtures in `tests/fixtures/ipx/` and `tests/fixtures/basiclu/`.
- HiPO: `--test t_w4_hipo --test t_w4_orderings`.
- Trace emitters: `RUSOLVE_IPX_TRACE`, `RUSOLVE_HIPO_TRACE`, `RUSOLVE_HIPO_PROF`.
- Benchmarks: `tools/bench/ipm_bench.py`.

## Change recipes
1. A trace diverges: the test reports the first record; find the C++ function from the doc line, compare reduction order and constants.
2. Port a missing orderings or BLAS piece: `docs/skills/contrib-port-function/SKILL.md`; oracle patches for these traces live in `tools/highs-oracle/highs-oracle-ipx-hooks.patch` and `highs-oracle-w4-hooks.patch`.
3. Add an IPM option: `docs/skills/contrib-add-option/SKILL.md`; read it in `hipo_wrapper/` or `ipx_wrapper.rs`.

## Pitfalls
- IPX results feed the MIP analytic centre: a small numerical change here can show up as a node-count change in identity, not as an LP failure.
- The IPM failure fallback (`lp_ipm_failure_policy`: cold IPX restart after a HiPO failure, crossover on a near-optimal IPX stall, simplex on IPM error) lives in `ipx_wrapper.rs`, `hipo_wrapper/` and `lp_data/highs_solve.rs`.
