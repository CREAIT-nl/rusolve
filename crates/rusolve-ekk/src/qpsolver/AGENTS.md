# qpsolver

The HiGHS convex QP solver (Quass): active-set method with Devex, Dantzig or steepest-edge pricing, a basis with Cholesky-style factor, ratio test, perturbation and scaling.

## What lives here
- `quass/` (`solve.rs`, `mod.rs`): the main loop; `a_quass.rs`, `a_asm.rs`: the `a_*` ports of `a_quass.cpp` and `a_asm.cpp`.
- `basis/`, `factor.rs`, `matrix.rs`, `instance.rs`, `qpvector.rs`, `gradient.rs`, `reducedcosts.rs`, `reducedgradient.rs`.
- Pricing: `pricing.rs`, `devexpricing.rs`, `dantzigpricing.rs`, `steepestedgepricing.rs`; `ratiotest.rs`, `perturbation.rs`, `scaling.rs`.
- `feasibility_bounded.rs`, `feasibility_highs.rs`, `settings.rs`, `runtime.rs`, `statistics.rs`, `eventhandler.rs`, `qpconst.rs`, `snippets.rs`, `oracle_trace.rs`.

## HiGHS files mirrored
`qpsolver/*.hpp` and `*.cpp` (`basis`, `quass`, `ratiotest`, `perturbation`, `a_asm`, `a_quass`, pricing headers). Table in `docs/agent/codemap.md`.

## Invariants
- Same iteration path as HiGHS: pricing choice, ratio test tie-breaks and perturbation draws must match; `perturbation.rs` notes that `std::default_random_engine` is libstdc++ `minstd_rand0`.
- File names stay lowercase like the C++ headers so the codemap and the doc lines line up.
- Convexity is decided before this solver runs (router in `crates/rusolve-cli/src/qp.rs`); nonconvex input must never reach it.
- Not in the default MIP search; identity unaffected but run the gate.

## How to test
- `cargo test --release -p rusolve-ekk --test t_w1_qp`; fixtures come from `tools/highs-oracle` (hooks patch `tools/highs-oracle/highs-oracle-w1-hooks.patch`).
- Trace: `RUSOLVE_QP_TRACE` (see `oracle_trace.rs`).
- Unit tests: `cargo test --release -p rusolve-ekk --lib qpsolver`; CLI: `crates/rusolve-cli/tests/qp.rs`.

## Change recipes
1. QP differs from HiGHS: run `t_w1_qp`, read the first diverging iteration, open the matching `.cpp`/`.hpp` function.
2. Add a QP option: `docs/skills/contrib-add-option/SKILL.md`; QP options such as `qp_regularization_value` are read in `src/lp_data/highs_qp.rs`; strategy enums are in `settings.rs`.
3. Port a missing piece: `docs/skills/contrib-port-function/SKILL.md`.

## Pitfalls
- File names are lowercase and match the C++ headers, which differs from the rest of the crate; do not rename them.
