# model

The HiGHS model container: `HighsModel` (an LP plus an optional Hessian) and `HighsHessian` with its normalisation utilities.

## What lives here
- `highs_model.rs`: `HighsModel` (`model/HighsModel.h`, `.cpp`).
- `highs_hessian.rs`: `HighsHessian` (CSC Hessian, `kTriangular` format by default).
- `highs_hessian_utils/`: `HighsHessianUtils.cpp` (`normalise.rs` and the rest of the utilities).
- `mod.rs`: module doc and re-exports.

## HiGHS files mirrored
`model/HighsModel.h`, `model/HighsModel.cpp`, `model/HighsHessian.h`, `model/HighsHessian.cpp`, `model/HighsHessianUtils.cpp`. Table in `docs/agent/codemap.md`.

## Invariants
- Hessians are stored in HiGHS's triangular CSC format (`Default` is the cleared state, `dim_ = 0`); assessment and normalisation produce the same ordering and the same error codes.
- Same names as C++ (`HighsModel`, `HighsHessian`); no FMA.
- `crates/rusolve-model` builds a `HighsLp`; this directory is the container the QP path consumes, not the builder.
- No effect on the MIP search; run the gate anyway.

## How to test
- Unit tests live in the files: `cargo test --release -p rusolve-ekk --lib model`.
- Through the driver: `tests/t_w1_qp.rs` (QP models), and `crates/rusolve-cli/tests/qp.rs`.

## Change recipes
1. Port a missing Hessian utility: `docs/skills/contrib-port-function/SKILL.md`; put it in `highs_hessian_utils/` and give it a doc line naming the C++ function.
2. Add a model field: add it to `HighsModel`, update every clone and `Default`, and the `pass_model` paths in `src/lp_data/highs_pass_model_raw.rs`.
3. Convert from the Python builder: change `crates/rusolve-model/src/lib.rs`, not this directory.

## Pitfalls
- Hessian handling is shared with `src/qpsolver` and `src/lp_data/highs_qp.rs`; a change in normalisation changes the QP input.
