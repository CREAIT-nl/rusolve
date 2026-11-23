# highs-oracle fixture format

Every fixture is **JSON Lines**: one JSON object per line, UTF-8, `\n`
terminated, no trailing whitespace. `serde_json::from_str::<serde_json::Value>`
parses each line. Field order inside an object is fixed by the generator (the
files are byte-deterministic) but readers must not rely on it.

## Scalars

| C++ type | JSON encoding |
|---|---|
| `double` | string of **16 lowercase hex digits** = the IEEE-754 binary64 bit pattern, big-endian digit order (`f64::to_bits` printed with `{:016x}`). Decode with `f64::from_bits(u64::from_str_radix(s, 16)?)`. `-0.0` is `"8000000000000000"`, `+inf` is `"7ff0000000000000"`, `kHighsZero` (1e-50) is `"358dee7a4ad4b81f"`. |
| `HighsInt` (int32 in this build), `int8_t`, enum, `int` | JSON number |
| `bool` | JSON number `0` / `1` |
| hash | string of 16 lowercase hex digits (u64), see *Hashes* |
| string | JSON string |

Vectors are JSON arrays of the element encoding (`vf` = array of hex strings,
`vi` / `vi8` = array of numbers).

## Hashes

Fields named `h_*` (and every field of a `hash` object) are **FNV-1a 64**
(offset `0xcbf29ce484222325`, prime `0x100000001b3`) over the little-endian
bytes of the array: each `f64` as the 8 bytes of `to_bits().to_le_bytes()`,
each `HighsInt` as 4 bytes (`i32::to_le_bytes`), each `int8_t` as 1 byte. The
empty array hashes to the offset basis. Hashes exist to localise the first
diverging iteration cheaply; full arrays are dumped at solve boundaries.

## Composite objects

**HVector** (`hv`):

```json
{"size":27,"count":3,"index":[13,0,9],"synthetic_tick":"0000000000000000",
 "nz_i":[0,9,13],"nz_v":["3ff0000000000000","…","…"],
 "packFlag":1,"packCount":3,"packIndex":[…],"packValue":[…]}
```

- `index` = `index[0..count)` in HiGHS' exact order (order matters: HiGHS loops
  over it). Absent when `count < 0` (dense mode, indices unknown).
- `nz_i` / `nz_v` = every position `i` of `array` whose **bit pattern** is
  non-zero, ascending. So `-0.0` and `1e-50` entries appear; all other
  positions are `+0.0`. `array.len() == size`.
- `synthetic_tick` is part of the state: HEkk sums vector ticks to decide
  reinversion (`kRebuildReasonSyntheticClockSaysInvert`).
- `pack*` fields only where the record says `with pack` (FTRAN/BTRAN results
  used by `HFactor::update`). `packIndex`/`packValue` are `[0..packCount)`.

**HighsSparseMatrix** (`mat`): `format` (1 colwise, 2 rowwise, 3 rowwise
partitioned), `num_col`, `num_row`, `start`, `p_end` (partitioned only),
`index`, `value`.

**HighsLp** (`lp` in `ekk_solve_enter`): `num_col`, `num_row`, `sense`
(1 min, -1 max), `offset`, `col_cost`, `col_lower`, `col_upper`, `row_lower`,
`row_upper`, `a` (matrix), `is_scaled`, `is_moved`, `scale` {`strategy`,
`has_scaling`, `num_col`, `num_row`, `cost`, `col`, `row`}.

**basis**: `basicIndex` (variables `0..num_col` are structurals,
`num_col + i` is the logical of row `i`), `nonbasicFlag`, `nonbasicMove`.

## Record types (field `t`)

### `random.jsonl` (HighsRandom)

`random`: `seed`, then the outputs of one generator drawn **in this order**:
`integer` (16 × `integer()`), `sup` / `integer_sup` (64 × `integer(sup)`),
`fraction` (16 × `fraction()`), `closedFraction` (16), `real_m3_5`
(16 × `real(-3, 5)`), `bit` (16), `shuffle20` (`shuffle` of `0..20`).

### `unit/<lp>.jsonl` (HighsSparseMatrix, HVector)

Inputs are the unscaled column-wise matrix of the LP as read from MPS.

| `t` | fields |
|---|---|
| `lp` | `name`, `num_col`, `num_row`, `a` |
| `createRowwise` | `out` (from `a`) |
| `createRowwisePartitioned` | `null_partition` (1 = `in_partition == NULL`), `in_partition` (per column), `out` |
| `priceByColumn` | `in_kind`, `quad`, `in`, `out` — `a.priceByColumn(quad, out, in)` into a cleared `out` |
| `priceByRow` | `matrix` (`rowwise` / `partitioned`), `in_kind`, `quad`, `in`, `out` |
| `priceByRowWithSwitch` | `matrix`, `in_kind`, `quad`, `expected_density`, `from_index`, `switch_density`, `in`, `out` |
| `computeDot` | `in_kind`, `in`, `use_col` (list), `out` (list, same order) |
| `collectAj` | `size`, `steps`: list of {`col`, `multiplier`, `after`} applied in sequence to one vector that starts cleared |
| `updatePartitioned` | `in_partition`, `before` (partitioned matrix), `steps`: list of {`var_in`, `var_out`, `after`} applied in sequence (`update(var_in, var_out, a)`) |
| `hvector_tight` / `hvector_reIndex` / `hvector_norm2` | `in`, `out` |
| `hvector_pack` | `in`, `out` (with pack; `packFlag` set before `pack()`) |
| `hvector_saxpy` | `pivotX`, `in`, `pivot`, `out` = `in.saxpy(pivotX, &pivot)` |

The partitioned matrix's `in_partition[j] = 1` means column `j` lies in
`[start, p_end)` — the role `nonbasicFlag` plays in the simplex.

### `factor/<lp>.jsonl` (HFactor)

All bases factor the **scaled** matrix HEkk factorises (`a` in `factor_lp`).
`HFactor::setupGeneral(num_col, num_row, num_row, …, pivot_threshold = 0.1,
pivot_tolerance = 1e-10, debug = 0, log = NULL, use_original_HFactor_logic =
true, update_method)` — exactly what `HSimplexNla::setup` does, except for the
explicit `update_method`.

| `t` | fields |
|---|---|
| `factor_lp` | `name`, `num_col`, `num_row`, `optimal_iteration_count`, `mid_iteration_count`, `pivot_threshold`, `pivot_tolerance`, `a` |
| `build` | `basis` (kind), `basic_index` (input), `result`, `solves` |
| `refactor` | `basis`, `basic_index` (input), `refactor_pivot_row/var/type`, `refactor_build_synthetic_tick` (the `RefactorInfo` given with `use = true`), `result`, `solves` |
| `update_sequence_start` | `lp`, `basis`, `method`, `update_method`, `basic_index`, `rank_deficiency`, `basic_index_after` |
| `update_step` | `k`, `variable_in`, `row_out`, `variable_out`, `aq_density`, `ep_density`, `aq` / `ep` (with pack, inputs to `update`), `hint`, `row_after`, `aq_after` / `ep_after` (with pack, after `update`), `solves` |
| `update_sequence_end` | `basic_index`, `invert` |

`result` = {`rank_deficiency`, `basic_index_after` (HFactor permutes
`basic_index` in place), `row_with_no_pivot`, `col_with_no_pivot`,
`var_with_no_pivot`, `build_synthetic_tick`, `basis_matrix_num_el`,
`invert_num_el`, `kernel_dim`, `kernel_num_el`, `refactor_pivot_row`,
`refactor_pivot_var`, `refactor_pivot_type` (the pivot sequence recorded in
`refactor_info_`), `refactor_build_synthetic_tick`, `invert`}.

`invert` = every array of `HFactor::getInvert()` (`InvertibleRepresentation`)
at full vector length, including any reserved tail: `l_pivot_index`,
`l_pivot_lookup`, `l_start`, `l_index`, `l_value`, `lr_start`, `lr_index`,
`lr_value`, `u_pivot_lookup`, `u_pivot_index`, `u_pivot_value`, `u_start`,
`u_last_p`, `u_index`, `u_value`, `ur_start`, `ur_lastp`, `ur_space`,
`ur_index`, `ur_value`, `pf_start`, `pf_index`, `pf_value`, `pf_pivot_index`,
`pf_pivot_value`.

`solves` = list of {`op` (`ftran`/`btran` = `ftranCall`/`btranCall(HVector&,
expected_density)`; `ftran_vec`/`btran_vec` = the `std::vector<double>&`
overloads), `rhs_kind`, `expected_density` (HVector ops only), `in`, `out`}.
After `build`/`refactor`: RHS {unit row 0, unit m/2, unit m-1, column a_{n/2},
sparse, dense} × densities {1.0, 0.1, 0.01} × {ftran, btran} plus the two
vector overloads on the dense RHS. After each `update_step`: unit m/2 at 0.01
and dense at 1.0, ftran and btran.

Update-step construction (deterministic, reproducible from the record alone):
`aq.clear(); aq.packFlag = true; a.collectAj(aq, variable_in, 1);
ftranCall(aq, aq_density)`; `row_out` = first index of max `|aq|`;
`ep.clear(); ep.count = 1; ep.index[0] = row_out; ep.array[row_out] = 1;
ep.packFlag = true; btranCall(ep, ep_density)`; `hint = 0;
update(&aq, &ep, &row, &hint)`; then the caller sets
`basic_index[row_out] = variable_in` (as `HEkk::updatePivots` does). A
sequence stops after a non-zero `hint`.

### `simplex/<lp>.<dual|primal>.jsonl` (HEkk end-to-end)

Records in emission order. `solve_call` counts `HEkk::solve` calls in the run.

| `t` | when | fields |
|---|---|---|
| `trace_header` | first line | `name`, `simplex_strategy`, `random_seed`, `highs_version`, `highs_githash` |
| `ekk_solve_enter` | `HEkk::solve`, before `initialiseForSolve` | `solve_call`, `iteration_count`, `status`, `options` (every HighsOptions record except file paths), `cost_scale`, `lp` (the scaled LP HEkk received), `basis` (only if `status.has_basis`) |
| `invert` | end of every `HEkk::computeFactor` that ran INVERT | `iteration_count`, `rank_deficiency`, `factor_build_synthetic_tick`, `build_synthetic_tick`, `total_synthetic_tick`, `basis_matrix_num_el`, `invert_num_el`, `kernel_dim`, `kernel_num_el`, `has_invert`, `h_basicIndex` (after INVERT's permutation) |
| `ekk_solve_initialised` | after `initialiseForSolve` | `status`, `basis`, `scalars`, `numTotRandomValue`, `numTotPermutation`, `numColPermutation`, strategy fields, perturbation multipliers, `factor_pivot_threshold`, `update_limit`, `arrays`, `ar_matrix` (partitioned row-wise matrix) |
| `dual_start` | `HEkkDual::solve`, before the major loop | `dual_call`, `solve_phase`, `edge_weight_mode` (0 Dantzig, 1 Devex, 2 steepest edge), `initial_basis_is_logical`, `allow_dual_steepest_edge_to_devex_switch` (only when `edge_weight_mode` = 2: otherwise HEkkDual leaves it uninitialised and never reads it), `Tp`, `Td`, `objective_bound`, `force_phase2`, `scalars`, `arrays` |
| `primal_start` | `HEkkPrimal::solve`, before the major loop | `primal_call`, `state`, `scalars`, `arrays` |
| `rebuild` | after every `rebuild()` | `algorithm`, `reason` (`rebuild_reason` on entry: 0 = first rebuild of a phase, 1 update limit, 2 synthetic clock, 3 possibly optimal, 4 possibly phase-1 feasible, 5 possibly primal unbounded, 6 possibly dual unbounded, 7 possibly singular, 8 primal infeasible in primal, 9 CHUZC fail, 10 force refactor, 11 excessive primal value, -1 cleanup), `solve_phase`, `rebuild_reason_after`, dual only: `dualInfeasCount`, `edge_weight_mode`; primal: `state`; `scalars`, `hash` |
| `iter` | after an `iterate()` that performed a basis change or bound flip (`iteration_count_` advanced) | see below |
| `iter_none` | after an `iterate()` that returned without iterating (CHUZR found nothing, bad basis change, …) | `algorithm`, `iteration`, `solve_phase`, `rebuild_reason`, `row_out` (dual) / `p` (primal; before the first iteration of an `HEkkPrimal::solve` call `p` omits the pivot fields `move_in`, `variable_out`, `move_out`, `theta_*`, `value_in`, `alpha_*`, `numericalTrouble`, which are uninitialised then) |
| `cleanup` | end of `cleanup()` (perturbation / shift removal) | `algorithm`, `solve_phase`, `scalars`, `hash` (+ `state` for primal) |
| `phase_end` | after `solvePhase1()` / `solvePhase2()` returns | `algorithm`, `solve_phase` (the phase decided next; 0 = exit), `rebuild_reason`, `model_status`, `iteration` |
| `ekk_solve_exit` | start of `HEkk::returnFromEkkSolve` | `model_status` (HighsModelStatus as int, 7 = optimal), `exit_algorithm`, `solve_bailout`, `scalars`, `status`, `basis`, `arrays` |
| `result` | last line, after `Highs::run` | `run_status`, `model_status`, `model_status_string`, `objective_function_value`, `simplex_iteration_count`, `primal_solution_status`, `dual_solution_status`, `num_primal_infeasibilities`, `num_dual_infeasibilities`, `basis_valid`, `col_status`, `row_status` (HighsBasisStatus as int), `col_value`, `col_dual`, `row_value`, `row_dual` (unscaled) |

**Dual `iter`**: `algorithm`, `iteration` (`iteration_count_` after the
iteration), `solve_phase`, `rebuild_reason`, `row_out`, `variable_out`,
`variable_in`, `move_out`, `delta_primal`, `theta_dual`, `theta_primal`,
`alpha_col`, `alpha_row`, `numericalTrouble`, `computed_edge_weight` (DSE
weight of `row_out` computed in CHUZR, the weight actually used),
`edge_weight_mode`, `num_flips` / `flips` (BFRT bound flips: `dualRow.workCount`
and the flipped variables `dualRow.workData[0..workCount).first`, in order),
`update_count`, `updated_dual_objective_value`, `total_synthetic_tick`,
`row_ep_density`, `row_ap_density`, `col_aq_density`, `row_DSE_density`,
`edge_weight_row_out_after` (`dual_edge_weight_[row_out]` after the update,
i.e. the weight of the entering variable), `h_basicIndex`, `h_workDual`,
`h_baseValue`, `h_dual_edge_weight`.

**Primal `iter`**: `algorithm`, `iteration`, `p` {`solve_phase`,
`rebuild_reason`, `variable_in`, `move_in`, `row_out` (-1 = bound flip),
`variable_out`, `move_out`, `theta_dual`, `theta_primal`, `value_in`,
`alpha_col`, `alpha_row`, `numericalTrouble`, `num_flip_since_rebuild`,
`edge_weight_mode`, `h_edge_weight`, `edge_weight_in`}, `update_count`,
`updated_primal_objective_value`, `total_synthetic_tick`, `h_basicIndex`,
`h_workDual`, `h_baseValue`, `h_workValue`.

`status` = the `HighsSimplexStatus` flags. `scalars` = {`iteration_count`,
`update_count`, `updated_dual_objective_value`,
`updated_primal_objective_value`, `dual_objective_value`,
`primal_objective_value`, `num_/sum_primal_infeasibilities`,
`num_/sum_dual_infeasibilities`, `costs_perturbed`, `costs_shifted`,
`bounds_perturbed`, `bounds_shifted`, `cost_perturbation_base`,
`cost_perturbation_max_abs_cost`, `build_synthetic_tick`,
`total_synthetic_tick`, the `*_density` running averages,
`dual_simplex_cleanup_level`, `simplex_in_scaled_space`}. `arrays` = full
`workCost`, `workDual`, `workShift`, `workLower`, `workUpper`, `workRange`,
`workValue`, `baseLower`, `baseUpper`, `baseValue`, `dual_edge_weight`.
`hash` = FNV-1a of `basicIndex`, `nonbasicMove`, and the same arrays (minus
`workRange`). Primal `state` = {`solve_phase`, `rebuild_reason`,
`edge_weight_mode`, `primal_feasibility_tolerance`,
`dual_feasibility_tolerance`, `num_free_col`, `num_flip_since_rebuild`,
`h_edge_weight`}.

Some scalars are only meaningful at certain points (e.g.
`num_dual_infeasibilities` is -1 until computed); they are dumped as HiGHS
holds them, uninterpreted.

## IPX records (`ipx/model`, `ipx/linalg`, `ipx/trace`, `ipx/run`)

The `ipx-model` and `ipx-linalg` commands (`src/ipx.cpp`) drive the IPX
classes of HiGHS 1.15.1 directly, below `LpSolver`. Private members are read
through `#define private public`, in that translation unit only. The
`ipx-trace` records are written by the second patch
(`highs-oracle-ipx-hooks.patch`) during a full `Highs::run` with
`solver=ipx`.

**Common objects:**
- **`Int`:** IPX's `ipxint` (= `HighsInt`, int32 in this build), written as a JSON number.
- **`ipx_mat`:** an `ipx::SparseMatrix` with `rows`, `cols`, `colptr` (`cols+1` entries), and `rowidx` and `values` (the first `colptr[cols]` entries; the capacity beyond that is not dumped).
- **`info`:** every `ipx::Info` field except the wall-clock `time_*` fields, in `ipx_info.h` order.
- **`iterate`:** an `ipx::Iterate`, with:
  - `tag`;
  - `x`, `xl`, `xu`, `y`, `zl`, `zu`;
  - `state` (`variable_state_` as int: 0 barrier lb, 1 barrier ub, 2 barrier box, 3 free, 4 fixed, 5 implied lb, 6 implied ub, 7 implied eq);
  - `feasibility_tol`, `optimality_tol`, `start_crossover_tol`;
  - `bounds_measure`, `costs_measure`;
  - `postprocessed`;
  - the evaluated quantities (the getters run `Evaluate()`): `pobjective`, `dobjective`, `pobjective_after_postproc`, `dobjective_after_postproc`, `presidual`, `dresidual`, `complementarity`, `mu`, `mu_min`, `mu_max`, `offset`, `feasible`, `optimal`, `term_crit_reached`;
  - the residual vectors `rb`, `rl`, `ru`, `rc`;
  - `scaling` (`ScalingFactor(j)` for every j).
- **`basis`:** an `ipx::Basis`, with:
  - `basis`, `map2basis`;
  - `fresh` (`factorization_is_fresh_`);
  - `factorizations`, `updates`, `ftran`, `btran`, `ftran_sparse`, `btran_sparse`;
  - `fill_factors`, `sum_ftran_density`, `sum_btran_density`;
  - `mean_fill`, `max_fill`;
  - `pivottol`;
  - `lu`, which is `"basiclu"` or `"forrest_tomlin"`:
    - BasicLu adds the BASICLU memory sizes after the `Reallocate` loop (`lmem`, `umem`, `wmem`), `istore_size`, `xstore_size`, `Li_size`, `Ui_size`, `Wi_size` and `fill_factor`;
    - ForrestTomlin adds `ft_updates` (the number of eta columns) and `fill_factor`;
  - `factors`, only when `fresh`: `GetLuFactors` → `L`, `U` (`ipx_mat`), `rowperm` and `colperm`.

**Pseudo-random inputs.** The vectors the oracle generates come from `HighsRandom`: seed 17 in `ipx-model`, seed 23 in `ipx-linalg`. Every record carries the vectors it used, so a replay never regenerates them.

### `ipx/model/<lp>.<cfg>.jsonl` (`ipx-model`)

`<cfg>` is `default` or `dual` (`dualize=1`).

| `t` | fields |
|---|---|
| `params` | `cfg` (the override string), then every `ipx::Parameters` field as the oracle set it. The base values are IpxWrapper's for the HiGHS defaults: feasibility tol 1e-7, optimality tol 1e-8, start_crossover 1e-8, crossover on, silent. |
| `input` | the `Model::Load` arguments, exactly as IpxWrapper's `fillInIpxData` builds them: `name`, `num_var`, `num_constr`, `offset`, `obj`, `lb`, `ub`, `Ap`, `Ai`, `Ax`, `rhs`, `constr_type` (chars as ints: `=` 61, `<` 60, `>` 62) |
| `model` | `errflag`, `dualized`, `num_rows`, `num_cols`, `num_dense_cols`, `nz_dense`, `AI`, `AIt`, `b`, `c`, `lb`, `ub`, `norm_bounds`, `norm_c`, `num_constr`, `num_eqconstr`, `num_var`, `num_free_var`, `num_entries`, `boxed_vars`, `constr_type`, `norm_obj`, `norm_rhs`, `offset`, `scaled_obj`, `scaled_rhs`, `scaled_lbuser`, `scaled_ubuser`, `A`, `flipped_vars`, `colscale`, `rowscale`, `filippo` (`filippoDualizationTest()`), `info` (`GetInfo`) |
| `presolve_starting_point` | inputs `x_user`, `slack_user`, `y_user`, `z_user`; outputs `x`, `y`, `z` (solver space) |
| `presolve_ipm_starting_point` | inputs `x_user`, `xl_user`, `xu_user`, `slack_user`, `y_user`, `zl_user`, `zu_user` (`xl`/`xu` are `inf` and `zl`/`zu` are 0 where the user bound is infinite, using the user bounds of flipped variables; slack and `y` signs follow the constraint type); `errflag`; outputs `x`, `xl`, `xu`, `y`, `zl`, `zu` |
| `postsolve_interior` | inputs `x_solver`, `xl_solver`, `xu_solver`, `y_solver`, `zl_solver`, `zu_solver` (infinite solver bounds have `xl`/`xu` = `inf` and `zl`/`zu` = 0); outputs `x`, `xl`, `xu`, `slack`, `y`, `zl`, `zu` (`PostsolveInteriorSolution`); `info` (`EvaluateInteriorSolution`) |

The three transform records are present only when `errflag` is 0.

### `ipx/linalg/<lp>.<cfg>.jsonl` (`ipx-linalg`)

`<cfg>` is one of:
- `default`;
- `lu1` (`lu_kernel=1,update_heuristic=0`: the ForrestTomlin kernel);
- `dual` (`dualize=1`);
- `slack` (`crash_basis=0`).

**Driver.** The real C++ IPM moves the iterate between stages, following `LpSolver::InteriorPointSolve`:
1. `StartingPoint`;
2. `RunInitialIPM` with KKTSolverDiag;
3. `BuildStartingBasis`;
4. at most 4 main-IPM iterations with KKTSolverBasis.

The Rust replay does not rerun the IPM. Instead, it reloads the dumped iterate at each stage and replays every linear-algebra component on it. For both KKTSolverDiag instances, `kkt_maxiter = min(500, 10 + m/20)`.

**Records, in order.** `params`, `input` and `model` are as in `ipx/model`. The records that follow are:

| `t` | fields |
|---|---|
| `iterate` | `iterate` object, top-level, with `tag` `constructed` (`Iterate(model)` plus the LpSolver tolerances), `starting_point`, `initial_ipm` or `main_ipm` |
| `ipm_starting_point`, `ipm_initial`, `ipm_main` | `info` after that IPM stage. `ipm_starting_point` is the last record when `status_ipm` is no longer `not_run`. |
| `kkt_diag` | `tag` (`starting_point_null`: `Factorize(nullptr)`; `starting_point`; `initial_ipm`), `with_iterate`, `info_in`, `info_factorize`, `basis_changes`, `iterate` (after `Factorize`, when `with_iterate`), then the solve: `a`, `b`, `tol`, `x`, `y`, `iter`, `info_solve`. The solve fields are absent when `info_factorize.errflag` is non-zero. |
| `starting_basis` | `info_in`, `info`, `basis`, `iterate` (after `StartingBasis`: crash basis and dependency handling) |
| `basis_ops` | `tag` (`starting_basis` / `basis_load`), `rhs`, `lhs_n` / `lhs_t` (`SolveDense` with `'N'` / `'T'`), `lhs_alias` (`SolveDense(rhs, rhs, 'N')`, aliased as in starting_basis.cc:76), `min_singular_value`, `rowcounts`, `colcounts` (`SymbolicInvert`), `density_inverse`, `random_permute` (symbolic_invert.cc `RandomPermute` of `basis`, the same code against the same libc++ `default_random_engine`) |
| `kkt_basis` | as `kkt_diag`, plus `basis` (after `Factorize`: `DropPrimal`/`DropDual`, maxvolume, basis updates). `tag` `starting_basis` uses the starting-basis iterate; `basis_load` uses the `main_ipm` iterate and the reloaded basis. |
| `basis_load` | `status` (`StatusOf(j)` of the post-main-IPM basis), `errflag` (`Basis::Load` into a fresh Basis), `basis` |
| `postsolve` | `info` and `iterate` (`postprocessed`), after `Postprocess` and `EvaluatePostsolved`, as `LpSolver::InteriorPointSolve` does; `x`, `xl`, `xu`, `slack`, `y`, `zl`, `zu` (`PostsolveInteriorSolution` of that iterate); `x_crossover`, `y_crossover`, `z_crossover` (`DropToComplementarity`); `x_basic`, `y_basic`, `z_basic` (`ComputeBasicSolution` on the loaded basis); `basic_status`; `xb_user`, `slackb_user`, `yb_user`, `zb_user` (`PostsolveBasicSolution`); `info_basic` (`EvaluateBasicSolution`); `cbasis`, `vbasis` (`PostsolveBasis`) |

**`basic_status`** is synthesised, in this order:
1. `IPX_basic` if `IsBasic(j)`.
2. Otherwise, `nonbasic_lb` / `nonbasic_ub` when `x_crossover[j]` equals that bound.
3. Otherwise, `superbasic` for a structural (j < n).
4. Otherwise, the slack's finite bound (`nonbasic_lb` before `nonbasic_ub`), or `IPX_nonbasic` when both are infinite.

Slacks are never superbasic, because `Model::DualizeBackBasis` asserts that.

**Dualized models.** When the model is dualized, `info.dependent_rows`/`dependent_cols` and `rows_inconsistent`/`cols_inconsistent` are swapped after `starting_basis`, as `LpSolver::BuildStartingBasis` does.

A run stops early, after the record that failed, on:
- a non-zero `errflag` from `Load` or `StartingBasis`;
- a non-zero `errflag` from `kkt_basis` on the starting basis;
- a non-zero `errflag` from `basis_load`.

### `ipx/trace/<lp>.<cfg>.jsonl` (`ipx-trace`, oracle (b))

The run is a full `Highs::run` with:
- `solver=ipx`;
- presolve off, serial, and the base oracle options (`random_seed=0`);
- integrality cleared.

`<cfg>` is `default` or `dual` (`ipx_dualize_strategy=1`).

The patch-2 hooks run only when `HIGHS_ORACLE_IPX_DUMP` names a file. They write the `ipm` and `info` lines, and the oracle frames them with a header and a result:

| `t` | when | fields |
|---|---|---|
| `trace_header` | first line | `name`, `cfg`, `solver`, `highs_version`, `highs_githash` |
| `ipm` | after every `IPM::PrintOutput` | `call` (0-based over all hook records), `kind` (`start`: after `ComputeStartingPoint`; `iter`: each `Driver` iteration; `centring`: each centring step), `iter`, `mu`, `step_primal`, `step_dual`, `presidual`, `dresidual`, `pobjective`, `dobjective`, `pobjective_after_postproc`, `dobjective_after_postproc`, `complementarity`, `mu_min`, `mu_max`, `kkt_iter`, `basis_changes` (of the active KKT solver), `kktiter1`, `kktiter2`, `updates_ipm`, `primal_dropped`, `dual_dropped`, `errflag`, `h_x`, `h_xl`, `h_xu`, `h_y`, `h_zl`, `h_zu` (FNV-1a of the f64 bits) |
| `info` | `kind` `driver_end` at the end of every `IPM::Driver`; `solve_end` at the end of `LpSolver::Solve` (after crossover) | `call`, `kind`, `info` fields inline (not nested) |
| `result` | last line | `run_status`, `model_status`, `model_status_string`, `objective_function_value`, `ipm_iteration_count`, `crossover_iteration_count`, `simplex_iteration_count`, `primal_solution_status`, `dual_solution_status`, `basis_valid`, `col_status`, `row_status`, `col_value`, `col_dual`, `row_value`, `row_dual` |

**Calls.** `call` orders all hook records of the run. Every committed and target trace holds exactly one `LpSolver::Solve`, so it has one `start` … `solve_end` group, and the test checks that. A run that re-solved would show several groups.

### `ipx/run/<lp>.jsonl` (`ipx`, oracle (b) + (c))

End-to-end IPX through `Highs::passModel` + `Highs::run` (integrality cleared), with the
base oracle options (presolve off, serial, `random_seed=0`) and `solver=ipx`, in four
configurations, each on a fresh `Highs`:

| `cfg` | overrides |
|---|---|
| `xon` | `run_crossover=on` |
| `xoff` | `run_crossover=off` |
| `ipm` | `solver=ipm,run_crossover=on` (routed by `useHipo`, see `hipo_available`) |
| `dual` | `ipx_dualize_strategy=1,run_crossover=on` |

| `t` | when | fields |
|---|---|---|
| `ipx_header` | first line | `name`, `highs_version`, `highs_githash`, `hipo_available` (`HighsExternalApi::isAvailable<HighsExtras::hipo>()`, 0 in the oracle build), `lp` (as read) |
| `run_header` | before each configuration | `cfg`, `overrides` |
| `ipm`, `info` | the patch-2 hook lines of that run (counter reset per run), exactly as in `ipx/trace` | see `ipx/trace` |
| `run_result` | after each configuration | `cfg`, `run_status`, `model_status`, `info` (every LP `HighsInfo` field), `value_valid`, `dual_valid`, `col_value`, `col_dual`, `row_value`, `row_dual`, `basis_valid`, `basis_alien`, `col_status`, `row_status` |

The Rust test (`tests/ipx_run.rs`) captures the port's own hook lines
(`ipm/ipx/oracle_trace.rs`) and compares them **byte for byte** with the `ipm`/`info`
lines, then compares the `run_result` bit for bit.

## MIP IPM records (`mip/<instance>.ac`)

Written by `mip_trace` when `HIGHS_ORACLE_AC_DUMP=<path>` is set
(`gen_fixtures.sh --mip` runs it into `<target fixture dir>/mip/`; the committed
copies are not overwritten). Each line is one IPM (IPX) result that the MIP solver uses and that the
port cannot reproduce itself. The file is truncated by the first record. Nothing
is written to stdout, so the `.trace` is unchanged. Command:
`HIGHS_ORACLE_AC_DUMP=<out>.ac mip_trace <X>.mps [presolve=off ...]`.

Every record has `t`, `call` (0-based index over all records of the process, in
emission order; the key is `(submip_level, call)`), `submip` (0/1),
`submip_level` (`HighsMipSolver::submip_level`, 0 at the root), and
`model_status` (`HighsModelStatus` as an integer, 7 = optimal, 15 = unknown).
Vectors are JSON arrays of hexfloat strings (`%a`, exact; `-0.0` keeps its
sign) unless noted.

| `t` | when | fields |
|---|---|---|
| `ac` | `HighsMipSolverData::startAnalyticCenterComputation`, after the IPM solve and before the size check (so both success and failure are recorded) | `stored` (1 when the size check passed, so `analyticCenterStatus` and `analyticCenter` were assigned; 0 means the IPM failed and the centre stays empty), `centring_status` (1.16.0 oracle only, see below), `centre` (every entry of the IPM's `col_value`, length = number of columns when `stored` = 1) |
| `lp_ipm` | `HighsLpRelaxation::run`, IPM fallback on `kIterationLimit` (HighsLpRelaxation.cpp ~1358), after the IPM solve and any HiPO-to-IPX retry, before `setBasis` | `ipm_solver` (`"ipx"` or `"hipo"`, read back from the `solver` option), `basis_valid` (0/1), `col_status` and `row_status` (`HighsBasisStatus` as integers, the basis after crossover), `col_value` and `row_value` (the IPM solution, unscaled) |

`centring_status` (written between `stored` and `centre`, by the 1.16.0
patch set `patches-1.16.0/` only) is `ipm.getInfo().centring_status` as an
integer: 0 = `kCentringStatusNone`, 1 = `kCentringStatusFailed`, 2 =
`kCentringStatusSuccess`. `model_status` stays the raw IPM status: HiGHS 1.16
then stores `analyticCenterStatus = kUnknown` whenever `centring_status` is
not 2, so a replay derives the stored status from both fields. Readers ignore
unknown keys (`ipm_replay.rs`), so 1.15.1 and 1.16.0 records share one parser.

The IPM runs in the analytic-centre and LP-fallback code, so `submip_level` is
the depth of the submip that made the call. Submip records appear when the
heuristics run a sub-MIP (for example `air05`, `call` 1-6).

Committed fixtures: `crates/rusolve-ekk/tests/fixtures/mip/{sp150x300d,air05,pg,enlight_hard}.ac`
(the traces at default options), and `sp150x300d.presolve_off.{trace,ac}`. The
`enlight_hard.ac` file contains an `lp_ipm` record (its root LP hits the iteration
limit and the IPM fallback runs). The other three contain only `ac` records.

## W1 records (`w1/qp`, `w1/ranging`, `w1/iis`, `w1/icrash`)

Written by the W1 oracle commands (`src/w1.cpp`: `w1-qp`, `w1-ranging`,
`w1-iis`, `w1-icrash`). Each
drives `Highs::passModel` (+ `passHessian`) and `Highs::run` on a fresh
`Highs` with the base oracle options (output off, `threads=1`,
`parallel=off`, `presolve=off`, `random_seed=0`) plus the record's
`overrides` (`k=v,k=v`, or `default`). The QUASS hook lines are written by the
third patch (`highs-oracle-w1-hooks.patch`, `highs/qpsolver/HighsOracleQpDump.h`)
when `HIGHS_ORACLE_QP_DUMP=<path>` is set; `w1-qp` writes the `w1_header`
line, closes the file, sets the variable to the fixture itself and runs; the
hook opens the file in append mode on its first record, so the hook lines
appear in place. All W1 fixtures are `target`-only. Instances named after
`reference/HiGHS/check/instances` are read from `$HIGHS_REF/check/instances`
(not copied into `lps/`).

Synthetic LPs live in `lps/w1/` (not `lps/`, so the IPX loop over
`lps/*.mps` does not pick them up).

**Common records:**

| `t` | when | fields |
|---|---|---|
| `w1_header` | first line | `name`, `highs_version`, `highs_githash`, `overrides`, `lp` (the `HighsLp` passed to `passModel`, encoding as `lp` above, unscaled), `hessian` (w1/qp only: `dim`, `format` (1 triangular, 2 square), `start`, `index`, `value`, as passed to `passHessian`) |
| `w1_result` | after `Highs::run` | `run_status`, `model_status`, `objective_function_value`, `simplex_iteration_count`, `qp_iteration_count`, `primal_solution_status`, `dual_solution_status`, `num_primal_infeasibilities`, `max_primal_infeasibility`, `sum_primal_infeasibilities`, `num_dual_infeasibilities`, `max_dual_infeasibility`, `sum_dual_infeasibilities`, `basis_valid`, `col_status`, `row_status` (HighsBasisStatus as int), `col_value`, `col_dual`, `row_value`, `row_dual` |

### `w1/qp/<inst>.jsonl` (`w1-qp`, QUASS)

Instances: `qjh`, `qjh_uncon`, `qjh_quadobj`, `qjh_qmatrix`, `2821`, `2894`,
`atwood0`, `primal1` (`reference/HiGHS/check/instances`), `overrides=default`
(the QP path is selected by the non-empty Hessian). Records: `w1_header`,
then the hook lines, then `w1_result`.

Hook lines (`call` = 0-based counter over the hook lines of the run):

| `t` | when | fields |
|---|---|---|
| `qp_iter` | end of every pass of the `while (true)` loop of `Quass::solve` (quass.cpp 351-524): after `redcosts.update()` on a step, after `atfsep = true` on a null step | `call`, `iter` (`runtime.statistics.num_iterations`), `atfsep` (0/1, after the pass), `minidx` (the pricing result of a major pass, -1 on a minor pass), `limiting` (`stepres.limitingconstraint`, -2 on a null step), `alpha` (`stepres.alpha`, `0` on a null step), `num_active`, `num_inactive` (`basis.getnumactive/getnuminactive()`), `h_primal`, `h_rowactivity`, `h_gradient` (the stored `Gradient::gradient` vector, read through the patch's `Gradient::oracleGradient()` accessor: `getGradient()` may recompute it, which would perturb the run), `h_p` (the search direction) |
| `qp_end` | after the post-loop `loginformation` of `Quass::solve` (quass.cpp 526) | `call`, `iter`, `status` (`QpModelStatus` as int, 2 = optimal), `h_primal`, `h_rowactivity` |

The `h_*` fields hash `QpVector::value[0..dim)` (the dense array; `num_nz`
and `index` are not hashed). The Rust test (`tests/t_w1_qp.rs`) captures the
port's lines (`qpsolver/oracle_trace.rs`, env `RUSOLVE_QP_TRACE`) and diffs
them **byte for byte**, then compares `w1_result` bit for bit.

Notes on the trace shape (all faithful HiGHS 1.15.1 behaviour):

- A major pass increments `num_iterations` twice, so `iter` advances by 2.
- A QP without constraints (`qjh_uncon`) is solved by
  `computeStartingPointBounded` (a_quass.cpp) and never enters
  `Quass::solve`: its fixture has no hook lines at all.
- The early `return`s of `Quass::solve` (unbounded, error, nonconvex, large
  nullspace) skip the post-loop `loginformation`, so they write no `qp_end`.
  `2894` ends that way (10 `qp_iter` lines, no `qp_end`; `run_status=-1`,
  `model_status=4` kSolveError), as `check/TestQpSolver.cpp` expects.

### `w1/ranging/<lp>.jsonl` (`w1-ranging`, getRanging)

LPs (netlib): `afiro`, `sc50a`, `kb2`, `adlittle`, `blend`, `share2b`,
`sc105`, `scagr7`, `overrides=default` (simplex). Records:

| `t` | when | fields |
|---|---|---|
| `w1_header` | first line | as above |
| `w1_result` | after `run` | as above |
| `ranging_status` | after `Highs::getRanging` | `status` (HighsStatus as int), `valid` |
| `ranging` | one per array, in the order `col_cost_up`, `col_cost_dn`, `col_bound_up`, `col_bound_dn`, `row_bound_up`, `row_bound_dn` | `which`, `value`, `objective`, `in_var`, `ou_var` (`HighsRangingRecord::{value_,objective_,in_var_,ou_var_}`) |

### `w1/iis/<lp>.jsonl` (`w1-iis`, getIis)

LPs: `galenet`, `forest6`, `woodinfe`, `klein1`
(`reference/HiGHS/check/instances`), `iis_rows`, `iis_bnds` (`lps/w1/`).
Records: `w1_header` (`overrides=default`), then one `iis_result` per
configuration, each on a fresh `Highs`, `run` then `getIis`:

| `cfg` | overrides |
|---|---|
| `light` | `iis_strategy=0` |
| `lp` | `iis_strategy=2` |
| `ray_irr` | `iis_strategy=5` |
| `lp_irr` | `iis_strategy=6` |

| `t` | fields |
|---|---|
| `iis_result` | `cfg`, `overrides`, `run_status`, `model_status`, `get_iis_status`, `valid`, `status` (`kIisModelStatus*`), `strategy`, `col_index`, `row_index`, `col_bound`, `row_bound` (`kIisBoundStatus*`), `col_status`, `row_status` (`kIisStatus*`), `info` {`num_lp_solved`, `sum_simplex_iteration_counts`, `min_simplex_iteration_count`, `max_simplex_iteration_count`} (the times are not dumped) |

### `w1/icrash/<lp>.jsonl` (`w1-icrash`, iCrash + crossover)

LPs (netlib): `afiro`, `sc50a`, `adlittle`; `overrides=icrash=true,simplex_strategy=4`
(the iCrash block of `Highs::run` asserts `kSimplexStrategyPrimal`). Records:

| `t` | when | fields |
|---|---|---|
| `w1_header` | first line | as above |
| `icrash_iter` | one per `Highs::icrash_info_.details` entry, after `run` | `num`, `weight`, `lambda_norm_2`, `lp_objective`, `quadratic_objective`, `residual_norm_2` (no `time`) |
| `icrash_info` | after the `icrash_iter` lines | `num_iterations`, `final_lp_objective`, `final_quadratic_objective`, `final_residual_norm_2`, `starting_weight`, `final_weight`, `x_values` (no `total_time`) |
| `w1_result` | last line | as above |

`adlittle` diverges inside iCrash: from `icrash_iter` 1 on, `lp_objective`,
`quadratic_objective`, `residual_norm_2` and (from 2) `lambda_norm_2` are NaN,
as are `final_lp_objective`, `final_quadratic_objective`,
`final_residual_norm_2` and `x_values`; crossover still reaches the optimum
(`w1_result` objective 225494.963...). The NaNs are the reference values (bit
pattern `7ff8000000000000`), not a writer fault; the port must reproduce them.
For `afiro` and `sc50a` `quadratic_objective` exceeds 1e12 once the weight
falls to 1e-13.

## W3 records (`w3/cupdlp`, `w3/hipdlp`)

The
Rust writer (`crates/rusolve-ekk/src/pdlp/oracle_trace.rs`, env
`RUSOLVE_PDLP_TRACE`) emits exactly these lines; change this section, the
patch and `oracle_trace.rs` together.

Written by the W3 oracle commands (`src/w3.cpp`: `w3-cupdlp`, `w3-hipdlp`).
Each drives `Highs::passModel` and `Highs::run` on a fresh `Highs` with the
base oracle options (output off, `threads=1`, `parallel=off`, `presolve=off`,
`random_seed=0`) plus the record's `overrides` (`k=v,k=v`). The hook lines are
written by the fourth patch (`highs-oracle-w3-hooks.patch`,
`highs/pdlp/HighsOracleW3Dump.h`, C-compatible because the `cupdlp_*.c` files
compile as C) when `HIGHS_ORACLE_PDLP_DUMP=<path>` is set; `w3-*` writes the
`w3_header` line, closes the file, sets the variable to the fixture itself and
runs; the hook appends, so the hook lines appear in place. All W3 fixtures are
`target`-only.

Integers are JSON numbers; doubles are 16 hex digits of the bit pattern; an
`h_*` field is the FNV-1a 64 hash (as in "Hashes") of the named vector.

LPs (netlib, `benches/data/netlib`): `afiro`, `sc50a`, `kb2`, `adlittle`,
`blend`, `share2b`, `sc105`, `scagr7`.

| fixture | command | overrides | hook records |
|---|---|---|---|
| `w3/cupdlp/<lp>.jsonl` | `w3-cupdlp` | `solver=pdlp,pdlp_iteration_limit=2000` | `pdlp_check`, `pdlp_restart`, `pdlp_end` |
| `w3/hipdlp/<lp>.jsonl` | `w3-hipdlp` | `solver=hipdlp,pdlp_iteration_limit=2000` | `hipdlp_check`, `hipdlp_restart`, `hipdlp_end` |

**Common records:**

| `t` | when | fields |
|---|---|---|
| `w3_header` | first line | `name`, `highs_version`, `highs_githash`, `overrides`, `lp` (the `HighsLp` passed to `passModel`, encoding as `lp` above, unscaled) |
| `w3_result` | last line, after `Highs::run` | `run_status`, `model_status`, `objective_function_value`, `pdlp_iteration_count`, `primal_solution_status`, `dual_solution_status`, `num_primal_infeasibilities`, `max_primal_infeasibility`, `sum_primal_infeasibilities`, `num_dual_infeasibilities`, `max_dual_infeasibility`, `sum_dual_infeasibilities`, `basis_valid`, `col_value`, `col_dual`, `row_value`, `row_dual` |

**Hook lines** (each starts `{"t":"<name>","call":N`, `call` = 0-based counter
over the hook lines of the run; the fields follow in this order):

| `t` | when | fields |
|---|---|---|
| `pdlp_check` | `cupdlp_solver.c::PDHG_Solve`, right after `PDHG_Compute_Infeas_Residuals(pdhg)` inside `if (bool_checking)` | `iter` (`timers->nIter`), `p_step`, `d_step`, `beta`, `ss_iter` (`stepsize->{dPrimalStep, dDualStep, dBeta, nStepSizeIter}`), `last_restart` (`iterates->iLastRestartIter`), `pobj`, `dobj`, `gap`, `compl`, `pfeas`, `dfeas`, `relgap` (`resobj->{dPrimalObj, dDualObj, dDualityGap, dComplementarity, dPrimalFeas, dDualFeas, dRelObjGap}`), `pobj_a`, `dobj_a`, `gap_a`, `compl_a`, `pfeas_a`, `dfeas_a`, `relgap_a` (the `*Average` fields), `pinf_obj`, `dinf_obj`, `pinf_res`, `dinf_res` (`resobj->{dPrimalInfeasObj, dDualInfeasObj, dPrimalInfeasRes, dDualInfeasRes}`), `h_x`, `h_y` (`iterates->x[nIter % 2]` over `nCols`, `y[nIter % 2]` over `nRows`), `h_xa`, `h_ya` (`xAverage`, `yAverage`) |
| `pdlp_restart` | `cupdlp_proj.c::PDHG_Restart_Iterate_GPU`, last statement (only when a restart happened) | `iter`, `choice` (1 current, 2 average), `p_step`, `d_step`, `beta`, `pfeas_lr`, `dfeas_lr`, `gap_lr` (`resobj->{dPrimalFeasLastRestart, dDualFeasLastRestart, dDualityGapLastRestart}`), `h_x`, `h_y` |
| `pdlp_end` | `cupdlp_solver.c::PDHG_Solve`, right after the `for` loop | `iter`, `term_code` (`resobj->termCode`), `term_iterate` (`resobj->termIterate`), `h_x`, `h_y`, `h_xa`, `h_ya` |
| `hipdlp_check` | `pdhg.cc::PDLPSolver::solve`, right after each `runConvergenceCheck` (initial and in-loop) | `iter` (`final_iter_count_`), `halpern` (`halpern_iteration_`), `stop` (0/1, the call's result), `p_step`, `d_step`, `beta` (`stepsize_.{primal_step, dual_step, beta}`), `omega`, `eta` (`params_.{omega, eta}`), `pw` (`primal_weight_`), `fpe` (`fpe_`), `ifpe` (`initial_fpe_`), `pobj`, `dobj`, `gap`, `pfeas`, `dfeas`, `relgap` (`results_.{primal_obj, dual_obj, duality_gap, primal_feasibility, dual_feasibility, relative_obj_gap}`), `h_x`, `h_y`, `h_xn`, `h_yn` (`x_current_`, `y_current_`, `x_next_`, `y_next_`, full length) |
| `hipdlp_restart` | `pdhg.cc::PDLPSolver::solve`, end of the `if (do_restart)` block | `iter`, `p_step`, `d_step`, `beta`, `omega`, `pw`, `best_pw` (`best_primal_weight_`), `pw_err_sum` (`primal_weight_error_sum_`), `pw_last_err` (`primal_weight_last_error_`), `h_x`, `h_y`, `h_xn`, `h_yn` |
| `hipdlp_end` | `pdhg.cc::PDLPSolver::solveReturn`, after `results_.term_code = term_code` | `iter`, `term_code` (`TerminationStatus`: NOTSET 0, OPTIMAL 1, INFEASIBLE 2, UNBOUNDED 3, MAXITER 4, TIMEOUT 5, ERROR 6), `h_x`, `h_y`, `h_xn`, `h_yn` |

The HiPDLP hooks sit inside `PDLPSolver` member functions and pass the member
values explicitly. The Rust tests (`tests/t_w3_{cupdlp,hipdlp}.rs`) capture
the port's lines and diff them **byte for byte**, then compare `w3_result`
bit for bit.

## W4 records (`w4/orderings`, `w4/hipo`)

The
Rust writer (`crates/rusolve-ekk/src/ipm/hipo/oracle_trace.rs`, env
`RUSOLVE_HIPO_TRACE`) emits exactly these lines; change this section, the
patch and `oracle_trace.rs` together.

Written by the W4 oracle commands (`src/w4.cpp`: `w4-orderings`, `w4-hipo`).
Each drives `Highs::passModel` and `Highs::run` on a fresh `Highs` with the
base oracle options (output off, `threads=1`, `parallel=off`, `presolve=off`,
`random_seed=0`) plus the record's `overrides` (`k=v,k=v`). They come from a
second oracle configuration, `target/highs-oracle/build-hipo`
(`-DORACLE_HIPO=ON`; the default `build/` stays `HIPO=OFF`), in which
`CMakeLists.txt` sets `HIPO=ON` and presets HiGHS' `HIGHS_BLAS_*` variables
(`HIGHS_BLAS_CONFIGURED`, `HIGHS_BLAS_TARGET=refblas`, vendor
`rusolve refblas`), so `FindHipoDeps` runs no BLAS search and HiPO links the
reference BLAS `tools/highs-oracle/refblas/refblas.c` (netlib loop order,
`-O2 -ffp-contract=off`, MIT, ours). (`-DBLAS_LIBRARIES` does not pin it:
`FindBLAS` clears that variable.) `gen_fixtures.sh` checks the configure log
for the pin and for the absence of a BLAS search, and runs the refblas unit
checks (`test_refblas`); `w4-*` exits 3 unless
`HighsExternalApi::isAvailable<HighsExtras::hipo>()` holds and the BLAS
provider is `rusolve refblas`. So the BLAS the oracle's HiPO calls is the one
`crates/rusolve-ekk/src/ipm/hipo/blas.rs` ports. The hook lines are written
by the fifth patch (`highs-oracle-w4-hooks.patch`,
`highs/ipm/hipo/HighsOracleW4Dump.h`) when `HIGHS_ORACLE_HIPO_DUMP=<path>` is
set; `w4-*` writes the `w4_header` line, closes the file, sets the variable to
the fixture itself and runs; the hook appends, so the hook lines appear in
place. All W4 fixtures are `target`-only.

Integers are JSON numbers; doubles are 16 hex digits of the bit pattern; an
`h_*` field is the FNV-1a 64 hash (as in "Hashes") of the named vector;
`ptr`, `rows` and `perm` are integer arrays written in full.

LPs (netlib, `benches/data/netlib`): `afiro`, `sc50a`, `kb2`, `adlittle`,
`blend`, `share2b`, `sc105`, `scagr7`.

| fixture | command | overrides | hook records |
|---|---|---|---|
| `w4/orderings/<lp>.<ordering>.jsonl` (`<ordering>` = `metis`, `amd`, `rcm`) | `w4-orderings` | `solver=hipo,hipo_ordering=<ordering>` | `hipo_ord`, `hipo_fact`, `hipo_iter`, `hipo_end` |
| `w4/hipo/<lp>.jsonl` | `w4-hipo` | `solver=hipo` (`hipo_system=choose`, `hipo_ordering=choose`) | `hipo_ord`, `hipo_fact`, `hipo_iter`, `hipo_end` |

Every fixture holds one run (one `w4_header`, its hook lines, one
`w4_result`); the Rust discovery (`tests/port/w4.rs::w4_fixtures`) takes every
`*.jsonl` of the directory.

**Common records:**

| `t` | when | fields |
|---|---|---|
| `w4_header` | first line | `name`, `highs_version`, `highs_githash`, `overrides`, `lp` (the `HighsLp` passed to `passModel`, encoding as `lp` above, unscaled) |
| `w4_result` | last line, after `Highs::run` | `run_status`, `model_status`, `objective_function_value`, `ipm_iteration_count`, `crossover_iteration_count`, `primal_solution_status`, `dual_solution_status`, `num_primal_infeasibilities`, `max_primal_infeasibility`, `sum_primal_infeasibilities`, `num_dual_infeasibilities`, `max_dual_infeasibility`, `sum_dual_infeasibilities`, `basis_valid`, `col_value`, `col_dual`, `row_value`, `row_dual` |

**Hook lines** (each starts `{"t":"<name>","call":N`, `call` = 0-based counter
over the hook lines of the run; the fields follow in this order):

| `t` | when | fields |
|---|---|---|
| `hipo_ord` | `FactorHiGHSSolver.cpp::chooseOrdering`: the lambda `run_ordering_and_analyse` keeps `permutations[i]` and `failure[i]` (right after `logger_.printInfo("Finished %s for %s\n", ...)`, before the `if (failure[i])` test); the W5 patch writes one line per ordering after the `for_each`, from the calling thread, for `i = 0, 1, ...` (`choose` tries `amd`, then `metis`). With `parallel=off` this is the order and content the W4 patch wrote from inside the lambda | `nla` (the `nla` argument, `"AS"` / `"NE"`), `ordering` (`orderings_to_try[i]`), `n`, `failure` (0/1, `failure[i]` after the ordering call), `ptr`, `rows` (`full_ptr`, `full_rows`: the full symmetric pattern without diagonal from `fullFromLower`, 0-based CSC), `perm` (`permutations[i]`) |
| `hipo_fact` | `FactorHiGHSSolver.cpp::setup`, right after `setParallel()` | `nla` (`"AS"` for `kOptionNlaAugmented`, `"NE"` for `kOptionNlaNormEq`, `options_.nla` after `setNla`), `ordering` (`ordering_AS_` / `ordering_NE_` of that system), `flops`, `spops`, `nz` (`S_.flops()`, `S_.spops()`, `(double)S_.nz()`) |
| `hipo_iter` | `Solver.cpp::Solver::runIpm`, right after each `makeStep()` | `iter` (`iter_`), `pobj`, `dobj`, `pinf`, `dinf`, `pdgap`, `mu` (`it_->{pobj, dobj, pinf, dinf, pdgap, mu}`), `alpha_p`, `alpha_d`, `sigma` (`alpha_primal_`, `alpha_dual_`, `sigma_`), `correctors` (`info_.correctors`), `h_x`, `h_xl`, `h_xu`, `h_y`, `h_zl`, `h_zu` (`it_->{x, xl, xu, y, zl, zu}`, full length) |
| `hipo_end` | `IpxWrapper.cpp::solveHipo`, right after `const hipo::Info hipo_info = hipo.getInfo()` | `status`, `ipm_iter`, `ipx_used` (0/1), `ipx_iter` (`ipx_info.iter`), `status_crossover`, `updates_crossover` (`ipx_info.*`), `m_solver`, `n_solver`, then the doubles `p_obj`, `d_obj`, `pd_gap`, `p_res_abs`, `p_res_rel`, `d_res_abs`, `d_res_rel` |

The `hipo_ord`, `hipo_fact` and `hipo_iter` hooks sit inside member functions
(`FactorHiGHSSolver`, `Solver`) and pass the member values explicitly. The
Rust tests capture the port's lines and diff them **byte for byte**
(`tests/t_w4_hipo.rs`), then compare `w4_result` bit for bit;
`tests/t_w4_orderings.rs` also replays each `hipo_ord` record standalone
(the recorded `n` / `ptr` / `rows` through the ported METIS / AMD / RCM with
the arguments of `chooseOrdering`) and checks `perm` and `failure`.
## W5 records (`w5/simplex`, `w5/hipo`, `w5/mip`)

These are the existing simplex, W4 HiPO and MIP traces, run with an
explicit `threads` = N and `parallel=on`. The record formats are unchanged:
the W1 simplex records ("Record types": `simplex/<lp>.<dual|primal>.jsonl`), the
W4 records ("W4 records") and the `mip_trace` text format
(`tests/fixtures/mip/<instance>.trace`). Only the header records gain keys.
There is no new oracle configuration: `w5-simplex` runs on the default
`target/highs-oracle/build`, `w5-hipo` on `build-hipo`, and the MIP traces use
`build/mip_trace`. Both builds have `HIGHS_NO_DEFAULT_THREADS`, so the oracle
always sets `threads`. Each fixture is a separate oracle process, so every run
starts the HiGHS scheduler fresh with N workers. All W5 fixtures are
`target`-only.

LPs (netlib, `benches/data/netlib`): `afiro`, `sc50a`, `kb2`, `adlittle`,
`blend`, `share2b`, `sc105`, `scagr7` (the W1 ranging / W3 / W4 set). MIP
instances (MIPLIB, `benches/data/miplib`): `air05`, `enlight_hard`, `pg`,
`sp150x300d` (the `tests/fixtures/mip` set).

| fixture | command | options on top of the base | records |
|---|---|---|---|
| `w5/simplex/<lp>_s<k>_t<N>.jsonl`, `k` in {2, 3}, N in {1, 2, 4, 8} | `w5-simplex <mps> <name> <k> <N> <out>` | `simplex_strategy=k` (2 = SIP, 3 = PAMI), `threads=N`, `parallel=on`; `simplex_max_concurrency` at its default (8) | the W1 simplex records |
| `w5/hipo/<lp>_<type>_<sys>_t<N>.jsonl`, `type` in {`tree`, `node`, `both`}, `sys` in {`choose` (N = 1 only), `as`, `ne`}, N in {1, 2, 4, 8} | `w5-hipo <mps> <name> <type> <choose\|augmented\|normaleq> <N> <out>` | `solver=hipo`, `threads=N`, `parallel=on`, `hipo_parallel_type=type`, `hipo_system` = `choose` / `augmented` / `normaleq` | the W4 records |
| `w5/mip/<instance>.t<N>.trace`, N in {2, 4, 8} | `mip_trace <mps> threads=N parallel=on` (no `HIGHS_ORACLE_AC_DUMP`) | `threads=N`, `parallel=on` | the `mip_trace` lines |

`w5-hipo` rejects `hipo_system=choose` at N > 1 (exit 64): `chooseNla` runs
the AS analyse and the NE structure as two tasks, which race at N > 1. The `as` / `ne` files exist at every N so
that each `(lp, type, sys)` can be compared across thread counts. The MIP
traces at threads = 1 are the committed `tests/fixtures/mip/<instance>.trace`.
The `.ac` dump (`HighsOracleAcDump.h`, a static `FILE*` and counter) is not
thread-safe and is off at N > 1.

**Header keys.** `trace_header` (simplex) is `name`, `simplex_strategy`,
`threads`, `parallel` (`"on"`), `random_seed`, `highs_version`,
`highs_githash`. `w4_header` (HiPO) is `name`, `highs_version`,
`highs_githash`, `overrides` (the full `k=v,...` list above, including
`threads=N` and `parallel=on`), `threads`, `parallel` (`"on"`), `lp`. W1 and
W4 fixtures have neither key (their writers pass threads 0 and the old header
is written unchanged). The `options` object in `ekk_solve_enter` already
records `threads` and `parallel`.

**`hipo_ord` order.** `highs-oracle-w5-hooks.patch` (sixth patch, applied
after the W4 patch with its own stamp `highs-src.w5.stamp`) moves the
`hipo_ord` write out of the `run_ordering_and_analyse` lambda, which
`highs::parallel::for_each` may run on several workers at once: the lambda
saves `permutations[i]` and `failure[i]`, and the records are written after
the `for_each` in index order from the calling thread. At `parallel=off` (all
W4 fixtures) the bytes are those of the W4 patch.

**Determinism.** `gen_fixtures.sh --determinism` writes `w5/` three times
(`target/highs-fixtures`, `det-a`, `det-b`) and fails on any difference.
Identity of a `(lp, strategy)` or `(lp, type, sys)` trace across thread
counts is a measured property, not a format rule: compare the files without
the header line (whose `threads` differs) and without the `threads` value in
`ekk_solve_enter.options` / the `w4_header` overrides.

## G1 records (`g1/<dir>/<name>.<scenario>.jsonl`)

The
consumer is `crates/rusolve-ekk/tests/port/g1.rs` (`g1_options`,
`apply_g1_call`, `check_g1_result`, `replay_g1`); change this section, `src/g1.cpp`
and `g1.rs` together.

Written by the G1 oracle commands (`src/g1.cpp`), all through the public
`Highs` API (no hook patch, so no new patch hash):

    highs-oracle g1-<x> <model> <name> <scenario> <out.jsonl>

Each command reads `<model>` (`.mps` or `.lp`) with `readModel` on a reader
`Highs`, copies the fields of the `lp` object below plus `integrality_` into
a fresh `HighsLp` (no names, column-wise), creates a fresh `Highs` with the
G1 base options (output off, `threads=1`, `parallel=off`, `random_seed=0`;
**presolve stays at its default** unless overridden), applies the scenario's
`overrides`, then calls `passModel` with that `HighsLp`. It writes the
`g1_header`, then performs the scenario's calls, writing one `g1_call` after
every API call and one `g1_result` after every `run`, `postsolve` and
`crossover`. The Rust replay does the same: options first, then
`pass_model(g1_lp)`, then the `g1_call`s in order. All G1 fixtures are
`target`-only.

Integers (`HighsInt`, enums, statuses) are JSON numbers, bools `0` / `1`,
doubles 16 hex digits of the bit pattern; arrays are written in full (no
hashes). `HighsStatus`: -1 kError, 0 kOk, 1 kWarning. `HighsVarType`: 0
continuous, 1 integer, 2 semi-continuous, 3 semi-integer, 4 implicit.

| `t` | when | fields |
|---|---|---|
| `g1_header` | first line | `name`, `highs_version`, `highs_githash`, `overrides` (`k=v,k=v` or `default`; keys as in `g1_options`, bools as `true` / `false`), `lp` (the `HighsLp` passed to `passModel`, encoding as `lp` above, unscaled, column-wise), `integrality` (only when non-empty) |
| `g1_call` | after each API call | `call` (0-based counter over the `g1_call` lines), `op`, `status` (the call's `HighsStatus`), then the op's inputs and outputs (below) |
| `g1_result` | after each `run` / `postsolve` / `crossover` call | `model_status`, `objective_function_value`, `primal_solution_status`, `dual_solution_status`, `simplex_iteration_count`, `max_primal_infeasibility`, `max_dual_infeasibility`, `value_valid`, `dual_valid`, `basis_valid` (`getBasis().valid`), `col_value`, `row_value`, `col_dual`, `row_dual` (`getSolution()`, as stored, possibly empty); when the incumbent is a MIP (`getLp().isMip()`): `mip_node_count`, `mip_dual_bound`, `mip_gap`, `max_integrality_violation`; when the basis is valid: `col_status`, `row_status` (`HighsBasisStatus` as int); always: `lp` + `integrality` (the incumbent `getLp()` after the call, column-wise copy) |
| `g1_file` | after the `run` of a `g1-mipfile` scenario | `text` (the whole `mip_improving_solution_file`, as a JSON string) |
| `clqn_result` | `g1-clique-normalise` only, after the call (1.16.0 oracle) | `result` (`HPresolve::Result` of `normaliseCliqueRows`), `num_reductions` (postsolve stack), `col_lower`, `col_upper`, `col_deleted`, `row_lower`, `row_upper`, `row_deleted` (the presolve model and flags after the call), `row_start` / `row_index` / `row_value` (every row's nonzeros in `getRowVector` order, none for a deleted row), `sum_lower`, `sum_upper`, `num_inf_sum_lower`, `num_inf_sum_upper` (`impliedRowBounds` of each row, 0 for a deleted row) |

`g1-clique-normalise` (`src/clqn.cpp`, HiGHS #3286 + #3318) does not run
`Highs::run`: after the `g1_header` it sets up a `HighsMipSolver` on the model as the
1.16.0 TestPresolveRules normalisation cases do (`mipdata_` init,
`setupDomainPropagation`), then `HPresolve::setInput(mipsolver, -1)`,
`okSetupPresolveDataStructures` and one `normaliseCliqueRows`, and writes one
`clqn_result`. Only scenario `default`; it fails with 64 on HiGHS < 1.16
(`gen_fixtures.sh` runs it only with `ORACLE_HIGHS=1.16.0`). Consumer:
`crates/rusolve-ekk/tests/t_g1_clique_normalise.rs`.

`g1_call` op fields (inputs, then outputs):

| `op` | inputs | outputs |
|---|---|---|
| `run` | — | — |
| `set_solution` | `col_value`, optional `row_value`, optional `col_dual` + `row_dual` | — |
| `set_solution_sparse` | `index`, `value` | — |
| `pass_linear_objectives` | `objectives`: list of `{weight, offset, coefficients, abs_tolerance, rel_tolerance, priority}` | — |
| `clear_linear_objectives` | — | — |
| `presolve` | — | `model_presolve_status` (`HighsPresolveStatus` as int: -1 kNotPresolved, 0 kNotReduced, 1 kInfeasible, 2 kUnboundedOrInfeasible, 3 kReduced, 4 kReducedToEmpty, 5 kTimeout, ...), `presolved` (`lp` + `integrality` of `getPresolvedLp()`; only when the status is not kError and the presolve status is kNotReduced, kReduced, kReducedToEmpty or kTimeout) |
| `postsolve`, `crossover` | the solution: `col_value`, `row_value`, optional `col_dual` + `row_dual` | — (a `g1_result` follows) |
| `get_standard_form_lp` | — | `num_col`, `num_row`, `num_nz`, `offset` (from the sizes-only call), `cost`, `rhs`, `start`, `index`, `value` (from the second call with arrays of those sizes; empty when the first call returned kError, and `status` is then the first call's) |
| `get_primal_ray`, `get_dual_ray`, `get_dual_unboundedness_direction` | `dim` (length of the zero-initialised array passed) | `has`, `value` (the array after the call, always written; compared only when `has`) |
| `get_fixed_lp` | — | `fixed` (`lp` + `integrality`; only when the status is not kError) |
| `get_ill_conditioning` | `constraint`, `method`, `ill_conditioning_bound` | `index`, `multiplier` (the `HighsIllConditioning::record` entries in order) |
| `get_saved_mip_solutions` | — | `objectives` (`HighsObjectiveSolution::objective` in order; `status` is 0) |
| `change_coeff` | `row`, `col`, `value` | — |
| `scale_col`, `scale_row` | `col` / `row`, `value` | — |
| `add_col` | `cost`, `lower`, `upper`, `index`, `value` | — |
| `add_row` | `lower`, `upper`, `index`, `value` | — |

A `HighsSolution` input (`set_solution`, `postsolve`, `crossover`) is built
exactly as `g1.rs::solution_from` builds it: `col_value` sets
`value_valid`, `row_value` is copied when present, `col_dual` + `row_dual`
set `dual_valid`; nothing else is set.

Models: netlib LPs from `benches/data/netlib/<name>.mps`; HiGHS check
instances from `$HIGHS_REF/check/instances/<name>.mps` (`CHECK_DIR`, not
copied); handcrafted models from `lps/g1/<name>.lp` or `.mps` (looked up in
that order, `gen_fixtures.sh` `g1_path`).

| dir / command | instances x scenarios | calls |
|---|---|---|
| `user_scale` / `g1-user-scale` | `afiro`, `sc50a`, `adlittle`, `blend` x `o3` (`user_objective_scale=3`), `b-2` (`user_bound_scale=-2`), `o-4b5` (`user_objective_scale=-4,user_bound_scale=5`), `o40` (`user_objective_scale=40`); `flugpl` x `o3b2` (`user_objective_scale=3,user_bound_scale=2`, MIP) | `run`, result |
| `inf_cost` / `g1-inf-cost` | `inf_cost_min`, `inf_cost_max`, `inf_cost_unbounded`, `inf_cost_infeasible`, `inf_cost_mip` x `default`, `ic1e15` (`infinite_cost=1e15`, so the 1e15 costs in the models are infinite too) | `run`, result |
| `semi` / `g1-semi` | `semi-continuous`, `semi-integer` (check), `semi_inf_upper`, `semi_inconsistent`, `semi_neg_lower` x `default`, `relax` (`solve_relaxation=true`) | `run`, result |
| `set_solution` / `g1-set-solution` | `afiro`, `adlittle`, `sc50a` x `full` (`set_solution` with `col_value`), `full-dual` (all four vectors), `sparse` (`set_solution_sparse` on columns 0, 2, 4, ...), `cold` (`use_warm_start=false`, `col_value`); `flugpl`, `egout`, `small_mip` x `mip-full` (`col_value`), `mip-partial` (`set_solution_sparse` on every other integer column), `mip-fractional` (`col_value` with the first integer column whose value + 0.5 <= upper moved by +0.5) | `set_solution*`, `run`, result |
| `multi_objective` / `g1-multi-objective` | `multiobj_lp`, `multiobj_mip`, `afiro` x `blend` (`blend_multi_objectives=true`), `lex` (`false`) | `pass_linear_objectives`, `run`, result, `clear_linear_objectives`, `run`, result |
| `presolve_api` / `g1-presolve-api` | `afiro`, `adlittle`, `sc50a`, `scagr7` x `presolve`, `postsolve`, `crossover`, `crossover-dual`; `egout`, `flugpl` x `presolve`, `postsolve` (MIP); `semi-continuous` x `presolve` (rejected: kError) | see below |
| `standard_form` / `g1-standard-form` | `afiro`, `adlittle`, `blend`, `sc50a`, `std_form` x `default` | `get_standard_form_lp`, `run`, result, `get_standard_form_lp` |
| `rays` / `g1-rays` | `ray_infeasible`, `ray_unbounded`, `ray_unbounded_dual` x `simplex` (`presolve=off`), `presolve` (`presolve=on`); `flugpl`, `egout`, `small_mip` x `fixed` | `simplex` / `presolve`: `run`, result, `get_dual_ray` (`dim` = rows), `get_dual_unboundedness_direction` (`dim` = columns), `get_primal_ray` (`dim` = columns); `fixed`: `run`, result, `get_fixed_lp` |
| `ill_conditioning` / `g1-ill-conditioning` | `afiro`, `sc50a`, `blend`, `ill_cond` x `m0c`, `m0v`, `m1c`, `m1v` (method 0 / 1, `constraint` 1 for `c`, 0 for `v`; bound 1e-4) | `run`, result, `get_ill_conditioning` |
| `modify` / `g1-modify` | `afiro`, `adlittle`, `sc50a` x `coeff`, `scale`, `add` | see below; every change is followed by `run`, result |
| `mipfile` / `g1-mipfile` | `flugpl`, `egout`, `small_mip` x `dense` (`mip_improving_solution_file=g1_improve.sol,mip_improving_solution_save=true`), `sparse` (plus `mip_improving_solution_report_sparse=true`) | `run`, result, `g1_file`, `get_saved_mip_solutions` |

Scenario details:
* **set_solution values** come from a first solve of the header model on a
  separate `Highs` (base options), so they are bit-identical inputs for both
  sides and the main `Highs` starts from `passModel` with no prior run (the
  same state as `run` then `clearSolver`); the record carries them in full.
* **multi_objective objectives.** `multiobj_lp`: (w 1, offset 0, c [1, 2, 0,
  1], abs 1e-6, rel 0, priority 2), (w -1, offset 1, c [0, 1, 1, 0], abs 0,
  rel 1e-4, priority 1), (w 0.5, offset -2, c [1, 0, 0, -1], abs 1e-6, rel
  1e-6, priority 0). `multiobj_mip`: (w 1, offset 0, c [-3, -2, -1], abs
  1e-6, rel 0, priority 1), (w 1, offset 0.5, c [0, 1, -1], abs 0, rel 1e-4,
  priority 0). `afiro`: objective 0 is the model cost (w 1, priority 1, abs
  1e-6, rel 0), objective 1 is all ones (w -0.5, priority 0, abs 0, rel
  1e-4), both offset 0. HiGHS leaves the incumbent costs as the blended cost
  (blend), or zeroes them and sets the offset to the last objective's with
  `objective_function_value` 0 and an invalid basis (lex,
  `returnFromLexicographicOptimization`), so the run after
  `clear_linear_objectives` solves that model (lex: objective = that
  offset, -2 for `multiobj_lp`, 0.5 for `multiobj_mip`); the `lp` of each
  result shows it.
* **presolve_api.** `presolve`: `presolve`, `run`, result. `postsolve`:
  `presolve`; a second `Highs` (base options, `presolve=off`) solves
  `getPresolvedLp()`; its `col_value` and `row_value` (plus `col_dual` and
  `row_dual` for an LP with `dual_valid`) go to `postsolve`, then result.
  `crossover`: run with `solver=ipm,run_crossover=off`, result, `crossover`
  with that run's `col_value` and `row_value` (IPX crossover from primal
  values only: kWarning, model status kUnknown), result. `crossover-dual`:
  the same with the duals as well (IPX `CrossoverFromStartingPoint` rejects
  them: kError, error flag 105).
* **modify.** `coeff`: `change_coeff` on column 0's first nonzero (value
  doubled), on column 1 at row r, the first r with column 1's r-th stored
  row index != r (value 0.5), on column 2's first nonzero (value 0: removal). `scale`: `scale_col`
  (0, 2), `scale_col` (1, -0.5), `scale_row` (0, -2). `add`: `add_col` (cost
  1, bounds [0, inf], rows 0, 1, 2 with 1, -1, 2), then `add_row` (bounds
  [-inf, 1000], columns 0, 1, 2 with 1, 1, 1).
* **mipfile file name.** The overrides say `g1_improve.sol`; the oracle
  writes it in a fresh `mkdtemp` directory, embeds the text and removes the
  file. The port writes under `std::env::temp_dir()` (`g1_temp_path`). The
  model has no names, so the file reads `NoName` for every column, as in the
  port.
* **ic1e15, not ic1e10.** `infinite_cost` has the range [1e15, inf], so
  `infinite_cost=1e10` is rejected by `setOptionValue`; the
  scenario uses the smallest legal value, which turns the models' 1e15 costs
  infinite.

Statuses HiGHS returns, kept as is (`status` / `model_status`):
`user_scale/{adlittle,blend,sc50a}.o40` run kError / kNotset (the scaled
bounds or costs leave the allowed range; `afiro.o40` solves);
`inf_cost/inf_cost_unbounded.ic1e15` run kError and
`inf_cost/inf_cost_infeasible.ic1e15` run kWarning, both kUnknown;
`presolve_api/*.crossover` crossover kWarning / kUnknown with the crossover
objective (afiro 440); `presolve_api/*.crossover-dual` crossover kError /
kNotset; `presolve_api/{egout,flugpl}.postsolve` (MIP) postsolve kWarning /
kUnknown; `presolve_api/semi-continuous.presolve` presolve kError with
`model_presolve_status` -1 (no `presolved`), the run then solves.

Non-finite values HiGHS reports, kept as is: `restoreInfCost` adds `value *
cost` for each infinite-cost column with a nonzero value, so
`objective_function_value` is `inf` for `inf_cost_min`, `inf_cost_max` and
`inf_cost_mip` under `ic1e15` (a 1e15-cost column fixed at a nonzero bound);
`assessCosts` stores those costs as `inf` in the incumbent `lp`.
`ill_conditioning_bound`, the bounds and the costs may be `inf` / `-inf` as
in any `lp`. No fixture holds a NaN.

## MANIFEST.tsv

`crates/rusolve-ekk/tests/fixtures/MANIFEST.tsv` (copied to
`target/highs-fixtures/`) lists every generated file: `file`, `bytes`,
`sha256`, `location` (`committed` = present in `tests/fixtures/`, `target` =
regenerate into `target/highs-fixtures/`). `#` lines carry the HiGHS commit and
the patch hashes (`patch_sha256`, `ipx_patch_sha256`, `w1_patch_sha256`, `w3_patch_sha256`, `w4_patch_sha256`, `w5_patch_sha256`). Tests that need a `target` fixture should skip (not fail) when
it is absent and may check its sha256 against the manifest. The W5 MIP
traces (`w5/mip/*.trace`) are not JSONL, and `tests/common/fixture.rs` loads
every row as JSONL, so they are listed as comment lines
`# trace <file> <bytes> <sha256>` (all `target`).
