# highs-oracle

A C++ golden-fixture generator for the `rusolve-ekk` faithful port of the
HiGHS simplex. It links an instrumented copy of HiGHS, calls HiGHS' internal
classes directly (`HighsSparseMatrix`, `HVector`, `HFactor`, `HighsRandom`,
`HEkk`, `HEkkDual`, `HEkkPrimal`, BASICLU and the `ipx::` classes), and writes
bit-exact JSONL fixtures. A Rust
translation matches HiGHS when it reproduces these files bit for bit. The
first diverging record shows where a translation went wrong.

HiGHS is MIT-licensed (© the HiGHS authors). This tool reads
`reference/HiGHS` but never modifies it. The instrumentation is a patch,
applied to a copy under `target/`.

- Pinned HiGHS: **1.15.1**, commit `73cac48c5340d775a477087198611862559be250`.
  `MANIFEST.tsv` records the commit and the six patch hashes.
- A second oracle on **HiGHS 1.16.0** (commit `f3f113d230429a8c0014b12ab0dd1a3847e58c55`)
  builds side by side for comparing the two versions (see *HiGHS 1.16.0 oracle*). It never writes the
  committed fixtures.
- Encoding: see [FORMAT.md](FORMAT.md).

## Build and regenerate

Run from the repository root (a worktree also works):

```sh
tools/highs-oracle/gen_fixtures.sh                     # build + regenerate everything
tools/highs-oracle/gen_fixtures.sh --determinism       # + regenerate twice more, byte-compare
tools/highs-oracle/gen_fixtures.sh --fp-contract-check # + rebuild with -ffp-contract=on, compare
tools/highs-oracle/gen_fixtures.sh --mip               # only the mip/ traces + .ac into target/highs-fixtures/mip
ORACLE_HIGHS=1.16.0 tools/highs-oracle/gen_fixtures.sh # the HiGHS 1.16.0 oracle (below)
```

**Requirements:**
- `cmake` ≥ 3.18
- a C++11 compiler (Apple clang and GCC 13 are both tested)
- `rsync`
- `patch`
- `shasum` or `sha256sum`
- `python3` (for `mkvariants.py`)

**Inputs:**
- `reference/HiGHS` and `benches/data/lp47` are gitignored. The script looks for them in the current checkout first, then in the main checkout (found through `git rev-parse --git-common-dir`).
- Override the locations with `HIGHS_REF=…` and `LP47_DIR=…`.
- Get the HiGHS tree with
  `git clone https://github.com/ERGO-Code/HiGHS reference/HiGHS && git -C reference/HiGHS checkout 73cac48c5340d775a477087198611862559be250`.
  `lp47` holds the LP relaxations of the pinned MIPLIB 2017 instances (`python3 tools/bench/fetch_corpus.py --set miplib-ref`; see `benches/data/ATTRIBUTION.md`).

**What the script does:**

1. **Patched copy.** It copies `reference/HiGHS` to `target/highs-oracle/highs-src` (without `.git`) and applies `highs-oracle-hooks.patch`, then `highs-oracle-ipx-hooks.patch`. Each patch has its own stamp file: `highs-src.stamp` holds the HiGHS commit and the first patch's sha256, and `highs-src.ipx.stamp` holds the second patch's sha256. The copy is refreshed when any of the three changes.
2. **Build.** It configures and builds `target/highs-oracle/build/highs-oracle`, plus the `mip_trace` and `api_seq` drivers in the same directory (`CMakeLists.txt` in this directory). The HiGHS build is:
   - static (`BUILD_SHARED_LIBS=OFF`);
   - serial (`HIGHS_NO_DEFAULT_THREADS=ON`);
   - `CMAKE_BUILD_TYPE=Release`, which uses HiGHS' own flags, `-O3 -DNDEBUG`;
   - IPO/LTO off;
   - `FAST_BUILD=ON`;
   - `ZLIB=OFF`;
   - no tests, examples, executable or extras;
   - **`-ffp-contract=off`** (see *Floating point*).

   It also configures a second build, `target/highs-oracle/build-hipo`, from the same patched copy with `-DORACLE_HIPO=ON`. Only the `w4/` fixtures come from it. In that build `CMakeLists.txt` sets `HIPO=ON` and builds `refblas/refblas.c`, a reference CBLAS written for rusolve (netlib loop order, `-O2 -ffp-contract=off`, MIT). It then presets HiGHS' `HIGHS_BLAS_*` variables (`HIGHS_BLAS_CONFIGURED`, `HIGHS_BLAS_TARGET=refblas`, vendor `rusolve refblas`), so `FindHipoDeps` searches for no BLAS and HiPO links refblas. `-DBLAS_LIBRARIES` would not pin it, because `FindBLAS` clears that variable. The script fails if the configure log lacks the pin or shows a BLAS search. It then builds and runs `test_refblas`, the hand-computed unit checks of refblas. METIS, AMD and RCM come from HiGHS' `extern/`. The default `build/` stays `HIPO=OFF`, so its fixtures do not change.
3. **Generate.** It writes every fixture to `target/highs-fixtures/`. The non-optimal variant LPs are first written to `target/highs-oracle/variants/` by `mkvariants.py`.
4. **Publish.** It copies the committed subset to `crates/rusolve-ekk/tests/fixtures/`. It also writes `MANIFEST.tsv` (file, bytes, sha256, `committed`/`target`) to both places.

A full run takes about 10 s on an M-series Mac after the one-off HiGHS build. On Linux aarch64 (6 jobs) the HiGHS build plus a full run takes about 45 s.

You can also run the oracle binary by hand:

```
highs-oracle unit    <mps> <name> <out.jsonl>
highs-oracle factor  <mps> <name> <out.jsonl>
highs-oracle basiclu <mps> <name> <out.jsonl>
highs-oracle basiclu-synth <dense70|rand120|neardep40> <out.jsonl>
highs-oracle simplex <mps> <name> <simplex_strategy> <out.jsonl>
highs-oracle random  <out.jsonl>
highs-oracle ipx-model  <mps> <name> <cfg> <out.jsonl>   # ipx::Model Load + transforms
highs-oracle ipx-linalg <mps> <name> <cfg> <out.jsonl>   # IPX linear-algebra components on real IPM iterates
highs-oracle ipx-trace  <mps> <name> <cfg> <out.jsonl>   # per-iteration IPX trace (patch 2)
highs-oracle ipx        <mps> <name> <out.jsonl>         # end-to-end IPX, 4 configurations (patch 2)
highs-oracle w1-qp      <model> <name> <out.jsonl>       # QUASS per-pass trace + result (patch 3)
highs-oracle w1-ranging <model> <name> <out.jsonl>       # getRanging after simplex
highs-oracle w1-iis     <model> <name> <out.jsonl>       # getIis, 4 iis_strategy configurations
highs-oracle w1-icrash  <model> <name> <out.jsonl>       # iCrash + crossover
highs-oracle w3-cupdlp  <mps> <name> <out.jsonl>         # cuPDLP-C trace + result (patch 4)
highs-oracle w3-hipdlp  <mps> <name> <out.jsonl>         # HiPDLP trace + result (patch 4)
highs-oracle w4-orderings <mps> <name> <metis|amd|rcm> <out.jsonl>  # HiPO, forced ordering (patch 5, build-hipo)
highs-oracle w4-hipo      <mps> <name> <out.jsonl>       # HiPO, its own choices (patch 5, build-hipo)
highs-oracle w5-simplex <mps> <name> <2|3> <threads> <out.jsonl>  # SIP / PAMI, threads=N, parallel=on
highs-oracle w5-hipo    <mps> <name> <tree|node|both> <choose|augmented|normaleq> <threads> <out.jsonl>
                                                          # HiPO, threads=N, parallel=on (patch 6, build-hipo)
highs-oracle g1-<x>     <model> <name> <scenario> <out.jsonl> # Highs API around run (no patch)
```

For `g1-<x>`, `<x>` is one of `user-scale`, `inf-cost`, `semi`, `set-solution`, `multi-objective`, `presolve-api`, `standard-form`, `rays`, `ill-conditioning`, `modify`, `mipfile`, `clique-normalise` (1.16 only); `<scenario>` selects the option overrides and API calls (FORMAT.md, *G1 records*).

For `ipx-model` and `ipx-linalg`, `<cfg>` is a comma-separated list of `ipx::Parameters` overrides (`key=value`), or `default`. For `ipx-trace`, it is a list of HiGHS option overrides.

## HiGHS 1.16.0 oracle (`ORACLE_HIGHS=1.16.0`)

`ORACLE_HIGHS=1.16.0` switches every path of the script; the default stays `1.15.1`.

| | 1.15.1 (default) | 1.16.0 |
|---|---|---|
| HiGHS tree (`HIGHS_REF`) | `reference/HiGHS` | `reference/HiGHS-1.16.0`, else `target/ref/HiGHS-1.16.0` (this checkout, then the main one) |
| patches | `highs-oracle-*.patch` | `patches-1.16.0/highs-oracle-*.patch` |
| build | `target/highs-oracle/` | `target/highs-oracle-1.16.0/` |
| fixtures + MANIFEST.tsv | `target/highs-fixtures/`, committed subset + MANIFEST in `crates/rusolve-ekk/tests/fixtures/` | `target/highs-fixtures-1.16.0/` only |

The 1.16 MANIFEST keeps the `committed`/`target` column of the 1.15.1 layout, so
`RUSOLVE_FIXTURES_DIR=target/highs-fixtures-1.16.0` (tests/common/fixture.rs) runs the
golden tests against it unchanged. `--mip` writes `<fixture dir>/mip/` (the committed
`mip/*.trace` / `.ac` set) and does not touch MANIFEST.tsv.

The 1.16 patch set is the 1.15.1 one re-cut against v1.16.0 (`git apply --check` clean,
no fuzz, in order), with these semantic changes:
- **hooks:** the `ac` record gains `centring_status` (FORMAT.md, *MIP IPM records*),
  written before 1.16's downgrade of a non-centred result to `kUnknown`.
- **w4:** the HiPO source is `FactorHighsSolver.cpp` (case-only rename). The
  `highsOracleW4Fact` hook follows `setParallelAfterSymbolic()` inside the new
  `if (kkt_.S.empty())` setup and reads `kkt_.S`; the iteration hook sits after
  `makeStep()` in the 1.16 `runIpm` loop.
- **w5:** the buffered `hipo_ord` records are written after the `chooseOrdering`
  block, which now has a parallel (`for_each`) and a serial branch.

Driver changes, all no-ops on 1.15.1 (its fixtures are unchanged):
- `src/ipx.cpp`: `ipx::Basis::Load` is gone in 1.16, so the `basis_load` record uses a
  verbatim copy of the 1.15.1 `Load` body (selected by `HIGHS_VERSION_MINOR`).
- `mip_trace.cpp` sets `parallel=off` (1.16 #3406 runs the MIP in parallel under
  `choose`); `key=value` arguments still override it (the `w5/mip` traces pass `parallel=on`).
- `lps/g1/inf_cost_*.lp` write infinite costs as `1e400`: the 1.16 LP reader rejects the
  `inf` keyword as a coefficient. 1.15.1 reads both forms to the same model (the `g1/inf_cost`
  fixtures are byte-identical).
- `CMakeLists.txt` compiles `cupdlp_solver.c` outside the unity C file (with `FAST_BUILD`
  the writer definitions of `HighsOracleW3Dump.h` were dropped by the include guard and
  `highs-oracle` failed to link on macOS); needs cmake 3.18.

## Options used

**Every oracle `Highs` instance (`oracleBaseOptions`):**
- `output_flag=false`
- `threads=1`
- `parallel=off`
- `presolve=off`
- `solver=simplex`
- `random_seed=0` (the HiGHS default, set explicitly)

The W5 commands and the `w5/mip` traces override `threads` and set `parallel=on` (FORMAT.md, *W5 records*).

**Simplex traces** add `simplex_strategy=1` (serial dual, `kSimplexStrategyDual`) or `simplex_strategy=4` (primal).
- Every other option keeps its HiGHS default.
- Each trace records the full option set as HEkk sees it, in `ekk_solve_enter.options`.
- Defaults that matter:
  - `simplex_scale_strategy=2`;
  - `simplex_dual_edge_weight_strategy=-1`, which picks steepest edge, `edge_weight_mode=2`;
  - `simplex_crash_strategy=0`;
  - `simplex_dualize_strategy=0` and `simplex_permute_strategy=0`;
  - `simplex_update_limit=5000`;
  - `time_limit=inf`;
  - `use_original_HFactor_logic=true`.

**lp47 instances** are MIP files. The oracle solves their LP relaxation: after `readModel` it calls `Highs::clearIntegrality()`.

**The factor fixtures' "mid" and "optimal" bases** come from separate dual runs:
- the optimal basis from a run without a limit;
- the mid basis from a run with `simplex_iteration_limit = optimal_iterations / 2`.

Both are snapshots of `HEkk::basis_.basicIndex_` taken when `HEkk::solve` returns.

## Fixture inventory

`MANIFEST.tsv` lists every file with its size and sha256. Committed files total about 7.9 MB; all files together total about 136 MB.

| file(s) | location | HiGHS internals exercised |
|---|---|---|
| `random.jsonl` | committed | `HighsRandom`: `initialise(seed)`, `integer()`, `integer(sup)`, `fraction()`, `closedFraction()`, `real(a,b)`, `bit()`, `shuffle` |
| `unit/{afiro,sc50a,kb2}.jsonl` | committed | `HighsSparseMatrix::createRowwise`, `createRowwisePartitioned` (NULL and mixed partition), `priceByColumn` (quad and double), `priceByRow` (row-wise and partitioned), `priceByRowWithSwitch` (hyper-sparse, dense-switch, `from_index>0`), `computeDot`, `collectAj` (with exact cancellation → `kHighsZero`), `update` (partitioned swap maintenance); `HVector::tight` (sparse and dense mode), `pack`, `reIndex`, `norm2`, `saxpy` |
| `unit/{adlittle,blend,share2b}.jsonl` | target | same |
| `factor/afiro.jsonl` | committed | `HFactor::setupGeneral` + `build` on slack, mid-solve, optimal and two rank-deficient bases: rank deficiency, `row/col/var_with_no_pivot`, permuted `basic_index`, pivot sequence (`refactor_info_`), every `getInvert()` array; `ftranCall` / `btranCall` (HVector and `std::vector` overloads) for unit, column, sparse and dense RHS at expected densities 1.0 / 0.1 / 0.01 (dense, sparse and hyper-sparse paths); `update` sequences (24 basis changes) under **FT, PF, MPF, APF** from slack and mid bases, each followed by FTRAN/BTRAN; `HFactor::rebuild` (HFactorRefactor) via `refactor_info_.use = true` |
| `factor/{sc50a,kb2,adlittle,blend,share2b}.jsonl` | target | same |
| `basiclu/sc50a.jsonl` | committed | BASICLU (the IPX LU kernel, `highs/ipm/basiclu`) through its object API on the slack, mid-solve, optimal and one rank-deficient basis of the scaled LP: `basiclu_obj_factorize` (status, reallocation sizes, the whole `xstore[0..512]` scalar block), `basiclu_obj_get_factors`, dense and sparse solves `'N'`/`'T'`; a Forrest-Tomlin sequence (10 updates: `solve_for_update` `'N'`/`'T'`, `basiclu_obj_update`, solves) and a refactorization of the updated basis on the used store; plus a long run from the slack basis (`slack_long`, m-1 updates, solves every 8th step) that reaches the L/U compression in `lu_update`. The object is zero-initialised (`basiclu_obj_initialize` mallocs; IPX's `BasicLu` and the Rust port use zeroed vectors) |
| `basiclu/{afiro,kb2,adlittle,blend,share2b}.jsonl` | target | same |
| `basiclu/{dense70,rand120}.jsonl` | committed | as `basiclu/sc50a.jsonl` on synthetic matrices (`basiclu-synth`): a dense 70 x 140 matrix (pivot columns with more than 64 off-diagonals take `lu_pivot_any`) and a random sparse 120 x 240 one (many file expansions and garbage collections); basis = the first m columns, 6 updates |
| `basiclu/neardep40.jsonl` | committed | as above on a dense 40 x 80 matrix whose basis columns 36..39 are columns 0..3 plus 1e-6 noise, factorised with the default and with IPX's strict absolute pivot tolerance 1e-3 (`lu_pivot`'s `lu_remove_col`, rank-deficient Markowitz exits, `BASICLU_WARNING_singular_matrix`); 3 updates each |
| `simplex/<lp>.{dual,primal}.jsonl` for afiro, sc50a, kb2, adlittle, blend, sc105, share2b | committed | end-to-end `Highs::run` → `solveLpSimplex` → `HEkk::solve` → `HEkkDual::solve` / `HEkkPrimal::solve` (see the trace list below) |
| `simplex/{afirocut,afiromax}.{dual,primal}.jsonl` | committed | as above, on non-optimal / maximisation variants (below) |
| `simplex/<lp>{cut,max}.{dual,primal}.jsonl` for adlittle, blend, kb2, sc105, sc50a, share2b | target | same |
| `ipx/model/{afiro,ipxmix}.{default,dual}.jsonl`, `ipx/model/ipxtall.default.jsonl` | committed | `ipx::Model::Load` with its own dualization decision, the AI / AIt / A matrices, b, c and bounds, norms, the scaled user model, flipped and boxed variables, equilibration, and `GetInfo`; also `PresolveStartingPoint`, `PresolveIPMStartingPoint`, `PostsolveInteriorSolution` and `EvaluateInteriorSolution` (FORMAT.md, *IPX records*) |
| `ipx/model/<lp>.default.jsonl` for every netlib LP and ipxdep; `ipx/model/{sc50a,kb2}.dual.jsonl` | target | same |
| `ipx/linalg/afiro.default`, `ipx/linalg/ipxmix.{default,lu1,dual}`, `ipx/linalg/ipxtall.default`, `ipx/linalg/ipxdep.{default,dual,slack}` (`.jsonl`) | committed | the IPX linear-algebra components on iterates of the real C++ IPM: Iterate (state, residuals, objectives, Postprocess); KKTSolverDiag (with NormalMatrix, DiagonalPrecond and ConjugateResiduals); StartingBasis and Basis (crash, dependent rows and columns, `SolveDense` aliased and not aliased, `MinSingularValue` / PowerMethod, `SymbolicInvert`, `RandomPermute`, `Load`, `GetLuFactors`); BasicLu and ForrestTomlin; KKTSolverBasis (SplittedNormalMatrix, maxvolume, DropPrimal/DropDual, basis updates); DropToComplementarity, `ComputeBasicSolution`, and the Model postsolve and evaluate paths |
| `ipx/linalg/<lp>.{default,lu1,dual,slack}.jsonl` for afiro, sc50a, kb2, adlittle, blend, share2b, ipxmix, ipxtall and ipxdep (the rest) | target | same |
| `w3/{cupdlp,hipdlp}/<lp>.jsonl` for afiro, sc50a, kb2, adlittle, blend, share2b, sc105, scagr7 | target | `Highs::run` with `solver=pdlp` (cuPDLP-C) and `solver=hipdlp` (HiPDLP), `pdlp_iteration_limit=2000`: every convergence check, restart and the termination, with the iterate hashes (patch 4; FORMAT.md, *W3 records*) |
| `w4/orderings/<lp>.{metis,amd,rcm}.jsonl` and `w4/hipo/<lp>.jsonl` for the same 8 LPs | target | `Highs::run` with `solver=hipo` on the `build-hipo` oracle (HiPO on refblas), once per forced `hipo_ordering` and once with `choose`: every ordering tried (with the full pattern and the permutation), the chosen factorisation, every IPM iteration with the iterate hashes, and `hipo::Info` (patch 5; FORMAT.md, *W4 records*) |
| `w5/simplex/<lp>_s{2,3}_t{1,2,4,8}.jsonl` for the same 8 LPs | target | the simplex trace with `simplex_strategy` 2 (SIP) or 3 (PAMI), `threads=N`, `parallel=on` (FORMAT.md, *W5 records*) |
| `w5/hipo/<lp>_{tree,node,both}_{as,ne}_t{1,2,4,8}.jsonl` and `w5/hipo/<lp>_{tree,node,both}_choose_t1.jsonl` for the same 8 LPs | target | the W4 HiPO trace on `build-hipo` with `threads=N`, `parallel=on`, `hipo_parallel_type` and `hipo_system` (augmented / normal equations; `choose` only at one thread) (patch 6) |
| `w5/mip/{air05,enlight_hard,pg,sp150x300d}.t{2,4,8}.trace` | target | `mip_trace <mps> threads=N parallel=on`, the `mip/<instance>.trace` format; listed in MANIFEST.tsv as `# trace` comment lines. See *Floating point* for the libstdc++ caveat |
| `g1/<dir>/<name>.<scenario>.jsonl`: `user_scale`, `inf_cost`, `semi`, `set_solution`, `multi_objective`, `presolve_api`, `standard_form`, `rays`, `ill_conditioning`, `modify`, `mipfile` (130 files) on netlib LPs, HiGHS check instances and the `lps/g1/` models | target | the Highs API around `run` through the public API (no patch): user objective / bound scaling, infinite costs, semi-variables, `setSolution`, multiple linear objectives, `presolve` / `postsolve` / `crossover`, `getStandardFormLp`, rays and `getFixedLp`, `getIllConditioning`, model modification, the improving-solution file; each call's status and outputs and the full incumbent state after every run (FORMAT.md, *G1 records*) |
| `ipx/trace/{afiro,ipxmix,ipxtall}.default.jsonl` | committed | oracle (b): per-iteration IPX trace of `Highs::run` with `solver=ipx`, written by the patch-2 hooks |
| `ipx/trace/<lp>.{default,dual}.jsonl` for every netlib LP and the synthetic LPs | target | same |
| `ipx/run/{afiro,ipxmix,ipxtall,ipxdep}.jsonl` | committed | oracle (b) + (c): `Highs::run` with IPX (crossover on and off, `solver=ipm`, dualized), the hook lines and the final HighsInfo / HighsSolution / HighsBasis (FORMAT.md, *ipx/run*) |
| `ipx/run/<lp>.jsonl` for every other netlib LP | target | same |
| `simplex/<lp>.{dual,primal}.jsonl` for the other 21 netlib instances (boeing2, bore3d, brandy, capri, e226, finnis, grow7, lotfi, scagr25, scagr7, scorpion, scrs8, scsd1, scsd6, sctap1, share1b, stair, standata, standmps, tuff, vtp.base) and the lp47 LP relaxations beasleyC3, binkar10_1, mc11, milo-v12-6-r2-40-1 | target | same |

Every netlib instance in `benches/data/netlib` (28) is traced with both strategies, and all 64 of those traces end `Optimal`.

**Synthetic IPX LPs** (`lps/`). These reach IPX paths that netlib does not:
- **`ipxmix`:** a free variable, boxed variables (one with a negative lower bound), a fixed variable, an upper-bound-only variable (which ipx::Model flips), a shifted lower bound, a ranged row, and a badly scaled row (power-of-two equilibration).
- **`ipxtall`:** 8 rows and 3 columns, so the default strategy dualizes it.
- **`ipxdep`:** two linearly dependent equality rows and two duplicate columns. It gives `dependent_rows` = 2 at default and `dependent_cols` = 2 with `dualize=1`.

**Variants** (`mkvariants.py`, from afiro, adlittle, sc50a, sc105, blend, kb2 and share2b; 28 more traces) reach the exits the stock instances never do:
- `<lp>cut` adds the row `c'x <= z* - (0.01|z*| + 1)` (z* = the published netlib optimum), so the LP is primal infeasible. The dual runs end `Infeasible` with a dual-ray proof. The primal runs end `Infeasible` in phase 1; HApp then re-solves the unscaled LP with dual simplex (forced Devex), so these traces have two `HEkk::solve` calls (`solve_call` 0 and 1).
- `<lp>max` adds `OBJSENSE MAX`. adlittle and blend become `Unbounded` (the dual run hands over to primal, `exit_algorithm` 1) and are re-solved unscaled with primal; the others stay `Optimal`.
- **Iteration counts:**
  - netlib: from 18 (afiro, primal) to 833 (scrs8, primal);
  - lp47: from 958 (beasleyC3, dual) to 6074 (milo, primal).
- **Dual BFRT:** traces include runs where it flips variables, for example grow7 (619 flips) and binkar10_1 (451).
- **Mixed run:** finnis's primal run ends in a dual cleanup, so one trace mixes both algorithms.
- **lp47 selection:** the four instances are the largest that solve their relaxation in under 2 s with presolve off. For comparison, the dual run takes about 0.2 s on beasleyC3 and 1.2 s on milo-v12-6-r2-40-1.

**What each simplex trace captures** (field lists in FORMAT.md):

- **`ekk_solve_enter`:**
  - the exact scaled LP HEkk receives: costs, bounds, matrix, and the scale factors `col`/`row`/`cost`;
  - every option value;
  - the simplex status flags;
  - the incoming basis, if there is one.
- **`ekk_solve_initialised`:**
  - the initial basis (logical here, because crash is off);
  - the random vectors `numTotRandomValue_`, `numTotPermutation_` and `numColPermutation_`;
  - the full work and base arrays after `initialiseForSolve`;
  - the partitioned row-wise matrix.
- **`dual_start` / `primal_start`:**
  - the phase chosen;
  - the edge-weight mode;
  - the tolerances;
  - the full arrays, including the initial DSE weights and any perturbed costs.
- **`invert`:** one record per INVERT, with rank deficiency, synthetic ticks, kernel statistics and a hash of the permuted `basicIndex_`.
- **`rebuild`:** one record per rebuild, with:
  - the **reason** (`rebuild_reason` captured on entry by a dedicated hook);
  - the phase;
  - the objective values;
  - the infeasibility counts;
  - the perturbation/shift flags;
  - synthetic ticks;
  - hashes of every work array.
- **`iter`:** one record per iteration (dual):
  - iteration count, phase, `rebuild_reason`;
  - `row_out`, `variable_out`, `variable_in`, `move_out`;
  - `delta_primal`, `theta_dual`, `theta_primal`, `alpha_col`, `alpha_row`, `numericalTrouble`;
  - the DSE weight of `row_out` (`computed_edge_weight`) and the weight stored after the update;
  - the BFRT flip count and the flipped variables;
  - `update_count`, `updated_dual_objective_value`, `total_synthetic_tick`;
  - the four density running averages;
  - hashes of `basicIndex_`, `workDual_`, `baseValue_` and `dual_edge_weight_`.

  Primal `iter` records carry the primal equivalents: `variable_in`, `move_in`, `row_out` (-1 = bound flip), `theta_*`, `value_in`, `alpha_*`, the edge weight of `variable_in`, and hashes.
- **`iter_none`:** an `iterate()` call that did not iterate, with the `rebuild_reason` it set.
- **`cleanup` / `phase_end`:** perturbation-removal cleanup, and phase switches with the phase decided next.
- **`ekk_solve_exit`:** model status, exit algorithm, the final basis, and the full final arrays.
- **`result`:** the `Highs::run` outcome:
  - status;
  - objective;
  - iterations;
  - the HighsBasis column and row statuses;
  - the unscaled primal and dual solution.

**Not captured, or only partly:**
- **Per-iteration full vectors.** `row_ep`, `col_aq`, `row_ap` and the full `workDual` are hashed per iteration, not dumped, to keep the files small. Full arrays appear only at solve, phase and rebuild boundaries (as hashes) and at solve start and exit (in full).
- **Dual-row internals.** The CHUZC candidate set, workGroup and BFRT group boundaries are not recorded. Only the final flips and the chosen pivot are.
- **Per-INVERT factor arrays inside the solve.** These would be large. The `factor/` fixtures cover HFactor in isolation, and each solve-time INVERT is identified by its tick and statistics.
- **Unexercised paths.** The second (unscaled) `HEkk::solve` is traced on 10 variant runs, but it never iterates there (it confirms the status from its start basis). The objective-bound third solve is never reached. The DualTasks and DualMulti strategies (2/3) are out of scope: they are parallel, and this build is serial.

## Patch summary (`highs-oracle-hooks.patch`)

The patch adds observer callbacks and nothing else. All hooks are null pointers unless the oracle installs them, and each one only reads solver state. Without hooks, the library behaves exactly like upstream HiGHS.

- **New header** `highs/simplex/HighsOracleHook.h`: a process-wide table of three function pointers (`ekk`, `dual`, `primal`), the event ids, and the macros `HIGHS_ORACLE_EKK`, `HIGHS_ORACLE_DUAL` and `HIGHS_ORACLE_PRIMAL`.
- **`HEkk.cpp`:**
  - hooks before and after `initialiseForSolve()` in `HEkk::solve`;
  - at the start of `returnFromEkkSolve`;
  - at the end of a successful `computeFactor`.
- **`HEkkDual.cpp`:**
  - before the major loop (solve start);
  - after `solvePhase1()` and `solvePhase2()`;
  - immediately before and after both `rebuild()` calls;
  - after both iterate-strategy switch blocks (that is, after every `iterate()`);
  - at the end of `cleanup()`.
- **`HEkkPrimal.cpp`:** the same events as the dual.
- **`mip/HighsOracleAcDump.h`, `HighsMipSolverData.cpp`, `HighsLpRelaxation.cpp`:** when
  `HIGHS_ORACLE_AC_DUMP` is set, append the analytic-centre and LP-fallback IPM results
  (see FORMAT.md, *MIP IPM records*). Nothing changes when it is unset.
- **`HEkkPrimal.h`, `HSimplexNla.h`:** `friend struct HighsOracleAccess;`, which gives read access to HEkkPrimal's private members and to `HSimplexNla::factor_`.

## Patch 2 (`highs-oracle-ipx-hooks.patch`, oracle (b))

The second patch is applied after the first and has its own stamp, `target/highs-oracle/highs-src.ipx.stamp`. MANIFEST.tsv records its hash in the `# ipx_patch_sha256` line. Like the first patch, it only observes: when `HIGHS_ORACLE_IPX_DUMP` is unset, the hooks return immediately.

- **New header** `highs/ipm/ipx/HighsOracleIpxDump.h`: inline writers for the `ipm` and `info` trace records (FORMAT.md, *ipx/trace*).
- **`ipm.cc`:**
  - after the first `PrintOutput` of `IPM::StartingPoint` (`start`);
  - after each `Driver` iteration (`iter`);
  - after each centring step (`centring`);
  - at the end of `Driver` (`driver_end`).
- **`lp_solver.cc`:** at the end of `LpSolver::Solve`, after crossover (`solve_end`).

The trace covers the whole IPX path that HiGHS runs, including the BasicLu-wrapper LU of KKTSolverBasis inside the IPM. The `ipm` records carry `basis_changes` and `updates_ipm`, and hash the iterate after every step.

## Patch 3 (`highs-oracle-w1-hooks.patch`)

Applied after patch 2, with its own stamp `target/highs-oracle/highs-src.w1.stamp`; MANIFEST.tsv records its hash in the `# w1_patch_sha256` line. When `HIGHS_ORACLE_QP_DUMP` is unset nothing is written.

- **New header** `highs/qpsolver/HighsOracleQpDump.h`: the `qp_iter` / `qp_end` writers (FORMAT.md, *W1 records*).
- **`gradient.hpp`:** a read-only `oracleGradient()` accessor, so the hook reads the stored gradient without the recompute in `getGradient()`.
- **`quass.cpp`:** one `qp_iter` at the end of every pass of the main loop, one `qp_end` after the post-loop `loginformation`.

## Patch 4 (`highs-oracle-w3-hooks.patch`)

Applied after patch 3, with its own stamp `target/highs-oracle/highs-src.w3.stamp`; MANIFEST.tsv records its hash in the `# w3_patch_sha256` line. When `HIGHS_ORACLE_PDLP_DUMP` is unset nothing is written. The oracle build compiles the cuPDLP-C and HiPDLP CPU paths only (`CUPDLP_GPU` and `HIGHS_GPU_LIB` forced off in `CMakeLists.txt`; `HConfig.h` defines `CUPDLP_CPU`).

- **New header** `highs/pdlp/HighsOracleW3Dump.h`: C-compatible writers for the `pdlp_*` and `hipdlp_*` records (FORMAT.md, *W3 records*). The dump `FILE*` and the `call` counter are `extern "C"` functions defined once, in `cupdlp_solver.c` (`HIGHS_ORACLE_W3_DUMP_IMPL`), so all translation units share them.
- **`cupdlp_solver.c`:** `pdlp_check` after `PDHG_Compute_Infeas_Residuals` in `PDHG_Solve`; `pdlp_end` right after the main loop.
- **`cupdlp_proj.c`:** `pdlp_restart` at the end of `PDHG_Restart_Iterate_GPU`.
- **`pdhg.cc`:** `hipdlp_check` after each `runConvergenceCheck`, `hipdlp_restart` at the end of the `do_restart` block, `hipdlp_end` in `solveReturn`.

## Patch 5 (`highs-oracle-w4-hooks.patch`)

Applied after patch 4, with its own stamp `target/highs-oracle/highs-src.w4.stamp`; MANIFEST.tsv records its hash in the `# w4_patch_sha256` line. When `HIGHS_ORACLE_HIPO_DUMP` is unset nothing is written, and HiPO is compiled only in the `build-hipo` configuration, so the default build's fixtures do not change.

- **New header** `highs/ipm/hipo/HighsOracleW4Dump.h`: the `hipo_*` writers (FORMAT.md, *W4 records*). HiPO is C++ only, so the dump `FILE*` and the `call` counter live in function-local statics of inline functions.
- **`FactorHiGHSSolver.cpp`:** `hipo_ord` after the "Finished %s for %s" log of each ordering in `chooseOrdering` (patch 6 moves the write after the `for_each`); `hipo_fact` after `setParallel()` in `setup`.
- **`Solver.cpp`:** `hipo_iter` after each `makeStep()` in `runIpm`.
- **`IpxWrapper.cpp`:** `hipo_end` after `hipo.getInfo()` in `solveHipo`.

## Patch 6 (`highs-oracle-w5-hooks.patch`)

Applied after patch 5, with its own stamp `target/highs-oracle/highs-src.w5.stamp`; MANIFEST.tsv records its hash in the `# w5_patch_sha256` line. It touches only the `hipo_ord` hook, so the default build is unchanged.

- **`FactorHiGHSSolver.cpp`:** in `chooseOrdering`, the `run_ordering_and_analyse` lambda can run on several workers at once (`highs::parallel::for_each` at `parallel=on`), so it no longer writes `hipo_ord`. It saves `permutations[i]` and `failure[i]` per ordering index, and after the `for_each` the calling thread writes the records in index order. With `parallel=off` this gives the bytes of patch 5, so the `w4/` fixtures are unchanged.
- **`HighsOracleW4Dump.h`:** comment only.

The simplex hooks (patch 1) already fire on the main thread after `iterate` / `iterateTasks` / `iterateMulti`. `HighsOracleAcDump.h` (the MIP `.ac` dump) is not thread-safe, so `gen_fixtures.sh` never sets `HIGHS_ORACLE_AC_DUMP` for the W5 MIP traces.

To change any patch:
1. Edit the files in `target/highs-oracle/highs-src`.
2. Regenerate the patch with `diff -u --label a/<f> --label b/<f> reference/HiGHS/<f> target/…/<f>` for each touched file (the new header is diffed against `/dev/null`).
3. Rerun the script. A changed patch hash invalidates the stamp and the fixtures' manifest header.

## Determinism

**Result:** `gen_fixtures.sh --determinism` generated all 247 files three times (the published run plus two more). All three runs were **byte-identical**. After patch 4 the same check, on Linux aarch64, covered all 319 files (including the 16 `w3/` traces): byte-identical again. After patch 5 it covered all 351 files, including the 32 `w4/` traces from `build-hipo`: byte-identical. The patch-5 hooks are inert in the default build, which reproduced every older MANIFEST sha unchanged. The G1 commands (no patch, default build) add the 130 `g1/` fixtures, which were checked byte-identical the same way. After patch 6 (Linux aarch64, GCC 13.3) it covered all 761 JSONL files plus the 12 `w5/mip` traces, including the 280 `w5/` JSONL traces at threads 1 to 8 with `parallel=on`: byte-identical. Every older MANIFEST row kept its sha256.

**Thread counts.** For every one of the 8 LPs, the SIP (`s2`) and PAMI (`s3`) traces at threads 2, 4 and 8 equal the threads = 1 trace record for record (ignoring the `threads` value in the header and in `ekk_solve_enter.options`). So do the HiPO traces for every `(lp, hipo_parallel_type, as|ne)`. The `choose_t1` HiPO traces equal `w4/hipo/<lp>.jsonl` after the header line. The MIP searches do depend on the thread count (`getMaxNumWorkers` grows with it): air05, enlight_hard and pg give different traces at 2, 4 and 8 threads (air05 at 4 and 8 is the same), and sp150x300d, solved at the root, gives the committed threads = 1 trace at every count.

**Sources of nondeterminism and how each is neutralised:**

| source | where | status |
|---|---|---|
| `HEkk::random_` (HighsRandom) | Seeded by `random_.initialise(options.random_seed)` in `HEkk::setSimplexOptions`. Used by `initialiseSimplexLpRandomVectors` (cost-perturbation values, CHUZC permutation `numTotPermutation_`, `numColPermutation_`), by `HEkkDualRHS::chooseNormal`/`chooseMultiGlobal` (random start of the CHUZR scan), and by `HEkkDual::correctDualInfeasibilities`. | Deterministic given the seed; `random_seed=0` is set explicitly. **Porting note:** `initialiseSimplexLpRandomVectors` runs **twice** on one continuing stream with no re-seed: once in `initialiseEkk` and again in `initialiseForSolve`. The vectors HEkk actually uses therefore come from the second draw. The traces dump those used values, and `random.jsonl` pins the generator. |
| `HEkkPrimal::random_` | Used only in PSE-weight debug checks. | Not on the solve path. |
| `HighsRandom random(1)` | `HEkk::factorSolveError`, `HSimplexNlaDebug`. | Fixed seed, so deterministic. |
| Wall-clock time limit | `HEkk::bailout` checks `timer_->read() > time_limit` only when `time_limit < inf`. `HFactor::build` checks a timeout only when `time_limit_` is finite, which only `HPresolve::removeDependentEquations` sets. | `time_limit=inf` (the default) and presolve is off, so neither check is ever armed. |
| Reinversion trigger | `kRebuildReasonSyntheticClockSaysInvert` compares `total_synthetic_tick_` against `build_synthetic_tick_`. | These are synthetic operation counts, not time, so deterministic. They are dumped (`synthetic_tick`, `*_synthetic_tick`) because the port must reproduce them to reinvert at the same iteration. |
| Threads | `HEkkDual::iterateTasks` / `iterateMulti` and the `HighsTaskExecutor`. | The library is built with `HIGHS_NO_DEFAULT_THREADS`, run with `threads=1` and `parallel=off`, and uses strategy 1 (serial dual). The W5 fixtures set `threads` and `parallel=on`: the simplex hooks fire on the main thread, patch 6 writes `hipo_ord` from the calling thread, `hipo_system=choose` (whose two tasks race) is used only at one thread, and the MIP `.ac` dump is off. |
| Debug-only switches | `debug_solve_report_` (keyed on `build_synthetic_tick` 445560), `analyse_*`, `highs_debug_level`. | Print only; off. |
| IPX | `ipx::Control` timers (`Elapsed`, `time_*` in `ipx::Info`), the `time_limit` check in `InterruptCheck`, and `RandomPermute`'s `std::default_random_engine` (symbolic_invert.cc). | `time_limit=inf`. The `time_*` fields are never dumped. The engine is default-seeded on every call, and `basis_ops.random_permute` pins its libc++ sequence. |
| Hash containers | `visited_basis_` (HighsHashTable over basis hashes). | Contents depend only on data, not addresses. |
| Uninitialised memory | Two dumped members can be uninitialised: `HEkkDual::allow_dual_steepest_edge_to_devex_switch` when the edge-weight mode is not steepest edge (finnis primal's dual cleanup, the Devex unscaled re-solves), and the `HEkkPrimal` pivot members (`variable_out`, `alpha_*`, …) in an `iter_none` before the first iteration of a primal solve (the unscaled re-solves). Neither is read by HiGHS. | The oracle omits them in those cases (FORMAT.md); before that, the max-variant traces differed between runs. Every dumped vector has a defined length. |

**Hash-equality check:** a Rust test can rerun `gen_fixtures.sh`, or simply compare its own output to the committed files, and check the sha256 values in `MANIFEST.tsv`.

## Floating point

**The oracle is built with `-ffp-contract=off`**, and this is required for bit-exactness. Apple clang (and clang ≥ 14 generally) defaults to `-ffp-contract=on`, which fuses `a*b + c` into one FMA within an expression on arm64 and on FMA-capable x86-64. GCC in GNU mode defaults to `fast`. GCC < 14 treats `on` as `off`, so `--fp-contract-check` passes `fast` to those compilers instead. Rust never contracts: `mul_add` is explicit. So a faithful translation reproduces only the non-contracted build.

**Measured effect.** `--fp-contract-check` builds the same oracle with `-ffp-contract=on`:
- 76 of 77 fixture files differ (only `random.jsonl`, which is integer-only, is unaffected);
- **35 of 64 simplex traces take a different pivot sequence**. The earliest divergence is at iteration 8 (boeing2, dual). Final iteration counts change, for example tuff primal 818 → 1052 and milo dual 5219 → 4636.

So contracted HiGHS (an ordinary `pip highspy`/Homebrew build on Apple Silicon) is **not** a bit-exact oracle for the port. Use this tool's fixtures, not a stock HiGHS build, for parity tests.

**Rust-side rules for parity:**
- No `mul_add`.
- No `-C target-cpu=native` assumptions that change FP semantics. Rust does not auto-fuse, so this is safe by default.
- Keep HiGHS' exact operation order: summation order, `(a*b)*c` versus `a*(b*c)`, and the `HighsCDouble` two-sum/two-product sequences used by the `quad_precision` PRICE.

**Other platform dependencies:**
- **libm.** Scaling (`HighsLpUtils.cpp`, `HSimplex.cpp`) uses `log`, `exp` and `pow`. The resulting scale factors are powers of two, and the exact scaled LP is dumped in `ekk_solve_enter`, so the simplex fixtures do not depend on the Rust port's libm. A port of *scaling itself* can be checked against `lp.scale` there.
- **Generation platform.** These fixtures were generated on macOS arm64 with Apple clang 21 and libc++. Regenerating on Linux aarch64 with GCC 13.3, libstdc++ and glibc gives byte-identical results for all 49 committed files. Two differences remain:
  - **`random_permute`.** `ipx/linalg` `basis_ops.random_permute` comes from `std::default_random_engine` and `uniform_int_distribution`. libstdc++ uses `minstd_rand0` and a different rejection scheme, so a stock libstdc++ build gives a different permutation. The Rust port mirrors libc++. When not built against libc++, `src/ipx.cpp` reproduces the libc++ sequence (`minstd_rand` with the independent-bits draw). Only the oracle's dump changes; the HiGHS library is untouched.
  - **`ipx/run/finnis.jsonl` (target-only).** The sha differs from the macOS one at the same size. Its xon, ipm and dual sections match the `ipx/trace` fixtures exactly. The likely cause is `info.mean_fill` (`std::pow(f, 1.0/num_factors)`, basis.cc): glibc `pow` is correctly rounded and Apple's is not always. This has not been confirmed against the macOS file. Committed files are unaffected, and the Rust test passes on Linux because it uses the platform `pow` too. The MANIFEST row for this file is platform-specific. So with the committed MANIFEST, a Linux-generated `target/highs-fixtures` makes `fixture_selftest` (`trace_hash_identities`, `decoders_and_builders_round_trip`) and `ipx_run` fail the sha check.
- **`mip_trace` (only the `w5/mip` traces are in `gen_fixtures.sh`).** Built against libstdc++, `mip_trace` reproduces `mip/sp150x300d.trace`. With presolve off it diverges from `mip/sp150x300d.presolve_off.trace` at line 17 (`cb nodes=0 lpit=217` versus 218 on macOS; final lp_iters 1057 versus 1126). Forcing correctly rounded libm results does not change this, so the cause is probably libstdc++ versus libc++ (`std::sort`/`nth_element` tie order, or unordered-container iteration). The Rust port matches the macOS trace on Linux, so regenerate the MIP traces on macOS only. At the default `threads=1`, the Linux (GCC 13.3, libstdc++) `mip_trace` built by `gen_fixtures.sh` also reproduces `mip/pg.trace`, but not `mip/air05.trace` (first difference at line 18; final nodes 19 against 21) or `mip/enlight_hard.trace` (line 19; nodes 2303 against 4012). The `w5/mip` traces of air05 and enlight_hard generated on Linux are therefore libstdc++ searches. Regenerate them on macOS (`gen_fixtures.sh --w5`) before using them as the port's oracle.

## Layout

```
tools/highs-oracle/
├── CMakeLists.txt              # builds patched HiGHS (static, serial) + the oracle
├── gen_fixtures.sh             # copy+patch, build, generate, publish, --determinism, --fp-contract-check, --w5, --mip; ORACLE_HIGHS
├── highs-oracle-hooks.patch    # observer hooks (applied to the copy, never to reference/)
├── highs-oracle-ipx-hooks.patch # patch 2: per-iteration IPX trace hooks (own stamp)
├── highs-oracle-w1-hooks.patch # patch 3: QUASS per-pass trace hooks (own stamp)
├── highs-oracle-w3-hooks.patch # patch 4: cuPDLP-C / HiPDLP trace hooks (own stamp)
├── highs-oracle-w4-hooks.patch # patch 5: HiPO trace hooks (own stamp)
├── highs-oracle-w5-hooks.patch # patch 6: hipo_ord written in ordering order after the parallel for_each (own stamp)
├── patches-1.16.0/             # the same six patches re-cut for HiGHS 1.16.0 (ORACLE_HIGHS=1.16.0)
├── refblas/                    # reference CBLAS for the HIPO=ON build (ours, MIT) + test_refblas.c
├── lps/                        # synthetic IPX LPs: ipxmix, ipxtall, ipxdep
│   ├── w1/                     # synthetic IIS LPs: iis_rows, iis_bnds
│   └── g1/                     # handcrafted G1 models: inf_cost_*, semi_*, std_form, multiobj_*, ray_*, ill_cond
├── FORMAT.md                   # fixture encoding
└── src/
    ├── main.cpp                # subcommand dispatch
    ├── oracle.h                # shared declarations
    ├── emit.{h,cpp}            # JSONL / hex-f64 / FNV-1a emitter
    ├── unit.cpp                # HighsSparseMatrix, HVector, HighsRandom fixtures
    ├── factor.cpp              # HFactor fixtures
    ├── basiclu.cpp             # BASICLU (IPX LU kernel) fixtures
    ├── ipx.cpp                 # ipx-model / ipx-linalg / ipx-trace
    ├── w1.cpp                  # w1-qp / w1-ranging / w1-iis / w1-icrash
    ├── w3.cpp                  # w3-cupdlp / w3-hipdlp
    ├── w4.cpp                  # w4-orderings / w4-hipo
    ├── w5.cpp                  # w5-simplex / w5-hipo (threads=N, parallel=on)
    ├── g1.cpp                  # g1-<x>: the Highs API around run
    └── simplex.cpp             # hooks + HEkk/HEkkDual/HEkkPrimal traces
```
