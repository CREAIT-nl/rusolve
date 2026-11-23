# LP repros

Standalone LP files (no integer markers) taken from MIP gate and perf runs.
Each one is a case for the LP lane.

## set_packing node-LP stall (from mip_perf, 2026-10-22)

The node LPs of `set_packing(m=150, n=400, density=0.03)` stall in
rusolve-lp. Every column is in [0,1]. Each file freezes one failing node's
column bounds, captured with a scratch-only `eprintln!` in bnb.rs. Seeds 1
and 2 then end `mip` with `Error: node LP failed: LP iteration limit`.

This is not a root-LP failure. The roots solve. The failing nodes are:

* s1: node 5, depth 4
* s2: node 17, depth 10

Binary used: e3e0b4f release `rusolve-cli`.

| file | `rusolve-cli solve f --csv --no-presolve` | HiGHS 1.15.1 |
|---|---|---|
| `set_packing_s1_node_stall.mps` (400×150, 131 cols fixed) | **IterationLimit** at p1 506 + p2 100,001 iters, 5.3 s, obj stuck at 738 | Optimal **746**, 103 iters (176 presolve off) |
| `set_packing_s2_node_stall.mps` (400×150, 165 cols fixed) | **IterationLimit** at p1 532 + p2 100,001 iters, 4.9 s, obj stuck at 654 | Optimal **660**, 101 iters (140 presolve off) |
| `set_packing_s1_root_slow.mps` (root LP, 400×150) | Optimal 758.333, p1 567 + p2 **10,021** iters, 0.9 s | Optimal 758.333, 220 iters (322 presolve off) |

* With presolve on, the results are the same: IterationLimit for both node
  LPs.
* The `mip` command runs node LPs with presolve off and otherwise uses
  `SimplexConfig::default()`. Today that config equals the CLI's.

The stderr log of each stall contains only these two lines before the limit:

```
rusolve simplex: Bland engaged after 25 stalled pivots
rusolve simplex: randomized deviation + drift correction after 0 stalled pivots
```

Lever probes, with the same command:

* `RUSOLVE_BLAND_AFTER=1000000000` (Bland never engages) solves every case:
  * s1 node: Optimal 746, p2 194 iters
  * s2 node: Optimal 660, p2 174 iters
  * s1 root: p2 313 iters instead of 10,021
* These levers still hit IterationLimit:
  * `RUSOLVE_NO_P2PERT=1`
  * `RUSOLVE_NO_SH=1`
  * `RUSOLVE_SCALE=1`

So once Bland mode engages on this degenerate packing LP, phase 2 never
leaves it. The probable cause is that Bland's rule interacts with the
deviation/drift-correction step. Regenerate:
`G.set_packing(m=150, n=400, density=0.03, seed=s)`, and make every column
continuous with the node bounds above.

## presolve false Optimal (from mip_diff, 2026-10-22)

`presolve_false_optimal_sp_s2.mps` is the root LP of `set_partition[infeasible]`
seed 2 (28 columns in [0,1], 11 equality rows). With presolve on,
`rusolve-cli solve` returns `Optimal 107` after 0 iterations: presolve removes
all 11 rows. With `--no-presolve` it returns `Infeasible`, which HiGHS
confirms. It still reproduces at e443fac.

**Fixed 2026-10-22 by 60d5395 (G16 phase-2 tie perturbation).** This is the
only crate commit between the repro binary (e3e0b4f) and master. With the same
command:

| file | result | phase-2 iters |
|---|---|---|
| s1 node | Optimal 746 | 302 |
| s2 node | Optimal 660 | 277 |
| s1 root | Optimal | 418 (was 10,021) |

`Bland engaged after 25 stalled pivots` still appears in the log, so Bland
mode is now escaped but its exit logic is unchanged. Iterations are 2.7–3×
HiGHS. Keep these files as regression cases.

## warm_reject_block_angular_s46.mps (C3, warm import rejects a round-tripped basis)

This is a MIP: `block_angular` seed 46 from `tools/diff/harness/generators.py`, with
15 rows, 8 of them equalities. Run it as
`rusolve-cli mip warm_reject_block_angular_s46.mps --no-presolve`. Every
cut-round re-solve comes back `rejected` (`reoptimize` returns `warm = false`
and solves cold):

    LPs: 0 warm ..., 1 cold, 4 fallbacks (0 unfit, 4 rejected, ...)

Root bound propagation fixes 3 columns (lb == ub). After that, `reoptimize`
also rejects the unchanged root LP's own exported optimal basis, and
`solve_with(Some(basis))` needs 15 phase-1 iterations. With the unpropagated
bounds every re-solve is warm.

## jobshop_wideM_s2 node LPs: Optimal with an infeasible x (C4)

Two node LPs of `jobshop_bigM[wideM]` seed 2: 61 columns, 100 model rows
plus 7 root cut rows (`CUT0`–`CUT6`), and node bounds frozen. They were
captured by a scratch-only hook in `RusolveEngine::solve_from` while running
`rusolve-cli mip <case> --no-warm` at 14d68ec + mip-rootinc. Rows are named
`R<i>` in row-form order (two-sided model rows split into L and G). Columns
are `C<j>`. There are no integer markers.

Check with
`rusolve-cli solve <file> --csv --no-presolve > x.csv`, then check the bound and
row violations of `x.csv` against `<file>`. The same result appears on the 14d68ec
binary and with presolve on.

| file | rusolve | worst bound viol | worst row viol | HiGHS 1.15 |
|---|---|---|---|---|
| `jobshop_wideM_s2_node_rowviol.mps` | Optimal 287.003315, p1 165 + p2 25 | 5.1e-5 (C33 = 1.000051, ub 1) | **42** (R12, a G big-M row) | Optimal 287.003315, 54 iters |
| `jobshop_wideM_s2_node_overshoot.mps` | Optimal 287.001443, p1 150 + p2 11 | **1.02e-4** (C44, ub 1) | 1.4e-14 | Optimal 287.001443, 48 iters |

* The objective matches HiGHS in both files, so only the reported x is
  wrong.
* With `RUSOLVE_NO_SH=1`, the rowviol file comes back feasible (row
  violation 1.4e-14) but still overshoots a bound by 6.0e-6. The overshoot
  file still overshoots by 5.8e-5.
* These levers change nothing: `RUSOLVE_NO_P2PERT`, `RUSOLVE_NO_P1PERT`,
  `RUSOLVE_BLAND_AFTER=1e9`.

In the whole `--no-warm` MIP run, 106 of about 1,100 node LPs came back
Optimal with a bound overshoot above 1e-6. Their values sit at
`1 + k·1e-6` (1.000051, 1.000063, …). With warm start on (the default)
none appear, because the warm path checks its point and these nodes solve
warm.

## warm_stall_market_split_s1: warm dual reopt stalls (C2)

The FIRST warm reopt that hits the MIP's iteration-limit floor on
`market_split` medium seed 1 (`m=3, n=20`; regenerate the MIP case with
`generators.market_split(m=3, n=20, seed=1)` →
`tools/diff/mip_cases/market_split_medium_s1.mps`). Captured at 9a63801 with a
scratch-only hook in `RusolveEngine::solve_warm`: the node LP (26 columns,
6 rows — 3 E, node bounds applied, no integer markers) as
`warm_stall_market_split_s1.mps`, plus the offered basis snapshot (and the
cutoff, `None` here) as `warm_stall_market_split_s1.basis`
(`C <j> <status>` / `R <i> <status>`, the `Basis` contract order).

Baseline run (release, 9a63801):

    rusolve-cli mip tools/diff/mip_cases/market_split_medium_s1.mps --no-presolve
    LPs: 8315 warm (492 cutoff, 857 infeasible re-checked, 0 overturned),
         1 cold, 64 fallbacks (0 unfit, 0 rejected, 64 iteration limit, 0 error);
         iterations 106779 warm / 8138 cold

The LP-level regression is `dual::tests::c2_warm_reopt_of_market_split_stall_lands_fast`
(`crates/rusolve-lp/src/dual.rs`): reoptimize of the captured LP+basis must
land Optimal WARM at the cold objective in < 100 iterations. At 9a63801 it
returns `IterationLimit` at the 200 cap.

## warm_stall_pg_cut{0,1,2}: warm reopt after root cut rounds (C7)

The first three warm reopts of `rusolve-cli mip benches/data/miplib/pg.mps`
(MIPLIB 2017 pg, presolve on) that returned `warm = false`. They were
captured at ff4725a with a scratch-only hook in `RusolveEngine::solve_warm`.
Each `.mps` is the engine's caller LP at that solve: the model rows
(two-sided rows split into L and G), followed by the cut rows of rounds
1–3. It has root bounds and no integer markers. Rows are named `R<i>` and
columns `C<j>`. Each `.basis` is the snapshot MIP offered, in the same
format as the market_split repro. Cut rows added since the snapshot are
already `Basic`.

All three are imported, then dropped at the dual engine's distress rung
(400 stalled visits) and re-solved cold. Mechanism (traced 2026-10-23): an
exact 4-visit cycle of flip-batch **farcross** events on boxed columns.
LP-level regression:
`dual::tests::c7_warm_reopt_of_pg_cut_rounds_lands_warm`.

