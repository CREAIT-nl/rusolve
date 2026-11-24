# Mixed-integer SOCP and convex MIQCQP by LP-based outer approximation

A model with integer columns and second-order cones, rotated second-order cones or convex
quadratic rows (class `mi_conic`) is refused by default. With the lever `misocp_oa=true` it is
solved by outer approximation over the existing MIP search and checked by `rusolve check`. This
extends the convex-MIQP route (`docs/design/miqp-oa.md`) from a quadratic objective to cone
constraints. Source files: `crates/rusolve-ekk/src/mip/oa.rs` (cone block type, cuts, separator,
oracle, polish hook), `crates/rusolve-ekk/src/mip/lazy.rs` (the exactness loop),
`crates/rusolve-cli/src/misocp_solve.rs` (model build, root, polish, answer), `route.rs`
(routing). Options: `misocp_oa`, `misocp_oa_disaggregate`, `misocp_oa_no_polish`,
`misocp_oa_lazy_only`, all default off.

## Method

* **Everything is a cone.** Each user cone, each convex quadratic row side and the objective
  Hessian becomes a constraint `||y(x)|| <= t(x)` with `y`, `t` affine. A rotated cone
  `2 x0 x1 >= ||w||^2` is `||(x0 - x1, sqrt2 w)|| <= x0 + x1`. A quadratic row side
  `1/2 x'Qx <= u(x)` with `Q = L L'` (the dense factor the convexity test already computes) is
  `||(u - 1, sqrt2 L'x)|| <= u + 1` (Lobo, Vandenberghe, Boyd and Lebret 1998, Linear Algebra and
  its Applications 284). A convex quadratic objective gets an epigraph column `eta` with the cone
  `1/2 x'Qx <= eta`. The LP keeps the linear rows, the bounds and the integrality; the cones exist
  only as rows added as cuts.
* **Gradient cuts.** For a point `zh` with `y(zh) != 0`, the row `u'y(x) <= t(x)` with
  `u = y(zh) / ||y(zh)||` is valid for every point of the cone, because `u'y <= ||y|| <= t`
  (Cauchy-Schwarz); it is the first-order outer approximation of Duran and Grossmann 1986
  (Math. Programming 36). It is tight at `zh` when `zh` is on the cone's boundary.
* **The apex.** At `y(zh) = 0` the norm is not differentiable. Any unit `u` is a subgradient, so
  the cut with `u = e_0` is used: it is valid (it is the case `|y_0| <= t`), and it cuts the point
  off exactly when `t(zh) < 0`. Every cone also starts with the 2 m rows for `u = +-e_i`
  (`|y_i| <= t`, at most 32 components), which bound `t` from below so the first LP is not
  unbounded. A cut with a coefficient above 1e10, or an empty one, is dropped and counted
  (`dropped_unsafe`); tiny coefficients are removed with the global column bounds, which keeps
  the row valid.
* **Root.** The continuous SOCP is solved by Clarabel (Goulart and Chen 2024, arXiv:2405.12762;
  the route of `conic_solve.rs`). The tangent of every cone at its optimum goes into the LP. A
  proved infeasible relaxation ends the run as `Infeasible`; an unbounded one is `Unknown`.
* **In the tree.** `OaSeparator` is a `MipHooks::separate` hook (Quesada and Grossmann 1992,
  Computers and Chemical Engineering 16): at each node LP point it adds the tangent of every
  violated cone, most violated first, up to 300 per call. Hooks work in the presolved column
  space, so tree separation forces `presolve=off`; `misocp_oa_lazy_only=true` keeps presolve and
  cuts only at the root and at candidate incumbents.
* **Exact incumbents.** A MIP solution whose cones are violated is not a solution. `OaOracle` is
  a lazy-constraint oracle: it adds the tangents at that solution and the search restarts. The
  loop has a round limit and the solve deadline.
* **Primal polish.** For each new integer assignment among the candidates, the integers are fixed
  and the continuous SOCP is solved by Clarabel. A result the checker accepts becomes an
  incumbent candidate and is handed to the next solve as a start. The assignment of the accepted
  OA point is polished once more at the end: the cuts leave the cones feasible only within
  `2 * mip_feasibility_tolerance`, which near an apex moves the objective by about its square
  root (a unit disc gave 7e-4 with the default tolerance), and the polished point is the exact
  one. The engine sets `mip_feasibility_tolerance=1e-9` for the same reason (`--opt` overrides
  it). `misocp_oa_no_polish=true` turns the polish off.
* **Disaggregation** (`misocp_oa_disaggregate`). A cone with three or more `y` components is
  written in the extended formulation of Vielma, Dunning, Huchette and Lubin 2017 (Math.
  Programming Computation 9): columns `z_i >= 0` with `sum z_i = t` and, for each `i`,
  `y_i^2 <= z_i t` written as the three-column cone `||(y_i, (z_i - t)/2)|| <= (z_i + t)/2`. The
  tangents of the small cones are sparse.

## Exactness

The answer is `Optimal` only when all of these hold: the last candidate was accepted by the
oracle; the reported point (the best of the polished points and the accepted candidate, with the
integers rounded) passes the same checker as `rusolve check` at its default tolerance, on the
original cones and rows; and the bound closes to the incumbent within `mip_rel_gap`. A bound that
lies beyond a checker-feasible point is wrong by definition: the status becomes `Unknown` and
`misocp_oa.bound_inconsistent` is true. A candidate the checker rejects is never reported. At a
time limit the status is `TimeLimit` with the best checked point and the bound.

## Scope and refusals

Cones of kind SOC and rotated SOC, and quadratic rows that are convex on each finite side, are
taken. A nonconvex row, a quadratic equality row, an indefinite objective Hessian, an exponential
or power cone, general constraints (indicator, SOS, piecewise-linear, min/max), and `--start`,
`--basis` and `--add-rows` stay refused or usage errors. With the lever off, behaviour is
unchanged: the typed refusal.

## Output

The envelope reports class `mi_conic`, engine `oa`, convexity, and a `misocp_oa` object: cone
counts, disaggregation, the Clarabel root, rounds, cuts by origin, the polish counters, the
origin of the incumbent, and the checker verdicts
(`docs/schema/rusolve-result-1.schema.json`).

## Known limits

The exactness loop re-solves the MIP from scratch each round, which dominates the time on
instances that need many rounds. The tree separator needs `presolve=off`. Cone cuts use a single
gradient per violated cone and node; there is no cut selection beyond the per-call limit.
