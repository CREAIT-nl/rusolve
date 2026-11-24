# Convex MIQP by LP-based outer approximation

A model with linear constraints, integer columns and a convex quadratic objective
`c'x + 1/2 x'Qx` (class `miqp`) is refused by default. With the lever `miqp_oa=true` it is solved
by outer approximation over the existing MIP search and checked by `rusolve check`. Source files:
`crates/rusolve-ekk/src/mip/oa.rs` (cuts, separator, oracle), `crates/rusolve-ekk/src/mip/lazy.rs`
(the exactness loop), `crates/rusolve-cli/src/miqp_solve.rs` (glue), `route.rs` (routing).

## Method

* **Epigraph.** The quadratic is replaced by one column `eta_k` per block of Q (objective cost 1;
  for a maximisation the sign flips and `eta_k <= 0`).
* **Gradient cuts.** For any point `xh`, `eta_k >= f_k(xh) + grad f_k(xh)'(x - xh)` is valid
  because `f_k` is convex (Duran and Grossmann 1986, Math. Programming 36). `xh` is projected on
  the bound box first, so that `g'xp = 2 f(xp)` and the right-hand side is exact. Tiny gradient
  coefficients are dropped using the global bounds (a valid relaxation); coefficients above 1e10
  are rejected; the violation threshold is twice `mip_feasibility_tolerance`.
* **Root tangent.** The continuous QP relaxation is solved first (HiPO, then QUASS). The tangent
  at its optimum makes the root LP value equal the QP value and keeps the LP bounded.
* **In the tree.** `OaSeparator` is a `MipHooks::separate` hook (Quesada and Grossmann 1992,
  Computers and Chemical Engineering 16): at each node LP point it adds the tangent of every block
  whose epigraph value is below the quadratic. Hooks work in the presolved column space, so the
  tree separator forces `presolve=off` (`miqp_oa_lazy_only=true` keeps presolve and relies on the
  lazy loop only).
* **Exact incumbents.** A MIP solution with `eta_k < f_k(x) - tol` is not a solution of the MIQP.
  `OaOracle` is a lazy-constraint oracle: it adds the tangent at that solution (plus up to 50
  remembered violated cuts) and the search is restarted. The accepted incumbent has the true
  objective up to `tol` per block. The loop has a round limit and the solve deadline; an optimal
  status without an accepted last candidate is reported as an iteration limit, never `Optimal`.
* **Disaggregation** (`miqp_oa_disaggregate`). One `eta_k` per connected component of Q's
  off-diagonal graph, so a separable Q gets sparse per-term tangents (the effect studied by
  Vielma, Dunning, Huchette and Lubin 2017, Math. Programming Computation 9).

## Scope and refusals

Convexity is decided by the existing Q analysis. A Q that is convex only after a diagonal shift,
a nonconvex Q, quadratic constraints, `--start`, `--basis` and `--add-rows` with the lever stay
refused (the cuts would be invalid or the epigraph columns would change the model).

## Output

The envelope reports class `miqp`, engine `oa`, convexity, and a `miqp_oa` object: blocks,
disaggregation, rounds, cuts by origin, the root QP, the exact objective and the largest epigraph
violation at the incumbent (`docs/schema/rusolve-result-1.schema.json`).

## Known limits

The exactness loop re-solves from scratch each round; there is no incumbent polish with the
continuous QP at fixed integers; no perspective cuts for semicontinuous structure; nonconvex
binary quadratics are not convexified.
