//! Primal-dual interior-point LP solver: infeasible-start
//! Mehrotra predictor-corrector with Gondzio correctors, normal equations by
//! supernodal sparse Cholesky with AMD, no crossover. Opt-in
//! (`rusolve solve --ipm` / `RUSOLVE_IPM=1`); the default simplex path does
//! not touch this module.
//!
//! Sources: Wright, *Primal-Dual Interior-Point Methods* (SIAM 1997);
//! Mehrotra, SIAM J. Optim. 2(4) (1992); Gondzio, COAP 6 (1996); HiGHS
//! HiPO and IPX (MIT). Each submodule header cites the sources it follows.

pub mod amd;
pub mod chol;
pub mod kernels;
pub mod model;
pub mod normal;
pub mod solver;

use rusolve_core::SparseLp;

use crate::simplex::{SimplexConfig, SimplexStatus};
pub use model::Verify;
pub use solver::{IpmConfig, IpmStats, IpmStatus};

/// IPM result in the caller's space.
#[derive(Debug, Clone)]
pub struct IpmSolution {
    pub status: IpmStatus,
    /// cᵀx on the caller's LP (NaN unless a point is reported; ±inf for a
    /// certified dual infeasibility).
    pub objective: f64,
    pub x: Vec<f64>,
    pub row_duals: Vec<f64>,
    pub reduced_costs: Vec<f64>,
    /// Duals are in the caller convention d = c − Aᵀπ (presolve-OFF path).
    pub duals_valid: bool,
    pub iterations: usize,
    pub stats: IpmStats,
    /// Original-LP check of the reported point.
    pub verify: Option<Verify>,
    /// (rows, cols, nnz) of the LP the IPM ran on (after presolve).
    pub solved_dims: (usize, usize, usize),
}

impl IpmSolution {
    fn bare(status: IpmStatus, maximize: bool) -> Self {
        Self {
            status,
            objective: match status {
                IpmStatus::DualInfeasible if maximize => f64::INFINITY,
                IpmStatus::DualInfeasible => f64::NEG_INFINITY,
                _ => f64::NAN,
            },
            x: vec![],
            row_duals: vec![],
            reduced_costs: vec![],
            duals_valid: false,
            iterations: 0,
            stats: IpmStats::default(),
            verify: None,
            solved_dims: (0, 0, 0),
        }
    }
}

/// IPM on `sp` as given (no presolve). Duals are reported.
pub fn solve_sparse(sp: &SparseLp, cfg: &IpmConfig) -> Result<IpmSolution, String> {
    let first = solve_sparse_once(sp, cfg)?;
    let Some(reg) = cfg.retry_reg else {
        return Ok(first);
    };
    if !matches!(first.status, IpmStatus::NoProgress | IpmStatus::Imprecise)
        || (cfg.rho >= reg && cfg.delta >= reg)
    {
        return Ok(first);
    }
    // One retry with stronger static regularization: Θ ≤ 1/ρ and M ⪰ δI
    // keep the Cholesky pivots above the cancellation floor that breaks the
    // default run (Saunders & Tomlin 1996 / Friedlander & Orban 2012
    // primal-dual regularization).
    if cfg.log {
        eprintln!(
            "rusolve ipm: {:?} with rho {:.0e} delta {:.0e}; retrying with rho = delta = {reg:.0e}",
            first.status, cfg.rho, cfg.delta
        );
    }
    let cfg2 = IpmConfig {
        rho: cfg.rho.max(reg),
        delta: cfg.delta.max(reg),
        ..cfg.clone()
    };
    let mut second = solve_sparse_once(sp, &cfg2)?;
    second.iterations += first.iterations;
    second.stats.iterations += first.stats.iterations;
    second.stats.t_factor += first.stats.t_factor;
    second.stats.t_solve += first.stats.t_solve;
    second.stats.t_ipm += first.stats.t_ipm;
    second.stats.retried = true;
    if second.status == IpmStatus::Optimal || first.status == IpmStatus::NoProgress {
        Ok(second)
    } else {
        Ok(first)
    }
}

fn solve_sparse_once(sp: &SparseLp, cfg: &IpmConfig) -> Result<IpmSolution, String> {
    let mut lp = match model::formulate(sp) {
        Ok(lp) => lp,
        Err(model::FormError::Infeasible(msg)) => {
            if cfg.log {
                eprintln!("rusolve ipm: infeasible at formulation: {msg}");
            }
            return Ok(IpmSolution::bare(IpmStatus::PrimalInfeasible, sp.maximize));
        }
        Err(model::FormError::Unsupported(msg)) => return Err(msg),
    };
    let dims = (lp.m, lp.n, lp.a.nnz());
    let core = solver::run_ipm(&mut lp, cfg);
    let mut sol = IpmSolution::bare(core.status, sp.maximize);
    sol.iterations = core.stats.iterations;
    sol.stats = core.stats;
    sol.solved_dims = dims;
    if core.status != IpmStatus::Optimal {
        return Ok(sol);
    }
    let (x, pi, d) = model::postprocess(sp, &lp, &core.point.x, &core.point.y);
    let v = model::verify(sp, &x, Some(&pi));
    if cfg.log {
        eprintln!(
            "rusolve ipm: original-LP check pinf {:.2e} dinf {:.2e} gap {:.2e} (tol {:.0e})",
            v.pinf, v.dinf, v.gap, cfg.verify_tol
        );
    }
    sol.objective = v.pobj;
    sol.status = if v.ok(cfg.verify_tol) {
        IpmStatus::Optimal
    } else {
        IpmStatus::Imprecise
    };
    sol.duals_valid = sol.status == IpmStatus::Optimal;
    sol.x = x;
    sol.row_duals = pi;
    sol.reduced_costs = d;
    sol.verify = Some(v);
    Ok(sol)
}

/// Presolve (with `ps`'s knobs, the same pipeline the simplex runs), IPM on
/// the reduced LP, primal postsolve, and the original-LP primal check.
/// Duals are not postsolved (`duals_valid = false`, as for the simplex).
pub fn solve_presolved(
    sp: &SparseLp,
    cfg: &IpmConfig,
    ps: &SimplexConfig,
) -> Result<IpmSolution, String> {
    let lp = sp.to_lp();
    match crate::presolve::run(
        &lp,
        ps.presolve_passes,
        ps.presolve_implied_eq,
        ps.presolve_parcol,
        ps.presolve_agg,
        ps.presolve_parcol_dom,
        ps.presolve_vub_eq,
    ) {
        crate::presolve::Presolved::Terminal { status } => {
            let st = match status {
                SimplexStatus::Infeasible => IpmStatus::PrimalInfeasible,
                SimplexStatus::Unbounded => IpmStatus::DualInfeasible,
                _ => IpmStatus::NoProgress,
            };
            Ok(IpmSolution::bare(st, sp.maximize))
        }
        crate::presolve::Presolved::Reduced {
            lp: reduced,
            live,
            postsolve,
        } => {
            if cfg.log {
                eprintln!(
                    "rusolve ipm presolve: rows {}->{} cols {}->{} nnz {}->{}",
                    lp.num_row,
                    reduced.num_row,
                    lp.num_col,
                    reduced.num_col,
                    lp.a.nnz(),
                    reduced.a.nnz()
                );
            }
            let rsp = SparseLp::from_lp(&reduced);
            let mut sol = solve_sparse(&rsp, cfg)?;
            if sol.status != IpmStatus::Optimal {
                // Reduced-LP point / duals are not the caller's.
                let mut out = IpmSolution::bare(sol.status, sp.maximize);
                out.iterations = sol.iterations;
                out.stats = sol.stats;
                out.solved_dims = sol.solved_dims;
                return Ok(out);
            }
            let (x, _obj) = crate::presolve::restore(&postsolve, &live, &sol.x, sol.objective);
            let v = model::verify(sp, &x, None);
            if cfg.log {
                eprintln!(
                    "rusolve ipm: original-LP primal check after postsolve pinf {:.2e} (tol {:.0e})",
                    v.pinf, cfg.verify_tol
                );
            }
            sol.status = if v.ok(cfg.verify_tol) {
                IpmStatus::Optimal
            } else {
                IpmStatus::Imprecise
            };
            sol.objective = v.pobj;
            sol.x = x;
            sol.row_duals = vec![];
            sol.reduced_costs = vec![];
            sol.duals_valid = false;
            sol.verify = Some(v);
            Ok(sol)
        }
    }
}

#[cfg(test)]
mod tests;
