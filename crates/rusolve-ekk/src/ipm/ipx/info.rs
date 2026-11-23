//! HiGHS source: ipm/ipx/info.h, info.cc
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use super::control::{fix2, sci2, sci8, textline};
use super::*;

/// HiGHS: ipx/info.cc::dump (lines 9-11)
///
/// `os << Textline(name) << value << '\n'` with default ostream state
/// (`Int` decimal, `double` `%g` precision 6, see `control::ostream_double`).
/// `value` is the already formatted text (`sci2(..)`, `fix2(..)` or an integer).
fn dump(os: &mut String, name: &str, value: &str) {
    os.push_str(&textline(&format!("info.{}", name)));
    os.push_str(value);
    os.push('\n');
}

/// HiGHS: ipx/info.cc::operator<<(std::ostream&, const Info&) (lines 13-98)
pub fn write_info(os: &mut String, info: &Info) {
    dump(os, "status", &info.status.to_string());
    dump(os, "status_ipm", &info.status_ipm.to_string());
    dump(os, "status_crossover", &info.status_crossover.to_string());
    dump(os, "errflag", &info.errflag.to_string());

    dump(os, "num_var", &info.num_var.to_string());
    dump(os, "num_constr", &info.num_constr.to_string());
    dump(os, "num_entries", &info.num_entries.to_string());

    dump(os, "num_rows_solver", &info.num_rows_solver.to_string());
    dump(os, "num_cols_solver", &info.num_cols_solver.to_string());
    dump(
        os,
        "num_entries_solver",
        &info.num_entries_solver.to_string(),
    );

    dump(os, "dualized", &info.dualized.to_string());
    dump(os, "dense_cols", &info.dense_cols.to_string());

    dump(os, "dependent_rows", &info.dependent_rows.to_string());
    dump(os, "dependent_cols", &info.dependent_cols.to_string());
    dump(os, "rows_inconsistent", &info.rows_inconsistent.to_string());
    dump(os, "cols_inconsistent", &info.cols_inconsistent.to_string());
    dump(os, "primal_dropped", &info.primal_dropped.to_string());
    dump(os, "dual_dropped", &info.dual_dropped.to_string());

    dump(os, "abs_presidual", &sci2(info.abs_presidual));
    dump(os, "abs_dresidual", &sci2(info.abs_dresidual));
    dump(os, "rel_presidual", &sci2(info.rel_presidual));
    dump(os, "rel_dresidual", &sci2(info.rel_dresidual));
    dump(os, "pobjval", &sci8(info.pobjval));
    dump(os, "dobjval", &sci8(info.dobjval));
    dump(os, "rel_objgap", &sci2(info.rel_objgap));
    dump(os, "complementarity", &sci2(info.complementarity));
    dump(os, "normx", &sci2(info.normx));
    dump(os, "normy", &sci2(info.normy));
    dump(os, "normz", &sci2(info.normz));

    dump(os, "objval", &sci8(info.objval));
    dump(os, "primal_infeas", &sci2(info.primal_infeas));
    dump(os, "dual_infeas", &sci2(info.dual_infeas));

    dump(os, "iter", &info.iter.to_string());
    dump(os, "kktiter1", &info.kktiter1.to_string());
    dump(os, "kktiter2", &info.kktiter2.to_string());
    dump(os, "basis_repairs", &info.basis_repairs.to_string());
    dump(os, "updates_start", &info.updates_start.to_string());
    dump(os, "updates_ipm", &info.updates_ipm.to_string());
    dump(os, "updates_crossover", &info.updates_crossover.to_string());

    dump(os, "time_total", &fix2(info.time_total));
    dump(os, "time_ipm1", &fix2(info.time_ipm1));
    dump(os, "time_ipm2", &fix2(info.time_ipm2));
    dump(os, "time_starting_basis", &fix2(info.time_starting_basis));
    dump(os, "time_crossover", &fix2(info.time_crossover));

    dump(os, "time_kkt_factorize", &fix2(info.time_kkt_factorize));
    dump(os, "time_kkt_solve", &fix2(info.time_kkt_solve));
    dump(os, "time_maxvol", &fix2(info.time_maxvol));
    dump(os, "time_cr1", &fix2(info.time_cr1));
    dump(os, "time_cr1_AAt", &fix2(info.time_cr1_AAt));
    dump(os, "time_cr1_pre", &fix2(info.time_cr1_pre));
    dump(os, "time_cr2", &fix2(info.time_cr2));
    dump(os, "time_cr2_NNt", &fix2(info.time_cr2_NNt));
    dump(os, "time_cr2_B", &fix2(info.time_cr2_B));
    dump(os, "time_cr2_Bt", &fix2(info.time_cr2_Bt));

    dump(os, "ftran_sparse", &fix2(info.ftran_sparse));
    dump(os, "btran_sparse", &fix2(info.btran_sparse));
    dump(os, "time_ftran", &fix2(info.time_ftran));
    dump(os, "time_btran", &fix2(info.time_btran));
    dump(os, "time_lu_invert", &fix2(info.time_lu_invert));
    dump(os, "time_lu_update", &fix2(info.time_lu_update));
    dump(os, "mean_fill", &fix2(info.mean_fill));
    dump(os, "max_fill", &fix2(info.max_fill));
    dump(os, "time_symb_invert", &fix2(info.time_symb_invert));

    dump(os, "maxvol_updates", &info.maxvol_updates.to_string());
    dump(os, "maxvol_skipped", &info.maxvol_skipped.to_string());
    dump(os, "maxvol_passes", &info.maxvol_passes.to_string());
    dump(os, "tbl_nnz", &info.tbl_nnz.to_string());
    dump(os, "tbl_max", &sci2(info.tbl_max));
    dump(os, "frobnorm_squared", &sci2(info.frobnorm_squared));
    dump(os, "lambdamax", &sci2(info.lambdamax));
    dump(os, "volume_increase", &sci2(info.volume_increase));
}

/// HiGHS: ipx/info.cc::StatusString (lines 100-122)
pub fn status_string(status: Int) -> String {
    let name = match status {
        IPX_STATUS_not_run => "not run",
        IPX_STATUS_solved => "solved",
        IPX_STATUS_stopped => "stopped",
        IPX_STATUS_no_model => "no model",
        IPX_STATUS_out_of_memory => "out of memory",
        IPX_STATUS_internal_error => "internal error",
        IPX_STATUS_optimal => "optimal",
        IPX_STATUS_imprecise => "imprecise",
        IPX_STATUS_primal_infeas => "primal infeas",
        IPX_STATUS_dual_infeas => "dual infeas",
        IPX_STATUS_time_limit => "time limit",
        IPX_STATUS_iter_limit => "iter limit",
        IPX_STATUS_no_progress => "no progress",
        IPX_STATUS_failed => "failed",
        IPX_STATUS_debug => "debug",
        _ => "unknown",
    };
    name.to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_status_and_dump_layout() {
        assert_eq!(status_string(IPX_STATUS_optimal), "optimal");
        assert_eq!(status_string(42), "unknown");
        let info = Info {
            status: IPX_STATUS_solved,
            ..Info::default()
        };
        let mut s = String::new();
        write_info(&mut s, &info);
        // Textline: 4 spaces, the name left-justified in 52 columns, then the value.
        let first = s.lines().next().unwrap();
        assert_eq!(
            first,
            format!("    {:<52}{}", "info.status", IPX_STATUS_solved)
        );
    }
}
