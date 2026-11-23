//! Translated from mip-heuristics (MIT License, Copyright (c) 2026 Simon Spoorendonk), commit 36f7d9d. See THIRD_PARTY_NOTICES.md.
//!
//! Ported from SCIP.

#[cfg(test)]
mod tests {
    use super::*;

    use crate::HighsInt;
    use crate::lp_data::hconst::HighsVarType;
    use crate::lp_data::hconst::kHighsInf;
    use crate::lp_data::highs_callback::HighsCallback;
    use crate::lp_data::highs_lp::HighsLp;
    use crate::lp_data::highs_options::HighsOptions;
    use crate::lp_data::hstruct::HighsSolution;
    use crate::lp_data::model_utils::util_model_status_to_string;
    use crate::mip::mip_solver::HighsMipSolver;

    /// Hand-built copy of tests/fixtures/s5/e2e_knap.mps (MPS file order, 20 columns, 3 rows, 60 nonzeros).
    fn s5_e2e_knap_lp() -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = 20;
        lp.num_row_ = 3;
        lp.col_cost_ = vec![
            -10.0, -44.0, -44.0, -24.0, -24.0, -46.0, -36.0, -36.0, -46.0, -56.0, -17.0, -24.0,
            -57.0, -14.0, -50.0, -10.0, -13.0, -23.0, -15.0, -23.0,
        ];
        lp.col_lower_ = vec![
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0,
        ];
        lp.col_upper_ = vec![
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0,
        ];
        lp.row_lower_ = vec![-kHighsInf, -kHighsInf, -kHighsInf];
        lp.row_upper_ = vec![219.0, 279.0, 261.0];
        lp.integrality_ = vec![
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
        ];
        lp.a_matrix_.num_col_ = 20;
        lp.a_matrix_.num_row_ = 3;
        lp.a_matrix_.start_ = vec![
            0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51, 54, 57, 60,
        ];
        lp.a_matrix_.index_ = vec![
            0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1,
            2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0,
            1, 2,
        ];
        lp.a_matrix_.value_ = vec![
            15.0, 39.0, 31.0, 28.0, 28.0, 32.0, 9.0, 29.0, 25.0, 10.0, 30.0, 38.0, 7.0, 27.0, 19.0,
            40.0, 40.0, 36.0, 21.0, 33.0, 9.0, 38.0, 26.0, 34.0, 19.0, 39.0, 27.0, 40.0, 12.0,
            32.0, 17.0, 25.0, 25.0, 18.0, 30.0, 38.0, 39.0, 27.0, 31.0, 12.0, 28.0, 28.0, 5.0,
            21.0, 5.0, 22.0, 38.0, 26.0, 19.0, 7.0, 7.0, 28.0, 40.0, 36.0, 29.0, 5.0, 9.0, 22.0,
            34.0, 34.0,
        ];
        lp
    }

    /// Hand-built copy of tests/fixtures/s5/e2e_setcover.mps (MPS file order, 25 columns, 15 rows, 105 nonzeros).
    fn s5_e2e_setcover_lp() -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = 25;
        lp.num_row_ = 15;
        lp.col_cost_ = vec![
            8.0, 9.0, 12.0, 3.0, 6.0, 9.0, 12.0, 3.0, 8.0, 3.0, 10.0, 11.0, 6.0, 5.0, 6.0, 11.0,
            4.0, 9.0, 12.0, 9.0, 8.0, 9.0, 12.0, 9.0, 10.0,
        ];
        lp.col_lower_ = vec![
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ];
        lp.col_upper_ = vec![
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        ];
        lp.row_lower_ = vec![
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        ];
        lp.row_upper_ = vec![
            kHighsInf, kHighsInf, kHighsInf, kHighsInf, kHighsInf, kHighsInf, kHighsInf, kHighsInf,
            kHighsInf, kHighsInf, kHighsInf, kHighsInf, kHighsInf, kHighsInf, kHighsInf,
        ];
        lp.integrality_ = vec![
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
        ];
        lp.a_matrix_.num_col_ = 25;
        lp.a_matrix_.num_row_ = 15;
        lp.a_matrix_.start_ = vec![
            0, 5, 10, 15, 19, 24, 28, 33, 37, 40, 43, 48, 52, 57, 62, 66, 70, 74, 78, 82, 85, 88,
            93, 97, 101, 105,
        ];
        lp.a_matrix_.index_ = vec![
            7, 8, 10, 11, 14, 0, 4, 5, 7, 11, 0, 1, 6, 7, 8, 3, 6, 8, 10, 1, 2, 5, 10, 11, 4, 7, 9,
            11, 1, 7, 12, 13, 14, 1, 5, 10, 13, 4, 5, 10, 0, 5, 11, 0, 1, 4, 10, 14, 0, 6, 8, 9, 1,
            2, 4, 5, 10, 1, 4, 7, 11, 14, 1, 3, 9, 12, 0, 1, 4, 13, 1, 3, 5, 14, 0, 1, 3, 9, 3, 5,
            9, 14, 1, 8, 11, 2, 3, 6, 3, 7, 10, 11, 14, 0, 3, 5, 9, 0, 3, 8, 13, 4, 5, 6, 7,
        ];
        lp.a_matrix_.value_ = vec![
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0,
        ];
        lp
    }

    /// Hand-built copy of tests/fixtures/s5/e2e_facility.mps (MPS file order, 28 columns, 10 rows, 52 nonzeros).
    fn s5_e2e_facility_lp() -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = 28;
        lp.num_row_ = 10;
        lp.col_cost_ = vec![
            38.0, 30.0, 58.0, 50.0, 48.0, 72.0, 75.0, 30.0, 10.0, 96.0, 12.0, 56.0, 75.0, 75.0,
            70.0, 36.0, 84.0, 40.0, 75.0, 120.0, 80.0, 108.0, 108.0, 56.0, 75.0, 75.0, 70.0, 84.0,
        ];
        lp.col_lower_ = vec![
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ];
        lp.col_upper_ = vec![
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        ];
        lp.row_lower_ = vec![
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -kHighsInf, -kHighsInf, -kHighsInf, -kHighsInf,
        ];
        lp.row_upper_ = vec![1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0];
        lp.integrality_ = vec![
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kContinuous,
            HighsVarType::kContinuous,
            HighsVarType::kContinuous,
            HighsVarType::kContinuous,
            HighsVarType::kContinuous,
            HighsVarType::kContinuous,
            HighsVarType::kContinuous,
            HighsVarType::kContinuous,
            HighsVarType::kContinuous,
            HighsVarType::kContinuous,
            HighsVarType::kContinuous,
            HighsVarType::kContinuous,
            HighsVarType::kContinuous,
            HighsVarType::kContinuous,
            HighsVarType::kContinuous,
            HighsVarType::kContinuous,
            HighsVarType::kContinuous,
            HighsVarType::kContinuous,
            HighsVarType::kContinuous,
            HighsVarType::kContinuous,
            HighsVarType::kContinuous,
            HighsVarType::kContinuous,
            HighsVarType::kContinuous,
            HighsVarType::kContinuous,
        ];
        lp.a_matrix_.num_col_ = 28;
        lp.a_matrix_.num_row_ = 10;
        lp.a_matrix_.start_ = vec![
            0, 1, 2, 3, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40,
            42, 44, 46, 48, 50, 52,
        ];
        lp.a_matrix_.index_ = vec![
            6, 7, 8, 9, 0, 6, 1, 6, 2, 6, 3, 6, 4, 6, 5, 6, 0, 7, 1, 7, 2, 7, 3, 7, 4, 7, 5, 7, 0,
            8, 1, 8, 2, 8, 3, 8, 4, 8, 5, 8, 0, 9, 1, 9, 2, 9, 3, 9, 4, 9, 5, 9,
        ];
        lp.a_matrix_.value_ = vec![
            -26.0, -27.0, -24.0, -33.0, 1.0, 12.0, 1.0, 8.0, 1.0, 15.0, 1.0, 15.0, 1.0, 10.0, 1.0,
            12.0, 1.0, 12.0, 1.0, 8.0, 1.0, 15.0, 1.0, 15.0, 1.0, 10.0, 1.0, 12.0, 1.0, 12.0, 1.0,
            8.0, 1.0, 15.0, 1.0, 15.0, 1.0, 10.0, 1.0, 12.0, 1.0, 12.0, 1.0, 8.0, 1.0, 15.0, 1.0,
            15.0, 1.0, 10.0, 1.0, 12.0,
        ];
        lp
    }

    // ---------------------------------------------------------------------------
    // Harness
    // ---------------------------------------------------------------------------

    /// One solve of the given LP with the given flags.
    /// Returns (objective, num_nodes, total_lp_iterations, report offers, report charges, heuristic_effort_used, modelstatus_string).
    #[allow(clippy::type_complexity)]
    fn s5_e2e_run(
        lp: &HighsLp,
        fpr: bool,
        lm: bool,
        seed: Option<HighsInt>,
    ) -> (
        f64,
        i64,
        i64,
        Vec<(&'static str, usize, f64)>,
        Vec<(&'static str, usize, bool)>,
        usize,
        String,
    ) {
        let mut options = HighsOptions::default();
        options.output_flag = false;
        options.log_options.output_flag = false;
        options.threads = 1;
        options.mip_heuristic_run_feasibility_jump = false;
        options.mip_heuristic_run_fpr = fpr;
        options.mip_heuristic_run_local_mip = lm;
        if let Some(s) = seed {
            options.random_seed = s;
        }
        let solution = HighsSolution::default();
        let mut callback = HighsCallback::default();
        let mut ms = HighsMipSolver::new(&mut callback, &options, lp, &solution, false, 0);
        ms.run();
        let obj = ms.solution_objective_;
        let nodes = ms.mipdata_.num_nodes;
        let lp_iters = ms.mipdata_.total_lp_iterations;
        let status_str = util_model_status_to_string(ms.modelstatus_);
        let (offers, charges) = if let Some(ref report) = ms.mipdata_.s5_report {
            let off: Vec<(&'static str, usize, f64)> = report
                .offers
                .iter()
                .map(|o| (o.name, o.effort_at, o.obj))
                .collect();
            let ch: Vec<(&'static str, usize, bool)> = report
                .heur
                .iter()
                .map(|h| (h.name, h.effort, h.found))
                .collect();
            (off, ch)
        } else {
            (Vec::new(), Vec::new())
        };
        let effort_used = ms.mipdata_.heuristic_effort_used;
        (
            obj,
            nodes,
            lp_iters,
            offers,
            charges,
            effort_used,
            status_str,
        )
    }

    // ---------------------------------------------------------------------------
    // Fixture parsing
    // ---------------------------------------------------------------------------

    struct ParsedHeurSol {
        name: String,
        effort_at: usize,
        obj: f64,
    }

    struct ParsedSequential {
        name: String,
        effort: usize,
        /// `[Heur] found` (effort_ledger.cpp:28-50): the chain entry had at
        /// least one accepted offer, i.e. a `[HeurSol]` line since the
        /// previous `[Sequential]` line.
        found: bool,
    }

    struct ParsedFinal {
        obj: f64,
        status: String,
    }

    struct ParsedFixture {
        heur_sols: Vec<ParsedHeurSol>,
        sequential: Vec<ParsedSequential>,
        final_line: ParsedFinal,
    }

    fn parse_fixture(text: &str) -> ParsedFixture {
        let mut heur_sols = Vec::new();
        let mut sequential = Vec::new();
        let mut final_line = ParsedFinal {
            obj: 0.0,
            status: String::new(),
        };
        let mut accepted_since_charge = false;
        for line in text.lines() {
            let line = line.trim();
            if line.is_empty() {
                continue;
            }
            if line.starts_with("[HeurSol]") {
                assert_eq!(extract_field(line, "accepted="), "1");
                accepted_since_charge = true;
                // [HeurSol] name=N ... effort_at=E ... obj=V accepted=1
                let name = extract_field(line, "name=").to_string();
                let effort_at: usize = extract_field(line, "effort_at=").parse().unwrap();
                let obj: f64 = extract_field(line, "obj=").parse().unwrap();
                heur_sols.push(ParsedHeurSol {
                    name,
                    effort_at,
                    obj,
                });
            } else if line.starts_with("[Sequential]") {
                // [Sequential] heur=N effort=E
                let name = extract_field(line, "heur=").to_string();
                let effort: usize = extract_field(line, "effort=").parse().unwrap();
                sequential.push(ParsedSequential {
                    name,
                    effort,
                    found: accepted_since_charge,
                });
                accepted_since_charge = false;
            } else if line.starts_with("final ") {
                // final obj=V status=S
                let obj: f64 = extract_field(line, "obj=").parse().unwrap();
                let status = extract_field(line, "status=").to_string();
                final_line = ParsedFinal { obj, status };
            }
        }
        ParsedFixture {
            heur_sols,
            sequential,
            final_line,
        }
    }

    fn extract_field<'a>(line: &'a str, key: &str) -> &'a str {
        let idx = line.find(key).unwrap();
        let rest = &line[idx + key.len()..];
        let end = rest.find(char::is_whitespace).unwrap_or(rest.len());
        &rest[..end]
    }

    // ---------------------------------------------------------------------------
    // Selfcheck parsing
    // ---------------------------------------------------------------------------

    struct SelfcheckBlock {
        model_name: String,
        rusolve_status: String,
        rusolve_obj: f64,
        rusolve_nodes: i64,
        rusolve_lp_iters: i64,
    }

    fn parse_selfcheck(text: &str) -> Vec<SelfcheckBlock> {
        let mut blocks = Vec::new();
        let mut current_model: Option<String> = None;
        for line in text.lines() {
            if line.trim().is_empty() {
                continue;
            }
            // Header line: "e2e_knap MATCH" or similar
            if !line.starts_with("  ") && !line.starts_with('\t') {
                let line = line.trim();
                // This is a header line; extract model name (first token)
                let first_token = line.split_whitespace().next().unwrap_or("");
                if first_token.starts_with("e2e_") {
                    current_model = Some(first_token.to_string());
                }
            } else if line.starts_with("  rusolve ") {
                // "  rusolve status=Optimal obj=<V> nodes=<N> lp_iters=<L>"
                if let Some(ref model) = current_model {
                    let status = extract_field(line, "status=").to_string();
                    let obj: f64 = extract_field(line, "obj=").parse().unwrap();
                    let nodes: i64 = extract_field(line, "nodes=").parse().unwrap();
                    let lp_iters: i64 = extract_field(line, "lp_iters=").parse().unwrap();
                    blocks.push(SelfcheckBlock {
                        model_name: model.clone(),
                        rusolve_status: status,
                        rusolve_obj: obj,
                        rusolve_nodes: nodes,
                        rusolve_lp_iters: lp_iters,
                    });
                    current_model = None;
                }
            }
        }
        blocks
    }

    fn find_selfcheck_block<'a>(
        blocks: &'a [SelfcheckBlock],
        model_prefix: &str,
    ) -> &'a SelfcheckBlock {
        blocks
            .iter()
            .find(|b| b.model_name.starts_with(model_prefix))
            .unwrap_or_else(|| panic!("selfcheck block not found for {}", model_prefix))
    }

    // ---------------------------------------------------------------------------
    // Fixture text constants
    // ---------------------------------------------------------------------------

    const FIXTURE_KNAP_FPR: &str = include_str!("../../../tests/fixtures/s5/e2e_knap.fpr.heur");
    const FIXTURE_KNAP_LM: &str = include_str!("../../../tests/fixtures/s5/e2e_knap.lm.heur");
    const FIXTURE_KNAP_BOTH: &str = include_str!("../../../tests/fixtures/s5/e2e_knap.both.heur");
    const FIXTURE_SETCOVER_FPR: &str =
        include_str!("../../../tests/fixtures/s5/e2e_setcover.fpr.heur");
    const FIXTURE_SETCOVER_LM: &str =
        include_str!("../../../tests/fixtures/s5/e2e_setcover.lm.heur");
    const FIXTURE_SETCOVER_BOTH: &str =
        include_str!("../../../tests/fixtures/s5/e2e_setcover.both.heur");
    const FIXTURE_FACILITY_FPR: &str =
        include_str!("../../../tests/fixtures/s5/e2e_facility.fpr.heur");
    const FIXTURE_FACILITY_LM: &str =
        include_str!("../../../tests/fixtures/s5/e2e_facility.lm.heur");
    const FIXTURE_FACILITY_BOTH: &str =
        include_str!("../../../tests/fixtures/s5/e2e_facility.both.heur");
    const SELFCHECK: &str = include_str!("../../../tests/fixtures/s5/selfcheck.txt");

    // ---------------------------------------------------------------------------
    // Tests
    // ---------------------------------------------------------------------------

    #[test]
    #[ignore = "s5_e2e slow"]
    fn s5_e2e_trace_matches_fixture() {
        let models: Vec<(&str, fn() -> HighsLp)> = vec![
            ("knap", s5_e2e_knap_lp),
            ("setcover", s5_e2e_setcover_lp),
            ("facility", s5_e2e_facility_lp),
        ];
        let configs: Vec<(&str, bool, bool, &str)> = vec![
            ("fpr", true, false, ""),
            ("lm", false, true, ""),
            ("both", true, true, ""),
        ];

        for (model_name, lp_fn) in &models {
            let lp = lp_fn();
            for (cfg_name, fpr, lm, _) in &configs {
                let fixture_text = match (*model_name, *cfg_name) {
                    ("knap", "fpr") => FIXTURE_KNAP_FPR,
                    ("knap", "lm") => FIXTURE_KNAP_LM,
                    ("knap", "both") => FIXTURE_KNAP_BOTH,
                    ("setcover", "fpr") => FIXTURE_SETCOVER_FPR,
                    ("setcover", "lm") => FIXTURE_SETCOVER_LM,
                    ("setcover", "both") => FIXTURE_SETCOVER_BOTH,
                    ("facility", "fpr") => FIXTURE_FACILITY_FPR,
                    ("facility", "lm") => FIXTURE_FACILITY_LM,
                    ("facility", "both") => FIXTURE_FACILITY_BOTH,
                    _ => unreachable!(),
                };
                let fixture = parse_fixture(fixture_text);
                let (_obj, _nodes, _lp_iters, offers, charges, _effort_used, status_str) =
                    s5_e2e_run(&lp, *fpr, *lm, None);

                // Compare offers: name, effort_at, obj.to_bits()
                assert_eq!(
                    offers.len(),
                    fixture.heur_sols.len(),
                    "offer count mismatch for {} {}",
                    model_name,
                    cfg_name
                );
                for (i, (run_offer, fix_offer)) in
                    offers.iter().zip(fixture.heur_sols.iter()).enumerate()
                {
                    assert_eq!(
                        run_offer.0, fix_offer.name,
                        "offer[{}] name mismatch for {} {}",
                        i, model_name, cfg_name
                    );
                    assert_eq!(
                        run_offer.1, fix_offer.effort_at,
                        "offer[{}] effort_at mismatch for {} {}",
                        i, model_name, cfg_name
                    );
                    assert_eq!(
                        run_offer.2.to_bits(),
                        fix_offer.obj.to_bits(),
                        "offer[{}] obj bits mismatch for {} {}",
                        i,
                        model_name,
                        cfg_name
                    );
                }

                // Compare charges: name, effort, in order
                assert_eq!(
                    charges.len(),
                    fixture.sequential.len(),
                    "charge count mismatch for {} {}",
                    model_name,
                    cfg_name
                );
                for (i, (run_charge, fix_charge)) in
                    charges.iter().zip(fixture.sequential.iter()).enumerate()
                {
                    assert_eq!(
                        run_charge.0, fix_charge.name,
                        "charge[{}] name mismatch for {} {}",
                        i, model_name, cfg_name
                    );
                    assert_eq!(
                        run_charge.1, fix_charge.effort,
                        "charge[{}] effort mismatch for {} {}",
                        i, model_name, cfg_name
                    );
                    assert_eq!(
                        run_charge.2, fix_charge.found,
                        "charge[{}] found mismatch for {} {}",
                        i, model_name, cfg_name
                    );
                }

                // Compare final obj and status
                assert_eq!(
                    _obj.to_bits(),
                    fixture.final_line.obj.to_bits(),
                    "final obj bits mismatch for {} {}",
                    model_name,
                    cfg_name
                );
                assert_eq!(
                    status_str, fixture.final_line.status,
                    "final status mismatch for {} {}",
                    model_name, cfg_name
                );
            }
        }
    }

    #[test]
    #[ignore = "s5_e2e slow"]
    fn s5_e2e_two_runs_identical() {
        let models: Vec<(&str, fn() -> HighsLp)> = vec![
            ("knap", s5_e2e_knap_lp),
            ("setcover", s5_e2e_setcover_lp),
            ("facility", s5_e2e_facility_lp),
        ];
        let configs: Vec<(&str, bool, bool)> = vec![
            ("fpr", true, false),
            ("lm", false, true),
            ("both", true, true),
        ];

        for (model_name, lp_fn) in &models {
            let lp = lp_fn();
            for (cfg_name, fpr, lm) in &configs {
                let r1 = s5_e2e_run(&lp, *fpr, *lm, None);
                let r2 = s5_e2e_run(&lp, *fpr, *lm, None);

                assert_eq!(
                    r1.0.to_bits(),
                    r2.0.to_bits(),
                    "obj bits differ for {} {}",
                    model_name,
                    cfg_name
                );
                assert_eq!(
                    r1.3.len(),
                    r2.3.len(),
                    "offer count differs for {} {}",
                    model_name,
                    cfg_name
                );
                for (i, (o1, o2)) in r1.3.iter().zip(r2.3.iter()).enumerate() {
                    assert_eq!(
                        o1.0, o2.0,
                        "offer[{}] name differs for {} {}",
                        i, model_name, cfg_name
                    );
                    assert_eq!(
                        o1.1, o2.1,
                        "offer[{}] effort_at differs for {} {}",
                        i, model_name, cfg_name
                    );
                    assert_eq!(
                        o1.2.to_bits(),
                        o2.2.to_bits(),
                        "offer[{}] obj bits differ for {} {}",
                        i,
                        model_name,
                        cfg_name
                    );
                }
                assert_eq!(
                    r1.4.len(),
                    r2.4.len(),
                    "charge count differs for {} {}",
                    model_name,
                    cfg_name
                );
                for (i, (c1, c2)) in r1.4.iter().zip(r2.4.iter()).enumerate() {
                    assert_eq!(
                        c1.0, c2.0,
                        "charge[{}] name differs for {} {}",
                        i, model_name, cfg_name
                    );
                    assert_eq!(
                        c1.1, c2.1,
                        "charge[{}] effort differs for {} {}",
                        i, model_name, cfg_name
                    );
                    assert_eq!(
                        c1.2, c2.2,
                        "charge[{}] found differs for {} {}",
                        i, model_name, cfg_name
                    );
                }
            }
        }
    }

    #[test]
    #[ignore = "s5_e2e slow"]
    fn s5_e2e_random_seed_1_optimum() {
        let models: Vec<(&str, fn() -> HighsLp)> = vec![
            ("knap", s5_e2e_knap_lp),
            ("setcover", s5_e2e_setcover_lp),
            ("facility", s5_e2e_facility_lp),
        ];

        for (model_name, lp_fn) in &models {
            let lp = lp_fn();
            // Off-path objective (both flags false)
            let (off_obj, _, _, _, _, _, _) = s5_e2e_run(&lp, false, false, None);
            // Both flags on, random_seed=1
            let (on_obj, _, _, _, _, _, _) = s5_e2e_run(&lp, true, true, Some(1));

            assert!(
                (on_obj - off_obj).abs() < 1e-6,
                "seed=1 obj {} differs from off-path obj {} for {}",
                on_obj,
                off_obj,
                model_name
            );
        }
    }

    #[test]
    #[ignore = "s5_e2e slow"]
    fn s5_e2e_off_matches_selfcheck() {
        let blocks = parse_selfcheck(SELFCHECK);
        let models: Vec<(&str, fn() -> HighsLp)> = vec![
            ("knap", s5_e2e_knap_lp),
            ("setcover", s5_e2e_setcover_lp),
            ("facility", s5_e2e_facility_lp),
        ];

        for (model_name, lp_fn) in &models {
            let lp = lp_fn();
            let (obj, nodes, lp_iters, _, _, effort_used, status_str) =
                s5_e2e_run(&lp, false, false, None);
            let block = find_selfcheck_block(&blocks, &format!("e2e_{}", model_name));

            assert_eq!(
                obj.to_bits(),
                block.rusolve_obj.to_bits(),
                "off-path obj bits mismatch for {}",
                model_name
            );
            assert_eq!(
                nodes, block.rusolve_nodes,
                "off-path nodes mismatch for {}",
                model_name
            );
            assert_eq!(
                lp_iters, block.rusolve_lp_iters,
                "off-path lp_iters mismatch for {}",
                model_name
            );
            assert_eq!(
                status_str, block.rusolve_status,
                "off-path status mismatch for {}",
                model_name
            );
            assert_eq!(
                effort_used, 0,
                "heuristic_effort_used should be 0 when off for {}",
                model_name
            );
        }
    }

    #[test]
    #[ignore = "s5_e2e slow"]
    fn s5_e2e_off_has_no_report() {
        let models: Vec<(&str, fn() -> HighsLp)> = vec![
            ("knap", s5_e2e_knap_lp),
            ("setcover", s5_e2e_setcover_lp),
            ("facility", s5_e2e_facility_lp),
        ];

        for (model_name, lp_fn) in &models {
            let lp = lp_fn();
            let mut options = HighsOptions::default();
            options.output_flag = false;
            options.log_options.output_flag = false;
            options.threads = 1;
            options.mip_heuristic_run_feasibility_jump = false;
            options.mip_heuristic_run_fpr = false;
            options.mip_heuristic_run_local_mip = false;
            let solution = HighsSolution::default();
            let mut callback = HighsCallback::default();
            let mut ms = HighsMipSolver::new(&mut callback, &options, &lp, &solution, false, 0);
            ms.run();

            assert!(
                ms.mipdata_.s5_report.is_none(),
                "s5_report should be None when off for {}",
                model_name
            );
            assert_eq!(
                ms.mipdata_.heuristic_effort_used, 0,
                "heuristic_effort_used should be 0 when off for {}",
                model_name
            );
        }
    }
}
