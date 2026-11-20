//! HiGHS source: mip/feasibilityjump.hh
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port.

use crate::prelude::*;

/// HiGHS: feasibilityjump.hh::FJStatus (lines 38-44)
#[derive(Clone, Debug)]
pub struct FJStatus<'a> {
    /// C++: `size_t totalEffort`
    pub total_effort: usize,
    /// C++: `size_t effortSinceLastImprovement`
    pub effort_since_last_improvement: usize,
    /// C++: `int numVars`
    pub num_vars: i32,
    /// C++: `double solutionObjectiveValue`
    pub solution_objective_value: f64,
    /// C++ `double* solution` (read-only view for the callback).
    /// PORT NOTE: C++ nullptr (periodic callback, no solution yet) is `None`.
    pub solution: Option<&'a [f64]>,
}

/// HiGHS: feasibilityjump.hh::IdxCoeff (lines 49-54)
#[derive(Clone, Copy, Default, Debug, PartialEq)]
pub struct IdxCoeff {
    /// C++: `uint32_t idx`
    pub idx: u32,
    /// C++: `double coeff`
    pub coeff: f64,
}

/// HiGHS: feasibilityjump.hh::Var (lines 56-62)
#[derive(Clone, Default, Debug)]
pub struct Var {
    /// C++: `VarType vartype`
    pub vartype: VarType,
    /// C++: `double lb`
    pub lb: f64,
    /// C++: `double ub`
    pub ub: f64,
    /// C++: `double objectiveCoeff`
    pub objective_coeff: f64,
    /// C++: `std::vector<IdxCoeff> coeffs`
    pub coeffs: Vec<IdxCoeff>,
}

/// HiGHS: feasibilityjump.hh::Constraint (lines 64-83)
#[derive(Clone, Default, Debug)]
pub struct Constraint {
    /// C++: `RowType sense`
    pub sense: RowType,
    /// C++: `double rhs`
    pub rhs: f64,
    /// C++: `std::vector<IdxCoeff> coeffs`
    pub coeffs: Vec<IdxCoeff>,
    /// C++: `double weight`
    pub weight: f64,
    /// C++: `double incumbentLhs`
    pub incumbent_lhs: f64,
    /// C++: `int32_t violatedIdx`
    pub violated_idx: i32,
}

/// HiGHS: feasibilityjump.hh::Move (lines 86-96)
#[derive(Clone, Copy, Default, Debug, PartialEq)]
pub struct Move {
    /// C++: `double value`
    pub value: f64,
    /// C++: `double score`
    pub score: f64,
}

/// HiGHS: feasibilityjump.hh::LhsModification (lines 101-107)
#[derive(Clone, Copy, Default, Debug, PartialEq)]
pub struct LhsModification {
    /// C++: `uint32_t varIdx`
    pub var_idx: u32,
    /// C++: `uint32_t constraintIdx`
    pub constraint_idx: u32,
    /// C++: `double coeff`
    pub coeff: f64,
    /// C++: `double oldLhs`
    pub old_lhs: f64,
    /// C++: `double newLhs`
    pub new_lhs: f64,
}

/// HiGHS: feasibilityjump.hh::Problem (lines 112-297)
#[derive(Clone, Debug)]
pub struct Problem {
    /// C++: `const double equalityTolerance`
    pub equality_tolerance: f64,
    /// C++: `const double violationTolerance`
    pub violation_tolerance: f64,
    /// C++: `std::vector<Var> vars`
    pub vars: Vec<Var>,
    /// C++: `std::vector<Constraint> constraints`
    pub constraints: Vec<Constraint>,
    /// C++: `std::vector<double> incumbentAssignment`
    pub incumbent_assignment: Vec<f64>,
    /// C++: `std::vector<uint32_t> violatedConstraints`
    pub violated_constraints: Vec<u32>,
    /// C++: `bool usedRelaxContinuous = false`
    pub used_relax_continuous: bool,
    /// C++: `size_t nNonzeros = 0`
    pub n_nonzeros: usize,
    /// C++: `double incumbentObjective = NAN`
    pub incumbent_objective: f64,
}

impl Default for Problem {
    fn default() -> Self {
        Self {
            equality_tolerance: Default::default(),
            violation_tolerance: Default::default(),
            vars: Default::default(),
            constraints: Default::default(),
            incumbent_assignment: Default::default(),
            violated_constraints: Default::default(),
            used_relax_continuous: false,
            n_nonzeros: 0,
            incumbent_objective: f64::NAN,
        }
    }
}

/// HiGHS: feasibilityjump.hh::JumpMove (lines 313-434)
#[derive(Clone, Default, Debug)]
pub struct JumpMove {
    /// C++: `const double equalityTolerance`
    pub equality_tolerance: f64,
    /// C++: `std::vector<Move> moves`
    pub moves: Vec<Move>,
    /// C++: `std::vector<std::pair<double, double>> bestShiftBuffer`
    pub best_shift_buffer: Vec<(f64, f64)>,
    /// HiGHS 1.16 #3310 (h16_fj): option flag (no C++ counterpart). Set from `options_mip_.h16_fj` in
    /// feasibility_jump.rs; `false` runs the HiGHS 1.15.1 statements.
    pub h16_fj: bool,
}

/// HiGHS: feasibilityjump.hh::FeasibilityJumpSolver (lines 436-798)
#[derive(Clone, Debug)]
pub struct FeasibilityJumpSolver {
    /// C++: `Problem problem`
    pub problem: Problem,
    /// C++: `JumpMove jumpMove`
    pub jump_move: JumpMove,
    /// C++: `std::vector<uint32_t> goodVarsSet`
    pub good_vars_set: Vec<u32>,
    /// C++: `std::vector<int32_t> goodVarsSetIdx`
    pub good_vars_set_idx: Vec<i32>,
    /// C++: `std::mt19937 rng`
    pub rng: Mt19937,
    /// C++: `double bestObjective = std :: numeric_limits < double > :: infinity ( )`
    pub best_objective: f64,
    /// C++: `double objectiveWeight = 0.0`
    pub objective_weight: f64,
    /// C++: `size_t bestViolationScore = SIZE_MAX`
    pub best_violation_score: usize,
    /// C++: `size_t effortAtLastCallback = 0`
    pub effort_at_last_callback: usize,
    /// C++: `size_t effortAtLastImprovement = 0`
    pub effort_at_last_improvement: usize,
    /// C++: `size_t effortAtLastLogging = 0`
    pub effort_at_last_logging: usize,
    /// C++: `size_t totalEffort = 0`
    pub total_effort: usize,
    /// C++: `double weightUpdateDecay`
    pub weight_update_decay: f64,
    /// C++: `double weightUpdateIncrement = 1.0`
    pub weight_update_increment: f64,
    /// C++: `size_t nBumps = 0`
    pub n_bumps: usize,
    /// C++: `bool need_logging_header = true`
    pub need_logging_header: bool,
    /// C++: `const double randomVarProbability = 0.001`
    pub random_var_probability: f64,
    /// C++: `const double randomCellProbability = 0.01`
    pub random_cell_probability: f64,
    /// C++: `const size_t maxMovesToEvaluate = 25`
    pub max_moves_to_evaluate: usize,
    /// C++: `const int kLoggingFrequency = 100000`
    pub k_logging_frequency: i32,
    /// C++: `const int kLoggingHeaderFrequency = 50`
    pub k_logging_header_frequency: i32,
    /// C++: `const size_t kMinEffortToLogging = 500 * kLoggingFrequency`
    pub k_min_effort_to_logging: usize,
    /// C++: `const double kMinRelativeObjectiveImprovement = 1e-4`
    pub k_min_relative_objective_improvement: f64,
    /// C++ `const HighsLogOptions& logOptions`: owned clone.
    pub log_options: HighsLogOptions,
    /// HiGHS 1.16 #3310 / #3269 (h16_fj): option flag (no C++ counterpart). Set from
    /// `options_mip_.h16_fj` in feasibility_jump.rs; `false` runs the HiGHS 1.15.1 statements.
    pub h16_fj: bool,
    /// rusolve throughput counter (fpr_local_mip/stats.rs):
    /// candidate variables whose best move was scored in `select_variable`. Never read by the search.
    pub ls_evals: u64,
    /// rusolve throughput counter: variable moves applied (`solve` steps). Never read by the search.
    pub ls_steps: u64,
}

impl Default for FeasibilityJumpSolver {
    fn default() -> Self {
        Self {
            problem: Default::default(),
            jump_move: Default::default(),
            good_vars_set: Default::default(),
            good_vars_set_idx: Default::default(),
            rng: Default::default(),
            best_objective: f64::INFINITY,
            objective_weight: 0.0,
            best_violation_score: usize::MAX,
            effort_at_last_callback: 0,
            effort_at_last_improvement: 0,
            effort_at_last_logging: 0,
            total_effort: 0,
            weight_update_decay: Default::default(),
            weight_update_increment: 1.0,
            n_bumps: 0,
            need_logging_header: true,
            random_var_probability: 0.001,
            random_cell_probability: 0.01,
            max_moves_to_evaluate: 25,
            k_logging_frequency: 100000,
            k_logging_header_frequency: 50,
            k_min_effort_to_logging: 500 * 100000, /* kLoggingFrequency */
            k_min_relative_objective_improvement: 1e-4,
            log_options: Default::default(),
            h16_fj: false,
            ls_evals: 0,
            ls_steps: 0,
        }
    }
}

/// HiGHS: feasibilityjump.hh::RowType (plain C enum: module consts).
pub type RowType = i32;
pub const Equal: i32 = 0;
pub const Lte: i32 = 1;
pub const Gte: i32 = 2;

/// HiGHS: feasibilityjump.hh::VarType (plain C enum: module consts).
pub type VarType = i32;
pub const Continuous: i32 = 0;
pub const Integer: i32 = 1;

/// HiGHS: feasibilityjump.hh::CallbackControlFlow (plain C enum: module consts).
pub type CallbackControlFlow = i32;
pub const Terminate: i32 = 0;
pub const Continue: i32 = 1;

/// HiGHS: feasibilityjump.hh::eq (lines 47-47)
/// C++: `static bool eq(double a, double b, double tol)`
pub fn eq(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() < tol
}

/// HiGHS: feasibilityjump.hh::modifyMove (lines 299-309)
/// C++: `static void modifyMove(LhsModification mod, Problem & problem, Move & move)`
pub fn modify_move(r#mod: LhsModification, problem: &mut Problem, r#move: &mut Move) {
    modify_move_in(
        r#mod,
        &problem.constraints[r#mod.constraint_idx as usize],
        &problem.incumbent_assignment,
        r#move,
    );
}

/// PORT NOTE (perf): the body of [`modify_move`] on the two `Problem` parts it reads
/// (`constraints[mod.constraint_idx]` and `incumbent_assignment`), so that
/// [`FeasibilityJumpSolver::do_variable_move`] can run it from inside `set_value` as the C++
/// lambda does. Same operations in the same order.
#[inline(always)]
fn modify_move_in(
    r#mod: LhsModification,
    cn: &Constraint,
    incumbent_assignment: &[f64],
    r#move: &mut Move,
) {
    let incumbent = incumbent_assignment[r#mod.var_idx as usize];
    let old_modified_lhs = r#mod.old_lhs + r#mod.coeff * (r#move.value - incumbent);
    let old_score_term = cn.weight * (cn.score(old_modified_lhs) - cn.score(r#mod.old_lhs));
    let new_modified_lhs = r#mod.new_lhs + r#mod.coeff * (r#move.value - incumbent);
    let new_score_term = cn.weight * (cn.score(new_modified_lhs) - cn.score(r#mod.new_lhs));
    r#move.score += new_score_term - old_score_term;
}

/// PORT NOTE (perf): the body of `FeasibilityJumpSolver::updateGoodMoves` (feasibilityjump.hh
/// lines 723-739) on the solver fields it touches, so it can run inside the `set_value`
/// callback. `bestMove(varIdx)` is inlined: `forEachMove` visits the single JumpMove move
/// `moves[varIdx]`, kept if its score beats `Move::undef().score`.
#[inline(always)]
fn update_good_moves_in(
    moves: &[Move],
    good_vars_set: &mut Vec<u32>,
    good_vars_set_idx: &mut [i32],
    var_idx: i32,
) {
    let mut best = Move::undef();
    let m = &moves[var_idx as usize];
    if m.score > best.score {
        best = *m;
    }
    let any_good_moves = best.score > 0.0;
    if any_good_moves && good_vars_set_idx[var_idx as usize] == -1 {
        // Became good, add to good set.
        good_vars_set_idx[var_idx as usize] = good_vars_set.len() as i32;
        good_vars_set.push(var_idx as u32);
    } else if !any_good_moves && good_vars_set_idx[var_idx as usize] != -1 {
        // Became bad, remove from good set.
        let last_set_idx = good_vars_set.len() - 1;
        let last_var_idx = good_vars_set[last_set_idx];
        let this_set_idx = good_vars_set_idx[var_idx as usize];
        good_vars_set.swap(this_set_idx as usize, last_set_idx);
        good_vars_set_idx[last_var_idx as usize] = this_set_idx;
        good_vars_set_idx[var_idx as usize] = -1;
        good_vars_set.pop();
    }
}

impl IdxCoeff {
    /// HiGHS: feasibilityjump.hh::IdxCoeff::IdxCoeff (lines 53-53)
    /// C++: `IdxCoeff(uint32_t idx, double coeff)`
    pub fn new(idx: u32, coeff: f64) -> Self {
        Self { idx, coeff }
    }
}

impl Constraint {
    /// HiGHS: feasibilityjump.hh::Constraint::score (lines 75-82)
    /// C++: `double score(double lhs)`
    pub fn score(&self, lhs: f64) -> f64 {
        if self.sense == Equal {
            -(lhs - self.rhs).abs()
        } else if self.sense == Lte {
            -cpp_max(0.0, lhs - self.rhs)
        } else {
            -cpp_max(0.0, self.rhs - lhs)
        }
    }
}

impl Move {
    /// HiGHS: feasibilityjump.hh::Move::undef (lines 90-95)
    /// C++: `static Move undef()`
    pub fn undef() -> Move {
        Move {
            value: f64::NAN,
            score: -f64::INFINITY,
        }
    }
}

impl Problem {
    /// HiGHS: feasibilityjump.hh::Problem::Problem (lines 117-119)
    /// C++: `Problem(double equalityTolerance, double violationTolerance)`
    pub fn new(equality_tolerance: f64, violation_tolerance: f64) -> Self {
        Self {
            equality_tolerance,
            violation_tolerance,
            ..Default::default()
        }
    }

    /// HiGHS: feasibilityjump.hh::Problem::addVar (lines 129-139)
    /// C++: `int addVar(VarType vartype, double lb, double ub, double objCoeff)`
    pub fn add_var(&mut self, vartype: VarType, lb: f64, ub: f64, obj_coeff: f64) -> i32 {
        let idx = self.vars.len();
        self.vars.push(Var {
            vartype,
            lb,
            ub,
            objective_coeff: obj_coeff,
            coeffs: Vec::new(),
        });
        self.incumbent_assignment.push(lb);
        idx as i32
    }

    /// HiGHS: feasibilityjump.hh::Problem::addConstraint (lines 141-206)
    /// C++: `int addConstraint(RowType sense, double rhs, int numCoeffs, int * rowVarIdxs, double * rowCoeffs, int relax_continuous)`
    pub fn add_constraint(
        &mut self,
        sense: RowType,
        rhs: f64,
        num_coeffs: i32,
        row_var_idxs: &mut [i32],
        row_coeffs: &mut [f64],
        relax_continuous: i32,
    ) -> i32 {
        let mut rhs = rhs;
        if relax_continuous != 0 {
            self.used_relax_continuous = true;
        }

        // If we are relaxing continuous variables, an equality needs to be split
        // into Gte and Lte.
        if relax_continuous > 0 && sense == Equal {
            let any_continuous = (0..num_coeffs as usize)
                .any(|i| self.vars[row_var_idxs[i] as usize].vartype == Continuous);
            if any_continuous {
                self.add_constraint(
                    Gte,
                    rhs,
                    num_coeffs,
                    row_var_idxs,
                    row_coeffs,
                    relax_continuous,
                );
                self.add_constraint(
                    Lte,
                    rhs,
                    num_coeffs,
                    row_var_idxs,
                    row_coeffs,
                    relax_continuous,
                );
                return i32::MAX;
            }
        }

        let mut coeffs: Vec<IdxCoeff> = Vec::new();
        for i in 0..num_coeffs as usize {
            if relax_continuous > 0 && self.vars[row_var_idxs[i] as usize].vartype == Continuous {
                if sense == Lte {
                    if row_coeffs[i] >= 0.0 {
                        rhs -= row_coeffs[i] * self.vars[row_var_idxs[i] as usize].lb;
                    } else {
                        rhs -= row_coeffs[i] * self.vars[row_var_idxs[i] as usize].ub;
                    }
                } else if sense == Gte {
                    if row_coeffs[i] >= 0.0 {
                        rhs -= row_coeffs[i] * self.vars[row_var_idxs[i] as usize].ub;
                    } else {
                        rhs -= row_coeffs[i] * self.vars[row_var_idxs[i] as usize].lb;
                    }
                } else {
                    return i32::MIN;
                }
            } else {
                coeffs.push(IdxCoeff::new(row_var_idxs[i] as u32, row_coeffs[i]));
            }
        }

        if coeffs.is_empty() {
            let ok = if sense == Lte {
                0.0 <= rhs + self.equality_tolerance
            } else if sense == Gte {
                0.0 + self.equality_tolerance >= rhs
            } else {
                eq(0.0, rhs, self.equality_tolerance)
            };
            return if ok { i32::MAX } else { i32::MIN };
        }

        let new_constraint_idx = self.constraints.len() as i32;
        for cc in coeffs.iter() {
            self.vars[cc.idx as usize]
                .coeffs
                .push(IdxCoeff::new(new_constraint_idx as u32, cc.coeff));
        }

        self.n_nonzeros += coeffs.len();
        let new_constraint = Constraint {
            sense,
            rhs,
            coeffs,
            weight: 1.0,
            incumbent_lhs: f64::NAN,
            violated_idx: -1,
        };
        self.constraints.push(new_constraint);
        new_constraint_idx
    }

    /// HiGHS: feasibilityjump.hh::Problem::resetIncumbent (lines 208-237)
    /// C++: `void resetIncumbent(double * initialValues)`
    pub fn reset_incumbent(&mut self, initial_values: &mut [f64]) {
        // PORT NOTE: C++ `if (initialValues)`: a null pointer is an empty slice here.
        if !initial_values.is_empty() {
            // Kept as the C++ index loop (its std::copy is commented out there).
            #[allow(clippy::manual_memcpy)]
            for i in 0..self.vars.len() {
                self.incumbent_assignment[i] = initial_values[i];
            }
        }

        // Reset the incumbent objective.
        self.incumbent_objective = 0.0;
        for i in 0..self.vars.len() {
            self.incumbent_objective += self.vars[i].objective_coeff * self.incumbent_assignment[i];
        }

        // Reset the constraint LHSs and the violatedConstraints list.
        self.violated_constraints.clear();
        for c_idx in 0..self.constraints.len() {
            let cstr = &mut self.constraints[c_idx];

            cstr.incumbent_lhs = 0.0;
            for vc in cstr.coeffs.iter() {
                cstr.incumbent_lhs += vc.coeff * self.incumbent_assignment[vc.idx as usize];
            }

            let lhs = cstr.incumbent_lhs;
            if cstr.score(lhs) < -self.violation_tolerance {
                cstr.violated_idx = self.violated_constraints.len() as i32;
                self.violated_constraints.push(c_idx as u32);
            } else {
                cstr.violated_idx = -1;
            }
        }
    }

    /// HiGHS: feasibilityjump.hh::Problem::setValue (lines 243-296)
    /// C++: `size_t setValue()`
    pub fn set_value(
        &mut self,
        var_idx: u32,
        new_value: f64,
        f: &mut dyn FnMut(LhsModification),
    ) -> usize {
        self.set_value_with(var_idx, new_value, |_, _, m| f(m))
    }

    /// PORT NOTE (perf): [`Problem::set_value`] as the C++ template (`F f` is inlined, not a
    /// `dyn` call). The callback also receives `constraints[m.constraint_idx]` (already holding
    /// the new `incumbent_lhs`) and `incumbent_assignment` (already holding the new value of
    /// `var_idx`), which is what C++ `f` reads through `problem`. Same operations in the same
    /// order as `set_value`.
    #[inline(always)]
    pub fn set_value_with<F: FnMut(&Constraint, &[f64], LhsModification)>(
        &mut self,
        var_idx: u32,
        new_value: f64,
        mut f: F,
    ) -> usize {
        let mut dt: usize = 0;
        let old_value = self.incumbent_assignment[var_idx as usize];
        let delta = new_value - old_value;
        self.incumbent_assignment[var_idx as usize] = new_value;
        self.incumbent_objective += self.vars[var_idx as usize].objective_coeff * delta;

        // Update the LHSs of all involved constraints.
        for cstr_coeff in self.vars[var_idx as usize].coeffs.iter() {
            let cidx = cstr_coeff.idx as usize;
            let old_lhs = self.constraints[cidx].incumbent_lhs;
            let new_lhs = old_lhs + cstr_coeff.coeff * delta;
            self.constraints[cidx].incumbent_lhs = new_lhs;
            let new_cost = self.constraints[cidx].score(new_lhs);

            // Add/remove from the violatedConstraints list.
            if new_cost < -self.violation_tolerance && self.constraints[cidx].violated_idx == -1 {
                // Became violated.
                self.constraints[cidx].violated_idx = self.violated_constraints.len() as i32;
                self.violated_constraints.push(cstr_coeff.idx);
            }
            if new_cost >= -self.violation_tolerance && self.constraints[cidx].violated_idx != -1 {
                // Became satisfied.
                let last_violated_idx = self.violated_constraints.len() - 1;
                let last_constraint_idx = self.violated_constraints[last_violated_idx];
                let this_violated_idx = self.constraints[cidx].violated_idx;
                self.violated_constraints
                    .swap(this_violated_idx as usize, last_violated_idx);
                self.constraints[last_constraint_idx as usize].violated_idx = this_violated_idx;
                self.constraints[cidx].violated_idx = -1;
                self.violated_constraints.pop();
            }

            // Now, report the changes in LHS for other variables.
            let cstr = &self.constraints[cidx];
            dt += cstr.coeffs.len();
            for var_coeff in cstr.coeffs.iter() {
                if var_coeff.idx != var_idx {
                    let m = LhsModification {
                        var_idx: var_coeff.idx,
                        constraint_idx: cstr_coeff.idx,
                        coeff: var_coeff.coeff,
                        old_lhs,
                        new_lhs,
                    };
                    f(cstr, &self.incumbent_assignment, m);
                }
            }
        }

        dt
    }
}

impl JumpMove {
    /// HiGHS: feasibilityjump.hh::JumpMove::JumpMove (lines 319-319)
    /// C++: `JumpMove(double equalityTolerance)`
    pub fn new(equality_tolerance: f64) -> Self {
        Self {
            equality_tolerance,
            moves: Vec::new(),
            best_shift_buffer: Vec::new(),
            h16_fj: false,
        }
    }

    /// HiGHS: feasibilityjump.hh::JumpMove::init (lines 320-320)
    /// C++: `void init(Problem & problem)`
    pub fn init(&mut self, problem: &mut Problem) {
        self.moves.resize(problem.vars.len(), Move::default());
    }

    /// HiGHS: feasibilityjump.hh::JumpMove::forEachVarMove (lines 322-325)
    /// C++: `void forEachVarMove()`
    #[inline(always)]
    pub fn for_each_var_move<F: FnMut(&mut Move)>(&mut self, var_idx: i32, mut f: F) {
        f(&mut self.moves[var_idx as usize]);
    }

    /// HiGHS: feasibilityjump.hh::JumpMove::updateValue (lines 327-433)
    /// HiGHS 1.16 (#3310): feasibilityjump.hh::JumpMove::updateValue (lines 327-441)
    /// C++: `void updateValue(Problem & problem, uint32_t varIdx)`
    pub fn update_value(&mut self, problem: &mut Problem, var_idx: u32) {
        self.best_shift_buffer.clear();
        let var_incumbent_value = problem.incumbent_assignment[var_idx as usize];
        let lower = problem.vars[var_idx as usize].lb;
        let upper = problem.vars[var_idx as usize].ub;
        let mut current_value = lower; // can be -inf!

        let mut current_slope = 0.0;

        // PORT NOTE (perf): the per-cell C++ `std::vector constraintBounds` (at most three
        // entries) is a stack array, and the cell-invariant `vartype`, `residualIncumbent` and
        // `1.0 / coeff` are evaluated once per cell (same operands, so the same values).
        let is_integer = problem.vars[var_idx as usize].vartype == Integer;
        for cell in problem.vars[var_idx as usize].coeffs.iter() {
            let constraint = &problem.constraints[cell.idx as usize];

            let mut constraint_bounds = [(0.0f64, 0.0f64); 3];
            let num_bounds = if constraint.sense == Lte {
                constraint_bounds[0] = (-f64::INFINITY, constraint.rhs);
                1
            } else if constraint.sense == Gte {
                constraint_bounds[0] = (constraint.rhs, f64::INFINITY);
                1
            } else {
                constraint_bounds[0] = (-f64::INFINITY, constraint.rhs);
                // HiGHS 1.16 #3310 (h16_fj): HiGHS 1.16 drops the (rhs, rhs) bound, feasibilityjump.hh lines 350-355 (1.16.0)
                if !self.h16_fj {
                    constraint_bounds[1] = (constraint.rhs, constraint.rhs);
                    constraint_bounds[2] = (constraint.rhs, f64::INFINITY);
                    3
                } else {
                    constraint_bounds[1] = (constraint.rhs, f64::INFINITY);
                    2
                }
            };

            let residual_incumbent = constraint.incumbent_lhs - cell.coeff * var_incumbent_value;
            let inv_coeff = 1.0 / cell.coeff;
            for bound in constraint_bounds[..num_bounds].iter() {
                // HiGHS 1.16 #3310 (h16_fj): the bounds are swapped when cell.coeff < 0, feasibilityjump.hh lines 361-372 (1.16.0)
                let mut valid_range: (f64, f64) = if self.h16_fj {
                    let (bound_for_valid_lower, bound_for_valid_upper) = if cell.coeff < 0.0 {
                        (bound.1, bound.0)
                    } else {
                        (bound.0, bound.1)
                    };
                    (
                        inv_coeff * (bound_for_valid_lower - residual_incumbent),
                        inv_coeff * (bound_for_valid_upper - residual_incumbent),
                    )
                } else {
                    (
                        inv_coeff * (bound.0 - residual_incumbent),
                        inv_coeff * (bound.1 - residual_incumbent),
                    )
                };

                if is_integer {
                    valid_range = (
                        (valid_range.0 - self.equality_tolerance).ceil(),
                        (valid_range.1 + self.equality_tolerance).floor(),
                    );
                }

                // HiGHS 1.16 #3310 (h16_fj): HiGHS 1.16 has no inverted-range skip, feasibilityjump.hh lines 369-378 (1.16.0)
                if !self.h16_fj && valid_range.0 > valid_range.1 {
                    continue;
                }

                if valid_range.0 > current_value {
                    current_slope -= constraint.weight;
                    if valid_range.0 < upper {
                        self.best_shift_buffer
                            .push((valid_range.0, constraint.weight));
                    }
                }

                if valid_range.1 <= current_value {
                    current_slope += constraint.weight;
                } else if valid_range.1 < upper {
                    self.best_shift_buffer
                        .push((valid_range.1, constraint.weight));
                }
            }
        }

        if lower.is_finite() {
            self.best_shift_buffer.push((lower, 0.0));
        }
        if upper.is_finite() {
            self.best_shift_buffer.push((upper, 0.0));
        }
        // PORT NOTE: C++ std::sort on (double, double) pairs. pdqsort_by with the pair operator<
        // gives the same sequence (manifest: any correct sort; ties differ only in -0.0 sign).
        pdqsort_by(
            &mut self.best_shift_buffer,
            |a: &(f64, f64), b: &(f64, f64)| a.0 < b.0 || (!(b.0 < a.0) && a.1 < b.1),
        );

        // We only care about relative scores, so can pick an arbitrary start level
        let mut current_score = 0.0;
        current_value = if self.best_shift_buffer.is_empty() {
            var_incumbent_value
        } else {
            self.best_shift_buffer[0].0
        };

        let mut best_score = current_score;
        let mut best_value = current_value;

        for item in self.best_shift_buffer.iter() {
            current_score += (item.0 - current_value) * current_slope;
            current_slope += item.1;
            current_value = item.0;

            if eq(best_value, var_incumbent_value, self.equality_tolerance)
                || (!eq(current_value, var_incumbent_value, self.equality_tolerance)
                    && current_score < best_score)
            {
                best_score = current_score;
                best_value = current_value;
            }

            // Slope is always increasing, so if we have a valid value, we can quit
            // as soon as the slope turns nonnegative, since we must already have
            // visited the minimum.
            if !eq(best_value, var_incumbent_value, self.equality_tolerance) && current_slope >= 0.0
            {
                break;
            }
        }

        self.moves[var_idx as usize].value = best_value;
    }
}

impl FeasibilityJumpSolver {
    /// HiGHS: feasibilityjump.hh::FeasibilityJumpSolver::FeasibilityJumpSolver (lines 480-489)
    /// C++: `FeasibilityJumpSolver(const HighsLogOptions & _logOptions, int seed = 0, double equalityTolerance = 1e-5, double violationTolerance = 1e-5, double _weightUpdateDecay = 1.0)`
    /// C++ default arguments: `seed = 0`, `equality_tolerance = 1e-5`, `violation_tolerance = 1e-5`, `_weight_update_decay = 1.0`
    pub fn new(
        _log_options: &HighsLogOptions,
        seed: i32,
        equality_tolerance: f64,
        violation_tolerance: f64,
        _weight_update_decay: f64,
    ) -> Self {
        Self {
            problem: Problem::new(equality_tolerance, violation_tolerance),
            jump_move: JumpMove::new(equality_tolerance),
            weight_update_decay: _weight_update_decay,
            // C++ `rng = std::mt19937(seed)`; int seed wraps to unsigned like C++.
            rng: Mt19937::new(seed as u32),
            log_options: _log_options.clone(),
            ..Default::default()
        }
    }

    /// HiGHS: feasibilityjump.hh::FeasibilityJumpSolver::addVar (lines 491-494)
    /// C++: `int addVar(VarType vartype, double lb, double ub, double objCoeff)`
    pub fn add_var(&mut self, vartype: VarType, lb: f64, ub: f64, obj_coeff: f64) -> i32 {
        self.good_vars_set_idx.push(-1);
        self.problem.add_var(vartype, lb, ub, obj_coeff)
    }

    /// HiGHS: feasibilityjump.hh::FeasibilityJumpSolver::addConstraint (lines 496-500)
    /// C++: `int addConstraint(RowType sense, double rhs, int numCoeffs, int * rowVarIdxs, double * rowCoeffs, int relax_continuous)`
    pub fn add_constraint(
        &mut self,
        sense: RowType,
        rhs: f64,
        num_coeffs: i32,
        row_var_idxs: &mut [i32],
        row_coeffs: &mut [f64],
        relax_continuous: i32,
    ) -> i32 {
        self.problem.add_constraint(
            sense,
            rhs,
            num_coeffs,
            row_var_idxs,
            row_coeffs,
            relax_continuous,
        )
    }

    /// HiGHS: feasibilityjump.hh::FeasibilityJumpSolver::solve (lines 502-556)
    /// HiGHS 1.16 (#3269): feasibilityjump.hh::FeasibilityJumpSolver::solve (lines 510-567)
    /// C++: `int solve(double * initialValues, std::function<CallbackControlFlow (FJStatus)> callback)`
    pub fn solve(
        &mut self,
        initial_values: &mut [f64],
        callback: &mut dyn FnMut(FJStatus<'_>) -> CallbackControlFlow,
    ) -> i32 {
        crate::highs_log_dev!(
            &self.log_options,
            HighsLogType::kInfo,
            "Feasibility Jump: starting solve. weightUpdateDecay={}, relaxContinuous={}  \n",
            self.weight_update_decay,
            self.problem.used_relax_continuous as i32
        );

        self.init(initial_values);

        // PORT NOTE: size_t wrap of `-kMinEffortToLogging` (enables step=0 logging).
        self.effort_at_last_logging = 0usize.wrapping_sub(self.k_min_effort_to_logging);
        let mut num_logging_lines_since_header: i32 = 0;
        let mut step: i32 = 0;
        while step < i32::MAX {
            if self.user_terminate(callback, None) {
                break;
            }

            // PORT NOTE: size_t wrapping add, paired with the wrapping sub above.
            if step % self.k_logging_frequency == 0
                && self.total_effort
                    > self
                        .effort_at_last_logging
                        .wrapping_add(self.k_min_effort_to_logging)
            {
                if self.need_logging_header {
                    self.logging(0, true);
                    num_logging_lines_since_header = 0;
                }
                self.logging(step, false);
                num_logging_lines_since_header += 1;
                self.need_logging_header =
                    num_logging_lines_since_header == self.k_logging_header_frequency;
            }

            if self.problem.violated_constraints.len() < self.best_violation_score {
                self.effort_at_last_improvement = self.total_effort;
                self.best_violation_score = self.problem.violated_constraints.len();
            }

            if self.problem.violated_constraints.is_empty()
                && self.problem.incumbent_objective < self.best_objective
            {
                let relative_objective_improvement = (self.best_objective
                    - self.problem.incumbent_objective)
                    / crate::util::highs_utils::cpp_max(
                        1.0,
                        self.problem.incumbent_objective.abs(),
                    );
                if relative_objective_improvement > self.k_min_relative_objective_improvement {
                    self.effort_at_last_improvement = self.total_effort;
                    self.best_objective = self.problem.incumbent_objective;

                    // PORT NOTE: split borrow. The incumbent is taken out for the call and put
                    // back before the break. An empty incumbent (libc++ data() == nullptr) is None.
                    let mut incumbent = std::mem::take(&mut self.problem.incumbent_assignment);
                    let solution = if incumbent.is_empty() {
                        None
                    } else {
                        Some(incumbent.as_mut_slice())
                    };
                    let terminate = self.user_terminate(callback, solution);
                    self.problem.incumbent_assignment = incumbent;
                    if terminate {
                        break;
                    }
                    // Repeat the header in case user callback logs new solution
                    self.need_logging_header = true;
                }
            }

            if self.problem.vars.is_empty() {
                break;
            }

            let var = self.select_variable();
            // HiGHS 1.16 #3269 (h16_fj): HiGHS 1.16 breaks on UINT_MAX, feasibilityjump.hh lines 559-562 (1.16.0)
            if self.h16_fj && var == u32::MAX {
                break;
            }
            self.do_variable_move(var);
            self.ls_steps += 1;
            step += 1;
        }

        0
    }

    /// HiGHS: feasibilityjump.hh::FeasibilityJumpSolver::logging (lines 559-574)
    /// C++: `void logging(const int step, const bool header = false)`
    /// C++ default arguments: `header = false`
    pub fn logging(&mut self, step: i32, header: bool) {
        let log_type = HighsLogType::kDetailed;
        if header {
            crate::highs_log_dev!(
                &self.log_options,
                log_type,
                "Feasibility Jump:        step  violations     good    bumps       effort (per step)          Objective\n"
            );
        } else {
            let effort_per_step = if step > 0 {
                self.total_effort / step as usize
            } else {
                0
            };
            crate::highs_log_dev!(
                &self.log_options,
                log_type,
                // C++ `%10.4g` -> `{:>10.4}` (io/highs_io.rs printf convention).
                " {:>10}    {:>8}   {:>6} {:>8} {:>12}    {:>6}          {:>10.4}\n",
                step,
                self.problem.violated_constraints.len(),
                self.good_vars_set.len(),
                self.n_bumps,
                self.total_effort,
                effort_per_step,
                self.problem.incumbent_objective
            );
            self.effort_at_last_logging = self.total_effort;
        }
    }

    /// HiGHS: feasibilityjump.hh::FeasibilityJumpSolver::init (lines 576-584)
    /// C++: `void init(double * initialValues)`
    pub fn init(&mut self, initial_values: &mut [f64]) {
        self.problem.reset_incumbent(initial_values);
        self.jump_move.init(&mut self.problem);
        self.total_effort += self.problem.n_nonzeros;

        // Reset the variable scores.
        self.good_vars_set.clear();
        for i in 0..self.problem.vars.len() {
            self.reset_moves(i as u32);
        }
    }

    /// HiGHS: feasibilityjump.hh::FeasibilityJumpSolver::selectVariable (lines 586-638)
    /// HiGHS 1.16 (#3269): feasibilityjump.hh::FeasibilityJumpSolver::selectVariable (lines 597-648)
    /// C++: `uint32_t selectVariable()`
    pub fn select_variable(&mut self) -> u32 {
        if !self.good_vars_set.is_empty() {
            if self.rng.uniform_real(0.0, 1.0) < self.random_var_probability {
                let idx = (self.rng.next_u32() as usize) % self.good_vars_set.len();
                return self.good_vars_set[idx];
            }

            let sample_size = cpp_min(self.max_moves_to_evaluate, self.good_vars_set.len());
            self.total_effort += sample_size;
            self.ls_evals += sample_size as u64;

            let mut best_score = f64::NEG_INFINITY;
            let mut best_var: u32 = u32::MAX;
            for _i in 0..sample_size {
                let setidx = (self.rng.next_u32() as usize) % self.good_vars_set.len();
                let var_idx = self.good_vars_set[setidx];
                let mv = self.best_move(var_idx);
                if mv.score > best_score {
                    best_score = mv.score;
                    best_var = var_idx;
                }
            }
            // HiGHS 1.16 #3269 (h16_fj): HiGHS 1.16 drops this assert, feasibilityjump.hh lines 606-617 (1.16.0)
            if !self.h16_fj {
                debug_assert!(best_var != u32::MAX);
            }
            return best_var;
        }

        // Local minimum, update weights.
        self.update_weights();

        if !self.problem.violated_constraints.is_empty() {
            let n_violated = self.problem.violated_constraints.len();
            let cstr_idx = self.problem.violated_constraints
                [(self.rng.next_u32() as usize) % n_violated] as usize;

            if self.rng.uniform_real(0.0, 1.0) < self.random_cell_probability {
                let n_cells = self.problem.constraints[cstr_idx].coeffs.len();
                let k = (self.rng.next_u32() as usize) % n_cells;
                return self.problem.constraints[cstr_idx].coeffs[k].idx;
            }

            let mut best_score = f64::NEG_INFINITY;
            let mut best_var_idx: u32 = u32::MAX;

            let n_cells = self.problem.constraints[cstr_idx].coeffs.len();
            self.ls_evals += n_cells as u64;
            for k in 0..n_cells {
                let cell_idx = self.problem.constraints[cstr_idx].coeffs[k].idx;
                let mv = self.best_move(cell_idx);
                if mv.score > best_score {
                    best_score = mv.score;
                    best_var_idx = cell_idx;
                }
            }
            return best_var_idx;
        }

        // Fallback to random choice.
        ((self.rng.next_u32() as usize) % self.problem.vars.len()) as u32
    }

    /// HiGHS: feasibilityjump.hh::FeasibilityJumpSolver::updateWeights (lines 640-694)
    /// HiGHS 1.16 (#3310): feasibilityjump.hh::FeasibilityJumpSolver::updateWeights (lines 650-705)
    /// C++: `void updateWeights()`
    pub fn update_weights(&mut self) {
        crate::highs_log_dev!(
            &self.log_options,
            HighsLogType::kVerbose,
            "Feasibility Jump: Reached a local minimum.\n"
        );
        self.n_bumps += 1;
        let mut rescale_all_weights = false;
        let mut dt: usize = 0;

        if self.problem.violated_constraints.is_empty() {
            self.objective_weight += self.weight_update_increment;
            if self.objective_weight > 1.0e20 {
                rescale_all_weights = true;
            }

            dt += self.problem.vars.len();
            // HiGHS 1.16 #3310 (h16_fj): read the flag once, before the move loop.
            let h16_fj = self.h16_fj;
            for var_idx in 0..self.problem.vars.len() {
                self.jump_move
                    .for_each_var_move(var_idx as i32, &mut |m: &mut Move| {
                        if h16_fj {
                            // HiGHS 1.16 #3310 (h16_fj): `-=` to align objective minimization with score maximization, feasibilityjump.hh lines 664-667 (1.16.0)
                            m.score -= self.weight_update_increment
                                * self.problem.vars[var_idx].objective_coeff
                                * (m.value - self.problem.incumbent_assignment[var_idx]);
                        } else {
                            m.score += self.weight_update_increment
                                * self.problem.vars[var_idx].objective_coeff
                                * (m.value - self.problem.incumbent_assignment[var_idx]);
                        }
                    });
            }
        } else {
            for k in 0..self.problem.violated_constraints.len() {
                let c_idx = self.problem.violated_constraints[k] as usize;
                self.problem.constraints[c_idx].weight += self.weight_update_increment;
                if self.problem.constraints[c_idx].weight > 1.0e20 {
                    rescale_all_weights = true;
                }

                dt += self.problem.constraints[c_idx].coeffs.len();
                for j in 0..self.problem.constraints[c_idx].coeffs.len() {
                    let cell = self.problem.constraints[c_idx].coeffs[j];
                    self.jump_move
                        .for_each_var_move(cell.idx as i32, &mut |m: &mut Move| {
                            let constraint = &mut self.problem.constraints[c_idx];
                            let inc_lhs = constraint.incumbent_lhs;
                            let candidate_lhs = inc_lhs
                                + cell.coeff
                                    * (m.value
                                        - self.problem.incumbent_assignment[cell.idx as usize]);
                            let diff = self.weight_update_increment
                                * (constraint.score(candidate_lhs) - constraint.score(inc_lhs));
                            m.score += diff;
                        });

                    self.update_good_moves(cell.idx as i32);
                }
            }
        }

        self.weight_update_increment /= self.weight_update_decay;
        if rescale_all_weights {
            self.weight_update_increment *= 1.0e-20;
            self.objective_weight *= 1.0e-20;

            for cc in self.problem.constraints.iter_mut() {
                cc.weight *= 1.0e-20;
            }
            dt += self.problem.constraints.len();

            for i in 0..self.problem.vars.len() {
                self.reset_moves(i as u32);
            }
        }

        self.total_effort += dt;
    }

    /// HiGHS: feasibilityjump.hh::FeasibilityJumpSolver::bestMove (lines 696-703)
    /// C++: `Move bestMove(uint32_t varIdx)`
    pub fn best_move(&mut self, var_idx: u32) -> Move {
        let mut best = Move::undef();
        self.jump_move
            .for_each_var_move(var_idx as i32, &mut |m: &mut Move| {
                if m.score > best.score {
                    best = *m;
                }
            });
        best
    }

    /// HiGHS: feasibilityjump.hh::FeasibilityJumpSolver::doVariableMove (lines 705-721)
    /// C++: `void doVariableMove(uint32_t varIdx)`
    pub fn do_variable_move(&mut self, var_idx: u32) {
        // First, we get the best move for the variable;
        let m = self.best_move(var_idx);
        let new_value = m.value;

        // PORT NOTE (perf): C++ calls modifyMove(mod, problem, m) and updateGoodMoves(mod.varIdx)
        // from inside setValue's callback; so does this, on split borrows of the solver:
        // `set_value_with` hands the callback the constraint and the incumbent that modifyMove
        // reads, and the callback owns the moves and the good-variable set. Same calls in the
        // same order as C++ (the buffered replay this replaces gave the same values but grew a
        // Vec on every move).
        let moves = &mut self.jump_move.moves;
        let good_vars_set = &mut self.good_vars_set;
        let good_vars_set_idx = &mut self.good_vars_set_idx;
        let dt = self.problem.set_value_with(
            var_idx,
            new_value,
            |cn: &Constraint, incumbent_assignment: &[f64], lm: LhsModification| {
                // forEachMove(mod.varIdx, modifyMove): the single JumpMove move of the variable
                modify_move_in(
                    lm,
                    cn,
                    incumbent_assignment,
                    &mut moves[lm.var_idx as usize],
                );
                update_good_moves_in(moves, good_vars_set, good_vars_set_idx, lm.var_idx as i32);
            },
        );
        self.total_effort += dt;

        self.reset_moves(var_idx);
    }

    /// HiGHS: feasibilityjump.hh::FeasibilityJumpSolver::updateGoodMoves (lines 723-739)
    /// C++: `void updateGoodMoves(int32_t varIdx)`
    pub fn update_good_moves(&mut self, var_idx: i32) {
        update_good_moves_in(
            &self.jump_move.moves,
            &mut self.good_vars_set,
            &mut self.good_vars_set_idx,
            var_idx,
        );
    }

    /// HiGHS: feasibilityjump.hh::FeasibilityJumpSolver::forEachMove (lines 741-747)
    /// C++: `void forEachMove()`
    #[inline(always)]
    pub fn for_each_move<F: FnMut(&mut Move)>(&mut self, var_idx: i32, f: F) {
        self.jump_move.for_each_var_move(var_idx, f);
    }

    /// HiGHS: feasibilityjump.hh::FeasibilityJumpSolver::resetMoves (lines 749-770)
    /// HiGHS 1.16 (#3310): feasibilityjump.hh::FeasibilityJumpSolver::resetMoves (lines 760-782)
    /// C++: `void resetMoves(uint32_t varIdx)`
    pub fn reset_moves(&mut self, var_idx: u32) {
        self.total_effort += self.problem.vars[var_idx as usize].coeffs.len();
        self.jump_move.update_value(&mut self.problem, var_idx);

        let vi = var_idx as usize;
        // HiGHS 1.16 #3310 (h16_fj): read the flag once, before the move closure.
        let h16_fj = self.h16_fj;
        self.jump_move
            .for_each_var_move(var_idx as i32, &mut |mv: &mut Move| {
                mv.score = 0.0;
                if h16_fj {
                    // HiGHS 1.16 #3310 (h16_fj): `-=` to align objective minimization with score maximization, feasibilityjump.hh lines 764-768 (1.16.0)
                    mv.score -= self.objective_weight
                        * self.problem.vars[vi].objective_coeff
                        * (mv.value - self.problem.incumbent_assignment[vi]);
                } else {
                    mv.score += self.objective_weight
                        * self.problem.vars[vi].objective_coeff
                        * (mv.value - self.problem.incumbent_assignment[vi]);
                }

                for cell in self.problem.vars[vi].coeffs.iter() {
                    let constraint = &mut self.problem.constraints[cell.idx as usize];
                    let inc_lhs = constraint.incumbent_lhs;
                    let candidate_lhs =
                        inc_lhs + cell.coeff * (mv.value - self.problem.incumbent_assignment[vi]);
                    mv.score += constraint.weight
                        * (constraint.score(candidate_lhs) - constraint.score(inc_lhs));
                }
            });

        self.update_good_moves(var_idx as i32);
    }

    /// HiGHS: feasibilityjump.hh::FeasibilityJumpSolver::user_terminate (lines 772-797)
    /// C++: `bool user_terminate(std::function<CallbackControlFlow (FJStatus)> callback, double * solution)`
    pub fn user_terminate(
        &mut self,
        callback: &mut dyn FnMut(FJStatus<'_>) -> CallbackControlFlow,
        solution: Option<&mut [f64]>,
    ) -> bool {
        const CALLBACK_EFFORT: usize = 500000; // Originally 500000
        if solution.is_some()
            || self.total_effort.wrapping_sub(self.effort_at_last_callback) > CALLBACK_EFFORT
        {
            crate::highs_log_dev!(
                &self.log_options,
                HighsLogType::kVerbose,
                "Feasibility Jump: calling user termination.\n"
            );
            self.effort_at_last_callback = self.total_effort;

            let status = FJStatus {
                total_effort: self.total_effort,
                effort_since_last_improvement: self
                    .total_effort
                    .wrapping_sub(self.effort_at_last_improvement),
                solution: solution.as_deref(),
                num_vars: self.problem.vars.len() as i32,
                solution_objective_value: self.problem.incumbent_objective,
            };

            let result = callback(status);
            if result == Terminate {
                crate::highs_log_dev!(
                    &self.log_options,
                    HighsLogType::kVerbose,
                    "Feasibility Jump: quitting.\n"
                );
                return true;
            }
        }
        false
    }
}

#[cfg(test)]
mod h16_fj_tests {
    use super::*;

    const TOL: f64 = 1e-5;

    // Expected values are from the C++ driver built against highs/mip/feasibilityjump.hh of
    // HiGHS 1.16.0 (on) and 1.15.1 (off), with a single column and a single row.

    /// One column (index 0) and one row with coefficient `coeff`, then one `updateValue`.
    fn one_var_jump(
        vartype: VarType,
        lb: f64,
        ub: f64,
        sense: RowType,
        rhs: f64,
        coeff: f64,
        init: f64,
        flag: bool,
    ) -> JumpMove {
        let mut problem = Problem::new(TOL, TOL);
        problem.add_var(vartype, lb, ub, 0.0);
        let mut idx = [0i32];
        let mut cf = [coeff];
        problem.add_constraint(sense, rhs, 1, &mut idx, &mut cf, 0);
        problem.reset_incumbent(&mut [init]);
        let mut jm = JumpMove::new(TOL);
        jm.h16_fj = flag;
        jm.init(&mut problem);
        jm.update_value(&mut problem, 0);
        jm
    }

    /// Solver with one column (objective coeff 2, [0, 1]) whose moves are reset once.
    fn objective_solver(flag: bool) -> FeasibilityJumpSolver {
        let mut s = FeasibilityJumpSolver::new(&HighsLogOptions::default(), 0, TOL, TOL, 1.0);
        s.h16_fj = flag;
        s.jump_move.h16_fj = flag;
        s.add_var(Integer, 0.0, 1.0, 2.0);
        s.problem.reset_incumbent(&mut [0.0]);
        s.jump_move.init(&mut s.problem);
        s.objective_weight = 1.0;
        s.reset_moves(0);
        s
    }

    /// Runs `solve` on one column with objective `objc` and returns the first reported solution.
    fn first_solution(
        flag: bool,
        objc: f64,
        ub: f64,
        sense: RowType,
        rhs: f64,
        coeff: f64,
    ) -> Vec<f64> {
        let mut s = FeasibilityJumpSolver::new(&HighsLogOptions::default(), 0, TOL, TOL, 1.0);
        s.h16_fj = flag;
        s.jump_move.h16_fj = flag;
        s.add_var(Integer, 0.0, ub, objc);
        let mut idx = [0i32];
        let mut cf = [coeff];
        s.add_constraint(sense, rhs, 1, &mut idx, &mut cf, 0);
        let mut init = [0.0f64];
        let mut found: Vec<f64> = Vec::new();
        s.solve(&mut init, &mut |st: FJStatus<'_>| {
            if let Some(sol) = st.solution {
                found = sol[..st.num_vars as usize].to_vec();
                Terminate
            } else {
                Continue
            }
        });
        found
    }

    #[test]
    fn h16_fj_on_lte_negative_coeff_jumps_to_feasible_value() {
        // 1.16.0 C++: move=1, buffer (0,0) (1,1) (5,0)
        let jm = one_var_jump(Integer, 0.0, 5.0, Lte, -1.0, -1.0, 0.0, true);
        assert_eq!(jm.moves[0].value, 1.0);
        assert_eq!(
            jm.best_shift_buffer,
            vec![(0.0, 0.0), (1.0, 1.0), (5.0, 0.0)]
        );
    }

    #[test]
    fn h16_fj_off_lte_negative_coeff_keeps_1_15_1_jump() {
        // 1.15.1 C++: move=5, buffer (0,0) (5,0): the inverted range is skipped
        let jm = one_var_jump(Integer, 0.0, 5.0, Lte, -1.0, -1.0, 0.0, false);
        assert_eq!(jm.moves[0].value, 5.0);
        assert_eq!(jm.best_shift_buffer, vec![(0.0, 0.0), (5.0, 0.0)]);
    }

    #[test]
    fn h16_fj_on_equality_negative_coeff_value() {
        // 1.16.0 C++: move=1, buffer (0,0) (1,0)
        let jm = one_var_jump(Integer, 0.0, 1.0, Equal, -1.0, -1.0, 0.0, true);
        assert_eq!(jm.moves[0].value, 1.0);
        assert_eq!(jm.best_shift_buffer, vec![(0.0, 0.0), (1.0, 0.0)]);
    }

    #[test]
    fn h16_fj_off_equality_negative_coeff_value() {
        // 1.15.1 C++: move=1, buffer (0,0) (1,0)
        let jm = one_var_jump(Integer, 0.0, 1.0, Equal, -1.0, -1.0, 0.0, false);
        assert_eq!(jm.moves[0].value, 1.0);
        assert_eq!(jm.best_shift_buffer, vec![(0.0, 0.0), (1.0, 0.0)]);
    }

    #[test]
    fn h16_fj_on_equality_positive_coeff_drops_duplicate_bound() {
        // 1.16.0 C++: move=2, buffer (0,0) (2,1) (2,1) (5,0)
        let jm = one_var_jump(Integer, 0.0, 5.0, Equal, 2.0, 1.0, 0.0, true);
        assert_eq!(jm.moves[0].value, 2.0);
        assert_eq!(
            jm.best_shift_buffer,
            vec![(0.0, 0.0), (2.0, 1.0), (2.0, 1.0), (5.0, 0.0)]
        );
    }

    #[test]
    fn h16_fj_off_equality_positive_coeff_three_bounds() {
        // 1.15.1 C++: move=2, buffer (0,0) (2,1) x4 (5,0)
        let jm = one_var_jump(Integer, 0.0, 5.0, Equal, 2.0, 1.0, 0.0, false);
        assert_eq!(jm.moves[0].value, 2.0);
        assert_eq!(
            jm.best_shift_buffer,
            vec![
                (0.0, 0.0),
                (2.0, 1.0),
                (2.0, 1.0),
                (2.0, 1.0),
                (2.0, 1.0),
                (5.0, 0.0)
            ]
        );
    }

    #[test]
    fn h16_fj_on_reset_moves_objective_sign() {
        // 1.16.0 C++ resetMoves: `score -= objectiveWeight * objCoeff * (value - incumbent)`
        let s = objective_solver(true);
        assert_eq!(s.jump_move.moves[0].value, 1.0);
        assert_eq!(s.jump_move.moves[0].score, -2.0);
    }

    #[test]
    fn h16_fj_off_reset_moves_objective_sign() {
        // 1.15.1 C++ resetMoves: `score += ...`
        let s = objective_solver(false);
        assert_eq!(s.jump_move.moves[0].value, 1.0);
        assert_eq!(s.jump_move.moves[0].score, 2.0);
    }

    #[test]
    fn h16_fj_on_update_weights_objective_sign() {
        // 1.16.0 C++ updateWeights (no violated rows): `score -= 1 * 2 * (1 - 0)` after resetMoves
        let mut s = objective_solver(true);
        s.update_weights();
        assert_eq!(s.objective_weight, 2.0);
        assert_eq!(s.jump_move.moves[0].score, -4.0);
    }

    #[test]
    fn h16_fj_off_update_weights_objective_sign() {
        // 1.15.1 C++ updateWeights: `score += 1 * 2 * (1 - 0)` after resetMoves
        let mut s = objective_solver(false);
        s.update_weights();
        assert_eq!(s.objective_weight, 2.0);
        assert_eq!(s.jump_move.moves[0].score, 4.0);
    }

    #[test]
    fn h16_fj_on_solve_finds_feasible_point() {
        // 1.16.0 C++ solve: first reported solution is 1
        let found = first_solution(true, 1.0, 5.0, Lte, -1.0, -1.0);
        assert_eq!(found, vec![1.0]);
    }

    #[test]
    fn h16_fj_off_solve_keeps_1_15_1_point() {
        // 1.15.1 C++ solve: first reported solution is 5 (also feasible)
        let found = first_solution(false, 1.0, 5.0, Lte, -1.0, -1.0);
        assert_eq!(found, vec![5.0]);
    }
}
