//! HiGHS source: lp_data/HighsRanging.cpp (lines 70-601)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! `getRangingData` is a single 532-line function; it is kept whole
//! rather than split mid-function. Use `cpp_max`/`cpp_min` for every
//! `max`/`min` (the C++ uses `using std::max/min`), `inf_product` /
//! `poss_inf_product` from the parent module, and the `HEkk` passed in for
//! every `ekk_instance_` access (`solver_object.ekk_instance_` in C++).

use super::*;
use crate::prelude::*;
use crate::simplex::nla::NlaRefs;

/// HiGHS: HighsRanging.cpp::getRangingData (lines 70-601)
///
/// PORT NOTE: C++ `solver_object.ekk_instance_` is the `ekk` parameter.
pub fn get_ranging_data(
    ranging: &mut HighsRanging,
    solver_object: &mut HighsLpSolverObject,
    ekk: &mut HEkk,
) -> HighsStatus {
    ranging.clear();
    if solver_object.model_status_ != HighsModelStatus::kOptimal {
        crate::highs_log_user!(
            &solver_object.options_.log_options,
            HighsLogType::kError,
            "Cannot get ranging without an optimal solution\n"
        );
        return HighsStatus::kError;
    }
    if !ekk.status_.initialised_for_solve {
        crate::highs_log_user!(
            &solver_object.options_.log_options,
            HighsLogType::kError,
            "Cannot get ranging without a valid Simplex instance\n"
        );
        return HighsStatus::kError;
    }
    let use_lp = &*solver_object.lp_;

    // Unscale the simplex data if the LP has been solved in the scaled space
    ekk.unscale_simplex(use_lp);

    // Aliases
    let ekk: &HEkk = ekk;
    let simplex_info = &ekk.info_;
    let simplex_basis = &ekk.basis_;
    let value_ = &simplex_info.work_value_;
    let dual_ = &simplex_info.work_dual_;
    let cost_ = &simplex_info.work_cost_;
    let lower_ = &simplex_info.work_lower_;
    let upper_ = &simplex_info.work_upper_;
    let bvalue_ = &simplex_info.base_value_;
    let blower_ = &simplex_info.base_lower_;
    let bupper_ = &simplex_info.base_upper_;
    let nflag_ = &simplex_basis.nonbasic_flag_;
    let nmove_ = &simplex_basis.nonbasic_move_;
    let bindex_ = &simplex_basis.basic_index_;
    let matrix = &use_lp.a_matrix_;

    // Local copies of scalars

    let num_row = use_lp.num_row_;
    let num_col = use_lp.num_col_;
    let num_total = num_col + num_row;
    let h_tt: f64 = 1e-13;
    let h_inf: f64 = kHighsInf;
    let objective = solver_object.highs_info_.objective_function_value;

    // Code written for minimization problems. Maximization problems are
    // solved by using negated costs in the simplex solver and
    // minimizing. Thus dual information in the simplex solver is
    // negated for maximization problems. The objective has the right
    // sign, though. Maximization problems are, thus, accommodated by
    // applying the sign multiplier to dual information.
    let mut sense: HighsInt = 1;
    if use_lp.sense_ == ObjSense::kMaximize {
        sense = -1;
    }

    let mut i_work_ = vec![0 as HighsInt; num_total as usize];
    let mut d_work_ = vec![0.0f64; num_total as usize];
    let mut column = HVector::default();
    column.setup(num_row);

    let mut xi = bvalue_.clone();
    for i in 0..num_row {
        xi[i as usize] = cpp_max(xi[i as usize], blower_[i as usize]);
        xi[i as usize] = cpp_min(xi[i as usize], bupper_[i as usize]);
    }
    let mut dj = dual_.clone();
    for j in 0..num_total {
        if nflag_[j as usize] != 0 && (lower_[j as usize] != upper_[j as usize]) {
            if value_[j as usize] == lower_[j as usize] {
                dj[j as usize] = cpp_max(dj[j as usize], 0.0);
            }
            if value_[j as usize] == upper_[j as usize] {
                dj[j as usize] = cpp_min(dj[j as usize], 0.0);
            }
            if lower_[j as usize] == -h_inf && upper_[j as usize] == h_inf {
                dj[j as usize] = 0.0;
            }
        }
    }
    //
    // Ranging 1.2. prepare "delta" space
    //
    let mut dxi_inc = vec![0.0f64; num_row as usize];
    let mut dxi_dec = vec![0.0f64; num_row as usize];
    for i in 0..num_row {
        dxi_inc[i as usize] = bupper_[i as usize] - xi[i as usize];
        dxi_dec[i as usize] = blower_[i as usize] - xi[i as usize];
    }

    let mut ddj_inc = vec![0.0f64; num_total as usize];
    let mut ddj_dec = vec![0.0f64; num_total as usize];
    for j in 0..num_total {
        if nflag_[j as usize] != 0 {
            ddj_inc[j as usize] = if value_[j as usize] == lower_[j as usize] {
                h_inf
            } else {
                -dj[j as usize]
            };
            ddj_dec[j as usize] = if value_[j as usize] == upper_[j as usize] {
                -h_inf
            } else {
                -dj[j as usize]
            };
        }
    }

    //
    // Ranging 1.3. prepare "theta" space
    //
    let tol_a: f64 = 1e-9;
    let theta_inf: f64 = h_inf / 1e40;

    let mut txj_inc = vec![theta_inf; num_total as usize]; // theta
    let mut axj_inc = vec![0.0f64; num_total as usize]; // alpha
    let mut ixj_inc = vec![-1 as HighsInt; num_total as usize]; // i-out
    let mut wxj_inc = vec![0 as HighsInt; num_total as usize]; // which bound is limiting
    let mut jxj_inc = vec![-1 as HighsInt; num_total as usize]; // j = n(i), (with bound flip)

    let mut txj_dec = vec![-theta_inf; num_total as usize];
    let mut axj_dec = vec![0.0f64; num_total as usize];
    let mut ixj_dec = vec![-1 as HighsInt; num_total as usize];
    let mut wxj_dec = vec![0 as HighsInt; num_total as usize];
    let mut jxj_dec = vec![-1 as HighsInt; num_total as usize];

    let mut tci_inc = vec![theta_inf; num_row as usize]; // theta
    let mut aci_inc = vec![0.0f64; num_row as usize]; // alpha
    let mut jci_inc = vec![-1 as HighsInt; num_row as usize]; // column index

    let mut tci_dec = vec![-theta_inf; num_row as usize];
    let mut aci_dec = vec![0.0f64; num_row as usize];
    let mut jci_dec = vec![-1 as HighsInt; num_row as usize];

    // Major "theta" loop
    for j in 0..num_total {
        // Skip basic column
        if nflag_[j as usize] == 0 {
            continue;
        }

        // Form updated column
        column.clear();
        matrix.collect_aj(&mut column, j, 1.0);

        let expected_density = ekk.info_.col_aq_density;
        // PORT NOTE: HEkk::ftran (HEkk.cpp 402-405) is
        // `simplex_nla_.ftran(rhs, expected_density)`. After the solve, C++
        // `simplex_nla_` points at the incumbent LP (HApp.h 499
        // `setNlaPointersForLpAndScale(incumbent_lp)`), i.e. `use_lp`, whose
        // scale_ is applied; `ekk.lp_` has been moved back out, so
        // `ekk.ftran` would skip the scaling (lib.rs rule (b), highs_api.rs).
        debug_assert!(ekk.status_.has_nla);
        let refs = NlaRefs::new(use_lp, &ekk.basis_.basic_index_, &ekk.options_);
        ekk.simplex_nla_
            .ftran(&refs, &mut column, expected_density, None);
        let mut n_work: HighsInt = 0;
        for k in 0..column.count {
            let i_row = column.index[k as usize];
            let alpha = column.array[i_row as usize];
            if alpha.abs() > tol_a {
                i_work_[n_work as usize] = i_row;
                d_work_[n_work as usize] = alpha;
                n_work += 1;
            }
        }
        // Standard primal ratio test
        let mut myt_inc: f64 = theta_inf;
        let mut myt_dec: f64 = -theta_inf;
        let mut myk_inc: HighsInt = -1;
        let mut myk_dec: HighsInt = -1;
        for k in 0..n_work {
            let i = i_work_[k as usize];
            let alpha = d_work_[k as usize];
            let theta_inc = if alpha < 0.0 {
                dxi_inc[i as usize]
            } else {
                dxi_dec[i as usize]
            } / -alpha;
            let theta_dec = if alpha > 0.0 {
                dxi_inc[i as usize]
            } else {
                dxi_dec[i as usize]
            } / -alpha;
            if myt_inc > theta_inc {
                myt_inc = theta_inc;
                myk_inc = k;
            }
            if myt_dec < theta_dec {
                myt_dec = theta_dec;
                myk_dec = k;
            }
        }

        if myk_inc != -1 {
            let i = i_work_[myk_inc as usize];
            let alpha = d_work_[myk_inc as usize];
            ixj_inc[j as usize] = i;
            axj_inc[j as usize] = alpha;
            txj_inc[j as usize] = if alpha < 0.0 {
                dxi_inc[i as usize]
            } else {
                dxi_dec[i as usize]
            } / -alpha;
            wxj_inc[j as usize] = if alpha < 0.0 { 1 } else { -1 };
        }

        if myk_dec != -1 {
            let i = i_work_[myk_dec as usize];
            let alpha = d_work_[myk_dec as usize];
            ixj_dec[j as usize] = i;
            axj_dec[j as usize] = alpha;
            txj_dec[j as usize] = if alpha > 0.0 {
                dxi_inc[i as usize]
            } else {
                dxi_dec[i as usize]
            } / -alpha;
            wxj_dec[j as usize] = if alpha > 0.0 { 1 } else { -1 };
        }

        // Accumulated dual ratio test
        let myd_inc = ddj_inc[j as usize];
        let myd_dec = ddj_dec[j as usize];
        for k in 0..n_work {
            let i = i_work_[k as usize];
            let alpha = d_work_[k as usize];
            let theta_inc = if alpha < 0.0 { myd_inc } else { myd_dec } / -alpha;
            let theta_dec = if alpha > 0.0 { myd_inc } else { myd_dec } / -alpha;
            if tci_inc[i as usize] > theta_inc {
                tci_inc[i as usize] = theta_inc;
                aci_inc[i as usize] = alpha;
                jci_inc[i as usize] = j;
            }
            if tci_dec[i as usize] < theta_dec {
                tci_dec[i as usize] = theta_dec;
                aci_dec[i as usize] = alpha;
                jci_dec[i as usize] = j;
            }
        }
    }

    // Additional j-out for primal ratio test (considering bound flip)
    for j in 0..num_total {
        if nflag_[j as usize] != 0 {
            // J-out for x_j = l_j
            if nmove_[j as usize] == 1 {
                let value = value_[j as usize] + txj_inc[j as usize];
                if ixj_inc[j as usize] != -1 && value <= upper_[j as usize] {
                    jxj_inc[j as usize] = bindex_[ixj_inc[j as usize] as usize];
                } else if value > upper_[j as usize] {
                    jxj_inc[j as usize] = j;
                }
            }
            // J-out for x_j = u_j
            if nmove_[j as usize] == -1 {
                let value = value_[j as usize] + txj_dec[j as usize];
                if ixj_dec[j as usize] != -1 && value >= lower_[j as usize] {
                    jxj_dec[j as usize] = bindex_[ixj_dec[j as usize] as usize];
                } else if value < lower_[j as usize] {
                    jxj_dec[j as usize] = j;
                }
            }
            // J-out for free variable
            if lower_[j as usize] == -h_inf && upper_[j as usize] == h_inf {
                // C++ chained assignment `jxj_inc[j] = jxj_dec[j] = Bindex_[..]`:
                // both receive Bindex_[..].
                if ixj_inc[j as usize] != -1 {
                    jxj_dec[j as usize] = bindex_[ixj_inc[j as usize] as usize];
                    jxj_inc[j as usize] = jxj_dec[j as usize];
                }
                if ixj_dec[j as usize] != -1 {
                    jxj_dec[j as usize] = bindex_[ixj_dec[j as usize] as usize];
                    jxj_inc[j as usize] = jxj_dec[j as usize];
                }
            }
        }
    }

    //
    // Ranging 2. cost ranging
    //
    let mut c_up_c = vec![0.0f64; num_total as usize];
    let mut c_dn_c = vec![0.0f64; num_total as usize];
    let mut c_up_f = vec![0.0f64; num_total as usize];
    let mut c_dn_f = vec![0.0f64; num_total as usize];
    let mut c_up_e = vec![0 as HighsInt; num_total as usize];
    let mut c_dn_e = vec![0 as HighsInt; num_total as usize];
    let mut c_up_l = vec![0 as HighsInt; num_total as usize];
    let mut c_dn_l = vec![0 as HighsInt; num_total as usize];

    //
    // Ranging 2.1. nonbasic cost ranging
    //
    //  const HighsInt check_col = 2951;
    for j in 0..num_col {
        if nflag_[j as usize] != 0 {
            // Primal value and its sign
            let value = value_[j as usize];
            let vsign: f64 = if value > 0.0 {
                1.0
            } else if value < 0.0 {
                -1.0
            } else {
                0.0
            };

            // Increase c_j
            if ddj_inc[j as usize] != h_inf {
                c_up_c[j as usize] = cost_[j as usize] + ddj_inc[j as usize];
                c_up_f[j as usize] =
                    objective + sense as f64 * poss_inf_product(ddj_inc[j as usize], value);
                c_up_e[j as usize] = j;
                c_up_l[j as usize] = jxj_dec[j as usize];
            } else {
                c_up_c[j as usize] = h_inf;
                c_up_f[j as usize] = objective + sense as f64 * inf_product(vsign);
                c_up_e[j as usize] = -1;
                c_up_l[j as usize] = -1;
            }

            // Decrease c_j
            if ddj_dec[j as usize] != h_inf {
                c_dn_c[j as usize] = cost_[j as usize] + ddj_dec[j as usize];
                c_dn_f[j as usize] =
                    objective + sense as f64 * poss_inf_product(ddj_dec[j as usize], value);
                c_dn_e[j as usize] = j;
                c_dn_l[j as usize] = jxj_inc[j as usize];
            } else {
                c_up_c[j as usize] = -h_inf;
                c_up_f[j as usize] = objective - sense as f64 * inf_product(vsign);
                c_up_e[j as usize] = -1;
                c_up_l[j as usize] = -1;
            }
        }
    }

    //
    // Ranging 2.2. basic cost ranging
    //

    for i in 0..num_row {
        if bindex_[i as usize] < num_col {
            // Primal variable and its sign
            let j = bindex_[i as usize];
            let mut je: HighsInt;
            let value = xi[i as usize];
            let vsign: f64 = if value > 0.0 {
                1.0
            } else if value < 0.0 {
                -1.0
            } else {
                0.0
            };

            // Increase c_i
            if jci_inc[i as usize] != -1 {
                c_up_c[j as usize] = cost_[j as usize] + tci_inc[i as usize];
                c_up_f[j as usize] =
                    objective + sense as f64 * poss_inf_product(tci_inc[i as usize], value);
                je = jci_inc[i as usize];
                c_up_e[j as usize] = je;
                c_up_l[j as usize] = if nmove_[je as usize] > 0 {
                    jxj_inc[je as usize]
                } else {
                    jxj_dec[je as usize]
                };
            } else {
                c_up_c[j as usize] = h_inf;
                c_up_f[j as usize] = objective + sense as f64 * inf_product(vsign);
                c_up_e[j as usize] = -1;
                c_up_l[j as usize] = -1;
            }

            // Decrease c_i
            if jci_dec[i as usize] != -1 {
                c_dn_c[j as usize] = cost_[j as usize] + tci_dec[i as usize];
                c_dn_f[j as usize] =
                    objective + sense as f64 * poss_inf_product(tci_dec[i as usize], value);
                je = jci_dec[i as usize];
                c_dn_e[j as usize] = je;
                c_dn_l[j as usize] = if nmove_[je as usize] > 0 {
                    jxj_inc[je as usize]
                } else {
                    jxj_dec[je as usize]
                };
            } else {
                c_dn_c[j as usize] = -h_inf;
                c_dn_f[j as usize] = objective - sense as f64 * inf_product(vsign);
                c_dn_e[j as usize] = -1;
                c_dn_l[j as usize] = -1;
            }
        }
    }

    //
    // Ranging 3. bounds ranging
    //
    let mut b_up_b = vec![0.0f64; num_total as usize];
    let mut b_dn_b = vec![0.0f64; num_total as usize];
    let mut b_up_f = vec![0.0f64; num_total as usize];
    let mut b_dn_f = vec![0.0f64; num_total as usize];
    let mut b_up_e = vec![0 as HighsInt; num_total as usize];
    let mut b_dn_e = vec![0 as HighsInt; num_total as usize];
    let mut b_up_l = vec![0 as HighsInt; num_total as usize];
    let mut b_dn_l = vec![0 as HighsInt; num_total as usize];

    //
    // Ranging 3.1. nonbasic bounds ranging
    //
    for j in 0..num_total {
        if nflag_[j as usize] != 0 {
            // FREE variable
            if lower_[j as usize] == -h_inf && upper_[j as usize] == h_inf {
                b_up_b[j as usize] = h_inf;
                b_up_f[j as usize] = objective;
                b_up_e[j as usize] = -1;
                b_up_l[j as usize] = -1;
                b_dn_b[j as usize] = -h_inf;
                b_dn_f[j as usize] = objective;
                b_dn_e[j as usize] = -1;
                b_dn_l[j as usize] = -1;
                continue;
            }

            // Dual value and its sign
            let dualv = dj[j as usize];
            let dsign: f64 = if dualv > 0.0 {
                1.0
            } else if dualv < 0.0 {
                -1.0
            } else {
                0.0
            };

            // Increase x_j
            if ixj_inc[j as usize] != -1 {
                let i = ixj_inc[j as usize];
                b_up_b[j as usize] = value_[j as usize] + txj_inc[j as usize];
                b_up_f[j as usize] =
                    objective + sense as f64 * poss_inf_product(txj_inc[j as usize], dualv);
                b_up_e[j as usize] = if wxj_inc[j as usize] > 0 {
                    jci_inc[i as usize]
                } else {
                    jci_dec[i as usize]
                };
                b_up_l[j as usize] = bindex_[i as usize];
            } else {
                b_up_b[j as usize] = h_inf;
                b_up_f[j as usize] = objective + sense as f64 * inf_product(dsign);
                b_up_e[j as usize] = -1;
                b_up_l[j as usize] = -1;
            }

            // Check if b_up_b > upper
            if value_[j as usize] != upper_[j as usize] && b_up_b[j as usize] > upper_[j as usize] {
                b_up_b[j as usize] = upper_[j as usize];
                debug_assert!(lower_[j as usize] > -kHighsInf);
                b_up_f[j as usize] =
                    objective + sense as f64 * (upper_[j as usize] - lower_[j as usize]) * dualv;
                b_up_e[j as usize] = j;
                b_up_l[j as usize] = j;
            }

            // Decrease x_j
            if ixj_dec[j as usize] != -1 {
                let i = ixj_dec[j as usize];
                b_dn_b[j as usize] = value_[j as usize] + txj_dec[j as usize];
                b_dn_f[j as usize] =
                    objective + sense as f64 * poss_inf_product(txj_dec[j as usize], dualv);
                b_dn_e[j as usize] = if wxj_dec[j as usize] > 0 {
                    jci_inc[i as usize]
                } else {
                    jci_dec[i as usize]
                };
                b_dn_l[j as usize] = bindex_[i as usize];
            } else {
                b_dn_b[j as usize] = -h_inf;
                b_dn_f[j as usize] = objective - sense as f64 * inf_product(dsign);
                b_dn_e[j as usize] = -1;
                b_dn_l[j as usize] = -1;
            }

            // Check if b_dn_b < lower
            if value_[j as usize] != lower_[j as usize] && b_dn_b[j as usize] < lower_[j as usize] {
                b_dn_b[j as usize] = lower_[j as usize];
                debug_assert!(upper_[j as usize] < kHighsInf);
                b_dn_f[j as usize] =
                    objective + sense as f64 * (lower_[j as usize] - upper_[j as usize]) * dualv;
                b_dn_e[j as usize] = j;
                b_dn_l[j as usize] = j;
            }
        }
    }

    //
    // Ranging 3.2. basic bounds ranging
    //
    for i in 0..num_row {
        for dir in [-1, 1] {
            let j = bindex_[i as usize];
            let newx: &mut f64 = if dir == -1 {
                &mut b_dn_b[j as usize]
            } else {
                &mut b_up_b[j as usize]
            };
            let newf: &mut f64 = if dir == -1 {
                &mut b_dn_f[j as usize]
            } else {
                &mut b_up_f[j as usize]
            };
            let j_enter: &mut HighsInt = if dir == -1 {
                &mut b_dn_e[j as usize]
            } else {
                &mut b_up_e[j as usize]
            };
            let j_leave: &mut HighsInt = if dir == -1 {
                &mut b_dn_l[j as usize]
            } else {
                &mut b_up_l[j as usize]
            };

            let j_in: HighsInt = if dir == -1 {
                jci_inc[i as usize]
            } else {
                jci_dec[i as usize]
            };
            let a_in: f64 = if dir == -1 {
                aci_inc[i as usize]
            } else {
                aci_dec[i as usize]
            };
            if j_in != -1 {
                let jmove = nmove_[j_in as usize];
                let i_out: HighsInt = if jmove > 0 {
                    ixj_inc[j_in as usize]
                } else {
                    ixj_dec[j_in as usize]
                };
                let j_out: HighsInt = if jmove > 0 {
                    jxj_inc[j_in as usize]
                } else {
                    jxj_dec[j_in as usize]
                };
                let w_out: HighsInt = if jmove > 0 {
                    wxj_inc[j_in as usize]
                } else {
                    wxj_dec[j_in as usize]
                };
                let tt: f64 = if jmove > 0 {
                    txj_inc[j_in as usize]
                } else {
                    txj_dec[j_in as usize]
                };
                if j_out == j_in {
                    // Bound flip
                    let delta = jmove as f64 * (upper_[j_in as usize] - lower_[j_in as usize]);
                    *newx = xi[i as usize] - delta * a_in;
                    *newf = objective + sense as f64 * delta * dual_[j_in as usize];
                    *j_enter = j_in;
                    *j_leave = j_out;
                } else if j_out != -1 {
                    // Regular
                    let delta: f64 = if w_out > 0 {
                        dxi_inc[i_out as usize]
                    } else {
                        dxi_dec[i_out as usize]
                    };
                    let a_out: f64 = if jmove > 0 {
                        axj_inc[j_in as usize]
                    } else {
                        axj_dec[j_in as usize]
                    };
                    *newx = xi[i as usize] + delta * a_in / a_out;
                    *newf = objective + sense as f64 * tt * dual_[j_in as usize];
                    *j_enter = j_in;
                    *j_leave = j_out;
                } else {
                    // Primal ratio test failed - change unlimitedly
                    //
                    // While still limited by its own bounds
                    //
                    // Its own bounds could just be inf
                    *newx = if dir == -1 {
                        lower_[j as usize]
                    } else {
                        upper_[j as usize]
                    };
                    *newf = objective;
                    *j_enter = -1;
                    *j_leave = -1;
                }
            } else {
                // Dual ratio test failed - just stay
                *newx = xi[i as usize];
                *newf = objective;
                *j_enter = -1;
                *j_leave = -1;
            }
        }
    }

    //
    // Ranging 4.1.1 Trim small value to zero
    //

    for j in 0..num_col {
        if c_up_c[j as usize].abs() < h_tt {
            c_up_c[j as usize] = 0.0;
        }
        if c_dn_c[j as usize].abs() < h_tt {
            c_dn_c[j as usize] = 0.0;
        }
        if b_up_b[j as usize].abs() < h_tt {
            b_up_b[j as usize] = 0.0;
        }
        if b_dn_b[j as usize].abs() < h_tt {
            b_dn_b[j as usize] = 0.0;
        }
    }
    {
        let mut i: HighsInt = 0;
        let mut j: HighsInt = num_col;
        while i < num_row {
            if b_up_b[j as usize].abs() < h_tt {
                b_up_b[j as usize] = 0.0;
            }
            if b_dn_b[j as usize].abs() < h_tt {
                b_dn_b[j as usize] = 0.0;
            }
            i += 1;
            j += 1;
        }
    }

    //
    // Ranging 4.2. Put to output buffer
    //

    if sense > 0 {
        ranging.col_cost_up.value_ = c_up_c;
        ranging.col_cost_dn.value_ = c_dn_c;
        ranging.col_cost_up.objective_ = c_up_f;
        ranging.col_cost_dn.objective_ = c_dn_f;
        ranging.col_cost_up.in_var_ = c_up_e;
        ranging.col_cost_dn.in_var_ = c_dn_e;
        ranging.col_cost_up.ou_var_ = c_up_l;
        ranging.col_cost_dn.ou_var_ = c_dn_l;
    } else {
        // For maximization problems, flip data and negate the cost values
        ranging.col_cost_up.value_ = c_dn_c;
        ranging.col_cost_dn.value_ = c_up_c;
        ranging.col_cost_up.objective_ = c_dn_f;
        ranging.col_cost_dn.objective_ = c_up_f;
        ranging.col_cost_up.in_var_ = c_dn_e;
        ranging.col_cost_dn.in_var_ = c_up_e;
        ranging.col_cost_up.ou_var_ = c_dn_l;
        ranging.col_cost_dn.ou_var_ = c_up_l;
        for v in ranging.col_cost_up.value_.iter_mut() {
            *v = -*v;
        }
        for v in ranging.col_cost_dn.value_.iter_mut() {
            *v = -*v;
        }
    }

    ranging.col_bound_up.value_ = b_up_b[..num_col as usize].to_vec();
    ranging.col_bound_dn.value_ = b_dn_b[..num_col as usize].to_vec();
    ranging.col_bound_up.objective_ = b_up_f[..num_col as usize].to_vec();
    ranging.col_bound_dn.objective_ = b_dn_f[..num_col as usize].to_vec();
    ranging.col_bound_up.in_var_ = b_up_e[..num_col as usize].to_vec();
    ranging.col_bound_dn.in_var_ = b_dn_e[..num_col as usize].to_vec();
    ranging.col_bound_up.ou_var_ = b_up_l[..num_col as usize].to_vec();
    ranging.col_bound_dn.ou_var_ = b_dn_l[..num_col as usize].to_vec();

    // Flip all data and negate the row bound values
    ranging.row_bound_up.value_ = b_dn_b[num_col as usize..num_total as usize].to_vec();
    for v in ranging.row_bound_up.value_.iter_mut() {
        *v = -*v;
    }

    ranging.row_bound_dn.value_ = b_up_b[num_col as usize..num_total as usize].to_vec();
    for v in ranging.row_bound_dn.value_.iter_mut() {
        *v = -*v;
    }

    ranging.row_bound_up.objective_ = b_dn_f[num_col as usize..num_total as usize].to_vec();
    ranging.row_bound_dn.objective_ = b_up_f[num_col as usize..num_total as usize].to_vec();
    ranging.row_bound_up.in_var_ = b_dn_e[num_col as usize..num_total as usize].to_vec();
    ranging.row_bound_dn.in_var_ = b_up_e[num_col as usize..num_total as usize].to_vec();
    ranging.row_bound_up.ou_var_ = b_dn_l[num_col as usize..num_total as usize].to_vec();
    ranging.row_bound_dn.ou_var_ = b_up_l[num_col as usize..num_total as usize].to_vec();

    ranging.valid = true;
    // writeRangingFile (HighsRanging.cpp 595-599) is file-writer code and is
    // not ported.
    HighsStatus::kOk
}
