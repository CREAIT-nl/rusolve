//! rusolve debug infrastructure for mip_symmetry_dynamic.
//!
//! No SCIP counterpart exists. Semantics model: HiGHS HighsDebugSol.cpp
//! 27-30 (activate gate), 70-74 (objective), 125 (deactivate), 167-177
//! (checkCut violation). See THIRD_PARTY_NOTICES.md.

use super::*;
use crate::prelude::*;
use std::collections::{HashMap, HashSet};
use std::sync::OnceLock;
use std::sync::atomic::{AtomicBool, Ordering};

pub(crate) const S3_DEBUG_SOL_ENV: &str = "RUSOLVE_S3_DEBUG_SOL";
pub(crate) const S3_DEBUG_ALIVE_CAP: usize = 100_000;
pub(crate) const S3_DEBUG_TOL: f64 = 1e-6;
/// Memory guard of the LexRed BFS (debug only): stored image values
/// (images times support size) above this report "unchecked".
pub(crate) const S3_DEBUG_ALIVE_MAX_VALUES: usize = 50_000_000;

static DEBUG_SOL: OnceLock<Option<Vec<f64>>> = OnceLock::new();
static ANNOUNCED: AtomicBool = AtomicBool::new(false);

/// Read the debug solution from the environment variable once per process.
/// Returns None when the variable is unset, the file cannot be read, or no
/// values parse. Never panics.
pub(crate) fn debugsol_values() -> Option<&'static [f64]> {
    DEBUG_SOL
        .get_or_init(|| {
            let path = std::env::var(S3_DEBUG_SOL_ENV).ok()?;
            let content = std::fs::read_to_string(&path).ok()?;
            let vals: Vec<f64> = content
                .split_whitespace()
                .filter_map(|s| s.parse::<f64>().ok())
                .collect();
            if vals.is_empty() { None } else { Some(vals) }
        })
        .as_deref()
}

/// True only when all activation conditions hold.
pub(crate) fn debugsol_active(ms: &HighsMipSolver, x: &[f64]) -> bool {
    if ms.submip {
        return false;
    }
    if ms.options_mip_.presolve != kHighsOffString {
        return false;
    }
    if ms.mipdata_.num_restarts != 0 {
        return false;
    }
    if ms.has_multiple_workers() {
        return false;
    }
    if x.len() != ms.orig_model_.num_col_ as usize {
        return false;
    }
    // presolve off: the search domain has the original columns (indexing guard)
    if x.len() != mip_model!(ms).num_col_ as usize {
        return false;
    }
    let obj = debugsol_objective(ms, x);
    obj <= ms.mipdata_.upper_limit
}

/// HighsCDouble objective of x on orig_model_ (HighsDebugSol.cpp:70-74).
pub(crate) fn debugsol_objective(ms: &HighsMipSolver, x: &[f64]) -> f64 {
    let mut debugsolobj = HighsCDouble::from(0.0);
    for i in 0..ms.orig_model_.num_col_ as usize {
        // C++ `col_cost_[i] * HighsCDouble(x[i])` is `HighsCDouble(x[i]) * col_cost_[i]`.
        debugsolobj += HighsCDouble::from(x[i]) * ms.orig_model_.col_cost_[i];
    }
    f64::from(debugsolobj + ms.orig_model_.offset_)
}

/// Bipartite matching (Kuhn augmenting paths) for a single orbitope.
/// Alive iff a perfect matching exists from x's column values to domain-allowed
/// columns.
pub(crate) fn alive_orbitope(
    orb: &S3OrbitopeData,
    x: &[f64],
    lb: &[f64],
    ub: &[f64],
    tol: f64,
) -> bool {
    let ncols = orb.ncols as usize;
    let nrows = orb.nrows as usize;
    if ncols == 0 {
        return true;
    }

    // match_col[k] = which source column j is matched to target column k
    let mut match_col: Vec<usize> = vec![usize::MAX; ncols];

    // Build adjacency: for each source column j, which target columns k are allowed?
    // Edge (j,k) iff for every row i: lb[vars[i][k]] - tol <= x[vars[i][j]] <= ub[vars[i][k]] + tol
    let adj: Vec<Vec<usize>> = (0..ncols)
        .map(|j| {
            (0..ncols)
                .filter(|&k| {
                    (0..nrows).all(|i| {
                        let src_col = orb.vars[i][j] as usize;
                        let tgt_col = orb.vars[i][k] as usize;
                        let val = x[src_col];
                        lb[tgt_col] - tol <= val && val <= ub[tgt_col] + tol
                    })
                })
                .collect()
        })
        .collect();

    // Kuhn's algorithm
    for j in 0..ncols {
        let mut visited = vec![false; ncols];
        if !kuhn_dfs(j, &adj, &mut match_col, &mut visited) {
            return false;
        }
    }
    true
}

fn kuhn_dfs(j: usize, adj: &[Vec<usize>], match_col: &mut [usize], visited: &mut [bool]) -> bool {
    for &k in &adj[j] {
        if visited[k] {
            continue;
        }
        visited[k] = true;
        if match_col[k] == usize::MAX || kuhn_dfs(match_col[k], adj, match_col, visited) {
            match_col[k] = j;
            return true;
        }
    }
    false
}

/// BFS over images of x under all lexdata generators.
/// Returns Some(true) if an image lies in the domain, Some(false) if BFS
/// exhausts without finding one, None if the cap is exceeded or lexdatas is empty.
pub(crate) fn alive_lexdatas(
    lexdas: &[S3LexData],
    x: &[f64],
    lb: &[f64],
    ub: &[f64],
    tol: f64,
    cap: usize,
) -> Option<bool> {
    if lexdas.is_empty() {
        return None;
    }

    // The union of all columns touched by any lexdata; images are stored on these
    // columns only (the permutations fix every other column).
    let mut all_cols: Vec<HighsInt> = Vec::new();
    let mut pos: HashMap<HighsInt, usize> = HashMap::new();
    for ld in lexdas {
        for &col in &ld.vars {
            if let std::collections::hash_map::Entry::Vacant(e) = pos.entry(col) {
                e.insert(all_cols.len());
                all_cols.push(col);
            }
        }
    }
    if all_cols.is_empty() {
        return None;
    }
    let m = all_cols.len();

    // (source, destination) positions of each generator
    let gens: Vec<Vec<(usize, usize)>> = lexdas
        .iter()
        .map(|ld| {
            (0..ld.vars.len())
                .map(|i| (pos[&ld.vars[i]], pos[&ld.vars[ld.perm[i] as usize]]))
                .collect()
        })
        .collect();

    let in_domain = |img: &[f64]| -> bool {
        all_cols.iter().zip(img).all(|(&col, &v)| {
            let ci = col as usize;
            lb[ci] - tol <= v && v <= ub[ci] + tol
        })
    };

    let start: Vec<f64> = all_cols.iter().map(|&col| x[col as usize]).collect();
    if in_domain(&start) {
        return Some(true);
    }

    // BFS
    let key = |v: &[f64]| -> Vec<u64> { v.iter().map(|f| f.to_bits()).collect() };
    let mut visited: HashSet<Vec<u64>> = HashSet::new();
    let mut queue: Vec<Vec<f64>> = Vec::new();
    visited.insert(key(&start));
    queue.push(start);

    let mut head = 0;
    while head < queue.len() {
        if visited.len() > cap || visited.len().saturating_mul(m) > S3_DEBUG_ALIVE_MAX_VALUES {
            return None;
        }
        let cur = std::mem::take(&mut queue[head]);
        head += 1;

        for g in &gens {
            if g.is_empty() {
                continue;
            }
            let mut img = cur.clone();
            for &(src, dst) in g {
                img[dst] = cur[src];
            }
            if !visited.insert(key(&img)) {
                continue;
            }
            if in_domain(&img) {
                return Some(true);
            }
            queue.push(img);
        }
    }

    Some(false)
}

/// True iff no column appears both in an orbitope (dynamic or static) and in a
/// lexdata.
pub(crate) fn sym_components_disjoint(data: &S3SymData, static_orbs: &[S3OrbitopeData]) -> bool {
    let mut orb_cols: HashSet<HighsInt> = HashSet::new();
    for orb in data.orbitopes.iter().chain(static_orbs) {
        for row in &orb.vars {
            for &col in row {
                orb_cols.insert(col);
            }
        }
    }
    for ld in &data.lexdatas {
        for &col in &ld.vars {
            if orb_cols.contains(&col) {
                return false;
            }
        }
    }
    true
}

/// The HiGHS static orbitopes kept by D6 (`data.static_orbitopes`, reduced by
/// `orbital_fixing` in D8 step 1) as column matrices `vars[row][col]`, so that the
/// alive check also accounts for their column symmetry.
pub(crate) fn debug_static_orbitopes(ms: &HighsMipSolver, data: &S3SymData) -> Vec<S3OrbitopeData> {
    data.static_orbitopes
        .iter()
        .map(|&idx| {
            let o = &ms.mipdata_.symmetries.orbitopes[idx];
            let nrows = o.num_rows;
            let ncols = o.row_length;
            let vars: Vec<Vec<HighsInt>> = (0..nrows)
                .map(|i| {
                    (0..ncols)
                        .map(|j| o.matrix[(i + j * nrows) as usize])
                        .collect()
                })
                .collect();
            S3OrbitopeData {
                vars,
                nrows,
                ncols,
                ..S3OrbitopeData::default()
            }
        })
        .collect()
}

/// Tree-global alive check: some symmetric image of x* lies in [lb, ub].
/// `static_orbs` are the HiGHS static orbitopes ([`debug_static_orbitopes`]).
/// Returns None if components are not disjoint or any component is unchecked.
pub(crate) fn dom_alive(
    data: &S3SymData,
    static_orbs: &[S3OrbitopeData],
    x: &[f64],
    lb: &[f64],
    ub: &[f64],
    tol: f64,
) -> Option<bool> {
    if !sym_components_disjoint(data, static_orbs) {
        return None;
    }

    // Collect all symmetric columns
    let mut sym_cols: HashSet<HighsInt> = HashSet::new();
    for orb in data.orbitopes.iter().chain(static_orbs) {
        for row in &orb.vars {
            for &col in row {
                sym_cols.insert(col);
            }
        }
    }
    for ld in &data.lexdatas {
        for &col in &ld.vars {
            sym_cols.insert(col);
        }
    }

    // Non-symmetric columns must satisfy bounds
    for col in 0..x.len() {
        if !sym_cols.contains(&(col as HighsInt)) {
            if !(lb[col] - tol <= x[col] && x[col] <= ub[col] + tol) {
                return Some(false);
            }
        }
    }

    // Every orbitope must be alive
    for orb in data.orbitopes.iter().chain(static_orbs) {
        if !alive_orbitope(orb, x, lb, ub, tol) {
            return Some(false);
        }
    }

    // Lexdatas via BFS
    if !data.lexdatas.is_empty() {
        match alive_lexdatas(&data.lexdatas, x, lb, ub, tol, S3_DEBUG_ALIVE_CAP) {
            Some(true) => {}
            Some(false) => return Some(false),
            None => return None,
        }
    }

    Some(true)
}

/// Scan a list of open-node domains. Returns (checked, any_alive, any_unchecked).
pub(crate) fn scan_open_domains(
    data: &S3SymData,
    static_orbs: &[S3OrbitopeData],
    x: &[f64],
    domains: &[(Vec<f64>, Vec<f64>)],
    tol: f64,
) -> (usize, bool, bool) {
    let checked = domains.len();
    let mut any_alive = false;
    let mut any_unchecked = false;
    for (lb, ub) in domains {
        match dom_alive(data, static_orbs, x, lb, ub, tol) {
            Some(true) => any_alive = true,
            Some(false) => {}
            None => any_unchecked = true,
        }
    }
    (checked, any_alive, any_unchecked)
}

/// True iff no open node is alive, none unchecked, and incumbent exceeds debug obj.
pub(crate) fn violation(
    any_alive: bool,
    any_unchecked: bool,
    incumbent: f64,
    debug_obj: f64,
) -> bool {
    !any_alive && !any_unchecked && incumbent > debug_obj + S3_DEBUG_TOL
}

/// Applies `chgs` to `lb`/`ub` as tightenings.
fn apply_domain_changes(lb: &mut [f64], ub: &mut [f64], chgs: &[HighsDomainChange]) {
    for dc in chgs {
        let col = dc.column as usize;
        match dc.boundtype {
            HighsBoundType::kLower => {
                lb[col] = cpp_max(lb[col], dc.boundval);
            }
            HighsBoundType::kUpper => {
                ub[col] = cpp_min(ub[col], dc.boundval);
            }
        }
    }
}

/// Scan the open nodes for alive images of x*. Prints checked count and VIOLATION
/// if appropriate.
///
/// Open nodes: every node of the global queue (its domchgstack over the global
/// domain), and the sibling of every branching on `localdom`'s stack (the stack
/// prefix before the branching plus `localdom.flip` of it, as HighsSearch::backtrack
/// builds it). HighsSearch.nodestack is not reachable from s3_propagate, so siblings
/// that were already explored are scanned too: this can only suppress a VIOLATION,
/// never report a false one.
pub(crate) fn kill_scan(
    ms: &HighsMipSolver,
    data: &S3SymData,
    static_orbs: &[S3OrbitopeData],
    x: &[f64],
    localdom: &HighsDomain,
) {
    let feastol = ms.options_mip_.mip_feasibility_tolerance;
    let global_lb = &ms.mipdata_.domains[0].col_lower_;
    let global_ub = &ms.mipdata_.domains[0].col_upper_;

    let nodequeue = &*ms.mipdata_.nodequeue;
    let free_set: HashSet<i64> = nodequeue.freeslots.iter().map(|r| r.0).collect();

    let mut domains: Vec<(Vec<f64>, Vec<f64>)> = Vec::new();

    for (idx, node) in nodequeue.nodes.iter().enumerate() {
        if free_set.contains(&(idx as i64)) {
            continue;
        }
        let mut lb = global_lb.clone();
        let mut ub = global_ub.clone();
        apply_domain_changes(&mut lb, &mut ub, &node.domchgstack);
        domains.push((lb, ub));
    }

    let stack = localdom.get_domain_change_stack();
    for &pos in localdom.get_branching_positions() {
        let pos = pos as usize;
        let mut lb = global_lb.clone();
        let mut ub = global_ub.clone();
        apply_domain_changes(&mut lb, &mut ub, &stack[..pos]);
        let sibling = localdom.flip(ms, &stack[pos]);
        apply_domain_changes(&mut lb, &mut ub, std::slice::from_ref(&sibling));
        domains.push((lb, ub));
    }

    let (checked, any_alive, any_unchecked) =
        scan_open_domains(data, static_orbs, x, &domains, feastol);

    eprintln!("s3-debugsol checked n={}", checked);

    let incumbent = ms.mipdata_.upper_bound;
    let debug_obj = debugsol_objective(ms, x);

    if violation(any_alive, any_unchecked, incumbent, debug_obj) {
        eprintln!(
            "s3-debugsol VIOLATION obj={} incumbent={}",
            debug_obj, incumbent
        );
    }
}

/// Print the activation line once per process.
pub(crate) fn announce_once(obj: f64) {
    if !ANNOUNCED.swap(true, Ordering::Relaxed) {
        eprintln!("s3-debugsol active obj={}", obj);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mip::symmetry_dynamic::lexred::lexicographic_reduction_add_permutation;

    #[test]
    fn s3_dbg_alive_orbitope_matching() {
        // Orbitope 2x2: vars = [[0,1],[2,3]], nrows=2, ncols=2
        let orb = S3OrbitopeData {
            vars: vec![vec![0, 1], vec![2, 3]],
            nrows: 2,
            ncols: 2,
            nbranchrows: 0,
            rowindexmap: vec![0, 0, 1, 1],
            colindexmap: vec![0, 1, 0, 1],
            columnordering: 0,
            rowordering: 1,
        };
        let x = [1.0, 0.0, 0.0, 1.0];
        let tol = 1e-9;

        // Domain A: lb = ub = [0,1,1,0] -> column swap maps x into domain
        let lb_a = [0.0, 1.0, 1.0, 0.0];
        let ub_a = [0.0, 1.0, 1.0, 0.0];
        assert!(alive_orbitope(&orb, &x, &lb_a, &ub_a, tol));

        // Domain B: lb = ub = [1,1,1,1] -> no image of x fits
        let lb_b = [1.0, 1.0, 1.0, 1.0];
        let ub_b = [1.0, 1.0, 1.0, 1.0];
        assert!(!alive_orbitope(&orb, &x, &lb_b, &ub_b, tol));
    }

    #[test]
    fn s3_dbg_alive_lexred_swap() {
        let mut lexdatas: Vec<S3LexData> = Vec::new();
        let mut symvarmap: std::collections::HashMap<HighsInt, HighsInt> =
            std::collections::HashMap::new();
        let mut nsymvars: HighsInt = 0;
        lexicographic_reduction_add_permutation(
            &mut lexdatas,
            &mut symvarmap,
            &mut nsymvars,
            &[0, 1],
            &[1, 0],
            false,
        );

        let x = [1.0, 0.0];
        let tol = 1e-9;

        // Domain lb=[0,1], ub=[0,1]: x=[1,0] is not in domain, but image [0,1] is
        let lb1 = [0.0, 1.0];
        let ub1 = [0.0, 1.0];
        let result = alive_lexdatas(&lexdatas, &x, &lb1, &ub1, tol, S3_DEBUG_ALIVE_CAP);
        assert_eq!(result, Some(true));

        // Domain lb=[1,1], ub=[1,1]: neither x=[1,0] nor image [0,1] is in domain
        let lb2 = [1.0, 1.0];
        let ub2 = [1.0, 1.0];
        let result = alive_lexdatas(&lexdatas, &x, &lb2, &ub2, tol, S3_DEBUG_ALIVE_CAP);
        assert_eq!(result, Some(false));
    }

    #[test]
    fn s3_dbg_violation_needs_dead_open_nodes() {
        let mut lexdatas: Vec<S3LexData> = Vec::new();
        let mut symvarmap: std::collections::HashMap<HighsInt, HighsInt> =
            std::collections::HashMap::new();
        let mut nsymvars: HighsInt = 0;
        lexicographic_reduction_add_permutation(
            &mut lexdatas,
            &mut symvarmap,
            &mut nsymvars,
            &[0, 1],
            &[1, 0],
            false,
        );

        let data = S3SymData {
            orbitopes: Vec::new(),
            static_orbitopes: Vec::new(),
            lexdatas,
        };

        let x = [1.0, 0.0];
        let tol = 1e-9;

        // Case 1: one open domain with the image in it -> alive, no violation
        let domains_alive = vec![(vec![0.0, 1.0], vec![0.0, 1.0])];
        let (checked, any_alive, any_unchecked) =
            scan_open_domains(&data, &[], &x, &domains_alive, tol);
        assert_eq!(checked, 1);
        assert!(any_alive);
        assert!(!any_unchecked);
        assert!(!violation(any_alive, any_unchecked, 10.0, 0.0));

        // Case 2: only open domain lb=ub=[1,1] -> dead, violation
        let domains_dead = vec![(vec![1.0, 1.0], vec![1.0, 1.0])];
        let (checked, any_alive, any_unchecked) =
            scan_open_domains(&data, &[], &x, &domains_dead, tol);
        assert_eq!(checked, 1);
        assert!(!any_alive);
        assert!(!any_unchecked);
        assert!(violation(false, false, 10.0, 0.0));

        // Unchecked suppresses violation
        assert!(!violation(false, true, 10.0, 0.0));
    }

    #[test]
    fn s3_dbg_static_orbitope_counts_as_symmetry() {
        // 1-row static orbitope over columns 0,1 (D6 keeps it for orbital_fixing)
        let stat = S3OrbitopeData {
            vars: vec![vec![0, 1]],
            nrows: 1,
            ncols: 2,
            ..S3OrbitopeData::default()
        };
        let data = S3SymData::default();
        let x = [1.0, 0.0, 1.0];
        let tol = 1e-9;
        // x is outside, its column swap [0,1,1] is inside
        let lb = [0.0, 1.0, 1.0];
        let ub = [0.0, 1.0, 1.0];
        assert_eq!(
            dom_alive(&data, std::slice::from_ref(&stat), &x, &lb, &ub, tol),
            Some(true)
        );
        // without the static orbitope, columns 0,1 are non-symmetric
        assert_eq!(dom_alive(&data, &[], &x, &lb, &ub, tol), Some(false));
        // non-symmetric column 2 outside its bounds
        let ub2 = [0.0, 1.0, 0.0];
        assert_eq!(dom_alive(&data, &[stat], &x, &lb, &ub2, tol), Some(false));
    }
}
