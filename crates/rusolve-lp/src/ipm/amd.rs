//! Approximate minimum degree (AMD) fill-reducing ordering.
//!
//! Provenance: Amestoy, Davis & Duff, "An approximate minimum degree ordering
//! algorithm", SIAM J. Matrix Anal. Appl. 17(4):886–905 (1996), and Davis,
//! *Direct Methods for Sparse Linear Systems* (SIAM 2006), ch. 7. Implemented
//! from the paper's description (quotient graph, approximate external degree
//! `d_i = min(n−k, d_i_old + |Lp∖i|, |Ai∖i| + |Lp∖i| + Σ|Le∖Lp|)`, element
//! absorption, aggressive absorption, mass elimination, indistinguishable
//! supervariable detection by hashing, dense-row postponement at
//! `max(16, 10√n)`). No SuiteSparse code was read; the in-place `iw` arena
//! of the reference implementation is replaced by per-node `Vec` lists
//! (simpler, safe Rust; quotient-graph size never exceeds the input).
//!
//! Correctness contract: the output is ALWAYS a permutation of `0..n`
//! (fill quality is the only thing a bug here could hurt; the unit tests
//! check both permutation validity and fill against natural order).

const NONE: u32 = u32::MAX;

/// AMD ordering of the symmetric pattern given as full adjacency lists
/// (`adj[ptr[i]..ptr[i+1]]` = neighbours of `i`, both directions present, no
/// self loops, duplicates tolerated). Returns `perm` with `perm[k]` = the node
/// eliminated `k`-th.
pub fn amd_order(n: usize, ptr: &[usize], adj: &[u32]) -> Vec<usize> {
    if n == 0 {
        return vec![];
    }
    let dense_thresh = 16usize.max((10.0 * (n as f64).sqrt()) as usize);
    // Variable lists (adjacent variables), element lists (adjacent
    // elements) and element member lists.
    let mut vars: Vec<Vec<u32>> = Vec::with_capacity(n);
    let mut elems: Vec<Vec<u32>> = vec![Vec::new(); n];
    let mut lelem: Vec<Vec<u32>> = vec![Vec::new(); n];
    let mut nv: Vec<i64> = vec![1; n];
    let mut is_elem = vec![false; n];
    let mut dead_elem = vec![false; n];
    let mut degree: Vec<i64> = vec![0; n];
    let mut w: Vec<i64> = vec![0; n];
    let mut parent_of: Vec<u32> = vec![NONE; n];
    let mut pivot_step: Vec<u32> = vec![NONE; n];
    let mut mark: Vec<u32> = vec![0; n];
    let mut mark_stamp: u32 = 0;

    // Dedup the input lists (marker) and drop self loops.
    for i in 0..n {
        mark_stamp += 1;
        mark[i] = mark_stamp;
        let mut v = Vec::with_capacity(ptr[i + 1] - ptr[i]);
        for &j in &adj[ptr[i]..ptr[i + 1]] {
            if mark[j as usize] != mark_stamp {
                mark[j as usize] = mark_stamp;
                v.push(j);
            }
        }
        vars.push(v);
    }
    // Dense rows: postponed to the end, removed from the graph (nv = 0 makes
    // every scan skip them).
    let mut dense: Vec<usize> = Vec::new();
    for i in 0..n {
        if vars[i].len() > dense_thresh {
            dense.push(i);
            nv[i] = 0;
        }
    }
    for &d in &dense {
        vars[d] = Vec::new();
    }
    let n_eff = (n - dense.len()) as i64;

    // Degree lists.
    let mut head: Vec<u32> = vec![NONE; n + 1];
    let mut next: Vec<u32> = vec![NONE; n];
    let mut prev: Vec<u32> = vec![NONE; n];
    let insert = |head: &mut [u32], next: &mut [u32], prev: &mut [u32], i: usize, d: usize| {
        let h = head[d];
        next[i] = h;
        prev[i] = NONE;
        if h != NONE {
            prev[h as usize] = i as u32;
        }
        head[d] = i as u32;
    };
    let remove = |head: &mut [u32], next: &mut [u32], prev: &mut [u32], i: usize, d: usize| {
        let (nx, pv) = (next[i], prev[i]);
        if nx != NONE {
            prev[nx as usize] = pv;
        }
        if pv != NONE {
            next[pv as usize] = nx;
        } else {
            head[d] = nx;
        }
    };
    for i in 0..n {
        if nv[i] == 0 {
            continue;
        }
        let d = vars[i].iter().filter(|&&j| nv[j as usize] > 0).count();
        degree[i] = d as i64;
        insert(&mut head, &mut next, &mut prev, i, d);
    }

    let mut mindeg = 0usize;
    let mut nel: i64 = 0;
    let mut step: u32 = 0;
    let mut wflg: i64 = 2;
    let mut lme: Vec<u32> = Vec::new();
    let mut hashes: Vec<(u64, u32)> = Vec::new();

    while nel < n_eff {
        // 1. Pivot: minimum approximate degree.
        while head[mindeg] == NONE {
            mindeg += 1;
        }
        let me = head[mindeg] as usize;
        remove(&mut head, &mut next, &mut prev, me, mindeg);
        pivot_step[me] = step;
        step += 1;
        let mut nvpiv = nv[me];
        nel += nvpiv;
        nv[me] = -nvpiv;

        // 2. New element Lme = (A_me ∪ ⋃ Le) ∖ me; absorb the Le.
        lme.clear();
        let mut degme: i64 = 0;
        let me_vars = std::mem::take(&mut vars[me]);
        let me_elems = std::mem::take(&mut elems[me]);
        for &i in &me_vars {
            let iu = i as usize;
            if nv[iu] > 0 && !is_elem[iu] {
                degme += nv[iu];
                remove(&mut head, &mut next, &mut prev, iu, degree[iu] as usize);
                nv[iu] = -nv[iu];
                lme.push(i);
            }
        }
        for &e in &me_elems {
            let eu = e as usize;
            if dead_elem[eu] {
                continue;
            }
            let members = std::mem::take(&mut lelem[eu]);
            for &i in &members {
                let iu = i as usize;
                if nv[iu] > 0 && !is_elem[iu] {
                    degme += nv[iu];
                    remove(&mut head, &mut next, &mut prev, iu, degree[iu] as usize);
                    nv[iu] = -nv[iu];
                    lme.push(i);
                }
            }
            dead_elem[eu] = true;
            w[eu] = 0;
        }
        is_elem[me] = true;

        // 3. |Le ∖ Lme| for every element adjacent to Lme: w[e] − wflg.
        for &i in &lme {
            let iu = i as usize;
            let nvi = -nv[iu];
            for &e in &elems[iu] {
                let eu = e as usize;
                if dead_elem[eu] {
                    continue;
                }
                let we = w[eu];
                w[eu] = if we >= wflg {
                    we - nvi
                } else {
                    degree[eu] + wflg - nvi
                };
            }
        }

        // 4. Degree update, pruning, aggressive absorption, mass elimination.
        hashes.clear();
        for idx in 0..lme.len() {
            let i = lme[idx] as usize;
            let nvi = -nv[i];
            let mut deg: i64 = 0;
            let mut el = std::mem::take(&mut elems[i]);
            let mut h: u64 = me as u64;
            el.retain(|&e| {
                let eu = e as usize;
                if dead_elem[eu] {
                    return false;
                }
                let dext = w[eu] - wflg;
                if dext > 0 {
                    deg += dext;
                    h = h.wrapping_add(e as u64);
                    true
                } else {
                    // Aggressive absorption: Le ⊆ Lme.
                    dead_elem[eu] = true;
                    w[eu] = 0;
                    lelem[eu] = Vec::new();
                    false
                }
            });
            let mut vl = std::mem::take(&mut vars[i]);
            vl.retain(|&j| {
                let ju = j as usize;
                if nv[ju] > 0 && !is_elem[ju] {
                    deg += nv[ju];
                    h = h.wrapping_add(j as u64);
                    true
                } else {
                    false
                }
            });
            if el.is_empty() && vl.is_empty() {
                // Mass elimination: i is adjacent to me only.
                parent_of[i] = me as u32;
                nv[i] = 0;
                degme -= nvi;
                nvpiv += nvi;
                nel += nvi;
            } else {
                degree[i] = degree[i].min(deg);
                el.insert(0, me as u32);
                elems[i] = el;
                vars[i] = vl;
                hashes.push((h, i as u32));
            }
        }
        degree[me] = degme;

        // 5. Supervariable detection among Lme (identical pruned lists).
        hashes.sort_unstable();
        let mut g = 0;
        while g < hashes.len() {
            let mut g_end = g + 1;
            while g_end < hashes.len() && hashes[g_end].0 == hashes[g].0 {
                g_end += 1;
            }
            if g_end - g > 1 {
                for a in g..g_end {
                    let i = hashes[a].1 as usize;
                    if nv[i] == 0 {
                        continue;
                    }
                    mark_stamp = mark_stamp.wrapping_add(1);
                    if mark_stamp == 0 {
                        mark.iter_mut().for_each(|m| *m = 0);
                        mark_stamp = 1;
                    }
                    for &e in &elems[i] {
                        mark[e as usize] = mark_stamp;
                    }
                    for &v in &vars[i] {
                        mark[v as usize] = mark_stamp;
                    }
                    for b in a + 1..g_end {
                        let j = hashes[b].1 as usize;
                        if nv[j] == 0
                            || elems[j].len() != elems[i].len()
                            || vars[j].len() != vars[i].len()
                        {
                            continue;
                        }
                        let same = elems[j].iter().all(|&e| mark[e as usize] == mark_stamp)
                            && vars[j].iter().all(|&v| mark[v as usize] == mark_stamp);
                        if same {
                            // Both are flagged negative inside Lme.
                            nv[i] += nv[j];
                            nv[j] = 0;
                            parent_of[j] = i as u32;
                            elems[j] = Vec::new();
                            vars[j] = Vec::new();
                        }
                    }
                }
            }
            g = g_end;
        }

        // 6. Final approximate degrees; restore flags; reinsert.
        let nleft = n_eff - nel;
        let mut kept: Vec<u32> = Vec::with_capacity(lme.len());
        for &i in &lme {
            let iu = i as usize;
            if nv[iu] == 0 {
                continue;
            }
            let nvi = -nv[iu];
            let d = (degree[iu] + degme - nvi).min(nleft - nvi).max(0);
            degree[iu] = d;
            nv[iu] = nvi;
            insert(&mut head, &mut next, &mut prev, iu, d as usize);
            if (d as usize) < mindeg {
                mindeg = d as usize;
            }
            kept.push(i);
        }
        lelem[me] = kept;
        nv[me] = 0;
        let _ = nvpiv;
        // Every w[e] set this step is < wflg + n: move past it.
        wflg += n as i64 + 1;
    }

    // Ordering: pivots in selection order, each followed by the variables it
    // absorbed (supervariable members / mass-eliminated), dense rows last.
    let mut root_step: Vec<u32> = vec![NONE; n];
    for v in 0..n {
        if pivot_step[v] != NONE {
            root_step[v] = pivot_step[v];
        }
    }
    let mut stack: Vec<usize> = Vec::new();
    for v in 0..n {
        if root_step[v] != NONE || parent_of[v] == NONE {
            continue;
        }
        let mut u = v;
        while root_step[u] == NONE {
            stack.push(u);
            u = parent_of[u] as usize;
        }
        let r = root_step[u];
        for &s in &stack {
            root_step[s] = r;
        }
        stack.clear();
    }
    let mut keyed: Vec<(u32, u8, u32)> = Vec::with_capacity(n);
    for v in 0..n {
        if root_step[v] != NONE {
            let follower = u8::from(pivot_step[v] == NONE);
            keyed.push((root_step[v], follower, v as u32));
        }
    }
    keyed.sort_unstable();
    let mut perm: Vec<usize> = keyed.iter().map(|&(_, _, v)| v as usize).collect();
    perm.extend(dense.iter().copied());
    // Defensive: the output must be a permutation whatever happened above
    // (duplicates dropped, missing nodes appended last).
    if perm.len() != n {
        let mut seen = vec![false; n];
        perm.retain(|&p| !std::mem::replace(&mut seen[p], true));
        perm.extend((0..n).filter(|&v| !seen[v]));
    }
    perm
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::XorShift64;

    fn sym_lists(n: usize, edges: &[(usize, usize)]) -> (Vec<usize>, Vec<u32>) {
        let mut l: Vec<Vec<u32>> = vec![Vec::new(); n];
        for &(a, b) in edges {
            if a != b {
                l[a].push(b as u32);
                l[b].push(a as u32);
            }
        }
        let mut ptr = vec![0];
        let mut adj = Vec::new();
        for v in l {
            adj.extend(v);
            ptr.push(adj.len());
        }
        (ptr, adj)
    }

    /// Symbolic fill of elimination order `perm` (dense, small n).
    fn fill(n: usize, edges: &[(usize, usize)], perm: &[usize]) -> usize {
        let mut pinv = vec![0; n];
        for (k, &p) in perm.iter().enumerate() {
            pinv[p] = k;
        }
        let mut g = vec![vec![false; n]; n];
        for &(a, b) in edges {
            let (a, b) = (pinv[a], pinv[b]);
            g[a][b] = true;
            g[b][a] = true;
        }
        let mut nnz = 0;
        for k in 0..n {
            let nb: Vec<usize> = (k + 1..n).filter(|&j| g[k][j]).collect();
            nnz += nb.len();
            for &a in &nb {
                for &b in &nb {
                    if a != b {
                        g[a][b] = true;
                    }
                }
            }
        }
        nnz
    }

    fn is_perm(p: &[usize], n: usize) -> bool {
        let mut s = vec![false; n];
        p.len() == n
            && p.iter()
                .all(|&v| v < n && !std::mem::replace(&mut s[v], true))
    }

    #[test]
    fn arrow_matrix_orders_hub_last() {
        // Star: hub 0 connected to all. Natural order fills completely;
        // minimum degree eliminates the leaves first (zero fill).
        let n = 12;
        let edges: Vec<(usize, usize)> = (1..n).map(|j| (0, j)).collect();
        let (ptr, adj) = sym_lists(n, &edges);
        let p = amd_order(n, &ptr, &adj);
        assert!(is_perm(&p, n));
        assert_eq!(fill(n, &edges, &p), n - 1);
        let natural: Vec<usize> = (0..n).collect();
        assert!(fill(n, &edges, &natural) > fill(n, &edges, &p));
    }

    #[test]
    fn grid_fill_beats_natural_and_is_permutation() {
        let k = 15;
        let n = k * k;
        let mut edges = Vec::new();
        for r in 0..k {
            for c in 0..k {
                let v = r * k + c;
                if c + 1 < k {
                    edges.push((v, v + 1));
                }
                if r + 1 < k {
                    edges.push((v, v + k));
                }
            }
        }
        let (ptr, adj) = sym_lists(n, &edges);
        let p = amd_order(n, &ptr, &adj);
        assert!(is_perm(&p, n));
        let natural: Vec<usize> = (0..n).collect();
        let (fa, fn_) = (fill(n, &edges, &p), fill(n, &edges, &natural));
        assert!(fa < fn_, "amd fill {fa} vs natural {fn_}");
    }

    #[test]
    fn random_graphs_always_permutation() {
        let mut rng = XorShift64::new(7);
        for trial in 0..40 {
            let n = 1 + (trial * 7) % 90;
            let m = n * 3;
            let mut edges = Vec::new();
            for _ in 0..m {
                let a = (rng.uniform(0.0, n as f64) as usize).min(n - 1);
                let b = (rng.uniform(0.0, n as f64) as usize).min(n - 1);
                edges.push((a, b));
            }
            // a few dense hubs
            if n > 40 {
                for j in 0..n {
                    edges.push((0, j));
                }
            }
            let (ptr, adj) = sym_lists(n, &edges);
            let p = amd_order(n, &ptr, &adj);
            assert!(is_perm(&p, n), "trial {trial}");
            let natural: Vec<usize> = (0..n).collect();
            // Not a strict guarantee, but AMD should never be wildly worse.
            assert!(fill(n, &edges, &p) <= fill(n, &edges, &natural) * 2 + n);
        }
    }
}
