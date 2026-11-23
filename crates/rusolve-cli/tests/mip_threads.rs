//! `rusolve mip --threads N` runs the 1.16 multi-worker tree search:
//! every thread count reaches the brute-force optimum of a small 3-constraint knapsack that
//! needs branching, and a repeated run at the same thread count and seed repeats the node
//! count (HiGHS 1.16's parallel search is deterministic by design).

use std::process::Command;

const N: usize = 22;
const ROWS: usize = 3;

/// Deterministic correlated data: values close to the mean weight make the LP bound loose.
fn data() -> (Vec<i64>, Vec<Vec<i64>>, Vec<i64>) {
    let mut state = 12345u64;
    let mut next = |m: i64| {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((state >> 33) as i64) % m
    };
    let w: Vec<Vec<i64>> = (0..ROWS)
        .map(|_| (0..N).map(|_| 10 + next(40)).collect())
        .collect();
    let v: Vec<i64> = (0..N)
        .map(|j| (w[0][j] + w[1][j] + w[2][j]) / 3 + next(5))
        .collect();
    let cap: Vec<i64> = w.iter().map(|r| r.iter().sum::<i64>() * 2 / 5).collect();
    (v, w, cap)
}

fn lp_text() -> String {
    let (v, w, cap) = data();
    let obj: Vec<String> = v
        .iter()
        .enumerate()
        .map(|(j, c)| format!("{c} x{j}"))
        .collect();
    let mut s = format!("Maximize\n obj: {}\nSubject To\n", obj.join(" + "));
    for (r, row) in w.iter().enumerate() {
        let t: Vec<String> = row
            .iter()
            .enumerate()
            .map(|(j, c)| format!("{c} x{j}"))
            .collect();
        s += &format!(" c{r}: {} <= {}\n", t.join(" + "), cap[r]);
    }
    s += "Binary\n";
    for j in 0..N {
        s += &format!(" x{j}");
    }
    s + "\nEnd\n"
}

fn brute_force() -> i64 {
    let (v, w, cap) = data();
    let mut best = 0;
    for mask in 0u32..(1 << N) {
        let mut load = [0i64; ROWS];
        let mut val = 0;
        for j in 0..N {
            if mask >> j & 1 == 1 {
                val += v[j];
                for r in 0..ROWS {
                    load[r] += w[r][j];
                }
            }
        }
        if (0..ROWS).all(|r| load[r] <= cap[r]) && val > best {
            best = val;
        }
    }
    best
}

fn solve(path: &std::path::Path, threads: u32) -> serde_json::Value {
    let out = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
        .env("RUSOLVE_MIP_HIGHS", "1")
        .arg("mip")
        .arg(path)
        .args(["--threads", &threads.to_string(), "--format", "json"])
        .args(["--opt", "presolve=off", "--opt", "random_seed=3"])
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    serde_json::from_slice(&out.stdout).unwrap()
}

#[test]
fn threads_reach_the_optimum_and_repeat() {
    let dir = std::env::temp_dir().join(format!("rusolve-mip-threads-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("mk.lp");
    std::fs::write(&path, lp_text()).unwrap();
    let expected = brute_force() as f64;
    for threads in [1, 2, 4] {
        let r = solve(&path, threads);
        assert_eq!(r["status"], "Optimal", "threads {threads}: {r}");
        let obj = r["objective"].as_f64().unwrap();
        assert!(
            (obj - expected).abs() < 1e-6,
            "threads {threads}: {obj} != {expected}"
        );
        assert_eq!(r["run"]["threads"], threads, "{r}");
        assert!(
            r["stats"]["nodes"].as_i64().unwrap() > 1,
            "the model must branch to reach the worker loop: {r}"
        );
    }
    let a = solve(&path, 4);
    let b = solve(&path, 4);
    assert_eq!(a["stats"]["nodes"], b["stats"]["nodes"], "{a} {b}");
    assert_eq!(a["objective"], b["objective"]);
    let _ = std::fs::remove_dir_all(&dir);
}
