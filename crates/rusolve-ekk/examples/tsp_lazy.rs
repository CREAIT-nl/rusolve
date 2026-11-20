//! Example: TSP with lazy subtour-elimination constraints.
//! `cargo run --release -p rusolve-ekk --example tsp_lazy [n] [seed]`
//! Only the degree rows are in the model; `SubtourOracle` returns one row per subtour of each
//! candidate (Dantzig, Fulkerson, Johnson 1954) and `run_with_lazy_constraints` re-solves.

#[path = "common/mod.rs"]
mod common;

use common::*;

fn main() {
    let mut args = std::env::args().skip(1);
    let n: usize = args.next().and_then(|a| a.parse().ok()).unwrap_or(9);
    let seed: u64 = args.next().and_then(|a| a.parse().ok()).unwrap_or(4);
    let pts = cities(n, seed);
    let mut h = quiet(tsp_model(&pts, false));
    let mut oracle = SubtourOracle::new(n);
    let (_, rep) = h.run_with_lazy_constraints(&mut oracle, 100);
    let obj = h.get_info().objective_function_value;
    println!("cities {n}  lazy tour length {obj:.6}  {rep:?}");
    if n <= 9 {
        println!("brute force optimum {:.6}", brute_force_tour(&pts));
    }
}
