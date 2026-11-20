//! `cargo run --release -p rusolve-mip --example node_lp_repro -- model.mps [out.mps]`
//!
//! Solves a MIP with node-LP presolve on and off; when a node LP returns an
//! infeasible "optimal" point, writes that node's LP (integrality dropped,
//! node bounds applied) to `out.mps` as a standalone LP repro.
use rusolve_mip::{MipConfig, MipError, MipProblem, solve};
fn main() {
    let path = std::env::args().nth(1).unwrap();
    let out = std::env::args().nth(2);
    let mut m = rusolve_io::read_mps(&std::fs::read_to_string(path).unwrap()).unwrap();
    let p = MipProblem::from_model(&m);
    for presolve in [true, false] {
        let mut cfg = MipConfig::default();
        cfg.lp.presolve = presolve;
        match solve(&p, &cfg) {
            Ok(s) => println!(
                "lp presolve={presolve}: {:?} obj={} nodes={}",
                s.status, s.objective, s.nodes
            ),
            Err(MipError::InvalidLpSolution {
                violation,
                lower,
                upper,
            }) => {
                println!("lp presolve={presolve}: invalid LP point, violation {violation}");
                if let Some(out) = &out {
                    for (j, b) in m.bounds.iter_mut().enumerate() {
                        b.lower = lower[j];
                        b.upper = upper[j];
                        b.integer = false;
                    }
                    std::fs::write(out, rusolve_io::write_mps(&m)).unwrap();
                    println!("node LP written to {out}");
                }
            }
            Err(e) => println!("lp presolve={presolve}: error {e}"),
        }
    }
}
