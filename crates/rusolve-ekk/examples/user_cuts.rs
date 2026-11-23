//! Example: user cuts equal the same cuts given in the model.
//! `cargo run --release -p rusolve-ekk --example user_cuts`
//! A 14-item knapsack; its minimal cover inequalities are (a) rows of the model and (b) returned
//! by a `MipHooks::separate` hook. Both reach the plain optimum.

#[path = "common/mod.rs"]
mod common;

use common::*;
use rusolve_ekk::mip::hooks::*;

struct Covers(Vec<Vec<usize>>);
impl MipHooks for Covers {
    fn separate(&mut self, _: &SeparationInfo<'_>) -> Vec<UserCut> {
        cover_cuts(&self.0)
    }
}

fn main() {
    let (v, w, cap) = knapsack_data();
    let covers = minimal_covers(&w, cap);
    let solve = |cuts: &[Vec<usize>], hook: bool| {
        let mut h = quiet_with(knapsack_model(&v, &w, cap, cuts), false);
        if hook {
            h.set_mip_hooks(Box::new(Covers(covers.clone())));
        }
        h.run();
        (h.get_info().objective_function_value, h.mip_hook_stats())
    };
    println!("{} minimal covers", covers.len());
    println!("plain          {:?}", solve(&[], false).0);
    println!("cuts in model  {:?}", solve(&covers, false).0);
    let (o, st) = solve(&[], true);
    println!(
        "user separator {o:?}  accepted {} cuts",
        st.user_cuts_accepted
    );
}
