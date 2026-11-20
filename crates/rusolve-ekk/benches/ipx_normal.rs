//! IPX normal-equations kernels: `NormalMatrix::_Apply` (A*W*A' times a vector)
//! and a fixed-iteration preconditioned CR solve with `DiagonalPrecond`, the
//! inner loop of every `KKTSolverDiag` Newton solve. Synthetic LP with a
//! deterministic sparse pattern.

use std::hint::black_box;
use std::rc::Rc;

use criterion::{Criterion, criterion_group, criterion_main};
use rusolve_ekk::ipm::ipx::conjugate_residuals::ConjugateResiduals;
use rusolve_ekk::ipm::ipx::control::Control;
use rusolve_ekk::ipm::ipx::diagonal_precond::DiagonalPrecond;
use rusolve_ekk::ipm::ipx::linear_operator::LinearOperator;
use rusolve_ekk::ipm::ipx::model::Model;
use rusolve_ekk::ipm::ipx::normal_matrix::NormalMatrix;
use rusolve_ekk::ipm::ipx::{Info, Vector};

/// `m` equality rows, `n` columns with `per_col` nonzeros each (LCG pattern).
fn synthetic_model(control: &Control, m: i32, n: i32, per_col: i32) -> Model {
    let mut state: u64 = 0x9E37_79B9_7F4A_7C15;
    let mut next = || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (state >> 33) as u32
    };
    let (mut ap, mut ai, mut ax) = (vec![0i32], Vec::new(), Vec::new());
    for _ in 0..n {
        let mut rows: Vec<i32> = (0..per_col).map(|_| (next() % m as u32) as i32).collect();
        rows.sort_unstable();
        rows.dedup();
        for r in rows {
            ai.push(r);
            ax.push(1.0 + (next() % 1000) as f64 / 250.0);
        }
        ap.push(ai.len() as i32);
    }
    let rhs = vec![1.0; m as usize];
    let ct = vec![b'='; m as usize];
    let obj: Vec<f64> = (0..n).map(|j| 1.0 + (j % 7) as f64).collect();
    let lb = vec![0.0; n as usize];
    let ub = vec![f64::INFINITY; n as usize];
    let mut model = Model::new();
    let err = model.load(
        control,
        m,
        n,
        Some(&ap),
        Some(&ai),
        Some(&ax),
        Some(&rhs),
        Some(&ct),
        0.0,
        Some(&obj),
        Some(&lb),
        Some(&ub),
    );
    assert_eq!(err, 0);
    model
}

fn bench(c: &mut Criterion) {
    let control = Control::new();
    let model = Rc::new(synthetic_model(&control, 4000, 12000, 6));
    let m = model.rows();
    let nm = model.cols() + m;
    let w: Vec<f64> = (0..nm).map(|i| 0.5 + (i % 13) as f64 * 0.1).collect();
    let rhs = Vector((0..m).map(|i| 1.0 + (i % 5) as f64).collect());

    let mut normal = NormalMatrix::new(model.clone());
    normal.prepare(Some(&w));
    let mut lhs = Vector::new(m);
    c.bench_function("ipx_normal_matrix_apply", |b| {
        b.iter(|| {
            let mut d = 0.0;
            normal.apply(black_box(&rhs), &mut lhs, Some(&mut d));
            black_box(d)
        })
    });

    let mut precond = DiagonalPrecond::new(model.clone());
    let mut info = Info::default();
    precond.factorize(Some(&w), &mut info);
    let mut cr = ConjugateResiduals::new(Rc::new(Control::new()));
    c.bench_function("ipx_cr_precond_20iter", |b| {
        b.iter(|| {
            let mut x = Vector::new(m);
            cr.solve_precond(
                &mut normal,
                &mut precond,
                black_box(&rhs),
                0.0,
                None,
                20,
                &mut x,
            );
            black_box(x[0])
        })
    });
}

criterion_group!(benches, bench);
criterion_main!(benches);
