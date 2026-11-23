use super::*;
use rusolve_ekk::lp_data::hconst::HighsModelStatus;
use rusolve_ekk::lp_data::highs_options::HighsOptions;
use rusolve_ekk::lp_data::highs_run::Highs;
use rusolve_ekk::lp_data::highs_status::HighsStatus;

const INF: f64 = f64::INFINITY;

fn cols(b: &mut ModelBuilder, n: usize) -> Range<usize> {
    b.add_cols(
        n,
        Fill::Scalar(0.0),
        Fill::Scalar(0.0),
        Fill::Scalar(INF),
        FillBool::Scalar(false),
    )
    .unwrap()
}

fn row_of(b: &ModelBuilder, r: usize) -> Vec<(i32, f64)> {
    let (s, e) = (b.row_start[r], b.row_start[r + 1]);
    let mut v: Vec<_> = (s..e).map(|p| (b.row_index[p], b.row_value[p])).collect();
    v.sort_by_key(|x| x.0);
    v
}

fn batch_rows(b: &mut ModelBuilder) -> RowBatch {
    let mut rb = b.take().rows;
    // canonical order inside rows, to compare paths that order entries differently
    for r in 0..rb.lower.len() {
        let (s, e) = (rb.start[r] as usize, rb.start[r + 1] as usize);
        let mut v: Vec<_> = (s..e).map(|p| (rb.index[p], rb.value[p])).collect();
        v.sort_by_key(|x| x.0);
        for (k, (j, x)) in v.into_iter().enumerate() {
            rb.index[s + k] = j;
            rb.value[s + k] = x;
        }
    }
    rb
}

#[test]
fn duplicates_merged_zero_dropped() {
    let mut b = ModelBuilder::new();
    cols(&mut b, 5);
    b.push(1, 2.0).unwrap();
    b.push(3, 1.0).unwrap();
    b.push(1, 3.0).unwrap();
    b.push(4, 7.0).unwrap();
    b.push(4, -7.0).unwrap(); // a + (-a)
    b.push(0, 0.0).unwrap(); // explicit zero
    assert_eq!(b.finish_row(-INF, 4.0).unwrap(), 0);
    assert_eq!(row_of(&b, 0), vec![(1, 5.0), (3, 1.0)]);
    assert_eq!(b.num_pending_nonzeros(), 2);
    // marks are clear: a second row works the same
    b.push(1, 1.0).unwrap();
    assert_eq!(b.finish_row(0.0, 0.0).unwrap(), 1);
    assert_eq!(row_of(&b, 1), vec![(1, 1.0)]);
    assert!(b.mark.iter().all(|&m| m == 0));
}

#[test]
fn empty_and_ranged_rows() {
    let mut b = ModelBuilder::new();
    cols(&mut b, 2);
    assert_eq!(b.finish_row(-INF, INF).unwrap(), 0);
    b.push(0, 1.0).unwrap();
    b.push(1, 1.0).unwrap();
    assert_eq!(b.finish_row(1.0, 4.0).unwrap(), 1);
    let t = b.take();
    assert_eq!(t.rows.lower, vec![-INF, 1.0]);
    assert_eq!(t.rows.upper, vec![INF, 4.0]);
    assert_eq!(t.rows.start, vec![0, 0, 2]);
}

#[test]
fn push_errors_name_row_and_col() {
    let mut b = ModelBuilder::new();
    cols(&mut b, 8);
    for _ in 0..3 {
        b.finish_row(0.0, 1.0).unwrap();
    }
    let e = b.push(7, f64::NAN).unwrap_err();
    assert_eq!(e.to_string(), "coefficient of column 7 in row 3 is NaN");
    let e = b.push(2, f64::INFINITY).unwrap_err();
    assert_eq!(
        e.to_string(),
        "coefficient of column 2 in row 3 is infinite"
    );
    assert_eq!(
        b.push(8, 1.0).unwrap_err(),
        BuildError::ColumnOutOfRange {
            row: Some(3),
            col: 8,
            num_cols: 8
        }
    );
    assert!(matches!(
        b.push(-1, 1.0),
        Err(BuildError::ColumnOutOfRange { col: -1, .. })
    ));
    assert!(b.push(8, 1.0).unwrap_err().to_string().contains("row 3"));
    // nothing was pushed
    b.finish_row(0.0, 0.0).unwrap();
    assert_eq!(b.num_pending_nonzeros(), 0);
}

#[test]
fn finish_row_nan_bounds_abort_row() {
    let mut b = ModelBuilder::new();
    cols(&mut b, 2);
    b.push(0, 1.0).unwrap();
    let e = b.finish_row(f64::NAN, 1.0).unwrap_err();
    assert_eq!(e.to_string(), "row lower bound of row 0 is NaN");
    assert_eq!(b.num_rows(), 0);
    assert_eq!(b.row_index.len(), 0);
    b.push(0, 1.0).unwrap();
    assert!(b.finish_row(0.0, f64::NAN).is_err());
    b.push(1, 1.0).unwrap();
    b.abort_row();
    b.finish_row(0.0, 1.0).unwrap();
    assert_eq!(b.num_pending_nonzeros(), 0);
}

#[test]
fn merged_overflow_to_inf_is_error() {
    // inf can't be pushed, but finite duplicates can sum to inf: that is an infinite
    // coefficient (python-api.md, Semantics), so the row is refused and nothing is added.
    let mut b = ModelBuilder::new();
    cols(&mut b, 1);
    b.push(0, f64::MAX).unwrap();
    b.push(0, f64::MAX).unwrap();
    let e = b.finish_row(0.0, 1.0).unwrap_err();
    assert!(
        matches!(
            e,
            BuildError::Infinite {
                row: Some(0),
                col: Some(0),
                ..
            }
        ),
        "{e:?}"
    );
    assert_eq!(b.num_rows(), 0);
    assert_eq!(b.num_pending_nonzeros(), 0);
    // the bulk path refuses it too, all or nothing
    let e = b
        .add_rows_csr(
            &[0i32, 2],
            &[0i32, 0],
            &[f64::MAX, f64::MAX],
            None,
            Fill::Scalar(0.0),
            Fill::Scalar(1.0),
        )
        .unwrap_err();
    assert!(matches!(e, BuildError::Infinite { .. }), "{e:?}");
    assert_eq!(b.num_rows(), 0);
}

#[test]
fn add_cols_checks() {
    let mut b = ModelBuilder::new();
    let e = b
        .add_cols(
            3,
            Fill::Slice(&[1.0, f64::NAN, 2.0]),
            Fill::Scalar(0.0),
            Fill::Scalar(1.0),
            FillBool::Scalar(false),
        )
        .unwrap_err();
    assert_eq!(e.to_string(), "cost of column 1 is NaN");
    let e = b
        .add_cols(
            1,
            Fill::Scalar(f64::INFINITY),
            Fill::Scalar(0.0),
            Fill::Scalar(1.0),
            FillBool::Scalar(false),
        )
        .unwrap_err();
    assert_eq!(e.to_string(), "cost of column 0 is infinite");
    let e = b
        .add_cols(
            2,
            Fill::Scalar(0.0),
            Fill::Scalar(f64::NAN),
            Fill::Scalar(1.0),
            FillBool::Scalar(false),
        )
        .unwrap_err();
    assert_eq!(e.to_string(), "lower bound of column 0 is NaN");
    let e = b
        .add_cols(
            2,
            Fill::Scalar(0.0),
            Fill::Scalar(0.0),
            Fill::Slice(&[1.0]),
            FillBool::Scalar(false),
        )
        .unwrap_err();
    assert_eq!(e.to_string(), "upper bound has length 1, expected 2");
    let e = b
        .add_cols(
            2,
            Fill::Scalar(0.0),
            Fill::Scalar(0.0),
            Fill::Scalar(1.0),
            FillBool::Slice(&[true]),
        )
        .unwrap_err();
    assert!(matches!(e, BuildError::Length { .. }));
    assert_eq!(b.num_cols(), 0);
    // infinite bounds and lower > upper are fine
    b.add_cols(
        2,
        Fill::Scalar(1.0),
        Fill::Slice(&[-INF, 5.0]),
        Fill::Slice(&[INF, 1.0]),
        FillBool::Scalar(false),
    )
    .unwrap();
    assert_eq!(b.num_cols(), 2);
}

/// A 3x4 matrix with a duplicate-free layout, as dense / CSR / CSC / COO.
fn sample() -> (Vec<f64>, Vec<i32>, Vec<i32>, Vec<f64>) {
    // [1 0 2 0; 0 3 0 0; 4 0 0 5]
    let dense = vec![1., 0., 2., 0., 0., 3., 0., 0., 4., 0., 0., 5.];
    let indptr = vec![0, 2, 3, 5];
    let indices = vec![0, 2, 1, 0, 3];
    let data = vec![1., 2., 3., 4., 5.];
    (dense, indptr, indices, data)
}

fn fresh(n: usize) -> ModelBuilder {
    let mut b = ModelBuilder::new();
    cols(&mut b, n);
    b
}

#[test]
fn csr_csc_coo_dense_agree() {
    let (dense, indptr, indices, data) = sample();
    let lo = [-INF, 0.0, 1.0];
    let hi = [1.0, 0.0, INF];
    let (lo, hi) = (Fill::Slice(&lo), Fill::Slice(&hi));

    let mut b = fresh(4);
    assert_eq!(
        b.add_rows_csr(&indptr, &indices, &data, None, lo, hi)
            .unwrap(),
        0..3
    );
    let want = batch_rows(&mut b);
    assert_eq!(want.start, vec![0, 2, 3, 5]);

    // CSC
    let cptr = vec![0i64, 2, 3, 4, 5];
    let cidx = vec![0i64, 2, 1, 0, 2];
    let cdat = vec![1., 4., 3., 2., 5.];
    let mut b = fresh(4);
    b.add_rows_csc(3, &cptr, &cidx, &cdat, None, lo, hi)
        .unwrap();
    assert_eq!(batch_rows(&mut b), want);

    // COO, scrambled order
    let r = [2i32, 0, 1, 2, 0];
    let c = [3i64, 2, 1, 0, 0];
    let d = [5., 2., 3., 4., 1.];
    let mut b = fresh(4);
    b.add_rows_coo(3, &r, &c, &d, None, lo, hi).unwrap();
    assert_eq!(batch_rows(&mut b), want);

    // dense
    let mut b = fresh(4);
    b.add_rows_dense(3, 4, &dense, None, lo, hi).unwrap();
    assert_eq!(batch_rows(&mut b), want);
}

#[test]
fn coo_duplicates_and_zeros() {
    let r = [0i32, 0, 1, 1, 1, 2];
    let c = [1i32, 1, 0, 0, 2, 3];
    let d = [1.0, 2.0, 5.0, -5.0, 0.0, 9.0];
    let mut b = fresh(4);
    b.add_rows_coo(3, &r, &c, &d, None, Fill::Scalar(0.0), Fill::Scalar(1.0))
        .unwrap();
    let t = b.take().rows;
    assert_eq!(t.start, vec![0, 1, 1, 2]);
    assert_eq!(t.index, vec![1, 3]);
    assert_eq!(t.value, vec![3.0, 9.0]);
}

#[test]
fn csr_duplicates_zero_fast_and_slow_path() {
    // row 0 clean, row 1 has duplicate, row 2 has explicit zero, row 3 cancels, row 4 empty
    let indptr = [0i32, 2, 5, 7, 9, 9];
    let indices = [0i32, 1, 2, 2, 0, 1, 0, 3, 3];
    let data = [1., 2., 1., 4., 3., 0.0, 6., 2., -2.];
    let mut b = fresh(4);
    b.add_rows_csr(
        &indptr,
        &indices,
        &data,
        None,
        Fill::Scalar(0.0),
        Fill::Scalar(0.0),
    )
    .unwrap();
    assert_eq!(row_of(&b, 0), vec![(0, 1.0), (1, 2.0)]);
    assert_eq!(row_of(&b, 1), vec![(0, 3.0), (2, 5.0)]);
    assert_eq!(row_of(&b, 2), vec![(0, 6.0)]);
    assert_eq!(row_of(&b, 3), vec![]);
    assert_eq!(row_of(&b, 4), vec![]);
    assert!(b.mark.iter().all(|&m| m == 0));
}

#[test]
fn bulk_paths_are_all_or_nothing() {
    let mut b = fresh(3);
    b.push(0, 1.0).unwrap();
    b.finish_row(0.0, 1.0).unwrap();
    let before = (
        b.row_index.clone(),
        b.row_start.clone(),
        b.row_lower.clone(),
    );
    let f = Fill::Scalar(0.0);
    // CSR: bad column in last row, NaN in last row
    let e = b
        .add_rows_csr(&[0i32, 1, 2], &[0i32, 5], &[1.0, 1.0], None, f, f)
        .unwrap_err();
    assert_eq!(
        e,
        BuildError::ColumnOutOfRange {
            row: Some(2),
            col: 5,
            num_cols: 3
        }
    );
    let e = b
        .add_rows_csr(&[0i32, 1, 2], &[0i32, 1], &[1.0, f64::NAN], None, f, f)
        .unwrap_err();
    assert_eq!(e.to_string(), "coefficient of column 1 in row 2 is NaN");
    let e = b
        .add_rows_csr(
            &[0i32, 1, 2],
            &[0i32, 1],
            &[1.0, 1.0],
            None,
            Fill::Slice(&[0.0, f64::NAN]),
            f,
        )
        .unwrap_err();
    assert_eq!(e.to_string(), "row lower bound of row 2 is NaN");
    // CSC / COO / dense
    assert!(
        b.add_rows_csc(1, &[0i32, 1, 2], &[0i32, 3], &[1., 1.], None, f, f)
            .is_err()
    );
    assert!(
        b.add_rows_coo(2, &[0i32, 1], &[0i32, 9], &[1., 1.], None, f, f)
            .is_err()
    );
    assert!(
        b.add_rows_coo(2, &[0i32, 2], &[0i32, 1], &[1., 1.], None, f, f)
            .is_err()
    );
    assert!(
        b.add_rows_dense(2, 3, &[1., 0., 0., 0., f64::INFINITY, 0.], None, f, f)
            .is_err()
    );
    assert!(b.add_rows_dense(2, 3, &[1.; 5], None, f, f).is_err());
    assert!(
        b.add_rows_csr(&[0i32, 2, 1], &[0i32, 1], &[1., 1.], None, f, f)
            .is_err()
    );
    assert!(
        b.add_rows_csr(&[1i32, 2], &[0i32, 1], &[1., 1.], None, f, f)
            .is_err()
    );
    assert_eq!(
        before,
        (
            b.row_index.clone(),
            b.row_start.clone(),
            b.row_lower.clone()
        )
    );
    assert_eq!(b.row_value.len(), 1);
    assert_eq!(b.row_upper.len(), 1);
    assert!(b.mark.iter().all(|&m| m == 0));
    // and the builder still works
    b.add_rows_csr(&[0i32, 1], &[2i32], &[4.0], None, f, f)
        .unwrap();
    assert_eq!(b.num_rows(), 2);
}

#[test]
fn open_row_blocks_bulk() {
    let mut b = fresh(2);
    b.push(0, 1.0).unwrap();
    let f = Fill::Scalar(0.0);
    assert!(
        b.add_rows_csr(&[0i32, 1], &[0i32], &[1.0], None, f, f)
            .is_err()
    );
}

#[test]
fn col_map_all_paths() {
    // matrix columns 0,1,2 -> model columns 4,0,2
    let map = [4i32, 0, 2];
    let f = Fill::Scalar(0.0);
    let mut want = None;
    for path in 0..4 {
        let mut b = fresh(5);
        match path {
            0 => b
                .add_rows_csr(
                    &[0i32, 2, 3],
                    &[0i32, 2, 1],
                    &[1., 2., 3.],
                    Some(&map),
                    f,
                    f,
                )
                .unwrap(),
            1 => b
                .add_rows_csc(
                    2,
                    &[0i32, 1, 2, 3],
                    &[0i32, 1, 0],
                    &[1., 3., 2.],
                    Some(&map),
                    f,
                    f,
                )
                .unwrap(),
            2 => b
                .add_rows_coo(
                    2,
                    &[0i32, 1, 0],
                    &[0i32, 1, 2],
                    &[1., 3., 2.],
                    Some(&map),
                    f,
                    f,
                )
                .unwrap(),
            _ => b
                .add_rows_dense(2, 3, &[1., 0., 2., 0., 3., 0.], Some(&map), f, f)
                .unwrap(),
        };
        let rb = batch_rows(&mut b);
        assert_eq!(rb.start, vec![0, 2, 3]);
        assert_eq!(&rb.index[..2], &[2, 4]);
        assert_eq!(rb.index[2], 0);
        assert_eq!(rb.value, vec![2., 1., 3.]);
        if let Some(w) = &want {
            assert_eq!(&rb, w);
        }
        want = Some(rb);
    }
    // length mismatches and bad targets
    let mut b = fresh(5);
    assert!(matches!(
        b.add_rows_csc(
            1,
            &[0i32, 1, 2, 3],
            &[0i32; 3],
            &[1.; 3],
            Some(&[0, 1]),
            f,
            f
        ),
        Err(BuildError::Length { .. })
    ));
    assert!(matches!(
        b.add_rows_dense(1, 3, &[1.; 3], Some(&[0, 1]), f, f),
        Err(BuildError::Length { .. })
    ));
    // CSR / COO: index beyond col_map
    assert!(matches!(
        b.add_rows_csr(&[0i32, 1], &[2i32], &[1.], Some(&[0, 1]), f, f),
        Err(BuildError::ColumnOutOfRange {
            col: 2,
            num_cols: 2,
            ..
        })
    ));
    assert!(matches!(
        b.add_rows_coo(1, &[0i32], &[2i32], &[1.], Some(&[0, 1]), f, f),
        Err(BuildError::ColumnOutOfRange { .. })
    ));
    // col_map entry beyond the model
    assert!(matches!(
        b.add_rows_csr(&[0i32, 1], &[0i32], &[1.], Some(&[9]), f, f),
        Err(BuildError::ColumnOutOfRange {
            col: 9,
            num_cols: 5,
            ..
        })
    ));
    assert_eq!(b.num_rows(), 0);
}

#[test]
fn col_map_duplicates_merge() {
    // two matrix columns mapped to the same model column sum up
    let f = Fill::Scalar(0.0);
    let mut b = fresh(2);
    b.add_rows_dense(1, 3, &[1., 2., 3.], Some(&[1, 1, 0]), f, f)
        .unwrap();
    assert_eq!(row_of(&b, 0), vec![(0, 3.0), (1, 3.0)]);
}

#[test]
fn take_twice_persists_indices() {
    let mut b = fresh(3);
    b.add_row(&[0, 2], &[1.0, 2.0], 0.0, 5.0).unwrap();
    let t1 = b.take();
    assert_eq!((t1.cols.first, t1.rows.first), (0, 0));
    assert_eq!(t1.cols.cost.len(), 3);
    assert_eq!((b.synced_cols(), b.synced_rows()), (3, 1));
    assert_eq!((b.num_cols(), b.num_rows()), (3, 1));
    // a second take with nothing new is empty
    let e = b.take();
    assert_eq!((e.cols.first, e.rows.first), (3, 1));
    assert!(e.cols.cost.is_empty() && e.rows.lower.is_empty());
    assert_eq!(e.rows.start, vec![0]);
    // new column and rows reference old and new columns by global index
    let r = cols(&mut b, 2);
    assert_eq!(r, 3..5);
    b.add_row(&[1, 4], &[1.0, 1.0], -INF, 7.0).unwrap();
    assert_eq!(b.add_row(&[3], &[1.0], 0.0, 0.0).unwrap(), 2);
    assert_eq!(b.num_pending_nonzeros(), 3);
    let t2 = b.take();
    assert_eq!((t2.cols.first, t2.rows.first), (3, 1));
    assert_eq!(t2.cols.cost.len(), 2);
    assert_eq!(t2.rows.index, vec![1, 4, 3]);
    assert_eq!(t2.rows.start, vec![0, 2, 3]);
    // old entities are synced
    assert_eq!(
        b.set_col_cost(0, 1.0),
        Err(BuildError::Synced {
            what: "column",
            index: 0
        })
    );
    assert_eq!(
        b.set_row_bounds(0, 0.0, 1.0),
        Err(BuildError::Synced {
            what: "row",
            index: 0
        })
    );
}

#[test]
fn names_and_defaults() {
    let mut b = fresh(2);
    assert!(b.take().cols.names.is_empty());
    let r = cols(&mut b, 3); // 2..5
    b.set_col_names(r.start + 1..r.end, Names::Prefix("x".into()))
        .unwrap();
    assert_eq!(b.col_name(3), Some("x[0]"));
    assert_eq!(b.col_name(2), None);
    b.add_row(&[2], &[1.0], 0.0, 1.0).unwrap();
    b.add_row(&[3], &[1.0], 0.0, 1.0).unwrap();
    b.set_row_names(1..2, Names::List(vec!["cap".into()]))
        .unwrap();
    assert_eq!(b.row_name(1), Some("cap"));
    assert_eq!(b.row_name(0), None);
    b.add_row(&[4], &[1.0], 0.0, 1.0).unwrap(); // added after naming
    let t = b.take();
    assert_eq!(t.cols.names, vec!["c2", "x[0]", "x[1]"]);
    assert_eq!(t.rows.names, vec!["r0", "cap", "r2"]);
    // columns added after naming also get defaults
    cols(&mut b, 1);
    assert_eq!(b.take().cols.names, vec!["c5"]);
    // errors
    cols(&mut b, 1);
    assert!(matches!(
        b.set_col_names(6..7, Names::List(vec![])),
        Err(BuildError::Length { .. })
    ));
    assert!(matches!(
        b.set_col_names(0..1, Names::Prefix("y".into())),
        Err(BuildError::Synced { .. })
    ));
    assert!(matches!(
        b.set_col_names(6..9, Names::Prefix("y".into())),
        Err(BuildError::ColumnOutOfRange { .. })
    ));
    assert!(matches!(
        b.set_row_names(0..1, Names::Prefix("y".into())),
        Err(BuildError::Synced { .. })
    ));
    assert!(matches!(
        b.set_row_names(3..4, Names::Prefix("y".into())),
        Err(BuildError::RowOutOfRange { .. })
    ));
}

#[test]
fn integrality() {
    let mut b = ModelBuilder::new();
    b.add_cols(
        2,
        Fill::Scalar(0.0),
        Fill::Scalar(0.0),
        Fill::Scalar(1.0),
        FillBool::Scalar(false),
    )
    .unwrap();
    let t = b.take();
    assert!(t.cols.integrality.is_empty());
    b.add_cols(
        3,
        Fill::Scalar(0.0),
        Fill::Scalar(0.0),
        Fill::Scalar(1.0),
        FillBool::Slice(&[false, true, false]),
    )
    .unwrap();
    b.set_col_integer(4, true).unwrap();
    let t = b.take();
    use HighsVarType::*;
    assert_eq!(t.cols.integrality, vec![kContinuous, kInteger, kInteger]);
    // once ever integer, always filled
    cols(&mut b, 1);
    assert_eq!(b.take().cols.integrality, vec![kContinuous]);
}

#[test]
fn semi_columns_override_integrality() {
    let mut b = ModelBuilder::new();
    cols(&mut b, 4);
    b.set_col_integer(1, true).unwrap();
    b.set_col_semi(0, false).unwrap();
    b.set_col_semi(1, true).unwrap();
    b.set_col_semi(2, true).unwrap();
    // a later set_col_integer replaces the semi kind
    b.set_col_integer(2, false).unwrap();
    use HighsVarType::*;
    assert_eq!(
        b.take().cols.integrality,
        vec![kSemiContinuous, kSemiInteger, kContinuous, kContinuous]
    );
}

#[test]
fn into_lp_matches_hand_built() {
    let mut b = ModelBuilder::new();
    b.add_cols(
        3,
        Fill::Slice(&[1.0, -2.0, 0.0]),
        Fill::Scalar(0.0),
        Fill::Slice(&[4.0, INF, 1.0]),
        FillBool::Slice(&[false, false, true]),
    )
    .unwrap();
    b.add_row(&[0, 1], &[1.0, 2.0], -INF, 8.0).unwrap();
    b.add_row(&[2, 1, 0], &[3.0, 4.0, 5.0], 1.0, 9.0).unwrap();
    b.add_row(&[], &[], 0.0, 0.0).unwrap();
    b.set_sense(Sense::Maximize);
    b.set_offset(2.5);
    b.set_name("m".into());
    b.set_col_names(0..3, Names::Prefix("v".into())).unwrap();
    let lp = b.take().into_lp();
    let want = HighsLp {
        num_col_: 3,
        num_row_: 3,
        col_cost_: vec![1.0, -2.0, 0.0],
        col_lower_: vec![0.0; 3],
        col_upper_: vec![4.0, INF, 1.0],
        row_lower_: vec![-INF, 1.0, 0.0],
        row_upper_: vec![8.0, 9.0, 0.0],
        a_matrix_: rusolve_ekk::util::sparse_matrix::HighsSparseMatrix {
            format_: MatrixFormat::kColwise,
            num_col_: 3,
            num_row_: 3,
            start_: vec![0, 2, 4, 5],
            p_end_: vec![],
            index_: vec![0, 1, 0, 1, 1],
            value_: vec![1.0, 5.0, 2.0, 4.0, 3.0],
        },
        sense_: ObjSense::kMaximize,
        offset_: 2.5,
        model_name_: "m".into(),
        col_names_: vec!["v[0]".into(), "v[1]".into(), "v[2]".into()],
        row_names_: vec![],
        integrality_: vec![
            HighsVarType::kContinuous,
            HighsVarType::kContinuous,
            HighsVarType::kInteger,
        ],
        ..Default::default()
    };
    assert_eq!(format!("{lp:?}"), format!("{want:?}"));
}

#[test]
fn set_coef_insert_update_remove() {
    let mut b = fresh(4);
    b.add_row(&[0, 1], &[1.0, 2.0], 0.0, 1.0).unwrap();
    b.add_row(&[2], &[5.0], 0.0, 1.0).unwrap();
    b.set_coef(0, 1, 9.0).unwrap(); // update
    assert_eq!(row_of(&b, 0), vec![(0, 1.0), (1, 9.0)]);
    b.set_coef(0, 3, 4.0).unwrap(); // insert
    assert_eq!(row_of(&b, 0), vec![(0, 1.0), (1, 9.0), (3, 4.0)]);
    assert_eq!(row_of(&b, 1), vec![(2, 5.0)]);
    b.set_coef(0, 0, 0.0).unwrap(); // remove
    assert_eq!(row_of(&b, 0), vec![(1, 9.0), (3, 4.0)]);
    assert_eq!(row_of(&b, 1), vec![(2, 5.0)]);
    b.set_coef(0, 2, 0.0).unwrap(); // remove absent
    assert_eq!(b.num_pending_nonzeros(), 3);
    assert!(matches!(
        b.set_coef(0, 0, f64::NAN),
        Err(BuildError::Nan { .. })
    ));
    assert!(matches!(
        b.set_coef(0, 4, 1.0),
        Err(BuildError::ColumnOutOfRange { .. })
    ));
    assert!(matches!(
        b.set_coef(2, 0, 1.0),
        Err(BuildError::RowOutOfRange { .. })
    ));
    b.take();
    assert_eq!(
        b.set_coef(0, 0, 1.0),
        Err(BuildError::Synced {
            what: "row",
            index: 0
        })
    );
}

#[test]
fn set_bounds_cost_clear() {
    let mut b = fresh(2);
    b.set_col_bounds(1, -1.0, 3.0).unwrap();
    b.set_col_cost(0, 2.0).unwrap();
    assert!(b.set_col_cost(0, f64::NAN).is_err());
    assert!(b.set_col_bounds(0, f64::NAN, 1.0).is_err());
    assert!(b.set_col_cost(2, 1.0).is_err());
    b.clear_unsynced_costs();
    let t = b.take();
    assert_eq!(t.cols.cost, vec![0.0, 0.0]);
    assert_eq!((t.cols.lower[1], t.cols.upper[1]), (-1.0, 3.0));
    assert_eq!(
        b.set_col_bounds(0, 0.0, 1.0),
        Err(BuildError::Synced {
            what: "column",
            index: 0
        })
    );
    assert_eq!(
        b.set_col_integer(1, true),
        Err(BuildError::Synced {
            what: "column",
            index: 1
        })
    );
}

#[test]
fn display_sentences() {
    let e = BuildError::Nan {
        what: "cost",
        row: None,
        col: Some(2),
    };
    assert_eq!(e.to_string(), "cost of column 2 is NaN");
    let e = BuildError::RowOutOfRange {
        row: 4,
        num_rows: 3,
    };
    assert_eq!(
        e.to_string(),
        "row index 4 is out of range: the model has 3 rows"
    );
    let e = BuildError::Synced {
        what: "row",
        index: 0,
    };
    assert!(e.to_string().starts_with("row 0 was already"));
}

fn solve(lp: HighsLp) -> (Highs, f64) {
    let mut o = HighsOptions::recommended();
    o.output_flag = false;
    o.log_to_console = false;
    o.sync_log_options();
    let mut h = Highs::with_options(o);
    assert_eq!(h.pass_model(lp), HighsStatus::kOk);
    h.run();
    assert_eq!(h.get_model_status(), HighsModelStatus::kOptimal);
    let obj = h.get_info().objective_function_value;
    (h, obj)
}

/// Knapsack: items (value, weight) a(10,5) b(13,7) c(7,4) d(4,3), capacity 14.
#[test]
fn end_to_end_knapsack_and_incremental() {
    let mut b = ModelBuilder::new();
    b.add_cols(
        4,
        Fill::Slice(&[10.0, 13.0, 7.0, 4.0]),
        Fill::Scalar(0.0),
        Fill::Scalar(1.0),
        FillBool::Scalar(true),
    )
    .unwrap();
    b.set_sense(Sense::Maximize);
    b.add_row(&[0, 1, 2, 3], &[5.0, 7.0, 4.0, 3.0], -INF, 14.0)
        .unwrap();
    // best set is b+c+d: 13+7+4 = 24 (weight 14)
    let (mut h, obj) = solve(b.take().into_lp());
    assert!((obj - 24.0).abs() < 1e-6, "obj {obj}");

    // second take: b + c <= 1; then the best is a+b = 23 (weight 12)
    b.add_row(&[1, 2], &[1.0, 1.0], -INF, 1.0).unwrap();
    let t = b.take();
    assert_eq!((t.rows.first, t.cols.first), (1, 4));
    let st = h.add_rows(
        1,
        &t.rows.lower,
        &t.rows.upper,
        t.rows.index.len() as HighsInt,
        &t.rows.start[..1],
        &t.rows.index,
        &t.rows.value,
    );
    assert_eq!(st, HighsStatus::kOk);
    h.run();
    assert_eq!(h.get_model_status(), HighsModelStatus::kOptimal);
    let obj = h.get_info().objective_function_value;
    assert!((obj - 23.0).abs() < 1e-6, "obj {obj}");
}
