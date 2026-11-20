//! `ModelBuilder::truncate_rows`: drop unsynced rows added after a point, so a batch call
//! (`Model.add_constrs` in crates/rusolve-py) that fails half-way leaves the model as it was
//! (docs/design/python-api.md, "Semantics": a failed call changes nothing).

use super::ModelBuilder;

impl ModelBuilder {
    /// Drop every row with global index `>= rows` (and the open row's entries). Rows already
    /// taken are never dropped: `rows` below `synced_rows()` truncates to `synced_rows()`.
    pub fn truncate_rows(&mut self, rows: usize) {
        self.abort_row();
        let keep = rows.max(self.synced_rows).min(self.num_rows()) - self.synced_rows;
        let nz = self.row_start[keep];
        self.row_index.truncate(nz);
        self.row_value.truncate(nz);
        self.row_start.truncate(keep + 1);
        self.row_lower.truncate(keep);
        self.row_upper.truncate(keep);
        self.row_names.truncate(keep);
    }
}

#[cfg(test)]
mod tests {
    use super::super::{Fill, FillBool, ModelBuilder};

    #[test]
    fn truncate_rows_restores_the_row_count() {
        let mut b = ModelBuilder::new();
        b.add_cols(
            3,
            Fill::Scalar(0.0),
            Fill::Scalar(0.0),
            Fill::Scalar(1.0),
            FillBool::Scalar(false),
        )
        .unwrap();
        b.add_row(&[0, 1], &[1.0, 2.0], 0.0, 1.0).unwrap();
        b.add_row(&[2], &[3.0], 0.0, 1.0).unwrap();
        b.push(1, 5.0).unwrap();
        b.truncate_rows(1);
        assert_eq!(b.num_rows(), 1);
        assert_eq!(b.num_pending_nonzeros(), 2);
        b.add_row(&[2], &[4.0], 0.0, 2.0).unwrap();
        let batch = b.take();
        assert_eq!(batch.rows.start, vec![0, 2, 3]);
        assert_eq!(batch.rows.value, vec![1.0, 2.0, 4.0]);
        b.truncate_rows(0); // synced rows stay
        assert_eq!(b.num_rows(), 2);
    }
}
