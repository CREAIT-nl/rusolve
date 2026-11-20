//! Typed decoders for the composite objects of FORMAT.md. Field names follow
//! FORMAT.md / HiGHS (snake_cased only where Rust requires it); every field
//! doc says which JSON key it comes from.

use super::fixture::Rec;

/// `pack*` fields of an HVector (FORMAT.md: "only where the record says
/// `with pack`").
#[derive(Clone, Debug, PartialEq)]
pub struct PackRec {
    /// `packFlag`
    pub pack_flag: bool,
    /// `packCount`
    pub pack_count: i32,
    /// `packIndex[0..packCount)`
    pub pack_index: Vec<i32>,
    /// `packValue[0..packCount)`
    pub pack_value: Vec<f64>,
}

/// An `HVector` as dumped by the oracle (`Obj::hv` in emit.cpp).
#[derive(Clone, Debug, PartialEq)]
pub struct HVectorRec {
    /// `size`
    pub size: i32,
    /// `count`; negative = dense mode (indices unknown).
    pub count: i32,
    /// `index[0..count)` in HiGHS' exact order; `None` when `count < 0`.
    pub index: Option<Vec<i32>>,
    /// `synthetic_tick`: part of the state, summed by HEkk to decide reinversion.
    pub synthetic_tick: f64,
    /// `nz_i`: every position whose bit pattern is non-zero, ascending.
    pub nz_i: Vec<i32>,
    /// `nz_v`: the values at `nz_i` (may include -0.0 and kHighsZero).
    pub nz_v: Vec<f64>,
    /// `packFlag`/`packCount`/`packIndex`/`packValue`, when dumped.
    pub pack: Option<PackRec>,
}

impl HVectorRec {
    pub fn from_rec(r: &Rec) -> Self {
        let count = r.int("count");
        let index = r.has("index").then(|| r.ints("index"));
        if count >= 0 {
            assert!(index.is_some(), "{}: count >= 0 but no index", r.ctx);
        }
        let nz_i = r.ints("nz_i");
        let nz_v = r.f64s("nz_v");
        assert_eq!(
            nz_i.len(),
            nz_v.len(),
            "{}: nz_i/nz_v length mismatch",
            r.ctx
        );
        let pack = r.has("packFlag").then(|| PackRec {
            pack_flag: r.bool("packFlag"),
            pack_count: r.int("packCount"),
            pack_index: r.ints("packIndex"),
            pack_value: r.f64s("packValue"),
        });
        HVectorRec {
            size: r.int("size"),
            count,
            index,
            synthetic_tick: r.f64("synthetic_tick"),
            nz_i,
            nz_v,
            pack,
        }
    }

    /// The full `array` (length `size`): `nz_v` at `nz_i`, `+0.0` elsewhere.
    pub fn array(&self) -> Vec<f64> {
        let mut a = vec![0.0; self.size as usize];
        for (&i, &v) in self.nz_i.iter().zip(&self.nz_v) {
            a[i as usize] = v;
        }
        a
    }

    /// `index[0..count)` or an empty slice in dense mode.
    pub fn index_slice(&self) -> &[i32] {
        self.index.as_deref().unwrap_or(&[])
    }

    /// An `index` vector of full length `size` (as `HVector::setup` sizes it),
    /// with `index[0..count)` filled and zeros after: use it to build an
    /// input HVector for a replayed call.
    pub fn index_padded(&self) -> Vec<i32> {
        let mut idx = vec![0; self.size as usize];
        let src = self.index_slice();
        idx[..src.len()].copy_from_slice(src);
        idx
    }
}

/// `HighsSparseMatrix::format_` values (`MatrixFormat`).
pub const MATRIX_FORMAT_COLWISE: i32 = 1;
pub const MATRIX_FORMAT_ROWWISE: i32 = 2;
pub const MATRIX_FORMAT_ROWWISE_PARTITIONED: i32 = 3;

/// A `HighsSparseMatrix` (`Obj::mat` in emit.cpp).
#[derive(Clone, Debug, PartialEq)]
pub struct MatrixRec {
    /// `format`: 1 colwise, 2 rowwise, 3 rowwise partitioned.
    pub format: i32,
    pub num_col: i32,
    pub num_row: i32,
    /// `start`
    pub start: Vec<i32>,
    /// `p_end`, only for the partitioned format.
    pub p_end: Option<Vec<i32>>,
    pub index: Vec<i32>,
    pub value: Vec<f64>,
}

impl MatrixRec {
    pub fn from_rec(r: &Rec) -> Self {
        let format = r.int("format");
        let p_end = r.has("p_end").then(|| r.ints("p_end"));
        assert_eq!(
            p_end.is_some(),
            format == MATRIX_FORMAT_ROWWISE_PARTITIONED,
            "{}: p_end present iff format == 3",
            r.ctx
        );
        MatrixRec {
            format,
            num_col: r.int("num_col"),
            num_row: r.int("num_row"),
            start: r.ints("start"),
            p_end,
            index: r.ints("index"),
            value: r.f64s("value"),
        }
    }
}

/// A `SimplexBasis` (`basisObj` in simplex.cpp).
#[derive(Clone, Debug, PartialEq)]
pub struct BasisRec {
    /// `basicIndex`: `0..num_col` structurals, `num_col + i` logical of row i.
    pub basic_index: Vec<i32>,
    /// `nonbasicFlag`
    pub nonbasic_flag: Vec<i8>,
    /// `nonbasicMove`
    pub nonbasic_move: Vec<i8>,
}

impl BasisRec {
    pub fn from_rec(r: &Rec) -> Self {
        BasisRec {
            basic_index: r.ints("basicIndex"),
            nonbasic_flag: r.i8s("nonbasicFlag"),
            nonbasic_move: r.i8s("nonbasicMove"),
        }
    }
}

/// The JSON keys of an `arrays` object, in emission order (`workArrays` in
/// simplex.cpp). `hash` objects carry the same minus `workRange`, plus
/// `basicIndex` and `nonbasicMove`.
pub const WORK_ARRAY_KEYS: [&str; 11] = [
    "workCost",
    "workDual",
    "workShift",
    "workLower",
    "workUpper",
    "workRange",
    "workValue",
    "baseLower",
    "baseUpper",
    "baseValue",
    "dual_edge_weight",
];

/// The JSON keys of a `hash` object, in emission order (`hashObj`).
pub const HASH_KEYS: [&str; 12] = [
    "basicIndex",
    "nonbasicMove",
    "workCost",
    "workDual",
    "workShift",
    "workLower",
    "workUpper",
    "workValue",
    "baseValue",
    "baseLower",
    "baseUpper",
    "dual_edge_weight",
];

/// The 11 HEkk work/base arrays (`HighsSimplexInfo::work*_`, `base*_`,
/// `HEkk::dual_edge_weight_`).
#[derive(Clone, Debug, PartialEq)]
pub struct WorkArraysRec {
    pub work_cost: Vec<f64>,
    pub work_dual: Vec<f64>,
    pub work_shift: Vec<f64>,
    pub work_lower: Vec<f64>,
    pub work_upper: Vec<f64>,
    pub work_range: Vec<f64>,
    pub work_value: Vec<f64>,
    pub base_lower: Vec<f64>,
    pub base_upper: Vec<f64>,
    pub base_value: Vec<f64>,
    pub dual_edge_weight: Vec<f64>,
}

impl WorkArraysRec {
    pub fn from_rec(r: &Rec) -> Self {
        WorkArraysRec {
            work_cost: r.f64s("workCost"),
            work_dual: r.f64s("workDual"),
            work_shift: r.f64s("workShift"),
            work_lower: r.f64s("workLower"),
            work_upper: r.f64s("workUpper"),
            work_range: r.f64s("workRange"),
            work_value: r.f64s("workValue"),
            base_lower: r.f64s("baseLower"),
            base_upper: r.f64s("baseUpper"),
            base_value: r.f64s("baseValue"),
            dual_edge_weight: r.f64s("dual_edge_weight"),
        }
    }

    /// `(json key, values)` pairs in emission order.
    pub fn named(&self) -> [(&'static str, &[f64]); 11] {
        [
            ("workCost", &self.work_cost),
            ("workDual", &self.work_dual),
            ("workShift", &self.work_shift),
            ("workLower", &self.work_lower),
            ("workUpper", &self.work_upper),
            ("workRange", &self.work_range),
            ("workValue", &self.work_value),
            ("baseLower", &self.base_lower),
            ("baseUpper", &self.base_upper),
            ("baseValue", &self.base_value),
            ("dual_edge_weight", &self.dual_edge_weight),
        ]
    }
}

/// `HighsScale` inside `lp`.
#[derive(Clone, Debug, PartialEq)]
pub struct ScaleRec {
    pub strategy: i32,
    pub has_scaling: bool,
    pub num_col: i32,
    pub num_row: i32,
    pub cost: f64,
    pub col: Vec<f64>,
    pub row: Vec<f64>,
}

/// The scaled `HighsLp` HEkk receives (`lp` in `ekk_solve_enter`).
#[derive(Clone, Debug, PartialEq)]
pub struct LpRec {
    pub num_col: i32,
    pub num_row: i32,
    /// `sense`: 1 minimise, -1 maximise.
    pub sense: i32,
    pub offset: f64,
    pub col_cost: Vec<f64>,
    pub col_lower: Vec<f64>,
    pub col_upper: Vec<f64>,
    pub row_lower: Vec<f64>,
    pub row_upper: Vec<f64>,
    /// `a`: column-wise constraint matrix (scaled).
    pub a: MatrixRec,
    pub is_scaled: bool,
    pub is_moved: bool,
    pub scale: ScaleRec,
}

impl LpRec {
    pub fn from_rec(r: &Rec) -> Self {
        let s = r.obj("scale");
        LpRec {
            num_col: r.int("num_col"),
            num_row: r.int("num_row"),
            sense: r.int("sense"),
            offset: r.f64("offset"),
            col_cost: r.f64s("col_cost"),
            col_lower: r.f64s("col_lower"),
            col_upper: r.f64s("col_upper"),
            row_lower: r.f64s("row_lower"),
            row_upper: r.f64s("row_upper"),
            a: r.mat("a"),
            is_scaled: r.bool("is_scaled"),
            is_moved: r.bool("is_moved"),
            scale: ScaleRec {
                strategy: s.int("strategy"),
                has_scaling: s.bool("has_scaling"),
                num_col: s.int("num_col"),
                num_row: s.int("num_row"),
                cost: s.f64("cost"),
                col: s.f64s("col"),
                row: s.f64s("row"),
            },
        }
    }
}

/// Element type of a dumped array.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ArrKind {
    /// `vf`: doubles (hex bit patterns)
    F64,
    /// `vi`: HighsInt
    Int,
    /// `vi8`: int8_t
    I8,
}

/// The arrays of `InvertibleRepresentation` (`invertObj` in factor.cpp), in
/// emission order.
pub const INVERT_KEYS: [(&str, ArrKind); 25] = [
    ("l_pivot_index", ArrKind::Int),
    ("l_pivot_lookup", ArrKind::Int),
    ("l_start", ArrKind::Int),
    ("l_index", ArrKind::Int),
    ("l_value", ArrKind::F64),
    ("lr_start", ArrKind::Int),
    ("lr_index", ArrKind::Int),
    ("lr_value", ArrKind::F64),
    ("u_pivot_lookup", ArrKind::Int),
    ("u_pivot_index", ArrKind::Int),
    ("u_pivot_value", ArrKind::F64),
    ("u_start", ArrKind::Int),
    ("u_last_p", ArrKind::Int),
    ("u_index", ArrKind::Int),
    ("u_value", ArrKind::F64),
    ("ur_start", ArrKind::Int),
    ("ur_lastp", ArrKind::Int),
    ("ur_space", ArrKind::Int),
    ("ur_index", ArrKind::Int),
    ("ur_value", ArrKind::F64),
    ("pf_start", ArrKind::Int),
    ("pf_index", ArrKind::Int),
    ("pf_value", ArrKind::F64),
    ("pf_pivot_index", ArrKind::Int),
    ("pf_pivot_value", ArrKind::F64),
];

/// Every array of `HFactor::getInvert()` at full vector length, including
/// any reserved tail (FORMAT.md, `invert`).
#[derive(Clone, Debug, PartialEq)]
pub struct InvertRec {
    pub l_pivot_index: Vec<i32>,
    pub l_pivot_lookup: Vec<i32>,
    pub l_start: Vec<i32>,
    pub l_index: Vec<i32>,
    pub l_value: Vec<f64>,
    pub lr_start: Vec<i32>,
    pub lr_index: Vec<i32>,
    pub lr_value: Vec<f64>,
    pub u_pivot_lookup: Vec<i32>,
    pub u_pivot_index: Vec<i32>,
    pub u_pivot_value: Vec<f64>,
    pub u_start: Vec<i32>,
    pub u_last_p: Vec<i32>,
    pub u_index: Vec<i32>,
    pub u_value: Vec<f64>,
    pub ur_start: Vec<i32>,
    pub ur_lastp: Vec<i32>,
    pub ur_space: Vec<i32>,
    pub ur_index: Vec<i32>,
    pub ur_value: Vec<f64>,
    pub pf_start: Vec<i32>,
    pub pf_index: Vec<i32>,
    pub pf_value: Vec<f64>,
    pub pf_pivot_index: Vec<i32>,
    pub pf_pivot_value: Vec<f64>,
}

impl InvertRec {
    pub fn from_rec(r: &Rec) -> Self {
        InvertRec {
            l_pivot_index: r.ints("l_pivot_index"),
            l_pivot_lookup: r.ints("l_pivot_lookup"),
            l_start: r.ints("l_start"),
            l_index: r.ints("l_index"),
            l_value: r.f64s("l_value"),
            lr_start: r.ints("lr_start"),
            lr_index: r.ints("lr_index"),
            lr_value: r.f64s("lr_value"),
            u_pivot_lookup: r.ints("u_pivot_lookup"),
            u_pivot_index: r.ints("u_pivot_index"),
            u_pivot_value: r.f64s("u_pivot_value"),
            u_start: r.ints("u_start"),
            u_last_p: r.ints("u_last_p"),
            u_index: r.ints("u_index"),
            u_value: r.f64s("u_value"),
            ur_start: r.ints("ur_start"),
            ur_lastp: r.ints("ur_lastp"),
            ur_space: r.ints("ur_space"),
            ur_index: r.ints("ur_index"),
            ur_value: r.f64s("ur_value"),
            pf_start: r.ints("pf_start"),
            pf_index: r.ints("pf_index"),
            pf_value: r.f64s("pf_value"),
            pf_pivot_index: r.ints("pf_pivot_index"),
            pf_pivot_value: r.f64s("pf_pivot_value"),
        }
    }
}

/// `result` of a factor `build` / `refactor` record.
#[derive(Clone, Debug, PartialEq)]
pub struct BuildResultRec {
    /// return value of `HFactor::build()`
    pub rank_deficiency: i32,
    /// `basic_index` after build (HFactor permutes it in place)
    pub basic_index_after: Vec<i32>,
    pub row_with_no_pivot: Vec<i32>,
    pub col_with_no_pivot: Vec<i32>,
    pub var_with_no_pivot: Vec<i32>,
    /// `HFactor::build_synthetic_tick`
    pub build_synthetic_tick: f64,
    pub basis_matrix_num_el: i32,
    pub invert_num_el: i32,
    pub kernel_dim: i32,
    pub kernel_num_el: i32,
    /// `refactor_info_.pivot_row`
    pub refactor_pivot_row: Vec<i32>,
    /// `refactor_info_.pivot_var`
    pub refactor_pivot_var: Vec<i32>,
    /// `refactor_info_.pivot_type`
    pub refactor_pivot_type: Vec<i8>,
    /// `refactor_info_.build_synthetic_tick`
    pub refactor_build_synthetic_tick: f64,
    pub invert: InvertRec,
}

impl BuildResultRec {
    pub fn from_rec(r: &Rec) -> Self {
        BuildResultRec {
            rank_deficiency: r.int("rank_deficiency"),
            basic_index_after: r.ints("basic_index_after"),
            row_with_no_pivot: r.ints("row_with_no_pivot"),
            col_with_no_pivot: r.ints("col_with_no_pivot"),
            var_with_no_pivot: r.ints("var_with_no_pivot"),
            build_synthetic_tick: r.f64("build_synthetic_tick"),
            basis_matrix_num_el: r.int("basis_matrix_num_el"),
            invert_num_el: r.int("invert_num_el"),
            kernel_dim: r.int("kernel_dim"),
            kernel_num_el: r.int("kernel_num_el"),
            refactor_pivot_row: r.ints("refactor_pivot_row"),
            refactor_pivot_var: r.ints("refactor_pivot_var"),
            refactor_pivot_type: r.i8s("refactor_pivot_type"),
            refactor_build_synthetic_tick: r.f64("refactor_build_synthetic_tick"),
            invert: r.invert("invert"),
        }
    }
}

/// Input/output of one solve: an HVector (`ftran`/`btran`) or a plain
/// `std::vector<double>` (`ftran_vec`/`btran_vec`).
#[derive(Clone, Debug, PartialEq)]
pub enum SolveIo {
    Hv(HVectorRec),
    Vec(Vec<f64>),
}

impl SolveIo {
    pub fn hv(&self) -> &HVectorRec {
        match self {
            SolveIo::Hv(h) => h,
            SolveIo::Vec(_) => panic!("solve io is a plain vector, not an HVector"),
        }
    }
    pub fn vec(&self) -> &[f64] {
        match self {
            SolveIo::Vec(v) => v,
            SolveIo::Hv(_) => panic!("solve io is an HVector, not a plain vector"),
        }
    }
}

/// One entry of a factor record's `solves` list.
#[derive(Clone, Debug, PartialEq)]
pub struct SolveRec {
    /// Location of this solve, for messages.
    pub ctx: String,
    /// `ftran` / `btran` = `ftranCall` / `btranCall(HVector&, expected_density)`;
    /// `ftran_vec` / `btran_vec` = the `std::vector<double>&` overloads.
    pub op: String,
    /// `unit` / `column` / `sparse` / `dense`.
    pub rhs_kind: String,
    /// HVector ops only.
    pub expected_density: Option<f64>,
    pub input: SolveIo,
    pub out: SolveIo,
}

impl SolveRec {
    pub fn from_rec(r: &Rec) -> Self {
        let op = r.str("op").to_string();
        let (input, out) = if op.ends_with("_vec") {
            (SolveIo::Vec(r.f64s("in")), SolveIo::Vec(r.f64s("out")))
        } else {
            (SolveIo::Hv(r.hv("in")), SolveIo::Hv(r.hv("out")))
        };
        SolveRec {
            ctx: r.ctx.clone(),
            op,
            rhs_kind: r.str("rhs_kind").to_string(),
            expected_density: r.opt_f64("expected_density"),
            input,
            out,
        }
    }

    pub fn is_ftran(&self) -> bool {
        self.op.starts_with("ftran")
    }
}

/// One `random` record of `random.jsonl` (HighsRandom draws in this order).
#[derive(Clone, Debug, PartialEq)]
pub struct RandomRec {
    /// `seed` (HighsUInt)
    pub seed: u32,
    /// 16 x `integer()`
    pub integer: Vec<i32>,
    /// the `sup` argument of each `integer(sup)` call
    pub sup: Vec<i32>,
    /// 64 x `integer(sup)`
    pub integer_sup: Vec<i32>,
    /// 16 x `fraction()`
    pub fraction: Vec<f64>,
    /// 16 x `closedFraction()`
    pub closed_fraction: Vec<f64>,
    /// 16 x `real(-3, 5)`
    pub real_m3_5: Vec<f64>,
    /// 16 x `bit()`
    pub bit: Vec<i32>,
    /// `shuffle` of `0..20`
    pub shuffle20: Vec<i32>,
}

impl RandomRec {
    pub fn from_rec(r: &Rec) -> Self {
        let seed = r.i64("seed");
        RandomRec {
            seed: u32::try_from(seed)
                .unwrap_or_else(|_| panic!("{}: seed {seed} not a u32", r.ctx)),
            integer: r.ints("integer"),
            sup: r.ints("sup"),
            integer_sup: r.ints("integer_sup"),
            fraction: r.f64s("fraction"),
            closed_fraction: r.f64s("closedFraction"),
            real_m3_5: r.f64s("real_m3_5"),
            bit: r.ints("bit"),
            shuffle20: r.ints("shuffle20"),
        }
    }
}
