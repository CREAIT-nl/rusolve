//! HiGHS source: util/HighsUtils.h, util/HighsUtils.cpp (subset)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Index collections plus the small numeric helpers used on the simplex path.
//! `analyseVectorValues`, `analyseMatrixSparsity`, value distributions and scatter data are
//! analysis-only (out of scope): `HighsValueDistribution` is kept as plain data so
//! `HighsSimplexAnalysis` can hold it.

use crate::prelude::*;

/// HiGHS: HighsUtils.h::HighsIndexCollection
#[derive(Clone, Debug, PartialEq)]
pub struct HighsIndexCollection {
    pub dimension_: HighsInt,
    pub is_interval_: bool,
    pub from_: HighsInt,
    pub to_: HighsInt,
    pub is_set_: bool,
    pub set_num_entries_: HighsInt,
    pub set_: Vec<HighsInt>,
    pub is_mask_: bool,
    pub mask_: Vec<HighsInt>,
}

impl Default for HighsIndexCollection {
    fn default() -> Self {
        Self {
            dimension_: -1,
            is_interval_: false,
            from_: -1,
            to_: -2,
            is_set_: false,
            set_num_entries_: -1,
            set_: Default::default(),
            is_mask_: false,
            mask_: Default::default(),
        }
    }
}

/// HiGHS: HighsUtils.h::HighsValueDistribution
#[derive(Clone, Debug, Default)]
pub struct HighsValueDistribution {
    pub distribution_name_: String,
    pub value_name_: String,
    pub num_count_: HighsInt,
    pub num_zero_: HighsInt,
    pub num_one_: HighsInt,
    pub min_value_: f64,
    pub max_value_: f64,
    pub limit_: Vec<f64>,
    pub count_: Vec<HighsInt>,
    pub sum_count_: HighsInt,
}

/// HiGHS: HighsUtils.cpp::create (lines 21-30, 10 lines)
/// Overload `create(index_collection, from_col, to_col, dimension)` (interval).
pub fn create_interval(
    index_collection: &mut HighsIndexCollection,
    from_col: HighsInt,
    to_col: HighsInt,
    dimension: HighsInt,
) -> HighsInt {
    // HighsUtils.h: kIndexCollectionCreateOk = 0, kIndexCollectionCreateIllegalInterval = 1
    const kIndexCollectionCreateOk: HighsInt = 0;
    const kIndexCollectionCreateIllegalInterval: HighsInt = 1;
    if from_col < 0 {
        return kIndexCollectionCreateIllegalInterval;
    }
    if to_col >= dimension {
        return kIndexCollectionCreateIllegalInterval;
    }
    index_collection.dimension_ = dimension;
    index_collection.is_interval_ = true;
    index_collection.from_ = from_col;
    index_collection.to_ = to_col;
    kIndexCollectionCreateOk
}

/// HiGHS: HighsUtils.cpp::create (lines 32-48, 17 lines)
/// Overload `create(index_collection, num_set_entries, set, dimension)` (set).
pub fn create_set(
    index_collection: &mut HighsIndexCollection,
    num_set_entries: HighsInt,
    set: &[HighsInt],
    dimension: HighsInt,
) -> HighsInt {
    // HighsUtils.h: kIndexCollectionCreateOk = 0, kIndexCollectionCreateIllegalSetSize = 1,
    // kIndexCollectionCreateIllegalSetDimension = 2, kIndexCollectionCreateIllegalSetOrder = 3
    const kIndexCollectionCreateOk: HighsInt = 0;
    const kIndexCollectionCreateIllegalSetSize: HighsInt = 1;
    const kIndexCollectionCreateIllegalSetDimension: HighsInt = 2;
    const kIndexCollectionCreateIllegalSetOrder: HighsInt = 3;
    // Create an index collection for the given set - so long as it is strictly
    // ordered
    if num_set_entries < 0 {
        return kIndexCollectionCreateIllegalSetSize;
    }
    if dimension < 0 {
        return kIndexCollectionCreateIllegalSetDimension;
    }
    index_collection.dimension_ = dimension;
    index_collection.is_set_ = true;
    index_collection.set_ = set[..num_set_entries as usize].to_vec();
    index_collection.set_num_entries_ = num_set_entries;
    // HighsSort.cpp::increasingSetOk(const vector<HighsInt>&, HighsInt, HighsInt, bool)
    // (lines 261-287), inlined: the integer overload is not part of this port's HighsSort.
    let increasing_set_ok = |set: &[HighsInt],
                             set_entry_lower: HighsInt,
                             set_entry_upper: HighsInt,
                             strict: bool|
     -> bool {
        let set_num_entries = set.len() as HighsInt;
        let check_bounds = set_entry_lower <= set_entry_upper;
        let mut previous_entry = if check_bounds {
            if strict {
                set_entry_lower - 1
            } else {
                set_entry_lower
            }
        } else {
            -kHighsIInf
        };
        for k in 0..set_num_entries {
            let entry = set[k as usize];
            if strict {
                if entry <= previous_entry {
                    return false;
                }
            } else {
                if entry < previous_entry {
                    return false;
                }
            }
            if check_bounds && entry > set_entry_upper {
                return false;
            }
            previous_entry = entry;
        }
        true
    };
    if !increasing_set_ok(&index_collection.set_, 1, 0, true) {
        return kIndexCollectionCreateIllegalSetOrder;
    }
    for ix in 0..num_set_entries {
        if set[ix as usize] < 0 || set[ix as usize] >= dimension {
            return -(ix + 1);
        }
    }
    kIndexCollectionCreateOk
}

/// HiGHS: HighsUtils.cpp::create (lines 50-58, 9 lines)
/// Overload `create(index_collection, mask, dimension)` (mask).
pub fn create_mask(
    index_collection: &mut HighsIndexCollection,
    mask: &[HighsInt],
    dimension: HighsInt,
) -> HighsInt {
    // HighsUtils.h: kIndexCollectionCreateOk = 0, kIndexCollectionCreateIllegalMaskSize = 1
    const kIndexCollectionCreateOk: HighsInt = 0;
    const kIndexCollectionCreateIllegalMaskSize: HighsInt = 1;
    // Create an index collection for the given mask
    if dimension < 0 {
        return kIndexCollectionCreateIllegalMaskSize;
    }
    index_collection.dimension_ = dimension;
    index_collection.is_mask_ = true;
    index_collection.mask_ = mask[..dimension as usize].to_vec();
    kIndexCollectionCreateOk
}

/// HiGHS: HighsUtils.cpp::ok (lines 90-174, 85 lines)
pub fn ok(index_collection: &HighsIndexCollection) -> bool {
    // HighsSort.cpp::increasingSetOk(const vector<HighsInt>&, HighsInt, HighsInt, bool)
    // (lines 261-287), inlined; used only by the debug assert below.
    let increasing_set_ok = |set: &[HighsInt],
                             set_entry_lower: HighsInt,
                             set_entry_upper: HighsInt,
                             strict: bool|
     -> bool {
        let set_num_entries = set.len() as HighsInt;
        let check_bounds = set_entry_lower <= set_entry_upper;
        let mut previous_entry = if check_bounds {
            if strict {
                set_entry_lower - 1
            } else {
                set_entry_lower
            }
        } else {
            -kHighsIInf
        };
        for k in 0..set_num_entries {
            let entry = set[k as usize];
            if strict {
                if entry <= previous_entry {
                    return false;
                }
            } else {
                if entry < previous_entry {
                    return false;
                }
            }
            if check_bounds && entry > set_entry_upper {
                return false;
            }
            previous_entry = entry;
        }
        true
    };
    // Check parameter for each technique of defining an index collection
    if index_collection.is_interval_ {
        // Changing by interval: check the parameters and that check set and mask
        // are false
        if index_collection.is_set_ {
            println!("Index collection is both interval and set");
            return false;
        }
        if index_collection.is_mask_ {
            println!("Index collection is both interval and mask");
            return false;
        }
        if index_collection.from_ < 0 {
            println!(
                "Index interval lower limit is {} < 0",
                index_collection.from_
            );
            return false;
        }
        if index_collection.to_ > index_collection.dimension_ - 1 {
            println!(
                "Index interval upper limit is {} > {}",
                index_collection.to_,
                index_collection.dimension_ - 1
            );
            return false;
        }
    } else if index_collection.is_set_ {
        // Changing by set: check the parameters and check that interval and mask
        // are false
        if index_collection.is_interval_ {
            println!("Index collection is both set and interval");
            return false;
        }
        if index_collection.is_mask_ {
            println!("Index collection is both set and mask");
            return false;
        }
        if index_collection.set_.is_empty() {
            println!("Index set is NULL");
            return false;
        }
        // Check that the values in the vector of integers are ascending
        let set = &index_collection.set_;
        let num_entries = index_collection.set_num_entries_;
        let entry_upper = index_collection.dimension_ - 1;
        let mut prev_set_entry: HighsInt = -1;
        for k in 0..num_entries {
            let kk = k as usize;
            if set[kk] < 0 || set[kk] > entry_upper {
                println!(
                    "Index set entry set[{}] = {} is out of bounds [0, {}]",
                    k, set[kk], entry_upper
                );
                return false;
            }
            if set[kk] <= prev_set_entry {
                println!(
                    "Index set entry set[{}] = {} is not greater than previous entry {}",
                    k, set[kk], prev_set_entry
                );
                return false;
            }
            prev_set_entry = set[kk];
        }
        // This was the old check done independently, and should be
        // equivalent.
        debug_assert!(increasing_set_ok(set, 0, entry_upper, true));
    } else if index_collection.is_mask_ {
        // Changing by mask: check the parameters and check that set and interval
        // are false
        if index_collection.mask_.is_empty() {
            println!("Index mask is NULL");
            return false;
        }
        if index_collection.is_interval_ {
            println!("Index collection is both mask and interval");
            return false;
        }
        if index_collection.is_set_ {
            println!("Index collection is both mask and set");
            return false;
        }
    } else {
        // No method defined
        println!("Undefined index collection");
        return false;
    }
    true
}

/// HiGHS: HighsUtils.cpp::limits (lines 176-190, 15 lines)
pub fn limits(index_collection: &HighsIndexCollection, from_k: &mut HighsInt, to_k: &mut HighsInt) {
    if index_collection.is_interval_ {
        *from_k = index_collection.from_;
        *to_k = index_collection.to_;
    } else if index_collection.is_set_ {
        *from_k = 0;
        *to_k = index_collection.set_num_entries_ - 1;
    } else if index_collection.is_mask_ {
        *from_k = 0;
        *to_k = index_collection.dimension_ - 1;
    } else {
        // C++ assert(1 == 0): active only in debug builds.
        if cfg!(debug_assertions) {
            panic!("assert(1 == 0) failed");
        }
    }
}

/// HiGHS: HighsUtils.cpp::updateOutInIndex (lines 192-238, 47 lines)
pub fn update_out_in_index(
    index_collection: &HighsIndexCollection,
    out_from_ix: &mut HighsInt,
    out_to_ix: &mut HighsInt,
    in_from_ix: &mut HighsInt,
    in_to_ix: &mut HighsInt,
    current_set_entry: &mut HighsInt,
) {
    if index_collection.is_interval_ {
        *out_from_ix = index_collection.from_;
        *out_to_ix = index_collection.to_;
        *in_from_ix = index_collection.to_ + 1;
        *in_to_ix = index_collection.dimension_ - 1;
    } else if index_collection.is_set_ {
        *out_from_ix = index_collection.set_[*current_set_entry as usize];
        *out_to_ix = *out_from_ix;
        *current_set_entry += 1;
        let current_set_entry0 = *current_set_entry;
        for set_entry in current_set_entry0..index_collection.set_num_entries_ {
            let ix = index_collection.set_[set_entry as usize];
            if ix > *out_to_ix + 1 {
                break;
            }
            *out_to_ix = index_collection.set_[*current_set_entry as usize];
            *current_set_entry += 1;
        }
        *in_from_ix = *out_to_ix + 1;
        if *current_set_entry < index_collection.set_num_entries_ {
            *in_to_ix = index_collection.set_[*current_set_entry as usize] - 1;
        } else {
            // Account for getting to the end of the set
            *in_to_ix = index_collection.dimension_ - 1;
        }
    } else {
        *out_from_ix = *in_to_ix + 1;
        *out_to_ix = index_collection.dimension_ - 1;
        for ix in (*in_to_ix + 1)..index_collection.dimension_ {
            if index_collection.mask_[ix as usize] == 0 {
                *out_to_ix = ix - 1;
                break;
            }
        }
        *in_from_ix = *out_to_ix + 1;
        *in_to_ix = index_collection.dimension_ - 1;
        for ix in (*out_to_ix + 1)..index_collection.dimension_ {
            if index_collection.mask_[ix as usize] != 0 {
                *in_to_ix = ix - 1;
                break;
            }
        }
    }
}

/// HiGHS: HighsUtils.cpp::dataSize (lines 240-250, 11 lines)
pub fn data_size(index_collection: &HighsIndexCollection) -> HighsInt {
    if index_collection.is_set_ {
        index_collection.set_num_entries_
    } else {
        if index_collection.is_interval_ {
            index_collection.to_ - index_collection.from_ + 1
        } else {
            index_collection.dimension_
        }
    }
}

/// HiGHS: HighsUtils.cpp::getNorm2 (lines 289-294, 6 lines)
pub fn get_norm2(values: &[f64]) -> f64 {
    let mut sum: f64 = 0.0;
    let values_size = values.len() as HighsInt;
    for i in 0..values_size {
        sum += values[i as usize] * values[i as usize];
    }
    sum
}

/// HiGHS: HighsUtils.cpp::highs_isInfinity (lines 296-299, 4 lines)
pub fn highs_is_infinity(val: f64) -> bool {
    val >= kHighsInf
}

/// HiGHS: HighsUtils.cpp::highsRelativeDifference (lines 301-303, 3 lines)
pub fn highs_relative_difference(v0: f64, v1: f64) -> f64 {
    (v0 - v1).abs() / cpp_max(v0, cpp_max(v1, 1.0))
}

/// HiGHS: HighsUtils.cpp::nearestPowerOfTwoScale (lines 1211-1232, 22 lines)
pub fn nearest_power_of_two_scale(value: f64) -> f64 {
    // std::frexp(value, &exp_scale), inlined (Rust std has no frexp). Returns x in
    // [0.5, 1) (sign kept) with value == x * 2^exp; zero, inf and nan are returned as is
    // with exp 0. Subnormals are first scaled by 2^64 (exact) and the exponent corrected.
    let subnormal = value != 0.0 && value.is_finite() && value.abs() < f64::MIN_POSITIVE;
    let (frexp_src, subnormal_adj): (f64, HighsInt) = if subnormal {
        (value * 2f64.powf(64.0), 64)
    } else {
        (value, 0)
    };
    let mut exp_scale: HighsInt;
    let mut check_x: f64;
    if frexp_src == 0.0 || !frexp_src.is_finite() {
        check_x = value;
        exp_scale = 0;
    } else {
        let bits = frexp_src.to_bits();
        let biased = ((bits >> 52) & 0x7ff) as HighsInt;
        exp_scale = biased - 1022 - subnormal_adj;
        check_x = f64::from_bits((bits & !(0x7ffu64 << 52)) | (1022u64 << 52));
    }
    // Decompose value into a normalized fraction and an integral power
    // of two.
    //
    // If arg is zero, returns zero and stores zero in *exp. Otherwise
    // (if arg is not zero), if no errors occur, returns the value x in
    // the range (-1;-0.5], [0.5; 1) and stores an integer value in *exp
    // such that x×2(*exp)=arg
    if check_x.abs() == 0.5 {
        check_x *= 2.0;
        exp_scale -= 1;
    }
    let check_value = check_x * 2f64.powf(exp_scale as f64);
    debug_assert!(check_value == value);
    exp_scale = -exp_scale;
    // Multiply a floating point value x(=1) by the number 2 raised to
    // the exp power
    // std::ldexp(1, exp_scale), as 2^exp_scale (exact for every representable power of two).
    2f64.powf(exp_scale as f64)
}

/// HiGHS: HighsUtils.cpp::highsAssert (lines 1234-1246, 13 lines)
/// C++ default arguments: message = ""
/// In C++ this logs and calls `assert`; translate as a `panic!` with the message when `!assertion`.
pub fn highs_assert(assert_condition: bool, message: &str) {
    if assert_condition {
        return;
    }
    println!("Failing highsAssert(\"{}\")", message);
    panic!("assert({}) failed", message);
}

/// HiGHS: HighsUtils.cpp::analyseVectorValues(const HighsLogOptions*, const std::string&,
/// HighsInt, const std::vector<double>&, bool, std::string)
/// C++ default arguments: analyse_value_list = false, model_name = "Unknown"
/// Value-distribution report (called only from `HSimplexNla::report*`): out of scope;
/// stub (no-op). `log_options = None` is the C++ `nullptr` (printf) case.
pub fn analyse_vector_values(
    log_options: Option<&HighsLogOptions>,
    message: &str,
    vec_dim: HighsInt,
    vec: &[f64],
    analyse_value_list: bool,
    model_name: &str,
) {
}

/// C++ `std::max(a, b)`: `(a < b) ? b : a`.
///
/// Translate every `std::max`/`max` with this, never `f64::max`: they differ
/// for NaN and for signed zeros (`std::max(-0.0, 0.0)` is `-0.0`), and either
/// difference can change a pivot.
pub fn cpp_max<T: PartialOrd>(a: T, b: T) -> T {
    if a < b { b } else { a }
}

/// C++ `std::min(a, b)`: `(b < a) ? b : a`. See [`cpp_max`]; never `f64::min`.
pub fn cpp_min<T: PartialOrd>(a: T, b: T) -> T {
    if b < a { b } else { a }
}

/// HiGHS: HighsUtils.h::fractionality (lines 210-217)
/// C++: `template <typename T> T fractionality(T input, T* intval = nullptr)` (`T = double`
/// at every MIP call site). `std::round` rounds half away from zero, as `f64::round` does.
pub fn fractionality(input: f64, intval: Option<&mut f64>) -> f64 {
    let val = input.round();
    if let Some(intval) = intval {
        *intval = val;
    }
    (input - val).abs()
}
