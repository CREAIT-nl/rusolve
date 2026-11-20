//! HiGHS source: util/ (subset)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! HiGHS `util/` subset: sparse vectors/matrices, LU factorization, small helpers,
//! and the presolve leaves (data stack, hash table, linear sum bounds, matrix
//! slices, pdqsort, splay), and the MIP leaves (disjoint sets, hash tree,
//! integers, red-black tree) and the libc++ algorithms they depend on
//! (cpp_std).

pub mod c_printf;
pub mod cpp_std;
pub mod factor_timer;
pub mod hfactor;
pub mod hfactor_const;
pub mod highs_cdouble;
pub mod highs_data_stack;
pub mod highs_disjoint_sets;
pub mod highs_hash;
pub mod highs_hash_table;
pub mod highs_hash_tree;
pub mod highs_integers;
pub mod highs_linear_sum_bounds;
pub mod highs_matrix_slice;
pub mod highs_matrix_utils;
pub mod highs_pdqsort;
pub mod highs_random;
pub mod highs_rb_tree;
pub mod highs_sort;
pub mod highs_sort_set;
pub mod highs_sparse_vector_sum;
pub mod highs_splay;
pub mod highs_timer;
pub mod highs_utils;
pub mod hset;
pub mod hvector;
pub mod sparse_matrix;
