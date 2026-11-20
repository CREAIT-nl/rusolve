;; Example user separator (ABI v1 `separate`, see ../src/lib.rs): the cover inequality
;; x0 + x1 + x2 <= 2 of a knapsack row in which items 0, 1, 2 together do not fit, emitted when the
;; node LP solution violates it. Result: one record `nnz, rhs, (col, val)*nnz` = 3, 2, 0,1, 1,1, 2,1
;; (8 f64), returned as the number of f64 written; 0 = no cut. Hand-written WAT, as the others.
(module
  (memory (export "memory") 1)
  ;; separate(n, in_ptr, out_ptr, out_cap) -> f64 slots written, -1 = no opinion
  (func (export "separate") (param $n i32) (param $in i32) (param $out i32) (param $cap i32) (result i32)
    (if (i32.lt_u (local.get $n) (i32.const 3)) (then (return (i32.const -1))))
    (if (i32.lt_u (local.get $cap) (i32.const 8)) (then (return (i32.const 0))))
    (if (f64.le
          (f64.add
            (f64.add (f64.load (local.get $in))
                     (f64.load offset=8 (local.get $in)))
            (f64.load offset=16 (local.get $in)))
          (f64.const 2.0001))
      (then (return (i32.const 0))))
    (f64.store (local.get $out) (f64.const 3))
    (f64.store offset=8 (local.get $out) (f64.const 2))
    (f64.store offset=16 (local.get $out) (f64.const 0))
    (f64.store offset=24 (local.get $out) (f64.const 1))
    (f64.store offset=32 (local.get $out) (f64.const 1))
    (f64.store offset=40 (local.get $out) (f64.const 1))
    (f64.store offset=48 (local.get $out) (f64.const 2))
    (f64.store offset=56 (local.get $out) (f64.const 1))
    (i32.const 8)))
