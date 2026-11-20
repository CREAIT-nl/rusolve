;; Example component (ABI v1, see ../src/lib.rs): branch on the most fractional candidate and
;; rank pool cuts by raw violation. Hand-written WAT because the wasm32-unknown-unknown target is
;; not installed; a Rust `#![no_std]` cdylib exporting the same two functions is equivalent.
(module
  (memory (export "memory") 1)
  ;; branch_score(depth, n, in_ptr, out_ptr) -> 0
  (func (export "branch_score") (param $depth i32) (param $n i32) (param $in i32) (param $out i32) (result i32)
    (local $i i32) (local $v f64) (local $frac f64)
    (block $done
      (loop $next
        (br_if $done (i32.ge_u (local.get $i) (local.get $n)))
        (local.set $v
          (f64.load offset=8 (i32.add (local.get $in) (i32.shl (local.get $i) (i32.const 5)))))
        (local.set $frac (f64.sub (local.get $v) (f64.floor (local.get $v))))
        (f64.store
          (i32.add (local.get $out) (i32.shl (local.get $i) (i32.const 3)))
          (f64.sub (f64.const 0.5) (f64.abs (f64.sub (local.get $frac) (f64.const 0.5)))))
        (local.set $i (i32.add (local.get $i) (i32.const 1)))
        (br $next)))
    (i32.const 0))
  ;; cut_score(violation, nnz, active_nnz, active_norm, default_score) -> score
  (func (export "cut_score") (param f64 i32 i32 f64 f64) (result f64)
    (local.get 0)))
