;; Example component (ABI v1, see ../src/lib.rs): node selection depth-first, ties by the
;; lowest lower bound (a dive-first rule many MIP solvers mix into best-bound). Records at
;; in_ptr are four f64 each: lower_bound, estimate, depth, num_domain_changes. Hand-written WAT
;; because the wasm32-unknown-unknown target is not installed; a Rust `#![no_std]` cdylib
;; exporting the same function is equivalent.
(module
  (memory (export "memory") 1)
  ;; select_node(n, in_ptr) -> chosen index in [0, n), or -1 for "no opinion"
  (func (export "select_node") (param $n i32) (param $in i32) (result i32)
    (local $i i32) (local $best i32)
    (local $best_depth f64) (local $best_lb f64)
    (local $depth f64) (local $lb f64)
    (local.set $best (i32.const -1))
    (local.set $best_depth (f64.const -inf))
    (local.set $best_lb (f64.const inf))
    (block $done
      (loop $next
        (br_if $done (i32.ge_u (local.get $i) (local.get $n)))
        ;; record base: in_ptr + 32*i; lower_bound at +0, depth at +16
        (local.set $lb (f64.load (i32.add (local.get $in) (i32.shl (local.get $i) (i32.const 5)))))
        (local.set $depth
          (f64.load offset=16 (i32.add (local.get $in) (i32.shl (local.get $i) (i32.const 5)))))
        (if
          (i32.or
            (f64.gt (local.get $depth) (local.get $best_depth))
            (i32.and
              (f64.eq (local.get $depth) (local.get $best_depth))
              (f64.lt (local.get $lb) (local.get $best_lb))))
          (then
            (local.set $best (local.get $i))
            (local.set $best_depth (local.get $depth))
            (local.set $best_lb (local.get $lb))))
        (local.set $i (i32.add (local.get $i) (i32.const 1)))
        (br $next)))
    (local.get $best)))
