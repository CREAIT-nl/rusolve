;; Selfcheck fixture for the extend_wasm eval (evals/graders.py), NOT part of the task prompt.
;;
;; The reference rule: depth-first node selection. Of the open nodes the host offers, dive into
;; the deepest one, breaking ties by the better local estimate. Depth-first (dive-first) node
;; choice instead of pure best-bound is textbook branch-and-bound (Wolsey, "Mixed Integer
;; Programming", Wiley 2020, ch. 11 "Branch-and-bound: depth-first search"; most MIP solvers
;; mix a diving phase into best-bound for exactly this reason). On this instance it reaches the
;; optimum in 534 nodes against 912 for the best-bound pop under the HiGHS child rule
;; (mip_child_selection=-1, the setting the grader uses): a dive finds a good
;; incumbent early, which prunes the rest of the tree.
;;
;; ABI v1 `select_node(n, in_ptr) -> i32`: the host writes n records of four f64 at in_ptr
;; (lower_bound +0, estimate +8, depth +16, num_domain_changes +24, 32 bytes per record) and the
;; component returns the index in [0, n) of the node to explore next, or -1 for "no opinion".
(module
  (memory (export "memory") 1)
  (func (export "select_node") (param $n i32) (param $in i32) (result i32)
    (local $i i32) (local $best i32)
    (local $best_depth f64) (local $best_est f64)
    (local $depth f64) (local $est f64)
    (local.set $best (i32.const -1))
    (local.set $best_depth (f64.const -inf))
    (local.set $best_est (f64.const inf))
    (block $done
      (loop $next
        (br_if $done (i32.ge_u (local.get $i) (local.get $n)))
        ;; record i starts at in_ptr + 32*i: estimate at +8, depth at +16
        (local.set $est
          (f64.load offset=8 (i32.add (local.get $in) (i32.shl (local.get $i) (i32.const 5)))))
        (local.set $depth
          (f64.load offset=16 (i32.add (local.get $in) (i32.shl (local.get $i) (i32.const 5)))))
        (if
          (i32.or
            (f64.gt (local.get $depth) (local.get $best_depth))
            (i32.and
              (f64.eq (local.get $depth) (local.get $best_depth))
              (f64.lt (local.get $est) (local.get $best_est))))
          (then
            (local.set $best (local.get $i))
            (local.set $best_depth (local.get $depth))
            (local.set $best_est (local.get $est))))
        (local.set $i (i32.add (local.get $i) (i32.const 1)))
        (br $next)))
    (local.get $best)))
