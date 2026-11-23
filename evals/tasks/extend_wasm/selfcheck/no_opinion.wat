;; Selfcheck fixture for the extend_wasm eval (evals/graders.py), NOT part of the task prompt.
;;
;; The degenerate component: it exports `select_node` and always answers -1, "no opinion". Every
;; call falls back to the HiGHS best-bound pop, so the tree is exactly the default one. The host
;; counts a -1 answer as `declined`, never `ok`, so this component fails the grader twice: `ok`
;; is 0 (the rule made no decision) and the node count is the default one.
(module
  (memory (export "memory") 1)
  (func (export "select_node") (param $n i32) (param $in i32) (result i32)
    (i32.const -1)))
