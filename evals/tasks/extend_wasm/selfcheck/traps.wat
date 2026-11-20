;; Selfcheck fixture for the extend_wasm eval (evals/graders.py), NOT part of the task prompt.
;;
;; A component that traps on every call. The sandbox survives (the solve still ends Optimal, at
;; the default node count), so the grader must fail it on its counters, `traps > 0` and `ok == 0`,
;; not on the objective.
(module
  (memory (export "memory") 1)
  (func (export "select_node") (param $n i32) (param $in i32) (result i32)
    (unreachable)))
