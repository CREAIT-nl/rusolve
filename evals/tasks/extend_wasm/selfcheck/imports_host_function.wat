;; Selfcheck fixture for the extend_wasm eval (evals/graders.py), NOT part of the task prompt.
;;
;; A component that cannot be loaded at all: the sandbox has no host functions, so a module that
;; imports one is rejected before the solve (the CLI exits 3 naming the file, the MCP `solve`
;; tool answers a tool error naming the file).
(module
  (import "env" "now" (func $now (result i64)))
  (memory (export "memory") 1)
  (func (export "select_node") (param $n i32) (param $in i32) (result i32)
    (drop (call $now))
    (i32.const 0)))
