"""One highspy MIP solve for tools/bench/abtest.py (arm kind `highspy`).

    python tools/bench/abtest_highs.py <model.mps[.gz]> <time_limit> <seed> [<node_limit>] [k=v ...]

Needs highspy (`pip install highspy`). Prints the CSV header and row that
`rusolve mip --csv` prints (without x), so abtest.py parses both solvers the
same way. Options: threads=1, the given random_seed and time_limit, HiGHS
default gaps, plus any trailing k=v pairs.
"""

import sys
import time

import highspy


def main():
    path, tl, seed = sys.argv[1], float(sys.argv[2]), int(sys.argv[3])
    rest = sys.argv[4:]
    node_limit = None
    if rest and "=" not in rest[0]:
        node_limit = int(rest[0])
        rest = rest[1:]
    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    h.setOptionValue("threads", 1)
    h.setOptionValue("random_seed", seed)
    h.setOptionValue("time_limit", tl)
    if node_limit is not None:
        h.setOptionValue("mip_max_nodes", node_limit)
    for kv in rest:
        k, v = kv.split("=", 1)
        for cast in (int, float):
            try:
                v = cast(v)
                break
            except ValueError:
                pass
        if v in ("true", "false"):
            v = v == "true"
        h.setOptionValue(k, v)
    h.readModel(path)
    t0 = time.perf_counter()
    h.run()
    secs = time.perf_counter() - t0
    st = h.modelStatusToString(h.getModelStatus())
    status = {"Optimal": "Optimal", "Infeasible": "Infeasible", "Time limit reached": "TimeLimit",
              "Solution limit reached": "NodeLimit"}.get(st, st)
    info = h.getInfo()
    obj = info.objective_function_value if info.primal_solution_status else float("nan")
    print("status,objective,bound,nodes,lp_iters,time_s,x")
    print(f"{status},{obj!r},{info.mip_dual_bound!r},{info.mip_node_count},"
          f"{info.simplex_iteration_count},{secs},[]")


if __name__ == "__main__":
    main()
