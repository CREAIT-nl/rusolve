#!/usr/bin/env python3
"""Runs the agent-task evals: a real agent drives rusolve, a deterministic grader scores it.

    python3 evals/run.py --tasks solve_report,debug_infeasible --model haiku
    python3 evals/run.py --all --model haiku --mode cli
    python3 evals/run.py --mode python --model haiku           # word problems -> rusolve.Model code
    python3 evals/run.py --mode python --model sonnet          # a second model: any `claude --model` value
    HIGHS_MCP_CMD="npx highs-mcp" python3 evals/run.py --mode highs-mcp   # same word problems, highs-mcp as the solver

The agent is Claude Code headless (`claude -p`). `--mode mcp` (default) gives it only the
`rusolve-cli mcp` server (docs/design/mcp.md); `--mode cli` gives it only Bash on the rusolve-cli
binary; `--mode python` gives it Bash on a Python interpreter that has the `rusolve` package
(`RSV_PY`, default `python3`) plus Write, Edit and Read, for the word-problem tasks whose `modes` list `python`
(it writes `rusolve.Model` code); `--mode highs-mcp` runs the same word problems against the
`highs-mcp` server started by `$HIGHS_MCP_CMD` (for comparison; the server is not part of this repo).
For each task it records success (graders.py, ground truth from rusolve itself), the
number of tool calls, tool errors, turns, tokens and cost, and writes
`evals/results/<UTC stamp>-<model>-<mode>.json`. One run per task per invocation: repeat the
command for a sample, there is deliberately no sweep driver.
"""
import argparse
import datetime
import json
import os
import shlex
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "lib"))
import graders  # noqa: E402
import rusolve as rs  # noqa: E402

TASKS = ["solve_report", "debug_infeasible", "verify_formulation", "tune_family", "add_cuts", "extend_wasm",
         "class_router", "class_flow", "class_cover", "class_scheduling", "class_lotsizing", "class_facility",
         "class_routing", "class_pwl_semi", "class_qp",
         "diet_plan", "project_selection", "fixed_charge_supply", "shift_cover", "mcp_build_model"]
MODES = ["mcp", "cli", "python", "highs-mcp"]

MCP_NOTE = ("\n\nYou are connected to the rusolve solver through MCP tools (load_model, solve, get_solution, "
            "build_model, modify, check_solution, explain, lint, ranging, feasrelax, free_model). Use them; do not "
            "estimate by hand what the solver can compute.")
CLI_NOTE = ("\n\nUse the rusolve command line through Bash: the binary is {bin}. `{bin} mip <model.mps> "
            "--format json --solution` solves (set RUSOLVE_MIP_HIGHS=1 in the environment; `--opt key=value` "
            "sets a solver option); `check`, `explain`, `feasrelax`, `lint`, `ranging` and `options --json` are "
            "subcommands. Do not estimate by hand what the solver can compute.")
PY_NOTE = ("\n\nWrite the model with the Python package `rusolve` (`import rusolve as rs; m = rs.Model()`); its "
           "quickstart is {doc} (read it first). Save your code in the current directory with Write (a relative "
           "name such as model.py) and run it as exactly `{py} model.py`, with no cd, pipe or redirect. Do not "
           "work the answer out by hand: build the model, solve it, print what the solver returns, and report that.")
HIGHS_NOTE = ("\n\nSolve it with the HiGHS MCP tools you are connected to. Do not work the answer out by hand: "
              "build the model, solve it, and report what the solver returns.")


def task_modes(task):
    """Word-problem tasks name their modes; the file-based tasks run through MCP or the CLI."""
    return task.get("modes", ["mcp", "cli"])


def load_task(tid):
    with open(os.path.join(HERE, "tasks", tid, "task.json")) as f:
        t = json.load(f)
    t["prompt"] = t["prompt"].replace("{dir}", os.path.join(HERE, "tasks", tid))
    return t


def agent_command(prompt, model, mode, workdir, max_budget, needs_write=False):
    if mode == "python":
        py = os.environ.get("RSV_PY", "python3")
        return ["claude", "-p", prompt, "--model", model, "--output-format", "stream-json", "--verbose",
                "--permission-mode", "dontAsk", "--no-session-persistence", "--setting-sources", "",
                "--max-budget-usd", str(max_budget), "--strict-mcp-config",
                "--allowedTools", f"Bash({py} *)", "Write", "Edit", "Read"]
    if mode == "highs-mcp":
        server = shlex.split(os.environ.get("HIGHS_MCP_CMD", ""))
        if not server:
            raise SystemExit("--mode highs-mcp needs HIGHS_MCP_CMD, the command that starts the highs-mcp server "
                             "(npm package highs-mcp, e.g. HIGHS_MCP_CMD='npx highs-mcp'); see evals/README.md")
        cfg = os.path.join(workdir, "mcp.json")
        with open(cfg, "w") as f:
            json.dump({"mcpServers": {"highs": {"command": server[0], "args": server[1:]}}}, f)
        return ["claude", "-p", prompt, "--model", model, "--output-format", "stream-json", "--verbose",
                "--permission-mode", "dontAsk", "--no-session-persistence", "--setting-sources", "",
                "--max-budget-usd", str(max_budget), "--mcp-config", cfg, "--strict-mcp-config",
                "--allowedTools", "mcp__highs__*"]
    binary = rs.binary()
    cmd = ["claude", "-p", prompt, "--model", model, "--output-format", "stream-json", "--verbose",
           "--permission-mode", "dontAsk", "--no-session-persistence", "--setting-sources", "",
           "--max-budget-usd", str(max_budget)]
    tools = ["mcp__rusolve__*"] if mode == "mcp" else [f"Bash({binary} *)", f"Bash(RUSOLVE_MIP_HIGHS=1 {binary} *)"]
    if needs_write:  # extend_wasm: the agent must save its Wasm component to a file it then names
        tools.append("Write")
    if mode == "mcp":
        cfg = os.path.join(workdir, "mcp.json")
        with open(cfg, "w") as f:
            json.dump({"mcpServers": {"rusolve": {"command": binary, "args": ["mcp"]}}}, f)
        cmd += ["--mcp-config", cfg, "--strict-mcp-config", "--allowedTools", *tools]
    else:
        cmd += ["--strict-mcp-config", "--allowedTools", *tools]
    return cmd


def parse_stream(text):
    """Tool calls, tool errors, final text, usage and cost from claude's stream-json output."""
    calls, errors, final, result = [], 0, "", {}
    for line in text.splitlines():
        try:
            ev = json.loads(line)
        except ValueError:
            continue
        if ev.get("type") == "assistant":
            for blk in ev.get("message", {}).get("content", []):
                if blk.get("type") == "tool_use":
                    calls.append(blk.get("name"))
        elif ev.get("type") == "user":
            content = ev.get("message", {}).get("content", [])
            if isinstance(content, list):
                errors += sum(1 for b in content if b.get("type") == "tool_result" and b.get("is_error"))
        elif ev.get("type") == "result":
            result = ev
            final = ev.get("result") or ""
    usage = result.get("usage", {})
    return {
        "tool_calls": len(calls),
        "tool_names": calls,
        "tool_errors": errors,
        "turns": result.get("num_turns"),
        "input_tokens": usage.get("input_tokens", 0),
        "output_tokens": usage.get("output_tokens", 0),
        "cache_read_tokens": usage.get("cache_read_input_tokens", 0),
        "cache_creation_tokens": usage.get("cache_creation_input_tokens", 0),
        "cost_usd": result.get("total_cost_usd"),
        "agent_error": bool(result.get("is_error")) or not result,
    }, final


def run_task(tid, model, mode, max_budget, timeout):
    task = load_task(tid)
    if mode not in task_modes(task):
        raise SystemExit(f"task {tid} does not run in mode {mode} (its modes: {task_modes(task)})")
    note = {"mcp": MCP_NOTE, "python": PY_NOTE.replace("{py}", os.environ.get("RSV_PY", "python3")).replace("{doc}", os.path.join(rs.ROOT, "crates", "rusolve-py", "README.md")),
            "highs-mcp": HIGHS_NOTE}.get(mode) or CLI_NOTE.replace("{bin}", rs.binary())
    prompt = task["prompt"] + note
    with tempfile.TemporaryDirectory() as work:  # empty cwd: no project instructions leak in
        cmd = agent_command(prompt, model, mode, work, max_budget, bool(task.get("needs_write")))
        t0 = time.time()
        try:
            p = subprocess.run(cmd, capture_output=True, text=True, cwd=work, timeout=timeout,
                               env=dict(os.environ, RUSOLVE_MIP_HIGHS="1"))
            out, err = p.stdout, p.stderr
        except subprocess.TimeoutExpired as e:
            out, err = (e.stdout or b"").decode() if isinstance(e.stdout, bytes) else (e.stdout or ""), "timeout"
        wall = time.time() - t0
        rec, final = parse_stream(out)
        # Grade from the agent's working directory: an answer may name a file the agent wrote
        # there (extend_wasm names the component file it saved), and the temp dir is gone outside.
        home = os.getcwd()
        os.chdir(work)
        try:
            ok, detail = graders.grade(tid, final)
        finally:
            os.chdir(home)
    rec.update({"task": tid, "success": ok, "detail": detail, "wall_s": round(wall, 1), "final_text": final[-1500:]})
    if err.strip():
        rec["stderr_tail"] = err.strip()[-300:]
    return rec


def git_commit():
    """The checkout's commit, so results are tracked per commit."""
    try:
        return subprocess.run(["git", "-C", HERE, "describe", "--always", "--dirty", "--abbrev=8"],
                              capture_output=True, text=True, timeout=10).stdout.strip() or None
    except Exception:
        return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tasks", help="comma-separated task ids")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--model", default="haiku", help="any `claude --model` value; run a second model by repeating "
                    "the command with another (results files are named per model)")
    ap.add_argument("--mode", choices=MODES, default="mcp")
    ap.add_argument("--max-budget-usd", type=float, default=0.5, help="per task")
    ap.add_argument("--timeout", type=int, default=600, help="seconds per task")
    ap.add_argument("--out-dir", default=os.path.join(HERE, "results"))
    a = ap.parse_args()
    ids = TASKS if a.all else (a.tasks.split(",") if a.tasks else [])
    if a.mode in ("python", "highs-mcp") and not a.tasks:  # default: every word problem
        ids = ids or TASKS
    if a.all or (not a.tasks and a.mode in ("python", "highs-mcp")):  # keep the tasks that fit the mode
        ids = [t for t in ids if a.mode in task_modes(load_task(t))]
    if not ids or set(ids) - set(TASKS):
        ap.error(f"give --all or --tasks from {TASKS}")
    # every task's mode is checked before the first (paid) run, not when the loop reaches it
    wrong = [f"{t} (modes {task_modes(load_task(t))})" for t in ids if a.mode not in task_modes(load_task(t))]
    if wrong:
        ap.error(f"--mode {a.mode} does not fit: {', '.join(wrong)}")
    records = []
    for tid in ids:
        r = run_task(tid, a.model, a.mode, a.max_budget_usd, a.timeout)
        records.append(r)
        print(f"{tid:20s} {'PASS' if r['success'] else 'FAIL'}  calls={r['tool_calls']} errors={r['tool_errors']} "
              f"turns={r['turns']} out_tok={r['output_tokens']} cost=${r['cost_usd']}  {r['detail']}", flush=True)
    summary = {
        "tasks": len(records), "passed": sum(r["success"] for r in records),
        "tool_calls": sum(r["tool_calls"] for r in records),
        "output_tokens": sum(r["output_tokens"] for r in records),
        "input_tokens": sum(r["input_tokens"] + r["cache_read_tokens"] + r["cache_creation_tokens"] for r in records),
        "cost_usd": round(sum(r["cost_usd"] or 0 for r in records), 4),
    }
    solver = rs.mip(os.path.join(HERE, "tasks", "solve_report", "model.mps")).get("solver", {})
    doc = {"schema": "rusolve.evals/1", "model": a.model, "mode": a.mode, "commit": git_commit(),
           "utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
           "solver": solver, "summary": summary, "records": records}
    os.makedirs(a.out_dir, exist_ok=True)
    path = os.path.join(a.out_dir, f"{datetime.datetime.now(datetime.timezone.utc):%Y%m%dT%H%M%SZ}-{a.model}-{a.mode}.json")
    with open(path, "w") as f:
        json.dump(doc, f, indent=1)
    print(f"{summary['passed']}/{summary['tasks']} passed; {summary['tool_calls']} tool calls; "
          f"{summary['output_tokens']} output / {summary['input_tokens']} input tokens; ${summary['cost_usd']}")
    print(f"wrote {path}")
    return 0 if summary["passed"] == summary["tasks"] else 1


if __name__ == "__main__":
    sys.exit(main())
