"""Ground truth from rusolve itself: the CLI for one-shot runs, a stdio MCP client for edits.

Graders use this so the expected answer is whatever the solver under test computes, never a
number typed into a task file. Standard library only.
"""
import json
import os
import subprocess

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def binary():
    """The rusolve-cli binary: $RSV_CLI, else target/release/rusolve-cli of this checkout."""
    p = os.environ.get("RSV_CLI") or os.path.join(ROOT, "target", "release", "rusolve-cli")
    if not os.path.exists(p):
        raise SystemExit(f"rusolve-cli not found at {p}: cargo build --release -p rusolve-cli or set RSV_CLI")
    return p


def cli(args, extra_env=None):
    env = dict(os.environ, RUSOLVE_MIP_HIGHS="1", **(extra_env or {}))
    return subprocess.run([binary(), *map(str, args)], capture_output=True, text=True, env=env, timeout=300)


def cli_json(args):
    r = cli(args)
    if r.returncode != 0:
        raise RuntimeError(f"rusolve-cli {args} exited {r.returncode}: {r.stderr.strip()[:400]}")
    return json.loads(r.stdout)


def mip(model, options=None, time_limit=None, wasm_hook=None):
    """`rusolve mip --format json` envelope (the ported HiGHS MIP path)."""
    args = ["mip", model, "--format", "json", "--solution", "--seed", "0"]
    if time_limit is not None:
        args += ["--time-limit", time_limit]
    if wasm_hook is not None:
        args += ["--wasm-hook", wasm_hook]
    for k, v in (options or {}).items():
        args += ["--opt", f"{k}={str(v).lower() if isinstance(v, bool) else v}"]
    return cli_json(args)


class Mcp:
    """Minimal newline-delimited JSON-RPC client for `rusolve-cli mcp` (docs/design/mcp.md)."""

    def __init__(self):
        self.p = subprocess.Popen([binary(), "mcp"], stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                  stderr=subprocess.DEVNULL, text=True)
        self.n = 0
        self._rpc("initialize", {"protocolVersion": "2025-11-25", "capabilities": {},
                                 "clientInfo": {"name": "rusolve-evals", "version": "1"}})
        self._send({"jsonrpc": "2.0", "method": "notifications/initialized"})

    def _send(self, msg):
        self.p.stdin.write(json.dumps(msg) + "\n")
        self.p.stdin.flush()

    def _rpc(self, method, params):
        self.n += 1
        self._send({"jsonrpc": "2.0", "id": self.n, "method": method, "params": params})
        while True:
            line = self.p.stdout.readline()
            if not line:
                raise RuntimeError("rusolve mcp closed stdout")
            msg = json.loads(line)
            if msg.get("id") == self.n:
                return msg

    def call(self, tool, **args):
        """The tool's structuredContent; raises RuntimeError on a tool error."""
        msg = self._rpc("tools/call", {"name": tool, "arguments": args})
        if "error" in msg:
            raise RuntimeError(f"{tool}: {msg['error']}")
        res = msg["result"]
        if res.get("isError"):
            raise RuntimeError(f"{tool}: {res['content'][0]['text']}")
        return res.get("structuredContent")

    def close(self):
        try:
            self.p.stdin.close()
            self.p.wait(timeout=5)
        except Exception:
            self.p.kill()

    def __enter__(self):
        return self

    def __exit__(self, *a):
        self.close()
