#!/usr/bin/env python3
"""Wheel smoke test: install the wheel in a clean venv; `rusolve --help` and an MCP initialize
over `python -m rusolve mcp` must answer. Usage: python tools/smoke_wheel.py <wheel-dir>"""
import glob
import json
import os
import subprocess
import sys
import tempfile
import venv

wheels = glob.glob(os.path.join(sys.argv[1], "*.whl"))
assert len(wheels) == 1, f"expected one wheel in {sys.argv[1]}, found {wheels}"
with tempfile.TemporaryDirectory() as tmp:
    venv.create(tmp, with_pip=True)
    bindir = os.path.join(tmp, "Scripts" if os.name == "nt" else "bin")
    py = os.path.join(bindir, "python")
    subprocess.check_call([py, "-m", "pip", "install", "-q", wheels[0]])
    out = subprocess.run([os.path.join(bindir, "rusolve"), "--help"], capture_output=True, text=True)
    assert "usage: rusolve" in out.stdout + out.stderr, out
    init = {"jsonrpc": "2.0", "id": 1, "method": "initialize", "params": {
        "protocolVersion": "2025-06-18", "capabilities": {}, "clientInfo": {"name": "smoke", "version": "0"}}}
    r = subprocess.run([py, "-m", "rusolve", "mcp"], input=json.dumps(init) + "\n",
                       capture_output=True, text=True, timeout=60)
    reply = json.loads(r.stdout.splitlines()[0])
    assert reply["result"]["serverInfo"]["name"] == "rusolve", reply
    subprocess.check_call([py, "-c", "import rusolve; rusolve.Model()"])
    print("smoke ok:", os.path.basename(wheels[0]))
