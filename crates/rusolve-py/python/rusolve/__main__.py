"""`rusolve` console script and `python -m rusolve`: run the bundled CLI.

The wheel carries the `rusolve-cli` binary in `rusolve/_bin/` (placed there by
`tools/build_wheel.sh` before `maturin build`; maturin packs everything under the
python-source package directory, https://www.maturin.rs/project_layout.html). This shim hands
the arguments to it, so `rusolve mcp`, `python -m rusolve mcp` and `uvx --from <wheel> rusolve mcp`
all start the same stdio MCP server (docs/design/mcp.md).
"""
import os
import shutil
import stat
import subprocess
import sys
from pathlib import Path

_NAME = "rusolve-cli" + (".exe" if os.name == "nt" else "")


def find_binary() -> str:
    bundled = Path(__file__).resolve().parent / "_bin" / _NAME
    if bundled.is_file():
        if os.name != "nt" and not os.access(bundled, os.X_OK):
            try:  # some installers drop the mode bit
                bundled.chmod(bundled.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
            except OSError:
                pass
        return str(bundled)
    found = shutil.which("rusolve-cli")
    if found:
        return found
    sys.stderr.write(
        "rusolve: the bundled rusolve-cli binary is missing from this install "
        "(a source checkout without tools/build_wheel.sh?) and none is on PATH.\n"
    )
    raise SystemExit(4)


def main() -> None:
    exe = find_binary()
    args = [exe, *sys.argv[1:]]
    if os.name == "nt":  # no real exec on Windows
        raise SystemExit(subprocess.call(args))
    sys.stdout.flush()
    os.execv(exe, args)


if __name__ == "__main__":
    main()
