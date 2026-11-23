"""The stubs match the compiled module, and the documented examples run and type-check.

* every public name of `rusolve._rusolve` (and each class's members) is in `_rusolve.pyi`, and
  nothing is in the stub that the module lacks;
* the `>>>` examples in the docstrings run (`doctest`);
* `mypy --strict` passes on the README's Python blocks and on the docstring examples.
"""
import ast
import doctest
import inspect
import pathlib
import re
import shutil
import subprocess
import sys

import pytest

import rusolve
from rusolve import _rusolve

PKG = pathlib.Path(rusolve.__file__).resolve().parent
README = pathlib.Path(__file__).resolve().parents[1] / "README.md"


def stub_members(path):
    tree = ast.parse(path.read_text())
    top = {}
    for node in tree.body:
        if isinstance(node, ast.ClassDef):
            members = set()
            for b in node.body:
                if isinstance(b, (ast.FunctionDef, ast.AsyncFunctionDef)):
                    members.add(b.name)
                elif isinstance(b, ast.AnnAssign) and isinstance(b.target, ast.Name):
                    members.add(b.target.id)
            top[node.name] = members
        elif isinstance(node, ast.FunctionDef):
            top[node.name] = None
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            top[node.target.id] = None
    return top


DUNDER_OK = {"__eq__", "__ge__", "__le__", "__hash__", "__init__", "__repr__", "__doc__"}


def test_stub_covers_the_module():
    stub = stub_members(PKG / "_rusolve.pyi")
    public = {n for n in dir(_rusolve) if not n.startswith("_")}
    missing = sorted(n for n in public if n not in stub)
    assert not missing, f"names of rusolve._rusolve missing from _rusolve.pyi: {missing}"
    for name in sorted(public):
        obj = getattr(_rusolve, name)
        if not inspect.isclass(obj):
            continue
        have = {m for m in dir(obj) if not m.startswith("_") or m == "_sync"}
        have |= {m for m in dir(obj) if m in {"__getitem__", "__len__", "__iter__", "__bool__", "__neg__", "__pos__"}}
        have |= {m for m in dir(obj) if m.startswith("__") and m.endswith("__") and m not in DUNDER_OK
                 and m in stub[name]}
        lacking = sorted(m for m in have if m not in stub[name] and m not in DUNDER_OK)
        assert not lacking, f"{name}: members missing from the stub: {lacking}"


def test_stub_has_nothing_the_module_lacks():
    stub = stub_members(PKG / "_rusolve.pyi")
    private_helpers = {"_Enum", "Floats", "Ints", "VarLike", "VarsLike", "Row", "Objective", "Sense", "Diagnostic"}
    for name, members in stub.items():
        if name in private_helpers:
            continue
        assert hasattr(_rusolve, name), f"_rusolve.pyi declares {name}, the module has no such name"
        if members:
            obj = getattr(_rusolve, name)
            for m in members:
                if m in DUNDER_OK or m == "__init__":
                    continue
                if name in {"HighsStatus", "HighsModelStatus", "HighsVarType", "ObjSense", "MatrixFormat"}:
                    pass
                assert hasattr(obj, m), f"_rusolve.pyi declares {name}.{m}, the module has no such member"


def test_package_stub_exports_match_all():
    stub = stub_members(PKG / "__init__.pyi")
    text = (PKG / "__init__.pyi").read_text()
    exported = set(re.findall(r"import (\w+) as \1", text)) | {"inf"}
    assert exported == set(rusolve.__all__), sorted(exported ^ set(rusolve.__all__))
    assert (PKG / "py.typed").exists()
    assert stub["inf"] is None


def test_scipy_stub_matches_the_module():
    import rusolve.scipy as s

    stub = stub_members(PKG / "scipy.pyi")
    assert set(stub) == set(s.__all__)


# ------------------------------------------------------------------ docstring examples


def docstring_sources():
    """(name, docstring) of every documented class and method of the public classes."""
    out = []
    for cname in sorted(n for n in dir(_rusolve) if not n.startswith("_")):
        cls = getattr(_rusolve, cname)
        if not inspect.isclass(cls):
            continue
        if cls.__doc__:
            out.append((cname, cls.__doc__))
        for mname, member in vars(cls).items():
            doc = getattr(member, "__doc__", None)
            if doc and ">>>" in doc:
                out.append((f"{cname}.{mname}", doc))
    return out


def test_docstring_examples_run():
    parser = doctest.DocTestParser()
    runner = doctest.DocTestRunner(optionflags=doctest.ELLIPSIS | doctest.NORMALIZE_WHITESPACE)
    n = 0
    for name, doc in docstring_sources():
        test = parser.get_doctest(doc, {}, name, None, 0)
        n += len(test.examples)
        runner.run(test)
    assert n >= 15, f"only {n} docstring examples"
    assert runner.failures == 0, f"{runner.failures} docstring examples failed (see the output above)"


def test_every_new_method_documents_an_example():
    needed = ["check", "explain", "lint", "ranging", "feasrelax", "write", "diagnose", "callback_stats"]
    docs = dict(docstring_sources())
    for m in needed:
        assert f"Model.{m}" in docs, f"Model.{m} has no runnable example"
    assert "Solution.to_json" in docs and "Solution.to_dict" in docs


# ------------------------------------------------------------------ mypy --strict

PRELUDE = "from typing import Any\n"
MYPY = [sys.executable, "-m", "mypy", "--strict", "--python-version", "3.11"]


def run_mypy(tmp_path, name, source):
    f = tmp_path / f"{name}.py"
    f.write_text(source)
    proc = subprocess.run(MYPY + [str(f)], capture_output=True, text=True, cwd=tmp_path)
    return proc.returncode, proc.stdout + proc.stderr


def readme_blocks():
    text = README.read_text()
    return re.findall(r"^```python\n(.*?)^```", text, flags=re.S | re.M)


def declare_free_names(tmp_path, name, source):
    """README fragments use data they do not define (`n`, `ub_array`, ...): declare each as Any,
    then require everything else to type-check."""
    declared = []
    for _ in range(4):
        pre = PRELUDE + "".join(f"{n}: Any = None\n" for n in declared)
        code, out = run_mypy(tmp_path, name, pre + source)
        missing = re.findall(r'Name "(\w+)" is not defined', out)
        new = [n for n in dict.fromkeys(missing) if n not in declared]
        if not new:
            return code, out
        declared += new
    return code, out


@pytest.mark.skipif(shutil.which("mypy") is None and not pathlib.Path(sys.prefix, "bin/mypy").exists(),
                    reason="mypy is not installed")
def test_mypy_strict_on_readme_examples(tmp_path):
    blocks = readme_blocks()
    assert len(blocks) >= 6
    for i, block in enumerate(blocks):
        code, out = declare_free_names(tmp_path, f"readme_{i}", block)
        assert code == 0, f"README python block {i} fails mypy --strict:\n{block}\n{out}"


@pytest.mark.skipif(shutil.which("mypy") is None and not pathlib.Path(sys.prefix, "bin/mypy").exists(),
                    reason="mypy is not installed")
def test_mypy_strict_on_docstring_examples(tmp_path):
    parser = doctest.DocTestParser()
    for i, (name, doc) in enumerate(docstring_sources()):
        src = "".join(ex.source for ex in parser.get_examples(doc))
        if not src.strip():
            continue
        # a bare expression statement is fine; its value is what the docstring shows
        code, out = declare_free_names(tmp_path, f"doc_{i}", src)
        assert code == 0, f"docstring example of {name} fails mypy --strict:\n{src}\n{out}"
