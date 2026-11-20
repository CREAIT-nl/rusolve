"""Typed exceptions of `rusolve` (docs/design/diagnostics.md).

Every exception carries `.diagnostics`: the `diagnostics` array of the `rusolve.diagnostics/1`
document (docs/schema/rusolve-diagnostics-1.schema.json), one dict per diagnostic with `code`
(`RS4001`...), `severity`, `message`, `location`, `notes` and `help`. `rusolve --explain RS4001`
(or `rusolve.explain_code`'s source, `docs/agent/diagnostics.md`) says what each code means.
The layout follows rustc's diagnostics (https://rustc-dev-guide.rust-lang.org/diagnostics.html).

    >>> import rusolve as rs
    >>> m = rs.Model()
    >>> try:
    ...     m.add_var(lb=float("nan"), name="x")
    ... except rs.ModelError as e:
    ...     e.code
    'RS4001'

Hierarchy (the model errors are also `ValueError`/`IndexError`, as before these classes existed):

    RusolveError
    +-- ModelError (ValueError)        RS4001 invalid argument, RS4003 mixed models
    |   +-- ModelIndexError (IndexError)   RS4002
    +-- SolveError
        +-- InfeasibleError            RS4004
        +-- UnboundedError             RS4005
        +-- NoSolutionError (ValueError)   RS4007 no primal values to read
        +-- LimitError                 RS4006 stopped without an optimal solution
    +-- CallbackError                  RS4008 a hook of `Model.solve(callback=...)` raised
"""
from __future__ import annotations

import json
import linecache
import sys
import types
from typing import Any, Optional

__all__ = [
    "RusolveError",
    "ModelError",
    "ModelIndexError",
    "SolveError",
    "InfeasibleError",
    "UnboundedError",
    "NoSolutionError",
    "LimitError",
    "CallbackError",
]


class RusolveError(Exception):
    """Base of every error `rusolve` raises on purpose. `.diagnostics` is a list of dicts."""

    diagnostics: list[dict[str, Any]]

    def __init__(self, message: str = "", diagnostics: Any = ()) -> None:
        super().__init__(message)
        self.diagnostics = list(diagnostics)

    @classmethod
    def _from_json(cls, message: str, doc: str) -> "RusolveError":
        return cls(message, json.loads(doc).get("diagnostics", []))

    def __reduce__(self) -> Any:
        return (type(self), (str(self), self.diagnostics))

    @property
    def code(self) -> Optional[str]:
        """The code of the first diagnostic (`"RS4001"`), or None."""
        return self.diagnostics[0].get("code") if self.diagnostics else None

    @property
    def location(self) -> Optional[dict[str, Any]]:
        """The location of the first diagnostic: `{"kind": "python", "path", "line", "source_line"}`
        for the line of Python that made the call, or None."""
        return self.diagnostics[0].get("location") if self.diagnostics else None

    @property
    def help(self) -> list[str]:
        """The help lines (concrete fixes) of the first diagnostic."""
        return list(self.diagnostics[0].get("help", [])) if self.diagnostics else []


class ModelError(RusolveError, ValueError):
    """A model-building call got an argument it cannot use (bounds, sizes, names, mixed Models)."""


class ModelIndexError(ModelError, IndexError):
    """A row or column index is outside the model."""


class SolveError(RusolveError):
    """A solve did not end with an optimal solution; raised by `Solution.raise_for_status()`."""


class InfeasibleError(SolveError):
    """The model is infeasible. `Model.explain()` returns an irreducible infeasible subsystem;
    `.diagnostics` also carries the post-solve RS2010 (the conflicting rows)."""


class UnboundedError(SolveError):
    """The model is unbounded (or infeasible or unbounded); `.diagnostics` also carries the ray (RS2011)."""


class LimitError(SolveError):
    """The solve stopped at a limit or without a verdict (time, iterations, nodes, error)."""


class NoSolutionError(SolveError, ValueError):
    """A primal value was read from a Solution that has none (`Solution.ok` is False)."""


class CallbackError(RusolveError):
    """A hook of the callback given to `Model.solve` raised (or returned the wrong type).

    `__cause__` is the original exception; `.solution` is the partial `Solution` of the stopped
    solve (None if it could not be built). `Model.callback_stats` counts the calls and errors.
    """

    solution: Any = None


def _caller() -> Optional[tuple[str, int, str]]:
    """`(path, line, source line)` of the innermost Python frame outside `rusolve` itself."""
    try:
        f: Optional[types.FrameType] = sys._getframe(1)
        while f is not None and f.f_globals.get("__name__", "").split(".")[0] == "rusolve":
            f = f.f_back
        if f is None:
            return None
        path = f.f_code.co_filename
        line = f.f_lineno
        return path, line, linecache.getline(path, line, f.f_globals).strip()
    except Exception:  # never let the error path raise a different error
        return None
