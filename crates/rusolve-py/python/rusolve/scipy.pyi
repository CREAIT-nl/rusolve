"""Type stubs of `rusolve.scipy`."""
from collections.abc import Callable, Mapping, Sequence
from typing import Any

class OptimizeResult(dict[str, Any]):
    def __getattr__(self, name: str) -> Any: ...
    def __setattr__(self, name: str, value: Any) -> None: ...

def milp(
    c: Any,
    *,
    integrality: Any = None,
    bounds: Any = None,
    constraints: Any = None,
    options: Mapping[str, Any] | None = None,
) -> OptimizeResult: ...
def linprog(
    c: Any,
    A_ub: Any = None,
    b_ub: Any = None,
    A_eq: Any = None,
    b_eq: Any = None,
    bounds: Any = None,
    method: str = "highs",
    callback: Callable[..., Any] | None = None,
    options: Mapping[str, Any] | None = None,
    x0: Any = None,
    integrality: Any = None,
) -> OptimizeResult: ...
