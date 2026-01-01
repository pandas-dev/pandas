from __future__ import annotations

from typing import (
    TYPE_CHECKING,
    Any,
    TypeVar,
)

if TYPE_CHECKING:
    from collections.abc import Callable

from pandas.diagnostics.collector import DiagnosticsCollector
from pandas.diagnostics.report import (
    DiagnosticsReport,
    attach_report,
    explain_attached,
)

T = TypeVar("T")


def report(
    func: Callable[..., T], *args: Any, memory: bool = False, **kwargs: Any
) -> T:
    """
    Run a callable with pandas diagnostics enabled for the duration of the call.
    The resulting DiagnosticsReport is attached to the returned object (when possible).
    Retrieve it via obj.explain().
    """
    collector = DiagnosticsCollector(track_memory=memory)

    with collector.activate():
        result = collector.run(func, *args, **kwargs)

    attach_report(result, collector.to_report())
    return result


def explain(obj, *, format: str = "text"):
    """
    Format the diagnostics report attached to an object
    returned by pd.diagnostics.report(...).
    """
    return explain_attached(obj, format=format)


__all__ = ["DiagnosticsReport", "explain", "report"]
