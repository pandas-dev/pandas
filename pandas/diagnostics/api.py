from __future__ import annotations

from typing import (
    TYPE_CHECKING,
    Any,
    TypeVar,
)

if TYPE_CHECKING:
    from collections.abc import Callable

from pandas.diagnostics.collector import DiagnosticsCollector
from pandas.diagnostics.report import attach_report

T = TypeVar("T")


def report(
    func: Callable[..., T], *args: Any, memory: bool = False, **kwargs: Any
) -> T:
    collector = DiagnosticsCollector(track_memory=memory)
    with collector.activate():
        result = collector.run(func, *args, **kwargs)
    attach_report(result, collector.to_report())
    return result
