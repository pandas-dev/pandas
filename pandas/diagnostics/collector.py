from __future__ import annotations

import contextlib
import contextvars
from dataclasses import dataclass
import time
import tracemalloc
from typing import (
    TYPE_CHECKING,
    Any,
    TypeVar,
)

from pandas.diagnostics.report import (
    DiagnosticsReport,
    PhaseRecord,
)

if TYPE_CHECKING:
    from collections.abc import (
        Callable,
        Iterator,
    )

T = TypeVar("T")

_ACTIVE_COLLECTOR: contextvars.ContextVar[DiagnosticsCollector | None] = (
    contextvars.ContextVar(
        "pandas_diagnostics_active_collector",
        default=None,
    )
)


def collector() -> DiagnosticsCollector | None:
    """
    Internal helper for hot paths.

    Returns the active DiagnosticsCollector if diagnostics are enabled for the current
    context; otherwise returns None.
    """
    return _ACTIVE_COLLECTOR.get()


@contextlib.contextmanager
def phase(
    collector: DiagnosticsCollector | None,
    name: str,
    **meta: Any,
) -> Iterator[None]:
    """
    Cheap phase context manager: no-op if collector is None.
    """
    if collector is None:
        yield
        return
    with collector.phase(name, **meta):
        yield


@dataclass(frozen=True)
class _PhaseDone:
    name: str
    seconds: float
    mem_bytes_net: int | None = None


class DiagnosticsCollector:
    """
    Captures coarse timing phases (and optional tracemalloc deltas) for one operation.
    """

    def __init__(
        self, *, track_memory: bool = False, track_phase_memory: bool = False
    ) -> None:
        self._track_memory = track_memory
        self._track_phase_memory = track_phase_memory

        self._t0: float | None = None
        self._total_s: float | None = None

        self._phase_stack: list[tuple[str, float]] = []
        self._phases_done: list[_PhaseDone] = []

        self._counters: dict[str, int] = {}
        self._meta: dict[str, Any] = {}

        self._mem_bytes_net: int | None = None

    @contextlib.contextmanager
    def activate(self) -> Iterator[None]:
        token = _ACTIVE_COLLECTOR.set(self)
        try:
            yield
        finally:
            _ACTIVE_COLLECTOR.reset(token)

    def set_meta(self, **meta: Any) -> None:
        self._meta.update(meta)

    def inc(self, name: str, value: int = 1) -> None:
        self._counters[name] = self._counters.get(name, 0) + int(value)

    @contextlib.contextmanager
    def phase(self, name: str, **meta: Any) -> Iterator[None]:
        # meta currently unused; kept for future extensions
        snap_before = None
        if self._track_memory and self._track_phase_memory and tracemalloc.is_tracing():
            snap_before = tracemalloc.take_snapshot()

        start = time.perf_counter()
        self._phase_stack.append((name, start))
        try:
            yield
        finally:
            n, s = self._phase_stack.pop()
            end = time.perf_counter()

            mem_delta = None
            if snap_before is not None:
                snap_after = tracemalloc.take_snapshot()
                stats = snap_after.compare_to(snap_before, "filename")
                mem_delta = int(sum(st.size_diff for st in stats))

            self._phases_done.append(
                _PhaseDone(name=n, seconds=end - s, mem_bytes_net=mem_delta)
            )

    def run(self, func: Callable[..., T], *args: Any, **kwargs: Any) -> T:
        self._t0 = time.perf_counter()

        mem_was_on = tracemalloc.is_tracing()
        snap_before = None

        if self._track_memory:
            if not mem_was_on:
                tracemalloc.start()
            snap_before = tracemalloc.take_snapshot()

        try:
            return func(*args, **kwargs)
        finally:
            t0 = self._t0 if self._t0 is not None else time.perf_counter()
            self._total_s = time.perf_counter() - t0

            if self._track_memory and snap_before is not None:
                snap_after = tracemalloc.take_snapshot()
                stats = snap_after.compare_to(snap_before, "filename")
                self._mem_bytes_net = int(sum(st.size_diff for st in stats))

                if not mem_was_on:
                    tracemalloc.stop()

    def to_report(self) -> DiagnosticsReport:
        return DiagnosticsReport(
            total_seconds=float(self._total_s or 0.0),
            phases=[
                PhaseRecord(p.name, p.seconds, mem_bytes_net=p.mem_bytes_net)
                for p in self._phases_done
            ],
            counters=dict(self._counters),
            metadata=dict(self._meta),
            memory_bytes_net=self._mem_bytes_net,
        )
